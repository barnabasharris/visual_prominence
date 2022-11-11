library(terra)
library(here)
library(purrr)
library(snow)
library(vapour)
library(glue)
library(stringr)

if (str_detect(here(),'lustre')) {
  env <- 'KATHLEEN'
} else env <- 'LOCAL'

wd <- here()

if (dir.exists('tmp')) {
  unlink('tmp',recursive=T)
}

dir.create('tmp')

if (dir.exists('outputs')) {
  unlink('outputs',recursive=T)
}

dir.create('outputs')

if (dir.exists('logs')) {
  unlink('logs',recursive=T)
}

dir.create('logs')

# set system tmp
Sys.setenv(TMPDIR=file.path(here(),'tmp'))

print('system tempdir is...')
sysTmpDir <- Sys.getenv("TMPDIR")
print(sysTmpDir)

print('main R tempdir is...')
print(tempdir())

# set permissions for python script to allow executable
system(glue('chmod -R +x {wd}/python/'))

# pre-process date for analysis -----
gridRes <- 6500 # decided by maximum job length allowed by kathleen, see  below
demFile <- 'bigdata/os_50m_masked.tif'

# load terrain raster
r <- rast(demFile)
# make polygon outline of non-na cells
r.bin <- r
r.bin[!is.nan(r.bin)] <- 1
r.bin.pol <- as.polygons(r.bin)

# optional: how many vs to compute?
r.bin[is.nan(r.bin)] <- 0
global(r.bin, 'sum')

# copy original r to template
template <- terra::rast(r)
# set new res of template
res(template) <- gridRes
# add arbitrary values
values(template) <- 1:ncell(template)

# convert to polygons
template.pols <- as.polygons(template)
# clip grid by outline
template.pols.m <- terra::mask(template.pols,r.bin.pol)

# visualize
if (env == 'LOCAL') {
  plot(r)
  plot(template.pols.m, add=T)
}

# x <- 100
# loop through each polygon, load intersecting part of overall raster
r.tiles <- 
  1:length(template.pols.m) %>%
  map(.f = function(x) {
    print(x)
    dm <- gridRes / res(r)
    crsObj <- crs(r)
    p <- template.pols[x]
    exSpat <- ext(p)
    exNum <- c(exSpat$xmin,exSpat$xmax,exSpat$ymin,exSpat$ymax)
    names(exNum) <- NULL

    vals <- vapour::vapour_warp_raster(demFile, extent = exNum, dimension = dm, projection = crsObj)
    # can't set CRS here, as for some reason Kathleen throws an error when unwrapping the rast
    r <- setValues(rast(extent=exSpat, nrows = dm[2], ncols = dm[1]), vals[[1]])
    
    if (all(is.nan(r[]))) {
      return(NULL)
    } else {
      return(wrap(r))
    }
  }) %>% 
  compact()


if (env == 'LOCAL') {
  r.tiles.length <- r.tiles %>% 
    map_int(.f = function(x) {
      length(as.points(rast(x)))
    })
  
  # assumptions for estimates
  kathleenSpeed <- 9.43 # seconds per viewshed
  agisSpeed <- 0.5 # seconds per viewshed
  kathleenCores <- 250
  agisCores <- 7
  
  # split into reasonably sized chunks based on length
  tileDf <- data.frame(tileID = 1:length(r.tiles.length),
                       tilePoints = r.tiles.length,
                       tilePointsGrp = cut(r.tiles.length,5,labels=LETTERS[5:1])
  )
  
  m <- tileDf %>% 
    dplyr::filter(tilePoints == 16900)
  dput(m$tileID[1:160])
  
  tileDf.sum <- tileDf %>% 
    dplyr::group_by(tilePointsGrp) %>% 
    dplyr::summarise(maxPoints = max(tilePoints),
                     grpNum = length(tilePoints)) %>% 
    dplyr::mutate(JobRunTimeKathleen = ((maxPoints * kathleenSpeed) / 60) / 60,
                  numJobsKathleen = ceiling(grpNum / kathleenCores),
                  totalRunTimeKathleen = numJobsKathleen * JobRunTimeKathleen,
                  JobRunTimeAGIS = ((maxPoints * agisSpeed) / 60) / 60,
                  numJobsAGIS = ceiling(grpNum / agisCores),
                  totalRunTimeAGIS = numJobsAGIS * JobRunTimeAGIS
    )
  
}


##•┣ prepare grass env ----
grassloc <- glue('{sysTmpDir}/grassdb')
grassPermanent <- file.path(grassloc,'PERMANENT')

# create location if required
if (dir.exists(grassloc)) unlink(grassloc,recursive = T)
system(glue('grass -c EPSG:27700 {grassloc} -e'))

sharedDemName <- 'dem'
# add raster to PERMANENT mapset for shared access across nodes
system(glue('grass {grassPermanent} --exec r.in.gdal {file.path(here(),demFile)} output={sharedDemName}'))

# multi-thread analysis ----

##•┣ define function ----
# x <- 1
# x <- 221
viewpointAnalysis <- function(x, cva=F) {
  st <- Sys.time()
  library(terra)
  library(glue)
  setwd(wd)
  # extract tile
  if (cva) {
    sink(glue('logs/analysis_cva_sinkout_{x}.txt'))
  } else sink(glue('logs/analysis_sinkout_{x}.txt'))
  
  print('extracting raster tile...')
  r <- rast(r.tiles[[x]])
  crs(r) <- 'EPSG:27700'
  
  # make into points
  print('converting to points...')
  
  r.points <- as.points(r)
  
  # write points to disk, ready for grass import
  print('writing points to disk...')
  pointsLoc <- glue('{sysTmpDir}/viewpoints_{x}.gpkg')
  if (file.exists(pointsLoc)) file.remove(pointsLoc)
  writeVector(r.points,pointsLoc,overwrite=T)
  
  # create grass mapset for node
  print('creating mapset...')
  grassMapset <- glue('{grassloc}/mapset_{x}')
  if (dir.exists(grassMapset)) unlink(grassMapset,recursive = T)
  system(glue('grass -c {grassloc}/mapset_{x} -e'))
  
  dem <- sharedDemName
  viewdist <- 10000
  viewobserver <-  1.75
  viewtarget <- 1.75
  pointsloc <- pointsLoc
  output_loc <- file.path(wd,'outputs')
  output_it <- x
  
  print('node preparation complete')
  diff <- Sys.time() -  st
  print(diff)
  # execute python script
  st <- Sys.time()
  # pick function
  if (cva) {
    pyLoc = file.path(wd,'python','viz_viewshed_cva_debug.py')
    logStem = 'viz_out_cva'
  } else {
    pyLoc = file.path(wd,'python','viz_viewshed_debug.py')
    logStem = 'viz_out'
  }
  
  print(paste0('running python script at...',pyLoc))
  system2('grass',
          paste(grassMapset,
                '--exec',
                pyLoc,
                dem,
                viewdist,
                viewobserver,
                viewtarget,
                pointsloc,
                output_loc,
                output_it
                ),
          stderr = paste0(getwd(),'/logs/',logStem,'_e',output_it,'.txt'),
          stdout = paste0(getwd(),'/logs/',logStem,'_o',output_it,'.txt')
          )
  print('viewshed analysis complete')
  diff <- Sys.time() -  st
  print(diff)
  # remove the grass mapset
  # unlink(grassMapset,recursive = T)
  sink()
  return(print(glue('node with job {x} finished')))
}

varsToExport <- c('r.tiles','grassloc','sysTmpDir','wd','sharedDemName')

##•┣ {parallel} version for testing ----
if (env == 'LOCAL') {
  cl <- parallel::makePSOCKcluster(6)
  wd <- here::here()
  parallel::clusterExport(cl, varlist = varsToExport)
  datOut <- parallel::clusterApply(cl, c(4,7,8,11,12,13), viewpointAnalysis, cva=T)
}

##•┣ {snow} version for kathleen ----
if (env == 'KATHLEEN') {
  #!! if offSet = T the use pd$tiles$pol !!
  print('getting MPI cluster...')
  cl <- snow::getMPIcluster()
  print('done!')
  
  # Display info about each process in the cluster
  print(clusterCall(cl, function() Sys.info()))
  print('exporting vars to nodes...')
  snow::clusterExport(cl, varsToExport)
  print(Sys.time())
  print('running analysis...')
  datOut <- snow::clusterApply(cl, c(17L, 25L, 84L, 96L, 97L, 99L, 108L, 115L, 116L, 117L, 118L, 
                                     119L, 120L, 121L, 122L, 123L, 128L, 133L, 136L, 138L, 139L, 140L, 
                                     141L, 142L, 143L, 144L, 145L, 146L, 147L, 148L, 155L, 156L, 157L, 
                                     163L, 164L, 165L, 166L, 167L, 168L, 169L, 170L, 171L, 172L, 173L, 
                                     174L, 175L, 176L, 177L, 178L, 184L, 185L, 186L, 192L, 193L, 194L, 
                                     195L, 196L, 197L, 198L, 199L, 200L, 201L, 202L, 203L, 204L, 205L, 
                                     206L, 207L, 214L, 215L, 226L, 227L, 228L, 229L, 230L, 231L, 232L, 
                                     233L, 234L, 235L, 236L, 237L, 242L, 244L, 245L, 251L, 252L, 255L, 
                                     256L, 257L, 258L, 259L, 260L, 261L, 262L, 263L, 264L, 265L, 266L, 
                                     271L, 279L, 280L, 281L, 282L, 283L, 284L, 285L, 286L, 287L, 288L, 
                                     289L, 290L, 291L, 298L, 305L, 306L, 307L, 308L, 309L, 310L, 311L, 
                                     312L, 313L, 314L, 315L, 316L, 317L, 330L, 331L, 332L, 333L, 334L, 
                                     335L, 336L, 337L, 338L, 339L, 340L, 341L, 358L, 359L, 360L, 361L, 
                                     362L, 363L, 364L, 365L, 366L, 386L, 387L, 388L, 389L, 390L, 393L, 
                                     409L, 410L, 411L, 412L, 413L, 414L),
                               viewpointAnalysis, cva=F)
  # datOut <- snow::clusterApplyLB(cl, 1:50, viewpointAnalysis)
}

# remove grass location
unlink(grassloc,recursive = T)
print(Sys.time())
print('done!')
# Clean up the cluster and release the relevant resources.
stopCluster(cl)
mpi.quit()




