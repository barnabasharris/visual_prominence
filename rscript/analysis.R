# !diagnostics off
library(terra)
library(here)
library(purrr)
# library(Rmpi)
library(snow)
library(vapour)
library(glue)

if (here() == "/media/mal/git/visual_prominence") {
  env <- 'LOCAL'
}

if (here() == "/lustre/scratch/tcrnbgh/visual_prominence") {
  env <- 'KATHLEEN'
}

wd <- here()

if (env == 'KATHLEEN') Sys.setenv(TMPDIR='/home/tcrnbgh/Scratch/tmp')
if (env == 'LOCAL') Sys.setenv(TMPDIR='/tmp')

print('system tempdir is...')
sysTmpDir <- Sys.getenv("TMPDIR")
print(sysTmpDir)

print('main R tempdir is...')
print(tempdir())

if (dir.exists('logs')) {
  unlink('logs',recursive=T)
}

if (dir.exists('outputs')) {
  unlink('outputs',recursive=T)
}

dir.create('logs')
dir.create('outputs')

system(glue('chmod -R +x {wd}/python/'))

# pre-process date for analysis -----
gridRes <- 10000
demFile <- 'bigdata/os_50m_masked.tif'

# load terrain raster
r <- terra::rast(demFile)
# make polygon outline of non-na cells
r.bin <- r
r.bin[!is.nan(r.bin)] <- 1
r.bin.pol <- as.polygons(r.bin)
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
  # 1:500 %>%
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

##?????? prepare grass env ----
grassloc <- glue('{sysTmpDir}/grassdb')
grassPermanent <- file.path(grassloc,'PERMANENT')

# create location
if (dir.exists(grassloc)) unlink(grassloc,recursive = T)
system(glue('grass -c EPSG:27700 {grassloc} -e'))

sharedDemName <- 'dem'
# add raster to PERMANENT mapset for shared access across nodes
system(glue('grass {grassPermanent} --exec r.in.gdal {file.path(here(),demFile)} output={sharedDemName}'))

# multi-thread analysis ----

##?????? define function ----
# x <- 1
# x <- 221
viewpointAnalysis <- function(x) {
  st <- Sys.time()
  library(terra)
  library(glue)
  setwd(wd)
  # extract tile
  
  sink(glue('logs/analysis_sinkout_{x}.txt'))
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
  print('running python script...')
  system2('grass',
          paste(grassMapset,
                '--exec',
                file.path(wd,'python','viz_viewshed.py'),
                dem,
                viewdist,
                viewobserver,
                viewtarget,
                pointsloc,
                output_loc,
                output_it
                ),
          stderr = paste0(getwd(),'/logs/viz_out_e',output_it,'.txt'),
          stdout = paste0(getwd(),'/logs/viz_out_o',output_it,'.txt')
          )
  print('viewshed analysis complete')
  diff <- Sys.time() -  st
  print(diff)
  # remove the grass mapset
  unlink(grassMapset,recursive = T)
  sink()
  return(print(glue('node with job {x} finished')))
}

varsToExport <- c('r.tiles','grassloc','sysTmpDir','wd','sharedDemName')

##?????? {parallel} version for testing ----
if (env == 'LOCAL') {
  cl <- parallel::makePSOCKcluster(6)
  wd <- here::here()
  parallel::clusterExport(cl, varlist = varsToExport)
  datOut <- parallel::clusterApply(cl, 1:6, viewpointAnalysis)
}


##?????? {snow} version for kathleen ----
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
  datOut <- snow::clusterApply(cl, 1:230, viewpointAnalysis)
}

# remove grass location
unlink(grassloc,recursive = T)
print(Sys.time())
print('done!')
# Clean up the cluster and release the relevant resources.
stopCluster(cl)
mpi.quit()



