# !diagnostics off
library(terra)
library(here)
library(purrr)
# library(Rmpi)
library(snow)
library(vapour)
library(glue)

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
gridRes <- 10000
demFile <- 'bigdata/os_50m_masked.tif'

# load terrain raster
r <- terra::rast(demFile)
# make polygon outline of non-na cells
r.bin <- r
r.bin[!is.nan(r.bin)] <- 1
r.bin.pol <- as.polygons(r.bin)
# how many vs to compute?
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
}

# calculate threads / batches 
if (env == 'LOCAL') {
  length(r.tiles)
  plot(rast(r.tiles[[122]]))
  r.p <- as.points(rast(r.tiles[[122]]))
  secs <- nrow(r.p) * 3 # avg time per viewshed in seconds
  mins <- secs / 60
  hours <- mins / 60
  length(r.tiles)/4
  # 4 batches of 230 threads running for 34 hours 
}

# work out speed
if (env == 'LOCAL') {
  r <- rast(r.tiles[[50]])
  plot(r)
  200*200
  p.geom <- terra::geom(as.points(rast(r.tiles[[24]])))
  
  p.geom.coords <- 1:nrow(p.geom) %>% 
    map_chr(.f = function(x) {
      paste0(p.geom[x,c('x','y')],collapse='|')
    })
  
  which(p.geom.coords == '322075|998025') / nrow(p.geom)
  st <- lubridate::as_datetime('2022-11-09 14:18:51.726620')
  now <- lubridate::as_datetime('2022-11-10 10:12:00.851957')
  dur <- difftime(now, st, units='secs') / which(p.geom.coords == '322075|998025')
}


##•┣ prepare grass env ----
grassloc <- glue('{sysTmpDir}/grassdb')
grassPermanent <- file.path(grassloc,'PERMANENT')

# create location
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

# map storage size in megabytes
mapMB <- mean(c(7.8,7.6,6.8,5.0,8.8,5.3)) / 100
# memory required for all max points per tiles (40000)
40000 * mapMB
# max scratch memory is 25000
# max num of cores if all max points
floor(250000 / (40000 * mapMB))

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
  datOut <- snow::clusterApply(cl, 1:50, viewpointAnalysis, cva=T)
  # datOut <- snow::clusterApplyLB(cl, 1:50, viewpointAnalysis)
}

# remove grass location
unlink(grassloc,recursive = T)
print(Sys.time())
print('done!')
# Clean up the cluster and release the relevant resources.
stopCluster(cl)
mpi.quit()




