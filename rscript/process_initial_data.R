# !diagnostics off
library(terra)
library(here)
library(parallel)
library(furrr)
library(purrr)

hgts <- list.files('/media/mal/srtm',pattern='*hgt$',
                   full.names = T)

plan(multisession, workers=6)

hgts %>% 
  future_walk(.f = function(x) {
    print(x)
    srtm <- rast(x)
    srtm.osgb <- terra::project(srtm, st_crs(27700,parameters=T)$wkt)
    srtm.osgb.50m <- terra::aggregate(srtm.osgb, 2)
    srtm.osgb.50m[srtm.osgb.50m == 0] <- NA
    srtm.osgb.50m.trimmed <- terra::trim(srtm.osgb.50m)
    srtm.osgb.50m.trimmed.exact <- srtm.osgb.50m.trimmed
    terra::res(srtm.osgb.50m.trimmed.exact) <- 50
    srtm.osgb.final <- terra::resample(srtm.osgb.50m.trimmed, srtm.osgb.50m.trimmed.exact)
    terra::writeRaster(srtm.osgb.final,
                       paste0('/media/mal/srtm_osgb/',str_replace(basename(x),'hgt','tif'))
    )
  })

tifs <- list.files('/media/mal/srtm_osgb/',pattern='*tif$',
                   full.names = T)


library(future)
library(furrr)
plan(multisession, workers=6)
1500*1024^2

options(future.globals.maxSize = 1572864000)

# 1:6 %>% 
1:length(r.tiles) %>%
  future_walk(.f = function(x) {
    grassloc <- '/tmp/grassdb'
    wd <- '/media/mal/git/visual_prominence'
    
    r <- rast(r.tiles[[x]])
    r.tmp <- r
    r.tmp[!is.nan(r.tmp)] <- 1
    r.poly <- as.polygons(r.tmp)
    r.poly.filled <- fillHoles(r.poly)
    
    r.tmp <- rast(r)
    r.filled <- rasterize(r.poly.filled, r.tmp)
    
    writeRaster(r, glue('outputs/rast_{x}.tif'))
    writeRaster(r.filled, glue('outputs/rast_mask_{x}.tif'))
    
    grassMapset <- glue('{grassloc}/mapset_{x}')
    if (dir.exists(grassMapset)) unlink(grassMapset,recursive = T)
    system(glue('grass -c {grassloc}/mapset_{x} -e'))
    # interpolate
    system(glue('grass {grassloc}/mapset_{x} --exec r.in.gdal input={wd}/outputs/rast_{x}.tif output=dem_{x}'))
    system(glue('grass {grassloc}/mapset_{x} --exec r.in.gdal input={wd}/outputs/rast_mask_{x}.tif output=dem_mask_{x}'))
    system(glue('grass {grassloc}/mapset_{x} --exec g.region rast=dem_mask_{x}'))
    system(glue('grass {grassloc}/mapset_{x} --exec r.mask rast=dem_mask_{x}'))
    system(glue('grass {grassloc}/mapset_{x} --exec r.fillnulls input=dem_{x} output=dem_int_{x}'))
    system(glue('grass {grassloc}/mapset_{x} --exec r.out.gdal input=dem_int_{x} output={wd}/outputs/rast_int_{x}.tif'))
    
  })



