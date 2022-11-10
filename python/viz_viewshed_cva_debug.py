#!/usr/bin/env python3
import os
import grass.script as gscript
import random
import sys
import datetime
import itertools

def my_print(text):
    currentDT = datetime.datetime.now()
    sys.stdout.write('\n')
    sys.stdout.write(str(currentDT))
    sys.stdout.write('\n' + str(text))
    sys.stdout.flush()

# dummy vars for debugging
dem = 'dem'
viewdist = 10000
viewobserver =  1.75
viewtarget = 1.75
pointsloc = '/tmp/viewpoints_221.gpkg'
output_loc = "/media/mal/git/visual_prominence/outputs"
output_it = '1'

def calcViz(dem, viewdist, viewobserver, viewtarget, pointsloc, output_loc, output_it):
    # set region and mask
    gscript.run_command('g.region',raster=dem)
    gscript.run_command('r.mask', raster=dem)
    
    # import points
    my_print('importing data...')
    gscript.run_command('v.in.ogr',
                    input=pointsloc,
                    output='points',
                    overwrite=True)
    
    # set region to just points
    gscript.run_command('g.region',vector='points')
    
    # calculate buffer in cells required around points to accommadate viewdistance
    nsres = float(gscript.read_command('g.region', flags='pg').split('\n')[6].split('nsres=')[1])
    growCells = str(float(viewdist) / nsres)
    # apply buffer
    gscript.run_command('g.region',grow=growCells)
    
    # DEBUG: take a random sample of points
    ## Create a list with point coords
    pointinfo = gscript.read_command('v.info',
                                    map='points',
                                    flags='t')
    pointnum = pointinfo.split('\n')[1].split('points=')[1]
    pointnum = int(pointnum)
    
    smpnum = 1000
    my_print('sampling ' + str(smpnum) +' points from data...')
    smp = random.sample(range(0,pointnum),smpnum)
    smp_strings = [str(x) for x in smp]
    sqlquery = 'cat IN (' + ','.join(smp_strings) +')'
    
    my_print('extracting ' + str(smpnum) + ' points from data...')
    gscript.run_command('v.extract',
                    input='points',
                    where=sqlquery,
                    output='points_smp',
                    overwrite = True)
    
    my_print('running viewshed.cva...')
    gscript.run_command('r.viewshed.cva',
                        input=dem,
                        vector='points_smp',
                        output='visual_prominence',
                        max_distance=viewdist,
                        observer_elevation=viewobserver,
                        target_elevation=viewtarget,
                        memory='16000',
                        overwrite = True)
    
    my_print('exporting map...')
    gscript.run_command('r.out.gdal',
                    input='visual_prominence',
                    output = output_loc + '/visual_prominence_cva_' + str(output_it) + '.tif',
                    overwrite = True)
    my_print('complete...')


# parse system arguments
dem = sys.argv[1] # 'dem'
viewdist = sys.argv[2] # field of view distance 10000
viewobserver = sys.argv[3] # height of the viewer 1.75
viewtarget = sys.argv[4] # height of the target 1.75
pointsloc = sys.argv[5] # '/tmp/viewpoints_1.gpkg'
output_loc = sys.argv[6] # '/media/mal/git/visual_prominence/outputs'
output_it = sys.argv[7] # '1'

# run function
calcViz(dem, viewdist, viewobserver, viewtarget, pointsloc, output_loc, output_it)

