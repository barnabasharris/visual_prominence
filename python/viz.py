#!/usr/bin/env python3
import os
import grass.script as gscript
import random
import sys
import datetime

def my_print(text):
    currentDT = datetime.datetime.now()
    sys.stdout.write('\n')
    sys.stdout.write(str(currentDT))
    sys.stdout.write('\n' + str(text))
    sys.stdout.flush()


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
    
    ## Create a list with point names from viewpoints
    pointlist1 = gscript.read_command('v.db.select',
                                      map='points',
                                      columns='cat',
                                      flags='c',
                                      overwrite = True)
    
    pointlist2 = pointlist1.split()
    pointlist=list(set(pointlist2)) # remove duplicates
    
    # set up empty visual prominence map
    my_print('creating empty visual prominence map...')
    gscript.raster.mapcalc('visual_prominence =' + dem + ' * 0',
                           overwrite=True)
    
    # DEBUG: take a random sample of points
    my_print('sampling 10 points from data...')
    smp = random.sample(range(0,len(pointlist)),10)                 
    pointlist = [pointlist[i] for i in smp]
    
    ## Reiterate through viewpoints
    my_print('iterating through points...')
    for point in pointlist:
        my_print('point number ' + point + '...')
        gscript.run_command('v.extract',
                            input='points',
                            where='cat=' + """'""" + point + """'""",
                            output='point',
                            overwrite = True)
        
        gscript.run_command('r.viewshed.cva',
                            input=dem,
                            vector='point',
                            output='point_view',
                            max_distance=viewdist,
                            observer_elevation=viewobserver,
                            target_elevation=viewtarget,
                            memory='16000',
                            overwrite = True)
        
        gscript.raster.mapcalc('visual_prominence = visual_prominence + point_view',
                               overwrite=True)
        
        gscript.run_command('g.remove',
                            flags='f',
                            type='raster,vector',
                            name='point_view,point')
    
    
    # export prominence map
    my_print('exporting prominence map...')
    gscript.run_command('r.out.gdal',
                        input='visual_prominence',
                        output = output_loc + '/visual_prominence_' + output_it + '.tif',
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

