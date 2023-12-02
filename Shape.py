"""
made to create shape files of all counties in the state and store here
"""

import glob
import os
from osgeo import gdal
from osgeo.gdalconst import GA_ReadOnly
import numpy as np
import rasterio
import rioxarray as rxr
from shapely.geometry import Polygon
import geopandas as gpd
from fiona.crs import from_epsg
import time
study_area="D:\\Study_Areas\\Lonoke\\20190906"
#all_input_path = glob.glob('/media/16TB-HDD1/old_tiles backup/MS_055Issaquena/')
all_input_path = glob.glob(study_area+"*\\*m\\*")
print(all_input_path)

def subdivide(in_path, out_path):
    all_input_path = glob.glob(in_path+'\\*.tif')
    input_filename = all_input_path[0].split('\\')[-1]
    output_filename = str(input_filename.split('.')[0]+'_')
    temp = "D:\\Study_Areas\\temp\\"
 #   com_string1 = "gdalwarp -overwrite -t_srs EPSG:26915 -r near -of GTiff " + str(in_path) + str(
 #       input_filename) + " " + str(temp) + str(input_filename)
 #   os.system(com_string1)
    print(output_filename)
    print(all_input_path)
    split_input_path = in_path.split('\\')
    # Desired num pixels for each subdivided TIF in x direction
    tile_size_x = 5000
    # Desired num pixels for each subdivided TIF in y direction
    tile_size_y = 5000

    # Open the input TIF
    ds = gdal.Open('{}'.format(all_input_path[0]))
    # Get the first band (used for calculating size,
    # so this assumes all bands have same size)
    band = ds.GetRasterBand(1)
    # Set size of TIF in x direction
    xsize = band.XSize
    # Set size of TIF in y direction
    ysize = band.YSize

    # Tile index, used to name the tiles
    tileIndex = 0
    if not os.path.exists(out_path):
        os.mkdir(out_path)
        print(f'{out_path} created')

    # i & j here correspond to the x-value and y-value, respectively, of the
    # top-left pixel index for each of the newly created tiles (e.g. the first
    # tile will have (i, j) = (0,0), the second will have
    # (i, j) = (0, tile_size_x), etc.)

    if False:
        for i in range(0, xsize, tile_size_x):
            for j in range(0, ysize, tile_size_y):
                # Index the tileIndex used for labeling; this is done before the first
                # tile is created so that the first index label is 001 rather than 000
                tileIndex = tileIndex + 1
                # Convert the tile index to a string with leading zeros (if index is <100)
                tileIndexStr = "%03i" % (tileIndex)
                # Set the string command to run.
                # This uses gdal's gdal_translate function to create the tiles.
                com_string = "gdal_translate -of GTIFF -srcwin " + str(i)+ ", " + str(j) + ", " + str(tile_size_x) + ", " + str(tile_size_y) + " " + str(in_path) +'\\'+ str(input_filename) + " " + str(out_path) + str(output_filename) + str(tileIndexStr) + ".tif"
                # Run the command
                os.system(com_string)
    else:
        shape_file_path = split_input_path[0]+'\\'+split_input_path[1]+'\\'+split_input_path[2]+'\\'+split_input_path[3]+'\\Shape\\{}.shp'.format(split_input_path[2])
        print(shape_file_path)
        shapefile = gpd.read_file(shape_file_path)
        for i in range(len(shapefile)):
            minx, miny, maxx, maxy = shapefile['geometry'].iloc[i].bounds
            print('{},{},{},{}'.format(minx,miny,maxx,maxy))
            #com_string = "gdal_translate -of GTIFF -r bilinear -outsize 5000 5000 -projwin_srs epsg:4326 -projwin " + str(minx) + " " + str(maxy) + " " + str(maxx) + " " + str(miny) + " " + str(in_path) + str(input_filename) + " " + str(out_path) + str(shapefile['title'].iloc[i]) + ".tif"
            com_string = "gdal_translate -of GTIFF -projwin " + str(minx) + ' ' + str(maxy) + ' ' + str(maxx) + ' ' + str(miny) + ' -projwin_srs EPSG:32615 ' +str(in_path)+'\\'+ str(input_filename) + " " + str(out_path) + str(output_filename) + str(shapefile['tile'].iloc[i][-3:]) + ".tif"
            os.system(com_string)

    # Set the the new input path (where the tiles were exported to)
    inPath = out_path + "*.tif"

    # Get the individual paths to tiles
    tileFilePaths = glob.glob(inPath)
    # Sort the tiles so that they are in alphabetical order
    tileFilePaths.sort()

    # Create array to hold the null tile paths
    nullTilePaths = []

    # Check each tile in the tile paths to determine if they are null
    for tilePath in tileFilePaths:
        print(tilePath)
        # Open the current tile
        with rasterio.open(tilePath) as raster:
            if tilePath.split('\\')[5] == 'BandTCI':
                # Open the RGB bands and get their data
                rBand = raster.read(1)
                gBand = raster.read(2)
                bBand = raster.read(3)
                # If all data in all three bands is null, append the corresponding
                # tile path to the nullTilePaths array
                if np.all(rBand == 0) and np.all(gBand == 0) and np.all(bBand == 0):
                    nullTilePaths.append(tilePath)
            else:
                Band = raster.read(1)
                #if np.all(Band == 0):
                #    nullTilePaths.append(tilePath)

    for tilePath in nullTilePaths:
        # Make sure that the file exists
        if os.path.exists(tilePath):
            # If the file exists, delete it
            os.remove(tilePath)
        # If the file path doesn't exist, print that out
        else:
            print("File does not exist: " + tilePath)
    # Get the paths to all tiles now that the null tiles have been removed
    tileFilePaths = glob.glob(inPath)
    # Sort the tile paths in alphabetical order
    tileFilePaths.sort()

    # For each of the tiles,
    # rename the file to account for the now removed tiles
    for tilePathIndex in range(len(tileFilePaths)):
        # Create the buffered string (with leading zeros) to be used
        # for the new file names
        tileIndexStr = "%03i" % (tilePathIndex + 1)

        # Get the tile directory and tile name from the current file path
        tileDir, tileName = os.path.split(tileFilePaths[tilePathIndex])

        # Get the beginning of the file name (before the index string)
        fileNameBeg = tileName.split('.')[0]
        fileNameBeg = fileNameBeg[:-4]
        # Create the new file name by taking the beginning of the file name
        # and appending the new index to the end (as well as the extension)
        newFileName = fileNameBeg + "_" + tileIndexStr + ".tif"

        # Finally, rename the file with the new file name
        os.rename(tileFilePaths[tilePathIndex], tileDir + "\\" + newFileName)
    # Get paths to individual TIF files
    if False:
        inFiles = glob.glob(inPath)
        # Sort the paths
        inFiles.sort()
        print(inFiles)
        # Open first file with rasterio to get the crs (assumes
        # crs same for all tiles)
        test = rxr.open_rasterio(inFiles[0]).squeeze()
        tileCrs = test.rio.crs
        print(tileCrs)
        # Create dictionary to hold bounding box polygon data
        boundingBoxData = {'tile': [], 'geometry' : []}

        # For each tile, create a bounding box and append it to
        # the dictionary we just declared
        for tileIndex in range(len(inFiles)):

            # Open the raster
            data = gdal.Open(inFiles[tileIndex], GA_ReadOnly)
            # Get the geo transform of the raster
            geoTransform = data.GetGeoTransform()
            # Get the minimum x value (min longitude)
            minx = geoTransform[0]
            # Get the maximum y value (max latitude)
            maxy = geoTransform[3]
            # Get the max x value (max longitude)
            maxx = minx + geoTransform[1] * data.RasterXSize
            # Get the min y value (min latitude)
            miny = maxy + geoTransform[5] * data.RasterYSize

            # Set the list of latitudes for bounding box
            lat_point_list = [maxy, maxy, miny, miny, maxy]
            # Set the list of longitudes for bounding box
            lon_point_list = [minx, maxx, maxx, minx, minx]

            # Get the tile name
            tileName = os.path.split(inFiles[tileIndex])[1]

            # Remove the extension from the tile name and save to variable
            tileNameNoExt = tileName.split('.')[0]

            boundingBoxData['tile'].append(tileNameNoExt)
            boundingBoxData['geometry'].append(Polygon(zip(lon_point_list, lat_point_list)))
        # Create a geopandas data frame using the previously created polygon and crs
        polygon = gpd.GeoDataFrame(boundingBoxData, crs=tileCrs)
        print(output_filename)
        # Save the polygon to shapefile
        polygon.to_file(filename=out_path+'\\'+output_filename[:-1] + '.shp', driver="ESRI Shapefile")

for in_path in all_input_path:
#in_path = "D:\\Study_Areas\\SA2\\20190831\\10m\\Band02\\"
#out_path = "D:\\Study_Areas\\SA2\\20190831\\10m\\Band02\\subdivided\\"
    out_path = in_path+'\\subdivided\\'
    #projection = {'EPSG:32652': 52}
    subdivide(in_path,out_path)
 