'''Subdivides the tiles'''
import rasterio
import glob
import numpy as np
import os, gdal
from gdalconst import GA_ReadOnly
import rioxarray as rxr
from shapely.geometry import Polygon
import geopandas as gpd


# Input path for directory to TIF file(merge file)
    #in_path = '/mnt/passport/Irrigation/'
in_path = '/media/18TB-HDD1/Levee_mapping_project/TempMerge/'
# Name of TIF file
input_filename = 'AR_2015_Randolph.tif'
#input_filename = '{}.tif'.format(merged)

# Output path for directory for subdivided TIFs.(where we want the output)
#out_path = '/mnt/passport/Irrigation/MS_163Yazoo/'
out_path = '/media/18TB-HDD1/Aniv_test/AR_2015_Randolph/'
# File name for subdivided TIFs (indexes will be
# appended to end later)
output_filename = 'Randolph_'

# Desired num pixels for each subdivided TIF in x direction
tile_size_x = 5000
# Desired num pixels for each subdivided TIF in y direction
tile_size_y = 5000

# Open the input TIF
ds = gdal.Open(in_path + input_filename)
#ds=image
# Get the first band (used for calculating size,
# so this assumes all bands have same size)
band = ds.GetRasterBand(1)
# Set size of TIF in x direction
xsize = band.XSize
# Set size of TIF in y direction
ysize = band.YSize

# Tile index, used to name the tiles
tileIndex = 0

# i & j here correspond to the x-value and y-value, respectively, of the
# top-left pixel index for each of the newly created tiles (e.g. the first
# tile will have (i, j) = (0,0), the second will have
# (i, j) = (0, tile_size_x), etc.)
for i in range(0, ysize, tile_size_y):
    for j in range(0, xsize, tile_size_x):
        # Index the tileIndex used for labeling; this is done before the first
        # tile is created so that the first index label is 001 rather than 000
        tileIndex = tileIndex + 1
        # Convert the tile index to a string with leading zeros (if index is <100)
        tileIndexStr = "%03i" % (tileIndex)
        # Set the string command to run.
        # This uses gdal's gdal_translate function to create the tiles.
        com_string = "gdal_translate -of GTIFF -b 1 -b 2 -b 3 -projwin_srs 'EPSG:26915' -srcwin " + str(j)+ ", " + str(i) + ", " + str(tile_size_x) + ", " + str(tile_size_y) + " " + str(in_path) + str(input_filename) + " " + str(out_path) + str(output_filename) + str(tileIndexStr) + ".tif"

        # Run the command
        os.system(com_string)


# ## Find Null TIFs
# We don't want any TIFs in our dataset that have no relevant information
# as this increases processing time & storage requirements unnecessarily.
#
# Here, we will find these empty TIF files.gdal_merge.py -ot Byte -of GTiff -co COMPRESS=LZW -co BIGTIFF=YES -o /media/18TB-HDD1/Levee\ mapping\ project/TempMerge/FL-2015-StLucia/merged.tif /media/18TB-HDD1/Levee\ mapping\ project/4\ RawData/FL-2015-StLucia/*.tif


# In[3]:


# Set the new input path (where the tiles were exported to)
inPath = out_path + "*.tif"

# Get the individual paths to tiles
tileFilePaths = glob.glob(inPath)
# Sort the tiles so that they are in alphabetical order
tileFilePaths.sort()

# Create array to hold the null tile paths
nullTilePaths = []

# Check each tile in the tile paths to determine if they are null
for tilePath in tileFilePaths:
    # Open the current tile
    with rasterio.open(tilePath) as raster:
        # Open the RGB bands and get their data
        rBand = raster.read(1)
        gBand = raster.read(2)
        bBand = raster.read(3)
        # If all data in all three bands is null, append the corresponding
        # tile path to the nullTilePaths array
        if np.all(rBand == 0) and np.all(gBand == 0) and np.all(bBand == 0):
            nullTilePaths.append(tilePath)


# ## Remove Null TIFs
# Here, we remove the blank TIFs we previously found.

# In[4]:


# For each tile in the list of null tiles,
# delete it
for tilePath in nullTilePaths:
    # Make sure that the file exists
    if os.path.exists(tilePath):
        # If the file exists, delete it
        os.remove(tilePath)
    # If the file path doesn't exist, print that out
    else:
        print("File does not exist: " + tilePath)


    # ## Rename TIFs
# Now that the null TIFs have been removed, the index in the file names needs to be re-indexed.

# In[5]:


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
    fileNameBeg = tileName.split('_')[0]

    # Create the new file name by taking the beginning of the file name
    # and appending the new index to the end (as well as the extension)
    newFileName = fileNameBeg + "_" + tileIndexStr + ".tif"

    # Finally, rename the file with the new file name
    os.rename(tileFilePaths[tilePathIndex], tileDir + "/" + newFileName)


# ## Create Bounding boxes
# This section creates bounding boxes (an outline of the tile perimeter)
# for the tiles we previously created. This is saved as a single shapefile,
# where each of the tiles' bounding boxes are individual polygon features.
#
# This can be commented out or removed if not needed. It's your code now,
# do what you want.

# In[6]:


# Get paths to individual TIF files
inFiles = glob.glob(inPath)
# Sort the paths
inFiles.sort()

# Open first file with rasterio to get the crs (assumes
# crs same for all tiles)
test = rxr.open_rasterio(inFiles[0]).squeeze()
tileCrs = test.rio.crs

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
    tileNameNoExt = tileNameNoExt.split('_')[-1]

    boundingBoxData['tile'].append(tileNameNoExt)
    boundingBoxData['geometry'].append(Polygon(zip(lon_point_list, lat_point_list)))

# Create a geopandas data frame using the previously created polygon and crs
polygon = gpd.GeoDataFrame(boundingBoxData, crs=tileCrs)

# Save the polygon to shapefile
polygon.to_file(filename=output_filename + '.shp', driver="ESRI Shapefile")