import argparse
from datetime import datetime
import getopt
import sys

import ee
import ee.mapclient
import folium
from folium import plugins
from IPython.display import HTML, display
import geopandas
from matplotlib import pyplot as plt
import numpy as np
import pandas


# ee.Authenticate()
ee.Initialize()


def maskL8sr(image):
    # Bits 3 and 5 are cloud shadow and cloud
    cloudShadowBitMask = (1 << 3)
    cloudsBitMask = (1 << 5)
    # Get pixel QA band
    qa = image.select('BQA')
    # Both flags should be zero, indicating clear conditions
    mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0).And(
        qa.bitwiseAnd(cloudsBitMask).eq(0)
    )

    return image.updateMask(mask)


def getLandsatCollection():
    """
    merge landsat 5, 7, 8 collection 1 
    tier 1 SR imageCollections and standardize band names
    """
    ## standardize band names
    bn8 = ['B1', 'B2', 'B3', 'B4', 'B6', 'pixel_qa', 'B5', 'B7']
    bn7 = ['B1', 'B1', 'B2', 'B3', 'B5', 'pixel_qa', 'B4', 'B7']
    bn5 = ['B1', 'B1', 'B2', 'B3', 'B5', 'pixel_qa', 'B4', 'B7']
    bns = ['uBlue', 'Blue', 'Green', 'Red', 'Swir1', 'BQA', 'Nir', 'Swir2']

    # create a merged collection from landsat 5, 7, and 8
    ls5 = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR").select(bn5, bns)

    ls7 = (ee.ImageCollection("LANDSAT/LE07/C01/T1_SR")
           .filterDate('1999-04-15', '2003-05-30')
           .select(bn7, bns))

    ls8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR").select(bn8, bns)

    merged = ls5.merge(ls7).merge(ls8)

    return(merged)


def getImage(year_month, polygon):
    # Get begining and end
    begin = year_month + '-01'
    end = year_month + '-28'

    # Filter image collection by
    return allLandsat.map(
        maskL8sr
    ).filterDate(
        begin, end 
    ).median().clip(
        polygon
    )


def calcMean(img, bands, scale):
    # gets the mean NDVI for the area in this img
    return img.reduceRegion(
        ee.Reducer.median(), 
        img.geometry().getInfo(), 
        scale 
    )


def add_ee_layer(self, ee_object, vis_params, name):
    
    try:    
        # display ee.Image()
        if isinstance(ee_object, ee.image.Image):    
            map_id_dict = ee.Image(ee_object).getMapId(vis_params)
            folium.raster_layers.TileLayer(
            tiles = map_id_dict['tile_fetcher'].url_format,
            attr = 'Google Earth Engine',
            name = name,
            overlay = True,
            control = True
            ).add_to(self) 
        # display ee.ImageCollection()
        elif isinstance(ee_object, ee.imagecollection.ImageCollection):    
            ee_object_new = ee_object.mosaic()
            map_id_dict = ee.Image(ee_object_new).getMapId(vis_params)
            folium.raster_layers.TileLayer(
            tiles = map_id_dict['tile_fetcher'].url_format,
            attr = 'Google Earth Engine',
            name = name,
            overlay = True,
            control = True
            ).add_to(self)
        # display ee.Geometry()
        elif isinstance(ee_object, ee.geometry.Geometry):    
            folium.GeoJson(
            data = ee_object.getInfo(),
            name = name,
            overlay = True,
            control = True
        ).add_to(self)
        # display ee.FeatureCollection()
        elif isinstance(ee_object, ee.featurecollection.FeatureCollection):  
            ee_object_new = ee.Image().paint(ee_object, 0, 2)
            map_id_dict = ee.Image(ee_object_new).getMapId(vis_params)
            folium.raster_layers.TileLayer(
            tiles = map_id_dict['tile_fetcher'].url_format,
            attr = 'Google Earth Engine',
            name = name,
            overlay = True,
            control = True
        ).add_to(self)
        else:
            folium.GeoJson(
                data = v,
                name = k
            ).add_to(self)
    except:
        print("Could not display {}".format(name))


# Folium Init
basemaps = {
    'Google Maps': folium.TileLayer(
        tiles = 'https://mt1.google.com/vt/lyrs=m&x={x}&y={y}&z={z}',
        attr = 'Google',
        name = 'Google Maps',
        overlay = True,
        control = True
    ),
    'Google Satellite': folium.TileLayer(
        tiles = 'https://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}',
        attr = 'Google',
        name = 'Google Satellite',
        overlay = True,
        control = True
    ),
    'Google Terrain': folium.TileLayer(
        tiles = 'https://mt1.google.com/vt/lyrs=p&x={x}&y={y}&z={z}',
        attr = 'Google',
        name = 'Google Terrain',
        overlay = True,
        control = True
    ),
    'Google Satellite Hybrid': folium.TileLayer(
        tiles = 'https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}',
        attr = 'Google',
        name = 'Google Satellite',
        overlay = True,
        control = True
    ),
    'Esri Satellite': folium.TileLayer(
        tiles = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
        attr = 'Esri',
        name = 'Esri Satellite',
        overlay = True,
        control = True
    )
}

# Add EE drawing method to folium.
folium.Map.add_ee_layer = add_ee_layer


# Get collection of all landsat images
allLandsat = ee.ImageCollection(getLandsatCollection())

# Folium Vis
lanVis = {
    'bands': ['Red', 'Green', 'Blue'],
    'min': 0,
    'max': 3000,
#    'palette':['225ea8','41b6c4','a1dab4','ffffcc']
}

# Make test polygon
coords = [
    (-122.824519, 38.683394),
    (-122.822585, 38.681848),
    (-122.820949, 38.683173),
    (-122.822873, 38.684691),
]
polygon = ee.Geometry.Polygon(coords)

# Set up how to sample the points
years = [i for i in range(1985, 2021)]
months = [i for i in range(1, 13)]
bands = ['uBlue', 'Blue', 'Green', 'Red', 'Swir1', 'Nir', 'Swir2']

# Get collection
data = {band: [] for band in bands}
data['dt'] = []
for year in years:
    for month in months:
        year_month = str(year) + '-' + str(month)
        image = getImage(year_month, polygon)
        mean = calcMean(image, bands, 30).getInfo()
        print(mean)
        if mean:
            data['dt'].append(datetime(year, month, 1))
            for band in bands:
                data[band].append(mean[band])

df = pandas.DataFrame(data).dropna(how='any')

df.to_csv('landsat_data.csv')

# Load in climate data
climate = pandas.read_csv('SantaRosaClimate.csv')
climate = climate[['DATE', 'PRCP', 'TAVG','TMAX', 'TMIN']]
climate['dt'] = pandas.to_datetime(climate['DATE'])

septs = [str(year) + '-09-01' for year in years]
fig, ax = plt.subplots(1, 1)
for band in bands:
    ax.plot(df['dt'], df[band], label=band)
for sept in septs:
    ax.axvline(x=sept, linewidth=0.5, color='black', linestyle='dashed')

#ax2 = ax.twinx()
#ax2.plot(climate['dt'], climate['PRCP'], label='Precipitation')
#
#ax3 = ax.twinx()
#ax3.plot(climate['dt'], climate['TAVG'], label='Tavg')
#ax3.plot(climate['dt'], climate['TMAX'], label='Tmax')
#ax3.plot(climate['dt'], climate['TMIN'], label='Tmin')

ax.legend()
plt.show()

# Visualization
my_map = folium.Map(
    location=(coords[0][1], coords[0][0]), 
    zoom_start=50, 
    height=500,
)

# Add custom basemaps
basemaps['Google Maps'].add_to(my_map)
basemaps['Google Satellite Hybrid'].add_to(my_map)

# Add to mao
my_map.add_ee_layer(image, lanVis, 'Landsat')
my_map.add_ee_layer(polygon, {}, 'Area')

# Add a layer control panel to the map.
my_map.add_child(folium.LayerControl())

# Add fullscreen button
plugins.Fullscreen().add_to(my_map)

# Display the map.
display(my_map)

my_map.save('index.html')
