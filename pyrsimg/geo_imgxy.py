### ----- 
# author: luo xin, 
# creat: 2021.6.15, modify: 2024.6.8
# des: image location transform between different coordinate system. 
# -----

import numpy as np
from osgeo import osr, ogr
import math

def get_utm_zone(lon):
  '''
  des: get utm zone from the given wgs84 coordinates.
  lon: the given longitute, should be in the range of [-180, 180].
  return: utm_zone number.
  '''
  utm_zone = np.floor(lon/6)+31
  return int(utm_zone)

def coor2coor(srs_from, srs_to, x, y):
    """
    Transform coordinates from srs_from to srs_to
    input:
        srs_from and srs_to: EPSG number, (e.g., 4326, 3031)
        x and y are x-coord and y-coord corresponding to srs_from and srs_to    
    return:
        x-coord and y-coord in srs_to 
    """
    sr_in = osr.SpatialReference(); sr_in.ImportFromEPSG(int(srs_from))    
    sr_out = osr.SpatialReference(); sr_out.ImportFromEPSG(int(srs_to))     
    if int(srs_from) == 4326:
        sr_in.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    if int(srs_to) == 4326:
        sr_out.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(float(x), float(y))
    point.AssignSpatialReference(sr_in)
    point.TransformTo(sr_out)
    return (point.GetX(), point.GetY())

def geo2imagexy(x, y, gdal_trans, rsimg_array, integer=True):
    '''
    des: from georeferenced location (i.e., lon, lat) to image location(col,row).
    note: the coordinate system should be same between x/y and gdal_trans.
    input:
        gdal_trans: obtained by gdal.Open() and .GetGeoTransform(), or by geotif_io.readTiff()['geotrans']
        x: project or georeferenced x, i.e.,lon
        y: project or georeferenced y, i.e., lat
        rsimg_array: np.array(), remote sensing image array, shape = (row, col)/(row, col, bands)
    return: 
        image col and row corresponding to the georeferenced location.
    '''
    a = np.array([[gdal_trans[1], gdal_trans[2]], [gdal_trans[4], gdal_trans[5]]])
    b = np.array([x - gdal_trans[0], y - gdal_trans[3]])
    col_img, row_img = np.linalg.solve(a, b)
    if integer:
        row_img, col_img = np.floor(row_img).astype('int'), np.floor(col_img).astype('int')
    ## Mask out the points outside the image.
    ids_out = np.where((row_img>=rsimg_array.shape[0]) | (col_img>=rsimg_array.shape[1]))[0]
    while len(ids_out) > 0:
        raise IndexError('The x and y out of image range')
    return row_img, col_img


def imagexy2geo(row, col, gdal_trans):
    '''
    input: 
        img_gdal: GDAL data (read by gdal.Open()
        row and col are corresponding to input image (dataset)
    :return:  
        geographical coordinates (left up of pixel)
    '''
    x = gdal_trans[0] + col * gdal_trans[1] + row * gdal_trans[2]
    y = gdal_trans[3] + col * gdal_trans[4] + row * gdal_trans[5]
    return x, y


def deg2meter_resolution(degree_res, center_lat=0):
    """
    des: convert degree resolution to meter resolution
    params:
        degree_res (float): resolution in degrees
        center_lat (float): center latitude (default is equator)
    return:
        tuple: (resolution in meters along longitude, and for latitude)
    """
    R = 6371000  # mean radius of the Earth in meters    
    lat_res_m = degree_res * (math.pi * R / 180)  # convert latitude resolution (constant)
    lon_res_m = degree_res * (math.pi * R / 180) * math.cos(math.radians(center_lat))  # convert longitude resolution (varies with latitude)
    return (lon_res_m, lat_res_m)


def meter2deg_resolution(meter_res, center_lat=0):
    """
    des: convert meter resolution to degree resolution
    params:
        meter_res (float): resolution in meters
        center_lat (float): reference latitude (default is equator)    
    return:
        tuple: (resolution in degrees for longitude, and for latitude)
    """
    # radius of the Earth in meters
    R = 6371000
    lat_res_deg = (meter_res * 180) / (math.pi * R)  # convert latitude resolution (constant)
    lon_res_deg = (meter_res * 180) / (math.pi * R * math.cos(math.radians(center_lat)))  # convert longitude resolution (varies with latitude)
    return (lon_res_deg, lat_res_deg)