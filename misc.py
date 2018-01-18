import subprocess
import glob
import os
import gdal
from shapely import geometry, wkb
import ogr
import numpy as np
import xml.etree.ElementTree
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import cv2
import osr
import itertools
import paramiko
import socket
from stat import S_ISDIR
from fmask import landsatangles
import rios
from osgeo import gdalnumeric
from PIL import Image, ImageDraw
from scipy.stats import norm
import matplotlib.mlab as mlab
import datetime
import warnings


def get_granule_info(search_folder, band_s2=3):

    """
    Get basic info on the input image
    
    :param search_folder: path to the granule or root folder of the input dataset (S-2 or L-8)
    :param band_s2: band number of the S2 band to be used for metadata and crs info, 1=green, 2=green, 3=red
    :return: Dictionary with satellite type, paths to the bands, and acquisition time, date and metadata
    
    # DEBUG
    search_folder = subfolder_list[0]
    search_folder =  '/home/stumpf/Data/Nepal/Durham/coseismic/LC081410412015060101T1-SC20171120065236'
    
    """

    # allocate output dictionary
    image_info = {}

    # check for L-8 granule metadata
    metadata_path = os.path.join(search_folder, '*MTL.txt')
    metadata_path = glob.glob(metadata_path)

    if metadata_path:

        print('Found L-8 granule...')
        image_info['satellite'] = 'L8'
        image_info['granule_path'] = search_folder
        image_info['metadata'] = metadata_path[0]

        # get band paths
        file_names = []
        lines = open(metadata_path[0], "r")
        for line in lines:
            if "FILE_NAME_BAND" in line: file_names.append(line.split("\"")[1])
            if "DATE_ACQUIRED" in line: image_info['date'] = line.split(" ")[-1].split('\n')[0]
            if "SCENE_CENTER_TIME" in line: image_info['time'] = line.split(" ")[-1].split('\n')[0]
        lines.close()
        # filter out quality band
        file_names = [file_name for file_name in file_names if not file_name.endswith("BQA.TIF")]
        bands_paths = []
        for i in file_names:
            bands_paths.append(os.path.join(search_folder, i))
        image_info['band_paths'] = bands_paths
        image_info['match_band'] = bands_paths[7]

        # get, geotransform, bounding boxes and SRS
        ds = gdal.Open(image_info['match_band'])
        gt = ds.GetGeoTransform()
        proj_wkt = ds.GetProjection()
        image_info['gt'] = gt
        image_info['proj'] = proj_wkt

        # get aproximate azimuth angle
        imgInfo = rios.fileinfo.ImageInfo(bands_paths[7])
        corners = landsatangles.findImgCorners(bands_paths[7], imgInfo)
        nadirLine = landsatangles.findNadirLine(corners)
        satAzimuth = landsatangles.satAzLeftRight(nadirLine)
        angle = abs(90 - 180 / np.pi * satAzimuth[0])
        image_info['azimuth'] = angle

    elif os.path.basename(search_folder).startswith('S2') or os.path.basename(search_folder).startswith('L1C'):
        image_info['satellite'] = 'S2'

        bands_paths = glob.glob(os.path.join(search_folder, 'IMG_DATA', "*.jp2"))

        if not bands_paths:
            granule_path = glob.glob(os.path.join(search_folder, 'GRANULE', '*'))[0]
            if not granule_path:
                raise Exception('Could not find any granule in ' + search_folder)
            if len(glob.glob(os.path.join(search_folder, 'GRANULE', '*'))) > 1:
                raise Exception('Found multiple granule in ' + search_folder +
                                'Cannot figure out which one to use by myself. Please specify full granule path')
            bands_paths = glob.glob(os.path.join(granule_path, 'IMG_DATA', "*.jp2"))
        else:
            granule_path = search_folder

        print('Found S-2 granule...')
        image_info['granule_path'] = granule_path
        bands_paths.sort()
        band_order = [0, 1, 2, 3, 4, 5, 6, 7, 12, 8, 9, 10, 11]
        bands_paths = [bands_paths[j] for j in band_order]
        image_info['band_paths'] = bands_paths
        image_info['match_band'] = bands_paths[band_s2]

        # get, geotransform, bounding boxes and SRS
        ds = gdal.Open(image_info['match_band'])
        gt = ds.GetGeoTransform()
        proj_wkt = ds.GetProjection()
        image_info['gt'] = gt
        image_info['proj'] = proj_wkt

        # get date
        xml_path = glob.glob(os.path.join(granule_path, '*.xml'))
        image_info['metadata'] = xml_path[0]
        e = xml.etree.ElementTree.parse(xml_path[0]).getroot()
        acquisition_time = e[0][3].text
        image_info['date'] = (acquisition_time[0:10])
        image_info['time'] = (acquisition_time[12:])

        # get gml path
        qi_path = os.path.join(granule_path, 'QI_DATA')
        gml_path = glob.glob(os.path.join(qi_path, '*DETFOO*B0' + str(band_s2 + 1) + '*.gml'))[0]
        image_info['azimuth'] = gml_path

    else:
        raise Exception('Could not identify a valid S-2 or L-8 product.')

    if not image_info:
        raise Exception(('Unable to gather granule info for folder ' + search_folder))

    return image_info


def get_all_granule_info(input_folder, band_s2=3):

    """
    Wrapper to get_granule_info to iterate over all subfolders in an input folder
    Args:
        input_folder: full path to the folder containting the image datasets
        band_s2: band number of the S2 band to be used for metadata and crs info, 1=green, 2=green, 3=red

    Returns: list for granule info dictionaries

    """

    # get granule infos
    subfolder_list = [d for d in os.listdir(input_folder) if os.path.isdir(os.path.join(input_folder, d)) and
                      (d.endswith('SAFE') or d.startswith('LC8'))]
    subfolder_list = [os.path.join(input_folder, subfolder) for subfolder in subfolder_list]

    granule_info_list = []
    for subfolder in subfolder_list:
        granule_info_list.append(get_granule_info(subfolder, band_s2=band_s2))

    return granule_info_list


def pixel2geo(dx, dy, transform):
    x = dx * transform[1] + transform[0]
    y = dy * transform[5] + transform[3]
    return x, y


def world_to_pixel(geo_matrix, x, y):

    """
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
    the pixel location of a geospatial coordinate; from:
    http://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html#clip-a-geotiff-with-shapefile
    """

    ulX = geo_matrix[0]
    ulY = geo_matrix[3]
    xDist = geo_matrix[1]
    pixel = int((x - ulX) / xDist)
    line = int((ulY - y) / xDist)
    return pixel, line


def granule2vrt(granule_info, pixel_size_x, pixel_size_y,
                select_bands=["B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B09", "B10", "B11", "B12"],
                output_folder=None):

    """
    Wrapper to gdalbuildvrt to construct a virtual raster for a specified S-2 Granule

    :param granule_info: dictionary with granule infos
    :param pixel_size_x: pixel size of the output vrt
    :param pixel_size_y: pixel size of the output vrt
    :param select_bands: list of strings with bands to select
    :param output_folder: optional to specify the output folder to which the VRT will be written
    :return: virtual raster on disk with specified resolution, and path to it

    """

    band_paths_subset = []
    for file in granule_info['band_paths']:
        do_use = [band_name for band_name in select_bands if band_name in file]
        if do_use:
            band_paths_subset.append(file)

    # build virtual raster
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)
        vrt_path = os.path.join(output_folder, 'bands.vrt')
    else:
        vrt_path = os.path.join(os.path.dirname(granule_info['band_paths'][0]), 'bands.vrt')

    cmd = ['gdalbuildvrt',
           '-resolution', 'user',
           '-tr', str(pixel_size_x), str(pixel_size_y),
           '-separate', vrt_path]
    for j in band_paths_subset:
        cmd.append(j)
    subprocess.call(cmd)

    return vrt_path


def make_pair(granule_info_master, granule_info_slave):

    """
    Combine two granule info dictionaries into a master-slave pair
    
    :param granule_info_master: granule info dictionary of the master granule
    :param granule_info_slave: granule info dictionary of the slave granule
    :return: master slave dictionary
    """

    master_slave_info = {}
    master_slave_info['master'] = granule_info_master
    master_slave_info['slave'] = granule_info_slave

    return master_slave_info


def order_by_dates(granule_info_list):

    """
    Takes a list of granule info dictionaries and orders them by their sensing date
    Args:
        granule_info_list: list of granule info dictionaries

    Returns: date ordered list of granule info dictionaries

    """

    sensing_dates = []
    for granule_info in granule_info_list:
        sensing_dates.append(granule_info['date'])

    import numpy as np

    sensing_dates_order = np.argsort(sensing_dates)
    granule_info_list = [granule_info_list[i] for i in sensing_dates_order]

    return granule_info_list


def find_master_by_date(granule_info_list, master_date):
    """
    Finds index of master according to the specified date
    Args:
        granule_info_list: list of granule info dictionaries
        master_date: string of the master date in the form of 'yyyy-mm-dd'

    Returns: an integer for the master index

    """

    master_index = []
    for i, granule_info in enumerate(granule_info_list):

        if granule_info['date'] == master_date:
            master_index.append(i)

    if not master_index:
        raise Exception('Could not find the master with the specified sensing date ' + master_date)
    elif len(master_index) > 1:
        raise Exception('Found more than one granule with the specified sensing date ' + master_date)
    else:
        master_index = master_index[0]

    return master_index


def make_all_pairs_with_master(granule_info_list, master_index):
    """
    Pairs all granules in the list of granule info dictionaries with the granule specified through the master index
    Args:
        granule_info_list: list of granule info dictionaries
        master_index: an integer for the master index

    Returns: a list of master slave dictionaries

    """

    master_slave_info_list = []
    for i, granule_info_slave in enumerate(granule_info_list):

        if i == master_index:
            print('Skip master with master pairing')
            continue

        master_slave_info_list.append(make_pair(granule_info_list[master_index], granule_info_slave))

    return master_slave_info_list


def make_all_pairs_with_master_by_date(granule_info_list, master_date):

    """
    Generates a list of master slave dictionaries with the master specified by sensing date
    Args:
        granule_info_list: list of granule info dictionaries
        master_date: string of the master date in the form of 'yyyy-mm-dd'

    Returns: a list of master slave dictionaries

    """

    # order by dates
    granule_info_list = order_by_dates(granule_info_list)
    master_index = find_master_by_date(granule_info_list, master_date)
    master_slave_info_list = make_all_pairs_with_master(granule_info_list, master_index)

    return master_slave_info_list


def make_pair_sequence(granule_info_list, match_range = 1, backward_matching = True):

    """
    Make a list of master slave dictionaries matching granules in sequential temporal order
    and according to the defined steps size
    Args:
        granule_info_list: list of granule info dictionaries
        match_range: maximum temporal step size for forming pairs
        backward_matching: if True also inverse matchig will be applied

    Returns: list of master slave dictionaries

    """

    granule_info_list = order_by_dates(granule_info_list)
    master_slave_info_list = []
    for i, granule_info in enumerate(granule_info_list):
        # print(granule_info['date'])

        for j in range(match_range):
            # print(j)
            try:
                master_slave_info_list.append(make_pair(granule_info,granule_info_list[i+j+1]))
            except:
                pass
            if backward_matching:
                try:
                    master_slave_info_list.append(make_pair(granule_info_list[i + j + 1], granule_info))
                except:
                    pass

    return master_slave_info_list


def make_pairs_before_after(granule_info_list, split_date, backward_matching=False):

    """
    Make a list of master slave dictionaries matching granules in sequential temporal order
    and according to the defined steps size.
    Args:
        granule_info_list: list of granule info dictionaries
        split_date: date at which the time series will be split to form pairs bridging before and after
        backward_matching: if true also inverse pairs will be formed

    Returns:

    """


    granule_set_before = []
    granule_set_after = []

    for granule_info in granule_info_list:

        if (datetime.datetime.strptime(split_date, "%Y-%m-%d") - datetime.datetime.strptime(granule_info['date'],
                                                                                            "%Y-%m-%d")).days > 0:
            granule_set_before.append(granule_info)
        else:
            granule_set_after.append(granule_info)

    pairs_forward = list(itertools.product(granule_set_before, granule_set_after))
    pairs_backward = list(itertools.product(granule_set_after, granule_set_before))

    if backward_matching:
        pairs = pairs_forward + pairs_backward
    else:
        pairs = pairs_forward

    master_slave_info_list = []
    for pair in pairs:
        master_slave_info = {}
        master_slave_info['master'] = pair[0]
        master_slave_info['slave'] = pair[1]
        master_slave_info_list.append(master_slave_info)

    return master_slave_info_list


def generate_pan(master_slave_info, resampling_method='cubic', output_folder=None):

    """
    Generates bands for matching S-2/S-2 or L-8/S-2, handles differences in projection and resolution
    :param master_slave: a master-slave dictionary as generated by misc.get_satellite and misc.make_pair
    :param resampling_method: options for resampling as in http://www.gdal.org/gdalwarp.html
    :param optinal output folder to which the resampled bands will be written
    :return: master slave dictionary including the full pathes to the matching bands
    """

    # S-2/ S-2 case
    if master_slave_info['master']['satellite'] == master_slave_info['slave']['satellite'] == 'S2':

        master_in = master_slave_info['master']['match_band']
        slave_in = master_slave_info['slave']['match_band']

        if output_folder:
            pair_folder = os.path.basename(master_slave_info['master']['granule_path'][:-7]) + '_' +\
                          os.path.basename(master_slave_info['slave']['granule_path'][:-7])
            pair_folder = os.path.join(output_folder, pair_folder)
            os.makedirs(pair_folder, exist_ok=True)

        # if jp2 convert to tif
        if master_in.endswith('.jp2'):
            print('Converting to JP2 to TIFF... ')

            if output_folder:
                pan_master = os.path.join(pair_folder, os.path.splitext(os.path.basename(master_in))[0] + '.tif')
            else:
                pan_master = os.path.splitext(master_in)[0] + '.tif'

            subprocess.call(['gdal_translate',
                             '-of', 'GTiff',
                             '-ot', 'UInt16',
                             master_in,
                             pan_master])
        else:
            pan_master = master_in

        # if jp2 convert to tif
        if slave_in.endswith('.jp2'):
            print('Converting to JP2 to TIFF... ')

            if output_folder:
                pan_slave = os.path.join(pair_folder, os.path.splitext(os.path.basename(slave_in))[0] + '.tif')
            else:
                pan_slave = os.path.splitext(slave_in)[0] + '.tif'

            subprocess.call(['gdal_translate',
                             '-of', 'GTiff',
                             '-ot', 'UInt16',
                             slave_in,
                             pan_slave])
        else:
            pan_slave = slave_in

        # reproject if the geotransform does not match
        if master_slave_info['master']['gt'] == master_slave_info['slave']['gt']:
            print('Both grids already have the same extent and resolution. No warping performed...')
            master_slave_info['pan_master'] = pan_master
            master_slave_info['pan_slave'] = pan_slave

        else:

            # compute rest the output coordinates coordinates
            print('Warping slave image to master...')
            xmin = master_slave_info['master']['gt'][0]
            ymax = master_slave_info['master']['gt'][3]
            xres = master_slave_info['master']['gt'][1]
            yres = master_slave_info['master']['gt'][5]
            ds = gdal.Open(pan_master)
            n_x = ds.RasterXSize
            n_y = ds.RasterYSize
            ds = None
            xmax = xmin + n_x * xres
            ymin = ymax + n_y * yres

            pan_slave_out = os.path.splitext(pan_slave)[0] + '_aligned.tif'

            subprocess.call(['gdalwarp',
                             pan_slave,
                             '-t_srs', master_slave_info['master']['proj'],
                             '-te', str(xmin), str(ymin), str(xmax), str(ymax),
                             '-tr', str(xres), str(abs(yres)),
                             '-overwrite',
                             '-r', resampling_method,
                             pan_slave_out])

            master_slave_info['pan_master'] = pan_master
            master_slave_info['pan_slave'] = pan_slave_out

    if master_slave_info['master']['satellite'] == master_slave_info['slave']['satellite'] == 'L8':
        # raise Exception('...Co-registration L-8 to L-8 is currently not implemented.')

        master_in = master_slave_info['master']['match_band']
        slave_in = master_slave_info['slave']['match_band']

        if output_folder:
            pair_folder = os.path.basename(master_slave_info['master']['granule_path']) + '_' +\
                          os.path.basename(master_slave_info['slave']['granule_path'])
            pair_folder = os.path.join(output_folder, pair_folder)
            os.makedirs(pair_folder, exist_ok=True)

        # reproject if the geotransform does not match
        if master_slave_info['master']['gt'] == master_slave_info['slave']['gt']:
            print('Both grids already have the same extent and resolution. No warping performed...')
            master_slave_info['pan_master'] = master_in
            master_slave_info['pan_slave'] = slave_in

        else:

            # compute rest the output coordinates coordinates
            print('Warping slave image to master...')
            xmin = master_slave_info['master']['gt'][0]
            ymax = master_slave_info['master']['gt'][3]
            xres = master_slave_info['master']['gt'][1]
            yres = master_slave_info['master']['gt'][5]
            ds = gdal.Open(master_in)
            n_x = ds.RasterXSize
            n_y = ds.RasterYSize
            ds = None
            xmax = xmin + n_x * xres
            ymin = ymax + n_y * yres

            if output_folder:
                pan_slave_out = os.path.join(pair_folder, os.path.splitext(os.path.basename(slave_in))[0] + '_aligned.tif')
            else:
                pan_slave_out = os.path.splitext(slave_in)[0] + '_aligned.tif'

            subprocess.call(['gdalwarp',
                             slave_in,
                             '-t_srs', master_slave_info['master']['proj'],
                             '-te', str(xmin), str(ymin), str(xmax), str(ymax),
                             '-tr', str(xres), str(abs(yres)),
                             '-overwrite',
                             '-r', resampling_method,
                             pan_slave_out])

            master_slave_info['pan_master'] = master_in
            master_slave_info['pan_slave'] = pan_slave_out

    # S-2/L8 case, L8 is used as master but processing is performed in the S-2 reference sytem
    if master_slave_info['master']['satellite'] == 'L8' and master_slave_info['slave']['satellite'] == 'S2':
        # raise Exception('...Co-registration S-2 to L-8 is currently not implemented.')

        master_in = master_slave_info['master']['match_band']
        slave_in = master_slave_info['slave']['match_band']

        if output_folder:
            pair_folder = os.path.basename(master_slave_info['master']['granule_path']) + '_' +\
                          os.path.basename(master_slave_info['slave']['granule_path'][:-7])
            pair_folder = os.path.join(output_folder, pair_folder)
            os.makedirs(pair_folder, exist_ok=True)


        # generate panchromatic band for the slave image
        for i in (1, 2, 3):
            input_band = master_slave_info['slave']['band_paths'][i]
            ds = gdal.Open(input_band)
            print('Reading S2 band ' + input_band)
            band = ds.GetRasterBand(1)

            if i == 1:
                band_array = band.ReadAsArray()
                datatype = band.DataType
            else:
                band_array = np.dstack((band_array, band.ReadAsArray()))

        band_array = band_array[:, :, 0] * 0.2 + band_array[:, :, 1] * 0.4 + band_array[:, :, 2] * 0.4
        band_array = np.rint(band_array)
        band_array = band_array.astype(int)

        if output_folder:
            pan_slave = os.path.join(pair_folder, os.path.splitext(os.path.basename(slave_in))[0][:-3] + 'syn_pan.tif')
        else:
            pan_slave = os.path.splitext(slave_in)[0][:-3] + 'syn_pan.tif'

        # write new raster to disk
        print('Writing panchromatic output band to ' + pan_slave)
        driver = gdal.GetDriverByName('GTiff')
        out_raster = driver.Create(pan_slave, ds.RasterYSize, ds.RasterXSize, 1, datatype)
        out_raster.SetGeoTransform(master_slave_info['slave']['gt'])
        out_raster.SetProjection(ds.GetProjection())
        outband = out_raster.GetRasterBand(1)

        if hasattr(band_array, 'mask'):
            outband.WriteArray(band_array.filled())
        else:
            outband.WriteArray(band_array)
        outband.FlushCache()

        out_raster = None

        # project slave image to extent of the master
        if output_folder:
            pan_master = os.path.join(pair_folder, os.path.splitext(os.path.basename(master_in))[0] + '_pan.tif')
        else:
            pan_master = os.path.splitext(master_in)[0] + '_pan.tif'

        xmin = master_slave_info['slave']['gt'][0]
        ymax = master_slave_info['slave']['gt'][3]
        xres = master_slave_info['slave']['gt'][1]
        yres = master_slave_info['slave']['gt'][5]
        xmax = xmin + ds.RasterXSize * xres
        ymin = ymax + ds.RasterYSize * yres

        ds = None

        # TODO: check if transformation to reflectance makes any difference https://github.com/mapbox/rio-toa/blob/master/rio_toa/reflectance.py

        subprocess.call(['gdalwarp',
                         master_in,
                         '-t_srs', master_slave_info['master']['proj'],
                         '-te', str(xmin), str(ymin), str(xmax), str(ymax),
                         '-tr', str(xres), str(abs(yres)),
                         '-overwrite',
                         '-r', resampling_method,
                         pan_master])

        master_slave_info['pan_master'] = pan_master
        master_slave_info['pan_slave'] = pan_slave

    # L8/S2 case
    if master_slave_info['master']['satellite'] == 'S2' and master_slave_info['slave']['satellite'] == 'L8':
        print('....generating correlation bands for co-registration L-8 to S-2...')

        master_in = master_slave_info['master']['match_band']
        slave_in = master_slave_info['slave']['match_band']

        if output_folder:
            pair_folder = os.path.basename(master_slave_info['master']['granule_path'][:-7]) + '_' +\
                          os.path.basename(master_slave_info['slave']['granule_path'])
            pair_folder = os.path.join(output_folder, pair_folder)
            os.makedirs(pair_folder, exist_ok=True)

        # generate panchromatic band for the master image
        for i in (1, 2, 3):
            input_band = master_slave_info['master']['band_paths'][i]
            ds = gdal.Open(input_band)
            print('Reading S2 band ' + input_band)
            band = ds.GetRasterBand(1)

            if i == 1:
                band_array = band.ReadAsArray()
                datatype = band.DataType
            else:
                band_array = np.dstack((band_array, band.ReadAsArray()))

        band_array = band_array[:,:,0] * 0.2 + band_array[:,:,1] * 0.4 +  band_array[:,:,2] * 0.4
        band_array = np.rint(band_array)
        band_array = band_array.astype(int)

        if output_folder:
            pan_master = os.path.join(pair_folder, os.path.splitext(os.path.basename(master_in))[0][:-3] + 'syn_pan.tif')
        else:
            pan_master = os.path.splitext(master_in)[0][:-3] + 'syn_pan.tif'

        # write new raster to disk
        print('Writing panchromatic output band to ' + pan_master)
        driver = gdal.GetDriverByName('GTiff')
        out_raster = driver.Create(pan_master, ds.RasterYSize, ds.RasterXSize, 1, datatype)
        out_raster.SetGeoTransform(master_slave_info['master']['gt'])
        out_raster.SetProjection(ds.GetProjection())
        outband = out_raster.GetRasterBand(1)

        if hasattr(band_array, 'mask'):
            outband.WriteArray(band_array.filled())
        else:
            outband.WriteArray(band_array)
        outband.FlushCache()

        out_raster = None

        # project slave image to extent of the master
        if output_folder:
            pan_slave = os.path.join(pair_folder, os.path.splitext(os.path.basename(slave_in))[0] + '_pan.tif')
        else:
            pan_slave = os.path.splitext(slave_in)[0] + '_pan.tif'

        xmin = master_slave_info['master']['gt'][0]
        ymax = master_slave_info['master']['gt'][3]
        xres = master_slave_info['master']['gt'][1]
        yres = master_slave_info['master']['gt'][5]
        xmax = xmin + ds.RasterXSize * xres
        ymin = ymax + ds.RasterYSize * yres

        ds = None

        # TODO: check if transformation to reflectance makes any difference https://github.com/mapbox/rio-toa/blob/master/rio_toa/reflectance.py

        subprocess.call(['gdalwarp',
                         slave_in,
                         '-t_srs', master_slave_info['master']['proj'],
                         '-te', str(xmin), str(ymin), str(xmax), str(ymax),
                         '-tr', str(xres), str(abs(yres)),
                         '-overwrite',
                         '-r', resampling_method,
                         pan_slave])

        master_slave_info['pan_master'] = pan_master
        master_slave_info['pan_slave'] = pan_slave

    return master_slave_info


def generate_pan_all(master_slave_info_list, resampling_method='cubic', output_folder=None):
    """
    Iterates over a list of master slave dictionaries to generate bands for matching S-2/S-2 or L-8/S-2,
     handles differences in projection and resolution
    :param master_slave: a master-slave dictionary as generated by misc.get_satellite and misc.make_pair
    :param resampling_method: options for resampling as in http://www.gdal.org/gdalwarp.html
    :param optinal output folder to which the resampled bands will be written
    :return: list of master slave dictionary including the full pathes to the matching bands
    """

    master_slave_info_list_updated = []
    for master_slave_info in master_slave_info_list:
        master_slave_info_list_updated.append(generate_pan(master_slave_info,
                                                                resampling_method=resampling_method,
                                                                output_folder=output_folder))
    return master_slave_info_list_updated


# TODO replace check granules by something shorter based on the output of get_granule_info
def check_input(granules, band_s2=3):

    """
    Function to check the consistency of the extents of the input grids and to determine the target spatial reference
    system (i.e. the most frequent projection of the input files)

    Args:
        granules: List with full paths to the S-2 granules

    Returns:
        granules:           list of intersecting granules
        master_granule:     full path to the master granule which defines the spatial reference
        bands3:             list of full pathes to the bands3 which will be used for matching
        dates:              list of the sensing dates
        date_order:         index to order the granules according the their sensing dates
        band_s2:            flag to change the S-2 band for S-2 to S-2 registration, blue=2, green=3, red=4

    """

    # lexicographical order
    granules.sort
    for granule in granules:
        if (os.path.isdir(granule)):
            print('Found the folder ' + os.path.basename(granule))
        else:
            raise Exception('Could not find the folder ' + os.path.basename(granule))

    # get dates
    dates = []
    bands = []
    for granule in granules:
        bands_path = glob.glob(os.path.join(granule, 'IMG_DATA', "*.jp2"))
        bands_path.sort()
        bands.append(bands_path[band_s2-1])
        xml_path = glob.glob(os.path.join(granule, '*.xml'))
        e = xml.etree.ElementTree.parse(xml_path[0]).getroot()
        acquisition_time = e[0][3].text
        dates.append(acquisition_time[0:10])

    # get tile projections
    tile_proj = []
    for granule in granules:
        xml_path = glob.glob(os.path.join(granule, '*.xml'))
        e = xml.etree.ElementTree.parse(xml_path[0]).getroot()
        tile_proj.append(e[1][0][0].text)
    # find the most common tile projection
    master_proj = max(set(tile_proj), key=tile_proj.count)
    print('Master projection is ' + master_proj)

    # get tile geotransforms
    geotransforms = []
    for band in bands:
        ds = gdal.Open(band)
        geotransforms.append(ds.GetGeoTransform())
    # find the most common geotransform
    master_gt = max(set(geotransforms), key=geotransforms.count)
    print('Master geotransform is ' + str(master_gt))

    # get the index of the master
    satisfy_proj = [i for i, x in enumerate(tile_proj) if x == master_proj]
    satisfy_gt = [i for i, x in enumerate(geotransforms) if x == master_gt]
    satisfy_both = list(set(satisfy_proj) & set(satisfy_gt))
    master_index = satisfy_both[0]
    master_granule = granules[master_index]
    print('The master tile is ...')
    print(granules[master_index])

    # remove all datasets which do not intersect the master extent

    do_intersect = []
    for band in bands:

        # link data
        image1_ds = gdal.Open(bands[master_index])
        image2_ds = gdal.Open(band)

        # get geoinfos for first
        gt1 = image1_ds.GetGeoTransform()
        proj_1_wkt = image1_ds.GetProjectionRef()

        # get geoinfos for second
        gt2 = image2_ds.GetGeoTransform()
        proj_2_wkt = image2_ds.GetProjectionRef()

        # find each image's bounding box
        # r1 has left, top, right, bottom of dataset's bounds in geospatial coordinates.
        r1 = [gt1[0], gt1[3], gt1[0] + (gt1[1] * image1_ds.RasterXSize), gt1[3] + (gt1[5] * image1_ds.RasterYSize)]
        r2 = [gt2[0], gt2[3], gt2[0] + (gt2[1] * image2_ds.RasterXSize), gt2[3] + (gt2[5] * image2_ds.RasterYSize)]
        print('\t1 bounding box: %s' % str(r1))
        print('\t2 bounding box: %s' % str(r2))

        # if not the same projections > reproject
        if not proj_1_wkt == proj_2_wkt:
            print('\t The two grids have different projections. Reprojecting bounding box to ' + master_proj + ' ... ')

            proj_1 = osr.SpatialReference()
            proj_1.ImportFromWkt(proj_1_wkt)

            proj_2 = osr.SpatialReference()
            proj_2.ImportFromWkt(proj_2_wkt)

            # define transformation
            transform = osr.CoordinateTransformation(proj_2, proj_1)

            # apply transform
            ulx, uly, _ = transform.TransformPoint(r2[0], r2[1])
            urx, ury, _ = transform.TransformPoint(r2[2], r2[1])
            lrx, lry, _ = transform.TransformPoint(r2[2], r2[3])
            llx, lly, _ = transform.TransformPoint(r2[0], r2[3])

            # generate polygon
            poly2 = geometry.Polygon([(ulx, uly), (urx, ury), (lrx, lry), (llx, lly)])

        else:

            # generate polygon
            poly2 = geometry.Polygon([(r2[0], r2[1]), (r2[2], r2[1]), (r2[2], r2[3]), (r2[0], r2[3])])

        # generate polygon and intersect
        poly1 = geometry.Polygon([(r1[0], r1[1]), (r1[2], r1[1]), (r1[2], r1[3]), (r1[0], r1[3])])
        poly_intersect = poly1.intersection(poly2)

        # # plot to check
        # from matplotlib import pyplot as plt
        # from descartes import PolygonPatch
        # p = PolygonPatch(poly_intersect, facecolor='blue', alpha=1)
        # fig, ax = plt.subplots(figsize=(8, 8))
        # ax.add_patch(p)
        # ax.relim()
        # ax.set_aspect('equal', 'datalim')
        # ax.autoscale()

        do_intersect.append(poly_intersect.bounds is not None)

        if poly_intersect.bounds is None:
            print(band)
            print('\t...does not intersect the master tile and is removed.')

    # filter the input granules
    granules = list(itertools.compress(granules, do_intersect))
    bands = list(itertools.compress(bands, do_intersect))
    dates = list(itertools.compress(dates, do_intersect))
    # get order by date
    date_order = np.argsort(dates)

    return granules, master_granule, bands, dates, date_order


def find_minimum_bounding_box(vrt_1_path, vrt_2_path, aoi=None):

    """
    Finds minimum bounding box of the input images and the optional area of interest

    Args:
        vrt_1_path: full path to the first image
        vrt_2_path: full path to the second image
        aoi: optional, full path to the shapefile holding an AOI polygon, only the first polygon in the shp is considered

    Returns: tuple of bounding box coordinates
        x_left, y_down, x_right, y_up
    """

    # load data
    image1_ds = gdal.Open(vrt_1_path)
    image2_ds = gdal.Open(vrt_2_path)

    gt1 = image1_ds.GetGeoTransform()
    gt2 = image2_ds.GetGeoTransform()

    # find each image's bounding box
    # r1 has left, top, right, bottom of dataset's bounds in geospatial coordinates.
    r1 = [gt1[0], gt1[3], gt1[0] + (gt1[1] * image1_ds.RasterXSize), gt1[3] + (gt1[5] * image1_ds.RasterYSize)]
    r2 = [gt2[0], gt2[3], gt2[0] + (gt2[1] * image2_ds.RasterXSize), gt2[3] + (gt2[5] * image2_ds.RasterYSize)]
    print('\t1 bounding box: %s' % str(r1))
    print('\t2 bounding box: %s' % str(r2))

    poly1 = geometry.box(r1[0], r1[1], r1[2], r1[3], ccw=True)
    poly2 = geometry.box(r2[0], r2[1], r2[2], r2[3], ccw=True)

    poly_intersect = poly1.intersection(poly2)

    if aoi:
        source = ogr.Open(aoi)
        layer = source.GetLayer()
        feat_count = layer.GetFeatureCount()
        if feat_count > 1:
            print('Found more than 1 feature in the shapefile but only the first will be considered ...')
        feat = layer.GetNextFeature()
        geom = wkb.loads(feat.GetGeometryRef().ExportToWkb())
        poly_intersect = poly_intersect.intersection(geom)

        # plot to check
        # from matplotlib import pyplot as plt
        # from descartes import PolygonPatch
        # p = PolygonPatch(poly_intersect, facecolor='blue', alpha=1)
        # fig, ax = plt.subplots(figsize=(8, 8))
        # ax.add_patch(p)
        # ax.relim()
        # ax.autoscale()

    return poly_intersect.bounds


def align_first2second(raster_path_1, raster_path_2, resampling_method='cubic', output_folder=None):

    """
    Align first image to the second image. Account for differences in extent, pixel grid, and projection. NOT shift,
    translation, etc.
    :param raster_path_1: full path to the raster that should be aligned
    :param raster_path_2: full path of the raster used as reference
    :param resampling_method: see option for see http://www.gdal.org/gdalwarp.html
    :param output_folder: folder to which the output rasters will be written
    :return: full path to the aligned raster

    DEBUG
    raster_path_1 = cropped_files[0]
    raster_path_2 = cropped_files[2]
    output_folder = work_folder
    resampling_method='cubic'
    """

    # open first and get geoinfos
    ds_1 = gdal.Open(raster_path_1)
    gt1 = ds_1.GetGeoTransform()
    proj_slave = ds_1.GetProjection()

    # open second and get geoinfos
    ds_2 = gdal.Open(raster_path_2)
    gt2 = ds_2.GetGeoTransform()
    proj_master = ds_2.GetProjection()

    if gt1 == gt2:
        print('Both grids already have the same extent and resolution. No warping performed...')
        return raster_path_1

    xmin = gt2[0]
    ymax = gt2[3]
    xres = gt2[1]
    yres = gt2[5]
    n_x = ds_2.RasterXSize
    n_y = ds_2.RasterYSize

    # compute rest of the coordinates
    xmax = xmin + n_x * xres
    ymin = ymax + n_y * yres

    if output_folder:
        out_path = os.path.join(output_folder, os.path.basename(os.path.splitext(raster_path_1)[0]) + '_aligned.tif')
    else:
        out_path = os.path.splitext(raster_path_1)[0] + '_aligned.tif'

    subprocess.call(['gdalwarp',
                     raster_path_1,
                     '-t_srs', proj_master,
                     '-te', str(xmin), str(ymin), str(xmax), str(ymax),
                     '-te_srs', proj_master,
                     '-tr', str(xres), str(abs(yres)),
                     '-overwrite',
                     '-r', resampling_method,
                     out_path])

    return out_path


def plot_layer(lyr, title=''):

    # lyr = layer3
    # Get extent and calculate buffer size
    ext = lyr.GetExtent()
    xoff = (ext[1]-ext[0])/50
    yoff = (ext[3]-ext[2])/50

    # Prepare figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(title, fontsize=16)
    ax.set_xlim(ext[0]-xoff,ext[1]+xoff)
    ax.set_ylim(ext[2]-yoff,ext[3]+yoff)

    paths = []
    lyr.ResetReading()

    # Read all features in layer and store as paths
    for feat in lyr:
        geom = feat.geometry()
        codes = []
        all_x = []
        all_y = []
        for i in range(geom.GetGeometryCount()):
            # Read ring geometry and create path
            r = geom.GetGeometryRef(i)
            x = [r.GetX(j) for j in range(r.GetPointCount())]
            y = [r.GetY(j) for j in range(r.GetPointCount())]
            # skip boundary between individual rings
            codes += [mpath.Path.MOVETO] + \
                         (len(x)-1)*[mpath.Path.LINETO]
            all_x += x
            all_y += y
        path = mpath.Path(np.column_stack((all_x,all_y)), codes)
        paths.append(path)

    # Add paths as patches to axes
    for path in paths:
        patch = mpatches.PathPatch(path, \
                facecolor='blue', edgecolor='black')
        ax.add_patch(patch)

    ax.set_aspect(1.0)
    plt.show()


def set_geospatial_reference(input_raster_paths, ref):

    # add spatial reference to micmac outputs
    ds_ref = gdal.Open(ref)
    gt_ref = ds_ref.GetGeoTransform()
    proj_ref = ds_ref.GetProjection()

    for file_path in input_raster_paths:
        ds_px = gdal.Open(file_path, gdal.GA_Update)
        ds_px.SetProjection(proj_ref)
        ds_px.SetGeoTransform(gt_ref)
        ds_px = None


def getTranslationMatrix2d(dx, dy):

    """
    Returns a numpy affine transformation matrix for a 2D translation of
    (dx, dy)
    """
    return np.matrix([[1, 0, dx], [0, 1, dy], [0, 0, 1]])


def rotateImage(image, angle):

    """
    image: 'input image'
    angle: 'rotation angle'
    Rotates the given image about it's centre

    Result:
    The rotated image
    The coordinates of the 4 corners of the original image in the rotated reference frame
    If the input image is a masked array an additional array holding the rotated mask will be returned

    """

    image_size = (image.shape[1], image.shape[0])
    image_center = tuple(np.array(image_size) / 2)

    rot_mat = np.vstack([cv2.getRotationMatrix2D(image_center, angle, 1.0), [0, 0, 1]])
    trans_mat = np.identity(3)

    w2 = image_size[0] * 0.5
    h2 = image_size[1] * 0.5

    rot_mat_notranslate = np.matrix(rot_mat[0:2, 0:2])

    tl = (np.array([-w2, h2]) * rot_mat_notranslate).A[0]
    tr = (np.array([w2, h2]) * rot_mat_notranslate).A[0]
    bl = (np.array([-w2, -h2]) * rot_mat_notranslate).A[0]
    br = (np.array([w2, -h2]) * rot_mat_notranslate).A[0]

    x_coords = [pt[0] for pt in [tl, tr, bl, br]]
    x_pos = [x for x in x_coords if x > 0]
    x_neg = [x for x in x_coords if x < 0]

    y_coords = [pt[1] for pt in [tl, tr, bl, br]]
    y_pos = [y for y in y_coords if y > 0]
    y_neg = [y for y in y_coords if y < 0]

    right_bound = max(x_pos)
    left_bound = min(x_neg)
    top_bound = max(y_pos)
    bot_bound = min(y_neg)

    new_w = int(abs(right_bound - left_bound))
    new_h = int(abs(top_bound - bot_bound))
    new_image_size = (new_w, new_h)

    new_midx = new_w * 0.5
    new_midy = new_h * 0.5

    tl_orig = [new_w - new_midx + tl[0], new_h - new_midy - tl[1]]
    tr_orig = [new_w - new_midx + tr[0], new_h - new_midy - tr[1]]
    br_orig = [new_w - new_midx + br[0], new_h - new_midy - br[1]]
    bl_orig = [new_w - new_midx + bl[0], new_h - new_midy - bl[1]]

    extent_orig = np.vstack([tl_orig, tr_orig, br_orig, bl_orig])

    dx = int(new_midx - w2)
    dy = int(new_midy - h2)

    trans_mat = getTranslationMatrix2d(dx, dy)
    affine_mat = (np.matrix(trans_mat) * np.matrix(rot_mat))[0:2, :]
    result = cv2.warpAffine(image, affine_mat, new_image_size, flags=cv2.INTER_NEAREST, borderValue=-255)

    # transfer mask
    if hasattr(image, 'mask'):
        mask_rot = cv2.warpAffine(np.uint8(image.mask), affine_mat, new_image_size, flags=cv2.INTER_NEAREST, borderValue=-255)

        return result, extent_orig, mask_rot
    else:
        return result, extent_orig


def ordered_pairs(input_list, range=2):
    iter_range = []

    # pairwise combinations according to range
    for i, current_image in enumerate(input_list):

        # get the the following in the sequence according to range
        if i + 1 > len(input_list):
            pass  # skip last image

        if i + 1 + range > len(input_list):
            follow_images = input_list[(i + 1):]  # if range after the current image exceeds list length take all remaining
            for paired_image in follow_images:
                iter_range.append([current_image, paired_image])
        else:
            follow_images = input_list[(i + 1):(i + 1 + range)]
            for paired_image in follow_images:
                iter_range.append([current_image, paired_image])

    return iter_range


class SSHSession(object):
    # Usage:
    # Detects DSA or RSA from key_file, either as a string filename or a
    # file object.  Password auth is possible, but I will judge you for
    # using it. So:
    # ssh=SSHSession('targetserver.com','root',key_file=open('mykey.pem','r'))
    # ssh=SSHSession('targetserver.com','root',key_file='/home/me/mykey.pem')
    # ssh=SSHSession('targetserver.com','root','mypassword')
    # ssh.put('filename','/remote/file/destination/path')
    # ssh.put_all('/path/to/local/source/dir','/path/to/remote/destination')
    # ssh.get_all('/path/to/remote/source/dir','/path/to/local/destination')
    # ssh.command('echo "Command to execute"')

    def __init__(self, hostname, username='root', key_file=None, password=None):
        #
        #  Accepts a file-like object (anything with a readlines() function)
        #  in either dss_key or rsa_key with a private key.  Since I don't
        #  ever intend to leave a server open to a password auth.
        #
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.connect((hostname,22))
        self.t = paramiko.Transport(self.sock)
        self.t.start_client()
        keys = paramiko.util.load_host_keys(os.path.expanduser('~/.ssh/known_hosts'))
        key = self.t.get_remote_server_key()
        # supposed to check for key in keys, but I don't much care right now to find the right notation
        if key_file is not None:
            if isinstance(key,str):
                key_file=open(key,'r')
            key_head=key_file.readline()
            key_file.seek(0)
            if 'DSA' in key_head:
                keytype=paramiko.DSSKey
            elif 'RSA' in key_head:
                keytype=paramiko.RSAKey
            else:
                raise Exception("Can't identify key type")
            pkey=keytype.from_private_key(key_file)
            self.t.auth_publickey(username, pkey)
        else:
            if password is not None:
                self.t.auth_password(username,password,fallback=False)
            else: raise Exception('Must supply either key_file or password')
        self.sftp=paramiko.SFTPClient.from_transport(self.t)

    def command(self,cmd):
        #  Breaks the command by lines, sends and receives
        #  each line and its output separately
        #
        #  Returns the server response text as a string

        chan = self.t.open_session()
        chan.get_pty()
        chan.invoke_shell()
        chan.settimeout(20.0)
        ret=''
        try:
            ret+=chan.recv(1024)
        except:
            chan.send('\n')
            ret+=chan.recv(1024)
        for line in cmd.split('\n'):
            chan.send(line.strip() + '\n')
            ret+=chan.recv(1024)
        return ret

    def put(self,localfile,remotefile):
        #  Copy localfile to remotefile, overwriting or creating as needed.
        self.sftp.put(localfile,remotefile)

    def put_all(self,localpath,remotepath):
        #  recursively upload a full directory
        os.chdir(os.path.split(localpath)[0])
        parent=os.path.split(localpath)[1]
        for walker in os.walk(parent):
            try:
                self.sftp.mkdir(os.path.join(remotepath,walker[0]))
            except:
                pass
            for file in walker[2]:
                self.put(os.path.join(walker[0],file),os.path.join(remotepath,walker[0],file))

    def get(self,remotefile,localfile):
        #  Copy remotefile to localfile, overwriting or creating as needed.
        self.sftp.get(remotefile,localfile)

    def sftp_walk(self,remotepath):
        # Kind of a stripped down  version of os.walk, implemented for
        # sftp.  Tried running it flat without the yields, but it really
        # chokes on big directories.
        path=remotepath
        files=[]
        folders=[]
        for f in self.sftp.listdir_attr(remotepath):
            if S_ISDIR(f.st_mode):
                folders.append(f.filename)
            else:
                files.append(f.filename)
        print (path,folders,files)
        yield path,folders,files
        for folder in folders:
            new_path=os.path.join(remotepath,folder)
            for x in self.sftp_walk(new_path):
                yield x

    def get_all(self,remotepath,localpath):
        #  recursively download a full directory
        #  Harder than it sounded at first, since paramiko won't walk
        #
        # For the record, something like this would gennerally be faster:
        # ssh user@host 'tar -cz /source/folder' | tar -xz

        self.sftp.chdir(os.path.split(remotepath)[0])
        parent=os.path.split(remotepath)[1]
        try:
            os.mkdir(localpath)
        except:
            pass
        for walker in self.sftp_walk(parent):
            try:
                os.mkdir(os.path.join(localpath,walker[0]))
            except:
                pass
            for file in walker[2]:
                self.get(os.path.join(walker[0],file),os.path.join(localpath,walker[0],file))

    def write_command(self,text,remotefile):
        #  Writes text to remotefile, and makes remotefile executable.
        #  This is perhaps a bit niche, but I was thinking I needed it.
        #  For the record, I was incorrect.
        self.sftp.open(remotefile,'w').write(text)
        self.sftp.chmod(remotefile,755)


def clip_raster(rast, features_path, gt=None, nodata=-9999):

    """
    Clips a raster (given as either a gdal.Dataset or as a numpy.array
    instance) to a polygon layer provided by a Shapefile (or other vector
    layer). If a numpy.array is given, a "GeoTransform" must be provided
    (via dataset.GetGeoTransform() in GDAL). Returns an array. Clip features
    must be a dissolved, single-part geometry (not multi-part). Modified from:

    http://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html
    #clip-a-geotiff-with-shapefile

    Arguments:
        rast            A gdal.Dataset or a NumPy array
        features_path   The path to the clipping features
        gt              An optional GDAL GeoTransform to use instead
        nodata          The NoData value; defaults to -9999.
    """

    def array_to_image(a):
        """
        Converts a gdalnumeric array to a Python Imaging Library (PIL) Image.
        """
        i = Image.fromstring('L',(a.shape[1], a.shape[0]),
            (a.astype('b')).tostring())
        return i

    def image_to_array(i):
        """
        Converts a Python Imaging Library (PIL) array to a gdalnumeric image.
        """
        a = gdalnumeric.fromstring(i.tobytes(), 'b')
        a.shape = i.im.size[1], i.im.size[0]
        return a

    def world_to_pixel(geo_matrix, x, y):
        """
        Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
        the pixel location of a geospatial coordinate; from:
        http://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html#clip-a-geotiff-with-shapefile
        """
        ulX = geo_matrix[0]
        ulY = geo_matrix[3]
        xDist = geo_matrix[1]
        yDist = geo_matrix[5]
        rtnX = geo_matrix[2]
        rtnY = geo_matrix[4]
        pixel = int((x - ulX) / xDist)
        line = int((ulY - y) / xDist)
        return (pixel, line)

    # Can accept either a gdal.Dataset or numpy.array instance
    if not isinstance(rast, np.ndarray):
        gt = rast.GetGeoTransform()
        rast = rast.ReadAsArray()

    # Create an OGR layer from a boundary shapefile
    features = ogr.Open(features_path)
    if features.GetDriver().GetName() == 'ESRI Shapefile':
        lyr = features.GetLayer(os.path.split(os.path.splitext(features_path)[0])[1])

    else:
        lyr = features.GetLayer()

    # Get the first feature
    poly = lyr.GetNextFeature()

    # Convert the layer extent to image pixel coordinates
    minX, maxX, minY, maxY = lyr.GetExtent()
    # round to nearest pixel
    rounding_x =1/gt[1]
    rounding_y = 1/abs(gt[5])
    minX = round((minX * rounding_x) / rounding_x)
    maxX = round((maxX * rounding_x) / rounding_x)
    minY = round((minY * rounding_y) / rounding_y)
    maxY = round((maxY * rounding_y) / rounding_y)
    ulX, ulY = world_to_pixel(gt, minX, maxY)
    lrX, lrY = world_to_pixel(gt, maxX, minY)

    # Calculate the pixel size of the new image
    pxWidth = int(lrX - ulX)
    pxHeight = int(lrY - ulY)

    # If the clipping features extend out-of-bounds and ABOVE the raster...
    if gt[3] < maxY:
        # In such a case... ulY ends up being negative--can't have that!
        iY = ulY
        ulY = 0

    # Multi-band image?
    try:
        clip = rast[:, ulY:lrY, ulX:lrX]

    except IndexError:
        clip = rast[ulY:lrY, ulX:lrX]

    # Create a new geomatrix for the image
    gt2 = list(gt)
    gt2[0] = minX
    gt2[3] = maxY

    # Map points to pixels for drawing the boundary on a blank 8-bit,
    #   black and white, mask image.
    points = []
    pixels = []
    geom = poly.GetGeometryRef()
    pts = geom.GetGeometryRef(0)

    for p in range(pts.GetPointCount()):
        points.append((pts.GetX(p), pts.GetY(p)))

    for p in points:
        pixels.append(world_to_pixel(gt2, p[0], p[1]))

    raster_poly = Image.new('L', (pxWidth, pxHeight), 1)
    rasterize = ImageDraw.Draw(raster_poly)
    rasterize.polygon(pixels, 0) # Fill with zeroes

    # If the clipping features extend out-of-bounds and ABOVE the raster...
    if gt[3] < maxY:
        # The clip features were "pushed down" to match the bounds of the
        #   raster; this step "pulls" them back up
        premask = image_to_array(raster_poly)
        # We slice out the piece of our clip features that are "off the map"
        mask = np.ndarray((premask.shape[-2] - abs(iY), premask.shape[-1]), premask.dtype)
        mask[:] = premask[abs(iY):, :]
        mask.resize(premask.shape) # Then fill in from the bottom

        # Most importantly, push the clipped piece down
        gt2[3] = maxY - (maxY - gt[3])

    else:
        mask = image_to_array(raster_poly)

    # Clip the image using the mask
    try:
        clip = gdalnumeric.choose(mask, (clip, nodata))

    # If the clipping features extend out-of-bounds and BELOW the raster...
    except ValueError:
        # We have to cut the clipping features to the raster!
        rshp = list(mask.shape)
        if mask.shape[-2] != clip.shape[-2]:
            rshp[0] = clip.shape[-2]

        if mask.shape[-1] != clip.shape[-1]:
            rshp[1] = clip.shape[-1]

        mask.resize(*rshp, refcheck=False)

        clip = gdalnumeric.choose(mask, (clip, nodata))

    return (clip, ulX, ulY, gt2)


def iterative_gaussian_fit(sample, max_iterations=100, min_change_sigma=0.005, diagnose=False, sigma_level=2.575,
                           title='Sample distribution', xlim=[-1, 1], ylim=[0, 3]):

    """
    Iterative fits the sample population with a Gaussian and removes outliers outside the defined confidence interval
    Args:
        sample: Array with the samples
        max_iterations: Maximum number of iterations
        min_change_sigma: Minimum change in sigma to continue fitting
        diagnose: if True a plot will be shown
        sigma_level: sigma_level*sigma > defines threshold for outliers, 2.575 corresponds to a 99% confidence level
        title: title of the upper figure
        xlim: lower limit, upper limit of the x - axis
        ylim: lower limit, upper limit of the y - axis
        
    Returns:
        1D array of the sample population with outliers and the RMSE

    """

    # remove mask and array to vector
    if isinstance(sample, np.ma.MaskedArray):
        sample = sample.compressed()
    if sample.ndim > 1:
        sample = sample.flatten()

    # estimate initial sigma and RMSE
    (mu, sigma) = norm.fit(sample)
    RMSE = '%.3f' % (np.ma.sqrt(np.ma.mean(np.square(sample))))


    if diagnose: # plot
        lower_bin = mu-3*sigma
        upper_bin = mu+3*sigma
        fig, axes = plt.subplots(nrows=2)
        hist, bins = np.histogram(sample, range=[lower_bin, upper_bin], density=False, bins=100)
        bars = axes[0].bar((bins[:-1] + bins[1:]) / 2, hist, align='center', width=(bins[1] - bins[0]))
        axes[0].set_title(title)
        axes[1].set_ylim(ylim[0], ylim[1])
        axes[0].set_xlim(xlim[0], xlim[1])
        text = axes[0].text(0.9, 0.9, 'RMSE=' + str(RMSE),
                            horizontalalignment='right',
                            verticalalignment='top',
                            transform=axes[0].transAxes)
        # text = axes[0].text(-1., 1000000., 'RMSE=' + str(RMSE), fontsize=12)
        y = mlab.normpdf( bins, mu, sigma)
        axes[1].plot(bins, y, 'r--', linewidth=2)
        axes[1].set_title('Gaussian fit', fontsize=12)
        axes[1].set_ylim(ylim[0], ylim[1])
        axes[1].set_xlim(xlim[0], xlim[1])
        axes[1].set_title('Gaussian fit')

    for i in range(max_iterations):

        sigma_start = sigma # update starting sigma

        # remove outliers
        sample = np.extract(sample < mu + sigma_level * sigma, sample)  # remove upper tail
        sample = np.extract(sample > mu - sigma_level * sigma, sample)  # remove lower tail

        # re-estimate Gaussian and RMSE
        (mu, sigma) = norm.fit(sample)
        RMSE = '%.3f' % (np.ma.sqrt(np.ma.mean(np.square(sample))))

        if diagnose: # update plot
            hist, bins = np.histogram(sample, range=[lower_bin, upper_bin], density=False, bins=100)
            y = mlab.normpdf(bins, mu, sigma)
            for rect, h in zip(bars, hist):
                rect.set_height(h)
            for line in axes[1].lines:
                line.set_xdata(bins)
                line.set_ydata(y)
            text.set_text('RMSE=' + str(RMSE))
            plt.pause(0.1)

        if sigma_start - sigma < min_change_sigma:

            return sample, mu, sigma, RMSE

    return sample, mu, sigma, RMSE


def compute_rmse_xy(px1, px2, cc_field, cc_thresh=0.3, zoom=0.25, diagnose=True,
                    titles=['Residuals X before correction', 'Residuals Y before correction'],
                    xlim=[-1, 1], ylim=[0, 3]):

    """
    Compute the rmse_xy from the two input offsets fields after masking and Gaussian fit

    Args:
        px1: full path to the displacement in X
        px2: full path to the displacement in Y
        cc_field: full path to the correlation coefficient raster
        cc_thresh: threshold on the correlation coefficient below which observation will not be considered
        zoom: downsampling factor, e.g. 0.25 means computation on a quarter of the original pixel size
        diagnose: plot the Gaussian fitting
        titles: title of the upper figure in the two plots
        xlim: lower limit, upper limit of the x - axis
        ylim: lower limit, upper limit of the y - axis

    Returns: RMSE_XY

    """

    # px1 = '/home/stumpf/Data/CoregisSites/34UFU/correlation/S2A_OPER_MSI_L1C_TL_MPS__20160716T131808_A005567_T34UFU_B03_S2A_OPER_MSI_L1C_TL_MPS__20160726T131755_A005710_T34UFU_B03/MEC/Px1_Num6_DeZoom1_LeChantier.tif'
    # px2 = '/home/stumpf/Data/CoregisSites/34UFU/correlation/S2A_OPER_MSI_L1C_TL_MPS__20160716T131808_A005567_T34UFU_B03_S2A_OPER_MSI_L1C_TL_MPS__20160726T131755_A005710_T34UFU_B03/MEC/Px2_Num6_DeZoom1_LeChantier.tif'
    # cc_field = '/home/stumpf/Data/CoregisSites/34UFU/correlation/S2A_OPER_MSI_L1C_TL_MPS__20160716T131808_A005567_T34UFU_B03_S2A_OPER_MSI_L1C_TL_MPS__20160726T131755_A005710_T34UFU_B03/MEC/Correl_LeChantier_Num_5.tif'
    # cc_thresh = 0.5
    # zoom = 0.25

    import gdal
    import numpy as np
    import cv2
    import misc

    # import offset
    displacement_field_link = gdal.Open(px1)
    px_band = displacement_field_link.GetRasterBand(1)
    px1_surface = np.array(px_band.ReadAsArray())

    displacement_field_link = gdal.Open(px2)
    px_band = displacement_field_link.GetRasterBand(1)
    px2_surface = np.array(px_band.ReadAsArray())

    # downsample
    if not zoom == 1:
        px1_surface = cv2.resize(px1_surface, None, fx=zoom, fy=zoom, interpolation=cv2.INTER_AREA)
        px2_surface = cv2.resize(px2_surface, None, fx=zoom, fy=zoom, interpolation=cv2.INTER_AREA)

    # import correlation coefficient
    cc_link = gdal.Open(cc_field)
    px_band = cc_link.GetRasterBand(1)
    in_cc = np.array(px_band.ReadAsArray())

    # downsample
    if not zoom == 1:
        in_cc = cv2.resize(in_cc, None, fx=zoom, fy=zoom, interpolation=cv2.INTER_AREA)
    # create mask
    in_cc = (in_cc - 127) / 128 < cc_thresh
    # apply mask
    px1_surface = np.ma.masked_array(px1_surface, mask=in_cc)
    px2_surface = np.ma.masked_array(px2_surface, mask=in_cc)

    # filter outliers
    _, mean_gauss, std_gauss, rmse_gauss = misc.iterative_gaussian_fit(px1_surface, diagnose=diagnose, title=titles[0],
                                                                       xlim=xlim, ylim=ylim)
    upper_threshold = mean_gauss + std_gauss * 2.575
    lower_threshold = mean_gauss - std_gauss * 2.575
    px1_surface = np.ma.masked_where(px1_surface >= upper_threshold, px1_surface)  # remove upper tail
    px1_surface = np.ma.masked_where(px1_surface <= lower_threshold, px1_surface)  # remove lower tail

    _, mean_gauss, std_gauss, rmse_gauss = misc.iterative_gaussian_fit(px2_surface, diagnose=diagnose, title=titles[1],
                                                                       xlim=xlim, ylim=ylim)
    upper_threshold =  mean_gauss + std_gauss * 2.575
    lower_threshold =  mean_gauss - std_gauss * 2.575
    px2_surface = np.ma.masked_where(px2_surface >= upper_threshold, px2_surface) # remove upper tail
    px2_surface = np.ma.masked_where(px2_surface <= lower_threshold, px2_surface) # remove lower tail

    rmse_xy = np.ma.sqrt(np.ma.mean(np.square(px1_surface) + np.square(px2_surface)))

    return rmse_xy


def add_georef(in_raster, ref_raster):

    # add extent and projection of the reference raster to the input raster

    # Opens source dataset
    src_ds = gdal.Open(in_raster, 1)
    format = "GTiff"
    driver = gdal.GetDriverByName(format)

    # Open reference dataset
    ref_ds = gdal.Open(ref_raster, 1)
    ref_gt = ref_ds.GetGeoTransform()
    ref_proj = ref_ds.GetProjection()

    # Set location
    src_ds.SetGeoTransform(ref_gt)

    # Set projection
    src_ds.SetProjection(ref_proj)

    # Close file
    src_ds = None

    return


def check_queue(host, user, password):

    """
    Wrapper to SLURMs squeue
    :param host: hostname
    :param user: username
    :param password: user password
    :return: prints current queue
    """

    cmd = 'squeue -u' + user
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.load_system_host_keys()
    ssh.connect(host, username=user, password=password)
    ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command(cmd)
    s = ssh_stdout.read()
    s = s.decode("utf-8")
    print(s)


def launch_remote(master_slave_info, host, user, password, partition, grant, mail, hpc_workfolder, job_number,
                  slurm_template='slurm_template_1'):

    """
    :param master_slave_info: master slave dictionary
    :param host: remote host
    :param user: user name
    :param password: user password
    :param partition: grant partition
    :param grant: grant name
    :param mail: user email
    :param hpc_workfolder: name of the general work folder which will be created as a sub folder in the remote home
    :param job_number: number that will be attached the name of the SLURM script
    :param slurm_template: slurm template to be used
    :return: master slave dictionary including infos on the correlation folder and the job file
    
    remote_home: user home on the remote server, will be queried automatically
    hpc_workfolder: name of the general work folder which will be created as a sub folder in the remote home
    hpc_job_folder: hpc_workfolder + a subfolder for the current job
    remote_correl_folder: remote_home + hpc_job_folder
    
    """

    hpc_job_folder = os.path.join(hpc_workfolder, os.path.basename(master_slave_info['correl_folder']))

    # get home folder
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(host, username=user, password=password)
    stdin, stdout, stderr = client.exec_command('pwd')
    remote_home = stdout.read().rstrip().decode("utf-8")
    client.close()

    # Write slurm job file
    job_file = os.path.join(os.path.dirname(master_slave_info['correl_folder']),
                            'runMicMacJob' + str(job_number) + '.sh')

    # Read in the template
    with open(slurm_template, 'r') as file:
        filedata = file.read()

    # Replace the target strings
    filedata = filedata.replace('yourmail', mail)
    filedata = filedata.replace('yourpartion', partition)
    filedata = filedata.replace('yourgrant', grant)

    # Write the job file
    with open(job_file, 'w') as file:
        file.write(filedata)
        file.write('cd ' + os.path.join(remote_home, hpc_job_folder) + ' \n')
        file.write('MICMAC ' + os.path.basename(master_slave_info['param_xml']) + '\n')

    # create work folder under remote home
    ssh = SSHSession(host, username=user, password=password)
    remote_folder = os.path.join(remote_home, hpc_workfolder)
    try:
        ssh.sftp.stat(remote_folder)
        print('Remote base folder exists ...')
    except IOError:
        print('Remote base folder is beeing created ..')
        ssh.sftp.mkdir(os.path.join(remote_home, hpc_workfolder))

    # copy all files needed to the work folder
    remote_correl_folder = os.path.join(remote_home, hpc_job_folder)
    try:
        ssh.sftp.stat(remote_correl_folder)
        print('Remote work folder exists ...')
    except IOError:
        print('Remote work folder is being created ..')
        ssh.sftp.mkdir(remote_correl_folder)

    print('Copying files to remote host...')
    ssh.sftp.put(master_slave_info['pan_master'],
                 os.path.join(remote_correl_folder, os.path.basename(master_slave_info['pan_master'])))
    ssh.sftp.put(master_slave_info['pan_slave'],
                 os.path.join(remote_correl_folder, os.path.basename(master_slave_info['pan_slave'])))
    ssh.sftp.put(master_slave_info['cloud_mask'],
                 os.path.join(remote_correl_folder, os.path.basename(master_slave_info['cloud_mask'])))
    ssh.sftp.put(master_slave_info['cloud_mask_xml'],
                 os.path.join(remote_correl_folder, os.path.basename(master_slave_info['cloud_mask_xml'])))
    ssh.sftp.put(master_slave_info['param_xml'],
                 os.path.join(remote_correl_folder, os.path.basename(master_slave_info['param_xml'])))
    ssh.sftp.put(job_file, os.path.join(remote_folder, os.path.basename(job_file)))

    master_slave_info['remote_job_folder'] = remote_correl_folder
    master_slave_info['slurm_script'] = job_file

    # launch via sbatch
    cmd = 'cd ' + os.path.join(remote_home, hpc_workfolder) + '\n sbatch ' + os.path.basename(job_file)
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.load_system_host_keys()
    ssh.connect(host, username=user, password=password)
    ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command(cmd)
    stdout = ssh_stdout.read().rstrip().decode("utf-8")
    stderr = ssh_stderr.read().rstrip().decode("utf-8")
    print(stdout)
    if stderr:
        raise Exception(stderr)

    check_queue(host, user, password)

    return master_slave_info


def launch_all_remote(master_slave_info_list,
                      host,
                      user,
                      password,
                      partition,
                      grant,
                      mail,
                      hpc_workfolder,
                      slurm_template='slurm_template_1'):

    """
    Wrapper to launch_remote to submit multiple jobs from a list of master slave dictonaries
    Args:
        master_slave_info_list: list of master slave dictionaries
        host: remote host
        user: user name
        password: user password
        partition: grant partition
        grant: grant name
        mail: user email
        hpc_workfolder: name of the general work folder which will be created as a sub folder in the remote home
        slurm_template: slurm template to be used

    Returns:

    """

    master_slave_info_list_updated =[]
    total_jobs = len(master_slave_info_list)

    for job_number, master_slave_info in enumerate(master_slave_info_list):

        # launch on remote host
        master_slave_info = launch_remote(master_slave_info, host, user, password, partition, grant, mail, hpc_workfolder,
                                               job_number=job_number, slurm_template=slurm_template)
        print('Submitted job number ' + str(job_number+1) + ' out of ' + str(total_jobs))

        master_slave_info_list_updated.append(master_slave_info)

    return master_slave_info_list_updated


def download_result(master_slave_info, host, user, password, last_step=6):
    """
    Wrapper around ssh.get to query and download the results of matching for a master slave pair
    :param master_slave_info: master slave dictionary
    :param host: host name
    :param user: user name
    :param password: user password
    :param last_step: last matching step
    :return: master slave dictionary including paths to px1, px2 and cc plus the downloaded file on disk
    """

    cmd = 'cd ' + master_slave_info['remote_job_folder'] + '\n find . -type f -name ' + '*Px1_Num' + str(
        last_step) + '*.tif' + \
          ' -o -name *Px2_Num' + str(last_step) + '*.tif' + ' -o -name *Correl_LeChantier_Num_' + str(
        last_step - 1) + '*.tif'
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.load_system_host_keys()
    ssh.connect(host, username=user, password=password)
    ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command(cmd)
    s = ssh_stdout.read().decode("utf-8").split('\n')
    s = s[:-1]
    if len(s) == 3:
        print('Found 3 output files...')
        for i in s: print(i)
    else:
        warnings.warn('Could not find three output files on the remote host. Something went wrong :-(')
        return master_slave_info

    ssh = SSHSession(host, username=user, password=password)
    for i, file in enumerate(s):

        if file[0:2] == './':
            file = file[2:]
        remote_file = os.path.join(master_slave_info['remote_job_folder'], file)
        local_file = os.path.join(master_slave_info['correl_folder'], file)
        os.makedirs(os.path.dirname(local_file), exist_ok=True)
        print('Copying... \n', remote_file, '\n to \n', local_file)
        ssh.get(remote_file, local_file)

        if 'Correl' in os.path.basename(local_file):
            master_slave_info['cc'] = local_file
        if 'Px1' in os.path.basename(local_file):
            master_slave_info['px1'] = local_file
        if 'Px2' in os.path.basename(local_file):
            master_slave_info['px2'] = local_file

    return master_slave_info


def download_all_result(master_slave_info_list, host, user, password, last_step=6):

    """
    Wrapper to download all results via download_result according to a list of master slave dictonaries
    Args:
        master_slave_info_list: list of master slave dictionaries
        host: remote host
        user: user name
        password: user password
        last_step: last matching step
    Returns: list of master slave dictionaries including paths to px1, px2 and cc plus the downloaded file on disk

    """

    master_slave_info_list_updated =[]

    for master_slave_info in master_slave_info_list:

        # launch on remote host
        master_slave_info = download_result(master_slave_info, host, user, password, last_step=last_step)

        master_slave_info_list_updated.append(master_slave_info)

    return master_slave_info_list_updated