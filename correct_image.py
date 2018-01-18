import subprocess
import os
import pprint
import gc
import array2geotiff
import osr
import ogr
import misc
from footprints import get_detector_footprint
from osgeo import ogr
from mask import geom2mask
import matplotlib.pyplot as plt
import cv2
import numpy as np
import gdal


def warp(input_grid, warping_parameters, output_file_path, target_srs, res_ref=1,
         resampling_method='cubic'):

    xres_new = warping_parameters[0]
    yres_new = warping_parameters[1]

    # call to gdal
    subprocess.call(['gdalwarp',
                     input_grid,
                     '-te', str(warping_parameters[2]), str(warping_parameters[3]), str(warping_parameters[4]),
                     str(warping_parameters[5]),
                     '-t_srs', target_srs,
                     '-tr', str(xres_new), str(abs(yres_new)),
                     '-overwrite',
                     '-r', resampling_method,
                     output_file_path])

    # adjust displacement values in pixels according to grid resolution
    ds = gdal.Open(output_file_path, gdal.GA_Update)
    band_link = ds.GetRasterBand(1)
    band_array = band_link.ReadAsArray()
    band_array = band_array / (xres_new / res_ref)
    band_link.WriteArray(band_array)
    ds = None


def warp_correctiongrid(px1_correction_grid,
                        px2_correction_grid,
                        band_paths,
                        satellite='L8',
                        resampling_method = 'cubic'):

    """
    warp correction grids to the extents and pixel grids of the input bands
    :param px1_correction_grid: full path to px1 correction grid
    :param px2_correction_grid: full path to px2 correction grid
    :param band_paths: list of the band paths
    :param satellite: of the slave 'S2' for Sentinel-2 and 'L8' for Landsat-8
    :param resampling_method: see http://www.gdal.org/gdalwarp.html
    :return: full paths to warped correction grids as a dictionary, e.g. of the form
    
    {'L8': {15: {'px1_path': '/.../px1_correction_grid_L8_15m_bands.tif',
                   'px2_path': '/.../px2_correction_grid_L8_15m_bands.tif'},
            30: {'px1_path': '/.../px1_correction_grid_L8_30m_bands.tif',
                   'px2_path': '/.../px2_correction_grid_L8_30m_bands.tif'}}}
    """

    # get geoinformation of the correction grid
    ds_ref = gdal.Open(px1_correction_grid)
    cols_ref = ds_ref.RasterXSize
    rows_ref = ds_ref.RasterYSize
    gt_ref = ds_ref.GetGeoTransform()
    prj_ref = ds_ref.GetProjection()
    xres_ref = gt_ref[1]
    yres_ref = gt_ref[5]

    warping_parameters = []

    for raster_file in band_paths:

        print('Reading band ' + raster_file)
        # link data source
        ds = gdal.Open(raster_file)

        # get geotransform
        gt = ds.GetGeoTransform()
        prj = ds.GetProjection()

        # handle cases where slave image is in a different SRS
        if not prj == prj_ref:
            source = osr.SpatialReference()
            source.ImportFromWkt(prj_ref)

            target = osr.SpatialReference()
            target.ImportFromWkt(prj)

            transform = osr.CoordinateTransformation(source, target)

            upper_left_ref = ogr.Geometry(ogr.wkbPoint)
            upper_left_ref.AddPoint(gt_ref[0], gt_ref[3])

            upper_left_ref.Transform(transform)

            xmin_ref = upper_left_ref.GetX()
            ymax_ref = upper_left_ref.GetY()

            target_srs = prj

        else:
            target_srs = prj_ref
            xmin_ref = gt_ref[0]
            ymax_ref = gt_ref[3]

        # get target resolution
        xres = gt[1]
        yres = gt[5]

        # find closest pixel coordinates of the corners of the correction grid in the grid of the slave image
        px_min_x = int(round((xmin_ref - gt[0]) / gt[1])) # left x = metric offset / pixel size
        px_max_y = int(round((ymax_ref - gt[3]) / gt[5])) # upper_y = metric offset / pixel size
        cols = int(round(cols_ref / (gt[1] / gt_ref[1]))) # n pixels x = n pixels x / pixel size ratio
        rows = int(round(rows_ref / (gt[5] / gt_ref[5]))) # n pixels y = n pixels y / pixel size ratio
        px_max_x = px_min_x + cols  # right x
        px_min_y = px_max_y + rows  # bottom y

        # translate pixel coordinates to map coordinates
        xmin_new, ymax_new = misc.pixel2geo(px_min_x, px_max_y, gt)
        xmax_new, ymin_new = misc.pixel2geo(px_max_x, px_min_y, gt)

        warping_parameters.append([xres, yres, xmin_new, ymin_new, xmax_new, ymax_new])

        ds = None

    # remove duplicates
    warping_parameters = np.array(warping_parameters)
    warping_parameters = np.vstack({tuple(row) for row in warping_parameters})

    # generate output dictionary
    out_files = {}

    if satellite == 'L8':

        # note in the output dictionary
        out_files['L8'] = {}

        # get warping parameters for 30 m bands
        warp_param = [row for row in warping_parameters if 30 in row[0:1]]

        if not warp_param:
            raise Exception('Could not find warping parameters for Landsat-8 30 m bands')
        else: warp_param = warp_param[0]

        # generate and store output path in the output dictionary
        px1_30m = os.path.splitext(px1_correction_grid)[0] + '_' + satellite + '_30m_bands.tif'
        px2_30m = os.path.splitext(px2_correction_grid)[0] + '_' + satellite + '_30m_bands.tif'
        out_files['L8'][30] = {}
        out_files['L8'][30]['px1_path'] = px1_30m
        out_files['L8'][30]['px2_path'] = px2_30m

        # warp px1
        warp(px1_correction_grid, warp_param, px1_30m, target_srs, res_ref=xres_ref, resampling_method=resampling_method)

        # warp px2
        warp(px2_correction_grid, warp_param, px2_30m, target_srs, res_ref=xres_ref, resampling_method=resampling_method)

       # get warping parameters for 15 m bands
        warp_param = [row for row in warping_parameters if 15 in row[0:1]]
        if not warp_param:
            raise Exception('Could not find warping parameters for Landsat-8 15 m bands')
        else: warp_param = warp_param[0]

        # generate and store output path in the output dictionary
        px1_15m = os.path.splitext(px1_correction_grid)[0] + '_' + satellite + '_15m_bands.tif'
        px2_15m = os.path.splitext(px2_correction_grid)[0] + '_' + satellite + '_15m_bands.tif'
        out_files['L8'][15] = {}
        out_files['L8'][15]['px1_path'] = px1_15m
        out_files['L8'][15]['px2_path'] = px2_15m

        # warp px1
        warp(px1_correction_grid, warp_param, px1_15m, target_srs, res_ref=xres_ref, resampling_method=resampling_method)

        # warp px2
        warp(px2_correction_grid, warp_param, px2_15m, target_srs, res_ref=xres_ref, resampling_method=resampling_method)

    elif satellite == "S2":

        out_files['S2'] = {}

        # get warping parameters for 10 m bands
        warp_param = [row for row in warping_parameters if 10 in row[0:1]]
        if not warp_param:
            raise Exception('Could not find warping parameters for Sentinel-2 10 m bands')
        else:
            warp_param = warp_param[0]

        # check if input geotransform is identical to destination geotransform
        if gt_ref == (warp_param[2], warp_param[0], 0.0, warp_param[5], 0.0,  warp_param[1]):
            print('Input and output transformation are identical. Keeping original grid.')
            px1_10m = px1_correction_grid
            px2_10m = px2_correction_grid
        else:
            # generate output paths
            px1_10m = os.path.splitext(px1_correction_grid)[0] + '_' + satellite + '_10m_bands.tif'
            px2_10m = os.path.splitext(px2_correction_grid)[0] + '_' + satellite + '_10m_bands.tif'

            # warp px1
            warp(px1_correction_grid, warp_param, px1_10m, target_srs, res_ref=xres_ref, resampling_method=resampling_method)
            # warp px2
            warp(px2_correction_grid, warp_param, px2_10m, target_srs, res_ref=xres_ref, resampling_method=resampling_method)

        # store output paths in the output dictionary
        out_files['S2'][10] = {}
        out_files['S2'][10]['px1_path'] = px1_10m
        out_files['S2'][10]['px2_path'] = px2_10m

        # get warping parameters for 20 m bands
        warp_param = [row for row in warping_parameters if 20 in row[0:1]]
        if not warp_param:
            raise Exception('Could not find warping parameters for Sentinel-2 20 m bands')
        else:
            warp_param = warp_param[0]

        # check if input geotransform is identical to destination geotransform
        if gt_ref == (warp_param[2], warp_param[0], 0.0, warp_param[5], 0.0, warp_param[1]):
            print('Input and output transformation are identical. Keeping original grid.')
            px1_20m = px1_correction_grid
            px2_20m = px2_correction_grid
        else:
            # generate output paths
            px1_20m = os.path.splitext(px1_correction_grid)[0] + '_' + satellite + '_20m_bands.tif'
            px2_20m = os.path.splitext(px2_correction_grid)[0] + '_' + satellite + '_20m_bands.tif'

            # warp px1
            warp(px1_correction_grid, warp_param, px1_20m, target_srs, res_ref=xres_ref, resampling_method=resampling_method)
            # warp px2
            warp(px2_correction_grid, warp_param, px2_20m, target_srs, res_ref=xres_ref, resampling_method=resampling_method)

        # store output paths in the output dictionary
        out_files['S2'][20] = {}
        out_files['S2'][20]['px1_path'] = px1_20m
        out_files['S2'][20]['px2_path'] = px2_20m

        # get warping parameters for 60 m bands
        warp_param = [row for row in warping_parameters if 60 in row[0:1]]
        if not warp_param:
            raise Exception('Could not find warping parameters for Sentinel-2 60 m bands')
        else:
            warp_param = warp_param[0]

        # check if input geotransform is identical to destination geotransform
        if gt_ref == (warp_param[2], warp_param[0], 0.0, warp_param[5], 0.0, warp_param[1]):
            print('Input and output transformation are identical. Keeping original grid.')
            px1_60m = px1_correction_grid
            px2_60m = px2_correction_grid
        else:
            # generate output paths
            px1_60m = os.path.splitext(px1_correction_grid)[0] + '_' + satellite + '_60m_bands.tif'
            px2_60m = os.path.splitext(px2_correction_grid)[0] + '_' + satellite + '_60m_bands.tif'

            # warp px1
            warp(px1_correction_grid, warp_param, px1_60m, target_srs, res_ref=xres_ref, resampling_method=resampling_method)

            # warp px2
            warp(px2_correction_grid, warp_param, px2_60m, target_srs, res_ref=xres_ref, resampling_method=resampling_method)

        # store output paths in the output dictionary
        out_files['S2'][60] = {}
        out_files['S2'][60]['px1_path'] = px1_60m
        out_files['S2'][60]['px2_path'] = px2_60m

    else:
        raise Exception(satellite + ' is not a suported option. Choose S2 of L8.')

    return out_files


def apply_correction(band_paths, corr_grid_dict, satellite):

    """
    A function to apply the correction grids obtained from warp_correctiongrid to all bands (on at a time)
    of the slave input image.
    :param band_paths: list of the slave band paths
    :param corr_grid_dict: dictionary with the paths of the warping grids as returned by warp_correctiongrid
    :param satellite: satellite name of the slave being 'S2' or 'L8'
    :return: list to the full paths of the corrected bands
    """

    # assemble dictionary with band infos
    band_info = {}
    band_info[satellite] = {}
    for index, band in enumerate(band_paths):

        ds = gdal.Open(band)
        gt = ds.GetGeoTransform()
        res = set([abs(gt[i]) for i in [1, 5]])
        if len(res) > 1:
            raise Exception("None squared pixels are not supported...")
        else: res = int(list(res)[0])

        # band_info[sat][res] = {}
        if res not in band_info[satellite]:
            band_info[satellite][res] = {}
        band_info[satellite][res][index] = {}
        band_info[satellite][res][index] = band


    print('Found the following bands ...')
    pprint.pprint(band_info)

    print('Using the following correction grids ...')
    pprint.pprint(corr_grid_dict)

    corrected_files = []

    # loop over all different resolutions
    res_keys = list(band_info[satellite].keys())
    for res in res_keys:

        # get all bands at this resolution
        current_band_paths = list(band_info[satellite][res].values())

        # get corresponding correction grids
        corr_grid_path_px1 = corr_grid_dict[satellite][res]['px1_path']
        corr_grid_path_px2 = corr_grid_dict[satellite][res]['px2_path']

        ds_px1 = gdal.Open(corr_grid_path_px1)
        ds_px2 = gdal.Open(corr_grid_path_px2)
        gt_corr_grid = ds_px1.GetGeoTransform()

        # get common bounding box at this resolution
        bb = misc.find_minimum_bounding_box(corr_grid_dict[satellite][res]['px1_path'], current_band_paths[0])
        ulX, ulY = misc.world_to_pixel(gt_corr_grid, bb[0], bb[3])
        lrX, lrY = misc.world_to_pixel(gt_corr_grid, bb[2], bb[1])
        cols = (lrX - ulX)
        rows = (lrY - ulY)

        # read
        band_link_px1 = ds_px1.GetRasterBand(1)
        band_link_px2 = ds_px2.GetRasterBand(1)
        px1_corr_grid = band_link_px1.ReadAsArray(ulX, ulY, cols, rows)
        px2_corr_grid = band_link_px2.ReadAsArray(ulX, ulY, cols, rows)

        ds_px1 = None
        ds_px2 = None

        # correction grids to cv2 map grids taken into account the image dimensions
        ny, nx = px1_corr_grid.shape
        x = np.arange(0, nx)
        y = np.arange(0, ny)
        xv, yv = np.meshgrid(x, y)
        px1_corr_grid += xv
        px2_corr_grid += yv
        px1_corr_grid = px1_corr_grid.astype('float32')
        px2_corr_grid = px2_corr_grid.astype('float32')
        del nx, ny, x, y, xv, yv
        gc.collect()

        # iterate over all corresponding bands
        for index, band in enumerate(current_band_paths):

            print('Reading band ' + band)
            ds = gdal.Open(band)

            # figure out subset in the slave image
            gt = ds.GetGeoTransform()
            prj = ds.GetProjection()
            px_min_x, px_max_y = misc.world_to_pixel(gt, bb[0], bb[3])

            if not px_min_x == int(px_min_x):
                raise Exception('Computed output resolution is not an integer, '
                                'indicates alignment problem between image and correction grids')
            if not px_max_y == int(px_max_y):
                raise Exception('Computed output resolution is not aninteger, '
                                'indicates alignment problem between image and correction grids')

            # get subset of band
            band_link = ds.GetRasterBand(1)
            band_array = band_link.ReadAsArray(px_min_x, px_max_y, cols, rows)
            # band_array = band_array.astype('float32')

            if band.endswith('BQA.TIF'):
                # Landsat-8 quality band should be resampled with nearest neighbour
                band_corrected_array = cv2.remap(band_array,
                                                 px1_corr_grid,
                                                 px2_corr_grid,
                                                 interpolation=cv2.INTER_NEAREST)

            else:
                # all others should be resampled with cubic interpolation
                band_corrected_array = cv2.remap(band_array,
                                                 px1_corr_grid,
                                                 px2_corr_grid,
                                                 interpolation=cv2.INTER_CUBIC)

            # save corrected band to disk
            band_corrected_array = np.rint(band_corrected_array)
            dst_filename = os.path.join(os.path.dirname(corr_grid_path_px1),
                                        os.path.splitext(os.path.basename(band))[0] + '_corrected.tif')
            try:
                os.remove(dst_filename)
            except OSError:
                pass
            print('Writing ' + dst_filename)
            out_gt = (bb[0], gt[1], gt[2], bb[3], gt[4], gt[5])
            array2geotiff.array2geotiff(band_corrected_array, dst_filename, (out_gt, prj))

            corrected_files.append(dst_filename)

    return corrected_files


def correct_image(slave_image_path,
                  slave_gml_path,
                  master_path,
                  correction_surface_x,
                  correction_surface_y,
                  diagnose=False):

    """
    A function to apply the correction grids obtained from refineCoregistration() to all bands (on at a time)
     of the slave input image (second date).

    Args:
        slave_image_path (string): full path to the slave image (band) that should be corrected
        slave_gml_path (string): full path to the corresponding GML holding the detector footprint
        master_path (string): full path to the master image (band) with which the slave was correlated
        correction_surface_x (numpy array): correction surface for the x coordinates
        correction_surface_y (numpy array): correction surface for the y coordinates
        diagnose (bool): flag to enable plotting of some diagnostic plots

    Result:
        slave_corrected (numpy array): band corrected for offsets in x and y

    """

    # make a copy to avoid that arrays get changed outside the function as well
    corr_surface_x = np.ma.copy(correction_surface_x)
    corr_surface_y = np.ma.copy(correction_surface_y)

    # get georeference infos for the slave image
    slave_link = gdal.Open(slave_image_path)
    slave_band = slave_link.GetRasterBand(1)
    slave_array = np.array(slave_band.ReadAsArray())
    slave_projection = slave_link.GetProjection()
    slave_georeference = slave_link.GetGeoTransform()
    max_x_slave = slave_georeference[0] + slave_georeference[1] * slave_link.RasterXSize
    min_y_slave = slave_georeference[3] + slave_georeference[5] * slave_link.RasterYSize

    # get georeference infos for the master image
    master_link = gdal.Open(master_path)
    master_projection = master_link.GetProjection()
    master_georeference = master_link.GetGeoTransform()
    max_x_master = master_georeference[0] + master_georeference[1] * master_link.RasterXSize
    min_y_master = master_georeference[3] + master_georeference[5] * master_link.RasterYSize

    # consistency check
    if not slave_projection == master_projection:
        raise Exception('Projection of the slave and master images are not identical')

    if not (slave_georeference[0], max_x_slave, min_y_slave, slave_georeference[3]) == (master_georeference[0],
                                                                                        max_x_master,
                                                                                        min_y_master,
                                                                                        master_georeference[3]):
        raise Exception('Extent of the slave and master images are not identical.')

    # if input band has a resolution that differs from the master resample the correction grid and its mask
    if not corr_surface_x.shape == (slave_link.RasterXSize, slave_link.RasterYSize):
        print("Correction surface and slave image size differ. The correction surface is resampled to fit...")

        if hasattr(corr_surface_x, 'mask'):
            mask = cv2.resize(corr_surface_x.mask.astype(float),
                              slave_array.shape,
                              interpolation=cv2.INTER_AREA)
            mask = mask > 0.5

        corr_surface_x = cv2.resize(corr_surface_x,
                                    slave_array.shape,
                                    interpolation=cv2.INTER_AREA)

        corr_surface_y = cv2.resize(corr_surface_y,
                                    slave_array.shape,
                                    interpolation=cv2.INTER_AREA)

        if 'mask' in locals():
            corr_surface_x = np.ma.masked_array(corr_surface_x, mask=mask)
            corr_surface_y = np.ma.masked_array(corr_surface_y, mask=mask)

    # check if the extent of the input band is fully covered by the corrections surface
    # get extent of the area covered in the slave image from the detector footprints
    slave_detector_footrprint = get_detector_footprint(slave_gml_path, diagn=diagnose)

    union = None
    for i in range(len(slave_detector_footrprint) - 1):

        if union is None:

            ' '.join([str(j) for j in slave_detector_footrprint[i]])

            wkt1 = "POLYGON((" + ', '.join([str(j).strip("[]") for j in slave_detector_footrprint[i]]) + "))"
            wkt2 = "POLYGON((" + ', '.join([str(j).strip("[]") for j in slave_detector_footrprint[i + 1]]) + "))"
            poly1 = ogr.CreateGeometryFromWkt(wkt1)
            poly2 = ogr.CreateGeometryFromWkt(wkt2)

            union = poly1.Union(poly2)

        else:
            wkt2 = "POLYGON((" + ', '.join([str(j).strip("[]") for j in slave_detector_footrprint[i + 1]]) + "))"
            poly2 = ogr.CreateGeometryFromWkt(wkt2)
            union = union.Union(poly2)

    slave_mask, ulX, ulY, gt2 = geom2mask(slave_array, union, gt=slave_georeference)

    mask_difference = slave_mask - corr_surface_x.mask

    if diagnose:
        fig, axes = plt.subplots(ncols=3)
        axes[0].imshow(slave_mask)
        axes[0].set_title('Extent slave image', fontsize=14)
        axes[0].tick_params(labelsize=14)
        axes[1].imshow(corr_surface_x.mask)
        axes[1].set_title('Extent correction surface', fontsize=14)
        axes[1].tick_params(labelsize=14)
        axes[2].imshow(mask_difference)
        axes[2].set_title('Overlap: -1 -> not covered', fontsize=14)
        axes[2].tick_params(labelsize=14)
        plt.show()

    if (mask_difference == -1).any():
        raise Exception('Correction surface does not fully cover the extent of the input image. Not permitted!')

    # convert correction surface to format required by cv2.remap
    nx, ny = corr_surface_x.shape
    x = np.arange(0, nx)
    y = np.arange(0, ny)
    xv, yv = np.meshgrid(x, y)
    corr_surface_x += xv
    corr_surface_x = corr_surface_x.astype(np.float32)
    corr_surface_y += yv
    corr_surface_y = corr_surface_y.astype(np.float32)

    slave_corrected = cv2.remap(slave_array, corr_surface_x, corr_surface_y, interpolation=cv2.INTER_CUBIC)

    return slave_corrected