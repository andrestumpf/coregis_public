import subprocess
import os
import gdal
from array2geotiff import array2geotiff
import numpy as np
from osgeo import gdalnumeric
from PIL import Image, ImageDraw
import misc
import warnings
import itertools
import coregis_filter
import ogr


def create_cloud_mask(granule_info,
                      fmask_path,
                      pixel_size_x=20,
                      pixel_size_y=20,
                      output_folder=None,
                      overwrite=False):

    """
    Wrapper to FMASK for Landsat-8 and Sentinel-2 cloud masking
    
    :param granule_info: dictionary with granule info
    :param fmask_path: path to the python-fmask folder
    :param pixel_size_x: pixel size for the cloud mask in y
    :param pixel_size_y: pixel size for the cloud mask in y
    :param output_folder: full path to the output folder
    :param overwrite: if True existing cloud masks will be overwritten
    :return: dictionary with the full path to cloud mask image on disk, if output_folder is provided the files will be
             written to a sub folder named after the granule, otherwise written to the input folder
    
    # DEBUG
    granule_info = granule_info_master
    
    """

    # figure out which script to use
    if granule_info['satellite'] == 'S2':

        # create output folder
        if output_folder:
            granule_folder = os.path.join(output_folder, os.path.basename(granule_info['granule_path'])[:-7])
            os.makedirs(granule_folder, exist_ok=True)
        else:
            granule_folder = os.path.dirname(granule_info['band_paths'][0])

        # check if mask already exists
        cloud_mask_path = os.path.join(granule_folder, 'cloud_mask.img')
        if os.path.exists(cloud_mask_path):
            print('Cloud mask for this datasets alredy exists in ' + cloud_mask_path)
            if not overwrite:
                print('Existing cloud mask will be used.')
                granule_info['cloud_mask'] = cloud_mask_path
                return granule_info
            else:
                print('Cloud mask will be recomputed...')

        print('Using Fmask pipeline for Sentinel-2 imagery...')
        vrt_path = misc.granule2vrt(granule_info,
                                    pixel_size_x=pixel_size_x,
                                    pixel_size_y=pixel_size_y,
                                    output_folder=granule_folder)

        # generate sun angle grids
        angle_path = os.path.join(granule_folder, 'sun_angle.img')
        cmd = ['python3', os.path.join(fmask_path, 'fmask_sentinel2makeAnglesImage.py'), '-i']
        cmd.append(granule_info['metadata'])
        cmd = cmd + ['-o', angle_path]
        subprocess.call(cmd)

        # create the cloud mask output image
        cmd = ['python3', os.path.join(fmask_path, 'fmask_sentinel2Stacked.py'),
               '-a', vrt_path,
               '-z', angle_path,
               '-o', cloud_mask_path]

        print('Generating cloud mask...')
        subprocess.call(cmd)

    elif granule_info['satellite'] == 'L8':

        print('Using Fmask pipeline for Landsat-8 imagery...')

        # create output folder
        if output_folder:
            granule_folder = os.path.join(output_folder, os.path.basename(granule_info['granule_path']))
            os.makedirs(granule_folder, exist_ok=True)
        else:
            granule_folder = granule_info['granule_path']

        # check if mask already exists
        cloud_mask_path = os.path.join(granule_folder, 'cloud_mask.img')
        if os.path.exists(cloud_mask_path):
            print('Cloud mask for this datasets alredy exists in ' + cloud_mask_path)
            if not overwrite:
                print('Existing cloud mask will be used.')
                granule_info['cloud_mask'] = cloud_mask_path
                return granule_info
            else:
                print('Cloud mask will be recomputed...')

        # concatenate Landsat-8 multi-band image
        ref_img_path = os.path.join(granule_folder, 'ref.img')
        bands_paths = granule_info['band_paths']
        print('Concatenating Landsat 8 image...')
        subprocess.call(['gdal_merge.py',
                         '-separate',
                         '-of', 'HFA',
                         '-co', 'COMPRESSED=YES',
                         '-o', ref_img_path,
                         bands_paths[0],
                         bands_paths[1],
                         bands_paths[2],
                         bands_paths[3],
                         bands_paths[4],
                         bands_paths[5],
                         bands_paths[6],
                         bands_paths[8]
                        ])

        # concatenate Landsat-8 thermal multi-band image
        therm_img_path = os.path.join(granule_folder, 'thermal.img')
        print('Concatenating Landsat 8 image thermal bands...')
        subprocess.call(['gdal_merge.py',
                         '-separate',
                         '-of', 'HFA',
                         '-co', 'COMPRESSED=YES',
                         '-o', therm_img_path,
                         bands_paths[9],
                         bands_paths[10]
                         ])

        # generate sun angle image
        angles_path = os.path.join(granule_folder, 'sun_angle.img')
        meta_data_path = granule_info['metadata']
        print('Computing image of incidence angles...')
        print(meta_data_path)
        subprocess.call(['python3', os.path.join(fmask_path, 'fmask_usgsLandsatMakeAnglesImage.py'),
                         '-m', meta_data_path,
                         '-t', ref_img_path,
                         '-o', angles_path])

        # generate saturation mask
        saturationmask_path = os.path.join(granule_folder, 'saturationmask.img')
        print('Computing saturation mask image...')
        subprocess.call(['python3', os.path.join(fmask_path, 'fmask_usgsLandsatSaturationMask.py'),
                         '-i', ref_img_path,
                         '-m', meta_data_path,
                         '-o', saturationmask_path])

        # compute TOA
        toa_path = os.path.join(granule_folder, 'toa.img')
        print('Computing top-of-atmosphere reflectance...')
        subprocess.call(['python3', os.path.join(fmask_path, 'fmask_usgsLandsatTOA.py'),
                         '-i', ref_img_path,
                         '-m', meta_data_path,
                         '-z', angles_path,
                         '-o', toa_path])

        # compute cloud mask
        cloud_mask_path = os.path.join(granule_folder, 'cloud_mask.img')
        print('Computing cloud mask...')
        subprocess.call(['python3', os.path.join(fmask_path, 'fmask_usgsLandsatStacked.py'),
                         '-t', therm_img_path,
                         '-a', toa_path,
                         '-m', meta_data_path,
                         '-z', angles_path,
                         '-s', saturationmask_path,
                         '-o', cloud_mask_path])

    else:
        raise Exception(granule_info['satellite'] + ' is not a suported option. Choose S2 of L8.')

    granule_info['cloud_mask'] = cloud_mask_path

    return granule_info


def create_all_cloud_masks(granule_info_list,
                           fmask_path,
                           pixel_size_x=20,
                           pixel_size_y=20,
                           output_folder=None,
                           overwrite=False):
    """
    Wrapper to create_cloud_mask to iterate over a list of granule info dictionaries
    Args:
        granule_info_list: list of granule info dictionaries as returned by get_granule_info
        fmask_path: full path to the python-fmask folder
        pixel_size_x:  pixel size for the cloud mask in y
        pixel_size_y: pixel size for the cloud mask in y
        output_folder: full path to the output folder
        overwrite: if True existing cloud masks will be overwritten

    Returns: list of dictionaries with the full path to cloud mask image on disk, if output_folder is provided the files
             will be written to a sub folder named after the granule, otherwise written to the input folder

    """

    # generate cloud masks
    granule_info_list_updated = []
    for granule_info in granule_info_list:

        # generate cloud masks
        granule_info_list_updated.append(create_cloud_mask(granule_info,
                                                           fmask_path,
                                                           pixel_size_x=pixel_size_x,
                                                           pixel_size_y=pixel_size_y,
                                                           output_folder=output_folder,
                                                           overwrite=overwrite))
    return granule_info_list_updated


def create_cloud_mask_by_threshold(band_path,
                                   pixel_size_x=20,
                                   pixel_size_y=20,
                                   threshold=2000):
    """
    Wrapper to FMASK for Landsat-8 and Sentinel-2 cloud masking

    Args:
        band_path: full path to the input tile
        fmask_path:  path to the python-fmask folder
        pixel_size_x: pixel size for the cloud mask in y
        pixel_size_y: pixel size for the cloud mask in y
        threshold: all pixels exceeding this value will be marked as masked

    Returns: Full path to cloud mask image on disk

    """
    # read input
    ds = gdal.Open(band_path)
    gt = ds.GetGeoTransform()
    band = ds.GetRasterBand(1)
    band_array = band.ReadAsArray()

    # threshold
    band_array[band_array <= threshold] = 1
    band_array[band_array > threshold] = 2

    # write to disk
    dst_filename = os.path.join(os.path.dirname(band_path), 'cloud_mask.tif')
    array2geotiff(band_array, dst_filename, band_path, out_datatype=None)

    # check if input pixel size deviates from output pixel size
    if pixel_size_x is not gt[1] or pixel_size_y is not abs(gt[5]):

        out_filename = os.path.splitext(dst_filename)[0] + '.img'

        # reopen
        ds_1 = gdal.Open(dst_filename)
        gt1 = ds_1.GetGeoTransform()

        # figure projection paramters
        min_x = gt1[0]
        max_x = gt1[0] + (gt1[1] * ds_1.RasterXSize)
        min_y = gt1[3] + (gt1[5] * ds_1.RasterYSize)
        max_y = gt1[3]

        # warp
        subprocess.call(['gdalwarp',
                         dst_filename,
                         '-of', 'HFA',
                         '-te', str(min_x), str(min_y), str(max_x), str(max_y),
                         '-tr', str(pixel_size_x), str(pixel_size_y),
                         '-overwrite',
                         out_filename])

        # remove intermediate file
        os.remove(dst_filename)

        return out_filename

    else:
        # return intermediate file
        return dst_filename


def check_cloud_mask(granule_info, cloud_values=[2, 4, 5], maximum_covered=0.7):
    """
    Checks if the cloud mask indicates less coverage than the defined fraction
    Args:
        granule_info: granule info dictionary
        cloud_values: values for which recoding returns 0 in the mask, 0=null, 1=clear, 2=cloud, 3=shadow, 4=snow,
                         5=water
        maximum_covered: maximum allowed fraction of 0s

    Returns: False or True according to if the test passed or not

    """

    ds_cm = gdal.Open(granule_info['cloud_mask'])
    band_cm = ds_cm.GetRasterBand(1)
    band_array_cm = band_cm.ReadAsArray()
    classes = np.unique(band_array_cm)
    ds_cm = None

    none_zeros = np.count_nonzero(band_array_cm)

    for class_i in classes:
        # print(class_i)
        if class_i in cloud_values:
            band_array_cm[band_array_cm == class_i] = 0
        else:
            band_array_cm[band_array_cm == class_i] = 1

    zeros = np.count_nonzero(band_array_cm == 0)
    if zeros / none_zeros > maximum_covered:
        warnings.warn('Cloud coverage fraction of ' + str(zeros / none_zeros) + ' exceeds the maximum fraction of ' +
                      str(maximum_covered))
        warnings.warn('The dataset +' + os.path.dirname(granule_info['band_paths'][0]) + ' will not be processed...')
        return False
    else:
        print('Cloud coverage fraction for ' + os.path.dirname(granule_info['band_paths'][0]) + 'is ' + str(
            zeros / none_zeros))
        print('Passed cloud mask check!')
        return True


def filter_by_cloud_masks(granule_info_list, cloud_values=[2, 4, 5], maximum_covered=0.7):
    """
    Filter out scenes with a too large cloud coverage
    Args:
        granule_info_list: list of granule info dictionaries
        cloud_values: values for which recoding returns 0 in the mask, 0=null, 1=clear, 2=cloud, 3=shadow, 4=snow,
                         5=water
        maximum_covered: maximum allowed fraction of 0s

    Returns: list of filterd granule info dictionaries

    """

    passed = []
    for granule_info in granule_info_list:
        passed.append(check_cloud_mask(granule_info, cloud_values=cloud_values, maximum_covered=maximum_covered))

    granule_info_list = list(itertools.compress(granule_info_list, passed))

    return granule_info_list


def recode_mask(cloud_mask_path):
    """
    Recode fmask cloud mask to binary raster where only clouds are marked with 1
    Args:
        cloud_mask_path: full path to the input cloud mask

    Returns: full path to the binary output cloud mask

    """

    input_image_link = gdal.Open(cloud_mask_path)
    px_band = input_image_link.GetRasterBand(1)
    in_mask = np.array(px_band.ReadAsArray())

    out_path_woext = os.path.splitext(cloud_mask_path)[0]
    out_path = out_path_woext + '_bin.tif'

    in_mask[in_mask == 1] = 0
    in_mask[in_mask == 2] = 1
    in_mask[in_mask == 3] = 0
    in_mask[in_mask == 4] = 1
    in_mask[in_mask == 5] = 1

    array2geotiff(in_mask, out_path, cloud_mask_path)

    input_image_link = None

    return out_path


def generate_all_combined_masks(master_slave_info_list,
                                cloud_values=[0, 2, 4, 5],
                                close_radius=100,
                                open_radius=200,
                                otbpath=None,
                                output_folder=None):
    """
    Generate combined cloud masks for all master slave pairs in the list including filtering and alignement
    Args:
        master_slave_info_list: list of master slave dictionaries
        cloud_values: values for which recoding returns 0 in the mask, 0=null, 1=clear, 2=cloud, 3=shadow, 4=snow,
                         5=water
        close_radius:  radius of the structuring element (ball) for closing in meters
        open_radius: radius of the structuring element (ball) for open in meters
        otbpath: path to the OTB binaries including otbcli_BinaryMorphologicalOperation
        output_folder:  folder to which the output rasters will be written, if not specified it defaults to
                        os.path.dirname(master_slave_info['pan_master'])
    Returns: list of master slave dictionaries including the paths to the combined cloud masks

    """

    master_slave_info_list_updated = []
    for master_slave_info in master_slave_info_list:

        if not output_folder:
            out_path = os.path.dirname(master_slave_info['pan_master'])
        else:
            out_path = output_folder

        # combine cloud masks
        master_slave_info = combine_masks(master_slave_info, cloud_values=cloud_values,
                                               output_folder=out_path)

        # morphological filter on combined cloud mask
        # close areas which will be correlated, than open (i.e. delete) areas which are smaller than open_radius
        master_slave_info['cloud_mask'] = coregis_filter.morphological_filter(master_slave_info['cloud_mask'],
                                                                              close_radius=close_radius,
                                                                              open_radius=open_radius,
                                                                              otbpath=otbpath,
                                                                              foreground_value=255,
                                                                              background_value=0,
                                                                              output_folder=os.path.dirname(
                                                                          master_slave_info['cloud_mask']))

        # upsample mask to full resolution and align with master
        master_slave_info['cloud_mask'] = misc.align_first2second(master_slave_info['cloud_mask'],
                                                                  master_slave_info['pan_master'],
                                                                  resampling_method='near',
                                                                  output_folder=os.path.dirname(
                                                                      master_slave_info['cloud_mask']))

        master_slave_info_list_updated.append(master_slave_info)

    return master_slave_info_list_updated


def check_no_data(raster_1_path, raster_2_path, no_data_value=0):

    """
    Function to generate a mask of no data areas in the two input grids
    Args:
        raster_1_path: full path to the first input image
        raster_2_path: full path to the second input image
        no_data_value: value signaling no data, e.g. for S-2

    Returns: boolean numpy array with True where at least one of the two raster contains a no data value

    """
    # DEBUG
    # raster_1_path = pan_1_path
    # raster_2_path = pan_2_aligned_path

    ds_1 = gdal.Open(raster_1_path)
    band_link_1 = ds_1.GetRasterBand(1)
    band_array_1 = band_link_1.ReadAsArray()

    ds_2 = gdal.Open(raster_2_path)
    band_link_2 = ds_2.GetRasterBand(1)
    band_array_2 = band_link_2.ReadAsArray()

    band_array_1 = band_array_1 == no_data_value
    band_array_2 = band_array_2 == no_data_value

    if not np.any(band_array_1) and not np.any(band_array_1):
        print('Data does not contain no data areas...')
        no_data_mask =  None

    else:
        no_data_mask = np.logical_or(band_array_1, band_array_2)

    return no_data_mask


def combine_masks(master_slave_info,
                  cloud_values=[0, 2, 4, 5],
                  output_folder=None):
    """
    Combine cloud mask from an image pair
    :param master_slave_info: master slave dictionary as returned by make_pair > generate_pan
    :param cloud_values: values for which recoding returns 0 in the mask, 0=null, 1=clear, 2=cloud, 3=shadow, 4=snow, 
                         5=water
    :param output_folder: path to the output folder
    :return: 
    """

    cloud_mask_master = master_slave_info['master']['cloud_mask']
    cloud_mask_slave = master_slave_info['slave']['cloud_mask']

    # open first and get geoinfos
    ds_master = gdal.Open(cloud_mask_master)
    gt_master = ds_master.GetGeoTransform()
    proj_master = ds_master.GetProjection()

    # open second and get geoinfos
    ds_slave = gdal.Open(cloud_mask_slave)
    gt_slave = ds_slave.GetGeoTransform()
    proj_slave = ds_slave.GetProjection()

    # check if in the same reference frame and warp if necessary
    if gt_slave == gt_master and proj_slave == proj_master:
        print('First tile is already in the master geometry. No warping performed...')
        cloud_mask_slave_prj = cloud_mask_slave

    else:

        print('First tile is beeing reprojected to the geometry of the master...')
        cloud_mask_slave_prj = os.path.splitext(cloud_mask_slave)[0] + '_reproj.tif'

        min_x = gt_master[0]
        max_x = gt_master[0] + (gt_master[1] * ds_master.RasterXSize)
        min_y = gt_master[3] + (gt_master[5] * ds_master.RasterYSize)
        max_y = gt_master[3]

        subprocess.call(['gdalwarp',
                         cloud_mask_slave,
                         '-t_srs', proj_master,
                         '-te', str(min_x), str(min_y), str(max_x), str(max_y),
                         '-tr', str(abs(gt_master[1])), str(abs(gt_master[5])),
                         '-overwrite',
                         cloud_mask_slave_prj])

    # open raster grids and compute union
    print('Combining cloud masks ' + cloud_mask_master + ' and ' + cloud_mask_slave_prj)
    cm1_ds = gdal.Open(cloud_mask_slave_prj)
    cm1_band = cm1_ds.GetRasterBand(1)
    cm1_array = cm1_band.ReadAsArray()

    cm2_ds = gdal.Open(cloud_mask_master)
    cm2_band = cm2_ds.GetRasterBand(1)
    cm2_array = cm2_band.ReadAsArray()

    classes_cm1 = np.unique(cm1_array)
    classes_cm2 = np.unique(cm2_array)
    classes = np.unique(np.concatenate([classes_cm1, classes_cm2]))

    for class_i in classes:
        # print(class_i)
        if class_i in cloud_values:
            cm1_array[cm1_array == class_i] = 0
            cm2_array[cm2_array == class_i] = 0
        else:
            cm1_array[cm1_array == class_i] = 1
            cm2_array[cm2_array == class_i] = 1

    combined_mask = cm1_array + cm2_array
    combined_mask[combined_mask > 1] = 255
    combined_mask[combined_mask == 1] = 0

    if output_folder:
        mask_name = 'cloud_masks_assembled_' + master_slave_info['master']['date'] + \
                    master_slave_info['master']['time'] + '_' + master_slave_info['slave']['date'] + \
                    master_slave_info['slave']['time'] + '.tif'
        dst_filename = os.path.join(output_folder, mask_name)
    else:
        dst_filename = os.path.splitext(cloud_mask_master)[0] + '_assembled.tif'

    array2geotiff(combined_mask, dst_filename, cloud_mask_master, out_datatype=gdal.GDT_Byte)

    master_slave_info['cloud_mask'] = dst_filename

    return master_slave_info


def geom2mask(rast, polygon, gt=None):

    """
    Generates a mask by rasterizing the input geom polygon on the provided
    array according to its geoTransform. Returns a mask array.

    Arguments:
        rast            NumPy array used to rasterize the input geometry
        polygon         The path to the clipping features
        gt              An optional GDAL GeoTransform to use instead
    """

    def image_to_array(i):

        """
        Converts a Python Imaging Library (PIL) array to a gdalnumeric image.
        """

        a = gdalnumeric.fromstring(i.tobytes(), 'b')
        a.shape = i.im.size[1], i.im.size[0]
        return a

    # Only accepts numpy.array instance
    if not isinstance(rast, np.ndarray):
        raise TypeError('rast has to be a 2D numpy array.')

    # Convert the layer extent to image pixel coordinates
    minX, maxX, minY, maxY = polygon.GetEnvelope()
    ulX, ulY = misc.world_to_pixel(gt, minX, maxY)
    lrX, lrY = misc.world_to_pixel(gt, maxX, minY)

    # Calculate the pixel size of the new image
    # pxWidth = int(lrX - ulX)
    # pxHeight = int(lrY - ulY)

    if lrX > rast.shape[0] or lrY > rast.shape[1]:
        raise ValueError('The extent of the input geometry seems to exceed the raster extent.')

    pxWidth = rast.shape[1]
    pxHeight = rast.shape[0]

    # Create a new geomatrix for the image
    gt2 = list(gt)
    gt2[0] = minX
    gt2[3] = maxY

    # Map points to pixels for drawing the boundary on a blank 8-bit,
    #   black and white, mask image.
    points = []
    pixels = []
    # geom = poly.GetGeometryRef()
    pts = polygon.GetGeometryRef(0)

    for p in range(pts.GetPointCount()):
        points.append((pts.GetX(p), pts.GetY(p)))

    for p in points:
        pixels.append(misc.world_to_pixel(gt, p[0], p[1]))

    raster_poly = Image.new('L', (pxWidth, pxHeight), 1)
    rasterize = ImageDraw.Draw(raster_poly)
    rasterize.polygon(pixels, 0)  # Fill with zeroes

    mask = image_to_array(raster_poly)

    return mask, ulX, ulY, gt2


def add2mask(mask_path, shp_path, add_value=255):

    """
    Set values in a raster within the border of polyon to defined value
    Args:
        mask_path: full path to the mask to be modified
        shp_path: full path to the shapefile holding a single polygon, if multiple only the first is used
        add_value: value to be added

    Returns: modified raster on disk

    """
    ds_mask = gdal.Open(mask_path, 1)
    band_mask = ds_mask.GetRasterBand(1)
    band_array = band_mask.ReadAsArray()

    ds_shape = ogr.Open(shp_path)
    layer = ds_shape.GetLayer()
    feat = layer.GetNextFeature()
    geom = feat.GetGeometryRef()

    geotransform = ds_mask.GetGeoTransform()
    output_mask, ulX, ulY, gt2 = geom2mask(band_array, geom, gt=geotransform)
    band_array[np.where(output_mask==0)] = add_value

    band_mask.WriteArray(band_array)
    ds_mask = None