import numpy as np
import gdal
from osgeo import gdal_array


def array2geotiff(band_array, dst_filename, ref, out_datatype=None):

    """
    Wrapper to GDAL to write an array to disk in GeoTiff format.

    Args:
        band_array: array that will be written to disk
        dst_filename: full path to the output GeoTiff
        ref: reference image with defined geotransform and SRS or a tuple of the form (geotransform, SRS_wkt)
        out_datatype: GDAL datatype for the output, if None the function tries to figure out the corresponding
                        type automatically

    Returns: full path to the GeoTiff on disk.

    """

    # DEBUG
    # band_array = band_corrected_array
    # dst_filename = dst_filename
    # ref = (out_gt, prj)

    # if the input array is a masked array this will set the masked values to zero
    np.ma.set_fill_value(band_array, 0)

    # figure out the output datatype according to the type of the input array
    if not out_datatype:
        out_datatype = gdal_array.NumericTypeCodeToGDALTypeCode(band_array.dtype)

    # from tuple
    if isinstance(ref, tuple):
        geotrans = ref[0]
        proj_info = ref[1]
        cols = band_array.shape[1]
        rows = band_array.shape[0]
    # or provided reference image
    else:
        ref_image = gdal.Open(ref)
        geotrans = ref_image.GetGeoTransform()
        proj_info = ref_image.GetProjection()
        cols = round(ref_image.RasterXSize)
        rows = round(ref_image.RasterYSize)

    driver = gdal.GetDriverByName('GTiff')
    driver.Create(dst_filename, cols, rows, 1, gdal.GDT_Byte)
    out_raster = driver.Create(dst_filename, cols, rows, 1, eType=out_datatype)
    out_raster.SetGeoTransform(geotrans)
    out_raster.SetProjection(proj_info)
    outband = out_raster.GetRasterBand(1)

    if hasattr(band_array, 'mask'):
        outband.WriteArray(band_array.filled())
    else:
        outband.WriteArray(band_array)
    outband.FlushCache()

    out_raster = None
    ref_image = None

    return dst_filename

