import mask
import os
import sys
import gen
import misc
import array2geotiff

# set paths
fmask_path = os.path.dirname(sys.executable)
otbpath = os.path.expanduser('/usr/bin')
micmac_path = os.path.expanduser('~/micmac')
work_folder = os.path.expanduser('~/Documents/Morbihan/output')

dsm_t1 = os.path.expanduser('~/Documents/Morbihan/v1/GEOTIFF-20180814T190122Z-001/20171010_resample_1m1.tif')
dsm_t2 = os.path.expanduser('~/Documents/Morbihan/v1/GEOTIFF-20180814T190122Z-001/20171109_resample_1m1.tif')
dates = ['2017-10-10', '2017-11-09']

SzW=6
CorrelMin=0.5
Inc=10
RegulBase=0.5

# assure existance of output folder
if work_folder:
    os.makedirs(work_folder, exist_ok=True)

# align
dsm_t2 = misc.align_first2second(dsm_t2, dsm_t1, resampling_method='cubic', output_folder=work_folder)

# no data mask
mask_path = mask.make_no_data_mask(dsm_t1, dsm_t2, work_folder, no_data_value=-99999)
mask_xml = gen.generate_mask_xml(mask_path)

# generate param_xml
_, param_xml_path = gen.generate_param_xml(mask_path,
                                           dsm_t1,
                                           dsm_t2,
                                           work_folder,
                                           SzW=SzW,
                                           CorrelMin=CorrelMin,
                                           Inc=Inc,
                                           RegulBase=RegulBase,
                                           mode='bathy')

# run correlation
px1, px2, cc = gen.do_micmac(dsm_t1,
                             dsm_t2,
                             mask_path,
                             mask_xml,
                             param_xml_path,
                             work_folder,
                             micmac_path,
                             run=True,
                             last_step=6)

# transfer CRS to results
misc.add_georef(px1, dsm_t1)
misc.add_georef(px2, dsm_t1)
misc.add_georef(cc, dsm_t1)