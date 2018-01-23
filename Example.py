# Example 1: Co-register two S-2 tiles (local processing)
import misc
import mask
import pprint
import os
import pickle
import coregis_filter
import gen
import correct
import correct_image
import gdal
import matplotlib.pyplot as plt
import cv2
import numpy as np
import statsmodels.api as sm

# set paths
fmask_path = './python-fmask-0.4.3/build/scripts.linux-x86_64-3.5'
otbpath = os.path.expanduser('~/OTB/install/bin')
micmac_path = os.path.expanduser('~/micmac')
work_folder = os.path.expanduser('~/Documents/DebreSina/out')

# set parameters for mask generation
pixel_size_x_cloud_mask = 30
pixel_size_y_cloud_mask = 30

# set parameters for matching
SzW=4
CorrelMin=0.3
Inc=2
RegulBase=0.2
mode='coregis'

# set parameters for corrections
iterations=3
zoom=0.25
cc_thresh=0.33

# set the input folders
master_granule = '/home/stumpf/Documents/DebreSina/S2A_OPER_PRD_MSIL1C_PDMC_20160128T013226_R092_V20160125T080606_20160125T080606.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_SGS__20160125T113126_A003092_T37PEL_N02.01'
slave_granule = '/home/stumpf/Documents/DebreSina/S2A_MSIL1C_20170119T074231_N0204_R092_T37PEL_20170119T075325.SAFE/GRANULE/L1C_T37PEL_A008240_20170119T075325'

# assure existance of output folder
if work_folder:
    os.makedirs(work_folder, exist_ok=True)

granule_info_master = misc.get_granule_info(master_granule)
granule_info_slave = misc.get_granule_info(slave_granule)
pprint.pprint(granule_info_master)
pprint.pprint(granule_info_slave)

# generate cloud mask for master
granule_info_master = mask.create_cloud_mask(granule_info_master,
                                             fmask_path,
                                             pixel_size_x=pixel_size_x_cloud_mask,
                                             pixel_size_y=pixel_size_y_cloud_mask,
                                             output_folder=work_folder)

# generate cloud mask for slave
granule_info_slave = mask.create_cloud_mask(granule_info_slave,
                                            fmask_path,
                                            pixel_size_x=pixel_size_x_cloud_mask,
                                            pixel_size_y=pixel_size_y_cloud_mask,
                                            output_folder=work_folder)

# check cloud mask
cloud_mask_ok_master = mask.check_cloud_mask(granule_info_master, cloud_values=[2, 4, 5], maximum_covered=0.7)
cloud_mask_ok_slave = mask.check_cloud_mask(granule_info_slave, cloud_values=[2, 4, 5], maximum_covered=0.7)

# make simple pair
master_slave_info = misc.make_pair(granule_info_master, granule_info_slave)
pprint.pprint(master_slave_info)

# generate bands for matching
master_slave_info = misc.generate_pan(master_slave_info, resampling_method='cubic', output_folder=work_folder)
pprint.pprint(master_slave_info)

# combine cloud masks
master_slave_info = mask.combine_masks(master_slave_info, cloud_values=[0,2,4,5],
                                       output_folder=os.path.dirname(master_slave_info['pan_master']))

# morphological filter on combined cloud mask
# close areas which will be correlated, than open (i.e. delete) areas which are smaller than open_radius
master_slave_info['cloud_mask'] = coregis_filter.morphological_filter(master_slave_info['cloud_mask'],
                                                                      close_radius=100,
                                                                      open_radius=200,
                                                                      otbpath=otbpath,
                                                                      foreground_value=255,
                                                                      background_value=0,
                                                                      output_folder=os.path.dirname(master_slave_info['cloud_mask']))

# upsample mask to full resolution and align with master
master_slave_info['cloud_mask'] = misc.align_first2second(master_slave_info['cloud_mask'],
                                                          master_slave_info['pan_master'],
                                                          resampling_method='near',
                                                          output_folder=os.path.dirname(master_slave_info['cloud_mask']))

# generate mask xml
master_slave_info['cloud_mask_xml'] = gen.generate_mask_xml(master_slave_info['cloud_mask'])

# generate MicMac parameter file
master_slave_info['correl_folder'],\
master_slave_info['param_xml'] = gen.generate_param_xml(master_slave_info['cloud_mask'],
                                                        master_slave_info['pan_master'],
                                                        master_slave_info['pan_slave'],
                                                        work_folder=os.path.dirname(master_slave_info['cloud_mask']),
                                                        SzW=SzW,
                                                        CorrelMin=CorrelMin,
                                                        Inc=Inc,
                                                        RegulBase=RegulBase,
                                                        mode='coregis')

# OPTION 1: launch the correlation locally
master_slave_info['px1'],\
master_slave_info['px2'],\
master_slave_info['cc'] = gen.do_micmac(master_slave_info['pan_master'],
                          master_slave_info['pan_slave'],
                          master_slave_info['cloud_mask'],
                          master_slave_info['cloud_mask_xml'],
                          master_slave_info['param_xml'],
                          master_slave_info['correl_folder'],
                          micmac_path,
                          run=True)

# run estimation of correction parameters
master_slave_info = correct.get_correction(master_slave_info,
                                           iterations=iterations,
                                           cc_thresh=cc_thresh,
                                           zoom=zoom, diagnose=True,
                                           full_resolution=True,
                                           do_destripe=True,
                                           horizontal_correction=True)

# warp correction grids to grid of the slave granule
master_slave_info['out']['warp_grid_dict'] = correct_image.warp_correctiongrid(master_slave_info['out']['correction_grid_px1'],
                                                                               master_slave_info['out']['correction_grid_px2'],
                                                                               master_slave_info['slave']['band_paths'],
                                                                               satellite=master_slave_info['slave']['satellite'],
                                                                               resampling_method='cubic')

# apply corrections to slave
corrected_bands = correct_image.apply_correction(master_slave_info['slave']['band_paths'],
                                                 master_slave_info['out']['warp_grid_dict'],
                                                 master_slave_info['slave']['satellite'])
master_slave_info['out']['corrected_bands'] = corrected_bands

# save and load project state
with open(os.path.join(work_folder, "info.txt"), "wb") as fp:
    pickle.dump(master_slave_info, fp)
with open(os.path.join(work_folder, "info.txt"), "rb") as fp:
    master_slave_info = pickle.load(fp)

# master image
ds = gdal.Open(master_slave_info['pan_master'])
band =ds.GetRasterBand(1)
band_array = band.ReadAsArray()
band_array = cv2.resize(band_array, None, fx=zoom, fy=zoom, interpolation=cv2.INTER_AREA)
band_array = band_array/10000
ulx, xres, xskew, uly, yskew, yres = ds.GetGeoTransform()
lrx = ulx + (ds.RasterXSize * xres)
lry = uly + (ds.RasterYSize * yres)
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ns_plot = plt.imshow(band_array, extent=(ulx, lrx, lry, uly), cmap='gray', vmin=0.05, vmax=0.3)
ax1.tick_params(axis='both', labelsize=8)
cbar = plt.colorbar(ns_plot)
cbar.set_label('TOA reflectance')
plt.title('Master image')
ax1.ticklabel_format(useOffset=False, style='plain')

# slave image
ds = gdal.Open(master_slave_info['pan_slave'])
band =ds.GetRasterBand(1)
band_array = band.ReadAsArray()
band_array = cv2.resize(band_array, None, fx=zoom, fy=zoom, interpolation=cv2.INTER_AREA)
band_array = band_array/10000

ulx, xres, xskew, uly, yskew, yres = ds.GetGeoTransform()
lrx = ulx + (ds.RasterXSize * xres)
lry = uly + (ds.RasterYSize * yres)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ns_plot = plt.imshow(band_array, extent=(ulx, lrx, lry, uly), cmap='gray', vmin=0.05, vmax=0.3)
ax1.tick_params(axis='both', labelsize=8)
cbar = plt.colorbar(ns_plot)
cbar.set_label('TOA reflectance')
plt.title('Slave image')
ax1.ticklabel_format(useOffset=False, style='plain')

# cloud mask
ds = gdal.Open(master_slave_info['cloud_mask'])
band =ds.GetRasterBand(1)
band_array = band.ReadAsArray()
band_array = cv2.resize(band_array, None, fx=zoom, fy=zoom, interpolation=cv2.INTER_AREA)
ulx, xres, xskew, uly, yskew, yres = ds.GetGeoTransform()
lrx = ulx + (ds.RasterXSize * xres)
lry = uly + (ds.RasterYSize * yres)
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ns_plot = plt.imshow(band_array, extent=(ulx, lrx, lry, uly), cmap='gray', vmin=0.0, vmax=255)
ax1.tick_params(axis='both', labelsize=8)
cbar = plt.colorbar(ns_plot)
cbar.set_label('Cloud mask')
ax1.ticklabel_format(useOffset=False, style='plain')

# correlation image
ds = gdal.Open(master_slave_info['out']['cc'])
band =ds.GetRasterBand(1)
band_array = band.ReadAsArray()
band_array = cv2.resize(band_array, None, fx=zoom, fy=zoom, interpolation=cv2.INTER_AREA)
band_array = (band_array - 127) / 128
ulx, xres, xskew, uly, yskew, yres = ds.GetGeoTransform()
lrx = ulx + (ds.RasterXSize * xres)
lry = uly + (ds.RasterYSize * yres)
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ns_plot = plt.imshow(band_array, extent=(ulx, lrx, lry, uly), cmap='gray', vmin=0.0, vmax=0.99)
ax1.tick_params(axis='both', labelsize=8)
cbar = plt.colorbar(ns_plot)
cbar.set_label('Correlation coefficient')
ax1.ticklabel_format(useOffset=False, style='plain')

# correction grids
ds = gdal.Open(master_slave_info['out']['correction_grid_px1'])
band =ds.GetRasterBand(1)
band_array = band.ReadAsArray()
band_array = cv2.resize(band_array, None, fx=zoom, fy=zoom, interpolation=cv2.INTER_AREA)
ulx, xres, xskew, uly, yskew, yres = ds.GetGeoTransform()
lrx = ulx + (ds.RasterXSize * xres)
lry = uly + (ds.RasterYSize * yres)
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ns_plot = plt.imshow(band_array, extent=(ulx, lrx, lry, uly), cmap='gray', vmin=-0.3, vmax=0.3)
ax1.tick_params(axis='both', labelsize=8)
cbar = plt.colorbar(ns_plot)
cbar.set_label('Offset [pixel]')
ax1.ticklabel_format(useOffset=False, style='plain')

ds = gdal.Open(master_slave_info['out']['correction_grid_px2'])
band =ds.GetRasterBand(1)
band_array = band.ReadAsArray()
band_array = cv2.resize(band_array, None, fx=zoom, fy=zoom, interpolation=cv2.INTER_AREA)
ulx, xres, xskew, uly, yskew, yres = ds.GetGeoTransform()
lrx = ulx + (ds.RasterXSize * xres)
lry = uly + (ds.RasterYSize * yres)
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ns_plot = plt.imshow(band_array, extent=(ulx, lrx, lry, uly), cmap='gray', vmin=-0.2, vmax=0.8)
ax1.tick_params(axis='both', labelsize=8)
cbar = plt.colorbar(ns_plot)
cbar.set_label('Offset [pixel]')
ax1.ticklabel_format(useOffset=False, style='plain')


# regression analysis before after
ds = gdal.Open(master_slave_info['pan_master'])
band =ds.GetRasterBand(1)
band_slave = band.ReadAsArray(0,0,500,500)
band_slave = band_slave/10000
band_slave = band_slave.flatten()
band_slave = band_slave.reshape(-1, 1)

ds = gdal.Open(master_slave_info['pan_slave'])
band =ds.GetRasterBand(1)
band_master = band.ReadAsArray(0,0,500,500)
band_master = band_master/10000
band_master = band_master.flatten()
band_master  = band_master.reshape(-1, 1)

results = sm.OLS(band_slave,sm.add_constant(band_master)).fit()
r2 = '%.3f' % results.rsquared
fig, ax = plt.subplots()
ax.scatter(band_master, band_slave, s=10, alpha=0.01)
X_plot = np.linspace(0,1,100)
ax.plot(band_master, results.fittedvalues, color='red', alpha=0.5, linewidth=2)
plt.xlabel('TOA reflectance master')
plt.ylabel('TOA reflectance slave')
ax.text(0.95, 0.95, r'$R^2$=' + r2,
        horizontalalignment='right',
        verticalalignment='top',
         transform=ax.transAxes)
ax.set_xlim((0.05,0.35))
ax.set_ylim((0.05,0.35))
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
ax.set_aspect(abs(x1-x0)/abs(y1-y0))


ds = gdal.Open(master_slave_info['out']['corrected_bands'][8])
band =ds.GetRasterBand(1)
band_slave = band.ReadAsArray(0,0,500,500)
band_slave = band_slave/10000
band_slave = band_slave.flatten()
band_slave = band_slave.reshape(-1, 1)

results = sm.OLS(band_slave,sm.add_constant(band_master)).fit()
r2 = '%.3f' % results.rsquared
fig, ax = plt.subplots()
ax.scatter(band_master, band_slave, s=10, alpha=0.01)
X_plot = np.linspace(0,1,100)
ax.plot(band_master, results.fittedvalues, color='red', alpha=0.5, linewidth=2)
plt.xlabel('TOA reflectance master')
plt.ylabel('TOA reflectance slave')
ax.text(0.95, 0.95, r'$R^2$=' + r2,
        horizontalalignment='right',
        verticalalignment='top',
         transform=ax.transAxes)
ax.set_xlim((0.05,0.35))
ax.set_ylim((0.05,0.35))
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
ax.set_aspect(abs(x1-x0)/abs(y1-y0))

plt.close('all')