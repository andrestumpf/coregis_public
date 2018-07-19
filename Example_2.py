#imports
import misc
import mask
import pprint
import os
import pickle
import gen
import correct
import coregis_filter
import sys

# set paths
fmask_path = os.path.dirname(sys.executable)
otbpath = os.path.expanduser('/usr/bin')
micmac_path = os.path.expanduser('~/micmac')
mpic_path = "/home/stumpf/Data/a2s-mpic/bin/MPIC" # TO BE ADAPTED
work_folder = '/home/astumpf/Data/Sawnla'
input_folder = '/home/astumpf/Data/Sawnla'
hpc_workfolder = 'sinard_time_series'

# get granule infos
granule_info_list = misc.get_all_granule_info(input_folder)

# HPC connection
host = 'xxx'
user = 'xxx'
password = 'xxx'
hpc_home = 'xxx'
mail = 'xxx'
grant = 'xxx'
partition = 'xxx'

# assure existance of output folder
if work_folder:
    os.makedirs(work_folder, exist_ok=True)

# set parameters for mask generation
pixel_size_x_cloud_mask = 30
pixel_size_y_cloud_mask = 30

# set parameters for matching
SzW=3
CorrelMin=0.2
Inc=7
RegulBase=0.2
mode='displacement'

# set parameters for correction
iterations=3
zoom=0.25
cc_thresh=0.33

# create cloud masks
granule_info_list = mask.create_all_cloud_masks(granule_info_list,
                                                fmask_path,
                                                pixel_size_x=pixel_size_x_cloud_mask,
                                                pixel_size_y=pixel_size_y_cloud_mask,
                                                output_folder=work_folder)

# create master slave pair sequence
master_slave_info_list = misc.make_pair_sequence(granule_info_list, match_range = 3, backward_matching = True)

# generate bands for matching
master_slave_info_list = misc.generate_pan_all(master_slave_info_list,
                                               resampling_method='cubic',
                                               output_folder=work_folder)

# generate combined cloud masks at full resolution
master_slave_info_list = mask.generate_all_combined_masks(master_slave_info_list,
                                                          cloud_values=[0,2,4,5],
                                                          close_radius=100,
                                                          open_radius=200,
                                                          otbpath=otbpath)

# generate all mask xmls
master_slave_info_list = gen.generate_all_mask_xml(master_slave_info_list)

# generate MicMac parameter file
master_slave_info_list = gen.generate_all_param_xml(master_slave_info_list,
                                                    SzW=SzW,
                                                    CorrelMin=CorrelMin,
                                                    Inc=Inc,
                                                    RegulBase=RegulBase,
                                                    mode=mode)

# launch on remote host
master_slave_info_list = misc.launch_all_remote(master_slave_info_list, host, user, password, partition, grant, mail,
                                                hpc_workfolder, slurm_template='slurm_template_1')

# download results
master_slave_info_list = misc.download_all_result(master_slave_info_list, host, user, password)

# save and load project state
with open(os.path.join(work_folder, "info.txt"), "wb") as fp:
    pickle.dump(master_slave_info_list, fp)
with open(os.path.join(work_folder, "info.txt"), "rb") as fp:
    master_slave_info_list = pickle.load(fp)

master_slave_info_list_updated = []
for i, master_slave_info in enumerate(master_slave_info_list):
    print('Correcting correlogram number ' + str(i))
    master_slave_info_updated = correct.get_correction(master_slave_info,
                                                       iterations=iterations,
                                                       cc_thresh=cc_thresh,
                                                       zoom=zoom,
                                                       diagnose=False,
                                                       full_resolution=True,
                                                       mask_full_resolution_output=False,
                                                       do_destripe=True,
                                                       horizontal_correction=False,
                                                       xlim=[-1, 1], ylim=[0, 3],
                                                       save_correction_grids=False)
    master_slave_info_list_updated.append(master_slave_info_updated)

with open(os.path.join(work_folder, "info_final.txt"), "wb") as fp:
    pickle.dump(master_slave_info_list_updated, fp)
with open(os.path.join(work_folder, "info_final.txt"), "rb") as fp:
    master_slave_info_list_updated = pickle.load(fp)


# compute mean per time step
allready_processed = []
master_slave_filtered_list = []
for i, master_slave_info in enumerate(master_slave_info_list_updated):

    print('Processing ' + str(i+1) + ' out of ' + str(len(master_slave_info_list_updated)))
    master_slave_current = {}

    # master_slave_info = master_slave_info_list_updated[12]
    if i in allready_processed:
        continue

    start_date = master_slave_info['master']['date']
    end_date = master_slave_info['slave']['date']
    master_slave_current['start_date'] = start_date
    master_slave_current['end_date'] = end_date
    master_slave_current['rmse_xy_output_forward'] = master_slave_info['corrections']['rmse_xy_output']

    # search inverted pair in the list
    for j, msi in enumerate(master_slave_info_list_updated):
        if msi['master']['date'] == end_date and msi['slave']['date'] == start_date:
            inverted_pair = msi
            allready_processed.append(j)

    master_slave_current['rmse_xy_output_backward'] = inverted_pair['corrections']['rmse_xy_output']
    cc_list = [master_slave_info['out']['cc'], inverted_pair['out']['cc']]
    master_slave_current['cc_forward'] = master_slave_info['out']['cc']
    master_slave_current['cc_backward'] = inverted_pair['out']['cc']
    forward_list = [1,-1]

    px_list = [master_slave_info['out']['px1_corrected'], inverted_pair['out']['px1_corrected']]
    outfile = os.path.join(work_folder, start_date + '_' + end_date + '_' + 'px1_mean_filter.tif')
    master_slave_current['px1_filtered'] = outfile

    coregis_filter.call_filter_block(px_list,
                                     cc_list,
                                     outfile,
                                     forward_list=forward_list,
                                     filter_type='mean',
                                     window_size=3,
                                     mask=True,
                                     cc_thresh=0.22,
                                     max_displacement=6,
                                     block_size=1024,
                                     n_threads=16)

    px_list = [master_slave_info['out']['px2_corrected'], inverted_pair['out']['px2_corrected']]
    outfile = os.path.join(work_folder, start_date + '_' + end_date + '_' + 'px2_mean_filter.tif')
    master_slave_current['px2_filtered'] = outfile

    coregis_filter.call_filter_block(px_list,
                                     cc_list,
                                     outfile,
                                     forward_list=forward_list,
                                     filter_type='mean',
                                     window_size=3,
                                     mask=True,
                                     cc_thresh=0.22,
                                     max_displacement=6,
                                     block_size=1024,
                                     n_threads=16)

    master_slave_filtered_list.append(master_slave_current)

with open(os.path.join(work_folder, "filtered_results.txt"), "wb") as fp:
    pickle.dump(master_slave_filtered_list, fp)
with open(os.path.join(work_folder, "filtered_results.txt"), "rb") as fp:
    master_slave_filtered_list = pickle.load(fp)

# compute correlation fields
master_slave_filtered_list_update = []
for master_slave_filtered in master_slave_filtered_list:

    start_date = master_slave_filtered['start_date']
    end_date = master_slave_filtered['end_date']
    cc_list = [master_slave_filtered['cc_forward'], master_slave_filtered['cc_backward']]
    cc_list_dummy = cc_list

    outfile = os.path.join(work_folder, start_date + '_' + end_date + '_' + 'cc_mean_filter.tif')
    master_slave_filtered['cc_filtered'] = outfile

    coregis_filter.call_filter_block(cc_list,
                                     cc_list,
                                     outfile,
                                     numpy_outtype='uint8',
                                     forward_list=None,
                                     filter_type='mean',
                                     window_size=3,
                                     mask=False,
                                     block_size=1024,
                                     n_threads=16)

    master_slave_filtered_list_update.append(master_slave_filtered)

master_slave_filtered_list = master_slave_filtered_list_update
with open(os.path.join(work_folder, "filtered_results.txt"), "wb") as fp:
    pickle.dump(master_slave_filtered_list, fp)
with open(os.path.join(work_folder, "filtered_results.txt"), "rb") as fp:
    master_slave_filtered_list = pickle.load(fp)


# call MPIC
import gdal
import datetime
import subprocess

def call_mpic(master_slave_filtered_list,
              outfolder,
              mpic_path=None,
              xres=None,
              yres=None,
              subpx_u=0.1,
              subpx_v=0.1,
              ref=None,
              window_size=3,
              min_disp=0.0,
              ccmin=0.2,
              gb_ram=3,
              min_frac=0.4,
              compute_mean_displacement=True,
              compute_mean_velocity=True):

    # parse arguments
    if ref is None:
        print('No reference defined. Using first displacement field as reference grid...')
        ref = master_slave_filtered_list[0]['px1_filtered']

    if xres is None:
        print('No pixel size X set. Using pixel size of the reference grid')
        ds_ref = gdal.Open(ref)
        xres = ds_ref.GetGeoTransform()[1]
        ds_ref = None

    if yres is None:
        print('No pixel size X set. Using pixel size of the reference grid')
        ds_ref = gdal.Open(ref)
        yres = abs(ds_ref.GetGeoTransform()[5])
        ds_ref = None

    start_date_list = []
    end_date_list = []
    duration_days_list = []
    for master_slave in master_slave_filtered_list:
        start_date = (datetime.datetime.strptime(master_slave['start_date'], '%Y-%m-%d'))
        end_date = (datetime.datetime.strptime(master_slave['end_date'], '%Y-%m-%d'))
        start_date_list.append(start_date)
        end_date_list.append(end_date)
        duration_days_list.append((end_date-start_date).days)

    start_series = min(start_date_list).strftime("%Y-%m-%d")
    end_series = max(end_date_list).strftime("%Y-%m-%d")

    # outfile = os.path.join(outfolder, start_series + '_' + end_series + '_' + 'VC' + '_' + str(window_size)  + '.tif')

    time_series_length = len(master_slave_filtered_list)
    cmd = [mpic_path,
           "-x", str(xres),
           "-y", str(yres),
           "-u", str(subpx_u),
           "-v", str(subpx_v),
           "-r", ref,
           "-n", str(window_size),
           "-d", str(min_disp),
           "-c", str(ccmin)]

    if compute_mean_displacement:
        cmd.append("-e")

    cmd.append("-g")
    cmd.append(str(gb_ram))
    cmd.append("-l")
    cmd.append(str(time_series_length))
    cmd.append("-m")
    cmd.append(str(min_frac))

    if compute_mean_velocity:
        cmd.append("-V")

    cmd.append("-o")
    cmd.append(outfolder)

    for master_slave, duration in zip(master_slave_filtered_list, duration_days_list):
        cmd.append(master_slave['px1_filtered'])
        cmd.append(master_slave['px2_filtered'])
        cmd.append(master_slave['cc_filtered'])
        cmd.append(str(duration))

    subprocess.run(cmd)

    return cmd


cmd = call_mpic(master_slave_filtered_list,
                work_folder,
                mpic_path=mpic_path,
                xres=None,
                yres=None,
                subpx_u=0.1,
                subpx_v=0.1,
                ref=None,
                window_size=3,
                min_disp=0.0,
                ccmin=0.2,
                gb_ram=3,
                min_frac=0.4,
                compute_mean_displacement=True,
                compute_mean_velocity=True)