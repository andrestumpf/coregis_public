import gdal
import subprocess
import os
import shutil
import glob


def generate_mask_xml(mask):

    """
    Wrapper to gen_xml_mask.pl to generate an XML for a mask compatible with MicMac
    :param mask: Full path to the mask raster
    :return: full path to the generated XML file in the same folder as the mask
    """

    ds = gdal.Open(mask)
    x_size = ds.RasterXSize
    y_size = ds.RasterYSize
    mask_xml = os.path.splitext(mask)[0] + '.xml'
    subprocess.call(['perl', './perl/gen_xml_mask.pl', str(x_size), str(y_size), mask, mask_xml])
    ds = None

    return mask_xml


def generate_all_mask_xml(master_slave_info_list):

    """
    Wrapper to generate_mask_xml to iterate over a list of master slave dictionaries
    Args:
        master_slave_info_list: list of master slave dictionaries

    Returns: list of master slave dictionaries including the path to the mask XML

    """

    master_slave_info_list_updated = []
    for master_slave_info in master_slave_info_list:

        master_slave_info['cloud_mask_xml'] = generate_mask_xml(master_slave_info['cloud_mask'])
        master_slave_info_list_updated.append(master_slave_info)

    return master_slave_info_list_updated


def generate_param_xml(mask_path, pan_1_path, pan_2_path, work_folder, sub_folder=None, SzW=4, CorrelMin=0.3, Inc=1,
                       RegulBase=0.5,
                       mode='coregis'):
    """
    Wrapper to gen_xml_params_displacement.pl to generate a MicMac parameter file

    :param mask_path: full path to the mask raster
    :param pan_1_path: full path to the first image to be correlated
    :param pan_2_path: full path to the second image to be correlated
    :param work_folder: full path to the work folder where processing is performed
    :param sub_folder: optional name of a sub folder within the work folder where correlation is performed
    :param SzW: size of the correlation window
    :param CorrelMin: minimum correlation threshold
    :param Inc: size of the search area (actual search area is ceiling(Inc*0.8)+2
    :param RegulBase: regularization parameter
    :param mode: 'coregis' or 'displacement' or 'bathy', displacement yields finer resolution but considerably increases
                 the runtime, 'bathy' is adapted for large displacements of sumarine dunes
    :return: full paths to the current work folder and the generated MicMac parameter file
    """

    mask_name = os.path.splitext(os.path.basename(mask_path))[0]
    pan_1_name = os.path.basename(pan_1_path)
    pan_2_name = os.path.basename(pan_2_path)

    if sub_folder:
        correl_folder = os.path.join(work_folder, sub_folder)
        os.makedirs(correl_folder, exist_ok=True)
        correl_folder_current = os.path.join(correl_folder,
                                             os.path.splitext(pan_1_name)[0] + '_' + os.path.splitext(pan_2_name)[0])
        os.makedirs(correl_folder_current, exist_ok=True)
    else:
        correl_folder_current = work_folder

    out_xml = os.path.join(correl_folder_current, 'param.xml')

    # generate MicMac XML file
    if mode == 'coregis':
        script = './perl/gen_xml_params_coregis.pl'
    if mode == 'displacement':
        script = './perl/gen_xml_params_displacement.pl'
    if mode == 'bathy':
        script = './perl/gen_xml_params_bathy.pl'

    subprocess.call(['perl', script,
                     '-s', str(SzW),
                     '-c', str(CorrelMin),
                     '-i', str(Inc),
                     '-r', str(RegulBase),
                     '-m', mask_name,
                     '-o', out_xml,
                     pan_1_name,
                     pan_2_name])

    return correl_folder_current, out_xml


def generate_all_param_xml(master_slave_info_list, SzW=4, CorrelMin=0.3, Inc=2, RegulBase=0.2, mode='coregis'):

    """
    Wrapper to generate_param_xml to iterate over a list of master slave dictionaries
    Args:
        master_slave_info_list: list of master slave dictionaries
        SzW: size of the correlation window
        CorrelMin: minimum correlation threshold
        Inc: size of the search area (actual search area is ceiling(Inc*0.8)+2
        RegulBase: regularization parameter
        mode: 'coregis' or 'displacement', displacement yields finer resolution but considerably increases the runtime

    Returns: list of master slave dictonaries including the path to the parameter xml

    """

    master_slave_info_list_updated = []

    for master_slave_info in master_slave_info_list:
        master_slave_info['correl_folder'],\
        master_slave_info['param_xml'] = generate_param_xml(master_slave_info['cloud_mask'],
                                                            master_slave_info['pan_master'],
                                                            master_slave_info['pan_slave'],
                                                            work_folder=os.path.dirname(master_slave_info['cloud_mask']),
                                                            SzW=SzW,
                                                            CorrelMin=CorrelMin,
                                                            Inc=Inc,
                                                            RegulBase=RegulBase,
                                                            mode=mode)

        master_slave_info_list_updated.append(master_slave_info)

    return master_slave_info_list_updated


def do_micmac(raster1, raster2, mask, mask_xml, param_xml, correl_folder, micmac_path, last_step=6, run=True):

    """
    Wrapper to launch MicMac or prepare MicMac 2D correlation
    :param raster1: full path to first input image
    :param raster2: full path to second input image
    :param mask: full path to micmac compatible mask
    :param mask_xml: full path to micmac compatible mask XML
    :param param_xml: full path to the MicMac parameter XML
    :param correl_folder: full path to work folder for execution
    :param micmac_path:  full path to micmac install folder
    :param last_step: last step of the correlation
    :param run: if True correlation will be launched, set False to prepare only the work folder
    :return: full paths to Px1, Px2, and CC raster outputs, or full paths to the work folder if run=False

    DEBUG
    correl_folder = master_slave_info['correl_folder']
    raster1 = master_slave_info['pan_master']
    """

    # copy necessary files to work folder
    print('Copying all necessary files to the correlation work folder...')
    if os.path.exists(os.path.join(correl_folder, os.path.basename(raster1))):
        print(raster1 + ' is already in the correlation folder')
    else:
        try:
            shutil.copy(raster1, correl_folder)
        except IOError:
            print('Unable to copy file ' + raster1)
    if os.path.exists(os.path.join(correl_folder, os.path.basename(raster2))):
        print(raster2 + ' is already in the correlation folder')
    else:
        try:
            shutil.copy(raster2, correl_folder)
        except IOError:
            print('Unable to copy file ' + raster2)
    if os.path.exists(os.path.join(correl_folder, os.path.basename(mask))):
        print(mask + ' is already in the correlation folder')
    else:
        try:
            shutil.copy(mask, correl_folder)
        except IOError:
            print('Unable to copy file ' + mask)
    if os.path.exists(os.path.join(correl_folder, os.path.basename(mask_xml))):
        print(mask_xml + ' is already in the correlation folder')
    else:
        try:
            shutil.copy(mask_xml, correl_folder)
        except IOError:
            print('Unable to copy file ' + mask_xml)
    if os.path.exists(os.path.join(correl_folder, os.path.basename(param_xml))):
        print(param_xml + ' is already in the correlation folder')
    else:
        try:
            shutil.copy(param_xml, correl_folder)
        except IOError:
            print('Unable to copy file ' + param_xml)

    # launch correlation
    if run:
        print('Running image correlation...')
        mm3d_path = os.path.join(micmac_path, 'bin/mm3d')
        param_xml_name = os.path.basename(param_xml)
        subprocess.call([mm3d_path, 'MICMAC', param_xml_name], cwd=correl_folder)

        # get paths to MicMac's outputs
        print('Query outputs ... ')
        px1 = glob.glob(correl_folder + '/**/*Px1_Num' + str(last_step) + '*.tif', recursive=True)[0]
        px2 = glob.glob(correl_folder + '/**/*Px2_Num' + str(last_step) + '*.tif', recursive=True)[0]
        cc = glob.glob(correl_folder + '/**/*Correl_LeChantier_Num_' + str(int(last_step) - 1) + '*.tif',
                       recursive=True)[0]

    else:
        px1 = 'not yet computed'
        px2 = 'not yet computed'
        cc = 'not yet computed'

    return px1, px2, cc
