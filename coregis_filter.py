import subprocess
import os
import gdal
import cv2
import numpy as np
import warnings
from rios import applier
from rios import cuiprogress
from rios import parallel


def morphological_filter(input_path, close_radius=None, open_radius=200, otbpath='',
                         foreground_value=1, background_value=0, output_folder=None):

    """
    Wrapper to OTB mathematical morphology filter command line application. First runs a closing
    and subsequently an opening operation.
    :param input_path: full path to the input image to be filtered
    :param close_radius: radius of the structuring element (ball) for closing in meters
    :param open_radius: radius of the structuring element (ball) for open in meters
    :param otbpath: path to the OTB binaries including otbcli_BinaryMorphologicalOperation
    :param foreground_value: Foreground Value
    :param background_value: Background Value
    :param output_folder: folder to which the output rasters will be written
    :return: path to the raster after binary morphological filtering
    """

    # figure out thieving sizes in pixels
    image_link = gdal.Open(input_path)
    gt = image_link.GetGeoTransform()
    x_res = gt[1]
    y_res = abs(gt[5])

    if close_radius is not None:

        # figure out thieving sizes in pixels
        close_radius_pixel_x = round(close_radius / x_res )
        close_radius_pixel_y = round(close_radius / y_res)

        # concatenate output name
        file_path_woext = os.path.splitext(input_path)[0]
        if output_folder:
            outputimage_close = os.path.join(output_folder, os.path.basename(file_path_woext) + '_close.tif')
        else:
            outputimage_close = file_path_woext + '_close.tif'

        # call otbcli_BinaryMorphologicalOperation closing
        subprocess.call([os.path.join(otbpath, "otbcli_BinaryMorphologicalOperation"),
                         "-in", input_path,
                         "-out", outputimage_close, "uint8",
                         "-structype.ball.xradius", str(close_radius_pixel_x),
                         "-structype.ball.yradius", str(close_radius_pixel_y),
                         "-filter", "closing",
                         "-filter.dilate.foreval", str(foreground_value),
                         "-filter.dilate.backval", str(background_value),
                         "-filter.erode.foreval", str(foreground_value),
                         "-filter.erode.backval", str(background_value),
                         "-filter.opening.foreval", str(foreground_value),
                         "-filter.opening.backval", str(background_value),
                         "-filter.closing.foreval", str(foreground_value)])

    else:
        outputimage_close = input_path

    if open_radius is not None:

        # figure out thieving sizes in pixels
        open_radius_pixel_x = round(open_radius / x_res )
        open_radius_pixel_y = round(open_radius / y_res)

        # concatenate output name
        file_path_woext = os.path.splitext(outputimage_close)[0]
        if output_folder:
            outputimage_open = os.path.join(output_folder, os.path.basename(file_path_woext) + '_open.tif')
        else:
            outputimage_open = file_path_woext + '_open.tif'

        # call otbcli_BinaryMorphologicalOperation closing
        subprocess.call([os.path.join(otbpath, "otbcli_BinaryMorphologicalOperation"),
                         "-in", outputimage_close,
                         "-out", outputimage_open, "uint8",
                         "-structype.ball.xradius", str(open_radius_pixel_x),
                         "-structype.ball.yradius", str(open_radius_pixel_y),
                         "-filter", "opening",
                         "-filter.dilate.foreval", str(foreground_value),
                         "-filter.dilate.backval", str(background_value),
                         "-filter.erode.foreval", str(foreground_value),
                         "-filter.erode.backval", str(background_value),
                         "-filter.opening.foreval", str(foreground_value),
                         "-filter.opening.backval", str(background_value),
                         "-filter.closing.foreval", str(foreground_value)])

        return outputimage_open

    elif close_radius is not None:

        return outputimage_close

    else:
        raise Exception("Both open and close radius are None.")


def filter_block(info, inputs, outputs, otherargs):

    # print('Current block is')
    # print(info.getBlockCount())
    # print('out of')
    # print(info.getTotalBlocks())

    if not otherargs.filter in ['NLM', 'median', 'mean']:
       raise SyntaxError(otherargs.filter + ' is not a valid filter. Use NLM, median, or mean')
    if (otherargs.window_size % 2) == 0:
       raise SyntaxError('window size must be an odd integer')
    if not 0 <= otherargs.cc_thresh < 1:
       raise SyntaxError('correlation threshold must be  >=0 and <1')
    if otherargs.mask:
        if len(inputs.offsets) == inputs.cc:
            raise SyntaxError('if mask=True the length of the input lists must be identical')
    if otherargs.filter == 'NLM':
        nlm_h = otherargs.nlm_h
        nlm_search_size = otherargs.nlm_search_size

    filter = otherargs.filter
    window_size = otherargs.window_size
    px_block = inputs.offsets
    forward_list = otherargs.forward_list

    # DEBUG
    # from rios.imagereader import ImageReader
    # from rios.readerinfo import ReaderInfo
    # reader = ImageReader(inputs.offsets, overlap=otherargs.window_size,  windowxsize=block_size, windowysize=block_size)
    # px_block = reader.readBlock(0)
    # block_coordinates_x, block_coordinates_y = px_block[0].getBlockCoordArrays()
    # px_block = px_block[1]
    # print(block_coordinates_x[0,0])
    # print(block_coordinates_y[0,0])


    # masking if requested
    if otherargs.mask:
        cc_block = inputs.cc

        # DEBUG
        # reader = ImageReader(inputs.cc, overlap=otherargs.window_size, windowxsize=block_size, windowysize=block_size)
        # cc_block = reader.readBlock(0)
        # cc_block = cc_block[1]

        # print('This is cc_block')
        # print(cc_block)
        # print(len(cc_block))
        cc_thresh_real = otherargs.cc_thresh * 128 + 127
        cc_bool = []
        for cc_patch in cc_block:
            # print('Shape of cc_patch is..')
            # print(cc_patch.shape)
            cc_patch = cc_patch[0]
            cc_bool.append(cc_patch < cc_thresh_real)

        # print('This is cc_bool is..')
        # print(cc_bool)
        # print(len(cc_bool))
        px_patch_masked = []
        for px_patch, cc_patch, forward in zip(px_block, cc_bool, forward_list):
            # print('Shape of px_patch is..')
            # print(px_patch.shape)
            px_patch  = px_patch[0]*forward
            px_patch = np.ma.masked_array(px_patch, mask=cc_patch)
            px_patch_masked.append(np.ma.masked_where(np.absolute(px_patch) > abs(otherargs.max_displacement), px_patch))
    else:
        px_patch_masked = []
        for px_patch, forward in zip(px_block, forward_list):
            px_patch = px_patch[0] * forward
            px_patch_masked.append(px_patch)

        # plt.imshow(cc_patch)
        # plt.imshow(cc_bool[1])
        # plt.imshow(px_patch_masked[1])

    if filter == 'median' or filter == 'mean':

        px_patch_stack = np.ma.dstack(px_patch_masked)
        # print(px_patch_stack.dtype)
        if otherargs.mask:
            px_patch_stack = px_patch_stack.filled(np.nan)
        # print('Shape of px_patch_stack is..')
        # print(px_patch_stack.shape)
        w, h, d = px_patch_stack.shape
        # array_mask = image.mask
        strided_image = np.lib.stride_tricks.as_strided(px_patch_stack,
                                                        shape=[w - window_size + 1,
                                                               h - window_size + 1,
                                                               1,
                                                               window_size,
                                                               window_size,
                                                               d],
                                                        strides=px_patch_stack.strides + px_patch_stack.strides)
        strided_image = np.squeeze(strided_image)

        # print('Shape of strided stack is..')
        # print(strided_image.shape)

        if filter == 'median':
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                average_image = np.nanmedian(strided_image, axis=(2, 3, 4))

        if filter == 'mean':
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                average_image = np.nanmean(strided_image, axis=(2,3,4))

        # padding to get shape expected by the rios writer
        average_image = np.lib.pad(average_image, int((window_size-1)/2), 'constant', constant_values=0)
        # plt.imshow(average_image)

        average_image = np.expand_dims(average_image, axis=0)
        # print('Shape of average_image is..')
        # print(average_image.shape)

    if filter =='NLM':

        # stretch input to 0,255 range required by the OpenCV implementation
        min_max = [func(l) for l in px_patch_masked for func in (np.min, np.max)]
        min_px1 = np.min(min_max)
        max_px1 = np.max(min_max)
        noisy = [(px1_patch - min_px1)/(max_px1-min_px1)*255 for px1_patch in px_patch_masked]
        noisy = [np.uint8(i) for i in noisy]

        # fig, axs = plt.subplots(3,3)
        # for ax, image in zip(axs.reshape(-1), noisy):
        #     ax = ax.imshow(image, vmin=0, vmax=255)

        # Denoise with NLM
        target_img = int(len(noisy)/2)
        temp_window = (len(noisy) - target_img) // 2 * 2 + 1

        average_image = cv2.fastNlMeansDenoisingMulti(noisy,
                                                     target_img,
                                                     temp_window,
                                                     None,
                                                     nlm_h,
                                                     window_size,
                                                     nlm_search_size)
        # convert to pixel
        average_image = average_image/255*(max_px1-min_px1) + min_px1

        # fig, axs = plt.subplots(1,2)
        # axs[0].imshow(average_image, vmin=-2, vmax=2)
        # axs[1].imshow(px_block[1][0][0, :, :], vmin=-2, vmax=2)


        # np.array_equal(px1_median_filtered[ :, :, 0], px1_median_filtered[ :, :, 1])
        # plt.imshow(px1_median_filtered[ :, :, 0] - px1_median_filtered[ :, :, 1])
        #
        #
        # fig, axs = plt.subplots(2,3)
        # axs[0,0].imshow(px1_median_filtered[ :, :, 0], vmin=-2, vmax=2)
        # axs[0,0].set_title("Median filter 21*21")
        # axs[0,1].imshow(px1_median_filtered[ :, :, 1], vmin=-2, vmax=2)
        # axs[0,1].set_title("Median filter 21*21")
        # axs[0,2].imshow(px1_median_filtered[ :, :, 6], vmin=-2, vmax=2)
        # axs[0,2].set_title("Median filter 21*21")
        # axs[1,0].imshow(px1_block[1][0][0, :, :], vmin=-2, vmax=2)
        # axs[1,0].set_title("Input")
        # axs[1,1].imshow(px1_denoised, vmin=-2, vmax=2)
        # axs[1,1].set_title("NLM_seq 3, 31, 51")
        # axs[1,2].imshow(px1_median_filtered_mean, vmin=-2, vmax=2)
        # axs[1,2].set_title("Median of the medians")

        # plt.close('all')

    if otherargs.numpy_outtype:
        average_image = average_image.astype(otherargs.numpy_outtype)

    outputs.outimage = average_image


def call_filter_block(px_list,
                      cc_list,
                      outfile,
                      numpy_outtype=None,
                      forward_list=None,
                      filter_type='mean',
                      window_size=11,
                      nlm_h=2,
                      nlm_search_size=21,
                      mask=True,
                      cc_thresh=0.23,
                      max_displacement=4,
                      block_size=512,
                      n_threads = 4):

    """
    Caller for filter_block which enables to apply mean, median or NLM filter on list of redundant displacement fields
    Args:
        px_list: list of file paths to the displacement fields
        cc_list: list of the corresponding correlation coefficients
        outfile: full path to the output file which will be written in GeoTiff format
        numpy_outtype: provide a numpy dtype as string to gain control over the output type
        forward_list: list of indicating if the pair is in temporal (1) or inverse order (-1), if not provided temporal
                      order is assumed for all pairs
        filter_type: valid options are 'median', 'mean', 'NLM',
                     the non local mean filter makes use of the OpenCV implementation fastNlMeansDenoisingMulti,
                     the image in the middle
        window_size:     window size for filtering
        nlm_h:           the h parameter for NLM
        nlm_search_size: the search window for NLM
        mask:            if set true max_displacement and cc_thresh will be used for masking
        cc_thresh:       threshold on the correlation coefficientn, set to zero to filter only by max_displacement values
        max_displacement: offset values above this value will be filtered out
        block_size:       block size for RIOS block wise processing
        n_threads:        number of paralled threads used with RIOS

    Returns: filtered GeoTiff on disk according to the specified output path

    """

    # inputs
    inputs = applier.FilenameAssociations()
    inputs.offsets = px_list
    inputs.cc = cc_list

    # parameters
    otherargs = applier.OtherInputs()
    otherargs.numpy_outtype = numpy_outtype
    otherargs.filter = filter_type
    otherargs.window_size = window_size
    otherargs.mask = mask
    otherargs.cc_thresh = cc_thresh
    otherargs.max_displacement = max_displacement
    otherargs.nlm_h=nlm_h
    otherargs.nlm_search_size=nlm_search_size
    if forward_list:
        otherargs.forward_list = forward_list
    else:
        otherargs.forward_list = [1] * len(px_list)

    # controls
    controls = applier.ApplierControls()
    controls.setOverlap(otherargs.window_size)
    controls.setWindowXsize(block_size)
    controls.setWindowYsize(block_size)
    controls.progress = cuiprogress.CUIProgressBar()
    controls.setOutputDriverName("GTiff")
    controls.setCreationOptions(["COMPRESS=NONE"])
    controls.setCreationOptions(["INTERLEAVE=PIXEL"])

    # parallelization
    if n_threads > 1:
        controls.setNumThreads(n_threads)
        controls.setJobManagerType('subproc')

    #output file
    outfiles = applier.FilenameAssociations()
    outfiles.outimage = outfile

    applier.apply(filter_block, inputs, outfiles, otherargs, controls=controls)


def min_max_filter(info, inputs, outputs, otherargs):

    print('Current block is')
    print(info.getBlockCount())
    print('out of')
    print(info.getTotalBlocks())

    in_array = inputs.image
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        in_array[in_array < otherargs.min] = np.nan
        in_array[in_array > otherargs.max] = np.nan

    outputs.outimage = in_array


def call_min_max_filter(infile, outfile, min=-1, max=1, block_size=512, n_threads=5):

    """
    Function caller to min_max_filter which set values below the minimum or above the maximum to nan
    Args:
        infile: full path to the input file
        outfile: full path to the output file
        min: minimum allowed value
        max: maximum allowed value
        block_size: block size for RIOS block wise processing
        n_threads: number of threads

    Returns: GeoTiff with extreme values set to nan

    """

    # inputs
    inputs = applier.FilenameAssociations()
    inputs.image = infile

    # parameters
    otherargs = applier.OtherInputs()
    otherargs.min = min
    otherargs.max = max

    # controls
    controls = applier.ApplierControls()
    controls.setWindowXsize(block_size)
    controls.setWindowYsize(block_size)
    controls.progress = cuiprogress.CUIProgressBar()
    controls.setOutputDriverName("GTiff")
    controls.setCreationOptions(["COMPRESS=NONE"])
    controls.setCreationOptions(["INTERLEAVE=PIXEL"])

    # parallelization currently not working due to issues with pickeling
    if n_threads > 1:
        # parallel.jobmanager.getAvailableJobManagerTypes()
        controls.setNumThreads(n_threads)
        controls.setJobManagerType('subproc')

    # output file
    outfiles = applier.FilenameAssociations()
    outfiles.outimage = outfile

    applier.apply(min_max_filter, inputs, outfiles, otherargs, controls=controls)


def call_compute_variance(px_list,
                          cc_list,
                          outfile,
                          forward_list=None,
                          mask=True,
                          cc_thresh=0.33,
                          max_displacement=1.5,
                          block_size=512,
                          n_threads=4):

    """
    RIOS style caller to compute_variance, which computes the variance per pixel along the time axis

    Args:
        px_list: list of file paths to the displacement fields
        cc_list: list of the corresponding correlation coefficients
        outfile: full path to the output file which will be written in GeoTiff format
        forward_list: list of indicating if the pair is in temporal (1) or inverse order (-1), if not provided temporal
                      order is assumed for all pairs
        mask:         if set true max_displacement and cc_thresh will be used for masking
        cc_thresh:    threshold on the correlation coefficientn, set to zero to filter only by max_displacement values
        max_displacement: offset values above this value will be filtered out
        block_size:   block size for RIOS block wise processing
        n_threads:    number of paralled threads used with RIOS

    Returns: filtered GeoTiff on disk according to the specified output path

    """

    # inputs
    inputs = applier.FilenameAssociations()
    inputs.offsets = px_list
    inputs.cc = cc_list

    # parameters
    otherargs = applier.OtherInputs()
    otherargs.mask = mask
    otherargs.cc_thresh = cc_thresh
    otherargs.max_displacement = max_displacement
    if forward_list:
        otherargs.forward_list = forward_list
    else:
        otherargs.forward_list = [1] * len(px_list)

    # controls
    controls = applier.ApplierControls()
    controls.setWindowXsize(block_size)
    controls.setWindowYsize(block_size)
    controls.progress = cuiprogress.CUIProgressBar()
    controls.setOutputDriverName("GTiff")
    controls.setCreationOptions(["COMPRESS=NONE"])
    controls.setCreationOptions(["INTERLEAVE=PIXEL"])

    # parallelization
    if n_threads > 1:
        controls.setNumThreads(n_threads)
        controls.setJobManagerType('subproc')

    #output file
    outfiles = applier.FilenameAssociations()
    outfiles.outimage = outfile

    applier.apply(compute_std, inputs, outfiles, otherargs, controls=controls)


def compute_std(info, inputs, outputs, otherargs):

    if not 0 <= otherargs.cc_thresh < 1:
       raise SyntaxError('correlation threshold must be  >=0 and <1')
    if otherargs.mask:
        if len(inputs.offsets) == inputs.cc:
            raise SyntaxError('if mask=True the length of the input lists must be identical')

    px_block = inputs.offsets
    forward_list = otherargs.forward_list

    # masking if requested
    if otherargs.mask:
        cc_block = inputs.cc

        # compute cc mask
        cc_thresh_real = otherargs.cc_thresh * 128 + 127
        cc_bool = []
        for cc_patch in cc_block:
            cc_patch = cc_patch[0]
            cc_bool.append(cc_patch < cc_thresh_real)

        # apply cc mask and max_displacment threshold
        px_patch_masked = []
        for px_patch, cc_patch, forward in zip(px_block, cc_bool, forward_list):
            # print('Shape of px_patch is..')
            # print(px_patch.shape)
            px_patch = px_patch[0] * forward
            px_patch = np.ma.masked_array(px_patch, mask=cc_patch)
            px_patch_masked.append(
                np.ma.masked_where(np.absolute(px_patch) > abs(otherargs.max_displacement), px_patch))

    else:
        px_patch_masked = []
        for px_patch, forward in zip(px_block, forward_list):
            px_patch = px_patch[0] * forward
            px_patch_masked.append(px_patch)

    px_patch_stack = np.ma.dstack(px_patch_masked)
    # print('Shape of px_patch_stack is..')
    # print(px_patch_stack.shape)
    if otherargs.mask:
        px_patch_stack = px_patch_stack.filled(np.nan)

    # compute standard deviation for each pixel
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        sigma_px = np.nanstd(px_patch_stack, axis=2)

    # allocate to output and add 3d dimension before
    outputs.outimage = np.expand_dims(sigma_px, axis=0)