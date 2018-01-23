import numpy as np
import gdal
import warnings
import cv2
import matplotlib.pyplot as plt
import gc
import misc
import os
import ogr
import footprints
import mask
from scipy import stats
import array2geotiff
import correct_image
import destripe


def fit_plane_iter(in_surface, tune=4.685, iterations=10, tol=1e-5):

    """
    Fit a plane two a 2D array using robust least square (L1-norm)

    Args:
        in_surface:  2D array to be fitted, if masked array the mask is not taken into account
        tune:        Tuckey tuning constant (default 4.685, currently set to common default)
        iterations:  Maximum number of iterations (default 10)
        tol:         Stopping criteria. The fiting will stop before iter is reached if the residuals are below
                     this threshold.

    Returns: fitted plane as 2D numpy array, MAD of the fits during each iteration

    """

    # create a coordinate matrix
    x_size = in_surface.shape[1]
    y_size = in_surface.shape[0]
    nx = np.linspace(-1, 1, x_size)
    ny = np.linspace(-1, 1, y_size)
    x, y = np.meshgrid(nx, ny)

    # construct design matrix
    x_fl = x.flatten()
    y_fl = y.flatten()
    z_ones = np.ones([x.size, 1])

    # transfer mask
    if hasattr(in_surface, 'mask'):
        mask = in_surface.mask.flatten()
        x_fl = x_fl[np.invert(mask)]
        y_fl = y_fl[np.invert(mask)]
        z_ones = z_ones[np.invert(mask)]

    X = np.hstack((np.reshape(x_fl, ([len(x_fl), 1])), np.reshape(y_fl, ([len(y_fl), 1])), z_ones))

    # construct response matrix
    Z = np.zeros(in_surface.shape)
    Z[:] = in_surface
    Z_fl = Z.flatten()
    # transfer mask
    if hasattr(in_surface, 'mask'):
        Z_fl = Z_fl[np.invert(mask)]
    Z = np.reshape(Z_fl, ([len(Z_fl), 1]))

    # initiate with standard least square
    A_lsq = np.linalg.lstsq(X, Z)[0]

    # compute absolute value of residuals (fit minus data)
    abs_resid = abs(np.dot(X, A_lsq) - Z)

    # store residuals for output
    iter_mad = np.median(abs_resid)

    # iterate till the fit converges

    for i in range(iterations):

        print('Running iteration', i)

        # compute the scaling factor for the standardization of residuals
        # using the median absolute deviation of the residuals
        # 6.9460 is a tuning constant (4.685/0.6745)
        # equivalent to  http://users.stat.umn.edu/~sandy/courses/8053/handouts/robust.pdf

        if np.median(abs_resid) < tol:
            print('MAR=', np.median(abs_resid))
            print('Convergence reached after', i, 'iteration(s)')

            if i == 0:

                print('This may indicate an issue with the input data as for example constant 0 grids.')
                print('Output will be empty.')
                out_surface = None
                return out_surface

            else:
                break

        abs_res_scale = 6.9460 * np.median(abs_resid)

        # standardize residuals
        w = abs_resid / abs_res_scale

        # compute the robust bisquare weights excluding outliers
        outliers = (w > 1) * 1
        w[outliers.nonzero()] = 0

        good_values = (w != 0) * 1

        # calculate robust weights for 'good' points
        tmp = 1 - np.power(w[good_values.nonzero()], 2)
        w[good_values.nonzero()] = np.power(tmp, 2)

        # get weighted X'es
        XW = np.tile(w, (1, 3)) * X

        a = np.dot(XW.T, X)
        b = np.dot(XW.T, Z)

        # get the least-squares solution to a linear matrix equation
        A_robust = np.linalg.lstsq(a, b)[0]

        # recompute absolute value of residuals (fit minus data)
        abs_resid = abs(np.dot(X, A_robust) - Z)

        iter_mad = np.append(iter_mad, np.median(abs_resid))
        print('MAD =', np.median(abs_resid))

    # get values of the aproximated plane

    # reconstruct unmasked design matrix to get an output values for all input elements (including masked values)
    x_fl = x.flatten()
    y_fl = y.flatten()
    z_ones = np.ones([x.size, 1])
    X = np.hstack((np.reshape(x_fl, ([len(x_fl), 1])), np.reshape(y_fl, ([len(y_fl), 1])), z_ones))

    out_values = np.dot(X, A_robust)

    # # bring it into the format of the input matrix

    # if hasattr(in_surface, 'mask'):
    #     out_surface = np.zeros(in_surface.shape)
    #     out_surface[np.invert(in_surface.mask)] = out_values.flatten()
    # else:

    out_surface = np.reshape(out_values, in_surface.shape)

    return out_surface, iter_mad.tolist()


def mask_large_values(in_surface, max_displacement):
    """Returns a masked array if max_displacement is exeeded anywhere on the map"""
    mask_large_values = np.absolute(in_surface) > abs(max_displacement)
    if np.any(mask_large_values):
        return np.ma.masked_where(np.absolute(in_surface) > abs(max_displacement), in_surface)
    else:
        return in_surface


def get_full_resolution_output(displacement_field, correction_surface, cc_field, mask_full_resolution_output,
                               cc_thresh, max_displacement):

    """

    Args:
        displacement_field: full path to the displacement field
        correction_surface:  numpy array holding the correction surface
        cc_field:  full path to the correlation coefficient
        mask_full_resolution_output: if set False the masks (cc, max_displacement, will not be applied on the output)
        cc_thresh: threshold for masking values according to the correlation coefficient
        max_displacement: all displacement values whose absolute magnitude exceeds this value will be masked out for the correction

    Returns: corrected, resampled (and optionally masked displacement field) at the full resolution

    """

    # re-import displacement field at original resolution
    displacement_field_link = gdal.Open(displacement_field)
    px_band = displacement_field_link.GetRasterBand(1)
    in_surface = np.array(px_band.ReadAsArray())

    # resample correction surface to the resolution of ht displacement field
    correction_surface = cv2.resize(correction_surface,
                                    (in_surface.shape[1], in_surface.shape[0]),
                                    interpolation=cv2.INTER_CUBIC)

    # apply correction
    in_surface -= correction_surface

    if mask_full_resolution_output:

        if cc_field is not None:

            # re-import correlation coefficient at original resolution
            cc_link = gdal.Open(cc_field)
            px_band = cc_link.GetRasterBand(1)
            in_cc = np.array(px_band.ReadAsArray())

            # create mask
            in_cc = (in_cc - 127) / 128 < cc_thresh

            # apply mask
            in_surface = np.ma.masked_array(in_surface, mask=in_cc)

    # mask large displacements
    in_surface = mask_large_values(in_surface, max_displacement)

    return in_surface, correction_surface


def get_correction_plane_ms(displacement_field,
                            master_slave,
                            ref_path,
                            cc_field=None,
                            cc_thresh=0.333,
                            iterations=3,
                            zoom=1.0,
                            diagnose=False,
                            full_resolution=False,
                            mask_full_resolution_output=True,
                            do_destripe=False,
                            horizontal_correction=False,
                            destriping_statistics='median',
                            max_displacement=-9999,
                            xlim=[-1, 1], ylim=[0, 3]):

    """
    Approximate the displacement field with a plane, destripe and return the the plane that can be used as a correction
     grid.

    Args:
        displacement_field: full path to the displacement fields that should be corrected
        master_slave: master slave dictionary
        ref_path: full path to a reference image, if the displacement field has no georeference the reference of the
                    reference image will be used
        cc_field: full path to the raster holding the correlation coefficients
        cc_thresh: threshold on the correlation coefficient, pixels below this value will be ignored during the
                     correction
        iterations: number of iterations for iteratively reweighted plane fitting
        zoom: optional downsampling, mainly for reducing memory size and processing time
        diagnose: if True plotting is activated to illustrate the corrections
        full_resolution: if True the final grid will be returned at the image resolution
        mask_full_resolution_output: if set False the masks (cc, max_displacement, will not be applied on the output)
        do_destripe: if True a destriping routine will applied
        horizontal_correction: if True a second destriping routine will be applied to also correct horizontal artifacts,
                                only if do_destripe=True
        destriping_statistics: statistic type which will be used to compute the line average, currently only median is implemented
        max_displacement: all displacement values whose absolute magnitude exceeds this value will be masked out for the correction
        xlim: x axis limits of the histogram plots
        ylim: y axis limits of the histogram plots

    Returns: the corrected grid, a dictionary with the log of the corrections and the correction grid

    """

    # DEBUG
    # displacement_field = master_slave_info_2009_2015['px1']
    # ref_path = master_slave_info_2009_2015['pan_master']
    # cc_field=None
    # cc_thresh=0.333
    # iterations = 3
    # zoom = 0.25
    # diagnose = True
    # full_resolution = True
    # do_destripe = False
    # get_correction_grid = True
    # max_displacement = -9999
    # xlim=[-5, -10]
    # ylim=[0, 3]

    # allocate output dictionary
    out_dict = {}

    # color map for diagnostic plots
    cmap = 'seismic'

    # this allows GDAL to throw Python Exceptions
    gdal.UseExceptions()

    # import offset
    displacement_field_link = gdal.Open(displacement_field)
    px_band = displacement_field_link.GetRasterBand(1)
    in_surface = np.array(px_band.ReadAsArray())

    # figure out the georeference of the displacement field
    projection = displacement_field_link.GetProjection()
    if projection == '':
        print('Input field has no georeference. Trying master image instead...')
        try:
            master_link = gdal.Open(ref_path)
            geotransform = master_link.GetGeoTransform()
        except:
            raise IOError('Cannot read georeference from ' + ref_path)
    else:
        geotransform = displacement_field_link.GetGeoTransform()

    # downsample
    if not zoom == 1:
        in_surface = cv2.resize(in_surface, None, fx=zoom, fy=zoom, interpolation=cv2.INTER_AREA)

    if cc_field is not None:
        print('Masking displacement field with correlation coefficients')

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
        in_surface = np.ma.masked_array(in_surface, mask=in_cc)

        # check if there are "enough" valid pixels
        if in_cc.sum() / in_cc.size > 0.5:
            warnings.warn("Less than 50% of the pixels are valid.", UserWarning)
            warnings.warn("The reason might be low correlation or small overlap among the input images.", UserWarning)

        del in_cc

    # mask large displacements
    in_surface = mask_large_values(in_surface, max_displacement)

    # get infos on the input
    val_mean = np.ma.mean(in_surface)
    val_std = np.ma.std(in_surface)
    val_min = val_mean - val_std
    val_max = val_mean + val_std
    out_dict['input'] = {}
    out_dict['input']['name'] = displacement_field
    out_dict['input']['mean'] = val_mean
    out_dict['input']['std'] = val_std
    out_dict['input']['rmse'] = np.sqrt(np.ma.mean(np.square(in_surface)))
    out_dict['input']['max'] = np.ma.max(in_surface)
    out_dict['input']['min'] = np.ma.min(in_surface)
    if isinstance(in_surface,np.ma.MaskedArray):
        out_dict['input']['n_samples'] = in_surface.count()
    else:
        out_dict['input']['n_samples'] = in_surface.size
    _, mean_gauss, std_gauss, rmse_gauss = misc.iterative_gaussian_fit(in_surface, diagnose=diagnose,
                                                                       title='Residuals before corrections',
                                                                       xlim=xlim, ylim=ylim)
    out_dict['input']['rmse_gauss'] = rmse_gauss
    out_dict['input']['std_gauss'] = std_gauss
    out_dict['input']['mean_gauss'] = mean_gauss

    # fit plane
    correction_surface, iter_mad = fit_plane_iter(in_surface, iterations=iterations)

    if diagnose:
        # plot
        vis_max = max(abs(val_min), abs(val_max))
        vis_min = -vis_max
        fig, axes = plt.subplots(ncols=3)
        axes[0].imshow(in_surface, cmap=cmap, vmin=vis_min, vmax=vis_max)
        axes[0].set_title('Input', fontsize=16)
        axes[0].tick_params(labelsize=14)
        plt.show()
        plt.draw()
        plt.pause(0.001)
    else:
        vis_max=0
    # correction with plane
    in_surface -= correction_surface

    if diagnose:
        # plot
        axes[1].imshow(correction_surface, cmap=cmap, vmin=vis_min, vmax=vis_max)
        axes[1].set_title('Fitted plane', fontsize=16)
        axes[1].tick_params(labelsize=14)
        im = axes[2].imshow(in_surface, cmap=cmap, vmin=vis_min, vmax=vis_max)
        axes[2].set_title('After plane correction', fontsize=16)
        axes[2].tick_params(labelsize=14)
        fig.subplots_adjust(right=0.91)
        cbar_ax = fig.add_axes([0.92, 0.35, 0.02, 0.3])
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.ax.tick_params(labelsize=16)

    # get infos after plane fitting
    out_dict['plane_corrected'] = {}
    out_dict['plane_corrected']['iter_mad'] = iter_mad
    out_dict['plane_corrected']['mean'] = np.ma.mean(in_surface)
    out_dict['plane_corrected']['std'] = np.ma.std(in_surface)
    out_dict['plane_corrected']['rmse'] = np.sqrt(np.ma.mean(np.square(in_surface)))
    out_dict['plane_corrected']['max'] = np.ma.max(in_surface)
    out_dict['plane_corrected']['min'] = np.ma.min(in_surface)
    _, mean_gauss, std_gauss, rmse_gauss = misc.iterative_gaussian_fit(in_surface, diagnose=diagnose,
                                                                       title='Residuals after deramping',
                                                                       xlim=xlim, ylim=ylim)
    out_dict['plane_corrected']['rmse_gauss'] = rmse_gauss
    out_dict['plane_corrected']['std_gauss'] = std_gauss
    out_dict['plane_corrected']['mean_gauss'] = mean_gauss

    if do_destripe:

        if master_slave['master']['satellite'] == master_slave['slave']['satellite'] == 'S2':
            raise Exception('...for co-registration S-2 to S-2 call get_correction_plane_S2.')

        if master_slave['master']['satellite'] == master_slave['slave']['satellite'] == 'L8':
            raise Exception('...Co-registration L-8 to L-8 still needs to be implemented.')

        if master_slave['master']['satellite'] == 'L8' and master_slave['slave']['satellite'] == 'S2':
            print('...Generating corrections bands for S-2 to L-8 scenario...')
            # get detector footprints
            detector_footprint = footprints.get_detector_footprint(master_slave['slave']['azimuth'], diagn=diagnose)

        if master_slave['master']['satellite'] == 'S2' and master_slave['slave']['satellite'] == 'L8':
            print('...Generating corrections bands for L-8 to S-2 scenario...')
            # get detector footprints
            detector_footprint = footprints.get_detector_footprint(master_slave['master']['azimuth'], diagn=diagnose)

        # no intersection needed
        detector_intersections_path = None
        detector_offsets = None
        detector_stds = None

        # get angle
        angle, angle_std = footprints.getAngle([detector_footprint])

        # angle
        in_surface, destripe_grid = destripe.destripe(in_surface,
                                                       angle=angle,
                                                       horizontal_correction=horizontal_correction,
                                                       diagnose=diagnose,
                                                       cmap=cmap,
                                                       color_max=vis_max)

        correction_surface += destripe_grid

        # get infos after along track destriping
        out_dict['destripe'] = {}
        out_dict['destripe']['detfoo1'] = master_slave['master']['azimuth']
        out_dict['destripe']['detfoo2'] = master_slave['slave']['azimuth']
        out_dict['destripe']['detfoo_intersect'] = detector_intersections_path
        out_dict['destripe']['detfoo_means'] = detector_offsets
        out_dict['destripe']['detfoo_stds'] = detector_stds
        out_dict['destripe']['mean'] = np.ma.mean(in_surface)
        out_dict['destripe']['std'] = np.ma.std(in_surface)
        out_dict['destripe']['rmse'] = np.sqrt(np.ma.mean(np.square(in_surface)))
        out_dict['destripe']['max'] = np.ma.max(in_surface)
        out_dict['destripe']['min'] = np.ma.min(in_surface)
        _, mean_gauss, std_gauss, rmse_gauss = misc.iterative_gaussian_fit(in_surface, diagnose=diagnose,
                                                                           title='Residuals after destriping',
                                                                           xlim=xlim, ylim=ylim)
        out_dict['destripe']['rmse_gauss'] = rmse_gauss
        out_dict['destripe']['std_gauss'] = std_gauss
        out_dict['destripe']['mean_gauss'] = mean_gauss

        del destripe_grid
        gc.collect()

    # if the correction grid is computed at lower resolution but to be applied on the full resolution input
    if full_resolution:
        if not zoom == 1:
            del in_surface
            gc.collect()
            in_surface, correction_surface = get_full_resolution_output(displacement_field,
                                                                        correction_surface,
                                                                        cc_field,
                                                                        mask_full_resolution_output,
                                                                        cc_thresh,
                                                                        max_displacement)

    return in_surface, out_dict, correction_surface


def get_correction_plane_S2(displacement_field,
                            master_slave,
                            ref_path,
                            cc_field=None,
                            current_work_folder=None,
                            cc_thresh=0.333,
                            iterations=3,
                            zoom=1.0,
                            diagnose=False,
                            full_resolution=False,
                            mask_full_resolution_output=True,
                            do_destripe=False,
                            horizontal_correction=False,
                            max_detector_offset=None,
                            destriping_statistics='mean',
                            max_displacement=-9999,
                            xlim=[-1, 1], ylim=[0, 3]):


    """
    Approximate the displacement field with a plane, destripe and return the corrected grid, a dictionary with the log of
    the corrections and the correction grid that can be used to transform the slave image.
    :param displacement_field:  full path to the displacement fields that should be corrected
    :param master_slave:        master slave dictionary
    :param ref_path:            full path to a reference image, if the displacement field has no georeference the reference
                                of the reference image will be used
    :param cc_field:                full path to the raster holding the correlation coefficients
    :param current_work_folder:     work folder where detector shapefiles will be written, if not specified they will be
                                    written to os.path.dirname(os.path.dirname(Px1))
    :param cc_thresh:               threshold on the correlation coefficient, pixels below this value will be ignored during the
                                    correction
    :param iterations:              number of iterations for iteratively reweighted plane fitting
    :param zoom:                    optional downsampling, mainly for reducing memory size and processing time
    :param diagnose:                if True plotting is activated to illustrate the corrections
    :param full_resolution:         if True the final grid will be returned at the image resolution
    :param mask_full_resolution_output: if set False the masks (cc, max_displacement, will not be applied on the output)
    :param do_destripe:             if True a destriping routine will applied
    :param horizontal_correction:   if True across-track destriping will be performed
    :param max_detector_offset:     optional threshold for maximum values which will be considered when computing detector offsets
    :param destriping_statistics:   statistic type which will be used to compute the detector average, options are 'mean', 'median', 'mode'
    :param max_displacement:        all displacement values whose absolute magnitude exceeds this value will be masked out for the correction
    :param xlim                     x axis limits of the histogram plots
    :param ylim                     y axis limits of the histogram plots
    :return:                        the corrected grid, a dictionary with the log of the corrections and the correction grid

    """

    # allocate output dictionary
    out_dict = {}

    cmap = 'seismic'

    # this allows GDAL to throw Python Exceptions
    gdal.UseExceptions()

    # import offset
    displacement_field_link = gdal.Open(displacement_field)
    px_band = displacement_field_link.GetRasterBand(1)
    in_surface = np.array(px_band.ReadAsArray())

    # figure out the georeference of the displacement field
    projection = displacement_field_link.GetProjection()
    if projection == '':
        print('Input field has no georeference. Trying reference image instead...')
        try:
            master_link = gdal.Open(ref_path)
            geotransform = master_link.GetGeoTransform()
        except:
            raise IOError('Cannot read georeference from ' + ref_path)
    else:
        geotransform = displacement_field_link.GetGeoTransform()

    # downsample
    if not zoom == 1:
        in_surface = cv2.resize(in_surface, None, fx=zoom, fy=zoom, interpolation=cv2.INTER_AREA)

    if cc_field is not None:
        print('Masking displacement field with correlation coefficients')

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
        in_surface = np.ma.masked_array(in_surface, mask=in_cc)

        # check if there are "enough" valid pixels
        if in_cc.sum() / in_cc.size > 0.5:
            warnings.warn("Less than 50% of the pixels are valid.", UserWarning)
            warnings.warn("The reason might be low correlation or small overlap among the input images.", UserWarning)

        del in_cc

    # mask large displacements
    in_surface = mask_large_values(in_surface, max_displacement)

    # get infos on the input
    val_mean = np.ma.mean(in_surface)
    val_std = np.ma.std(in_surface)
    val_min = val_mean - val_std
    val_max = val_mean + val_std
    out_dict['input'] = {}
    out_dict['input']['name'] = displacement_field
    out_dict['input']['mean'] = val_mean
    out_dict['input']['std'] = val_std
    out_dict['input']['rmse'] = np.sqrt(np.ma.mean(np.square(in_surface)))
    out_dict['input']['max'] = np.ma.max(in_surface)
    out_dict['input']['min'] = np.ma.min(in_surface)
    if isinstance(in_surface,np.ma.MaskedArray):
        out_dict['input']['n_samples'] = in_surface.count()
    else:
        out_dict['input']['n_samples'] = in_surface.size
    _, mean_gauss, std_gauss, rmse_gauss = misc.iterative_gaussian_fit(in_surface, diagnose=diagnose,
                                                                       title='Residuals before corrections',
                                                                       xlim=xlim, ylim=ylim)
    out_dict['input']['rmse_gauss'] = rmse_gauss
    out_dict['input']['std_gauss'] = std_gauss
    out_dict['input']['mean_gauss'] = mean_gauss

    # fit plane
    correction_surface, iter_mad = fit_plane_iter(in_surface, iterations=iterations)

    if diagnose:
        # plot
        fig, axes = plt.subplots(ncols=3)
        vis_max = max(abs(val_min), abs(val_max))
        vis_min = -vis_max
        axes[0].imshow(in_surface, cmap=cmap, vmin=vis_min, vmax=vis_max)
        axes[0].set_title('Input', fontsize=16)
        axes[0].tick_params(labelsize=14)
        plt.show()
        plt.draw()
        plt.pause(0.001)

    # correction with plane
    in_surface -= correction_surface

    if diagnose:
        # plot
        axes[1].imshow(correction_surface, cmap=cmap, vmin=vis_min, vmax=vis_max)
        axes[1].set_title('Fitted plane', fontsize=16)
        axes[1].tick_params(labelsize=14)
        im = axes[2].imshow(in_surface, cmap=cmap, vmin=vis_min, vmax=vis_max)
        axes[2].set_title('After plane correction', fontsize=16)
        axes[2].tick_params(labelsize=14)
        fig.subplots_adjust(right=0.91)
        cbar_ax = fig.add_axes([0.92, 0.35, 0.02, 0.3])
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.ax.tick_params(labelsize=16)

    # get infos after plane fitting
    out_dict['plane_corrected'] = {}
    out_dict['plane_corrected']['iter_mad'] = iter_mad
    out_dict['plane_corrected']['mean'] = np.ma.mean(in_surface)
    out_dict['plane_corrected']['std'] = np.ma.std(in_surface)
    out_dict['plane_corrected']['rmse'] = np.sqrt(np.ma.mean(np.square(in_surface)))
    out_dict['plane_corrected']['max'] = np.ma.max(in_surface)
    out_dict['plane_corrected']['min'] = np.ma.min(in_surface)
    _, mean_gauss, std_gauss, rmse_gauss = misc.iterative_gaussian_fit(in_surface, diagnose=diagnose,
                                                                       title='Residuals after deramping',
                                                                       xlim=xlim, ylim=ylim)
    out_dict['plane_corrected']['rmse_gauss'] = rmse_gauss
    out_dict['plane_corrected']['std_gauss'] = std_gauss
    out_dict['plane_corrected']['mean_gauss'] = mean_gauss

    if do_destripe:

        # get detector footprint intersections, includes check if bother master and slave are S2 images
        if current_work_folder is None:
            current_work_folder = os.path.dirname(os.path.dirname(displacement_field))
        detector_intersections_path = footprints.get_S2_detector_intersections(master_slave,
                                                                               current_work_folder,
                                                                               diagnose=diagnose)
        destripe_grid = np.zeros_like(in_surface)

        # adjust geotransform according to zoom resolution
        gt_adj = tuple([geotransform[0], geotransform[1] / zoom, geotransform[2], geotransform[3], geotransform[4],
                        geotransform[5] / zoom])

        ds_detectors = ogr.Open(detector_intersections_path)
        layer_detectors = ds_detectors.GetLayer()
        layer_detectors.GetFeatureCount()

        detector_stds = []
        detector_offsets = []
        count = 0

        # layer_detectors.ResetReading()
        for feat in layer_detectors:

            # feat = layer_detectors.GetNextFeature()
            count += 1
            print('Computing average of detector element ' + str(count))

            # get mask from geometry
            geom = feat.GetGeometryRef()
            output_mask, ulX, ulY, gt2 = mask.geom2mask(in_surface, geom, gt=gt_adj)

            # get values
            val = in_surface[output_mask == 0]
            if max_detector_offset is not None:  # mask values exceeding the expected detector offset pattern
                val = np.ma.masked_greater(val, max_detector_offset)
                val = np.ma.masked_less(val, -max_detector_offset)

            # compute detector average
            if destriping_statistics == 'mode':
                val_round = np.round(val, 3)  # escape the floating point hell
                detector_average = stats.mode(val_round, axis=None)
                detector_average = detector_average[0]
            elif destriping_statistics == 'median':
                detector_average = np.ma.median(val)
            elif destriping_statistics == 'mean':
                detector_average = np.ma.mean(val)
            else:
                raise Exception(destriping_statistics + ' is not an accepted option. Use mode, mean or median')

            # compute standard deviation of the detector offset
            detector_stds.append(np.ma.std(val))

            # assign detector average
            destripe_grid[np.where(output_mask == 0)] = detector_average
            detector_offsets.append(detector_average)

        if diagnose:
            fig, axes = plt.subplots(ncols=3)

            axes[0].imshow(in_surface, cmap=cmap, vmin=vis_min, vmax=vis_max)
            axes[0].set_title('Input', fontsize=16)
            axes[0].tick_params(labelsize=14)
            plt.show()
            plt.draw()
            plt.pause(0.001)

            in_surface -= destripe_grid

            # plot
            axes[1].imshow(destripe_grid, cmap=cmap, vmin=vis_min, vmax=vis_max)
            axes[1].set_title('Destriping grid', fontsize=16)
            axes[1].tick_params(labelsize=14)
            im = axes[2].imshow(in_surface, cmap=cmap, vmin=vis_min, vmax=vis_max)
            axes[2].set_title('After destriping', fontsize=16)
            axes[2].tick_params(labelsize=14)
            fig.subplots_adjust(right=0.91)
            cbar_ax = fig.add_axes([0.92, 0.35, 0.02, 0.3])
            cb = fig.colorbar(im, cax=cbar_ax)
            cb.ax.tick_params(labelsize=16)

        correction_surface += destripe_grid

        del destripe_grid
        gc.collect()

        # get infos after along track destriping
        out_dict['destripe'] = {}
        out_dict['destripe']['detfoo1'] = master_slave['master']['azimuth']
        out_dict['destripe']['detfoo2'] = master_slave['slave']['azimuth']
        out_dict['destripe']['detfoo_intersect'] = detector_intersections_path
        out_dict['destripe']['detfoo_means'] = detector_offsets
        out_dict['destripe']['detfoo_stds'] = detector_stds
        out_dict['destripe']['mean'] = np.ma.mean(in_surface)
        out_dict['destripe']['std'] = np.ma.std(in_surface)
        out_dict['destripe']['rmse'] = np.sqrt(np.ma.mean(np.square(in_surface)))
        out_dict['destripe']['max'] = np.ma.max(in_surface)
        out_dict['destripe']['min'] = np.ma.min(in_surface)
        _, mean_gauss, std_gauss, rmse_gauss = misc.iterative_gaussian_fit(in_surface, diagnose=diagnose,
                                                                           title='Residuals after along-track destriping',
                                                                           xlim=xlim, ylim=ylim)
        out_dict['destripe']['rmse_gauss'] = rmse_gauss
        out_dict['destripe']['std_gauss'] = std_gauss
        out_dict['destripe']['mean_gauss'] = mean_gauss

        if horizontal_correction:

            # transform shapes to polygon lists
            ds = ogr.Open(detector_intersections_path)
            layer = ds.GetLayer()
            n_detectors = layer.GetFeatureCount()
            detector_ids = list()
            detector_footprint2 = []
            for i in range(0, n_detectors):
                feat = layer.GetFeature(i)
                detector_ids.append(list(feat.items().values()))
                geom = feat.GetGeometryRef()
                geom.GetGeometryName()
                poly = geom.GetGeometryRef(0)
                coords = np.asarray(poly.GetPoints())
                detector_footprint2.append(coords)

            # get angle
            angle, angle_std = footprints.getAngle(detector_footprint2, tol=0.5)

            # rotate
            if hasattr(in_surface, 'mask'):
                in_surface_rot, orig, mask_rot = misc.rotateImage(in_surface, angle)
            else:
                in_surface_rot, orig = misc.rotateImage(in_surface, angle)

                # apply mask
            if 'mask_rot' in locals():
                in_surface_rot = np.ma.masked_where(mask_rot > 0, in_surface_rot)

            # remove pixels not covered in the original image
            in_surface_rot = np.ma.masked_where(in_surface_rot == -255, in_surface_rot)

            if diagnose:
                fig = plt.figure()
                fig.suptitle('Rotated to along track direction', fontsize=16)
                plt.imshow(in_surface_rot)
                step_size = 100
                for i in np.arange(step_size, in_surface_rot.shape[1], step_size):
                    plt.axvline(i)

            # compute average of each profile
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                if destriping_statistics == 'mode':
                    val_round = np.round(val, 3)  # escape the floating point hell
                    row_average = stats.mode(val_round, axis=1)
                    detector_average = detector_average[0]
                elif destriping_statistics == 'median':
                    row_average = np.ma.median(in_surface_rot, axis=1)
                elif destriping_statistics == 'mean':
                    row_average = np.ma.mean(in_surface_rot, axis=1)
                else:
                    raise Exception(destriping_statistics + ' is not an accepted option. Use mode, mean or median')

            # generate and rotate stripe correction grid
            destripe_grid_horizontal = np.zeros((in_surface_rot.shape[0], in_surface_rot.shape[1]))
            for i in range(in_surface_rot.shape[1]):
                destripe_grid_horizontal[:, i] = row_average

            rows_in, cols_in = destripe_grid_horizontal.shape
            M = cv2.getRotationMatrix2D((cols_in / 2, rows_in / 2), -angle, 1)
            rows_out, cols_out = in_surface.shape
            M[0, 2] += (rows_out - rows_in) / 2
            M[1, 2] += (cols_out - cols_in) / 2
            destripe_grid_horizontal = cv2.warpAffine(destripe_grid_horizontal, M, (cols_out, rows_out),
                                                      flags=cv2.INTER_NEAREST)

            if diagnose:
                val_mean = np.mean(in_surface)
                val_std = np.std(in_surface)
                val_min = val_mean - val_std
                val_max = val_mean + val_std
                vis_max = max(abs(val_min), abs(val_max))
                vis_min = -vis_max
                fig, axes = plt.subplots(ncols=3)
                axes[0].imshow(in_surface, cmap=cmap, vmin=vis_min, vmax=vis_max)
                axes[0].set_title('Before across-track destriping', fontsize=16)
                axes[0].tick_params(labelsize=14)

            # destripe horizontal
            in_surface -= destripe_grid_horizontal

            if diagnose:
                axes[1].imshow(destripe_grid_horizontal, cmap=cmap, vmin=vis_min, vmax=vis_max)
                axes[1].set_title('Across-track destriping grid', fontsize=16)
                axes[1].tick_params(labelsize=14)
                im = axes[2].imshow(in_surface, cmap=cmap, vmin=vis_min, vmax=vis_max)
                axes[2].set_title('After across-track destriping', fontsize=16)
                axes[2].tick_params(labelsize=14)
                fig.subplots_adjust(right=0.91)
                cbar_ax = fig.add_axes([0.92, 0.35, 0.02, 0.3])
                cb = fig.colorbar(im, cax=cbar_ax)
                cb.ax.tick_params(labelsize=16)

            correction_surface += destripe_grid_horizontal

            # get infos after across-track destriping
            out_dict['destripe_act'] = {}
            out_dict['destripe_act']['mean'] = np.ma.mean(in_surface)
            out_dict['destripe_act']['std'] = np.ma.std(in_surface)
            out_dict['destripe_act']['rmse'] = np.sqrt(np.ma.mean(np.square(in_surface)))
            out_dict['destripe_act']['max'] = np.ma.max(in_surface)
            out_dict['destripe_act']['min'] = np.ma.min(in_surface)
            _, mean_gauss, std_gauss, rmse_gauss = misc.iterative_gaussian_fit(in_surface, diagnose=diagnose,
                                                                               title='Residuals after across-track destriping')
            out_dict['destripe_act']['rmse_gauss'] = rmse_gauss
            out_dict['destripe_act']['std_gauss'] = std_gauss
            out_dict['destripe_act']['mean_gauss'] = mean_gauss

    else:
        out_dict['destripe'] = 'no destriping performed'

    # if the correction grid is computed at lower resolution but to be applied on the full resolution input
    if full_resolution:
        if not zoom == 1:
            del in_surface
            gc.collect()
            in_surface, correction_surface = get_full_resolution_output(displacement_field,
                                                                        correction_surface,
                                                                        cc_field,
                                                                        mask_full_resolution_output,
                                                                        cc_thresh,
                                                                        max_displacement)

    return in_surface, out_dict, correction_surface


def get_correction(master_slave_info, iterations=3, cc_thresh=0.33, zoom=0.25, diagnose=False, full_resolution=True,
                   mask_full_resolution_output=True, do_destripe=False, horizontal_correction=False,
                   save_correction_grids=True, xlim=[-1, 1], ylim=[0, 3]):

    """
    Wrapper to get_correction_plane functions, writes correction grids and corrected grids to disk, the output folder is
    master_slave_info['correl_folder']
    :param master_slave_info: master_slave_info: master slave dictionary
    :param iterations: iterations: number of iterations for iteratively reweighted plane fitting
    :param cc_thresh: threshold on the correlation coefficient, pixels below this value will be ignored during the
                      correction
    :param zoom: optional downsampling, mainly for reducing memory size and processing time
    :param diagnose: if True plotting is activated to illustrate the corrections
    :param full_resolution:  if True the final grid will be returned at the image resolution
    :param do_destripe: if True a destriping routine will applied
    :param horizontal_correction: if True across-track destriping will be performed
    :param save_correction_grids: if True the correction grids will be stored to disk in GeoTiff format
    :param xlim: x axis limits of the histogram plots
    :param ylim: y axis limits of the histogram plots
    :return: master slave dictionary including a log of the corrections and full paths to the output files on disk
             being the correction grids, the corrected grids and the georeferenced correlation grid 
    """

    # init new dictionary entries
    master_slave_info['corrections'] = {}
    master_slave_info['out'] = {}

    # compute input RMSE xy
    rmse_xy_before = misc.compute_rmse_xy(master_slave_info['px1'], master_slave_info['px2'], master_slave_info['cc'],
                                          cc_thresh=cc_thresh, zoom=zoom, diagnose=diagnose,
                                          titles=['Residuals X before correction', 'Residuals Y before correction'],
                                          xlim=xlim, ylim=ylim)
    master_slave_info['corrections']['rmse_xy_input'] = rmse_xy_before

    # correct X component
    print('Correcting X component ' + master_slave_info['px1'])
    if master_slave_info['master']['satellite'] == master_slave_info['slave']['satellite'] == 'S2':
        out_surface_px1, out_dict_px1, px1_correction_grid = get_correction_plane_S2(master_slave_info['px1'],
                                                                                     master_slave_info,
                                                                                     master_slave_info['pan_master'],
                                                                                     cc_field=master_slave_info['cc'],
                                                                                     cc_thresh=cc_thresh,
                                                                                     iterations=iterations,
                                                                                     zoom=zoom,
                                                                                     diagnose=diagnose,
                                                                                     full_resolution=full_resolution,
                                                                                     mask_full_resolution_output=mask_full_resolution_output,
                                                                                     do_destripe=do_destripe,
                                                                                     horizontal_correction=horizontal_correction,
                                                                                     xlim=xlim, ylim=ylim)

    if not master_slave_info['master']['satellite'] == master_slave_info['slave']['satellite']:
        out_surface_px1, out_dict_px1, px1_correction_grid = get_correction_plane_ms(master_slave_info['px1'],
                                                                                     master_slave_info,
                                                                                     master_slave_info['pan_master'],
                                                                                     cc_field=master_slave_info['cc'],
                                                                                     cc_thresh=cc_thresh,
                                                                                     iterations=iterations,
                                                                                     zoom=zoom,
                                                                                     diagnose=diagnose,
                                                                                     full_resolution=full_resolution,
                                                                                     mask_full_resolution_output=mask_full_resolution_output,
                                                                                     do_destripe=do_destripe,
                                                                                     horizontal_correction=horizontal_correction,
                                                                                     xlim=xlim, ylim=ylim)
    if save_correction_grids:
        # save correction grid for X to disk
        dst_filename_px1_cg = os.path.join(master_slave_info['correl_folder'], 'px1_correction_grid.tif')
        array2geotiff.array2geotiff(px1_correction_grid, dst_filename_px1_cg, master_slave_info['pan_master'])
        master_slave_info['out']['correction_grid_px1'] = dst_filename_px1_cg
    del px1_correction_grid
    gc.collect()

    # correct Y component
    print('Correcting Y component ' + master_slave_info['px2'])
    if master_slave_info['master']['satellite'] == master_slave_info['slave']['satellite'] == 'S2':
        out_surface_px2, out_dict_px2, px2_correction_grid = get_correction_plane_S2(master_slave_info['px2'],
                                                                                     master_slave_info,
                                                                                     master_slave_info['pan_master'],
                                                                                     cc_field=master_slave_info['cc'],
                                                                                     cc_thresh=cc_thresh,
                                                                                     iterations=iterations,
                                                                                     zoom=zoom,
                                                                                     diagnose=diagnose,
                                                                                     full_resolution=full_resolution,
                                                                                     mask_full_resolution_output=mask_full_resolution_output,
                                                                                     do_destripe=do_destripe,
                                                                                     horizontal_correction=horizontal_correction,
                                                                                     xlim=xlim, ylim=ylim)

    if not master_slave_info['master']['satellite'] == master_slave_info['slave']['satellite']:
        out_surface_px2, out_dict_px2, px2_correction_grid = get_correction_plane_ms(master_slave_info['px2'],
                                                                                     master_slave_info,
                                                                                     master_slave_info['pan_master'],
                                                                                     cc_field=master_slave_info['cc'],
                                                                                     cc_thresh=cc_thresh,
                                                                                     iterations=iterations,
                                                                                     zoom=zoom,
                                                                                     diagnose=diagnose,
                                                                                     full_resolution=full_resolution,
                                                                                     mask_full_resolution_output=mask_full_resolution_output,
                                                                                     do_destripe=do_destripe,
                                                                                     horizontal_correction=horizontal_correction,
                                                                                     xlim=xlim, ylim=ylim)
    if save_correction_grids:
        dst_filename_px2_cg = os.path.join(master_slave_info['correl_folder'], 'px2_correction_grid.tif')
        array2geotiff.array2geotiff(px2_correction_grid, dst_filename_px2_cg, master_slave_info['pan_master'])
        master_slave_info['out']['correction_grid_px2'] = dst_filename_px2_cg
    del px2_correction_grid
    gc.collect()

    # save corrected grids to disk
    dst_filename_px1 = os.path.join(master_slave_info['correl_folder'], 'px1_corrected_georef.tif')
    array2geotiff.array2geotiff(out_surface_px1, dst_filename_px1, master_slave_info['pan_master'])
    dst_filename_px2 = os.path.join(master_slave_info['correl_folder'], 'px2_corrected_georef.tif')
    array2geotiff.array2geotiff(out_surface_px2, dst_filename_px2, master_slave_info['pan_master'])
    misc.add_georef(master_slave_info['cc'], master_slave_info['pan_master'])
    master_slave_info['out']['px1_corrected'] = dst_filename_px1
    master_slave_info['out']['px2_corrected'] = dst_filename_px2
    master_slave_info['out']['cc'] = master_slave_info['cc']

    # remove outliers to compute output rmse xy
    outlier_px1 = out_dict_px1['destripe']['std_gauss'] * 2.575
    outlier_px2 = out_dict_px2['destripe']['std_gauss'] * 2.575
    out_surface_px1 = np.ma.masked_where(out_surface_px1 >= outlier_px1, out_surface_px1)
    out_surface_px1 = np.ma.masked_where(out_surface_px1 <= -outlier_px1, out_surface_px1)
    out_surface_px2 = np.ma.masked_where(out_surface_px2 >= outlier_px2, out_surface_px2)
    out_surface_px2 = np.ma.masked_where(out_surface_px2 <= -outlier_px2, out_surface_px2)
    rmse_xy_after = np.ma.sqrt(np.ma.mean(np.square(out_surface_px1) + np.square(out_surface_px2)))
    master_slave_info['corrections']['rmse_xy_output'] = rmse_xy_after
    master_slave_info['corrections']['px1_log'] = out_dict_px1
    master_slave_info['corrections']['px2_log'] = out_dict_px1

    plt.close("all")

    return master_slave_info


def correct_all(master_slave_info_list,
                iterations=3,
                diagnose=True,
                zoom=0.25,
                cc_thresh=0.33,
                full_resolution=True,
                do_destripe=True,
                horizontal_correction=False,
                xlim=[-1, 1], ylim=[0, 3]):

    """
    Wrapper to get_correction, warp_correctiongrid and apply_correction to process a list of master slave dictonaries
    Args:
        master_slave_info_list: list of master slave dictionaries
        iterations: iterations: number of iterations for iteratively reweighted plane fitting
        diagnose: if True plotting is activated to illustrate the corrections
        zoom: optional downsampling, mainly for reducing memory size and processing time
        cc_thresh: threshold on the correlation coefficient, pixels below this value will be ignored during the
                   correction
        full_resolution: if True the final grid will be returned at the image resolution
        do_destripe: if True a destriping routine will applied
        horizontal_correction: if True across-track destriping will be performed
        xlim: x axis limits of the histogram plots
        ylim: y axis limits of the histogram plots

    Returns: master slave dictionary including a log of the corrections and full paths to the output files on disk
             being the correction grids, the corrected grids, the georeferenced correlation grid and the corrected
             image bands

    """

    master_slave_info_list_updated = []

    for master_slave_info in master_slave_info_list:

        # run estimation of correction parameters
        master_slave_info = get_correction(master_slave_info,
                                           iterations=iterations,
                                           cc_thresh=cc_thresh,
                                           zoom=zoom,
                                           diagnose=diagnose,
                                           full_resolution=full_resolution,
                                           do_destripe=do_destripe,
                                           horizontal_correction=horizontal_correction,
                                           xlim=xlim, ylim=ylim)

        # warp correction grids to grid of the slave granule
        master_slave_info['out']['warp_grid_dict'] = correct_image.warp_correctiongrid(
            master_slave_info['out']['correction_grid_px1'],
            master_slave_info['out']['correction_grid_px2'],
            master_slave_info['slave']['band_paths'],
            satellite=master_slave_info['slave']['satellite'],
            resampling_method='cubic')

        # apply corrections to slave
        corrected_bands = correct_image.apply_correction(master_slave_info['slave']['band_paths'],
                                                         master_slave_info['out']['warp_grid_dict'],
                                                         master_slave_info['slave']['satellite'])

        master_slave_info['out']['corrected_bands'] = corrected_bands

        master_slave_info_list_updated.append(master_slave_info)

    return master_slave_info_list_updated


def deramp(displacement_field,
           cc_field=None,
           cc_thresh=0.333,
           iterations=3,
           zoom=1.0,
           diagnose=False,
           full_resolution=False,
           mask_full_resolution_output=True,
           get_correction_grid=False,
           max_displacement=-9999,
           xlim=[-1, 1], ylim=[0, 3]):
    """
    Basically a sipmlified version of get_correction_plane_S2 which takes any offset field and performs
    a simple deramping.
    Args:
        displacement_field: full path to the displacement fields that should be corrected
        cc_field:           full path to the raster holding the correlation coefficients
        cc_thresh:          threshold for masking values according to the correlation coefficient, pixels below this value
                            will be ignored during the correction
        iterations:         number of iterations for iteratively reweighted plane fitting
        zoom:               optional downsampling, mainly for reducing memory size
        diagnose:           if true plotting is activated to illustrate the corrections
        full_resolution:    upsample correction grid to original resolution before correction
        mask_full_resolution_output: if set False the masks (cc, max_displacement, will not be applied on the output)
        get_correction_grid: get correction grid in addition to the corrected surface, upsampled to full resolution
        max_displacement: all displacement values whose absolute magnitude exceeds this value will be masked out for the correction
        xlim: x axis limits of the histogram plots
        ylim: y axis limits of the histogram plots

    Returns: corrected grid, a dictionary with the log of the corrections and optionally the correction grid

# DEBUG
displacement_field = master_slave_info_2009_2015['px1']
cc_field=master_slave_info_2009_2015['cc']
cc_thresh=0.7
iterations=3
zoom=0.25
diagnose=True
full_resolution=True
mask_full_resolution_output=True
get_correction_grid=True
max_displacement=-9999
xlim=[-1, 1]
ylim=[0, 3]


    """

    # allocate output dictionary
    out_dict = {}

    # color map for diagnostic plots
    cmap = 'seismic'

    # this allows GDAL to throw Python Exceptions
    gdal.UseExceptions()

    # import offset
    displacement_field_link = gdal.Open(displacement_field)
    px_band = displacement_field_link.GetRasterBand(1)
    in_surface = np.array(px_band.ReadAsArray())

    # figure out the georeference of the displacement field
    projection = displacement_field_link.GetProjection()

    # downsample
    if not zoom == 1:
        in_surface = cv2.resize(in_surface, None, fx=zoom, fy=zoom, interpolation=cv2.INTER_AREA)

    if cc_field is not None:
        print('Masking displacement field with correlation coefficients')

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
        in_surface = np.ma.masked_array(in_surface, mask=in_cc)

        # check if there are "enough" valid pixels
        if in_cc.sum() / in_cc.size > 0.5:
            warnings.warn("Less than 50% of the pixels are valid.", UserWarning)
            warnings.warn("The reason might be low correlation or small overlap among the input images.", UserWarning)

        del in_cc

    # mask large displacements
    in_surface = mask_large_values(in_surface, max_displacement)

    # get infos on the input
    val_mean = np.ma.mean(in_surface)
    val_std = np.ma.std(in_surface)
    val_min = val_mean - val_std
    val_max = val_mean + val_std
    out_dict['input'] = {}
    out_dict['input']['name'] = displacement_field
    out_dict['input']['mean'] = val_mean
    out_dict['input']['std'] = val_std
    out_dict['input']['rmse'] = np.sqrt(np.ma.mean(np.square(in_surface)))
    out_dict['input']['max'] = np.ma.max(in_surface)
    out_dict['input']['min'] = np.ma.min(in_surface)
    if isinstance(in_surface,np.ma.MaskedArray):
        out_dict['input']['n_samples'] = in_surface.count()
    else:
        out_dict['input']['n_samples'] = in_surface.size
    _, mean_gauss, std_gauss, rmse_gauss = misc.iterative_gaussian_fit(in_surface, diagnose=diagnose,
                                                                       title='Residuals before corrections',
                                                                       xlim=xlim, ylim=ylim)
    out_dict['input']['rmse_gauss'] = rmse_gauss
    out_dict['input']['std_gauss'] = std_gauss
    out_dict['input']['mean_gauss'] = mean_gauss

    # fit plane
    correction_surface, iter_mad = fit_plane_iter(in_surface, iterations=iterations)

    if diagnose:
        # plot
        vis_max = max(abs(val_min), abs(val_max))
        vis_min = -vis_max
        fig, axes = plt.subplots(ncols=3)
        axes[0].imshow(in_surface, cmap=cmap, vmin=vis_min, vmax=vis_max)
        axes[0].set_title('Input', fontsize=16)
        axes[0].tick_params(labelsize=14)
        plt.show()
        plt.draw()
        plt.pause(0.001)

    # correction with plane
    in_surface -= correction_surface

    if diagnose:
        # plot
        axes[1].imshow(correction_surface, cmap=cmap, vmin=vis_min, vmax=vis_max)
        axes[1].set_title('Fitted plane', fontsize=16)
        axes[1].tick_params(labelsize=14)
        im = axes[2].imshow(in_surface, cmap=cmap, vmin=vis_min, vmax=vis_max)
        axes[2].set_title('After plane correction', fontsize=16)
        axes[2].tick_params(labelsize=14)
        fig.subplots_adjust(right=0.91)
        cbar_ax = fig.add_axes([0.92, 0.35, 0.02, 0.3])
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.ax.tick_params(labelsize=16)

    # get infos after plane fitting
    out_dict['plane_corrected'] = {}
    out_dict['plane_corrected']['iter_mad'] = iter_mad
    out_dict['plane_corrected']['mean'] = np.ma.mean(in_surface)
    out_dict['plane_corrected']['std'] = np.ma.std(in_surface)
    out_dict['plane_corrected']['rmse'] = np.sqrt(np.ma.mean(np.square(in_surface)))
    out_dict['plane_corrected']['max'] = np.ma.max(in_surface)
    out_dict['plane_corrected']['min'] = np.ma.min(in_surface)
    _, mean_gauss, std_gauss, rmse_gauss = misc.iterative_gaussian_fit(in_surface, diagnose=diagnose,
                                                                       title='Residuals after deramping',
                                                                       xlim=xlim, ylim=ylim)
    out_dict['plane_corrected']['rmse_gauss'] = rmse_gauss
    out_dict['plane_corrected']['std_gauss'] = std_gauss
    out_dict['plane_corrected']['mean_gauss'] = mean_gauss

    # if the correction grid is computed at lower resolution but to be applied on the full resolution input
    if full_resolution:
        if not zoom == 1:
            del in_surface
            gc.collect()
            in_surface, correction_surface = get_full_resolution_output(displacement_field,
                                                                        correction_surface,
                                                                        cc_field,
                                                                        mask_full_resolution_output,
                                                                        cc_thresh,
                                                                        max_displacement)

    if get_correction_grid:

        return in_surface, out_dict, correction_surface

    else:

        return in_surface, out_dict
