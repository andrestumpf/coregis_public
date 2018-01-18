import numpy as np
import matplotlib.pyplot as plt
import cv2
import warnings



def getTranslationMatrix2d(dx, dy):
    """
    Returns a numpy affine transformation matrix for a 2D translation of
    (dx, dy)
    """
    return np.matrix([[1, 0, dx], [0, 1, dy], [0, 0, 1]])


def rotateImage(image,
                angle):
    """
    :param image: 'input image'
    :param angle: 'rotation angle'
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


def destripe(in_surface,
             angle,
             in_mask=None,
             horizontal_correction=False,
             diagnose=False,
             cmap='seismic',
             color_max=0.4):

    """Casts profiles across the input array and stacks them
    to quantify and correct striped artifacts.

    Keyword arguments:
    in_surface -- 2D array to be fitted
    angle -- angle of the striped arifacts in the in_surface in degree 0=vertical to 90=horizontal
    in_mask -- optional binary mask where 1 marks the area of interest (typically the space of valid image pixels)
               this will take preference over the in_surface.mask
    diagn -- if true plotting is ativated to illustrate the corrections

    """

    # DEBUG in_surface_rot = result
    # DEBUG orig = extent_orig

    # rotate
    if hasattr(in_surface, 'mask'):
        in_surface_rot, orig, mask_rot = rotateImage(in_surface, angle)
    else:
        in_surface_rot, orig = rotateImage(in_surface, angle)

    # rotate in_mask if available, will overwrite mask_rot derived from in_surface.mask
    if in_mask is not None:
        mask_rot, orig = rotateImage(in_mask, angle)[0:1]

    # apply mask
    if 'mask_rot' in locals():
        in_surface_rot = np.ma.masked_where(mask_rot > 0, in_surface_rot)

    # remove pixels not covered in the original image
    in_surface_rot = np.ma.masked_where(in_surface_rot == -255, in_surface_rot)


    if diagnose:
        fig = plt.figure()
        plt.imshow(in_surface_rot)
        step_size = 100
        for i in np.arange(step_size, in_surface_rot.shape[1], step_size):
            plt.axvline(i)


    # compute median of all profiles
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        col_median = np.ma.median(in_surface_rot, axis=0)

    if diagnose:
        # -- Plot...
        fig, axes = plt.subplots(3)
        axes[0].imshow(in_surface)
        axes[1].imshow(in_surface_rot)
        axes[1].autoscale(False)
        axes[1].plot(orig[:, 0], orig[:, 1], 'ro')
        axes[2].plot(col_median)
        plt.show()

    # generate and rotate stripe correction grid
    destripe_grid = np.zeros((in_surface_rot.shape[0], in_surface_rot.shape[1]))
    for i in range(in_surface_rot.shape[0]):
        destripe_grid[i, :] = col_median

    rows_in, cols_in = destripe_grid.shape
    M = cv2.getRotationMatrix2D((cols_in / 2, rows_in / 2), -angle, 1)
    rows_out, cols_out = in_surface.shape
    M[0, 2] += (rows_out - rows_in) / 2
    M[1, 2] += (cols_out - cols_in) / 2
    destripe_grid = cv2.warpAffine(destripe_grid, M, (cols_out, rows_out), flags=cv2.INTER_NEAREST)

    # destripe
    out_surface = in_surface - destripe_grid

    if diagnose:
        fig, axes = plt.subplots(ncols=3)
        axes[0].imshow(in_surface, cmap=cmap, vmin=-color_max, vmax=color_max)
        axes[0].set_title('Before destriping', fontsize=16)
        axes[0].tick_params(labelsize=14)
        axes[1].imshow(destripe_grid, cmap=cmap, vmin=-color_max, vmax=color_max)
        axes[1].set_title('Destriping grid', fontsize=16)
        axes[1].tick_params(labelsize=14)
        im = axes[2].imshow(out_surface, cmap=cmap, vmin=-color_max, vmax=color_max)
        axes[2].set_title('After destriping', fontsize=16)
        axes[2].tick_params(labelsize=14)
        fig.subplots_adjust(right=0.91)
        cbar_ax = fig.add_axes([0.92, 0.35, 0.02, 0.3])
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.ax.tick_params(labelsize=16)

    # horizontal_correction = True
    if (horizontal_correction):

        # rotate
        if hasattr(in_surface, 'mask'):
            in_surface_rot, orig, mask_rot = rotateImage(in_surface, angle)
        else:
            in_surface_rot, orig = rotateImage(in_surface, angle)

        # rotate in_mask if available, will overwrite mask_rot derived from in_surface.mask
        if in_mask is not None:
            mask_rot, orig = rotateImage(in_mask, angle)[0:1]

        # apply mask
        if 'mask_rot' in locals():
            in_surface_rot = np.ma.masked_where(mask_rot > 0, in_surface_rot)

        # remove pixels not covered in the original image
        in_surface_rot = np.ma.masked_where(in_surface_rot == -255, in_surface_rot)

        if diagnose:
            fig = plt.figure()
            plt.imshow(in_surface_rot)
            step_size = 100
            for i in np.arange(step_size, in_surface_rot.shape[1], step_size):
                plt.axvline(i)

        # compute median of all profiles
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            row_median = np.ma.median(in_surface_rot, axis=1)

        # generate and rotate stripe correction grid
        destripe_grid_horizontal = np.zeros((in_surface_rot.shape[0], in_surface_rot.shape[1]))
        for i in range(in_surface_rot.shape[1]):
            destripe_grid_horizontal[:,i] = row_median

        rows_in, cols_in = destripe_grid_horizontal.shape
        M = cv2.getRotationMatrix2D((cols_in / 2, rows_in / 2), -angle, 1)
        rows_out, cols_out = in_surface.shape
        M[0, 2] += (rows_out - rows_in) / 2
        M[1, 2] += (cols_out - cols_in) / 2
        destripe_grid_horizontal = cv2.warpAffine(destripe_grid_horizontal, M, (cols_out, rows_out), flags=cv2.INTER_NEAREST)

        if diagnose:
            fig, axes = plt.subplots(ncols=3)
            axes[0].imshow(out_surface, cmap=cmap, vmin=-color_max, vmax=color_max)
            axes[0].set_title('Before horizontal destriping', fontsize=16)
            axes[0].tick_params(labelsize=14)

        # destripe horizontal
        out_surface = out_surface - destripe_grid_horizontal

        if diagnose:
            axes[1].imshow(destripe_grid_horizontal, cmap=cmap, vmin=-color_max, vmax=color_max)
            axes[1].set_title('Horizontal Destriping grid', fontsize=16)
            axes[1].tick_params(labelsize=14)
            im = axes[2].imshow(out_surface, cmap=cmap, vmin=-color_max, vmax=color_max)
            axes[2].set_title('After horizontal destriping', fontsize=16)
            axes[2].tick_params(labelsize=14)
            fig.subplots_adjust(right=0.91)
            cbar_ax = fig.add_axes([0.92, 0.35, 0.02, 0.3])
            cb = fig.colorbar(im, cax=cbar_ax)
            cb.ax.tick_params(labelsize=16)

        destripe_grid += destripe_grid_horizontal

    return out_surface, destripe_grid