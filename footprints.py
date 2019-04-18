import numpy as np
import ogr
import osr
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import glob
import os
import misc
import gdal
import math
from xml.dom import minidom

# date when ESA changed the format of the detector footprints
DETFOO_CHANGE_DATE = "20181106"


def get_sensing_time(gml_path):

    """
    Reads the sensing time from string gml:id="" in the detfoo GML
    :param gml_path: full path of the GML file
    :return: string of the format YYYYMMDD
    """

    xml = minidom.parse(gml_path)
    root_elem = xml.getElementsByTagName("eop:Mask")[0]
    root_items = list(root_elem.attributes.items())
    for key, val in root_items:
        if key == 'gml:id':
            gml_id = val
    return gml_id.split('_')[-5][:8]


def get_detector_footprint(gml_path, diagn=False):

    """
    Reads a S-2 L1C detector footprint gml and finds the border lines between neighbouring detectors
    Returns a list of arrays with the footprint outlines. WARNING: This is very ugly code and you may get blind if you
    look at it too long.

    Keyword arguments:
    gml_path -- full path of the GML file
    diagn -- if True a series of plots will be shown, mainly for debugging

    Example:
        gml_path = "S2A_OPER_MSK_DETFOO_SGS__20160125T113126_A003092_T37PEL_B02_MSIL1C.gml"
    """
    # DEBUG
    # gml_path = '/home/stumpf/Data/Aiguelles/S2A_OPER_PRD_MSIL1C_PDMC_20160805T220539_R108_V20160803T103724_20160803T103724.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_SGS__20160803T173235_A005825_T32TLQ_N02.04/QI_DATA/S2A_OPER_MSK_DETFOO_SGS__20160803T173235_A005825_T32TLQ_B03_MSIL1C.gml'
    # gml_path = '/home/stumpf/Data/Aiguelles/S2A_OPER_PRD_MSIL1C_PDMC_20160425T185203_R108_V20160425T103025_20160425T103025.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_SGS__20160425T144616_A004395_T32TLQ_N02.01/QI_DATA/S2A_OPER_MSK_DETFOO_SGS__20160425T144616_A004395_T32TLQ_B03_MSIL1C.gml'
    # gml_path = '/home/stumpf/Data/CapeTown/S2A_OPER_PRD_MSIL1C_PDMC_20151218T181852_R121_V20151218T085036_20151218T085036.SAFE_2/GRANULE/S2A_OPER_MSI_L1C_TL_SGS__20151218T140504_A002549_T34HBJ_N02.01/QI_DATA/S2A_OPER_MSK_DETFOO_SGS__20151218T140504_A002549_T34HBJ_B03_MSIL1C.gml'
    # gml_path = '/home/stumpf/Data/CapeTown/S2A_OPER_PRD_MSIL1C_PDMC_20170530T013237_R121_V20151208T082332_20151208T082332.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_EPA__20170529T222118_A002406_T34HBH_N02.04/QI_DATA/S2A_OPER_MSK_DETFOO_EPA__20170529T222118_A002406_T34HBH_B03_MSIL1C.gml'
    #
    # gml_path = '/home/stumpf/Data/Aiguelles/S2A_MSIL1C_20170420T103021_N0204_R108_T32TLQ_20170420T103454.SAFE/GRANULE/L1C_T32TLQ_A009543_20170420T103454/QI_DATA/MSK_DETFOO_B03.gml'
    # gml_path = '/home/stumpf/Data/Aiguelles/S2A_MSIL1C_20170510T103031_N0205_R108_T32TLQ_20170510T103025.SAFE/GRANULE/L1C_T32TLQ_A009829_20170510T103025/QI_DATA/MSK_DETFOO_B03.gml'
    #
    # gml_path = '/home/stumpf/Data/CoregisSites/36UUU/S2A_MSIL1C_20160617T090020_N0204_R007_T36UUU_20160618T074301.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_SGS__20160617T124623_A005152_T36UUU_N02.04/QI_DATA/S2A_OPER_MSK_DETFOO_SGS__20160617T124623_A005152_T36UUU_B03_MSIL1C.gml'
    # gml_path = '/home/stumpf/Data/CoregisSites/36UUU/S2A_MSIL1C_20160405T085012_N0201_R107_T36UUU_20160407T042952.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_SGS__20160405T124601_A004108_T36UUU_N02.01/QI_DATA/S2A_OPER_MSK_DETFOO_SGS__20160405T124601_A004108_T36UUU_B03_MSIL1C.gml'
    #
    # gml_path = '/home/stumpf/Data/CoregisSites/14SKF/S2A_MSIL1C_20160112T174916_N0201_R055_T14SKF_20160113T013403.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_MTI__20160112T205753_A002912_T14SKF_N02.01/QI_DATA/S2A_OPER_MSK_DETFOO_MTI__20160112T205753_A002912_T14SKF_B04_MSIL1C.gml'
    # gml_path = '/home/stumpf/Data/CoregisSites/14SKF/S2A_MSIL1C_20160511T174344_N0202_R055_T14SKF_20160512T065256.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_MTI__20160511T223327_A004628_T14SKF_N02.02/QI_DATA/S2A_OPER_MSK_DETFOO_MTI__20160511T223327_A004628_T14SKF_B04_MSIL1C.gml'
    #
    # gml_path = '/home/stumpf/Data/CoregisSites/20HNH/S2A_MSIL1C_20160313T141151_N0201_R067_T20HNH_20160313T230452.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_SGS__20160313T211303_A003782_T20HNH_N02.01/QI_DATA/S2A_OPER_MSK_DETFOO_SGS__20160313T211303_A003782_T20HNH_B04_MSIL1C.gml'


    # gml_path = '/home/stumpf/Data/CoregisSites/20HNH/S2A_MSIL1C_20160313T140130_N0201_R067_T20HNH_20160313T222921.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_MTI__20160313T202751_A003782_T20HNH_N02.01/QI_DATA/S2A_OPER_MSK_DETFOO_MTI__20160313T202751_A003782_T20HNH_B04_MSIL1C.gml'
    # gml_path = '/home/stumpf/Data/Sinard/full_series/S2A_MSIL1C_20160717T104026_N0204_R008_T31TGK_20160718T061224.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_SGS__20160717T161142_A005582_T31TGK_N02.04/QI_DATA/S2A_OPER_MSK_DETFOO_SGS__20160717T161142_A005582_T31TGK_B04_MSIL1C.gml'
    # diagn = True


    # test_list = ['/home/stumpf/Data/Aiguelles/S2A_OPER_PRD_MSIL1C_PDMC_20160805T220539_R108_V20160803T103724_20160803T103724.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_SGS__20160803T173235_A005825_T32TLQ_N02.04/QI_DATA/S2A_OPER_MSK_DETFOO_SGS__20160803T173235_A005825_T32TLQ_B03_MSIL1C.gml',
    #             '/home/stumpf/Data/Aiguelles/S2A_OPER_PRD_MSIL1C_PDMC_20160425T185203_R108_V20160425T103025_20160425T103025.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_SGS__20160425T144616_A004395_T32TLQ_N02.01/QI_DATA/S2A_OPER_MSK_DETFOO_SGS__20160425T144616_A004395_T32TLQ_B03_MSIL1C.gml',
    #             '/home/stumpf/Data/CapeTown/S2A_OPER_PRD_MSIL1C_PDMC_20151218T181852_R121_V20151218T085036_20151218T085036.SAFE_2/GRANULE/S2A_OPER_MSI_L1C_TL_SGS__20151218T140504_A002549_T34HBJ_N02.01/QI_DATA/S2A_OPER_MSK_DETFOO_SGS__20151218T140504_A002549_T34HBJ_B03_MSIL1C.gml',
    #             '/home/stumpf/Data/CapeTown/S2A_OPER_PRD_MSIL1C_PDMC_20170530T013237_R121_V20151208T082332_20151208T082332.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_EPA__20170529T222118_A002406_T34HBH_N02.04/QI_DATA/S2A_OPER_MSK_DETFOO_EPA__20170529T222118_A002406_T34HBH_B03_MSIL1C.gml',
    #             '/home/stumpf/Data/Aiguelles/S2A_MSIL1C_20170420T103021_N0204_R108_T32TLQ_20170420T103454.SAFE/GRANULE/L1C_T32TLQ_A009543_20170420T103454/QI_DATA/MSK_DETFOO_B03.gml',
    #             '/home/stumpf/Data/Aiguelles/S2A_MSIL1C_20170510T103031_N0205_R108_T32TLQ_20170510T103025.SAFE/GRANULE/L1C_T32TLQ_A009829_20170510T103025/QI_DATA/MSK_DETFOO_B03.gml',
    #             '/home/stumpf/Data/CoregisSites/36UUU/S2A_MSIL1C_20160617T090020_N0204_R007_T36UUU_20160618T074301.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_SGS__20160617T124623_A005152_T36UUU_N02.04/QI_DATA/S2A_OPER_MSK_DETFOO_SGS__20160617T124623_A005152_T36UUU_B03_MSIL1C.gml',
    #             '/home/stumpf/Data/CoregisSites/36UUU/S2A_MSIL1C_20160405T085012_N0201_R107_T36UUU_20160407T042952.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_SGS__20160405T124601_A004108_T36UUU_N02.01/QI_DATA/S2A_OPER_MSK_DETFOO_SGS__20160405T124601_A004108_T36UUU_B03_MSIL1C.gml',
    #             '/home/stumpf/Data/CoregisSites/14SKF/S2A_MSIL1C_20160112T174916_N0201_R055_T14SKF_20160113T013403.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_MTI__20160112T205753_A002912_T14SKF_N02.01/QI_DATA/S2A_OPER_MSK_DETFOO_MTI__20160112T205753_A002912_T14SKF_B04_MSIL1C.gml',
    #             '/home/stumpf/Data/CoregisSites/14SKF/S2A_MSIL1C_20160511T174344_N0202_R055_T14SKF_20160512T065256.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_MTI__20160511T223327_A004628_T14SKF_N02.02/QI_DATA/S2A_OPER_MSK_DETFOO_MTI__20160511T223327_A004628_T14SKF_B04_MSIL1C.gml',
    #             '/home/stumpf/Data/CoregisSites/20HNH/S2A_MSIL1C_20160313T141151_N0201_R067_T20HNH_20160313T230452.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_SGS__20160313T211303_A003782_T20HNH_N02.01/QI_DATA/S2A_OPER_MSK_DETFOO_SGS__20160313T211303_A003782_T20HNH_B04_MSIL1C.gml',
    #             '/home/stumpf/Data/CoregisSites/20HNH/S2A_MSIL1C_20160313T140130_N0201_R067_T20HNH_20160313T222921.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_MTI__20160313T202751_A003782_T20HNH_N02.01/QI_DATA/S2A_OPER_MSK_DETFOO_MTI__20160313T202751_A003782_T20HNH_B04_MSIL1C.gml',
    #             '/home/stumpf/Data/Sinard/full_series/S2A_MSIL1C_20160717T104026_N0204_R008_T31TGK_20160718T061224.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_SGS__20160717T161142_A005582_T31TGK_N02.04/QI_DATA/S2A_OPER_MSK_DETFOO_SGS__20160717T161142_A005582_T31TGK_B04_MSIL1C.gml']
    #
    # for gml_path in test_list:
    #     get_detector_footprint(gml_path, diagn=True)

    def line(p1, p2):
        A = (p1[1] - p2[1])
        B = (p2[0] - p1[0])
        C = (p1[0] * p2[1] - p2[0] * p1[1])
        return A, B, -C


    def intersection(L1, L2):
        D = L1[0] * L2[1] - L1[1] * L2[0]
        Dx = L1[2] * L2[1] - L1[1] * L2[2]
        Dy = L1[0] * L2[2] - L1[2] * L2[0]
        if D != 0:
            x = Dx / D
            y = Dy / D
            return x, y
        else:
            return False


    def get_border(poly_coords, side='right'):

        """
        Finds and iterates over one of the two sides of the polygon marking the detector coordinates
        and returns the vertices of the segments along the polyline
        """

        if side == 'right':
            # first find UPPER RIGHT points w
            # ind_start = np.argmax(poly_coords, axis=0)[0]
            # find points on the rigth > highest x coordinates
            ind_r = np.where(poly_coords[:, 0] == max(poly_coords[:, 0]))
            # from those find the one with the lowest y coordinate
            ind_start = ind_r[0][np.argmin(poly_coords[ind_r, 1])]

            # print(ind_start)

        if side == 'left':  # UPPER LEFT
            # find points on the left > smallest x coordinates
            ind_l = np.where(poly_coords[:, 0] == min(poly_coords[:, 0]))
            # from those find the one with the highest y coordinate
            ind_start = ind_l[0][np.argmax(poly_coords[ind_l, 1])]

            # print(ind_start)

        # iterate over points as long as the change in x is negative > RIGHT LINE or positive > LEFT LINE

        i = 0
        line = []
        ind = []

        if side == 'right':  # start from UPPER RIGHT and move down the line until horizontal
            delta_y = -1
            while delta_y < 0:
                delta = poly_coords[ind_start + i + 1, :] - poly_coords[ind_start + i]
                delta_y = delta[1]
                line.append(poly_coords[ind_start + i])
                ind.append(ind_start + i)
                i += 1

        if side == 'left':  # start from UPPER LEFT and move up the line until horizontal
            delta_y = 1
            while delta_y > 0:
                delta = poly_coords[ind_start + i + 1, :] - poly_coords[ind_start + i]
                delta_y = delta[1]
                line.append(poly_coords[ind_start + i])
                ind.append(ind_start + i)
                i += 1

        line = np.array(line)
        ind = np.array(ind)

        if side == 'left':  # to have them always starting from high y to low y
            line = np.flipud(line)
            ind = np.array(ind)

        return line, ind


    def get_mid_points(line_iter, line_project, diagnose=False):

        """
        Projects the vertices from the first line on the second line and finds the
        the mid points between the original and the projected points
        """

        # DEBUG
        # line_iter = points_right
        # line_project = points_left

        def project_point(x, y, x1, y1, x2, y2):
            """
            x, y is the source point and x1, y1 to x2, y2 defines the line segment.
            Returns the closest points on the line relative to the source point
            """

            A = x - x1
            B = y - y1
            C = x2 - x1
            D = y2 - y1

            dot = A * C + B * D
            len_sq = C * C + D * D

            if len_sq != 0:
                param = dot / len_sq

                xx = x1 + param * C
                yy = y1 + param * D

            else:
                raise Exception('Undertermined solution when trying to project point on detector border')

            return xx, yy

        def mid_point(p1, p2):
            """
            get the point half way between the two points
            """
            p3 = p2 - (p2 - p1) / 2
            return p3

        def line_intersection(line1, line2):

            xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
            ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

            def det(a, b):
                return a[0] * b[1] - a[1] * b[0]

            div = det(xdiff, ydiff)
            if div == 0:
                raise Exception('lines do not intersect')

            d = (det(*line1), det(*line2))
            x = det(d, xdiff) / div
            y = det(d, ydiff) / div

            return x, y

        line_mid = []

        for i in range(0, len(line_iter)):

            # find the two closest points on the other line
            dist = np.sqrt(np.sum((line_project - line_iter[i]) ** 2, axis=1))
            ind_nn = dist.argsort()[:2]

            # get point coordinates
            x = line_iter[i][0]
            y = line_iter[i][1]
            x1 = line_project[ind_nn][0, 0]
            x2 = line_project[ind_nn][1, 0]
            y1 = line_project[ind_nn][0, 1]
            y2 = line_project[ind_nn][1, 1]

            # project current point on the line
            x_new, y_new = project_point(x, y, x1, y1, x2, y2)

            # get mid point
            p1 = np.array([x, y])
            p2 = np.array([x_new, y_new])
            p_new = mid_point(p1, p2)

            # check if points falls outside the tile and if yes project it back to the tile border
            line_xmin = np.min(np.concatenate((line_iter[:, 0], line_project[:, 0])))
            line_xmax = np.max(np.concatenate((line_iter[:, 0], line_project[:, 0])))
            line_ymin = np.min(np.concatenate((line_iter[:, 1], line_project[:, 1])))
            line_ymax = np.max(np.concatenate((line_iter[:, 1], line_project[:, 1])))

            if not (line_xmin <= p_new[0] <= line_xmax and line_ymin <= p_new[1] <= line_ymax):

                border_flag = True
                # cast a normal form the midpoints
                #  if we define dx=x2-x1 and dy=y2-y1, then the normals are (-dy, dx) and (dy, -dx)
                dx = p_new[0] - x
                dy = p_new[1] - y
                x_new1 = p_new[0] - dy
                y_new1 = p_new[1] + dx

                # and intersect it with the image border
                line1 = ((x, y), (x1, y1))
                line2 = ((p_new[0], p_new[1]), (x_new1, y_new1))
                x_border, y_border = line_intersection(line1, line2)
                p_border = np.array([x_border, y_border])

                line_mid.append(p_border)

            else:
                border_flag = False
                line_mid.append(p_new)

            if diagnose:
                plt_mid = np.array(line_mid)
                plt_min = np.min(np.array(line_mid), axis=0)
                plt_max = np.max(np.array(line_mid), axis=0)
                if i == 0:
                    plt.figure()
                plt.plot(line_iter[:, 0], line_iter[:, 1], 'ro-')
                plt.plot(line_project[:, 0], line_project[:, 1], 'co-')
                plt.gca().set_aspect(1.0)
                plt.axis([plt_min[0] - 5000, plt_max[0] + 5000, plt_min[1] - 5000, plt_max[1] + 5000])
                plt.plot([x, x_new], [y, y_new], 'bo-')
                if border_flag:
                    plt.plot([p_new[0], x_new1], [p_new[1], y_new1], 'mo-')
                    plt.plot([p_border[0]], [p_border[1]], 'mo')
                plt.plot(plt_mid[:, 0], plt_mid[:, 1], 'go-')

        line_mid = np.array(line_mid)

        return line_mid

    # handles new detfoo metadata format since 06/11/2018
    sensing_date = get_sensing_time(gml_path)
    if sensing_date >= DETFOO_CHANGE_DATE:
        use_new_df = True
    else:
        use_new_df = False

    # get detector positions from GML metadata files
    inSource = ogr.Open(gml_path)
    inLayer = inSource.GetLayer()
    scene_extent = inLayer.GetExtent()

    n_detectors = inLayer.GetFeatureCount()
    detector_ids = list()
    detector_coord = []

    if diagn:
        jet = plt.get_cmap('jet')
        colors = iter(jet(np.linspace(0, 1, 10)))
        fig1 = plt.figure()

    for i in range(0, n_detectors):
        feat = inLayer.GetFeature(i)
        detector_ids.append(list(feat.items().values()))
        geom = feat.GetGeometryRef()
        geom.GetGeometryName()
        poly = geom.GetGeometryRef(0)
        coords = np.delete(np.asarray(poly.GetPoints()), -1, 1)
        detector_coord.append(coords)

    # Plot detector layout
    if diagn:
        ax = plt.gca()
        ax.set_title('Detector footprints', fontsize=16)
        ax.set_xlim(scene_extent[0], scene_extent[1])
        ax.set_ylim(scene_extent[2], scene_extent[3])

        for i in range(0, n_detectors):
            # create paths for plotting
            codes = []
            codes += [Path.MOVETO] + (len(detector_coord[i]) - 2) * [Path.LINETO] + [Path.CLOSEPOLY]
            path = Path(detector_coord[i], codes)

            patch = patches.PathPatch(path, facecolor='0.8', edgecolor='black', alpha=0.1)
            ax.add_patch(patch)

            plt.gca().add_patch(patch)
            plt.gca().add_artist(patch)
            ax.set_aspect(1.0)

        plt.show()

    if use_new_df:  # we can safe all the hustle
        return detector_coord

    new_detector_coord = []
    final_iter = len(detector_coord)
    start_iter = 0

    # check if there if there is only one detector element
    if len(detector_coord) == 1:
        new_detector_coord.append(detector_coord[0])
        if diagn:
            plt.plot(new_detector_coord[0][:, 0], new_detector_coord[0][:, 1], 'o-', color=next(colors))
        return new_detector_coord

    # test if the first polygon is so small that the second still starts at the upper left
    points_left, ind_left = get_border(detector_coord[start_iter + 1], side='left')
    if len(points_left) == 1:
        # happens if the first detector element is so small that the second still starts in the upper left
        start_iter += 1
        # David fix: test again after having skipped the first detector
        # check if there is only one detector element
        if (final_iter - start_iter) == 1:
            new_detector_coord.append(detector_coord[start_iter])
            if diagn:
                plt.plot(new_detector_coord[start_iter][:, 0], new_detector_coord[start_iter][:, 1], 'o-', color=next(colors))
            return new_detector_coord

    for i in range(start_iter, final_iter):

        if i == final_iter - 1:  # in the final iteration only the left side has to be fixed

            mask = np.ones(len(detector_coord[i]), np.bool)
            mask[ind_left_previous] = 0
            # the last coordinate is always a duplicate
            mask[-1] = 0
            # left side of the polygon, y coordinates must be ascending
            insert_left_points = mid_points_previous[mid_points_previous[:, 1].argsort()]

            # put them all together
            y_ll = scene_extent[3] - 10900 * 10
            if y_ll < detector_coord[i][:, 1].min() and detector_coord[i][mask, 0].min() == detector_coord[i][mask, 0].max(): # case that lower left comprises a triangle of no data
                # i.e. if the min detector extent does not reach the theoretical scene extent and the left border is straight
                new_detector_coord.append(np.vstack((mid_points_previous, detector_coord[i][mask], mid_points_previous[0])))
            else:
                new_detector_coord.append(np.vstack((insert_left_points, detector_coord[i][mask], insert_left_points[0])))

        else:

            # get borders on both sides of the overlap
            points_right, ind_right = get_border(detector_coord[i], side='right')
            points_left, ind_left = get_border(detector_coord[i + 1], side='left')

            # get center points for each vertices along the borders
            mid_points_right = get_mid_points(points_right, points_left, diagnose=False)
            if np.all(points_right == [scene_extent[1], scene_extent[2]]): # when arriving at the right side and the right side of the overlap is out of bounds
                mid_points = mid_points_right
            else:
                mid_points_left = get_mid_points(points_left, points_right, diagnose=False)

                # stack and sort them in descending order
                mid_points = np.vstack((mid_points_right, mid_points_left))
                mid_points = mid_points[mid_points[:, 1].argsort()[::-1]]

            if i == start_iter: # first detector element

                # invert index and use it to get the polygon coordinates which are not along the detector border
                # right side of the current polygon
                mask = np.ones(len(detector_coord[i]), np.bool)
                mask[ind_right] = 0
                # the last coordinate is always a duplicate
                mask[-1] = 0

                # check if the detector element covers the upper left corner
                contains_upper_left_corner = any(
                    [([scene_extent[0], scene_extent[3]] == pair).all() for pair in detector_coord[i]])
                # check if the detector element covers the lower scene extent
                contains_lower_scene_extent = scene_extent[3] - 10900 * 10 == detector_coord[i][:, 1].min()

                # check if the intersection area hits the scene corner
                intersection_at_scene_corner = (mid_points[-1][0] - scene_extent[0] > 0.0001 and mid_points[-1][1] - scene_extent[2] > 0.0001)

                if intersection_at_scene_corner and not contains_lower_scene_extent:

                        # # intersect the extension of the mid line with the outlines of the detectors
                        # L1 = line(mid_points[-2], mid_points[-1])
                        # L2 = line(detector_coord[i][-2], detector_coord[i][-1])
                        # R1 = intersection(L1, L2)
                        #
                        # L1 = line(mid_points[-2], mid_points[-1])
                        # L2 = line(detector_coord[i][-3], detector_coord[i][-2])
                        # R2 = intersection(L1, L2)
                        #
                        # # add the intersection point with the shortest distance
                        # dist1 = np.linalg.norm(mid_points[-1] - R1)
                        # dist2 = np.linalg.norm(mid_points[-1] - R2)
                        #
                        # if dist1 < dist2:
                        #     mid_points = np.vstack((mid_points, R1))
                        # else:
                        #     mid_points = np.vstack((mid_points, R2))
                        #
                        # new_detector_coord.append(np.vstack((mid_points,
                        #                                      detector_coord[i][-2],
                        #                                      mid_points[0])))

                        # intersect the extension of the mid line with the outlines of the detectors
                        L1 = line(mid_points[-2], mid_points[-1])
                        L2 = line(detector_coord[i][-2], detector_coord[i][-1])
                        R1 = intersection(L1, L2)

                        L1 = line(mid_points[-2], mid_points[-1])
                        L2 = line(detector_coord[i][-3], detector_coord[i][-2])
                        R2 = intersection(L1, L2)

                        L1 = line(mid_points[-2], mid_points[-1])
                        L2 = line(detector_coord[i][-4], detector_coord[i][-3])
                        R3 = intersection(L1, L2)

                        # add the intersection point with the shortest distance
                        dist1 = np.linalg.norm(mid_points[-1] - R1)
                        dist2 = np.linalg.norm(mid_points[-1] - R2)
                        dist3 = np.linalg.norm(mid_points[-1] - R3)

                        # find closest intersection
                        index_min = np.argmin([dist1, dist2, dist3])

                        if index_min == 0:  # horizontal tile border
                            mid_points = np.vstack((mid_points, R1))
                            mid_points = np.vstack((mid_points, [scene_extent[0], scene_extent[2]]))
                            new_detector_coord.append(np.vstack((mid_points,
                                                                 detector_coord[i][-2],
                                                                 mid_points[0])))
                        elif index_min == 1:  # vertical tile border
                            mid_points = np.vstack((mid_points, R2))
                            new_detector_coord.append(np.vstack((mid_points,
                                                                 detector_coord[i][-2],
                                                                 mid_points[0])))
                        elif index_min == 2:  # the detector element itself
                            new_detector_coord.append(np.vstack((mid_points,
                                                                 detector_coord[i][-3],
                                                                 detector_coord[i][-2],
                                                                 mid_points[0])))

                elif contains_upper_left_corner and contains_lower_scene_extent: # case that image starts at the upper left

                    # merge and add the duplicate back to close the polygon
                    new_detector_coord.append(np.vstack((mid_points, detector_coord[i][mask][::-1], mid_points[0])))

                elif contains_upper_left_corner and not contains_lower_scene_extent:
                    # case that upper left comprise no data AND lower left comprises no data
                    # i.e. if the lower detector extent is not reach the theoretical scene extent

                    # if there is a duplicate point in the two arrays storing the left original points of the polygon
                    if any((detector_coord[i][:min(ind_right)] == x).all() for x in detector_coord[i][mask]):
                        new_detector_coord.append(np.vstack((mid_points,
                                                             detector_coord[i][mask][-1],
                                                             detector_coord[i][:min(ind_right)],
                                                             mid_points[0])))
                    else:

                        new_detector_coord.append(np.vstack((mid_points,
                                                             detector_coord[i][mask],
                                                             detector_coord[i][:min(ind_right)],
                                                             mid_points[0])))

                else:
                    # case that only the upper left comprises no data
                    # TODO maybe use this notation in the other cases for clarity
                    new_detector_coord.append(np.vstack((mid_points,
                                                         detector_coord[i][max(ind_right)+1:],
                                                         detector_coord[i][:min(ind_right)],
                                                         mid_points[0])))

                # Prepare left side of the next polygon
                mid_points_previous = mid_points
                ind_left_previous = ind_left

                if diagn:
                    plt.plot(new_detector_coord[i-start_iter][:, 0], new_detector_coord[i-start_iter][:, 1],
                             'o-', color=next(colors))
                    plt.gca().set_aspect(1.0)
                    plt.axis([scene_extent[0] - 1000, scene_extent[1] + 1000, scene_extent[2] - 1000,
                              scene_extent[3] + 1000])

            else: # not the first detector element

                # fix both sides of the current polygon
                mask = np.ones(len(detector_coord[i]), np.bool)
                ind_all = np.concatenate((ind_left_previous, ind_right))
                mask[ind_all] = 0
                # the last coordinate is always a duplicate
                mask[-1] = 0
                # left side of the polygon, y coordinates must be in ascending order
                insert_left_points = mid_points_previous[mid_points_previous[:, 1].argsort()]

                # put them all together

                # if the intersection area hits the scene corner mid_points do not comprise a proper point on the border
                if (mid_points[-1][0] - scene_extent[0] > 0.0001 and mid_points[-1][1] - scene_extent[2] > 0.0001):

                    # This condition is also true if the lower part comprises no data
                    if scene_extent[3] - 10900 * 10 < detector_coord[i][:, 1].min():
                        # case that lower left comprises no data
                        # i.e. if the lower detector extent does not reach the theoretical scene extent
                        new_detector_coord.append(np.vstack((insert_left_points, mid_points, detector_coord[i][mask], insert_left_points[0])))

                    # ...or well the intersection hits the border
                    elif len(mid_points) > 1: # if there are several mid_points determine the intersection from a line through the mid points

                        # intersect the extension of the mid line with the outlines of the detectors
                        L1 = line(mid_points[-2], mid_points[-1])
                        L2 = line(detector_coord[i][-2], detector_coord[i][-1])
                        R1 = intersection(L1, L2)

                        L1 = line(mid_points[-2], mid_points[-1])
                        L2 = line(detector_coord[i][-3], detector_coord[i][-2])
                        R2 = intersection(L1, L2)

                        # add the intersection point with the shortest distance
                        dist1 = np.linalg.norm(mid_points[-1] - R1)
                        dist2 = np.linalg.norm(mid_points[-1] - R2)

                        if dist1 < dist2:
                            mid_points = np.vstack((mid_points, R1))
                            new_detector_coord.append(np.vstack((insert_left_points, mid_points, insert_left_points[0])))
                            mid_points = np.vstack((mid_points, [scene_extent[0], scene_extent[2]] ))
                        else:
                            mid_points = np.vstack((mid_points, R2))
                            new_detector_coord.append(np.vstack((insert_left_points, mid_points, detector_coord[i][-2], insert_left_points[0])))

                    else: # if there is only one mid_point determine the intersection from the detector limited shifted to the mid_point
                        # plt.plot(points_left[:,0], points_left[:,1], 'or')
                        # plt.plot(mid_points[0][0], mid_points[0][1], 'or')
                        L_left = line(points_left[-2], points_left[-1])
                        C_parallel = L_left[0]*mid_points[0][0] + L_left[1]*mid_points[0][1]
                        L_parallel = (L_left[0], L_left[1], C_parallel)

                        # plt.plot(scene_extent[0], scene_extent[2], 'og')
                        # plt.plot(scene_extent[1], scene_extent[2], 'og')
                        lower_limit = line([scene_extent[0], scene_extent[2]], [scene_extent[1], scene_extent[2]])
                        lower_point = intersection(L_parallel, lower_limit)
                        # plt.plot(lower_point[0], lower_point[1], 'oy')

                        # plt.plot(scene_extent[1], scene_extent[2], 'ob')
                        # plt.plot(scene_extent[1], scene_extent[3], 'ob')
                        right_limit = line([scene_extent[1], scene_extent[2]], [scene_extent[1], scene_extent[3]])
                        right_point = intersection(L_parallel, right_limit)
                        # plt.plot(right_point[0], right_point[1], 'oy')

                        mid_points = np.vstack((right_point, mid_points, lower_point))

                        new_detector_coord.append(np.vstack((insert_left_points, mid_points, insert_left_points[0])))

                elif not detector_coord[i][mask].size or detector_coord[i][mask][0][0] == scene_extent[0]:  # if on the left side or original points
                    # new_detector_coord.append(np.vstack((insert_left_points, mid_points, insert_left_points[0])))
                    new_detector_coord.append(np.vstack((insert_left_points, mid_points, detector_coord[i][mask], insert_left_points[0])))

                # when arriving at the right side of the scene change the order of concatenation
                # elif detector_coord[i][mask][0][0] == scene_extent[1]:
                #     new_detector_coord.append(
                #         np.vstack((insert_left_points, detector_coord[i][mask], mid_points, insert_left_points[0])))
                # when arriving at the right side of the scene change the order of concatenation
                elif scene_extent[1] in detector_coord[i][:, 0]:
                    # check if the element's left side comprises the theoretical lower scene extent
                    # i.e. there is no no data area at the lower left of this element
                    lower_left_on_scene_extent = any(
                        [y_coordinate <= scene_extent[3] - 10900 * 10 for y_coordinate in insert_left_points[:, 1]])
                    if lower_left_on_scene_extent:
                        new_detector_coord.append(
                            np.vstack((insert_left_points, detector_coord[i][mask], mid_points, insert_left_points[0])))
                    else:
                        new_detector_coord.append(np.vstack((insert_left_points,
                                                             detector_coord[i][mask][1],
                                                             mid_points,
                                                             detector_coord[i][mask][0],
                                                             insert_left_points[0])))

                elif detector_coord[i][mask].size == 2: # Not yet on the right side but nodata area at the upper right
                    new_detector_coord.append(np.vstack((insert_left_points, detector_coord[i][mask], mid_points, insert_left_points[0])))

                # David: added this correction, but I'm not absolutely sure about correctness ... experience will tell us <<<
                elif detector_coord[i][mask].size == 4 and detector_coord[i][mask][1][0] == scene_extent[0] and \
                                detector_coord[i][mask][1][1] == scene_extent[2]:
                    # Not yet on the right side but nodata area at the upper right and lower left corner in detector
                    new_detector_coord.append(np.vstack((insert_left_points, detector_coord[i][mask][0],
                                                         mid_points, detector_coord[i][mask][1],
                                                         insert_left_points[0])))
                    # David:  end >>>

                # the element is on the left extent of the scene but does not contain the lower left corner
                # i.e. no data area at the lower left
                elif not any(detector_coord[i][detector_coord[i][:, 0] == scene_extent[0], 1] == scene_extent[2]):
                    new_detector_coord.append(np.vstack((insert_left_points, mid_points, detector_coord[i][mask], insert_left_points[0])))

                else:
                    raise AssertionError("An unkown error occured when trying to generate the detector footprints")

                # Prepare left side of the next polygon
                mid_points_previous = mid_points
                ind_left_previous = ind_left

        if diagn:
            plt.plot(new_detector_coord[i-start_iter][:, 0], new_detector_coord[i-start_iter][:, 1],
                     'o-', color=next(colors))

    return new_detector_coord


def plot_ogr_polygon(ogrGeometry, plot_title=""):
    scene_extent = ogrGeometry.GetEnvelope()
    plt.figure()
    ax = plt.gca()
    ax.set_xlim(scene_extent[0], scene_extent[1])
    ax.set_ylim(scene_extent[2], scene_extent[3])
    plt.title(plot_title)

    if ogrGeometry.GetGeometryType() == ogr.wkbPolygon:
        paths = []
        codes = []
        all_x = []
        all_y = []
        for i in range(ogrGeometry.GetGeometryCount()):
            # Read ring geometry and create path
            r = ogrGeometry.GetGeometryRef(i)
            x = [r.GetX(j) for j in range(r.GetPointCount())]
            y = [r.GetY(j) for j in range(r.GetPointCount())]
            # skip boundary between individual rings
            codes += [mpath.Path.MOVETO] + (len(x) - 1) * [mpath.Path.LINETO]
            all_x += x
            all_y += y
        path = mpath.Path(np.column_stack((all_x, all_y)), codes)
        paths.append(path)

        # Add paths as patches to axes
        for path in paths:
            patch = mpatches.PathPatch(path, facecolor='0.8', edgecolor='black')
            ax.add_patch(patch)
            plt.gca().add_patch(patch)
            plt.gca().add_artist(patch)
            ax.set_aspect(1.0)

        plt.show()

    else:
        raise TypeError('The input is not of the type osgeo.ogr.Geometry Polygon')


def intersect_footprint_outlines(detector_footprint1, detector_footprint2, diagnose=False):

    # get union of the first footprints
    union1 = None

    for i in range(len(detector_footprint1) - 1):

        if union1 is None:

            ' '.join([str(j) for j in detector_footprint1[i]])

            wkt1 = "POLYGON((" + ', '.join([str(j).strip("[]") for j in detector_footprint1[i]]) + "))"
            wkt2 = "POLYGON((" + ', '.join([str(j).strip("[]") for j in detector_footprint1[i + 1]]) + "))"
            poly1 = ogr.CreateGeometryFromWkt(wkt1)
            poly2 = ogr.CreateGeometryFromWkt(wkt2)

            union1 = poly1.Union(poly2)

        else:
            wkt2 = "POLYGON((" + ', '.join([str(j).strip("[]") for j in detector_footprint1[i + 1]]) + "))"
            poly2 = ogr.CreateGeometryFromWkt(wkt2)
            union1 = union1.Union(poly2)

    # get union of the second footprints
    union2 = None

    for i in range(len(detector_footprint2) - 1):

        if union2 is None:

            ' '.join([str(j) for j in detector_footprint2[i]])

            wkt1 = "POLYGON((" + ', '.join([str(j).strip("[]") for j in detector_footprint2[i]]) + "))"
            wkt2 = "POLYGON((" + ', '.join([str(j).strip("[]") for j in detector_footprint2[i + 1]]) + "))"
            poly1 = ogr.CreateGeometryFromWkt(wkt1)
            poly2 = ogr.CreateGeometryFromWkt(wkt2)

            union2 = poly1.Union(poly2)

        else:
            wkt2 = "POLYGON((" + ', '.join([str(j).strip("[]") for j in detector_footprint2[i + 1]]) + "))"
            poly2 = ogr.CreateGeometryFromWkt(wkt2)
            union2 = union2.Union(poly2)

    # intersect
    footprint_intersection = union1.Intersection(union2)

    if diagnose:
        plot_ogr_polygon(footprint_intersection, plot_title ='Intersection of tiles')

    return footprint_intersection


def get_S2_detector_intersections(master_slave, current_work_folder, diagnose=False):

    """
    Generate a SHP file of the intersection of two detector footprints, performs an internal check to avoid
    multipolygons in the output shapefile.
    :param master_slave: master slave dictionary
    :param current_work_folder: folder to which the outputs will be written
    :param diagnose: if True diangostic plots are provided
    :return: path to SHP file on disk
    """

    if not master_slave['master']['satellite'] == master_slave['slave']['satellite'] == 'S2':
        raise Exception('Both granules should be S2. Found ' + master_slave['master'] + ' and ' + master_slave['slave'])

    # query image bands to for geospatial referencing
    band_path1 = master_slave['master']['match_band']
    if not band_path1:
        raise Exception("No JP2 files found in the specified folder")
    band_path2 = master_slave['slave']['match_band']
    if not band_path2:
        raise Exception("No JP2 files found in the specified folder")

    # get footprints by figuring out the middle between neighboring detector elements
    detector_footprint1 = get_detector_footprint(master_slave['master']['azimuth'], diagn=diagnose)
    detector_footprint2 = get_detector_footprint(master_slave['slave']['azimuth'], diagn=diagnose)

    # write detector footprints
    out_shape_1 = os.path.join(current_work_folder, 'detector_footprint1.shp')
    out_shape_2 = os.path.join(current_work_folder, 'detector_footprint2.shp')
    polylist2shp(detector_footprint1, out_shape_1, band_path1)
    polylist2shp(detector_footprint2, out_shape_2, band_path2)

    # read back the detector footprints
    # out_shape_1 = '/home/stumpf/Data/CoregisSites/Kaikoura/correlation/S2A_OPER_MSI_L1C_TL_SGS__20160101T001139_A002743_T59GQP_B04_S2A_OPER_MSI_L1C_TL_MPS__20161206T171745_A007605_T59GQP_B04/detector_footprint1.shp'
    # out_shape_2 = '/home/stumpf/Data/CoregisSites/Kaikoura/correlation/S2A_OPER_MSI_L1C_TL_SGS__20160101T001139_A002743_T59GQP_B04_S2A_OPER_MSI_L1C_TL_MPS__20161206T171745_A007605_T59GQP_B04/detector_footprint2.shp'

    ds1 = ogr.Open(out_shape_1)
    layer1 = ds1.GetLayer()
    ds2 = ogr.Open(out_shape_2)
    layer2 = ds2.GetLayer()

    # intersect and write to disk
    out_shape3 = os.path.join(current_work_folder, 'detector_intersections.shp')
    # check if output file already exists
    if os.path.isfile(out_shape3):
        os.remove(out_shape3)
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds3 = driver.CreateDataSource(out_shape3)
    srs = layer1.GetSpatialRef()
    layer3 = ds3.CreateLayer('layer3', srs, ogr.wkbPolygon)
    layer1.Intersection(layer2, layer3)

    # get all polygons
    polygons = []
    for feat in layer3:
        geom = feat.geometry()
        if geom.GetGeometryName() == 'MULTIPOLYGON':
            for polygon in geom:
                # print(polygon.GetGeometryName())
                for i in range(polygon.GetGeometryCount()):
                    r = polygon.GetGeometryRef(i)
                    x = [r.GetX(j) for j in range(r.GetPointCount())]
                    y = [r.GetY(j) for j in range(r.GetPointCount())]
                    polygons.append((x, y))
        else:
            for i in range(geom.GetGeometryCount()):
                r = geom.GetGeometryRef(i)
                x = [r.GetX(j) for j in range(r.GetPointCount())]
                y = [r.GetY(j) for j in range(r.GetPointCount())]
                polygons.append((x, y))

    # if the intersection shape contains multipolygons
    if not layer3.GetFeatureCount() == len(polygons):

        # write to disk
        out_shape4 = os.path.join(current_work_folder, 'detector_intersections_clean.shp')

        # check if output file already exists
        if os.path.isfile(out_shape4):
            os.remove(out_shape4)

        # create output shapefile
        driver = ogr.GetDriverByName('ESRI Shapefile')
        ds4 = driver.CreateDataSource(out_shape4)
        srs = layer1.GetSpatialRef()
        layer4 = ds4.CreateLayer('layer', srs, ogr.wkbPolygon)
        layer_defn = layer4.GetLayerDefn()  # gets parameters of the current shapefile

        # write polygons
        for i, polygon in enumerate(polygons):
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for x, y in zip(polygon[0], polygon[1]):
                ring.AddPoint(x, y)
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)

            feature = ogr.Feature(layer_defn)
            feature.SetGeometry(poly)
            feature.SetFID(i)
            layer4.CreateFeature(feature)
            feature = None

        outshape = out_shape4
        print('Resulting shapefile contains ' + str(layer4.GetFeatureCount()) + ' features.')
        if diagnose:
            misc.plot_layer(layer4)

    else:
        outshape = out_shape3
        print('Resulting shapefile contains ' + str(layer3.GetFeatureCount()) + ' features.')
        if diagnose:
            misc.plot_layer(layer3)

    ds1 = None
    ds2 = None
    ds3 = None
    ds4 = None

    return outshape


def polylist2shp(polylist, outfile, ref, overwrite=True):
    """Takes a list of numpy arrays which hold the ordered vertices coordinates of polygons and write them
    to a shapefile.

    Keyword arguments:
    polylist -- list of 2D numpy arrays with two columns holding the x and y coordinates
                Example for two polygons. Note that the start point must be duplicated!
                detector_footprint1[0]
                array([[  604894.83558543,   990239.99999997],
                       [  604931.8479541 ,   990407.32603429],
                       [  607384.45994306,  1001462.00991814],
                       [  607416.88491911,  1001608.47026906],
                       [  609019.91180068,  1008827.94445133],
                       [  609779.99999999,  1012244.43344573],
                       [  609780.        ,   990240.        ],
                       [  604894.83558543,   990239.99999997]])
                detector_footprint1[1]
                array([[  507151.11350485,  1100040.        ],
                       [  507114.46128252,  1099871.55147799],
                       [  505522.20413839,  1092599.60161843],
                       [  505494.83529568,  1092474.59333894],
                       [  501968.06849405,  1076162.25267087],
                       [  501936.88815993,  1076020.7818152 ],
                       [  500752.08407593,  1070645.14637844],
                       [  499980.        ,  1067136.11626875],
                       [  499980.        ,  1100040.        ],
                       [  507151.11350485,  1100040.        ]])

    outfile -- full path of the output shape file

    ref -- path two a geospatial raster from which the CRS code will be extracted

    overwrite -- set to False to avoid that existing files will be overwritten

    """

    ########################################## side note on problems when reading GML ##################################
    # inSource = ogr.Open(gml_path1)
    # inLayer = inSource.GetLayer()
    # sr = inLayer.GetSpatialRef()              # it seems the GML reader does not pass on the EPSG code

    # inSource = ogr.Open(gml_path1)            # also this does not work
    # inLayer = inSource.GetLayer()
    # f = inLayer.GetNextFeature()
    # g = f.GetGeometryRef()
    # r = g.GetSpatialReference()
    # str(r.GetAuthorityCode("PROJCS"))

    # while the conversion to SHP with ogr2ogr seems to pass on the CRS and it is read properly from the shape file
    # inSource = ogr.Open("S2A_OPER_MSK_DETFOO_SGS__20160125T113126_A003092_T37PEL_B02_MSIL1C/MaskFeature.shp")
    # inLayer = inSource.GetLayer()
    # sr = inLayer.GetSpatialRef()
    # sr.ExportToProj4()
    ####################################################################################################################

    # test
    # polylist = detector_footprint1

    # get CRS from reference raster
    ref_rast = gdal.Open(ref)
    spatialReference = osr.SpatialReference()
    spatialReference.ImportFromWkt(ref_rast.GetProjectionRef())
    # spatialReference.ExportToProj4()

    # set shapefile driver
    driver = ogr.GetDriverByName('ESRI Shapefile')

    # check if target folder exists if not create
    if not os.path.isdir(os.path.dirname(outfile)):
        os.makedirs(os.path.dirname(outfile))

    # check if output file already exists
    if os.path.isfile(outfile):
        if overwrite:
            os.remove(outfile)
        else:
            raise IOError('File already exists and overwrite was set to False')

    # create shapefile
    shapeData = driver.CreateDataSource(outfile)

    # create layer within the shape file
    options = []
    options.append('OVERWRITE=YES')
    layer = shapeData.CreateLayer('detectorFootprint', spatialReference,
                                  ogr.wkbPolygon)  # this will create a corresponding layer for our data with given spatial information.
    layer_defn = layer.GetLayerDefn()  # gets parameters of the current shapefile

    # write the polygons into the layer
    for i in range(0, len(polylist)):

        ring = ogr.Geometry(ogr.wkbLinearRing)

        for j in range(0, len(polylist[i])):
            ring.AddPoint(polylist[i][j,][0], polylist[i][j,][1])
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)

        feature = ogr.Feature(layer_defn)
        feature.SetGeometry(poly)
        feature.SetFID(i)
        layer.CreateFeature(feature)

    # free the shapefile
    shapeData.Destroy()

    return


def getAngle(detector_footprints, tol=0.5):

    """
    Takes a list of arrays with x,y point coordinates and computes the mean angle
    of the lines spanning between the points. Horizontal and vertical angles are excluded.

    Keyword arguments:
    detector_footprints -- list of arrays or list of list of arrays with x,y coordinates
    tol -- +- deviation for not considering horizontal and vertical angles in degrees
    
    DEBUG:
    detector_footprints = [detector_footprint1, detector_footprint2]
    
    """
    # TODO test one of the GMLs used for testing get_detector_footprints testing

    angles = []

    if type(detector_footprints[0]) is list:

        for i in range(0, len(detector_footprints)):
            current_footprint = detector_footprints[i]

            for j in range(0, len(current_footprint)):

                current_detector = current_footprint[j]

                for k in range(0, len(current_detector)-1):
                    dx, dy = current_detector[k+1] - current_detector[k]
                    rads = math.atan2(abs(dy), abs(dx))
                    degs = 90-math.degrees(rads)
                    angles = np.append(angles, degs)
    else:

        for j in range(0, len(detector_footprints)):

            current_detector = detector_footprints[j]

            for k in range(0, len(current_detector) - 1):
                dx, dy = current_detector[k + 1] - current_detector[k]
                rads = math.atan2(abs(dy), abs(dx))
                degs = 90 - math.degrees(rads)
                angles = np.append(angles, degs)


    ind_horizontal = np.where(np.logical_and((90 - tol) <= angles, angles <= (90 + tol)))
    ind_vertical = np.where(np.logical_and((0 - tol) <= angles, angles <= (0 + tol)))
    ind_all = np.append(ind_horizontal, ind_vertical)
    angles = np.delete(angles, ind_all)

    mean_angle = np.mean(angles)
    angle_variance = np.std(angles)

    return  mean_angle, angle_variance


def reproject_shapefile(shapefile_path_in, target_srs_wkt, shapefile_path_out):


    # DEBUG
    # shapefile_path_in = out_shape
    # target_srs_wkt = proj_band
    # shapefile_path_out = shapefile_path_out

    driver = ogr.GetDriverByName('ESRI Shapefile')

    # input SpatialReference
    inDataSet = driver.Open(shapefile_path_in)
    inLayer = inDataSet.GetLayer()
    inSpatialRef = inLayer.GetSpatialRef()

    # output SpatialReference
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromWkt(target_srs_wkt)

    # create the CoordinateTransformation
    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

    # get the input layer
    inDataSet = driver.Open(shapefile_path_in)
    inLayer = inDataSet.GetLayer()

    # create the output layer
    outputShapefile = shapefile_path_out
    if os.path.exists(outputShapefile):
        driver.DeleteDataSource(outputShapefile)
    outDataSet = driver.CreateDataSource(outputShapefile)
    layer_name = inLayer.GetName()
    geom_type = inLayer.GetGeomType()
    outLayer = outDataSet.CreateLayer(layer_name, outSpatialRef, geom_type=geom_type)

    # add fields
    inLayerDefn = inLayer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)

    # get the output layer's feature definition
    outLayerDefn = outLayer.GetLayerDefn()

    # loop through the input features
    inFeature = inLayer.GetNextFeature()
    while inFeature:
        # get the input geometry
        geom = inFeature.GetGeometryRef()
        # reproject the geometry
        geom.Transform(coordTrans)
        # create a new feature
        outFeature = ogr.Feature(outLayerDefn)
        # set the geometry and attribute
        outFeature.SetGeometry(geom)
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
        # add the feature to the shapefile
        outLayer.CreateFeature(outFeature)
        # dereference the features and get the next input feature
        outFeature = None
        inFeature = inLayer.GetNextFeature()

    # Save and close the shapefiles
    inDataSet = None
    outDataSet = None


def align_detector_footprints(s2_bands, aoi_path):

    # get reference info
    ds_ref = ogr.Open(aoi_path)
    layer_ref = ds_ref.GetLayer()
    srs_ref = layer_ref.GetSpatialRef()
    proj_ref = srs_ref.ExportToWkt()

    # DEBUG
    # s2_bands = bands3_tiff_ord
    # aoi_path= aoi_path

    gml_list = []
    for band3 in s2_bands:

        # DEBUG
        # print(band3)
        # band3 = '/home/stumpf/Data/Aiguelles/S2A_OPER_PRD_MSIL1C_PDMC_20161101T190224_R108_V20161101T103202_20161101T103202.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_SGS__20161101T155732_A007112_T31TGK_N02.04/IMG_DATA/S2A_OPER_MSI_L1C_TL_SGS__20161101T155732_A007112_T31TGK_B03.tif'
        # band3 = '/home/stumpf/Data/Aiguelles/S2A_OPER_PRD_MSIL1C_PDMC_20160805T220539_R108_V20160803T103724_20160803T103724.SAFE/GRANULE/S2A_OPER_MSI_L1C_TL_SGS__20160803T173235_A005825_T32TLQ_N02.04/IMG_DATA/S2A_OPER_MSI_L1C_TL_SGS__20160803T173235_A005825_T32TLQ_B03.tif'

        # get image GML
        qi_path = os.path.join(os.path.dirname(os.path.dirname(band3)), 'QI_DATA')
        gml_path = glob.glob(os.path.join(qi_path, '*DETFOO*B03*.gml'))[0]

        # get detector footprint and store as shapefile
        detector_footprint = get_detector_footprint(gml_path, diagn=False)
        out_shape = os.path.join(os.path.splitext(gml_path)[0] + '.shp')
        polylist2shp(detector_footprint, out_shape, band3)

        # get infos on detector footprints
        ds_detfoo = ogr.Open(out_shape)
        layer_detfoo = ds_detfoo.GetLayer()
        srs_layer_detfoo = layer_detfoo.GetSpatialRef()
        proj_detfoo = srs_layer_detfoo.ExportToWkt()
        ds_detfoo = None

        # if not in the reference projection do reproject
        if not proj_detfoo == proj_ref:
            print('Projections do not match. Reprojecting and clipping to AOI...')
            shapefile_reprojected = os.path.splitext(out_shape)[0] + '_aligned.shp'
            reproject_shapefile(out_shape, proj_ref, shapefile_reprojected)

            # read back the detector footprints
            ds_detfoo = ogr.Open(shapefile_reprojected)
            layer_detfoo = ds_detfoo.GetLayer()

            # read aoi
            ds_aoi = ogr.Open(aoi_path)
            layer_aoi = ds_aoi.GetLayer()

            # create output file
            driver = ogr.GetDriverByName('ESRI Shapefile')
            shapefile_clipped = os.path.splitext(shapefile_reprojected)[0] + '_clipped.shp'
            if os.path.exists(shapefile_clipped):
                driver.DeleteDataSource(shapefile_clipped)
            ds = driver.CreateDataSource(shapefile_clipped)
            input_layer_defn = layer_detfoo.GetLayerDefn()
            out_layer = ds.CreateLayer(input_layer_defn.GetName(), layer_detfoo.GetSpatialRef())

            layer_detfoo.Clip(layer_aoi, out_layer)

            gml_list.append(shapefile_clipped)

        else:

            gml_list.append(out_shape)

    return gml_list
