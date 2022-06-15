# encoding: utf-8

import os
import numpy as np

# 矩阵上下翻转
def flip_array(in_array):
    """
    :param in_array:2维
    :return:
    """
    shp = in_array.shape
    rows = shp[0]
    for i in range(rows//2):
        temp = in_array[i, :].copy()
        in_array[i, :] = in_array[rows - 1 -i , :].copy()
        in_array[rows - 1 -i , :] = temp
    return in_array


# 将二维矩阵保存成tif文件
def write_tif(dic_var, rect, cell, outname, folder="D:\\GRACE_ini_datasets\\Final_data\\linear\\GIS_tif_v2"):
    """

    :param dic_var:
    :param low_left_point:
    :return: 返回tif文件
    """
    import arcpy
    # 设置工作环境
    arcpy.env.workspace = "D:\\temp_202203"
    arcpy.env.overwriteOutput = True
    arcpy.env.outCoordinateSystem = arcpy.SpatialReference(4326)

    sr = arcpy.SpatialReference(4326)

    box = outname.split(".")

    # 把结果保存成raster
    for key, value in dic_var.items():
        array = flip_array(value)
        raster = arcpy.NumPyArrayToRaster(in_array=array, lower_left_corner=arcpy.Point(rect[0]-cell/2, rect[2]-cell/2), x_cell_size=cell,  y_cell_size=cell, value_to_nodata=-9999.0)
        arcpy.DefineProjection_management(raster, sr)
        # 输出的名字
        box_name = box[0] + "_" +key + ".tif"
        filename = os.path.join(folder, box_name)
        raster.save(filename)

    print("Save as tiff file: successful!")
    return 0