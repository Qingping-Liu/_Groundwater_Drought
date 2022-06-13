# encoding: utf-8

"""
this scripts used to clip the groundwater storage variation from April 2002 to December 2021 according to the boundary of NCP.

gloas:
    -- clip the dataset of groundwater storage anamolies according to the shape file of NCP
input:
    -- GLOBAL_Terrestial_Monthly_Groundwater_Storage_Anomalies.nc
    -- NCP.shp

output
    -- NCP_Terrestial_Monthly_Groundwater_Storage_Anomalies.nc
"""

# import modules
import netCDF4 as nc
from netCDF4 import date2num, num2date
import numpy as np
import datetime
import os

from write_nc import  write_nc

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


# 矩阵平移， 确保栅格与栅格对应
def tranform_maxtrix(in_array):
    """
    :param in_array: 2维
    :return:
    """
    shp = in_array.shape
    columns = shp[1]  # 一共多少列
    col = columns // 2

    temp = in_array[:, 0: col].copy()
    in_array[:, :col] = in_array[:, col:]  # 180°-360°的部分平移到-180°-0°
    in_array[:, col:] = temp

    return in_array*10


# 根据NCP的boundary裁剪数据
def clip_data(nc_path="", bound="", keyword="groundwater storage anomalies", left_low_point=(-179.875,59.875), region_clip=False):
    """_summary_

    Args:
        nc_path (str): the file of groundwater storage anomalies.
        boudary (str): the boundary of NCP
    """
    # import module
    import arcpy

    # 设置工作空间
    path = "D:\\temp_202203"
    arcpy.env.workspace = "D:\\temp_202203"
    # 设置坐标系
    sr = arcpy.SpatialReference(4326)
    arcpy.env.outCoordinateSystem = sr
    # 允许覆盖
    arcpy.env.overwriteOutput = True


    # 读取nc数据
    nc_dataset = nc.Dataset(nc_path, mode="r")
    # 读取相关变量
    lats = nc_dataset.variables["lat"][:]
    lons = nc_dataset.variables["lon"][:]
    # 把时间处理成datetime格式
    times = nc_dataset.variables["time"][:]
    unit = nc_dataset.variables["time"].units
    calend = nc_dataset.variables["time"].calendar
    times = num2date(times, unit, calend)

    #创建一个空的三维数组
    dic_var = {}
    with arcpy.da.SearchCursor(bound, ["SHAPE@", "Name"]) as rows:
        for row in rows:
            shape = row[0]  # 几何对象类型
            name = row[1]  # shape的名字

            new_vars = []  # 这里的值可能有问题; 数据类型
            # 根据要素逐个clip
            for i in range(len(times)):
                print("region: %s\tstart to process: %s" % (name, times[i].strftime("%Y-%m")))
                vars = np.array(nc_dataset.variables[keyword][i, :, :])

                # # ！！！！！！！！GRACE 要平移
                # matrix = tranform_maxtrix(flip_array(vars))   # 翻转

                matrix = flip_array(vars)

                yearMonth = times[i].strftime("%Y-%m") + ".tif"
                raster = arcpy.NumPyArrayToRaster(matrix, arcpy.Point(left_low_point[0], left_low_point[1]), x_cell_size=0.25, y_cell_size=0.25, value_to_nodata=-9999.0)      # 角点非常重要！！！
                arcpy.DefineProjection_management(raster, sr)
                try:
                    arcpy.Clip_management(in_raster=raster, out_raster=yearMonth, in_template_dataset=shape, clipping_geometry=True)

                except:
                    print("%s have been meeting a problem." % name)
                    Flag = False
                    break
                else:
                    Flag = True
                    # 获取裁剪后栅格的角标
                    des = arcpy.Describe(yearMonth)
                    xmin, xmax = des.Extent.XMin, des.Extent.XMax
                    ymin, ymax = des.Extent.YMin, des.Extent.YMax

                    rect = [xmin, xmax - 0.001, ymin, ymax - 0.001]

                    new_matrix = arcpy.RasterToNumPyArray(os.path.join(path, yearMonth), nodata_to_value=-9999.0)
                    new_vars.append(flip_array(new_matrix))  # 翻转回来


            if region_clip and Flag:
                outname = "D:\\GRACE_ini_datasets\\Final_data\\Regional_DSI" + "\\" + name + ".nc"
                # 以nc形式保存
                dic_var["DSI"] = np.asarray(new_vars)
                write_nc(dic_var, rect, outname, scope="NORTH CHINA PLAIN", sp_re=0.25, iunits=" ")

            del new_vars
    # 定义一个字典
    dic_var = {}
    dic_var[keyword] = new_vars

    return dic_var, rect



def sum_gldas(nc_path):
    # 读取nc数据
    nc_dataset = nc.Dataset(nc_path, mode="r")
    # 读取相关变量
    lats = nc_dataset.variables["lat"][:]
    lons = nc_dataset.variables["lon"][:]
    # 把时间处理成datetime格式
    times = nc_dataset.variables["time"][:]
    unit = nc_dataset.variables["time"].units
    calend = nc_dataset.variables["time"].calendar
    times = num2date(times, unit, calend)

    # 一个空的字典
    vars = np.zeros(shape=(len(times), len(lats), len(lons)), dtype=np.float32)



    keys = ["canopy surface water anomalies", "snow depth water anomalies", "soil moisture content anomalies", "surface runoff anomalies"]
    for i in range(len(times)):
        c = nc_dataset.variables[keys[0]][i, :, :]
        sn = np.array(nc_dataset.variables[keys[1]][i, :, :])
        so = np.array(nc_dataset.variables[keys[2]][i, :, :])
        su = np.array(nc_dataset.variables[keys[3]][i, :, :])

        # 替换填充值
        c[c[:, :] == -9999.0] = np.nan
        sn[sn[:, :] == -9999.0] = np.nan
        so[so[:, :] == -9999.0] = np.nan
        su[su[:, :] == -9999.0] = np.nan

        vars[i, :, :] = c + sn + so + su

    # 替换填充值
    vars[np.isnan(vars)] = -9999.0


    dic_var = {}
    dic_var["terrestrial water storage anamolies"] = vars

    # 数据范围
    xmin, xmax = float(nc_dataset.variables["lon"].valid_min), float(nc_dataset.variables["lon"].valid_max)
    ymin, ymax = float(nc_dataset.variables["lat"].valid_min), float(nc_dataset.variables["lat"].valid_max)

    rect = [xmin, xmax, ymin, ymax]

    return dic_var, rect

if __name__ == "__main__":

    # top_path = "D:\\GRACE_ini_datasets\\mediate_data"
    # filename = "GLOBAL_Terrestial_Monthly_Groundwater_Storage_Anomalies_v2.nc"
    # nc_path = os.path.join(top_path, filename)
    #
    # # 华北平原
    # path = r"C:\Users\11072\Desktop\GRACE_dataset\huebei\111.gdb"
    # name = "NCP_boundary"
    # bound = path + "\\" + name
    #
    # dic_var, rect = clip_data(nc_path, bound)
    # 文件名
    # outname = "NCP_Terrestial_Monthly_Groundwater_Storage_Anomalies.nc"
    # write_nc(dic_var, rect, outname, sp_re=0.25)


    # #全球地下水水库   428 1106
    #path = r"C:\Users\11072\Desktop\GRACE_dataset\huebei\111.gdb"
    #global_name = "Global_Groundwater_Reservoir"
    #global_bound = path + "\\" + global_name
    #dic_var, rect = clip_data(nc_path, global_bound)
    # 文件名
    #outname1 = "Global_Groundwater_Reservoir_Monthly_Groundwater_Storage_Anomalies.nc"
    #write_nc(dic_var, rect, outname1)
    
    
    # clip MODIS in NCP
    #path = "C:\\Users\\11072\\Desktop\\GRACE_dataset\\huebei\\111.gdb"
    #name = "NCP_boundary"
    #bound_path= os.path.join(path, name)

    #keyword="CMG 0.25 Deg Monthly EVI"
    #dict_var, rect = clip_data(nc_path, bound_path, keyword)

    # 输出的文件名
    #outname = "NCP_Terrestial_Monthly_EVI.nc"
    #write_nc(dict_vars=dict_var, rectangle=rect, out_name=outname, scope="NORTH CHINA PLAIN", iunits="--", sp_re=0.25, dtype=np.int16)
    



    # # clip MODIS in GLOBAL GROUNDWATER RESERVOIRS
    # path = "C:\\Users\\11072\\Desktop\\GRACE_dataset\\huebei\\111.gdb"
    # name = "Global_Groundwater_Reservoir"
    # global_bound = os.path.join(path, name)
    #
    # keyword = "CMG 0.25 Deg Monthly EVI"
    # dict_var, rect = clip_data(nc_path, global_bound, keyword)
    #
    # # 输出文件
    # outname = "Global_Groundwater_Reservoir_Terrestial_Monthly_EVI.nc"
    # write_nc(dict_var, rect, outname, "", "GLOBAL GROUNDWATER RESERVOIRS", iunits="--", sp_re=0.25)



    # # 提取GLDAS的数据
    # from write_nc import write_nc
    #
    # top_path = r"D:\GRACE_ini_datasets\initial_INPUT_data\_Global"
    # name = "GLDAS_Terrestrial_Water_Storage_Anamolies_v2.nc"
    # nc_path = os.path.join(top_path, name)
    #
    # # # 先计算gladas的和
    # # dic_var, rect = sum_gldas(nc_path)
    # # outname = "GLDAS_Terrestrial_Water_Storage_Anamolies_v2.nc"
    #
    # # 华北平原
    # path = r"C:\Users\11072\Desktop\GRACE_dataset\huebei\111.gdb"
    # name = "NCP_boundary"
    # bound = path + "\\" + name
    #
    # keyword = "terrestrial water storage anamolies"
    # dic_var, rect = clip_data(nc_path, bound, keyword)
    #
    # outname = "NCP_GLDAS_Terrestial_Storage_Anomalies.nc"
    # write_nc(dic_var, rect, outname, sp_re=0.25)



    # # # 提取GRACE的数据
    # from write_nc import write_nc
    # top_path = r"D:\GRACE_ini_datasets\initial_INPUT_data\_Global"
    # name = "BNU_GRACE_GRACE-FO_RL06_Mascons_Linear_interpolation_v02.nc"
    # nc_path = os.path.join(top_path, name)
    #
    # # 华北平原
    # path = r"C:\Users\11072\Desktop\GRACE_dataset\huebei\111.gdb"
    # name = "NCP_boundary"
    # bound = path + "\\" + name
    #
    # keyword = "lwe_thickness"
    # dic_var, rect = clip_data(nc_path, bound, keyword)
    #
    # outname = "NCP_GRACE_Terrestial_Storage_Anomalies.nc"
    # write_nc(dic_var, rect, outname, sp_re=0.25)



    # 按省裁剪要素
    top_path = r"D:\GRACE_ini_datasets\Final_data"
    filename = "NCP_GWSA-DSI_v2.nc"
    nc_path = os.path.join(top_path, filename)
    keyword = "groundwater storage anomalies DSI"

    # 华北平原
    path = r"C:\Users\11072\Desktop\GRACE_dataset\huebei\111.gdb"
    name = "NCP_Provinces"
    bound = path + "\\" + name
    dic_var, rect = clip_data(nc_path, bound, keyword, left_low_point=(112.5, 32.0),region_clip=True)