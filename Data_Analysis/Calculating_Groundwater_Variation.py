# encoding: utf-8

"""
this scripts used to compute the groundwater storage variation from April 2002 to December 2021 in the NCP.

gloas:
    -- calculating the groundwater storage variation according to the following expression:
                     GWSA = TWSA * 10 - SMSA - SWEA - CWSA - SWSA
    -- where:
            GWSA: groundwater storage anomalies [kg m-2].
            TWSA: terrestrial water storage anomalies [cm].
            SMSA: soil moisture storage anomalies [kg m-2].
            SWEA: snow water equivalent anomalies [kg m-2].
            CWSA: canopy water storage anomalies [kg m-2]
            SWSA: surface water storage anomalies [kg m-2].

input:
    -- BNU_GRACE_GRACE-FO_RL06_Mascons_interpolation_v02.nc
    -- GLDAS_GLOBAL_PRO_PROCESS_ANOMALIES.nc

output
    -- GLOBAL_Terrestial_Monthly_Groundwater_Storage_Anomalies.nc

note:
    -- before computing the values of GWSA, we need to make sure the units of related variables is consistent.
    -- the latitude of gldas varies from -60 to 90, doesn't consider the Antarctica, which isn't consistent with the dataset of GRACE.
    
"""


# import modules
import netCDF4 as nc
from netCDF4 import date2num, num2date
import numpy as np
import datetime
import os 



# 矩阵平移， 确保栅格与栅格对应
def tranform_maxtrix(in_array):
    shp = in_array.shape
    columns = shp[2]        # 一共多少列
    col = columns //2 
    
    temp = in_array[:, :, 0: col].copy()    
    in_array[:, :, :col] = in_array[:, :, col:] # 180°-360°的部分平移到-180°-0°
    in_array[:, :, col: ] = temp
    
    return in_array



def computing_GWSA(grace_path="", gldas_path=""):
    """_summary_

    Args:
        grace_path (str): the file provides the dataset of TWSA
        gldas_path (str): the file provides the datasets of SMSA, SWEA, CWSA, SWSA.
    """
    # 先读取grace的数据
    nc_grace = nc.Dataset(grace_path, mode="r")
    grace_lwe = np.array(nc_grace.variables["lwe_thickness"][: , 120:, :])      # note: 注意范围
    # 将grace的数据180度到360度之间平移到前面
    TWSA = tranform_maxtrix(grace_lwe)
    
    # 再读取gldas数据
    nc_gldas = nc.Dataset(gldas_path, mode="r")
    lats = nc_gldas.variables["lat"][:]
    lons = nc_gldas.variables["lon"][:]
    
    times = nc_gldas.variables["time"][:]
    unit = nc_gldas.variables["time"].units
    calend = nc_gldas.variables["time"].calendar
    # 把时间处理成datetime格式
    times = num2date(times, unit, calend)
    # 定义几个空的数组
    SMSA = np.empty(shape = (len(times), len(lats), len(lons)))

    vars = nc_gldas.variables.keys()
    print(vars)

    SMSA = np.array(nc_gldas.variables["soil moisture content anomalies"][:, :, :])
    CWSA = np.array(nc_gldas.variables["canopy surface water anomalies"][:, :, :])
    SWEA = np.array(nc_gldas.variables["snow depth water anomalies"][:, :, :])
    SWSA = np.array(nc_gldas.variables["surface runoff anomalies"][:, :, :])
    
    # 处理缺失值
    SMSA[SMSA[:, :, :] == -9999.0] = np.nan
    CWSA[CWSA[:, :, :] == -9999.0] = np.nan
    SWEA[SWEA[:, :, :] == -9999.0] = np.nan
    SWSA[SWSA[:, :, :] == -9999.0] = np.nan
    
    # 根据表达式计算
    GWSA = np.empty(shape=(len(times), len(lats), len(lons)))
    for i in range(len(times)):
        year_month = times[i].strftime("%Y-%m")
        print("the file of %s is processing..." % year_month)
        GWSA[i, :, :] = TWSA[i, :, :] * 10 - SMSA[i, :, :]  - CWSA[i, :, :]  - SWEA[i, :, :]  - SWSA[i, :, :] 
    
    # 填充nan值
    GWSA[np.isnan(GWSA)] = -9999.0
    
    
    # 定义一个空的字典
    dic_var = {}
    dic_var["groundwater storage anomalies"] = GWSA
    
    return dic_var


def write_nc(dict_vars, rectangle, out_name, in_time=""):
    """
    :param in_array:
    :param keyword:
    :param rectangle: x_min, x_max, y_min, y_max, consider 0.125
    :in_time: default 2002-04 to 2021-12
    :return:
    默认像元大小： 0.25
    """
    # create a new dataset
    print("start to create new nc file.")
    nc_data = nc.Dataset(out_name, mode="w")

    xmin = float(rectangle[0])
    xmax = float(rectangle[1])

    ymin = float(rectangle[2])
    ymax = float(rectangle[3])

    len_lon = len(np.arange(xmin, xmax+0.001, 0.25))
    len_lat = len(np.arange(ymin, ymax+0.001, 0.25))

    # 新建维度
    __ = nc_data.createDimension("lon", size=len_lon)
    __ = nc_data.createDimension("lat", size=len_lat)
    __ = nc_data.createDimension("time", None)

    # 创建变量
    lats = nc_data.createVariable("lat", np.float32, ["lat", ])
    lons = nc_data.createVariable("lon", np.float32, ["lon", ])
    times = nc_data.createVariable("time", np.float32, ["time", ])
    # 创建变量
    vars = {}       # 定义1个空字典，存储所有的变量 键为关键字，值为对应的变量
    for key in dict_vars.keys():
        vars[key] = nc_data.createVariable(key, np.float32, ["time", "lat", "lon"])

    # 设置全局属性
    filename = out_name.split(".")[0]
    nc_data.filename = str(filename)
    nc_data.institution =  "COLLEGE OF WATER SICENCE, BEIJING NORMAL UNIVERSITY"
    nc_data.creator_name = "QingPing Liu"

    nc_data.time_coverage_start = "2000-04-01 00:00:00"
    nc_data.time_coverage_end = "2021-12-31 00:00:00"
    nc_data.geospatial_lat_min = str(ymin)
    nc_data.geospatial_lat_max = str(ymax)
    nc_data.geospatial_lat_units = "degrees_north"
    nc_data.geospatial_lon_min = str(xmin)
    nc_data.geospatial_lon_max = str(xmax)
    nc_data.geospatial_lon_units = "degrees_east"
    nc_data.geospatial_resolution = "0.25"
    nc_data.mask = "Gobal"

    # 设置局部属性
    lats.valid_min = str(ymin)
    lats.valid_max = str(ymax)
    lats.units = "degrees_north"
    lats.long_name = "Latitude"
    lats.standard_name = "Latitude"

    # 设置局部属性
    lons.valid_min = str(xmin)
    lons.valid_max = str(xmax)
    lons.units = "degrees_east"
    lons.long_name = "Longitude"
    lons.standard_name = "Longitude"

    # 设置时间属性
    times.standard_name = "Time"
    times.long_name = "Time"
    times.calendar = "standard"
    times.units = "days since 2000-01-01 00:00:00"

    # 给变量赋值
    lons[:] = np.arange(xmin, xmax+0.001, 0.25).tolist()
    lats[:] = np.arange(ymin, ymax+0.001, 0.25).tolist()
    
    
    #默认时间范围是2002年4月至2021年12月
    date = []
    year = 2002
    while year <= 2021:
        month = 1
        while month <= 12:
            if year*100 + month >= 200204:
                year_month = str(year)+str(month)
                date.append(datetime.datetime.strptime(year_month, "%Y%m"))
            month += 1  # 结束条件
        year += 1 # 结束条件
    # 把时间格式处理为num格式
    num_date = date2num(date, units=times.units, calendar=times.calendar)
    
    # 给时间赋值
    times[:] = num_date[:]


    # 给每个变量设置属性
    for key, value in vars.items():
        # 自动添加属性
        box = key.split(" ")
        standard_name = ""
        for i in range(len(box)-1):
            standard_name = standard_name + str(box[i]) + "_" 
        standard_name = standard_name + str(box[-1])
        long_name = key 
        value.standard_name = standard_name
        value.long_name = long_name
        value.units = "mm"
        value.Units = "mm"
        value.fill_value = -9999.0
        value.missing_value = -9999.0
        value[:, :, :] = dict_vars[key].copy()      # 给变量赋值
            
    
    # 保存
    print("start to save file")
    nc_data.close()
    print("DONE!")

    return 0


if __name__ == "__main__":
    top_path = "D:\\GRACE_ini_datasets\\initial_INPUT_data\\_Global"
    gracename = "BNU_GRACE_GRACE-FO_RL06_Mascons_Linear_interpolation_v02.nc"
    gldasname = "GLDAS_GLOBAL_PRE_PROCESS_ANOMALIES_v2.nc"
    
    grace_path = os.path.join(top_path, gracename)
    gldas_path = os.path.join(top_path, gldasname)
    
    dic_var = computing_GWSA(grace_path, gldas_path)
    
    # 处理的范围
    # gldas 的数据范围
    rect = [-179.875, 179.875, -59.875, 89.875]     # xmin, xmax, ymin, ymax
    
    
    # 文件名
    outname = "GLOBAL_Terrestial_Monthly_Groundwater_Storage_Anomalies.nc"
    
    write_nc(dic_var, rect, outname)