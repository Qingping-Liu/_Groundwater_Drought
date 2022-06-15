# encoding: utf-8

"""
this script used to pre-processing the dataset of GLDAS.

goals:
    -- organize the scattered datas from April 2002 to Dec 2021 of GLDAS into a complete, clean NC file. named "GLDAS_GLOBAL_PRE_PROCESS.nc".
    -- calculate the monthly variation of soil water, canopy water, snow equivalent water, surface runoff depth relative to the average of 2004 to 2009. 
    which including : 
        -- soil moisture storage anomaly(SMSA).
        -- snow water equivalent anomaly(SWEA).
        -- canopy water storage anomaly(CWSA).
        -- surface water storage anomaly(SWSA).
    the file named "GLDAS_GLOBAL_PRO_PROCESS_ANOMALIES.nc"

output:
    -- GLDAS_GLOBAL_PRE_PROCESS.nc
    -- GLDAS_GLOBAL_PRO_PROCESS_ANOMALIES.nc
"""

# import modules
from cmath import nan
import os
import netCDF4 as nc
from netCDF4 import date2num, num2date
import numpy as np
import datetime


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
        if key == "SWE_inst":
            value.standard_name = "surface_snow_amount"
            value.long_name = "snow depth water equivalent"
            value.units = "kg m-2"
            value.fill_value = -9999.0
            value.missing_value = -9999.0
            value[:, :, :] = dict_vars[key].copy()       # 给变量赋值
            
        elif key == "SoilMoi0_10cm_inst":
            value.standard_name = "soil_moisture"
            value.long_name = "soil moisture content"
            value.units = "kg m-2"
            value.fill_value = -9999.0
            value.missing_value = -9999.0
            value[:, :, :] = dict_vars[key].copy()       # 给变量赋值
            
        elif key == "SoilMoi10_40cm_inst":
            value.standard_name = "soil_moisture"
            value.long_name = "soil moisture content"
            value.units = "kg m-2"
            value.fill_value = -9999.0
            value.missing_value = -9999.0
            value[:, :, :] = dict_vars[key].copy()       # 给变量赋值
            
        elif key == "SoilMoi40_100cm_inst":
            value.standard_name = "soil_moisture"
            value.long_name ="soil moisture content"
            value.units = "kg m-2"
            value.fill_value = -9999.0
            value.missing_value = -9999.0
            value[:, :, :] = dict_vars[key].copy()       # 给变量赋值
            
        elif key == "SoilMoi100_200cm_inst":
            value.standard_name = "soil_moisture"
            value.long_name ="soil moisture content"
            value.units = "kg m-2"
            value.fill_value = -9999.0
            value.missing_value = -9999.0
            value[:, :, :] = dict_vars[key].copy()       # 给变量赋值
            
        elif key == "CanopInt_inst":
            value.standard_name = "canopy_water"
            value.long_name = "plant canopy surface water"
            value.units = "kg m-2"
            value.fill_value = -9999.0
            value.missing_value = -9999.0
            value[:, :, :] = dict_vars[key].copy()       # 给变量赋值
            
        elif key == "Qs_acc":
            value.standard_name = "storm_surface_runoff"
            value.long_name = "storm surface runoff"
            value.units = "kg m-2 per 3h"
            value.fill_value = -9999.0
            value.missing_value = -9999.0
            value[:, :, :] = dict_vars[key].copy()       # 给变量赋值
            
        elif key == "Rainf_f_tavg":
            value.standard_name = "total_precipitation"
            value.long_name = "total precipitation rate"
            value.units = "kg m-2 s-1"
            value.fill_value = -9999.0
            value.missing_value = -9999.0
            value[:, :, :] = dict_vars[key].copy()       # 给变量赋值
            
        elif key == "Evap_tavg":
            value.standard_name = "total_evapotranspiration"
            value.long_name = "total evapotranspiration rate"
            value.units = "kg m-2 s-1"
            value.fill_value = -9999.0
            value.missing_value = -9999.0
            value[:, :, :] = dict_vars[key].copy()     # 给变量赋值
    
        else:
            # 自动添加属性
            box = key.split(" ")
            standard_name = ""
            for i in range(len(box)-1):
                standard_name = standard_name + str(box[i]) + "_" 
            standard_name = standard_name + str(box[-1])
            long_name = key 
            value.standard_name = standard_name
            value.long_name = long_name
            value.units = "kg m-3"
            value.fill_value = -9999.0
            value.missing_value = -9999.0
            value[:, :, :] = dict_vars[key].copy()      # 给变量赋值
            
    
    
    
    # 保存
    print("start to save file")
    nc_data.close()
    print("DONE!")

    return 0

# 对glads数据进行整合
def organize_gldas_file(top_path: str):
    """_summary_
    Args:
        top_path (str): the file that gldas datasets exist.
        name_date (str): in order to make sure the correction of output, we need to jugle the date.
    return: default "GLDAS_GLOBAL_PRE_PROCESS.nc"
    """
    
    # 遍历top_path的所有文件
    walk = os.walk(top_path)
    for dirpaths, dirnames, filenames in walk:
        print("dirpath:", dirpaths)
        print("dirnames:", dirnames)
        print("filenames", filenames)
        gldas_date = []     # 读取文件中每个gldas的时间，并存储到list中
        for filename in filenames:
            box = filename.split(sep=".", maxsplit=5)
            str_date = box[1].split("A", 2)[1]
            gldas_date.append(datetime.datetime.strptime(str_date, "%Y%m"))
        
    # 检查是否有文件缺失
    print("lentime:", len(gldas_date))
    print("gldastime:", gldas_date)


    # 从2002-04月开始读取数据,
    start_time = datetime.datetime.strptime("2004-01", "%Y-%m")
    end_time = datetime.datetime.strptime("2009-12", "%Y-%m")


    # 9个变量：土壤水（4层）， 冠层水，雪水当量，径流深，降水、潜在蒸散发
    snow_datas = np.zeros(shape=(237, 600, 1440), dtype=np.float32)
    soil_10_datas = np.zeros(shape=(237, 600, 1440), dtype=np.float32)
    soil_40_datas = np.zeros(shape=(237, 600, 1440), dtype=np.float32)
    soil_100_datas = np.zeros(shape=(237, 600, 1440), dtype=np.float32)
    soil_200_datas = np.zeros(shape=(237, 600, 1440), dtype=np.float32)
    canop_datas = np.zeros(shape=(237, 600, 1440), dtype=np.float32)
    surwater_datas = np.zeros(shape=(237, 600, 1440), dtype=np.float32)
    precipitation_datas = np.zeros(shape=(237, 600, 1440), dtype=np.float32)
    evapotranspiration_datas = np.zeros(shape=(237,600, 1440), dtype=np.float32)
    
    
    # 文件命名格式
    temple_head = "GLDAS_NOAH025_M.A"
    temple_tail = ".021.nc4.SUB.nc4"
    
    count = 0
    for time in gldas_date:
        str_time = datetime.datetime.strftime(time, "%Y%m")
        name = temple_head + str_time + temple_tail
        filepath = os.path.join(top_path, name)
        nc_datas = nc.Dataset(filepath, mode="r")
        snow_datas[count, :, :] = nc_datas.variables["SWE_inst"][:, :]      # 雪水当量
        soil_10_datas[count, :, :] = nc_datas.variables["SoilMoi0_10cm_inst"][:, :]     # 土壤水0-10
        soil_40_datas[count, :, :] = nc_datas.variables["SoilMoi10_40cm_inst"][:, :]    # 土壤水当量10-40
        soil_100_datas[count, :, :] = nc_datas.variables["SoilMoi40_100cm_inst"][:, :]   # 土壤水40-100
        soil_200_datas[count, :, :] = nc_datas.variables["SoilMoi100_200cm_inst"][:, :]  # 土壤水100-200
        canop_datas[count, :, :] = nc_datas.variables["CanopInt_inst"][:, :]    # 植物冠层水
        surwater_datas[count, :, :] = nc_datas.variables["Qs_acc"][:, :]    # 地表径流深
        precipitation_datas[count, :, :] = nc_datas.variables["Rainf_f_tavg"][:, :]  # 降雨深
        evapotranspiration_datas[count, :, :] = nc_datas.variables["Evap_tavg"][:, :]     # 蒸散发量
        count += 1 
    

    # 字典的值
    keys = ["SWE_inst","SoilMoi0_10cm_inst","SoilMoi10_40cm_inst","SoilMoi40_100cm_inst","SoilMoi100_200cm_inst","CanopInt_inst","Qs_acc","Rainf_f_tavg","Evap_tavg"]
    # 将结果保存在字典中
    dict_vars = {}
    dict_vars["SWE_inst"] = snow_datas
    dict_vars["SoilMoi0_10cm_inst"] = soil_10_datas
    dict_vars["SoilMoi10_40cm_inst"] = soil_40_datas
    dict_vars["SoilMoi40_100cm_inst"] = soil_100_datas
    dict_vars["SoilMoi100_200cm_inst"] = soil_200_datas
    dict_vars["CanopInt_inst"] = canop_datas
    dict_vars["Qs_acc"] = surwater_datas
    dict_vars["Rainf_f_tavg"] = precipitation_datas
    dict_vars["Evap_tavg"] = evapotranspiration_datas
    
    return dict_vars


# 根据处理好的gldas文件，计算2004年-2009年的均值，并计算相对变化量
def compute_anomalies(nc_path: str, start_time: str, end_time: str):
    """_summary_

    Args:
        nc_path (str): the file format is NETCDF.
        start_time (str): to calculate the relative variations of each variable.
        end_time (str): to calculate the average of each variable.
        return:
            --"GLDAS_GLOBAL_PRO_PROCESS_ANOMALIES.nc"
    """
    # 读取整个数据集
    nc_datasets = nc.Dataset(nc_path, mode="r")
    # 读取变量
    lons = nc_datasets.variables["lon"][:]
    lats = nc_datasets.variables["lat"][:]
    times = nc_datasets.variables["time"][:]
    calend = nc_datasets.variables["time"].calendar # 日历
    unit = nc_datasets.variables["time"].units  #起始时间        
    # 把时间格式转换成datetime格式
    start_time = datetime.datetime.strptime(start_time, "%Y-%m")
    end_time = datetime.datetime.strptime(end_time, "%Y-%m")
    times = num2date(times, unit, calend)
    

    
    # 计算均值
    average_soil = np.zeros(shape=(len(lats),len(lons)))
    average_canop = np.zeros(shape=(len(lats),len(lons)))
    average_surface = np.zeros(shape=(len(lats),len(lons)))
    average_snow = np.zeros(shape=(len(lats),len(lons)))
    count = 0
    for i in range(len(times)):
        if start_time <= times[i] and times[i] <= end_time:
            # 开始处理哪一年的数据
            print("%s have been executed" % times[i].strftime("%Y-%m"))
            # 处理缺失值
            soil_10 = np.array(nc_datasets.variables["SoilMoi0_10cm_inst"][i, :, :])          # 0 - 10 cm土壤水
            soil_40 = np.array(nc_datasets.variables["SoilMoi10_40cm_inst"][i, :, :])         # 10 - 40 cm土壤水
            soil_100 = np.array(nc_datasets.variables["SoilMoi40_100cm_inst"][i, :, :])        # 40 - 100 cm土壤水
            soil_200 = np.array(nc_datasets.variables["SoilMoi100_200cm_inst"][i, :, :])      # 100-200 cm土壤水
            snow = np.array(nc_datasets.variables["SWE_inst"][i, :, :])                       # 雪水当量
            canop = np.array(nc_datasets.variables["CanopInt_inst"][i, :, :])                 # 冠层水当量
            surface = np.array(nc_datasets.variables["Qs_acc"][i, :, :])                      # 径流深     
         
            # -9999.0填充为 np.nan
            soil_10[soil_10[:, :] == -9999.0] = np.nan
            soil_40[soil_40[:, :] == -9999.0] = np.nan
            soil_100[soil_100[:, :] == -9999.0] = np.nan
            soil_200[soil_200[:, :] == -9999.0] = np.nan
            
            snow[snow[:, :] == -9999.0] = np.nan
            canop[canop[:, :] == -9999.0] = np.nan
            surface[surface[:, :] == -9999.0] = np.nan
            
            # 累加土壤水、 冠层水、雪水、径流深
            average_soil = average_soil + soil_10 + soil_40 + soil_100 + soil_200
            average_canop = average_canop + canop
            average_snow = average_snow + snow
            average_surface = average_surface + surface * 8 * 30
            count += 1
            
            
    # 查看是否有问题
    print(np.count_nonzero(np.isnan(average_soil)))
    
    # 计算多年平均值, 此时的缺失值是np.nan
    average_soil = average_soil / count
    average_canop = average_canop / count
    average_snow = average_snow / count
    average_surface = average_surface / count
    
    # 准备存放gldas相对变化量的array
    anomalies_soil = np.empty(shape=(len(times), len(lats), len(lons)))
    anomalies_canop = np.empty(shape=(len(times), len(lats), len(lons)))
    anomalies_snow = np.empty(shape=(len(times), len(lats), len(lons)))
    anomalies_surface = np.empty(shape=(len(times), len(lats), len(lons)))
    
    for i in range(len(times)):
        print("the file of %s have been executed." % times[i].strftime("%Y-%m"))
        soil_10 = np.array(nc_datasets.variables["SoilMoi0_10cm_inst"][i, :, :])          # 0 - 10 cm土壤水
        soil_40 = np.array(nc_datasets.variables["SoilMoi10_40cm_inst"][i, :, :])         # 10 - 40 cm土壤水
        soil_100 = np.array(nc_datasets.variables["SoilMoi40_100cm_inst"][i, :, :])        # 40 - 100 cm土壤水
        soil_200 = np.array(nc_datasets.variables["SoilMoi100_200cm_inst"][i, :, :])      # 100-200 cm土壤水
        snow = np.array(nc_datasets.variables["SWE_inst"][i, :, :])                       # 雪水当量
        canop = np.array(nc_datasets.variables["CanopInt_inst"][i, :, :])                 # 冠层水当量
        surface = np.array(nc_datasets.variables["Qs_acc"][i, :, :])                      # 径流深 

        # 计算相对变化量
        anomalies_soil[i, :, :] = soil_10 + soil_40 + soil_100 + soil_200 - average_soil
        anomalies_canop[i, :, :] = canop - average_canop
        anomalies_snow[i, :, :] = snow - average_snow
        anomalies_surface[i, :, :] = surface - average_surface
    
    
    # 计算完成后，np.nan值填充为-9999.0
    anomalies_soil[np.isnan(anomalies_soil)] = -9999.0
    anomalies_canop[np.isnan(anomalies_canop)] = -9999.0
    anomalies_snow[np.isnan(anomalies_snow)] = -9999.0
    anomalies_surface[np.isnan(anomalies_surface)] = -9999.0
        
    # 保存在字典中
    anoma_dic = {}
    anoma_dic["soil moisture content anomalies"] = anomalies_soil
    anoma_dic["canopy surface water anomalies"] = anomalies_canop
    anoma_dic["snow depth water anomalies"] = anomalies_snow
    anoma_dic["surface runoff anomalies"] = anomalies_surface
    
    
    return anoma_dic
    


if __name__ == "__main__":
    nc_path = "D:\\GRACE_ini_datasets\\initial_INPUT_data\\_Global\\GLDAS_GLOBAL_PRE_PROCESS.nc"
    start_time = "2004-01"  
    end_time = "2009-12"
    out_name = "D:\\GRACE_ini_datasets\\initial_INPUT_data\\_Global\\GLDAS_GLOBAL_PRE_PROCESS_ANOMALIES_v2.nc"
    
    # gldas 的数据范围
    rect = [-179.875, 179.875, -59.875, 89.875]     # xmin, xmax, ymin, ymax
    
    # 所有数据保存在字典里
    dict_vars = compute_anomalies(nc_path, start_time, end_time)
    
    # 将数据保存成一个nc文件
    write_nc(dict_vars, rect, out_name)