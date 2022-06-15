# encoding: utf-8

"""
this script used to pre-processing the dataset of MODIS.

goals:
    -- organize the scattered datas from April 2002 to Dec 2021 of MODIS into a complete, clean NC file. named "MODIS_GLOBAL_PRE_PROCESS.nc".
    which including : 
        -- normalized differenced vegetation index (NDVI).
        -- enhanced vegetation index (EVI).
    the file named "MODIS_GLOBAL_PRE_PROCESS.nc"

output:
    -- MODIS_EVI_GLOBAL_PRE_PROCESS.nc
    -- MODIS_NDVI_GLOBAL_PRE_PROCESS.nc
"""

# import modules
from itertools import count
import os
import netCDF4 as nc
from netCDF4 import date2num
from pyhdf.SD import SD, SDC
import numpy as np
import datetime


# 将结果保存成NC文件
def write_nc(dict_vars, rectangle, out_name, in_time="", scope="Gobal", sp_re=0.05, iunits="mm"):
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

    len_lon = len(np.arange(xmin, xmax+0.001, sp_re))
    len_lat = len(np.arange(ymin, ymax+0.001, sp_re))

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
        vars[key] = nc_data.createVariable(key, np.int16, ["time", "lat", "lon"])

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
    nc_data.geospatial_resolution = "0.05"
    nc_data.mask = scope

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
    lons[:] = np.arange(xmin, xmax+0.001, sp_re).tolist()
    lats[:] = np.arange(ymin, ymax+0.001, sp_re).tolist()
    print(lats)
    
    
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
        value.units = iunits    # 默认单位是mm
        value.Units = iunits
        value.fill_value = -9999.0
        value.missing_value = -9999.0
        value[:, :, :] = dict_vars[key].copy()      # 给变量赋值

            
    # 保存
    print("start to save file")
    nc_data.close()
    print("DONE!")
    
    return 0


# 矩阵上下颠倒
def flip(in_array):
    shp = in_array.shape
    rows = shp[0]
    for row in range(rows//2):
        temp = in_array[row, :].copy()
        in_array[row, :] = in_array[rows-1-row, :]
        in_array[rows-1-row, :] = temp

    return in_array


# 对modis数据进行整合
def organize_modis_file(top_path:str, txt_path:str):
    """_summary_

    Args:
        top_path (str): the file that modis datasets exist.
        txt_path (str): in order to make sure the correction of output, we need to jugle the date.
    return:
        --MODIS_GLOBAL_PRE_PROCESS.nc
    """
    # 从txt文件中提取date，并以键值对的形式保存在字典中
    dict_date = {}   # 读取每个文件所属的时间
    with open(txt_path, "r") as f:
        readlines = f.readlines()
        # 倒序堆栈
        readlines.reverse()
        for line in readlines:
            print(line[69:107])
            print(line[58:65])
            dict_date[line[69:107]] = datetime.datetime.strptime(line[58:65], "%Y.%m")
    

    # 遍历文件夹中所有文件
    walk = os.walk(top_path)
    for dirpath, dirnames, filenames in walk:
        print("dirpath:", dirpath)
        print("dirnames:", dirnames)
        print("filenames:", filenames)
        
        modis_name = [] # 读取文件中的每个modis文件，并存储到list中
        for filename in filenames:
            modis_name.append(filename)
    
    # 开辟存放evi和ndvi的空间
    evi_matrix = np.zeros(shape=(237, 3600, 7200), dtype=np.int16)
    
    count = 0
    for key in dict_date.keys():
        if key == modis_name[count]:
            # 说明该时间所对应的文件找到了
            path = os.path.join(top_path, modis_name[count])
            try:
                hdf_data = SD(path, mode=SDC.READ)      # 从2002年开始读取文件
            # 读取evi数据集中的数据
                evi = hdf_data.select("CMG 0.05 Deg Monthly EVI")
                evi_matrix[count, :, :] = flip(evi.get())
                print(count)
                count += 1
            except:
                print("*"*20)
                print("something wrong:", path)
                print("*"*20)
    
    # 处理缺失值
    evi_matrix[evi_matrix[:, :, :] == -3000] = -9999
    
    # 将结果保存在字典中
    dict_vars = {}
    dict_vars["CMG_0.05_Deg_Monthly_EVI"] = evi_matrix
    
    return dict_vars
    
# 对modis数据进行整合
def organize_modis_file(top_path:str, txt_path:str):
    """_summary_

    Args:
        top_path (str): the file that modis datasets exist.
        txt_path (str): in order to make sure the correction of output, we need to jugle the date.
    return:
        --MODIS_GLOBAL_PRE_PROCESS.nc
    """
    # 从txt文件中提取date，并以键值对的形式保存在字典中
    dict_date = {}   # 读取每个文件所属的时间
    with open(txt_path, "r") as f:
        readlines = f.readlines()
        # 倒序堆栈
        readlines.reverse()
        for line in readlines:
            print(line[69:107])
            print(line[58:65])
            dict_date[line[69:107]] = datetime.datetime.strptime(line[58:65], "%Y.%m")
    

    # 遍历文件夹中所有文件
    walk = os.walk(top_path)
    for dirpath, dirnames, filenames in walk:
        print("dirpath:", dirpath)
        print("dirnames:", dirnames)
        print("filenames:", filenames)
        
        modis_name = [] # 读取文件中的每个modis文件，并存储到list中
        for filename in filenames:
            modis_name.append(filename)
    
    # 开辟存放evi和ndvi的空间
    ndvi_matrix = np.zeros(shape=(237, 3600, 7200), dtype=np.int16)
    
    count = 0
    for key in dict_date.keys():
        if key == modis_name[count]:
            # 说明该时间所对应的文件找到了
            path = os.path.join(top_path, modis_name[count])
            try:
                hdf_data = SD(path, mode=SDC.READ)      # 从2002年开始读取文件
            # 读取evi数据集中的数据
                ndvi = hdf_data.select("CMG 0.05 Deg Monthly NDVI")
                ndvi_matrix[count, :, :] = flip(ndvi.get())
                print(count)
                count += 1
            except:
                print("*"*20)
                print("something wrong:", path)
                print("*"*20)
    
    # 处理缺失值
    ndvi_matrix[ndvi_matrix[:, :, :] == -3000] = -9999
    
    # 将结果保存在字典中
    dict_vars = {}
    dict_vars["CMG_0.05_Deg_Monthly_NDVI"] = ndvi_matrix
    
    return dict_vars    
    
    
            
if __name__=="__main__":
    top_path = "D:\\GRACE_ini_datasets\\MODIS_datasets"
    txt_path = "D:\\wget\\01download.txt"
    # 输出文件的名字
    
    out_name1 = "D:\\GRACE_ini_datasets\\MODIS_datasets\\MODIS_EVI_GLOBAL_PRE_PROCESS.nc"  # evi
    
    out_name2 = "D:\\GRACE_ini_datasets\\MODIS_datasets\\MODIS_NDVI_GLOBAL_PRE_PROCESS.nc"  # ndvi
    # 全球的范围
    rect = [-179.975, 179.975, -89.975, 89.975]
    
    dict_vars = organize_modis_file(top_path, txt_path)
    
    # 保存
    write_nc(dict_vars, rect, out_name2)
    
    
    

    