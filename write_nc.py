# encoding: utf-8

import os
import netCDF4 as nc
from netCDF4 import date2num
import numpy as np
import datetime


# 将结果保存成NC文件
def write_nc(dict_vars, rectangle, out_name, in_time="", scope="Gobal", sp_re=0.05, iunits="mm", dtype=np.float32):
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

    len_lon = len(np.arange(xmin, xmax + 0.001, sp_re))
    len_lat = len(np.arange(ymin, ymax + 0.001, sp_re))

    # 新建维度
    __ = nc_data.createDimension("lon", size=len_lon)
    __ = nc_data.createDimension("lat", size=len_lat)
    __ = nc_data.createDimension("time", None)

    # 创建变量
    lats = nc_data.createVariable("lat", np.float32, ["lat", ])
    lons = nc_data.createVariable("lon", np.float32, ["lon", ])
    times = nc_data.createVariable("time", np.float32, ["time", ])
    # 创建变量
    vars = {}  # 定义1个空字典，存储所有的变量 键为关键字，值为对应的变量
    for key in dict_vars.keys():
        vars[key] = nc_data.createVariable(key, dtype, ["time", "lat", "lon"])

    # 设置全局属性
    filename = out_name.split(".")[0]
    nc_data.filename = str(filename)
    nc_data.institution = "COLLEGE OF WATER SICENCE, BEIJING NORMAL UNIVERSITY"
    nc_data.creator_name = "QingPing Liu"

    nc_data.time_coverage_start = "2000-04-01 00:00:00"
    nc_data.time_coverage_end = "2021-12-31 00:00:00"
    nc_data.geospatial_lat_min = str(ymin)
    nc_data.geospatial_lat_max = str(ymax)
    nc_data.geospatial_lat_units = "degrees_north"
    nc_data.geospatial_lon_min = str(xmin)
    nc_data.geospatial_lon_max = str(xmax)
    nc_data.geospatial_lon_units = "degrees_east"
    nc_data.geospatial_resolution = sp_re
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
    lons[:] = np.arange(xmin, xmax + 0.001, sp_re).tolist()
    lats[:] = np.arange(ymin, ymax + 0.001, sp_re).tolist()
    print(lats)
    
    if in_time == "":
    # 默认时间范围是2002年4月至2021年12月
        date = []
        year = 2002
        while year <= 2021:
            month = 1
            while month <= 12:
                if year * 100 + month >= 200204:
                    year_month = str(year) + str(month)
                    date.append(datetime.datetime.strptime(year_month, "%Y%m"))
                month += 1  # 结束条件
            year += 1  # 结束条件
        # 把时间格式处理为num格式
        num_date = date2num(date, units=times.units, calendar=times.calendar)
    
    elif in_time == "2002-2021":
        date = []
        year = 2002
        while year <= 2021:
            _year = str(year)
            date.append(datetime.datetime.strptime(_year, "%Y"))
            year += 1   # 结束条件
        # 把时间处理成num格式
        num_date = date2num(date, units=times.units, calendar=times.calendar)

    # 给时间赋值
    times[:] = num_date[:]

    # 给每个变量设置属性
    for key, value in vars.items():
        # 自动添加属性
        box = key.split(" ")
        standard_name = ""
        for i in range(len(box) - 1):
            standard_name = standard_name + str(box[i]) + "_"
        standard_name = standard_name + str(box[-1])
        long_name = key
        value.standard_name = standard_name
        value.long_name = long_name
        value.units = iunits  # 默认单位是mm
        value.Units = iunits
        value.fill_value = -9999.0
        value.missing_value = -9999
        value[:, :, :] = dict_vars[key]  # 给变量赋值

    # 保存
    print("start to save file")
    nc_data.close()
    print("DONE!")

    return 0



