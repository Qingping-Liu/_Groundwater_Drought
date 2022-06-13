# encoding:utf-8

import netCDF4 as nc
from netCDF4 import date2num, num2date
import datetime
import numpy as np
from scipy import interpolate as inter
import os
from write_nc import write_nc




# 分段线性插值
def interpolate_method(total_time: list[float], obser_time: list[float], obser_value: np.array([float])):
    """

    :param totall_time:  including observation_time and missing_time
    :param obser_time:  merely have the time of observation
    :param obser_value:
    :return:
    """
    func = inter.interp1d(obser_time, obser_value, kind="linear")
    totall_value = func(total_time)

    return totall_value


def interpolate_GRACE(filepath, keyword):
    print("Start processing..")
    # 读取数据集
    data_nc = nc.Dataset(filepath, mode="r", format="NETCDF4")  # 读取nc文件，返回一个对象
    
    # 获取经度、维度
    lons = data_nc.variables["lon"][:]
    lats = data_nc.variables["lat"][:]
    
    # 获取时间戳
    unit = data_nc.variables["time"].Units
    # 获取日历的方式
    calen = data_nc.variables["time"].calendar
    obser_time = list(data_nc.variables["time"][:])  # nc格式的时间
    print("obser_time:", obser_time)
    # 获取缺失的时间
    missing_str = data_nc.months_missing.split(sep=";")
    missing_1 = []
    for i in missing_str:
        date = i.split(sep=",")
        if len(date) < 2:
            if i not in missing_1:      # 查重
                missing_1.append(date[0])
        else:
            for j in date:
                if j not in missing_1:  # 查重
                    missing_1.append(j)
    print("*"*20)
    print("date_missing: ", missing_1)
    missing = [datetime.datetime.strptime(i, "%Y-%m") for i in missing_1]  # 注意 "%Y-%m"格式

    # 把时间转换成nums
    print("datetime format:", missing)
    missing_time = date2num(missing, units="days since 2002-01-01T00:00:00Z", calendar=calen).tolist()  # 把它转换成list
    print("missing_time:", missing_time)
    total_time = obser_time + missing_time
    total_time.sort()  # 按时间先后排序


    vars = np.array(data_nc.variables[keyword])

    # 新建一个字典
    dict_var = {}
    dict_var[keyword] = np.empty(shape=(len(total_time), len(lats), len(lons)))

    for i in range(len(lats)):
        for j in range(len(lons)):
            obser_value = vars[:, i, j]
            total_value = interpolate_method(total_time, obser_time, obser_value)
            dict_var[keyword][:, i, j] = total_value



    return dict_var 

if __name__ == "__main__":
    top_path = "C:\\Users\\11072\\Desktop\\python_code\\py_dataExtraction"
    name = "CSR_GRACE_GRACE-FO_RL06_Mascons_all-corrections_v02.nc"
    
    filepath = os.path.join(top_path, name)
    keyword = "lwe_thickness"
    # 范围
    rect = [0.125, 359.875, -89.875, 89.875]
    
    out_name = "BNU_GRACE_GRACE-FO_RL06_Mascons_Linear_interpolation_v02.nc"
    
    
    
    
    dict_var = interpolate_GRACE(filepath, keyword)
    
    
    write_nc(dict_var, rect, out_name, sp_re=0.25, iunits="cm")
    
    