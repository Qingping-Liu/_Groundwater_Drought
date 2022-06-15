# encoding:utf-8

"""
this script used to calculate the Groundwater Storage Anomalies - Drought Severity Index from April 2002 to December 2021 in the NCP.

goals:
    -- calculating the Drought Severity Index according to the following expression:
            average_GWSA(month) = sum(GWSA(month, year)) / n(month)     n = 1, 2, 3, 4 ...
            DSI(month, year) = GWSA(month, year) - average_GWSA(month)      month = 1,2,3,...,12

input:
    -- NCP_Terrestial_Monthly_Groundwater_Storage_Anomalies.nc

output:
    -- NCP_GWSA-DSI.nc
    -- The GWSA-DSI changes in the North China Plain from 2002 to 2021.txt
note:
    -- 归一化后的指数不能直接参与运算
    -- 计算不同时间尺度的dsi时，应先计算不同时间尺度的地下水储量变化量
"""

# import modules
import netCDF4 as nc
from netCDF4 import num2date, date2num
import numpy as np
import datetime
import os

# 引入随机模块
import random
from write_txt import write_txt
from write_nc import write_nc
from write_tif import write_tif

# 重力卫星干旱指数计算
def drought_index(array):
    """
    :param array: 一维的时间序列数组,数据类型为float型
    :return:
    """
    # 计算均值
    array = np.asarray(array)   # 可能是个list
    n = len(array)  # 数组的长度
    if array[0] != -9999.0:     # 判断是否是缺失值
        average = np.sum(array)/n   # 数据的均值
        std = pow(np.sum((array - average)**2)/n, 0.5)
        DSI = (array - average) / std
    else:
        DSI = array       # 以原来的缺失值填充
    return DSI


# 计算均值和标准差
def statistic_index(array):
    """
    :param array: 一维的时间序列数组,数据类型为float型
    :return:
    """
    # 计算均值

    array = np.asarray(array)
    n = len(array)  # 数组的长度
    _max = np.max(array)    # 计算最大值
    average = np.sum(array)/n   # 数据的均值
    std = pow(np.sum((array - average)**2)/n, 0.5)  # 计算标准差
    _min = np.min(array)        # 计算最小值

    return average, std, _min, _max



# 按月份提取数据并单独保存，用每一个像元的时间的时间序列来计算地下水干旱指数
# 计算不同时间尺度的 dsi 的时间序列  ---->应先参与运算再归一化
def compute_raster_dsi(nc_path, keyword, fre_period="monthly"):
    """_summary_

    Args:
        nc_path (str): the nc file of groundwater water anomalies
    """

    # 读取数据集
    dataset = nc.Dataset(nc_path, mode="r")
    # 读取变量
    lons = dataset.variables["lon"][:]
    lats = dataset.variables["lat"][:]
    # 读取时间
    times = dataset.variables["time"][:]
    unit = dataset.variables["time"].units
    calend = dataset.variables["time"].calendar
    # 将时间转换为datetime格式
    times = num2date(times, unit, calend)
    month = [i.month for i in times]
    # 按月份提取重复的次数
    unique, tp = np.unique(month, return_counts=True)
    
    if fre_period == "monthly":
        # 定义一个字典
        dic_matrix = {}
        for i in range(len(unique)):     # 1,2,3表示Jan, Feb, March
            mon = str(unique[i])
            dic_matrix[mon] = np.empty(shape=(tp[i], len(lats), len(lons)))
        
        # 按月份提取数据
        for i in range(len(times)):
            mon = str(times[i].month)
            j = i // 12     # 向下取整
            dic_matrix[mon][j, :, :] = dataset.variables[keyword][i, :, :]     
        
        # 定义一个存放dsi的三维数组
        dic_dsi = {}
        for i in range(len(unique)):
            mon = str(unique[i])
            dic_dsi[mon] = np.empty(shape=(tp[i], len(lats), len(lons)))
        
        # 根据栅格计算每一个像元的dsi
        for i in range(len(lats)):      # 纬度
            for j in range(len(lons)):  # 经度
                for mon, value in dic_matrix.items():
                    dic_dsi[mon][:, i, j] = drought_index(value[:, i, j])
        
        # 计算完成后，按时间顺序存放
        DSI = np.empty(shape=(len(times), len(lats), len(lons)))
        for i in range(1 + len(times)//12):
            for mon, value in dic_dsi.items():
                if int(mon) >= 1 and int(mon)<=3 and i < 19:
                    j = int(mon) + i*12 + 8    # 判断索引
                    DSI[j, :, :] = value[i, :, :]      # 从1开始
                elif int(mon) >3 :
                    j = int(mon) + (i-1)*12 + 8    # 判断索引
                    DSI[j, :, :] = value[i, :, :]  
                                
    dic_vars = {}
    dic_vars["groundwater storage anomalies DSI"] = DSI
    
    return dic_vars




# 计算不同时间尺度的 dsi 的时间序列  ---->应先参与运算再归一化
# bootstrap重采样300次，采样率90%
def compute_timeseries_dsi(nc_path, keyword, fre_period="monthly", bootNum=30, bootPercent=0.75):
    
    """
    :param nc_path: 输入地下水储量变化量nc文件
    :param fre_period: ["monthly", "season":["spring", "summer", "autumn", "winter"], "year"]
    :return:
    """
    
    # 读取数据集
    dataset = nc.Dataset(nc_path, mode="r")
    # 读取变量
    lons = dataset.variables["lon"][:]
    lats = dataset.variables["lat"][:]
    # 读取时间
    times = dataset.variables["time"][:]
    unit = dataset.variables["time"].units
    calend = dataset.variables["time"].calendar
    # 将时间转换为datetime格式
    times = num2date(times, unit, calend)
    
    vars = np.asarray(dataset.variables[keyword][:, :, :])
    shp = vars.shape
    # 缺失值替换
    vars[vars[:, :, :] == -9999.0] = np.nan

    # 新建一个空的字典
    dic_var = {}
    # 时间键值对 先创建是为了后续写入txt有顺序
    keys = ["date", "Average", "Std", "Min", "Max"]
    dic_var["date"] = times[:]

    if fre_period=="monthly":
        
        pcc = ["month"]
        # 重采样
        boot = []
        for i in range(30):
            bootname = "bootstrap_" + str(i+1)
            pcc.append(bootname)
            boot.append(bootname)
            
        # 新建一个空的字典
        for i in range(len(pcc)):
                dic_var[pcc[i]] = np.empty(shape=(1, len(times)))

        for key in pcc:
            teplist = []
            dsi = []
            if key == "month":
                for i in range(shp[0]):
                    traster = vars[i, :, :]
                    excep_nan = traster[~np.isnan(traster)].flatten()     # 除nan值像元
                    num = len(excep_nan)
                    value = np.sum(excep_nan)/num        # 计算面地下水储量变化量
                    teplist.append(value)       # 某一个月的所有栅格的累积地下水储量变化量
                    
            else:
                for i in range(shp[0]): 
                    traster = vars[i, :, :]
                    excep_nan = traster[~np.isnan(traster)].flatten()     # 除nan值像元
                    # 随机抽样
                    num = int(len(excep_nan) * bootPercent)    # 按总体75%随机抽样
                    np.random.shuffle(excep_nan)    # 随机打乱
                    bvalue = excep_nan[:num]  # 
                    value = np.sum(bvalue)/num
                    teplist.append(value)
                    
            # 计算每个月的dsi值
            dsi = drought_index(teplist)
            dic_var[key] = np.asarray(dsi)
        
        del teplist, dsi    
        
        # 额外功能，用来计算统计指标：均值，标准差，最小值， 最大值
        Average = []
        Std = []
        Min = []
        Max = []
        for i in range(len(times)):
            templist = []
            _ave = 0
            _std = 0
            _min = 0
            _max = 0
            for key, value in dic_var.items():
                if key == "month":
                    pass
                elif key in boot:
                    templist.append(value[i])
                
            _ave, _std, _min, _max = statistic_index(templist)
            
            Average.append(_ave)
            Std.append(_std)
            Min.append(_min)
            Max.append(_max)
    
    elif fre_period == "seasonly":
        
        
        

        dic_var["Average"] = np.asarray(Average)
        dic_var["Std"] = np.asarray(Std)
        dic_var["MIn"] = np.asarray(Max)
        dic_var["Max"] = np.asarray(Min)    
    
    
    return dic_var


# 
    
    
                    
        
        
    




if __name__== "__main__":
    
    
    # 地下水储量干旱指数dsi的时间序列
    top_path = "D:\\GRACE_ini_datasets\\mediate_data"
    filename1 = "NCP_Terrestial_Monthly_Groundwater_Storage_Anomalies_v2.nc"
    keyword = "groundwater storage anomalies"
    
    nc_path1 = os.path.join(top_path, filename1)

    # 计算dsi
    rect1 = [112.625, 122.375, 32.125, 40.375]        # NCP的范围
    dic_vars = compute_timeseries_dsi(nc_path1, keyword)
    
    outname = "The GWSA_DSI changes in the NCP from 2002 to 2021_v3.txt"
    # 存储
    write_txt(dic_vars, outname)
    
    
    
    
    '''     
    # 全球地下水水库相对变化量
    top_path = "D:\\GRACE_ini_datasets\\mediate_data\\linear"
    filename2= "Global_Groundwater_Reservoir_Monthly_Groundwater_Storage_Anomalies.nc"
    nc_path2 = os.path.join(top_path, filename2)
    outname2 = "GLOBAL_GWSA-DSI.nc"
    # 范围
    rect2 = [-124.875, 151.375, -33.875, 72.875]
    # 计算
    dic_vars2 = compute_dsi(nc_path2)
    # 存储
    write_nc(dic_vars2, rect2, outname2)
    '''