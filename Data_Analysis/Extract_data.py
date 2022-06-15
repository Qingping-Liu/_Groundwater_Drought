# encoding: utf-8

"""

this script used to extract data from nc file with providng frequence period:
    -- monthly
    -- seasonly
    -- yearly


return:
    -- a dictionary 
"""

# import modules
import numpy as np
import netCDF4 as nc
from netCDF4 import num2date
import os


# 根据年季月来提取数据
def extract_raster_series(nc_path, keyword, fre_period="monthly"):
    """_summary_

    Args:
        nc_path (_type_): the file of nc
        period (str, optional): 
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

    # 读取所有值
    vars = np.array(dataset.variables[keyword][:, :, :], dtype=np.float32)

    # 将缺失值替换为 -9999.0
    # vars[vars.mask == True] = np.nan
    # else:
    vars[vars[:, :, :] == -9999.0] = np.nan

    # 定义一个存储数据位置的字典
    dic_matrix = {}
    Flag = {}

    # 按月分别提取数据
    if fre_period == "monthly":
        month = [i.month for i in times]
        # 按月份提取重复的次数
        unique, tm = np.unique(month, return_counts=True)
        
        for i in range(len(unique)):  # 1,2,3表示Jan, Feb, March
            mon = str(unique[i])
            dic_matrix[mon] = np.empty(shape=(tm[i], len(lats), len(lons)))     # 二维矩阵
            Flag[mon] = np.zeros(tm[i])

        # 按月份提取数据
        for i in range(len(times)):
            mon = str(times[i].month)
            j = i // 12  # 向下取整
            dic_matrix[mon][j, :, :] = vars[i, :, :]
            Flag[mon][j] = i
            
        
        
        
            
            
    # 按季节分别提取数据
    elif fre_period == "seasonly":
        season = ["spring", "summer", "autumn", "winter"]
        # ！！！！！！
        tp = 20  # 一共20年

        # 定义一个字典
        dic_matrix = {}
        for i in range(len(season)):
            dic_matrix[season[i]] = np.zeros(shape=(tp, len(lats), len(lons)))

        i = 0
        while i < tp:
            if i == 0:
                j = 4
            else:
                j = 1
            while j <= 12:
                date = i * 12 + j - 4
                if j == 3 or j == 4 or j == 5:
                    dic_matrix["spring"][i, :, :] += vars[date, :, :]

                elif j == 6 or j == 7 or j == 8:
                    dic_matrix["summer"][i, :, :] += vars[date, :, :]

                elif j == 9 or j == 10 or j == 11:
                    dic_matrix["autumn"][i, :, :] += vars[date, :, :]

                elif j == 12 and i < 19:
                    dic_matrix["winter"][i, :, :] = vars[date, :, :] + vars[date + 1, :, :] + vars[date + 2, :, :]

                elif j == 12 and i == 19:
                    dic_matrix["winter"][i, :, :] = vars[date]
                j += 1

            for key, value in dic_matrix.items():
                if i == 0 and key == "spring":  # 第一年春季
                    dic_matrix[key][i, :, :] = value[i, :, :] / 2
                elif i == 19 and key == "winter":  # 21年冬季
                    dic_matrix[key][i, :, :] = value[i, :, :]
                else:
                    dic_matrix[key][i, :, :] = value[i, :, :] / 3
            i += 1

    # 按生长季节提取数据
    elif fre_period == "grow_period":
        grow = ["grow_period"]
        # ！！！！！！
        tp = 20  # 一共20年

        # 定义一个字典
        dic_matrix = {}
        for i in range(len(grow)):
            dic_matrix[grow[i]] = np.zeros(shape=(tp, len(lats), len(lons)))

        i = 0
        while i < tp:
            if i == 0:
                j = 4
            else:
                j = 1
            while j <= 12:
                date = i * 12 + j - 4
                if j>=4 and j<=10:
                    dic_matrix["grow_period"][i, :, :] += vars[date, :, :]
                j += 1

            for key, value in dic_matrix.items():
                dic_matrix[key][i, :, :] = value[i, :, :] / 7

            i += 1

    # 按年平均提取数据
    elif fre_period == "yearly":
        year = ["year"]
        tp = 20
        # 定义一个字典
        dic_matrix = {}
        for i in range(len(year)):
            dic_matrix[year[i]] = np.zeros(shape=(20, len(lats), len(lons)))

        i = 0
        while i < tp:
            # 定义一个空的数据
            arr = np.zeros(shape=(len(lats), len(lons)))
            if i == 0:
                j = 4
            else:
                j = 1
            while j <= 12:
                date = i * 12 + j - 4
                arr += vars[date, :, :]
                j += 1

            if i == 0:
                arr = arr / 9
            else:
                arr = arr / 12    
                    
            #缺失值
            arr[np.isnan(arr)] = -9999
            dic_matrix["year"][i, :, :] = arr
            i += 1
    
        # 数据范围
    xmin, xmax = float(dataset.variables["lon"].valid_min), float(dataset.variables["lon"].valid_max)
    ymin, ymax = float(dataset.variables["lat"].valid_min), float(dataset.variables["lat"].valid_max)

    rect = [xmin, xmax, ymin, ymax]
      

    return dic_matrix, Flag



# 提取时间序列
def extract_time_series(nc_path, keyword, fre_period="monthly"):
    """

    :param nc_path:
    :param keyword:
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

    # 读取变量
    vars = np.array(dataset.variables[keyword][:, :, :], dtype=np.float32)
    
    # 把缺失的变量用0填充
    vars[vars[:, :, :] == -9999.0] = 0
    
    # 把每一时刻面上所有的地下水储量变化量累加起来转换成时间序列
    lis_var = np.asarray([np.sum(vars[i, :, :]) for i in range(len(times))])
    
    # 定义一个字典
    dic_matrix = {}
    Flag = {}   # 记录数据的原始位置， 方便后续按时间顺序堆放
    
    if fre_period == "momthly":
        month = [i.month for i in times]
        # 按月份提取重复的次数
        unique, tm = np.unique(month, return_counts=True)
        
        for i in range(len(unique)):  # 1,2,3表示Jan, Feb, March
            mon = str(unique[i])
            dic_matrix[mon] = np.zeros(tm[i])
            Flag[mon] = np.zeros(tm[i])

        # 按月份提取数据
        for i in range(len(times)):
            mon = str(times[i].month)
            j = i // 12  # 向下取整
            dic_matrix[mon][j] = lis_var[i]
            Flag[mon][j] = i        # 记录数据的原始位置， 方便后续还原
    
    # 按季节分别提取数据
    elif fre_period == "seasonly":
        season = ["spring", "summer", "autumn", "winter"]
        # ！！！！！！
        tp = 20  # 一共20年

        # 定义一个空字典
        for key in season:
            dic_matrix[key] = np.zeros(tp)

        i = 0
        while i < tp:
            if i == 0:
                j = 4
            else:
                j = 1
            while j <= 12:
                date = i * 12 + j - 4
                if j == 3 or j == 4 or j == 5:
                    dic_matrix["spring"][i] += lis_var[date]

                elif j == 6 or j == 7 or j == 8:
                    dic_matrix["summer"][i] += lis_var[date]

                elif j == 9 or j == 10 or j == 11:
                    dic_matrix["autumn"][i] += lis_var[date]

                elif j == 12 and i < 19:
                    dic_matrix["winter"][i] = lis_var[date] + lis_var[date+1] + lis_var[date+2]

                elif j == 12 and i == 19:
                    dic_matrix["winter"][i] = lis_var[date]
                    
                    
                j += 1

            for key, value in dic_matrix.items():
                if i == 0 and key == "spring":  # 第一年春季
                    dic_matrix[key][i] = value[i] / 2
                elif i == 19 and key == "winter":  # 21年冬季
                    dic_matrix[key][i] = value[i]
                else:
                    dic_matrix[key][i] = value[i] / 3
            i += 1
            
    # 按生长季节提取数据
    elif fre_period == "grow_period":
        grow = ["grow_period"]
        # ！！！！！！
        tp = 20  # 一共20年
        
        # 定义一个字典
        for key in grow:
            dic_matrix[key] = np.zeros(tp)

        i = 0
        while i < tp:
            if i == 0:
                j = 4
            else:
                j = 1
            while j <= 12:
                date = i * 12 + j - 4
                if j>=4 and j<=10:
                    dic_matrix["grow_period"][i] += lis_var[date]
                j += 1

            for key, value in dic_matrix.items():
                dic_matrix[key][i] = value[i] / 7
            i += 1
        
        # 按年平均提取数据
    elif fre_period == "yearly":
        years = ["year"]
        tp = 20
        # 定义一个字典
        for key in years:
            dic_matrix[key] = np.zeros(tp)

        i = 0
        while i < tp:
            if i == 0:
                j = 4
            else:
                j = 1
            # 定义一个临时变量
            temp = 0
            while j <= 12:
                date = i * 12 + j - 4
                temp += lis_var[date]
                j += 1

            if i == 0:
                temp = temp / 9
            else:
                temp = temp / 12    
            dic_matrix["year"][i] = temp
            i += 1


    return dic_matrix, Flag