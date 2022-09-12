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
def extract_data(nc_path, keyword, fre_period="monthly"):
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

    if fre_period == "monthly":
        month = [i.month for i in times]
        # 按月份提取重复的次数
        unique, tm = np.unique(month, return_counts=True)

        # 定义一个字典
        dic_matrix = {}
        for i in range(len(unique)):  # 1,2,3表示Jan, Feb, March
            mon = str(unique[i])
            dic_matrix[mon] = np.empty(shape=(tm[i], len(lats), len(lons)))

        # 按月份提取数据
        for i in range(len(times)):
            mon = str(times[i].month)
            j = i // 12  # 向下取整
            dic_matrix[mon][j, :, :] = vars[i, :, :]


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
      

    return dic_matrix, rect


# 计算华北平原地区内多年平均值
def extrac_average_area(nc_path, keyword):
    from GWSA_DSI import drought_index
    """
    :param nc_path: dsi为归一化的结果，不能做简单累加再平均，故输入数据应该为NCP地下水储量变化量的数据
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
    vars[vars[:, :, :] == -9999.0] = np.nan

    # 创建一个字典
    dic_var = {}
    keys = ["multiyears_average"]
    for key in keys:
        dic_var[key] = np.empty(shape=(len(lats), len(lons)))


    temp = np.zeros(shape=(len(lats), len(lons)))
    for i in range(len(lats)):
        for j in range(len(lons)):
            temp[i, j] = np.sum(vars[:, i, j]) / len(times)

        print(temp[i, :])

    # 坦化成一维
    shp = temp.shape
    in_f = temp.flatten()
    Flag = np.isnan(in_f).tolist()
    in_arr = in_f[~np.isnan(in_f)]
    # 归一化
    out_dsi = drought_index(in_arr)

    # 还原成二维
    count = 0
    for i in range(shp[0]):
        for j in range(shp[1]):
            k = i * shp[1] + j
            if Flag[k] == True:
                temp[i, j] = -9999.0
                count += 1
            else:
                m = k-count
                temp[i, j] = out_dsi[m]
        print(temp[i, :])

    # 保存到字典
    dic_var[keys[0]] = temp.copy()

    # 数据范围
    xmin, xmax = float(dataset.variables["lon"].valid_min), float(dataset.variables["lon"].valid_max)
    ymin, ymax = float(dataset.variables["lat"].valid_min), float(dataset.variables["lat"].valid_max)

    rect = [xmin, xmax, ymin, ymax]



    return dic_var, rect


# 提取时间序列
def extract_time_series(nc_path, keyword):
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
    vars[vars[:, :, :] == -9999.0] = np.nan

    # 创建一个字典
    dic_var = {}
    keys = ["date", "average", "std", "_min", "_max"]
    for key in keys:
        if key == "date":
            dic_var[key] = [i for i in range(len(times))]
        else:
            dic_var[key] = np.empty(shape=(len(times)))

    for i in range(len(times)):
        li = []
        for j in range(len(lats)):
            for k in range(len(lons)):
                if np.isnan(vars[i, j, k]):
                    pass
                else:
                    li.append(vars[i, j, k])
        # 计算统计值
        average, std, min, max = statistic_index(li)
        dic_var["date"][i] = times[i]
        dic_var["average"][i] = average
        dic_var["std"][i] = std
        dic_var["_min"][i] = min
        dic_var["_max"][i] = max

    return dic_var


# 计算均值,标准差,最大值,最小值
def statistic_index(array):
    """
    :param array: 一维的时间序列数组,数据类型为float型
    :return:
    """
    # 计算均值

    array = np.asarray(array)
    n = len(array)  # 数组的长度
    _max = np.max(array)  # 计算最大值
    average = np.sum(array) / n  # 数据的均值
    std = pow(np.sum((array - average) ** 2) / n, 0.5)  # 计算标准差
    _min = np.min(array)  # 计算最小值

    return average, std, _min, _max


def main():
    from write_txt import write_txt
    from write_nc import write_nc
    from write_tif import write_tif

    """
    #计算dsi的时间数据
    # 测试，提取DSI的时间序列数据
    nc_path = "D:\\GRACE_ini_datasets\\Final_data"
    name = "NCP_GWSA-DSI_v2.nc"
    filepath = os.path.join(nc_path, name)
    keyword = "groundwater storage anomalies DSI"

    # 提取时间序列
    dict_var = extract_time_series(filepath, keyword)
    # ASCII码保存

    outname = "The GWSA_DSI changes in the NCP from 2002 to 2021.txt"
    write_txt(dict_var, outname)
    print("Process Done!")
    """



    """
    # 计算DSI的年内平均值
    nc_path = "D:\\GRACE_ini_datasets\\Final_data"
    name = "NCP_GWSA-DSI_v2.nc"
    filepath = os.path.join(nc_path, name)
    keyword = "groundwater storage anomalies DSI"
    dic_var, rect = extract_data(filepath, keyword, fre_period="yearly" )
    
    outname = "The Yearly changes of GWSA_DSI in the NCP from 2002 to 2021.nc"
    write_nc(dic_var, rect, outname, in_time="2002-2021" ,scope="North China Plain", sp_re=0.25)
    """

    # 计算多年DSI均值在空间上的多年平均值
    top_path = "D:\\GRACE_ini_datasets\\mediate_data"
    name = "NCP_Terrestial_Monthly_Groundwater_Storage_Anomalies_v2.nc"
    filepath = os.path.join(top_path, name)
    keyword = "groundwater storage anomalies"
    dic_var, rect = extrac_average_area(filepath, keyword)

    outname = "The spatial distribution of multi_year_average.tif"
    write_tif(dic_var, rect, 0.25, outname, folder="D:\\GRACE_ini_datasets\\Final_data\\spatial_distribution_of_DSI")

    
    # # 计算EVI的年内平均值
    # nc_path = "D:\\GRACE_ini_datasets\\Final_data"
    # name = "NCP_Terrestial_Monthly_EVI.nc"
    # filepath = os.path.join(nc_path, name)
    # keyword = "CMG 0.25 Deg Monthly EVI"
    # dic_var, rect = extract_data(filepath, keyword, fre_period="yearly" )
    #
    # outname = "The Yearly changes of EVI in the NCP from 2002 to 2021.nc"
    # write_nc(dic_var, rect, outname, in_time="2002-2021" ,scope="North China Plain", sp_re=0.25)
    
    
    
    """
    # 计算GLDAS的时间序列数据
    nc_path = "C:\\Users\\11072\\Desktop\\python_code\\data_pre_process"
    name = "NCP_GLDAS_Terrestial_Storage_Anomalies.nc"
    filepath = os.path.join(nc_path, name)
    keyword = "terrestrial water storage anamolies"
    
    # 提取时间序列
    dic_var = extract_time_series(filepath, keyword)
    # 保存成txt
    outname = "The_GLDAS_Terrestial_Storage_Anomalies in the NCP from 2002 to 2021.txt"
    write_txt(dic_var, outname)
    """
    
    """
    # 计算GRACE的时间序列数据
    nc_path = "C:\\Users\\11072\\Desktop\\python_code\\data_pre_process"
    name = "NCP_GRACE_Terrestial_Storage_Anomalies.nc"
    filepath = os.path.join(nc_path, name)
    keyword = "lwe_thickness"
    
    # 提取时间序列
    dic_var = extract_time_series(filepath, keyword)
    # 保存成txt
    outname = "The GRACE Terrestrial Storage Anomalies in the NCP from 2002 to 2021.txt"
    write_txt(dic_var, outname)
    """
    
    """
    # 提取evi的时间序列数据
    nc_path = "C:\\Users\\11072\\Desktop\\python_code\\data_pre_process"
    name = "NCP_Terrestial_Monthly_EVI.nc"
    filepath = os.path.join(nc_path, name)
    keyword = "CMG 0.25 Deg Monthly EVI"
    
    # 提取时间序列
    dic_var = extract_time_series(filepath, keyword)
    outname = "The changes of EVI in the NCP from 2002 to 2021.txt"
    write_txt(dic_var, outname)
    """

if __name__ == "__main__":
    main()
