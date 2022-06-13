# encoding:utf-8


# import modules
import netCDF4 as nc
from netCDF4 import num2date
import numpy as np
import os
from collections import OrderedDict

# 导入计算干旱指标的函数
from dsi import drought_index

# 导入计算统计指标的函数
from statistic_index import statistic_index

# 导入提取数据函数
from Extract_data import extract_time_series

# 导入保存为txt文档
from data_pre_process.write_txt import write_txt


# 计算不同时间尺度的 dsi 的时间序列  ---->应先参与运算再归一化
# bootstrap重采样300次，采样率90%
def compute_timeseries_dsi(nc_path, keyword, fre_period="monthly", bootstrap=False, bootNum=30, bootPercent=0.75):
    
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
    dic_var = OrderedDict()
    
    # 采用bootstrap采样技术
    if bootstrap:
        # 时间键值对 先创建是为了后续写入txt有顺序
        keys = ["Date", "Average", "Std", "Min", "Max"]

        # 创建键值对
        for key in keys:
            if key != "Date":
                dic_var[key] = np.zeros(len(times))
            elif key == "Date":
                dic_var["Date"] = times[:]

      # 开始计算dsi
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
                    dic_var[pcc[i]] = np.zeros(len(times))
                    
            # 下一步，要排除年内周期性变化的影响
            temp_dic = {}
            temp_dsi ={}
            month = [i.month for i in times]
            # 按月份提取重复的次数
            unique, tm = np.unique(month, return_counts=True)
            for key in unique:
                temp_dic[str(key)] = np.zeros(tm[np.argwhere(unique==key)])
                temp_dsi[str(key)] = np.zeros(tm[np.argwhere(unique==key)])
            
            # 存放按时间序列排放的dsi结果
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
                        teplist.append(value)       # 按时间顺序排列的地下水储量
                        
                # 按月份提取数据
                for i in range(shp[0]):
                    mon = str(times[i].month)
                    j = i // 12  # 向下取整
                    temp_dic[mon][j] = teplist[i]

                # 计算每个月的dsi
                for mon, value in temp_dic.items():
                    temp_dsi[mon] = drought_index(value)
                    # 按时间顺序还原
                    for i in range(value.shape[0]):
                        if mon == "1" or mon == "2" or mon == "3":
                            j = int(i * 12 + 8 + int(mon))
                            dic_var[key][j] = temp_dsi[mon][i]
                        else:
                            j = int(i*12 + int(mon) - 4)
                            dic_var[key][j] = temp_dsi[mon][i]

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

            dic_var["Average"] = np.asarray(Average)
            dic_var["Std"] = np.asarray(Std)
            dic_var["MIn"] = np.asarray(Max)
            dic_var["Max"] = np.asarray(Min)



            # 不采用bootstrap技术
    else:
        pass


    return dic_var



if __name__== "__main__":
    
    
    # 地下水储量干旱指数dsi的时间序列
    top_path = "D:\\GRACE_ini_datasets\\mediate_data"
    filename1 = "NCP_Terrestial_Monthly_Groundwater_Storage_Anomalies_v2.nc"
    keyword = "groundwater storage anomalies"
    
    nc_path1 = os.path.join(top_path, filename1)

    # 计算dsi
    dic_vars= compute_timeseries_dsi(nc_path1, keyword, fre_period="monthly", bootstrap=True)
    
    outname = "The GWSA_DSI changes in the NCP from 2002 to 2021_v4.txt"
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