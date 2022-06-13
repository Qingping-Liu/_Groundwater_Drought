# encoding:utf-8


# import modules
import netCDF4 as nc
from netCDF4 import num2date
import numpy as np
import os


# 导入计算干旱指标的函数
from dsi import drought_index

# 导入提取文件的数据
from Extract_data import extract_raster_series

# 导入保存为nc文件的函数
from write_nc import write_nc


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
        # 定义一个字典, 存放各个月份的结果
        dic_matrix = {}
        for i in range(len(unique)):     # 1,2,3表示Jan, Feb, March
            mon = str(unique[i])
            dic_matrix[mon] = np.empty(shape=(tp[i], len(lats), len(lons)))
        
        # 按月份提取数据
        dic_matrix, Flag = extract_raster_series(nc_path, keyword, fre_period="monthly")
        
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
        
        # 
        for mon, value in dic_dsi.items():
            for i in range(value.shape[0]):
                t = Flag[mon][i]
                DSI[t, :, :] = value[t, :, :]
                                
    dic_vars = {}
    dic_vars["groundwater storage anomalies DSI"] = DSI
    
    # 数据范围
    xmin, xmax = float(dataset.variables["lon"].valid_min), float(dataset.variables["lon"].valid_max)
    ymin, ymax = float(dataset.variables["lat"].valid_min), float(dataset.variables["lat"].valid_max)
    rect = [xmin, xmax, ymin, ymax]
    
    return dic_vars, rect



if __name__== "__main__":
    
    
    # 地下水储量干旱指数dsi的时间序列
    top_path = "D:\\GRACE_ini_datasets\\mediate_data"
    filename1 = "NCP_Terrestial_Monthly_Groundwater_Storage_Anomalies_v2.nc"
    keyword = "groundwater storage anomalies"
    
    nc_path1 = os.path.join(top_path, filename1)

    # 计算dsi
    dic_vars, rect = compute_raster_dsi(nc_path1, keyword)
    
    outname = "NCP_GWSA-DSI_v2.nc"
    # 存储
    write_nc(dic_vars, rect, outname, scope="North China Plain", sp_re=0.25, iunits=" ")
    
    
    
    
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