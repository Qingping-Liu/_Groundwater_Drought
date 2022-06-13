# encoding:utf-8

"""
this script used to resampling EVI to make sure the spatial resolution(0.25) is consistent with the GLDAS and GRACE.

goals:
    --use the average value composite method to regridded 0.25*0.25 spatial resolution of EVI.

input:
    -- MODIS_EVI_GLOBAL_PRE_PROCESS.nc

output:
    -- MODIS_GLOBAL_EVI_Resample_0.25.nc
"""

# import modules
import netCDF4 as nc
from netCDF4 import num2date, date2num
import numpy as np
import datetime 
import os

# 从本地文件导入函数
from MODIS_Assimilation import write_nc

# 定义重采样函数
def resample_modis(nc_path:str, spatial_resolution: float, chunk_size=12):
    """_summary_
    Args:
        nc_path (str): the file of modis
        spatial_resolution (float): _description_
        chunk_size: 默认分块读取
    """
    # 读取数据集
    dataset = nc.Dataset(nc_path, mode="r")
    # 读取相关变量
    lats = dataset.variables["lat"][:]
    lons = dataset.variables["lon"][:]
    
    unit = dataset.variables["time"].units
    calend = dataset.variables["time"].calendar
    times = dataset.variables["time"][:]
    # 把时间格式转换成datetime格式
    times = num2date(times, unit, calend)
    # 读取其空间分辨率
    sp_re = float(dataset.geospatial_resolution)
    
    n_skip= int(spatial_resolution / sp_re)
    
    # 读取其值
    shp = (len(times), len(lats), len(lons))        # 重采样前的行列数
    
    new_rows = int(shp[1]/n_skip)         # 重采样以后的行数
    new_colums = int(shp[2]/n_skip)       # 重采样以后的列数
    
    # 定义一个空的三维数组
    matrix = np.empty(shape=(len(times), new_rows, new_colums), dtype=np.int16)
    
    '''     # 开始重采样
    for t in range(chunk_size):
        print(t)
        for i in range(new_rows):
            for j in range(new_colums):
                step_x = 0
                value = 0   
                count = 0       # 计有效值的个数
                while step_x < n_skip:
                    step_y = 0
                    while step_y < n_skip:
                        if vars[t, i*n_skip + step_x, j*n_skip + step_y] is not np.ma.masked:
                            num = vars[t, i*n_skip + step_x, j*n_skip + step_y]
                            value = value + num
                            count += 1 
                        else:
                            pass
                        step_y += 1 
                    step_x += 1
                if value != 0:
                    matrix[t, i, j] = value / count
                else:
                    matrix[t, i, j] = -9999 
    '''
    
    # 重采样算法2
    start = 0
    chunk = [50,100,150,200,237]
    for end in chunk:
        vars = dataset.variables["CMG_0.05_Deg_Monthly_EVI"][start:end, :, :]
        for t in np.arange(end-start):
            print(times[start+t])
            for i in range(new_rows):
                for j in range(new_colums):     
                    arr = np.array(vars[t, i*n_skip:i*n_skip+5, j*n_skip:j*n_skip+5])
                    mx = np.ma.masked_values(arr, -9999)
                    value = mx.sum()
                    count = 25 - np.count_nonzero(mx.mask) 
                    if count != 0:
                        matrix[start+t, i, j] = value / count
                    else:
                        matrix[start+t, i, j] = -9999.0

        start = end
        
    # 重采样算法3

    
    
    
    # 定义一个空字典
    dict_vars = {}
    dict_vars["CMG 0.25 Deg Monthly EVI"] = matrix
    
    return dict_vars


if __name__ == "__main__":
    top_path = "D:\\GRACE_ini_datasets\\initial_INPUT_data\\_Global"
    filename = "MODIS_EVI_GLOBAL_PRE_PROCESS.nc"
    
    nc_path = os.path.join(top_path, filename)
    
    # 范围
    rect = [-179.875, 179.875, -89.875, 89.875]
    
    # 重采样的空间分辨率
    sp = 0.25
    dict_var = resample_modis(nc_path, 0.25)
    
    outname = "MODIS_GLOBAL_EVI_Resample_025.nc"
    # 存储
    write_nc(dict_var, rect, outname, sp_re=0.25, iunits="--")