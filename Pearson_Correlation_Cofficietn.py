# encoding: utf-8

"""
this script used to calculating the pearson corellation cofficient with gaven nc file.


"""

# import modules
from ast import keyword
import numpy as np
import netCDF4 as nc
from netCDF4 import date2num, num2date
import datetime
import os
import sys
from scipy import stats

from Draw import r2_rmse
from write_tif import write_tif
from EXTRACT_DATA import extract_data


def func(x, k, b):
    return k*x+b



def pearson_co_cofficient(x, y):
    """_summary_

    Args:
        x (np.array): time series of independent variable.
        y (np.array): time series of independent variable.
    """
    # 判断长度是否一致
    assert len(x) == len(y)
    # 计算相关系数和p
    rou, p = stats.pearsonr(x, y)
    return rou, p


# nc文件和txt文件的相关性分析
def pcoco(nc_path1, keyword1, nc_path2, keyword2, format=True, fre_period="normal"):
    """_summary_

    Args:
        path1 (str): _description_
        path2 (str): _description_
        format (bool, optional): True represent the data format is nc, while False represent the data format is txt
    """
    if format:
        # 读取数据集
        dataset1 = nc.Dataset(nc_path1, mode="r")
        dataset2 = nc.Dataset(nc_path2, mode="r")

        # 读取相关变量
        lats = dataset1.variables["lat"][:]
        lons = dataset2.variables["lon"][:]

        times = dataset1.variables["time"][:]
        unit = dataset1.variables["time"].units
        calen = dataset1.variables["time"].calendar

        times = num2date(times, units=unit, calendar=calen)
        # 把所有数据读入缓存
        vars_1 = np.array(dataset1.variables[keyword1][:, :, :], dtype=np.float32)
        vars_2 = np.array(dataset2.variables[keyword2][:, :, :], dtype=np.float32)

        # 新建一个相关系数字典, 用于计算植被相关系数
        co_matrix = {}
        p_matrix = {}

        if fre_period == "normal":
            pcc = ["pearson_correlation_cofficient"]
            pks = ["p_month"]

            for i in range(len(pcc)):
                co_matrix[pcc[i]] = np.empty(shape=(len(lats), len(lons)))
                p_matrix[pks[i]] = np.empty(shape=(len(lats), len(lons)))

            # 处理缺失值
            vars_1[vars_1 == -9999.0] = np.nan
            vars_2[vars_2 == -9999.0] = np.nan
            vars1_dic = vars_1.copy()
            vars2_dic = vars_2.copy()

        elif fre_period == "monthly":
            pcc = ["Jan", "Feb", "Mar", "Ari", "May", "Jun", "Jul", "Aut", "Sep", "Oct", "Nov", "Dec"]
            pks = ["p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10", "p11", "p12"]

            for i in range(len(pcc)):
                co_matrix[pcc[i]] = np.empty(shape=(len(lats), len(lons)))
                p_matrix[pks[i]] = np.empty(shape=(len(lats), len(lons)))

            # 按季节提取数据
            var1_dic = extract_data(nc_path1, keyword1, "monthly")
            var2_dic = extract_data(nc_path2, keyword2, "monthly")

        # 按季节提取数据
        elif fre_period == "seasonly":
            pcc = ["spring", "summer", "autumn", "winter"]
            pks = ["p_spring", "p_summer", "p_autumn", "p_winter"]

            for i in range(len(pcc)):
                co_matrix[pcc[i]] = np.empty(shape=(len(lats), len(lons)))
                p_matrix[pks[i]] = np.empty(shape=(len(lats), len(lons)))

            # 按季节提取数据
            var1_dic = extract_data(nc_path1, keyword1, "seasonly")[0]
            var2_dic = extract_data(nc_path2, keyword2, "seasonly")[0]

        elif fre_period == "grow_period":
            pcc = ["grow_period"]
            pks = ["p_grow"]

            for i in range(len(pcc)):
                co_matrix[pcc[i]] = np.empty(shape=(len(lats), len(lons)))
                p_matrix[pks[i]] = np.empty(shape=(len(lats), len(lons)))

            # 按季节提取数据
            var1_dic = extract_data(nc_path1, keyword1, "grow_period")[0]
            var2_dic = extract_data(nc_path2, keyword2, "grow_period")[0]


        # 按年内平均计算相关系数和
        elif fre_period == "yearly":
            pcc = ["year"]
            pks = ["p_year"]

            for i in range(len(pcc)):
                co_matrix[pcc[i]] = np.empty(shape=(len(lats), len(lons)))
                p_matrix[pks[i]] = np.empty(shape=(len(lats), len(lons)))
            
            # 按季节提取数据
            var1_dic = extract_data(nc_path1, keyword1, "yearly")[0]
            var2_dic = extract_data(nc_path2, keyword2, "yearly")[0]

        del vars_1, vars_2

        for k in range(len(pcc)):
            vars_1 = var1_dic[pcc[k]][:, :, :]
            vars_2 = var2_dic[pcc[k]][:, :, :]
            for i in range(len(lats)):
                for j in range(len(lons)):
                    # 判断是否是缺失值
                    flag1 = np.count_nonzero(np.isnan(vars_1[:, i, j]))
                    flag2 = np.count_nonzero(np.isnan(vars_2[:, i, j]))
                    if not flag1 or not flag2:
                        value_x = vars_1[:, i, j]
                        value_y = vars_2[:, i, j] / 10000

                        #分母有可能为0
                        rou, p_value = pearson_co_cofficient(value_x, value_y)
                        co_matrix[pcc[k]][i, j] = rou
                        p_matrix[pks[k]][i, j] = p_value
                    else:
                        co_matrix[pcc[k]][i, j] = -9999.0
                        p_matrix[pks[k]][i, j] = -9999.0

        # 数据范围
        xmin, xmax = float(dataset2.variables["lon"].valid_min), float(dataset2.variables["lon"].valid_max)
        ymin, ymax = float(dataset1.variables["lat"].valid_min), float(dataset1.variables["lat"].valid_max)

        rect = [xmin, xmax, ymin, ymax]

        return co_matrix, p_matrix, rect
        
# 未准备好
def k_value(nc_path1, keyword1, nc_path2, keyword2, format=True, fre_period="normal"):
    """
    :param nc_path:
    :param keyword:
    :param format:
    :param fre_period:
    :return:
    """
    if format:
        # 读取数据集
        dataset1 = nc.Dataset(nc_path1, mode="r")
        dataset2 = nc.Dataset(nc_path2, mode="r")

        # 读取相关变量
        lats = dataset1.variables["lat"][:]
        lons = dataset2.variables["lon"][:]

        times = dataset1.variables["time"][:]
        unit = dataset1.variables["time"].units
        calen = dataset1.variables["time"].calendar

        times = num2date(times, units=unit, calendar=calen)

        # 把所有数据读入缓存
        vars_1 = np.array(dataset1.variables[keyword1][:, :, :], dtype=np.float32)
        vars_2 = np.array(dataset2.variables[keyword2][:, :, :], dtype=np.float32)

        co_matrix = {}
        # 用于分析非植被数据
        if fre_period == "R2_RMSE":
            from scipy import optimize

            # 把所有数据读入缓存
            vars_1 = vars_1[:117, :, :]
            vars_2 = vars_2[:117, :, :]
            print(times[117].strftime("%Y-%m"))

            keys = ["R2", "RMSE"]
            for key in keys:
                co_matrix[key] = np.empty(shape=(len(lats), len(lons)))

            for i in range(len(lats)):
                for j in range(len(lons)):
                    # 判断是否是缺失值
                    if vars_1[0, i, j] != -9999.0 or vars_2[0, i, j] != -9999.0:
                        value_x = vars_1[:, i, j]
                        value_y = vars_2[:, i, j]
                        k, b = optimize.curve_fit(func, value_x, value_y)[0]
                        co_matrix["R2"][i, j], co_matrix["RMSE"][i, j] = r2_rmse(value_x, value_y, k, b)
                    else:
                        co_matrix["R2"][i, j], co_matrix["RMSE"][i, j] = -9999.0, -9999.0

    return 0

if __name__ == "__main__":

    # 文件所在目录
    top_path = r"D:\GRACE_ini_datasets\Final_data"

    dsi_name = "NCP_GWSA-DSI_v2.nc"
    nc_path1 = os.path.join(top_path, dsi_name)
    keyword1 = "groundwater storage anomalies DSI"

    evi_name = "NCP_Terrestial_Monthly_EVI.nc"
    nc_path2 = os.path.join(top_path, evi_name)
    keyword2 = "CMG 0.25 Deg Monthly EVI"

    # 计算相关系数
    co_matrix, p_matrix, rect = pcoco(nc_path1, keyword1, nc_path2, keyword2, fre_period="yearly")

    outname1 = "DSI-EVI_Pearson_correlation_cofficient.tif"
    outname2 = "signigicant.tif"

    # 保存成nc文件
    #write_tif(dic_var=co_matrix, rect=rect, cell=0.25, outname=outname1)
    write_tif(dic_var=p_matrix, rect=rect, cell=0.25, outname=outname2)


    # # 文件所在目录
    # top_path = "C:\\Users\\11072\\Desktop\\python_code"
    #
    # dsi_name = "The Yearly changes of GWSA_DSI in the NCP from 2002 to 2021.nc"
    # nc_path1 = os.path.join(top_path, dsi_name)
    # keyword1 = "year"
    #
    # evi_name = "The Yearly changes of EVI in the NCP from 2002 to 2021.nc"
    # nc_path2 = os.path.join(top_path, evi_name)
    # keyword2 = "year"
    #
    # # 计算相关系数
    # co_matrix, p_matrix, rect = pcoco(nc_path1, keyword1, nc_path2, keyword2, fre_period="normal")
    #
    # outname = "Yearly_DSI-EVI_Pearson_correlation_cofficient.tif"
    # # 保存成nc文件
    # write_tif(dic_var=co_matrix, rect=rect, cell=0.25, outname=outname)
    # write_tif(dic_var=p_matrix, rect=rect, cell=0.25, outname=outname)



    # # GRACE 数据验证
    # top_path = "C:\\Users\\11072\\Desktop\\python_code\\data_pre_process"
    # filename1 = "NCP_GLDAS_Terrestial_Storage_Anomalies.nc"
    # filepath1 = os.path.join(top_path, filename1)
    # keyword1 = "terrestrial water storage anamolies"
    #
    # filename2 = "NCP_GRACE_Terrestial_Storage_Anomalies.nc"
    # filepath2 = os.path.join(top_path, filename2)
    # keyword2 = "lwe_thickness"
    #
    # # 计算相关系数
    # dic_var, rect = pcoco(filepath1, keyword1, filepath2, keyword2, fre_period="R2_RMSE")
    #
    # outname = "monthly pearson correlation between GLDAS and GRACE.tif"
    # write_tif(dic_var=dic_var, rect=rect, cell=0.25, outname=outname)
    