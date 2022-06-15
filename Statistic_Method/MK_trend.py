# encoding:utf-8

# from __future__ import division
import numpy as np
from scipy.stats import norm
import netCDF4 as nc
from netCDF4 import num2date
import os

from write_tif import write_tif


# mk检验
def mk_test(x, alpha=0.05):
    """
    This function is derived from code originally posted by Sat Kumar Tomer
    (satkumartomer@gmail.com)
    See also: http://vsp.pnnl.gov/help/Vsample/Design_Trend_Mann_Kendall.htm

    The purpose of the Mann-Kendall (MK) test (Mann 1945, Kendall 1975, Gilbert
    1987) is to statistically assess if there is a monotonic upward or downward
    trend of the variable of interest over time. A monotonic upward (downward)
    trend means that the variable consistently increases (decreases) through
    time, but the trend may or may not be linear. The MK test can be used in
    place of a parametric linear regression analysis, which can be used to test
    if the slope of the estimated linear regression line is different from
    zero. The regression analysis requires that the residuals from the fitted
    regression line be normally distributed; an assumption not required by the
    MK test, that is, the MK test is a non-parametric (distribution-free) test.
    Hirsch, Slack and Smith (1982, page 107) indicate that the MK test is best
    viewed as an exploratory analysis and is most appropriately used to
    identify stations where changes are significant or of large magnitude and
    to quantify these findings.

    Input:
        x:   a vector of data
        alpha: significance level (0.05 default)

    Output:
        trend: tells the trend (increasing, decreasing or no trend)
        h: True (if trend is present) or False (if trend is absence)
        p: p value of the significance test
        z: normalized test statistics

    Examples
    --------

    """
    if x[0] != -9999.0:
        # 获取list的长度
        n = len(x)

        # calculate S
        s = 0
        for k in range(n - 1):
            for j in range(k + 1, n):
                s += np.sign(x[j] - x[k])

        # calculate the unique data
        unique_x, tp = np.unique(x, return_counts=True)  # 返回去重后的元素，并且返回去重元素的数量
        g = len(unique_x)

        # calculate the var(s)
        if n == g:  # there is no tie
            var_s = (n * (n - 1) * (2 * n + 5)) / 18  # 如果每个元素唯一，则计算直接计算统计量s方差
        else:  # there are some ties in data    # 如果有重复的元素
            var_s = (n * (n - 1) * (2 * n + 5) - np.sum(tp * (tp - 1) * (2 * tp + 5))) / 18

        # 计算z值
        if s > 0:
            z = (s - 1) / np.sqrt(var_s)
        elif s < 0:
            z = (s + 1) / np.sqrt(var_s)
        else:  # s == 0:
            z = 0

        # calculate the p_value
        p = 2 * (1 - norm.cdf(abs(z)))  # two tail test 两端检验
        if abs(z) > norm.ppf(1 - alpha / 2):
            h = 1  # 1 表示接受假设
        else:
            h = 0  # 0 表示拒绝假设

        if (z < 0) and h:  # 接受
            trend = -1  # -1 表示下降趋势
        elif (z > 0) and h:  # 接受
            trend = 1  # 1 表示上升趋势
        else:
            trend = 0  # 0 表示无趋势

    else:
        h = -9999.0
        trend = -9999.0
        p = -9999.0
        z = -9999.0

    return z, h


# 趋势检验
def trend_test(nc_path, keyword, fre_period="normal", region_test=True):
    """
    :param nc_path:
    :param keyword:
    :param fre_period:
    :param format:
    :return: 返回的趋势性检验的结果
    """
    from EXTRACT_DATA import extract_data

    if not region_test:
        # 读取数据集
        dataset = nc.Dataset(nc_path, mode="r")

        # 读取相关变量
        lats = dataset.variables["lat"][:]
        lons = dataset.variables["lon"][:]

        times = dataset.variables["time"][:]
        unit = dataset.variables["time"].units
        calen = dataset.variables["time"].calendar

        times = num2date(times, units=unit, calendar=calen)
        # 把所有数据读入缓存
        vars = np.array(dataset.variables[keyword][:, :, :], dtype=np.float32)

        # 新建一个相关系数字典, 用于计算植被相关系数
        test_matrix = {}

        if fre_period == "normal":
            pcc = ["normal"]
            keys = ["z", "h"]

            for i in range(len(pcc)):
                for j in range(len(keys)):
                    key = pcc[i] + "_" + keys[j]
                    test_matrix[key] = np.empty(shape=(len(lats), len(lons)))

            # 处理缺失值
            var_dic = {}
            vars[vars == -9999.0] = np.nan
            var_dic["normal"] = vars.copy()


        # 按季节提取数据
        elif fre_period == "seasonly":
            pcc = ["spring", "summer", "autumn", "winter"]
            keys = ["z", "h"]

            for i in range(len(pcc)):
                for j in range(len(keys)):
                    key = pcc[i] + "_" + keys[j]
                    test_matrix[key] = np.empty(shape=(len(lats), len(lons)))

            # 按季节提取数据
            var_dic = extract_data(nc_path, keyword, "seasonly")[0]



        elif fre_period == "grow_period":
            pcc = ["grow_period"]
            keys = ["z", "h"]
            for i in range(len(pcc)):
                for j in range(len(keys)):
                    key = pcc[i] + "_" + keys[j]
                    test_matrix[key] = np.empty(shape=(len(lats), len(lons)))

            # 按季节提取数据
            var_dic = extract_data(nc_path, keyword, "grow_period")[0]



        # 按年内平均计算相关系数和
        elif fre_period == "yearly":
            pcc = ["year"]
            keys = ["z", "h"]

            for i in range(len(pcc)):
                for j in range(len(keys)):
                    key = pcc[i] + "_" + keys[j]
                    test_matrix[key] = np.empty(shape=(len(lats), len(lons)))

            # 按季节提取数据
            var_dic = extract_data(nc_path, keyword, "yearly")[0]

        del vars

        for k in range(len(pcc)):
            vars = var_dic[pcc[k]][:, :, :]
            for i in range(len(lats)):
                for j in range(len(lons)):
                    flag = np.count_nonzero(np.isnan(vars[:, i, j]))
                    if not flag:
                        value_x = vars[:, i, j]
                        # 分母有可能为0
                        z, h = mk_test(value_x)
                        test_matrix[pcc[k] + "_z"][i, j] = z
                        test_matrix[pcc[k] + "_h"][i, j] = h
                    else:
                        test_matrix[pcc[k] + "_z"][i, j] = -9999.0
                        test_matrix[pcc[k] + "_h"][i, j] = -9999.0

        # 数据范围
        xmin, xmax = float(dataset.variables["lon"].valid_min), float(dataset.variables["lon"].valid_max)
        ymin, ymax = float(dataset.variables["lat"].valid_min), float(dataset.variables["lat"].valid_max)

        rect = [xmin, xmax, ymin, ymax]

        return test_matrix, rect


# 热力图检验
def region_trend(top_path, keyword="DSI"):
    from EXTRACT_DATA import extract_data

    walk = os.walk(top_path)

    # 空数组
    dic_z = {}
    dic_h = {}
    pcc = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "spring", "summer", "autumn", "winter", "year"]
    city = ["BEIJING", "TIANJIN", "HENAN", "ANHUI", "JIANGSU", "SHANDONG", "SHANXI", "HEBEI"]
    for i in city:
        dic_z[i] = []
        dic_h[i] = []

    for root, dirs, filenames in walk:
        for name in filenames:
            ncpath = os.path.join(root, name)
            period = ["monthly", "seasonly", "yearly"]
            # 获取处理的名字
            region = name.split(".")[0]

            for p in period:
                dic = extract_data(ncpath, keyword, p)[0]
                # 保证按顺序读取并存放
                if p == "monthly":
                    freq = pcc[:12]
                elif p == "seasonly":
                    freq = pcc[12:16]
                elif p == "yearly":
                    freq = [pcc[-1]]
                for key in freq:
                    value = dic[key][:, :, :]
                    shp = value.shape[0]
                    sequence = []
                    for i in range(shp):
                        data = value[i, :, :].copy().flatten()
                        excepnan = data[~np.isnan(data)]
                        if excepnan.size !=0:
                            mean = excepnan.mean()
                            sequence.append(mean)
                        else:
                            sequence = [-9999.0]
                            break
                    # 计算z和h
                    z, h = mk_test(sequence)
                    dic_z[region].append(z)
                    dic_h[region].append(h)

                    del sequence

    return dic_z, dic_h

def main():
    # # 文件所在目录
    # top_path = r"D:\GRACE_ini_datasets\Final_data"
    #
    # dsi_name = "NCP_GWSA-DSI_v2.nc"
    # nc_path1 = os.path.join(top_path, dsi_name)
    # keyword1 = "groundwater storage anomalies DSI"
    #
    # # 计算相关系数
    # dic_matrix, rect = trend_test(nc_path1, keyword1, fre_period="normal")
    # outname1 = "DSI-Trend_test.tif"
    # # 保存成nc文件
    # write_tif(dic_var=dic_matrix, rect=rect, cell=0.25, outname=outname1, folder=r"D:\GRACE_ini_datasets\Final_data\linear\trend_tif")

    # # 文件所在目录
    # top_path = r"D:\GRACE_ini_datasets\Final_data"
    # evi_name = "NCP_Terrestial_Monthly_EVI.nc"
    # nc_path2 = os.path.join(top_path, evi_name)
    # keyword2 = "CMG 0.25 Deg Monthly EVI"
    # # # mk趋势检验
    # dic_matrix, rect = trend_test(nc_path2, keyword2, fre_period="yearly")
    # outname1 = "EVI-Trend_test.tif"
    #  # 保存成nc文件
    # write_tif(dic_var=dic_matrix, rect=rect, cell=0.25, outname=outname1, folder=r"D:\GRACE_ini_datasets\Final_data\linear\trend_tif")

    # # DSI趋势检验
    # top_path = "D:\\GRACE_ini_datasets\\Final_data"
    # dsi_name = "NCP_GWSA-DSI_v2.nc"
    # nc_path = os.path.join(top_path, dsi_name)
    # keyword = "groundwater storage anomalies DSI"
    # dic_matrix, rect = trend_test(nc_path, keyword, fre_period="grow_period")
    # outname = "DSI-Trend_test.tif"
    # write_tif(dic_var=dic_matrix, rect=rect, cell=0.25, outname=outname,
    #           folder=r"D:\GRACE_ini_datasets\Final_data\linear\trend_tif\DSI_trend")

    from write_txt import write_txt
    top_path = r"D:\GRACE_ini_datasets\Final_data\Regional_DSI"
    dic_z, dic_h = region_trend(top_path)
    write_txt(dic_z, "DSI_Trend_z.txt")
    write_txt(dic_h, "DSI_Trend_h.txt")

if __name__ == "__main__":
    main()
