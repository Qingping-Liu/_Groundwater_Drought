# encoding: utf-8

"""
this script used to draw Line chart, Scatter chart, 2-D GIS figure.

"""

import datetime
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import collections  # 线段集合模块
from scipy import optimize

# BG分割算法
from B_G_segement import *



# 定义常量
# 共33个月数据缺失
START_INDEX = [2, 14, 105, 110, 121, 126, 131, 136, 143, 148, 153, 159, 163, 170, 174, 179, 201]
END_INDEX = [4, 15, 106, 111, 122, 127, 132, 138, 144, 149, 154, 160, 165, 171, 175, 198, 203]

# 干旱阈值
THRESHOLD = -0.8

# 定义一个线性的函数
def func(x, k, b):
    return k * x + b


# 计算决定系数R，RMSE
def r2_rmse(x, y, k, b):
    """_summary_

    Args:
        x (_type_): 1维array
        y (_type_): 1维array
        k (_type_): 斜率
        b (_type_): 截距

    Returns:
        _type_: R2, RMSE
    """
    assert len(x) == len(y)
    y_mean = np.mean(y)

    y_pre = x * k + b
    R2 = np.sum((y_pre - y_mean) ** 2) / np.sum((y - y_mean) ** 2)
    RMSE = pow(np.sum((y_pre - y) ** 2) / len(x), 0.5)

    return round(R2, 2), round(RMSE, 2)


# 游程理论
def fish(in_array):
    """
    :param in_array:
    :return:
    """
    start = []
    end = []
    left = 0
    while left < len(in_array) - 1:
        if in_array[left] < THRESHOLD:
            start.append(left)
            right = left
            while right < len(in_array) - 1:
                if in_array[right] < THRESHOLD:
                    right += 1
                else:
                    end.append(right)
                    left = right + 1
                    break
            else:
                left = right
                end.append(right)
        else:
            left = left + 1
    print(start)
    print(end)
    return start, end


# 绘制干旱指标DSI的折线图
def draw_line_pic(in_txt, outname, error_bar=True, show_missing=False):
    """_summary_

    Args:
        in_txt (_type_): _description_
        outname (_type_): _description_
        error_bar (bool, optional): _description_. Defaults to True.
    """
    date = []
    ave = []
    std = []
    # 将数据导入
    with open(in_txt, "r") as f:
        readlines = f.readlines()[4:]
        # 跳过前5行
        for read in readlines:
            read = read.split("\t")
            date.append(datetime.datetime.strptime(read[0], "%Y-%m"))
            ave.append(float(read[1]))
            std.append(float(read[2]))

    f.close()

    ave = np.array(ave)
    std = np.array(std)

    # 确定图框的大小
    y_max = round(np.max(ave + std), 0)
    y_min = round(np.min(ave - std), 0)
    # 统计出干旱的频率
    starts, ends = fish(ave)

    # B-G分割算法,识别时间序列突变点
    segements = BG_test(ave, 0.95)


    # 创建图框
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot()
    
    # 把插值的区域高亮显示出来
    if show_missing:
        miss_start = START_INDEX
        miss_end = END_INDEX
        for i in range(len(miss_start)):
            ax.fill_betweenx([y_min, y_max], date[miss_start[i]], date[miss_end[i]], color="gray", alpha=0.10)
    
    # 识别发生干旱的时间段，并高亮显示
    for i in range(len(starts)):
        ax.fill_betweenx([y_min, y_max], date[starts[i]], date[ends[i]], color="orange", alpha=0.75)
    
    # 折线图
    yerr = np.zeros([2, len(ave)])
    yerr[0, :] = std
    yerr[1, :] = std
    ax.errorbar(date, ave, yerr=yerr[:, :], fmt="none", ecolor="gray", elinewidth=1, capsize=2, capthick=1, alpha=0.5)
    ax.bar(date, ave, width=31, color="green")

    # B-G分割算法结果
    ax.plot(date, segements[0, :], lw=1.5)




    # 设置标题
    # ax.set_title('The GWSA-DSI changes in the North Plain from 2002 to 2021', fontsize=15)
    ax.set_xlabel('Year', fontsize=14)
    ax.set_ylabel('GWSA-DSI', fontsize=14)

    # 设置横轴为时间坐标轴
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

    # 坐标轴的范围
    time_1 = datetime.datetime.strptime("2002", "%Y")
    time_2 = datetime.datetime.strptime("2022", "%Y")
    ax.set_xlim((time_1, time_2))
    # y坐标轴的范围
    ax.set_ylim(y_min, y_max)
    # 设置坐标轴的刻度
    ax.set_xticks = [datetime.datetime.strptime(str(i), "%Y") for i in np.arange(2002, 2022, 2)]

    # ax.set_yticks = [-3, -2, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0]

    ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(.5))

    # 绘制横线标识出干旱的区域
    ycateg = [-0.8, -1.3, -1.6, -2.0]
    for i in range(len(ycateg)):
         ax.axhline(ycateg[i], lw=1, color="r", linestyle="--")


    ax.yaxis.grid(True, which='minor', linewidth=0.25, linestyle='-', color='orange')
    ax.yaxis.grid(True, which='major', linewidth=0.75, linestyle='-', color='orange')

    plt.savefig(outname, dpi=300)
    plt.show()

    return 0


# 验证陆面水储量数据的折线图
def drwa_plot(in_txt1, in_txt2, outname, error_bar=True, show_missing=True):
    date = []
    ave1 = []
    std1 = []
    # 将数据导入
    with open(in_txt1, "r") as f:
        readlines = f.readlines()[4:]
        # 跳过前5行
        for read in readlines:
            read = read.split("\t")
            # print(read)
            date.append(datetime.datetime.strptime(read[0], "%Y-%m"))
            ave1.append(float(read[1]))
            std1.append(float(read[2]))

        f.close()

    ave2 = []
    std2 = []
    with open(in_txt2, "r") as f:
        readlines = f.readlines()[4:]
        # 跳过前5行
        for read in readlines:
            read = read.split("\t")

            ave2.append(float(read[1]))
            std2.append(float(read[2]))

        f.close()

    ave = np.array([ave1, ave2])
    std = np.array([std1, std2])

    # 确定图框的大小
    y_max = round(max(np.max(ave[0] + std[0]), np.max(ave[1] + std[1])), 0)
    y_min = round(min(np.min(ave[0] - std[0]), np.min(ave[1] - std[1])), 0)

    # 创建图框
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot()

    # 折线图
    yerr_1 = np.zeros([2, len(date)])
    yerr_1[0, :] = std[0]
    yerr_1[1, :] = std[0]

    yerr_2 = np.zeros([2, len(date)])
    yerr_2[0, :] = std[1]
    yerr_2[1, :] = std[1]

    ax.plot(date, ave[0], color="g", marker="v", markersize=3, lw=1)
    ax.plot(date, ave[1], color="orange", marker="s", markersize=3, lw=1)
    if error_bar:
        ax.errorbar(date, ave[0], yerr=yerr_1[:, :], fmt="none", ecolor="g", elinewidth=1, capsize=2, capthick=1,
                    alpha=0.5)
        ax.errorbar(date, ave[1], yerr=yerr_2[:, :], fmt="none", ecolor="orange", elinewidth=1, capsize=2, capthick=1,
                    alpha=0.5)

    # 把插值的区域高亮显示出来
    if show_missing:
        miss_start = START_INDEX
        miss_end = END_INDEX
        for i in range(len(miss_start)):
            ax.fill_betweenx([y_min, y_max], date[miss_start[i]], date[miss_end[i]], color="gray", alpha=0.20)

    # 设置标题
    ax.set_title('The TWSA of GRACE and GLDAS-NOAH2.1 in the NCP from 2002 to 2021', fontsize=15)
    ax.set_xlabel('Year', fontsize=14)
    ax.set_ylabel('TWSA [mm]', fontsize=14)

    # 设置横轴为时间坐标轴
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

    # 坐标轴的范围
    time_1 = datetime.datetime.strptime("2002", "%Y")
    time_2 = datetime.datetime.strptime("2021", "%Y")
    ax.set_xlim((time_1, time_2))
    # y坐标轴的范围
    ax.set_ylim(y_min, y_max)
    # 设置坐标轴的刻度
    ax.set_xticks = [datetime.datetime.strptime(str(i), "%Y") for i in np.arange(2002, 2022, 2)]

    plt.show()
    fig.savefig(outname, dpi=300)


# 数据验证的散点图
def draw_scatter(in_txt1, in_txt2, outname, mark_missing=True):
    date = []
    ave1 = []
    # 将数据导入
    with open(in_txt1, "r") as f:
        readlines = f.readlines()[4:]
        # 跳过前5行
        for read in readlines:
            read = read.split("\t")
            date.append(datetime.datetime.strptime(read[0], "%Y-%m"))
            ave1.append(float(read[1]))
        f.close()

    ave2 = []
    with open(in_txt2, "r") as f:
        readlines = f.readlines()[4:]
        # 跳过前5行
        for read in readlines:
            read = read.split("\t")
            ave2.append(float(read[1]))
        f.close()

    ave = np.array([ave1, ave2])

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot()

    x_group = ave[0, 0:117]
    y_group = ave[1, 0:117]

    print(date[117])
    # 散点图
    ax.scatter(ave[0, :], ave[1, :], s=30, color="green", alpha=0.75)
    ax.scatter(x_group, y_group, s=30, label="Observation Values(2002-2012)")

    # 线性拟合拟合的参数
    k, b = optimize.curve_fit(func, x_group, y_group)[0]
    R2, RMSE = r2_rmse(x_group, y_group, k, b)

    # 确定x的大小
    x_max = int(np.max(x_group))
    x_min = int(np.min(x_group))

    x = np.arange(x_min, x_max, 0.01)
    y = k * x + b
    # 拟合曲线
    ax.plot(x, y, lw=2, color="r", label="Linear Regression Curve")

    # 文本框
    ax.text(x=0.05, y=0.85, s="R2=" + str(R2) + "0", fontdict=dict(fontsize=10, color="black"), transform=ax.transAxes)
    ax.text(x=0.05, y=0.80, s="RMSE=" + str(RMSE) + "mm", fontdict=dict(fontsize=10, color="black"),
            transform=ax.transAxes)

    # 坐标轴名称
    ax.set_title("The correlation between DSI and GLDAS from 2002 to 2012", fontsize=14)
    ax.set_xlabel("TWSA-GRACE [mm]", fontsize=12)
    ax.set_ylabel("TWSA-GLDAS [mm]", fontsize=12)

    # 设置坐标轴的刻度
    ax.set_xticks = [datetime.datetime.strptime(str(i), "%Y") for i in np.arange(2002, 2022, 2)]

    plt.legend()
    # plt.show()
    plt.savefig(outname, dpi=300)

    return 0


# 绘制EVI和DSI的双坐标图，分析evi和dsi的相关性
def draw_EVI_DSI(in_txt1, in_txt2, outname, show_missing=True):
    from Pearson_Correlation_Cofficietn import pearson_co_cofficient

    date = []
    ave1 = []
    # 将数据导入
    with open(in_txt1, "r") as f:
        readlines = f.readlines()[4:]
        # 跳过前5行
        for read in readlines:
            read = read.split("\t")
            date.append(datetime.datetime.strptime(read[0], "%Y-%m"))
            ave1.append(float(read[1]))
        f.close()

    ave2 = []
    with open(in_txt2, "r") as f:
        readlines = f.readlines()[4:]
        # 跳过前5行
        for read in readlines:
            read = read.split("\t")
            ave2.append(float(read[1]) / 10000)  # ave2的为evi时，要除以尺度因子
        f.close()
    # 相关系数
    cov, p = pearson_co_cofficient(ave1, ave2)

    ave = np.array([ave1, ave2])

    x_group = [i for i in range(len(date))]
    y1_group = ave[0, :]
    y2_group = ave[1, :]
    # 线性拟合拟合的参数
    k1, b1 = optimize.curve_fit(func, x_group, y1_group)[0]
    k2, b2 = optimize.curve_fit(func, x_group, y2_group)[0]

    fig = plt.figure(figsize=(16, 8))
    ax1 = fig.add_subplot(111)

    p1, = ax1.plot(date, ave[0, :], lw=2, color="b")
    ax2 = ax1.twinx()
    p2, = ax2.plot(date, ave[1, :], lw=2, color="green")

    # 拟合曲线
    y1 = [k1 * float(i) + b1 for i in x_group]
    y2 = [k2 * float(i) + b2 for i in x_group]

    p3, = ax1.plot(date, y1, lw=2, linestyle=":", color="b")
    p4, = ax2.plot(date, y2, lw=2, linestyle=":", color="green")

    # 确定图框的大小
    y_max = 1.5
    y_min = -2.0

    # 把插值的区域高亮显示出来
    if show_missing:
        miss_start = np.asarray(START_INDEX) - 8
        miss_end = np.asarray(END_INDEX) - 8
        for i in range(len(miss_start)):
            ax1.fill_betweenx([y_min, y_max], date[miss_start[i]], date[miss_end[i]], color="gray", alpha=0.20)

    # 设置坐标轴范围
    ax1.set_xticks = [datetime.datetime.strptime(str(i), "%Y") for i in np.arange(2002, 2022, 2)]
    ax2.set_ylim(0, 0.5)
    # 设置坐标轴字体大小
    ax1.tick_params(axis="x", labelsize=14)
    ax2.tick_params(axis="y", labelsize=14)
    ax1.tick_params(axis="y", labelsize=14)

    # 设置坐标轴和图例
    # ax1.set_title("the relationship between GWSA-DSI and EVI", fontsize=24)
    ax1.set_xlabel("Year", fontsize=18)
    ax1.set_ylabel("GWSA-DSI", fontsize=18)
    ax2.set_ylabel("EVI", fontsize=18)

    plt.legend([p1, p2, p3, p4], ["DWSA-DSI", "EVI", "EVI trend", "GWSA-DSI trend"], loc=3, frameon=False)
    plt.savefig(outname, dpi=300)


# 绘制多年月变化情况
def draw_seasonly_cycle(nc_path1, keyword1, nc_path2, keyword2, outname):
    from EXTRACT_DATA import extract_data
    months = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

    # 植被
    dic_evi = extract_data(nc_path1, keyword1)[0]
    # dsi
    dic_dsi = extract_data(nc_path2, keyword2)[0]

    # 计算每月的均值
    y_values = np.empty((2, len(months)))
    estd = []
    dstd = []
    for key in months:
        evis = dic_evi[str(key)]
        dsis = dic_dsi[str(key)]
        evi_ave = evis[~np.isnan(evis)].mean() /10000  # 检查是否是乱序的
        evi_std = evis[~np.isnan(evis)].std() / 10000
        estd.append(evi_std)

        dsi_ave = dsis[~np.isnan(dsis)].mean()
        dsi_std = dsis[~np.isnan(evis)].std()
        dstd.append(dsi_std)
        y_values[0, months.index(key)] = evi_ave
        y_values[1, months.index(key)] = dsi_ave

    fig = plt.figure(figsize=(8, 8))
    ax1 = fig.add_subplot(111)
    # ax2 = ax1.twinx()

    # evi的误差棒
    evi_yerr = np.zeros([2, len(months)])
    evi_yerr[0, :] = estd[:]
    evi_yerr[1, :] = estd[:]
    ax1.errorbar(months, y_values[0, :], yerr=evi_yerr[:, :], fmt="none", ecolor="gray", elinewidth=1, capsize=2, capthick=1, alpha=0.5)

    # # dsi的误差棒
    # dsi_yerr = np.zeros([2, len(months)])
    # dsi_yerr[0, :] = dstd[:]
    # dsi_yerr[1, :] = dstd[:]
    # ax2.errorbar(months, y_values[1, :], yerr=dsi_yerr[:, :], fmt="none", ecolor="gray", elinewidth=1, capsize=2, capthick=1, alpha=0.5)

    ax1.plot(months, y_values[0, :], lw=2, color="g", label="EVI")
    # ax2.plot(months, y_values[1, :], lw=2, color="b")

    # 设置xy坐标
    ax1.set_xticks(months)
    ax1.set_xticklabels(["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"], fontsize=14)
    ax1.tick_params(axis="y", labelsize=16)

    ax1.set_xlabel("Month", fontsize=16)
    ax1.set_ylabel("EVI", fontsize=16)

    plt.legend(fontsize=16)
    plt.show()


# 绘制带标注的热力图
def anno_hotmap(in_txt1, in_txt2, outname):
    from matplotlib import colors

    matrix_z = np.empty(shape=(17, 8))

    with open(in_txt1, "r") as f:
        readlines = f.readlines()[4:]
        # 跳过前5行
        for i in range(len(readlines)):
            read = readlines[i].split("\t")[:-1]
            matrix_z[i] = np.array(read)
        f.close()

    matrix_h = np.empty(shape=(17, 8), dtype=int)
    with open(in_txt2, "r") as f:
        readlines = f.readlines()[4:]
        # 跳过前5行
        for i in range(len(readlines)):
            read = readlines[i].split("\t")[:-1]
            matrix_h[i] = np.array(read)
        f.close()

    # 矩阵转置
    matrix_h = matrix_h.T
    matrix_z = matrix_z.T

    # 无趋势提取出来
    for i in range(matrix_h.shape[0]):
        for j in range(matrix_h.shape[1]):
            if matrix_h[i, j] == 0:
                matrix_z[i, j] = 0  # 无趋势

    pcc = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Spr", "Sum", "Aut",
           "Win", "Year"]
    city = ["BEIJING", "TIANJIN", "HENAN", "ANHUI", "JIANGSU", "SHANDONG", "SHANXI", "HEBEI"]

    # 绘制热力图
    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot(111)

    # 设置colorbar的颜色
    color = ["darkorange", "gold", "lawngreen", "lightseagreen", "deepskyblue"]
    # cmap = colors.ListedColormap(["darkorange", "gold", "lawngreen", "lightseagreen", "deepskyblue"])
    cmap = colors.ListedColormap(color)
    # 根据z值的大小分配颜色
    z_discret = np.array([np.min(matrix_z), -2.58, -1.96, -1.65, -1.28, 0.1])
    norm = colors.BoundaryNorm(boundaries=z_discret, ncolors=5)

    # 绘制热力图
    im = ax.imshow(matrix_z, cmap=cmap, norm=norm)

    # colorbar
    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, cmap=cmap, orientation='vertical', shrink=0.70)

    # label的位置标注在中心
    clab=[np.min(matrix_z), 0]
    # for i in range(len(z_discret)-1):
    #     clab.append((z_discret[i] + z_discret[i+1])/2)
    cbar.ax.set_yticks(z_discret)
    cbar.ax.set_yticklabels(labels=["-5.80", "-2.58", "-1.96", "-1.65", "-1.28", "0.1"], fontdict={"family":"Times New Roman",'weight': 'normal',"size":12})
    cbar.ax.set_ylabel("Z Values of GWSA-DSI", rotation=-90, va="bottom",fontsize=16)


    # Loop over data dimensions and create text annotations.
    for i in range(len(city)):
        for j in range(len(pcc)):
            if matrix_z[i, j] < -2.58 or matrix_z[i, j] > 2.58:
                text = ax.text(j, i, "***", ha="center", va="center", color="w")
            elif -2.58 <= matrix_z[i, j] < -1.96 or 2.58 > matrix_z[i, j] >= 1.96:
                text = ax.text(j, i, "**", ha="center", va="center", color="w")
            elif -1.96 <= matrix_z[i, j] < -1.68 or 1.96 > matrix_z[i, j] >= 1.68:
                text = ax.text(j, i, "*", ha="center", va="center", color="w")
            else:
                pass



    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(pcc)))
    ax.set_xticklabels(pcc, fontdict={'family': 'Times New Roman','weight': 'normal',"size":12})
    ax.set_yticks(np.arange(len(city)))
    ax.set_yticklabels(city, fontdict={'family': 'Times New Roman','weight': 'bold', "size":12})

    # x坐标轴标签旋转90度
    plt.setp(ax.get_xticklabels(), rotation=0, ha="center",
             rotation_mode="anchor")

    # 关闭骨架并创建白色网格
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    # 设置白色网格
    ax.set_xticks(np.arange(len(pcc) + 1) - .5, minor=True)
    ax.set_yticks(np.arange(len(city) + 1) - .5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)

    fig.tight_layout()
    plt.show()


def main():
    # # DSI 的时间序列
    path = "D:\\GRACE_ini_datasets\\Final_data\\NCP_DSI_v3"
    name = "The GWSA_DSI changes in the NCP from 2002 to 2021_v4.txt"
    filepath = os.path.join(path, name)
    
    outname = "The GWSA_DSI changes in the NCP from 2002 to 2021_v4—2.jpg"
    
    draw_line_pic(filepath, outname, show_missing=True)

    # # 陆面水储量验证折线图
    # path = "C:\\Users\\11072\\Desktop\\python_code"
    # name1 = "The GRACE Terrestrial Storage Anomalies in the NCP from 2002 to 2021.txt"
    # filepath1 = os.path.join(path, name1)
    # name2 = "The_GLDAS_Terrestial_Storage_Anomalies in the NCP from 2002 to 2021.txt"
    # filepath2 = os.path.join(path, name2)
    #
    # outname = "Time series of TWSA of CSR and Noah in NCP from April 2002 to Dec 2021.jpg"
    #
    # drwa_plot(filepath1, filepath2, outname)

    # # DSI_EVI散点图
    # path = "D:\\GRACE_ini_datasets\\Final_data\\fig"
    # name1 = "The GRACE Terrestrial Storage Anomalies in the NCP from 2002 to 2021.txt"
    # filepath1 = os.path.join(path, name1)
    # name2 = "The_GLDAS_Terrestial_Storage_Anomalies in the NCP from 2002 to 2021.txt"
    # filepath2 = os.path.join(path, name2)
    #
    # outname = "The correlation between DSI and GLDAS_all.jpg"
    #
    # draw_scatter(filepath1, filepath2, outname)

    # EVI和GWSA-DSI的双坐标图
    

    # path = "C:\\Users\\11072\\Desktop\\python_code"
    # name1 = "GWSA-DSI_sliding_window_scalar_8_months.txt"
    # filepath1 = os.path.join(path, name1)
    # name2 ="EVI_sliding_window_scalar_8_months.txt"
    # filepath2 = os.path.join(path, name2)
    #
    # outname = "the_relationship_of_EVI_and_GWSA-DSI_in_the_NCP.png"
    #
    # draw_EVI_DSI(filepath1, filepath2, outname)

    # # EVI NC文件
    # top_path = "D:\GRACE_ini_datasets\Final_data"
    # filename1 = "NCP_Terrestial_Monthly_EVI.nc"
    # filepath1 = os.path.join(top_path, filename1)
    # keyword1="CMG 0.25 Deg Monthly EVI"
    #
    # filename2 = "NCP_GWSA-DSI_v2.nc"
    # keyword2 = "groundwater storage anomalies DSI"
    # filepath2 = os.path.join(top_path, filename2)
    # outname = "Seasonly cycles of EVI and DSI.png"
    # draw_seasonly_cycle(filepath1, keyword1, filepath2, keyword2, outname)

    # ## 趋势检验热力图
    # top_path = r"C:\Users\11072\Desktop\python_code\data_pre_process"
    # txt1 = "DSI_Trend_z.txt"
    # path_txt1 = os.path.join(top_path, txt1)
    # txt2 = "DSI_Trend_h.txt"
    # path_txt2 = os.path.join(top_path, txt2)
    # outname = "DSI_hotmap.jpg"
    # anno_hotmap(path_txt1, path_txt2, outname)


if __name__ == "__main__":
    main()
