# encoding: utf-8
from scipy import stats
import numpy as np
from matplotlib import pyplot as plt
import os
import datetime


def sk(data):
    n = len(data)
    Sk = [0]
    UFk = [0]
    s = 0
    E = [0]
    Var = [0]
    for i in range(1, n):
        for j in range(i):
            if data[i] > data[j]:
                s = s + 1
            else:
                s = s + 0
        Sk.append(s)
        E.append((i + 1) * (i + 2) / 4)  # Sk[i]的均值
        Var.append((i + 1) * i * (2 * (i + 1) + 5) / 72)  # Sk[i]的方差
        UFk.append((Sk[i] - E[i]) / np.sqrt(Var[i]))
    UFk = np.array(UFk)
    return UFk


# a为置信度
def MK(data, a):
    ufk = sk(data)  # 顺序列
    ubk1 = sk(data[::-1])  # 逆序列
    ubk = -ubk1[::-1]  # 逆转逆序列

    # 输出突变点的位置
    p = []
    u = ufk - ubk
    for i in range(1, len(ufk)):
        if u[i - 1] * u[i] < 0:
            p.append(i)
    if p:
        print("突变点位置：", p)
    else:
        print("未检测到突变点")

    # 画图
    conf_intveral = stats.norm.interval(a, loc=0, scale=1)  # 获取置信区间

    plt.figure(figsize=(10, 5))
    plt.plot(range(len(data)), ufk, label='UFk', color='r')
    plt.plot(range(len(data)), ubk, label='UBk', color='b')
    plt.ylabel('UFk-UBk', fontsize=25)
    x_lim = plt.xlim()
    plt.ylim([-6, 7])
    plt.plot(x_lim, [conf_intveral[0], conf_intveral[0]], '--', color='r', label='95%显著区间')
    plt.plot(x_lim, [conf_intveral[1], conf_intveral[1]], '--', color='r')
    plt.axhline(0, ls="--", c="k")
    plt.legend(loc='upper center', frameon=False, ncol=3, fontsize=20)  # 图例
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.show()

    return p

# 突变点检测
def Pettitt(data):
    data = np.array(data)
    n = len(data)
    sk = [0]

    for i in range(1, n):
        s = 0
        for j in range(i):
            if data[i] > data[j]:
                s = s + 1
            if data[i] < data[j]:
                s = s - 1
            else:
                s = s + 0
        sk.append(s)

    k = np.max(sk)
    kt = sk.index(max(sk))
    print(k, kt)
    p = 2 * np.exp((-6 * (k ** 2)) / (n ** 3 + n ** 2))

    if p <= 0.05:
        a = '显著'
    else:
        a = '不显著'
    print('突变点位置:%s, %s' % (kt, a))

    # 画图
    plt.plot(range(len(data)), data)
    plt.plot([0, kt], [np.mean(data[0:kt]), np.mean(data[0:kt])], 'm--', color='r')
    plt.plot([kt, len(data)], [np.mean(data[kt::]), np.mean(data[kt::])], 'm--', color='r')
    plt.axvline(x=kt, ls="--", c="g")
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    return kt  # ,Pettitt_result


# 滑动T检验
def MTT(data, step):
    n = len(data)
    v = step + step - 2  # 自由度
    t = np.zeros((n - step - step + 1))
    ss = np.sqrt(1 / step + 1 / step)

    ttest = 3.4663  # step=5,alpha=0.05，这个需要根据需要查表做改动

    for i in range(len(t)):
        n1 = data[i:i + step]
        n2 = data[i + step:i + step + step]
        x1 = np.mean(n1)  # 平均值
        x2 = np.mean(n2)
        s1 = np.std(n1)  # 方差
        s2 = np.std(n2)
        s = np.sqrt((step * s1 * s1 + step * s2 * s2) / v)
        t[i] = (x1 - x2) / (s * ss)

    plt.plot(t, "m.-", label="Move T value(step=%s)" % step)
    plt.axhline(0, ls="--", c="k")
    plt.axhline(ttest, ls="--", c="r", label='95% significant')
    plt.axhline(-ttest, ls="--", c="r")
    plt.legend(loc='upper center', frameon=False, ncol=2, fontsize=14)  # 图例

    return t






import scipy.special
import math, sys
import random

np.random.seed(0)

def B_G(data):
      X = list(data["X"])

      P0 = 0.95
      MIN_SUB_LENGTH = 100
      N = len(X)
      flag = [0] * N
      results = {}

      T, Tmax, Tmax_idx, PTmax = FindTmax(X)
      print(PTmax)
      if PTmax < P0:
        print("No find!")
        sys.exit()

      flag[Tmax_idx] = 1
      break_idx = [0, Tmax_idx, len(X)-1]
      results[Tmax_idx] = {"Tmax":Tmax, "T":T, "PTmax":PTmax, "start_idx":0, "end_idx":len(X)-1, "break_order":1}





if __name__ == "__main__":

    top_path = "D:\\GRACE_ini_datasets\\Final_data\\NCP_DSI_v3"
    txt = "The GWSA_DSI changes in the NCP from 2002 to 2021_v4.txt"

    file_path = os.path.join(top_path, txt)

    ave = []
    date =[]
    std = []
    # 将数据导入
    with open(file_path, "r") as f:
        readlines = f.readlines()[4:]
        # 跳过前5行
        for read in readlines:
            read = read.split("\t")
            # print(read)
            date.append(datetime.datetime.strptime(read[0], "%Y-%m"))
            ave.append(float(read[1]))
            std.append(float(read[2]))

    # # 输入数据和置信度即可
    tp_1 = MK(ave, 0.95)
    for i in range(len(tp_1)):
        print(date[tp_1[i]])

    # tp_2 = Pettitt(ave)
    # print(date[tp_2])

    # tp_3 = MTT(ave, 30)
    # for i in range(len(tp_3)):
    #     print(date[tp_3[i]])