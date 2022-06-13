# encoding:utf-8

# B-G分割算法

from scipy.special import betainc
import numpy as np
import sys
import os
import datetime



print(betainc(0.2, 3.5, 1))

TMAX_id = []
P_value = []


def tmax(data, left, right):
    data = np.asarray(data[left:right])
    N = right - left
    if N <= 30:
        print("can't split anymore!")
        return 0
    # Tmax
    Tmax_id = int(FindTmax(data)[2])
    P_Tmax = FindTmax(data)[3]

    if right - Tmax_id >= 30 and Tmax_id - left >= 30:
        middle = Tmax_id

        tmax(data, left, middle)
        tmax(data, middle + 1, right)

        print(Tmax_id)
        print(P_Tmax)

        return TMAX_id.append(Tmax_id), P_value.append(P_Tmax)


def FindTmax(data):
    #np.random.seed(0)
    data = np.asarray(data)
    N = len(data)

    T = np.empty(N - 1)
    for i in range(N - 1):
        Nl = len(data[:i + 1])
        Nr = len(data[i + 1:])
        U_l = np.mean(data[:i + 1])
        U_r = np.mean(data[i + 1:])
        S_l = np.var(data[:i + 1])
        S_r = np.var(data[i + 1:])
        S_D = pow((S_l ** 2 + S_r ** 2) / (Nl + Nr - 2), 0.5) * pow((1 / Nl + 1 / Nr), 0.5)
        T[i] = (U_l - U_r) / S_D

    # Tmax
    Tmax = np.max(T)
    Tmax_id = np.where(T == Tmax)[0]
    Xt = data[Tmax_id]

    # beta参数
    theta = 0.4
    alpha = 4.19 * np.log(N) - 11.54

    # 计算 p
    mu = len(data) - 2
    x = mu / (mu + Tmax)
    PTmax = betainc(theta, alpha, x)

    return Xt, Tmax, Tmax_id, PTmax

# 假设性检验，通过假设性检验的结果才可靠
def BG_test(data, alpha):
    """
    :param data: 一维时间序列
    :param alpha: 置信度
    :return:

    """
    N = len(data)

    tmax(data, 0, N)


    assert len(TMAX_id) == len(P_value)

    outp_id = []
    outp_id.append(0)
    for i in range(len(P_value)):
        if P_value[i] > alpha:
            outp_id.append(TMAX_id[i])

    outp_id.append(N)

    segements= np.zeros((1, N), dtype=np.float32)
    for i in range(len(outp_id)-1):

        start = outp_id[i]
        end = outp_id[i+1]
        mea = np.mean(data[start: end])
        segements[0, start: end] = mea.copy()

    print(outp_id) 
        

    # 返回的时序列的均值
    return segements

if __name__ == "__main__":
    
    # # DSI 的时间序列
    path = "D:\\GRACE_ini_datasets\\Final_data\\NCP_DSI_v3"
    name = "The GWSA_DSI changes in the NCP from 2002 to 2021_v4.txt"
    filepath = os.path.join(path, name)
    
    ave = []
    date =[]
    std = []
    # 将数据导入
    with open(filepath, "r") as f:
        readlines = f.readlines()[4:]
        # 跳过前5行
        for read in readlines:
            read = read.split("\t")
            # print(read)
            date.append(datetime.datetime.strptime(read[0], "%Y-%m"))
            ave.append(float(read[1]))
            std.append(float(read[2]))

    oupid = BG_test(ave, 0.95)



    print(oupid)
    

