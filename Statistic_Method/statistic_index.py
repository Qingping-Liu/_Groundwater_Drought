# encoding:utf-8

import numpy as np


# 计算均值、标准差、极大值、极小值
def statistic_index(array):
    """
    :param array: 一维的时间序列数组,数据类型为float型
    :return:
    """
    # 计算均值

    array = np.asarray(array)
    n = len(array)  # 数组的长度
    _max = np.max(array)    # 计算最大值
    average = np.sum(array)/n   # 数据的均值
    std = pow(np.sum((array - average)**2)/n, 0.5)  # 计算标准差
    _min = np.min(array)        # 计算最小值

    return average, std, _min, _max