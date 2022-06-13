# encoding: utf-8

import numpy as np

# 重力卫星干旱指数计算
def drought_index(array):
    """
    :param array: 一维的时间序列数组,数据类型为float型
    :return:
    """
    # 计算均值
    array = np.asarray(array)   # 可能是个list
    n = len(array)  # 数组的长度
    if array[0] != -9999.0:     # 判断是否是缺失值
        average = np.sum(array)/n   # 数据的均值
        std = pow(np.sum((array - average)**2)/n, 0.5)
        DSI = (array - average) / std
    else:
        DSI = array       # 以原来的缺失值填充
    return DSI