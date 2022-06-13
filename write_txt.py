# encoding: utf-8
import datetime
import os
import numpy as np

def write_txt(dic_var, outname):
    """

    :param dic_var: 为时间序列的dic
    :return:
    """
    with open(outname, "w") as f:
        # 先写大标题
        f.write(outname.split(".")[0])
        f.write("\n")
        f.write("author: Qingping Liu\n")
        now_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write("created on:  " + now_time)
        f.write("\n")
        # 先写标题
        for key in dic_var.keys():
            f.write(key)
            f.write("\t")
        f.write("\n")
        # 再写入数据
        shp = np.asarray(dic_var[key]).shape
        for i in range(shp[0]):
            for key, value in dic_var.items():
                if key == "Date":
                    print(value[i])
                    f.write(value[i].strftime("%Y-%m"))
                    f.write("\t")
                else:
                    f.write(str(value[i]))
                    f.write("\t")
            f.write("\n")
        f.close()
    
    print("save as txtfile successfully!")
    return 0

