# encoding:utf-8

# import modules
from datetime import datetime
import numpy as np
import pandas as pd
import sys

sys.path.append("C:\\Users\\11072\\Desktop\\python_code\\data_pre_process")
from write_txt import write_txt

# 根据city_id匹配来提取数据
def water_use(filepath, li=""):
    f = pd.ExcelFile(filepath)
    data={}
    for i in f.sheet_names:
        # 默认第一行为表头
        data[i] = pd.read_excel(filepath, sheet_name=i)
    # 删除sheet D1工作表中第二行的数据
    data["D1"].drop(0, inplace=True)
    
    datasets = data["D1"].loc[:16709, ["City_ID", "Year", "IRR(km3 yr-1)", "IND(km3 yr-1)", "URB", "RUR", "Total water use"]]
    datasets["Year"] = datasets["Year"].apply(int)
    datasets["Year"] = datasets["Year"].apply(str)
    datasets["Year"] = pd.to_datetime(datasets["Year"]).dt.year
    lookupID = data["NCP"].loc[:, ["Perfecture", "Province_n"]]
    
    print(lookupID)
    
    
    pcc = ["Henan", "Anhui", "Hebei", "Shandong", "Jiangsu", "Henan", "Tianjin", "Shanxi"]
        
    # 新建一个存放结果的dataframes
    water_use = pd.DataFrame(datasets["Year"][:49], columns=["Year"])
    # 将时间设置为行索引
    water_use.set_index("Year", inplace=True)
    water_use["Irrigation"] = 0
    water_use["Industry"] =0
    water_use["Urban"] = 0
    water_use["Rural"] = 0
    water_use["Total_wateruse"] = 0
        
    

    print(water_use)
    
    
    for i in range(len(lookupID["Perfecture"].values)):
        cityid = lookupID.iloc[i, 0]
        extracts = datasets[datasets["City_ID"] == cityid]
        extracts.set_index("Year", inplace=True)
        print(extracts)

        # 省份
        province = lookupID.iloc[i, 1]
    
        if province =="Jiangsu":
            water_use["Irrigation"] = water_use["Irrigation"] + extracts["IRR(km3 yr-1)"]
            water_use["Industry"] = water_use["Industry"] + extracts["IND(km3 yr-1)"]  
            water_use["Urban"] = water_use["Urban"] + extracts["URB"]
            water_use["Rural"] = water_use["Rural"] + extracts["RUR"]
            water_use["Total_wateruse"] = water_use["Total_wateruse"] + extracts["Total water use"]
        

        

        

        print(water_use)
    # 把索引还原
    water_use.reset_index(inplace=True)
    print(water_use)
    water_use.rename(columns={'Year':'Date'}, inplace=True) 
    #  转换成字典
    dic_var = water_use.to_dict('list')
    #  把字典的year处理处理成datetime格式
    dic_var["Date"] = [datetime.strptime(str(i), "%Y") for i in dic_var["Date"]]
    
    return dic_var
    
    

if __name__=="__main__":
    # 写成txt文件


    filepath = "C:\\Users\\11072\\Desktop\\GRACE_dataset\\water_cons\\Zhou et al_2020_PNAS_dataset.xlsx"
    watercon = water_use(filepath)
    
    outname = "the water consumption from 1965 to 2013_jiangsu.txt"
    # 保存成txt文件
    write_txt(watercon, outname)

