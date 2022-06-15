# encoding:utf-8

# import modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# 读取excel文件



filepath = "C:\\Users\\11072\\Desktop\\GRACE_dataset\\water_cons\\Water_Supply.xlsx"

# 读取excel文件
f = pd.ExcelFile(filepath)
data = {}
for i in f.sheet_names:
    # 默认第一行为表头
    data[i] = pd.read_excel(filepath, sheet_name=i)
    

beijing = data["anhui"]
beijing["Year"] = pd.to_datetime(beijing["Year"]).dt.year

x_group =[str(i) for i in beijing["Year"].values] 
surface = beijing["Surface water supply"].values
groundwater = beijing["Groundwater supply"].values
reclaimed = beijing["Reclaimed water supply"].values
diverted = beijing["Diverted water supply"].values

# 计算bottom值
b1 = surface
b2 = b1 + groundwater
b3 = b2 + reclaimed


print(x_group)
print(beijing)

# 绘图
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot()

# 对每一个供水来源绘制一个bar
width = 0.6
p1 = ax.bar(x_group, surface,  width, label='Surface',color='tab:red')
p2 = ax.bar(x_group, groundwater, width, bottom=b1,  label='Ground',color='tab:blue')
p3 = ax.bar(x_group, reclaimed, width, bottom=b2, label='Reclaimed',color='tab:green')
p4 = ax.bar(x_group, diverted, width, bottom=b3, label='Diverted',color='tab:orange')


# 绘制折线图
ax2 = ax.twinx()
surf_p= beijing["Surface percentage"].values
ground_p = beijing["Groundwater percentage"].values
reclaim_p = beijing["Reclaimed percentage"].values
divert_p = beijing["Diverted percentage"].values

p21 = ax2.plot(x_group, surf_p, label="Fraction", color="r", lw=2.5)
p22 = ax2.plot(x_group, ground_p, label="Fraction", color="b", lw=2.5)
p23 = ax2.plot(x_group, reclaim_p, label="Fraction", color="g", lw=2.5)
p24 = ax2.plot(x_group, divert_p, label="Fraction", color="orange", lw=2.5)

ax2.set_ylim(0, 1)

ax.set_title("Water Supply Structure in Anhui")
ax.set_xlabel("Year")
ax.set_ylabel("Annual Water Use[bn m3]")

ax2.set_ylabel("Fraction")

plt.legend([p1, p2, p3, p4], ["Surface", "Ground", "Reclaimed", "Diverted"])
plt.savefig("Water Supply Structure in Anhui.jpg", dpi=300)
print("Done!")