# encoding:utf-8

#  Import modules
import pdfplumber
import pandas as pd
import os 



print("start...")

top_path = "D:\\GRACE_ini_datasets\\initial_INPUT_data\\water_use\\hebei"
filename = "hebei2010.pdf"
filepath = os.path.join(top_path, filename)
# 读取pdf文件
pdf = pdfplumber.open(filepath)

#  访问供水和用水数据
page = pdf.pages[11]

#  自动读取表格信息，返回列表
table = page.extract_table()


text = page.extract_text()

row = len(table)
col = len(table[-1])

print(len(pdf.pages))
print(table)

outname = filename.split(sep=".")[0] + ".txt"
with open(outname, "w") as f:
    for i in range(row):
        li = table[i]
        for j in range(col): 
            try:
                f.write(li[j])
            except:
                pass
            f.write("\t")
        f.write("\n")
    
    f.close()
    
print("Done!")