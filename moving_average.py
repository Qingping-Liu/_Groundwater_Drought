#encoding:utf-8

import datetime
import os
from itertools import islice
from collections import deque
import numpy as np


# 滑动平均法
def movingaverage(data, subset_size, data_is_list = None, avoid_fp_drift = False):
	'''Return the moving averages of the data, with a window size of
	`subset_size`.  `subset_size` must be an integer greater than 0 and
	less than the length of the input data, or a ValueError will be raised.

	`data_is_list` can be used to tune the algorithm for list or iteratable
	as an input.  The default value, `None` will auto-detect this.
	The algorithm used if `data` is a list is almost twice as fast as if
	it is an iteratable.

	`avoid_fp_drift`, if True (the default) sums every sub-set rather than
	keeping a "rolling sum" (which may be subject to floating-point drift).
	While more correct, it is also dramatically slower for subset sizes
	much larger than 20.

	NOTE: You really should consider setting `avoid_fp_drift = False` unless
	you are dealing with very small numbers (say, far smaller than 0.00001)
	or require extreme accuracy at the cost of execution time.  For
	`subset_size` < 20, the performance difference is very small.
	'''
 
	if subset_size < 1:
		raise ValueError('subset_size must be 1 or larger')

	if data_is_list is None:
		data_is_list = hasattr(data, '__getslice__')

	divisor = float(subset_size)
	if data_is_list:
		#  This only works if we can re-access old elements, but is much faster.
		#  In other words, it can't be just an iterable, it needs to be a list.

		if subset_size > len(data):
			raise ValueError('subset_size must be smaller than data set size')

		if avoid_fp_drift:
			for x in range(subset_size, len(data) + 1):
				yield sum(data[x - subset_size:x]) / divisor
		else:
			cur = sum(data[0:subset_size])
            # 定义了一个迭代器
			yield cur / divisor
			for x in range(subset_size, len(data)):
				cur += data[x] - data[x - subset_size]
				yield cur / divisor
	else:
		#  Based on the recipe at:
		#     http://docs.python.org/library/collections.html#deque-recipes
        # 优雅的分块读取
		it = iter(data)     # 使其变成迭代器
		d = deque(islice(it, subset_size))

		if subset_size > len(d):
			raise ValueError('subset_size must be smaller than data set size')

		if avoid_fp_drift:
			yield sum(d) / divisor
			for elem in it:
				d.popleft()
				d.append(elem)
				yield sum(d) / divisor
		else:
			s = sum(d)
			yield s / divisor
			for elem in it:
				s += elem - d.popleft()
				d.append(elem)
				yield s / divisor



def main():
	from write_txt import write_txt

	# # 处理EVI滑动平均
    # top_path = "D:\\GRACE_ini_datasets\\Final_data\\fig"
    # name = "The changes of EVI in the NCP from 2002 to 2021.txt"
    # filepath = os.path.join(top_path, name)
    #
    # date = []
    # ave = []
    # with open(filepath, "r") as f:
    #     readlines = f.readlines()[4:]
    #     # 跳过前5行
    #     for read in readlines:
    #         read = read.split("\t")
    #         # print(read)
    #         date.append(datetime.datetime.strptime(read[0], "%Y-%m"))
    #         ave.append(float(read[1]))
    #
    #     f.close()
    # new_ave = movingaverage(ave, 8, data_is_list=True)
    # # 创建字典
    # dic_var = {}
    # dic_var["date"] = np.asarray(date[:])
    # dic_var["ave"] = np.asarray([i for i in new_ave])
    #
    # new_name = "EVI_sliding_window_scalar_8_months.txt"
    # write_txt(dic_var, new_name)

    # 处理dsi
	toppath = "D:\\GRACE_ini_datasets\\Final_data\\fig"
	name = "The GWSA_DSI changes in the NCP from 2002 to 2021.txt"
	filepath = os.path.join(toppath, name)
	with open(filepath, "r"):
		date = []
		ave = []
		with open(filepath, "r") as f:
			readlines = f.readlines()[4:]
			# 跳过前5行
			for read in readlines:
				read = read.split("\t")
				# print(read)
				date.append(datetime.datetime.strptime(read[0], "%Y-%m"))
				ave.append(float(read[1]))

		f.close()
	new_ave = movingaverage(ave, 8, data_is_list=True)
	# 创建字典
	dic_var = {}
	dic_var["date"] = np.asarray(date[:])
	dic_var["ave"] = np.asarray([i for i in new_ave])

	new_name = "GWSA-DSI_sliding_window_scalar_8_months.txt"
	write_txt(dic_var, new_name)

	return 0

##########################
if __name__ == '__main__':

	main()

