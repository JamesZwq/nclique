import matplotlib.pyplot as plt
import numpy as np

# 读取文件中的数据（假设每行只有一个数字）
data = []
dataFile = "/Users/zhangwenqian/UNSW/pivoter/adj_list.txt"
data = np.loadtxt(dataFile, dtype=int)
# remove if data = 0
data = data[data != 0]
maxData = np.max(data)
minData = np.min(data)
print("maxData: ", maxData)
print("minData: ", minData)
plt.figure(figsize=(20, 6))
# bins可以根据数据特点调整，比如这里分100个区间
plt.hist(data, bins=maxData, edgecolor='black')
# log scale
plt.yscale('log')
plt.savefig(dataFile + ".png", dpi=300)
#  print the number of date for each value