# 文件名
import os

from comparative_experiments.CVF.sort import Sort
from services import read_dem

# 将工作目录改变到脚本所在目录
os.chdir(os.path.dirname(os.path.abspath(__file__)))

filename = '../../output/csv/cvf_KM_mountain_feature.txt'

# 创建一个空列表来存储二维列表
clusters = []

# 固定半径为600 划分区域分60*60 分辨率为30
radius = 200
# 读取地形文件
dem = read_dem.Dem("../../data/plain_1000.tif", 5)

# 视点距地面高程
h_stand = 2
# 打开文件进行读取
with open(filename, 'r') as file:
    # 逐行读取文件
    for line in file:
        # 去掉每行的换行符并以逗号分割，然后将每个元素转换为整数
        int_list = [int(x) for x in line.strip().split(',')]
        # 将转换后的整数列表添加到列表a中
        clusters.append(int_list)


def main():
    clusters_after_sort = Sort.descend_sort(clusters, dem)
    return clusters_after_sort

# # 存储排序聚类后的簇
# with open('../../output/csv/cvf_sorted_km_mountain_feature.txt', 'w') as file:
#     # 遍历列表中的每个子列表
#     for sublist in clusters:
#         # 将子列表中的每个元素转换为字符串，并用逗号连接
#         line = ','.join(map(str, sublist))
#         # 将整理好的字符串写入文件，并在每行的末尾添加换行符
#         file.write(line + '\n')
