import os
from datetime import datetime

import numpy as np
import random
import math

import DPDERL
from cul_wcoms_sorted_by_heap import grid_point_to_lon_lat, judge_extent_beyond_area, transform_matrix
from services import read_dem

# 将工作目录改变到脚本所在目录
os.chdir(os.path.dirname(os.path.abspath(__file__)))
# 定义全局变量
M, N = 1000, 1000

# 保留点数量
number = 45
T = 100  # 起始温度
Tmin = 1e-5  # 终止温度
k = 10  # 迭代的次数
features = []
finalresult = np.zeros((number, 2), dtype=int)
eur = 2000
re = 5
delta = 0.9
radius = 200
# 读取地形文件
dem = read_dem.Dem("../../data/plain_1000.tif", 5)
demdata = dem.height
# 视点距地面高程
h_stand = 2


# 提取地形特征点
def extract_features():
    global features
    # 假设这个函数以某种方式填充 features
    datafilepath = "../../output/csv/sa_mountain_feature.txt"
    features = np.loadtxt(fname=datafilepath, delimiter=',')
    # 使用布尔索引过滤掉需要删除的点
    features = features[~np.array([judge_extent_beyond_area(x, y) for x, y in features])]


# 模拟退火算法
def simulated_annealing():
    global finalresult
    coverage = compute_all_view(finalresult)
    temp_location = np.copy(finalresult)
    print(coverage)
    t = T
    cnt = 0
    while t > Tmin:
        date = datetime.now()
        print(f"#{cnt}:{date} {t}℃--->{Tmin}℃")
        cnt += 1
        for _ in range(k):
            for j in range(len(finalresult)):
                flag = 0
                while flag != 1:
                    r = random.randint(0, len(features) - 1)
                    tempx, tempy = features[r]
                    tempx = int(tempx)
                    tempy = int(tempy)
                    # 使用 distance 函数计算新视点与当前视点的距离，如果距离小于阈值 eur，则接受这个新视点，并更新 temp_location。
                    if distance(tempx, tempy, demdata[tempx][tempy], temp_location[j][0], temp_location[j][1],
                                demdata[temp_location[j][0]][temp_location[j][1]]) < eur:
                        temp_location[j] = [tempx, tempy]
                        flag = 1

            new_coverage = compute_all_view(temp_location)
            if new_coverage > coverage:
                finalresult[:] = temp_location[:]
                coverage = new_coverage
                print(coverage)
            else:
                p = 1 / (1 + math.exp(-(new_coverage - coverage) / T))
                if random.random() < p:
                    finalresult[:] = temp_location[:]
                    coverage = new_coverage
        t *= delta


def compute_all_view(temp):
    """
       计算给定视点集合的总视域覆盖区域。

       参数:
       - temp: 包含视点坐标的二维数组，每行为一个视点坐标 [x, y]

       返回:
       - percent: 所有视点的视域覆盖区域与总区域的比例
       """
    allview = np.zeros(demdata.shape)  # 创建一个二维数组来累积所有视点的视域

    # 遍历每个视点
    for i in range(len(temp)):
        x, y = temp[i]
        lon, lat = grid_point_to_lon_lat(x, y)

        # 计算视点的可视性指数
        tempview, _ = DPDERL.analysis_by_spderl_simplified(lon, lat, demdata[x][y],
                                                           radius, h_stand)
        tempview = transform_matrix(tempview, [x, y])
        # 将当前视点的视域累加到总视域中
        allview += tempview

    # 计算总视域中被视点覆盖的单元格数量
    count = np.sum(allview)

    # 计算覆盖比例
    percent = count / (M * N)
    return percent


def distance(x1, y1, h1, x2, y2, h2):
    dis = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (h1 - h2) ** 2)
    return dis * re


# 主函数
if __name__ == "__main__":
    last_time = datetime.now()
    extract_features()
    # 随机生成初始解
    for i in range(number):
        rand_index = random.randint(0, len(features) - 1)
        finalresult[i] = features[rand_index]
    simulated_annealing()
    print("成功保存模拟退火筛选视点坐标")
    now_time = datetime.now()
    sa_time = str(now_time - last_time).replace(":", "-")
    print(f"SA算法用时:{sa_time}")
    # 保存这个数组到文本文件中
    np.savetxt(f'../../output/csv/sa_result_{sa_time}.txt', finalresult, fmt='%d', delimiter=',')
    print(f"保留点：sa_result_{sa_time}.txt")
