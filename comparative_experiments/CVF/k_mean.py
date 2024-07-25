import numpy as np
import random
import math
from services import read_dem

datafilepath = "../../output/csv/cvf_mountain_feature.txt"
point = np.loadtxt(fname=datafilepath, delimiter=',')

# 分辨率
re = 5

dem_data = read_dem.Dem("../../data/plain_1000.tif", re)
demdata = dem_data.height

# 簇量
number = 45

# 簇中心
centers = []

# 簇
clusters = []





def KM(demdata, points, number=30):
    # 随机初始化簇中心
    global centers
    for _ in range(number):
        random_index = random.randint(0, len(points) - 1)
        p = points[random_index]
        center = [p[0], p[1], demdata[int(p[0])][int(p[1])]]
        centers.append(center)

    flag = 0
    while flag != 1:
        clusters.clear()
        KMeans(demdata, points)

        for i, cluster in enumerate(clusters):
            if cluster[0] == 0 and cluster[1] == 0:
                break
            if i + 1 == len(clusters):
                flag = 1

    return clusters


def KMeans(demdata, points):
    # 初始化簇
    global clusters
    cluster = [[] for _ in range(len(centers))]
    flag = 0
    while flag < 100:
        print(f"flag:{flag}")
        cluster = [[] for _ in range(len(centers))]

        # 将离中心点最近的候选点赋值到对应簇中
        for point in points:
            row, col = point
            row = int(row)
            col = int(col)
            min_distance = com_distance(row, col, demdata[row][col], centers[0][0], centers[0][1], centers[0][2])
            tempj = 0
            for j in range(len(centers)):
                distance = com_distance(row, col, demdata[row][col], centers[j][0], centers[j][1], centers[j][2])
                if distance < min_distance:
                    min_distance = distance
                    tempj = j
            cluster[tempj].append(row)
            cluster[tempj].append(col)

        # 重新分配簇中心
        for i in range(len(cluster)):
            sum_row = sum(cluster[i][j] for j in range(2, len(cluster[i]), 2))
            sum_col = sum(cluster[i][j] for j in range(3, len(cluster[i]), 2))
            sum_height = sum(demdata[cluster[i][j]][cluster[i][j + 1]] for j in range(2, len(cluster[i]), 2))
            num_points = len(cluster[i]) // 2 - 1
            if num_points > 0:
                cen_row = sum_row / num_points
                cen_col = sum_col / num_points
                cen_height = sum_height / num_points
                new_center = [cen_row, cen_col, cen_height]
                centers[i] = new_center

        flag += 1

    clusters = []
    for i in range(len(cluster)):
        # 为这个簇创建一个新的列表
        clusters.append([])

        # 添加簇中心坐标，如果需要的话，转换为整数（假设center包含浮点值）
        clusters[i].append(int(centers[i][0]))
        clusters[i].append(int(centers[i][1]))

        # 将当前簇的所有点添加到新列表中
        for point in cluster[i]:
            clusters[i].append(point)


def com_distance(row, col, h1, i, j, h2):
    distance = math.sqrt((row - i) ** 2 * re ** 2 + (col - j) ** 2 * re ** 2 + (h1 - h2) ** 2)
    return distance


KM(demdata, point, number)

# 存储聚类后的簇
with open('../../output/csv/cvf_KM_mountain_feature.txt', 'w') as file:
    # 遍历列表中的每个子列表
    for sublist in clusters:
        # 将子列表中的每个元素转换为字符串，并用逗号连接
        line = ','.join(map(str, sublist))
        # 将整理好的字符串写入文件，并在每行的末尾添加换行符
        file.write(line + '\n')
print(f"已保存好{number}个簇")
