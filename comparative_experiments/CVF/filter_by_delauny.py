import math
import os
from datetime import datetime

import numpy as np

import DPDERL
from comparative_experiments.CVF.sort import Sort
from cul_wcoms_sorted_by_heap import grid_point_to_lon_lat, get_overlap, transform_matrix
from services import read_dem

# 将工作目录改变到脚本所在目录
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# 初始化变量
initial_size = 10000  # 初始三角形的大小
triangles = []  # 存储每个三角形的三个点坐标
edge_buffer = []  # 存储边的列表，每个边是一个包含两个点坐标的列表
point_to_triangles = []  # 存储各簇中心点所连接的三角形
threshold = 1  # 阈值，筛选后的结果除了簇中心坐标外，只保留一个边界点
d = 0  # 当前簇索引
clusters = []  # 簇列表

# 固定半径为600 划分区域分60*60 分辨率为30
radius = 200
# 读取地形文件
demdata = read_dem.Dem("../../data/plain_1000.tif", 5)
# 视点距地面高程
h_stand = 2
clusters_cover = {}
start_time = datetime.now()


def read_clusters():
    global clusters
    # 读取排序过后的簇
    with open("../../output/csv/cvf_sorted_km_mountain_feature.txt", "r") as file:
        for line in file:
            points = list(map(int, line.strip().split(",")))
            clusters.append(points)
    # with open('../../data/clusters_points35.txt', 'r') as file:
    #     for line in file:
    #         clean_line = line.strip().rstrip(',')  # 去除行尾的逗号
    #         oneline = clean_line.split(',')
    #         cluster = [int(item) for item in oneline if item]  # 使用列表推导式并检查空字符串
    #         clusters.append(cluster)


def manage_clusters_cover(key, cluster):
    # 管理一个总覆盖率数组，键为数字编号，值为NumPy二维数组。
    global clusters_cover
    if key in clusters_cover:
        # 如果键已存在，不做操作
        pass
    else:
        # 如果键不存在，创建新覆盖数组
        cover = Sort.cul_all_cover(cluster)
        clusters_cover[key] = cover


def deal():
    global clusters, clusters_cover
    """
    根据各个簇的视点数据进行处理，移除不满足条件的视点并重新调整簇。
    """
    flag = 0
    cnt = 0
    global d  # 全局变量，表示当前处理的簇的索引

    while flag == 0:  # 如果每个簇未满足阈值要求则继续处理
        if len(clusters[d]) // 2 > (threshold + 1):
            date = datetime.now()
            print(f"#{cnt}: {date}:当前分析簇为{d} 簇内剩余点数为{len(clusters[d]) // 2}")
            cnt += 1
            # 簇中叶子节点的位置
            m = len(clusters[d]) - 2
            n = len(clusters[d]) - 1
            # 取出簇中叶子节点的元素值
            a = int(clusters[d][m])
            b = int(clusters[d][n])
            # 将簇d中的叶子节点从该簇中移除
            clusters[d].pop()
            clusters[d].pop()

            lon, lat = grid_point_to_lon_lat(a, b)
            # 计算视点的可视性指数
            coverage = Sort.compute_v(a, b)


            j = 1
            done = False
            while j <= len(clusters) - 1:
                l = (d + j) % len(clusters)  # l代表与当前簇d所比较的簇，判断l簇是否为相邻簇
                flag0 = 0
                # 遍历当前簇d在三角形中的相邻簇
                for p in range(2, len(point_to_triangles[d]), 2):
                    if point_to_triangles[d][p] == clusters[l][0] and point_to_triangles[d][p + 1] == clusters[l][1]:
                        break
                    if p + 2 == len(point_to_triangles[d]):
                        flag0 = 1

                if flag0 == 1:
                    j += 1
                    continue

                current_overlap = Sort.compute_o(a, b, clusters[l])

                if current_overlap == 0:
                    j += 1
                    continue

                current_indicator = current_overlap / (coverage ** 2)
                if current_indicator < 0:
                    print(f"IVC<0:{current_indicator}, {clusters_cover[l]},{coverage} ")
                print(f"当前点在簇{d}中的IVC:{current_indicator}")

                count = 0
                k = 2
                while k <= len(clusters[l]) - 1 and count < 5:
                    row = int(clusters[l][k])
                    col = int(clusters[l][k + 1])
                    lon, lat = grid_point_to_lon_lat(row, col)
                    # 计算视点的可视性指数
                    temp_viewshed, _ = DPDERL.analysis_by_spderl_simplified(lon, lat, demdata.height[row][col],
                                                                            radius, h_stand)
                    temp_coverage = int(np.sum(temp_viewshed))
                    temp_overlap = Sort.compute_o(row, col, clusters[l])
                    temp_indicator = temp_overlap / (temp_coverage ** 2)

                    print(f"簇{d}的相邻簇{l}中点的IVC:{temp_indicator}")

                    if current_indicator < temp_indicator:
                        clusters[l].insert(k, a)
                        clusters[l].insert(k + 1, b)
                        d = l
                        done = True
                        break
                    k += 2
                    count += 1
                if done is True:
                    break
                j += 1
        else:
            d = (d + 1) % len(clusters)

        # 检查是否所有簇都满足阈值要求
        for i in range(len(clusters)):
            if len(clusters[i]) // 2 > (threshold + 1):
                break
            if i == len(clusters) - 1:
                flag = 1


def compute_cr(x1, y1, x2, y2, x3, y3):
    info = [0.0, 0.0, 0.0]
    denominator1 = 2 * ((x3 - x1) * (y2 - y1) - (x2 - x1) * (y3 - y1))
    denominator2 = 2 * ((y3 - y1) * (x2 - x1) - (y2 - y1) * (x3 - x1))

    if denominator1 == 0 or denominator2 == 0:
        return info

    a = ((y2 - y1) * (y3 * y3 - y1 * y1 + x3 * x3 - x1 * x1) - (y3 - y1) * (
            y2 * y2 - y1 * y1 + x2 * x2 - x1 * x1)) / denominator1
    b = ((x2 - x1) * (x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1) - (x3 - x1) * (
            x2 * x2 - x1 * x1 + y2 * y2 - y1 * y1)) / denominator2
    r = math.sqrt((x1 - a) ** 2 + (y1 - b) ** 2)

    info[0] = a
    info[1] = b
    info[2] = r
    return info


def delaunay(points):
    triangles = []  # 初始化存储所有三角形的二维列表
    edge_buffer = []  # 边缘缓存列表

    # 添加一个“超级三角形”包含了所有其他点
    triangles.append([
        -initial_size, -initial_size,
        initial_size, -initial_size,
        0, initial_size
    ])
    triangle_count = 0

    # 遍历每个点以构建三角形
    for x, y in points:
        edge_buffer.clear()
        j = 0
        while j < len(triangles):
            # 计算当前三角形的外接圆半径和中心
            circle_radius = compute_cr(*triangles[j])
            if circle_radius[2] == 0:  # 如果半径为0，则跳过此三角形
                triangles.pop(j)
                triangle_count -= 1
                continue

            # 如果点在外接圆内，则不是Delaunay三角形
            if ((x - circle_radius[0]) ** 2 + (y - circle_radius[1]) ** 2) < circle_radius[2] ** 2:
                x1, y1, x2, y2, x3, y3 = triangles[j]
                edge_buffer.extend([(x1, y1, x2, y2), (x1, y1, x3, y3), (x2, y2, x3, y3)])
                triangles.pop(j)
                triangle_count -= 1
            else:
                j += 1

        # 去重边
        l = 0
        while l < len(edge_buffer):
            removed = False  # Flag to check if we have removed elements
            h = 0
            while h < len(edge_buffer):
                if h != l:
                    el = edge_buffer[l]
                    eh = edge_buffer[h]

                    if ((el[0] == eh[0] and el[1] == eh[1] or el[0] == eh[2] and el[1] == eh[3]) and
                            (el[2] == eh[0] and el[3] == eh[1] or el[2] == eh[2] and el[3] == eh[3])):

                        if h > l:
                            edge_buffer.pop(h)
                            edge_buffer.pop(l)
                        else:
                            edge_buffer.pop(l)
                            edge_buffer.pop(h)
                        removed = True
                        break
                if not removed:
                    h += 1

            if not removed:
                l += 1

        print(edge_buffer)

        # 使用唯一的边和当前点创建新的三角形
        for edge in edge_buffer:
            triangles.append([int(x), int(y), edge[0], edge[1], edge[2], edge[3]])
            triangle_count += 1

    # 打印所有三角形
    count = 0
    for triangle in triangles:
        print(" ".join(map(str, triangle)))
        count += 1
    print(count)

    return triangles


def main():
    global point_to_triangles
    read_clusters()

    # 计算视点的可视信息
    count = sum(len(cluster) - 2 for cluster in clusters) // 2
    print(count)

    # 获取聚类中心
    points = np.array([[int(cluster[0]), int(cluster[1])] for cluster in clusters], dtype=int)
    print("开始构建德劳内三角形")
    triangles = delaunay(points)
    print("构建德劳内三角形完毕")

    # 假设 `points` 是一个列表，其中包含坐标点的子列表，如 [x, y]。
    # 假设 `triangles` 是一个列表，每个子列表表示一个三角形，包含三个坐标点的六个值，如 [x1, y1, x2, y2, x3, y3]。

    point_to_triangles = []  # 用于存储每个点及其对应的三角形中的其他点

    # 遍历所有点
    for i in range(len(points)):
        tempx = points[i][0]
        tempy = points[i][1]
        point_to_triangles.append([int(tempx), int(tempy)])  # 初始化时，先将当前点添加到列表中

        # 遍历所有三角形，寻找包含当前点的三角形
        for j in range(len(triangles)):
            # 遍历当前三角形的每个点（每两个数是一个坐标点）
            for k in range(0, len(triangles[j]), 2):
                # 检查当前点是否为三角形的一个顶点
                if tempx == triangles[j][k] and tempy == triangles[j][k + 1]:
                    # 如果当前点是三角形的一部分，则添加除当前点外的其他两个点
                    for l in range(0, len(triangles[j]), 2):
                        if l != k:  # 确保不重复添加当前点
                            point_to_triangles[i].append(triangles[j][l])
                            point_to_triangles[i].append(triangles[j][l + 1])

    # 输出结果，查看每个点对应的其他三角形点
    print(point_to_triangles)

    print("根据德劳内三角形找出相邻簇并记录")
    # 过滤和重分类c
    deal()
    end_time = datetime.now()
    cvf_time = str(end_time-start_time).replace(":", '-')
    print(cvf_time)

    # 写入文件
    output_path = f"../../output/csv/filtered_cvf_mountain_feature_{cvf_time}.csv"
    with open(output_path, 'w') as file:
        print("把筛选后的簇写入文本文件中")
        for cluster in clusters:
            file.write(",".join(map(str, cluster)) + "\n")
    print("写入文本结束")


if __name__ == "__main__":
    main()
