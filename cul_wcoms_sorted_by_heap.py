import copy
import gc
import heapq
import json
from datetime import datetime

import numpy as np

import DPDERL
from services import read_dem
import os

# 将工作目录改变到脚本所在目录
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# 固定半径为600 划分区域分60*60 分辨率为30
radius = 200
# 读取地形文件
dem = read_dem.Dem("data/plain_1000.tif", 5)
# 视点距地面高程
h_stand = 2

# 公式权重
#w1 0.47 w2 0.53
#w1 0.3 w2 0.7
w1 = 0.5
w2 = 0.5

# 簇的总覆盖率
cover_clusters = []


# 加权覆盖重叠度量公式：F=w1×Coverage−w2×Overlap
def cul_wcom(coverage, overlap):
    return w1 * coverage - w2 * overlap


# 将格网点坐标转为经纬度
def grid_point_to_lon_lat(x, y):
    lat = dem.start_y + y * dem.dy
    lon = dem.start_x + x * dem.dx
    return lon, lat


# 重叠率计算有误，这里的viewshed只是以观察点为中心的相对可视矩阵(2*x_grid,2*y_grid)并没有结合到地形3600*3600,需进行矩阵转化
def transform_matrix(a, point):
    """
    将矩阵a放置到一个新的3600x3600矩阵b中，使a的中心与b中的p1点对齐。

    参数:
    a (np.array): 输入的二维矩阵a。
    p1 (tuple): 矩阵b中的目标中心坐标点，例如 (1800, 1800)。

    返回:
    np.array: 调整后的矩阵b。
    """
    # 创建目标矩阵b
    b = np.zeros(dem.height.shape, np.uint16)

    # 计算a的中心坐标
    center_a = (a.shape[0] // 2, a.shape[1] // 2)

    # 计算在b中a应该开始填充的起始索引
    start_index = (int(point[0]) - center_a[0], int(point[1]) - center_a[1])

    # 计算a在b中的填充区域的终止索引
    end_index = (start_index[0] + a.shape[0], start_index[1] + a.shape[1])

    # 在b中填充a
    b[start_index[0]:end_index[0], start_index[1]:end_index[1]] = a

    return b


# 计算整个簇的重叠率 值代表该网格点由多少个视点所覆盖
def cul_cluster_cover(cluster):
    cluster_cover = np.array([])
    to_remove = []
    for point in cluster['points']:
        if judge_extent_beyond_area(*point) is True:
            # 收集需要删除的元素
            to_remove = [point for point in cluster['points'] if judge_extent_beyond_area(*point)]
            continue
        lon, lat = grid_point_to_lon_lat(*point)
        print(f"当前分析点为{point[0], point[1]}")
        # if point[0] == 2.0 and point[1] == 719.0:
        #     c = 1
        viewshed, size = DPDERL.analysis_by_spderl_simplified(lon, lat, dem.height[int(point[0])][int(point[1])],
                                                              radius, h_stand)
        viewshed = transform_matrix(viewshed, point)
        # 如果cluster_cover是空的，直接赋值
        if cluster_cover.size == 0:
            cluster_cover = np.array(viewshed)
        else:
            # 否则，确保尺寸匹配并累加
            cluster_cover += np.array(viewshed)

    # 实际删除元素 通过list操作删除元素原始列表也会修改
    for item in to_remove:
        cluster['points'].remove(item)

    return cluster_cover


def get_overlap(cur_view, cluster_cover, coverage):
    # 确保是numpy数组操作
    cur_view = np.array(cur_view)
    cluster_cover = np.array(cluster_cover)
    # 选出当前点可视且其可视网格点同时其他点也可视的索引
    mask = (cur_view == 1) & (cluster_cover >= 1)
    # 根据索引获取对应的值
    selected_elements = cluster_cover[mask]
    # 值求和，为重叠率
    overlap = np.sum(selected_elements)
    # 减去自身的可视区域
    return overlap - coverage


# 计算簇的wcom
def cul_cluster_wcoms(cluster, i):
    cluster_cover = cul_cluster_cover(cluster)
    cover_clusters[i] = cluster_cover
    for point in cluster['points']:
        lon, lat = grid_point_to_lon_lat(*point)
        cur_view, size = DPDERL.analysis_by_spderl_simplified(lon, lat, dem.height[int(point[0])][int(point[1])],
                                                              radius, h_stand)
        global_view = transform_matrix(cur_view, point)

        # 找出当前可视数组为1且总视域数组中相应元素也为1的点的坐标，即CC覆盖贡献点
        cc_mask = (global_view == 1) & (cluster_cover == 1)
        cc = np.array(np.where(cc_mask)).T  # 转置以获得[[x1, y1], [x2, y2], ...]形式

        # 找出当前可视数组为1且总视域数组中相应元素大于1的点的坐标，即OC重叠贡献点
        oc_mask = (global_view == 1) & (cluster_cover > 1)
        oc = np.array(np.where(oc_mask)).T  # 转置以获得[[x1, y1], [x2, y2], ...]形式

        # 覆盖点个数(包含重叠点)
        coverage = np.sum(cur_view == 1)

        # 重叠点个数
        overlap = len(oc)
        wcom = cul_wcom(coverage, overlap)
        cluster['wcoms'].append(wcom)
        cluster['cc'].append(cc.tolist())
        cluster['oc'].append(oc.tolist())


def judge_extent_beyond_area(x, y):
    r_grid = radius / dem.rdx
    x_grid_count = int(r_grid / dem.dx)
    y_grid_count = int(r_grid / dem.dy)

    min_x: int = x - (x_grid_count - 1)
    max_x: int = x + x_grid_count
    min_y: int = y - (y_grid_count - 1)
    max_y: int = y + y_grid_count
    if min_x <= 0 or min_y <= 0 or max_x >= dem.x_size or max_y >= dem.y_size:
        return True
    return False


def main():
    global cover_clusters
    # 从txt文件读取字典
    with open('output/csv/our_KM_clusters_data.txt', 'r') as file:
        clusters = json.load(file)

    # 将所有字符串键转换为整数键
    clusters = {int(k): v for k, v in clusters.items()}
    start_time = datetime.now()
    print(f"成功获取聚类后数据-{start_time}")

    # 使用clusters
    print(clusters)

    cover_clusters = np.zeros((len(clusters), *dem.height.shape), dtype=np.uint16)
    # 计算内存占用
    num_elements = np.prod(cover_clusters.shape)  # 总的元素数量
    bytes_per_element = cover_clusters.itemsize  # 每个元素的字节数
    total_bytes = num_elements * bytes_per_element / (1024 * 1024)  # 总的字节数
    total_mb = total_bytes / 1024  # 转换为兆字节（MB）
    print(f"Total memory usage for cover_clusters: {total_mb:.2f} GB")

    # 根据公式计算簇内点的WCOM 并基于WCOM排序
    for i in range(len(clusters)):
        print(f"当前分析簇为：簇{i}")
        if (i == 7):
            a = 2
        cul_cluster_wcoms(clusters[i], i)
        # 对 wcoms 数组进行排序的同时，保持 points 列表中相应点的对应关系
        # 使用浅拷贝
        original_wcoms = np.array(clusters[i]['wcoms'])
        original_points = np.array(clusters[i]['points'])
        original_cc = list(clusters[i]['cc'])
        original_oc = list(clusters[i]['oc'])
        # original_wcoms = np.array(copy.deepcopy(clusters[i]['wcoms']))
        # original_points = np.array(clusters[i]['points'])
        # original_cc = copy.deepcopy(clusters[i]['cc'])
        # original_oc = copy.deepcopy(clusters[i]['oc'])
        clusters[i]['points'] = []
        clusters[i]['cc'] = []
        clusters[i]['oc'] = []
        heapq.heapify(clusters[i]['wcoms'])

        # 确定元素从旧位置到新位置的映射关系 保持wcoms与point的同步
        new_indices = np.argsort(original_wcoms)
        # new_indices = np.array([np.where(original_wcoms == item)[0][0] for item in clusters[i]['wcoms']])
        # 展为一维
        # new_indices = new_indices.flatten()
        # 应用映射关系重排point
        new_points = original_points[new_indices]
        # 列表无法直接索引
        new_oc = [original_oc[index] for index in new_indices]
        new_cc = [original_cc[index] for index in new_indices]
        clusters[i]['points'] = new_points.tolist()
        clusters[i]['oc'] = new_oc
        clusters[i]['cc'] = new_cc

    # # 存储总视域覆盖矩阵
    # # 保存这个三维数组到.npy文件
    # np.save('output/csv/our_cover_clusters.npy', cover_clusters)
    #
    # # 计算WCOM后及以堆结构存储后的结果
    # with open('output/csv/our_heapq_clusters_data.txt', 'w') as file:
    #     json.dump(clusters, file, indent=4)
    print("成功计算WCOM并已堆存储")
    return clusters, cover_clusters


# 如果模块是被自身直接运行的，则main()被运行，如果模块是被导入的，则main()不被运行。
if __name__ == "__main__":
    main()
