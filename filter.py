import copy
import gc
import heapq
import json
import math
import time
from datetime import datetime

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import DPDERL
import cul_wcoms_sorted_by_heap
from cul_wcoms_sorted_by_heap import get_overlap, cul_wcom, transform_matrix, grid_point_to_lon_lat
from preprocessing.extract_filiter_point import ours_extract_point
from services import read_dem

# 固定半径为200 划分区域分40*40 分辨率为5
radius = 200
# 读取地形文件
dem = read_dem.Dem("data/plain_1000.tif", 5)
# 视点距地面高程
h_stand = 2

# 当前分析簇数
d = 0

clusters, cover_clusters = cul_wcoms_sorted_by_heap.main()


# # 读取总视域数组
# cover_clusters = np.load('output/csv/our_cover_clusters.npy')
#
#
# print("开始读取簇")
# # 从txt文件读取字典
# with open('output/csv/our_heapq_clusters_data.txt', 'r') as file:
#     clusters = json.load(file)
#
# print("文件读取成功")
# # 将所有字符串键转换为整数键
# clusters = {int(k): v for k, v in clusters.items()}
#
# # 遍历字典，并将每个 'wcoms' 列表转换为堆
# for key, value in clusters.items():
#     if 'wcoms' in value and isinstance(value['wcoms'], list):
#         heapq.heapify(value['wcoms'])


# 是否结束过滤，若每个簇仅剩余一个点则结束过滤
def judge_end_filter(clusters):
    is_end_filter = False
    for i in range(len(clusters)):
        if len(clusters[i]['points']) > 1:
            return is_end_filter
    is_end_filter = True
    return is_end_filter


# 更新簇内的wcoms和cover实现细节
def update_detail(cluster):
    for i in range(len(cluster['points'])):
        point = cluster['points'][i]
        lon, lat = grid_point_to_lon_lat(*point)
        cur_view, size = DPDERL.analysis_by_spderl_simplified(lon, lat, dem.height[int(point[0])][int(point[1])],
                                                              radius, h_stand)
        global_view = transform_matrix(cur_view, cluster['points'][i])
        # 覆盖率
        coverage = np.sum(cur_view == 1)
        # 重叠率
        overlap = get_overlap(global_view, cluster['cover'], coverage)
        wcom = cul_wcom(coverage, overlap)
        cluster['wcoms'].append(wcom)


# 某簇删除或添加点需要重新更新该簇的wcom
# cul_cover操作点的可视矩阵，cluster为点所对应的操作簇，opt为0代表移除操作，1代表添加操作
def update_wcom(cur_cover, cluster, oc, opt):
    cur_cover = np.array(cur_cover)
    # 移除该点
    # 如果移除某点导致总视域矩阵变为1则这些点应从重叠贡献点转为覆盖贡献点
    if opt == 0:
        cover_clusters[d] -= cur_cover
        # 如果不存在重叠贡献点则不对簇内其他店产生影响
        if len(oc) != 0:
            cluster['wcoms'] = []
            rows, cols = zip(*oc)  # 解压坐标列表到两个元组

            # 使用解压后的行和列作为索引来获取对应位置的值
            values = cover_clusters[d][rows, cols]

            # 查找值为 1 的位置
            indices_of_ones = tuple(tuple(oc[i]) for i, value in enumerate(values) if value == 1)

            for i in range(len(cluster['oc'])):
                # 转换每个内部列表为元组
                tuple_of_cc = tuple(tuple(x) for x in cluster['cc'][i])
                tuple_of_oc = tuple(tuple(x) for x in cluster['oc'][i])

                # 转换为集合并使用并差集操作
                cc_set = set(tuple_of_cc)
                oc_set = set(tuple_of_oc)
                indices_set = set(indices_of_ones)

                # 求出重叠贡献点与弹出点的存在交集
                insec = oc_set.intersection(indices_set)

                # 这些点从重叠贡献点集中删除
                oc_set.difference_update(insec)

                # 转为覆盖贡献点
                cc_set.update(insec)

                list_of_cc = [list(inner_tuple) for inner_tuple in tuple_of_cc]
                list_of_oc = [list(inner_tuple) for inner_tuple in tuple_of_oc]

                # 转换回列表
                cluster['oc'][i] = list_of_oc
                cluster['cc'][i] = list_of_cc

                # 覆盖点个数(包含重叠点)
                coverage = len(cluster['oc'][i]) + len(cluster['cc'][i])
                # 重叠点个数
                overlap = len(cluster['oc'][i])
                wcom = cul_wcom(coverage, overlap)
                cluster['wcoms'].append(wcom)

    # 插入操作
    # 如果插入点导致总视域矩阵从0到1 则点为覆盖贡献点 从1到2则这些点为重叠贡献点
    elif opt == 1 and len(oc) != 0:
        cover_clusters[d] += cur_cover
        if len(oc) > 0:
            cluster['wcoms'] = []
            # 插入操作该点，由于其与原先簇脱离则覆盖、重叠失去意义，全视为可见点
            # 转换坐标列表为 NumPy 索引
            rows, cols = zip(*oc)  # 解压坐标列表到两个元组

            # 使用解压后的行和列作为索引来获取对应位置的值
            values = cover_clusters[d][rows, cols]

            # 值变为重叠贡献点的位置
            indices_of_ones = tuple(tuple(oc[i]) for i, value in enumerate(values) if value == 2)

            for i in range(len(cluster['oc'])):
                # 转换每个内部列表为元组
                tuple_of_cc = tuple(tuple(x) for x in cluster['cc'][i])
                tuple_of_oc = tuple(tuple(x) for x in cluster['oc'][i])

                # 转换为集合并使用并差集操作
                cc_set = set(tuple_of_cc)
                oc_set = set(tuple_of_oc)
                indices_set = set(indices_of_ones)

                # 求出覆盖贡献点与弹出点的存在交集
                insec = cc_set.intersection(indices_set)

                # 这些点从覆盖贡献点集中删除
                cc_set.difference_update(insec)

                # 转为重叠贡献点
                oc_set.update(insec)

                list_of_cc = [list(inner_tuple) for inner_tuple in tuple_of_cc]
                list_of_oc = [list(inner_tuple) for inner_tuple in tuple_of_oc]

                # 转换回列表
                cluster['oc'][i] = list_of_oc
                cluster['cc'][i] = list_of_cc

                # 覆盖点个数(包含重叠点)
                coverage = len(cluster['oc'][i]) + len(cluster['cc'][i])

                # 重叠点个数
                overlap = len(cluster['oc'][i])
                wcom = cul_wcom(coverage, overlap)
                cluster['wcoms'].append(wcom)

    # 对 wcoms 数组进行排序的同时，保持 points 列表中相应点的对应关系
    # 浅拷贝
    original_wcoms = np.array(cluster['wcoms'])
    original_points = np.array(cluster['points'])
    original_oc = list(cluster['oc'])
    original_cc = list(cluster['cc'])
    # original_wcoms = np.array(copy.deepcopy(cluster['wcoms']))
    # original_points = np.array(cluster['points'])
    # original_oc = copy.deepcopy(cluster['oc'])
    # original_cc = copy.deepcopy(cluster['cc'])

    cluster['points'] = []
    cluster['oc'] = []
    cluster['cc'] = []
    heapq.heapify(cluster['wcoms'])

    # 确定元素从旧位置到新位置的映射关系 保持wcoms与point的同步
    new_indices = np.argsort(original_wcoms)
    # new_indices = np.array([np.where(original_wcoms == item)[0][0] for item in cluster['wcoms']])
    #
    # # 展为一维
    # new_indices = new_indices.flatten()

    # 应用映射关系重排point
    new_points = original_points[new_indices]
    new_oc = [original_oc[index] for index in new_indices]
    new_cc = [original_cc[index] for index in new_indices]
    cluster['points'] = new_points.tolist()
    cluster['oc'] = new_oc
    cluster['cc'] = new_cc

    # 删除不再需要的深拷贝数据
    # del original_wcoms
    # del original_points
    # del original_cc
    # del original_oc

    # 手动触发垃圾回收 会增加内存消耗，大幅度降低运行时间
    # gc.collect()


# 符合条件的相邻簇
k = 0
# 优于前%20 才视为更有价值点
prop = 0.3

# 记录循环次数
cnt = 0
# 过滤算法
'''
目的是从三个簇找出视域贡献最低的点，但是只对簇内进行了排序，簇间比较后还需要才簇间重新计算WCOM值并排序
优化目标：在每个簇中选择一个唯一的元素，这个元素在本簇内的值是最小的，并且在任何相邻簇中，这个元素的值都不是最大的，元素的值会被相邻簇影响
将相邻簇对应的视域贡献程度，用堆来存储，每次弹出当前簇最小的值元素，计算该元素在相邻簇的值，并判断是否大于前n(设置阈值)个，若大于则插入堆中，并在该堆继续，小于则继续弹出在当前堆最小值依次循环
'''
last_time = datetime.now()

results = []  # 用于记录实验结果

while not judge_end_filter(clusters):
    for tmp in range(len(clusters)):
        start_time = time.time()
        date = datetime.now()
        num_point = len(clusters[d]['points'])
        print(f"#{cnt}: {date}:当前分析簇为{d},剩余点数为{num_point}")
        cnt += 1
        # 只有当簇点数大于1时才需要过滤
        if len(clusters[d]['points']) > 1:
            # 移除并弹出最小WCOM值和对应的点
            cur_wcom = heapq.heappop(clusters[d]['wcoms'])
            cur_point = clusters[d]['points'].pop(0)
            cur_cc = clusters[d]['cc'].pop(0)
            cur_oc = clusters[d]['oc'].pop(0)
            lon, lat = grid_point_to_lon_lat(*cur_point)
            cur_view, size = DPDERL.analysis_by_spderl_simplified(lon, lat,
                                                                  dem.height[int(cur_point[0])][int(cur_point[1])],
                                                                  radius, h_stand)
            global_view = transform_matrix(cur_view, cur_point)
            if len(cur_oc) == 0:
                sad = 1
            update_wcom(global_view, clusters[d], cur_oc, 0)

            # 与相邻簇最小值作比较，若小于所有则进行弹出当前簇最小值，若大于某一簇的最小值，则插入至该簇
            for i in clusters[d]['neighbour']:
                if cur_wcom > clusters[i]['wcoms'][0]:
                    k = i

                    cluster_cover = cover_clusters[i]
                    tmp_cluster_cover = cluster_cover + global_view

                    # 覆盖率
                    coverage = np.sum(cur_view == 1)

                    # 找出当前可视数组为1且总视域数组中相应元素也为1的点的坐标，即CC覆盖贡献点
                    cc_mask = (global_view == 1) & (tmp_cluster_cover == 1)
                    tmp_cc = np.array(np.where(cc_mask)).T  # 转置以获得[[x1, y1], [x2, y2], ...]形式

                    # 找出当前可视数组为1且总视域数组中相应元素大于1的点的坐标，即OC重叠贡献点
                    oc_mask = (global_view == 1) & (tmp_cluster_cover > 1)
                    tmp_oc = np.array(np.where(oc_mask)).T  # 转置以获得[[x1, y1], [x2, y2], ...]形式

                    # 重叠点个数
                    overlap = len(tmp_oc)
                    # 如果重叠率为0 无法视为相邻簇
                    if overlap == 0:
                        continue

                    # 相邻簇的wcom
                    wcom = cul_wcom(coverage, overlap)

                    # 根据求出的wcom判断是否插入相邻簇，设置阈值当优于一定的百分比时才视其更有价值，并插入到簇中
                    # 从堆中取wcoms前n大的值，向上取整,防止n为0的情况
                    n = math.ceil(len(clusters[k]['wcoms']) * prop)
                    nth_large = heapq.nlargest(n, clusters[k]['wcoms']).pop()
                    if wcom > nth_large:
                        # 将元素添加到堆
                        heapq.heappush(clusters[k]['wcoms'], wcom)
                        # 获取该元素位于堆中的索引
                        index = clusters[k]['wcoms'].index(wcom)
                        # 在同样位置插入对应点和可视性数组
                        clusters[k]['points'].insert(index, cur_point)
                        clusters[k]['oc'].insert(index, tmp_oc.tolist())
                        clusters[k]['cc'].insert(index, tmp_cc.tolist())
                        # 更新簇
                        update_wcom(global_view, clusters[k], tmp_oc, 1)
                        d = k
                        break
        else:
            d = (d + 1) % len(clusters)
        # 记录结果迭代次数
        # end_time = time.time()
        # results.append({
        #     'iteration': cnt,
        #     'computation_time': end_time - start_time,
        # })
now_time = datetime.now()
our_time = str(now_time - last_time)
print(f"本发明算法用时：{our_time}")
safe_file_name = our_time.replace(':', '-')
input_file_path = f'output/csv/our_filtered_clusters_data_{safe_file_name}.txt'
output_file_path = f'output/csv/our_filtered_clusters_data_point_{safe_file_name}.txt'
# 过滤后的结果
with open(input_file_path, 'w') as file:
    json.dump(clusters, file, indent=4)
print("过滤完成")

ours_extract_point(input_file_path, output_file_path)
print(f"保留点：{output_file_path}")

# # 生成统计图
# iterations = [result['iteration'] for result in results]
# computation_times = [result['computation_time'] for result in results]
#w
# plt.figure(figsize=(10, 5))
#
# # 计算时间图
# plt.plot(iterations, computation_times, marker='o')
# plt.title('Computation Time per Iteration')
# plt.xlabel('Iteration')
# plt.ylabel('Computation Time (seconds)')
# plt.grid(True)
#
# plt.tight_layout()
# plt.show()
#
# # 保存实验结果到Excel文件
# results_df = pd.DataFrame(results)
# results_df.to_excel('output/csv/experiment_results.xlsx', index=False)
# print("实验结果已保存至 'output/csv/experiment_results.xlsx'")
