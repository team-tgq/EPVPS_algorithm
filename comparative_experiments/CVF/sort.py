import numpy as np

import DPDERL
from cul_wcoms_sorted_by_heap import grid_point_to_lon_lat, judge_extent_beyond_area, get_overlap, transform_matrix
from services import read_dem

# 固定半径为600 划分区域分60*60 分辨率为30
radius = 200
# 读取地形文件
demdata = read_dem.Dem("data/plain_1000.tif", 5)
# 视点距地面高程
h_stand = 2

d = 0


class Sort:
    @staticmethod
    def descend_sort(clusters, demdata):
        # 根据给定的评判指数升序排序簇
        for i in range(len(clusters)):
            global d
            d = i
            indices_to_remove = []
            print(f"当前簇为:{i}")

            # 计算视点信息
            for j in range(2, len(clusters[i]), 2):
                if judge_extent_beyond_area(clusters[i][j], clusters[i][j + 1]):
                    indices_to_remove.append(j)
                    continue

            # 删除视点
            for m in sorted(indices_to_remove, reverse=True):
                del clusters[i][m:m + 2]  # 删除点及其下一个元素

            # 重新计算视点信息
            viewshed_num = [0] * (len(clusters[i]) // 2 - 1)
            overlap_num = [0] * (len(clusters[i]) // 2 - 1)
            fitness = [0.0] * (len(clusters[i]) // 2 - 1)

            all_cover = Sort.cul_all_cover(clusters[i])

            print(f"计算簇{i}的IVC")
            for j in range(2, len(clusters[i]), 2):
                print(f"计算簇{i}的第{(j - 2) // 2}个点IVC")
                a = int(clusters[i][j])
                b = int(clusters[i][j + 1])
                lon, lat = grid_point_to_lon_lat(a, b)
                # 计算视点的可视性指数
                current_viewshed, _ = DPDERL.analysis_by_spderl_simplified(lon, lat, demdata.height[int(a)][int(b)],
                                                                           radius, h_stand)
                current_viewshed = transform_matrix(current_viewshed, [a, b])
                coverage = int(np.sum(current_viewshed))
                current_overlap = get_overlap(current_viewshed, all_cover, coverage)
                viewshed_num[(j - 2) // 2] = coverage
                overlap_num[(j - 2) // 2] = current_overlap

            for k in range(2, len(clusters[i]), 2):
                fitness[(k - 2) // 2] = overlap_num[(k - 2) // 2] / (viewshed_num[(k - 2) // 2] ** 2)

            print(f"按IVC升序排序簇{i}")
            for s in range(len(fitness) - 1):
                for t in range(len(fitness) - 1 - s):
                    if fitness[t] > fitness[t + 1]:
                        Sort.swap(clusters[i], fitness, t, t + 1)

        return clusters

    @staticmethod
    def swap(list_, fitness, s, t):
        # 交换列表中两个元素的位置
        temp = fitness[s]
        fitness[s] = fitness[t]
        fitness[t] = temp

        row = list_[s * 2 + 2]
        col = list_[s * 2 + 3]

        list_[s * 2 + 2] = list_[t * 2 + 2]
        list_[s * 2 + 3] = list_[t * 2 + 3]

        list_[t * 2 + 2] = row
        list_[t * 2 + 3] = col

    @staticmethod
    def compute_o(i, j, cluster):
        # 计算可视域重叠
        lon, lat = grid_point_to_lon_lat(i, j)
        # 计算视点的可视性指数
        view, _ = DPDERL.analysis_by_spderl_simplified(lon, lat, demdata.height[i][j],
                                                       radius, h_stand)
        current_view = transform_matrix(view, [i, j])
        count = 0

        for k in range(2, len(cluster), 2):
            # print(f"当前簇为:{d}")
            a = cluster[k]
            b = cluster[k + 1]
            if a != i or b != j:
                lon, lat = grid_point_to_lon_lat(a, b)
                if judge_extent_beyond_area(a, b) is True:
                    continue
                # 计算视点的可视性指数
                temp_view, _ = DPDERL.analysis_by_spderl_simplified(lon, lat, demdata.height[int(a)][int(b)],
                                                                    radius, h_stand)
                temp_view = transform_matrix(temp_view, [a, b])
                for x in range(len(current_view)):
                    for y in range(len(current_view[0])):
                        if temp_view[x][y] == 1 and current_view[x][y] == 1:
                            count += 1

        return count

    @staticmethod
    def compute_v(i, j):
        lon, lat = grid_point_to_lon_lat(i, j)

        # 计算视点的可视性指数
        view, _ = DPDERL.analysis_by_spderl_simplified(lon, lat, demdata.height[int(i)][int(j)],
                                                       radius, h_stand)
        count = 0

        for a in range(len(view)):
            for b in range(len(view[0])):
                if view[a][b] == 1:
                    count += 1

        return count

    @staticmethod
    def cul_all_cover(cluster):
        cluster_cover = np.zeros(demdata.height.shape)
        for k in range(2, len(cluster), 2):
            # print(f"当前簇为:{d}")
            a = cluster[k]
            b = cluster[k + 1]
            lon, lat = grid_point_to_lon_lat(a, b)
            # 计算视点的可视性指数
            viewshed, _ = DPDERL.analysis_by_spderl_simplified(lon, lat, demdata.height[int(a)][int(b)],
                                                               radius, h_stand)
            viewshed = transform_matrix(viewshed, [a, b])
            # 如果cluster_cover是空的，直接赋值
            if cluster_cover.size == 0:
                cluster_cover = np.array(viewshed)
            else:
                # 否则，确保尺寸匹配并累加
                cluster_cover += np.array(viewshed)
        return cluster_cover
