import math

import numpy as np
from matplotlib import pyplot as plt

from ECKM import assign_cluster_center
from scipy.spatial import distance, voronoi_plot_2d

from scipy.spatial import Voronoi
import heapq
import json

# 防止中文乱码
plt.rcParams['font.sans-serif'] = [u'SimHei']
plt.rcParams['axes.unicode_minus'] = False
k = 0
pred = []
# 初始化聚类中心变化量
center_shift = np.inf
iter_count = 0
file_name = ""
data_type = ""


# 将数据点与聚类中心可视化
def visualising_points_and_clustering_centres(X, clusters, is_end=False):
    labels = [f'C{i}' for i in range(k)]
    fig, ax = plt.subplots()
    # 使用轴对象绘图
    ax.scatter(X[:, 0], X[:, 1], c=pred, s=10)
    for i in clusters:
        center = clusters[i]['center']
        ax.scatter(center[0], center[1], marker='^', color='red', s=5)
        ax.text(center[0], center[1], labels[i], fontsize=10, ha='right', va='bottom')
    # 显示图形
    # plt.show()
    # 设置坐标轴的比例为相等
    # 可以手动设置坐标轴的范围来确保起始刻度相同
    # min_tick = min(X[:, 0].min(), X[:, 1].min()) - 1
    # max_tick = max(X[:, 0].max(), X[:, 1].max()) + 1
    # ax.set_xlim([0, max_tick])
    # ax.set_ylim([0, max_tick])

    fig.gca().set_aspect(1)

    if iter_count % 10 == 0 or is_end is True:
        # fig.savefig(f"output/image/{data_type}/mountain_feature_clustering_{iter_count}", dpi=300)
        fig.savefig(f"output/image/{data_type}/mountain_feature_clustering_{iter_count}")


# 定义欧几里得距离
def distance_eu(p1, p2):
    return np.sqrt(np.sum((p1 - p2) ** 2))


# 将数据点分配给最近的聚类中心
def assign_clusters(X, clusters):
    global pred
    pred = []
    for idx in range(X.shape[0]):
        dist = []

        curr_x = X[idx]

        # 计算当前数据点与聚类中心的距离
        for i in range(k):
            dis = distance_eu(curr_x, clusters[i]['center'])
            dist.append(dis)

        # 求出当前点与所有簇中心最小距离对应的簇索引
        curr_cluster = np.argmin(dist)
        pred.append(curr_cluster)
        # 将该点加入簇中
        clusters[curr_cluster]['points'].append(curr_x.tolist())
    visualising_points_and_clustering_centres(X, clusters)
    return clusters


# 根据 K-means 聚类中分配点的平均值更新聚类中心
def update_clusters(X, clusters):
    global center_shift
    max_shift = 0
    for i in range(k):
        points = np.array(clusters[i]['points'])
        if points.shape[0] > 0:
            # 对x,y取其均值
            new_center = points.mean(axis=0)
            shift = distance_eu(new_center, clusters[i]['center'])
            clusters[i]['center'] = new_center.tolist()
            max_shift = max(max_shift, shift)
    center_shift = max_shift
    return clusters


# 对最终的簇中心基于空圆的思想确定相邻簇,未使用凸包，因为簇中心构造的维诺圆顶点范围应不会超出
def determine_neighbouring_clusters_by_empty_circle(clusters):
    clusters_center_list = []
    for i in range(k):
        clusters_center_list.append(clusters[i]['center'])
    clusters_center = np.array(clusters_center_list)
    # 构建维诺图
    vor = Voronoi(clusters_center)
    vor_array = np.array(vor.vertices)

    fig = voronoi_plot_2d(vor, show_vertices=False, line_colors='orange', line_width=1,
                          line_alpha=0.8, alpha=0.8, s=10)
    fig.savefig(f"output/image/{data_type}/voronoi_mountain_feature")

    # 计算内部维诺顶点与原始数据点之间的欧式距离 rad_dist[i][j]表示in_vor中第i个点与points第j个点的距离
    rad_dist = (distance.cdist(vor_array, clusters_center, metric='euclidean'))
    radlist = []
    # 计算每个维诺顶点与数据点的最小距离
    for ii in range(len(rad_dist)):
        radlist.append(min(rad_dist[ii]))
    radlist = np.array(radlist)

    eps = 0.001
    nsmall = math.ceil(k / 3.0)

    # 记录相邻簇字典
    neighbour_cluster = {}

    # LECs and SECs cif为确定空圆的数量 以最大距离对应的维诺点作为圆心
    circlefind = list(range(k))
    # 循环的目的在于将初始中心点数量大于等于给定的k值
    for loop in circlefind:
        rad_small = np.array(heapq.nsmallest(nsmall, radlist))

        sm_cen_idx = heapq.nsmallest(nsmall, range(len(radlist)), radlist.take)

        # 获取最大及最小距离对应的维诺点坐标
        sm_cen = vor_array[sm_cen_idx]

        # 具有最小距离维诺顶点与原始数据点的欧式距离
        sm_dst = np.array(distance.cdist(sm_cen, clusters_center, 'euclidean'))

        ncranges = np.arange(0, len(rad_small), 1)

        # ncrangel的长度生成一个新的列表，每个元素都是一个空列表
        sidx = [[] for i in ncranges]
        # 记录相邻簇字典
        neighbour_cluster = {}

        # 将计算得到的差值四舍五入到小数点后三位并于eps判断 返回满足条件的元素的列索引并转为列表
        # dst[i]表示最大距离的维诺点坐标P与其余原始数据集距离列表，rad_lrg[i]表示最大距离值
        # (np.where((np.round([abs(dst[ncount] - rad_lrg[ncount])], 3)) < eps)该表达式是找寻圆心为P半径为rad_lrg[ncount]的圆上圆周点，该圆周点同时也属于原始数据集上的候选点
        for ncount in ncranges:
            sidx[ncount] = ((np.where((np.round([abs(sm_dst[ncount] - rad_small[ncount])], 3)) < eps))[1].tolist())
            for element in sidx[ncount]:
                # 使用列表推导式获取除当前元素外的其他所有元素
                neighbour_elements = [x for x in sidx[ncount] if x != element]
                # 先出现的优先级高
                if element in neighbour_cluster.keys():
                    continue
                # 由于之前初始化 因此在此处索引的顺序即代表簇数
                neighbour_cluster[element] = neighbour_elements

        if len(neighbour_cluster) > k:
            break
        else:
            nsmall += 1
    for i in range(k):
        if i not in neighbour_cluster:
            neighbour_cluster[i] = []
        clusters[i]['neighbour'] = neighbour_cluster[i]
    print(f"优化后相邻簇：{neighbour_cluster}")
    return clusters


def main():
    global k, file_name, data_type
    clusters, X, noc, fn, dt = assign_cluster_center()
    k = noc
    file_name = fn
    data_type = dt

    global iter_count
    max_iter = 1000
    tol = 0.0001  # 阈值
    clusters = assign_clusters(X, clusters)
    while center_shift > tol and iter_count < max_iter:
        for i in range(k):
            clusters[i]['points'] = []
        clusters = assign_clusters(X, clusters)
        clusters = update_clusters(X, clusters)
        iter_count += 1

    # 防止空簇
    for i in range(k):
        if len(clusters[i]['points']) == 0:
            clusters[i]['points'].append(clusters[i]['center'])
    visualising_points_and_clustering_centres(X, clusters, True)
    clusters = determine_neighbouring_clusters_by_empty_circle(clusters)
    print(f"center_shift:{center_shift},tol:{tol}\niter_count:{iter_count},max_iter:{max_iter}")
    # 将字典保存到txt文件
    with open('output/csv/our_KM_clusters_data.txt', 'w') as file:
        json.dump(clusters, file, indent=4)
    print(f"成功保存{k}个聚类后数据")


# 如果模块是被自身直接运行的，则main()被运行，如果模块是被导入的，则main()不被运行。
if __name__ == "__main__":
    main()
