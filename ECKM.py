import math
import os

import numpy as np
import sys
import heapq

from matplotlib import pyplot as plt
from scipy.spatial import distance, voronoi_plot_2d

from scipy.spatial import ConvexHull
from scipy.spatial import Voronoi

# 簇数量
noc = 45

plt.rcParams.update({'figure.max_open_warning': 0})
# datafilepath = "data/compound_data.txt"
# read_data = np.loadtxt(fname=datafilepath)
datafilepath = "output/csv/our_mountain_feature.txt"
read_data = np.loadtxt(fname=datafilepath, delimiter=',')
# 使用os.path.basename获取路径中的最后一部分（文件名）
full_file_name = os.path.basename(datafilepath)
# 使用os.path.splitext去除文件名的扩展名
file_name, _ = os.path.splitext(full_file_name)
data_type = ""
if file_name == "our_mountain_feature":
    data_type = "my_data"
else:
    data_type = "test_data"
main_data = np.array(read_data)
fig, ax = plt.subplots()
ax.scatter(main_data[:, 0], main_data[:, 1], alpha=0.8, s=10)
ax.set(xlabel='X', ylabel='Y', title='Dataset')
fig.savefig(f"output/image/{data_type}/data_{file_name}")
points = main_data
# 构建维诺图
vor = Voronoi(points)
# 构架凸包，使得分析区域在凸包内，返回值为凸包顶点的索引
hull = ConvexHull(points)
# 获取维诺图的维诺顶点，维诺顶点可能在凸包的外部
myList1 = vor.vertices
# 获取凸包顶点
myList = points[hull.vertices]
plt.plot(points[:, 0], points[:, 1], 'o', markersize=3)
for simplex in hull.simplices:
    plt.plot(points[simplex, 0], points[simplex, 1], 'k-')
plt.savefig(f'output/image/{data_type}/convex_hull_{file_name}.png')  # 保存图像为PNG文件

nat = points.shape[1]
insideVor = []
b = myList

# 遍历维诺顶点并添加至凸包中，判断凸包结构是否变化，若发生变化则表明生成的维诺顶点在凸包外部则舍去该点
for i in range(0, len(myList1)):
    x = (np.append(b, myList1[i]))
    x = np.array(np.reshape(np.array([x]), (-1, nat)))
    hull2 = ConvexHull(x)
    # oldHull = np.array(points[hull.vertices])
    # newHull = np.array(x[hull2.vertices])
    # if np.array_equal(oldHull, newHull):
    #     insideVor.append(myList1[i])
    # 顺序会乱但值相同 用set更合适
    oldHull = set(tuple(point) for point in points[hull.vertices])
    newHull = set(tuple(point) for point in x[hull2.vertices])
    if oldHull == newHull:
        insideVor.append(myList1[i])
    else:
        continue
# 内部点
in_vor = np.array(insideVor)

fig = voronoi_plot_2d(vor, show_vertices=False, line_colors='orange', line_width=1,
                      line_alpha=0.8, alpha=0.8, s=10)
fig.savefig(f"output/image/{data_type}/voronoi_{file_name}")

# 外部点
out_vor = np.array([elem for elem in myList1 if elem not in in_vor])

# 计算内部维诺顶点与原始数据点之间的欧式距离 rad_dist[i][j]表示in_vor中第i个点与points第j个点的距离
rad_dist = (distance.cdist(in_vor, points, metric='euclidean'))
radlist = []
# 计算每个维诺顶点与数据点的最小距离
for ii in range(len(rad_dist)):
    radlist.append(min(rad_dist[ii]))
radlist = np.array(radlist)


cfi = noc
eps = 0.001
nlarge = math.ceil(noc / 3.0)
nsmall = math.ceil(noc / 3.0)

# 记录相邻簇字典
neighbour_cluster = {}

# LECs and SECs cif为确定空圆的数量 以最大距离对应的维诺点作为圆心
circlefind = list(range(1, cfi))
lar_cen = []
res = []
rad_lrg = []
# 循环的目的在于将初始中心点数量大于等于给定的k值
for loop in circlefind:
    # 通过堆排 返回nlarge个最大最小距离(圆半径)的元素 返回nsmall个最小的元素
    rad_lrg = np.array(heapq.nlargest(nlarge, radlist))

    # 包含radlist中最大最小距离的nlarge个元素的索引 numpy.take() 函数沿索引从数组中返回元素
    cen_idx = heapq.nlargest(nlarge, range(len(radlist)), radlist.take)
    sm_cen_idx = heapq.nsmallest(nsmall, range(len(radlist)), radlist.take)

    # 获取最大及最小距离对应的维诺点坐标
    lar_cen = in_vor[cen_idx]
    sm_cen = in_vor[sm_cen_idx]

    # 具有最大距离、最小距离维诺顶点与原始数据点的欧式距离
    dst = np.array(distance.cdist(lar_cen, points, 'euclidean'))
    sm_dst = np.array(distance.cdist(sm_cen, points, 'euclidean'))

    ncrangel = np.arange(0, len(rad_lrg), 1)

    # ncrangel的长度生成一个新的列表，每个元素都是一个空列表
    idx = [[] for i in ncrangel]
    # 记录相邻簇字典
    neighbour_cluster = {}

    # 按照圆周点对应圆的半径为第一比较关系，圆周点对应索引未第二比较关系确定簇。即第一大的圆对应索引:a<b<c圆周点，则a为簇0
    # 将计算得到的差值四舍五入到小数点后三位并于eps判断 返回满足条件的元素的列索引并转为列表
    # dst[i]表示最大距离的维诺点坐标P与其余原始数据集距离列表，rad_lrg[i]表示最大距离值
    # (np.where((np.round([abs(dst[ncount] - rad_lrg[ncount])], 3)) < eps)该表达式是找寻圆心为P半径为rad_lrg[ncount]的圆上圆周点，该圆周点同时也属于原始数据集上的候选点
    for ncount in ncrangel:
        idx[ncount] = ((np.where((np.round([abs(dst[ncount] - rad_lrg[ncount])], 3)) < eps))[1].tolist())
        # 创建一个空集合用于跟踪已经见过的元素
        seen = set()
        # 使用列表推导式进行去重，同时保持原有顺序
        flattened_list = [item for sublist in idx for item in sublist if not (item in seen or seen.add(item))]

        for element in idx[ncount]:
            # 使用列表推导式获取除当前元素外的其他所有元素
            neighbour_elements = [x for x in idx[ncount] if x != element]
            neighbour_indexes = [flattened_list.index(item) for item in neighbour_elements if item in flattened_list]
            # 使用列表推导式去除所有大于30的元素
            filtered_neighbour_indexes = [index for index in neighbour_indexes if index < noc]
            neighbour_cluster[flattened_list.index(element)] = filtered_neighbour_indexes

    # Indexing Part
    all_idx = []
    sm_all_idx = []
    xidx = [i for i in idx if i != []]
    for ncount in ncrangel:
        all_idx = all_idx + xidx[ncount]

    # Feasibility Checking
    res = []
    for i in all_idx:
        if i not in res:
            res.append(i)
    sm_res = []
    for i in sm_all_idx:
        if i not in sm_res:
            sm_res.append(i)

    if len(res) > noc:
        break
    else:
        nlarge += 1

if len(res) < noc:
    sys.exit("Sorry, Error is there, increase the the length of circlefind list")

print(res)
circumpts_all = points[res]

# LEC and SEC Method (M1)
strt_l_m1 = np.array(circumpts_all[:noc])

k5_cen = circumpts_all
c5 = lar_cen
vorxx = c5[:, 0]
voryy = c5[:, 1]
r5 = np.array(rad_lrg)
fig, ax = plt.subplots()
for simplex in hull.simplices:
    plt.plot(points[simplex, 0], points[simplex, 1], 'k-')
plt.scatter(k5_cen[:, 0], k5_cen[:, 1], marker="x", color='r', s=100, alpha=1)
plt.scatter(points[:, 0], points[:, 1], s=10)
for i in range(len(c5)):
    circle5 = plt.Circle((vorxx[i], voryy[i]), r5[i], fill=False, color='g')
    ax.add_artist(circle5)
voronoi_plot_2d(vor, ax, show_vertices=False, line_colors='orange', line_width=1,
                line_alpha=0.8, show_points=False)
fig.savefig(f'output/image/{data_type}/circle_{file_name}.png')

# 算法名称 ECHVF
# 定义簇结构
# Weighted Coverage Overlap Metric:加权覆盖重叠度量 F=w1×Coverage−w2×Overlap
"""
neighbour:相邻簇
center:簇中心
points:簇包含的点
wcoms:加权覆盖重叠度量
cover:簇总视域矩阵（数据过大，另外用一个文本存储）
Coverage contribution:覆盖贡献点
overlap contribution:重叠贡献点
"""


def assign_cluster_center():
    clusters = {}
    for idx in range(noc):
        cluster = {
            'neighbour': neighbour_cluster[idx],
            'center': circumpts_all[idx].tolist(),
            'points': [],
            'wcoms': [],
            'cc': [],
            'oc': []
        }
        clusters[idx] = cluster
    print(f"初始相邻簇：{neighbour_cluster}")
    return clusters, main_data, noc, file_name, data_type
