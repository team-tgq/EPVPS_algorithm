import numpy as np

from services import read_dem


# 假设这是你的DEM高程数组
dem_data = read_dem.Dem("../../data/plain_1000.tif", 5)
M, N = dem_data.height.shape
VR = np.zeros((M, N), dtype=int)
#分析半径对应的网格数
grid_num = 40
# 与Ours特征点数保持一致
size_m = 30
size_n = 30
feature_points = []


# 寻找并返回指定子区块中高程值最高的点的坐标。
def compute(i, j, s_m, s_n):
    max_value = -float('inf')
    max_x, max_y = -1, -1

    # 确保索引不越界
    end_i = min(i + s_m, M - grid_num - 1)
    end_j = min(j + s_n, N - grid_num - 1)

    # 寻找最大高程值的坐标
    for x in range(i, end_i):
        for y in range(j, end_j):
            if dem_data.height[x][y] > max_value:
                max_value = dem_data.height[x][y]
                max_x = x
                max_y = y

    # 返回最高点的坐标
    return [max_x, max_y]


# 这个函数的目的是在二维数组 dem_data.height 的一个特定子区块内找到具有最高高程值且对应的 VR 值不为零的点的坐标，并将其存储在列表中
def partition_find(i, j, s_m, s_n):
    current_list = []
    while len(current_list) < 1:
        coords = compute(i, j, s_m, s_n)
        if not current_list:
            current_list.append(coords)
            VR[coords[0]][coords[1]] = 0
    return current_list


"""
局部极大点（高峰）和局部极小点（低谷）被标记为 3 和 -2，分别表示地形的高点和低点。
特定的视线相关性标记为 1 或 -1，表示在东西方向上的视线通过性与阻碍。
数组 VR 用于存储这些特征值，可以用于进一步的分析，如可视性分析、路径规划等。
"""


def feature_point():
    # 遍历数组，跳过边缘以便访问所有邻居
    for i in range(1, dem_data.height.shape[0] - 1):
        for j in range(1, dem_data.height.shape[1] - 1):
            center = dem_data.height[i, j]
            N = dem_data.height[i - 1, j]
            S = dem_data.height[i + 1, j]
            W = dem_data.height[i, j - 1]
            E = dem_data.height[i, j + 1]
            NE = dem_data.height[i - 1, j + 1]
            SE = dem_data.height[i + 1, j + 1]
            WN = dem_data.height[i - 1, j - 1]
            WS = dem_data.height[i + 1, j - 1]

            if (W - center) * (E - center) > 0:
                if E > center:
                    VR[i, j] = -1
                if E < center:
                    VR[i, j] = 1
            elif (S - center) * (N - center) > 0:
                if N > center:
                    VR[i, j] = -1
                if N < center:
                    VR[i, j] = 1

            # 检测局部极大点
            if center > W and center > N and center > E and center > S and center > NE and center > SE and center > WN and center > WS:
                VR[i, j] = 3

            # 检测局部极小点
            if center < W and center < N and center < E and center < S and center < NE and center < SE and center < WN and center < WS:
                VR[i, j] = -2


# 与边界区域留出分析半径的距离
def extract():
    global feature_points
    feature_point()
    for i in range(1, (M // size_m)):
        for j in range(1, (N // size_n)):
            #  if i * size < 40 or i * size >= 960 or j * size < 40 or j * size >= 960:
            if i * size_m < grid_num or i * size_m >= M - grid_num or j * size_n < grid_num or j * size_n >= N - grid_num:
                continue
            block_list = partition_find(i * size_m, j * size_n, size_m, size_n)
            feature_points.extend(block_list)
    return feature_points


# 主函数调用
features = extract()
# 指定要写入的文件路径
file_path = "../../output/csv/cvf_mountain_feature.txt"

# 将数组写入文件
with open(file_path, 'w') as f:
    for row in feature_points:
        # 将每行转换为以逗号分隔的字符串，并写入文件
        line = ','.join(map(str, row))
        f.write(line + '\n')
