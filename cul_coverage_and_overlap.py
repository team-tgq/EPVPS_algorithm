import numpy as np
import rasterio
from affine import Affine
from osgeo import gdal

import DPDERL
from cul_wcoms_sorted_by_heap import grid_point_to_lon_lat, transform_matrix
from services import read_dem

datafilepath = "output/csv/plain_961to_48_96/cvf_filtered_plain_feature_point_0-16-56.886992.txt"
points = np.loadtxt(fname=datafilepath, delimiter=',')

# 固定半径为600 划分区域分60*60 分辨率为30
radius = 200
# 读取地形文件
dem = read_dem.Dem("data/plain_1000.tif", 5)
N, M = dem.height.shape
# 视点距地面高程
h_stand = 2

sum_viewshed = np.array([])

sum_point = N * M


def vis(sum_viewshed):
    sum_viewshed = np.rot90(sum_viewshed)
    # 可视化覆盖数组
    dem_data = gdal.Open("data/plain_1000.tif")
    geo_transform = dem_data.GetGeoTransform()
    geo_list = list(geo_transform)
    crs = dem_data.GetProjection()
    geo_transform = dem_data.GetGeoTransform()
    geo_list = list(geo_transform)
    # 创建 Affine 对象
    transform = Affine(geo_list[1], geo_list[2], geo_list[0],
                       geo_list[4], geo_list[5], geo_list[3])
    # 保存栅格数据为GeoTIFF
    with rasterio.open(
            "output/gis/our_viewed_mountain_210.tif", 'w', driver='GTiff',
            height=sum_viewshed.shape[0], width=sum_viewshed.shape[1],
            count=1, dtype=sum_viewshed.dtype,
            crs=crs, transform=transform
    ) as dst:
        dst.write(sum_viewshed, 1)

    print("可视化完成")


def count_elements_ge(array, threshold):
    """
    统计二维数组中大于等于某个阈值的元素个数。

    参数:
    array (np.array): 输入的二维数组。
    threshold (float): 要比较的阈值。

    返回:
    int: 大于等于阈值的元素个数。
    """
    # 生成布尔数组，其中元素大于等于阈值的位置为True
    greater_equal_mask = array >= threshold

    # 计算True的个数，即大于等于阈值的元素个数
    count = np.sum(greater_equal_mask)

    return count


for point in points:
    lon, lat = grid_point_to_lon_lat(*point)
    print(f"当前分析点为{point[0], point[1]}")
    if (point[0], point[1]) == (801.0, 959.0):
        a = 3
    viewshed, size = DPDERL.analysis_by_spderl_simplified(lon, lat, dem.height[int(point[0])][int(point[1])],
                                                          radius, h_stand)
    viewshed = transform_matrix(viewshed, [int(point[0]), int(point[1])])
    a = np.sum(viewshed)
    # 如果final_cover是空的，直接赋值
    if sum_viewshed.size == 0:
        sum_viewshed = np.array(viewshed)
    else:
        # 否则，确保尺寸匹配并累加
        sum_viewshed += np.array(viewshed)
    b = np.sum(sum_viewshed)

# 覆盖率：总视域数组中网格点值为大于等于1的点个数总和，除以总视域数组点数
coverage = count_elements_ge(sum_viewshed, 1) / sum_point
# 重叠率1：总视域数组中网格点值为大于等于2的点个数总和，除以总视域数组点数
overlap = count_elements_ge(sum_viewshed, 2) / sum_point
# 重叠率2：总视域数组中网格点值为大于等于2的元素值的总和再减去数量总和(覆盖)，除以总视域数组点数
# 找出数组中所有大于等于2的元素
greater_than_two = sum_viewshed[sum_viewshed >= 2]

# 计算这些元素的总和
sum_greater_than_two = np.sum(greater_than_two)

# 最后，计算得到的总和减去这些元素的数量：

# overlap = (sum_greater_than_two - count_elements_ge(sum_viewshed, 2)) / sum_point
print(f"覆盖率:{coverage} 重叠率:{overlap}")

# 可视化
# vis(sum_viewshed)
