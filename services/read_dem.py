from osgeo import gdal
import numpy as np


class Dem(object):
    def __init__(self, file_path, real_distance):
        # DEM文件路径
        self.file_path = file_path
        # DEM分辨率
        self.real_distance = real_distance
        # 左下角原点的X值
        self.start_x = 0.0
        # 左下角原点的Y值
        self.start_y = 0.0
        # 横间隔
        self.dx = 0
        # 横间隔实地距离
        self.rdx = 1.0
        # 纵间隔
        self.dy = 0
        # 纵间隔实际距离
        self.rdy = 1.0
        # dem行列值
        self.x_size = 0
        self.y_size = 0

        # 打开图像
        dataset = gdal.Open(file_path, gdal.GA_Update)

        # 读入DEM图像行和宽
        self.x_size = dataset.RasterXSize
        self.y_size = dataset.RasterYSize

        # 存储的高程列表
        # self.height[X][Y]的第一位（即Height[0][0]）代表区域左下角的高程,X,Y代表格网点序号
        # 在北半球经度越大,X越大,纬度越大,Y越大
        self.height = np.zeros((self.y_size, self.x_size), dtype=np.float64)

        # 获取栅格数据的点的数量
        count = dataset.RasterCount

        # 获取第一个band
        demband = dataset.GetRasterBand(1)

        # 获取屏幕坐标转换到实际地理坐标的参数 包含横纵网格间距

        # 输出地理转换参数
        # "左上角 x 坐标:", geotransform[0]
        # "水平分辨率:", geotransform[1]
        # "旋转参数:", geotransform[2]
        # "左上角 y 坐标:", geotransform[3]
        # "旋转参数:", geotransform[4]
        # "垂直分辨率:", geotransform[5]
        gt = dataset.GetGeoTransform()

        # 读取DEM数据集的栅格值 以左上角为原点 在使用band.ReadAsArray函数读取高程数据时,返回的NumPy数组的维度顺序是 (y, x),而不是 (x, y).这意味着数组的第一个维度是行数（y）,第二个维度是列数（x）.
        dem_array = demband.ReadAsArray()
        # 以左下角为原点的方式存储 DEM 数组
        dem_array = dem_array[::-1, :].T
        # 正值表示像素大小与地理坐标的方向一致,负值表示相反方向,一致则翻转
        if gt[5] > 0:
            dem_array = dem_array.flipud(dem_array)
            self.start_y = gt[3]
            self.dy = gt[5]
        else:
            self.start_y = gt[3] + self.y_size * gt[5]
            self.dy = -gt[5]

        # 异常值处理
        dem_array[dem_array < 0] = 0
        dem_array[dem_array > 20000] = 20000

        self.height = dem_array

        self.start_x = gt[0]
        self.dx = gt[1]

        if abs(gt[1]) < 1:
            # 108000
            self.rdx = real_distance / abs(self.dx)
            self.rdy = real_distance / abs(self.dy)
        self.max_lon = self.start_x + self.x_size * self.dx
        self.max_lat = self.start_y + self.y_size * self.dy

        # 提取的3600个特征点
        self.feature_point = []

        # 分辨率为5m 1000*1000
        # 3250 ->16
        # 2116 ->20
        # 900  ->30

        # 分辨率为5m 3600*3600
        # 3250 ->56
        # 2116 ->69
        # 900  ->105

        # 分辨率为1 3000,3000
        # 3250 ->56
        # 2116 ->68
        # 900  ->87
        self.feature_point_extraction(87, 87)

    # 获取观察者网格位置

    """ 
    x_grid_center:观察者处与x网格点位置
    y_grid_center:观察者处与y网格点位置
    y_grid_inner:实际位置与x格网点差距
    y_grid_inner:实际位置与y格网点差距
    """

    def get_point_location(self, lon, lat):
        # 判断观察点是否位于格网点上
        is_in_x_grid_point = False
        is_in_y_grid_point = False
        # 不一定刚好落在格网点上
        x_grid_observe = (lon - self.start_x) / self.dx
        y_grid_observe = (lat - self.start_y) / self.dy

        # 最邻近格网点
        # 中心格网点对应的x,y坐标 是要小于观察者点对应的x,y坐标
        x_grid_center = int(x_grid_observe)
        y_grid_center = int(y_grid_observe)

        # 只要靠近格网点都算
        if x_grid_observe - x_grid_center <= 1e-7 or 1 - x_grid_observe - x_grid_center <= 1e-7:
            is_in_x_grid_point = True
        if y_grid_observe - y_grid_center <= 1e-7 or 1 - y_grid_observe - y_grid_center <= 1e-7:
            is_in_y_grid_point = True

        # 差距值
        x_grid_inner = self.dx * (x_grid_observe - x_grid_center)
        y_grid_inner = self.dy * (y_grid_observe - y_grid_center)

        return x_grid_observe, y_grid_observe, x_grid_center, y_grid_center, x_grid_inner, y_grid_inner, is_in_x_grid_point, is_in_y_grid_point

    def feature_point_extraction(self, h, w):
        rows, cols = self.height.shape
        self.feature_point = np.zeros(((rows - 2 * h) * (cols - 2 * w) // (h * w), 2), dtype=np.int64)
        # 初始化一个空列表来存储每个方块中最大值的坐标
        max_points = []

        # 遍历每个1000个20*20的方块
        # 与边界框留出分析半径距离
        # 20 -->rows - 40
        # 16、30 --> rows - 40 - h不能整除
        # 分析网格数
        grid_num = 200
        for i in range(grid_num, rows - grid_num - h, h):
            for j in range(grid_num, cols - grid_num - h, w):
                # 获取当前3x3方块
                block = self.height[i:i + h, j:j + w]

                # 找到当前方块中最大值的索引
                local_max_idx = np.unravel_index(block.argmax(), block.shape)

                # 计算最大值的全局坐标
                global_max_idx = (i + local_max_idx[0], j + local_max_idx[1])

                # 添加最大值的坐标到列表
                max_points.append(global_max_idx)

        self.feature_point = np.array(max_points)
