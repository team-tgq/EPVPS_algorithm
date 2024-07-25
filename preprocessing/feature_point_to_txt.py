from services import read_dem

# DEM相对高程数组
dem_data = read_dem.Dem("../data/CW_relative_3000_1.tif", 1)

feature_point_array = dem_data.feature_point.copy()

# 指定要写入的文件路径
file_path = "../output/csv/our_mountain_feature.txt"

# 将数组写入文件
with open(file_path, 'w') as f:
    for row in feature_point_array:
        # 将每行转换为以逗号分隔的字符串，并写入文件
        line = ','.join(map(str, row))
        f.write(line + '\n')

print(f"文件已保存到 {file_path}")
