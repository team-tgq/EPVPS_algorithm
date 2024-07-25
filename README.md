# Language

- [English](#english)
- [中文](#中文)

---

## English
## Project Introduction
In this Project, we introduce an algorithm named XXX, which provides an important reference for multi-viewpoint planning. Firstly, we determined the candidate viewpoints based on relative elevation, and then we used ECKM algorithm to determine the initial centre of the clustering of candidate viewpoints. In the clustering process, we used the ECKM idea to determine the neighbouring clusters. Secondly, we proposed the viewpoint evaluation metric WCOM, and in the iterative process, the viewable domain of the candidate points was decomposed into the set of coverage contributing points and the set of overlapping contributing points, based on which WCOM was computed and stored as the minimum heap in the order of its value. Finally, the filtering algorithm adjusts the set of covering and overlapping contributing points to complete the operation of adding and deleting viewpoints in the cluster until there is only one most valuable viewpoint left in each cluster and then the filtering is finished. To evaluate the effectiveness of the Ours algorithm, we compared it with the Candidate Viewpoint Filtering (CVF) algorithm and Simulated Annealing (SA) algorithm. The experimental results show that the Ours algorithm performs best in terms of computational efficiency, coverage and overlap.

## Code Introduction
### Preprocessing
The perprocessing folder contains pre-processing methods based on the GDAL library for extracting candidate points from terrain files, resampling and so on.
### The main programme of the algorithm
keam_clustering implements the whole process of clustering, including initialising the cluster centre, determining adjacent clusters, etc.; cul_coverage_and_overlap calculates the WCOM within each cluster and stores it as a minimum heap; filter implements the filtering algorithm.
### Comparative experiments
The comparative_experiments folder contains specific implementations of the SA, CVF algorithms
### Reference line classes
The reference line class structure and corresponding methods are defined under models.
### Implementation of viewshed algorithms
The services folder contains the two steps of viewshed algorithm implementation, read_dem reads the input terrain file data and converts it to structured data, and partition_optimal_reference_lineal_gorithm implements the core part of DPDERL.
### Data
The /data/ file directory contains all the DSM data used in this study.

---

## 中文

## 项目介绍
本研究实现了一种名为XXX的算法，为多视点规划提供了重要参考。首先，我们基于相对高程确定了候选视点，随后利用ECKM算法确定了候选视点聚类的初始中心，并对候选视点进行了聚类。在聚类过程中，我们采用了ECKM思想确定相邻簇。其次，我们提出了视点评估指标WCOM，在迭代过程中，将候选点的可视域分解为覆盖贡献点集和重叠贡献点集，并基于此计算WCOM并按照其值排序存储为最小堆。最后，通过过滤算法调整覆盖贡献点集和重叠贡献点集，完成了簇内新增和删除视点的操作，直至每个簇仅剩一个最有价值视点则过滤结束。为了评估Ours算法的有效性，我们将其与候选视点过滤（CVF）算法和模拟退火 (SA) 算法进行了比较。实验结果表明，Ours算法在计算效率、覆盖率和重叠率等性能方面表现最佳。

## 代码简介
### 预处理
perprocessing文件夹包含了基于GDAL库提取地形文件候选点，重采样等预处理方法
### 本算法主程序
keam_clustering实现了聚类的整个过程，包含初始化簇中心，确定相邻簇等；cul_coverage_and_overlap计算各簇内WCOM并存储为最小堆；filter实现了过滤算法。
### 对比试验
comparative_experiments文件夹包含了SA、CVF算法的具体实现
### 参考线类
models下定义了参考线类结构及相应方法
### 视域算法具体实现
services下包含了视域算法实现的两个步骤，read_dem读取输入的地形文件数据并转为结构化数据，partition_optimal_reference_lineal_gorithm实现了DPDERL核心部分。
### 数据
在/data/文件目录下包含了本研究所用的所有DSM数据


## 版本控制

该项目使用Git进行版本管理。您可以在repository参看当前可用版本。
