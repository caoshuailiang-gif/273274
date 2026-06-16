"""单细胞 RNA-seq 分析工作流。

模块划分：
    io_utils   —— 数据读写
    qc         —— 质量控制与过滤
    preprocess —— 标准化、高变基因、缩放、PCA
    cluster    —— 邻接图、UMAP、Leiden 聚类、批次整合
    markers    —— 标记基因 / 差异表达
    plotting   —— 图表输出
    pipeline   —— 端到端流程编排
"""

__version__ = "0.1.0"
