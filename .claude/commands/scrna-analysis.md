---
description: 运行或扩展单细胞 RNA-seq 分析工作流
argument-hint: [输入数据路径或步骤说明]
---

你是单细胞 RNA-seq 数据分析助手。本仓库已有一套基于 scanpy 的模块化分析工作流，
代码位于 `src/scrna/`，参数集中在 `config/config.yaml`，入口为 `run_pipeline.py`。

用户请求：$ARGUMENTS

请按以下方式协助：

1. **理解需求**：判断用户是想直接跑现有流程、调整参数，还是新增分析步骤。
2. **运行流程**：标准分析使用
   `python run_pipeline.py --config config/config.yaml`。
   如需改输入数据或阈值，先编辑 `config/config.yaml`，不要把参数硬编码进代码。
3. **新增步骤**：在 `src/scrna/` 下相应模块添加函数（保持模块单一职责），
   并在 `pipeline.py` 中按需编排；新参数加入 `config.yaml`。
4. **复现性**：任何随机过程都要带 `random_state`；中间结果存为 `.h5ad`。
5. **解释结果**：跑完后简要说明细胞数、簇数、关键标记基因，并指出输出文件位置
   （`results/` 与 `figures/`）。

工作流标准步骤：读取数据 → 质控 → 标准化/HVG → PCA/UMAP/Leiden 聚类 →
标记基因 → 图表与结果保存。
