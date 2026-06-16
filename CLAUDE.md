# CLAUDE.md

本仓库是一套基于 Python (scanpy) 的**单细胞 RNA-seq 分析工作流**。本文件为 Claude Code 提供项目约定。

## 项目结构

```
config/config.yaml        所有分析参数（输入输出、QC、预处理、聚类、标记基因）
run_pipeline.py           命令行入口
src/scrna/
  io_utils.py             数据读写（10x mtx/h5、h5ad、csv/tsv 自动识别）
  qc.py                   质量控制：QC 指标、过滤、双细胞检测
  preprocess.py           归一化、log1p、高变基因、缩放、PCA
  cluster.py              批次整合、邻接图、UMAP、Leiden 聚类
  markers.py              标记基因 / 差异表达
  plotting.py             图表输出
  pipeline.py             端到端流程编排
```

## 运行

```bash
pip install -r requirements.txt
python run_pipeline.py --config config/config.yaml
```

输出：`results/`（`.h5ad` 结果、标记基因 CSV）与 `figures/`（QC、聚类、标记基因图）。

## 开发约定

- **配置驱动**：参数放 `config/config.yaml`，不要硬编码进函数。每个模块的 `run_*` 函数接收配置字典。
- **模块单一职责**：新增分析能力时，加到对应模块并在 `pipeline.py` 编排。
- **复现性**：所有随机过程传 `random_state`（默认取自 `config.cluster.random_state`）。
- **AnnData 约定**：obs=细胞、var=基因；原始计数存于 `layers['counts']`，完整表达谱存于 `.raw`，缩放只作用于 HVG 子集。
- **可选依赖优雅降级**：scrublet（双细胞）、harmonypy（批次整合）缺失时应跳过并提示，不中断流程。
- **无界面绘图**：plotting 使用 `matplotlib Agg` 后端保存图片，不弹窗。

## 数据约定

`data/` 与 `results/`、`figures/` 已被 `.gitignore` 忽略，不要提交大体积数据/结果文件。
