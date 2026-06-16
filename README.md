# 单细胞 RNA-seq 分析工作流

一套基于 Python（[scanpy](https://scanpy.readthedocs.io/)）的模块化单细胞 RNA-seq 分析流程，
配置驱动、可复现，并集成 Claude Code 工作流（`/scrna-analysis` 命令）。

## 功能

从原始表达矩阵到聚类与标记基因的完整端到端分析：

| 步骤 | 内容 |
|------|------|
| 1. 读取数据 | 自动识别 10x mtx 目录 / 10x h5 / `.h5ad` / csv·tsv |
| 2. 质量控制 | QC 指标、按基因数/线粒体比例过滤、可选双细胞检测（scrublet） |
| 3. 标准化 | 文库归一化、log1p、高变基因、缩放、PCA |
| 4. 降维聚类 | 邻接图、UMAP、Leiden 聚类，可选 Harmony 批次整合 |
| 5. 标记基因 | 各簇差异表达（Wilcoxon 等），导出 CSV |
| 6. 输出 | 保存 `.h5ad` 结果与 QC/聚类/标记基因图表 |

## 安装

```bash
pip install -r requirements.txt
```

建议 Python ≥ 3.10。`scrublet`（双细胞检测）与 `harmonypy`（批次整合）为可选依赖，缺失时自动跳过。

## 快速开始

1. 把数据放到 `data/raw/`（例如 10x 的 `pbmc3k` 目录）。
2. 按需修改 `config/config.yaml`（输入路径、过滤阈值、聚类分辨率等）。
3. 运行：

```bash
python run_pipeline.py --config config/config.yaml
```

结果输出到 `results/`，图表输出到 `figures/`。

### 用 PBMC3k 示例数据试跑

```python
import scanpy as sc
adata = sc.datasets.pbmc3k()            # 下载示例数据
adata.write_h5ad("data/raw/pbmc3k.h5ad")
```

然后把 `config.yaml` 的 `io.input` 改为 `data/raw/pbmc3k.h5ad` 即可。

## 配置说明

所有参数集中在 `config/config.yaml`，分为 `io` / `qc` / `preprocess` / `cluster` / `markers`
五组，每个键都有中文注释。修改参数无需改动代码。

## 在 Claude Code 中使用

本仓库带有 `CLAUDE.md`（项目约定）和自定义斜杠命令 `/scrna-analysis`。
在 Claude Code 中输入：

```
/scrna-analysis 用 config/config.yaml 跑一遍分析，并解释聚类结果
```

即可让 Claude 运行流程、调整参数或新增分析步骤。

## 项目结构

```
config/config.yaml     分析参数
run_pipeline.py        命令行入口
src/scrna/             分析模块（io_utils / qc / preprocess / cluster / markers / plotting / pipeline）
.claude/commands/      Claude Code 斜杠命令
CLAUDE.md              项目约定
```

## 扩展

新增分析能力时，在 `src/scrna/` 对应模块添加函数，并在 `pipeline.py` 中编排，
新参数加入 `config.yaml`。详见 `CLAUDE.md` 的开发约定。
