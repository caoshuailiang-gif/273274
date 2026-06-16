"""端到端流程编排：从原始数据到聚类、标记基因与图表。"""
from __future__ import annotations

import yaml

import scanpy as sc

from . import cluster, io_utils, markers, plotting, preprocess, qc


def load_config(config_path: str) -> dict:
    """读取 YAML 配置。"""
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def run_pipeline(config_path: str):
    """按配置执行完整单细胞分析流程，返回最终 AnnData。"""
    cfg = load_config(config_path)
    sc.settings.verbosity = 1

    io_cfg = cfg["io"]
    output_dir = io_cfg["output_dir"]
    figures_dir = io_cfg["figures_dir"]
    sample = io_cfg["sample_name"]

    print("=" * 60)
    print("步骤 1/6  读取数据")
    adata = io_utils.load_data(io_cfg["input"], io_cfg.get("transpose", False))

    print("=" * 60)
    print("步骤 2/6  质量控制")
    adata = qc.run_qc(adata, cfg["qc"])

    print("=" * 60)
    print("步骤 3/6  标准化与特征选择")
    adata = preprocess.run_preprocess(adata, cfg["preprocess"],
                                      random_state=cfg["cluster"]["random_state"])

    print("=" * 60)
    print("步骤 4/6  降维与聚类")
    adata = cluster.run_cluster(adata, cfg["cluster"])

    print("=" * 60)
    print("步骤 5/6  标记基因分析")
    adata = markers.run_markers(adata, cfg["markers"], output_dir, sample)

    print("=" * 60)
    print("步骤 6/6  生成图表并保存结果")
    plotting.run_plots(adata, figures_dir)
    io_utils.save_adata(adata, output_dir, sample, "processed")

    print("=" * 60)
    print(f"完成！结果目录: {output_dir}/  图表目录: {figures_dir}/")
    return adata
