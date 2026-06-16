#!/usr/bin/env python
"""单细胞 RNA-seq 分析工作流命令行入口。

用法：
    python run_pipeline.py --config config/config.yaml

也可单独运行某些步骤——见 src/scrna 下各模块。
"""
import argparse
import sys
from pathlib import Path

# 允许从源码目录直接运行
sys.path.insert(0, str(Path(__file__).parent / "src"))

from scrna.pipeline import run_pipeline


def main():
    parser = argparse.ArgumentParser(description="单细胞 RNA-seq 分析工作流")
    parser.add_argument("--config", "-c", default="config/config.yaml",
                        help="YAML 配置文件路径 (默认: config/config.yaml)")
    args = parser.parse_args()
    run_pipeline(args.config)


if __name__ == "__main__":
    main()
