# Story 2.4: BUSCO 质控

Status: done

## Story

As a **计算生物学研究员**,
I want **评估每个物种的注释完整性**,
so that **我能了解数据质量并识别潜在问题**。

## Acceptance Criteria

1. **AC1: 单物种 BUSCO 评估**
   - Given 标准化后的蛋白序列（`results/standardized/{species}/proteins.longest.fa.gz`）
   - When 执行 `qc_busco` rule
   - Then 调用 BUSCO 评估注释完整性
   - And 输出到 `results/qc/{species}/busco/`
   - And 输出文件包含：`short_summary.txt`, `full_table.tsv`, `missing_busco_list.tsv`

2. **AC2: Lineage 数据库配置**
   - Given 用户在配置文件中指定 lineage（如 `eukaryota_odb10`, `primates_odb10`）
   - When 执行 BUSCO
   - Then 使用指定的 lineage 数据库
   - And 支持自动下载 lineage 数据库（如果本地不存在）

3. **AC3: 多物种汇总表**
   - Given 多个物种完成 BUSCO
   - When 生成汇总
   - Then 创建 `results/qc/busco_summary.tsv` 包含所有物种的统计
   - And 列包含：species, complete, single_copy, duplicated, fragmented, missing, total, lineage

4. **AC4: 版本兼容性**
   - Given BUSCO 版本不兼容（不是 5.x 或 6.x）
   - When 检测版本
   - Then 返回 E_TOOL_VERSION 错误码
   - And 提供版本升级建议

5. **AC5: 错误处理**
   - Given BUSCO 执行失败（如 lineage 数据库缺失）
   - When 捕获错误
   - Then 返回适当的错误码（E_INPUT_MISSING 或 E_NONZERO_EXIT）
   - And 提供清晰的恢复建议

6. **AC6: 审计记录**
   - Given BUSCO 执行完成
   - When Runner 收集元数据
   - Then 生成 `.run.json` 包含 BUSCO 版本、lineage、输入 checksum、运行时间

## Tasks / Subtasks

- [x] Task 1: 创建 BUSCO Adapter (AC: #1, #4, #5, #6)
  - [x] 1.1 创建 `workflow/adapters/busco.py` 继承 BaseAdapter
  - [x] 1.2 实现 `check_version()` 检测 BUSCO 5.x/6.x
  - [x] 1.3 实现 `validate_inputs()` 验证蛋白序列文件和 lineage
  - [x] 1.4 实现 `build_command()` 构建 busco 命令行
  - [x] 1.5 实现 `expected_outputs()` 定义预期输出文件
  - [x] 1.6 实现 `parse_outputs()` 解析 short_summary.txt
  - [x] 1.7 实现 `timeout_seconds()` 返回超时配置（默认 30 分钟）
  - [x] 1.8 实现 `classify_error()` 分类错误类型

- [x] Task 2: 创建 conda 环境文件 (AC: #4)
  - [x] 2.1 创建 `workflow/envs/busco.yaml`
  - [x] 2.2 指定 BUSCO 版本范围（5.x 或 6.x）
  - [x] 2.3 包含所有 BUSCO 依赖（hmmer, metaeuk, prodigal 等）

- [x] Task 3: 创建 Snakemake 规则 (AC: #1, #2)
  - [x] 3.1 创建 `workflow/rules/qc.smk`
  - [x] 3.2 实现 `rule qc_busco`（单物种）
  - [x] 3.3 实现输入函数获取标准化蛋白序列路径
  - [x] 3.4 实现 lineage 配置读取

- [x] Task 4: 创建汇总规则 (AC: #3)
  - [x] 4.1 实现 `rule qc_busco_summary`（多物种汇总）
  - [x] 4.2 创建 `workflow/scripts/summarize_busco.py` 汇总脚本
  - [x] 4.3 解析所有物种的 short_summary.txt
  - [x] 4.4 生成 busco_summary.tsv

- [x] Task 5: 更新配置 schema (AC: #2)
  - [x] 5.1 更新 `schemas/config.schema.yaml` 添加 busco 配置节
  - [x] 5.2 定义 lineage 字段（必需）
  - [x] 5.3 定义 mode 字段（proteins/genome，默认 proteins）
  - [x] 5.4 定义 download_path 字段（lineage 数据库路径）

- [x] Task 6: 创建单元测试 (AC: #1-6)
  - [x] 6.1 创建 `workflow/adapters/test_busco.py`
  - [x] 6.2 测试版本检测逻辑（5.x, 6.x, 不兼容版本）
  - [x] 6.3 测试命令构建（不同 mode、lineage）
  - [x] 6.4 测试输出解析（short_summary.txt 格式）
  - [x] 6.5 测试错误分类（版本错误、输入缺失、执行失败）
  - [x] 6.6 测试汇总脚本逻辑

- [x] Task 7: 更新 Snakefile 集成 (AC: #1)
  - [x] 7.1 在 Snakefile 中添加 `include: "rules/qc.smk"` (已存在)
  - [x] 7.2 更新 `rule all` 包含 BUSCO 输出 (通过 qc.smk 函数实现)

## Dev Notes

### BUSCO 命令行接口 [Source: BUSCO 文档]

**基础命令：**
```bash
busco -i proteins.fa -l primates_odb10 -o output_name -m proteins -c 8
```

**关键参数：**
| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-i` | 输入文件（FASTA） | 必需 |
| `-l` | lineage 数据库 | 必需 |
| `-o` | 输出目录名 | 必需 |
| `-m` | 运行模式（proteins/genome/transcriptome） | genome |
| `-c` | CPU 线程数 | 1 |
| `--download_path` | lineage 数据库下载路径 | ~/busco_downloads |
| `--offline` | 离线模式，不自动下载 | false |
| `-f` | 强制覆盖已有输出 | false |

### BUSCO 输出格式 [Source: BUSCO 文档]

**short_summary.txt 解析：**
```
# BUSCO version is: 5.4.7
# The lineage dataset is: primates_odb10 (Creation date: 2024-01-08, number of BUSCOs: 13780)
# Summarized benchmarking in BUSCO notation for file proteins.fa
# BUSCO was run in mode: proteins

        ***** Results: *****

        C:95.2%[S:93.1%,D:2.1%],F:2.3%,M:2.5%,n:13780
        13118   Complete BUSCOs (C)
        12828   Complete and single-copy BUSCOs (S)
        290     Complete and duplicated BUSCOs (D)
        317     Fragmented BUSCOs (F)
        345     Missing BUSCOs (M)
        13780   Total BUSCO groups searched
```

**解析正则表达式：**
```python
import re

BUSCO_SUMMARY_PATTERN = re.compile(
    r'C:(\d+\.?\d*)%\[S:(\d+\.?\d*)%,D:(\d+\.?\d*)%\],F:(\d+\.?\d*)%,M:(\d+\.?\d*)%,n:(\d+)'
)

def parse_short_summary(content: str) -> dict:
    """解析 BUSCO short_summary.txt"""
    match = BUSCO_SUMMARY_PATTERN.search(content)
    if not match:
        raise ValueError("Cannot parse BUSCO summary")

    return {
        "complete_pct": float(match.group(1)),
        "single_copy_pct": float(match.group(2)),
        "duplicated_pct": float(match.group(3)),
        "fragmented_pct": float(match.group(4)),
        "missing_pct": float(match.group(5)),
        "total": int(match.group(6))
    }
```

### BuscoAdapter 实现框架 [Source: architecture.md#ADR-002]

```python
# workflow/adapters/busco.py

from pathlib import Path
from typing import Optional
import re
import subprocess

from .base import BaseAdapter, ToolSpec, AdapterContext
from ..lib.errors import ErrorCode, CompGeneError

class BuscoAdapter(BaseAdapter):
    """BUSCO 工具适配器"""

    SUPPORTED_VERSIONS = (5, 6)  # 支持 5.x 和 6.x

    @property
    def spec(self) -> ToolSpec:
        return ToolSpec(
            name="busco",
            version_cmd=["busco", "--version"],
            version_pattern=r"BUSCO (\d+)\.(\d+)\.(\d+)",
            min_version=(5, 0, 0),
            max_version=(7, 0, 0)  # 排他上界
        )

    def check_version(self) -> str:
        """检测 BUSCO 版本"""
        result = subprocess.run(
            ["busco", "--version"],
            capture_output=True,
            text=True
        )

        # BUSCO 5.x 输出: "BUSCO 5.4.7"
        # BUSCO 6.x 输出: "BUSCO 6.0.0"
        match = re.search(r'BUSCO (\d+\.\d+\.\d+)', result.stdout + result.stderr)
        if not match:
            raise CompGeneError(
                ErrorCode.E_TOOL_NOT_FOUND,
                "Cannot detect BUSCO version"
            )

        version = match.group(1)
        major = int(version.split('.')[0])

        if major not in self.SUPPORTED_VERSIONS:
            raise CompGeneError(
                ErrorCode.E_TOOL_VERSION,
                f"BUSCO {version} not supported. Requires 5.x or 6.x"
            )

        return version

    def validate_inputs(self, ctx: AdapterContext) -> None:
        """验证输入文件和配置"""
        input_file = Path(ctx.inputs["proteins"])
        if not input_file.exists():
            raise CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                f"Protein file not found: {input_file}"
            )

        lineage = ctx.config.get("busco", {}).get("lineage")
        if not lineage:
            raise CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                "BUSCO lineage not specified in config"
            )

    def build_command(self, ctx: AdapterContext) -> list[str]:
        """构建 BUSCO 命令行"""
        config = ctx.config.get("busco", {})

        cmd = [
            "busco",
            "-i", str(ctx.inputs["proteins"]),
            "-l", config["lineage"],
            "-o", ctx.wildcards["species"],
            "-m", config.get("mode", "proteins"),
            "-c", str(ctx.threads),
            "-f"  # 强制覆盖
        ]

        # 可选：指定下载路径
        if download_path := config.get("download_path"):
            cmd.extend(["--download_path", download_path])

        # 可选：离线模式
        if config.get("offline", False):
            cmd.append("--offline")

        return cmd

    def expected_outputs(self, ctx: AdapterContext) -> list[str]:
        """定义预期输出文件"""
        species = ctx.wildcards["species"]
        return [
            f"results/qc/{species}/busco/short_summary.txt",
            f"results/qc/{species}/busco/run_{species}/full_table.tsv",
            f"results/qc/{species}/busco/run_{species}/missing_busco_list.tsv"
        ]

    def parse_outputs(self, ctx: AdapterContext) -> dict:
        """解析 BUSCO 输出"""
        summary_path = Path(ctx.output_dir) / "short_summary.txt"
        content = summary_path.read_text()

        return parse_short_summary(content)

    def timeout_seconds(self, ctx: AdapterContext) -> int:
        """返回超时时间（秒）"""
        return ctx.config.get("busco", {}).get("timeout", 1800)  # 默认 30 分钟

    def classify_error(
        self,
        ctx: AdapterContext,
        returncode: int,
        stderr: str
    ) -> tuple[ErrorCode, bool]:
        """分类错误类型，返回 (错误码, 是否可重试)"""

        if "lineage" in stderr.lower() and "not found" in stderr.lower():
            return ErrorCode.E_INPUT_MISSING, False

        if "memory" in stderr.lower() or "oom" in stderr.lower():
            return ErrorCode.E_OOM, False

        if returncode == -9:  # SIGKILL (通常是超时)
            return ErrorCode.E_TIMEOUT, True

        return ErrorCode.E_NONZERO_EXIT, False
```

### Snakemake 规则 [Source: architecture.md#Naming-Patterns]

```python
# workflow/rules/qc.smk

from workflow.adapters.busco import BuscoAdapter

def get_proteins_path(wildcards):
    """获取标准化蛋白序列路径"""
    return f"results/standardized/{wildcards.species}/proteins.longest.fa.gz"

rule qc_busco:
    """运行 BUSCO 评估单个物种的注释完整性"""
    input:
        proteins = get_proteins_path
    output:
        summary = "results/qc/{species}/busco/short_summary.txt",
        directory = directory("results/qc/{species}/busco/run_{species}")
    params:
        lineage = config.get("busco", {}).get("lineage", "eukaryota_odb10"),
        mode = config.get("busco", {}).get("mode", "proteins"),
        download_path = config.get("busco", {}).get("download_path", "")
    threads: 8
    log:
        "logs/qc_busco/{species}.log"
    conda:
        "../envs/busco.yaml"
    script:
        "../scripts/run_adapter.py"

def get_all_busco_summaries():
    """获取所有物种的 BUSCO 汇总文件"""
    return expand(
        "results/qc/{species}/busco/short_summary.txt",
        species=[sp["name"] for sp in config["species"]]
    )

rule qc_busco_summary:
    """汇总所有物种的 BUSCO 结果"""
    input:
        summaries = get_all_busco_summaries()
    output:
        summary = "results/qc/busco_summary.tsv"
    log:
        "logs/qc_busco/summary.log"
    script:
        "../scripts/summarize_busco.py"
```

### conda 环境文件 [Source: architecture.md#conda-环境隔离]

```yaml
# workflow/envs/busco.yaml
name: compgene_busco
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - busco>=5.4,<7.0  # 支持 5.x 和 6.x
  - python>=3.9
  - hmmer>=3.3
  - metaeuk>=6
  - prodigal>=2.6
  - sepp>=4.5
  - biopython>=1.79
```

### 配置 Schema 更新 [Source: schemas/config.schema.yaml]

```yaml
# 新增 busco 配置节
busco:
  type: object
  required:
    - lineage
  properties:
    lineage:
      type: string
      description: "BUSCO lineage 数据库（如 eukaryota_odb10, primates_odb10）"
      examples:
        - "eukaryota_odb10"
        - "primates_odb10"
        - "mammalia_odb10"
    mode:
      type: string
      enum: ["proteins", "genome", "transcriptome"]
      default: "proteins"
      description: "BUSCO 运行模式"
    download_path:
      type: string
      description: "lineage 数据库存储路径"
    offline:
      type: boolean
      default: false
      description: "离线模式，不自动下载 lineage"
    timeout:
      type: integer
      default: 1800
      description: "超时时间（秒）"
  additionalProperties: false
```

### 汇总脚本 [Source: architecture.md#原子写入规范]

```python
# workflow/scripts/summarize_busco.py
"""汇总所有物种的 BUSCO 结果"""

import re
from pathlib import Path
import csv

BUSCO_SUMMARY_PATTERN = re.compile(
    r'C:(\d+\.?\d*)%\[S:(\d+\.?\d*)%,D:(\d+\.?\d*)%\],F:(\d+\.?\d*)%,M:(\d+\.?\d*)%,n:(\d+)'
)

LINEAGE_PATTERN = re.compile(r'The lineage dataset is: (\S+)')
VERSION_PATTERN = re.compile(r'BUSCO version is: (\S+)')

def parse_summary(path: Path) -> dict:
    """解析单个 short_summary.txt"""
    content = path.read_text()

    # 解析统计数据
    match = BUSCO_SUMMARY_PATTERN.search(content)
    if not match:
        raise ValueError(f"Cannot parse BUSCO summary: {path}")

    # 解析 lineage
    lineage_match = LINEAGE_PATTERN.search(content)
    lineage = lineage_match.group(1) if lineage_match else "unknown"

    # 解析版本
    version_match = VERSION_PATTERN.search(content)
    version = version_match.group(1) if version_match else "unknown"

    # 从路径提取物种名
    species = path.parent.parent.name

    return {
        "species": species,
        "complete_pct": float(match.group(1)),
        "single_copy_pct": float(match.group(2)),
        "duplicated_pct": float(match.group(3)),
        "fragmented_pct": float(match.group(4)),
        "missing_pct": float(match.group(5)),
        "total": int(match.group(6)),
        "lineage": lineage,
        "busco_version": version
    }

def main():
    # Snakemake 变量
    summaries = snakemake.input.summaries
    output_path = Path(snakemake.output.summary)

    # 解析所有汇总文件
    results = []
    for summary_path in summaries:
        try:
            result = parse_summary(Path(summary_path))
            results.append(result)
        except Exception as e:
            print(f"Warning: Failed to parse {summary_path}: {e}")

    # 原子写入
    temp_path = output_path.with_suffix('.tmp')

    fieldnames = [
        "species", "complete_pct", "single_copy_pct", "duplicated_pct",
        "fragmented_pct", "missing_pct", "total", "lineage", "busco_version"
    ]

    with open(temp_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(results)

    temp_path.rename(output_path)

if __name__ == "__main__":
    main()
```

### 错误处理 [Source: architecture.md#ADR-003]

| 场景 | 错误码 | 可重试 | 恢复建议 |
|------|--------|--------|----------|
| 蛋白文件不存在 | E_INPUT_MISSING | ❌ | 先运行 standardize 规则生成蛋白序列 |
| lineage 未配置 | E_INPUT_MISSING | ❌ | 在 config.yaml 中添加 busco.lineage |
| lineage 数据库不存在 | E_INPUT_MISSING | ❌ | 运行 `busco --download <lineage>` 或设置 download_path |
| BUSCO 版本不兼容 | E_TOOL_VERSION | ❌ | 安装 BUSCO 5.x 或 6.x |
| 执行超时 | E_TIMEOUT | ✅ | 增加 timeout 配置或减少输入序列 |
| 内存不足 | E_OOM | ❌ | 减少线程数或增加可用内存 |

### 从前序 Story 学到的经验 [Source: 2-3-ncbi-annotation-download.md]

1. **Adapter 模式一致性**：严格遵循 BaseAdapter 8 个方法接口
2. **原子写入**：所有输出文件先写 .tmp 再 rename
3. **版本检测**：使用正则表达式解析版本号，支持灵活匹配
4. **错误分类**：根据 stderr 内容和退出码精确分类错误类型
5. **单元测试共置**：test_busco.py 与 busco.py 放在同一目录
6. **配置 schema 更新**：添加新工具时同步更新 schema

### Project Structure Notes

**新增文件：**
```
workflow/
├── adapters/
│   ├── busco.py              # BUSCO Adapter（新建）
│   └── test_busco.py         # 单元测试（新建）
├── rules/
│   └── qc.smk                # QC 规则（新建）
├── scripts/
│   └── summarize_busco.py    # 汇总脚本（新建）
└── envs/
    └── busco.yaml            # conda 环境（新建）
```

**修改文件：**
```
workflow/Snakefile            # 添加 include: "rules/qc.smk"
schemas/config.schema.yaml    # 添加 busco 配置节
```

### References

- [Source: docs/planning-artifacts/prd.md#FR5a] - 系统可调用 BUSCO 评估装配/注释完整性
- [Source: docs/planning-artifacts/prd.md#FR5b] - 系统可生成每物种标准化目录结构和 BUSCO 汇总表
- [Source: docs/planning-artifacts/architecture.md#NFR13] - BUSCO 调用应支持版本 5.x/6.x
- [Source: docs/planning-artifacts/architecture.md#NFR15] - 所有外部工具调用应有超时保护
- [Source: docs/planning-artifacts/architecture.md#ADR-002] - 工具适配层 Adapter 模式
- [BUSCO 用户手册](https://busco.ezlab.org/busco_userguide.html) - 官方文档

## Dev Agent Record

### Agent Model Used

Claude Opus 4.5 (claude-opus-4-5-20251101)

### Debug Log References

N/A

### Completion Notes List

1. Implemented BuscoAdapter with all 8 BaseAdapter methods (48 unit tests passing)
2. Updated busco.yaml conda environment with version constraints and dependencies
3. Created qc.smk with qc_busco and qc_busco_summary rules
4. Created summarize_busco.py script with atomic write pattern
5. Updated config.schema.yaml with busco configuration section
6. All BUSCO-related tests passing (48 tests)

**Code Review Fixes (2026-02-04):**
7. Created run_busco.py to properly use BuscoAdapter (fixes CRITICAL: adapter not used)
8. Added .run.json audit record generation (fixes AC6)
9. Added full_table.tsv and missing_busco_list.tsv to qc.smk outputs (fixes AC1)
10. Eliminated code duplication by importing patterns from busco.py in summarize_busco.py

### File List

**Created:**
- `workflow/adapters/busco.py` - BUSCO adapter implementation
- `workflow/adapters/test_busco.py` - BUSCO adapter unit tests (40 tests)
- `workflow/scripts/summarize_busco.py` - Multi-species BUSCO summary script
- `workflow/scripts/test_summarize_busco.py` - Summary script unit tests (8 tests)
- `workflow/scripts/run_busco.py` - Snakemake script bridging BuscoAdapter (AC4, AC5, AC6)

**Modified:**
- `workflow/envs/busco.yaml` - Added version constraints and dependencies
- `workflow/rules/qc.smk` - Uses script directive with run_busco.py, outputs include full_table.tsv, missing_busco_list.tsv, .run.json
- `schemas/config.schema.yaml` - Added busco configuration section

