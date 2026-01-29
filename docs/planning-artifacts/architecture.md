---
stepsCompleted: [1, 2, 3, 4, 5, 6, 7, 8]
inputDocuments:
  - "prd.md"
  - "product-brief-CompGene-2026-01-20.md"
workflowType: 'architecture'
project_name: 'CompGene'
user_name: 'Yuchenghuang'
date: '2026-01-20'
---

# Architecture Decision Document

_This document builds collaboratively through step-by-step discovery. Sections are appended as we work through each architectural decision together._

## Project Context Analysis

### Requirements Overview

**Functional Requirements:**
- 49 条 FR 覆盖 10 个能力领域
- 核心是"工具编排"而非"算法实现"
- 主要工作是数据格式转换和流程串联

**Non-Functional Requirements:**
- 性能：上游步骤（标准化/QC）单物种 < 30 分钟；OrthoFinder 按物种集合规模计时
- 可靠性：100% 可复现，断点恢复
- 集成：NCBI 速率限制，多工具版本兼容
- 可维护性：模块低耦合，配置驱动

**Scale & Complexity:**
- Primary domain: CLI Pipeline / Bioinformatics
- Complexity level: Medium
- Estimated architectural components: 10 模块

### Scope Boundaries (Critical Clarifications)

#### 增量运行定义

| 增量类型 | 能力 | 说明 |
|----------|------|------|
| **增量数据**（同物种更新） | ✅ 支持 | checksum 变化触发下游重算 |
| **增量物种**（新增物种） | ⚠️ 部分支持 | 上游缓存可复用，OrthoFinder 需全量重算 |

> OrthoFinder 的正交分组是全物种集合计算，新增物种 = 重算 orthogroups。增量仅指复用"标准化/蛋白提取/BUSCO"等上游产物。

#### 性能指标现实化

| 步骤 | 性能目标 | 条件 |
|------|----------|------|
| 上游（标准化/BUSCO/Liftoff） | 单物种 < 30 分钟 | 标准服务器 |
| OrthoFinder | 按物种集合规模 | 3 物种 ~1-2h，7 物种 ~数小时 |
| 缓存命中 | < 1 秒 | 本地已有产物，不含网络/解压 |

#### 输出契约（Data Products）

```
{output_dir}/
├── standardized/{species}/
│   ├── genome.fa.gz
│   ├── annotation.gff3.gz
│   └── proteins.longest.fa.gz
├── orthology/
│   ├── orthogroups.tsv
│   ├── presence_absence.tsv
│   └── copy_number.tsv
├── annotation/eggnog/
│   └── annotations.tsv
├── validation/liftoff/
│   └── true_absence_candidates.tsv
├── qc/
│   └── busco_summary.tsv
├── expression/
│   ├── DE_results.tsv
│   └── normalized_counts.tsv
├── reports/
│   ├── summary.html
│   └── summary.md
└── logs/
    ├── run.log
    └── run.json
```

> 输出契约稳定 = 下游不崩。换工具/换参数只要契约不变，集成方无感知。

### Technical Constraints & Dependencies

- **外部工具依赖**：OrthoFinder 2.5.x, eggNOG-mapper 2.1.x, Liftoff 1.6.x, BUSCO 5.x/6.x, DESeq2 (R/Bioconductor)
- **运行环境**：本地服务器，Linux/macOS，Python 3.x
- **数据格式**：GFF3/GTF 输入，TSV 输出
- **网络依赖**：NCBI API（速率限制 3 req/s）
- **工作流引擎**：Snakemake（推荐）或 Nextflow - 不自建调度器

### Architecture Module Breakdown (10 Modules)

| # | 模块 | 职责 |
|---|------|------|
| 1 | **Config & Validation** | YAML 解析、schema 验证、参数覆盖 |
| 2 | **Tool Adapter Layer** | 每个外部工具一个 adapter（输入检查、命令拼装、版本检测） |
| 3 | **Data Standardization** | FASTA/GFF 清洗、ID 映射、代表转录本选择 |
| 4 | **QC** | BUSCO + 基础统计 |
| 5 | **Orthology Core** | OrthoFinder 调用与结果解析 |
| 6 | **Functional Annotation** | eggNOG-mapper 调用与汇总 |
| 7 | **Presence/Absence Derivation** | 矩阵生成、阈值规则、缺失候选集 |
| 8 | **Validation & Cross-check** | Liftoff 核验、真缺失筛选 |
| 9 | **Reporting** | 表格汇总、图、可追溯 metadata |
| 10 | **Runtime Services** | 日志/审计、缓存、状态、错误恢复 |

### Cross-Cutting Concerns

1. **日志与审计** - 统一框架，贯穿所有模块
2. **缓存与状态** - 文件级依赖，checksum 驱动失效
3. **配置管理** - YAML 解析，验证，参数覆盖
4. **错误处理** - 统一错误码，重试逻辑，恢复建议
5. **版本追踪** - 工具版本、数据库版本、分析参数、输入快照

### Key Architecture Decision Required

**ADR-001: 工作流引擎选型**
- **Snakemake**（推荐）：Python 生态友好，bioinformatics 主流，天然支持文件级依赖/断点/conda
- **Nextflow**：容器/云/HPC 更强，适合未来扩展
- **自建**：❌ 不推荐 - 会把范围从"中等"拉到"重型平台工程"

## Starter Template Evaluation

### Primary Technology Domain

CLI Pipeline / Bioinformatics - 基于 Python，本地服务器运行

### Starter Options Considered

| 选项 | 评估 |
|------|------|
| **Snakemake 9.x** | ✅ 推荐 - Python 原生、文件依赖驱动、conda 集成 |
| **Nextflow 25.x** | ⚪ 备选 - 云/HPC/平台生态更强，未来扩展时考虑 |
| **自建调度器** | ❌ 排除 - 范围膨胀风险 |

### Selected Starter: Snakemake 9.14.8

**Rationale for Selection:**
- Python 技术栈一致性
- 文件依赖模型匹配生物信息学工作流
- 断点续跑、conda 集成开箱即用
- 社区活跃（>11 citations/week）
- 本地服务器运行，不需要 Seqera/nf-core 的平台化能力

**Nextflow 备选条件（未来切换触发点）：**
- 需要上云/HPC 集群
- 团队已普遍使用 Nextflow/nf-core
- 需要 Seqera Platform 的监控/调度能力

### Architectural Decisions Provided by Starter

**Language & Runtime:**
- Python 3.11+
- Snakemake 9.x DSL
- 插件点：logger plugin（自定义日志）、config schema validation（PEP/eido）

**Code Organization:**
```
compgene/
├── workflow/
│   ├── Snakefile              # 主入口
│   ├── rules/                 # 模块化规则
│   │   ├── standardize.smk
│   │   ├── qc.smk
│   │   ├── orthology.smk
│   │   ├── annotation.smk
│   │   ├── validation.smk
│   │   ├── expression.smk
│   │   └── reporting.smk
│   ├── scripts/               # Python 辅助脚本
│   └── envs/                  # 每模块 conda 环境 YAML
│       ├── orthofinder.yaml
│       ├── eggnog.yaml
│       ├── liftoff.yaml
│       ├── busco.yaml
│       └── deseq2.yaml
├── config/
│   └── config.yaml            # 用户配置
├── schemas/
│   └── config.schema.yaml     # 配置校验 schema
├── profiles/
│   └── local/
│       └── config.yaml        # 本地资源配置（cores、tempdir）
├── resources/                 # 静态资源/数据库
├── results/
│   └── meta/                  # 审计元数据
│       ├── tool_versions.tsv
│       └── input_checksums.tsv
└── logs/                      # 统一日志路径
    └── {rule}/{species}/
```

**可复现性三层产物：**
1. `results/meta/tool_versions.tsv` - 每次 run 自动采集 `--version`
2. `workflow/envs/*.yaml` - 环境文件纳入 git，版本锁定
3. `results/meta/input_checksums.tsv` - 输入数据 checksum

**断点续跑团队约定：**
- 大文件规则：先写临时文件再原子重命名（避免 incomplete 半成品）
- 日志路径：`logs/{rule}/{species}.log`（方便审计）
- 运行命令：`snakemake --rerun-incomplete`

**缓存边界定义：**

| 缓存类型 | 规则 | 说明 |
|----------|------|------|
| **强缓存** | standardize, proteins, BUSCO, eggNOG | 同输入+同版本=同输出 |
| **半缓存** | OrthoFinder | 固定线程、记录命令与版本、一般可复现 |

**配置校验落地：**
- `schemas/config.schema.yaml` 定义必需字段和类型
- Snakefile 启动时调用 `validate(config, schema="schemas/config.schema.yaml")`

**Initialization Command:**

```bash
# 创建项目结构
mkdir -p compgene/{workflow/{rules,scripts,envs},config,schemas,profiles/local,resources,results/meta,logs}

# 安装 Snakemake
conda create -n compgene python=3.11 snakemake=9.14 -c conda-forge -c bioconda

# 初始化主 Snakefile
cat > compgene/workflow/Snakefile << 'EOF'
from snakemake.utils import min_version, validate
min_version("9.0")

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

include: "rules/standardize.smk"
include: "rules/qc.smk"
include: "rules/orthology.smk"
include: "rules/annotation.smk"
include: "rules/validation.smk"
include: "rules/expression.smk"
include: "rules/reporting.smk"

rule all:
    input:
        "results/reports/summary.html"
EOF
```

**Development Experience:**
- `snakemake -n` / `snakemake --dry-run` 预览执行计划
- `snakemake --dag | dot -Tsvg > dag.svg` 可视化 DAG
- `snakemake --rerun-incomplete` 断点续跑
- `snakemake --lint` 检查规则质量

**DAG 工程规范（避免超大/动态 DAG 问题）：**
- 尽量用静态文件依赖
- 减少 checkpoints 使用
- 动态部分前移到生成 manifest（species.tsv）

## Core Architectural Decisions

### Decision Priority Analysis

**Critical Decisions (Block Implementation):**
1. 数据架构：文件流 + SQLite 元数据层
2. 工具适配层：Adapter 模块 + 统一 Runner
3. 错误处理：详细错误码 + 恢复建议 + 混合重试
4. 日志审计：双输出 + 规则级审计

**Deferred Decisions (Post-MVP):**
- 云/HPC 部署（触发条件：需要 Nextflow 切换时）
- Web 报告服务（触发条件：团队协作需求上升时）

### ADR-001: 数据架构

**Decision:** 文件流 + SQLite 元数据层（混合模式）

**Context:**
- Snakemake 依赖文件产物进行缓存和断点续跑
- 跨步骤 join、追溯、报告汇总需要快速查询

**Decision:**
- **文件产物（Source of Truth）**：FASTA/GFF/TSV，Snakemake 依赖跟踪
- **SQLite（Index of Record）**：`results/meta/compgene.db`
  - ID 映射（gene ↔ transcript ↔ protein ↔ orthogroup）
  - 运行元数据（版本、参数、checksum、时间戳）
  - 报告/统计的快速 join

**Key Tables:**
- `species`, `gene`, `transcript`, `protein`
- `orthogroup`, `orthogroup_member`
- `qc_busco`
- `run`, `artifact`, `tool_version`

**Write Strategy:**
- 单写入 rule 模式：各模块输出 TSV/JSON，`rule ingest_db` 统一写入
- 避免多 rule 并发写 SQLite

**ID Mapping Strategy:**
- 三层 ID：raw（原始）+ norm（标准化）+ uid（内部稳定）
- 代表蛋白选择可追溯（is_representative 标记）

### ADR-002: 工具适配层

**Decision:** Adapter 模块 + 统一 Runner

**Architecture:**
```
workflow/
├── adapters/           # 纯 Python 模块（可单测）
│   ├── base.py         # BaseAdapter, ToolSpec, RunResult
│   ├── orthofinder.py
│   ├── eggnog.py
│   ├── liftoff.py
│   ├── busco.py
│   └── deseq2.py
├── scripts/
│   └── run_adapter.py  # 统一 Runner
└── lib/
    ├── errors.py       # 错误分类
    └── audit.py        # 审计元数据
```

**Adapter Interface:**
```python
class BaseAdapter(ABC):
    spec: ToolSpec
    def check_version(self) -> str
    def validate_inputs(self, ctx) -> None
    def build_command(self, ctx) -> list[str]
    def expected_outputs(self, ctx) -> list[str]
    def parse_outputs(self, ctx) -> dict
    def timeout_seconds(self, ctx) -> int
    def classify_error(self, ctx, returncode, stderr) -> tuple[str, bool]
```

**Runner Responsibilities:**
- 版本检测、输入校验
- 命令执行（subprocess，非 shell）
- 超时 kill（SIGTERM → grace → SIGKILL）
- stdout/stderr 捕获
- 输出解析 → summary.json
- 审计元数据 → run.json

### ADR-003: 错误处理与恢复

**Decision:** 详细错误码 + 恢复建议 + 混合重试

**Error Code Schema:**
| 错误码 | 含义 | 可重试 | 恢复建议 |
|--------|------|--------|----------|
| E_INPUT_MISSING | 输入缺失 | ❌ | 检查输入文件路径 |
| E_INPUT_FORMAT | 格式错误 | ❌ | 验证 GFF/FASTA 格式 |
| E_TOOL_NOT_FOUND | 工具缺失 | ❌ | 激活 conda 环境 |
| E_TOOL_VERSION | 版本不符 | ❌ | 更新工具版本 |
| E_TIMEOUT | 超时 | ✅ | 增加超时或减少输入 |
| E_OOM | 内存不足 | ❌ | 减少 threads 或增加 memory |
| E_NET_RATE_LIMIT | 网络限流 | ✅ | 等待后重试 |
| E_DISK_FULL | 磁盘满 | ❌ | 清理磁盘空间 |
| E_NONZERO_EXIT | 其他错误 | ❌ | 查看日志 |

**Retry Strategy (Hybrid):**
- Snakemake 层：`retries: 1` 兜底
- Runner 层：网络错误 3 次 exponential backoff

**Snakemake 运行规范:**
| 场景 | 命令 |
|------|------|
| 正常运行 | `snakemake --cores 8` |
| 断点续跑 | `snakemake --rerun-incomplete` |
| 部分失败后继续 | `snakemake --keep-going` |
| 调试单规则 | `snakemake {target} --forcerun` |

### ADR-004: 日志与审计

**Decision:** 双输出日志 + 规则级审计

**Log Format (Dual Output):**
- `logs/{rule}/{wildcards}.log` - 人类可读文本
- `logs/{rule}/{wildcards}.jsonl` - JSON Lines（jq/pandas 可查询）

**Audit Granularity (Rule-level):**
每个 rule 执行生成：
- `results/meta/{rule}/{wildcards}.run.json`
  ```json
  {
    "rule": "orthofinder",
    "wildcards": {"species_set": "lemur_macaque"},
    "cmd": ["orthofinder", "-f", "..."],
    "tool_version": "2.5.5",
    "input_checksums": {"proteins/": "sha256:..."},
    "threads": 8,
    "runtime_seconds": 3600,
    "exit_code": 0,
    "timestamp": "2026-01-20T10:15:32Z"
  }
  ```

**Audit Ingestion:**
- `rule ingest_db` 收集所有 `.run.json` 写入 SQLite
- 支持追溯：任意输出 → 生成它的命令/版本/输入

### Decision Impact Analysis

**Implementation Sequence:**
1. 项目结构初始化（目录、Snakefile 骨架）
2. `lib/errors.py` + `lib/audit.py`（基础设施）
3. `adapters/base.py` + `run_adapter.py`（框架）
4. 各工具 Adapter 实现（按 DAG 顺序）
5. SQLite schema + `rule ingest_db`
6. 报告生成（查询 SQLite）

**Cross-Component Dependencies:**
- 所有 Adapter 依赖 `lib/errors.py` 的错误码
- 所有 rule 输出的 `.run.json` 依赖统一 schema
- 报告依赖 SQLite 完整性

## Implementation Patterns & Consistency Rules

### Pattern Categories Defined

**Critical Conflict Points Addressed:** 10 个潜在冲突点

### Naming Patterns

#### Snakemake Naming Conventions

**Rule Naming:** `{module}_{action}`
```python
# ✅ Good
rule orthology_infer:
rule qc_busco:
rule standardize_proteins:
rule annotation_eggnog:
rule validation_liftoff:

# ❌ Bad
rule run_orthofinder:      # 动词在前
rule BUSCO:                # 全大写
rule doEggnog:             # camelCase
```

**Wildcard Naming:** 描述性、snake_case
```python
# ✅ Good
"{species}/proteins.fa"
"{species_set}/orthogroups.tsv"
"{rule_name}/{species}.log"

# ❌ Bad
"{sp}/proteins.fa"         # 太简短
"{speciesSet}/..."         # camelCase
```

**Output Path Pattern:** 层级式
```
results/
├── standardized/{species}/
├── qc/{species}/busco/
├── orthology/{species_set}/
├── annotation/eggnog/{species}/
├── validation/liftoff/{source}_to_{target}/
├── expression/
├── matrices/
├── reports/
├── meta/{rule}/{wildcards}/
└── logs/{rule}/{wildcards}/
```

#### Python Naming Conventions

**Functions & Variables:** snake_case
```python
# ✅ Good
def parse_orthogroups(input_path: Path) -> dict:
    species_count = len(species_list)
    gene_ids = extract_gene_ids(gff_path)

# ❌ Bad
def parseOrthogroups():    # camelCase
def ParseOrthogroups():    # PascalCase
```

**Classes:** PascalCase
```python
# ✅ Good
class OrthoFinderAdapter(BaseAdapter):
class BuscoResult:
class ToolSpec:

# ❌ Bad
class orthofinder_adapter:  # snake_case
class BUSCO_RESULT:         # UPPER_CASE
```

**Constants:** UPPER_SNAKE_CASE
```python
# ✅ Good
DEFAULT_TIMEOUT = 3600
MAX_RETRY_COUNT = 3
NCBI_RATE_LIMIT = 3  # requests per second

# ❌ Bad
defaultTimeout = 3600
default_timeout = 3600
```

### Structure Patterns

#### Project Organization

```
compgene/
├── workflow/
│   ├── Snakefile                    # 主入口
│   ├── rules/                       # Snakemake rules
│   │   ├── standardize.smk
│   │   ├── qc.smk
│   │   └── ...
│   ├── scripts/                     # 被 rule 调用的脚本
│   │   └── run_adapter.py
│   ├── adapters/                    # 工具适配器
│   │   ├── __init__.py
│   │   ├── base.py
│   │   ├── orthofinder.py
│   │   ├── test_orthofinder.py     # 单元测试共置
│   │   └── ...
│   ├── lib/                         # 共享库
│   │   ├── errors.py
│   │   ├── audit.py
│   │   └── io.py
│   └── envs/                        # conda 环境
│       └── *.yaml
├── config/
│   └── config.yaml
├── schemas/
│   └── config.schema.yaml
├── profiles/
│   └── local/config.yaml
├── tests/                           # 集成测试集中
│   ├── integration/
│   │   ├── test_full_pipeline.py
│   │   └── test_data/
│   └── conftest.py
├── resources/
├── results/
└── logs/
```

#### Test Organization (Hybrid)

**Unit Tests:** 与模块共置
```
adapters/
├── orthofinder.py
├── test_orthofinder.py    # pytest 自动发现
├── eggnog.py
└── test_eggnog.py
```

**Integration Tests:** 集中在 tests/
```
tests/
├── integration/
│   ├── test_standardize_to_orthology.py
│   └── test_full_pipeline.py
└── conftest.py            # 共享 fixtures
```

### Format Patterns

#### TSV Output Format

**Column Naming:** snake_case
```tsv
gene_id	species_id	orthogroup_id	is_representative
mmur:GENE001	mmur	OG0001234	1
lcat:GENE042	lcat	OG0001234	1
```

**Standard Columns Across Files:**
| 列名 | 类型 | 说明 |
|------|------|------|
| `gene_id` | string | `{species}:{raw_id}` 格式 |
| `species_id` | string | 物种短码（mmur, lcat, mfa） |
| `orthogroup_id` | string | OG 编号 |
| `timestamp` | string | ISO 8601 |

#### JSON Output Format

**Field Naming:** snake_case
```json
{
  "rule_name": "orthology_infer",
  "species_set": "lemur_macaque",
  "tool_version": "2.5.5",
  "input_checksums": {
    "proteins_dir": "sha256:abc123..."
  },
  "runtime_seconds": 3600,
  "exit_code": 0,
  "timestamp": "2026-01-20T10:15:32Z"
}
```

#### Log Format

**Human-readable (.log):**
```
2026-01-20 10:15:32 [INFO] orthology_infer: Starting with 3 species
2026-01-20 10:15:33 [INFO] orthology_infer: Input validation passed
2026-01-20 11:15:32 [INFO] orthology_infer: Completed in 3600s
```

**Machine-readable (.jsonl):**
```json
{"ts":"2026-01-20T10:15:32Z","level":"INFO","rule":"orthology_infer","msg":"Starting","species_count":3}
```

### Process Patterns

#### Input Validation Timing

**Rule:** 在 Adapter.validate_inputs() 中执行，Runner 调用
- 文件存在性检查
- 格式基础验证（header 检查）
- 版本兼容性检查

#### File Write Pattern (Atomic)

**Rule:** 所有大文件输出使用原子写入
```python
# ✅ Good - 原子写入
def write_output(data, output_path: Path):
    temp_path = output_path.with_suffix('.tmp')
    write_to_file(data, temp_path)
    temp_path.rename(output_path)  # 原子操作

# ❌ Bad - 直接写入（可能产生不完整文件）
def write_output(data, output_path: Path):
    write_to_file(data, output_path)
```

### Enforcement Guidelines

**All AI Agents MUST:**

1. 使用 `{module}_{action}` 格式命名 Snakemake rules
2. 使用 snake_case 命名所有 Python 函数、变量、TSV 列、JSON 字段
3. 使用层级式输出路径 `results/{category}/{wildcards}/`
4. 所有大文件写入使用原子模式（先写 .tmp 再 rename）
5. 所有日期时间使用 ISO 8601 格式
6. 单元测试与被测模块共置，集成测试放 tests/integration/

**Pattern Verification:**
- `snakemake --lint` 检查 rule 质量
- `ruff check` / `flake8` 检查 Python 代码风格
- CI 中运行 `pytest` 确保测试通过

### Anti-Patterns to Avoid

```python
# ❌ 避免：混合命名风格
rule runOrthoFinder:           # camelCase rule name
    output: "results/OrthoFinder_output.tsv"  # 混合大小写

# ❌ 避免：扁平输出结构
rule orthology:
    output: "orthogroups.tsv"  # 无层级

# ❌ 避免：非原子写入
rule report:
    shell: "python script.py > {output}"  # 直接重定向

# ✅ 正确做法
rule orthology_infer:
    output: "results/orthology/{species_set}/orthogroups.tsv"
    script: "scripts/run_adapter.py"  # 内部使用原子写入
```

## Project Structure & Boundaries

### Complete Project Directory Structure

```
compgene/
├── README.md
├── LICENSE
├── pyproject.toml                    # Python 项目配置（含 ruff/pytest）
├── setup.py                          # 可选：pip install -e .
├── .gitignore
├── .github/
│   └── workflows/
│       ├── ci.yml                    # 测试 + lint
│       └── release.yml               # 版本发布
│
├── workflow/
│   ├── Snakefile                     # 主入口
│   │
│   ├── rules/                        # Snakemake 规则（按模块）
│   │   ├── common.smk                # 共享函数、配置加载
│   │   ├── standardize.smk           # 数据标准化
│   │   ├── qc.smk                    # BUSCO QC
│   │   ├── orthology.smk             # OrthoFinder
│   │   ├── annotation.smk            # eggNOG-mapper
│   │   ├── validation.smk            # Liftoff 核验
│   │   ├── matrices.smk              # 存在/缺失矩阵
│   │   ├── expression.smk            # DESeq2
│   │   ├── reporting.smk             # 报告生成
│   │   └── ingest.smk                # SQLite 写入
│   │
│   ├── scripts/                      # 被 rule 调用的脚本
│   │   ├── run_adapter.py            # 统一 Runner
│   │   ├── ingest_db.py              # SQLite 写入脚本
│   │   └── generate_report.py        # 报告生成
│   │
│   ├── adapters/                     # 工具适配器
│   │   ├── __init__.py
│   │   ├── base.py                   # BaseAdapter, ToolSpec, RunResult
│   │   ├── orthofinder.py
│   │   ├── test_orthofinder.py
│   │   ├── eggnog.py
│   │   ├── test_eggnog.py
│   │   ├── liftoff.py
│   │   ├── test_liftoff.py
│   │   ├── busco.py
│   │   ├── test_busco.py
│   │   ├── deseq2.py
│   │   └── test_deseq2.py
│   │
│   ├── lib/                          # 共享库
│   │   ├── __init__.py
│   │   ├── errors.py                 # ErrorCode enum + ERROR_RECOVERY
│   │   ├── audit.py                  # 审计元数据收集
│   │   ├── io.py                     # 原子写入、路径工具
│   │   ├── gff.py                    # GFF/GTF 解析
│   │   ├── fasta.py                  # FASTA 处理
│   │   ├── presence_absence.py       # 矩阵生成逻辑
│   │   └── report.py                 # 报告模板渲染
│   │
│   └── envs/                         # conda 环境定义
│       ├── orthofinder.yaml
│       ├── eggnog.yaml
│       ├── liftoff.yaml
│       ├── busco.yaml
│       ├── deseq2.yaml
│       └── base.yaml                 # Python + 基础依赖
│
├── config/
│   ├── config.yaml                   # 用户配置模板
│   └── species.tsv                   # 物种清单（可选）
│
├── schemas/
│   ├── config.schema.yaml            # 配置校验 schema
│   ├── run.schema.json               # run.json 审计 schema
│   └── summary.schema.json           # summary.json schema
│
├── profiles/
│   └── local/
│       └── config.yaml               # 本地资源配置
│
├── tests/
│   ├── conftest.py                   # pytest fixtures
│   ├── integration/
│   │   ├── test_standardize_flow.py
│   │   ├── test_orthology_flow.py
│   │   ├── test_full_pipeline.py
│   │   └── test_data/
│   │       ├── mini_genome.fa
│   │       ├── mini_annotation.gff3
│   │       └── expected_outputs/
│   └── fixtures/
│       └── sample_config.yaml
│
├── resources/                        # 静态资源（可选）
│   └── report_template.html
│
├── results/                          # 运行时输出（gitignore）
│   ├── standardized/{species}/
│   ├── qc/{species}/busco/
│   ├── orthology/{species_set}/
│   ├── annotation/eggnog/{species}/
│   ├── validation/liftoff/
│   ├── matrices/
│   ├── expression/
│   ├── reports/
│   └── meta/
│       ├── compgene.db               # SQLite
│       ├── {rule}/{wildcards}.run.json
│       └── {rule}/{wildcards}.summary.json
│
└── logs/                             # 运行日志（gitignore）
    └── {rule}/{wildcards}.log
    └── {rule}/{wildcards}.jsonl
```

### Architectural Boundaries

**数据流边界：**
```
输入文件 (GFF/FASTA)
    ↓
[standardize] → standardized/{species}/ (Source of Truth)
    ↓
[qc_busco] → qc/{species}/busco/
    ↓
[orthology_infer] → orthology/{species_set}/ (需要所有物种)
    ↓
[annotation_eggnog] → annotation/eggnog/{species}/
    ↓
[matrices_generate] → matrices/*.tsv
    ↓
[validation_liftoff] → validation/liftoff/
    ↓
[ingest_db] → meta/compgene.db (Index of Record)
    ↓
[reporting_summary] → reports/summary.html
```

**模块边界（依赖方向）：**
```
lib/ ← adapters/ ← scripts/ ← rules/
  ↑
schemas/
```

- `lib/` 不依赖任何其他模块
- `adapters/` 只依赖 `lib/`
- `scripts/` 只依赖 `adapters/` 和 `lib/`
- `rules/` 只依赖 `scripts/`（通过 `script:` 指令）

**SQLite 写入边界：**
- 只有 `ingest.smk` 中的 `rule ingest_db` 写入 SQLite
- 其他 rule 只输出 `.run.json` 和 `.summary.json`

### Integration Points

**内部通信：** 纯文件（Snakemake 依赖跟踪）

**外部集成：**
| 工具 | 集成方式 | 输入 | 输出 |
|------|----------|------|------|
| OrthoFinder | subprocess | `proteins/` 目录 | `OrthoFinder/Results_*` |
| eggNOG-mapper | subprocess | FASTA | TSV |
| Liftoff | subprocess | GFF + FASTA | GFF |
| BUSCO | subprocess | FASTA | summary.txt |
| DESeq2 | Rscript | counts.tsv + metadata.tsv | DE_results.tsv |

**网络依赖：**
- NCBI API（可选，下载注释）→ `lib/ncbi.py` 处理速率限制

### File Organization Patterns

**配置文件层级：**
1. `config/config.yaml` - 用户配置
2. `profiles/local/config.yaml` - 资源覆盖
3. 命令行参数 - 最高优先级

**环境隔离：**
- 每个外部工具一个 conda env
- `workflow/envs/*.yaml` 版本锁定
- Snakemake `--use-conda` 自动激活

**测试数据：**
- `tests/test_data/` - 最小化测试数据集
- `tests/fixtures/` - 配置 fixtures
- `tests/integration/expected_outputs/` - 预期输出快照

### Requirements to Structure Mapping

| FR 类别 | 对应文件/目录 |
|---------|---------------|
| 数据准备 (FR1-5b) | `rules/standardize.smk`, `adapters/busco.py`, `lib/gff.py`, `lib/fasta.py` |
| 直系同源 (FR6-10) | `rules/orthology.smk`, `adapters/orthofinder.py` |
| 功能注释 (FR11-15) | `rules/annotation.smk`, `adapters/eggnog.py` |
| 缺失核验 (FR16-19) | `rules/validation.smk`, `adapters/liftoff.py` |
| 存在/缺失 (FR20-25) | `rules/matrices.smk`, `lib/presence_absence.py` |
| RNA-Seq (FR26-31) | `rules/expression.smk`, `adapters/deseq2.py` |
| 工作流编排 (FR32-37) | `Snakefile`, `scripts/run_adapter.py` |
| 配置管理 (FR38-41) | `schemas/config.schema.yaml`, `config/`, `profiles/` |
| 输出报告 (FR42-46) | `rules/reporting.smk`, `lib/report.py`, `resources/report_template.html` |
| 错误处理 (FR47-49) | `lib/errors.py`, `lib/audit.py` |

## Architecture Validation

### Coherence Validation

| 检查项 | 状态 | 说明 |
|--------|------|------|
| **ADR 互相兼容** | ✅ | 文件 Source of Truth + SQLite Index of Record 与 Snakemake 文件依赖模型一致 |
| **Adapter 模式与 Error 模式** | ✅ | `classify_error()` 返回 `lib/errors.py` 定义的错误码 |
| **日志模式与审计模式** | ✅ | `.run.json` 由 Runner 生成，`ingest_db` 统一写入 SQLite |
| **命名规范一致性** | ✅ | 全栈 snake_case（rule、Python、TSV、JSON） |
| **原子写入规范** | ✅ | 写入模式已定义，所有大文件先 .tmp 再 rename |

### Requirements Coverage

#### Functional Requirements（49 条）

| FR 范围 | 覆盖状态 | 架构映射 |
|---------|----------|----------|
| FR1-5b（数据准备/QC） | ✅ | `rules/standardize.smk`, `adapters/busco.py`, `lib/gff.py`, `lib/fasta.py` |
| FR6-10（直系同源） | ✅ | `rules/orthology.smk`, `adapters/orthofinder.py` |
| FR11-15（功能注释） | ✅ | `rules/annotation.smk`, `adapters/eggnog.py` |
| FR16-19（缺失核验） | ✅ | `rules/validation.smk`, `adapters/liftoff.py` |
| FR20-25（存在/缺失） | ✅ | `rules/matrices.smk`, `lib/presence_absence.py` |
| FR26-31（RNA-Seq） | ✅ | `rules/expression.smk`, `adapters/deseq2.py` |
| FR32-37（工作流编排） | ✅ | `Snakefile`, `scripts/run_adapter.py`, Snakemake 内建 |
| FR38-41（配置管理） | ✅ | `schemas/config.schema.yaml`, `config/`, `profiles/` |
| FR42-46（输出报告） | ✅ | `rules/reporting.smk`, `lib/report.py` |
| FR47-49（错误处理） | ✅ | `lib/errors.py`, ADR-003 重试策略 |

**FR 覆盖率：49/49 = 100%**

#### Non-Functional Requirements（18 条）

| NFR 范围 | 覆盖状态 | 架构支持 |
|----------|----------|----------|
| NFR1-4（性能） | ✅ | Snakemake 并行、缓存边界定义、大文件分块策略 |
| NFR5-8（可靠性） | ✅ | 三层版本追踪、断点续跑、ADR-003 重试策略 |
| NFR9-15（集成） | ✅ | Adapter 层版本检查、超时保护、速率限制 |
| NFR16-18（可维护性） | ✅ | 模块边界、日志级别配置、YAML 驱动 |

**NFR 覆盖率：18/18 = 100%**

### Implementation Readiness

| 维度 | 状态 | 说明 |
|------|------|------|
| **目录结构定义** | ✅ | 完整的项目结构已定义 |
| **依赖边界清晰** | ✅ | `lib/ ← adapters/ ← scripts/ ← rules/` |
| **数据流边界清晰** | ✅ | 输出契约明确，SQLite 写入单点 |
| **接口定义完整** | ✅ | `BaseAdapter` 8 个方法、`ErrorCode` enum |
| **命名规范可执行** | ✅ | 示例代码覆盖所有场景 |
| **测试策略定义** | ✅ | Hybrid（单元共置 + 集成集中） |

### Gap Analysis

| 潜在缺口 | 严重度 | 建议处理 |
|----------|--------|----------|
| CLI 入口点未明确 | 低 | PRD 定义了 `compgene run`，实现阶段补充 click/typer 包装 |
| NCBI 速率限制具体实现 | 低 | `lib/ncbi.py` 预留位置已有，实现阶段添加 |
| DESeq2 R 脚本调用模式 | 低 | `adapters/deseq2.py` 调用 Rscript，模式与其他 Adapter 一致 |

**结论：无阻塞性缺口，所有 Gap 可在实现阶段解决。**

### Validation Summary

| 维度 | 结果 |
|------|------|
| ADR 一致性 | ✅ PASS |
| FR 覆盖率 | ✅ 100% |
| NFR 覆盖率 | ✅ 100% |
| 实现就绪 | ✅ PASS |
| 阻塞缺口 | ✅ 无 |

## Quick Reference

### Architecture at a Glance

```
┌─────────────────────────────────────────────────────────────────┐
│                        CompGene Pipeline                        │
├─────────────────────────────────────────────────────────────────┤
│  Workflow Engine: Snakemake 9.14.8                              │
│  Language: Python 3.11+                                         │
│  Data Store: Files (SoT) + SQLite (Index)                       │
└─────────────────────────────────────────────────────────────────┘

数据流：
  输入 → [standardize] → [qc_busco] → [orthology_infer]
                                            ↓
  [reporting] ← [ingest_db] ← [matrices] ← [annotation_eggnog]
                                            ↓
                                     [validation_liftoff]

外部工具：
  OrthoFinder 2.5.x │ eggNOG-mapper 2.1.x │ Liftoff 1.6.x
  BUSCO 5.x/6.x     │ DESeq2 (R)
```

### Key Decisions Summary

| 决策 | 选择 | 原因 |
|------|------|------|
| 工作流引擎 | Snakemake 9.14.8 | Python 原生，文件依赖，bioinformatics 主流 |
| 数据架构 | 混合（文件 + SQLite） | 文件=SoT for Snakemake，SQLite=快速查询/追溯 |
| 工具集成 | Adapter 模式 | 统一接口，可测试，易扩展 |
| 错误处理 | 详细错误码 + 混合重试 | 用户友好，自动恢复 |
| 日志审计 | 双格式 + 规则级 | 人类可读 + 机器可查 |

### Critical Files Reference

| 文件 | 职责 |
|------|------|
| `workflow/Snakefile` | 主入口，include 所有 rules |
| `workflow/adapters/base.py` | BaseAdapter 接口定义 |
| `workflow/lib/errors.py` | ErrorCode enum + ERROR_RECOVERY |
| `workflow/scripts/run_adapter.py` | 统一 Runner |
| `schemas/config.schema.yaml` | 配置校验 schema |
| `results/meta/compgene.db` | SQLite 元数据库 |

### Command Quick Reference

```bash
# 安装
conda create -n compgene python=3.11 snakemake=9.14 -c conda-forge -c bioconda

# 运行
snakemake --cores 8 --use-conda

# 断点续跑
snakemake --rerun-incomplete

# 部分失败后继续
snakemake --keep-going

# 预览 DAG
snakemake --dag | dot -Tsvg > dag.svg

# 代码检查
snakemake --lint
ruff check workflow/
pytest workflow/adapters/
```

## Implementation Guidance

### Recommended Implementation Order

**Phase 1: 基础设施（Week 1）**
1. 初始化项目结构（目录、pyproject.toml、.gitignore）
2. `lib/errors.py` - ErrorCode enum + ERROR_RECOVERY dict
3. `lib/io.py` - 原子写入、路径工具
4. `lib/audit.py` - 审计元数据收集
5. `adapters/base.py` - BaseAdapter 抽象类
6. `scripts/run_adapter.py` - 统一 Runner

**Phase 2: 数据准备（Week 2）**
7. `lib/gff.py`, `lib/fasta.py` - 格式解析
8. `adapters/busco.py` + 单元测试
9. `rules/standardize.smk`, `rules/qc.smk`
10. `schemas/config.schema.yaml`

**Phase 3: 核心分析（Week 3-4）**
11. `adapters/orthofinder.py` + 单元测试
12. `adapters/eggnog.py` + 单元测试
13. `adapters/liftoff.py` + 单元测试
14. `rules/orthology.smk`, `rules/annotation.smk`, `rules/validation.smk`
15. `lib/presence_absence.py`, `rules/matrices.smk`

**Phase 4: 表达分析 + 报告（Week 5）**
16. `adapters/deseq2.py` + 单元测试
17. `rules/expression.smk`
18. SQLite schema + `rules/ingest.smk`
19. `lib/report.py`, `rules/reporting.smk`

**Phase 5: 集成测试（Week 6）**
20. `tests/integration/` - 端到端测试
21. 文档完善
22. CI/CD 配置

### First File to Create

```python
# workflow/lib/errors.py
from enum import Enum

class ErrorCode(str, Enum):
    E_INPUT_MISSING = "E_INPUT_MISSING"
    E_INPUT_FORMAT = "E_INPUT_FORMAT"
    E_TOOL_NOT_FOUND = "E_TOOL_NOT_FOUND"
    E_TOOL_VERSION = "E_TOOL_VERSION"
    E_TIMEOUT = "E_TIMEOUT"
    E_OOM = "E_OOM"
    E_NET_RATE_LIMIT = "E_NET_RATE_LIMIT"
    E_DISK_FULL = "E_DISK_FULL"
    E_NONZERO_EXIT = "E_NONZERO_EXIT"

ERROR_RECOVERY: dict[ErrorCode, tuple[bool, str]] = {
    ErrorCode.E_INPUT_MISSING: (False, "检查输入文件路径是否正确"),
    ErrorCode.E_INPUT_FORMAT: (False, "验证 GFF/FASTA 格式是否符合规范"),
    ErrorCode.E_TOOL_NOT_FOUND: (False, "激活对应的 conda 环境"),
    ErrorCode.E_TOOL_VERSION: (False, "更新工具版本至要求范围"),
    ErrorCode.E_TIMEOUT: (True, "增加超时时间或减少输入数据规模"),
    ErrorCode.E_OOM: (False, "减少 threads 或增加可用内存"),
    ErrorCode.E_NET_RATE_LIMIT: (True, "等待后自动重试"),
    ErrorCode.E_DISK_FULL: (False, "清理磁盘空间后重试"),
    ErrorCode.E_NONZERO_EXIT: (False, "查看日志了解详细错误"),
}

def get_recovery(code: ErrorCode) -> tuple[bool, str]:
    """返回 (是否可重试, 恢复建议)"""
    return ERROR_RECOVERY.get(code, (False, "未知错误，请查看日志"))
```

### Architecture Compliance Checklist

实现时务必遵循：

- [ ] Rule 命名：`{module}_{action}` 格式
- [ ] Python 命名：snake_case（函数/变量）、PascalCase（类）
- [ ] 输出路径：`results/{category}/{wildcards}/`
- [ ] 大文件写入：先 .tmp 再 rename
- [ ] 日期时间：ISO 8601 格式
- [ ] 单元测试：与被测模块共置（`test_*.py`）
- [ ] SQLite 写入：仅通过 `rule ingest_db`
- [ ] 错误码：使用 `lib/errors.py` 定义的枚举

---

_Architecture Document Complete_
_Generated: 2026-01-20_
_Author: Yuchenghuang + Claude Architect_
