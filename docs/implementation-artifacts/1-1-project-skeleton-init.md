# Story 1.1: 项目骨架初始化

Status: done

## Story

As a **计算生物学研究员**,
I want **通过标准命令初始化 CompGene 项目结构**,
so that **我可以在一致的目录布局中开始配置和运行分析**.

## Acceptance Criteria

1. **Given** 用户在空目录中
   **When** 执行项目初始化
   **Then** 生成完整目录结构（workflow/, config/, schemas/, profiles/, results/, logs/）

2. **Given** 项目初始化完成
   **When** 检查 Snakefile
   **Then** 包含 `min_version("9.0")` 和 `validate()` 调用

3. **Given** 项目初始化完成
   **When** 检查 pyproject.toml
   **Then** 定义项目元数据和依赖

## Tasks / Subtasks

- [x] Task 1: 创建完整目录结构 (AC: #1)
  - [x] 1.1 创建 `workflow/` 及子目录 (rules/, scripts/, adapters/, lib/, envs/)
  - [x] 1.2 创建 `config/` 目录
  - [x] 1.3 创建 `schemas/` 目录
  - [x] 1.4 创建 `profiles/local/` 目录
  - [x] 1.5 创建 `results/meta/` 目录
  - [x] 1.6 创建 `logs/` 目录
  - [x] 1.7 创建 `tests/integration/` 和 `tests/fixtures/` 目录
  - [x] 1.8 创建 `resources/` 目录

- [x] Task 2: 创建主 Snakefile (AC: #2)
  - [x] 2.1 添加 `min_version("9.0")` 版本检查
  - [x] 2.2 添加 `configfile` 指令
  - [x] 2.3 添加 `validate(config, schema=...)` 调用
  - [x] 2.4 添加所有 rules 的 `include` 语句
  - [x] 2.5 创建 `rule all` 定义最终目标

- [x] Task 3: 创建规则骨架文件 (AC: #1, #2)
  - [x] 3.1 创建 `rules/common.smk` (共享函数、配置加载)
  - [x] 3.2 创建空的 rules 文件: standardize.smk, qc.smk, orthology.smk, annotation.smk, validation.smk, matrices.smk, expression.smk, reporting.smk, ingest.smk

- [x] Task 4: 创建 pyproject.toml (AC: #3)
  - [x] 4.1 定义项目元数据 (name, version, description)
  - [x] 4.2 定义 Python 依赖 (snakemake, pyyaml, pandas 等)
  - [x] 4.3 配置开发工具 (ruff, pytest)

- [x] Task 5: 创建配置模板文件
  - [x] 5.1 创建 `config/config.yaml` 用户配置模板
  - [x] 5.2 创建 `profiles/local/config.yaml` 本地资源配置

- [x] Task 6: 创建 Python 包初始化
  - [x] 6.1 创建 `workflow/adapters/__init__.py`
  - [x] 6.2 创建 `workflow/lib/__init__.py`

- [x] Task 7: 创建项目基础文件
  - [x] 7.1 创建 `.gitignore` (忽略 results/, logs/, *.pyc 等)
  - [x] 7.2 创建 README.md 基本框架

## Dev Notes

### 技术栈要求 [Source: docs/planning-artifacts/architecture.md#Starter-Template-Evaluation]

- **工作流引擎**: Snakemake 9.14.8 (严格版本)
- **Python**: 3.11+
- **包管理**: conda + pip (pyproject.toml)

### 目录结构规范 [Source: docs/planning-artifacts/architecture.md#Complete-Project-Directory-Structure]

```
compgene/
├── workflow/
│   ├── Snakefile              # 主入口
│   ├── rules/                 # Snakemake rules
│   │   ├── common.smk         # 共享函数
│   │   ├── standardize.smk
│   │   ├── qc.smk
│   │   ├── orthology.smk
│   │   ├── annotation.smk
│   │   ├── validation.smk
│   │   ├── matrices.smk
│   │   ├── expression.smk
│   │   ├── reporting.smk
│   │   └── ingest.smk
│   ├── scripts/               # 被 rule 调用的脚本
│   ├── adapters/              # 工具适配器
│   ├── lib/                   # 共享库
│   └── envs/                  # conda 环境
├── config/
│   └── config.yaml
├── schemas/
│   └── config.schema.yaml     # Story 1.2 创建
├── profiles/
│   └── local/config.yaml
├── tests/
│   ├── integration/
│   └── fixtures/
├── resources/
├── results/
│   └── meta/
└── logs/
```

### Snakefile 模板 [Source: docs/planning-artifacts/architecture.md#Initialization-Command]

```python
from snakemake.utils import min_version, validate
min_version("9.0")

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

include: "rules/common.smk"
include: "rules/standardize.smk"
include: "rules/qc.smk"
include: "rules/orthology.smk"
include: "rules/annotation.smk"
include: "rules/validation.smk"
include: "rules/matrices.smk"
include: "rules/expression.smk"
include: "rules/reporting.smk"
include: "rules/ingest.smk"

rule all:
    input:
        "results/reports/summary.html"
```

### pyproject.toml 关键内容

```toml
[project]
name = "compgene"
version = "0.1.0"
description = "Comparative genomics pipeline for ortholog analysis"
requires-python = ">=3.11"
dependencies = [
    "snakemake>=9.0,<10",
    "pyyaml>=6.0",
    "pandas>=2.0",
]

[project.optional-dependencies]
dev = [
    "pytest>=7.0",
    "ruff>=0.1",
]

[tool.ruff]
line-length = 100
select = ["E", "F", "I", "N", "W"]

[tool.pytest.ini_options]
testpaths = ["tests", "workflow/adapters"]
python_files = ["test_*.py"]
```

### .gitignore 关键条目

```
# Results and logs (runtime output)
results/
logs/

# Python
__pycache__/
*.pyc
*.pyo
.pytest_cache/

# Snakemake
.snakemake/

# IDE
.vscode/
.idea/

# Conda
.conda/

# OS
.DS_Store
```

### 命名规范 [Source: docs/planning-artifacts/architecture.md#Naming-Patterns]

- **Rule 命名**: `{module}_{action}` (e.g., `orthology_infer`)
- **Python 函数/变量**: snake_case
- **Python 类**: PascalCase
- **输出路径**: `results/{category}/{wildcards}/`

### 关键约束

1. **Snakefile 位置**: 必须在 `workflow/` 目录下
2. **configfile 路径**: 相对于 Snakefile 位置 (`config/config.yaml`)
3. **schema 路径**: 相对于 Snakefile 位置 (`../schemas/config.schema.yaml`)
4. **空的 rules 文件**: 创建骨架避免 include 报错

### Project Structure Notes

- 目录结构完全遵循 Architecture 文档定义
- `results/` 和 `logs/` 需要存在但被 gitignore
- 所有 Python 包目录需要 `__init__.py`

### References

- [Source: docs/planning-artifacts/architecture.md#Complete-Project-Directory-Structure]
- [Source: docs/planning-artifacts/architecture.md#Starter-Template-Evaluation]
- [Source: docs/planning-artifacts/architecture.md#Naming-Patterns]
- [Source: docs/planning-artifacts/architecture.md#Initialization-Command]
- [Source: docs/planning-artifacts/epics.md#Story-1.1]

## Dev Agent Record

### Agent Model Used

Claude Opus 4.5 (claude-opus-4-5-20251101)

### Debug Log References

N/A - No debug issues encountered

### Completion Notes List

- Created complete directory structure matching Architecture specification
- Snakefile includes min_version("9.0") and validate() calls as required
- pyproject.toml defines project metadata with snakemake>=9.0, pyyaml, pandas dependencies
- All 10 rule skeleton files created (common.smk + 9 module rules)
- Configuration templates created for user config and local profile
- Python package __init__.py files created for adapters and lib
- .gitignore and README.md created for project basics
- Placeholder config.schema.yaml created (full schema in Story 1.2)

### File List

- compgene/.gitignore
- compgene/README.md
- compgene/pyproject.toml
- compgene/config/config.yaml
- compgene/profiles/local/config.yaml
- compgene/schemas/config.schema.yaml
- compgene/workflow/Snakefile
- compgene/workflow/adapters/__init__.py
- compgene/workflow/lib/__init__.py
- compgene/workflow/rules/common.smk
- compgene/workflow/rules/standardize.smk
- compgene/workflow/rules/qc.smk
- compgene/workflow/rules/orthology.smk
- compgene/workflow/rules/annotation.smk
- compgene/workflow/rules/validation.smk
- compgene/workflow/rules/matrices.smk
- compgene/workflow/rules/expression.smk
- compgene/workflow/rules/reporting.smk
- compgene/workflow/rules/ingest.smk
- compgene/results/.gitkeep (code review fix)
- compgene/logs/.gitkeep (code review fix)
- compgene/tests/integration/.gitkeep (code review fix)
- compgene/tests/fixtures/.gitkeep (code review fix)
- compgene/resources/.gitkeep (code review fix)
- compgene/workflow/scripts/.gitkeep (code review fix)
- compgene/workflow/envs/.gitkeep (code review fix)

## Code Review Record

### Review Model Used

Claude Opus 4.5 (claude-opus-4-5-20251101)

### Issues Found: 7 (1 HIGH, 4 MEDIUM, 2 LOW)

| ID | Severity | File | Issue | Resolution |
|----|----------|------|-------|------------|
| H1 | HIGH | workflow/Snakefile:13 | configfile path incorrect - used "config/config.yaml" instead of "../config/config.yaml" | Fixed |
| M1 | MEDIUM | pyproject.toml:55-63 | ruff.lint.select added extra rules not in spec | Fixed - reverted to spec |
| M2 | MEDIUM | pyproject.toml:34 | snakemake version range >=9.0,<10 instead of strict 9.14.8 | Fixed - pinned to ==9.14.8 |
| M3 | MEDIUM | workflow/rules/common.smk | get_species_list() silently returns [] for missing config | Fixed - added ValueError |
| M4 | MEDIUM | Multiple directories | Missing .gitkeep files for git tracking of empty dirs | Fixed - added 7 .gitkeep files |
| L1 | LOW | README.md:47-48 | Quick Start instructions assumed wrong working directory | Fixed - clarified commands |
| L2 | LOW | profiles/local/config.yaml:21 | runtime=60 missing unit annotation | Fixed - added "# minutes" comment |

### Review Outcome

PASS - All issues fixed and verified

## Change Log

- 2026-01-20: Initial implementation - created complete project skeleton with all directories, Snakefile, rules, pyproject.toml, configuration templates, and project files
- 2026-01-20: Code review completed - fixed 7 issues (1 HIGH, 4 MEDIUM, 2 LOW), added .gitkeep files, status changed to done
