# Story 3.2: Orthogroups 表生成

Status: done

## Story

As a **计算生物学研究员**,
I want **获得结构化的 orthogroups 表**,
so that **我可以查看跨物种的基因家族分组，为下游存在/缺失分析提供数据基础**。

## Acceptance Criteria

1. **AC1: Orthogroups 长格式表生成**
   - Given OrthoFinder 运行完成，`results/orthology/Orthogroups/Orthogroups.tsv` 存在
   - When 执行 `orthology_parse_orthogroups` rule
   - Then 生成 `results/orthology/orthogroups.tsv`（长格式）
   - And 包含列：`orthogroup_id`, `gene_id`, `species_id`
   - And 每个基因一行（OrthoFinder 原始格式为每行一个 orthogroup、逗号分隔基因）

2. **AC2: Orthogroup 统计表生成**
   - Given orthogroups 长格式表已生成
   - When 统计每 orthogroup 的物种分布
   - Then 生成 `results/orthology/orthogroup_stats.tsv`
   - And 包含列：`orthogroup_id`, `gene_count`, `species_count`, `is_single_copy`, `species_list`

3. **AC3: 物种重叠矩阵生成**
   - Given orthogroups 长格式表已生成
   - When 计算物种对之间共享的 orthogroup 数量
   - Then 生成 `results/orthology/species_overlap.tsv`
   - And 格式为方阵：行列均为物种名，值为共享 orthogroup 数

4. **AC4: 基因查询能力**
   - Given orthogroups.tsv 已生成
   - When 查询特定 gene_id
   - Then 可以找到其所属 orthogroup_id 以及该 orthogroup 中所有物种的基因

5. **AC5: 审计记录**
   - Given 表生成完成
   - When 收集元数据
   - Then 生成 `.run.json` 包含输入文件 checksum、输出统计、运行时间

## Tasks / Subtasks

- [x] Task 1: 创建 `workflow/lib/orthogroup_utils.py` (AC: #1, #2, #3, #4)
  - [x] 1.1 实现 `parse_orthogroups_tsv()` 解析 OrthoFinder Orthogroups.tsv → 长格式 list[dict]
  - [x] 1.2 实现 `write_orthogroups_long_format()` 写入 orthogroups.tsv（原子写入）
  - [x] 1.3 实现 `calculate_orthogroup_stats()` 统计每 orthogroup 的 gene_count、species_count、is_single_copy
  - [x] 1.4 实现 `write_orthogroup_stats()` 写入 orthogroup_stats.tsv
  - [x] 1.5 实现 `calculate_species_overlap()` 计算物种对间共享 orthogroup 方阵
  - [x] 1.6 实现 `write_species_overlap()` 写入 species_overlap.tsv

- [x] Task 2: 创建 `workflow/scripts/build_orthogroup_tables.py` (AC: #1, #2, #3, #5)
  - [x] 2.1 桥接 Snakemake → lib/orthogroup_utils.py
  - [x] 2.2 调用三个写入函数生成三个输出文件
  - [x] 2.3 生成 `.run.json` 审计记录

- [x] Task 3: 更新 `workflow/rules/orthology.smk` (AC: #1, #2, #3)
  - [x] 3.1 新增 `rule orthology_parse_orthogroups` 调用 build_orthogroup_tables.py
  - [x] 3.2 声明三个输出文件 + run_json

- [x] Task 4: 创建单元测试 `workflow/lib/test_orthogroup_utils.py` (AC: #1-5)
  - [x] 4.1 创建 fixture：模拟 Orthogroups.tsv（含多物种、逗号分隔基因、空单元格）
  - [x] 4.2 测试 `parse_orthogroups_tsv()` 解析逻辑（正常、空行、缺失物种）
  - [x] 4.3 测试 `calculate_orthogroup_stats()` 统计正确性
  - [x] 4.4 测试 `calculate_species_overlap()` 方阵对称性和值正确性
  - [x] 4.5 测试边界情况（单物种、空文件、大量 orthogroup）

## Dev Notes

### OrthoFinder Orthogroups.tsv 原始格式 [Source: 3-1-orthofinder-adapter.md]

```tsv
Orthogroup	species1	species2	species3
OG0000000	gene1, gene2	gene3	gene4, gene5, gene6
OG0000001	gene7	gene8
OG0000002		gene9	gene10
```

- 第一列：Orthogroup ID
- 后续列：该物种在该 orthogroup 中的基因，**逗号+空格分隔**（`, `）
- 空单元格表示该物种无此 orthogroup 中的基因
- 物种名来自 OrthoFinder 输入文件名（与 config 中的 species name 一致）

### 目标输出格式

**orthogroups.tsv（长格式）：**
```tsv
orthogroup_id	gene_id	species_id
OG0000000	gene1	species1
OG0000000	gene2	species1
OG0000000	gene3	species2
OG0000000	gene4	species3
OG0000000	gene5	species3
OG0000000	gene6	species3
OG0000001	gene7	species1
OG0000001	gene8	species2
```

**orthogroup_stats.tsv：**
```tsv
orthogroup_id	gene_count	species_count	is_single_copy	species_list
OG0000000	6	3	false	species1,species2,species3
OG0000001	2	2	true	species1,species2
OG0000002	2	2	true	species2,species3
```

**species_overlap.tsv（方阵）：**
```tsv
species	species1	species2	species3
species1	2	2	1
species2	2	3	2
species3	1	2	2
```
对角线值 = 该物种参与的 orthogroup 总数

### 已有代码模式参考 [Source: workflow/lib/absence_detection.py]

遵循 Epic 5 的 lib 模块模式：
- 纯函数设计，不依赖 Snakemake
- 使用 `logging.getLogger(__name__)` 记录日志
- 使用 `__version__ = "1.0.0"`
- 使用 `# ===...===` 分节注释
- 输入验证用 `CompGeneError` + `ErrorCode.E_INPUT_MISSING`
- 原子写入：写到 `.tmp` 再 `Path.rename()`
- 返回 dict 结构的统计数据

### Snakemake 规则模式 [Source: workflow/rules/orthology.smk]

新增规则命名遵循 `{module}_{action}` 格式：
```python
rule orthology_parse_orthogroups:
    input:
        orthogroups_raw="results/orthology/Orthogroups/Orthogroups.tsv",
    output:
        orthogroups="results/orthology/orthogroups.tsv",
        stats="results/orthology/orthogroup_stats.tsv",
        overlap="results/orthology/species_overlap.tsv",
        run_json="results/meta/orthology_parse_orthogroups/run.run.json",
    threads: 1
    log:
        "logs/orthology_parse_orthogroups/run.log",
    script:
        "../scripts/build_orthogroup_tables.py"
```

### 脚本桥接模式 [Source: workflow/scripts/classify_absence.py]

```python
# workflow/scripts/build_orthogroup_tables.py
from pathlib import Path
from workflow.lib.orthogroup_utils import (
    parse_orthogroups_tsv,
    write_orthogroups_long_format,
    calculate_orthogroup_stats,
    write_orthogroup_stats,
    calculate_species_overlap,
    write_species_overlap,
)
from workflow.lib.audit import create_and_write_audit
```

### 关键实现注意事项

1. **逗号分隔解析**：OrthoFinder 用 `, `（逗号+空格）分隔同一 orthogroup 内的基因。解析时需 `strip()` 每个 gene_id，处理尾部空白。

2. **空单元格处理**：某物种无此 orthogroup 的基因时，对应列为空字符串。解析时 `if cell.strip()` 过滤空值。

3. **原子写入**：所有输出文件先写 `.tmp` 再 `rename()`。[Source: architecture.md#断点续跑团队约定]

4. **不需要 Adapter**：此 Story 是对 OrthoFinder 输出的后处理，不调用外部工具，直接在 lib 模块中实现纯 Python 逻辑。

5. **输出路径契约**：`results/orthology/orthogroups.tsv` 是架构文档中定义的输出契约路径。[Source: architecture.md#输出契约]

6. **下游依赖**：此输出将被 Epic 6A（存在/缺失矩阵）和 Epic 4（功能注释汇总）使用，格式需稳定。

### 从前序 Story 3.1 学到的经验

1. **OrthoFinder 输出目录结构**：结果已由 `run_orthofinder.py` 移动到 `results/orthology/` 下，`Orthogroups/` 子目录包含所有需要的 TSV 文件。
2. **已有解析函数**：`orthofinder.py` 中的 `parse_gene_count_tsv()` 已解析 GeneCount.tsv，但 Story 3-2 需要解析的是 `Orthogroups.tsv`（包含基因 ID 而非计数）。
3. **避免 pipe buffer deadlock**：对文件操作无此问题，但保持日志重定向模式。
4. **物种名来自文件名**：OrthoFinder 用 `.fa` 文件名作物种标识，已在 `orthology_prepare_proteins` rule 中保证与 config 一致。

### Project Structure Notes

**新增文件：**
```
workflow/
├── lib/
│   ├── orthogroup_utils.py          # Orthogroups 解析与表生成（新建）
│   └── test_orthogroup_utils.py     # 单元测试（新建）
└── scripts/
    └── build_orthogroup_tables.py   # Snakemake 桥接脚本（新建）
```

**修改文件：**
```
workflow/rules/orthology.smk          # 新增 rule orthology_parse_orthogroups
```

### References

- [Source: docs/planning-artifacts/prd.md#FR7] - 系统可生成 orthogroups 表（跨物种基因家族分组）
- [Source: docs/planning-artifacts/architecture.md#ADR-001] - 数据架构：文件流 + SQLite 元数据层
- [Source: docs/planning-artifacts/architecture.md#输出契约] - `orthology/orthogroups.tsv`
- [Source: docs/planning-artifacts/architecture.md#Naming-Patterns] - Rule 命名 `{module}_{action}`
- [Source: docs/planning-artifacts/architecture.md#ADR-002] - 工具适配层模式
- [Source: docs/planning-artifacts/epics-and-stories.md#Story-3.2] - Orthogroups 表生成需求
- [Source: docs/implementation-artifacts/3-1-orthofinder-adapter.md] - 前序 Story 实现与经验
- [Source: workflow/lib/absence_detection.py] - lib 模块模式参考
- [Source: workflow/adapters/orthofinder.py#parse_gene_count_tsv] - 已有 GeneCount.tsv 解析逻辑

## Dev Agent Record

### Agent Model Used

Claude Opus 4.6 (claude-opus-4-6)

### Debug Log References

- 22 unit tests: all PASSED (0.12s)
- 771 total tests: 769 passed, 2 pre-existing failures (unrelated), zero regressions

### Completion Notes List

- ✅ Task 1: Created `orthogroup_utils.py` with 6 pure functions following absence_detection.py pattern
- ✅ Task 2: Created `build_orthogroup_tables.py` Snakemake bridge script with audit record generation
- ✅ Task 3: Added `rule orthology_parse_orthogroups` to orthology.smk with 3 outputs + run_json
- ✅ Task 4: Created 22 unit tests covering parsing, stats, overlap matrix, write functions, edge cases
- All AC satisfied: AC1 (long format), AC2 (stats), AC3 (overlap matrix), AC4 (gene query via TSV), AC5 (audit)

### Code Review Fixes

- 🔧 Finding 1 (HIGH): Fixed audit path mismatch — `meta_dir=run_json_output.parent` + `shutil.move` 兜底
- 🔧 Finding 2 (HIGH): Added `try/except CompGeneError/except Exception/finally` 错误处理模式
- 🔧 Finding 3 (MEDIUM): Added missing `shutil` and `CompGeneError` imports
- 🔧 Finding 4 (MEDIUM): Added `print()` statements for Snakemake log capture
- 🔧 Finding 5 (MEDIUM): Added `test_write_species_overlap_empty` test case
- ℹ️ Finding 6 (LOW): Noted — pure Python rule 不需 conda 指令，与项目约定一致

### File List

**New files:**
- `compgene/workflow/lib/orthogroup_utils.py` — Orthogroup parsing and table generation (6 functions)
- `compgene/workflow/lib/test_orthogroup_utils.py` — 22 unit tests
- `compgene/workflow/scripts/build_orthogroup_tables.py` — Snakemake bridge script

**Modified files:**
- `compgene/workflow/rules/orthology.smk` — Added `rule orthology_parse_orthogroups`
