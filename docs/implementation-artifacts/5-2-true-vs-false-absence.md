# Story 5.2: 映射质量分析与缺失检测

Status: done

## Story

As a **计算生物学研究员**,
I want **分析 Liftoff 映射质量并识别缺失基因**,
so that **我可以知道哪些基因在目标物种中缺失**。

## Acceptance Criteria

1. **AC1: 映射质量分类**
   - Given Liftoff 映射结果（lifted_annotation.gff3）
   - When 分析映射质量
   - Then 根据 coverage 和 identity 分类基因：
     | 映射状态 | coverage | identity | 分类 |
     |---------|----------|----------|------|
     | mapped | ≥50% | ≥50% | `present` (存在) |
     | partial | <50% | any | `uncertain` (不确定) |
     | unmapped | - | - | `missing` (缺失候选) |

2. **AC2: 缺失基因清单输出**
   - Given 分类完成
   - When 生成缺失清单
   - Then 输出 `results/liftoff/{reference}_to_{target}/missing_genes.tsv`：
     ```
     gene_id, gene_name, reference_species, target_species, liftoff_status, coverage, identity
     ```
   - And 输出 `results/liftoff/{reference}_to_{target}/gene_classification.tsv`（全部基因分类）

3. **AC3: 可配置阈值**
   - Given 用户需要调整判定严格程度
   - When 配置 config.yaml
   - Then 支持自定义阈值：
     ```yaml
     liftoff:
       min_coverage: 0.50
       min_identity: 0.50
     ```
   - And 阈值变化触发重新分类

4. **AC4: 批量比较支持**
   - Given 用户配置多对物种比较
   - When 执行分析
   - Then 为每对物种生成独立的分类结果
   - And 生成汇总表 `results/liftoff/absence_summary.tsv`

5. **AC5: 审计记录**
   - Given 分析完成
   - When Runner 收集元数据
   - Then 生成 `.run.json` 包含分类阈值、输入文件 checksum、统计摘要

## Tasks / Subtasks

- [x] Task 1: 创建缺失检测模块 (AC: #1, #3)
  - [x] 1.1 创建 `workflow/lib/absence_detection.py` 分析模块
  - [x] 1.2 实现 `classify_liftoff_genes()` 解析 GFF 并分类基因
  - [x] 1.3 实现 `extract_coverage_identity()` 从 GFF 属性提取 coverage/identity
  - [x] 1.4 实现阈值参数化（从 config 读取 min_coverage, min_identity）

- [x] Task 2: 创建 Snakemake 规则 (AC: #1, #2, #4)
  - [x] 2.1 创建 `rule liftoff_classify` 在 `workflow/rules/liftoff.smk` 中
  - [x] 2.2 输入：lifted_annotation.gff3, unmapped_features.txt
  - [x] 2.3 输出：missing_genes.tsv, gene_classification.tsv
  - [x] 2.4 实现批量比较的 expand 逻辑

- [x] Task 3: 创建汇总规则 (AC: #4)
  - [x] 3.1 创建 `rule liftoff_absence_summary`
  - [x] 3.2 聚合所有比较对的缺失统计
  - [x] 3.3 输出 `results/liftoff/absence_summary.tsv`

- [x] Task 4: 创建运行脚本 (AC: #5)
  - [x] 4.1 创建 `workflow/scripts/classify_absence.py` Snakemake 脚本
  - [x] 4.2 集成审计记录生成（使用 `create_and_write_audit()`）
  - [x] 4.3 记录分类统计（present/uncertain/missing 计数）

- [x] Task 5: 创建单元测试 (AC: #1-5)
  - [x] 5.1 创建 `workflow/lib/test_absence_detection.py`
  - [x] 5.2 测试 coverage/identity 提取（各种 GFF 属性格式）
  - [x] 5.3 测试基因分类逻辑（边界条件）
  - [x] 5.4 测试阈值参数化

## Dev Notes

### Liftoff GFF 输出属性 [Source: Story 5-1, Liftoff 文档]

**lifted_annotation.gff3 属性格式：**
```
chr1  liftoff  gene  1000  2000  .  +  .  ID=gene1;coverage=0.95;sequence_ID=NM_001234
chr1  liftoff  mRNA  1000  2000  .  +  .  ID=mrna1;Parent=gene1;coverage=0.95
```

**关键属性：**
- `coverage`: 映射覆盖度（0.0-1.0），表示参考基因有多少被映射
- `sequence_ID`: 原始参考序列 ID（用于追溯）
- 注意：identity 不是直接属性，需要从 minimap2 比对中推断或使用 coverage 近似

**unmapped_features.txt 格式：**
- 每行一个未映射的 gene ID
- 这些基因直接标记为 `missing`

### 分类逻辑 [Source: epics.md#Epic5-Story5.2]

```python
def classify_gene(coverage: float, identity: float,
                  min_coverage: float = 0.5, min_identity: float = 0.5) -> str:
    """
    分类基因状态。

    Args:
        coverage: Liftoff 报告的覆盖度 (0.0-1.0)
        identity: 序列相似度 (0.0-1.0)，如无则使用 coverage
        min_coverage: 最小覆盖度阈值
        min_identity: 最小相似度阈值

    Returns:
        'present': 成功映射
        'uncertain': 部分映射
        'missing': 未映射
    """
    if coverage >= min_coverage and identity >= min_identity:
        return 'present'
    elif coverage > 0:
        return 'uncertain'
    else:
        return 'missing'
```

### 输出文件格式

**gene_classification.tsv：**
```
gene_id	gene_name	reference_species	target_species	status	coverage	identity	source_gene
GENE001	BRCA1	human	mouse_lemur	present	0.95	0.92	NM_007294
GENE002	TP53	human	mouse_lemur	uncertain	0.42	0.38	NM_000546
GENE003	FOXP2	human	mouse_lemur	missing	0.0	0.0	NM_014491
```

**missing_genes.tsv：**
```
gene_id	gene_name	reference_species	target_species	liftoff_status	coverage	identity
GENE003	FOXP2	human	mouse_lemur	unmapped	0.0	0.0
GENE004	NOTCH1	human	mouse_lemur	partial	0.35	0.30
```

**absence_summary.tsv：**
```
reference	target	total_genes	present	uncertain	missing	present_rate	missing_rate
human	mouse_lemur	20000	18500	800	700	0.925	0.035
human	ring_tailed	20000	18200	900	900	0.910	0.045
```

### 从 Story 5-1 学到的经验

1. **统计需要区分 gene vs feature**：5-1 code review 发现统计所有 GFF 行而非 gene 行会导致误导，本 Story 必须只统计 gene 类型
2. **处理压缩文件**：输入 GFF 可能是 gzip 压缩的，使用 `decompress_if_gzipped()` 或直接读取
3. **审计记录必须生成**：使用 `create_and_write_audit()` 生成 `.run.json`
4. **字段命名一致性**：使用 `lifted_genes`/`unmapped_genes` 而非 `lifted_features`

### 依赖关系

**输入依赖（来自 Story 5-1）：**
- `results/liftoff/{reference}_to_{target}/lifted_annotation.gff3`
- `results/liftoff/{reference}_to_{target}/unmapped_features.txt`
- `results/liftoff/{reference}_to_{target}/liftoff_stats.tsv`

**配置依赖：**
```yaml
# config/config.yaml
liftoff:
  comparisons:
    - reference: human
      targets: [mouse_lemur, ring_tailed_lemur]
  min_coverage: 0.50
  min_identity: 0.50
```

### Project Structure Notes

**新增文件：**
```
workflow/
├── lib/
│   ├── absence_detection.py       # 缺失检测核心逻辑（新建）
│   └── test_absence_detection.py  # 单元测试（新建）
└── scripts/
    └── classify_absence.py        # Snakemake 脚本（新建）
```

**修改文件：**
```
workflow/rules/liftoff.smk         # 添加 liftoff_classify, liftoff_absence_summary 规则
```

### 架构合规性检查

| 要求 | 合规性 | 说明 |
|------|--------|------|
| 文件产物作为 Source of Truth | ✅ | TSV 输出，Snakemake 依赖跟踪 |
| 审计记录 | ✅ | 生成 .run.json |
| 错误处理 | ✅ | 使用 ErrorCode，CompGeneError |
| 命名规范 | ✅ | snake_case 函数，results/{category}/ 路径 |
| 单元测试 | ✅ | lib/ 模块有配套测试 |

### References

- [Source: docs/planning-artifacts/epics.md#Story-5.2] - 映射质量分析与缺失检测
- [Source: docs/planning-artifacts/prd.md#FR17] - 系统可识别"真缺失"vs"注释假缺失"候选
- [Source: docs/planning-artifacts/prd.md#FR18] - 系统可输出核验后的"真缺失候选清单"
- [Source: docs/planning-artifacts/architecture.md#Module-8] - Validation & Cross-check 模块
- [Source: docs/implementation-artifacts/5-1-liftoff-adapter.md] - Liftoff Adapter 实现和 code review 经验

## Dev Agent Record

### Agent Model Used

Claude Opus 4.5 (claude-opus-4-5-20251101)

### Debug Log References

N/A - No issues encountered during implementation

### Completion Notes List

1. **Task 1 Complete**: Created `workflow/lib/absence_detection.py` with:
   - `classify_gene()`: Classification function with configurable thresholds
   - `extract_coverage_identity()`: GFF attribute parser with percentage normalization
   - `classify_liftoff_genes()`: Main function to process Liftoff outputs
   - `calculate_classification_stats()`: Statistics aggregation
   - `filter_missing_genes()`: Filter for missing/uncertain genes
   - Follows red-green-refactor cycle (tests written first)

2. **Task 2 Complete**: Added Snakemake rules to `workflow/rules/liftoff.smk`:
   - `rule liftoff_classify`: Classifies genes for a single comparison
   - `rule liftoff_classify_all`: Triggers all comparisons
   - Added `get_classification_targets()` helper function
   - Updated module header to include Story 5.2 outputs

3. **Task 3 Complete**: Added `rule liftoff_absence_summary` in `workflow/rules/liftoff.smk`:
   - Aggregates classification stats across all comparisons
   - Outputs `results/liftoff/absence_summary.tsv` with present/missing rates

4. **Task 4 Complete**: Created `workflow/scripts/classify_absence.py`:
   - Snakemake script bridging rules to absence_detection module
   - Generates audit records using `create_and_write_audit()`
   - Writes both `gene_classification.tsv` and `missing_genes.tsv`

5. **Task 5 Complete**: Created `workflow/lib/test_absence_detection.py` with 26 tests:
   - 9 tests for `classify_gene()` function (boundary conditions, custom thresholds)
   - 5 tests for `extract_coverage_identity()` (various attribute formats)
   - 5 tests for `classify_liftoff_genes()` (GFF parsing, unmapped handling)
   - 2 tests for `calculate_classification_stats()`
   - 2 tests for `filter_missing_genes()`
   - 3 tests for edge cases (percentage coverage, extra attributes, gzip)

### File List

**New Files:**
- `workflow/lib/absence_detection.py` - Core absence detection module
- `workflow/lib/test_absence_detection.py` - Unit tests (26 tests, all pass)
- `workflow/scripts/classify_absence.py` - Snakemake script

**Modified Files:**
- `workflow/rules/liftoff.smk` - Added liftoff_classify, liftoff_classify_all, liftoff_absence_summary rules

## Change Log

- 2026-02-05: Story file created by create-story workflow
- 2026-02-05: Implementation complete - all 5 tasks done, 26 tests passing
- 2026-02-05: Code review complete - 8 issues fixed (4 HIGH, 4 MEDIUM), 29 tests passing

### Code Review Fixes Applied

**HIGH Issues Fixed:**
- H1: Removed unused `open_gff_file` import from absence_detection.py
- H2: Added warning logging when unmapped_file doesn't exist
- H3: Added `conda:` directive to `liftoff_classify` rule
- H4: Removed unused `audit_summary` variable, fixed command tracing in audit

**MEDIUM Issues Fixed:**
- M1: Added specific `CompGeneError` exception handling in classify_absence.py
- M2: Added test for missing unmapped file behavior (`test_missing_unmapped_file_logs_warning`)
- M3: Refactored `liftoff_absence_summary` to use new `read_classification_stats_from_tsv()` function
- M4: Replaced hardcoded tool_version with dynamic `absence_detection.__version__`

**LOW Issues Fixed:**
- L1: Removed unused `tempfile` import from test_absence_detection.py
