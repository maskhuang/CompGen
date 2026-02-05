# Story 5.3: 迁移注释输出

Status: done

## Story

As a **计算生物学研究员**,
I want **获得目标物种的迁移注释文件**,
so that **我可以用于下游分析**。

## Acceptance Criteria

1. **AC1: 增强的迁移注释文件**
   - Given Liftoff 成功映射的基因
   - When 生成迁移注释
   - Then 输出 `results/liftoff/{reference}_to_{target}/lifted_annotation_enhanced.gff3`
   - And 每个特征添加属性 `liftoff_coverage=X;liftoff_identity=Y;source_gene=Z`

2. **AC2: 来源追溯属性**
   - Given 需要追溯来源
   - When 查询迁移基因
   - Then 可通过 GFF3 属性追溯：
     - `source_gene`: 原始参考基因 ID
     - `reference_species`: 参考物种名称
     - `liftoff_coverage`: 映射覆盖度 (0.0-1.0)
     - `liftoff_identity`: 映射相似度 (0.0-1.0)
     - `liftoff_status`: 分类状态 (present/uncertain/missing)

3. **AC3: 仅包含成功映射基因**
   - Given 需要干净的注释用于下游分析
   - When 输出增强注释文件
   - Then 仅包含 `status == 'present'` 的基因（coverage ≥50%, identity ≥50%）
   - And 可配置是否包含 uncertain 基因

4. **AC4: 批量比较支持**
   - Given 用户配置多对物种比较
   - When 执行增强注释生成
   - Then 为每对物种生成独立的 `lifted_annotation_enhanced.gff3`

5. **AC5: 审计记录**
   - Given 增强注释生成完成
   - When Runner 收集元数据
   - Then 生成 `.run.json` 包含处理参数、输入文件 checksum、基因计数统计

## Tasks / Subtasks

- [x] Task 1: 创建注释增强模块 (AC: #1, #2)
  - [x] 1.1 创建 `workflow/lib/annotation_enhance.py` 模块
  - [x] 1.2 实现 `enhance_gff_with_liftoff_attrs()` 函数，添加 liftoff 属性
  - [x] 1.3 实现 `filter_by_status()` 函数，根据分类状态过滤
  - [x] 1.4 实现 `load_classification_data()` 关联分类结果与 GFF

- [x] Task 2: 创建 Snakemake 规则 (AC: #1, #3, #4)
  - [x] 2.1 在 `workflow/rules/liftoff.smk` 添加 `rule liftoff_enhance`
  - [x] 2.2 输入：lifted_annotation.gff3, gene_classification.tsv
  - [x] 2.3 输出：lifted_annotation_enhanced.gff3
  - [x] 2.4 添加 `include_uncertain` 参数支持

- [x] Task 3: 创建运行脚本 (AC: #5)
  - [x] 3.1 创建 `workflow/scripts/enhance_annotation.py` Snakemake 脚本
  - [x] 3.2 集成审计记录生成（使用 `create_and_write_audit()`）
  - [x] 3.3 记录处理统计（输入基因数、输出基因数、过滤基因数）

- [x] Task 4: 创建 all 规则 (AC: #4)
  - [x] 4.1 添加 `rule liftoff_enhance_all` 触发所有比较
  - [x] 4.2 添加 `get_enhance_targets()` 辅助函数

- [x] Task 5: 创建单元测试 (AC: #1-5)
  - [x] 5.1 创建 `workflow/lib/test_annotation_enhance.py`
  - [x] 5.2 测试 GFF 属性增强逻辑
  - [x] 5.3 测试状态过滤逻辑
  - [x] 5.4 测试边界条件（无匹配基因、空输入）

## Dev Notes

### 设计说明

本 Story 是 Epic 5 的最后一个 Story，完成 Liftoff 工作流的闭环。核心功能是将 Story 5.2 的分类结果（present/uncertain/missing）合并回 GFF3 文件，生成增强版的迁移注释，便于下游分析工具直接使用。

### 依赖关系

**输入依赖（来自 Story 5-1, 5-2）：**
- `results/liftoff/{reference}_to_{target}/lifted_annotation.gff3` (Story 5-1)
- `results/liftoff/{reference}_to_{target}/gene_classification.tsv` (Story 5-2)

**输出产物：**
- `results/liftoff/{reference}_to_{target}/lifted_annotation_enhanced.gff3`
- `results/meta/liftoff_enhance/reference={reference}_target={target}.run.json`

### GFF3 属性增强格式

**输入 GFF3（来自 Liftoff）：**
```
chr1  liftoff  gene  1000  2000  .  +  .  ID=gene1;coverage=0.95;sequence_ID=NM_001234
chr1  liftoff  mRNA  1000  2000  .  +  .  ID=mrna1;Parent=gene1;coverage=0.95
```

**输出 GFF3（增强后）：**
```
chr1  liftoff  gene  1000  2000  .  +  .  ID=gene1;coverage=0.95;sequence_ID=NM_001234;liftoff_coverage=0.95;liftoff_identity=0.95;liftoff_status=present;source_gene=NM_001234;reference_species=human
chr1  liftoff  mRNA  1000  2000  .  +  .  ID=mrna1;Parent=gene1;coverage=0.95;liftoff_coverage=0.95;liftoff_identity=0.95;liftoff_status=present;source_gene=NM_001234;reference_species=human
```

### 核心实现逻辑

```python
# workflow/lib/annotation_enhance.py

def enhance_gff_with_liftoff_attrs(
    gff_path: Path,
    classification_path: Path,
    output_path: Path,
    reference_species: str,
    include_uncertain: bool = False
) -> dict:
    """
    增强 GFF3 文件，添加 Liftoff 分类属性。

    Args:
        gff_path: 原始 lifted_annotation.gff3 路径
        classification_path: gene_classification.tsv 路径
        output_path: 输出 GFF3 路径
        reference_species: 参考物种名称
        include_uncertain: 是否包含 uncertain 状态的基因

    Returns:
        处理统计字典：
        - input_genes: 输入基因数
        - output_genes: 输出基因数
        - filtered_genes: 过滤掉的基因数
    """
    # 1. 读取分类结果，构建 gene_id -> classification 映射
    # 2. 遍历 GFF3，为每个 gene 特征添加 liftoff 属性
    # 3. 根据 include_uncertain 过滤非 present 基因
    # 4. 子特征（mRNA, CDS 等）继承父基因的属性
    # 5. 原子写入输出文件
```

### 配置支持

```yaml
# config/config.yaml
liftoff:
  comparisons:
    - reference: human
      targets: [mouse_lemur, ring_tailed_lemur]
  min_coverage: 0.50
  min_identity: 0.50
  include_uncertain: false  # 新增：是否包含 uncertain 基因
```

### 从 Story 5-1, 5-2 学到的经验

1. **统计需要区分 gene vs feature**：只对 gene 类型行计数和处理
2. **子特征继承父属性**：mRNA/CDS/exon 等子特征需要继承 gene 的 liftoff 属性
3. **原子写入**：大文件先写 .tmp 再 rename
4. **审计记录必须生成**：使用 `create_and_write_audit()` 生成 `.run.json`
5. **处理压缩文件**：输入 GFF 可能是 gzip 压缩的，使用 `parse_gff3()` 自动处理

### gene_classification.tsv 格式 [Source: Story 5-2]

```tsv
gene_id	gene_name	reference_species	target_species	status	coverage	identity	source_gene
GENE001	BRCA1	human	mouse_lemur	present	0.95	0.92	NM_007294
GENE002	TP53	human	mouse_lemur	uncertain	0.42	0.38	NM_000546
GENE003	FOXP2	human	mouse_lemur	missing	0.0	0.0	NM_014491
```

### Project Structure Notes

**新增文件：**
```
workflow/
├── lib/
│   ├── annotation_enhance.py       # 注释增强核心逻辑（新建）
│   └── test_annotation_enhance.py  # 单元测试（新建）
└── scripts/
    └── enhance_annotation.py       # Snakemake 脚本（新建）
```

**修改文件：**
```
workflow/rules/liftoff.smk         # 添加 liftoff_enhance, liftoff_enhance_all 规则
```

### Snakemake 规则设计

```python
# workflow/rules/liftoff.smk (新增)

rule liftoff_enhance:
    """
    Enhance lifted annotation with classification attributes.

    Adds liftoff_coverage, liftoff_identity, liftoff_status, source_gene,
    and reference_species attributes to the GFF3.

    Source: Story 5.3 - 迁移注释输出
    """
    input:
        lifted_gff="results/liftoff/{reference}_to_{target}/lifted_annotation.gff3",
        classification="results/liftoff/{reference}_to_{target}/gene_classification.tsv",
    output:
        enhanced_gff="results/liftoff/{reference}_to_{target}/lifted_annotation_enhanced.gff3",
        run_json="results/meta/liftoff_enhance/reference={reference}_target={target}.run.json",
    params:
        include_uncertain=lambda wildcards: config.get("liftoff", {}).get("include_uncertain", False),
    threads: 1
    log:
        "logs/liftoff_enhance/{reference}_to_{target}.log"
    conda:
        "../envs/liftoff.yaml"
    script:
        "../scripts/enhance_annotation.py"


rule liftoff_enhance_all:
    """
    Run enhancement for all Liftoff comparisons.
    """
    input:
        get_enhance_targets()
    output:
        touch("results/liftoff/.enhance_complete")
```

### 架构合规性检查

| 要求 | 合规性 | 说明 |
|------|--------|------|
| 文件产物作为 Source of Truth | ✅ | GFF3 输出，Snakemake 依赖跟踪 |
| 审计记录 | ✅ | 生成 .run.json |
| 错误处理 | ✅ | 使用 ErrorCode，CompGeneError |
| 命名规范 | ✅ | snake_case 函数，results/{category}/ 路径 |
| 单元测试 | ✅ | lib/ 模块有配套测试 |
| 原子写入 | ✅ | 先写 .tmp 再 rename |

### 测试策略

**单元测试（test_annotation_enhance.py）：**
1. 测试属性增强：验证 liftoff_* 属性正确添加
2. 测试状态过滤：只输出 present 基因（默认）
3. 测试 include_uncertain：包含 uncertain 状态基因
4. 测试子特征继承：mRNA/CDS 继承父 gene 属性
5. 测试无匹配基因：分类结果中无对应 gene
6. 测试空输入：空 GFF 或空分类文件

### References

- [Source: docs/planning-artifacts/epics.md#Story-5.3] - 迁移注释输出
- [Source: docs/planning-artifacts/prd.md#FR19] - 系统可生成注释补全后的对照注释文件
- [Source: docs/planning-artifacts/architecture.md#Module-8] - Validation & Cross-check 模块
- [Source: docs/implementation-artifacts/5-1-liftoff-adapter.md] - Liftoff Adapter 实现
- [Source: docs/implementation-artifacts/5-2-true-vs-false-absence.md] - 映射质量分析与缺失检测实现

## Dev Agent Record

### Agent Model Used

Claude Opus 4.5 (claude-opus-4-5-20251101)

### Debug Log References

N/A - No issues encountered during implementation

### Completion Notes List

1. **Task 1 Complete**: Created `workflow/lib/annotation_enhance.py` with:
   - `load_classification_data()`: Loads gene classification TSV into a dict
   - `format_gff3_attributes()`: Formats attributes dict to GFF3 string
   - `enhance_gff_line()`: Adds liftoff_* attributes to a GFF line
   - `filter_by_status()`: Determines allowed statuses for filtering
   - `enhance_gff_with_liftoff_attrs()`: Main function processing GFF with classification
   - Module version `__version__ = "1.0.0"` for audit records
   - Atomic file writing (write to .tmp, then rename)

2. **Task 5 Complete**: Created `workflow/lib/test_annotation_enhance.py` with 20 tests:
   - 3 tests for `load_classification_data()`
   - 2 tests for `format_gff3_attributes()`
   - 3 tests for `enhance_gff_line()`
   - 2 tests for `filter_by_status()`
   - 8 tests for `enhance_gff_with_liftoff_attrs()` (main function)
   - 2 tests for edge cases (percentage normalization, extra attributes)
   - Follows red-green-refactor cycle (tests written before implementation)

3. **Task 3 Complete**: Created `workflow/scripts/enhance_annotation.py`:
   - Snakemake script bridging rules to annotation_enhance module
   - Generates audit records using `create_and_write_audit()`
   - Records processing statistics (input/output/filtered genes)
   - Handles CompGeneError exceptions with proper error codes

4. **Task 2 & 4 Complete**: Added to `workflow/rules/liftoff.smk`:
   - `rule liftoff_enhance`: Enhances GFF with classification attributes
   - `rule liftoff_enhance_all`: Triggers all comparisons
   - `get_enhance_targets()` helper function
   - Updated module header to include Story 5.3 outputs

### File List

**New Files:**
- `workflow/lib/annotation_enhance.py` - Core annotation enhancement module
- `workflow/lib/test_annotation_enhance.py` - Unit tests (20 tests, all pass)
- `workflow/scripts/enhance_annotation.py` - Snakemake script

**Modified Files:**
- `workflow/rules/liftoff.smk` - Added liftoff_enhance, liftoff_enhance_all rules, get_enhance_targets() function

## Change Log

- 2026-02-05: Story file created by create-story workflow
- 2026-02-05: Implementation complete - all 5 tasks done, 20 tests passing
- 2026-02-05: Code review fixes applied - 22 tests passing:
  - H2: Removed unused `gene_children` variable (line 302)
  - H4: Added coverage/identity normalization for values >1.0 (percentage to decimal)
  - H1: Added gzip input test (`test_gzip_input_support`)
  - M3: Added special character test for source_gene (`test_special_characters_in_source_gene`)
  - Updated `test_percentage_coverage_normalization` to verify actual normalization
