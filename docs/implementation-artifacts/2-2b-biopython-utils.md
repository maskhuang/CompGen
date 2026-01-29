# Story 2.2b: Biopython 工具封装

Status: done

## Story

As a **计算生物学研究员**,
I want **使用 Biopython 处理复杂序列操作**,
so that **我可以处理非标准遗传密码、多格式 ID 解析和生成审计摘要**。

## Acceptance Criteria

1. **AC1: 多遗传密码表支持**
   - Given 需要翻译线粒体或其他非标准遗传密码序列
   - When 调用翻译函数
   - Then 使用 Biopython 的 CodonTable 支持 NCBI 定义的所有遗传密码表
   - And 支持常用密码表：Standard (1), Vertebrate Mitochondrial (2), Yeast Mitochondrial (3)

2. **AC2: 多格式 ID 解析**
   - Given 来自不同来源的 FASTA 文件（UniProt、NCBI、OrthoFinder 输出）
   - When 提取序列 ID
   - Then 正确解析各种 ID 格式：
     - UniProt: `sp|P12345|GENE_HUMAN` → `P12345`
     - NCBI: `gi|123|ref|NP_001234.1|` → `NP_001234.1`
     - OrthoFinder: `species|gene_id` → `gene_id`
     - Simple: `gene_id description` → `gene_id`

3. **AC3: ID 一致性检查（OrthoFinder/eggNOG）**
   - Given 标准化蛋白 FASTA 和 OrthoFinder/eggNOG 输出
   - When 执行一致性检查
   - Then 比较 ID 集合并报告：
     - 共有 ID 数量
     - 仅在输入中的 ID（未被处理）
     - 仅在输出中的 ID（可能格式变化）
     - ID 映射关系（输入 ID → 输出 ID）
   - And 输出 consistency_report.json 包含差异详情

4. **AC4: 审计摘要生成（summary.json）**
   - Given 完成数据标准化后的 FASTA 文件
   - When 执行 `generate_summary` 命令或 Snakemake rule
   - Then 在同目录生成 `summary.json` 包含：
     - sequence_count: 序列数量
     - total_length: 总长度
     - mean_length: 平均长度
     - min_length / max_length: 长度范围
     - n50: N50 值
     - gc_content: GC 含量（核酸序列，蛋白为 null）
     - n_content: N/X 碱基比例
     - length_histogram: 长度分布直方图数据
   - And 支持批量生成所有标准化 FASTA 的 summary.json

5. **AC5: 依赖检测与降级**
   - Given Biopython 未安装
   - When 调用 bio_utils 模块
   - Then 返回 E_TOOL_NOT_FOUND 并提供安装建议
   - And 核心功能可使用现有 fasta.py 模块降级运行

6. **AC6: 与现有模块集成**
   - Given 现有 fasta.py 和 standardize.py 模块
   - When 调用 bio_utils 功能
   - Then 无缝集成，不影响现有功能
   - And 共享 CompGeneError 错误处理体系

## Tasks / Subtasks

- [x] Task 1: 创建 lib/bio_utils.py 模块 (AC: #1, #2, #5, #6)
  - [x] 1.1 实现依赖检测 `check_biopython_available()`
  - [x] 1.2 实现多密码表翻译 `translate_with_table(sequence, table_id)`
  - [x] 1.3 实现 ID 解析器 `parse_fasta_id(header, format_hint=None)`
  - [x] 1.4 实现批量 ID 提取 `extract_ids(fasta_path, format_hint=None)` → `set[str]`

- [x] Task 2: 实现 ID 一致性检查功能 (AC: #3)
  - [x] 2.1 实现 `compare_id_sets(ids_a, ids_b)` 返回 ComparisonResult dataclass
  - [x] 2.2 实现 `check_orthofinder_consistency(input_fasta, orthofinder_output_dir)`
  - [x] 2.3 实现 `check_eggnog_consistency(input_fasta, eggnog_annotations)`
  - [x] 2.4 实现 `write_consistency_report(result, output_path)` 写入 JSON

- [x] Task 3: 实现审计摘要功能 (AC: #4)
  - [x] 3.1 实现 `generate_audit_summary(fasta_path)` 返回 dict
  - [x] 3.2 实现 `write_summary_json(fasta_path, output_path=None)` 写入 summary.json
  - [x] 3.3 支持核酸和蛋白序列的差异化统计（GC vs 氨基酸组成）
  - [x] 3.4 实现长度直方图计算 `compute_length_histogram(lengths, bins=None)`

- [x] Task 4: 创建 Snakemake 规则 (AC: #4)
  - [x] 4.1 创建 rules/audit.smk
  - [x] 4.2 实现 rule generate_fasta_summary（单文件）
  - [x] 4.3 实现 rule generate_all_summaries（批量）

- [x] Task 5: 创建单元测试 (AC: #1-6)
  - [x] 5.1 创建 lib/test_bio_utils.py
  - [x] 5.2 测试各遗传密码表翻译
  - [x] 5.3 测试各种 ID 格式解析
  - [x] 5.4 测试 ID 一致性检查
  - [x] 5.5 测试审计摘要生成
  - [x] 5.6 测试 Biopython 缺失时的降级行为

- [x] Task 6: 更新依赖配置 (AC: #5)
  - [x] 6.1 更新 pyproject.toml 添加 biopython>=1.81
  - [x] 6.2 验证安装和导入

## Dev Notes

### 模块设计 [Source: 现有架构模式]

```python
# workflow/lib/bio_utils.py

from pathlib import Path
from typing import Optional
import json

# 延迟导入 Biopython，支持降级
_BIOPYTHON_AVAILABLE: Optional[bool] = None

def check_biopython_available() -> bool:
    """检查 Biopython 是否可用"""
    global _BIOPYTHON_AVAILABLE
    if _BIOPYTHON_AVAILABLE is None:
        try:
            from Bio import SeqIO
            from Bio.Seq import Seq
            _BIOPYTHON_AVAILABLE = True
        except ImportError:
            _BIOPYTHON_AVAILABLE = False
    return _BIOPYTHON_AVAILABLE
```

### 遗传密码表映射 [Source: NCBI Genetic Codes]

```python
CODON_TABLE_NAMES = {
    1: "Standard",
    2: "Vertebrate Mitochondrial",
    3: "Yeast Mitochondrial",
    4: "Mold Mitochondrial",
    5: "Invertebrate Mitochondrial",
    6: "Ciliate Nuclear",
    9: "Echinoderm Mitochondrial",
    10: "Euplotid Nuclear",
    11: "Bacterial",
    12: "Alternative Yeast Nuclear",
    13: "Ascidian Mitochondrial",
    14: "Alternative Flatworm Mitochondrial",
    15: "Blepharisma Macronuclear",
    16: "Chlorophycean Mitochondrial",
    21: "Trematode Mitochondrial",
    22: "Scenedesmus obliquus Mitochondrial",
    23: "Thraustochytrium Mitochondrial",
}
```

### ID 格式解析规则

| 来源 | 格式示例 | 解析规则 |
|------|----------|----------|
| UniProt | `sp\|P12345\|GENE_HUMAN` | 分割 `\|`，取第2部分 |
| NCBI RefSeq | `ref\|NP_001234.1\|` | 分割 `\|`，取 ref 后部分 |
| NCBI GI | `gi\|123456\|ref\|NP_001\|` | 分割 `\|`，取 ref 后部分 |
| GenBank | `gb\|AAA12345.1\|` | 分割 `\|`，取第2部分 |
| OrthoFinder | `species\|gene_id` | 分割 `\|`，取最后部分 |
| Simple | `gene_id description text` | 分割空白，取第一部分 |

### ID 一致性检查 [核心用例 - OrthoFinder/eggNOG]

**使用场景：**
1. 标准化蛋白 FASTA → OrthoFinder 输入/输出 ID 对比
2. 标准化蛋白 FASTA → eggNOG-mapper 注释结果 ID 对比
3. 检测 ID 格式变化（如 `species|gene` vs `gene`）

**ComparisonResult dataclass：**
```python
@dataclass
class ComparisonResult:
    """ID 集合比较结果"""
    source_a: str              # 来源 A 描述
    source_b: str              # 来源 B 描述
    total_a: int               # A 中 ID 总数
    total_b: int               # B 中 ID 总数
    common: int                # 共有 ID 数量
    only_in_a: set[str]        # 仅在 A 中的 ID
    only_in_b: set[str]        # 仅在 B 中的 ID
    match_rate: float          # 匹配率 = common / total_a

    def is_consistent(self, threshold: float = 0.95) -> bool:
        """检查是否达到一致性阈值"""
        return self.match_rate >= threshold
```

**consistency_report.json 格式：**
```json
{
  "generated_at": "2026-01-29T12:00:00Z",
  "source_a": "results/standardized/mmur/proteins.longest.fa.gz",
  "source_b": "results/orthology/OrthoFinder/Orthogroups.tsv",
  "comparison": {
    "total_a": 12345,
    "total_b": 12340,
    "common": 12335,
    "only_in_a_count": 10,
    "only_in_b_count": 5,
    "match_rate": 0.9992
  },
  "only_in_a": ["gene_001", "gene_002", "..."],
  "only_in_b": ["gene_x", "..."],
  "status": "PASS",
  "threshold": 0.95
}
```

### summary.json 格式（审计摘要）

**输出位置：** 与 FASTA 同目录，如：
- `results/standardized/{species}/proteins.longest.fa.gz`
- `results/standardized/{species}/proteins.longest.summary.json`

```json
{
  "file_path": "proteins.longest.fa.gz",
  "file_size_bytes": 1234567,
  "generated_at": "2026-01-29T12:00:00Z",
  "sequence_type": "protein",
  "statistics": {
    "sequence_count": 12345,
    "total_length": 5678901,
    "mean_length": 460.1,
    "median_length": 423,
    "min_length": 50,
    "max_length": 8234,
    "n50": 512,
    "gc_content": null,
    "n_content": 0.0,
    "ambiguous_count": 0,
    "length_histogram": {
      "bins": [0, 100, 200, 500, 1000, 2000, 5000, 10000],
      "counts": [123, 456, 789, 234, 567, 89, 12]
    }
  },
  "id_format": {
    "detected": "species|gene_id",
    "sample_ids": ["mmur|MMUR_001", "mmur|MMUR_002", "mmur|MMUR_003"]
  }
}
```

### Snakemake 规则 [Source: architecture.md#Naming-Patterns]

```python
# rules/audit.smk

rule generate_fasta_summary:
    """为单个 FASTA 生成 summary.json"""
    input:
        fasta="{output_dir}/standardized/{species}/{file}.fa.gz"
    output:
        summary="{output_dir}/standardized/{species}/{file}.summary.json"
    run:
        from workflow.lib.bio_utils import write_summary_json
        write_summary_json(Path(input.fasta), Path(output.summary))

rule generate_all_summaries:
    """批量生成所有标准化 FASTA 的 summary.json"""
    input:
        expand(
            "{output_dir}/standardized/{species}/{file}.summary.json",
            output_dir=get_output_dir(),
            species=get_local_species(),
            file=["genome", "proteins.longest"]
        )
    output:
        touch("{output_dir}/standardized/.summaries_complete")
```

### 与现有模块关系

```
workflow/lib/
├── fasta.py           # 核心 FASTA 解析（保留）
├── standardize.py     # 标准化逻辑（保留）
├── bio_utils.py       # Biopython 封装（新建）
│   ├── 依赖: fasta.py（降级时使用）
│   └── 增强: 密码表、ID解析、审计、一致性检查
└── test_bio_utils.py  # 单元测试（新建）

workflow/rules/
├── standardize.smk    # 已有
└── audit.smk          # 新建 - summary.json 生成
```

### 降级策略

当 Biopython 不可用时：
- `translate_with_table()`: 对于标准密码表 (table_id=1) 使用 standardize.py 的 CODON_TABLE
- `parse_fasta_id()`: 使用正则表达式实现基本解析
- `generate_audit_summary()`: 使用 fasta.py 的 get_sequence_stats()

### References

- [Source: docs/planning-artifacts/architecture.md#Scope-Boundaries] - 模块边界
- [Source: docs/planning-artifacts/architecture.md#ADR-003] - 错误处理规范
- [Source: docs/implementation-artifacts/2-1-gff-fasta-parser.md] - FASTA 解析模块
- [Source: docs/implementation-artifacts/2-2-local-data-standardization.md] - 标准化模块
- [NCBI Genetic Codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) - 遗传密码表参考

## Dev Agent Record

### Agent Model Used

Claude Opus 4.5 (claude-opus-4-5-20251101)

### Debug Log References

N/A - No debug issues encountered

### Completion Notes List

- Task 1: Created bio_utils.py with 600+ lines of Biopython-enhanced utilities
- Task 2: Implemented ComparisonResult dataclass and consistency checking functions
- Task 3: Implemented audit summary generation with JSON output
- Task 4: Created audit.smk with rules for summary generation
- Task 5: Created comprehensive test suite with 57 tests (100% pass)
- Task 6: Biopython dependency added and verified

Key implementation details:
- Graceful degradation: Standard codon table (ID=1) works without Biopython
- ID format auto-detection: UniProt, NCBI (GI, RefSeq, GenBank), OrthoFinder, simple
- Consistency checking: compare_id_sets, check_orthofinder_consistency, check_eggnog_consistency
- Audit summary: sequence stats, length histogram, ID format detection
- All functions integrate with existing CompGeneError error handling

### File List

**New Files:**
- `workflow/lib/bio_utils.py` - Biopython 封装模块（ID解析、一致性检查、审计摘要）(~650 lines)
- `workflow/lib/test_bio_utils.py` - 单元测试 (57 tests)
- `workflow/rules/audit.smk` - summary.json 生成规则 (~180 lines)

**Modified Files:**
- `pyproject.toml` - 添加 biopython>=1.81 依赖，更新 pytest testpaths
- `workflow/rules/standardize.smk` - 集成 audit 规则引用
- `docs/planning-artifacts/epics.md` - 添加 Story 2.2b
- `docs/implementation-artifacts/sprint-status.yaml` - 添加 Story 状态

## Change Log

- 2026-01-29: Story file created
- 2026-01-29: Implementation completed - all 6 tasks done, 57 tests passing
- 2026-01-29: Code review completed - 6 issues fixed:
  - [H1] Added workflow/lib to pytest testpaths
  - [H2] Added truncation indicators to ComparisonResult.to_dict()
  - [M1] Updated File List to include standardize.smk
  - [M2] Removed redundant genome/proteins summary rules from audit.smk
  - [M3] Added warning/status fields for empty FASTA files
  - [M4/L1] Converted IDFormat to Enum for type safety
  - Added test_to_dict_truncation test (58 tests total)
