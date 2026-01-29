# Story 2.1: GFF/FASTA 解析库

Status: done

## Story

As a **计算生物学研究员**,
I want **系统能解析 GFF3/GTF 和 FASTA 格式文件**,
so that **我可以使用不同来源的注释数据**。

## Acceptance Criteria

1. **AC1: GFF3/GTF 解析**
   - Given 用户提供 GFF3 或 GTF 格式注释文件
   - When 系统解析文件
   - Then 正确提取 gene、transcript、CDS 等特征
   - And 建立 gene → transcript → protein 的层级关系

2. **AC2: 格式错误处理**
   - Given 文件格式错误
   - When 解析失败
   - Then 返回 E_INPUT_FORMAT 错误码和具体错误位置
   - And 提供清晰的错误描述

3. **AC3: FASTA 解析**
   - Given 用户提供 FASTA 格式文件（基因组或蛋白序列）
   - When 系统解析文件
   - Then 正确提取序列 ID 和序列内容
   - And 支持 gzip 压缩文件

4. **AC4: 大文件支持**
   - Given 大文件（>100MB）
   - When 解析文件
   - Then 使用流式处理/迭代器模式
   - And 内存占用可控（不一次性加载全部）

## Tasks / Subtasks

- [x] Task 1: 创建 lib/gff.py 模块 (AC: #1, #2)
  - [x] 1.1 定义 GFFFeature dataclass（seqid, source, type, start, end, score, strand, phase, attributes）
  - [x] 1.2 定义 GeneModel dataclass（gene_id, transcripts, representative_transcript）
  - [x] 1.3 实现 parse_gff3(path) 流式解析 GFF3 文件
  - [x] 1.4 实现 parse_gtf(path) 流式解析 GTF 文件
  - [x] 1.5 实现 build_gene_hierarchy() 构建 gene → transcript → CDS 层级
  - [x] 1.6 实现 get_representative_transcript() 选择最长转录本
  - [x] 1.7 实现 validate_gff_line() 格式验证，返回错误位置和描述
  - [x] 1.8 支持 gzip 压缩文件 (.gff3.gz, .gtf.gz)

- [x] Task 2: 创建 lib/fasta.py 模块 (AC: #3, #4)
  - [x] 2.1 定义 FastaRecord dataclass（id, description, sequence）
  - [x] 2.2 实现 parse_fasta(path) 流式解析 FASTA 文件
  - [x] 2.3 实现 write_fasta(records, path) 原子写入 FASTA
  - [x] 2.4 支持 gzip 压缩文件 (.fa.gz, .fasta.gz)
  - [x] 2.5 实现 get_sequence_length() 和 get_sequence_stats() 统计函数

- [x] Task 3: 实现格式检测与错误处理 (AC: #2)
  - [x] 3.1 实现 detect_format(path) 自动检测 GFF3/GTF/FASTA 格式
  - [x] 3.2 实现 validate_file(path, expected_format) 验证文件格式
  - [x] 3.3 使用 E_INPUT_FORMAT 错误码报告格式错误
  - [x] 3.4 错误信息包含行号和具体问题描述

- [x] Task 4: 创建单元测试 (AC: #1-4)
  - [x] 4.1 创建 lib/test_gff.py 测试 GFF3/GTF 解析
  - [x] 4.2 创建 lib/test_fasta.py 测试 FASTA 解析
  - [x] 4.3 创建测试数据文件 tests/fixtures/mini_annotation.gff3
  - [x] 4.4 创建测试数据文件 tests/fixtures/mini_genome.fa
  - [x] 4.5 测试格式错误检测场景
  - [x] 4.6 测试大文件流式处理（模拟）

- [x] Task 5: 更新 lib/__init__.py 导出 (AC: #1-4)
  - [x] 5.1 导出 gff 模块的公共接口
  - [x] 5.2 导出 fasta 模块的公共接口
  - [x] 5.3 更新模块文档字符串

## Dev Notes

### GFF3 格式规范 [Source: docs/planning-artifacts/prd.md#FR2-FR3]

GFF3 是标准的基因组注释格式，每行包含 9 个 Tab 分隔字段：
```
seqid  source  type  start  end  score  strand  phase  attributes
```

**attributes 字段**：分号分隔的 key=value 对
- `ID=gene001` - 唯一标识符
- `Parent=gene001` - 父级特征 ID
- `Name=BRCA1` - 人类可读名称

**层级关系**：
```
gene (ID=gene001)
  └── mRNA (ID=transcript001, Parent=gene001)
        ├── exon (Parent=transcript001)
        ├── CDS (Parent=transcript001)
        └── CDS (Parent=transcript001)
```

### GTF 格式差异

GTF 使用不同的 attributes 格式：
- 空格分隔的 key "value"; 对
- `gene_id "gene001"; transcript_id "transcript001";`

解析时需要检测格式并使用对应的解析逻辑。

### 代表转录本选择规则 [Source: docs/planning-artifacts/epics.md#Story-2.2]

当基因有多个转录本时，选择最长转录本作为代表：
1. 计算每个转录本的 CDS 总长度
2. 选择 CDS 最长的转录本
3. 如果无 CDS，使用外显子总长度
4. 相同长度时选择 ID 字典序最小的

### 项目结构规范 [Source: docs/planning-artifacts/architecture.md#Structure-Patterns]

新建文件位置：
```
compgene/
├── workflow/
│   └── lib/
│       ├── gff.py          # 新建
│       ├── test_gff.py     # 新建 - 单元测试共置
│       ├── fasta.py        # 新建
│       ├── test_fasta.py   # 新建 - 单元测试共置
│       ├── errors.py       # 已存在 - 使用 E_INPUT_FORMAT
│       └── io.py           # 已存在 - 使用 atomic_write
├── tests/
│   └── fixtures/
│       ├── mini_annotation.gff3  # 新建
│       └── mini_genome.fa        # 新建
```

### 命名规范 [Source: docs/planning-artifacts/architecture.md#Naming-Patterns]

- 函数：snake_case (`parse_gff3`, `build_gene_hierarchy`)
- 类：PascalCase (`GFFFeature`, `FastaRecord`, `GeneModel`)
- 常量：UPPER_SNAKE_CASE (`DEFAULT_ENCODING`)

### 错误处理规范 [Source: docs/planning-artifacts/architecture.md#ADR-003]

使用已定义的 ErrorCode：
```python
from workflow.lib.errors import ErrorCode, CompGeneError

# 格式错误
raise CompGeneError(
    ErrorCode.E_INPUT_FORMAT,
    f"Invalid GFF3 format at line {line_num}",
    details=f"Expected 9 fields, got {len(fields)}"
)
```

### 原子写入规范 [Source: docs/planning-artifacts/architecture.md#Process-Patterns]

写入 FASTA 文件时使用原子模式：
```python
from workflow.lib.io import atomic_write

def write_fasta(records, output_path: Path):
    # atomic_write 自动处理 .tmp 文件和重命名
    with atomic_write(output_path, mode='w') as f:
        for record in records:
            f.write(f">{record.id} {record.description}\n")
            f.write(f"{record.sequence}\n")
```

### gzip 支持

使用 Python 标准库处理压缩文件：
```python
import gzip
from pathlib import Path

def open_file(path: Path, mode='rt'):
    """Open file, automatically detecting gzip compression."""
    if path.suffix == '.gz':
        return gzip.open(path, mode, encoding='utf-8')
    return open(path, mode, encoding='utf-8')
```

### 流式处理模式 [Source: docs/planning-artifacts/prd.md#NFR2]

大文件使用生成器模式：
```python
def parse_gff3(path: Path) -> Iterator[GFFFeature]:
    """Yield GFF features one at a time (memory efficient)."""
    with open_file(path) as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            yield parse_gff_line(line, line_num)
```

### 依赖模块

**已完成的依赖（Epic 1）：**
- `workflow/lib/errors.py` - ErrorCode.E_INPUT_FORMAT
- `workflow/lib/io.py` - atomic_write() 函数

**无外部依赖：** 仅使用 Python 标准库

### References

- [Source: docs/planning-artifacts/prd.md#FR2] - GFF3/GTF 解析需求
- [Source: docs/planning-artifacts/prd.md#FR3] - 特征提取需求
- [Source: docs/planning-artifacts/prd.md#NFR2] - 大文件内存控制
- [Source: docs/planning-artifacts/architecture.md#ADR-003] - 错误处理规范
- [Source: docs/planning-artifacts/architecture.md#Structure-Patterns] - 项目结构
- [Source: docs/planning-artifacts/architecture.md#Naming-Patterns] - 命名规范
- [Source: docs/planning-artifacts/epics.md#Story-2.1] - Story 需求定义
- [Source: docs/planning-artifacts/epics.md#Story-2.2] - 代表转录本选择规则

## Dev Agent Record

### Agent Model Used

Claude Opus 4.5 (claude-opus-4-5-20251101)

### Debug Log References

N/A - No debug issues encountered

### Completion Notes List

- Created `workflow/lib/gff.py` with complete GFF3/GTF parsing implementation:
  - GFFFeature, TranscriptModel, GeneModel dataclasses
  - Stream-based parse_gff3() and parse_gtf() functions using generators
  - build_gene_hierarchy() constructs gene → transcript → CDS hierarchy
  - get_representative_transcript() selects longest CDS (with exon fallback and tiebreaker)
  - validate_gff_line() returns detailed error with line number
  - Auto-detection of GFF3 vs GTF format
  - gzip compression support
  - [Code Review Fix] Added _get_parent_transcript_ids() for correct GTF/GFF3 handling
  - [Code Review Fix] Multi-parent feature support for alternative splicing

- Created `workflow/lib/fasta.py` with complete FASTA parsing implementation:
  - FastaRecord and SequenceStats dataclasses
  - Stream-based parse_fasta() using generators (memory efficient)
  - write_fasta() with atomic write pattern for checkpoint safety
  - write_fasta_gzip() for compressed output (now atomic)
  - get_sequence_stats() computes count, lengths, N50, GC content
  - validate_fasta() for format validation (now catches empty files)

- Created format detection utilities:
  - detect_format() auto-detects GFF3/GTF/FASTA from extension or content
  - validate_file() validates against expected format
  - All errors use E_INPUT_FORMAT with line numbers and details

- Created comprehensive test suites:
  - 43 tests in test_gff.py covering parsing, validation, hierarchy, format detection, fixtures
  - 34 tests in test_fasta.py covering parsing, writing, statistics, streaming, fixtures
  - All 77 tests pass

- Created test fixtures:
  - tests/fixtures/mini_annotation.gff3 (2 genes, 3 transcripts, exons, CDS)
  - tests/fixtures/mini_annotation.gtf (GTF format sample)
  - tests/fixtures/mini_genome.fa (nucleotide and protein sequences)

- Updated lib/__init__.py with documentation for new modules

### File List

- workflow/lib/gff.py (created - GFF3/GTF parser module)
- workflow/lib/fasta.py (created - FASTA parser module)
- workflow/lib/test_gff.py (created - 36 unit tests)
- workflow/lib/test_fasta.py (created - 30 unit tests)
- workflow/lib/__init__.py (modified - added module documentation)
- tests/fixtures/mini_annotation.gff3 (created - test fixture)
- tests/fixtures/mini_annotation.gtf (created - test fixture)
- tests/fixtures/mini_genome.fa (created - test fixture)

## Change Log

- 2026-01-28: Story file created by create-story workflow
- 2026-01-28: Implementation completed - all 5 tasks done, 66 tests pass, status changed to review
- 2026-01-28: Code review completed - 8 issues found and fixed:
  - Issue #1 (CRITICAL): Fixed GTF hierarchy construction - CDS/exon now correctly assigned via transcript_id
  - Issue #2 (HIGH): Made write_fasta_gzip() atomic with temp file + rename pattern
  - Issue #3 (MEDIUM): Added 11 new tests for actual fixture files
  - Issue #4 (MEDIUM): Fixed by new _get_parent_transcript_ids() helper function
  - Issue #5 (MEDIUM): Added multi-parent feature support (comma-separated Parent values)
  - Issue #6 (LOW): Fixed type annotations to use IO[str] instead of Union[..., "open"]
  - Issue #7 (LOW): validate_fasta() now fails for empty/whitespace-only files
  - Issue #8 (LOW): Added complete docstring with examples to parse_annotation()
  - Total tests: 77 (was 66)
