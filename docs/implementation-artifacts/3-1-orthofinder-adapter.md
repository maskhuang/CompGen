# Story 3.1: OrthoFinder Adapter

Status: in-progress

## Story

As a **计算生物学研究员**,
I want **系统能调用 OrthoFinder 对多物种蛋白序列执行 orthogroup 推断**,
so that **我可以识别跨物种的基因家族关系，为下游分析提供直系同源分组**。

## Acceptance Criteria

1. **AC1: OrthoFinder 调用**
   - Given 多个物种的蛋白序列目录（`results/standardized/{species}/proteins.longest.fa.gz`）
   - When 执行 `orthology_infer` rule
   - Then 调用 OrthoFinder 执行 orthogroup 推断
   - And 输出到 `results/orthology/`

2. **AC2: 输出文件生成**
   - Given OrthoFinder 执行完成
   - When 解析输出
   - Then 生成 `Orthogroups/Orthogroups.tsv`（基因到 orthogroup 的映射）
   - And 生成 `Species_Tree/SpeciesTree_rooted.txt`（物种树）
   - And 生成 `Orthogroups/Orthogroups_UnassignedGenes.tsv`（未分配基因）

3. **AC3: 蛋白序列输入准备**
   - Given 多个物种的压缩蛋白序列（`.fa.gz`）
   - When 准备 OrthoFinder 输入
   - Then 解压到临时目录，每物种一个 `.fa` 文件
   - And 文件名为 `{species}.fa`（OrthoFinder 用文件名作为物种标识）

4. **AC4: 版本兼容性**
   - Given OrthoFinder 版本不兼容（不是 2.5.x）
   - When 检测版本
   - Then 返回 E_TOOL_VERSION 错误码
   - And 提供版本升级建议

5. **AC5: 错误处理**
   - Given OrthoFinder 执行失败
   - When 捕获错误
   - Then 返回适当的错误码（E_INPUT_MISSING, E_NONZERO_EXIT, E_OOM 等）
   - And 提供清晰的恢复建议

6. **AC6: 审计记录**
   - Given OrthoFinder 执行完成
   - When Runner 收集元数据
   - Then 生成 `.run.json` 包含 OrthoFinder 版本、物种列表、输入 checksum、运行时间

## Tasks / Subtasks

- [ ] Task 1: 创建 OrthoFinder Adapter (AC: #1, #4, #5, #6)
  - [ ] 1.1 创建 `workflow/adapters/orthofinder.py` 继承 BaseAdapter
  - [ ] 1.2 实现 `spec` 属性定义工具规格（OrthoFinder 2.5.x）
  - [ ] 1.3 实现 `check_version()` 检测 OrthoFinder 版本
  - [ ] 1.4 实现 `validate_inputs()` 验证蛋白序列目录和文件
  - [ ] 1.5 实现 `build_command()` 构建 orthofinder 命令行
  - [ ] 1.6 实现 `expected_outputs()` 定义预期输出文件
  - [ ] 1.7 实现 `parse_outputs()` 解析 orthogroup 统计
  - [ ] 1.8 实现 `timeout_seconds()` 返回超时配置（默认 4 小时）
  - [ ] 1.9 实现 `classify_error()` 分类错误类型

- [ ] Task 2: 创建 conda 环境文件 (AC: #4)
  - [ ] 2.1 创建 `workflow/envs/orthofinder.yaml`
  - [ ] 2.2 指定 OrthoFinder 版本范围（>=2.5.0,<3.0.0）
  - [ ] 2.3 包含所有 OrthoFinder 依赖（diamond, mafft, fasttree, mcl 等）

- [ ] Task 3: 创建 Snakemake 规则 (AC: #1, #2, #3)
  - [ ] 3.1 更新 `workflow/rules/orthology.smk`
  - [ ] 3.2 实现 `rule orthology_prepare_proteins`（解压蛋白序列到临时目录）
  - [ ] 3.3 实现 `rule orthology_infer`（调用 OrthoFinder）
  - [ ] 3.4 创建 `workflow/scripts/run_orthofinder.py` Snakemake 脚本
  - [ ] 3.5 实现输入函数获取所有物种的蛋白序列路径

- [ ] Task 4: 更新配置 schema (AC: #1)
  - [ ] 4.1 更新 `schemas/config.schema.yaml` 添加 orthofinder 配置节
  - [ ] 4.2 定义 orthofinder 参数（threads, search_method, timeout 等）

- [ ] Task 5: 创建单元测试 (AC: #1-6)
  - [ ] 5.1 创建 `workflow/adapters/test_orthofinder.py`
  - [ ] 5.2 测试版本检测逻辑（2.5.x 通过，其他拒绝）
  - [ ] 5.3 测试命令构建（不同参数组合）
  - [ ] 5.4 测试输出解析（orthogroup 统计数据提取）
  - [ ] 5.5 测试错误分类（版本错误、输入缺失、执行失败）

- [ ] Task 6: 更新 Snakefile 集成 (AC: #1)
  - [ ] 6.1 确认 Snakefile 已有 `include: "rules/orthology.smk"`（已存在）
  - [ ] 6.2 更新 orthology.smk 使其功能完整

## Dev Notes

### OrthoFinder 命令行接口 [Source: OrthoFinder GitHub]

**基础命令：**
```bash
# 从蛋白序列目录推断 orthogroups
orthofinder -f /path/to/proteins_dir -t 8 -a 1

# 使用 DIAMOND 搜索（默认）
orthofinder -f proteins_dir -S diamond

# 使用 BLAST 搜索
orthofinder -f proteins_dir -S blast
```

**关键参数：**
| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-f` | 蛋白序列目录（每物种一个 .fa 文件） | 必需 |
| `-t` | 序列搜索线程数 | 1 |
| `-a` | 分析线程数 | 1 |
| `-S` | 搜索方法（diamond, diamond_ultra_sens, blast） | diamond |
| `-M` | 基因树推断方法（dendroblast, msa） | dendroblast |
| `-n` | 结果目录后缀名 | 日期时间 |
| `-o` | 输出目录（覆盖默认路径） | proteins_dir/OrthoFinder/ |
| `-s` | 已有物种树文件（跳过物种树推断） | 无 |
| `-p` | 临时文件目录 | 无 |

**版本检测：**
```bash
orthofinder --version
# 输出: OrthoFinder version 2.5.5
```

**重要行为特性：**
- OrthoFinder 接收一个**目录**作为输入（`-f`），不是单个文件
- 目录中每个 `.fa`/`.fasta`/`.faa` 文件代表一个物种
- **文件名（不含扩展名）作为物种标识符**
- 输出默认写入 `{input_dir}/OrthoFinder/Results_{date}/`
- 使用 `-o` 参数可指定输出目录

### OrthoFinder 输出目录结构 [Source: OrthoFinder GitHub]

```
Results_<date>/
├── Orthogroups/
│   ├── Orthogroups.tsv                    # 主表：每行一个 orthogroup，列为物种
│   ├── Orthogroups.txt                    # 文本格式：每行一个 orthogroup 的所有基因
│   ├── Orthogroups_UnassignedGenes.tsv    # 未分配到 orthogroup 的基因
│   ├── Orthogroups_SingleCopyOrthologues.txt  # 单拷贝 orthogroup 列表
│   └── Orthogroups.GeneCount.tsv          # 每 orthogroup 每物种的基因数
├── Gene_Trees/                            # 每 orthogroup 的基因树（newick）
├── Resolved_Gene_Trees/                   # 解析后的基因树
├── Species_Tree/
│   ├── SpeciesTree_rooted.txt             # 根化物种树（newick）
│   └── SpeciesTree_rooted_node_labels.txt # 带节点标签
├── Phylogenetic_Hierarchical_Orthogroups/
│   └── N0.tsv                             # 层次正交组
├── Comparative_Genomics_Statistics/
│   ├── Statistics_Overall.tsv
│   ├── Statistics_PerSpecies.tsv
│   ├── OrthologuesStats_*.tsv
│   └── Duplications_per_Species_Tree_Node.tsv
└── WorkingDirectory/                      # 中间文件（可删除）
```

### OrthoFinder 输出格式

**Orthogroups.tsv：**
```tsv
Orthogroup	species1	species2	species3
OG0000000	gene1, gene2	gene3	gene4, gene5, gene6
OG0000001	gene7	gene8
OG0000002		gene9	gene10
```

**Orthogroups.GeneCount.tsv：**
```tsv
Orthogroup	species1	species2	species3	Total
OG0000000	2	1	3	6
OG0000001	1	1	0	2
```

### OrthoFinderAdapter 实现框架 [Source: architecture.md#ADR-002]

```python
# workflow/adapters/orthofinder.py

from pathlib import Path
from typing import Optional
import re
import subprocess

from .base import BaseAdapter, ToolSpec, RunResult, AdapterContext, classify_common_errors
from workflow.lib.errors import ErrorCode, CompGeneError


class OrthoFinderAdapter(BaseAdapter):
    """OrthoFinder 工具适配器 - 直系同源推断"""

    SUPPORTED_MAJOR_VERSION = 2
    SUPPORTED_MINOR_MIN = 5
    DEFAULT_TIMEOUT = 14400  # 4 hours

    @property
    def spec(self) -> ToolSpec:
        return ToolSpec(
            name="orthofinder",
            min_version="2.5.0",
            max_version="2.5.99",
            conda_env="orthofinder.yaml",
            description="OrthoFinder orthology inference tool",
        )

    def check_version(self) -> str:
        """检测 OrthoFinder 版本"""
        result = subprocess.run(
            ["orthofinder", "--version"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        # OrthoFinder version 2.5.5
        match = re.search(
            r'OrthoFinder version (\d+\.\d+\.\d+)',
            result.stdout + result.stderr
        )
        if not match:
            raise CompGeneError(
                ErrorCode.E_TOOL_NOT_FOUND,
                "Cannot detect OrthoFinder version"
            )

        version = match.group(1)
        parts = version.split('.')
        major, minor = int(parts[0]), int(parts[1])

        if major != self.SUPPORTED_MAJOR_VERSION or minor < self.SUPPORTED_MINOR_MIN:
            raise CompGeneError(
                ErrorCode.E_TOOL_VERSION,
                f"OrthoFinder {version} not supported. Requires 2.5.x"
            )

        return version

    def validate_inputs(self, ctx: AdapterContext) -> None:
        """验证输入: proteins_dir 目录包含至少 2 个 .fa 文件"""
        proteins_dir = ctx.inputs.get("proteins_dir")
        if not proteins_dir or not proteins_dir.is_dir():
            raise CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                f"Proteins directory not found: {proteins_dir}"
            )

        fa_files = list(proteins_dir.glob("*.fa"))
        if len(fa_files) < 2:
            raise CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                f"Need at least 2 species .fa files in {proteins_dir}, found {len(fa_files)}"
            )

    def build_command(self, ctx: AdapterContext) -> list[str]:
        """构建 OrthoFinder 命令行"""
        config = ctx.config.get("orthofinder", {})

        cmd = [
            "orthofinder",
            "-f", str(ctx.inputs["proteins_dir"]),
            "-t", str(ctx.threads),
            "-a", str(min(ctx.threads, 4)),  # 分析线程通常不需要太多
        ]

        # 输出目录
        if "output_dir" in ctx.outputs:
            cmd.extend(["-o", str(ctx.outputs["output_dir"])])

        # 搜索方法
        search_method = config.get("search_method", "diamond")
        cmd.extend(["-S", search_method])

        # 结果目录后缀名
        cmd.extend(["-n", "compgene"])

        return cmd

    def expected_outputs(self, ctx: AdapterContext) -> list[Path]:
        """定义预期输出文件"""
        results_dir = ctx.outputs["results_dir"]
        return [
            results_dir / "Orthogroups" / "Orthogroups.tsv",
            results_dir / "Orthogroups" / "Orthogroups.GeneCount.tsv",
            results_dir / "Species_Tree" / "SpeciesTree_rooted.txt",
        ]

    def parse_outputs(self, ctx: AdapterContext) -> RunResult:
        """解析 OrthoFinder 输出统计"""
        results_dir = ctx.outputs["results_dir"]
        gene_count_path = results_dir / "Orthogroups" / "Orthogroups.GeneCount.tsv"
        unassigned_path = results_dir / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv"

        # 统计 orthogroup 数量
        total_ogs = 0
        if gene_count_path.exists():
            with open(gene_count_path) as f:
                total_ogs = sum(1 for line in f) - 1  # 减去表头

        # 统计未分配基因
        unassigned_count = 0
        if unassigned_path.exists():
            with open(unassigned_path) as f:
                unassigned_count = sum(1 for line in f) - 1

        summary = {
            "total_orthogroups": total_ogs,
            "unassigned_genes": unassigned_count,
        }

        return RunResult(
            outputs={
                "orthogroups_tsv": results_dir / "Orthogroups" / "Orthogroups.tsv",
                "gene_count_tsv": gene_count_path,
                "species_tree": results_dir / "Species_Tree" / "SpeciesTree_rooted.txt",
            },
            summary=summary,
        )

    def timeout_seconds(self, ctx: AdapterContext) -> int:
        """返回超时时间 (OrthoFinder 对大物种集合可能很慢)"""
        return ctx.config.get("orthofinder", {}).get("timeout", self.DEFAULT_TIMEOUT)

    def classify_error(
        self, ctx: AdapterContext, returncode: int, stderr: str
    ) -> tuple[ErrorCode, bool]:
        """分类错误类型"""
        stderr_lower = stderr.lower()

        # DIAMOND/BLAST 数据库错误
        if "diamond" in stderr_lower and "error" in stderr_lower:
            return (ErrorCode.E_NONZERO_EXIT, False)

        # 输入文件问题
        if "no fasta files" in stderr_lower or "no sequences" in stderr_lower:
            return (ErrorCode.E_INPUT_MISSING, False)

        # MCL 错误
        if "mcl" in stderr_lower and "error" in stderr_lower:
            return (ErrorCode.E_NONZERO_EXIT, False)

        # 内存错误
        if "memory" in stderr_lower or "oom" in stderr_lower:
            return (ErrorCode.E_OOM, False)

        # SIGKILL 通常是超时
        if returncode == -9:
            return (ErrorCode.E_TIMEOUT, True)

        # 回退到通用分类
        return classify_common_errors(returncode, stderr)
```

### Snakemake 规则 [Source: architecture.md#Naming-Patterns]

```python
# workflow/rules/orthology.smk

rule orthology_prepare_proteins:
    """解压蛋白序列到临时目录，准备 OrthoFinder 输入"""
    input:
        proteins=lambda wildcards: expand(
            "results/standardized/{species}/proteins.longest.fa.gz",
            species=[sp["name"] for sp in config["species"]]
        )
    output:
        proteins_dir=directory("results/orthology/input_proteins"),
        done="results/orthology/input_proteins/.prepared"
    run:
        import gzip, shutil
        from pathlib import Path

        out_dir = Path(output.proteins_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        for proteins_gz in input.proteins:
            species = Path(proteins_gz).parent.name
            out_fa = out_dir / f"{species}.fa"
            with gzip.open(proteins_gz, 'rb') as f_in:
                with open(out_fa, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

        Path(output.done).touch()


rule orthology_infer:
    """调用 OrthoFinder 执行 orthogroup 推断"""
    input:
        proteins_dir="results/orthology/input_proteins",
        done="results/orthology/input_proteins/.prepared"
    output:
        orthogroups="results/orthology/Orthogroups/Orthogroups.tsv",
        gene_count="results/orthology/Orthogroups/Orthogroups.GeneCount.tsv",
        species_tree="results/orthology/Species_Tree/SpeciesTree_rooted.txt",
        run_json="results/meta/orthology_infer/run.run.json",
    params:
        search_method=lambda wildcards: config.get("orthofinder", {}).get("search_method", "diamond"),
        results_dir="results/orthology",
    threads: get_threads("orthology_infer")
    log:
        "logs/orthology_infer/run.log"
    conda:
        "../envs/orthofinder.yaml"
    script:
        "../scripts/run_orthofinder.py"
```

### conda 环境文件 [Source: architecture.md#conda-环境隔离]

```yaml
# workflow/envs/orthofinder.yaml
name: compgene_orthofinder
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - orthofinder>=2.5.0,<3.0.0
  - python>=3.9
  - diamond>=2.0
  - mafft>=7.0
  - fasttree>=2.1
  - mcl>=14.137
```

### 配置 Schema 更新 [Source: schemas/config.schema.yaml]

```yaml
# 新增 orthofinder 配置节
orthofinder:
  type: object
  properties:
    search_method:
      type: string
      enum: ["diamond", "diamond_ultra_sens", "blast"]
      default: "diamond"
      description: "序列搜索方法"
    tree_method:
      type: string
      enum: ["dendroblast", "msa"]
      default: "dendroblast"
      description: "基因树推断方法"
    timeout:
      type: integer
      minimum: 600
      maximum: 259200
      default: 14400
      description: "超时时间（秒），默认 4 小时"
  additionalProperties: false
```

### 错误处理 [Source: architecture.md#ADR-003]

| 场景 | 错误码 | 可重试 | 恢复建议 |
|------|--------|--------|----------|
| 蛋白序列目录不存在 | E_INPUT_MISSING | ❌ | 先运行 standardize 规则生成蛋白序列 |
| 物种数量不足（<2） | E_INPUT_MISSING | ❌ | 至少需要 2 个物种的蛋白序列 |
| OrthoFinder 版本不兼容 | E_TOOL_VERSION | ❌ | 安装 OrthoFinder 2.5.x |
| DIAMOND 数据库错误 | E_NONZERO_EXIT | ❌ | 检查 DIAMOND 安装 |
| 执行超时 | E_TIMEOUT | ✅ | 增加 timeout 或减少物种数 |
| 内存不足 | E_OOM | ❌ | 减少线程或增加内存 |

### run_orthofinder.py 脚本框架 [Source: 5-1-liftoff-adapter.md 模式]

```python
# workflow/scripts/run_orthofinder.py
"""
Run OrthoFinder through the Adapter framework.

Bridges Snakemake with OrthoFinderAdapter:
- Version checking (AC4)
- Error classification (AC5)
- Audit record generation (AC6)
- Input preparation (AC3): decompress .fa.gz to temp dir
- Output symlinking: link OrthoFinder results to expected paths
"""

import subprocess
import time
from pathlib import Path

from workflow.adapters.orthofinder import OrthoFinderAdapter
from workflow.adapters.base import AdapterContext
from workflow.lib.errors import CompGeneError, ErrorCode
from workflow.lib.audit import create_and_write_audit


def main():
    adapter = OrthoFinderAdapter()

    # 1. Version check
    version = adapter.check_version()

    # 2. Build context from snakemake variables
    ctx = AdapterContext(
        inputs={"proteins_dir": Path(snakemake.input.proteins_dir)},
        outputs={
            "results_dir": Path(snakemake.params.results_dir),
            "output_dir": Path(snakemake.params.results_dir),
        },
        config=snakemake.config,
        wildcards=dict(snakemake.wildcards),
        threads=snakemake.threads,
    )

    # 3. Validate inputs
    adapter.validate_inputs(ctx)

    # 4. Build and execute command
    cmd = adapter.build_command(ctx)
    start_time = time.time()

    try:
        timeout = adapter.timeout_seconds(ctx)
        # Redirect stdout/stderr to log file to avoid pipe buffer issues
        log_path = Path(snakemake.log[0])
        log_path.parent.mkdir(parents=True, exist_ok=True)

        with open(log_path, 'w') as log_fh:
            proc = subprocess.run(
                cmd,
                stdout=log_fh,
                stderr=log_fh,
                timeout=timeout,
            )

        runtime = time.time() - start_time

        if proc.returncode != 0:
            log_content = log_path.read_text() if log_path.exists() else ""
            error_code, retryable = adapter.classify_error(ctx, proc.returncode, log_content)
            raise CompGeneError(error_code, f"OrthoFinder failed (exit {proc.returncode})")

    except subprocess.TimeoutExpired:
        runtime = time.time() - start_time
        raise CompGeneError(ErrorCode.E_TIMEOUT, f"OrthoFinder timed out after {timeout}s")

    # 5. Handle OrthoFinder output directory naming
    # OrthoFinder writes to {input_dir}/OrthoFinder/Results_compgene/
    # or to -o directory if specified. Symlink/move results to expected paths.

    # 6. Parse outputs and generate audit
    result = adapter.parse_outputs(ctx)

    create_and_write_audit(
        rule="orthology_infer",
        wildcards=dict(snakemake.wildcards),
        cmd=cmd,
        tool_version=version,
        input_paths=[Path(snakemake.input.proteins_dir)],
        output_paths=[Path(p) for p in [
            snakemake.output.orthogroups,
            snakemake.output.gene_count,
            snakemake.output.species_tree,
        ]],
        runtime_seconds=runtime,
        exit_code=0,
        audit_path=Path(snakemake.output.run_json),
    )


if __name__ == "__main__":
    main()
```

### 从前序 Story 学到的经验

**从 Story 5.1（Liftoff Adapter）：**
1. **Adapter 需要配套 Snakemake 脚本**：创建 `run_orthofinder.py` 桥接 Snakemake 和 Adapter
2. **输出文件需显式声明**：在 Snakemake rule 中声明所有输出文件
3. **审计记录必须生成**：使用 `create_and_write_audit()` 生成 `.run.json`
4. **处理压缩文件**：输入蛋白序列是 gzip 压缩的，需要解压到临时目录
5. **版本检测正则**：根据实际输出格式调整正则表达式
6. **避免 capture_output=True**：大量输出时会导致 pipe buffer deadlock，应重定向到日志文件
7. **使用 .resolve() 处理路径**：subprocess 的 cwd 可能导致相对路径失败
8. **SGE 集群支持**：大基因组需要提交到集群节点（参考 profile/sge/config.yaml）

**从 Story 2.4（BUSCO QC）：**
1. **Adapter 模式一致性**：严格遵循 BaseAdapter 8 个方法接口
2. **原子写入**：所有输出文件先写 .tmp 再 rename
3. **配置 schema 更新**：添加新工具时同步更新 schema

### 关键实现注意事项

1. **OrthoFinder 输出目录命名**：OrthoFinder 会自动在输入目录下创建 `OrthoFinder/Results_<suffix>/` 子目录。使用 `-n compgene` 参数固定后缀名为 `compgene`，使用 `-o` 指定输出位置，避免路径不确定性。

2. **蛋白序列准备**：现有蛋白序列在 `results/standardized/{species}/proteins.longest.fa.gz`（由 Story 2.2 生成）。OrthoFinder 需要**未压缩的 .fa 文件**放在同一目录中。需要一个准备步骤解压到 `results/orthology/input_proteins/` 目录。

3. **物种标识符一致性**：OrthoFinder 使用**文件名**（不含扩展名）作为物种名。解压时文件命名为 `{species}.fa`，确保与 config 中的 species name 一致。

4. **OrthoFinder 是全物种集合计算**：不同于 BUSCO/Liftoff 可以按物种并行，OrthoFinder 必须一次处理所有物种。新增物种 = 全量重算。

5. **默认超时 4 小时**：3 物种约 1-2 小时，7+ 物种可能需要数小时。

### Project Structure Notes

**新增文件：**
```
workflow/
├── adapters/
│   ├── orthofinder.py           # OrthoFinder Adapter（新建）
│   └── test_orthofinder.py      # 单元测试（新建）
├── rules/
│   └── orthology.smk            # 直系同源规则（更新，替换骨架）
├── scripts/
│   └── run_orthofinder.py       # Snakemake 脚本（新建）
└── envs/
    └── orthofinder.yaml         # conda 环境（新建）
```

**修改文件：**
```
schemas/config.schema.yaml       # 添加 orthofinder 配置节
workflow/rules/orthology.smk     # 替换骨架为完整规则
```

### References

- [Source: docs/planning-artifacts/prd.md#FR6] - 系统可调用 OrthoFinder 执行 orthogroup 推断
- [Source: docs/planning-artifacts/prd.md#FR7] - 系统可生成 orthogroups 表
- [Source: docs/planning-artifacts/prd.md#FR8] - 系统可输出物种树和基因树
- [Source: docs/planning-artifacts/prd.md#NFR10] - OrthoFinder 调用应支持版本 2.5.x
- [Source: docs/planning-artifacts/prd.md#NFR15] - 所有外部工具调用应有超时保护
- [Source: docs/planning-artifacts/architecture.md#ADR-002] - 工具适配层 Adapter 模式
- [Source: docs/planning-artifacts/architecture.md#Integration-Points] - OrthoFinder 集成：subprocess, proteins/ 目录输入
- [Source: docs/planning-artifacts/epics.md#Story-3.1] - 蛋白序列提取（输入准备）
- [Source: docs/planning-artifacts/epics.md#Story-3.2] - OrthoFinder 运行（核心推断）

## Dev Agent Record

### Agent Model Used

Claude Opus 4.6 (claude-opus-4-6)

### Debug Log References

### Completion Notes List

### File List
