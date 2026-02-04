# Story 2.3: NCBI 注释下载

Status: done

## Story

As a **计算生物学研究员**,
I want **通过物种标识符自动从 NCBI 下载注释**,
so that **我不需要手动查找和下载数据**。

## Acceptance Criteria

1. **AC1: NCBI 标识符识别**
   - Given 配置文件中 annotation 为 NCBI 标识符（如 `NCBI:GCF_000165445`）
   - When 执行 standardize rule
   - Then 系统识别 NCBI 前缀并触发自动下载流程
   - And 支持 GCF_（RefSeq）和 GCA_（GenBank）两种 accession 格式

2. **AC2: 自动下载基因组和注释文件**
   - Given 有效的 NCBI accession
   - When 执行下载
   - Then 从 NCBI FTP/API 下载以下文件：
     - 基因组序列（*_genomic.fna.gz）
     - 注释文件（*_genomic.gff.gz）
   - And 文件保存到 `results/downloads/ncbi/{accession}/`
   - And 生成 download_manifest.json 记录下载详情

3. **AC3: 本地缓存机制**
   - Given 已下载过相同标识符的数据
   - When 再次请求下载
   - Then 使用本地缓存，不重复下载
   - And 通过文件存在性 + checksum 验证缓存有效性
   - And 支持 `--force-download` 标志强制重新下载

4. **AC4: NCBI 速率限制遵循**
   - Given NCBI API 调用频繁
   - When 超过速率限制（3 req/s）
   - Then 自动等待并重试，遵循 NCBI 官方限制
   - And 使用指数退避策略（1s, 2s, 4s...最大 60s）
   - And 支持 NCBI_API_KEY 环境变量提升速率限制

5. **AC5: 错误处理与恢复**
   - Given 网络错误或 NCBI 服务不可用
   - When 下载失败
   - Then 返回 E_NET_RATE_LIMIT 或 E_NET_ERROR 错误码
   - And 自动重试最多 3 次
   - And 提供清晰的错误信息和恢复建议

6. **AC6: 与标准化流程集成**
   - Given NCBI 下载完成
   - When 触发后续 standardize rule
   - Then 下载的文件作为 standardize 的输入
   - And 遵循现有标准化输出格式（genome.fa.gz, annotation.gff3.gz）

## Tasks / Subtasks

- [x] Task 1: 创建 lib/ncbi.py NCBI API 客户端模块 (AC: #1, #4, #5) ✅
  - [x] 1.1 实现 `NCBIClient` 类，封装 NCBI API 调用
  - [x] 1.2 实现 `parse_ncbi_accession(identifier)` 解析 NCBI:GCF_xxx 格式
  - [x] 1.3 实现速率限制器 `RateLimiter` 类，支持 3 req/s 限制
  - [x] 1.4 实现指数退避重试逻辑 `with_retry(func, max_retries=3)`
  - [x] 1.5 实现 NCBI_API_KEY 环境变量支持

- [x] Task 2: 实现 NCBI 文件下载功能 (AC: #2, #3) ✅
  - [x] 2.1 实现 `resolve_ncbi_urls(accession)` 获取 FTP URL 列表 → `resolve_ftp_directory()` + `get_assembly_files()`
  - [x] 2.2 实现 `download_file(url, dest_path)` 单文件下载
  - [x] 2.3 实现 `download_ncbi_dataset(accession, output_dir)` 完整数据集下载
  - [x] 2.4 实现下载进度显示和日志记录
  - [x] 2.5 实现 download_manifest.json 生成

- [x] Task 3: 实现缓存机制 (AC: #3) ✅
  - [x] 3.1 实现 `is_cached(accession, output_dir)` 缓存检查
  - [x] 3.2 实现 `verify_cache(accession, output_dir)` checksum 验证
  - [x] 3.3 实现缓存元数据存储（下载时间、来源 URL、checksum）
  - [x] 3.4 实现 `--force-download` 配置项支持 → `force` 参数

- [x] Task 4: 创建 Snakemake 规则 (AC: #2, #6) ✅
  - [x] 4.1 创建 rules/download.smk
  - [x] 4.2 实现 rule ncbi_download（单 accession）
  - [x] 4.3 实现 rule ncbi_download_all（批量下载）
  - [x] 4.4 更新 standardize.smk 集成 NCBI 下载输出 → 在 download.smk 中实现 standardize_ncbi_genome/annotation

- [x] Task 5: 更新配置 schema (AC: #1) ✅
  - [x] 5.1 更新 config.schema.yaml 支持 annotation: "NCBI:GCF_xxx" 格式
  - [x] 5.2 添加 ncbi 配置节（api_key、rate_limit、cache_dir）
  - [x] 5.3 添加 force_download 配置项

- [x] Task 6: 创建单元测试 (AC: #1-6) ✅
  - [x] 6.1 创建 lib/test_ncbi.py
  - [x] 6.2 测试 accession 解析（GCF/GCA 格式）
  - [x] 6.3 测试速率限制器行为
  - [x] 6.4 测试缓存机制（mock 文件系统）
  - [x] 6.5 测试错误处理和重试逻辑
  - [x] 6.6 测试与 standardize 的集成

## Dev Notes

### NCBI API 与 FTP 访问策略

**推荐方案：NCBI Datasets API + FTP 混合**

1. **NCBI Datasets API**（首选，需要 API key）
   - 端点: `https://api.ncbi.nlm.nih.gov/datasets/v2alpha/`
   - 优势: 结构化响应、元数据丰富
   - 限制: 需要 API key 以提升速率限制

2. **NCBI FTP**（备选，无需 key）
   - 基础 URL: `https://ftp.ncbi.nlm.nih.gov/genomes/all/`
   - 路径模式: `GCF/000/165/445/GCF_000165445.1_Mmur_1.0/`
   - 文件: `*_genomic.fna.gz`, `*_genomic.gff.gz`

### NCBI Accession 格式解析 [Source: NCBI RefSeq]

```python
# workflow/lib/ncbi.py

import re
from dataclasses import dataclass
from typing import Optional
from enum import Enum

class AccessionType(Enum):
    REFSEQ = "refseq"   # GCF_ prefix
    GENBANK = "genbank" # GCA_ prefix

@dataclass
class NCBIAccession:
    """解析后的 NCBI accession"""
    raw: str              # 原始输入
    prefix: str           # GCF or GCA
    numeric: str          # 000165445
    version: Optional[str] # 1 (可选)
    accession_type: AccessionType

    @property
    def full_accession(self) -> str:
        """返回完整 accession，如 GCF_000165445.1"""
        if self.version:
            return f"{self.prefix}_{self.numeric}.{self.version}"
        return f"{self.prefix}_{self.numeric}"

    @property
    def ftp_path_parts(self) -> tuple[str, str, str]:
        """返回 FTP 路径分段，如 ('000', '165', '445')"""
        return (self.numeric[0:3], self.numeric[3:6], self.numeric[6:9])

# 解析函数
NCBI_PATTERN = re.compile(
    r'^(?:NCBI:)?(GC[FA])_(\d{9})(?:\.(\d+))?$',
    re.IGNORECASE
)

def parse_ncbi_accession(identifier: str) -> NCBIAccession:
    """解析 NCBI accession 标识符"""
    match = NCBI_PATTERN.match(identifier.strip())
    if not match:
        raise ValueError(f"Invalid NCBI accession format: {identifier}")

    prefix = match.group(1).upper()
    numeric = match.group(2)
    version = match.group(3)

    acc_type = AccessionType.REFSEQ if prefix == "GCF" else AccessionType.GENBANK

    return NCBIAccession(
        raw=identifier,
        prefix=prefix,
        numeric=numeric,
        version=version,
        accession_type=acc_type
    )
```

### 速率限制器实现 [Source: architecture.md#NFR9]

```python
import time
import threading
from typing import TypeVar, Callable

T = TypeVar('T')

class RateLimiter:
    """令牌桶速率限制器"""

    def __init__(self, requests_per_second: float = 3.0):
        self.rate = requests_per_second
        self.tokens = requests_per_second
        self.last_update = time.monotonic()
        self.lock = threading.Lock()

    def acquire(self) -> float:
        """获取一个令牌，返回等待时间"""
        with self.lock:
            now = time.monotonic()
            elapsed = now - self.last_update
            self.tokens = min(self.rate, self.tokens + elapsed * self.rate)
            self.last_update = now

            if self.tokens >= 1:
                self.tokens -= 1
                return 0.0

            wait_time = (1 - self.tokens) / self.rate
            time.sleep(wait_time)
            self.tokens = 0
            self.last_update = time.monotonic()
            return wait_time

def with_retry(
    func: Callable[[], T],
    max_retries: int = 3,
    base_delay: float = 1.0,
    max_delay: float = 60.0,
    retryable_exceptions: tuple = (IOError, TimeoutError)
) -> T:
    """带指数退避的重试装饰器"""
    last_exception = None

    for attempt in range(max_retries + 1):
        try:
            return func()
        except retryable_exceptions as e:
            last_exception = e
            if attempt < max_retries:
                delay = min(base_delay * (2 ** attempt), max_delay)
                logger.warning(f"Retry {attempt + 1}/{max_retries} after {delay}s: {e}")
                time.sleep(delay)

    raise last_exception
```

### FTP URL 解析与下载 [Source: NCBI FTP 结构]

```python
from pathlib import Path
from typing import Optional
import urllib.request
import hashlib

NCBI_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"

def resolve_ncbi_urls(accession: NCBIAccession) -> dict[str, str]:
    """解析 NCBI accession 对应的 FTP 文件 URL"""
    p1, p2, p3 = accession.ftp_path_parts
    prefix = accession.prefix

    # FTP 目录路径
    # 例: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/165/445/GCF_000165445.1_Mmur_1.0/
    base_path = f"{NCBI_FTP_BASE}/{prefix}/{p1}/{p2}/{p3}"

    # 需要先获取目录列表以确定完整路径（包含 assembly 名称）
    # 这里返回模式，实际实现需要 HTTP 请求获取目录内容
    return {
        "base_path": base_path,
        "patterns": {
            "genome": "*_genomic.fna.gz",
            "annotation": "*_genomic.gff.gz",
            "report": "*_assembly_report.txt"
        }
    }

def download_file(
    url: str,
    dest_path: Path,
    rate_limiter: RateLimiter,
    chunk_size: int = 8192
) -> dict:
    """下载单个文件，返回元数据"""
    rate_limiter.acquire()

    dest_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = dest_path.with_suffix('.tmp')

    sha256 = hashlib.sha256()
    total_bytes = 0

    try:
        with urllib.request.urlopen(url) as response:
            with open(temp_path, 'wb') as f:
                while chunk := response.read(chunk_size):
                    f.write(chunk)
                    sha256.update(chunk)
                    total_bytes += len(chunk)

        # 原子重命名
        temp_path.rename(dest_path)

        return {
            "url": url,
            "path": str(dest_path),
            "size_bytes": total_bytes,
            "sha256": sha256.hexdigest()
        }
    except Exception:
        if temp_path.exists():
            temp_path.unlink()
        raise
```

### 缓存机制 [Source: architecture.md#ADR-001]

```python
import json
from datetime import datetime

CACHE_MANIFEST = "download_manifest.json"

def is_cached(accession: str, cache_dir: Path) -> bool:
    """检查数据是否已缓存"""
    manifest_path = cache_dir / accession / CACHE_MANIFEST
    if not manifest_path.exists():
        return False

    try:
        with open(manifest_path) as f:
            manifest = json.load(f)

        # 检查所有文件是否存在
        for file_info in manifest.get("files", []):
            file_path = Path(file_info["path"])
            if not file_path.exists():
                return False

        return True
    except (json.JSONDecodeError, KeyError):
        return False

def verify_cache(accession: str, cache_dir: Path) -> bool:
    """验证缓存的 checksum"""
    manifest_path = cache_dir / accession / CACHE_MANIFEST

    try:
        with open(manifest_path) as f:
            manifest = json.load(f)

        for file_info in manifest.get("files", []):
            file_path = Path(file_info["path"])
            expected_sha256 = file_info.get("sha256")

            if expected_sha256:
                actual_sha256 = compute_sha256(file_path)
                if actual_sha256 != expected_sha256:
                    return False

        return True
    except Exception:
        return False

def write_download_manifest(
    accession: str,
    output_dir: Path,
    files: list[dict],
    metadata: Optional[dict] = None
) -> Path:
    """写入下载清单"""
    manifest = {
        "accession": accession,
        "downloaded_at": datetime.utcnow().isoformat() + "Z",
        "files": files,
        "metadata": metadata or {}
    }

    manifest_path = output_dir / CACHE_MANIFEST
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    return manifest_path
```

### Snakemake 规则 [Source: architecture.md#Naming-Patterns]

```python
# rules/download.smk

from workflow.lib.ncbi import (
    parse_ncbi_accession,
    NCBIClient,
    is_cached,
    download_ncbi_dataset
)

def get_ncbi_species():
    """获取配置中使用 NCBI 标识符的物种"""
    return [
        sp for sp in config["species"]
        if str(sp.get("annotation", "")).upper().startswith("NCBI:")
    ]

rule ncbi_download:
    """下载单个 NCBI 数据集"""
    output:
        manifest = "{output_dir}/downloads/ncbi/{accession}/download_manifest.json",
        genome = "{output_dir}/downloads/ncbi/{accession}/genome.fna.gz",
        annotation = "{output_dir}/downloads/ncbi/{accession}/annotation.gff.gz"
    params:
        force = config.get("ncbi", {}).get("force_download", False),
        api_key = os.environ.get("NCBI_API_KEY")
    log:
        "{output_dir}/logs/download/{accession}.log"
    run:
        from workflow.lib.ncbi import download_ncbi_dataset

        output_dir = Path(output.manifest).parent
        download_ncbi_dataset(
            accession=wildcards.accession,
            output_dir=output_dir,
            force=params.force,
            api_key=params.api_key
        )

def get_ncbi_download_targets():
    """生成所有 NCBI 下载目标"""
    targets = []
    for sp in get_ncbi_species():
        accession = sp["annotation"].replace("NCBI:", "")
        targets.append(
            f"{config['output_dir']}/downloads/ncbi/{accession}/download_manifest.json"
        )
    return targets

rule ncbi_download_all:
    """批量下载所有 NCBI 数据集"""
    input:
        get_ncbi_download_targets()
    output:
        touch("{output_dir}/downloads/ncbi/.downloads_complete")
```

### 配置 Schema 更新 [Source: schemas/config.schema.yaml]

```yaml
# 新增 ncbi 配置节
ncbi:
  type: object
  properties:
    api_key:
      type: string
      description: "NCBI API key（可选，提升速率限制）"
    rate_limit:
      type: number
      default: 3.0
      description: "每秒请求数限制"
    cache_dir:
      type: string
      default: "results/downloads/ncbi"
      description: "下载缓存目录"
    force_download:
      type: boolean
      default: false
      description: "强制重新下载，忽略缓存"
  additionalProperties: false

# species.annotation 支持 NCBI 格式
species:
  items:
    properties:
      annotation:
        oneOf:
          - type: string
            pattern: "^/.*"  # 本地路径
          - type: string
            pattern: "^NCBI:GC[FA]_\\d{9}(\\.\\d+)?$"  # NCBI accession
```

### 与标准化流程集成

**数据流更新：**
```
配置文件 (annotation: "NCBI:GCF_xxx")
    ↓
[ncbi_download] → downloads/ncbi/{accession}/
    ↓                ├── genome.fna.gz
    ↓                ├── annotation.gff.gz
    ↓                └── download_manifest.json
    ↓
[standardize_species] → standardized/{species}/
                         ├── genome.fa.gz
                         └── annotation.gff3.gz
```

**standardize.smk 更新：**
```python
def get_annotation_input(wildcards):
    """获取物种的注释文件输入（本地或 NCBI）"""
    sp_config = get_species_config(wildcards.species)
    annotation = sp_config.get("annotation", "")

    if annotation.upper().startswith("NCBI:"):
        accession = annotation.replace("NCBI:", "")
        return f"{config['output_dir']}/downloads/ncbi/{accession}/annotation.gff.gz"
    else:
        return annotation  # 本地路径
```

### 错误处理 [Source: architecture.md#ADR-003]

| 场景 | 错误码 | 可重试 | 恢复建议 |
|------|--------|--------|----------|
| 无效 accession 格式 | E_INPUT_FORMAT | ❌ | 检查 accession 格式（GCF_/GCA_ + 9位数字） |
| 网络超时 | E_TIMEOUT | ✅ | 检查网络连接，稍后重试 |
| 速率限制 | E_NET_RATE_LIMIT | ✅ | 等待后自动重试，或使用 API key |
| 文件不存在 | E_INPUT_MISSING | ❌ | 检查 accession 是否存在于 NCBI |
| 下载中断 | E_NET_ERROR | ✅ | 自动断点续传重试 |

### Project Structure Notes

**新增文件：**
```
workflow/
├── lib/
│   ├── ncbi.py              # NCBI 客户端（新建）
│   └── test_ncbi.py         # 单元测试（新建）
├── rules/
│   └── download.smk         # 下载规则（新建）
```

**修改文件：**
```
config/
└── config.yaml              # 添加 ncbi 配置节示例
schemas/
└── config.schema.yaml       # 添加 ncbi 配置 schema
workflow/rules/
└── standardize.smk          # 集成 NCBI 下载输入
```

### 从前序 Story 学到的经验 [Source: 2-2b-biopython-utils.md]

1. **延迟导入模式**：对可选依赖（如 NCBI API 客户端库）使用延迟导入，支持降级
2. **dataclass 使用**：使用 dataclass 定义结构化数据（如 NCBIAccession、DownloadResult）
3. **错误处理一致性**：所有错误应使用 CompGeneError 体系
4. **单元测试共置**：test_ncbi.py 与 ncbi.py 放在同一目录
5. **原子写入**：下载时先写 .tmp 再 rename，避免不完整文件

### References

- [Source: docs/planning-artifacts/prd.md#FR1] - 用户可通过 NCBI 物种标识符自动下载基因组注释数据
- [Source: docs/planning-artifacts/prd.md#FR4] - 系统可缓存已下载的注释数据，避免重复下载
- [Source: docs/planning-artifacts/architecture.md#NFR9] - NCBI API 调用应遵循官方速率限制（每秒 3 请求）
- [Source: docs/planning-artifacts/architecture.md#ADR-003] - 错误处理与恢复策略
- [NCBI Datasets API](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/rest-api/) - REST API 参考
- [NCBI FTP 结构](https://ftp.ncbi.nlm.nih.gov/genomes/all/) - FTP 目录结构

## Dev Agent Record

### Agent Model Used

Claude Opus 4.5 (claude-opus-4-5-20251101)

### Debug Log References

None - all tests passed without issues.

### Completion Notes List

1. Implemented complete NCBI download module (`workflow/lib/ncbi.py`) with:
   - AccessionType enum and NCBIAccession dataclass for parsing GCF_/GCA_ accessions
   - RateLimiter class using token bucket algorithm (3 req/s default, 10 req/s with API key)
   - Exponential backoff retry logic via `with_retry()` function
   - FTP directory resolution and file URL construction
   - SHA256 checksum computation during download
   - Atomic file writes using .tmp + rename pattern
   - Cache mechanism with manifest-based validation
   - NCBIClient high-level class interface

2. Created comprehensive test suite (`workflow/lib/test_ncbi.py`) with 54 tests covering:
   - Accession parsing (valid/invalid formats, case sensitivity)
   - Rate limiter behavior (tokens, refill, thread safety)
   - Retry logic (exponential backoff, max delay cap)
   - API key handling (env var, URL building)
   - Cache management (existence check, checksum validation)
   - Integration tests with mocked network calls

3. Created Snakemake rules (`workflow/rules/download.smk`) with:
   - `ncbi_download` rule for single accession download
   - `ncbi_download_all` aggregate rule for batch downloads
   - `standardize_ncbi_genome` and `standardize_ncbi_annotation` for integration
   - Helper functions for NCBI species detection

4. Updated config schema (`schemas/config.schema.yaml`) with:
   - New `ncbi` configuration section (api_key, rate_limit, cache_dir, force_download)
   - Enhanced annotation field description for NCBI accession support

5. Updated main Snakefile to include `download.smk` module

### File List

**New Files:**
- `workflow/lib/ncbi.py` (1179 lines) - NCBI download module
- `workflow/lib/test_ncbi.py` (790 lines) - Unit tests
- `workflow/rules/download.smk` (277 lines) - Snakemake download rules

**Modified Files:**
- `workflow/Snakefile` - Added include for download.smk
- `schemas/config.schema.yaml` - Added ncbi configuration section

