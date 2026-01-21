# Story 1.4: 日志与审计基座

Status: ready-for-dev

## Story

As a **计算生物学研究员**,
I want **分析过程自动生成结构化日志和审计记录**,
so that **我可以追溯任何输出的生成过程**.

## Acceptance Criteria

1. **Given** 任意 rule 执行完成
   **When** Runner 收集元数据
   **Then** 生成双格式日志：`.log`（文本）+ `.jsonl`（JSON Lines）
   **And** 生成 `.run.json` 包含命令、版本、checksums、时间戳

## Tasks / Subtasks

- [ ] Task 1: 创建日志模块基础设施 (AC: #1)
  - [ ] 1.1 创建 `workflow/lib/logging.py` 模块
  - [ ] 1.2 实现 `DualLogger` 类，同时输出 `.log` 和 `.jsonl`
  - [ ] 1.3 实现日志级别配置 (DEBUG/INFO/WARNING/ERROR)
  - [ ] 1.4 实现日志路径生成函数 `get_log_path(rule, wildcards) -> tuple[Path, Path]`

- [ ] Task 2: 实现文本日志格式 (AC: #1)
  - [ ] 2.1 格式：`YYYY-MM-DD HH:MM:SS [LEVEL] {rule}: {message}`
  - [ ] 2.2 支持 ANSI 颜色输出（可配置禁用）
  - [ ] 2.3 日志文件路径：`logs/{rule}/{wildcards}.log`

- [ ] Task 3: 实现 JSON Lines 日志格式 (AC: #1)
  - [ ] 3.1 格式：每行一个 JSON 对象，包含 ts, level, rule, msg, 额外字段
  - [ ] 3.2 日志文件路径：`logs/{rule}/{wildcards}.jsonl`
  - [ ] 3.3 支持 jq/pandas 直接查询

- [ ] Task 4: 创建审计模块 (AC: #1)
  - [ ] 4.1 创建 `workflow/lib/audit.py` 模块
  - [ ] 4.2 定义 `RunMetadata` 数据类，包含所有审计字段
  - [ ] 4.3 实现 `collect_run_metadata()` 函数收集执行元数据
  - [ ] 4.4 实现 `write_run_json()` 函数生成 `.run.json` 文件

- [ ] Task 5: 实现 checksum 计算 (AC: #1)
  - [ ] 5.1 实现 `compute_file_checksum(path: Path, algorithm: str = "sha256") -> str`
  - [ ] 5.2 实现 `compute_dir_checksum(path: Path) -> str` 用于目录
  - [ ] 5.3 支持大文件分块计算（避免内存溢出）

- [ ] Task 6: 定义 .run.json Schema (AC: #1)
  - [ ] 6.1 创建 `schemas/run.schema.json` 定义审计格式
  - [ ] 6.2 必需字段：rule, wildcards, cmd, tool_version, input_checksums, threads, runtime_seconds, exit_code, timestamp
  - [ ] 6.3 可选字段：output_checksums, error_code, error_message

- [ ] Task 7: 创建原子写入工具 (AC: #1)
  - [ ] 7.1 创建 `workflow/lib/io.py` 模块
  - [ ] 7.2 实现 `atomic_write(path: Path, content: str)` 函数
  - [ ] 7.3 实现 `atomic_write_json(path: Path, data: dict)` 函数
  - [ ] 7.4 写入模式：先写 `.tmp` 文件，再原子 rename

- [ ] Task 8: 创建单元测试 (AC: #1)
  - [ ] 8.1 创建 `workflow/lib/test_logging.py`
  - [ ] 8.2 创建 `workflow/lib/test_audit.py`
  - [ ] 8.3 创建 `workflow/lib/test_io.py`
  - [ ] 8.4 测试日志格式、审计元数据收集、checksum 计算、原子写入

## Dev Notes

### 日志格式规范 [Source: docs/planning-artifacts/architecture.md#ADR-004]

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

### 审计文件格式 [Source: docs/planning-artifacts/architecture.md#ADR-004]

**审计文件路径：** `results/meta/{rule}/{wildcards}.run.json`

**示例 .run.json：**
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

### 命名规范 [Source: docs/planning-artifacts/architecture.md#Naming-Patterns]

- Python 命名：snake_case（函数/变量）、PascalCase（类）、UPPER_SNAKE_CASE（常量）
- 日志路径：`logs/{rule}/{wildcards}.log` 和 `logs/{rule}/{wildcards}.jsonl`
- 审计路径：`results/meta/{rule}/{wildcards}.run.json`
- 日期时间：ISO 8601 格式（UTC）

### 关键约束

1. **双格式日志**：每次日志调用同时写入 `.log` 和 `.jsonl`
2. **原子写入**：所有文件写入使用先 `.tmp` 再 rename 模式
3. **大文件处理**：checksum 计算使用分块读取，避免内存溢出
4. **日志级别**：支持 DEBUG/INFO/WARNING/ERROR，通过配置控制
5. **无外部依赖**：仅使用 Python 标准库（logging, json, hashlib, pathlib）

### 从 Story 1.3 继承的上下文

Story 1.3 已实现：
- `workflow/lib/errors.py` - ErrorCode enum 和 CompGeneError 异常类
- 错误码可用于审计记录中的 `error_code` 字段

本 Story 需要：
- 与 errors.py 集成，在审计记录中包含错误码信息
- 保持 `lib/` 模块边界：lib/ 不依赖其他项目模块（errors.py 除外）

### 与后续 Story 的关系

- **Story 1.5 (BaseAdapter)**: 将使用 audit 模块记录工具执行元数据
- **Story 1.6 (统一 Runner)**: 将调用 DualLogger 记录执行过程
- **Story 8.1 (SQLite 审计汇总)**: 将读取所有 `.run.json` 文件写入 SQLite

### Project Structure Notes

**新建文件：**
- `compgene/workflow/lib/logging.py` - 双格式日志模块
- `compgene/workflow/lib/audit.py` - 审计元数据收集
- `compgene/workflow/lib/io.py` - 原子写入工具
- `compgene/workflow/lib/test_logging.py` - 日志测试
- `compgene/workflow/lib/test_audit.py` - 审计测试
- `compgene/workflow/lib/test_io.py` - IO 测试
- `compgene/schemas/run.schema.json` - 审计文件 schema

**修改文件：**
- `compgene/workflow/lib/__init__.py` - 导出新模块

**模块边界：**
- `lib/` 只依赖 `lib/errors.py`（同层）
- 不依赖 `adapters/`, `scripts/`, `rules/`

### 代码骨架参考

**DualLogger 类骨架：**
```python
# workflow/lib/logging.py
import logging
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

class DualLogger:
    """同时输出文本和 JSON Lines 格式的日志器。"""

    def __init__(
        self,
        rule: str,
        wildcards: dict[str, str],
        log_dir: Path = Path("logs"),
        level: int = logging.INFO,
        use_color: bool = True
    ):
        self.rule = rule
        self.wildcards = wildcards
        self.log_path = self._get_log_path(log_dir, ".log")
        self.jsonl_path = self._get_log_path(log_dir, ".jsonl")
        # ... setup handlers

    def info(self, msg: str, **extra: Any) -> None: ...
    def warning(self, msg: str, **extra: Any) -> None: ...
    def error(self, msg: str, **extra: Any) -> None: ...
    def debug(self, msg: str, **extra: Any) -> None: ...
```

**RunMetadata 数据类骨架：**
```python
# workflow/lib/audit.py
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Optional

@dataclass
class RunMetadata:
    """rule 执行的审计元数据。"""
    rule: str
    wildcards: dict[str, str]
    cmd: list[str]
    tool_version: str
    input_checksums: dict[str, str]
    threads: int
    runtime_seconds: float
    exit_code: int
    timestamp: str  # ISO 8601
    output_checksums: Optional[dict[str, str]] = None
    error_code: Optional[str] = None
    error_message: Optional[str] = None

    def to_dict(self) -> dict:
        """转换为可 JSON 序列化的字典。"""
        return {k: v for k, v in asdict(self).items() if v is not None}
```

**原子写入骨架：**
```python
# workflow/lib/io.py
import json
from pathlib import Path
from typing import Any

def atomic_write(path: Path, content: str, encoding: str = "utf-8") -> None:
    """原子写入文本文件。"""
    temp_path = path.with_suffix(path.suffix + ".tmp")
    temp_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path.write_text(content, encoding=encoding)
    temp_path.rename(path)

def atomic_write_json(path: Path, data: Any, indent: int = 2) -> None:
    """原子写入 JSON 文件。"""
    content = json.dumps(data, ensure_ascii=False, indent=indent)
    atomic_write(path, content + "\n")
```

### References

- [Source: docs/planning-artifacts/architecture.md#ADR-004]
- [Source: docs/planning-artifacts/architecture.md#Log-Format]
- [Source: docs/planning-artifacts/architecture.md#Audit-Granularity]
- [Source: docs/planning-artifacts/epics.md#Story-1.4]
- [Source: docs/planning-artifacts/prd.md#FR45-46]
- [Source: docs/planning-artifacts/prd.md#NFR17] (日志级别可配置)

## Dev Agent Record

### Agent Model Used

### Debug Log References

### Completion Notes List

### File List
