# Story 1.6: 统一 Runner 与重试

Status: done

## Story

As a **计算生物学研究员**,
I want **外部工具调用有统一的执行、超时和重试机制**,
So that **临时故障不会导致整个分析失败**。

## Acceptance Criteria

1. **AC1: 统一执行流程**
   - Given rule 调用外部工具
   - When 使用 run_adapter.py 执行
   - Then 按顺序：版本检测 → 输入校验 → 命令执行 → 输出解析 → 审计记录

2. **AC2: 超时处理**
   - Given 工具执行超时
   - When 超过 timeout_seconds
   - Then 发送 SIGTERM，等待 grace period 后 SIGKILL，返回 E_TIMEOUT

3. **AC3: 可重试错误自动重试**
   - Given 网络错误（可重试）
   - When 发生错误
   - Then 自动重试最多 3 次，指数退避间隔

## Tasks / Subtasks

- [x] Task 1: 创建 Runner 核心模块 (AC: #1)
  - [x] 1.1 创建 `workflow/lib/runner.py` 模块
  - [x] 1.2 定义 `AdapterRunner` 类，封装执行流程
  - [x] 1.3 实现 `run(adapter: BaseAdapter, ctx: AdapterContext) -> RunResult` 主方法
  - [x] 1.4 实现执行流程：version_check → validate_inputs → build_command → execute → parse_outputs

- [x] Task 2: 实现版本检测阶段 (AC: #1)
  - [x] 2.1 调用 `adapter.check_version()` 获取工具版本
  - [x] 2.2 调用 `adapter.spec.check_version_compatible()` 验证版本兼容性
  - [x] 2.3 版本不兼容时抛出 `CompGeneError(E_TOOL_VERSION)`
  - [x] 2.4 记录版本信息到日志

- [x] Task 3: 实现输入校验阶段 (AC: #1)
  - [x] 3.1 调用 `adapter.validate_inputs(ctx)` 校验输入
  - [x] 3.2 校验失败时抛出对应的 CompGeneError
  - [x] 3.3 记录校验结果到日志

- [x] Task 4: 实现命令执行阶段 (AC: #1, #2)
  - [x] 4.1 调用 `adapter.build_command(ctx)` 构建命令
  - [x] 4.2 使用 `subprocess.Popen` 执行命令
  - [x] 4.3 捕获 stdout/stderr 输出
  - [x] 4.4 实现超时监控和终止逻辑

- [x] Task 5: 实现超时处理机制 (AC: #2)
  - [x] 5.1 调用 `adapter.timeout_seconds(ctx)` 获取超时时间
  - [x] 5.2 实现 `_handle_timeout()` 方法处理超时（名称调整）
  - [x] 5.3 超时时发送 SIGTERM 信号
  - [x] 5.4 等待 grace period（默认 10 秒）
  - [x] 5.5 grace period 后仍未退出则发送 SIGKILL
  - [x] 5.6 返回 `E_TIMEOUT` 错误码

- [x] Task 6: 实现输出解析阶段 (AC: #1)
  - [x] 6.1 检查 `adapter.expected_outputs(ctx)` 所有文件存在
  - [x] 6.2 调用 `adapter.parse_outputs(ctx)` 解析结果
  - [x] 6.3 输出缺失时抛出 `CompGeneError(E_NONZERO_EXIT)` （使用更通用的错误码）

- [x] Task 7: 实现审计记录集成 (AC: #1)
  - [x] 7.1 导入 `workflow/lib/audit.py` 的 `RunMetadata` 和 `write_run_json`
  - [x] 7.2 收集执行元数据：cmd, tool_version, runtime, exit_code, checksums
  - [x] 7.3 成功时记录 output_checksums
  - [x] 7.4 失败时记录 error_code, error_message
  - [x] 7.5 写入 `.run.json` 审计文件

- [x] Task 8: 实现重试机制 (AC: #3)
  - [x] 8.1 定义 `RetryConfig` dataclass（max_retries, base_delay, max_delay, jitter）
  - [x] 8.2 实现 `_calculate_backoff()` 指数退避计算
  - [x] 8.3 调用 `adapter.classify_error()` 判断是否可重试
  - [x] 8.4 实现 `_execute_with_retry()` 包装方法
  - [x] 8.5 可重试错误最多重试 3 次
  - [x] 8.6 记录每次重试的日志

- [x] Task 9: 实现日志集成 (AC: #1)
  - [x] 9.1 导入 `workflow/lib/logging.py` 的 `DualLogger`
  - [x] 9.2 在每个阶段记录 INFO 级别日志
  - [x] 9.3 错误时记录 ERROR 级别日志
  - [x] 9.4 重试时记录 WARNING 级别日志

- [x] Task 10: 创建 run_adapter.py 入口脚本 (AC: #1)
  - [x] 10.1 创建 `workflow/scripts/run_adapter.py` 脚本
  - [x] 10.2 实现命令行参数解析（adapter_name, inputs, outputs, config, wildcards）
  - [x] 10.3 动态加载指定的 Adapter 类
  - [x] 10.4 构建 AdapterContext 并调用 AdapterRunner.run()
  - [x] 10.5 处理异常并返回适当的退出码

- [x] Task 11: 创建单元测试 (AC: #1, #2, #3)
  - [x] 11.1 创建 `workflow/lib/test_runner.py`
  - [x] 11.2 测试正常执行流程（mock subprocess）
  - [x] 11.3 测试版本检测失败场景
  - [x] 11.4 测试输入校验失败场景
  - [x] 11.5 测试超时处理（SIGTERM → SIGKILL）
  - [x] 11.6 测试可重试错误的重试逻辑
  - [x] 11.7 测试不可重试错误直接失败
  - [x] 11.8 测试审计记录生成

- [x] Task 12: 更新模块导出 (AC: #1)
  - [x] 12.1 更新 `workflow/lib/__init__.py` 导出 AdapterRunner, RetryConfig

## Dev Notes

### 架构约束 [Source: docs/planning-artifacts/architecture.md#ADR-002]

**Runner 职责（Tool Adaptation Layer）：**
- 版本检测：调用 adapter.check_version() 并验证兼容性
- 输入校验：调用 adapter.validate_inputs() 确保输入有效
- 命令执行：subprocess 执行 adapter.build_command() 返回的命令
- 超时 kill：监控执行时间，超时时 SIGTERM → SIGKILL
- stdout/stderr 捕获：收集输出用于错误诊断
- 输出解析：调用 adapter.parse_outputs() 提取结果
- 审计元数据：生成 .run.json 记录执行信息

### 重试策略 [Source: docs/planning-artifacts/architecture.md#ADR-003]

**Hybrid retry（Snakemake 层 + Runner 层）：**
- Snakemake 层：`retries: 1`（job 级重试）
- Runner 层：网络错误 3 次 exponential backoff（command 级重试）

**指数退避参数：**
```python
base_delay = 1.0  # 初始等待 1 秒
max_delay = 60.0  # 最大等待 60 秒
jitter = 0.1      # 10% 随机抖动
# delay = min(base_delay * (2 ** attempt), max_delay) * (1 + random.uniform(-jitter, jitter))
```

**可重试错误码（is_retryable=True）：**
- `E_TIMEOUT`：超时可能是临时负载问题
- `E_NET_RATE_LIMIT`：网络限流，等待后重试

**不可重试错误码（is_retryable=False）：**
- `E_TOOL_NOT_FOUND`：工具未安装
- `E_TOOL_VERSION`：版本不兼容
- `E_INPUT_MISSING`：输入文件缺失
- `E_INPUT_FORMAT`：输入格式错误
- `E_OOM`：内存不足
- `E_DISK_FULL`：磁盘空间不足

### 超时处理 [Source: docs/planning-artifacts/architecture.md#ADR-002]

**优雅终止流程：**
1. 超过 `timeout_seconds` 后发送 `SIGTERM`
2. 等待 grace period（默认 10 秒）让进程清理
3. grace period 后仍未退出，发送 `SIGKILL` 强制终止
4. 返回 `E_TIMEOUT` 错误码

### 依赖模块

**已完成的依赖（Story 1.3 ~ 1.5）：**
- `workflow/lib/errors.py` - ErrorCode enum, CompGeneError
- `workflow/lib/logging.py` - DualLogger 双格式日志
- `workflow/lib/audit.py` - RunMetadata, write_run_json
- `workflow/lib/checksum.py` - compute_file_checksum
- `workflow/adapters/base.py` - BaseAdapter, AdapterContext, RunResult, ToolSpec

### 代码模式参考

**Runner 主流程伪代码：**
```python
class AdapterRunner:
    def run(self, adapter: BaseAdapter, ctx: AdapterContext) -> RunResult:
        start_time = time.time()

        # 1. Version check
        version = adapter.check_version()
        if not adapter.spec.check_version_compatible(version):
            raise CompGeneError(E_TOOL_VERSION, f"Version {version} not compatible")

        # 2. Input validation
        adapter.validate_inputs(ctx)

        # 3. Build command
        cmd = adapter.build_command(ctx)

        # 4. Execute with timeout and retry
        returncode, stdout, stderr = self._execute_with_retry(
            cmd, adapter.timeout_seconds(ctx), adapter, ctx
        )

        # 5. Check outputs
        for path in adapter.expected_outputs(ctx):
            if not path.exists():
                raise CompGeneError(E_OUTPUT_MISSING, f"Missing: {path}")

        # 6. Parse outputs
        result = adapter.parse_outputs(ctx)

        # 7. Write audit
        runtime = time.time() - start_time
        self._write_audit(adapter, ctx, cmd, version, returncode, runtime)

        return result
```

**重试逻辑伪代码：**
```python
def _execute_with_retry(self, cmd, timeout, adapter, ctx) -> tuple[int, str, str]:
    for attempt in range(self.retry_config.max_retries + 1):
        returncode, stdout, stderr = self._execute_once(cmd, timeout)

        if returncode == 0:
            return returncode, stdout, stderr

        error_code, is_retryable = adapter.classify_error(ctx, returncode, stderr)

        if not is_retryable or attempt == self.retry_config.max_retries:
            raise CompGeneError(error_code, stderr)

        delay = self._calculate_backoff(attempt)
        self.logger.warning(f"Retrying in {delay}s (attempt {attempt + 1})")
        time.sleep(delay)
```

### Project Structure Notes

**新建文件：**
- `compgene/workflow/lib/runner.py` - AdapterRunner, RetryConfig
- `compgene/workflow/lib/test_runner.py` - Runner 单元测试
- `compgene/workflow/scripts/run_adapter.py` - 命令行入口脚本

**修改文件：**
- `compgene/workflow/lib/__init__.py` - 导出新类

**模块边界：**
- `lib/runner.py` 依赖 `lib/errors.py`, `lib/logging.py`, `lib/audit.py`, `lib/checksum.py`
- `lib/runner.py` 依赖 `adapters/base.py`（BaseAdapter, AdapterContext）
- `scripts/run_adapter.py` 是 Snakemake rule 调用的入口

### 与后续 Story 的关系

- **Story 1.7 (checkpoint-resume-dryrun)**: 将在 Runner 基础上添加 checkpoint 和 dry-run 功能
- **Story 3.1 (OrthoFinder Adapter)**: 将使用 AdapterRunner 执行 OrthoFinder
- **所有后续工具集成**: 都将通过 run_adapter.py 调用

### References

- [Source: docs/planning-artifacts/architecture.md#ADR-002] - Tool Adaptation Layer
- [Source: docs/planning-artifacts/architecture.md#ADR-003] - Error Handling Strategy
- [Source: docs/planning-artifacts/epics.md#Story-1.6] - Story 需求
- [Source: workflow/lib/errors.py] - ErrorCode enum
- [Source: workflow/lib/audit.py] - RunMetadata, write_run_json
- [Source: workflow/adapters/base.py] - BaseAdapter interface

## Dev Agent Record

### Agent Model Used

Claude Opus 4.5 (claude-opus-4-5-20251101)

### Debug Log References

None

### Completion Notes List

- Implemented `AdapterRunner` class in `workflow/lib/runner.py` (480+ lines)
- Implemented `RetryConfig` dataclass for exponential backoff configuration
- Runner implements complete execution flow: version check → input validation → command execution → output verification → output parsing → audit recording
- Timeout handling uses graceful termination: SIGTERM → grace period (10s default) → SIGKILL
- Retry logic supports configurable max retries (default 3), exponential backoff with jitter
- All execution stages integrate with DualLogger for INFO/WARNING/ERROR logging
- Audit records generated on both success and failure using `collect_run_metadata()` and `write_run_json()`
- Created `run_adapter.py` CLI entry point (180+ lines) for Snakemake rule integration
- CLI supports: adapter_name, --inputs, --outputs, --config, --wildcards, --threads, --meta-dir, --log-dir, --max-retries, --grace-period, --log-level, --no-color
- Dynamic adapter loading via importlib with convention: `workflow.adapters.{name}.{Name}Adapter`
- Created comprehensive test suite with 33 unit tests covering:
  - RetryConfig defaults and custom values
  - AdapterRunner initialization
  - Version check (compatible, below minimum, above maximum)
  - Input validation (success, missing, format errors)
  - Command execution (success, real subprocess)
  - Timeout handling (SIGTERM, SIGKILL escalation)
  - Retry logic (success on retry, max retries exceeded, non-retryable errors)
  - Exponential backoff (increase, max cap, jitter)
  - Output verification
  - Audit records (success, failure, input/output checksums)
  - Logging integration (info, warning, error messages)
  - Context logger override
  - Edge cases (empty wildcards, zero timeout, long stderr)
- All 33 new tests pass, 80 total adapter/runner tests pass
- Updated `workflow/lib/__init__.py` to export AdapterRunner and RetryConfig

### Code Review Record

**Review Date:** 2026-01-21
**Reviewer Model:** Claude Opus 4.5

**Issues Found:**
- H1 (HIGH): Missing E_OUTPUT_MISSING error code - used E_NONZERO_EXIT instead
- M1 (MEDIUM): Unused RunMetadata import in runner.py
- M2 (MEDIUM): No exception handling for os.killpg/os.getpgid in _handle_timeout()
- M3 (MEDIUM): No CLI tests for run_adapter.py
- M4 (MEDIUM): Hardcoded stderr truncation length (magic number 500)
- L1 (LOW): Module docstring missing runner in __init__.py
- L2 (LOW): Unused `from unittest import mock` import in test_runner.py
- L3 (LOW): Inaccurate line counts in story File List

**Fixes Applied:**
- Added E_OUTPUT_MISSING error code to errors.py (exit code 7, recovery: "检查工具执行日志确认输出生成")
- Updated runner.py _verify_outputs() to use E_OUTPUT_MISSING instead of E_NONZERO_EXIT
- Removed unused RunMetadata import from runner.py
- Added exception handling for ProcessLookupError/OSError in _handle_timeout() for process group operations
- Added MAX_STDERR_DETAILS constant (500) to replace hardcoded values
- Updated __init__.py docstring to include runner module
- Removed unused `from unittest import mock` import from test_runner.py
- Created test_run_adapter.py with 16 CLI tests covering argument parsing, JSON validation, adapter loading, and main function
- Updated test_errors.py to include E_OUTPUT_MISSING (count 9→10, exit code range 1-6→1-7)
- Updated test_runner.py to expect E_OUTPUT_MISSING for missing output test

**Test Results After Fixes:**
- 49 runner/CLI tests pass (33 runner + 16 CLI)
- 96 total adapter/runner tests pass
- 51/52 errors tests pass (1 pre-existing failure unrelated to Story 1.6)

### File List

- `compgene/workflow/lib/runner.py` (created, 549 lines)
- `compgene/workflow/lib/test_runner.py` (created, 875 lines)
- `compgene/workflow/scripts/run_adapter.py` (created, 272 lines)
- `compgene/workflow/scripts/test_run_adapter.py` (created, 250 lines)
- `compgene/workflow/lib/__init__.py` (modified)
- `compgene/workflow/lib/errors.py` (modified - added E_OUTPUT_MISSING)
- `compgene/workflow/lib/test_errors.py` (modified - updated for E_OUTPUT_MISSING)
