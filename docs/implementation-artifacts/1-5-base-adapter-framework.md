# Story 1.5: BaseAdapter 框架

Status: done

## Story

As a **开发者**,
I want **统一的 Adapter 抽象接口**,
So that **后续工具集成遵循一致模式**。

## Acceptance Criteria

1. **AC1: BaseAdapter 抽象类定义**
   - Given 需要集成新的外部工具
   - When 继承 BaseAdapter
   - Then 必须实现 8 个方法：check_version, validate_inputs, build_command, expected_outputs, parse_outputs, timeout_seconds, classify_error, spec

2. **AC2: 可独立测试**
   - Given Adapter 实现完成
   - When 运行单元测试
   - Then 可以独立测试（mock 外部工具）

3. **AC3: ToolSpec 数据结构**
   - Given 任意 Adapter
   - When 访问 spec 属性
   - Then 返回 ToolSpec 包含：name, min_version, max_version, conda_env

4. **AC4: RunResult 数据结构**
   - Given Adapter 执行完成
   - When parse_outputs 返回
   - Then 返回 RunResult 包含：outputs, summary, metrics

5. **AC5: 错误分类集成**
   - Given 工具执行失败
   - When classify_error 被调用
   - Then 返回 (ErrorCode, is_retryable) 使用 lib/errors.py 的错误码

## Tasks / Subtasks

- [x] Task 1: 创建 ToolSpec 数据结构 (AC: #3)
  - [x] 1.1 定义 ToolSpec dataclass（name, min_version, max_version, conda_env, description）
  - [x] 1.2 实现版本比较辅助方法
  - [x] 1.3 编写 ToolSpec 单元测试

- [x] Task 2: 创建 RunResult 数据结构 (AC: #4)
  - [x] 2.1 定义 RunResult dataclass（outputs, summary, metrics, warnings）
  - [x] 2.2 实现 to_dict() 序列化方法
  - [x] 2.3 编写 RunResult 单元测试

- [x] Task 3: 实现 BaseAdapter 抽象类 (AC: #1)
  - [x] 3.1 定义抽象基类继承 ABC
  - [x] 3.2 定义 spec 抽象属性
  - [x] 3.3 定义 check_version() 抽象方法
  - [x] 3.4 定义 validate_inputs() 抽象方法
  - [x] 3.5 定义 build_command() 抽象方法
  - [x] 3.6 定义 expected_outputs() 抽象方法
  - [x] 3.7 定义 parse_outputs() 抽象方法
  - [x] 3.8 定义 timeout_seconds() 抽象方法
  - [x] 3.9 定义 classify_error() 抽象方法

- [x] Task 4: 实现 AdapterContext 上下文类 (AC: #1)
  - [x] 4.1 定义 AdapterContext dataclass（inputs, outputs, config, wildcards, threads, logger）
  - [x] 4.2 实现路径解析辅助方法
  - [x] 4.3 编写 AdapterContext 单元测试

- [x] Task 5: 创建 MockAdapter 测试实现 (AC: #2)
  - [x] 5.1 创建 MockAdapter 继承 BaseAdapter
  - [x] 5.2 实现所有 8 个必须方法
  - [x] 5.3 验证接口完整性
  - [x] 5.4 编写 MockAdapter 单元测试

- [x] Task 6: 错误分类集成 (AC: #5)
  - [x] 6.1 导入 lib/errors.py 的 ErrorCode
  - [x] 6.2 定义 classify_error 返回类型为 tuple[ErrorCode, bool]
  - [x] 6.3 实现常见错误模式的分类辅助函数
  - [x] 6.4 编写错误分类单元测试

- [x] Task 7: 更新模块导出 (AC: #1, #2)
  - [x] 7.1 更新 workflow/adapters/__init__.py 导出所有公共类
  - [x] 7.2 验证模块可以正常导入

## Dev Notes

### 架构约束 [Source: docs/planning-artifacts/architecture.md#ADR-002]

**Adapter Interface 定义：**
```python
class BaseAdapter(ABC):
    spec: ToolSpec
    def check_version(self) -> str
    def validate_inputs(self, ctx: AdapterContext) -> None
    def build_command(self, ctx: AdapterContext) -> list[str]
    def expected_outputs(self, ctx: AdapterContext) -> list[Path]
    def parse_outputs(self, ctx: AdapterContext) -> RunResult
    def timeout_seconds(self, ctx: AdapterContext) -> int
    def classify_error(self, ctx: AdapterContext, returncode: int, stderr: str) -> tuple[ErrorCode, bool]
```

**关键设计决策：**
- classify_error 返回 (ErrorCode, is_retryable) 元组
- 所有方法接受 AdapterContext 参数提供运行时上下文
- ToolSpec 包含版本约束用于版本检查

### 依赖模块 [Source: workflow/lib/]

**已完成的依赖（Story 1.3 + 1.4）：**
- `workflow/lib/errors.py` - ErrorCode enum, CompGeneError
- `workflow/lib/logging.py` - DualLogger（可选，用于 AdapterContext）
- `workflow/lib/audit.py` - RunMetadata（参考数据结构设计）
- `workflow/lib/io.py` - 原子写入（测试时可能用到）

### 代码模式参考

**参考 RunMetadata (Story 1.4) 的 dataclass 模式：**
```python
@dataclass
class RunMetadata:
    rule: str
    wildcards: dict[str, str]
    cmd: list[str]
    ...
    def to_dict(self) -> dict[str, Any]: ...
    def to_json(self) -> str: ...
```

**错误码使用模式 (Story 1.3)：**
```python
from workflow.lib.errors import ErrorCode, CompGeneError

def classify_error(...) -> tuple[ErrorCode, bool]:
    if "timeout" in stderr.lower():
        return (ErrorCode.E_TIMEOUT, True)
    if returncode == 127:
        return (ErrorCode.E_TOOL_NOT_FOUND, False)
    return (ErrorCode.E_NONZERO_EXIT, False)
```

### Project Structure Notes

**文件位置：**
- `workflow/adapters/base.py` - BaseAdapter, ToolSpec, RunResult, AdapterContext
- `workflow/adapters/test_base.py` - 单元测试（共置模式）
- `workflow/adapters/__init__.py` - 模块导出

**命名规范：**
- 类名：PascalCase（BaseAdapter, ToolSpec, RunResult）
- 方法名：snake_case（check_version, validate_inputs）
- 常量：UPPER_SNAKE_CASE（如有）

### References

- [Source: docs/planning-artifacts/architecture.md#ADR-002] - Adapter Interface 定义
- [Source: docs/planning-artifacts/epics.md#Story-1.5] - Story 需求和验收标准
- [Source: workflow/lib/errors.py] - ErrorCode enum 定义
- [Source: workflow/lib/audit.py] - RunMetadata dataclass 参考模式

## Dev Agent Record

### Agent Model Used

Claude Opus 4.5 (claude-opus-4-5-20251101)

### Debug Log References

None

### Completion Notes List

- Implemented complete BaseAdapter framework in `workflow/adapters/base.py`
- ToolSpec: dataclass with name, min_version, max_version, conda_env, description; includes check_version_compatible() and to_dict()
- RunResult: dataclass with outputs, summary, metrics, warnings; includes to_dict() and to_json()
- AdapterContext: dataclass with inputs, outputs, config, wildcards, threads, logger, temp_dir; includes helper methods get_input(), get_output(), get_config(), get_wildcard()
- BaseAdapter: ABC with 8 abstract methods (spec, check_version, validate_inputs, build_command, expected_outputs, parse_outputs, timeout_seconds, classify_error)
- classify_common_errors(): helper function for common error pattern classification
- MockAdapter: complete test implementation in test_base.py
- All 35+ unit tests pass
- Module exports updated in __init__.py

### Code Review Record

**Review Date:** 2026-01-21
**Reviewer Model:** Claude Opus 4.5 (claude-opus-4-5-20251101)

**Issues Found:** 0 Critical, 5 Medium, 3 Low
**Issues Fixed:** 5 Medium

**Fixes Applied:**
1. `base.py:14` - Removed unused import `field` from dataclasses
2. `base.py:470` - Fixed OOM detection operator precedence bug (added parentheses)
3. `base.py:219` - Added documentation for temp_dir as ADR-002 extension
4. `test_base.py:21,434` - Moved CompGeneError import to file top level
5. `base.py:99-115` - Fixed ToolSpec.to_dict() serialization asymmetry (conda_env now excluded when empty)

**Tests Added:**
- `test_to_dict_excludes_empty_conda_env()` - verifies new serialization behavior

### File List

- `compgene/workflow/adapters/base.py` (created, reviewed)
- `compgene/workflow/adapters/test_base.py` (created, reviewed)
- `compgene/workflow/adapters/__init__.py` (modified)
