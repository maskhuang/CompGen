# Story 6A.2: 多物种比较

Status: done

## Story

As a **计算生物学研究员**,
I want **比较多物种间的基因存在/缺失模式**,
so that **我可以识别共有和物种特异的基因，为下游分析提供分类统计**。

## Acceptance Criteria

1. **AC1: Orthogroup 共享分类**
   - Given `results/matrices/presence_absence.tsv`（Story 6A.1 输出）
   - When 执行 `matrices_compare` rule
   - Then 生成 `results/matrices/comparison.tsv`
   - And 每行一个 orthogroup，含分类标签：
     - `all_shared`：所有物种共有
     - `species_specific`：仅一个物种特有
     - `partial`：部分物种共有（>1 且 <全部）
   - And 含存在物种列表、缺失物种列表、物种计数

2. **AC2: 汇总统计**
   - Given comparison.tsv 已生成
   - When 生成汇总
   - Then 生成 `results/matrices/comparison_summary.tsv`
   - And 包含：total_orthogroups、all_shared 计数、partial_shared 计数、species_specific 总计数
   - And 包含每个物种的 species_specific 计数

3. **AC3: 边界条件处理**
   - Given 仅含单物种数据的 PA 矩阵
   - When 分类
   - Then 所有 orthogroups 分类为 `species_specific`
   - Given 所有 orthogroups 在所有物种中均存在
   - When 分类
   - Then 所有 orthogroups 分类为 `all_shared`

4. **AC4: 审计记录**
   - Given 比较分析完成
   - When 收集元数据
   - Then 生成 `results/meta/matrices_compare/run.run.json`
   - And 包含输入 checksum、orthogroup 总数、运行时间

## Tasks / Subtasks

- [x] Task 1: 在 `presence_absence.py` 新增 `read_presence_absence_matrix()` (AC: #1, #3)
  - [x] 1.1 读取 PA 矩阵 TSV，返回 `(matrix_dict, species_list)`
  - [x] 1.2 输入验证（文件不存在、表头格式错误）

- [x] Task 2: 在 `presence_absence.py` 新增 `classify_orthogroup_sharing()` (AC: #1, #3)
  - [x] 2.1 遍历矩阵每行，按物种存在数分类：all_shared / species_specific / partial
  - [x] 2.2 返回分类列表，含 orthogroup_id、category、present/absent species

- [x] Task 3: 在 `presence_absence.py` 新增 `summarize_sharing_counts()` (AC: #2)
  - [x] 3.1 从分类列表汇总计数
  - [x] 3.2 包含每物种的 species_specific 计数

- [x] Task 4: 在 `presence_absence.py` 新增 `write_comparison_tsv()` 和 `write_summary_tsv()` (AC: #1, #2)
  - [x] 4.1 原子写入 comparison.tsv
  - [x] 4.2 原子写入 comparison_summary.tsv

- [x] Task 5: 创建 `workflow/scripts/build_comparison.py` (AC: #1, #2, #4)
  - [x] 5.1 桥接 Snakemake → lib/presence_absence.py
  - [x] 5.2 try/except/finally 审计模式

- [x] Task 6: 在 `matrices.smk` 新增 `matrices_compare` rule (AC: #1, #2)
  - [x] 6.1 输入引用 Story 6A.1 输出
  - [x] 6.2 声明 comparison + summary + run_json 输出

- [x] Task 7: 在 `test_presence_absence.py` 新增全部测试 (AC: #1-4)
  - [x] 7.1 测试 `read_presence_absence_matrix()` 正常读取 + 边界 (5 tests)
  - [x] 7.2 测试 `classify_orthogroup_sharing()` 三种分类 (6 tests)
  - [x] 7.3 测试 `summarize_sharing_counts()` 汇总正确性 (3 tests)
  - [x] 7.4 测试 `write_comparison_tsv()` 和 `write_summary_tsv()` 输出格式 (5 tests)
  - [x] 7.5 端到端集成测试（读矩阵 → 分类 → 汇总 → 写文件 → 验证）(1 test)

- [x] Task 8: 版本号升至 1.1.0，端到端验证 (AC: #1-4)
  - [x] 8.1 更新 `presence_absence.py` 的 `__version__`
  - [x] 8.2 运行 pytest 确认全部通过

## Dev Notes

### 输入数据格式 [Source: 6a-1-presence-absence-matrix.md]

`results/matrices/presence_absence.tsv`（Story 6A.1 输出）：
```tsv
orthogroup_id	species1	species2	species3
OG0000000	1	1	1
OG0000001	1	1	0
OG0000002	0	1	1
```

### 目标输出格式

**comparison.tsv：**
```tsv
orthogroup_id	category	n_species_present	present_species	absent_species	specific_to
OG0000000	all_shared	3	species1,species2,species3
OG0000001	partial	2	species1,species2	species3
OG0000002	partial	2	species2,species3	species1
```

**comparison_summary.tsv：**
```tsv
metric	count
total_orthogroups	3
all_shared	1
partial_shared	2
species_specific_total	0
species_specific_species1	0
species_specific_species2	0
species_specific_species3	0
```

### 核心分类逻辑

```python
n_species = len(species_list)
n_present = sum(1 for sp in species_list if row[sp] == 1)

if n_present == n_species:
    category = "all_shared"
elif n_present == 1:
    category = "species_specific"
else:
    category = "partial"
```

### 复用已有代码 [Source: workflow/lib/presence_absence.py]

`write_matrix_tsv()` 的原子写入模式可复用于新的写入函数。
`read_orthogroups_long_format()` 的输入验证模式可复用。

### Bridge 脚本模式 [Source: workflow/scripts/build_matrices.py]

遵循 try/except CompGeneError / except Exception / finally 三段式。
finally 中始终写审计记录。

### 从前序 Story 学到的经验

1. audit 文件路径需用 `shutil.move()` 修正
2. bridge 脚本必须有 `try/except/finally` 错误处理
3. 纯 Python 脚本不需要 conda 指令
4. `print()` 用于 Snakemake 日志捕获

### Project Structure Notes

**修改文件：**
```
workflow/lib/presence_absence.py          # 新增 4 个函数
workflow/lib/test_presence_absence.py     # 新增测试类
workflow/rules/matrices.smk               # 新增 matrices_compare rule
```

**新建文件：**
```
workflow/scripts/build_comparison.py      # Snakemake 桥接脚本
```

### References

- [Source: docs/planning-artifacts/prd.md#FR21] - 系统可比较多物种间的基因存在/缺失模式
- [Source: docs/planning-artifacts/epics-and-stories.md#Story-6A.2] - 多物种比较需求
- [Source: docs/implementation-artifacts/6a-1-presence-absence-matrix.md] - 前序 Story 实现和 Code Review 经验
- [Source: workflow/lib/presence_absence.py] - 现有矩阵模块（复用写入模式）
- [Source: workflow/scripts/build_matrices.py] - bridge 脚本模式参考

## Dev Agent Record

### Agent Model Used

Claude Opus 4.6 (claude-opus-4-6)

### Debug Log References

- 20 new unit tests: all PASSED
- 49 total tests in test_presence_absence.py: 49 passed (0.34s)
- Code review found 1 logic bug (classification priority) and 1 test bug, both fixed before final pass

### Code Review Record

**Reviewer:** Dev Agent (self-review)

**Findings (3 issues):**

1. **[H1] Logic bug — Single-species classification**: `classify_orthogroup_sharing` checked `n_present == n_total` before `n_present == 1`. With single species (n_total=1, n_present=1), OGs were classified as `all_shared` instead of `species_specific` per AC3. **FIXED**: Reordered to check `species_specific` first.
2. **[L1] Test bug — per_species_order filter**: `startswith("species_specific_")` also matched `species_specific_total`. **FIXED**: Added exclusion for `species_specific_total`.
3. **[L2] Module docstring outdated**: Only mentioned Story 6A.1. **FIXED**: Updated to include Story 6A.2 features.

**Post-fix test results:** 49 passed (was 47 + 2 failures), zero regressions.

### Completion Notes List

- Task 1: `read_presence_absence_matrix()` reads PA TSV, returns (matrix_dict, species_list), validates header and file existence
- Task 2: `classify_orthogroup_sharing()` classifies OGs as all_shared/species_specific/partial, sorted by OG ID
- Task 3: `summarize_sharing_counts()` aggregates counts including per-species species_specific counts
- Task 4: `write_comparison_tsv()` + `write_summary_tsv()` with shared `_atomic_write_tsv()` helper
- Task 5: `build_comparison.py` bridge script with try/except/finally + shutil.move audit pattern
- Task 6: `matrices_compare` rule in matrices.smk with 3 outputs (comparison, summary, run_json)
- Task 7: 20 new tests across 6 test classes + 1 integration test
- Task 8: Version bumped to 1.1.0, all 49 tests pass
- All AC satisfied: AC1 (classification), AC2 (summary), AC3 (edge cases), AC4 (audit)

### File List

**New files:**
- `compgene/workflow/scripts/build_comparison.py` — Snakemake bridge script

**Modified files:**
- `compgene/workflow/lib/presence_absence.py` — Added 6 functions + 1 helper, version 1.0.0 → 1.1.0
- `compgene/workflow/lib/test_presence_absence.py` — Added 20 tests (29 → 49 total)
- `compgene/workflow/rules/matrices.smk` — Added `rule matrices_compare`
