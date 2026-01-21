---
stepsCompleted: [1, 2, 3, 4]
inputDocuments:
  - "prd.md"
  - "architecture.md"
---

# CompGene - Epic Breakdown

## Overview

This document provides the complete epic and story breakdown for CompGene, decomposing the requirements from the PRD, UX Design if it exists, and Architecture requirements into implementable stories.

## Requirements Inventory

### Functional Requirements

**数据准备与 QC 能力 (FR1-FR5b)**

- FR1: 用户可通过 NCBI 物种标识符自动下载基因组注释数据
- FR2: 系统可解析 GFF3 和 GTF 格式的注释文件
- FR3: 系统可提取基因、转录本、CDS 等特征信息并标准化存储
- FR4: 系统可缓存已下载的注释数据，避免重复下载
- FR5: 用户可指定本地注释文件路径替代自动下载
- FR5a: 系统可调用 BUSCO 评估装配/注释完整性
- FR5b: 系统可生成每物种标准化目录结构和 BUSCO 汇总表

**直系同源推断能力 - OrthoFinder (FR6-FR10)**

- FR6: 系统可调用 OrthoFinder 执行 orthogroup 推断
- FR7: 系统可生成 orthogroups 表（跨物种基因家族分组）
- FR8: 系统可输出物种树和基因树
- FR9: 系统可统计每物种的 orthogroup copy number
- FR10: 系统可识别物种特异性 orthogroups

**功能注释能力 - eggNOG-mapper (FR11-FR15)**

- FR11: 系统可调用 eggNOG-mapper 对基因/蛋白进行功能注释
- FR12: 用户可指定 eggNOG 参考数据库（如 eukaryota、vertebrates）
- FR13: 系统可提取 GO/KEGG/COG 功能标签
- FR14: 系统可将功能注释汇总到 orthogroup 层面
- ~~FR15: 系统可生成功能富集统计报告~~ → **Post-MVP**

**缺失核验能力 - Liftoff (FR16-FR19)**

- FR16: 系统可调用 Liftoff 将参考注释映射到目标 assembly
- FR17: 系统可识别"真缺失"vs"注释假缺失"候选
- FR18: 系统可输出核验后的"真缺失候选清单"
- FR19: 系统可生成注释补全后的对照注释文件

**基因存在/缺失分析能力 (FR20-FR25)**

- FR20: 系统可从 orthogroups 派生 presence/absence 矩阵（0/1 或 copy number）
- FR21: 系统可比较多物种间的基因存在/缺失模式
- FR22: 系统可按基因家族分类统计
- FR23: 系统可识别物种特异性 orthogroups
- FR24: 系统可整合 Liftoff 核验结果标记可信度
- FR25: 系统可生成比较分析汇总报告

**RNA-Seq 表达分析能力 (FR26-FR31)**

- FR26: 用户可提供 counts 矩阵和样本元数据作为输入
- FR27: 系统可调用 DESeq2 执行差异表达分析
- FR28: 系统可输出标准化后的表达矩阵
- FR29: 系统可输出差异表达结果（DESeq2 标准格式）
- FR30: 系统可将表达数据与 orthogroup 关联
- FR31: 系统可进行跨物种表达模式比较

**工作流编排能力 (FR32-FR37)**

- FR32: 用户可通过 YAML 配置文件定义完整分析流程
- FR33: 系统可自动解析模块依赖并按顺序执行
- FR34: 系统可缓存中间结果，支持断点续跑
- FR35: 用户可通过 --resume 标志从中断点继续
- FR36: 用户可通过 --incremental 标志仅处理新增数据
- FR37: 系统可支持 dry-run 模式预览执行计划

**配置管理能力 (FR38-FR41)**

- FR38: 用户可在配置文件中定义多个物种及其参数
- FR39: 系统可验证配置文件格式和必需字段
- FR40: 用户可通过命令行参数覆盖配置文件设置
- FR41: 系统可指定输出目录和资源限制（线程数、内存）

**结果输出与报告能力 (FR42-FR46)**

- FR42: 系统可输出所有数据结果为 TSV 格式
- FR43: 系统可生成 Markdown 格式汇总报告
- FR44: 系统可生成 HTML 格式可视化报告
- FR45: 系统可输出机器可解析的 JSON 日志
- FR46: 系统可输出人类可读的文本日志

**错误处理与恢复能力 (FR47-FR49)**

- FR47: 系统可在出错时提供清晰的错误信息和修复建议
- FR48: 系统可对网络请求实现指数退避重试
- FR49: 系统可返回标准化退出码（成功、配置错误、数据错误、网络错误、资源不足）

### NonFunctional Requirements

**性能要求 (NFR1-NFR4)**

- NFR1: 单物种完整分析流程应在 30 分钟内完成（标准服务器配置）
- NFR2: 大文件（>100MB）处理时内存占用应可控，支持分块读取
- NFR3: 多核并行处理时 CPU 利用率应达到 70% 以上
- NFR4: 中间结果缓存命中时，重复步骤跳过延迟 < 1 秒

**可靠性要求 (NFR5-NFR8)**

- NFR5: 分析结果 100% 可复现（相同输入 + 配置 = 相同输出）
- NFR6: 所有依赖工具版本和数据库版本必须记录在日志中
- NFR7: 支持从任意断点恢复，中间状态完整保存
- NFR8: 网络中断后自动重试，最多 3 次，指数退避间隔

**集成要求 (NFR9-NFR15)**

- NFR9: NCBI API 调用应遵循官方速率限制（每秒 3 请求）
- NFR10: OrthoFinder 调用应支持版本 2.5.x
- NFR11: eggNOG-mapper 调用应支持版本 2.1.x
- NFR12: Liftoff 调用应支持版本 1.6.x
- NFR13: BUSCO 调用应支持版本 5.x/6.x
- NFR14: DESeq2 输出格式应兼容 Bioconductor 标准
- NFR15: 所有外部工具调用应有超时保护（默认 30 分钟）

**可维护性要求 (NFR16-NFR18)**

- NFR16: 模块间低耦合，单个模块可独立测试
- NFR17: 日志级别可配置（DEBUG/INFO/WARNING/ERROR）
- NFR18: 配置文件变更无需修改代码

### Additional Requirements

**来自 Architecture 的技术需求：**

**Starter Template（关键 - 影响 Epic 1 Story 1）：**
- 使用 Snakemake 9.14.8 作为工作流引擎
- Python 3.11+ 运行时
- 项目结构已在 Architecture 中完整定义

**数据架构：**
- 文件产物作为 Source of Truth（Snakemake 依赖跟踪）
- SQLite 作为 Index of Record（ID 映射、审计、快速查询）
- 单写入 rule 模式避免 SQLite 并发问题
- 三层 ID 策略：raw + norm + uid

**工具适配层：**
- BaseAdapter 抽象类定义 8 个方法
- 每个外部工具一个 Adapter 模块
- 统一 Runner (run_adapter.py) 处理执行流程

**错误处理：**
- ErrorCode enum 定义 9 种错误码
- ERROR_RECOVERY dict 提供恢复建议
- 混合重试策略（Snakemake 层 + Runner 层）

**日志审计：**
- 双格式日志：.log（人类可读）+ .jsonl（机器可查）
- 规则级审计：每个 rule 生成 .run.json
- SQLite 审计表支持追溯

**命名规范：**
- Snakemake rule: `{module}_{action}`
- Python: snake_case（函数/变量）、PascalCase（类）
- 输出路径: `results/{category}/{wildcards}/`
- 所有大文件使用原子写入（先 .tmp 再 rename）

**conda 环境隔离：**
- 每个外部工具独立 conda 环境
- envs/*.yaml 版本锁定

### FR Coverage Map

| FR | Epic | 说明 |
|----|------|------|
| FR1 | Epic 2 | NCBI 注释下载 |
| FR2 | Epic 2 | GFF/GTF 解析 |
| FR3 | Epic 2 | 特征提取与标准化 |
| FR4 | Epic 2 | 缓存机制 |
| FR5 | Epic 2 | 本地文件支持 |
| FR5a | Epic 2 | BUSCO QC |
| FR5b | Epic 2 | 标准化目录 |
| FR6 | Epic 3 | OrthoFinder 调用 |
| FR7 | Epic 3 | Orthogroups 表 |
| FR8 | Epic 3 | 物种树/基因树 |
| FR9 | Epic 3 | Copy number 统计 |
| FR10 | Epic 3 | 物种特异 OG |
| FR11 | Epic 4 | eggNOG-mapper 调用 |
| FR12 | Epic 4 | 数据库选择 |
| FR13 | Epic 4 | GO/KEGG/COG 提取 |
| FR14 | Epic 4 | OG 层面汇总 |
| ~~FR15~~ | Post-MVP | 功能富集统计 |
| FR16 | Epic 5 | Liftoff 映射 |
| FR17 | Epic 5 | 真缺失识别 |
| FR18 | Epic 5 | 候选清单输出 |
| FR19 | Epic 5 | 补全注释文件 |
| FR20 | Epic 6A | Presence/absence 矩阵 |
| FR21 | Epic 6A | 多物种比较 |
| FR22 | Epic 6A | 基因家族统计 |
| FR23 | Epic 6B | 物种特异性识别 |
| FR24 | Epic 6B | Liftoff 可信度整合 |
| FR25 | Epic 6B | 比较汇总报告 |
| FR26 | Epic 7 | Counts 矩阵输入 |
| FR27 | Epic 7 | DESeq2 调用 |
| FR28 | Epic 7 | 标准化矩阵输出 |
| FR29 | Epic 7 | DE 结果输出 |
| FR30 | Epic 7 | OG 关联 |
| FR31 | Epic 7 | 跨物种比较 |
| FR32 | Epic 1 | YAML 配置 |
| FR33 | Epic 1 | 依赖解析 |
| FR34 | Epic 1 | 缓存与断点 |
| FR35 | Epic 1 | --resume 标志 |
| FR36 | Epic 1 | --incremental 标志 |
| FR37 | Epic 1 | dry-run 模式 |
| FR38 | Epic 1 | 多物种定义 |
| FR39 | Epic 1 | 配置验证 |
| FR40 | Epic 1 | 命令行覆盖 |
| FR41 | Epic 1 | 资源限制 |
| FR42 | Epic 8 | TSV 输出 |
| FR43 | Epic 8 | Markdown 报告 |
| FR44 | Epic 8 | HTML 报告 |
| FR45 | Epic 1 | JSON 日志 |
| FR46 | Epic 1 | 文本日志 |
| FR47 | Epic 1 | 错误信息 |
| FR48 | Epic 1 | 重试机制 |
| FR49 | Epic 1 | 退出码 |

## Epic List

### Epic 1: 项目基础设施、配置与可观测性

用户可以初始化 CompGene 项目、验证配置，所有后续 Epic 自动获得日志/审计能力。

**FRs 覆盖：** FR32-41, FR45-46, FR47-49

**实现内容：**
- Snakemake 骨架（Snakefile + rules/common.smk）
- 配置 schema（schemas/config.schema.yaml）
- BaseAdapter 框架（adapters/base.py）
- ErrorCode 系统（lib/errors.py）
- 审计基座（lib/audit.py + .run.json 生成）
- 双格式日志（.log + .jsonl）
- 统一 Runner（scripts/run_adapter.py）

---

### Epic 2: 数据标准化与质控

用户可以准备输入数据、验证格式，并评估注释完整性。

**FRs 覆盖：** FR1-5b

**实现内容：**
- rules/standardize.smk
- rules/qc.smk
- adapters/busco.py
- lib/gff.py, lib/fasta.py
- envs/busco.yaml

---

### Epic 3: 直系同源分析

用户可以识别跨物种的直系同源基因，获得 orthogroups 和系统发育树。

**FRs 覆盖：** FR6-10

**实现内容：**
- rules/orthology.smk
- adapters/orthofinder.py
- envs/orthofinder.yaml

---

### Epic 4: 功能注释

用户可以为基因添加 GO/KEGG/COG 功能标签，汇总到 orthogroup 层面。

**FRs 覆盖：** FR11-14

**实现内容：**
- rules/annotation.smk
- adapters/eggnog.py
- envs/eggnog.yaml

---

### Epic 5: 缺失核验

用户可以区分"真缺失"与"注释假缺失"，获得可信的缺失候选清单。

**FRs 覆盖：** FR16-19

**实现内容：**
- rules/validation.smk
- adapters/liftoff.py
- envs/liftoff.yaml

---

### Epic 6A: 存在/缺失矩阵（基础）

用户可以生成 presence/absence 矩阵，进行多物种比较和基因家族统计。

**FRs 覆盖：** FR20-22

**实现内容：**
- rules/matrices.smk（基础部分）
- lib/presence_absence.py

---

### Epic 6B: 存在/缺失分析（增强）

用户可以识别物种特异基因，整合 Liftoff 可信度标记，生成汇总报告。

**FRs 覆盖：** FR23-25

**依赖：** Epic 5 + Epic 6A

**实现内容：**
- rules/matrices.smk（增强部分）
- 可信度标记逻辑
- 比较汇总报告

---

### Epic 7: 表达分析

用户可以分析差异表达，并关联到 orthogroup 进行跨物种比较。

**FRs 覆盖：** FR26-31

**实现内容：**
- rules/expression.smk
- adapters/deseq2.py
- envs/deseq2.yaml

---

### Epic 8: 报告生成

用户可以生成综合报告，查询 SQLite 元数据进行追溯。

**FRs 覆盖：** FR42-44

**实现内容：**
- rules/reporting.smk
- rules/ingest.smk（SQLite 汇总写入）
- lib/report.py
- resources/report_template.html

---

---

## Epic 1: 项目基础设施、配置与可观测性

用户可以初始化 CompGene 项目、验证配置，所有后续 Epic 自动获得日志/审计能力。

### Story 1.1: 项目骨架初始化

As a **计算生物学研究员**,
I want **通过标准命令初始化 CompGene 项目结构**,
So that **我可以在一致的目录布局中开始配置和运行分析**。

**Acceptance Criteria:**

**Given** 用户在空目录中
**When** 执行项目初始化
**Then** 生成完整目录结构（workflow/, config/, schemas/, profiles/, results/, logs/）
**And** Snakefile 包含 `min_version("9.0")` 和 `validate()` 调用
**And** pyproject.toml 定义项目元数据和依赖

---

### Story 1.2: 配置 Schema 与验证

As a **计算生物学研究员**,
I want **系统在运行前验证我的配置文件格式**,
So that **我能在分析开始前发现配置错误**。

**Acceptance Criteria:**

**Given** 用户编写了 config/config.yaml
**When** 执行 `snakemake --dry-run`
**Then** 系统使用 schemas/config.schema.yaml 验证配置
**And** 缺少必需字段时返回清晰错误信息

**Given** 用户通过命令行传递 `--config key=value`
**When** 执行分析
**Then** 命令行参数覆盖配置文件中的对应值

---

### Story 1.3: 错误码系统

As a **计算生物学研究员**,
I want **系统返回标准化的错误码和修复建议**,
So that **我能快速定位问题并知道如何解决**。

**Acceptance Criteria:**

**Given** 分析过程中发生错误
**When** 系统捕获到错误
**Then** 返回 ErrorCode enum 中定义的错误码
**And** 输出对应的恢复建议

**Given** 分析成功完成
**When** 进程退出
**Then** 退出码为 0

**Given** 发生配置错误
**When** 进程退出
**Then** 退出码为 2

---

### Story 1.4: 日志与审计基座

As a **计算生物学研究员**,
I want **分析过程自动生成结构化日志和审计记录**,
So that **我可以追溯任何输出的生成过程**。

**Acceptance Criteria:**

**Given** 任意 rule 执行完成
**When** Runner 收集元数据
**Then** 生成双格式日志：`.log`（文本）+ `.jsonl`（JSON Lines）
**And** 生成 `.run.json` 包含命令、版本、checksums、时间戳

---

### Story 1.5: BaseAdapter 框架

As a **开发者**,
I want **统一的 Adapter 抽象接口**,
So that **后续工具集成遵循一致模式**。

**Acceptance Criteria:**

**Given** 需要集成新的外部工具
**When** 继承 BaseAdapter
**Then** 必须实现 8 个方法：check_version, validate_inputs, build_command, expected_outputs, parse_outputs, timeout_seconds, classify_error, spec

**Given** Adapter 实现完成
**When** 运行单元测试
**Then** 可以独立测试（mock 外部工具）

---

### Story 1.6: 统一 Runner 与重试

As a **计算生物学研究员**,
I want **外部工具调用有统一的执行、超时和重试机制**,
So that **临时故障不会导致整个分析失败**。

**Acceptance Criteria:**

**Given** rule 调用外部工具
**When** 使用 run_adapter.py 执行
**Then** 按顺序：版本检测 → 输入校验 → 命令执行 → 输出解析 → 审计记录

**Given** 工具执行超时
**When** 超过 timeout_seconds
**Then** 发送 SIGTERM，等待后 SIGKILL，返回 E_TIMEOUT

**Given** 网络错误（可重试）
**When** 发生错误
**Then** 自动重试最多 3 次，指数退避间隔

---

### Story 1.7: 断点续跑与 dry-run

As a **计算生物学研究员**,
I want **分析中断后可以从断点继续，并能预览执行计划**,
So that **我不需要重跑已完成的步骤**。

**Acceptance Criteria:**

**Given** 分析中途中断
**When** 执行 `snakemake --rerun-incomplete`
**Then** 仅重跑未完成的 rules

**Given** 用户想预览执行计划
**When** 执行 `snakemake --dry-run`
**Then** 显示将要执行的 rules，不实际执行

**Given** 所有输入文件未变化
**When** 重新执行
**Then** 所有 rules 被跳过，延迟 < 1 秒

---

## Epic 2: 数据标准化与质控

用户可以准备输入数据、验证格式，并评估注释完整性。

### Story 2.1: GFF/FASTA 解析库

As a **计算生物学研究员**,
I want **系统能解析 GFF3/GTF 和 FASTA 格式文件**,
So that **我可以使用不同来源的注释数据**。

**Acceptance Criteria:**

**Given** 用户提供 GFF3 或 GTF 格式注释文件
**When** 系统解析文件
**Then** 正确提取 gene、transcript、CDS 等特征
**And** 建立 gene → transcript → protein 的层级关系

**Given** 文件格式错误
**When** 解析失败
**Then** 返回 E_INPUT_FORMAT 错误码和具体错误位置

---

### Story 2.2: 本地数据标准化

As a **计算生物学研究员**,
I want **将本地注释文件标准化到统一目录结构**,
So that **后续分析可以使用一致的输入格式**。

**Acceptance Criteria:**

**Given** 用户在配置文件中指定本地注释路径
**When** 执行 standardize rule
**Then** 生成标准化目录：genome.fa.gz, annotation.gff3.gz, proteins.longest.fa.gz

**Given** 基因有多个转录本
**When** 提取蛋白序列
**Then** 选择最长转录本作为代表，记录 is_representative 标记

---

### Story 2.3: NCBI 注释下载

As a **计算生物学研究员**,
I want **通过物种标识符自动从 NCBI 下载注释**,
So that **我不需要手动查找和下载数据**。

**Acceptance Criteria:**

**Given** 配置文件中 annotation 为 NCBI 标识符
**When** 执行 standardize rule
**Then** 自动从 NCBI 下载基因组和注释文件

**Given** 已下载过相同标识符的数据
**When** 再次请求
**Then** 使用本地缓存，不重复下载

**Given** NCBI API 调用频繁
**When** 超过速率限制
**Then** 自动等待并重试，遵循 3 req/s 限制

---

### Story 2.4: BUSCO 质控

As a **计算生物学研究员**,
I want **评估每个物种的注释完整性**,
So that **我能了解数据质量并识别潜在问题**。

**Acceptance Criteria:**

**Given** 标准化后的蛋白序列
**When** 执行 qc_busco rule
**Then** 调用 BUSCO 评估注释完整性，输出到 results/qc/{species}/busco/

**Given** 多个物种完成 BUSCO
**When** 生成汇总
**Then** 创建 busco_summary.tsv 包含所有物种的统计

---

## Epic 3: 直系同源分析

用户可以识别跨物种的直系同源基因，获得 orthogroups 和系统发育树。

### Story 3.1: OrthoFinder Adapter

As a **计算生物学研究员**,
I want **系统能调用 OrthoFinder 进行直系同源推断**,
So that **我可以识别跨物种的基因家族关系**。

**Acceptance Criteria:**

**Given** 多个物种的标准化蛋白序列
**When** 执行 orthology_infer rule
**Then** 调用 OrthoFinder 执行 orthogroup 推断
**And** 输出到 results/orthology/{species_set}/

**Given** OrthoFinder 版本不兼容
**When** 检测版本
**Then** 返回 E_TOOL_VERSION 错误码

**Given** OrthoFinder 执行超时
**When** 超过配置的超时时间
**Then** 返回 E_TIMEOUT 并记录已完成进度

---

### Story 3.2: Orthogroups 表生成

As a **计算生物学研究员**,
I want **获得结构化的 orthogroups 表**,
So that **我可以查看跨物种的基因家族分组**。

**Acceptance Criteria:**

**Given** OrthoFinder 运行完成
**When** 解析输出
**Then** 生成 orthogroups.tsv 包含：orthogroup_id, gene_id, species_id

**Given** orthogroups 表生成完成
**When** 查询特定基因
**Then** 可以找到其所属的 orthogroup 和同组的其他物种基因

---

### Story 3.3: 物种树与基因树输出

As a **计算生物学研究员**,
I want **获得物种树和基因树**,
So that **我可以了解物种间的进化关系**。

**Acceptance Criteria:**

**Given** OrthoFinder 运行完成
**When** 提取树文件
**Then** 输出 species_tree.nwk（Newick 格式）
**And** 输出各 orthogroup 的基因树到 gene_trees/ 目录

---

### Story 3.4: Copy Number 统计

As a **计算生物学研究员**,
I want **统计每物种在各 orthogroup 中的基因拷贝数**,
So that **我可以识别基因家族扩张或收缩**。

**Acceptance Criteria:**

**Given** orthogroups 表生成完成
**When** 统计 copy number
**Then** 生成 copy_number.tsv：每行一个 orthogroup，每列一个物种，值为基因数量

**Given** copy number 矩阵
**When** 分析模式
**Then** 可以识别物种特异性 orthogroups（仅在单一物种中存在）

---

## Epic 4: 功能注释

用户可以为基因添加 GO/KEGG/COG 功能标签，汇总到 orthogroup 层面。

### Story 4.1: eggNOG-mapper Adapter

As a **计算生物学研究员**,
I want **系统能调用 eggNOG-mapper 进行功能注释**,
So that **我可以了解基因的功能分类**。

**Acceptance Criteria:**

**Given** 标准化的蛋白序列
**When** 执行 annotation_eggnog rule
**Then** 调用 eggNOG-mapper 进行功能注释
**And** 输出到 results/annotation/eggnog/{species}/

**Given** 用户指定参考数据库（如 eukaryota）
**When** 执行注释
**Then** 使用指定数据库进行比对

---

### Story 4.2: GO/KEGG/COG 提取

As a **计算生物学研究员**,
I want **提取结构化的功能标签**,
So that **我可以按功能分类分析基因**。

**Acceptance Criteria:**

**Given** eggNOG-mapper 输出
**When** 解析注释结果
**Then** 提取并结构化：GO terms, KEGG pathways, COG categories
**And** 生成 annotations.tsv：gene_id, go_terms, kegg_pathways, cog_category

---

### Story 4.3: Orthogroup 层面功能汇总

As a **计算生物学研究员**,
I want **将功能注释汇总到 orthogroup 层面**,
So that **我可以了解基因家族的整体功能**。

**Acceptance Criteria:**

**Given** 各物种的功能注释和 orthogroups 表
**When** 汇总到 orthogroup 层面
**Then** 生成 og_annotations.tsv：orthogroup_id, consensus_go, consensus_kegg, consensus_cog

**Given** orthogroup 内基因功能不一致
**When** 汇总时
**Then** 采用多数投票或全部保留策略（可配置）

---

## Epic 5: 缺失核验

用户可以区分"真缺失"与"注释假缺失"，获得可信的缺失候选清单。

### Story 5.1: Liftoff Adapter

As a **计算生物学研究员**,
I want **系统能调用 Liftoff 进行注释映射**,
So that **我可以验证基因缺失是真实的还是注释遗漏**。

**Acceptance Criteria:**

**Given** 参考物种注释和目标物种基因组
**When** 执行 validation_liftoff rule
**Then** 调用 Liftoff 将参考注释映射到目标 assembly
**And** 输出到 results/validation/liftoff/{source}_to_{target}/

---

### Story 5.2: 真缺失 vs 假缺失识别

As a **计算生物学研究员**,
I want **区分真正的基因缺失和注释遗漏**,
So that **我的比较分析结果更可靠**。

**Acceptance Criteria:**

**Given** Liftoff 映射结果
**When** 分析映射状态
**Then** 将基因分类为：
- 成功映射（注释存在）
- 部分映射（可能的假缺失）
- 无法映射（真缺失候选）

**Given** 分类完成
**When** 生成候选清单
**Then** 输出 true_absence_candidates.tsv：gene_id, source_species, target_species, confidence

---

### Story 5.3: 补全注释文件

As a **计算生物学研究员**,
I want **生成包含 Liftoff 补全的注释文件**,
So that **我可以使用更完整的注释进行下游分析**。

**Acceptance Criteria:**

**Given** Liftoff 成功映射的基因
**When** 生成补全注释
**Then** 输出 lifted_annotation.gff3：包含原始注释 + Liftoff 补全的基因

**Given** 补全注释与原始注释有冲突
**When** 合并时
**Then** 优先保留原始注释，Liftoff 结果作为补充

---

## Epic 6A: 存在/缺失矩阵（基础）

用户可以生成 presence/absence 矩阵，进行多物种比较和基因家族统计。

### Story 6A.1: Presence/Absence 矩阵生成

As a **计算生物学研究员**,
I want **从 orthogroups 生成 presence/absence 矩阵**,
So that **我可以直观看到各物种的基因分布**。

**Acceptance Criteria:**

**Given** orthogroups 表和 copy_number 表
**When** 执行 matrices_generate rule
**Then** 生成 presence_absence.tsv：
- 行：orthogroup_id
- 列：各物种
- 值：0/1（存在/缺失）或 copy number

---

### Story 6A.2: 多物种比较

As a **计算生物学研究员**,
I want **比较多物种间的基因存在/缺失模式**,
So that **我可以识别共有和物种特异的基因**。

**Acceptance Criteria:**

**Given** presence/absence 矩阵
**When** 分析模式
**Then** 统计：
- 所有物种共有的 orthogroups
- 各物种特异的 orthogroups
- 部分物种共有的 orthogroups

---

### Story 6A.3: 基因家族统计

As a **计算生物学研究员**,
I want **按基因家族分类统计**,
So that **我可以了解不同功能类别的分布**。

**Acceptance Criteria:**

**Given** presence/absence 矩阵和功能注释
**When** 按 COG/GO 分类统计
**Then** 生成 family_stats.tsv：functional_category, total_ogs, shared_ogs, species_specific_ogs

---

## Epic 6B: 存在/缺失分析（增强）

用户可以识别物种特异基因，整合 Liftoff 可信度标记，生成汇总报告。

### Story 6B.1: 物种特异性 Orthogroups 识别

As a **计算生物学研究员**,
I want **识别物种特异性的 orthogroups**,
So that **我可以研究物种独特的基因**。

**Acceptance Criteria:**

**Given** presence/absence 矩阵
**When** 筛选物种特异 OG
**Then** 输出 species_specific_ogs.tsv：orthogroup_id, species, gene_count, gene_ids

---

### Story 6B.2: Liftoff 可信度整合

As a **计算生物学研究员**,
I want **将 Liftoff 核验结果整合到比较矩阵中**,
So that **我可以知道哪些缺失是可信的**。

**Acceptance Criteria:**

**Given** presence/absence 矩阵和 Liftoff 真缺失候选
**When** 整合可信度
**Then** 更新矩阵添加 confidence 列：high（真缺失）, low（可能假缺失）, na（未验证）

---

### Story 6B.3: 比较分析汇总报告

As a **计算生物学研究员**,
I want **生成比较分析的汇总报告**,
So that **我可以快速了解分析结果**。

**Acceptance Criteria:**

**Given** 所有比较分析完成
**When** 生成汇总
**Then** 输出 comparison_summary.md 包含：
- 物种数量和名称
- Orthogroup 总数和分布
- 物种特异基因统计
- 可信度分布统计

---

## Epic 7: 表达分析

用户可以分析差异表达，并关联到 orthogroup 进行跨物种比较。

### Story 7.1: DESeq2 Adapter

As a **计算生物学研究员**,
I want **系统能调用 DESeq2 进行差异表达分析**,
So that **我可以识别差异表达的基因**。

**Acceptance Criteria:**

**Given** counts 矩阵和样本元数据
**When** 执行 expression_deseq2 rule
**Then** 调用 DESeq2（通过 Rscript）执行差异表达分析
**And** 输出到 results/expression/

---

### Story 7.2: 标准化矩阵输出

As a **计算生物学研究员**,
I want **获得标准化后的表达矩阵**,
So that **我可以进行跨样本比较**。

**Acceptance Criteria:**

**Given** DESeq2 分析完成
**When** 提取标准化结果
**Then** 输出 normalized_counts.tsv：gene_id × sample_id 矩阵

---

### Story 7.3: 差异表达结果输出

As a **计算生物学研究员**,
I want **获得 DESeq2 标准格式的差异表达结果**,
So that **我可以使用常规工具进行下游分析**。

**Acceptance Criteria:**

**Given** DESeq2 分析完成
**When** 提取 DE 结果
**Then** 输出 DE_results.tsv 包含：gene_id, baseMean, log2FoldChange, pvalue, padj

---

### Story 7.4: Orthogroup 表达关联

As a **计算生物学研究员**,
I want **将表达数据与 orthogroup 关联**,
So that **我可以在基因家族层面分析表达**。

**Acceptance Criteria:**

**Given** DE 结果和 orthogroups 表
**When** 关联数据
**Then** 输出 og_expression.tsv：orthogroup_id, species, gene_id, log2FC, padj

---

### Story 7.5: 跨物种表达比较

As a **计算生物学研究员**,
I want **进行跨物种的表达模式比较**,
So that **我可以识别保守和分化的表达模式**。

**Acceptance Criteria:**

**Given** 多物种的 og_expression 数据
**When** 比较表达模式
**Then** 识别：
- 表达保守的 orthogroups（各物种表达方向一致）
- 表达分化的 orthogroups（物种间表达方向相反）

---

## Epic 8: 报告生成

用户可以生成综合报告，查询 SQLite 元数据进行追溯。

### Story 8.1: SQLite 审计汇总

As a **计算生物学研究员**,
I want **所有审计记录汇总到 SQLite 数据库**,
So that **我可以快速查询和追溯分析历史**。

**Acceptance Criteria:**

**Given** 所有 rule 的 .run.json 文件
**When** 执行 ingest_db rule
**Then** 将审计记录写入 results/meta/compgene.db

**Given** SQLite 数据库
**When** 查询特定输出文件
**Then** 可以找到生成该文件的命令、参数、工具版本、输入 checksums

---

### Story 8.2: TSV 数据输出

As a **计算生物学研究员**,
I want **所有分析结果以 TSV 格式输出**,
So that **我可以用 Excel、R、Python 等工具处理**。

**Acceptance Criteria:**

**Given** 分析完成
**When** 检查输出
**Then** 所有数据结果为 TSV 格式（UTF-8，Tab 分隔）
**And** 列名使用 snake_case

---

### Story 8.3: Markdown 报告生成

As a **计算生物学研究员**,
I want **生成 Markdown 格式的汇总报告**,
So that **我可以在版本控制中追踪报告变化**。

**Acceptance Criteria:**

**Given** 所有分析完成
**When** 执行 reporting_summary rule
**Then** 生成 results/reports/summary.md 包含：
- 分析概览（物种、时间、参数）
- 数据质量摘要（BUSCO）
- 主要发现（orthogroup 统计、物种特异基因）
- 输出文件清单

---

### Story 8.4: HTML 报告生成

As a **计算生物学研究员**,
I want **生成可视化的 HTML 报告**,
So that **我可以直观地浏览分析结果**。

**Acceptance Criteria:**

**Given** 所有分析完成
**When** 执行 reporting_html rule
**Then** 生成 results/reports/summary.html 包含：
- 与 Markdown 报告相同的内容
- 交互式表格（可排序、筛选）
- 基础图表（BUSCO 柱状图、orthogroup 分布饼图）

---

## Post-MVP Features

- FR15: 功能富集统计报告
