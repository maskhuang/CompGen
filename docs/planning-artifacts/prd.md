---
stepsCompleted:
  - step-01-init
  - step-02-discovery
  - step-03-success
  - step-04-journeys
  - step-05-domain
  - step-06-innovation
  - step-07-project-type
  - step-08-scoping
  - step-09-functional
  - step-10-nonfunctional
  - step-11-polish
  - step-12-complete
inputDocuments:
  - "product-brief-CompGene-2026-01-20.md"
  - "Yucheng_project (1).docx"
workflowType: 'prd'
date: 2026-01-20
author: Yuchenghuang
documentCounts:
  briefs: 1
  research: 0
  brainstorming: 0
  projectDocs: 1
classification:
  projectType: cli_tool
  domain: scientific
  complexity: medium
  projectContext: greenfield
---

# Product Requirements Document - CompGene

**Author:** Yuchenghuang
**Date:** 2026-01-20

## Executive Summary

CompGene 是比较基因组学自动化分析管道，整合多种工具实现从数据准备到结果报告的一条命令运行。

**核心工具组合**：
- **OrthoFinder** - 构建 orthogroups（比较骨架）
- **eggNOG-mapper** - 功能注释（GO/KEGG/COG）
- **Liftoff** - 缺失核验（真缺失 vs 注释假缺失）
- **BUSCO** - 装配/注释完整性 QC
- **DESeq2** - RNA-Seq 差异表达分析

**目标物种**：Microcebus murinus、Lemur catta、猕猴（MFA-T2T）
**核心价值**：将数小时的手动分析缩短至 30 分钟，100% 可复现
**技术栈**：Python CLI + YAML 配置

## Success Criteria

### User Success

**效率提升**：
- 单物种分析时间：从数小时缩短至 30 分钟内
- 多物种批量处理：支持一条命令运行全部分析
- 新物种扩展：仅需添加输入文件，无需重新配置（< 5 分钟）

**质量保证**：
- 结果 100% 可复现：配置文件 + 完整日志记录
- 中间结果缓存：支持断点续跑，避免重复计算

**可用性**：
- 有 Python/R 经验的用户 15 分钟内可完成首次分析
- 清晰的错误提示和文档

### Business Success

**研究里程碑**：
1. **初始阶段**：完成 Microcebus murinus、Lemur catta 与猕猴的比较分析
2. **扩展阶段**：整合 Kinnex 长读长数据，支持 4 个额外狐猴物种
3. **进阶阶段**：扩展至猿类物种和 CHM13 人类参考组装

**产出目标**：
- 可复用的自动化分析管道
- 标准化的分析报告格式
- 为团队后续研究提供基础设施

### Technical Success

**性能要求**：
- 核心流程自动化率：100%
- 分析结果可复现率：100%

**可靠性要求**：
- 支持断点续跑，中断后可从上次位置继续
- 完整的错误处理和日志记录

**可维护性**：
- 模块化设计，便于添加新分析方法
- YAML 配置文件驱动，易于调整参数

### Measurable Outcomes

| KPI | 目标值 | 验证方式 |
|-----|--------|----------|
| 核心流程自动化率 | 100% | 端到端测试 |
| 分析结果可复现率 | 100% | 重复运行对比 |
| 新物种扩展配置时间 | < 5 分钟 | 实际操作计时 |
| 首次使用学习时间 | < 15 分钟 | 用户测试 |
| 单物种分析时间 | < 30 分钟 | 性能测试 |

## Product Scope

### MVP - Minimum Viable Product

**核心模块**：
1. **数据准备与 QC 模块** - NCBI 注释获取 + BUSCO 完整性检测
2. **直系同源推断模块** - OrthoFinder 构建 orthogroups（比较骨架）
3. **功能注释模块** - eggNOG-mapper 补充 GO/KEGG/COG 标签
4. **缺失核验模块** - Liftoff 区分真缺失 vs 注释假缺失
5. **基因存在/缺失分析** - 多物种比较矩阵和报告
6. **RNA-Seq 表达分析整合** - DESeq2 流程和跨物种比较
7. **工作流编排引擎** - YAML 配置、日志、缓存、断点续跑

**MVP 验收标准**：
- [ ] 成功完成三物种的基因存在/缺失分析
- [ ] 直系同源推断结果与手动分析一致
- [ ] RNA-Seq 分析流程可复现
- [ ] 全流程可通过配置文件一键运行

### Growth Features (Post-MVP)

- Kinnex 长读长测序数据支持（二月中旬数据到达后）
- 4 个额外狐猴物种支持
- 表达-直系同源综合分析增强

### Vision (Future)

- 扩展至猿类物种和 CHM13 人类参考组装
- 完整灵长类比较基因组学框架
- 可复用的比较基因组学分析平台
- 支持更广泛的物种和研究场景

## User Journeys

### Journey 1: 首次分析 - 成功路径

**用户**: Yucheng（计算生物学研究员）

**背景故事**：
Yucheng 刚加入比较基因组学项目，需要分析狐猴物种与猕猴的基因差异。之前他手动运行各种工具，每次分析需要数小时，而且步骤繁琐容易出错。

**开场**：
Yucheng 收到了三个物种的基因组注释数据，需要完成直系同源推断和基因存在/缺失分析。他打开终端，准备开始工作。

**旅程步骤**：
1. **安装配置**（5分钟）
   - 克隆 CompGene 仓库
   - 安装依赖（conda/pip）
   - 验证安装成功

2. **准备输入数据**（10分钟）
   - 创建项目目录
   - 复制/链接基因组注释文件
   - 编辑 YAML 配置文件，指定物种和参数

3. **运行分析**（一条命令）
   - `compgene run config.yaml`
   - 系统自动：下载 NCBI 注释 → 运行 eggNOG → 生成比较矩阵

4. **监控进度**
   - 查看实时日志输出
   - 了解当前步骤和预计剩余时间

5. **获取结果**
   - 标准化的直系同源对应表
   - 基因存在/缺失比较矩阵
   - 汇总报告（HTML/Markdown）

**高潮时刻**：
运行完成后，Yucheng 看到清晰的比较矩阵和报告，之前需要一整天的工作现在 30 分钟内完成。

**结局**：
Yucheng 可以将结果直接用于下游分析，与 Mihir 分享进行共线性可视化。

---

### Journey 2: 添加新物种 - 扩展场景

**用户**: Yucheng

**背景故事**：
初始分析完成后，Kinnex 数据到达，需要将 4 个新狐猴物种加入分析。

**开场**：
Yucheng 收到新物种的基因组数据，需要扩展现有分析。

**旅程步骤**：
1. **准备新数据**（2分钟）
   - 将新物种注释文件放入数据目录

2. **更新配置**（2分钟）
   - 在 YAML 配置中添加新物种条目
   - 指定新物种的参数

3. **增量运行**
   - `compgene run config.yaml --incremental`
   - 系统自动识别新物种，复用已有结果

4. **获取更新结果**
   - 扩展后的 7 物种比较矩阵
   - 更新的汇总报告

**高潮时刻**：
仅需 5 分钟配置，无需重跑已完成的分析。

**结局**：
分析扩展无缝完成，Yucheng 可以继续研究更广泛的物种比较。

---

### Journey 3: 错误恢复 - 边缘情况

**用户**: Yucheng

**背景故事**：
分析运行到一半时，服务器网络中断，导致 NCBI 下载失败。

**开场**：
Yucheng 发现分析中断，担心需要从头开始。

**旅程步骤**：
1. **发现错误**
   - 查看日志，看到清晰的错误信息
   - 错误提示网络问题和恢复建议

2. **解决问题**
   - 等待网络恢复

3. **断点续跑**
   - `compgene run config.yaml --resume`
   - 系统从中断点继续，跳过已完成步骤

4. **完成分析**
   - 分析从断点继续完成
   - 结果与无中断运行一致

**高潮时刻**：
无需重跑已完成的耗时步骤，节省大量时间。

**结局**：
Yucheng 对系统的可靠性建立信心，知道任何中断都可以恢复。

---

### Journey 4: 团队协作 - 未来场景

**用户**: Mihir（理论研究/可视化）

**背景故事**：
Mihir 需要使用 Yucheng 的直系同源分析结果进行共线性可视化。

**开场**：
Mihir 收到 Yucheng 分享的分析结果目录。

**旅程步骤**：
1. **获取结果**
   - 访问共享的结果目录
   - 查看汇总报告了解分析概况

2. **提取所需数据**
   - 直系同源对应表（标准 TSV 格式）
   - 基因位置信息

3. **导入可视化工具**
   - 将标准格式数据导入 AnchorWave/SVbyEye
   - 进行共线性可视化

**高潮时刻**：
标准化的输出格式使得数据交接无缝衔接。

**结局**：
团队协作效率提升，各自专注于擅长的领域。

---

### Journey Requirements Summary

| 旅程 | 揭示的能力需求 |
|------|----------------|
| 首次分析 | 安装向导、配置模板、进度显示、结果报告 |
| 添加新物种 | 增量运行、配置热更新、结果复用 |
| 错误恢复 | 清晰错误信息、断点续跑、中间结果缓存 |
| 团队协作 | 标准输出格式、汇总报告、结果目录组织 |

## Domain-Specific Requirements

### 科学可复现性

**版本控制**：
- 记录所有依赖工具版本（eggNOG-mapper、Python 包等）
- 记录使用的参考数据库版本（NCBI、eggNOG 数据库）
- 配置文件完整保存分析参数

**结果验证**：
- 提供结果校验和（checksum）
- 支持与手动分析结果对比
- 记录随机种子（如有随机过程）

### 数据完整性

**输入验证**：
- 验证基因组注释文件格式（GFF/GTF）
- 检查物种标识符一致性
- 验证必需字段完整性

**输出质量**：
- 生成质量控制报告
- 标记低置信度的直系同源推断
- 提供数据完整性统计

### 计算资源管理

**资源需求**：
- 估算并报告内存需求
- 支持多核并行处理
- 提供进度和资源使用监控

**容错处理**：
- 处理大文件时分块读取
- 网络请求重试机制
- 磁盘空间检查

### 数据来源追溯

**审计日志**：
- 记录所有外部数据下载时间和来源
- 保存 NCBI 查询参数和响应
- 记录完整的分析时间线

## CLI Tool Specific Requirements

### Project-Type Overview

CompGene 是一个命令行驱动的生物信息学分析管道，设计为脚本自动化模式运行，支持批处理和工作流集成。

### Command Structure

**主命令**：
```
compgene <command> [options] <config_file>
```

**核心命令**：

| 命令 | 说明 |
|------|------|
| `run` | 运行完整分析流程 |
| `run --resume` | 从断点继续运行 |
| `run --incremental` | 增量运行（仅处理新物种） |
| `validate` | 验证配置文件 |
| `status` | 查看运行状态 |
| `clean` | 清理中间文件 |

**全局选项**：

| 选项 | 说明 |
|------|------|
| `--config, -c` | 指定配置文件路径 |
| `--output, -o` | 指定输出目录 |
| `--threads, -t` | 并行线程数 |
| `--verbose, -v` | 详细日志输出 |
| `--quiet, -q` | 静默模式 |
| `--dry-run` | 模拟运行，不执行实际操作 |

### Output Formats

**数据输出**：

| 模块 | 输出文件 | 格式 |
|------|----------|------|
| 注释检索 | `annotations/{species}.tsv` | TSV |
| 直系同源推断 | `orthologs/ortholog_table.tsv` | TSV |
| 基因存在/缺失 | `comparison/presence_absence_matrix.tsv` | TSV |
| RNA-Seq 分析 | `expression/DE_results.tsv` | TSV（DESeq2 标准格式） |
| RNA-Seq 分析 | `expression/normalized_counts.tsv` | TSV |

**报告输出**：

| 文件 | 格式 | 说明 |
|------|------|------|
| `report/summary.md` | Markdown | 分析汇总，便于版本控制 |
| `report/summary.html` | HTML | 可视化报告，便于浏览 |
| `logs/run.log` | 文本 | 人类可读日志 |
| `logs/run.json` | JSON | 机器可解析日志 |

### Config Schema

**YAML 配置文件结构**：

```yaml
project:
  name: "my_analysis"
  output_dir: "./results"

species:
  - name: "macaque"
    assembly: "MFA-T2T"
    annotation: "/path/to/macaque.gff"
  - name: "microcebus"
    annotation: "NCBI:GCF_000165445"  # 自动下载

analysis:
  orthology:
    method: "eggnog"
    database: "eukaryota"
  expression:
    enabled: true
    counts_matrix: "/path/to/counts.tsv"
    metadata: "/path/to/samples.tsv"

resources:
  threads: 8
  memory: "16G"
```

### Scripting Support

**退出码**：

| 退出码 | 含义 |
|--------|------|
| 0 | 成功完成 |
| 1 | 一般错误 |
| 2 | 配置错误 |
| 3 | 输入数据错误 |
| 4 | 网络错误（可重试） |
| 5 | 资源不足 |

**脚本集成**：
- 支持 stdin/stdout 管道（适用时）
- JSON 格式状态输出（`--json` 选项）
- 非交互模式默认，无需用户确认

## Project Scoping & Phased Development

### MVP Strategy & Philosophy

**MVP 类型**：问题解决型 MVP
**核心目标**：将手动多工具流程自动化为一条命令
**资源需求**：单人开发，本地服务器环境

### MVP Feature Set (Phase 1)

**核心用户旅程支持**：
- ✅ 首次分析 - 成功路径
- ✅ 错误恢复 - 断点续跑
- ⚠️ 添加新物种 - 基础支持（完整增量运行在 Phase 2）

**必需能力**：
1. NCBI 注释数据自动获取 + BUSCO QC
2. OrthoFinder orthogroup 推断（比较骨架）
3. eggNOG-mapper 功能注释（GO/KEGG/COG）
4. Liftoff 缺失核验（真缺失 vs 注释假缺失）
5. 多物种 presence/absence 比较矩阵生成
6. DESeq2 RNA-Seq 差异表达分析流程
7. YAML 配置驱动的工作流编排
8. 日志记录和断点续跑支持

### Post-MVP Features

**Phase 2 - 扩展阶段**：
- Kinnex 长读长数据支持
- 完整增量运行（智能识别新物种）
- 4 个额外狐猴物种验证
- 表达-直系同源综合分析报告

**Phase 3 - 愿景阶段**：
- 猿类物种支持
- CHM13 人类参考组装整合
- 完整灵长类比较基因组学框架
- 可复用平台化

### Risk Mitigation Strategy

**技术风险缓解**：
- eggNOG：先用 subprocess 调用，后续可优化为 API
- NCBI：实现本地缓存 + 指数退避重试
- 内存：大文件分块处理

**资源风险缓解**：
- 严格聚焦 MVP 范围
- 模块化实现，可独立测试
- 优先实现核心流程，边缘功能后补

## Functional Requirements

### 数据准备与 QC 能力

- FR1: 用户可通过 NCBI 物种标识符自动下载基因组注释数据
- FR2: 系统可解析 GFF3 和 GTF 格式的注释文件
- FR3: 系统可提取基因、转录本、CDS 等特征信息并标准化存储
- FR4: 系统可缓存已下载的注释数据，避免重复下载
- FR5: 用户可指定本地注释文件路径替代自动下载
- FR5a: 系统可调用 BUSCO 评估装配/注释完整性
- FR5b: 系统可生成每物种标准化目录结构和 BUSCO 汇总表

### 直系同源推断能力（OrthoFinder）

- FR6: 系统可调用 OrthoFinder 执行 orthogroup 推断
- FR7: 系统可生成 orthogroups 表（跨物种基因家族分组）
- FR8: 系统可输出物种树和基因树
- FR9: 系统可统计每物种的 orthogroup copy number
- FR10: 系统可识别物种特异性 orthogroups

### 功能注释能力（eggNOG-mapper）

- FR11: 系统可调用 eggNOG-mapper 对基因/蛋白进行功能注释
- FR12: 用户可指定 eggNOG 参考数据库（如 eukaryota、vertebrates）
- FR13: 系统可提取 GO/KEGG/COG 功能标签
- FR14: 系统可将功能注释汇总到 orthogroup 层面
- FR15: 系统可生成功能富集统计报告

### 缺失核验能力（Liftoff）

- FR16: 系统可调用 Liftoff 将参考注释映射到目标 assembly
- FR17: 系统可识别"真缺失"vs"注释假缺失"候选
- FR18: 系统可输出核验后的"真缺失候选清单"
- FR19: 系统可生成注释补全后的对照注释文件

### 基因存在/缺失分析能力

- FR20: 系统可从 orthogroups 派生 presence/absence 矩阵（0/1 或 copy number）
- FR21: 系统可比较多物种间的基因存在/缺失模式
- FR22: 系统可按基因家族分类统计
- FR23: 系统可识别物种特异性 orthogroups
- FR24: 系统可整合 Liftoff 核验结果标记可信度
- FR25: 系统可生成比较分析汇总报告

### RNA-Seq 表达分析能力

- FR26: 用户可提供 counts 矩阵和样本元数据作为输入
- FR27: 系统可调用 DESeq2 执行差异表达分析
- FR28: 系统可输出标准化后的表达矩阵
- FR29: 系统可输出差异表达结果（DESeq2 标准格式）
- FR30: 系统可将表达数据与 orthogroup 关联
- FR31: 系统可进行跨物种表达模式比较

### 工作流编排能力

- FR32: 用户可通过 YAML 配置文件定义完整分析流程
- FR33: 系统可自动解析模块依赖并按顺序执行
- FR34: 系统可缓存中间结果，支持断点续跑
- FR35: 用户可通过 --resume 标志从中断点继续
- FR36: 用户可通过 --incremental 标志仅处理新增数据
- FR37: 系统可支持 dry-run 模式预览执行计划

### 配置管理能力

- FR38: 用户可在配置文件中定义多个物种及其参数
- FR39: 系统可验证配置文件格式和必需字段
- FR40: 用户可通过命令行参数覆盖配置文件设置
- FR41: 系统可指定输出目录和资源限制（线程数、内存）

### 结果输出与报告能力

- FR42: 系统可输出所有数据结果为 TSV 格式
- FR43: 系统可生成 Markdown 格式汇总报告
- FR44: 系统可生成 HTML 格式可视化报告
- FR45: 系统可输出机器可解析的 JSON 日志
- FR46: 系统可输出人类可读的文本日志

### 错误处理与恢复能力

- FR47: 系统可在出错时提供清晰的错误信息和修复建议
- FR48: 系统可对网络请求实现指数退避重试
- FR49: 系统可返回标准化退出码（成功、配置错误、数据错误、网络错误、资源不足）

## Non-Functional Requirements

### 性能要求

- NFR1: 单物种完整分析流程应在 30 分钟内完成（标准服务器配置）
- NFR2: 大文件（>100MB）处理时内存占用应可控，支持分块读取
- NFR3: 多核并行处理时 CPU 利用率应达到 70% 以上
- NFR4: 中间结果缓存命中时，重复步骤跳过延迟 < 1 秒

### 可靠性要求

- NFR5: 分析结果 100% 可复现（相同输入 + 配置 = 相同输出）
- NFR6: 所有依赖工具版本和数据库版本必须记录在日志中
- NFR7: 支持从任意断点恢复，中间状态完整保存
- NFR8: 网络中断后自动重试，最多 3 次，指数退避间隔

### 集成要求

- NFR9: NCBI API 调用应遵循官方速率限制（每秒 3 请求）
- NFR10: OrthoFinder 调用应支持版本 2.5.x
- NFR11: eggNOG-mapper 调用应支持版本 2.1.x
- NFR12: Liftoff 调用应支持版本 1.6.x
- NFR13: BUSCO 调用应支持版本 5.x/6.x
- NFR14: DESeq2 输出格式应兼容 Bioconductor 标准
- NFR15: 所有外部工具调用应有超时保护（默认 30 分钟）

### 可维护性要求

- NFR16: 模块间低耦合，单个模块可独立测试
- NFR17: 日志级别可配置（DEBUG/INFO/WARNING/ERROR）
- NFR18: 配置文件变更无需修改代码
