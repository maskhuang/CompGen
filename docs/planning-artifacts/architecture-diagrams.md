# CompGene Architecture Diagrams

本文档包含 CompGene 项目的 Mermaid 架构图。

## 1. 系统总览

```mermaid
graph TB
    subgraph "CompGene Pipeline"
        subgraph "Input Layer"
            CONFIG[("config.yaml")]
            GFF["GFF3/GTF<br/>注释文件"]
            FASTA["FASTA<br/>基因组/蛋白"]
            COUNTS["counts.tsv<br/>表达矩阵"]
        end

        subgraph "Workflow Engine"
            SMK["Snakemake 9.14.8"]
            RULES["Rules<br/>(.smk files)"]
            RUNNER["run_adapter.py<br/>统一 Runner"]
        end

        subgraph "Tool Adapter Layer"
            BUSCO_A["BuscoAdapter"]
            OF_A["OrthoFinderAdapter"]
            EGG_A["EggNOGAdapter"]
            LIFT_A["LiftoffAdapter"]
            DEQ_A["DESeq2Adapter"]
        end

        subgraph "External Tools"
            BUSCO["BUSCO 5.x/6.x"]
            OF["OrthoFinder 2.5.x"]
            EGG["eggNOG-mapper 2.1.x"]
            LIFT["Liftoff 1.6.x"]
            DEQ["DESeq2 (R)"]
        end

        subgraph "Core Libraries"
            ERRORS["lib/errors.py<br/>错误码系统"]
            AUDIT["lib/audit.py<br/>审计元数据"]
            GFF_LIB["lib/gff.py<br/>GFF解析"]
            FASTA_LIB["lib/fasta.py<br/>FASTA处理"]
            BIO_UTILS["lib/bio_utils.py<br/>Biopython封装"]
        end

        subgraph "Output Layer"
            FILES[("Files<br/>Source of Truth")]
            SQLITE[("SQLite<br/>Index of Record")]
            REPORTS["Reports<br/>HTML/Markdown"]
        end
    end

    CONFIG --> SMK
    GFF --> SMK
    FASTA --> SMK
    COUNTS --> SMK

    SMK --> RULES
    RULES --> RUNNER

    RUNNER --> BUSCO_A
    RUNNER --> OF_A
    RUNNER --> EGG_A
    RUNNER --> LIFT_A
    RUNNER --> DEQ_A

    BUSCO_A --> BUSCO
    OF_A --> OF
    EGG_A --> EGG
    LIFT_A --> LIFT
    DEQ_A --> DEQ

    BUSCO_A --> ERRORS
    OF_A --> ERRORS
    EGG_A --> ERRORS
    LIFT_A --> ERRORS
    DEQ_A --> ERRORS

    RUNNER --> AUDIT
    AUDIT --> FILES
    AUDIT --> SQLITE
    FILES --> REPORTS
    SQLITE --> REPORTS

    classDef input fill:#e1f5fe,stroke:#01579b
    classDef engine fill:#fff3e0,stroke:#e65100
    classDef adapter fill:#f3e5f5,stroke:#7b1fa2
    classDef tool fill:#e8f5e9,stroke:#2e7d32
    classDef lib fill:#fce4ec,stroke:#c2185b
    classDef output fill:#e0f2f1,stroke:#00695c

    class CONFIG,GFF,FASTA,COUNTS input
    class SMK,RULES,RUNNER engine
    class BUSCO_A,OF_A,EGG_A,LIFT_A,DEQ_A adapter
    class BUSCO,OF,EGG,LIFT,DEQ tool
    class ERRORS,AUDIT,GFF_LIB,FASTA_LIB,BIO_UTILS lib
    class FILES,SQLITE,REPORTS output
```

## 2. 数据流图 (Pipeline DAG)

```mermaid
flowchart LR
    subgraph "Phase 1: 数据准备"
        INPUT["输入文件<br/>GFF/FASTA"]
        STD["standardize<br/>数据标准化"]
        QC["qc_busco<br/>质量控制"]
    end

    subgraph "Phase 2: 直系同源分析"
        ORTHO["orthology_infer<br/>OrthoFinder"]
        OG_TABLE["orthogroups.tsv"]
        TREES["species_tree.nwk<br/>gene_trees/"]
    end

    subgraph "Phase 3: 功能注释"
        ANNOT["annotation_eggnog<br/>功能注释"]
        GO_KEGG["GO/KEGG/COG<br/>annotations.tsv"]
    end

    subgraph "Phase 4: 缺失核验"
        LIFTOFF["validation_liftoff<br/>Liftoff核验"]
        ABSENCE["true_absence<br/>_candidates.tsv"]
    end

    subgraph "Phase 5: 矩阵生成"
        MATRIX["matrices_generate<br/>存在/缺失矩阵"]
        PA["presence_absence.tsv"]
        CN["copy_number.tsv"]
    end

    subgraph "Phase 6: 表达分析"
        EXPR["expression_deseq2<br/>差异表达"]
        DE["DE_results.tsv"]
    end

    subgraph "Phase 7: 报告生成"
        INGEST["ingest_db<br/>SQLite写入"]
        REPORT["reporting_summary<br/>报告生成"]
        HTML["summary.html"]
    end

    INPUT --> STD
    STD --> QC
    STD --> ORTHO

    ORTHO --> OG_TABLE
    ORTHO --> TREES

    OG_TABLE --> ANNOT
    ANNOT --> GO_KEGG

    OG_TABLE --> MATRIX
    MATRIX --> PA
    MATRIX --> CN

    PA --> LIFTOFF
    LIFTOFF --> ABSENCE

    STD --> EXPR
    EXPR --> DE

    OG_TABLE --> INGEST
    GO_KEGG --> INGEST
    PA --> INGEST
    ABSENCE --> INGEST
    DE --> INGEST

    INGEST --> REPORT
    REPORT --> HTML

    style INPUT fill:#e3f2fd
    style HTML fill:#c8e6c9
```

## 3. 模块依赖图

```mermaid
graph BT
    subgraph "Layer 4: Rules"
        R_STD["rules/standardize.smk"]
        R_QC["rules/qc.smk"]
        R_ORTHO["rules/orthology.smk"]
        R_ANNOT["rules/annotation.smk"]
        R_VAL["rules/validation.smk"]
        R_MAT["rules/matrices.smk"]
        R_EXPR["rules/expression.smk"]
        R_REP["rules/reporting.smk"]
        R_AUDIT["rules/audit.smk"]
    end

    subgraph "Layer 3: Scripts"
        S_RUN["scripts/run_adapter.py"]
        S_INGEST["scripts/ingest_db.py"]
        S_REPORT["scripts/generate_report.py"]
    end

    subgraph "Layer 2: Adapters"
        A_BASE["adapters/base.py"]
        A_BUSCO["adapters/busco.py"]
        A_OF["adapters/orthofinder.py"]
        A_EGG["adapters/eggnog.py"]
        A_LIFT["adapters/liftoff.py"]
        A_DEQ["adapters/deseq2.py"]
    end

    subgraph "Layer 1: Libraries"
        L_ERR["lib/errors.py"]
        L_IO["lib/io.py"]
        L_AUDIT["lib/audit.py"]
        L_GFF["lib/gff.py"]
        L_FASTA["lib/fasta.py"]
        L_BIO["lib/bio_utils.py"]
        L_PA["lib/presence_absence.py"]
        L_REP["lib/report.py"]
    end

    %% Rules depend on Scripts
    R_STD --> S_RUN
    R_QC --> S_RUN
    R_ORTHO --> S_RUN
    R_ANNOT --> S_RUN
    R_VAL --> S_RUN
    R_EXPR --> S_RUN
    R_MAT --> S_RUN
    R_REP --> S_REPORT
    R_AUDIT --> L_BIO

    %% Scripts depend on Adapters
    S_RUN --> A_BASE
    S_RUN --> A_BUSCO
    S_RUN --> A_OF
    S_RUN --> A_EGG
    S_RUN --> A_LIFT
    S_RUN --> A_DEQ

    %% Adapters depend on base and libs
    A_BUSCO --> A_BASE
    A_OF --> A_BASE
    A_EGG --> A_BASE
    A_LIFT --> A_BASE
    A_DEQ --> A_BASE
    A_BASE --> L_ERR
    A_BASE --> L_AUDIT

    %% Libraries internal deps
    L_AUDIT --> L_ERR
    L_IO --> L_ERR
    L_GFF --> L_ERR
    L_FASTA --> L_ERR
    L_BIO --> L_FASTA
    L_BIO --> L_ERR
    L_PA --> L_ERR
    L_REP --> L_ERR

    classDef rules fill:#fff3e0,stroke:#e65100
    classDef scripts fill:#e8f5e9,stroke:#2e7d32
    classDef adapters fill:#f3e5f5,stroke:#7b1fa2
    classDef libs fill:#e3f2fd,stroke:#1565c0

    class R_STD,R_QC,R_ORTHO,R_ANNOT,R_VAL,R_MAT,R_EXPR,R_REP,R_AUDIT rules
    class S_RUN,S_INGEST,S_REPORT scripts
    class A_BASE,A_BUSCO,A_OF,A_EGG,A_LIFT,A_DEQ adapters
    class L_ERR,L_IO,L_AUDIT,L_GFF,L_FASTA,L_BIO,L_PA,L_REP libs
```

## 4. 外部工具集成

```mermaid
flowchart TB
    subgraph "CompGene"
        RUNNER["run_adapter.py"]

        subgraph "Adapters"
            A1["BuscoAdapter"]
            A2["OrthoFinderAdapter"]
            A3["EggNOGAdapter"]
            A4["LiftoffAdapter"]
            A5["DESeq2Adapter"]
        end
    end

    subgraph "Conda Environments"
        E1["envs/busco.yaml<br/>BUSCO 5.x/6.x"]
        E2["envs/orthofinder.yaml<br/>OrthoFinder 2.5.x"]
        E3["envs/eggnog.yaml<br/>eggNOG-mapper 2.1.x"]
        E4["envs/liftoff.yaml<br/>Liftoff 1.6.x"]
        E5["envs/deseq2.yaml<br/>R + DESeq2"]
    end

    subgraph "Tool Execution"
        T1["busco<br/>-i proteins.fa<br/>-l eukaryota"]
        T2["orthofinder<br/>-f proteins_dir/"]
        T3["emapper.py<br/>-i proteins.fa"]
        T4["liftoff<br/>-g ref.gff<br/>target.fa"]
        T5["Rscript<br/>deseq2_analysis.R"]
    end

    subgraph "Outputs"
        O1["busco_summary.txt"]
        O2["Orthogroups.tsv<br/>SpeciesTree_*.txt"]
        O3["*.emapper.annotations"]
        O4["lifted_annotation.gff3"]
        O5["DE_results.tsv"]
    end

    RUNNER --> A1 & A2 & A3 & A4 & A5

    A1 --> E1 --> T1 --> O1
    A2 --> E2 --> T2 --> O2
    A3 --> E3 --> T3 --> O3
    A4 --> E4 --> T4 --> O4
    A5 --> E5 --> T5 --> O5

    style RUNNER fill:#fff3e0
    style E1,E2,E3,E4,E5 fill:#e8f5e9
    style T1,T2,T3,T4,T5 fill:#e3f2fd
    style O1,O2,O3,O4,O5 fill:#f3e5f5
```

## 5. Epic 5: Liftoff 缺失核验工作流

```mermaid
flowchart TB
    subgraph "Input"
        OG["OrthoFinder 输出<br/>presence_absence.tsv"]
        REF_GFF["参考物种<br/>annotation.gff3"]
        TGT_FA["目标物种<br/>genome.fa"]
    end

    subgraph "Step 1: 识别缺失候选"
        IDENTIFY["识别缺失基因<br/>Species B 缺少 Gene X"]
        TARGETS["liftoff_targets.txt<br/>待核验基因列表"]
    end

    subgraph "Step 2: Liftoff 映射"
        LIFTOFF["Liftoff<br/>A → B 注释映射"]
        STATS["liftoff_stats.tsv<br/>coverage, identity"]
    end

    subgraph "Step 3: 判定规则"
        JUDGE{"判定逻辑"}

        MAPPED["mapped<br/>coverage ≥90%<br/>identity ≥85%"]
        PARTIAL_HIGH["partial<br/>coverage 50-90%"]
        PARTIAL_LOW["partial<br/>coverage <50%"]
        UNMAPPED["unmapped<br/>无法映射"]
    end

    subgraph "Step 4: 分类结果"
        FALSE_ABS["false_absence<br/>注释假缺失<br/>confidence: very_low"]
        LIKELY_FALSE["likely_false<br/>可能假缺失<br/>confidence: low"]
        LIKELY_TRUE["likely_true<br/>可能真缺失<br/>confidence: medium"]
        TRUE_ABS["true_absence<br/>真缺失候选<br/>confidence: high"]
    end

    subgraph "Output"
        CANDIDATES["true_absence_candidates.tsv"]
        LIFTED["lifted_annotation.gff3<br/>补全注释"]
        VERIFIED["presence_absence_verified.tsv<br/>含 confidence 列"]
    end

    OG --> IDENTIFY
    IDENTIFY --> TARGETS

    TARGETS --> LIFTOFF
    REF_GFF --> LIFTOFF
    TGT_FA --> LIFTOFF
    LIFTOFF --> STATS

    STATS --> JUDGE

    JUDGE --> MAPPED
    JUDGE --> PARTIAL_HIGH
    JUDGE --> PARTIAL_LOW
    JUDGE --> UNMAPPED

    MAPPED --> FALSE_ABS
    PARTIAL_HIGH --> LIKELY_FALSE
    PARTIAL_LOW --> LIKELY_TRUE
    UNMAPPED --> TRUE_ABS

    FALSE_ABS --> LIFTED
    LIKELY_FALSE --> CANDIDATES
    LIKELY_TRUE --> CANDIDATES
    TRUE_ABS --> CANDIDATES

    CANDIDATES --> VERIFIED
    LIFTED --> VERIFIED

    style OG fill:#e3f2fd
    style LIFTOFF fill:#fff3e0
    style TRUE_ABS fill:#ffcdd2
    style FALSE_ABS fill:#c8e6c9
    style CANDIDATES fill:#f3e5f5
    style VERIFIED fill:#e0f2f1
```

## 6. 数据存储架构

```mermaid
flowchart TB
    subgraph "Source of Truth (Files)"
        subgraph "standardized/{species}/"
            GENOME["genome.fa.gz"]
            ANNOT["annotation.gff3.gz"]
            PROTEIN["proteins.longest.fa.gz"]
            SUMMARY["*.summary.json"]
        end

        subgraph "orthology/{species_set}/"
            OG_TSV["orthogroups.tsv"]
            TREE["species_tree.nwk"]
            GENE_TREES["gene_trees/"]
        end

        subgraph "annotation/eggnog/{species}/"
            EGG_OUT["annotations.tsv"]
        end

        subgraph "validation/liftoff/"
            LIFTED["lifted_annotation.gff3"]
            TRUE_ABS["true_absence_candidates.tsv"]
        end

        subgraph "meta/{rule}/"
            RUN_JSON["*.run.json<br/>审计记录"]
        end
    end

    subgraph "Index of Record (SQLite)"
        DB[("compgene.db")]

        T_SPECIES["species"]
        T_GENE["gene"]
        T_TRANSCRIPT["transcript"]
        T_PROTEIN["protein"]
        T_OG["orthogroup"]
        T_OG_MEMBER["orthogroup_member"]
        T_QC["qc_busco"]
        T_RUN["run"]
        T_ARTIFACT["artifact"]
        T_VERSION["tool_version"]
    end

    subgraph "Reports"
        MD["summary.md"]
        HTML["summary.html"]
    end

    GENOME --> DB
    ANNOT --> DB
    PROTEIN --> DB
    OG_TSV --> DB
    EGG_OUT --> DB
    TRUE_ABS --> DB
    RUN_JSON --> DB

    DB --> T_SPECIES
    DB --> T_GENE
    DB --> T_TRANSCRIPT
    DB --> T_PROTEIN
    DB --> T_OG
    DB --> T_OG_MEMBER
    DB --> T_QC
    DB --> T_RUN
    DB --> T_ARTIFACT
    DB --> T_VERSION

    DB --> MD
    DB --> HTML

    style DB fill:#fff3e0,stroke:#e65100
    style HTML fill:#c8e6c9
```

## 7. 错误处理与重试策略

```mermaid
flowchart TB
    subgraph "Error Classification"
        E1["E_INPUT_MISSING<br/>输入缺失"]
        E2["E_INPUT_FORMAT<br/>格式错误"]
        E3["E_TOOL_NOT_FOUND<br/>工具缺失"]
        E4["E_TOOL_VERSION<br/>版本不符"]
        E5["E_TIMEOUT<br/>超时"]
        E6["E_OOM<br/>内存不足"]
        E7["E_NET_RATE_LIMIT<br/>网络限流"]
        E8["E_DISK_FULL<br/>磁盘满"]
        E9["E_NONZERO_EXIT<br/>其他错误"]
    end

    subgraph "Retry Decision"
        RETRY{"可重试?"}
        YES["自动重试<br/>指数退避"]
        NO["停止执行<br/>输出恢复建议"]
    end

    subgraph "Retry Strategy"
        RUNNER_RETRY["Runner 层重试<br/>网络错误 3次<br/>exponential backoff"]
        SMK_RETRY["Snakemake 层重试<br/>retries: 1"]
    end

    subgraph "Recovery Suggestions"
        R1["检查输入文件路径"]
        R2["验证 GFF/FASTA 格式"]
        R3["激活 conda 环境"]
        R4["更新工具版本"]
        R5["增加超时时间"]
        R6["减少 threads"]
        R7["等待后重试"]
        R8["清理磁盘空间"]
        R9["查看日志"]
    end

    E1 --> RETRY
    E2 --> RETRY
    E3 --> RETRY
    E4 --> RETRY
    E5 --> RETRY
    E6 --> RETRY
    E7 --> RETRY
    E8 --> RETRY
    E9 --> RETRY

    RETRY -->|Yes| YES
    RETRY -->|No| NO

    YES --> RUNNER_RETRY
    RUNNER_RETRY --> SMK_RETRY

    E1 --> R1
    E2 --> R2
    E3 --> R3
    E4 --> R4
    E5 --> R5
    E6 --> R6
    E7 --> R7
    E8 --> R8
    E9 --> R9

    style E5 fill:#c8e6c9
    style E7 fill:#c8e6c9
    style YES fill:#c8e6c9
    style NO fill:#ffcdd2
```

## 8. ID 一致性检查流程 (Story 2.2b)

```mermaid
flowchart LR
    subgraph "Input Sources"
        FASTA["proteins.longest.fa.gz<br/>标准化蛋白"]
        OF_OUT["Orthogroups.tsv<br/>OrthoFinder输出"]
        EGG_OUT["annotations.emapper<br/>eggNOG输出"]
    end

    subgraph "ID Extraction"
        EX1["extract_ids()<br/>从 FASTA"]
        EX2["extract_ids_from_tsv()<br/>从 OrthoFinder"]
        EX3["extract_ids_from_tsv()<br/>从 eggNOG"]
    end

    subgraph "ID Parsing"
        PARSE["parse_fasta_id()<br/>多格式解析"]

        FMT1["UniProt<br/>sp|P12345|GENE"]
        FMT2["NCBI<br/>ref|NP_001234|"]
        FMT3["OrthoFinder<br/>species|gene_id"]
        FMT4["Simple<br/>gene_id desc"]
    end

    subgraph "Comparison"
        CMP["compare_id_sets()"]
        RESULT["ComparisonResult"]

        COMMON["common: 共有ID"]
        ONLY_A["only_in_a: 仅输入"]
        ONLY_B["only_in_b: 仅输出"]
        RATE["match_rate: 匹配率"]
    end

    subgraph "Output"
        REPORT["consistency_report.json"]
        STATUS{"status"}
        PASS["PASS<br/>rate ≥ 95%"]
        FAIL["FAIL<br/>rate < 95%"]
    end

    FASTA --> EX1
    OF_OUT --> EX2
    EGG_OUT --> EX3

    EX1 --> PARSE
    PARSE --> FMT1 & FMT2 & FMT3 & FMT4

    EX1 --> CMP
    EX2 --> CMP
    EX3 --> CMP

    CMP --> RESULT
    RESULT --> COMMON & ONLY_A & ONLY_B & RATE

    RESULT --> REPORT
    REPORT --> STATUS
    STATUS -->|≥95%| PASS
    STATUS -->|<95%| FAIL

    style PASS fill:#c8e6c9
    style FAIL fill:#ffcdd2
```

## 9. 项目目录结构

```mermaid
graph LR
    subgraph "compgene/"
        ROOT["compgene/"]

        subgraph "workflow/"
            WF["workflow/"]
            SMK["Snakefile"]
            RULES["rules/*.smk"]
            SCRIPTS["scripts/*.py"]
            ADAPTERS["adapters/*.py"]
            LIB["lib/*.py"]
            ENVS["envs/*.yaml"]
        end

        subgraph "config/"
            CFG["config/"]
            CFG_YAML["config.yaml"]
            SPECIES["species.tsv"]
        end

        subgraph "schemas/"
            SCH["schemas/"]
            CFG_SCH["config.schema.yaml"]
            RUN_SCH["run.schema.json"]
        end

        subgraph "tests/"
            TST["tests/"]
            INT["integration/"]
            FIX["fixtures/"]
        end

        subgraph "results/"
            RES["results/"]
            STD["standardized/"]
            ORTHO["orthology/"]
            ANNOT["annotation/"]
            VAL["validation/"]
            META["meta/"]
            REP["reports/"]
        end

        subgraph "logs/"
            LOG["logs/"]
            LOG_RULE["<rule>/"]
            LOG_FILE["*.log, *.jsonl"]
        end
    end

    ROOT --> WF & CFG & SCH & TST & RES & LOG

    WF --> SMK & RULES & SCRIPTS & ADAPTERS & LIB & ENVS
    CFG --> CFG_YAML & SPECIES
    SCH --> CFG_SCH & RUN_SCH
    TST --> INT & FIX
    RES --> STD & ORTHO & ANNOT & VAL & META & REP
    LOG --> LOG_RULE --> LOG_FILE

    style ROOT fill:#e3f2fd
    style WF fill:#fff3e0
    style RES fill:#e8f5e9
```

---

_Generated: 2026-01-29_
_Source: architecture.md_
