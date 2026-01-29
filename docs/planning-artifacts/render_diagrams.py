#!/usr/bin/env python3
"""
Render Mermaid diagrams to PNG using mermaid.ink API.

Usage:
    python3 render_diagrams.py
"""

import base64
import re
import zlib
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import URLError

# Output directory
OUTPUT_DIR = Path(__file__).parent / "diagrams"
OUTPUT_DIR.mkdir(exist_ok=True)

# Mermaid diagrams - extracted and simplified for rendering
DIAGRAMS = {
    "01_system_overview": """
flowchart TB
    subgraph Input["Input Layer"]
        CONFIG[("config.yaml")]
        GFF["GFF3/GTF"]
        FASTA["FASTA"]
    end

    subgraph Engine["Workflow Engine"]
        SMK["Snakemake 9.14.8"]
        RULES["Rules (.smk)"]
        RUNNER["run_adapter.py"]
    end

    subgraph Adapters["Tool Adapter Layer"]
        BUSCO_A["BuscoAdapter"]
        OF_A["OrthoFinderAdapter"]
        EGG_A["EggNOGAdapter"]
        LIFT_A["LiftoffAdapter"]
    end

    subgraph Tools["External Tools"]
        BUSCO["BUSCO"]
        OF["OrthoFinder"]
        EGG["eggNOG-mapper"]
        LIFT["Liftoff"]
    end

    subgraph Libs["Core Libraries"]
        ERRORS["lib/errors.py"]
        AUDIT["lib/audit.py"]
        BIO["lib/bio_utils.py"]
    end

    subgraph Output["Output Layer"]
        FILES[("Files SoT")]
        SQLITE[("SQLite")]
        REPORTS["Reports"]
    end

    CONFIG --> SMK
    GFF --> SMK
    FASTA --> SMK
    SMK --> RULES --> RUNNER
    RUNNER --> BUSCO_A & OF_A & EGG_A & LIFT_A
    BUSCO_A --> BUSCO
    OF_A --> OF
    EGG_A --> EGG
    LIFT_A --> LIFT
    BUSCO_A & OF_A --> ERRORS
    RUNNER --> AUDIT --> FILES & SQLITE
    FILES --> REPORTS
""",

    "02_data_flow": """
flowchart LR
    subgraph P1["Phase 1: 数据准备"]
        INPUT["输入文件"]
        STD["standardize"]
        QC["qc_busco"]
    end

    subgraph P2["Phase 2: 直系同源"]
        ORTHO["orthology_infer"]
        OG["orthogroups.tsv"]
    end

    subgraph P3["Phase 3: 功能注释"]
        ANNOT["annotation_eggnog"]
        GO["GO/KEGG/COG"]
    end

    subgraph P4["Phase 4: 缺失核验"]
        LIFT["validation_liftoff"]
        ABS["true_absence"]
    end

    subgraph P5["Phase 5: 报告"]
        INGEST["ingest_db"]
        REPORT["reports"]
        HTML["summary.html"]
    end

    INPUT --> STD --> QC
    STD --> ORTHO --> OG
    OG --> ANNOT --> GO
    OG --> LIFT --> ABS
    OG & GO & ABS --> INGEST --> REPORT --> HTML
""",

    "03_module_deps": """
flowchart BT
    subgraph L4["Layer 4: Rules"]
        R1["standardize.smk"]
        R2["orthology.smk"]
        R3["validation.smk"]
    end

    subgraph L3["Layer 3: Scripts"]
        S1["run_adapter.py"]
    end

    subgraph L2["Layer 2: Adapters"]
        A0["base.py"]
        A1["orthofinder.py"]
        A2["liftoff.py"]
    end

    subgraph L1["Layer 1: Libraries"]
        LIB1["errors.py"]
        LIB2["audit.py"]
        LIB3["fasta.py"]
        LIB4["bio_utils.py"]
    end

    R1 & R2 & R3 --> S1
    S1 --> A0 & A1 & A2
    A1 & A2 --> A0
    A0 --> LIB1 & LIB2
    LIB4 --> LIB3
    LIB3 & LIB4 --> LIB1
""",

    "04_liftoff_workflow": """
flowchart TB
    subgraph Input["Input"]
        OG["presence_absence.tsv"]
        REF["参考注释 GFF3"]
        TGT["目标基因组 FA"]
    end

    subgraph Step1["Step 1: 识别缺失"]
        ID["识别缺失基因"]
        LIST["待核验列表"]
    end

    subgraph Step2["Step 2: Liftoff"]
        LIFT["Liftoff 映射"]
        STATS["coverage, identity"]
    end

    subgraph Step3["Step 3: 判定"]
        J{"判定规则"}
        M1["mapped ≥90%"]
        M2["partial 50-90%"]
        M3["partial <50%"]
        M4["unmapped"]
    end

    subgraph Step4["Step 4: 分类"]
        C1["false_absence"]
        C2["likely_false"]
        C3["likely_true"]
        C4["true_absence"]
    end

    subgraph Output["Output"]
        OUT1["true_absence_candidates.tsv"]
        OUT2["lifted_annotation.gff3"]
    end

    OG --> ID --> LIST
    LIST & REF & TGT --> LIFT --> STATS --> J
    J --> M1 & M2 & M3 & M4
    M1 --> C1 --> OUT2
    M2 --> C2 --> OUT1
    M3 --> C3 --> OUT1
    M4 --> C4 --> OUT1
""",

    "05_data_storage": """
flowchart TB
    subgraph Files["Source of Truth - Files"]
        F1["standardized/species/"]
        F2["orthology/species_set/"]
        F3["annotation/eggnog/"]
        F4["validation/liftoff/"]
        F5["meta/*.run.json"]
    end

    subgraph SQLite["Index of Record"]
        DB[("compgene.db")]
        T1["species"]
        T2["gene"]
        T3["orthogroup"]
        T4["run"]
    end

    subgraph Reports["Reports"]
        R1["summary.md"]
        R2["summary.html"]
    end

    F1 & F2 & F3 & F4 & F5 --> DB
    DB --> T1 & T2 & T3 & T4
    DB --> R1 & R2
""",

    "06_error_handling": """
flowchart TB
    subgraph Errors["Error Types"]
        E1["E_INPUT_MISSING"]
        E2["E_INPUT_FORMAT"]
        E3["E_TOOL_NOT_FOUND"]
        E4["E_TIMEOUT"]
        E5["E_NET_RATE_LIMIT"]
        E6["E_OOM"]
    end

    subgraph Decision["Retry Decision"]
        D{"可重试?"}
        YES["自动重试"]
        NO["停止+建议"]
    end

    subgraph Retry["Retry Strategy"]
        R1["Runner层 3次"]
        R2["Snakemake层 1次"]
    end

    E1 & E2 & E3 --> D
    E4 & E5 --> D
    E6 --> D
    D -->|Yes| YES --> R1 --> R2
    D -->|No| NO
""",

    "07_id_consistency": """
flowchart LR
    subgraph Input["Input Sources"]
        I1["proteins.fa.gz"]
        I2["Orthogroups.tsv"]
        I3["annotations.emapper"]
    end

    subgraph Extract["ID Extraction"]
        E1["extract_ids()"]
        E2["parse_fasta_id()"]
    end

    subgraph Compare["Comparison"]
        C1["compare_id_sets()"]
        C2["ComparisonResult"]
    end

    subgraph Output["Output"]
        O1["consistency_report.json"]
        O2{"status"}
        PASS["PASS ≥95%"]
        FAIL["FAIL <95%"]
    end

    I1 & I2 & I3 --> E1 --> E2 --> C1 --> C2 --> O1 --> O2
    O2 --> PASS & FAIL
""",

    "08_project_structure": """
flowchart LR
    ROOT["compgene/"]

    subgraph WF["workflow/"]
        W1["Snakefile"]
        W2["rules/"]
        W3["adapters/"]
        W4["lib/"]
        W5["envs/"]
    end

    subgraph CFG["config/"]
        C1["config.yaml"]
    end

    subgraph SCH["schemas/"]
        S1["config.schema.yaml"]
    end

    subgraph RES["results/"]
        R1["standardized/"]
        R2["orthology/"]
        R3["meta/"]
        R4["reports/"]
    end

    ROOT --> WF & CFG & SCH & RES
    WF --> W1 & W2 & W3 & W4 & W5
    CFG --> C1
    SCH --> S1
    RES --> R1 & R2 & R3 & R4
""",
}


def encode_mermaid(diagram: str) -> str:
    """Encode Mermaid diagram for mermaid.ink URL."""
    # Clean up the diagram
    diagram = diagram.strip()

    # Compress and encode
    compressed = zlib.compress(diagram.encode('utf-8'), 9)
    encoded = base64.urlsafe_b64encode(compressed).decode('ascii')

    return encoded


def render_diagram(name: str, diagram: str) -> bool:
    """Render a single diagram to PNG using mermaid.ink."""
    try:
        encoded = encode_mermaid(diagram)
        url = f"https://mermaid.ink/img/pako:{encoded}?type=png&bgColor=white"

        print(f"Rendering {name}...")

        req = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urlopen(req, timeout=30) as response:
            png_data = response.read()

        output_path = OUTPUT_DIR / f"{name}.png"
        output_path.write_bytes(png_data)
        print(f"  ✓ Saved to {output_path}")
        return True

    except URLError as e:
        print(f"  ✗ Network error: {e}")
        return False
    except Exception as e:
        print(f"  ✗ Error: {e}")
        return False


def main():
    """Render all diagrams."""
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Rendering {len(DIAGRAMS)} diagrams...\n")

    success = 0
    failed = 0

    for name, diagram in DIAGRAMS.items():
        if render_diagram(name, diagram):
            success += 1
        else:
            failed += 1

    print(f"\nDone: {success} succeeded, {failed} failed")

    if success > 0:
        print(f"\nPNG files saved to: {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
