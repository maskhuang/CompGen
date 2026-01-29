#!/usr/bin/env python3
"""
Render Mermaid diagrams to PNG using Kroki API.

Usage:
    python3 render_diagrams_kroki.py
"""

import base64
import zlib
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import URLError

# Output directory
OUTPUT_DIR = Path(__file__).parent / "diagrams"
OUTPUT_DIR.mkdir(exist_ok=True)

# Simplified Mermaid diagrams for better compatibility
DIAGRAMS = {
    "01_system_overview": """
graph TB
    CONFIG["config.yaml"] --> SMK["Snakemake"]
    GFF["GFF3/GTF"] --> SMK
    FASTA["FASTA"] --> SMK

    SMK --> RUNNER["run_adapter.py"]

    RUNNER --> A1["BuscoAdapter"]
    RUNNER --> A2["OrthoFinderAdapter"]
    RUNNER --> A3["EggNOGAdapter"]
    RUNNER --> A4["LiftoffAdapter"]

    A1 --> T1["BUSCO"]
    A2 --> T2["OrthoFinder"]
    A3 --> T3["eggNOG-mapper"]
    A4 --> T4["Liftoff"]

    RUNNER --> AUDIT["lib/audit.py"]
    AUDIT --> FILES["Files"]
    AUDIT --> SQLITE["SQLite"]
    FILES --> REPORTS["Reports"]
""",

    "02_data_flow": """
graph LR
    INPUT["输入文件"] --> STD["standardize"]
    STD --> QC["qc_busco"]
    STD --> ORTHO["orthology_infer"]
    ORTHO --> OG["orthogroups.tsv"]
    OG --> ANNOT["annotation_eggnog"]
    ANNOT --> GO["GO/KEGG/COG"]
    OG --> LIFT["validation_liftoff"]
    LIFT --> ABS["true_absence"]
    OG --> INGEST["ingest_db"]
    GO --> INGEST
    ABS --> INGEST
    INGEST --> REPORT["summary.html"]
""",

    "03_module_deps": """
graph BT
    R1["rules/*.smk"] --> S1["run_adapter.py"]
    S1 --> A0["adapters/base.py"]
    S1 --> A1["adapters/orthofinder.py"]
    S1 --> A2["adapters/liftoff.py"]
    A1 --> A0
    A2 --> A0
    A0 --> L1["lib/errors.py"]
    A0 --> L2["lib/audit.py"]
    L3["lib/bio_utils.py"] --> L4["lib/fasta.py"]
    L4 --> L1
""",

    "04_liftoff_workflow": """
graph TB
    OG["presence_absence.tsv"] --> ID["识别缺失基因"]
    ID --> LIST["待核验列表"]
    REF["参考注释"] --> LIFT["Liftoff映射"]
    TGT["目标基因组"] --> LIFT
    LIST --> LIFT
    LIFT --> STATS["coverage/identity"]
    STATS --> J{"判定规则"}
    J -->|"≥90%"| C1["false_absence"]
    J -->|"50-90%"| C2["likely_false"]
    J -->|"<50%"| C3["likely_true"]
    J -->|"unmapped"| C4["true_absence"]
    C1 --> OUT2["lifted_annotation.gff3"]
    C2 --> OUT1["true_absence_candidates.tsv"]
    C3 --> OUT1
    C4 --> OUT1
""",

    "05_data_storage": """
graph TB
    subgraph Files["Source of Truth"]
        F1["standardized/"]
        F2["orthology/"]
        F3["annotation/"]
        F4["meta/*.run.json"]
    end

    F1 --> DB["compgene.db"]
    F2 --> DB
    F3 --> DB
    F4 --> DB

    DB --> T1["species table"]
    DB --> T2["gene table"]
    DB --> T3["orthogroup table"]

    DB --> R1["summary.md"]
    DB --> R2["summary.html"]
""",

    "06_error_handling": """
graph TB
    E1["E_INPUT_MISSING"] --> D{"可重试?"}
    E2["E_TIMEOUT"] --> D
    E3["E_NET_RATE_LIMIT"] --> D
    E4["E_OOM"] --> D

    D -->|Yes| YES["自动重试"]
    D -->|No| NO["停止+建议"]

    YES --> R1["Runner层 3次"]
    R1 --> R2["Snakemake层 1次"]
""",

    "07_id_consistency": """
graph LR
    I1["proteins.fa.gz"] --> E1["extract_ids"]
    I2["Orthogroups.tsv"] --> E1
    E1 --> P["parse_fasta_id"]
    P --> C["compare_id_sets"]
    C --> R["ComparisonResult"]
    R --> O["consistency_report.json"]
    O --> S{"status"}
    S -->|"≥95%"| PASS["PASS"]
    S -->|"<95%"| FAIL["FAIL"]
""",

    "08_project_structure": """
graph LR
    ROOT["compgene/"] --> WF["workflow/"]
    ROOT --> CFG["config/"]
    ROOT --> SCH["schemas/"]
    ROOT --> RES["results/"]

    WF --> W1["Snakefile"]
    WF --> W2["rules/"]
    WF --> W3["adapters/"]
    WF --> W4["lib/"]

    RES --> R1["standardized/"]
    RES --> R2["orthology/"]
    RES --> R3["reports/"]
""",
}


def encode_kroki(diagram: str) -> str:
    """Encode diagram for Kroki URL."""
    compressed = zlib.compress(diagram.encode('utf-8'), 9)
    encoded = base64.urlsafe_b64encode(compressed).decode('ascii')
    return encoded


def render_with_kroki(name: str, diagram: str) -> bool:
    """Render diagram using Kroki service."""
    try:
        encoded = encode_kroki(diagram)
        url = f"https://kroki.io/mermaid/png/{encoded}"

        print(f"Rendering {name}...")

        req = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urlopen(req, timeout=60) as response:
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
    print(f"Rendering {len(DIAGRAMS)} diagrams using Kroki...\n")

    success = 0
    failed = 0

    for name, diagram in DIAGRAMS.items():
        if render_with_kroki(name, diagram):
            success += 1
        else:
            failed += 1

    print(f"\nDone: {success} succeeded, {failed} failed")

    if success > 0:
        print(f"\nPNG files saved to: {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
