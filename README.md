# synviz — Synteny visualization for any DNA sequences

A command-line pipeline for comparing two DNA sequences, generating:
- A **dual-axis synteny map** with ORF tracks, genomic island bands, key gene highlights, and curved light-blue ribbons connecting conserved (>90% identity) blocks
- A **standalone ribbon plot** showing syntenic connections between the two genomes
- A **per-window %identity plot** (500 bp sliding blastn windows)
- **Tabular data** with per-window blastn results and merged high-identity regions

![Example output](plasmid_comparison_map.svg)

---

## Software Requirements

| Software | Version | Purpose |
|----------|---------|---------|
| **Python** | ≥ 3.9 | Pipeline orchestration & plotting |
| **BLAST+** | ≥ 2.12 | `blastn` and `makeblastdb` for per-window identity |
| **matplotlib** | ≥ 3.5 | SVG figure generation |
| **numpy** | ≥ 1.21 | Numerical arrays |

All Python dependencies can be installed with:

```bash
pip install matplotlib numpy
```

BLAST+ must be installed separately. On macOS:

```bash
brew install blast
```

On Ubuntu/Debian:

```bash
sudo apt install ncbi-blast+
```

---

## Quick Start

```bash
# Show help (safe — no execution)
python scripts/run_synteny_pipeline.py

# Full run with all inputs
python scripts/run_synteny_pipeline.py --run \
    --ref-fasta ref.fasta --ref-gff ref.gff3 \
    --qry-fasta qry.fasta --qry-gff qry.gff3 \
    --islands islands.tsv --suffix my_comparison
```

Without `--run`, the script prints help and exits — it never runs by default.

---

## Usage

### Full example (all inputs)

```bash
python scripts/run_synteny_pipeline.py --run \
    --ref-fasta path/to/reference.fasta \
    --ref-gff   path/to/reference.gff3 \
    --qry-fasta path/to/query.fasta \
    --qry-gff   path/to/query.gff3 \
    --islands   path/to/islands.tsv \
    --suffix    my_comparison
```

This generates:

```
plasmid_blastn_identity_windows_my_comparison.tsv
plasmid_identity_plot_my_comparison.svg
plasmid_high_identity_regions_my_comparison.tsv
plasmid_comparison_map_my_comparison.svg
plasmid_ribbon_synteny_my_comparison.svg
```

### FASTA-only (no annotations)

**GFF files are optional.** If omitted, the comparison map draws genome bars and ribbons only (no ORF tracks, no island bands).

```bash
python scripts/run_synteny_pipeline.py --run \
    --ref-fasta ref.fasta \
    --qry-fasta qry.fasta \
    --suffix fasta_only
```

### Tagged outputs (preserve previous plots)

```bash
python scripts/run_synteny_pipeline.py --run \\
    --ref-fasta ref.fasta --qry-fasta qry.fasta \\
    --suffix comparison_A
```

All output files are tagged with `_comparison_A`:

```
plasmid_blastn_identity_windows_comparison_A.tsv
plasmid_identity_plot_comparison_A.svg
plasmid_high_identity_regions_comparison_A.tsv
plasmid_comparison_map_comparison_A.svg
plasmid_ribbon_synteny_comparison_A.svg
```

### Genomic island bands

Control the coloured annotation bands drawn behind the ORF tracks:

```bash
# Auto-generate from GFF annotations (default when GFF is present)
python scripts/run_synteny_pipeline.py --run \\
    --ref-fasta ref.fasta --ref-gff ref.gff3 \\
    --qry-fasta qry.fasta --qry-gff qry.gff3 \\
    --islands auto --suffix auto_islands

# No island bands
python scripts/run_synteny_pipeline.py --run \\
    --ref-fasta ref.fasta --qry-fasta qry.fasta \\
    --islands none --suffix plain

# Read island definitions from a TSV file
python scripts/run_synteny_pipeline.py --run \\
    --ref-fasta ref.fasta --ref-gff ref.gff3 \\
    --qry-fasta qry.fasta --qry-gff qry.gff3 \\
    --islands path/to/islands.tsv --suffix with_islands
```

**Island file format** (TSV):

```
# Comments start with #
# start  end     label                   color (optional)
0        10000   Rep / IS elements       #377EB8
35000    50000   DNA metabolism          #A6761D
70000    89000   T7SS / ESX-1 cluster    #984EA3
169000   177000  Conjugation cluster     #E41A1C
```

If the colour column is omitted, it is auto-assigned from the label text.

### Step control

| Flag | Effect |
|------|--------|
| `--skip blastn` | Re-use existing blastn data (saves ~1–2 min) |
| `--only map` | Regenerate only the comparison map (auto-resolves dependencies) |
| `--force` | Re-run steps even if outputs already exist |
| `--dry-run` | Preview what would run without executing |

---

## Pipeline Steps

The pipeline runs four steps in dependency order:

| Step | Script | Input | Output |
|------|--------|-------|--------|
| **1. blastn windows** | `blastn_identity_windows.py` | Ref + query FASTA | `plasmid_blastn_identity_windows.tsv`, `plasmid_identity_plot.svg` |
| **2. high-identity regions** | `high_identity_regions.py` | blastn TSV | `plasmid_high_identity_regions.tsv` |
| **3. comparison map** | `plasmid_comparison_map.py` | GFFs + regions TSV | `plasmid_comparison_map.svg` |
| **4. standalone ribbon** | `ribbon_synteny_plot.py` | regions TSV | `plasmid_ribbon_synteny.svg` |

### Step 1 — blastn identity windows

- Slides a 500 bp window across the query genome (250 bp step)
- Each window is blasted against the reference
- The best HSP (by bitscore) is recorded per window
- Output: TSV with query position, ref position, %identity, strand, bitscore
- Also generates a continuous %identity line plot (y-axis 60–100%)

### Step 2 — high-identity regions

- Reads the per-window TSV
- Merges consecutive windows with >90% identity into contiguous blocks
- Reports each block's coordinates in both genomes
- Output: 44 conserved blocks totalling ~111 kb (57% of query)

### Step 3 — comparison map (dual-axis synteny)

- Two horizontal tracks: reference (top) and query (bottom), each with independent x-axes
- ORF boxes coloured by functional category (conjugation, T7SS, IS elements, etc.)
- Genomic island bands (auto-generated or from file)
- Key gene labels with arrows
- Light-blue curved cubic-Bezier ribbons connecting conserved blocks
- Semi-transparent blue connector for the 1:1 aligned region
- Legend and summary box on the right

### Step 4 — standalone ribbon plot

- Clean two-genome diagram with only the curved light-blue ribbons
- No ORF tracks — focused on synteny structure

---

## Output Files

| File | Format | Description |
|------|--------|-------------|
| `plasmid_blastn_identity_windows*.tsv` | TSV (774 rows) | Per-window blastn results: query pos, ref pos, %id, strand, bitscore |
| `plasmid_identity_plot*.svg` | SVG | %identity vs query position (continuous line, 60–100% y-axis) |
| `plasmid_high_identity_regions*.tsv` | TSV (44 rows) | Merged >90% conserved blocks with coordinates in both genomes |
| `plasmid_comparison_map*.svg` | SVG | Dual-axis synteny map with ORFs, islands, key genes, and curved ribbons |
| `plasmid_ribbon_synteny*.svg` | SVG | Standalone light-blue curved-ribbon synteny diagram |

All `*` placeholders are replaced by the `--suffix` value (if provided).

---

## Architecture

```
scripts/
├── run_synteny_pipeline.py    # Master orchestrator (CLI entry point)
├── _pipeline_paths.py         # Shared path resolution (env-var overrides)
├── blastn_identity_windows.py # Step 1: per-window blastn
├── high_identity_regions.py   # Step 2: merge conserved blocks
├── plasmid_comparison_map.py  # Step 3: dual-axis synteny map + ribbons
└── ribbon_synteny_plot.py     # Step 4: standalone ribbon plot
```

Each step script can also be run independently using environment variables:

```bash
SYNTENY_REF_FASTA=ref.fasta \
SYNTENY_QRY_FASTA=qry.fasta \
SYNTENY_REF_GFF=ref.gff3 \
SYNTENY_QRY_GFF=qry.gff3 \
SYNTENY_SUFFIX=standalone \
python scripts/plasmid_comparison_map.py
```

---

## Environment Variables

For advanced / scripted use, the following environment variables override defaults:

| Variable | Purpose | Default |
|----------|---------|---------|
| `SYNTENY_REF_FASTA` | Reference FASTA path | `ref_annotation/JP-H-1_plasmid.fasta` |
| `SYNTENY_REF_GFF` | Reference GFF path | `ref_annotation/JP-H-1_plasmid.gff3` |
| `SYNTENY_QRY_FASTA` | Query FASTA path | `Final_annotation/IDR2500080001-01-01_plasmid.fasta` |
| `SYNTENY_QRY_GFF` | Query GFF path | `Final_annotation/IDR2500080001-01-01_plasmid.gff3` |
| `SYNTENY_SUFFIX` | Output filename suffix | (none) |
| `SYNTENY_ISLANDS` | Island mode: `auto`, `none`, or file path | `auto` |