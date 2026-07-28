#!/usr/bin/env python3
"""
=============================================================================
  Synteny Plot Generation Pipeline — Master Orchestrator
=============================================================================
Runs the complete comparison workflow between two sequences in the correct
order:

  Step 1 — blastn identity windows   (500 bp sliding window %identity)
  Step 2 — high-identity regions     (merge windows into conserved blocks above threshold)
  Step 3 — comparison map            (dual-axis ORF map + curved light-blue ribbons)
  Step 4 — standalone ribbon plot    (light-blue curved ribbons only)

If a --suffix is provided (e.g. "cmp1"), ALL output files are tagged so
previous plots are preserved:

  blastn_identity_windows_cmp1.tsv     identity_plot_cmp1.svg
  high_identity_regions_cmp1.tsv
  comparison_map_cmp1.svg              ribbon_synteny_cmp1.svg

Usage:
  # Default run
  conda run -n TB_plasmid python3 scripts/run_synteny_pipeline.py

  # With a named comparison (preserves previous outputs)
  conda run -n TB_plasmid python3 scripts/run_synteny_pipeline.py --suffix cmp1

  # Custom input files
  conda run -n TB_plasmid python3 scripts/run_synteny_pipeline.py \\
      --ref-fasta path/to/ref.fasta --ref-gff path/to/ref.gff3 \\
      --qry-fasta path/to/qry.fasta --qry-gff path/to/qry.gff3 \\
      --suffix comparison_A

  # Skip slow blastn step (re-use existing windows)
  conda run -n TB_plasmid python3 scripts/run_synteny_pipeline.py --skip blastn

  # Only regenerate the comparison map
  conda run -n TB_plasmid python3 scripts/run_synteny_pipeline.py --only map

Requirements:
  - Conda environment 'TB_plasmid' (matplotlib, numpy)
  - blastn 2.16+ at  /Users/Pascal/anaconda3/envs/TB/bin/
=============================================================================
"""

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

# ── Defaults (used when CLI args are not provided) ──────────────────────────
BASE    = Path(__file__).resolve().parent.parent
SCRIPTS = Path(__file__).resolve().parent
OUT_DIR = Path.cwd()

DEFAULT_REF_FASTA = BASE / "test_files/ref.fasta"
DEFAULT_REF_GFF   = BASE / "test_files/ref.gff3"
DEFAULT_QRY_FASTA = BASE / "test_files/qry.fasta"
DEFAULT_QRY_GFF   = BASE / "test_files/qry.gff3"

BLASTN      = "/Users/Pascal/anaconda3/envs/TB/bin/blastn"
MAKEBLASTDB = "/Users/Pascal/anaconda3/envs/TB/bin/makeblastdb"
PYTHON      = sys.executable
CONDA_ENV   = "TB_plasmid"

# ── Step definitions ────────────────────────────────────────────────────────

STEPS = [
    {
        "key":     "blastn",
        "name":    "blastn (identity windows)",
        "script":  SCRIPTS / "blastn_identity_windows.py",
        "outputs": lambda sfx: [
            OUT_DIR / f"blastn_identity_windows{sfx}.tsv",
            OUT_DIR / f"identity_plot{sfx}.svg",
        ],
        "desc":    "500 bp sliding-window blastn %identity vs reference",
        "slow":    True,
    },
    {
        "key":     "regions",
        "name":    "high-identity regions",
        "script":  SCRIPTS / "high_identity_regions.py",
        "outputs": lambda sfx: [OUT_DIR / f"high_identity_regions{sfx}.tsv"],
        "depends": "blastn",
        "desc":    "merge windows into conserved blocks above threshold",
    },
    {
        "key":     "map",
        "name":    "comparison map + ribbons",
        "script":  SCRIPTS / "comparison_map.py",
        "outputs": lambda sfx: [OUT_DIR / f"comparison_map{sfx}.svg"],
        "depends": "regions",
        "desc":    "dual-axis ORF synteny map with light-blue curved ribbon overlay",
    },
    {
        "key":     "ribbon",
        "name":    "standalone ribbon plot",
        "script":  SCRIPTS / "ribbon_synteny_plot.py",
        "outputs": lambda sfx: [OUT_DIR / f"ribbon_synteny{sfx}.svg"],
        "depends": "regions",
        "desc":    "clean light-blue curved-ribbon synteny diagram",
    },
]

# ── Helpers ─────────────────────────────────────────────────────────────────

def _read_fasta_length(path):
    """Return total length (bp) of a FASTA file without loading entire sequence."""
    total = 0
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def _hr(title="", char="=", width=68):
    if title:
        print(f"\n{char * 4}  {title}  {char * max(0, width - len(title) - 6)}")
    else:
        print(char * width)


def _resolve_suffix(raw):
    """Normalise a CLI suffix: empty → '', 'cmp1' → '_cmp1'."""
    if not raw:
        return ""
    s = raw.strip().replace(" ", "_")
    return s if s.startswith("_") else "_" + s


def _output_status(outputs):
    for p in outputs:
        if p.exists():
            size_kb = p.stat().st_size / 1024
            print(f"     ✓  {p.name:50s}  {size_kb:7.1f} KB")
        else:
            print(f"     ✗  {p.name:50s}  NOT GENERATED")


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Synteny plot generation pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --run \\
      --ref-fasta ref.fasta --ref-gff ref.gff3 \\
      --qry-fasta qry.fasta --qry-gff qry.gff3 \\
      --islands-ref ref_islands.tsv --islands-qry qry_islands.tsv \\
      --suffix my_comparison

  %(prog)s --run --ref-fasta R.fa --qry-fasta Q.fa --suffix fasta_only

  %(prog)s --run --ref-fasta R.fa --qry-fasta Q.fa --skip blastn

  %(prog)s --dry-run --ref-fasta R.fa --qry-fasta Q.fa
        """,
    )

    # ── Input / output arguments ─────────────────────────────────────────
    inp = parser.add_argument_group("Input files")
    inp.add_argument("--ref-fasta", type=Path, default=None,
                     help=f"Reference FASTA  [default: {DEFAULT_REF_FASTA.name}]")
    inp.add_argument("--ref-gff",   type=Path, default=None,
                     help=f"Reference GFF3 (optional)  [default: {DEFAULT_REF_GFF.name}]")
    inp.add_argument("--qry-fasta", type=Path, default=None,
                     help=f"Query FASTA      [default: {DEFAULT_QRY_FASTA.name}]")
    inp.add_argument("--qry-gff",   type=Path, default=None,
                     help=f"Query GFF3 (optional)       [default: {DEFAULT_QRY_GFF.name}]")

    out = parser.add_argument_group("Output")
    out.add_argument("--suffix", "-s", type=str, default="",
                     help="Tag appended to all output filenames (e.g. 'cmp1' → _cmp1)")
    out.add_argument("--islands-ref", type=str, default="auto",
                     help="Reference island bands: 'auto' (from GFF), 'none', or path to a TSV file")
    out.add_argument("--islands-qry", type=str, default="auto",
                     help="Query island bands: 'auto' (from GFF), 'none', or path to a TSV file")
    out.add_argument("--identity-threshold", type=float, default=90.0,
                     help="%% identity threshold for conserved blocks  [default: 90.0]")
    out.add_argument("--min-region-length", type=int, default=1000,
                     help="Minimum conserved block length (bp) to draw a ribbon  [default: 1000]")

    # ── Step selection ───────────────────────────────────────────────────
    ctl = parser.add_argument_group("Step control")
    step_keys = [s["key"] for s in STEPS]
    ctl.add_argument("--run", action="store_true",
                     help="Execute the pipeline (required; without it, help is shown)")
    ctl.add_argument("--skip",  nargs="*", choices=step_keys,
                     help="Steps to skip")
    ctl.add_argument("--only",  nargs="*", choices=step_keys,
                     help="Only run these steps (dependencies auto-resolved)")
    ctl.add_argument("--force", action="store_true",
                     help="Re-run steps even if outputs already exist")
    ctl.add_argument("--dry-run", action="store_true",
                     help="Show plan without executing")

    args = parser.parse_args()

    # Require --run (or --dry-run) to execute; otherwise show help
    if not args.run and not args.dry_run:
        parser.print_help()
        return

    # ── Resolve inputs ───────────────────────────────────────────────────
    ref_fasta = Path(args.ref_fasta or DEFAULT_REF_FASTA).resolve()
    ref_gff   = Path(args.ref_gff   or DEFAULT_REF_GFF).resolve()
    qry_fasta = Path(args.qry_fasta or DEFAULT_QRY_FASTA).resolve()
    qry_gff   = Path(args.qry_gff   or DEFAULT_QRY_GFF).resolve()
    suffix    = _resolve_suffix(args.suffix)

    # ── Determine which steps to run ─────────────────────────────────────
    if args.only:
        selected = set(args.only)
        changed = True
        while changed:
            changed = False
            for s in STEPS:
                dep = s.get("depends")
                if s["key"] in selected and dep and dep not in selected:
                    dep_step = next(x for x in STEPS if x["key"] == dep)
                    dep_outputs = dep_step["outputs"](suffix)
                    if not args.force and all(p.exists() for p in dep_outputs):
                        continue
                    selected.add(dep)
                    changed = True
        to_run = [s for s in STEPS if s["key"] in selected]
    else:
        to_run = [s for s in STEPS if args.skip is None or s["key"] not in (args.skip or [])]

    skipped_existing = []
    actual_run = []
    for s in to_run:
        outputs = s["outputs"](suffix)
        if not args.force and all(p.exists() for p in outputs):
            skipped_existing.append(s)
        else:
            actual_run.append(s)

    # ── Print plan ───────────────────────────────────────────────────────
    print()
    _hr("SYNTENY PLOT PIPELINE")
    print(f"\n  Project  : {BASE.name}")
    print(f"  Python   : {PYTHON}")
    print(f"  Ref      : {ref_fasta}")
    print(f"  Query    : {qry_fasta}")
    if suffix:
        print(f"  Suffix   : '{suffix}'  →  outputs tagged with '{suffix}'")

    if args.dry_run:
        print("\n  [DRY RUN — no actions taken]")
        if actual_run:
            print(f"\n  Would execute ({len(actual_run)} steps):")
            for s in actual_run:
                print(f"    • {s['name']}")
        print()
        return

    if not actual_run:
        print("\n  All outputs already exist — nothing to do.")
        if skipped_existing:
            print("\n  Existing outputs:")
            for s in skipped_existing:
                _output_status(s["outputs"](suffix))
        print("\n  Use --force to re-run.\n")
        return

    print(f"\n  Steps to execute ({len(actual_run)}):")
    for s in actual_run:
        tag = "⏳" if s.get("slow") else "  "
        print(f"    {tag}  {s['key']:10s}  →  {s['desc']}")

    # ── Check prerequisites ──────────────────────────────────────────────
    _hr("PREREQUISITES", "─")
    ok = True
    for label, p, required in [("Ref FASTA", ref_fasta, True),
                                ("Ref GFF",   ref_gff,   False),
                                ("Qry FASTA", qry_fasta, True),
                                ("Qry GFF",   qry_gff,   False)]:
        exists = p.exists()
        if exists:
            tag = "✓"
        elif required:
            tag = "✗ MISSING"
            ok = False
        else:
            tag = "—"
        print(f"  {tag}  {label:15s}  {p}")
    for label, p in [("blastn", BLASTN), ("makeblastdb", MAKEBLASTDB)]:
        exists = os.path.isfile(p) and os.access(p, os.X_OK)
        tag = "✓" if exists else "✗ NOT FOUND"
        print(f"  {tag}  {label:15s}  {p}")
        if not exists:
            ok = False
    env = os.environ.get("CONDA_DEFAULT_ENV", "")
    if env != CONDA_ENV:
        print(f"  ⚠  Conda env '{CONDA_ENV}' not active (current: '{env or 'none'}')."
              f"  Run: conda activate {CONDA_ENV}")

    if not ok:
        print("\n  Cannot proceed — missing prerequisites.\n")
        sys.exit(1)
    print()

    # ── Build environment for subprocesses ───────────────────────────────
    child_env = os.environ.copy()
    child_env["SYNTENY_REF_FASTA"] = str(ref_fasta)
    child_env["SYNTENY_REF_GFF"]   = str(ref_gff)
    child_env["SYNTENY_QRY_FASTA"] = str(qry_fasta)
    child_env["SYNTENY_QRY_GFF"]   = str(qry_gff)
    child_env["SYNTENY_SUFFIX"]             = suffix.lstrip("_") if suffix else ""
    child_env["SYNTENY_ISLANDS_REF"]        = args.islands_ref
    child_env["SYNTENY_ISLANDS_QRY"]        = args.islands_qry
    child_env["SYNTENY_IDENTITY_THRESHOLD"] = str(args.identity_threshold)
    child_env["SYNTENY_MIN_REGION_LENGTH"]  = str(args.min_region_length)
    child_env["SYNTENY_REF_LEN"] = str(_read_fasta_length(ref_fasta))
    child_env["SYNTENY_QRY_LEN"] = str(_read_fasta_length(qry_fasta))

    # ── Execute ──────────────────────────────────────────────────────────
    total = len(actual_run)
    pipeline_start = time.time()
    failures = 0

    for i, step in enumerate(actual_run, 1):
        _hr(f"STEP {i}/{total}  —  {step['name']}")
        print(f"  Script : {step['script'].name}")
        print(f"  Action : {step['desc']}")
        if step.get("slow"):
            print("  ⏳ This step may take several minutes (blastn search) …")
        outputs = step["outputs"](suffix)
        print(f"  Output : {', '.join(p.name for p in outputs)}")

        t0 = time.time()
        result = subprocess.run(
            [PYTHON, str(step["script"])],
            cwd=str(OUT_DIR),
            env=child_env,
            capture_output=False,
            text=True,
        )
        elapsed = time.time() - t0

        if result.returncode == 0:
            print(f"\n  ✓  Completed in {elapsed:.0f}s")
        else:
            print(f"\n  ✗  FAILED after {elapsed:.0f}s (exit code {result.returncode})")
            failures += 1
            for s in actual_run[i:]:
                if s.get("depends") == step["key"]:
                    print(f"\n  ⚠  Skipping '{s['name']}' (dependency '{step['key']}' failed)")
                    failures += 1
            break
        _output_status(outputs)

    # ── Summary ──────────────────────────────────────────────────────────
    pipeline_elapsed = time.time() - pipeline_start
    _hr("SUMMARY")
    print(f"\n  Total time : {pipeline_elapsed:.0f}s  ({pipeline_elapsed/60:.1f} min)")
    n_ok = total - failures
    print(f"  Result     : {n_ok}/{total} steps succeeded"
          + (f", {failures} failed" if failures else ""))
    if not failures:
        print(f"\n  Output files (in {OUT_DIR}/):")
        done_keys = {s["key"] for s in actual_run}
        for s in STEPS:
            if s["key"] in done_keys:
                for p in s["outputs"](suffix):
                    if p.exists():
                        print(f"    {p.name}")
        print()
    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
