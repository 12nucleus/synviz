#!/usr/bin/env python3
"""
Shared path resolution for the synteny pipeline.
All step scripts import this module to resolve input/output paths.

Command-line overrides are communicated via environment variables
set by the master orchestrator (run_synteny_pipeline.py):

  SYNTENY_REF_FASTA    — reference FASTA  (default: test_files/ref.fasta)
  SYNTENY_REF_GFF      — reference GFF    (default: test_files/ref.gff3)
  SYNTENY_QRY_FASTA    — query FASTA      (default: test_files/qry.fasta)
  SYNTENY_QRY_GFF      — query GFF        (default: test_files/qry.gff3)
  SYNTENY_SUFFIX       — output filename suffix (default: "" → no suffix)

When SYNTENY_SUFFIX is set (e.g. "cmp1"), outputs are named:
  blastn_identity_windows_cmp1.tsv
  identity_plot_cmp1.svg
  …
"""

import os
from pathlib import Path

# ── Project root ────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent

# ── Output-name suffix ──────────────────────────────────────────────────────
_SUFFIX = os.environ.get("SYNTENY_SUFFIX", "")
if _SUFFIX and not _SUFFIX.startswith("_"):
    _SUFFIX = "_" + _SUFFIX


def out_path(basename, ext=".tsv"):
    """Return a path in the current working directory with the optional suffix inserted before ext."""
    return Path.cwd() / f"{basename}{_SUFFIX}{ext}"


def base_name(basename, ext=".tsv"):
    """Return just the filename string (no directory)."""
    return f"{basename}{_SUFFIX}{ext}"


# ── Input files ─────────────────────────────────────────────────────────────
REF_FASTA = Path(os.environ.get(
    "SYNTENY_REF_FASTA",
    str(BASE / "test_files/ref.fasta")))
REF_GFF = Path(os.environ.get(
    "SYNTENY_REF_GFF",
    str(BASE / "test_files/ref.gff3")))
QRY_FASTA = Path(os.environ.get(
    "SYNTENY_QRY_FASTA",
    str(BASE / "test_files/qry.fasta")))
QRY_GFF = Path(os.environ.get(
    "SYNTENY_QRY_GFF",
    str(BASE / "test_files/qry.gff3")))

# ── Configurable parameters (overridable via environment variables) ──────────
IDENTITY_THRESHOLD = float(os.environ.get("SYNTENY_IDENTITY_THRESHOLD", "90.0"))

# ── Constants ───────────────────────────────────────────────────────────────
REF_LEN = 191151
QRY_LEN = 193713
REF_CUT = 177304   # 1-based rep-protein start for original → rep-oriented rotation
