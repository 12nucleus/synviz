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
import shutil
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

# ── BLAST executables (resolved from PATH) ────────────────────────────────────
# Use shutil.which() for resolution; if not found, fall back to bare command
# name so the subprocess runtime PATH (e.g. active conda env) is still used.
BLASTN      = shutil.which("blastn")      or "blastn"
MAKEBLASTDB = shutil.which("makeblastdb") or "makeblastdb"


# ── Configurable parameters (overridable via environment variables) ──────────
IDENTITY_THRESHOLD = float(os.environ.get("SYNTENY_IDENTITY_THRESHOLD", "90.0"))
MIN_REGION_LENGTH  = int(os.environ.get("SYNTENY_MIN_REGION_LENGTH", "1000"))
MAX_REF_GAP        = int(os.environ.get("SYNTENY_MAX_REF_GAP", "2000"))

# ── Genome lengths (overridable via environment variables) ──────────────────
def _read_fasta_length(path):
    """Return the total length (bp) of a (multi-line) FASTA file."""
    total = 0
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                total += len(line.strip())
    return total

REF_LEN = int(os.environ.get("SYNTENY_REF_LEN", "0")) or _read_fasta_length(REF_FASTA)
QRY_LEN = int(os.environ.get("SYNTENY_QRY_LEN", "0")) or _read_fasta_length(QRY_FASTA)
