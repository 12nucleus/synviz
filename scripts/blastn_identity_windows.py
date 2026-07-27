#!/usr/bin/env python3
"""
Per-window percent nucleotide identity between two sequences, computed with
blastn over a 500 bp sliding window (250 bp step) on the QUERY.

For each window the best blastn HSP vs the reference (in its native original
orientation) is recorded. The result table stores BOTH the query position and
the reference position of the match (original + rep-oriented), plus the %
identity.

Outputs:
  - blastn_identity_windows.tsv   (per-window table: query + ref positions, %id)
  - identity_plot.svg             (% identity vs query position, continuous)

Usage:
    conda run -n TB_plasmid python3 scripts/blastn_identity_windows.py
"""

import subprocess, tempfile, os, shutil, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Paths (overridable via SYNTENY_* environment variables) ──
import _pipeline_paths as _pp
BASE      = _pp.BASE
REF_FASTA = _pp.REF_FASTA
QRY_FASTA = _pp.QRY_FASTA
OUT_TSV   = _pp.out_path("blastn_identity_windows", ".tsv")
OUT_SVG   = _pp.out_path("identity_plot", ".svg")

BLASTN      = "/Users/Pascal/anaconda3/envs/TB/bin/blastn"
MAKEBLASTDB = "/Users/Pascal/anaconda3/envs/TB/bin/makeblastdb"

WINDOW = 500
STEP   = 250

# ── FASTA I/O ──────────────────────────────────────────
def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                continue
            seq.append(line.strip())
    return ''.join(seq).upper()

def write_fasta(path, header, seq):
    with open(path, 'w') as f:
        f.write(header + '\n')
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + '\n')

# ── Main ───────────────────────────────────────────────
def main():
    print("=" * 64)
    print("blastn per-window % identity (500 bp window, 250 bp step)")
    print("=" * 64)

    if not os.path.exists(BLASTN) or not os.path.exists(MAKEBLASTDB):
        print("ERROR: blastn / makeblastdb not found at expected path")
        sys.exit(1)

    print("[1/5] Reading sequences...")
    qseq = read_fasta(QRY_FASTA)
    rseq = read_fasta(REF_FASTA)
    QRY_LEN = len(qseq)
    print(f"  query : {QRY_LEN} bp (rep-oriented)")
    print(f"  ref   : {len(rseq)} bp (original orientation)")

    # Ref is used in its native (original) orientation for the blast, which
    # keeps the query->ref mapping collinear. The output table records the ref
    # coordinates directly so positions can be compared with the comparison map.
    tmp = tempfile.mkdtemp(prefix="blastn_win_")
    try:
        # Reference DB (rotated, rep-oriented)
        ref_db_fasta = os.path.join(tmp, "ref_orig.fasta")
        write_fasta(ref_db_fasta, ">reference", rseq)
        subprocess.run([MAKEBLASTDB, "-in", ref_db_fasta, "-dbtype", "nucl",
                        "-out", os.path.join(tmp, "refdb")],
                       check=True, capture_output=True, text=True)

        # Query windows multi-FASTA (500 bp, step 250)
        qwin_fasta = os.path.join(tmp, "query_windows.fasta")
        windows = []   # (win_id, qstart_1based, qend_1based)
        with open(qwin_fasta, 'w') as f:
            idx = 0
            start = 0
            while start + WINDOW <= QRY_LEN:
                end = start + WINDOW
                win = qseq[start:end]
                wid = f"win_{idx:04d}"
                f.write(f">{wid} qstart={start+1} qend={end}\n")
                for i in range(0, len(win), 80):
                    f.write(win[i:i+80] + '\n')
                windows.append((wid, start + 1, end))
                idx += 1
                start += STEP
            # trailing partial window
            if start < QRY_LEN:
                end = QRY_LEN
                win = qseq[start:end]
                wid = f"win_{idx:04d}"
                f.write(f">{wid} qstart={start+1} qend={end}\n")
                for i in range(0, len(win), 80):
                    f.write(win[i:i+80] + '\n')
                windows.append((wid, start + 1, end))
                idx += 1
        print(f"  windows: {len(windows)} (500 bp, step {STEP})")

        # blastn: query windows -> ref
        print("[2/5] Running blastn (query windows -> ref)...")
        db = os.path.join(tmp, "refdb")
        outfmt = "6 qseqid sseqid sstart send qstart qend pident length bitscore sstrand"
        res = subprocess.run(
            [BLASTN, "-query", qwin_fasta, "-db", db, "-outfmt", outfmt,
             "-task", "blastn", "-evalue", "1e-5", "-max_target_seqs", "5"],
            check=True, capture_output=True, text=True)
        lines = [l for l in res.stdout.strip().split('\n') if l]

        # Parse: keep the best HSP (by BITSCORE) for each window. Bitscore
        # balances identity and alignment length, so a short 100% micro-match
        # loses to the longer true ortholog. Every hit is reported with its
        # ref coordinates.
        print("[3/5] Parsing blastn output...")
        best = {}        # best HSP per window (by bitscore)
        for l in lines:
            p = l.split('\t')
            wid     = p[0]
            sstart, send = int(p[2]), int(p[3])
            qstart, qend = int(p[4]), int(p[5])
            pident  = float(p[6])
            alnlen  = int(p[7])
            bits    = float(p[8])
            sstrand = p[9]
            rec = dict(sstart=sstart, send=send, qstart=qstart, qend=qend,
                       pident=pident, alnlen=alnlen, bits=bits, sstrand=sstrand)
            if wid not in best or bits > best[wid]['bits']:
                best[wid] = rec

        # Build table
        print("[4/5] Building per-window table...")
        rows = []
        for wid, qa, qb in windows:
            r = best.get(wid)
            if r is not None:
                pid = round(r['pident'], 3)
                if r['sstrand'] == 'plus':        # ref coords (rep-oriented), start<end
                    rs, re = r['sstart'], r['send']
                else:
                    rs, re = r['send'], r['sstart']
                # Both orig and rep columns are the same — the REF FASTA is
                # already rep-oriented, so blastn coordinates ARE rep-oriented.
                rows.append((wid, qa, qb, pid, rs, re,
                             rs, re,
                             r['sstrand'], r['alnlen'], round(r['bits'], 1), 1))
            else:
                # No blastn hit at all: report 0 (not NaN) for every numeric field.
                rows.append((wid, qa, qb, 0, 0, 0, 0, 0, '', 0, 0, 0))

        with open(OUT_TSV, 'w') as f:
            f.write("window\tquery_start\tquery_end\tpident\tref_start_orig\tref_end_orig\t"
                    "ref_start_rep\tref_end_rep\tref_strand\taln_length\tbitscore\thas_hit\n")
            for r in rows:
                f.write("\t".join(str(x) for x in r) + "\n")
        print(f"  ✓ Saved: {OUT_TSV}  ({len(rows)} windows)")

        # ── Plot: % identity vs query position ──
        print("[5/5] Plotting...")
        xs = [(r[1] + r[2]) / 2 / 1000.0 for r in rows]   # query mid (kb)
        ys = [r[3] if r[11] == 1 else 0.0 for r in rows]   # continuous line; no-hit windows drop to 0
        fig, ax = plt.subplots(figsize=(22, 6))
        ax.plot(xs, ys, color='#2166AC', linewidth=0.7, zorder=3,
                label='best blastn hit per window')
        ax.fill_between(xs, ys, 0, color='#2166AC', alpha=0.12, zorder=2)
        ax.axhline(98, color='#E41A1C', ls='--', lw=0.5, alpha=0.4)
        ax.axhline(90, color='#FF7F00', ls='--', lw=0.5, alpha=0.4)
        ax.axhline(85, color='#999999', ls=':', lw=0.5, alpha=0.4)
        ax.set_ylim(60, 100.5)
        ax.set_xlim(0, QRY_LEN / 1000.0)
        ax.set_xlabel("Position on query — kb", fontsize=10)
        ax.set_ylabel("% Identity (blastn)", fontsize=10)
        ax.set_title("Per-window % nucleotide identity  (500 bp window, 250 bp step, blastn)",
                     fontsize=12, fontweight='bold')
        # Stats
        hits = [r[3] for r in rows if r[11] == 1]
        mean_id = float(np.mean(hits)) if hits else 0.0
        n_hit = sum(r[11] for r in rows)
        n_nohit = len(rows) - n_hit
        ax.legend(loc='lower left', fontsize=8, framealpha=0.9)
        ax.text(0.98, 0.05,
                f"Mean identity (all hits): {mean_id:.2f}%   "
                f"Windows with hit: {n_hit}/{len(rows)}   "
                f"No hit: {n_nohit}/{len(rows)}",
                transform=ax.transAxes, ha='right', va='bottom', fontsize=8.5,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        ax.text(0.02, 0.95, "98%", transform=ax.transAxes, fontsize=7,
                color='#E41A1C', alpha=0.7)
        ax.text(0.02, 0.78, "90%", transform=ax.transAxes, fontsize=7,
                color='#FF7F00', alpha=0.7)
        ax.text(0.02, 0.62, "85%", transform=ax.transAxes, fontsize=7,
                color='#999999', alpha=0.7)
        ax.grid(True, axis='y', alpha=0.3, lw=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        fig.savefig(str(OUT_SVG), dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  ✓ Saved: {OUT_SVG}")
        print(f"\n  Mean identity (all hits): {mean_id:.2f}%   "
              f"Windows with hit: {n_hit}/{len(rows)}   No hit: {n_nohit}/{len(rows)}")
        print("Done!")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

if __name__ == "__main__":
    main()
