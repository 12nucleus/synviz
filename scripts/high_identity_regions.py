#!/usr/bin/env python3
"""
Identify contiguous regions where the per-window blastn % identity is > 90%,
and report each region's start-end in BOTH genomes:
  - query coordinates (rep-oriented)
  - reference coordinates, original and rep-oriented

Reads blastn_identity_windows.tsv and writes high_identity_regions.tsv.

Usage:
    conda run -n TB_plasmid python3 scripts/high_identity_regions.py
"""

import csv
from pathlib import Path

import _pipeline_paths as _pp
BASE    = _pp.BASE
QRY_LEN = _pp.QRY_LEN
TSV     = _pp.out_path("blastn_identity_windows", ".tsv")
OUT     = _pp.out_path("high_identity_regions", ".tsv")
THRESH  = _pp.IDENTITY_THRESHOLD

rows = []
with open(TSV) as f:
    for d in csv.DictReader(f, delimiter='\t'):
        rows.append({
            'qstart': int(d['query_start']),
            'qend':   int(d['query_end']),
            'pident': float(d['pident']) if d['pident'] != '' else 0.0,
            'hit':    int(d['has_hit']),
            'ros':    int(d['ref_start_orig']) if d['ref_start_orig'] else 0,
            'roe':    int(d['ref_end_orig'])   if d['ref_end_orig']   else 0,
            'rrs':    int(d['ref_start_rep'])  if d['ref_start_rep']  else 0,
            'rre':    int(d['ref_end_rep'])    if d['ref_end_rep']    else 0,
        })

# A window qualifies if it has a hit and identity > THRESH
qual = [ (r['hit'] == 1 and r['pident'] > THRESH) for r in rows ]

# Merge consecutive qualifying windows into contiguous regions
regions = []
i, n = 0, len(rows)
while i < n:
    if qual[i]:
        j = i
        while j + 1 < n and qual[j + 1]:
            j += 1
        block = rows[i:j + 1]
        q_start, q_end = block[0]['qstart'], block[-1]['qend']
        ros, roe = min(b['ros'] for b in block), max(b['roe'] for b in block)
        rrs, rre = min(b['rrs'] for b in block), max(b['rre'] for b in block)
        mean_pid = sum(b['pident'] for b in block) / len(block)
        wraps = False  # ref coords are now directly rep-oriented (no orig↔rep conversion)
        regions.append((q_start, q_end, ros, roe, rrs, rre, mean_pid, len(block), wraps))
        i = j + 1
    else:
        i += 1

with open(OUT, 'w') as f:
    f.write("region\tquery_start\tquery_end\tref_start_orig\tref_end_orig\t"
            "ref_start_rep\tref_end_rep\tmean_pident\tn_windows\twraps_cut\n")
    for idx, (qs, qe, ros, roe, rrs, rre, mp, nw, w) in enumerate(regions, 1):
        f.write(f"{idx}\t{qs}\t{qe}\t{ros}\t{roe}\t{rrs}\t{rre}\t{mp:.2f}\t{nw}\t{int(w)}\n")

print(f"Found {len(regions)} contiguous region(s) with per-window identity > {THRESH}%")
total_q = sum(r[1] - r[0] for r in regions)
qry_total = QRY_LEN or 1
print(f"Total query span covered: {total_q:,} bp "
      f"({100*total_q/qry_total:.1f}% of query)\n")
print(f"{'region':>6}  {'query':>16}  {'ref_orig':>18}  {'ref_rep':>18}  {'mean':>7}  {'win':>4}")
for idx, (qs, qe, ros, roe, rrs, rre, mp, nw, w) in enumerate(regions, 1):
    wrapnote = "  [WRAPS rep cut]" if w else ""
    print(f"{idx:>6}  {qs:>7}-{qe:<7}  {ros:>8}-{roe:<8}  {rrs:>8}-{rre:<8}  "
          f"{mp:>6.2f}%  {nw:>4}{wrapnote}")
