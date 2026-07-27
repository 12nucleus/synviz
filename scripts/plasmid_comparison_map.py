#!/usr/bin/env python3
"""
Linear genome comparison map — JP-H-1 (rep-oriented, top) vs IDR2500080001-01-01 (bottom).
Dual independent x-axes: top = reference coordinates, bottom = query coordinates.
ORFs drawn as boxes above (+) or below (−) a central line.
Genomic islands shown as coloured bands. Key genes highlighted.
A shaded connector links the 1:1 aligned blocks across the two axes.

Usage:
    conda run -n TB_plasmid python3 scripts/plasmid_comparison_map.py
"""

import os
import numpy as np
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, PathPatch
from matplotlib.lines import Line2D
from matplotlib.path import Path as MPath

# ── Paths (overridable via SYNTENY_* environment variables) ──
import _pipeline_paths as _pp
BASE    = _pp.BASE
REF_GFF = _pp.REF_GFF
QRY_GFF = _pp.QRY_GFF
OUT_SVG = _pp.out_path("plasmid_comparison_map", ".svg")

REF_LEN = 191151
QRY_LEN = 193713

# ── Alignment mapping (from minimap2, rep-oriented ref) ─
# query[11239:183108] ↔ ref[20413:185539]
QRY_ALN_START, QRY_ALN_END = 11239, 183108
REF_ALN_START, REF_ALN_END = 20413, 185539

def qry_to_ref(qpos):
    """Map query position → reference position via alignment (kept for reference)."""
    if qpos < QRY_ALN_START or qpos > QRY_ALN_END:
        return None
    ratio = (REF_ALN_END - REF_ALN_START) / (QRY_ALN_END - QRY_ALN_START)
    return REF_ALN_START + (qpos - QRY_ALN_START) * ratio

# ── Feature categorisation ─────────────────────────────
def categorise(product):
    if not product or product == "hypothetical protein":
        return "hypothetical"
    p = product.lower()
    if any(k in p for k in ['conjugal transfer','trbl','virb','tcp','virb','aaa-like domain']):
        return "conjugation"
    if any(k in p for k in ['mmpl','efflux','antiporter','integral membrane']):
        return "efflux_transport"
    if any(k in p for k in ['ecc','mycp','esp','esx','t7s','wxg100','ppe','pe domain','esx-1','type vii']):
        return "t7ss_esx"
    if any(k in p for k in ['rep protein','putative rep','parb','dnab','replicative dna helicase']):
        return "replication"
    if any(k in p for k in ['transposase','ist','is21','is256','is607','is110','is3','is200','tnp','mutator transposase']):
        return "is_element"
    if any(k in p for k in ['toxin','antitoxin','vapc','vapb','hipa','hnh endonuclease','abiei','abigii','pin domain','xre family','res domain']):
        return "ta_system"
    if any(k in p for k in ['dna ligase','helicase','methyltransferase','methylase','recombinase','integrase','nuclease','recb','hera','whia','whib']):
        return "dna_metabolism"
    return "other"

CAT_COLORS = {
    "conjugation":      "#E41A1C",
    "efflux_transport": "#FF7F00",
    "t7ss_esx":         "#984EA3",
    "replication":      "#377EB8",
    "is_element":       "#4DAF4A",
    "ta_system":        "#FFD700",
    "dna_metabolism":   "#A6761D",
    "hypothetical":     "#D9D9D9",
    "other":            "#B3B3B3",
}

CAT_LABELS = {
    "conjugation":      "Conjugation (TrbL/VirB4/TcpC)",
    "efflux_transport": "Efflux / Transport (MmpL/antiporter)",
    "t7ss_esx":         "T7SS / ESX-1 secretion",
    "replication":      "Replication / Maintenance",
    "is_element":       "IS Elements / Transposases",
    "ta_system":        "Toxin-Antitoxin systems",
    "dna_metabolism":   "DNA metabolism",
    "hypothetical":     "Hypothetical",
    "other":            "Other",
}

# ── GFF parser ─────────────────────────────────────────
def parse_gff(path):
    feats = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.strip().split('\t')
            if len(p) < 9 or p[2] != "CDS":
                continue
            start, end = int(p[3]), int(p[4])
            strand = p[6]
            attrs = {}
            for a in p[8].split(';'):
                if '=' in a:
                    k, v = a.split('=', 1)
                    attrs[k] = v
            product = attrs.get('Name', attrs.get('product', ''))
            gene    = attrs.get('gene', '')
            locus   = attrs.get('locus_tag', '')
            if gene and gene != product:
                label = gene
            elif product and product != "hypothetical protein":
                label = product.split('%2C')[0].split(',')[0].strip()
                if len(label) > 30:
                    label = label[:28] + '..'
            else:
                label = locus
            cat = categorise(product)
            feats.append(dict(start=start-1, end=end, strand=strand,
                              label=label, product=product, category=cat,
                              locus_tag=locus))
    return feats

# ── Genomic islands (re-oriented reference coordinates) ─
REF_ISLANDS = [
    (0, 10000, "Rep/IS110/espK", "#377EB8"),
    (35000, 50000, "RecB/DNA mod/IS607", "#A6761D"),
    (60000, 65000, "dnaB/ParB", "#377EB8"),
    (70000, 89000, "T7SS/ESX-1 cluster", "#984EA3"),
    (90000, 100000, "TA systems", "#FFD700"),
    (169000, 177000, "Conjugation cluster", "#E41A1C"),
    (187000, 190000, "espK/IS110", "#4DAF4A"),
]

# Query islands — drawn directly in QUERY coordinates (no remapping)
QRY_ISLANDS = [
    (0, 5000, "Rep", "#377EB8"),
    (21000, 24000, "IS21 (istAB)", "#4DAF4A"),
    (60000, 65000, "MmpL13/Antiporter\n★ unique to query", "#FF7F00"),
    (77000, 83000, "dnaB/ParB/Int", "#377EB8"),
    (91000, 110000, "T7SS/ESX-1 cluster", "#984EA3"),
    (110000, 116000, "TA systems", "#FFD700"),
    (161000, 172000, "Conjugation cluster", "#E41A1C"),
    (182000, 187000, "IS21 (istAB)", "#4DAF4A"),
]

# ── Key genes (re-oriented reference coordinates) ──────
REF_KEY = [
    (2, 970, "Rep", "#377EB8"),
    (169818, 172037, "TrbL/VirB6", "#E41A1C"),
    (175295, 176263, "tcpC", "#E41A1C"),
    (172037, 174481, "VirB4-like AAA", "#E41A1C"),
    (80573, 84751, "FtsK/SpoIIIE", "#E41A1C"),
    (70142, 72106, "eccE", "#984EA3"),
    (72106, 73542, "mycP", "#984EA3"),
    (73539, 75098, "eccD", "#984EA3"),
    (84881, 86464, "eccB", "#984EA3"),
    (86479, 88329, "eccA1", "#984EA3"),
    (47050, 48510, "tnpB", "#4DAF4A"),
    (60326, 62464, "dnaB", "#377EB8"),
    (63505, 64368, "parB", "#377EB8"),
]

# Query key genes — drawn directly in QUERY coordinates
QRY_KEY = [
    (130, 921, "rep", "#377EB8"),
    (182570, 183718, "rep", "#377EB8"),
    (165918, 168137, "TrbL/VirB6", "#E41A1C"),
    (171396, 172364, "tcpC", "#E41A1C"),
    (168137, 170581, "VirB4-like AAA", "#E41A1C"),
    (101659, 105837, "FtsK/SpoIIIE", "#E41A1C"),
    (60893, 63154, "MmpL13", "#FF7F00"),
    (63376, 64620, "Ion antiporter", "#FF7F00"),
    (91215, 93182, "eccE", "#984EA3"),
    (93182, 94618, "mycP", "#984EA3"),
    (94615, 96174, "eccD", "#984EA3"),
    (105918, 107501, "eccB", "#984EA3"),
    (107516, 109366, "eccA1", "#984EA3"),
    (21223, 23504, "istAB", "#4DAF4A"),
    (183926, 186207, "istAB", "#4DAF4A"),
    (77124, 79247, "dnaB", "#377EB8"),
    (79674, 80537, "parB", "#377EB8"),
]

# ── Genomic island helpers ─────────────────────────────

def auto_islands(feats, track_len, min_gap=2000):
    """Auto-generate genomic island bands by grouping adjacent CDS of the
    same functional category.  Categories 'hypothetical' and 'other' are
    skipped.  Groups separated by ≤ *min_gap* bp are merged."""
    if not feats:
        return []
    sorted_feats = sorted(feats, key=lambda f: f['start'])
    groups = []
    cur_cat, cur_start, cur_end = None, None, None
    for f in sorted_feats:
        cat = f['category']
        if cat in ('hypothetical', 'other'):
            continue
        s, e = f['start'], f['end']
        if cur_cat == cat and s - cur_end <= min_gap:
            cur_end = max(cur_end, e)
        else:
            if cur_cat is not None:
                groups.append((cur_start, cur_end, cur_cat))
            cur_cat, cur_start, cur_end = cat, s, e
    if cur_cat is not None:
        groups.append((cur_start, cur_end, cur_cat))
    # Convert to (start, end, label, color) tuples
    islands = []
    for gs, ge, cat in groups:
        label = CAT_LABELS.get(cat, cat)
        color = CAT_COLORS.get(cat, "#B3B3B3")
        islands.append((gs, ge, label, color))
    return islands


def parse_island_file(path):
    """Read genomic island definitions from a TSV file.
    Format:  start<TAB>end<TAB>label[<TAB>color]
    Lines starting with '#' are comments.  If color is omitted it is
    auto-assigned from the label (matched against CAT_LABELS)."""
    islands = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            start, end, label = int(parts[0]), int(parts[1]), parts[2]
            if len(parts) >= 4:
                color = parts[3]
            else:
                # auto-assign colour from category
                color = "#B3B3B3"
                for cat, clabel in CAT_LABELS.items():
                    if clabel.lower() in label.lower() or cat.lower() in label.lower():
                        color = CAT_COLORS.get(cat, color)
                        break
            islands.append((start, end, label, color))
    return islands


# ── Main ───────────────────────────────────────────────
def main():
    print("Building linear plasmid comparison map (dual-axis)...")
    have_ref_gff = REF_GFF.exists()
    have_qry_gff = QRY_GFF.exists()
    ref_feats = parse_gff(REF_GFF) if have_ref_gff else []
    qry_feats = parse_gff(QRY_GFF) if have_qry_gff else []
    if have_ref_gff or have_qry_gff:
        print(f"  JP-H-1: {len(ref_feats)} CDS  |  IDR2500080001: {len(qry_feats)} CDS")
    else:
        print("  No GFF annotations — genome bars + ribbons only (no ORF tracks)")

    # ── Figure: main axes on left, legend on right ──
    # Top track = reference coordinates; bottom track = query coordinates (dual axis)
    fig = plt.figure(figsize=(28, 10), facecolor='white')
    gs = fig.add_gridspec(2, 1, left=0.05, right=0.78, top=0.94, bottom=0.12,
                          hspace=0.30, height_ratios=[1, 1])
    ax_top = fig.add_subplot(gs[0])
    ax_bot = fig.add_subplot(gs[1])   # independent x-axis (query coords)

    # ── Helper: draw one track (in its own coordinate space) ──
    def draw_track(ax, feats, islands, key_genes, title, track_len):
        """Draw a linear genome track with ORFs, islands, and key genes.
        Coordinates are the track's own (ref for top, query for bottom),
        so every CDS is shown regardless of alignment coverage."""
        ax.set_xlim(-track_len * 0.015, track_len * 1.015)
        ax.set_ylim(-1.6, 1.6)

        # ── Genomic island bands ──
        for istart, iend, ilabel, icolor in islands:
            rect = Rectangle((istart, -1.55), iend - istart, 3.1,
                             facecolor=icolor, alpha=0.07, edgecolor='none', zorder=0)
            ax.add_patch(rect)
            mid = (istart + iend) / 2
            ax.text(mid, 1.48, ilabel, fontsize=6.5, color=icolor,
                    fontweight='bold', ha='center', va='top',
                    bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                              edgecolor=icolor, alpha=0.6, lw=0.4), zorder=10)

        # ── Central axis line ──
        ax.axhline(y=0, color='#333333', linewidth=0.6, zorder=1)

        # ── ORF boxes (drawn directly in track coordinates) ──
        orf_h = 0.32
        for f in feats:
            s, e = f['start'], f['end']
            strand = f['strand']
            cat = f['category']
            color = CAT_COLORS.get(cat, "#B3B3B3")

            if strand == '+':
                y = 0.04
                h = orf_h
            else:
                y = -orf_h - 0.04
                h = orf_h

            rect = Rectangle((s, y), e - s, h,
                             facecolor=color, edgecolor=color,
                             linewidth=0.1, alpha=0.85, zorder=2)
            ax.add_patch(rect)

        # ── Key gene highlights ──
        for idx, (ks, ke, klabel, kcolor) in enumerate(key_genes):
            # Larger highlight box
            rect = Rectangle((ks, -orf_h - 0.12), ke - ks, orf_h * 2 + 0.24,
                             facecolor=kcolor, alpha=0.18, edgecolor=kcolor,
                             linewidth=1.0, zorder=5)
            ax.add_patch(rect)
            # Label — alternate above/below
            mid = (ks + ke) / 2
            label_y = 1.25 if idx % 2 == 0 else -1.25
            va = 'bottom' if label_y > 0 else 'top'
            ax.annotate(klabel, xy=(mid, 0), xytext=(mid, label_y),
                        fontsize=6, color=kcolor, fontweight='bold',
                        ha='center', va=va,
                        arrowprops=dict(arrowstyle='-', color=kcolor, lw=0.4, alpha=0.4),
                        zorder=11)

        # ── Title (inside plot area, top-left) ──
        ax.text(track_len * 0.003, 1.52, title, fontsize=10, fontweight='bold',
                color='#222', va='top', ha='left', zorder=20)

        # ── Strand labels ──
        ax.text(-track_len * 0.01, orf_h / 2 + 0.04, '+', fontsize=8, color='#888',
                ha='right', va='center')
        ax.text(-track_len * 0.01, -orf_h / 2 - 0.04, '−', fontsize=8, color='#888',
                ha='right', va='center')

        # ── Styling ──
        ax.set_yticks([])
        for spine in ['top', 'right', 'left']:
            ax.spines[spine].set_visible(False)
        ax.spines['bottom'].set_visible(False)

    # ── Resolve genomic islands ──────────────────────────
    islands_mode = os.environ.get("SYNTENY_ISLANDS", "auto" if (have_ref_gff or have_qry_gff) else "none")
    if islands_mode == "none":
        ref_islands, qry_islands = [], []
    elif islands_mode == "auto":
        ref_islands = auto_islands(ref_feats, REF_LEN) if have_ref_gff else []
        qry_islands = auto_islands(qry_feats, QRY_LEN) if have_qry_gff else []
    elif os.path.isfile(islands_mode):
        all_islands = parse_island_file(islands_mode)
        ref_islands = qry_islands = all_islands
    else:
        print(f"  ⚠ Unknown --islands value '{islands_mode}' — falling back to 'auto'")
        ref_islands = auto_islands(ref_feats, REF_LEN) if have_ref_gff else []
        qry_islands = auto_islands(qry_feats, QRY_LEN) if have_qry_gff else []

    # ── Draw tracks (each in its own coordinate space) ──
    if have_ref_gff:
        draw_track(ax_top, ref_feats, ref_islands, REF_KEY,
                   "JP-H-1  (M. avium subsp. hominissuis plasmid p1-JPH1)  —  Reference  (191,151 bp)", REF_LEN)
    else:
        draw_track(ax_top, [], ref_islands, [],
                   "Reference  (191,151 bp)  [FASTA only — no annotations]", REF_LEN)
        ax_top.text(REF_LEN / 2, 0, "No GFF annotations provided", ha='center', va='center',
                    fontsize=12, color='#999', style='italic', zorder=20)

    if have_qry_gff:
        draw_track(ax_bot, qry_feats, qry_islands, QRY_KEY,
                   "IDR2500080001-01-01  (M. tuberculosis clinical isolate)  —  Query  (193,713 bp)", QRY_LEN)
    else:
        draw_track(ax_bot, [], qry_islands, [],
                   "Query  (193,713 bp)  [FASTA only — no annotations]", QRY_LEN)
        ax_bot.text(QRY_LEN / 2, 0, "No GFF annotations provided", ha='center', va='center',
                    fontsize=12, color='#999', style='italic', zorder=20)

    # ── Aligned-block shading (in each track's own coordinates) ──
    ax_top.axvspan(REF_ALN_START, REF_ALN_END, alpha=0.05, color='#377EB8', zorder=0)
    ax_bot.axvspan(QRY_ALN_START, QRY_ALN_END, alpha=0.05, color='#377EB8', zorder=0)

    # ── Connector ribbon: link the 1:1 aligned blocks across the two axes ──
    overlay = fig.add_axes([0, 0, 1, 1], frameon=False)
    overlay.set_xlim(0, 1); overlay.set_ylim(0, 1); overlay.axis('off')

    def _to_fig(ax, x, y):
        disp = ax.transData.transform((x, y))
        return tuple(fig.transFigure.inverted().transform(disp))

    A = _to_fig(ax_top, REF_ALN_START, -1.6)
    B = _to_fig(ax_top, REF_ALN_END,   -1.6)
    C = _to_fig(ax_bot, QRY_ALN_END,    1.6)
    D = _to_fig(ax_bot, QRY_ALN_START,  1.6)
    overlay.add_patch(Polygon([A, B, C, D], closed=True,
                              facecolor=(0.2157, 0.4941, 0.7216, 0.07),
                              edgecolor=(0.2157, 0.4941, 0.7216, 0.4),
                              linewidth=0.6, zorder=0))
    mx = (A[0] + B[0] + C[0] + D[0]) / 4
    my = (A[1] + B[1] + C[1] + D[1]) / 4
    overlay.text(mx, my, "1:1 aligned block", fontsize=7, color='#377EB8',
                 ha='center', va='center', rotation=90, fontweight='bold',
                 transform=overlay.transData, zorder=1)

    # ── Conserved-block ribbons (light blue, curved) linking the two tracks ──
    RIBBON_TSV = _pp.out_path("plasmid_high_identity_regions", ".tsv")
    LIGHT_BLUE = (0.61, 0.78, 0.95, 0.40)   # light blue, semi-transparent
    RIBBON_EDGE = (0.30, 0.55, 0.80, 0.55)

    def _ribbon_path(tL, tR, bL, bR):
        """Closed polygon with cubic-Bezier (S-curve) left/right edges (fig coords)."""
        ym = 0.5 * (tL[1] + bL[1])
        verts = [tL, tR,
                 (tR[0], ym), (bR[0], ym), bR,
                 bL,
                 (bL[0], ym), (tL[0], ym), tL]
        codes = [MPath.MOVETO, MPath.LINETO,
                 MPath.CURVE4, MPath.CURVE4, MPath.CURVE4,
                 MPath.LINETO,
                 MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]
        return MPath(verts, codes)

    if RIBBON_TSV.exists():
        import csv as _csv
        n_rib = 0
        with open(RIBBON_TSV) as _f:
            for _r in _csv.DictReader(_f, delimiter='\t'):
                # top track = reference in REP-oriented coords (matches this map)
                ra, rb = int(_r['ref_start_rep']), int(_r['ref_end_rep'])
                # bottom track = query in rep-oriented coords
                qa, qb = int(_r['query_start']), int(_r['query_end'])
                # attach at bottom edge of top track / top edge of bottom track
                tL = _to_fig(ax_top, ra, -1.6)
                tR = _to_fig(ax_top, rb, -1.6)
                bL = _to_fig(ax_bot, qa,  1.6)
                bR = _to_fig(ax_bot, qb,  1.6)
                overlay.add_patch(PathPatch(_ribbon_path(tL, tR, bL, bR),
                                            facecolor=LIGHT_BLUE,
                                            edgecolor=RIBBON_EDGE,
                                            linewidth=0.3, zorder=1))
                n_rib += 1
        print(f"  ✓ Overlaid {n_rib} light-blue conserved-block ribbons")

    # ── X-axis ticks (independent for each track) ──
    for ax, L, name in [(ax_top, REF_LEN, 'reference'),
                        (ax_bot, QRY_LEN, 'query')]:
        ticks = np.arange(0, L + 1, 10000)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{int(t/1000)}kb" if t > 0 else "0" for t in ticks],
                           fontsize=7, rotation=90)
        ax.spines['bottom'].set_visible(True)
        ax.spines['bottom'].set_color('#999')
        ax.tick_params(axis='x', colors='#666', length=3, pad=2)
        ax.set_xlabel(f"Genome position ({name} coordinates)", fontsize=9,
                      color='#555', labelpad=8)

    # ── Legend (separate area on the right) ──
    leg_ax = fig.add_axes([0.80, 0.30, 0.18, 0.45])
    leg_ax.axis('off')
    leg_ax.text(0, 0.99, "Feature Categories", fontsize=11, fontweight='bold',
                transform=leg_ax.transAxes, va='top')

    handles = []
    for cat in ['conjugation','efflux_transport','t7ss_esx','replication',
                'is_element','ta_system','dna_metabolism','hypothetical']:
        handles.append(Line2D([0],[0], marker='s', color='w',
                              markerfacecolor=CAT_COLORS[cat], markersize=10,
                              label=CAT_LABELS[cat]))

    leg_ax.legend(handles=handles, loc='upper left', frameon=True,
                  fontsize=8, handletextpad=1.2, bbox_to_anchor=(0, 0.58))

    # Summary box
    if have_ref_gff and have_qry_gff:
        summary = (
            "Plasmid Comparison\n"
            "──────────────────\n"
            f"JP-H-1: {REF_LEN:,} bp, 66.3% GC\n"
            f"IDR2500080001: {QRY_LEN:,} bp, 65.9% GC\n"
            f"Alignment: 99.96% identity\n"
            f"over 87.4% of reference\n\n"
            "★ MmpL13 + antiporter\n"
            "  unique to query plasmid\n"
            "  → potential INH efflux\n\n"
            "Blue = 1:1 aligned block\n"
            "(connector links axes)\n\n"
            "Light-blue ribbons =\n"
            "conserved >90% blocks"
        )
    else:
        summary = (
            "Plasmid Comparison\n"
            "──────────────────\n"
            f"Ref: {REF_LEN:,} bp\n"
            f"Qry: {QRY_LEN:,} bp\n\n"
            "⚠ FASTA-only mode\n"
            "No GFF annotations\n"
            "provided — ORF tracks\n"
            "are not drawn.\n\n"
            "Blue = 1:1 aligned block\n"
            "(connector links axes)\n\n"
            "Light-blue ribbons =\n"
            "conserved >90% blocks"
        )
    leg_ax.text(0, 0.02, summary, fontsize=8.5, fontfamily='monospace',
                transform=leg_ax.transAxes, va='bottom',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow',
                          edgecolor='#ccc'))

    fig.savefig(str(OUT_SVG), dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  ✓ Saved: {OUT_SVG}")

if __name__ == "__main__":
    main()
