#!/usr/bin/env python3
"""
Linear genome comparison map — reference (top) vs query (bottom).
Dual independent x-axes: top = reference coordinates, bottom = query coordinates.
ORFs drawn as boxes above (+) or below (−) a central line.
Genomic islands shown as coloured bands. Key genes highlighted.
A shaded connector links the 1:1 aligned blocks across the two axes.

Usage:
    conda run -n TB_plasmid python3 scripts/comparison_map.py
"""

import os
import numpy as np
from pathlib import Path

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
OUT_SVG = _pp.out_path("comparison_map", ".svg")
ID_THR  = int(_pp.IDENTITY_THRESHOLD) if _pp.IDENTITY_THRESHOLD == int(_pp.IDENTITY_THRESHOLD) else _pp.IDENTITY_THRESHOLD

# Genome lengths — derived from the shared pipeline constants.
REF_LEN = _pp.REF_LEN
QRY_LEN = _pp.QRY_LEN

# ── Optional 1:1 aligned-block coords ──
# When set via the SYNTENY_REF_ALN_* / SYNTENY_QRY_ALN_* environment variables,
# a shaded connector is drawn between the two tracks. Leave unset to skip it.
def _aln_env(name):
    v = os.environ.get(name)
    return int(v) if v and v.strip() else None

REF_ALN_START = _aln_env("SYNTENY_REF_ALN_START")
REF_ALN_END   = _aln_env("SYNTENY_REF_ALN_END")
QRY_ALN_START = _aln_env("SYNTENY_QRY_ALN_START")
QRY_ALN_END   = _aln_env("SYNTENY_QRY_ALN_END")
HAS_ALN = all(x is not None for x in
              (REF_ALN_START, REF_ALN_END, QRY_ALN_START, QRY_ALN_END))

# ── Feature categorisation ─────────────────────────────
def categorise(product):
    if not product or product == "hypothetical protein":
        return "hypothetical"
    p = product.lower()

    # Structural / mobile elements
    if any(k in p for k in ['conjugal transfer','trbl','virb','tcp','aaa-like domain']):
        return "conjugation"
    if any(k in p for k in ['transposase','ist','is21','is256','is607','is110','is3','is200','tnp','mutator transposase']):
        return "is_element"
    if any(k in p for k in ['integrase', 'recombinase', 'resolvase', 'inversion']):
        return "mobile_element"
    if any(k in p for k in ['phage', 'prophage', 'capsid', 'tail fiber', 'lysin']):
        return "phage"

    # Secretion systems
    if any(k in p for k in ['ecc','mycp','esp','esx','t7s','wxg100','ppe','pe domain','esx-1','type vii']):
        return "t7ss_esx"
    if any(k in p for k in ['type iv', 'type i secretion', 'type ii secretion', 'type iii secretion',
                            'type vi secretion', 'tss', 'hcp', 'vgrg', 'dot', 'icm']):
        return "secretion"
    if any(k in p for k in ['flagell', 'pilin', 'fimbri', 'pilus']):
        return "motility_adhesion"

    # Transport / efflux
    if any(k in p for k in ['mmpl','efflux','antiporter','integral membrane']):
        return "efflux_transport"
    if any(k in p for k in ['abc transporter', 'permease', 'multidrug', 'mdr', 'transport']):
        return "transport"

    # Replication / maintenance
    if any(k in p for k in ['rep protein','putative rep','parb','dnab','replicative dna helicase']):
        return "replication"
    if any(k in p for k in ['parb', 'para', 'par', 'partition', 'centromere']):
        return "partition"

    # DNA metabolism
    if any(k in p for k in ['dna ligase','helicase','methyltransferase','methylase',
                            'nuclease','recb','hera','whia','whib',
                            'topoisomerase', 'gyrase', 'polymerase', 'dna primase',
                            'exonuclease', 'endonuclease']):
        return "dna_metabolism"
    if any(k in p for k in ['toxin','antitoxin','vapc','vapb','hipa',
                            'hnh endonuclease','abiei','abigii','pin domain',
                            'xre family','res domain']):
        return "ta_system"
    if any(k in p for k in ['restriction', 'modification', 'crispr', 'cas']):
        return "defence"

    # Metabolism / biosynthesis
    if any(k in p for k in ['dehydrogenase', 'oxidoreductase', 'reductase', 'oxidase',
                            'oxygenase', 'monooxygenase', 'dioxygenase']):
        return "metabolism"
    if any(k in p for k in ['synthase', 'synthetase', 'transferase', 'lyase',
                            'hydrolase', 'phosphatase', 'kinase', 'phosphorylase',
                            'isomerase', 'epimerase', 'racemase']):
        return "metabolism"
    if any(k in p for k in ['ribosomal', 'trna', 'srrna', 'translation',
                            'elongation factor', 'initiation factor']):
        return "translation"
    if any(k in p for k in ['rnase', 'ribonuclease', 'nuclease',
                            'transcription factor', 'sigma factor', 'rna polymerase']):
        return "transcription_regulation"
    if any(k in p for k in ['chaperone', 'heat shock', 'protease', 'peptidase',
                            'proteinase', 'proteasome', 'clpp', 'clpx', 'htpx', 'dnak', 'dnal', 'grpe']):
        return "protein_folding_turnover"

    # Regulation / signalling
    if any(k in p for k in ['regulator', 'transcriptional regulator', 'response regulator',
                            'two-component', 'histidine kinase', 'sensor',
                            'repressor', 'activator', 'tetr', 'gntr', 'lysr', 'laci',
                            'arac', 'xre', 'whib', 'sarp', 'luxr']):
        return "regulation"
    if any(k in p for k in ['sigma factor', 'anti-sigma', 'ecf']):
        return "regulation"

    # Stress / resistance
    if any(k in p for k in ['antibiotic resistance', 'beta-lactamase', 'drug resistance',
                            'chloramphenicol', 'tetracycline', 'kanamycin',
                            'streptomycin', 'rifampin', 'vancomycin',
                            'resistance protein', 'macrolide']):
        return "resistance"
    if any(k in p for k in ['stress protein', 'cold shock', 'osmotic', 'dna repair',
                            'uvr', 'uvra', 'uvrb', 'uvrc', 'muts', 'mutl', 'rec',
                            'formamidopyrimidine', 'sod', 'peroxidase', 'catalase',
                            'superoxide', 'glutaredoxin', 'thioredoxin']):
        return "stress_defence"

    return "other"


# 23 visually-distinct colours (hand-picked for hue + lightness separation).
# Grey reserved for hypothetical / other. No two adjacent-by-letter entries
# share the same hue family.
CAT_COLORS = {
    "conjugation":              "#E41A1C",  # crimson red
    "resistance":               "#FF6F00",  # bright orange
    "secretion":                "#FFC000",  # amber / gold
    "efflux_transport":         "#B45F06",  # burnt orange (sienna)
    "ta_system":                "#FFD700",  # pure gold (kept — strong contrast)
    "t7ss_esx":                 "#984EA3",  # purple
    "regulation":               "#6A3D9A",  # deep violet
    "transcription_regulation": "#CAB2D6",  # pale lavender
    "replication":              "#377EB8",  # strong blue
    "partition":                "#00BFC4",  # teal / cyan
    "stress_defence":           "#1F78B4",  # mid-blue (far enough from replication)
    "dna_metabolism":           "#7B5B00",  # dark mustard / olive-brown
    "is_element":               "#4DAF4A",  # medium green
    "metabolism":               "#33A02C",  # dark green
    "translation":              "#B2DF8A",  # pale green
    "motility_adhesion":        "#F781BF",  # magenta-pink
    "transport":                "#E78AC3",  # pink (lighter, distinct from FB9A99)
    "protein_folding_turnover": "#FB9A99",  # salmon
    "phage":                    "#A6D854",  # lime / yellow-green
    "mobile_element":           "#A65628",  # chocolate brown
    "defence":                  "#8B4513",  # saddle brown
    "hypothetical":             "#D9D9D9",  # light grey
    "other":                    "#B3B3B3",  # mid grey
}

CAT_LABELS = {
    "conjugation":              "Conjugation (TrbL/VirB4/TcpC)",
    "is_element":               "IS Elements / Transposases",
    "mobile_element":           "Integrases / Recombinases",
    "phage":                    "Phage / Prophage",
    "t7ss_esx":                 "T7SS / ESX-1 secretion",
    "secretion":                "Secretion (T1SS–T6SS)",
    "motility_adhesion":        "Motility / Adhesion",
    "efflux_transport":         "Efflux / MmpL / Antiporter",
    "transport":                "Transport / Permeases",
    "replication":              "Replication / Maintenance",
    "partition":                "Partition proteins",
    "dna_metabolism":           "DNA metabolism",
    "ta_system":                "Toxin-Antitoxin systems",
    "defence":                  "Restriction / CRISPR / Defence",
    "metabolism":               "Primary / secondary metabolism",
    "translation":              "Ribosomal / Translation",
    "transcription_regulation": "Transcription / RNA metabolism",
    "protein_folding_turnover": "Protein folding / turnover",
    "regulation":               "Regulatory / Signalling",
    "resistance":               "Antibiotic resistance",
    "stress_defence":           "Stress / DNA repair / ROS",
    "hypothetical":             "Hypothetical",
    "other":                    "Other",
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
    print("Building linear genome comparison map (dual-axis)...")
    have_ref_gff = REF_GFF.exists()
    have_qry_gff = QRY_GFF.exists()
    ref_feats = parse_gff(REF_GFF) if have_ref_gff else []
    qry_feats = parse_gff(QRY_GFF) if have_qry_gff else []
    if have_ref_gff or have_qry_gff:
        print(f"  Reference: {len(ref_feats)} CDS  |  Query: {len(qry_feats)} CDS")
    else:
        print("  No GFF annotations — genome bars + ribbons only (no ORF tracks)")

    # ── Figure: main axes on left, legend on right ──
    # Top track = reference coordinates; bottom track = query coordinates (dual axis)
    fig = plt.figure(figsize=(28, 14), facecolor='white')
    gs = fig.add_gridspec(2, 1, left=0.05, right=0.78, top=0.94, bottom=0.12,
                          hspace=0.45, height_ratios=[1, 1])
    ax_top = fig.add_subplot(gs[0])
    ax_bot = fig.add_subplot(gs[1])   # independent x-axis (query coords)

    # ── Helper: draw one track (in its own coordinate space) ──
    def draw_track(ax, feats, islands, key_genes, title, track_len,
                   label_islands=False, label_side="top"):
        """Draw a linear genome track with ORFs, islands, and key genes.
        Coordinates are the track's own (ref for top, query for bottom),
        so every CDS is shown regardless of alignment coverage.

        When *label_islands* is True (islands from a user TSV), each island
        band is annotated with its text label, staggered across tiers.
        *label_side*="top" places labels ABOVE the track (positive y, arrow
        down); *label_side*="bottom" places them BELOW (negative y, arrow up)."""
        # Extend the track to make room for the label tiers on the side where
        # they live.
        top_lim    = 3.6 if (label_islands and label_side == "top")    else 1.6
        bottom_lim = -3.6 if (label_islands and label_side == "bottom") else -1.6
        ax.set_xlim(-track_len * 0.015, track_len * 1.015)
        ax.set_ylim(bottom_lim, top_lim)

        N_TIERS   = 4          # number of staggered label tiers
        TIER_STEP = 0.35       # vertical gap between tiers
        if label_side == "top":
            TIER_BASE = 0.55   # first tier just above +strand ORFs
            band_edge = -1.55  # arrow anchors at the band's top edge
            va = 'bottom'
        else:
            TIER_BASE = -1.75  # first tier just below the band
            band_edge = 1.55   # arrow anchors at the band's bottom edge
            va = 'top'

        for idx, (istart, iend, ilabel, icolor) in enumerate(islands):
            rect = Rectangle((istart, -1.55), iend - istart, 3.1,
                             facecolor=icolor, alpha=0.07, edgecolor='none', zorder=0)
            ax.add_patch(rect)
            if label_islands:
                mid = (istart + iend) / 2
                tier = idx % N_TIERS
                if label_side == "top":
                    label_y = TIER_BASE + tier * TIER_STEP
                else:
                    label_y = TIER_BASE - tier * TIER_STEP
                lbl = ilabel if len(ilabel) <= 26 else ilabel[:24] + '…'
                ax.annotate(lbl, xy=(mid, band_edge), xytext=(mid, label_y),
                            fontsize=7, color=icolor, fontweight='bold',
                            ha='center', va=va, rotation=0,
                            arrowprops=dict(arrowstyle='-', color=icolor,
                                            lw=0.4, alpha=0.4),
                            bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                                      edgecolor=icolor, alpha=0.6, lw=0.4),
                            zorder=11)

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

        # Draw small key–gene highlights — only when key_genes list is non‑empty;
        # this is intended for low‑throughput manual use (supply via env or file).
        for idx, (ks, ke, klabel, kcolor) in enumerate(key_genes):
            rect = Rectangle((ks, -orf_h - 0.12), ke - ks, orf_h * 2 + 0.24,
                             facecolor=kcolor, alpha=0.18, edgecolor=kcolor,
                             linewidth=1.0, zorder=5)
            ax.add_patch(rect)
            mid = (ks + ke) / 2
            label_y = 1.25 if idx % 2 == 0 else -1.25
            va = 'bottom' if label_y > 0 else 'top'
            ax.annotate(klabel, xy=(mid, 0), xytext=(mid, label_y),
                        fontsize=6, color=kcolor, fontweight='bold',
                        ha='center', va=va,
                        arrowprops=dict(arrowstyle='-', color=kcolor, lw=0.4, alpha=0.4),
                        zorder=11)

        # ── Title (inside plot area, top-left) ──
        # Sits at the very top of the track; since island labels go BELOW the
        # central axis there is never any collision with them.
        ax.text(track_len * 0.003, 1.55, title, fontsize=10, fontweight='bold',
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

    # ── Resolve genomic islands (independent ref / qry modes) ──────────
    # Each side accepts: 'auto' | 'none' | <path-to-TSV>
    # ref and qry are resolved independently so the user can supply a TSV
    # for one genome and 'auto' for the other (or nothing at all).
    default_mode = "auto" if (have_ref_gff or have_qry_gff) else "none"
    ref_islands_mode = os.environ.get("SYNTENY_ISLANDS_REF", default_mode)
    qry_islands_mode = os.environ.get("SYNTENY_ISLANDS_QRY", default_mode)

    def _resolve_islands(mode, feats, track_len, label):
        """Return (islands_list, label_flag) for one side."""
        if mode == "none":
            return [], False
        if mode == "auto":
            return (auto_islands(feats, track_len) if feats else []), False
        if os.path.isfile(mode):
            return parse_island_file(mode), True
        print(f"  ⚠ Unknown --islands-{label} value '{mode}' — falling back to 'auto'")
        return (auto_islands(feats, track_len) if feats else []), False

    ref_islands, label_ref = _resolve_islands(ref_islands_mode, ref_feats, REF_LEN, "ref")
    qry_islands, label_qry = _resolve_islands(qry_islands_mode, qry_feats, QRY_LEN, "qry")

    # ── Draw tracks (each in its own coordinate space) ──
    # In-plot island labels are only drawn for sides whose islands came
    # from a user TSV file (auto / none modes keep the bands unlabelled).
    # Reference = top track → labels ABOVE the ORF track (positive y, arrow down)
    # Query     = bottom track → labels BELOW the ORF track (negative y, arrow up)
    if have_ref_gff:
        draw_track(ax_top, ref_feats, ref_islands, [],
                   f"Reference  ({REF_LEN:,} bp)", REF_LEN,
                   label_islands=label_ref, label_side="top")
    else:
        draw_track(ax_top, [], ref_islands, [],
                   f"Reference  ({REF_LEN:,} bp)  [FASTA only — no annotations]", REF_LEN,
                   label_islands=label_ref, label_side="top")
        ax_top.text(REF_LEN / 2, 0, "No GFF annotations provided", ha='center', va='center',
                    fontsize=12, color='#999', style='italic', zorder=20)

    if have_qry_gff:
        draw_track(ax_bot, qry_feats, qry_islands, [],
                   f"Query  ({QRY_LEN:,} bp)", QRY_LEN,
                   label_islands=label_qry, label_side="bottom")
    else:
        draw_track(ax_bot, [], qry_islands, [],
                   f"Query  ({QRY_LEN:,} bp)  [FASTA only — no annotations]", QRY_LEN,
                   label_islands=label_qry, label_side="bottom")
        ax_bot.text(QRY_LEN / 2, 0, "No GFF annotations provided", ha='center', va='center',
                    fontsize=12, color='#999', style='italic', zorder=20)

    # ── Aligned-block shading + connector ribbon ──
    RIBBON_TSV = _pp.out_path("high_identity_regions", ".tsv")
    overlay = fig.add_axes([0, 0, 1, 1], frameon=False)
    overlay.set_xlim(0, 1); overlay.set_ylim(0, 1); overlay.axis('off')

    def _to_fig(ax, x, y):
        disp = ax.transData.transform((x, y))
        return tuple(fig.transFigure.inverted().transform(disp))

    if HAS_ALN:
        # Shaded 1:1 aligned span on each track
        ax_top.axvspan(REF_ALN_START, REF_ALN_END, alpha=0.05, color='#377EB8', zorder=0)
        ax_bot.axvspan(QRY_ALN_START, QRY_ALN_END, alpha=0.05, color='#377EB8', zorder=0)

        # Connector polygon between the two axes
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
                # columns are rep-oriented reference and query coordinates
                ra, rb = int(_r['ref_start_rep']), int(_r['ref_end_rep'])
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

    # ── Right-hand panel: legend (top) + summary box (bottom) ──
    # Two stacked axes in a single right-hand column. The legend's height
    # is sized to its entry count so the summary box sits just below it.
    LEG_LEFT, LEG_W = 0.80, 0.18
    LEG_TOP, LEG_BOT = 0.94, 0.12
    TITLE_GAP = 0.012   # space reserved for the "Feature Categories" title
    PANEL_GAP = 0.015   # gap between legend and summary box

    # Show only categories that appear in the data (plus always hypothetical).
    cats_in_data = set(f['category'] for f in ref_feats + qry_feats)
    legend_cats = [c for c in CAT_LABELS
                   if c in cats_in_data or c == 'hypothetical']
    handles = []
    for cat in legend_cats:
        handles.append(Line2D([0],[0], marker='s', color='w',
                              markerfacecolor=CAT_COLORS[cat], markersize=10,
                              label=CAT_LABELS[cat]))

    # Estimate legend height: ~0.026 figure units per entry (fontsize 8,
    # markersize 10) plus a small fixed border/frame overhead.
    legend_h = min(LEG_TOP - LEG_BOT - 0.05,
                   0.026 * len(handles) + 0.035)

    legend_ax = fig.add_axes([LEG_LEFT, LEG_TOP - legend_h, LEG_W, legend_h])
    legend_ax.axis('off')
    legend = legend_ax.legend(handles=handles, loc='upper left', frameon=True,
                              fontsize=8, handletextpad=1.0, borderpad=0.6)
    # Title: placed just above the legend's top edge, anchored to the figure
    # (uses transAxes so it always sits a small gap above the frame).
    legend_ax.text(0, 1.0 + TITLE_GAP / max(legend_h, 1e-6),
                   "Feature Categories", fontsize=11, fontweight='bold',
                   transform=legend_ax.transAxes, va='bottom', ha='left')

    # ── Summary box (separate axes, packed just below the legend) ──
    summary_top = LEG_TOP - legend_h - PANEL_GAP
    summary_h = summary_top - LEG_BOT
    summary_ax = fig.add_axes([LEG_LEFT, LEG_BOT, LEG_W, summary_h])
    summary_ax.axis('off')

    summary_lines = [
        "Genome Comparison",
        "─────────────────",
        f"Reference: {REF_LEN:,} bp",
        f"Query:     {QRY_LEN:,} bp",
    ]
    if have_ref_gff and have_qry_gff:
        summary_lines += ["", "ORF tracks + islands", "from GFF annotations"]
    else:
        summary_lines += [
            "",
            "⚠ FASTA-only mode" if not (have_ref_gff or have_qry_gff) else "⚠ partial annotations",
            "ORF tracks not drawn",
            "where GFF is missing",
        ]
    if HAS_ALN:
        summary_lines += [
            "",
            "Blue = 1:1 aligned block",
            "(connector links axes)",
        ]
    else:
        summary_lines += ["", "No 1:1 aligned-block", "coords supplied"]
    summary_lines += [
        "",
        "Light-blue ribbons =",
        f"conserved >{ID_THR}% blocks",
    ]
    summary = "\n".join(summary_lines)

    summary_ax.text(0, 1.0, summary, fontsize=8.5, fontfamily='monospace',
                    transform=summary_ax.transAxes, va='top', ha='left',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow',
                              edgecolor='#ccc'))

    fig.savefig(str(OUT_SVG), dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  ✓ Saved: {OUT_SVG}")

if __name__ == "__main__":
    main()
