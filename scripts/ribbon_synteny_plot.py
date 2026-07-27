#!/usr/bin/env python3
"""Curved-ribbon synteny plot between reference and query genomes.

Draws two horizontal genome bars (reference on top, query on bottom) and
connects each conserved (>90% identity) block with a smooth cubic-Bezier
ribbon. Blocks are kept separate so the synteny structure is visible.
"""

import csv
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as mpatches
import numpy as np

# ----------------------------------------------------------------------------
# Inputs / outputs (overridable via SYNTENY_* env vars)
# ----------------------------------------------------------------------------
import _pipeline_paths as _pp

REGIONS_TSV = str(_pp.out_path("high_identity_regions", ".tsv"))
OUT_SVG     = str(_pp.out_path("ribbon_synteny", ".svg"))

REF_LEN = _pp.REF_LEN   # reference, original orientation
QRY_LEN = _pp.QRY_LEN   # query, rep-oriented

# Vertical layout of the two genome bars
Y_TOP_C = 1.0      # centre of reference bar
Y_BOT_C = 0.0      # centre of query bar
BAR_H = 0.18       # half-height of each bar
Y_TOP_EDGE = Y_TOP_C - BAR_H   # bottom edge of top bar (ribbon attaches here)
Y_BOT_EDGE = Y_BOT_C + BAR_H   # top edge of bottom bar (ribbon attaches here)

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def ribbon_path(x0_top, x1_top, x0_bot, x1_bot, y_top, y_bot):
    """Closed polygon with cubic-Bezier (S-curve) left/right edges."""
    ym = 0.5 * (y_top + y_bot)
    verts = [
        (x0_top, y_top),
        (x1_top, y_top),
        (x1_top, ym), (x1_bot, ym), (x1_bot, y_bot),   # right edge curve
        (x0_bot, y_bot),
        (x0_bot, ym), (x0_top, ym), (x0_top, y_top),   # left edge curve
    ]
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.LINETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
    ]
    return Path(verts, codes)


def fmt_kb(x):
    return f"{int(x)//1000}kb"


# ----------------------------------------------------------------------------
def main():
    regions = []
    with open(REGIONS_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            regions.append(row)
    print(f"Loaded {len(regions)} conserved blocks")

    fig, ax = plt.subplots(figsize=(20, 8))
    xmax = max(REF_LEN, QRY_LEN) * 1.02

    # --- genome background bars -------------------------------------------------
    ax.add_patch(mpatches.Rectangle((0, Y_TOP_C - BAR_H), REF_LEN, 2 * BAR_H,
                                    facecolor="#e8e8e8", edgecolor="#999999",
                                    linewidth=0.8, zorder=1))
    ax.add_patch(mpatches.Rectangle((0, Y_BOT_C - BAR_H), QRY_LEN, 2 * BAR_H,
                                    facecolor="#e8e8e8", edgecolor="#999999",
                                    linewidth=0.8, zorder=1))

    # --- all ribbons / blocks light blue (no %identity colouring) -------------
    LIGHT_BLUE = "#9ecbff"
    RIBBON_FILL = (0.62, 0.80, 0.96, 0.45)

    # --- draw each block + ribbon ---------------------------------------------
    for r in regions:
        qa, qb = int(r["query_start"]), int(r["query_end"])
        ra, rb = int(r["ref_start_orig"]), int(r["ref_end_orig"])
        col = LIGHT_BLUE

        # blocks on each bar
        ax.add_patch(mpatches.Rectangle((ra, Y_TOP_C - BAR_H), rb - ra, 2 * BAR_H,
                                        facecolor=col, edgecolor="none", zorder=2))
        ax.add_patch(mpatches.Rectangle((qa, Y_BOT_C - BAR_H), qb - qa, 2 * BAR_H,
                                        facecolor=col, edgecolor="none", zorder=2))

        # curved ribbon connecting the two blocks
        path = ribbon_path(ra, rb, qa, qb, Y_TOP_EDGE, Y_BOT_EDGE)
        ax.add_patch(mpatches.PathPatch(path, facecolor=RIBBON_FILL,
                                        edgecolor="none", zorder=0))

    # --- axes / ticks ----------------------------------------------------------
    ax.set_xlim(0, xmax)
    ax.set_ylim(Y_BOT_C - BAR_H - 0.35, Y_TOP_C + BAR_H + 0.35)
    ax.set_yticks([])

    # bottom x-axis = query
    q_ticks = list(range(0, int(QRY_LEN) + 1, 20000))
    ax.set_xticks(q_ticks)
    ax.set_xticklabels([fmt_kb(t) for t in q_ticks], fontsize=9)
    ax.set_xlabel("Query (rep-oriented) position", fontsize=12)

    # top x-axis = reference (twiny)
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    r_ticks = list(range(0, int(REF_LEN) + 1, 20000))
    ax2.set_xticks(r_ticks)
    ax2.set_xticklabels([fmt_kb(t) for t in r_ticks], fontsize=9)
    ax2.set_xlabel("Reference (original orientation) position", fontsize=12)
    ax2.xaxis.set_ticks_position("top")
    ax2.xaxis.set_label_position("top")

    # bar labels
    ax.text(REF_LEN / 2, Y_TOP_C + BAR_H + 0.12, "Reference (orig. orient.)",
            ha="center", va="bottom", fontsize=12, fontweight="bold")
    ax.text(QRY_LEN / 2, Y_BOT_C - BAR_H - 0.12, "Query (rep-oriented)",
            ha="center", va="top", fontsize=12, fontweight="bold")

    # colour bar removed — ribbons are uniformly light blue
    ax.set_title("Curved-ribbon synteny: reference  ↔  query\n"
                 f"{len(regions)} conserved (>90%) blocks, total "
                 f"{sum(int(r['query_end'])-int(r['query_start']) for r in regions):,} bp",
                 fontsize=14, fontweight="bold")

    plt.tight_layout()
    fig.savefig(OUT_SVG, dpi=150, bbox_inches="tight")
    print(f"Wrote {OUT_SVG}")


if __name__ == "__main__":
    main()
