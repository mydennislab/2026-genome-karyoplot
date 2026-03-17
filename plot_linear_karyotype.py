#!/usr/bin/env python3
"""
Linear karyotype ideogram with stacked repeat density tracks.
Requires: numpy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import argparse, re
from pathlib import Path

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8, 'figure.dpi': 300, 'savefig.dpi': 300,
})

# Colors
COL_LINE   = '#E64B35'
COL_SINE   = '#4DBBD5'
COL_LTR    = '#F39B7F'
COL_DNA_TE = '#8491B4'
COL_SIMPLE = '#91D1C2'
COL_CHR_FILL = '#E8E8E8'
COL_CHR_EDGE = '#424242'
COL_CENTRO   = '#424242'
COL_GAP      = '#D32F2F'
COL_TELO_YES = '#2E7D32'
COL_TELO_NO  = '#EF5350'
COL_SD       = '#7B1FA2'

REPEAT_CLASSES = ['Simple/Satellite', 'DNA', 'LTR', 'SINE', 'LINE']
CLASS_COLORS = {
    'LINE': COL_LINE, 'SINE': COL_SINE, 'LTR': COL_LTR,
    'DNA': COL_DNA_TE, 'Simple/Satellite': COL_SIMPLE,
}
CLASS_MAP = {
    'LINE': 'LINE', 'SINE': 'SINE', 'LTR': 'LTR', 'DNA': 'DNA',
    'Simple_repeat': 'Simple/Satellite',
    'Satellite': 'Simple/Satellite',
    'Low_complexity': 'Simple/Satellite',
}


# --- Parsers ---

def parse_chrom_sizes(filepath):
    """Parse .fai or 2-column TSV (chrom<TAB>size)."""
    sizes = {}
    with open(filepath) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                try:
                    sizes[parts[0]] = int(parts[1])
                except ValueError:
                    continue
    return sizes


def parse_bed(filepath):
    if filepath is None or not Path(filepath).exists():
        return {}
    regions = {}
    with open(filepath) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                try:
                    regions.setdefault(parts[0], []).append(
                        (int(parts[1]), int(parts[2])))
                except ValueError:
                    continue
    return regions


def parse_tidk(filepath):
    """Parse tidk search output — just needs first and last window per chrom."""
    if filepath is None or not Path(filepath).exists():
        return {}
    telomeres = {}
    header = None
    # Group rows by chrom, keep first and last window
    chrom_rows = {}
    with open(filepath) as f:
        for line in f:
            if header is None:
                header = line.strip().split('\t')
                continue
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            chrom = parts[0]
            fwd, rev = int(parts[2]), int(parts[3])
            total = fwd + rev
            window = int(parts[1])
            if chrom not in chrom_rows:
                chrom_rows[chrom] = {'first_win': window, 'first_val': total,
                                     'last_win': window, 'last_val': total}
            else:
                if window < chrom_rows[chrom]['first_win']:
                    chrom_rows[chrom]['first_win'] = window
                    chrom_rows[chrom]['first_val'] = total
                if window > chrom_rows[chrom]['last_win']:
                    chrom_rows[chrom]['last_win'] = window
                    chrom_rows[chrom]['last_val'] = total
    for chrom, d in chrom_rows.items():
        telomeres[chrom] = {'5prime': d['first_val'], '3prime': d['last_val']}
    return telomeres


def parse_segdups(filepath):
    """Parse BED or BEDPE segmental duplications.
    BEDPE (BISER/SEDEF format) with divergence in col 8: keeps <10% divergence and >=1 kb."""
    if filepath is None or not Path(filepath).exists():
        return {}
    regions = {}
    is_bedpe = None
    with open(filepath) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            # detect format from first data line
            if is_bedpe is None:
                is_bedpe = (len(parts) >= 6
                            and not parts[3].isdigit()
                            and bool(re.match(r'^[A-Za-z_]', parts[3])))
            if is_bedpe:
                try:
                    chrom1, s1, e1 = parts[0], int(parts[1]), int(parts[2])
                    chrom2, s2, e2 = parts[3], int(parts[4]), int(parts[5])
                except ValueError:
                    continue
                if len(parts) > 7:
                    try:
                        divergence = float(parts[7])
                        if divergence >= 10.0 or ((e1-s1)+(e2-s2))/2 < 1000:
                            continue
                    except ValueError:
                        pass
                for c, s, e in [(chrom1, s1, e1), (chrom2, s2, e2)]:
                    regions.setdefault(c, []).append((s, e))
            else:
                try:
                    regions.setdefault(parts[0], []).append(
                        (int(parts[1]), int(parts[2])))
                except ValueError:
                    continue
    return regions


def parse_repeatmasker(filepath):
    """Parse RepeatMasker .out, grouping into LINE/SINE/LTR/DNA/Simple+Satellite."""
    by_class = {c: {} for c in REPEAT_CLASSES}
    if filepath is None or not Path(filepath).exists():
        return by_class
    with open(filepath) as f:
        for i, line in enumerate(f):
            if i < 3:
                continue
            parts = line.split()
            if len(parts) < 11:
                continue
            try:
                chrom = parts[4]
                start, end = int(parts[5]), int(parts[6])
                class_family = parts[10]
            except (ValueError, IndexError):
                continue
            group = CLASS_MAP.get(class_family.split('/')[0])
            if group:
                by_class[group].setdefault(chrom, []).append((start, end))
    return by_class


# --- Helpers ---

def strip_prefix(name, prefix):
    if prefix and name.startswith(prefix):
        return name[len(prefix):]
    return name


def chrom_sort_key(chrom, prefix=""):
    name = strip_prefix(chrom, prefix)
    if name.startswith('chr'):
        name = name[3:]
    if name == 'X':  return (1, 100, chrom)
    if name == 'Y':  return (1, 101, chrom)
    if name in ('M', 'MT'): return (1, 102, chrom)
    try:
        return (0, int(name), chrom)
    except ValueError:
        m = re.match(r'(\d+)', name)
        return (0, int(m.group(1)), chrom) if m else (2, 0, chrom)


def compute_density(features, chrom_size, window_size):
    n_win = int(np.ceil(chrom_size / window_size))
    positions = np.minimum(
        np.arange(n_win) * window_size + window_size / 2, chrom_size - 1)
    starts = np.array([s for s, e in features])
    ends = np.array([e for s, e in features])
    lengths = ends - starts
    bins = np.clip(starts // window_size, 0, n_win - 1).astype(int)
    density = np.zeros(n_win)
    np.add.at(density, bins, lengths)
    density /= window_size
    np.minimum(density, 1.0, out=density)
    return positions, density


# --- Main ---

def plot_karyotype(chrom_sizes, centromeres=None, telomeres=None,
                   gaps=None, repeats_by_class=None, segdups=None,
                   output_prefix='assembly', title='Karyotype Ideogram',
                   window_size=500_000, chrom_prefix="", min_size=0):

    centromeres = centromeres or {}
    telomeres = telomeres or {}
    gaps = gaps or {}
    segdups = segdups or {}
    if repeats_by_class is None:
        repeats_by_class = {c: {} for c in REPEAT_CLASSES}

    # Filter and sort chromosomes
    chromosomes = sorted(
        [c for c in chrom_sizes
         if (not chrom_prefix or c.startswith(chrom_prefix))
         and chrom_sizes[c] >= min_size],
        key=lambda c: chrom_sort_key(c, chrom_prefix))

    if not chromosomes:
        print("ERROR: No chromosomes matched. Check --prefix and --min-size.")
        return
    max_size = max(chrom_sizes[c] for c in chromosomes)
    n_chrom = len(chromosomes)
    print(f"  {n_chrom} chromosomes, largest = {max_size / 1e6:.1f} Mb")

    # Compute repeat densities per window
    has_repeats = any(
        v for cls_dict in repeats_by_class.values() for v in cls_dict.values())
    all_densities = {}
    global_max_stacked = 0.0

    for chrom in chromosomes:
        cs = chrom_sizes[chrom]
        chrom_dens = {}
        for cls in REPEAT_CLASSES:
            features = repeats_by_class[cls].get(chrom, [])
            if features:
                pos, dens = compute_density(features, cs, window_size)
            else:
                n_win = int(np.ceil(cs / window_size))
                pos = np.arange(n_win) * window_size + window_size / 2
                pos = np.minimum(pos, cs - 1)
                dens = np.zeros(n_win)
            chrom_dens[cls] = (pos, dens)

        stacked = sum(chrom_dens[cls][1] for cls in REPEAT_CLASSES)
        global_max_stacked = max(global_max_stacked, stacked.max() if len(stacked) else 0)
        all_densities[chrom] = chrom_dens

    if global_max_stacked == 0:
        global_max_stacked = 1.0
    if has_repeats:
        print(f"  Global max stacked density: {global_max_stacked:.3f}")

    # Figure layout
    row_h = 0.6
    top_margin, bot_margin, legend_h = 0.8, 0.6, 0.8
    fig_h = n_chrom * row_h + top_margin + bot_margin + legend_h

    fig, ax = plt.subplots(figsize=(12, fig_h))
    x_max = max_size / 1e6 * 1.02
    y_total = n_chrom * row_h
    ax.set_xlim(-x_max * 0.08, x_max)
    ax.set_ylim(-bot_margin - legend_h, y_total + 0.2)
    ax.set_xlabel('Position (Mb)', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(left=False, labelleft=False)

    # Draw chromosomes
    print("  Drawing chromosomes...")
    chr_bar_h = row_h * 0.14
    area_h = row_h * 0.55
    sd_h = row_h * 0.12
    gap_between = row_h * 0.04

    for i, chrom in enumerate(chromosomes):
        cs_mb = chrom_sizes[chrom] / 1e6
        row_top = y_total - i * row_h
        chr_bar_y = row_top - area_h - gap_between - chr_bar_h
        area_base = chr_bar_y + chr_bar_h + gap_between
        area_top = row_top - row_h * 0.05
        sd_y = chr_bar_y - gap_between - sd_h

        ax.text(-x_max * 0.02, chr_bar_y + chr_bar_h / 2,
                strip_prefix(chrom, chrom_prefix),
                ha='right', va='center', fontsize=7, fontweight='bold')

        ax.add_patch(FancyBboxPatch(
            (0, chr_bar_y), cs_mb, chr_bar_h,
            boxstyle=f"round,pad=0,rounding_size={chr_bar_h * 0.4}",
            facecolor=COL_CHR_FILL, edgecolor=COL_CHR_EDGE,
            linewidth=0.6, zorder=2))

        for c_start, c_end in centromeres.get(chrom, []):
            ax.fill_between([c_start / 1e6, c_end / 1e6],
                            chr_bar_y, chr_bar_y + chr_bar_h,
                            color=COL_CENTRO, zorder=3)

        for gs, ge in gaps.get(chrom, []):
            g_mid = (gs + ge) / 2 / 1e6
            ax.plot([g_mid, g_mid], [chr_bar_y, chr_bar_y + chr_bar_h],
                    color=COL_GAP, linewidth=0.8, zorder=4)

        if telomeres:
            td = telomeres.get(chrom, {})
            for pos, key, marker in [(0, '5prime', '>'), (cs_mb, '3prime', '<')]:
                col = COL_TELO_YES if td.get(key, 0) >= 10 else COL_TELO_NO
                ax.scatter(pos, chr_bar_y + chr_bar_h / 2, marker=marker,
                           s=16, color=col, zorder=5, clip_on=False)

        if has_repeats:
            dens_data = all_densities[chrom]
            positions_mb = dens_data[REPEAT_CLASSES[0]][0] / 1e6
            cumulative = np.zeros_like(positions_mb)
            for cls in REPEAT_CLASSES:
                prev_y = area_base + cumulative / global_max_stacked * (area_top - area_base)
                cumulative = cumulative + dens_data[cls][1]
                curr_y = area_base + cumulative / global_max_stacked * (area_top - area_base)
                ax.fill_between(positions_mb, prev_y, curr_y,
                                facecolor=CLASS_COLORS[cls], alpha=0.85,
                                linewidth=0, zorder=1)

        for sd_s, sd_e in segdups.get(chrom, []):
            ax.add_patch(plt.Rectangle(
                (sd_s / 1e6, sd_y), (sd_e - sd_s) / 1e6, sd_h,
                facecolor=COL_SD, alpha=0.7, edgecolor='none', zorder=2))

    # Legend
    legend_elements = []
    if has_repeats:
        for cls in reversed(REPEAT_CLASSES):
            legend_elements.append(
                mpatches.Patch(facecolor=CLASS_COLORS[cls], alpha=0.85, label=cls))
    if segdups:
        legend_elements.append(mpatches.Patch(facecolor=COL_SD, alpha=0.7, label='Segmental dup.'))
    if centromeres:
        legend_elements.append(mpatches.Patch(facecolor=COL_CENTRO, label='Centromere'))
    if gaps:
        legend_elements.append(mpatches.Patch(facecolor=COL_GAP, label='Gap'))
    if telomeres:
        legend_elements.append(mpatches.Patch(facecolor=COL_TELO_YES, label='Telomere present'))
        legend_elements.append(mpatches.Patch(facecolor=COL_TELO_NO, label='Telomere absent'))
    if legend_elements:
        ax.legend(handles=legend_elements, loc='lower center',
                  bbox_to_anchor=(0.5, -0.02), ncol=5, fontsize=7,
                  frameon=True, edgecolor='#CCCCCC', handlelength=1.5,
                  handleheight=1.0, columnspacing=1.0,
                  bbox_transform=fig.transFigure)

    fig.suptitle(title, fontsize=12, fontweight='bold', y=0.98)
    for ext in ['png', 'pdf']:
        outpath = f"{output_prefix}_karyotype.{ext}"
        fig.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"Saved: {outpath}")
    plt.close()


if __name__ == '__main__':
    p = argparse.ArgumentParser(
        description='Linear karyotype ideogram with stacked repeat density')
    p.add_argument('--chrom-sizes', required=True,
                   help='.fai or 2-column TSV (chrom, size)')
    p.add_argument('--centromere', default=None, help='Centromere BED')
    p.add_argument('--tidk', default=None, help='tidk search TSV')
    p.add_argument('--gaps', default=None, help='Assembly gaps BED')
    p.add_argument('--repeatmasker', default=None, help='RepeatMasker .out')
    p.add_argument('--segdups', default=None, help='Segmental duplications BED/BEDPE')
    p.add_argument('-o', '--output', default='assembly', help='Output prefix')
    p.add_argument('-t', '--title', default='Karyotype Ideogram')
    p.add_argument('-w', '--window', type=int, default=500_000,
                   help='Density window size in bp (default: 500000)')
    p.add_argument('--prefix', default='', help='Chromosome name prefix filter')
    p.add_argument('--min-size', type=int, default=0, help='Min chromosome size (bp)')
    args = p.parse_args()

    print(f"Loading data for {args.title}...")
    chrom_sizes = parse_chrom_sizes(args.chrom_sizes)
    centromeres = parse_bed(args.centromere)
    gaps = parse_bed(args.gaps)
    telomeres = parse_tidk(args.tidk)
    segdups = parse_segdups(args.segdups)

    repeats_by_class = None
    if args.repeatmasker:
        print("  Parsing RepeatMasker (this may take a moment)...")
        repeats_by_class = parse_repeatmasker(args.repeatmasker)
        for cls in REPEAT_CLASSES:
            n = sum(len(v) for v in repeats_by_class[cls].values())
            print(f"    {cls}: {n:,} entries")

    plot_karyotype(
        chrom_sizes, centromeres=centromeres, telomeres=telomeres,
        gaps=gaps, repeats_by_class=repeats_by_class, segdups=segdups,
        output_prefix=args.output, title=args.title, window_size=args.window,
        chrom_prefix=args.prefix, min_size=args.min_size)
