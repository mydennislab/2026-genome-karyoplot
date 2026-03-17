#!/usr/bin/env python3
"""
Generate realistic synthetic test data for genome-karyotype.

Creates a small 6-chromosome genome (~180 Mb total) with:
  - Centromere calls
  - tidk telomere output
  - Assembly gaps
  - RepeatMasker .out (LINE/SINE/LTR/DNA/Simple)
  - Segmental duplications BEDPE

Repeat densities vary by chromosome region to produce
visually interesting karyotype figures.
"""

import random
import os

random.seed(42)

OUTDIR = os.path.dirname(os.path.abspath(__file__))

# ── Genome definition ─────────────────────────────────────────────────────
CHROMS = [
    ("chr1", 50_000_000),
    ("chr2", 42_000_000),
    ("chr3", 35_000_000),
    ("chr4", 28_000_000),
    ("chr5", 20_000_000),
    ("chrX", 32_000_000),
]

# Centromere positions (start, end) — roughly 1-4 Mb regions
CENTROMERES = {
    "chr1": (22_000_000, 24_500_000),
    "chr2": (18_000_000, 19_200_000),
    "chr3": (5_000_000,  6_800_000),
    "chr4": (12_000_000, 13_500_000),
    "chr5": (9_000_000,  10_000_000),
    "chrX": (14_000_000, 16_000_000),
}


def write_chrom_sizes():
    path = os.path.join(OUTDIR, "chrom.sizes")
    with open(path, "w") as f:
        for chrom, size in CHROMS:
            f.write(f"{chrom}\t{size}\n")
    print(f"  {path}")


def write_centromere_bed():
    path = os.path.join(OUTDIR, "centromeres.bed")
    with open(path, "w") as f:
        for chrom, _size in CHROMS:
            cs, ce = CENTROMERES[chrom]
            f.write(f"{chrom}\t{cs}\t{ce}\n")
    print(f"  {path}")


def write_tidk_tsv():
    """Generate tidk telomere output with 10 kb windows."""
    path = os.path.join(OUTDIR, "tidk.tsv")
    window = 10_000
    with open(path, "w") as f:
        f.write("id\twindow\tforward_repeat_number\treverse_repeat_number\n")
        for chrom, size in CHROMS:
            n_windows = size // window
            for i in range(n_windows):
                pos = i * window
                # Telomere signal: high at ends, zero in middle
                if i < 5:
                    fwd = random.randint(15, 80)
                    rev = random.randint(0, 3)
                elif i >= n_windows - 5:
                    fwd = random.randint(0, 3)
                    rev = random.randint(15, 80)
                else:
                    fwd = 0
                    rev = 0
                f.write(f"{chrom}\t{pos}\t{fwd}\t{rev}\n")
    print(f"  {path}")


def write_gaps_bed():
    """Generate a few assembly gaps per chromosome."""
    path = os.path.join(OUTDIR, "gaps.bed")
    with open(path, "w") as f:
        for chrom, size in CHROMS:
            # 0-3 gaps per chromosome, not near telomeres
            n_gaps = random.randint(0, 3)
            for _ in range(n_gaps):
                start = random.randint(size // 10, size - size // 10)
                gap_len = random.randint(100, 5000)
                f.write(f"{chrom}\t{start}\t{start + gap_len}\n")
    print(f"  {path}")


def region_repeat_rate(chrom, pos, size, repeat_class):
    """
    Return a probability of placing a repeat element at this position.
    Creates realistic density patterns:
      - Pericentromeric regions: high satellite/simple, high LINE
      - Subtelomeric: moderate SINE
      - Gene-rich regions (middle): lower repeats
    """
    cs, ce = CENTROMERES[chrom]
    frac = pos / size  # 0..1 along chromosome

    # Distance from centromere (0 = at centromere, 1 = far)
    centro_mid = (cs + ce) / 2 / size
    dist_centro = abs(frac - centro_mid)

    # Base rates per class (probability per 2 kb step)
    # Tuned so stacked density peaks ~50-60% with visible layers
    base = {
        "LINE":             0.025,
        "SINE":             0.022,
        "LTR":              0.018,
        "DNA":              0.015,
        "Simple_repeat":    0.008,
        "Satellite":        0.004,
    }

    rate = base.get(repeat_class, 0.01)

    # Boost near centromere for satellite/simple
    if repeat_class in ("Satellite", "Simple_repeat"):
        if dist_centro < 0.05:
            rate *= 15
        elif dist_centro < 0.10:
            rate *= 5

    # LINE enriched in pericentromeric and subtelomeric
    if repeat_class == "LINE":
        if dist_centro < 0.08:
            rate *= 3
        if frac < 0.05 or frac > 0.95:
            rate *= 2

    # SINE enriched in gene-rich (middle, away from centromere)
    if repeat_class == "SINE":
        if 0.2 < frac < 0.8 and dist_centro > 0.15:
            rate *= 2.5

    # LTR slightly enriched near centromere
    if repeat_class == "LTR":
        if dist_centro < 0.10:
            rate *= 2

    return rate


def write_repeatmasker_out():
    """Generate RepeatMasker .out with position-dependent density."""
    path = os.path.join(OUTDIR, "repeats.out")

    rm_classes = [
        ("LINE", "LINE/L1"),
        ("LINE", "LINE/L2"),
        ("SINE", "SINE/Alu"),
        ("SINE", "SINE/MIR"),
        ("LTR",  "LTR/ERV1"),
        ("LTR",  "LTR/ERVK"),
        ("DNA",  "DNA/hAT-Charlie"),
        ("DNA",  "DNA/TcMar-Tigger"),
        ("Simple_repeat", "Simple_repeat"),
        ("Satellite", "Satellite/centr"),
    ]

    step = 2000  # check every 2 kb
    elem_id = 1

    with open(path, "w") as f:
        # RepeatMasker header (3 lines)
        f.write("   SW   perc perc perc  query      position in query    matching          repeat          position in repeat\n")
        f.write("score   div. del. ins.  sequence   begin    end   (left) strand  repeat   class/family    begin   end  (left)  ID\n")
        f.write("\n")

        for chrom, size in CHROMS:
            for pos in range(0, size, step):
                for top_class, class_family in rm_classes:
                    rate = region_repeat_rate(chrom, pos, size, top_class)
                    if random.random() < rate:
                        # Random element length
                        if top_class == "LINE":
                            elem_len = random.randint(200, 6000)
                        elif top_class == "SINE":
                            elem_len = random.randint(100, 350)
                        elif top_class == "LTR":
                            elem_len = random.randint(300, 3000)
                        elif top_class == "DNA":
                            elem_len = random.randint(100, 1500)
                        else:  # Simple/Satellite
                            elem_len = random.randint(50, 2000)

                        start = pos + random.randint(0, step)
                        end = min(start + elem_len, size)
                        score = random.randint(200, 2000)
                        div = round(random.uniform(1, 35), 1)
                        name = class_family.split("/")[-1] if "/" in class_family else class_family

                        left = size - end
                        f.write(
                            f"  {score:>4}  {div:>5}  0.0  0.0  "
                            f"{chrom:<14} {start:>9} {end:>9} ({left}) + "
                            f"{name:<16} {class_family:<25} 1 {elem_len:>5} (0)  {elem_id}\n"
                        )
                        elem_id += 1

    print(f"  {path}  ({elem_id - 1:,} elements)")


def write_segdups_bedpe():
    """Generate segmental duplications: intra and inter-chromosomal."""
    path = os.path.join(OUTDIR, "segdups.bedpe")
    with open(path, "w") as f:
        f.write("#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tstrand\terror_rate\n")

        for chrom, size in CHROMS:
            # 2-6 intra-chromosomal SDs per chromosome
            n_sd = random.randint(2, 6)
            for _ in range(n_sd):
                len1 = random.randint(5_000, 80_000)
                s1 = random.randint(0, size - len1)
                e1 = s1 + len1
                # Second copy displaced by 0.5-5 Mb
                offset = random.randint(500_000, 5_000_000)
                s2 = min(s1 + offset, size - len1)
                e2 = s2 + len1
                error = round(random.uniform(1, 8), 2)
                f.write(f"{chrom}\t{s1}\t{e1}\t{chrom}\t{s2}\t{e2}\t+\t{error}\n")

        # A few inter-chromosomal SDs
        for _ in range(8):
            c1, s1_size = random.choice(CHROMS)
            c2, s2_size = random.choice(CHROMS)
            if c1 == c2:
                continue
            sd_len = random.randint(10_000, 50_000)
            s1 = random.randint(0, s1_size - sd_len)
            s2 = random.randint(0, s2_size - sd_len)
            error = round(random.uniform(1, 8), 2)
            f.write(f"{c1}\t{s1}\t{s1 + sd_len}\t{c2}\t{s2}\t{s2 + sd_len}\t+\t{error}\n")

    print(f"  {path}")


def main():
    print("Generating test data for genome-karyotype...")
    write_chrom_sizes()
    write_centromere_bed()
    write_tidk_tsv()
    write_gaps_bed()
    write_repeatmasker_out()
    write_segdups_bedpe()
    print("Done.")


if __name__ == "__main__":
    main()
