import re
import math
import matplotlib.pyplot as plt

# Not sure if this is normal: H5N1 segments 5 and 6 start with the NS sequence and H9N2 segment 3 has its first ~160 bases duplicated; sequences are used as they are in the files from NCBI.

GENOMES = {
    "H9N2 (2014)": "H9N2_2014.fna",
    "H5N1 (2022)": "H5N1_2022.fna",
}

ENZYMES = [
    {"name": "EcoRI", "site": "GAATTC", "cut": 1},
    {"name": "BamHI", "site": "GGATCC", "cut": 1},
    {"name": "HindIII", "site": "AAGCTT", "cut": 1},
]

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]", "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


def read_fasta(path):
    records = []
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header, seq = line[1:], []
            elif line:
                seq.append(line.upper())
    if header is not None:
        records.append((header, "".join(seq)))
    return records


SEGMENTS = {1: "PB2", 2: "PB1", 3: "PA", 4: "HA", 5: "NP", 6: "NA", 7: "M", 8: "NS"}


def segment_name(header):
    m = re.search(r"segment (\d+)", header)
    if m:
        return f"seg {m.group(1)} {SEGMENTS.get(int(m.group(1)), '')}"
    return header.split()[0]


def site_regex(site):
    return re.compile("(?=" + "".join(IUPAC[b] for b in site.upper()) + ")")


def digest(sequence, enzyme):
    cuts = [m.start() + enzyme["cut"] for m in site_regex(enzyme["site"]).finditer(sequence)]
    bounds = [0] + cuts + [len(sequence)]
    fragments = [bounds[i + 1] - bounds[i] for i in range(len(bounds) - 1)]
    return cuts, fragments


genomes = {name: read_fasta(path) for name, path in GENOMES.items()}

print("Restriction enzymes used:")
for e in ENZYMES:
    site = e["site"]
    print(f"  {e['name']:<8} recognition sequence {site:<8} cleavage site "
          f"{site[:e['cut']]}^{site[e['cut']:]}")

results = {}
for gname, segments in genomes.items():
    total = sum(len(s) for _, s in segments)
    print()
    print("=" * 78)
    print(f"{gname}: {len(segments)} segments, {total} nucleotides")
    for header, seq in segments:
        print(f"  {segment_name(header):<12}{len(seq):>6} nt   {header.split()[0]}")

    for e in ENZYMES:
        all_fragments = []
        total_cuts = 0
        print()
        print(f"  {e['name']} ({e['site'][:e['cut']]}^{e['site'][e['cut']:]})")
        for header, seq in segments:
            cuts, fragments = digest(seq, e)
            total_cuts += len(cuts)
            all_fragments.extend(fragments)
            if cuts:
                print(f"    {segment_name(header):<12} cleavages: {len(cuts)}   "
                      f"cut after nt: {cuts}   fragments: {fragments}")
            else:
                print(f"    {segment_name(header):<12} cleavages: 0   "
                      f"fragment: [{len(seq)}] (not cut)")
        print(f"    TOTAL cleavages: {total_cuts}   fragments: {len(all_fragments)}")
        results[(gname, e["name"])] = all_fragments

print()
print("=" * 78)
print("Comparison of the two versions:")
print(f"  {'Enzyme':<10}" + "".join(f"{g:>22}" for g in genomes))
for e in ENZYMES:
    row = f"  {e['name']:<10}"
    for g in genomes:
        cuts = len(results[(g, e['name'])]) - len(genomes[g])
        row += f"{str(cuts) + ' cleavages':>22}"
    print(row)


LADDER = [50, 100, 200, 300, 500, 700, 1000, 1500, 2000, 3000]
SMALL, LARGE = 50, 3000
D_SMALL, D_LARGE = 92.0, 10.0
GEL_LENGTH = 100


def distance(bp):
    bp = max(bp, 10)
    slope = (D_SMALL - D_LARGE) / (math.log10(LARGE) - math.log10(SMALL))
    return D_SMALL - slope * (math.log10(bp) - math.log10(SMALL))


lanes = [("Ladder", "", LADDER, "#f2c14e")]
for g, segments in genomes.items():
    lanes.append((g.split()[0], "uncut", [len(s) for _, s in segments], "#b8c4cc"))
colors = {list(genomes)[0]: "#7fe0ff", list(genomes)[1]: "#ff8fa3"}
for e in ENZYMES:
    for g in genomes:
        lanes.append((g.split()[0], e["name"], results[(g, e["name"])], colors[g]))

fig, ax = plt.subplots(figsize=(13, 8))
fig.patch.set_facecolor("#101418")
ax.set_facecolor("#1b2229")
w = 0.72
for x, (label, enzyme, sizes, color) in enumerate(lanes):
    ax.add_patch(plt.Rectangle((x - w / 2, -3), w, 2, color="#5a6570"))
    for bp in sizes:
        y = distance(bp)
        if y > GEL_LENGTH:
            continue
        alpha = 0.35 + 0.6 * min(1, bp / 1500)
        ax.add_patch(plt.Rectangle((x - w / 2, y - 0.5), w, 1.0,
                                   color=color, alpha=alpha))
    ax.text(x, -6.5, label, color="white", ha="center", fontsize=9)
    ax.text(x, -4.6, enzyme, color="#c9d1d9", ha="center", fontsize=8)

for bp in LADDER:
    ax.text(-0.6, distance(bp), f"{bp} bp", color="#f2c14e",
            ha="right", va="center", fontsize=8)

for i in range(len(ENZYMES)):
    x = 3 + 2 * i - 0.5
    ax.plot([x, x], [-3, GEL_LENGTH], color="#3a4550", linewidth=0.8)
ax.plot([2.5, 2.5], [-3, GEL_LENGTH], color="#3a4550", linewidth=0.8)

ax.set_xlim(-1.9, len(lanes) - 0.3)
ax.set_ylim(GEL_LENGTH, -9)
ax.axis("off")
ax.set_title("Restriction digest of two influenza A genomes (H9N2 2014 vs H5N1 2022)",
             color="white", fontsize=13)
plt.tight_layout()
plt.savefig("restriction_gel.png", dpi=150, facecolor=fig.get_facecolor())
plt.show()
