import matplotlib.pyplot as plt

S = "CGGACTGATCTATCT" + "A" * 27 + "CGTAGCATCTATCGATCTATCTAGCGATCTA"
WINDOW = 30
STEP = 1
PROMOTER_FILE = "gene_promoters_complete.fa"
PROMOTER_INDEX = 0
NUM_PROMOTERS = 1000


def cg_content(seq):
    seq = seq.upper()
    a, c, g, t = seq.count("A"), seq.count("C"), seq.count("G"), seq.count("T")
    return round(100 / (a + c + g + t) * (c + g), 2)


def information_content(seq):
    n = len(seq)
    total = 0
    for u in range(1, n):
        matches = sum(1 for i in range(n - u) if seq[i] == seq[i + u])
        total += matches / (n - u) * 100
    return round(total / (n - 1), 2)


def pattern(seq, window=WINDOW, step=STEP):
    points = []
    for i in range(0, len(seq) - window + 1, step):
        w = seq[i:i + window]
        points.append((i + 1, cg_content(w), information_content(w)))
    return points


def center_of_weight(points):
    xs = [p[1] for p in points]
    ys = [p[2] for p in points]
    return sum(xs) / len(xs), sum(ys) / len(ys)


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


def plot_pattern(points, title, filename):
    pos = [p[0] for p in points]
    cg = [p[1] for p in points]
    ic = [p[2] for p in points]
    cx, cy = center_of_weight(points)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    ax1.plot(pos, cg, color="#2a7f9e", label="C+G %")
    ax1.plot(pos, ic, color="#d1495b", label="Information content (IC)")
    ax1.set_xlabel("Sliding window start position")
    ax1.set_ylabel("Value (%)")
    ax1.set_title("Signals along the sequence")
    ax1.legend()
    ax1.grid(alpha=0.3)

    ax2.scatter(cg, ic, s=14, color="#1b2229", alpha=0.6, label="windows")
    ax2.scatter([cx], [cy], s=160, color="#d1495b", marker="X",
                label=f"center of weight ({cx:.2f}, {cy:.2f})")
    ax2.set_xlim(0, 100)
    ax2.set_ylim(0, 100)
    ax2.set_xlabel("C+G %")
    ax2.set_ylabel("Information content (IC)")
    ax2.set_title("DNA pattern (objective digital stain)")
    ax2.legend(loc="upper right")
    ax2.grid(alpha=0.3)

    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()


print("Test sequence S (length", len(S), "):")
print(S)
print()
print("Whole sequence:  CG =", cg_content(S), "  IC =", information_content(S))

points_s = pattern(S)
print(f"Sliding window {WINDOW} bases -> {len(points_s)} windows")
print()
print(f"{'Window':<8}{'Position':>9}{'CG %':>9}{'IC':>9}   Content")
for k, (p, cg, ic) in enumerate(points_s, start=1):
    print(f"{k:<8}{p:>9}{cg:>9.2f}{ic:>9.2f}   {S[p - 1:p - 1 + WINDOW]}")

cx, cy = center_of_weight(points_s)
print()
print(f"Center of weight of the pattern: CG = {cx:.2f}, IC = {cy:.2f}")
plot_pattern(points_s, "Test sequence S", "pattern_S.png")

promoters = read_fasta(PROMOTER_FILE)
selected = promoters[:NUM_PROMOTERS]
print()
print(f"Promoters in file: {len(promoters)}, analysed: {len(selected)}")

centers = []
for k, (header, seq) in enumerate(selected, start=1):
    centers.append(center_of_weight(pattern(seq)))
    if k % 100 == 0:
        print(f"  processed {k} promoters")

plt.figure(figsize=(7, 6))
plt.scatter([c[0] for c in centers], [c[1] for c in centers],
            s=10, color="#2a7f9e", alpha=0.5)
plt.xlim(0, 100)
plt.ylim(0, 100)
plt.xlabel("C+G % (center of weight)")
plt.ylabel("Information content (center of weight)")
plt.title(f"General distribution of pattern centers ({len(centers)} promoters)")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("centers_distribution.png", dpi=150)
plt.close()

header, seq = promoters[PROMOTER_INDEX]
points_p = pattern(seq)
px, py = center_of_weight(points_p)
print()
print("Promoter used for comparison with PromKappa:")
print(" ", header, f"({len(seq)} bp, {len(points_p)} windows)")
print(f"  center of weight: CG = {px:.2f}, IC = {py:.2f}")
plot_pattern(points_p, header, "promoter_pattern.png")

print()
print("Charts saved: pattern_S.png, centers_distribution.png, promoter_pattern.png")
