import math
import random
import matplotlib.pyplot as plt

FILE = "sequence.fna"
NUM_SAMPLES = 10
MIN_LEN, MAX_LEN = 100, 3000
GEL_LENGTH = 100
random.seed(5)


sequence = ""
with open(FILE) as f:
    for line in f:
        if line.startswith(">"):
            print(line.strip())
        else:
            sequence += line.strip().upper()
print("Sequence length:", len(sequence), "nucleotides")


max_len = min(MAX_LEN, len(sequence))
samples = []
for _ in range(NUM_SAMPLES):
    length = random.randint(MIN_LEN, max_len)
    start = random.randint(0, len(sequence) - length)
    samples.append(sequence[start:start + length])


SMALL, LARGE = 100, 3000
D_SMALL, D_LARGE = 0.90 * GEL_LENGTH, 0.15 * GEL_LENGTH


def migration_distance(bp):
    slope = (D_SMALL - D_LARGE) / (math.log10(LARGE) - math.log10(SMALL))
    return D_SMALL - slope * (math.log10(bp) - math.log10(SMALL))


print()
print(f"{'Sample':<8}{'Length (bp)':>12}{'Distance':>10}   First bases")
for i, s in enumerate(samples, start=1):
    print(f"{i:<8}{len(s):>12}{migration_distance(len(s)):>10.1f}   {s[:20]}...")

ranking = sorted(range(NUM_SAMPLES), key=lambda i: len(samples[i]))
print()
print("Order on the gel (fastest -> slowest):",
      ", ".join(f"S{i + 1} ({len(samples[i])} bp)" for i in ranking))


LADDER = [100, 200, 300, 500, 700, 1000, 1500, 2000, 3000]
lanes = [("Ladder", LADDER)] + \
        [(f"S{i + 1}", [len(s)]) for i, s in enumerate(samples)] + \
        [("All", [len(s) for s in samples])]

fig, ax = plt.subplots(figsize=(11, 7))
fig.patch.set_facecolor("#101418")
ax.set_facecolor("#1b2229")

lane_width = 0.7
for x, (name, sizes) in enumerate(lanes):
    ax.add_patch(plt.Rectangle((x - lane_width / 2, -3), lane_width, 2,
                               color="#5a6570"))
    for bp in sizes:
        y = migration_distance(bp)
        color = "#f2c14e" if name == "Ladder" else "#7fe0ff"
        ax.add_patch(plt.Rectangle((x - lane_width / 2, y - 0.6),
                                   lane_width, 1.2, color=color, alpha=0.9))
        if name not in ("Ladder", "All"):
            ax.text(x, y + 2.8, f"{bp}", color="white", ha="center", fontsize=8)
    ax.text(x, -5, name, color="white", ha="center", fontsize=9)

for bp in LADDER:
    ax.text(-0.6, migration_distance(bp), f"{bp} bp", color="#f2c14e",
            ha="right", va="center", fontsize=8)

ax.annotate("", xy=(len(lanes) - 0.3, GEL_LENGTH - 3),
            xytext=(len(lanes) - 0.3, 3),
            arrowprops=dict(arrowstyle="->", color="white"))
ax.text(len(lanes) - 0.1, GEL_LENGTH / 2, "migration  (-  to  +)",
        color="white", rotation=-90, va="center", fontsize=9)

ax.set_xlim(-1.8, len(lanes) + 0.3)
ax.set_ylim(GEL_LENGTH, -8)
ax.axis("off")
ax.set_title("Simulated gel electrophoresis of 10 random DNA samples",
             color="white", fontsize=13)
plt.tight_layout()
plt.savefig("gel_electrophoresis.png", dpi=150, facecolor=fig.get_facecolor())
plt.show()
