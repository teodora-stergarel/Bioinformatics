from itertools import product
import matplotlib.pyplot as plt

FILE = "genome.fna"

# 1. Read the genome (skip header lines that start with '>')
genome = ""
with open(FILE) as f:
    for line in f:
        if not line.startswith(">"):
            genome += line.strip().upper()

print("Genome length:", len(genome), "nucleotides")

# 2. Prepare all 64 codons with count 0
codons = {"".join(c): 0 for c in product("ACGT", repeat=3)}

# 3. Sliding window of 3 positions, moving 1 letter at a time
total = 0
for i in range(len(genome) - 2):
    window = genome[i:i + 3]
    if window in codons:          # skip windows with N or other letters
        codons[window] += 1
        total += 1

print("Total windows (codons) counted:", total)
print()

# 4. Sort by frequency and print all 64
ranking = sorted(codons.items(), key=lambda x: x[1], reverse=True)

print("Rank  Codon  Count  Frequency")
for rank, (codon, count) in enumerate(ranking, start=1):
    print(f"{rank:<5} {codon:<6} {count:<6} {count / total * 100:.2f}%")

# 5. Chart with the top 10 codons
top10 = ranking[:10]
names = [c for c, _ in top10]
values = [n / total * 100 for _, n in top10]

plt.figure(figsize=(10, 5.5))
bars = plt.bar(names, values, color="#2a7f9e")
for bar, (codon, count) in zip(bars, top10):
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
             f"{count}", ha="center", va="bottom", fontsize=9)
plt.title("Top 10 codons by frequency (sliding window of 3, step 1)")
plt.xlabel("Codon")
plt.ylabel("Frequency (%)")
plt.tight_layout()
plt.savefig("top10_codons.png", dpi=150)
plt.show()
