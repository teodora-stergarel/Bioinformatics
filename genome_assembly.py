"""
Bioinformatics Lab - Rebuilding a genome from random samples (shotgun assembly)

1. Read a whole viral genome from an NCBI .fna file.
2. Take 2000 random samples (reads) of 100-120 nucleotides and store them in a list.
3. Rebuild the genome only from the samples, by overlapping them (greedy assembly).
4. Compare the result with the original genome and show the main problem.
"""

import random
import matplotlib.pyplot as plt

FILE = "genome.fna"
NUM_SAMPLES = 2000
MIN_LEN, MAX_LEN = 100, 120
MIN_OVERLAP = 20          # minimum overlap accepted between two samples
random.seed(42)           # same random samples every run (remove for new ones)


# ---------------------------------------------------------------
# 1. Read the genome
# ---------------------------------------------------------------
genome = ""
with open(FILE) as f:
    for line in f:
        if not line.startswith(">"):
            genome += line.strip().upper()

print("Original genome length:", len(genome), "nucleotides")


# ---------------------------------------------------------------
# 2. Take 2000 random samples and store them in a list
# ---------------------------------------------------------------
samples = []
positions = []            # kept ONLY to check the result later, not used to rebuild
for _ in range(NUM_SAMPLES):
    length = random.randint(MIN_LEN, MAX_LEN)
    start = random.randint(0, len(genome) - length)
    samples.append(genome[start:start + length])
    positions.append((start, start + length))

total_bases = sum(len(s) for s in samples)
print("Samples taken:", len(samples))
print("Average coverage:", round(total_bases / len(genome), 2), "x")


# ---------------------------------------------------------------
# 3. Rebuild the genome from the samples (greedy overlap assembly)
# ---------------------------------------------------------------

# 3a. Remove duplicates and samples fully contained inside other samples
unique = sorted(set(samples), key=len, reverse=True)
reads = []
for s in unique:
    if not any(s in r for r in reads):
        reads.append(s)
print("Samples left after removing duplicates/contained ones:", len(reads))

# 3b. Index every sample by its first MIN_OVERLAP letters
prefix_index = {}
for j, r in enumerate(reads):
    prefix_index.setdefault(r[:MIN_OVERLAP], []).append(j)

# 3c. Find all overlaps: the end of sample A equals the beginning of sample B
edges = []                # (overlap length, A, B)
for i, a in enumerate(reads):
    for pos in range(1, len(a) - MIN_OVERLAP + 1):
        for j in prefix_index.get(a[pos:pos + MIN_OVERLAP], []):
            if j != i:
                ov = len(a) - pos
                if reads[j][:ov] == a[pos:]:
                    edges.append((ov, i, j))

# 3d. Greedy: join the pairs with the largest overlap first
edges.sort(reverse=True)
next_of, prev_of = {}, {}


def creates_cycle(a, b):
    node = b
    while node in next_of:
        node = next_of[node]
        if node == a:
            return True
    return node == a


overlap_of = {}
for ov, a, b in edges:
    if a not in next_of and b not in prev_of and not creates_cycle(a, b):
        next_of[a], prev_of[b] = b, a
        overlap_of[(a, b)] = ov

# 3e. Walk each chain and glue the samples together into contigs
contigs = []
for start_read in range(len(reads)):
    if start_read in prev_of:
        continue
    contig = reads[start_read]
    node = start_read
    while node in next_of:
        nxt = next_of[node]
        contig += reads[nxt][overlap_of[(node, nxt)]:]
        node = nxt
    contigs.append(contig)

contigs.sort(key=len, reverse=True)
rebuilt = contigs[0]


# ---------------------------------------------------------------
# 4. Compare with the original genome
# ---------------------------------------------------------------
print()
print("Number of pieces (contigs) rebuilt:", len(contigs))
print("Longest contig:", len(rebuilt), "nucleotides",
      f"({len(rebuilt) / len(genome) * 100:.1f}% of the genome)")
print("Contig lengths:", [len(c) for c in contigs])

correct = sum(1 for c in contigs if c in genome)
print("Contigs that match the original exactly:", correct, "of", len(contigs))
print("Genome rebuilt completely and correctly:", rebuilt == genome)

# Coverage: how many samples cover each position of the genome
coverage = [0] * len(genome)
for s, e in positions:
    for p in range(s, e):
        coverage[p] += 1
uncovered = [p for p, c in enumerate(coverage) if c == 0]

gaps = []
for p in uncovered:
    if gaps and p == gaps[-1][1] + 1:
        gaps[-1][1] = p
    else:
        gaps.append([p, p])

print()
print("Positions never sampled:", len(uncovered))
print("Gaps (regions with zero coverage):", len(gaps))
for s, e in gaps:
    print(f"  gap at positions {s + 1}-{e + 1} ({e - s + 1} nt)")

# Coverage chart
plt.figure(figsize=(12, 4.5))
plt.plot(range(1, len(genome) + 1), coverage, color="#2a7f9e", linewidth=0.7)
for s, e in gaps:
    plt.axvspan(s + 1, e + 1, color="#d1495b", alpha=0.6)
plt.title(f"Coverage of the genome by {NUM_SAMPLES} random samples "
          f"(red = gaps, never sampled)")
plt.xlabel("Position in genome")
plt.ylabel("Number of samples")
plt.tight_layout()
plt.savefig("coverage.png", dpi=150)
plt.show()
