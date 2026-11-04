import re
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

FASTA = "genome.fna"

IR_MIN = 10                  
INNER_MIN = 10               
INNER_MAX = 2500              
TSD_MIN, TSD_MAX = 3, 9       
STRONG_IR = 15                
CIRCULAR = False              

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def reverse_complement(seq):
    return "".join(COMPLEMENT.get(b, "N") for b in reversed(seq))


def pairs(a, b):
    return COMPLEMENT.get(a) == b


header = FASTA
genome = ""
with open(FASTA) as f:
    for line in f:
        if line.startswith(">"):
            header = line.strip()
        else:
            genome += line.strip().upper()
n = len(genome)
print(header)
print("Genome length:", n, "bases")


def base(i):
    """Base at position i (wraps around for a circular genome)."""
    if CIRCULAR:
        return genome[i % n]
    return genome[i] if 0 <= i < n else "N"


def region(i, length):
    return "".join(base(i + k) for k in range(length))

kmer_positions = {}
for i in range(n - IR_MIN + 1):
    kmer_positions.setdefault(genome[i:i + IR_MIN], []).append(i)


candidates = []
for i in range(n - IR_MIN + 1):
    left = genome[i:i + IR_MIN]
    for j in kmer_positions.get(reverse_complement(left), []):
        inner = j - (i + IR_MIN)
        if not (INNER_MIN <= inner <= INNER_MAX):
            continue
   
        if i > 0 and j + IR_MIN < n and pairs(genome[i - 1], genome[j + IR_MIN]):
            continue
        ir_len = IR_MIN
        while (i + ir_len < j - ir_len + IR_MIN - 1 - INNER_MIN
               and pairs(genome[i + ir_len], genome[j + IR_MIN - ir_len - 1])):
            ir_len += 1
        start, end = i, j + IR_MIN                       # [start, end)

        tsd = ""
        for t in range(TSD_MAX, TSD_MIN - 1, -1):
            if region(start - t, t) == region(end, t):
                tsd = region(end, t)
                break

        candidates.append({
            "start": start, "end": end, "length": end - start,
            "ir_len": ir_len, "ir": genome[i:i + ir_len],
            "inner": end - start - 2 * ir_len, "tsd": tsd,
        })

candidates.sort(key=lambda c: (c["ir_len"], len(c["tsd"])), reverse=True)

# expected number of random hits, to judge how meaningful the results are
expected_random = n * (INNER_MAX - INNER_MIN) / 4 ** IR_MIN
print(f"Candidates found: {len(candidates)} "
      f"(about {expected_random:.0f} expected just by chance)")


annotated = []
if GENBANK:
    with open(GENBANK) as f:
        for line in f:
            m = re.match(r"\s+(repeat_region|mobile_element)\s+(\d+)\.\.(\d+)", line)
            if m:
                annotated.append((int(m.group(2)) - 1, int(m.group(3))))
    print("Annotated repeat regions (NCBI):",
          [f"{s + 1}-{e}" for s, e in annotated])
else:
    print("No annotation file given - candidates cannot be checked against NCBI.")


def overlaps_annotation(c):
    return any(c["start"] < e and s < c["end"] for s, e in annotated)


def ends_match_annotation(c, tolerance=30):
    return any(abs(c["start"] - s) <= tolerance and abs(c["end"] - e) <= tolerance
               for s, e in annotated)


print()
print(f"{'Rank':<5}{'Start':>6}{'End':>7}{'Length':>8}{'IR':>4}  "
      f"{'Inverted repeat':<22}{'Direct rep.':<12}Note")
for rank, c in enumerate(candidates, start=1):
    note = "STRONG (too long to be chance)" if c["ir_len"] >= STRONG_IR else ""
    if ends_match_annotation(c):
        note = "<- matches annotated element"
    elif overlaps_annotation(c):
        note = "inside annotated region"
    print(f"{rank:<5}{c['start'] + 1:>6}{c['end']:>7}{c['length']:>8}"
          f"{c['ir_len']:>4}  {(c['ir'] if len(c['ir']) <= 20 else c['ir'][:17] + '...'):<22}{c['tsd'] or '-':<12}{note}")


fig, ax = plt.subplots(figsize=(12, 4))
ax.plot([1, n], [0, 0], color="black", linewidth=2)
for s, e in annotated:
    ax.add_patch(mpatches.Rectangle((s + 1, -0.6), e - s, 0.35,
                                    color="#e09f3e"))

for c in candidates:
    mid = (c["start"] + c["end"]) / 2
    width = c["end"] - c["start"]
    strong = ends_match_annotation(c) or c["ir_len"] >= STRONG_IR
    color = "#d1495b" if strong else "#2a7f9e"
    lw = 2.2 if strong else 0.9
    arc = mpatches.Arc((mid + 1, 0), width, width / n * 6, theta1=0,
                       theta2=180, color=color, linewidth=lw)
    ax.add_patch(arc)

ax.set_xlim(-100, n + 100)
ax.set_ylim(-0.8, 1.35)
ax.set_yticks([])
ax.set_xlabel("Position in genome")
ax.set_title(f"Inverted-repeat pairs (IR >= {IR_MIN} nt, {INNER_MIN}-{INNER_MAX} nt "
             f"between) in {header.lstrip('>').split(',')[0][:70]}")
handles = [
    mpatches.Patch(color="#2a7f9e", label=f"candidate (IR {IR_MIN}-{STRONG_IR - 1} nt, could be chance)"),
    mpatches.Patch(color="#d1495b", label=f"strong candidate (IR >= {STRONG_IR} nt)"),
]
if annotated:
    handles.append(mpatches.Patch(color="#e09f3e", label="annotated repeat element (NCBI)"))
ax.legend(handles=handles, loc="upper right", fontsize=8)
plt.tight_layout()
plt.savefig("transposon_candidates.png", dpi=150)
plt.show()
