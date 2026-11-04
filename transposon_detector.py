import random
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

random.seed(7)

IR = "GGTTCAAGCG"          
TSD_LEN = 5                
MIN_TE, MAX_TE = 30, 120   


def reverse_complement(seq):
    pairs = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(pairs[b] for b in reversed(seq))


IR_RIGHT = reverse_complement(IR)    


def random_dna(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def make_transposon(body_len):
    return IR + random_dna(body_len) + IR_RIGHT



host = random_dna(150)
while IR in host or IR_RIGHT in host:          
    host = random_dna(150)

dna = host
true_te = {}               

def insert(name, position, transposon):
    """Insert a transposon after `position`, duplicating the target site."""
    global dna
    tsd = dna[position:position + TSD_LEN]
    cut = position + TSD_LEN
    block = transposon + tsd                     
    dna = dna[:cut] + block + dna[cut:]

    for te in true_te.values():
        if te[0] >= cut:                         
            te[0] += len(block)
            te[1] += len(block)
        elif te[1] > cut:                        
            te[1] += len(block)
    true_te[name] = [cut, cut + len(transposon)]


insert("T1", 20, make_transposon(25))
insert("T2", 110, make_transposon(30))
insert("T3", true_te["T2"][0] + 20, make_transposon(15))
insert("T4", true_te["T1"][1] - 8, make_transposon(10))

print("Artificial DNA sequence length:", len(dna), "bases")
print(dna)
print()
print("Real positions (1-based, inclusive):")
for name, (s, e) in true_te.items():
    print(f"  {name}: start = {s + 1:>3}, end = {e:>3}, length = {e - s}")


def find_all(text, pattern):
    return [i for i in range(len(text) - len(pattern) + 1)
            if text[i:i + len(pattern)] == pattern]


left_irs = find_all(dna, IR)
right_irs = find_all(dna, IR_RIGHT)

detected = []        
unconfirmed = []     
for L in left_irs:
    for R in right_irs:
        start, end = L, R + len(IR_RIGHT)              
        if not (MIN_TE <= end - start <= MAX_TE):
            continue
        before = dna[start - TSD_LEN:start]
        after = dna[end:end + TSD_LEN]
        if len(before) == TSD_LEN and before == after:
            detected.append((start, end, before))
        else:
            unconfirmed.append((start, end))

paired_l = {s for s, _, _ in detected}
paired_r = {e - len(IR_RIGHT) for _, e, _ in detected}
orphan_l = [p for p in left_irs if p not in paired_l]
orphan_r = [p for p in right_irs if p not in paired_r]

print()
print("Left inverted repeats found at :", [p + 1 for p in left_irs])
print("Right inverted repeats found at:", [p + 1 for p in right_irs])
print()
print("DETECTED transposons (inverted repeats + direct repeats):")
for s, e, tsd in sorted(detected):
    match = [n for n, (ts, te) in true_te.items() if ts == s and te == e]
    label = match[0] if match else "FALSE POSITIVE"
    print(f"  start = {s + 1:>3}, end = {e:>3}, length = {e - s:>3}, "
          f"direct repeat = {tsd}  ->  {label}")

print()
print("Inverted-repeat pairs rejected (no matching direct repeats):")
for s, e in sorted(unconfirmed):
    print(f"  {s + 1}-{e}")
print("Left IRs with no partner :", [p + 1 for p in orphan_l])
print("Right IRs with no partner:", [p + 1 for p in orphan_r])

found_names = {n for n, (ts, te) in true_te.items()
               if any(ts == s and te == e for s, e, _ in detected)}
print()
for name in true_te:
    print(f"  {name}: {'found' if name in found_names else 'MISSED'}")


colors = {"T1": "#2a7f9e", "T2": "#e09f3e", "T3": "#9e2a2b", "T4": "#540b0e"}
fig, ax = plt.subplots(figsize=(12, 3.6))
ax.plot([1, len(dna)], [0, 0], color="black", linewidth=2)

for row, (name, (s, e)) in enumerate(true_te.items()):
    ax.add_patch(mpatches.Rectangle((s + 1, 1 + row * 0.5), e - s, 0.4,
                                    color=colors[name]))
    ax.text(e + 3, 1.2 + row * 0.5, name, va="center", fontsize=9)

for s, e, _ in detected:
    level = sum(1 for s2, e2, _ in detected if s2 <= s and e <= e2 and (s2, e2) != (s, e))
    ax.add_patch(mpatches.Rectangle((s + 1, -0.9 - level * 0.5), e - s, 0.4,
                                    color="#3a9e5c"))
for p in orphan_l + orphan_r:
    ax.plot(p + 1, -1.85, marker="x", color="#d1495b", markersize=8)

ax.text(-5, 1.9, "Real", ha="right", va="center")
ax.text(-5, -0.95, "Detected", ha="right", va="center")
ax.text(-5, -1.85, "Unpaired IR", ha="right", va="center", color="#d1495b")
ax.set_xlim(-40, len(dna) + 20)
ax.set_ylim(-2.3, 3.3)
ax.set_yticks([])
ax.set_xlabel("Position in DNA sequence")
ax.set_title("Transposable elements: real positions vs. detected positions")
plt.tight_layout()
plt.savefig("transposons.png", dpi=150)
plt.show()
