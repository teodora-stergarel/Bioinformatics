FILE = "sequence.fna"
UNIT_MIN, UNIT_MAX = 3, 5
TOTAL_MIN, TOTAL_MAX = 6, 10


def is_primitive(unit):
    for p in range(1, len(unit)):
        if len(unit) % p == 0 and unit[:p] * (len(unit) // p) == unit:
            return False
    return True


sequence = ""
with open(FILE) as f:
    for line in f:
        if line.startswith(">"):
            print(line.strip())
        else:
            sequence += line.strip().upper()
n = len(sequence)
print("Sequence length:", n, "nucleotides")


repeats = []
too_long = []

for k in range(UNIT_MIN, UNIT_MAX + 1):
    i = 0
    while i + 2 * k <= n:
        unit = sequence[i:i + k]
        copies = 1
        while sequence[i + copies * k:i + (copies + 1) * k] == unit:
            copies += 1

        is_start = i < k or sequence[i - k:i] != unit
        if copies >= 2 and is_start and is_primitive(unit):
            total = copies * k
            item = {"start": i, "end": i + total, "unit": unit,
                    "unit_len": k, "copies": copies, "total": total}
            if TOTAL_MIN <= total <= TOTAL_MAX:
                repeats.append(item)
            elif total > TOTAL_MAX:
                too_long.append(item)
            i += total - k + 1
        else:
            i += 1

repeats.sort(key=lambda r: r["start"])


print()
print(f"Repetitions found ({TOTAL_MIN}-{TOTAL_MAX} bases, unit of "
      f"{UNIT_MIN}-{UNIT_MAX} bases): {len(repeats)}")
print()
print(f"{'#':<4}{'Start':>6}{'End':>6}  {'Unit':<6}{'Copies':>7}{'Total':>7}   Context")
for idx, r in enumerate(repeats, start=1):
    s, e = r["start"], r["end"]
    context = (sequence[max(0, s - 5):s].lower() + sequence[s:e]
               + sequence[e:e + 5].lower())
    print(f"{idx:<4}{s + 1:>6}{e:>6}  {r['unit']:<6}{r['copies']:>7}{r['total']:>7}   {context}")

print()
print("Summary by unit length:")
for k in range(UNIT_MIN, UNIT_MAX + 1):
    group = [r for r in repeats if r["unit_len"] == k]
    print(f"  unit of {k} bases: {len(group)} repetitions")

counts = {}
for r in repeats:
    counts[r["unit"]] = counts.get(r["unit"], 0) + 1
common = sorted(counts.items(), key=lambda x: x[1], reverse=True)[:5]
print("Most frequent repeated units:", ", ".join(f"{u} ({c}x)" for u, c in common))

if too_long:
    print()
    print(f"Repetitions longer than {TOTAL_MAX} bases (not counted above):")
    for r in too_long:
        print(f"  {r['start'] + 1}-{r['end']}: {r['unit']} x {r['copies']} "
              f"= {r['total']} bases")
