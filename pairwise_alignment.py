import matplotlib.pyplot as plt

X = "GATTACAGCTTAGCCGTA"
Y = "CCGATTACATTAGCCTT"
MATCH = 1
MISMATCH = -1
GAP_COST = -2


def score(a, b):
    return MATCH if a == b else MISMATCH


def global_alignment(x, y, score, gap_cost):
    m, n = len(x), len(y)
    F = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        F[i][0] = F[i - 1][0] + gap_cost
    for j in range(1, n + 1):
        F[0][j] = F[0][j - 1] + gap_cost

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = F[i - 1][j - 1] + score(x[i - 1], y[j - 1])
            up = F[i - 1][j] + gap_cost
            left = F[i][j - 1] + gap_cost
            F[i][j] = max(diag, up, left)

    i, j = m, n
    a, b, path = "", "", [(i, j)]
    while i > 0 or j > 0:
        if i > 0 and j > 0 and F[i][j] == F[i - 1][j - 1] + score(x[i - 1], y[j - 1]):
            a, b = x[i - 1] + a, y[j - 1] + b
            i, j = i - 1, j - 1
        elif i > 0 and F[i][j] == F[i - 1][j] + gap_cost:
            a, b = x[i - 1] + a, "-" + b
            i -= 1
        else:
            a, b = "-" + a, y[j - 1] + b
            j -= 1
        path.append((i, j))
    return a, b, F[m][n], F, path


def local_alignment(x, y, score, gap_cost):
    m, n = len(x), len(y)
    H = [[0] * (n + 1) for _ in range(m + 1)]
    best_score, best_i, best_j = 0, 0, 0

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = H[i - 1][j - 1] + score(x[i - 1], y[j - 1])
            up = H[i - 1][j] + gap_cost
            left = H[i][j - 1] + gap_cost
            H[i][j] = max(0, diag, up, left)
            if H[i][j] > best_score:
                best_score, best_i, best_j = H[i][j], i, j

    i, j = best_i, best_j
    a, b, path = "", "", [(i, j)]
    while i > 0 and j > 0 and H[i][j] > 0:
        if H[i][j] == H[i - 1][j - 1] + score(x[i - 1], y[j - 1]):
            a, b = x[i - 1] + a, y[j - 1] + b
            i, j = i - 1, j - 1
        elif H[i][j] == H[i - 1][j] + gap_cost:
            a, b = x[i - 1] + a, "-" + b
            i -= 1
        else:
            a, b = "-" + a, y[j - 1] + b
            j -= 1
        path.append((i, j))
    return a, b, best_score, H, path


def match_line(a, b):
    return "".join("|" if p == q else (" " if "-" in (p, q) else ".")
                   for p, q in zip(a, b))


def print_matrix(M, x, y, name):
    print(f"Matrix {name}:")
    print("       " + "".join(f"{c:>4}" for c in "-" + y))
    for i, row in enumerate(M):
        label = "-" if i == 0 else x[i - 1]
        print(f"   {label:>2} " + "".join(f"{v:>4}" for v in row))
    print()


def print_alignment(title, a, b, best):
    identity = sum(1 for p, q in zip(a, b) if p == q)
    print(title)
    print("  " + a)
    print("  " + match_line(a, b))
    print("  " + b)
    print(f"  score = {best}, length = {len(a)}, identity = {identity}/{len(a)} "
          f"({identity / len(a) * 100:.1f}%)")
    print()


def plot_matrix(ax, M, x, y, path, title):
    im = ax.imshow(M, cmap="viridis")
    ax.plot([p[1] for p in path], [p[0] for p in path], color="#ff4d4d",
            linewidth=2, marker="o", markersize=3)
    ax.set_xticks(range(len(y) + 1))
    ax.set_xticklabels(["-"] + list(y))
    ax.set_yticks(range(len(x) + 1))
    ax.set_yticklabels(["-"] + list(x))
    ax.xaxis.tick_top()
    ax.set_title(title, pad=24)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


print(f"X = {X} ({len(X)})")
print(f"Y = {Y} ({len(Y)})")
print(f"match = {MATCH}, mismatch = {MISMATCH}, gap = {GAP_COST}")
print()

ga, gb, g_score, F, g_path = global_alignment(X, Y, score, GAP_COST)
la, lb, l_score, H, l_path = local_alignment(X, Y, score, GAP_COST)

print_matrix(F, X, Y, "F (global, Needleman-Wunsch)")
print_alignment("Global alignment (Needleman-Wunsch):", ga, gb, g_score)
print_matrix(H, X, Y, "H (local, Smith-Waterman)")
print_alignment("Local alignment (Smith-Waterman):", la, lb, l_score)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
plot_matrix(ax1, F, X, Y, g_path, f"Global (Needleman-Wunsch), score = {g_score}")
plot_matrix(ax2, H, X, Y, l_path, f"Local (Smith-Waterman), score = {l_score}")
plt.tight_layout()
plt.savefig("alignment_matrices.png", dpi=150)
plt.show()
