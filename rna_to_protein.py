GENETIC_CODE = {
    # U first
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop",
    "UGU": "Cys", "UGC": "Cys", "UGA": "Stop", "UGG": "Trp",
    # C first
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    # A first
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    # G first
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

RNA = "GCAUGGCUUCAAAGGGUUGGCACUAAGC"

start = RNA.find("AUG")
protein = []

print("RNA =", RNA)

if start == -1:
    print("No start codon (AUG) found.")
else:
    print("Start codon AUG found at position", start + 1)
    print()
    print("Codon  ->  Amino acid")
    for i in range(start, len(RNA) - 2, 3):
        codon = RNA[i:i + 3]
        amino_acid = GENETIC_CODE[codon]
        print(codon, "  -> ", amino_acid)
        if amino_acid == "Stop":
            break
        protein.append(amino_acid)

    print()
    print("Protein =", "-".join(protein))
