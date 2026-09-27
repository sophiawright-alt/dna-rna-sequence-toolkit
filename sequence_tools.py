

DNA_BASES = set("ACGTN")
RNA_BASES = set("ACGUN")


def clean(seq):
    return "".join(seq.split()).upper()


def validate(seq, alphabet=DNA_BASES):
    """Raise ValueError if seq contains a letter that isn't in alphabet."""
    bad = set(seq) - alphabet
    if bad:
        raise ValueError(f"Unexpected character(s) in sequence: {sorted(bad)}")


# GC calculator
def gc_content(seq):
    seq = clean(seq)
    validate(seq)
    if len(seq) == 0:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq) * 100



# Reverse complement
COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def reverse_complement(seq):
    seq = clean(seq)
    validate(seq)
    complemented = [COMPLEMENT[base] for base in seq]
    return "".join(reversed(complemented))


# DNA -> RNA (transcription)
def transcribe(seq):
    seq = clean(seq)
    validate(seq)
    return seq.replace("T", "U")


# DNA -> protein (translation)

CODON_TABLE = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(seq, stop_at_stop=True):
    rna = transcribe(seq)
    protein = []
    for i in range(0, len(rna) - 2, 3):        # -2 so a partial codon is skipped
        codon = rna[i:i + 3]
        amino_acid = CODON_TABLE.get(codon, "X")
        if amino_acid == "*":
            if stop_at_stop:
                break
            protein.append("*")
        else:
            protein.append(amino_acid)
    return "".join(protein)


# Motif finder
def find_motif(seq, motif):
    seq = clean(seq)
    motif = clean(motif)
    validate(seq)
    validate(motif)
    if motif == "":
        return []

    positions = []
    start = seq.find(motif)
    while start != -1:
        positions.append(start + 1)          # +1 to convert 0-based -> 1-based
        start = seq.find(motif, start + 1)   # step by 1 to allow overlaps
    return positions


# Codon search (in-frame)
STOP_CODONS = ("TAA", "TAG", "TGA")
START_CODON = "ATG"


def normalise_codon(codon):
    """Accept a DNA or RNA codon (e.g. 'uag' or 'TAG') and return it as DNA."""
    codon = clean(codon).replace("U", "T")
    validate(codon)
    if len(codon) != 3:
        raise ValueError(f"A codon must be exactly 3 bases, got '{codon}'")
    return codon

# find codon in frame
def find_codon(seq, codon, reverse=True):
    seq = clean(seq)
    validate(seq)
    codon = normalise_codon(codon)
    length = len(seq)
    hits = {}

    for offset in range(3):
        frame = offset + 1
        hits[frame] = [i + 1 for i in range(offset, length - 2, 3)
                       if seq[i:i + 3] == codon]

    if reverse:
        rc = reverse_complement(seq)
        for offset in range(3):
            frame = -(offset + 1)
            # a codon at rc[i:i+3] covers forward bases (length-i-3)..(length-i-1)
            hits[frame] = sorted(length - i - 2 for i in range(offset, length - 2, 3)
                                 if rc[i:i + 3] == codon)
    return hits


# find stop codon
def find_stop_codons(seq, reverse=True):
    results = {}
    for stop in STOP_CODONS:
        for frame, positions in find_codon(seq, stop, reverse).items():
            results.setdefault(frame, [])
            results[frame].extend((pos, stop) for pos in positions)
    for frame in results:
        results[frame].sort()
    return results


def print_frame_hits(hits, label):
    total = sum(len(v) for v in hits.values())
    print(f"   {label}: {total} in-frame hit(s)")
    for frame in sorted(hits, key=lambda f: (f < 0, abs(f))):
        found = hits[frame]
        if not found:
            shown = "-"
        elif isinstance(found[0], tuple):          # (position, codon) pairs
            shown = ", ".join(f"{pos} ({codon})" for pos, codon in found)
        else:
            shown = ", ".join(str(pos) for pos in found)
        print(f"      frame {frame:+d} : {shown}")


# checker
EXAMPLE = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"


def self_test():
    """Run every function on known inputs. Raises AssertionError if wrong."""
    assert gc_content("GGCC") == 100.0
    assert gc_content("ATAT") == 0.0
    assert gc_content("") == 0.0
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AAAA") == "TTTT"
    assert transcribe("ATGC") == "AUGC"
    assert translate("ATGTTTTAA") == "MF"                     # Met-Phe-Stop
    assert translate("ATGTTTTAA", stop_at_stop=False) == "MF*"
    assert find_motif("ATATA", "ATA") == [1, 3]               # overlapping
    assert find_motif("ACGT", "TTT") == []

    # codon search
    assert find_codon("ATGTTTTAA", "TAA", reverse=False) == {1: [7], 2: [], 3: []}
    assert find_codon("ATGTTTTAA", "uaa", reverse=False)[1] == [7]   # RNA input ok
    assert find_codon("CTAAG", "TAA", reverse=False) == {1: [], 2: [2], 3: []}
    # "TTA" on the forward strand is "TAA" on the reverse strand
    assert find_codon("TTA", "TAA")[-1] == [1]
    assert find_stop_codons("ATGTAGTGA", reverse=False)[1] == [(4, "TAG"), (7, "TGA")]
    try:
        normalise_codon("AT")
        assert False, "short codon should raise"
    except ValueError:
        pass


# the report
def report(seq):
    seq = clean(seq)
    validate(seq)
    print("=" * 55)
    print("SEQUENCE REPORT")
    print("Sequence           :", seq)
    print("Length             :", len(seq), "bp")
    print()

    print("Base composition:")
    for base in "ATGC":
        n = seq.count(base)
        pct = (n / len(seq) * 100) if seq else 0.0
        print(f"   {base} : {n:>4}   ({pct:5.1f} %)")
    other = len(seq) - sum(seq.count(b) for b in "ATGC")
    if other:
        print(f"   other (e.g. N) : {other}")
    print()

    print("GC content         : {:.2f} %".format(gc_content(seq)))
    print("AT content         : {:.2f} %".format(100 - gc_content(seq)))
    print()

    print("Reverse complement :", reverse_complement(seq))
    print("Transcribed (RNA)  :", transcribe(seq))
    print("Protein (to stop)  :", translate(seq))
    print("Protein (full)     :", translate(seq, stop_at_stop=False))
    print()

    print("Codons (positions are 1-based, forward strand):")
    print_frame_hits(find_codon(seq, START_CODON), "Start codons (ATG)")
    print_frame_hits(find_stop_codons(seq), "Stop codons (TAA/TAG/TGA)")
    print("=" * 55)


def motif_search(seq):
    """Ask the user for motifs and report where they occur. Blank line quits."""
    seq = clean(seq)
    print("\nMotif finder  -  type a motif and press Enter (blank to quit)")
    while True:
        query = input("   motif > ")
        if query.strip() == "":
            print("   done.")
            break
        try:
            hits = find_motif(seq, query)
        except ValueError as error:
            print("  ", error)
            continue
        motif = clean(query)
        if hits:
            print(f"   '{motif}' found {len(hits)}x at position(s): {hits}")
        else:
            print(f"   '{motif}' not found")


def codon_search(seq):
    """Ask the user for codons and report in-frame hits. Blank line quits."""
    seq = clean(seq)
    print("\nCodon finder  -  type a codon (e.g. TAG or UAG), 'stop' for all")
    print("stop codons, or press Enter to quit")
    while True:
        query = input("   codon > ").strip()
        if query == "":
            print("   done.")
            break
        try:
            if query.lower() == "stop":
                print_frame_hits(find_stop_codons(seq), "Stop codons")
            else:
                codon = normalise_codon(query)
                print_frame_hits(find_codon(seq, codon), f"'{codon}'")
        except ValueError as error:
            print("  ", error)


if __name__ == "__main__":
    self_test()
    print("All self-tests passed.\n")

    entered = input("DNA sequence: ")
    sequence = entered if entered.strip() else EXAMPLE

    report(sequence)
    motif_search(sequence)
    codon_search(sequence)
