# DNA/RNA Sequence Toolkit

A small Python toolkit for common DNA/RNA sequence analysis tasks, built from scratch using only the Python standard library (no external dependencies).

## Features
- **GC content** calculation
- **Reverse complement** generation
- **Transcription** (DNA → RNA)
- **Translation** (RNA → protein, using the standard codon table)
- **Motif searching**, including overlapping matches
- **Codon search** across all six reading frames (three forward, three on the reverse strand)
- **Stop codon finder** that locates every TAA, TAG and TGA in frame
- Built-in **self-tests** to verify correctness on every run
- Interactive command-line report: enter a sequence and get a full breakdown

## Usage
```bash
python sequence_tools.py
```
Enter a DNA sequence when prompted (or press Enter to use the built-in example). The script runs its self-tests, then prints a report — base composition, GC%, reverse complement, RNA transcript, and translated protein — and lets you search for motifs interactively.

## Codon search
Unlike the motif search, which matches at any position, the codon search only counts codons that line up with a reading frame's triplets. This makes it useful for spotting start and stop codons and identifying possible open reading frames.

- Codons can be entered as DNA or RNA (`TAG` or `UAG`)
- Type `stop` to list every stop codon
- Positions are 1-based and given on the forward strand (pointing at the codon's first base), so hits from different frames can be compared directly

Example output for the built-in sequence:
```
Stop codons (TAA/TAG/TGA): 3 in-frame hit(s)
   frame +1 : 22 (TGA), 37 (TAG)
   frame +2 : 11 (TAA)
   frame +3 : -
   frame -1 : -
   frame -2 : -
   frame -3 : -
```

## Why I built this
With a background in biology and completing an MSc in Bioinformatics, I wanted a project that turned concepts I already understood, like transcription, translation, and GC content, into working code, as a way to build confidence with Python and Git before moving on to larger pipelines. It's already proven useful for quickly breaking down a DNA sequence, checking its GC content, finding motifs, and generating a reverse complement, without needing to open a full bioinformatics suite. The same logic could be extended into a real pipeline step, for example as a quick pre-processing or quality-check stage before alignment or variant calling.
