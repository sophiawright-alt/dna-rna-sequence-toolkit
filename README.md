# DNA/RNA Sequence Toolkit

A small Python toolkit for common DNA/RNA sequence analysis tasks, built from scratch using only the Python standard library (no external dependencies).

## Features
- **GC content** calculation
- **Reverse complement** generation
- **Transcription** (DNA → RNA)
- **Translation** (RNA → protein, using the standard codon table)
- **Motif searching**, including overlapping matches
- Built-in **self-tests** to verify correctness on every run
- Interactive command-line report: enter a sequence and get a full breakdown

## Usage
```bash
python sequence_tools.py
```
Enter a DNA sequence when prompted (or press Enter to use the built-in example). The script runs its self-tests, then prints a report — base composition, GC%, reverse complement, RNA transcript, and translated protein — and lets you search for motifs interactively.

## Why I built this
With a background in biology and completing an MSc in Bioinformatics, I wanted a project that turned concepts I already understood, like transcription, translation, and GC content, into working code, as a way to build confidence with Python and Git before moving on to larger pipelines. It's already proven useful for quickly breaking down a DNA sequence, checking its GC content, finding motifs, and generating a reverse complement, without needing to open a full bioinformatics suite. The same logic could be extended into a real pipeline step, for example as a quick pre-processing or quality-check stage before alignment or variant calling.
