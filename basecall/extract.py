#!/usr/bin/env python3

# Extract telomere lengths from NCRF output files.
# The script reads forward-strand and reverse-strand NCRF results,
# extracts read ID and telomere length, and writes them into one TSV file.


with open("telomere_lengths.tsv", "w") as out:
    out.write("read_id\tlength\tstrand\n")

    # Parse one NCRF output file and label its strand.
    def parse(file, strand):
        with open(file) as f:
            for line in f:
                parts = line.strip().split()

                # Skip malformed or too-short lines.
                if len(parts) < 4:
                    continue

                # Extract telomere length from column 3.
                # The value may contain "bp", so remove it before conversion.
                try:
                    length = int(parts[2].replace("bp", ""))
                except ValueError:
                    continue
                out.write(f"{parts[0]}\t{length}\t{strand}\n")

    parse("telo_G.ncrf", "G")
    parse("telo_C.ncrf", "C")
