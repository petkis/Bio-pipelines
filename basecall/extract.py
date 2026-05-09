#!/usr/bin/env python3

with open("telomere_lengths.tsv", "w") as out:
    out.write("read_id\tlength\tstrand\n")

    def parse(file, strand):
        with open(file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) < 4:
                    continue
                try:
                    length = int(parts[2].replace("bp", ""))
                except ValueError:
                    continue
                out.write(f"{parts[0]}\t{length}\t{strand}\n")

    parse("telo_G.ncrf", "G")
    parse("telo_C.ncrf", "C")
