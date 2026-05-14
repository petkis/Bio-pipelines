import pandas as pd
from pathlib import Path
import os

# Combine per-chromosome alignment summaries into one CSV file.
# The script searches result folders, reads per_chrom_all.txt files,
# calculates the percentage of reads mapped to the expected same-side target,
# and saves the combined table.

ROOT = Path(".")
OUTFILE = Path("combined_results.csv")
OUTFILE_AGG = Path("combined_results_grouped.csv")

# Extract metadata from a result folder name.
# Expected information includes chromosome, haplotype side,
# reference type, number of reads, and read length.

def parse_folder_name(name: str):
    toks = name.split("_")

    chrom = next((t for t in toks if t.startswith("chr")), None)
    side = next((t for t in toks if t in ("MATERNAL", "PATERNAL")), None)
    
    # Extract reference subtype and simulation parameters from fixed tokens.
    chrom_type = toks[0]
    ref = toks[toks.index("to") + 1] if "to" in toks else None
    n_val = int(toks[toks.index("N") + 1])
    length = int(toks[toks.index("L") + 1])

    return chrom, side, ref, chrom_type, n_val, length


rows = []

# Search all processing result directories.
for proc_dir in ROOT.glob("*_PROC_RESULTS_*"):
    tool_name = proc_dir.name.split("_")[0]

    # Find every per_chrom_all.txt file inside the result directory.
    for txt in proc_dir.rglob("per_chrom_all.txt"):
        folder = txt.parent.name

        # Parse metadata from the parent folder name.
        chrom, side, ref, chrom_type, n_val, length = parse_folder_name(folder)

        # Skip folders where chromosome or haplotype side could not be detected.
        if not chrom or not side:
            continue

        row = {
            "tool": tool_name,
            "chrom": chrom,
            "side": side,
            "REF": ref,
            "chromType": chrom_type,
            "N": n_val,
            "Length": length,
        }

        # Read per-chromosome mapping counts and add them as columns.
        with txt.open("r") as f:
            for line in f:
                p = line.strip().split()
                if len(p) == 2:
                    row[p[0]] = int(p[1])

        target_key = f"{chrom}_{side}"

        correct_reads = row.get(target_key, 0)
        row["Percent"] = (correct_reads / row["N"]) * 100.0 if row["N"] > 0 else 0

        rows.append(row)

if rows:
    df = pd.DataFrame(rows)
    df.fillna(0).to_csv(OUTFILE, index=False)
