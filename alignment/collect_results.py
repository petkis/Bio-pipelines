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
        opposite_side = "PATERNAL" if side == "MATERNAL" else "MATERNAL"

        correct_chrom_correct_hap = 0
        correct_chrom_wrong_hap = 0
        wrong_chrom_correct_hap = 0
        wrong_chrom_wrong_hap = 0

        # Count columns like chr1_MATERNAL, chr1_PATERNAL, chrX_MATERNAL, etc.
        mapping_cols = [
            c for c in row.keys()
            if c.startswith("chr") and c not in {"chrom", "chromType"}
        ]

        for c in mapping_cols:
            count = row.get(c, 0)

            if c.endswith(f"_{side}"):
                same_hap = True
            else:
                same_hap = False

            if c.startswith(f"{chrom}_"):
                same_chrom = True
            else:
                same_chrom = False

            if same_chrom and same_hap:
                correct_chrom_correct_hap += count
            elif same_chrom and not same_hap:
                correct_chrom_wrong_hap += count
            elif not same_chrom and same_hap:
                wrong_chrom_correct_hap += count
            else:
                wrong_chrom_wrong_hap += count

        row["CorrectChromCorrectHap"] = correct_chrom_correct_hap
        row["CorrectChromWrongHap"] = correct_chrom_wrong_hap
        row["WrongChromCorrectHap"] = wrong_chrom_correct_hap
        row["WrongChromWrongHap"] = wrong_chrom_wrong_hap

        row["Percent_CorrectChromCorrectHap"] = (
            correct_chrom_correct_hap / row["N"]
        ) * 100.0 if row["N"] > 0 else 0

        row["Percent_CorrectChromWrongHap"] = (
            correct_chrom_wrong_hap / row["N"]
        ) * 100.0 if row["N"] > 0 else 0

        row["Percent_WrongChromCorrectHap"] = (
            wrong_chrom_correct_hap / row["N"]
        ) * 100.0 if row["N"] > 0 else 0

        row["Percent_WrongChromWrongHap"] = (
            wrong_chrom_wrong_hap / row["N"]
        ) * 100.0 if row["N"] > 0 else 0

        row["Percent"] = row["Percent_CorrectChromCorrectHap"]

        rows.append(row)

if rows:
    df = pd.DataFrame(rows)
    df.fillna(0).to_csv(OUTFILE, index=False)
else:
    print("XD")