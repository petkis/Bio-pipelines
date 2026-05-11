import pandas as pd
from pathlib import Path
import os

ROOT = Path(".")
OUTFILE = Path("combined_results.csv")
OUTFILE_AGG = Path("combined_results_grouped.csv")

MODE = os.environ.get("MODE", "same").lower()

if MODE not in ["same", "cross"]:
    raise ValueError("MODE must be 'same' or 'cross'")

print(f"Running in MODE = {MODE}")

def parse_folder_name(name: str):
    toks = name.split("_")

    chrom = next((t for t in toks if t.startswith("chr")), None)
    side = next((t for t in toks if t in ("MATERNAL", "PATERNAL")), None)

    chrom_type = toks[0]
    ref = toks[toks.index("to") + 1] if "to" in toks else None
    n_val = int(toks[toks.index("N") + 1])
    length = int(toks[toks.index("L") + 1])

    return chrom, side, ref, chrom_type, n_val, length


rows = []

for proc_dir in ROOT.glob("*_PROC_RESULTS_*"):
    tool_name = proc_dir.name.split("_")[0]

    for txt in proc_dir.rglob("per_chrom_all.txt"):
        folder = txt.parent.name

        chrom, side, ref, chrom_type, n_val, length = parse_folder_name(folder)

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

        with txt.open("r") as f:
            for line in f:
                p = line.strip().split()
                if len(p) == 2:
                    row[p[0]] = int(p[1])

        # ----------------------------
        # SWITCH LOGIC
        # ----------------------------
        if MODE == "cross":
            opposite = "MATERNAL" if side == "PATERNAL" else "PATERNAL"

            if chrom == "chrX":
                target_key = f"chrY_{opposite}"
            elif chrom == "chrY":
                target_key = f"chrX_{opposite}"
            else:
                target_key = f"{chrom}_{opposite}"

        else:
            target_key = f"{chrom}_{side}"

        correct_reads = row.get(target_key, 0)
        row["Percent"] = (correct_reads / row["N"]) * 100.0 if row["N"] > 0 else 0
        row["mode"] = MODE

        rows.append(row)

if rows:
    df = pd.DataFrame(rows)
    df.fillna(0).to_csv(OUTFILE, index=False)
