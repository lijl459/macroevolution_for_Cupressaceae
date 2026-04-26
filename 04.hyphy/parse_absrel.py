#!/usr/bin/env python3
import json
import glob
import os
import pandas as pd
from ete3 import Tree

ABSREL_DIR = "03.absrel"
TREE_DIR   = "01.tree_rr"
P_CUTOFF   = 0.05

def to_float(x):
    """Best-effort convert json value to float; return None if impossible."""
    if x is None:
        return None
    if isinstance(x, (int, float)):
        return float(x)
    if isinstance(x, str):
        try:
            return float(x)
        except Exception:
            return None
    # some HyPhy json may store p-values as dicts; try common keys
    if isinstance(x, dict):
        for k in ["value", "p", "P", "p-value", "P-value"]:
            if k in x:
                return to_float(x[k])
        return None
    return None

def get_pvalue(info: dict):
    """Try multiple possible p-value fields across HyPhy versions."""
    # Most common (Holm-Bonferroni corrected)
    for key in ["Corrected P-value", "Corrected p-value", "corrected p-value", "corrected_pvalue"]:
        if key in info:
            return to_float(info.get(key))

    # Sometimes uncorrected p-values exist under different names
    for key in ["P-value", "p-value", "pvalue", "p"]:
        if key in info:
            return to_float(info.get(key))

    # Sometimes nested under LRT test object
    # e.g., "LRT": {"p-value": 0.001}
    for key in ["LRT", "LRT test", "test", "Test"]:
        if key in info and isinstance(info[key], dict):
            for k2 in ["Corrected P-value", "P-value", "p-value", "p"]:
                if k2 in info[key]:
                    return to_float(info[key][k2])

    return None

def get_omega_list(info: dict):
    """Extract omega(s) from possible fields."""
    omegas = []

    # standard aBSREL key
    if "omega distribution" in info and isinstance(info["omega distribution"], list):
        for d in info["omega distribution"]:
            if isinstance(d, dict) and "omega" in d:
                o = to_float(d["omega"])
                if o is not None:
                    omegas.append(o)

    # fallback: sometimes stored under different key
    for k in ["omegas", "omega", "Omega", "ω"]:
        if k in info:
            val = info[k]
            if isinstance(val, list):
                for x in val:
                    o = to_float(x)
                    if o is not None:
                        omegas.append(o)
            else:
                o = to_float(val)
                if o is not None:
                    omegas.append(o)

    # de-duplicate (optional)
    return omegas

records = []

# -----------------------------
# 1) parse jsons
# -----------------------------
for jf in glob.glob(f"{ABSREL_DIR}/*.json"):
    gene = os.path.basename(jf).replace(".json", "")

    if os.path.getsize(jf) == 0:
        print(f"[WARN] Empty json file skipped: {jf}")
        continue

    try:
        with open(jf) as f:
            data = json.load(f)
    except Exception as e:
        print(f"[WARN] Failed to parse json: {jf} ({e})")
        continue

    branches = data.get("branch attributes", {})
    if not isinstance(branches, dict) or len(branches) == 0:
        print(f"[WARN] No branch attributes found: {jf}")
        continue

    for br, info in branches.items():
        if not isinstance(info, dict):
            continue

        p = get_pvalue(info)
        omegas = get_omega_list(info)

        records.append({
            "gene": gene,
            "branch": br,
            "corrected_p": p,              # may be None if not found
            "has_omega_gt1": any(o > 1 for o in omegas),
            "max_omega": max(omegas) if omegas else None
        })

if len(records) == 0:
    raise RuntimeError("No valid aBSREL json files were parsed. Check your inputs.")

df = pd.DataFrame(records)

# ensure numeric
df["corrected_p"] = pd.to_numeric(df["corrected_p"], errors="coerce")

# selected branch definition
df["is_selected"] = (df["corrected_p"].notna()) & (df["corrected_p"] < P_CUTOFF)

df.to_csv("absrel_branch_level.tsv", sep="\t", index=False)

# -----------------------------
# 2) per-gene summary
# -----------------------------
gene_sum = (
    df[df["is_selected"]]
    .groupby("gene")
    .size()
    .reset_index(name="n_selected_branches")
)
gene_sum.to_csv("absrel_gene_summary.tsv", sep="\t", index=False)

# -----------------------------
# 3) per-species evolutionary history summary
# -----------------------------
species_records = []

for gene, sub in df[df["is_selected"]].groupby("gene"):
    treefile = f"{TREE_DIR}/{gene}.treefile"
    if not os.path.exists(treefile):
        print(f"[WARN] Tree file missing: {treefile}")
        continue

    try:
        t = Tree(treefile, format=1)
    except Exception as e:
        print(f"[WARN] Failed to parse tree: {treefile} ({e})")
        continue

    selected_branches = set(sub["branch"])

    # NOTE: This relies on branch names in the Newick matching HyPhy branch labels.
    for leaf in t.iter_leaves():
        species = leaf.name
        count = 0
        node = leaf
        while node.up:
            if node.name in selected_branches:
                count += 1
            node = node.up

        species_records.append({"species": species, "n_selection_events": count})

if len(species_records) == 0:
    print("[WARN] No species-level records generated (branch names may not match between HyPhy JSON and tree).")
else:
    sp_df = (
        pd.DataFrame(species_records)
        .groupby("species", as_index=False)
        .agg({"n_selection_events": "sum"})
    )
    sp_df.to_csv("absrel_species_summary.tsv", sep="\t", index=False)

print("Done.")
