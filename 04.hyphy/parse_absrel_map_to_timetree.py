#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
import json
import math
import argparse
from typing import Dict, List, Optional

import pandas as pd
from ete3 import Tree


# -----------------------------
# numeric helpers
# -----------------------------
def to_float(x):
    if x is None:
        return None
    if isinstance(x, (int, float)):
        try:
            if math.isnan(float(x)):
                return None
        except Exception:
            pass
        return float(x)
    if isinstance(x, str):
        try:
            return float(x)
        except Exception:
            return None
    if isinstance(x, dict):
        for k in ["value", "p", "P", "p-value", "P-value"]:
            if k in x:
                return to_float(x[k])
    return None


def safe_omega(dn: Optional[float], ds: Optional[float], ds_eps: float) -> Optional[float]:
    if dn is None or ds is None:
        return None
    if ds <= ds_eps:
        return None
    return dn / ds


def clean_braces(s: str) -> str:
    return re.sub(r"\{[^{}]*\}", "", s)


def get_corrected_p(info: dict) -> Optional[float]:
    for key in ["Corrected P-value", "Corrected p-value", "corrected p-value", "corrected_pvalue"]:
        if key in info:
            return to_float(info.get(key))
    for key in ["Uncorrected P-value", "P-value", "p-value", "pvalue", "p"]:
        if key in info:
            return to_float(info.get(key))
    for key in ["LRT", "LRT test", "test", "Test"]:
        if key in info and isinstance(info[key], dict):
            for k2 in ["Corrected P-value", "Uncorrected P-value", "P-value", "p-value", "p"]:
                if k2 in info[key]:
                    return to_float(info[key][k2])
    return None


# -----------------------------
# log parsing (your rule)
# -----------------------------
def extract_annotated_tree(log_path: str) -> Optional[str]:
    """
    Count lines containing 'Annotated Tree' (case-insensitive).
    Use 2nd occurrence if exists, else 1st.
    Tree starts at NEXT line, concatenate until ';'.
    """
    try:
        with open(log_path, "r", errors="ignore") as f:
            lines = f.readlines()
    except Exception:
        return None

    pat = re.compile(r"annotated\s+tree", flags=re.IGNORECASE)
    hits = [i for i, line in enumerate(lines) if pat.search(line)]
    if not hits:
        return None

    header_idx = hits[1] if len(hits) >= 2 else hits[0]
    start = header_idx + 1
    if start >= len(lines):
        return None

    buf = []
    for j in range(start, len(lines)):
        s = lines[j].strip()
        if s == "":
            if buf and ";" in "".join(buf):
                break
            continue
        buf.append(s)
        if ";" in s:
            break

    if not buf:
        return None

    nwk = "".join(buf)
    if ";" not in nwk:
        return None

    nwk = nwk.split(";", 1)[0] + ";"
    nwk = clean_braces(nwk)
    nwk = re.sub(r"\s+", "", nwk)
    return nwk


# -----------------------------
# trees
# -----------------------------
def load_time_tree(path: str) -> Tree:
    tr = Tree(path, format=1)
    tips = [lf.name for lf in tr.iter_leaves()]
    if any((t is None or t == "") for t in tips):
        raise RuntimeError("Time tree has unnamed tips. Ensure tip labels are species names.")
    return tr


def build_gene_node_to_time_node_map(time_tree: Tree, annotated_newick: str) -> Dict[str, str]:
    try:
        gt = Tree(annotated_newick, format=1)
    except Exception:
        gt = Tree(annotated_newick, format=0)

    time_tips = set([lf.name for lf in time_tree.iter_leaves()])
    mapping: Dict[str, str] = {}

    for node in gt.traverse("postorder"):
        if node.is_leaf() or not node.name:
            continue

        taxa = sorted({lf.name for lf in node.iter_leaves() if lf.name in time_tips})
        if len(taxa) == 0:
            continue
        if len(taxa) == 1:
            mapping[node.name] = taxa[0]
            continue

        mrca = time_tree.get_common_ancestor(taxa)
        mapping[node.name] = mrca.name if mrca.name else taxa[0]

    return mapping


# -----------------------------
# main
# -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--absrel_dir", required=True)
    ap.add_argument("--time_tree", required=True)
    ap.add_argument("--p_cutoff", type=float, default=0.05)
    ap.add_argument("--ds_eps", type=float, default=1e-6)
    ap.add_argument("--out_prefix", default="absrel")
    ap.add_argument("--debug_show_tree", action="store_true")
    ap.add_argument("--max_debug_genes", type=int, default=5)
    args = ap.parse_args()

    time_tree = load_time_tree(args.time_tree)

    json_files = sorted(glob.glob(os.path.join(args.absrel_dir, "*.json")))
    if not json_files:
        raise SystemExit(f"[ERROR] No *.json found under: {args.absrel_dir}")

    node2bl: Dict[str, float] = {n.name: float(n.dist) for n in time_tree.traverse("preorder") if n.name}

    records: List[dict] = []
    node_map_rows: List[dict] = []
    problems: List[str] = []

    dbg_printed = 0

    for jf in json_files:
        gene = os.path.basename(jf).replace(".json", "")
        logf = os.path.join(args.absrel_dir, gene + ".log")

        if not os.path.exists(logf):
            problems.append(f"NOLOG\t{gene}\t{jf}\t{logf}\t")
            continue

        ann_nwk = extract_annotated_tree(logf)
        if not ann_nwk:
            problems.append(f"NO_ANNOTATED_TREE\t{gene}\t{jf}\t{logf}\t")
            continue

        if args.debug_show_tree and dbg_printed < args.max_debug_genes:
            print(f"[DEBUG] {gene}: {ann_nwk[:200]}")
            dbg_printed += 1

        try:
            gene_node2time = build_gene_node_to_time_node_map(time_tree, ann_nwk)
        except Exception as e:
            problems.append(f"ANNOTATED_TREE_PARSE_FAIL\t{gene}\t{jf}\t{logf}\t{e}")
            continue

        for gn, tn in gene_node2time.items():
            node_map_rows.append({"gene": gene, "gene_node": gn, "timeTree_node": tn})

        if os.path.getsize(jf) == 0:
            problems.append(f"EMPTYJSON\t{gene}\t{jf}\t{logf}\t")
            continue

        try:
            with open(jf, "r") as f:
                data = json.load(f)
        except Exception as e:
            problems.append(f"BADJSON\t{gene}\t{jf}\t{logf}\t{e}")
            continue

        # HARD FIX: your JSON is standard => take exactly this path
        try:
            branches = data["branch attributes"]["0"]
        except Exception as e:
            problems.append(f"NO_BRANCHATTR_0\t{gene}\t{jf}\t{logf}\t{e}")
            continue

        if not isinstance(branches, dict) or not branches:
            problems.append(f"EMPTY_BRANCHATTR_0\t{gene}\t{jf}\t{logf}\t")
            continue

        # parse branches
        for branch_key, info in branches.items():
            if branch_key == "attributes":
                continue
            if not isinstance(info, dict):
                continue

            branch_raw = clean_braces(str(branch_key))
            species_label = clean_braces(str(info.get("original name", branch_raw)))

            time_node = gene_node2time.get(branch_raw)
            if time_node is None:
                time_node = gene_node2time.get(species_label)
            if time_node is None:
                time_node = species_label

            dn = to_float(info.get("Full adaptive model (non-synonymous subs/site)"))
            ds = to_float(info.get("Full adaptive model (synonymous subs/site)"))
            omega = safe_omega(dn, ds, ds_eps=args.ds_eps)

            omega_baseline = to_float(info.get("Baseline MG94xREV omega ratio"))
            p = get_corrected_p(info)
            is_sel = (p is not None) and (p < args.p_cutoff)

            records.append({
                "gene": gene,
                "branch_raw": branch_raw,
                "species_label": species_label,
                "timeTree_node": time_node,
                "branch_length": node2bl.get(time_node, float("nan")),
                "dN": dn,
                "dS": ds,
                "omega": omega,
                "omega_baseline": omega_baseline,
                "corrected_p": p,
                "is_selected": is_sel,
            })

    # always write problems
    with open("absrel_parse_problems.tsv", "w") as w:
        w.write("type\tgene\tjson\tlog\textra\n")
        for p in problems:
            w.write(p + "\n")

    if not records:
        raise SystemExit("[ERROR] No records parsed. See absrel_parse_problems.tsv")

    df = pd.DataFrame(records)
    for c in ["dN", "dS", "omega", "omega_baseline", "corrected_p", "branch_length"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["is_selected"] = df["is_selected"].astype(bool)

    out_branch = f"{args.out_prefix}_branch_mapped.tsv"
    df.to_csv(out_branch, sep="\t", index=False)

    node_map_df = pd.DataFrame(node_map_rows).drop_duplicates()
    node_map_df.to_csv("node_map_timeTree.tsv", sep="\t", index=False)

    gnode = (
        df.groupby(["timeTree_node", "gene"], as_index=False)
          .agg(is_selected=("is_selected", "max"),
               branch_length=("branch_length", "first"))
    )
    summary = (
        gnode.groupby("timeTree_node", as_index=False)
             .agg(tested_genes=("gene", "nunique"),
                  events=("is_selected", "sum"),
                  branch_length=("branch_length", "first"))
    )
    summary["event_rate"] = summary.apply(
        lambda r: (r["events"] / r["tested_genes"]) if (r["tested_genes"] and r["tested_genes"] > 0) else float("nan"),
        axis=1
    )
    summary["density"] = summary.apply(
        lambda r: (r["events"] / (r["tested_genes"] * r["branch_length"]))
        if (r["tested_genes"] and r["tested_genes"] > 0 and pd.notna(r["branch_length"]) and r["branch_length"] > 0)
        else float("nan"),
        axis=1
    )
    summary.sort_values(["events", "tested_genes"], ascending=[False, False]).to_csv(
        "timeTree_node_summary.tsv", sep="\t", index=False
    )

    print("Done.")
    print(f"Parsed genes: {df['gene'].nunique()} | records: {len(df)}")
    print(f"Wrote: {out_branch}")
    print("Wrote: node_map_timeTree.tsv")
    print("Wrote: timeTree_node_summary.tsv")
    print("Wrote: absrel_parse_problems.tsv")


if __name__ == "__main__":
    main()
