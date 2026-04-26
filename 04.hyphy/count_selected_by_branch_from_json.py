#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Set, Tuple

import pandas as pd
from ete3 import Tree


# -------------------------
# generic helpers (keep your stable behavior)
# -------------------------
def safe_float(x: Any) -> float:
    try:
        if x is None:
            return float("nan")
        return float(x)
    except Exception:
        return float("nan")


def iter_json_files(in_dir: str) -> Iterable[Path]:
    p = Path(in_dir)
    if not p.exists():
        raise FileNotFoundError(f"Input directory not found: {in_dir}")
    for fp in sorted(p.glob("*.json")):
        yield fp


def get_gene_name(fp: Path, data: Dict[str, Any]) -> str:
    return fp.stem


def extract_branch_attributes(data: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """
    Return {branch_name: branch_stats_dict}
    Typical aBSREL:
      data["branch attributes"]["0"][branch_name] = {...}
    """
    ba = data.get("branch attributes")
    if not isinstance(ba, dict):
        return {}
    if "0" in ba and isinstance(ba["0"], dict):
        return ba["0"]
    if all(isinstance(v, dict) for v in ba.values()):
        return ba  # type: ignore
    return {}


def choose_pvalue(branch_stats: Dict[str, Any]) -> Tuple[float, str]:
    p_corr = safe_float(branch_stats.get("Corrected P-value"))
    if not math.isnan(p_corr):
        return p_corr, "corrected"
    p_raw = safe_float(branch_stats.get("Uncorrected P-value"))
    return p_raw, "uncorrected"


def is_selected_branch(branch_stats: Dict[str, Any], p_cutoff: float) -> bool:
    p, _ = choose_pvalue(branch_stats)
    if math.isnan(p):
        return False
    return p < p_cutoff


# -------------------------
# log parsing + node mapping
# -------------------------
ANNOTATED_TREE_TAG = "### Annotated Tree"

def find_second_annotated_tree_newick(log_path: Path) -> Optional[str]:
    """
    We need the SECOND occurrence of '### Annotated Tree' in the log,
    and the tree line(s) following it (the version WITHOUT {1} tags).

    Many logs have:
      ### Annotated Tree
      ((A{1},B{2})Node3{1},...)...
      ### Annotated Tree
      ((A,B)Node3,(...))...

    We'll extract the second tree robustly:
    - after the 2nd tag, collect subsequent lines that look like part of Newick
    - stop when parentheses are balanced and we hit end or a semicolon.
    """
    try:
        text = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except Exception:
        return None

    idxs = [i for i, line in enumerate(text) if line.strip().startswith(ANNOTATED_TREE_TAG)]
    if len(idxs) < 2:
        return None

    start = idxs[1] + 1
    # skip empty lines
    while start < len(text) and text[start].strip() == "":
        start += 1
    if start >= len(text):
        return None

    buf = []
    bal = 0
    started = False

    # a conservative “tree-line” detector:
    # lines that contain '(' or ')' or ',' and no obvious log prefix.
    for i in range(start, len(text)):
        line = text[i].strip()
        if not line:
            # allow blank lines inside wrapped newick; just skip
            continue

        # if we haven't started, require '(' to start
        if not started:
            if "(" not in line:
                continue
            started = True

        buf.append(line)

        # update bracket balance (ignore those inside quotes; usually none here)
        bal += line.count("(") - line.count(")")

        joined = "".join(buf)
        if ";" in joined:
            # cut at first ';'
            joined = joined.split(";", 1)[0] + ";"
            return joined

        # if balanced and seems like a complete newick, accept
        if started and bal == 0 and joined.endswith(")"):
            return joined + ";"

        # sometimes ends with tip (no trailing ')') but balanced==0
        if started and bal == 0 and (joined.endswith(",") is False):
            # still might be complete; add ';'
            return joined + ";"

    if not buf:
        return None
    joined = "".join(buf)
    if not joined.endswith(";"):
        joined += ";"
    return joined


def parse_tree_ete3(newick: str) -> Tree:
    """
    Try a few ete3 formats because node names like Node123 are internal,
    and sometimes there are odd characters.
    """
    # Try common formats first
    last_err = None
    for fmt in (1, 0, 2, 3, 5, 9):
        try:
            return Tree(newick, format=fmt, quoted_node_names=True)
        except Exception as e:
            last_err = e
    # final attempt without quoted_node_names
    for fmt in (1, 0, 2, 3, 5, 9):
        try:
            return Tree(newick, format=fmt)
        except Exception as e:
            last_err = e
    raise ValueError(str(last_err))


NODE_RE = re.compile(r"^Node\d+$")


def build_geneNode_to_speciesNode_map(
    gene_tree: Tree,
    species_tree: Tree,
) -> Dict[str, str]:
    """
    For each internal node in gene_tree named NodeX, map it to the MRCA node
    (named NodeY) in species_tree based on the leaf set under that gene-tree node.

    Returns dict: {gene_tree_node_name(NodeX): species_tree_node_name(NodeY)}
    """
    # Make sure species tree internal nodes have names
    # (you said Cup173.NodeNamed.nwk already does)
    sp_leaves = set(species_tree.get_leaf_names())

    mapping: Dict[str, str] = {}

    # iterate internal nodes with names Node###
    for n in gene_tree.traverse():
        if n.is_leaf():
            continue
        if not n.name or not NODE_RE.match(n.name):
            continue

        leaves = [x.name for x in n.get_leaves()]
        # keep only leaves existing in species tree
        leaves = [x for x in leaves if x in sp_leaves]

        if len(leaves) < 2:
            # MRCA undefined / not informative
            continue

        try:
            mrca = species_tree.get_common_ancestor(leaves)
        except Exception:
            continue

        if not mrca.name:
            # if unnamed, fallback to its node id string
            # but you requested unified NodeX, so we warn by leaving unmapped
            continue

        mapping[n.name] = mrca.name

    return mapping


def load_species_tree(path: str) -> Tree:
    nwk = Path(path).read_text(encoding="utf-8", errors="replace").strip()
    if not nwk.endswith(";"):
        nwk += ";"
    t = parse_tree_ete3(nwk)
    # ensure unique names; ete3 uses .name
    return t


# -------------------------
# main
# -------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Count number of selected genes per branch/species using HyPhy aBSREL JSON, "
                    "unifying internal NodeXX via species tree and per-gene annotated gene tree in .log."
    )
    ap.add_argument("-i", "--input_dir", required=True, help="Input directory, e.g. ../03.absrel")
    ap.add_argument("-t", "--species_tree", required=True, help="Species tree with Node-named internals, e.g. Cup173.NodeNamed.nwk")
    ap.add_argument("-o", "--output_prefix", required=True, help="Output prefix (dir/name prefix)")
    ap.add_argument("--alpha", type=float, default=0.05, help="P-value cutoff (default 0.05)")
    ap.add_argument("--use-raw-only", action="store_true", help="Use Uncorrected P-value only (ignore Corrected P-value)")
    args = ap.parse_args()

    in_dir = Path(args.input_dir)
    species_tree = load_species_tree(args.species_tree)

    tested_genes_by_branch: Dict[str, Set[str]] = {}
    selected_genes_by_branch: Dict[str, Set[str]] = {}

    bad_json = 0
    empty_json = 0
    log_missing = 0
    log_parse_fail = 0

    # cache: gene -> mapping dict
    mapping_cache: Dict[str, Dict[str, str]] = {}

    for json_fp in iter_json_files(str(in_dir)):
        # skip empty json
        if json_fp.stat().st_size == 0:
            empty_json += 1
            continue

        try:
            data = json.loads(json_fp.read_text(encoding="utf-8", errors="replace"))
        except Exception:
            bad_json += 1
            continue

        gene = get_gene_name(json_fp, data)

        # load mapping for this gene (from log)
        if gene not in mapping_cache:
            log_fp = in_dir / f"{gene}.log"
            if not log_fp.exists():
                mapping_cache[gene] = {}
                log_missing += 1
            else:
                newick = find_second_annotated_tree_newick(log_fp)
                if not newick:
                    mapping_cache[gene] = {}
                    log_parse_fail += 1
                else:
                    try:
                        gene_tree = parse_tree_ete3(newick)
                        mapping_cache[gene] = build_geneNode_to_speciesNode_map(gene_tree, species_tree)
                    except Exception:
                        mapping_cache[gene] = {}
                        log_parse_fail += 1

        node_map = mapping_cache[gene]

        branch_map = extract_branch_attributes(data)
        if not branch_map:
            continue

        for branch, stats in branch_map.items():
            if not isinstance(stats, dict):
                continue

            # unify branch name:
            # - tips keep as is
            # - Node### in gene tree -> mapped Node### in species tree (if present)
            unified = branch
            if isinstance(branch, str) and NODE_RE.match(branch):
                unified = node_map.get(branch, branch)

            tested_genes_by_branch.setdefault(unified, set()).add(gene)

            # decide selection
            if args.use_raw_only:
                p = safe_float(stats.get("Uncorrected P-value"))
                selected = (not math.isnan(p)) and (p < args.alpha)
            else:
                selected = is_selected_branch(stats, args.alpha)

            if selected:
                selected_genes_by_branch.setdefault(unified, set()).add(gene)

    # build output table
    rows = []
    all_branches = sorted(set(tested_genes_by_branch) | set(selected_genes_by_branch))
    for branch in all_branches:
        tested = tested_genes_by_branch.get(branch, set())
        selected = selected_genes_by_branch.get(branch, set())
        n_tested = len(tested)
        n_sel = len(selected)
        prop = (n_sel / n_tested) if n_tested > 0 else 0.0
        rows.append({
            "branch": branch,
            "n_tested_genes": n_tested,
            "n_selected_genes": n_sel,
            "prop_selected": prop,
        })

    df = pd.DataFrame(rows).sort_values(
        by=["n_selected_genes", "prop_selected", "n_tested_genes"],
        ascending=[False, False, False],
        kind="mergesort",
    )

    outdir = Path(args.output_prefix)
    # if output_prefix is a directory-like prefix, create its parent
    if outdir.suffix:
        # user gave file name; keep parent
        outdir.parent.mkdir(parents=True, exist_ok=True)
        out_path = outdir
    else:
        outdir.mkdir(parents=True, exist_ok=True)
        out_path = outdir / "branch_summary.tsv"

    df.to_csv(out_path, sep="\t", index=False)

    # optional: report logs
    report_path = out_path.parent / "run_report.txt"
    report = [
        f"input_dir={in_dir}",
        f"species_tree={args.species_tree}",
        f"alpha={args.alpha}",
        f"branches={len(df)}",
        f"empty_json={empty_json}",
        f"bad_json={bad_json}",
        f"log_missing={log_missing}",
        f"log_parse_fail={log_parse_fail}",
    ]
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")

    print(f"Done. Output: {out_path}")
    print("\n".join(report))


if __name__ == "__main__":
    main()
