#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Set, Tuple, List

from ete3 import Tree


NODE_RE = re.compile(r"^Node(\d+)$")
BRACE_RE = re.compile(r"\{[^}]*\}")


def safe_float(x: Any) -> float:
    try:
        if x is None:
            return float("nan")
        return float(x)
    except Exception:
        return float("nan")


def load_species_tree(tree_path: str) -> Tuple[Tree, Set[str], Set[str], int]:
    st = Tree(tree_path, format=1)
    tips = set(st.get_leaf_names())
    internal = set()
    max_node_id = 0
    for n in st.traverse():
        if n.is_leaf():
            continue
        if n.name:
            internal.add(n.name)
            m = NODE_RE.match(n.name)
            if m:
                max_node_id = max(max_node_id, int(m.group(1)))
    return st, tips, internal, max_node_id


def iter_json_files(in_dir: str):
    p = Path(in_dir)
    if not p.exists():
        raise FileNotFoundError(f"Input directory not found: {in_dir}")
    for fp in sorted(p.glob("*.json")):
        yield fp


def get_gene_name(fp: Path) -> str:
    return fp.stem


def extract_branch_attributes(data: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """
    Return {branch_name: stats_dict}
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


def choose_pvalue(stats: Dict[str, Any], use_raw_only: bool) -> float:
    if use_raw_only:
        return safe_float(stats.get("Uncorrected P-value"))

    p_corr = safe_float(stats.get("Corrected P-value"))
    if not math.isnan(p_corr):
        return p_corr
    return safe_float(stats.get("Uncorrected P-value"))


def parse_first_annotated_tree_next_line_from_log(log_path: Path) -> Optional[str]:
    """
    Your log format:
      ### Annotated Tree
      <ONE-LINE NEWICK WITHOUT TRAILING ';' POSSIBLY>
    We:
      - take the FIRST occurrence
      - take the next non-empty line as the tree line
      - remove {...} annotations
      - append ';' if missing
    """
    if not log_path.exists():
        return None

    try:
        lines = log_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return None

    for i, line in enumerate(lines):
        if line.strip().startswith("### Annotated Tree"):
            # next non-empty line
            j = i + 1
            while j < len(lines) and lines[j].strip() == "":
                j += 1
            if j >= len(lines):
                return None
            nwk = lines[j].strip()
            if not nwk:
                return None

            # strip brace annotations like {1} {23}
            nwk = BRACE_RE.sub("", nwk).strip()

            # ete3 is happier with semicolon
            if not nwk.endswith(";"):
                nwk += ";"
            return nwk

    return None


def load_gene_tree_from_log(log_path: Path) -> Optional[Tree]:
    nwk = parse_first_annotated_tree_next_line_from_log(log_path)
    if not nwk:
        return None

    for fmt in (1, 0, 9):
        try:
            return Tree(nwk, format=fmt)
        except Exception:
            pass

    return None


def map_branch_to_species_tree_node(
    branch: str,
    gene_tree: Tree,
    species_tree: Tree,
    species_tips: Set[str],
    max_species_node_id: int,
) -> Tuple[Optional[str], Optional[str]]:
    """
    Return (mapped, reason).
    NodeX must be mapped via gene_tree tips -> species_tree MRCA.
    """
    if branch in species_tips:
        return branch, None

    if not NODE_RE.match(branch):
        return None, "not_tip_and_not_NodeX"

    nodes = gene_tree.search_nodes(name=branch)
    if not nodes:
        return None, "NodeX_not_found_in_gene_tree"

    gt_node = nodes[0]
    leaves = gt_node.get_leaf_names()
    if not leaves:
        return None, "gene_tree_NodeX_has_no_leaves"

    leaves2 = [x for x in leaves if x in species_tips]
    if len(leaves2) == 0:
        return None, f"leaves_not_in_species_tree (n_leaves={len(leaves)})"
    if len(leaves2) == 1:
        return leaves2[0], None

    try:
        mrca = species_tree.get_common_ancestor(leaves2)
    except Exception as e:
        return None, f"species_tree_get_common_ancestor_failed: {type(e).__name__}"

    if not mrca.name:
        return None, "species_tree_MRCA_has_no_name"

    m2 = NODE_RE.match(mrca.name)
    if not m2:
        return None, "species_tree_MRCA_name_not_NodeX"

    if int(m2.group(1)) > max_species_node_id:
        return None, "species_tree_MRCA_NodeX_exceeds_max_id"

    return mrca.name, None


def main():
    ap = argparse.ArgumentParser(
        description="Count selected genes per species-tree branch for HyPhy aBSREL outputs, mapping gene-tree NodeX to species-tree NodeX by MRCA."
    )
    ap.add_argument("-i", "--input_dir", required=True, help="Input dir of aBSREL outputs, e.g. ../03.absrel")
    ap.add_argument("-t", "--tree", required=True, help="Species tree with internal nodes named, e.g. Cup173.NodeNamed.nwk")
    ap.add_argument("-o", "--output_prefix", required=True, help="Output prefix, e.g. absrel_branch_summary")
    ap.add_argument("--alpha", type=float, default=0.05, help="P-value cutoff (default 0.05)")
    ap.add_argument("--use-raw-only", action="store_true", help="Use Uncorrected P-value only (ignore Corrected P-value)")
    ap.add_argument("--write-gene-list", action="store_true", help="Write selected gene list per branch")

    ap.add_argument("--debug-map", action="store_true", help="Print mapping failures to stderr (limited by --debug-limit)")
    ap.add_argument("--debug-limit", type=int, default=200, help="Max number of debug lines to print (default 200)")
    args = ap.parse_args()

    species_tree, species_tips, species_internal, max_species_node_id = load_species_tree(args.tree)

    tested: Dict[str, Set[str]] = {}
    selected: Dict[str, Set[str]] = {}

    empty_json = 0
    bad_json = 0
    no_ba = 0

    warn_missing_log = 0
    warn_bad_gene_tree = 0
    warn_map_fail = 0

    debug_printed = 0

    gene_tree_cache: Dict[str, Optional[Tree]] = {}
    in_dir = Path(args.input_dir)

    def dprint(msg: str):
        nonlocal debug_printed
        if not args.debug_map:
            return
        if debug_printed >= args.debug_limit:
            return
        sys.stderr.write(msg.rstrip() + "\n")
        debug_printed += 1

    for fp in iter_json_files(args.input_dir):
        if fp.stat().st_size == 0:
            empty_json += 1
            continue

        try:
            with fp.open("r", encoding="utf-8") as f:
                data = json.load(f)
        except Exception:
            bad_json += 1
            continue

        gene = get_gene_name(fp)
        branch_map = extract_branch_attributes(data)
        if not branch_map:
            no_ba += 1
            continue

        # need gene tree if any NodeX exists (internal nodes must be mapped)
        need_gene_tree = any(NODE_RE.match(str(b)) is not None for b in branch_map.keys())

        gene_tree: Optional[Tree] = None
        if need_gene_tree:
            if gene in gene_tree_cache:
                gene_tree = gene_tree_cache[gene]
            else:
                log_path = in_dir / f"{gene}.log"
                if not log_path.exists():
                    warn_missing_log += 1
                    gene_tree_cache[gene] = None
                    gene_tree = None
                    dprint(f"[DEBUG] missing_log\tgene={gene}\tlog={log_path}")
                else:
                    gt = load_gene_tree_from_log(log_path)
                    if gt is None:
                        warn_bad_gene_tree += 1
                        dprint(f"[DEBUG] bad_gene_tree\tgene={gene}\tlog={log_path}")
                    gene_tree_cache[gene] = gt
                    gene_tree = gt

        for branch, stats in branch_map.items():
            if not isinstance(stats, dict):
                continue

            branch = str(branch)

            mapped: Optional[str] = None

            if branch in species_tips:
                mapped = branch
            else:
                if NODE_RE.match(branch):
                    if gene_tree is None:
                        warn_map_fail += 1
                        dprint(f"[DEBUG] map_fail_no_gene_tree\tgene={gene}\tbranch={branch}")
                        continue
                    mapped, reason = map_branch_to_species_tree_node(
                        branch=branch,
                        gene_tree=gene_tree,
                        species_tree=species_tree,
                        species_tips=species_tips,
                        max_species_node_id=max_species_node_id,
                    )
                    if mapped is None:
                        warn_map_fail += 1
                        dprint(f"[DEBUG] map_fail\tgene={gene}\tbranch={branch}\treason={reason}")
                        continue
                else:
                    continue

            if mapped is None:
                warn_map_fail += 1
                dprint(f"[DEBUG] map_fail_mapped_none\tgene={gene}\tbranch={branch}")
                continue

            tested.setdefault(mapped, set()).add(gene)

            p = choose_pvalue(stats, args.use_raw_only)
            if (not math.isnan(p)) and (p < args.alpha):
                selected.setdefault(mapped, set()).add(gene)

    all_branches = sorted(set(tested.keys()) | set(selected.keys()))
    rows = []
    for b in all_branches:
        tset = tested.get(b, set())
        sset = selected.get(b, set())
        n_tested = len(tset)
        n_sel = len(sset)
        prop = (n_sel / n_tested) if n_tested > 0 else 0.0
        row = {
            "branch": b,
            "n_tested_genes": n_tested,
            "n_selected_genes": n_sel,
            "prop_selected": prop,
        }
        if args.write_gene_list:
            row["selected_genes"] = ";".join(sorted(sset))
        rows.append(row)

    rows.sort(key=lambda r: (r["n_selected_genes"], r["prop_selected"], r["n_tested_genes"]), reverse=True)

    out_tsv = f"{args.output_prefix}.tsv"
    with open(out_tsv, "w", encoding="utf-8") as out:
        headers = list(rows[0].keys()) if rows else ["branch", "n_tested_genes", "n_selected_genes", "prop_selected"]
        out.write("\t".join(headers) + "\n")
        for r in rows:
            out.write("\t".join(str(r[h]) for h in headers) + "\n")

    max_out_node = 0
    for b in all_branches:
        m = NODE_RE.match(b)
        if m:
            max_out_node = max(max_out_node, int(m.group(1)))

    sys.stderr.write(
        "[INFO] Done.\n"
        f"[INFO] species_tree max Node = {max_species_node_id}\n"
        f"[INFO] output max Node = {max_out_node}\n"
        f"[INFO] branches_out = {len(all_branches)}\n"
        f"[INFO] empty_json = {empty_json}, bad_json = {bad_json}, no_branch_attributes = {no_ba}\n"
        f"[WARN] missing_log = {warn_missing_log}, bad_gene_tree = {warn_bad_gene_tree}, map_fail = {warn_map_fail}\n"
        f"[INFO] debug_printed = {debug_printed} (limit={args.debug_limit}, enabled={args.debug_map})\n"
        f"[INFO] wrote: {out_tsv}\n"
    )


if __name__ == "__main__":
    main()
