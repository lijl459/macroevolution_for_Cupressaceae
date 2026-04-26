#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, glob, json, math, hashlib
import pandas as pd
from ete3 import Tree

# ================== 需要你改的路径 ==================
SPECIES_TREE = "../Cup173.regular.tre"   # 物种树 newick（tip 名必须与基因树一致）
ABSREL_DIR   = "../03.absrel"          # aBSREL 输出目录（里面有 *.json 和 *.log）
P_CUTOFF     = 0.05                 # 显著性阈值
DS_EPS       = 10e-6                 # dS 太小就不算 omega（避免 1e-10 造成爆炸）
# ====================================================

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
    return None

def to_str(x):
    if x is None:
        return None
    if isinstance(x, str):
        return x
    if isinstance(x, (int, float, bool)):
        return str(x)
    if isinstance(x, (dict, list)):
        return json.dumps(x, sort_keys=True, ensure_ascii=False)
    return str(x)

def safe_omega(dn, ds):
    if dn is None or ds is None:
        return None
    if ds <= DS_EPS:
        return None
    return dn / ds

def leafset_id(leaves, prefix="ST"):
    """Stable ID for a species-tree node based on descendant leaf set."""
    s = "|".join(sorted(leaves))
    h = hashlib.md5(s.encode("utf-8")).hexdigest()[:12]
    return f"{prefix}_{h}"

def name_species_tree_internal_nodes(st: Tree):
    """Assign stable IDs (ST_xxx) to all internal nodes based on their descendant leaf set."""
    for node in st.traverse("postorder"):
        if node.is_leaf():
            continue
        leaves = [lf.name for lf in node.iter_leaves()]
        node.name = leafset_id(leaves, prefix="ST")

def extract_annotated_tree_newick(log_path):
    """
    Extract Newick string after '### Annotated Tree' from HyPhy log.
    Join wrapped lines until blank line or next '###'.
    """
    with open(log_path, "r", errors="ignore") as f:
        lines = f.readlines()

    start = None
    for i, line in enumerate(lines):
        if "### Annotated Tree" in line:
            start = i + 1
            break
    if start is None:
        return None

    buf = []
    for j in range(start, len(lines)):
        s = lines[j].strip()
        if s.startswith("###") and j > start:
            break
        if s == "" and buf:
            break
        if s != "":
            buf.append(s)

    if not buf:
        return None

    nwk = "".join(buf)
    if ";" in nwk:
        nwk = nwk.split(";", 1)[0] + ";"
    else:
        return None
    return nwk

def get_corrected_p(info: dict):
    """Prefer corrected p-value; fallback to uncorrected if needed."""
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

def get_branch_dict_from_absrel_json(data: dict):
    """
    你的输出常见结构：
      data["branch attributes"] = {"0": {branch->info}, "attributes": {...schema...}}
    取真正分支字典：["0"]；若没有则用顶层。
    """
    ba = data.get("branch attributes", {})
    if not isinstance(ba, dict) or len(ba) == 0:
        return None
    if "0" in ba and isinstance(ba["0"], dict):
        return ba["0"]
    return ba

def build_gene_node_to_unified_map(species_tree: Tree, annotated_tree_newick: str):
    """
    For each node in gene annotated tree:
      - if tip: unified = tip species name
      - if internal NodeX: unified = MRCA( descendant tips ) in species tree (node ST_xxx)
    Return dict: gene_node_name -> unified_name
    """
    # parse gene tree
    try:
        gt = Tree(annotated_tree_newick, format=1)
    except Exception:
        gt = Tree(annotated_tree_newick, format=0)

    # cache species-tree tip set
    st_tipset = set([lf.name for lf in species_tree.iter_leaves()])

    mapping = {}

    for node in gt.traverse("postorder"):
        if node.is_leaf():
            sp = node.name
            if sp:
                mapping[sp] = sp  # tip maps to itself
            continue

        # internal node name like Node3; sometimes may be empty
        gene_node = node.name
        if not gene_node:
            continue

        taxa = [lf.name for lf in node.iter_leaves() if lf.name in st_tipset]
        taxa = sorted(set(taxa))

        if len(taxa) == 0:
            continue
        if len(taxa) == 1:
            mapping[gene_node] = taxa[0]  # MRCA of one taxon = itself
            continue

        mrca = species_tree.get_common_ancestor(taxa)
        unified = mrca.name if (mrca.name and not mrca.is_leaf()) else mrca.name
        if mrca.is_leaf():
            unified = mrca.name
        mapping[gene_node] = unified

    return mapping

def main():
    # 1) load species tree and name internal nodes
    st = Tree(SPECIES_TREE, format=1)
    for lf in st.iter_leaves():
        if not lf.name:
            raise SystemExit("[ERROR] Species tree has unnamed tips. Please ensure all tip labels are species names.")
    name_species_tree_internal_nodes(st)
    print(f"[OK] Loaded species tree: {SPECIES_TREE}")

    # 2) collect absrel outputs
    json_files = sorted(glob.glob(os.path.join(ABSREL_DIR, "*.json")))
    if not json_files:
        raise SystemExit(f"[ERROR] No *.json found in {ABSREL_DIR}/")

    records = []
    problems = []

    for jf in json_files:
        gene = os.path.basename(jf).replace(".json", "")
        logf = os.path.join(ABSREL_DIR, gene + ".log")

        # read annotated tree (needed for NodeX->unified mapping)
        if not os.path.exists(logf):
            problems.append(f"NOLOG\t{jf}\t{logf}")
            continue

        ann_nwk = extract_annotated_tree_newick(logf)
        if not ann_nwk:
            problems.append(f"NO_ANNOTATED_TREE\t{jf}\t{logf}")
            continue

        try:
            node2uni = build_gene_node_to_unified_map(st, ann_nwk)
        except Exception as e:
            problems.append(f"ANNOTATED_TREE_PARSE_FAIL\t{jf}\t{logf}\t{e}")
            continue

        # read absrel json
        if os.path.getsize(jf) == 0:
            problems.append(f"EMPTYJSON\t{jf}")
            continue

        try:
            data = json.load(open(jf, "r"))
        except Exception as e:
            problems.append(f"BADJSON\t{jf}\t{e}")
            continue

        branches = get_branch_dict_from_absrel_json(data)
        if branches is None or not isinstance(branches, dict) or len(branches) == 0:
            problems.append(f"NOBRANCHATTR\t{jf}")
            continue

        # iterate branch keys (these should include tips and NodeX)
        for branch, info in branches.items():
            if branch == "attributes":
                continue
            if not isinstance(info, dict):
                continue

            branch_key = to_str(branch)

            # preferred species label from json, fallback branch key
            sp_label = to_str(info.get("original name", branch_key))

            # unified node name:
            # - if branch_key is NodeX, map via node2uni
            # - if branch_key is a tip name, maps to itself
            # - else fallback to sp_label
            branch_unified = node2uni.get(branch_key, node2uni.get(sp_label, sp_label))

            dn = to_float(info.get("Full adaptive model (non-synonymous subs/site)"))
            ds = to_float(info.get("Full adaptive model (synonymous subs/site)"))
            omega = safe_omega(dn, ds)

            omega_baseline = to_float(info.get("Baseline MG94xREV omega ratio"))
            p = get_corrected_p(info)
            is_sel = (p is not None) and (p < P_CUTOFF)

            records.append({
                "gene": gene,
                "branch_raw": branch_key,        # HyPhy branch key: tip or NodeX
                "species_label": sp_label,       # HyPhy original name (often tip name)
                "branch_unified": branch_unified,# MRCA projection on species tree (ST_xxx or tip)
                "dN": dn,
                "dS": ds,
                "omega": omega,
                "omega_baseline": omega_baseline,
                "corrected_p": p,
                "is_selected": is_sel
            })

    if not records:
        raise SystemExit("[ERROR] No records parsed. Check paths / log format / json format.")

    df = pd.DataFrame(records)

    # numeric columns
    for c in ["dN", "dS", "omega", "omega_baseline", "corrected_p"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # keys
    for c in ["gene", "branch_raw", "species_label", "branch_unified"]:
        df[c] = df[c].apply(to_str)

    # output 1: unified branch table (核心)
    df.to_csv("absrel_branch_unified.tsv", sep="\t", index=False)

    # output 2: per unified-branch summary across genes (便于跨基因比较同一物种树节点)
    ub = (
        df.groupby(["branch_unified"], as_index=False)
          .agg(
              n_genes=("gene", "nunique"),
              n_events=("is_selected", "sum"),
              min_p=("corrected_p", "min"),
              median_omega=("omega", "median"),
              mean_omega=("omega", "mean"),
              max_omega=("omega", "max"),
          )
          .sort_values(["n_events", "n_genes"], ascending=[False, False])
    )
    ub.to_csv("absrel_unified_branch_summary.tsv", sep="\t", index=False)

    # output 3: per species (tips) summary — 只统计 branch_unified == species_name 的 terminal 分支
    # 这样得到“每个物种 terminal branch 的事件/omega”等
    tips = set([lf.name for lf in st.iter_leaves()])
    df_tip = df[df["branch_unified"].isin(tips)].copy()

    sp = (
        df_tip.groupby("branch_unified", as_index=False)
              .agg(
                  n_genes_tested=("gene", "nunique"),
                  n_selection_events=("is_selected", "sum"),
                  median_omega=("omega", "median"),
                  mean_omega=("omega", "mean"),
                  max_omega=("omega", "max"),
                  mean_dN=("dN", "mean"),
                  mean_dS=("dS", "mean"),
              )
              .rename(columns={"branch_unified": "species"})
              .sort_values(["n_selection_events", "n_genes_tested"], ascending=[False, False])
    )
    sp.to_csv("absrel_species_summary.tsv", sep="\t", index=False)

    # output problems
    if problems:
        with open("absrel_parse_problems.tsv", "w") as w:
            w.write("type\tjson\tlog\textra\n")
            for p in problems:
                w.write(p + "\n")

    print("Done.")
    print(f"Parsed genes: {df['gene'].nunique()} | records: {len(df)}")
    print("Wrote:")
    print("  absrel_branch_unified.tsv")
    print("  absrel_unified_branch_summary.tsv")
    print("  absrel_species_summary.tsv")
    if problems:
        print("  absrel_parse_problems.tsv (missing log / missing annotated tree / bad json etc.)")

if __name__ == "__main__":
    main()
