#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd


def safe_float(x) -> float:
    try:
        if x is None:
            return float("nan")
        return float(x)
    except Exception:
        return float("nan")


def is_finite(x: float) -> bool:
    return isinstance(x, (int, float)) and (not math.isnan(x)) and math.isfinite(x)


def parse_rate_distributions(branch_obj: Dict[str, Any]) -> Tuple[float, float]:
    """
    Return:
      omega_max: max omega among rate classes
      w_pos: sum of weights where omega > 1
    """
    rd = branch_obj.get("Rate Distributions", None)
    if not isinstance(rd, list) or len(rd) == 0:
        return float("nan"), float("nan")

    omegas = []
    w_pos = 0.0
    for item in rd:
        if not (isinstance(item, (list, tuple)) and len(item) >= 2):
            continue
        omega = safe_float(item[0])
        w = safe_float(item[1])
        if is_finite(omega):
            omegas.append(omega)
        if is_finite(omega) and is_finite(w) and omega > 1.0:
            w_pos += w

    omega_max = max(omegas) if omegas else float("nan")
    return omega_max, (w_pos if omegas else float("nan"))


def get_branch_map(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    For your JSON structure:
      data["branch attributes"]["0"] is the branch->stats map.
    Fallback:
      try any key inside data["branch attributes"] that is a dict and looks like branch map.
    """
    ba = data.get("branch attributes", None)
    if not isinstance(ba, dict):
        return {}

    # Prefer "0"
    if "0" in ba and isinstance(ba["0"], dict):
        return ba["0"]

    # Sometimes numeric 0
    if 0 in ba and isinstance(ba[0], dict):
        return ba[0]

    # Fallback: pick the first dict value
    for k, v in ba.items():
        if isinstance(v, dict) and len(v) > 0:
            return v

    return {}


def derive_gene_name(json_path: Path) -> str:
    return json_path.stem


def records_from_json(gene: str, data: Dict[str, Any], p_cutoff: float, use: str) -> List[Dict[str, Any]]:
    branch_map = get_branch_map(data)
    if not branch_map:
        return []

    records: List[Dict[str, Any]] = []
    for branch, bobj in branch_map.items():
        if not isinstance(bobj, dict):
            continue

        lrt = safe_float(bobj.get("LRT", float("nan")))
        p_raw = safe_float(bobj.get("Uncorrected P-value", float("nan")))
        p_fdr = safe_float(bobj.get("Corrected P-value", float("nan")))

        dN = safe_float(bobj.get("Full adaptive model (non-synonymous subs/site)", float("nan")))
        dS = safe_float(bobj.get("Full adaptive model (synonymous subs/site)", float("nan")))

        omega_max, w_pos = parse_rate_distributions(bobj)

        # Determine is_selected
        if use == "raw":
            is_sel = is_finite(p_raw) and (p_raw < p_cutoff)
        elif use == "fdr":
            is_sel = is_finite(p_fdr) and (p_fdr < p_cutoff)
        else:  # auto
            if is_finite(p_fdr):
                is_sel = p_fdr < p_cutoff
            else:
                is_sel = is_finite(p_raw) and (p_raw < p_cutoff)

        records.append({
            "gene": gene,
            "branch": str(branch),
            "LRT": lrt,
            "p_raw": p_raw,
            "p_FDR": p_fdr,
            "is_selected": bool(is_sel),
            "omega_max": omega_max,
            "w_pos": w_pos,
            "dN": dN,
            "dS": dS,
        })

    return records


def main():
    ap = argparse.ArgumentParser(description="Summarize HyPhy aBSREL JSON into one branch-level table.")
    ap.add_argument("-i", "--input", required=True, help="Input directory OR one JSON file.")
    ap.add_argument("-o", "--output", required=True, help="Output .tsv or .csv")
    ap.add_argument("--recursive", action="store_true", help="Recursively scan input dir for *.json")
    ap.add_argument("--p-cutoff", type=float, default=0.05, help="Cutoff for is_selected (default 0.05)")
    ap.add_argument("--use", choices=["auto", "fdr", "raw"], default="auto",
                    help="Which p-value to use for is_selected (default auto: prefer FDR).")
    args = ap.parse_args()

    input_path = Path(args.input)

    if input_path.is_file():
        json_files = [input_path]
    else:
        pattern = "**/*.json" if args.recursive else "*.json"
        json_files = sorted(Path(args.input).glob(pattern))

    if not json_files:
        raise SystemExit("ERROR: no json files found.")

    all_records: List[Dict[str, Any]] = []
    skipped = 0
    no_branch_attr = 0

    for jp in json_files:
        try:
            with jp.open("r", encoding="utf-8") as f:
                data = json.load(f)
        except Exception:
            skipped += 1
            continue

        if not isinstance(data, dict) or not data:
            skipped += 1
            continue

        gene = derive_gene_name(jp)
        recs = records_from_json(gene, data, args.p_cutoff, args.use)
        if not recs:
            no_branch_attr += 1
            continue
        all_records.extend(recs)

    df = pd.DataFrame(all_records, columns=[
        "gene", "branch", "LRT", "p_raw", "p_FDR", "is_selected", "omega_max", "w_pos", "dN", "dS"
    ])

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.suffix.lower() == ".csv":
        df.to_csv(out, index=False)
    else:
        df.to_csv(out, sep="\t", index=False)

    print(f"Done. rows={len(df)}  out={out}")
    if skipped:
        print(f"  skipped unreadable/empty json: {skipped}")
    if no_branch_attr:
        print(f"  json missing branch attributes (or empty): {no_branch_attr}")


if __name__ == "__main__":
    main()
