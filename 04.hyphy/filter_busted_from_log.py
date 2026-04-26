#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filter positively selected genes from HyPhy BUSTED *.log files.

Key idea:
- Don't trust printed "p = 0.0000" (rounded).
- Recompute LRT and p-value from LogL(full) and LogL(null) found in the log.
  BUSTED (v2+): p uses a 0.5*(chi^2_0 + chi^2_2) mixture.
  For LRT > 0:  p = 0.5 * P(ChiSq(df=2) >= LRT) = 0.5 * exp(-LRT/2)

Outputs:
- busted_log_summary.tsv
- busted_positive_genes.txt
"""

import argparse
import glob
import math
import os
import re
from typing import Optional, Tuple, Dict, Any, List

RE_FULL_HEADER = re.compile(r"Performing the full .*branch-site model fit", re.IGNORECASE)
RE_NULL_HEADER = re.compile(r"Performing the constrained .*model fit", re.IGNORECASE)
RE_LOGL = re.compile(r"\*\s*Log\(L\)\s*=\s*([\-0-9.]+)")
RE_P_REPORTED = re.compile(r"Likelihood ratio test.*?\bp\s*=\s*([0-9.eE+\-]+)", re.IGNORECASE)

def gene_id_from_filename(path: str) -> str:
    base = os.path.basename(path)
    base = re.sub(r"\.log$", "", base, flags=re.IGNORECASE)
    return base

def busted_p_from_lrt(lrt: float) -> float:
    # Mixture: 0.5*chi2_0 + 0.5*chi2_2
    # For lrt > 0, p = 0.5 * exp(-lrt/2) (df=2 survival)
    if lrt <= 0:
        return 1.0
    return 0.5 * math.exp(-0.5 * lrt)

def parse_busted_log(log_path: str) -> Dict[str, Any]:
    """
    Return dict with:
      logl_full, logl_null, lrt, p_calc, p_reported
    """
    logl_full = None
    logl_null = None
    p_reported = None

    in_full_block = False
    in_null_block = False

    with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            # reported p
            m_p = RE_P_REPORTED.search(line)
            if m_p:
                try:
                    p_reported = float(m_p.group(1))
                except Exception:
                    pass

            # detect blocks
            if RE_FULL_HEADER.search(line):
                in_full_block = True
                in_null_block = False
                continue
            if RE_NULL_HEADER.search(line):
                in_null_block = True
                in_full_block = False
                continue

            # capture the first Log(L) after entering each block
            m_ll = RE_LOGL.search(line)
            if m_ll:
                try:
                    ll = float(m_ll.group(1))
                except Exception:
                    continue
                if in_full_block and logl_full is None:
                    logl_full = ll
                    in_full_block = False  # stop after first
                elif in_null_block and logl_null is None:
                    logl_null = ll
                    in_null_block = False

    lrt = None
    p_calc = None
    if logl_full is not None and logl_null is not None:
        lrt = 2.0 * (logl_full - logl_null)
        # numerical guard
        if lrt < 0 and abs(lrt) < 1e-8:
            lrt = 0.0
        p_calc = busted_p_from_lrt(lrt)

    return {
        "logl_full": logl_full,
        "logl_null": logl_null,
        "lrt": lrt,
        "p_calc": p_calc,
        "p_reported": p_reported
    }

def fmt(x: Optional[float]) -> str:
    if x is None:
        return ""
    # keep scientific notation for tiny p-values
    return f"{x:.12g}"

def main():
    ap = argparse.ArgumentParser(description="Select positively selected genes from HyPhy BUSTED .log outputs.")
    ap.add_argument("-i", "--input", required=True, help="Folder containing *.log (e.g., 02.busted)")
    ap.add_argument("-p", "--pvalue", type=float, default=0.05, help="p-value threshold (default 0.05)")
    ap.add_argument("-o", "--outdir", default=".", help="Output directory (default current)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    logs = sorted(glob.glob(os.path.join(args.input, "*.log")))
    if not logs:
        raise SystemExit(f"[ERROR] No .log files found in: {args.input}")

    summary_path = os.path.join(args.outdir, "busted_log_summary.tsv")
    pos_path = os.path.join(args.outdir, "busted_positive_genes.txt")

    records: List[Dict[str, Any]] = []
    pos_genes: List[str] = []

    for lp in logs:
        gene = gene_id_from_filename(lp)
        r = parse_busted_log(lp)
        r["gene"] = gene
        r["file"] = os.path.basename(lp)

        # decide which p to use for filtering: prefer p_calc
        p_use = r["p_calc"] if r["p_calc"] is not None else r["p_reported"]
        r["p_used"] = p_use

        records.append(r)
        if p_use is not None and p_use < args.pvalue:
            pos_genes.append(gene)

    # write summary
    with open(summary_path, "w", encoding="utf-8") as out:
        out.write("\t".join([
            "gene", "file",
            "logl_full", "logl_null", "lrt",
            "p_calc", "p_reported", "p_used"
        ]) + "\n")
        for r in sorted(records, key=lambda x: x["gene"]):
            out.write("\t".join([
                r["gene"], r["file"],
                fmt(r["logl_full"]), fmt(r["logl_null"]), fmt(r["lrt"]),
                fmt(r["p_calc"]), fmt(r["p_reported"]), fmt(r["p_used"])
            ]) + "\n")

    # write positives
    with open(pos_path, "w", encoding="utf-8") as out:
        for g in sorted(set(pos_genes)):
            out.write(g + "\n")

    print(f"[OK] Wrote summary: {summary_path}")
    print(f"[OK] Positive genes (p < {args.pvalue}): {pos_path}  (n={len(set(pos_genes))})")

if __name__ == "__main__":
    main()
