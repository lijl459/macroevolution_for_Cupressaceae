#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable, List, Set


DEFAULT_DONE_PATTERNS = [
    r"Adaptive branch site random effects likelihood test",
    r"Likelihood ratio test for episodic diversifying positive selection",
    r"branch site random effects likelihood",
]

DEFAULT_ERROR_PATTERNS = [
    r"\bERROR\b",
    r"\bFATAL\b",
    r"Traceback \(most recent call last\)",
    r"Segmentation fault",
    r"\bKilled\b",
    r"Out of memory",
    r"cannot allocate memory",
    r"AssertionError",
    r"Exception:",
]


def compile_any(patterns: Iterable[str]) -> re.Pattern:
    joined = "|".join(f"(?:{p})" for p in patterns)
    return re.compile(joined, flags=re.IGNORECASE)


def iter_gene_ids_from_tree_dir(tree_dir: Path, suffix: str) -> Set[str]:
    genes: Set[str] = set()
    for fp in tree_dir.glob(f"*{suffix}"):
        if fp.is_file():
            genes.add(fp.stem)  # OG0011726_inclade1_ortho1
    return genes


def tail_lines(path: Path, n: int = 2000) -> List[str]:
    """
    Read last n lines efficiently from a potentially large log.
    """
    try:
        with path.open("rb") as f:
            f.seek(0, 2)
            size = f.tell()
            block = 8192
            data = b""
            pos = size
            while pos > 0 and data.count(b"\n") <= n:
                step = block if pos >= block else pos
                pos -= step
                f.seek(pos)
                data = f.read(step) + data
            lines = data.decode("utf-8", errors="ignore").splitlines()
            return lines[-n:] if len(lines) > n else lines
    except Exception:
        return []


def log_is_finished(
    log_path: Path,
    done_re: re.Pattern,
    err_re: re.Pattern,
    tail_n: int,
) -> bool:
    """
    Finished if:
      - tail contains a done marker
      - and tail does NOT contain an error marker
    Otherwise unfinished.
    """
    lines = tail_lines(log_path, n=tail_n)
    if not lines:
        return False

    text = "\n".join(lines)

    if err_re.search(text):
        return False

    if done_re.search(text):
        return True

    return False


def gene_finished_in_any_dir(
    gene: str,
    absrel_dirs: List[Path],
    done_re: re.Pattern,
    err_re: re.Pattern,
    tail_n: int,
) -> bool:
    """
    Consider gene finished if ANY absrel_dir contains a finished log.
    """
    for d in absrel_dirs:
        log_path = d / f"{gene}.log"
        if not log_path.exists():
            continue
        if log_is_finished(log_path, done_re, err_re, tail_n):
            return True
    return False


def main():
    ap = argparse.ArgumentParser(
        description="Output unfinished genes for HyPhy aBSREL with multiple result dirs: missing log OR log not finished in ALL dirs."
    )
    ap.add_argument("-g", "--gene-tree-dir", required=True, help="Directory containing gene tree files")
    ap.add_argument(
        "-a", "--absrel-dir",
        required=True,
        action="append",
        help="aBSREL output directory (repeatable): expects <gene>.log inside",
    )
    ap.add_argument(
        "--tree-suffix",
        default=".treefile",
        help="Gene tree filename suffix (default: .treefile)",
    )
    ap.add_argument(
        "--tail",
        type=int,
        default=2000,
        help="Read last N lines of log for checks (default: 2000)",
    )
    ap.add_argument(
        "--done-pattern",
        action="append",
        default=[],
        help="Regex for completion marker (can repeat). Defaults include aBSREL summary lines.",
    )
    ap.add_argument(
        "--error-pattern",
        action="append",
        default=[],
        help="Regex for error marker (can repeat). Defaults include ERROR/FATAL/Traceback/Killed etc.",
    )
    ap.add_argument("-o", "--out", default="-", help="Output file (default: stdout)")
    args = ap.parse_args()

    tree_dir = Path(args.gene_tree_dir)
    if not tree_dir.exists():
        raise FileNotFoundError(f"Gene tree dir not found: {tree_dir}")

    absrel_dirs = [Path(x) for x in args.absrel_dir]
    for d in absrel_dirs:
        if not d.exists():
            raise FileNotFoundError(f"aBSREL dir not found: {d}")

    genes = iter_gene_ids_from_tree_dir(tree_dir, args.tree_suffix)
    if not genes:
        sys.stderr.write(f"[WARN] No gene trees found in {tree_dir} with suffix {args.tree_suffix}\n")

    done_patterns = args.done_pattern if args.done_pattern else DEFAULT_DONE_PATTERNS
    err_patterns = args.error_pattern if args.error_pattern else DEFAULT_ERROR_PATTERNS
    done_re = compile_any(done_patterns)
    err_re = compile_any(err_patterns)

    unfinished: List[str] = []
    for gene in sorted(genes):
        if not gene_finished_in_any_dir(gene, absrel_dirs, done_re, err_re, args.tail):
            unfinished.append(gene)

    out_fp = sys.stdout if args.out == "-" else open(args.out, "w", encoding="utf-8")
    try:
        for g in unfinished:
            out_fp.write(g + "\n")
    finally:
        if out_fp is not sys.stdout:
            out_fp.close()

    sys.stderr.write(
        f"[INFO] expected_genes={len(genes)} absrel_dirs={len(absrel_dirs)} unfinished={len(unfinished)}\n"
    )


if __name__ == "__main__":
    main()
