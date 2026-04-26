#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import sys
import shutil
from pathlib import Path
from typing import Iterable, List, Set

# 默认标志
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

def load_genes_from_list(list_file: Path) -> Set[str]:
    genes: Set[str] = set()
    with list_file.open("r", encoding="utf-8") as f:
        for line in f:
            gene = line.strip()
            if gene and not gene.startswith("#"):
                genes.add(gene)
    return genes

def tail_lines(path: Path, n: int = 2000) -> List[str]:
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

def log_is_finished(log_path: Path, done_re: re.Pattern, err_re: re.Pattern, tail_n: int) -> bool:
    lines = tail_lines(log_path, n=tail_n)
    if not lines:
        return False
    text = "\n".join(lines)
    if err_re.search(text): return False
    if done_re.search(text): return True
    return False

# --- 核心逻辑：使用 shutil.copy 进行复制 ---
def process_gene_files(
    gene: str, 
    absrel_dirs: List[Path], 
    dest_dir: Path, 
    done_re: re.Pattern, 
    err_re: re.Pattern, 
    tail_n: int
) -> bool:
    """
    检查基因是否完成。如果完成且指定了 dest_dir，则复制文件。
    """
    for d in absrel_dirs:
        log_path = d / f"{gene}.log"
        json_path = d / f"{gene}.json"
        
        if not log_path.exists():
            continue
            
        if log_is_finished(log_path, done_re, err_re, tail_n):
            if dest_dir:
                dest_dir.mkdir(parents=True, exist_ok=True)
                # 复制 log 文件
                if log_path.exists():
                    shutil.copy(str(log_path), str(dest_dir / log_path.name))
                # 复制 json 文件
                if json_path.exists():
                    shutil.copy(str(json_path), str(dest_dir / json_path.name))
            return True
    return False

def main():
    ap = argparse.ArgumentParser(description="Check aBSREL logs and copy finished files to a target directory.")
    ap.add_argument("-i", "--input-list", required=True, help="Text file with gene names (one per line)")
    ap.add_argument("-a", "--absrel-dir", required=True, action="append", help="aBSREL output directories (can repeat)")
    ap.add_argument("-d", "--dest-dir", help="Target directory to COPY finished .log and .json files")
    ap.add_argument("-o", "--out", default="-", help="Output for UNFINISHED genes (default: stdout)")
    ap.add_argument("--tail", type=int, default=2000, help="Lines to check at the end of log")
    
    args = ap.parse_args()

    input_list = Path(args.input_list)
    dest_path = Path(args.dest_dir) if args.dest_dir else None
    absrel_dirs = [Path(x) for x in args.absrel_dir]

    genes = load_genes_from_list(input_list)
    done_re = compile_any(DEFAULT_DONE_PATTERNS)
    err_re = compile_any(DEFAULT_ERROR_PATTERNS)

    unfinished: List[str] = []
    finished_count = 0

    for gene in sorted(genes):
        if process_gene_files(gene, absrel_dirs, dest_path, done_re, err_re, args.tail):
            finished_count += 1
        else:
            unfinished.append(gene)

    # 输出未完成的基因列表
    out_fp = sys.stdout if args.out == "-" else open(args.out, "w", encoding="utf-8")
    try:
        for g in unfinished:
            out_fp.write(g + "\n")
    finally:
        if out_fp is not sys.stdout:
            out_fp.close()

    sys.stderr.write(
        f"[INFO] Total_Input={len(genes)}, Finished_Found={finished_count}, Unfinished={len(unfinished)}\n"
    )
    if dest_path:
        sys.stderr.write(f"[INFO] Finished files copied to: {dest_path}\n")

if __name__ == "__main__":
    main()
