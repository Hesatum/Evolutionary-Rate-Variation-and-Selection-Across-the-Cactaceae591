#!/usr/bin/env python3
"""
build_reduced_supermatrix.py

Builds a reduced supermatrix + IQ-TREE partition file by removing one or
more orthogroups from the full 591-locus (18-taxon) Cactaceae591
supermatrix, so the species tree can be re-estimated without them
(manuscript Methods 2.3: "the species tree was re-estimated from a reduced
supermatrix under the same settings used for the full dataset").

Used to build all three sensitivity trees discussed in the manuscript and
its supplementary material:
  1. exons under positive selection removed (Table 2, 21 orthogroups)
  2. introns with sustained clade-specific acceleration removed (Table S2)
  3. both sets removed together (the tree reported in the manuscript,
     Robinson-Foulds distance = 0 from the full tree, Figure S2)

The full partitions.txt is expected in RAxML/IQ-TREE partition format:
    DNA, <locus_name> = <start>-<end>
one line per partition, coordinates 1-based inclusive, matching the
columns of the accompanying full-supermatrix FASTA (49 samples in the
source dataset; the 18 final taxa listed in Table 1 are extracted by
their sample codes).

Usage
-----
    python build_reduced_supermatrix.py \
        --partitions partitions_full.txt \
        --fasta supermatrix_full.fasta \
        --exclude-list exclude_lists/exons_under_selection_21.csv \
        --taxa-codes taxa_codes.txt \
        --out-prefix trees/exon_selection_removed/supermatrix_18tax

`--exclude-list` is a CSV with one orthogroup ID per line (a `reason`
column is optional and only used in the log, e.g. "exon_selection" or
"intron_acceleration"). A partition is excluded when its name starts with
one of the listed orthogroup IDs, since a single orthogroup may contribute
more than one partition (e.g. `OG0056812.onlyexons...` and
`OG0056812_allexon...`).

`--taxa-codes` is a text file with one sample code per line, in the order
they should appear in the output FASTA (see manuscript Table 1, e.g.
S102A10, S80F9, ...). Each code must match the end of exactly one FASTA
header in the source alignment.
"""

import argparse
import csv
import re
from pathlib import Path


def parse_partitions(path: Path) -> list[tuple[str, int, int]]:
    partitions = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            m = re.match(r"DNA,\s*(\S+)\s*=\s*(\d+)-(\d+)", line)
            if not m:
                print(f"WARN: could not parse line: {line}")
                continue
            name, start, end = m.group(1), int(m.group(2)), int(m.group(3))
            partitions.append((name, start, end))
    return partitions


def read_fasta(path: Path) -> dict[str, str]:
    seqs, cur_id, chunks = {}, None, []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(chunks)
                cur_id, chunks = line[1:].strip(), []
            else:
                chunks.append(line)
        if cur_id is not None:
            seqs[cur_id] = "".join(chunks)
    return seqs


def read_exclude_list(path: Path) -> dict[str, str]:
    """Returns {orthogroup_id: reason}. Accepts a bare list or a CSV with
    'og'/'reason' columns; a missing reason column defaults to the file stem."""
    ogs = {}
    with open(path, encoding="utf-8") as f:
        reader = csv.reader(f)
        rows = list(reader)
    header = [c.strip().lower() for c in rows[0]] if rows else []
    has_header = "og" in header
    data_rows = rows[1:] if has_header else rows
    og_col = header.index("og") if has_header else 0
    reason_col = header.index("reason") if has_header and "reason" in header else None
    for row in data_rows:
        if not row or not row[og_col].strip():
            continue
        og = row[og_col].strip()
        reason = row[reason_col].strip() if reason_col is not None and len(row) > reason_col else path.stem
        ogs[og] = reason
    return ogs


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--partitions", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--exclude-list", required=True, type=Path,
                         help="CSV/text file with orthogroup IDs to remove (see module docstring).")
    parser.add_argument("--taxa-codes", required=True, type=Path,
                         help="Text file, one sample code per line, in output order.")
    parser.add_argument("--out-prefix", required=True, type=Path)
    args = parser.parse_args()

    exclude = read_exclude_list(args.exclude_list)
    target_codes = [c.strip() for c in args.taxa_codes.read_text(encoding="utf-8").splitlines() if c.strip()]

    log = []

    def logp(msg):
        print(msg)
        log.append(msg)

    partitions = parse_partitions(args.partitions)
    logp(f"Total partitions parsed: {len(partitions)}")

    def excluding_og(name):
        for og in exclude:
            if name.startswith(og):
                return og
        return None

    kept, excluded = [], []
    for name, start, end in partitions:
        og = excluding_og(name)
        (excluded if og else kept).append((name, start, end, og) if og else (name, start, end))

    logp(f"Kept partitions: {len(kept)}")
    logp(f"Excluded partitions: {len(excluded)}")
    for name, start, end, og in excluded:
        logp(f"  EXCLUDED ({og}, {exclude[og]}): {name}  [{start}-{end}]  len={end - start + 1}")

    matched = {og for *_, og in excluded}
    missing = [og for og in exclude if og not in matched]
    if missing:
        logp(f"NOTE: these requested orthogroups were not found as partitions: {missing}")

    seqs = read_fasta(args.fasta)
    logp(f"Taxa read from fasta: {len(seqs)}")

    code_to_header = {}
    for code in target_codes:
        matches = [h for h in seqs if h.endswith(code)]
        if len(matches) != 1:
            logp(f"ERROR: code {code} matched {len(matches)} headers: {matches}")
        else:
            code_to_header[code] = matches[0]
    if len(code_to_header) != len(target_codes):
        raise SystemExit("ABORT: not all target codes matched uniquely (see log above).")

    new_seqs = {code: [] for code in target_codes}
    new_partitions = []
    cursor = 1
    for name, start, end in kept:
        length = end - start + 1
        for code in target_codes:
            new_seqs[code].append(seqs[code_to_header[code]][start - 1:end])
        new_partitions.append((name, cursor, cursor + length - 1))
        cursor += length

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    out_fasta = args.out_prefix.with_suffix(".fasta")
    out_part = Path(str(args.out_prefix) + "_partitions.txt")
    out_log = Path(str(args.out_prefix) + "_build_log.txt")

    with open(out_fasta, "w", encoding="utf-8") as f:
        for code in target_codes:
            seq = "".join(new_seqs[code])
            f.write(f">{code_to_header[code]}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")

    with open(out_part, "w", encoding="utf-8") as f:
        for name, start, end in new_partitions:
            f.write(f"DNA, {name} = {start}-{end}\n")

    out_log.write_text("\n".join(log), encoding="utf-8")
    logp(f"\n{len(target_codes)} taxa x {cursor - 1} sites, {len(new_partitions)} partitions kept "
         f"of {len(partitions)}")
    logp(f"Wrote: {out_fasta}, {out_part}, {out_log}")


if __name__ == "__main__":
    main()
