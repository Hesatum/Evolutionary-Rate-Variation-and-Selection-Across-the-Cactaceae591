#!/usr/bin/env python3
"""
reading_frame_validation.py

Six-frame translation and internal-stop-codon screen for codon alignments,
used before dN/dS analysis (manuscript Methods 2.2: "Each retained sequence
was then translated in all six reading frames with a custom Python script
with standard genetic code, and internal stop codons were counted per
frame. Sequences showing internal stop codons in the forward reading frame
(frame 0) were subsequently evaluated across alternative reading frames.").

For each alignment file, every sequence is translated in all six frames
(forward 0/1/2 and reverse-complement 0/1/2). A locus is:
  - OK               if the forward frame-0 translation is stop-free;
  - re-framed        if some other frame is stop-free (the sequence is
                      reverse-complemented and/or offset accordingly);
  - flagged          if no frame is stop-free (kept out of the corrected
                      set; candidates for MACSE re-alignment or exclusion,
                      see Methods 2.2).

Writes a combined TSV report and a folder with the OK + re-framed
sequences ready for the next step (frame validation via Blastx, see
frame_validation_blastx.py).

Usage
-----
    python reading_frame_validation.py --input-dirs alignments/ alignments/temp \
        --output-dir genes_finais/ --report frame_check_report.tsv
"""

import argparse
import shutil
from pathlib import Path

_COMP = str.maketrans(
    'ACGTacgtRYSWKMBDHVNrysWkmbdhvn-?',
    'TGCAtgcaYRSWMKVHDBNyrswmkvhdbN-?'
)


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


STOP_CODONS = frozenset({'TAA', 'TAG', 'TGA'})


def count_internal_stops(seq: str, frame: int = 0) -> int:
    """Count internal stop codons in `seq` starting at `frame`.

    Gaps/`?` are stripped before slicing into codons; a terminal stop codon
    is not counted (it is the legitimate end of the CDS), and codons
    containing 'N' are treated as ambiguous, not stop.
    """
    clean = seq.upper().replace('-', '').replace('?', '')
    codons = [clean[i:i + 3] for i in range(frame, len(clean) - 2, 3)]
    if not codons:
        return 0
    check = codons[:-1] if codons[-1] in STOP_CODONS else codons
    return sum(1 for c in check if len(c) == 3 and 'N' not in c and c in STOP_CODONS)


def read_fasta(path: Path) -> list[tuple[str, str]]:
    seqs, header, parts = [], None, []
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if header is not None:
                    seqs.append((header, ''.join(parts)))
                header, parts = line[1:], []
            elif line:
                parts.append(line)
    if header is not None:
        seqs.append((header, ''.join(parts)))
    return seqs


def write_fasta(path: Path, seqs: list[tuple[str, str]], width: int = 60) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', encoding='utf-8') as fh:
        for header, seq in seqs:
            fh.write(f'>{header}\n')
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + '\n')


FRAME_KEYS = [('fwd', 0), ('fwd', 1), ('fwd', 2), ('rc', 0), ('rc', 1), ('rc', 2)]
FRAME_LABEL = {
    ('fwd', 0): 'forward frame 0', ('fwd', 1): 'forward frame 1', ('fwd', 2): 'forward frame 2',
    ('rc', 0): 'reverse-complement frame 0', ('rc', 1): 'reverse-complement frame 1',
    ('rc', 2): 'reverse-complement frame 2',
}


def analyze(path: Path) -> dict | None:
    """Translate every sequence in `path` in all 6 frames and classify the locus."""
    seqs = read_fasta(path)
    if not seqs:
        return None

    ref_len = len(seqs[0][1].replace('-', '').replace('?', ''))

    max_stops: dict[tuple, int] = {}
    for strand, off in FRAME_KEYS:
        worst = 0
        for _, seq in seqs:
            s = reverse_complement(seq) if strand == 'rc' else seq
            n = count_internal_stops(s, off)
            worst = max(worst, n)
        max_stops[(strand, off)] = worst

    fwd0 = max_stops[('fwd', 0)]
    if fwd0 == 0:
        classification, clean_key, clean_label = 'OK', None, '-'
    else:
        clean_key = next((k for k in FRAME_KEYS[1:] if max_stops[k] == 0), None)
        if clean_key:
            clean_label, classification = FRAME_LABEL[clean_key], 'REFRAMED (strand/frame corrected)'
        else:
            clean_label, classification = 'none', 'FLAGGED (no clean ORF in any frame)'

    return dict(filepath=path, filename=path.name, n_seqs=len(seqs), ref_len=ref_len,
                stops_fwd0=fwd0, clean_key=clean_key, clean_label=clean_label,
                classification=classification, sequences=seqs)


def correct_seq(seq: str, strand: str, offset: int) -> str:
    if strand == 'rc':
        seq = reverse_complement(seq)
    if offset:
        seq = seq[offset:]
    seq = seq.replace('-', '').replace('?', '')
    extra = len(seq) % 3
    return seq[:-extra] if extra else seq


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-dirs", nargs="+", required=True, type=Path,
                         help="One or more directories with per-locus FASTA alignments (*.fas/*.fasta).")
    parser.add_argument("--output-dir", required=True, type=Path,
                         help="Directory for the OK + re-framed sequences ready for the next step.")
    parser.add_argument("--report", required=True, type=Path,
                         help="Path to the combined TSV report.")
    parser.add_argument("--skip", nargs="*", default=[],
                         help="File stems to skip (e.g. pre-existing compilation files).")
    args = parser.parse_args()

    all_results = []
    for src_dir in args.input_dirs:
        fas_files = sorted(src_dir.glob('*.fas')) + sorted(src_dir.glob('*.fasta'))
        fas_files = [f for f in fas_files if f.stem not in args.skip]
        print(f'[{src_dir}] {len(fas_files)} file(s) found')
        for fp in fas_files:
            r = analyze(fp)
            if r:
                r['source'] = str(src_dir)
                all_results.append(r)

    args.report.parent.mkdir(parents=True, exist_ok=True)
    with open(args.report, 'w', encoding='utf-8') as fh:
        fh.write('source\tfile\tn_seqs\tref_length\tstops_fwd_frame0\tclean_frame\tclassification\n')
        for r in all_results:
            fh.write(f"{r['source']}\t{r['filename']}\t{r['n_seqs']}\t{r['ref_len']}\t"
                      f"{r['stops_fwd0']}\t{r['clean_label']}\t{r['classification']}\n")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    copied_ok = corrected = skipped = 0
    for r in all_results:
        dest = args.output_dir / r['filename']
        if r['classification'] == 'OK':
            shutil.copy2(r['filepath'], dest)
            copied_ok += 1
        elif r['clean_key'] is not None:
            strand, off = r['clean_key']
            corrected_seqs = [(hdr, correct_seq(seq, strand, off)) for hdr, seq in r['sequences']]
            write_fasta(dest, corrected_seqs)
            corrected += 1
        else:
            skipped += 1

    print(f'\nTotal analyzed: {len(all_results)}')
    print(f'  OK               : {copied_ok}')
    print(f'  Re-framed        : {corrected}')
    print(f'  Flagged (no ORF) : {skipped}')
    print(f'\nCorrected/validated set: {args.output_dir}')
    print(f'Report: {args.report}')


if __name__ == '__main__':
    main()
