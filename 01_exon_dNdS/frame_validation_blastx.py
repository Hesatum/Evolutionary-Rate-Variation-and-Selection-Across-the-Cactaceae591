#!/usr/bin/env python3
"""
frame_validation_blastx.py

Remote NCBI blastx validation of reading-frame assignments by homology
(manuscript Methods 2.2: "Reading-frame assignments were validated by
Blastx searches against the NCBI non-redundant protein database restricted
to Viridiplantae. Loci whose best hit fell in a different frame were
re-examined, and re-aligned with MACSE v2 ... and retained, whereas loci
carrying internal frameshifts shared by all taxa in their homologous frame
were excluded as likely pseudogenes.").

Submits one representative sequence per locus (concatenated into a single
multi-FASTA query file) to the public NCBI blastx CGI, polls until the job
is ready, and reports the query frame of each locus' best hit — the
homology-based cross-check on top of the internal-stop-codon screen
(reading_frame_validation.py).

Usage
-----
    python frame_validation_blastx.py --query representatives.fasta \
        --email you@example.com --out-dir results/

Notes
-----
NCBI's usage guidelines ask for a valid contact e-mail and no more than one
concurrent job per user; this script polls at 60 s intervals and is meant
to be run once per locus set, not in a loop.
"""

import argparse
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
from xml.etree import ElementTree as ET

URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'


def post(params):
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(URL, data=data)
    return urllib.request.urlopen(req, timeout=120).read().decode('utf-8', 'replace')


def get(params):
    return urllib.request.urlopen(URL + '?' + urllib.parse.urlencode(params),
                                   timeout=180).read().decode('utf-8', 'replace')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--query", required=True, type=Path,
                         help="Multi-FASTA with one representative sequence per locus.")
    parser.add_argument("--email", required=True, help="Contact e-mail required by NCBI.")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--entrez-query", default="Viridiplantae[Organism]",
                         help="ENTREZ_QUERY filter (default restricts the search to green plants).")
    parser.add_argument("--timeout-min", type=int, default=120,
                         help="Max minutes to wait for the job to finish.")
    args = parser.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / 'frame_blast.log'
    report_path = out_dir / 'frame_blast_report.tsv'
    xml_path = out_dir / 'frame_blast_raw.xml'

    def log(msg):
        line = f'[{time.strftime("%H:%M:%S")}] {msg}'
        with open(log_path, 'a', encoding='utf-8') as fh:
            fh.write(line + '\n')
        print(line, flush=True)

    query = args.query.read_text(encoding='utf-8')
    log_path.write_text('', encoding='utf-8')
    log(f'Submitting blastx job to NCBI (nr, {args.entrez_query})...')
    resp = post({'CMD': 'Put', 'PROGRAM': 'blastx', 'DATABASE': 'nr',
                 'QUERY': query, 'ENTREZ_QUERY': args.entrez_query,
                 'HITLIST_SIZE': '1', 'email': args.email, 'tool': 'frame-validation'})
    m = re.search(r'RID = (\S+)', resp)
    rtoe = re.search(r'RTOE = (\d+)', resp)
    if not m:
        log('ERROR: no RID returned. First 500 chars:\n' + resp[:500])
        return
    rid = m.group(1)
    log(f'RID={rid}  estimated time={rtoe.group(1) if rtoe else "?"}s')

    deadline = time.time() + args.timeout_min * 60
    while time.time() < deadline:
        time.sleep(60)
        info = get({'CMD': 'Get', 'FORMAT_OBJECT': 'SearchInfo', 'RID': rid})
        st = re.search(r'Status=(\w+)', info)
        status = st.group(1) if st else '?'
        log(f'poll: Status={status}')
        if status == 'READY':
            break
        if status in ('FAILED', 'UNKNOWN'):
            log(f'ERROR: job {status}')
            return
    else:
        log('ERROR: timed out waiting for results')
        return

    log('Fetching XML results...')
    xml = get({'CMD': 'Get', 'RID': rid, 'FORMAT_TYPE': 'XML'})
    xml_path.write_text(xml, encoding='utf-8')
    try:
        root = ET.fromstring(xml)
    except Exception as e:
        log(f'ERROR parsing XML: {e}')
        return

    rows = [['orthogroup', 'best_hit_query_frame', 'pident', 'aln_len', 'evalue', 'verdict', 'top_hit']]
    counts = {'frame_OK(+1)': 0, 'offset(+2/+3)': 0, 'WRONG_STRAND(neg)': 0, 'no_hit': 0}
    for it in root.iter('Iteration'):
        qdef = it.findtext('Iteration_query-def') or '?'
        qid = qdef.split()[0]
        hsp = it.find('.//Hit/Hit_hsps/Hsp')
        if hsp is None:
            rows.append([qid, 'NA', 'NA', 'NA', 'NA', 'no_hit', '-'])
            counts['no_hit'] += 1
            continue
        frame = int(hsp.findtext('Hsp_query-frame'))
        idn = int(hsp.findtext('Hsp_identity'))
        alen = int(hsp.findtext('Hsp_align-len'))
        pid = round(idn / alen * 100, 1) if alen else 0
        ev = hsp.findtext('Hsp_evalue')
        hitdef = (it.findtext('.//Hit/Hit_def') or '')[:70]
        if frame == 1:
            verdict = 'frame_OK'
            counts['frame_OK(+1)'] += 1
        elif frame in (2, 3):
            verdict = 'OFFSET_frame'
            counts['offset(+2/+3)'] += 1
        else:
            verdict = 'WRONG_STRAND'
            counts['WRONG_STRAND(neg)'] += 1
        rows.append([qid, frame, pid, alen, ev, verdict, hitdef])

    report_path.write_text('\n'.join('\t'.join(map(str, r)) for r in rows), encoding='utf-8')
    log('=== SUMMARY ===')
    for k, v in counts.items():
        log(f'  {k}: {v}')
    flagged = [r[0] for r in rows[1:] if r[5] in ('WRONG_STRAND', 'OFFSET_frame')]
    log(f'FLAGGED (not forward frame +1, candidates for MACSE re-alignment): {flagged}')
    log(f'Report written: {report_path}')
    log('DONE')


if __name__ == '__main__':
    main()
