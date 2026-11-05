#!/usr/bin/env python3
"""sanibel_taxonomy.py — shared taxonomy / contamination primitives."""

import os
import re
from collections import namedtuple

SKANI_MIN_ANI = 95.0
SKANI_MIN_AF  = 50.0

# Minimum identity/length for a 16S hit to count as a contamination candidate.
CONTAM_PIDENT = 99.0
CONTAM_LENGTH = 1400

# Minimum identity/length for a 16S hit to count as a species-ID candidate
BLAST16S_MIN_LENGTH = 400
BLAST16S_MIN_PIDENT = 98.7
# Near full-length 16S; outranks any partial hit regardless of identity
BLAST16S_FULL_LENGTH = 1400

# RefSeq/GenBank accession pattern (GCF_*, GCA_*, NC_*, NZ_*)
ACCESSION_RE = re.compile(
    r'GC[FA]_[A-Za-z0-9]+(?:\.[0-9]+)?|N[CZ]_[A-Za-z0-9]+(?:\.[0-9]+)?'
)


def find_accession(text):
    """First RefSeq/GenBank accession in `text`, or None."""
    if not text:
        return None
    m = ACCESSION_RE.search(text)
    return m.group(0) if m else None


def parse_mash_ref(ref_name):
    segs      = ref_name.split('-.-')
    name_part = re.sub(r'\.fna.*', '', segs[-1])
    tokens    = [t for t in name_part.split('_') if t]
    genus     = tokens[0] if tokens else None
    species   = tokens[1] if len(tokens) > 1 else None
    accession = None
    for seg in segs[:-1]:
        acc = find_accession(seg)
        if acc:
            accession = acc
            break
    return genus, species, accession


# Genera indistinguishable by 16S; never flag each other as contamination.
SYNONYMOUS_PAIRS = {
    frozenset({'escherichia',   'shigella'}),
    frozenset({'klebsiella',    'enterobacter'}),
    frozenset({'klebsiella',    'raoultella'}),
    frozenset({'salmonella',    'citrobacter'}),
    frozenset({'yersinia',      'serratia'}),
    frozenset({'hafnia',        'escherichia'}),
    frozenset({'haemophilus',   'aggregatibacter'}),
    frozenset({'haemophilus',   'pasteurella'}),
    frozenset({'neisseria',     'kingella'}),
    frozenset({'neisseria',     'eikenella'}),
    frozenset({'streptococcus', 'lactococcus'}),
    frozenset({'streptococcus', 'enterococcus'}),
    frozenset({'staphylococcus', 'macrococcus'}),
    frozenset({'staphylococcus', 'mammaliicoccus'}),
    frozenset({'mycobacterium', 'mycobacteroides'}),
    frozenset({'mycobacterium', 'mycolicibacterium'}),
    frozenset({'campylobacter', 'arcobacter'}),
    frozenset({'campylobacter', 'aliarcobacter'}),
    frozenset({'campylobacter', 'helicobacter'}),
    frozenset({'listeria',      'brochothrix'}),
    frozenset({'bacillus',      'paenibacillus'}),
}


def genus_of(label):
    if not label:
        return None
    tokens = label.replace('_', ' ').split()
    return tokens[0] if tokens else None


def are_synonymous(genus_a, genus_b):
    if not genus_a or not genus_b:
        return False
    return frozenset({genus_a.lower(), genus_b.lower()}) in SYNONYMOUS_PAIRS


def ranges_overlap(s1, e1, s2, e2):
    lo1, hi1 = min(s1, e1), max(s1, e1)
    lo2, hi2 = min(s2, e2), max(s2, e2)
    return lo1 <= hi2 and lo2 <= hi1


# Parsed 16S BLAST hit
BlastHit = namedtuple(
    'BlastHit', 'qseqid pident length sstart send bitscore stitle genus species'
)


def iter_blast16s_rows(path):
    try:
        fh = open(path)
    except OSError:
        return
    with fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 13 or parts[0] == 'qseqid':
                continue
            try:
                pident   = float(parts[2])
                length   = int(parts[3])
                sstart   = int(parts[8])
                send     = int(parts[9])
                bitscore = float(parts[11])
            except (ValueError, IndexError):
                continue
            stitle  = parts[12].strip()
            words   = stitle.split()
            genus   = words[0] if words else None
            species = words[1] if len(words) > 1 else None
            yield BlastHit(parts[0], pident, length, sstart, send, bitscore,
                           stitle, genus, species)


SkaniRow = namedtuple('SkaniRow', 'ani align_fraction species reference')


def read_skani_rows(path):
    """skani dist rows, highest ANI first. References are named <species>__<accession>.fna."""
    rows = []
    try:
        fh = open(path)
    except OSError:
        return rows
    with fh:
        fh.readline()
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            try:
                ani = float(parts[2])
            except ValueError:
                continue
            name = os.path.basename(parts[0])
            if name.endswith('.fna'):
                name = name[:-4]
            species, _, reference = name.partition('__')
            if not reference and len(parts) > 5 and parts[5].strip():
                reference = parts[5].split()[0]
            rows.append(SkaniRow(ani, parts[4], species or None, reference or None))
    rows.sort(key=lambda r: r.ani, reverse=True)
    return rows


def iter_kraken_species_rows(path):
    try:
        fh = open(path)
    except OSError:
        return
    with fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6 or parts[3].strip() != 'S':
                continue
            try:
                pct   = float(parts[0].strip())
                reads = int(parts[1].strip())
            except ValueError:
                continue
            words   = parts[5].strip().split()
            genus   = words[0] if words else None
            species = words[1] if len(words) > 1 else None
            yield genus, species, pct, reads


def extract_contam_candidates(blast16s_tsv_path):
    rows = []
    for hit in iter_blast16s_rows(blast16s_tsv_path):
        if hit.pident < CONTAM_PIDENT or hit.length < CONTAM_LENGTH:
            continue
        if not hit.genus:
            continue
        rows.append((hit.pident, hit.length, hit.bitscore, hit.genus,
                     hit.qseqid, hit.sstart, hit.send))

    rows.sort(key=lambda r: (-r[0], abs(r[1] - 1500), -r[2]))
    candidates = {}
    for _pi, _ln, _bs, genus, qseqid, sstart, send in rows:
        if qseqid not in candidates:
            candidates[qseqid] = (genus, sstart, send)
    return candidates


def detect_contamination(contam_candidates, self_genus, suppress_synonymous=True):
    if not self_genus:
        return 'Not screened'
    if len(contam_candidates) < 2:
        return 'None'

    def is_foreign(genus):
        if genus.lower() == self_genus.lower():
            return False
        if suppress_synonymous and are_synonymous(genus, self_genus):
            return False
        return True

    contigs = list(contam_candidates.items())
    for i in range(len(contigs)):
        qid_a, (genus_a, s_a, e_a) = contigs[i]
        for j in range(i + 1, len(contigs)):
            qid_b, (genus_b, s_b, e_b) = contigs[j]
            if ranges_overlap(s_a, e_a, s_b, e_b):
                if is_foreign(genus_a):
                    return f"Possible contamination: {genus_a} detected on contig {qid_a}"
                if is_foreign(genus_b):
                    return f"Possible contamination: {genus_b} detected on contig {qid_b}"
    return 'None'


def blast16s_rank_key(hit):
    return (hit.length < BLAST16S_FULL_LENGTH, -hit.pident, -hit.bitscore)


def blast16s_top_hit(path):
    hits = sorted(
        (h for h in iter_blast16s_rows(path)
         if h.genus and h.length >= BLAST16S_MIN_LENGTH and h.pident >= BLAST16S_MIN_PIDENT),
        key=blast16s_rank_key,
    )
    if not hits:
        return None
    tier, neg_pident, _bs = blast16s_rank_key(hits[0])
    genera = []
    for h in hits:
        if blast16s_rank_key(h)[:2] != (tier, neg_pident):
            break
        if h.genus.lower() not in (g.lower() for g in genera):
            genera.append(h.genus)
    return genera, hits[0].pident


def format_16s(genera):
    return '/'.join(genera) + ' sp.'
