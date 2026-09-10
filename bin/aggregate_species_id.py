#!/usr/bin/env python3
"""aggregate_species_id.py - 2-of-3 species vote across Mash, Kraken2, and 16S BLAST.

Usage: aggregate_species_id.py <assembly_stats_csv> <kraken_report> <blast_16s_tsv>

Output (stdout, TSV):
  genus  species  confidence  evidence
  where evidence = mash:G S|kraken:G S|blast16s:G sp.

Confidence values: high (all 3 agree), medium (2 agree), low (otherwise).
16S reports genus only and agrees when the winning genus is among its tied top genera.
"""

import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from sanibel_taxonomy import blast16s_top_hit, format_16s, iter_kraken_species_rows


def parse_mash(assembly_stats):
    with open(assembly_stats) as fh:
        line = fh.readline().strip()
    if not line:
        return None, None
    fields = line.split(',')
    genus   = fields[0].strip() if len(fields) > 0 else None
    species = fields[1].strip() if len(fields) > 1 else None
    if not genus or genus == 'Unknown':
        return None, None
    return genus, species


def parse_kraken(kraken_report):
    for genus, species, _pct, _reads in iter_kraken_species_rows(kraken_report):
        if genus:
            return genus, species
    return None, None


def vote(mash, kraken, blast16s_genera):
    sources = {
        'mash':     mash,
        'kraken':   kraken,
        'blast16s': (blast16s_genera[0] if blast16s_genera else None, None),
    }

    genus_counts = Counter(
        g.lower() for g, _s in sources.values() if g is not None
    )
    if not genus_counts:
        return 'Unknown', 'unknown', 'low'

    top_genus_lower = genus_counts.most_common(1)[0][0]

    winner_genus = next(
        g for g, _s in sources.values()
        if g is not None and g.lower() == top_genus_lower
    )

    agreeing = [
        s for g, s in sources.values()
        if g is not None and g.lower() == top_genus_lower and s is not None
    ]

    if not agreeing:
        winner_species = 'unknown'
    else:
        top_sp_lower   = Counter(s.lower() for s in agreeing).most_common(1)[0][0]
        winner_species = next(s for s in agreeing if s.lower() == top_sp_lower)

    if winner_species in (None, 'unknown'):
        confidence = 'low'
    else:
        agree = sum(
            1 for g, s in sources.values()
            if g is not None and s is not None
            and g.lower() == winner_genus.lower() and s.lower() == winner_species.lower()
        )
        if top_genus_lower in (g.lower() for g in blast16s_genera):
            agree += 1
        confidence = 'high' if agree >= 3 else 'medium' if agree == 2 else 'low'

    return winner_genus, winner_species, confidence


def main():
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} <assembly_stats_csv> <kraken_report> <blast_16s_tsv>")

    assembly_stats = sys.argv[1]
    kraken_report  = sys.argv[2]
    blast_tsv      = sys.argv[3]

    mash     = parse_mash(assembly_stats)
    kraken   = parse_kraken(kraken_report)
    top16s   = blast16s_top_hit(blast_tsv)
    blast16s = top16s[0] if top16s else []

    genus, species, confidence = vote(mash, kraken, blast16s)

    if blast16s:
        token16s = format_16s(blast16s)
    else:
        token16s = 'Not detected' if os.path.isfile(blast_tsv) else 'No data'
    evidence = '|'.join(
        [f"{name}:{g} {s or 'sp.'}" for name, (g, s) in (('mash', mash), ('kraken', kraken)) if g]
        + [f"blast16s:{token16s}"]
    )

    disp_genus, disp_species = ('Inconclusive', 'unknown') if confidence == 'low' else (genus, species)

    print('genus\tspecies\tconfidence\tevidence')
    print(f"{disp_genus}\t{disp_species}\t{confidence}\t{evidence}")


if __name__ == '__main__':
    main()
