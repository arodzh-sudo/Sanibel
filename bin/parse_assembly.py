#!/usr/bin/env python3
"""parse_assembly.py - Mash top-hit for a sample from its distances file.

Usage: parse_assembly.py <distances_file>
Output (stdout, TSV with header):
  genus  species  distance  accession
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from sanibel_taxonomy import find_accession, parse_mash_ref


def parse_mash(distances_file):
    with open(distances_file) as fh:
        line = fh.readline().strip()

    fields = line.split('\t')
    ref_id = fields[0]
    dist   = fields[2] if len(fields) > 2 else 'NA'

    genus, species, accession = parse_mash_ref(ref_id)
    if accession is None and '-.-' not in ref_id:
        accession = find_accession(ref_id)

    return genus or 'Unknown', species or 'unknown', dist, accession or 'Unknown'


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <distances_file>")

    genus, species, dist, accession = parse_mash(sys.argv[1])
    print('genus\tspecies\tdistance\taccession')
    print(f"{genus}\t{species}\t{dist}\t{accession}")
