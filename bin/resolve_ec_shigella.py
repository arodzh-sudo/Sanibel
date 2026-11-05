#!/usr/bin/env python3
"""resolve_ec_shigella.py - Pick Escherichia or Shigella from the ShigaTyper call."""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from sanibel_taxonomy import genus_of, read_skani_rows, SKANI_MIN_ANI, SKANI_MIN_AF


def read_prediction(path):
    if not os.path.isfile(path):
        return ''
    found_header = False
    with open(path) as f:
        for line in f:
            if found_header:
                fields = line.rstrip('\n').split('\t')
                return fields[1].strip() if len(fields) >= 2 else ''
            if line.startswith('sample\t'):
                found_header = True
    return ''


def resolve_target(prediction):
    p = prediction.lower()
    if 'not shigella' in p:
        return 'Escherichia', None
    if 'shigella' in p:
        tokens = prediction.split()
        named = f'{tokens[0]}_{tokens[1]}' if len(tokens) > 1 and tokens[1].islower() else None
        return 'Shigella', named
    if 'eiec' in p:
        return 'Escherichia', None
    return None, None


def passing_species(path):
    """Species of the skani rows that clear the ANI and align fraction cutoffs, best first."""
    species = []
    for r in read_skani_rows(path):
        try:
            af = float(r.align_fraction)
        except ValueError:
            continue
        if r.species and r.ani >= SKANI_MIN_ANI and af >= SKANI_MIN_AF:
            species.append(r.species)
    return species


def main():
    ap = argparse.ArgumentParser(description='Resolve Escherichia vs Shigella from ShigaTyper.')
    ap.add_argument('--skani',      required=True)
    ap.add_argument('--shigatyper', required=True)
    ap.add_argument('--prefix',     required=True)
    args = ap.parse_args()

    rows       = passing_species(args.skani)
    prediction = read_prediction(args.shigatyper)
    genus, named_species = resolve_target(prediction)

    resolved = ''
    if not rows:
        print(f'resolve_ec_shigella: no skani row clears {SKANI_MIN_ANI}/{SKANI_MIN_AF}',
              file=sys.stderr)
    elif genus is None:
        print(f"resolve_ec_shigella: ShigaTyper prediction '{prediction}' is not usable, "
              f'keeping {rows[0]}', file=sys.stderr)
    else:
        pick = next((sp for sp in rows if sp == named_species), None) \
            or next((sp for sp in rows if genus_of(sp) == genus), None)
        if pick is None:
            print(f'resolve_ec_shigella: ShigaTyper says {genus} but no {genus} reference '
                  f'clears {SKANI_MIN_ANI}/{SKANI_MIN_AF}, keeping {rows[0]}', file=sys.stderr)
        else:
            resolved = pick
            if pick != rows[0]:
                print(f"resolve_ec_shigella: {rows[0]} -> {pick} on ShigaTyper '{prediction}'",
                      file=sys.stderr)

    with open(f'{args.prefix}_species_resolved.txt', 'w') as fh:
        fh.write(f'{resolved}\n' if resolved else '')


if __name__ == '__main__':
    main()
