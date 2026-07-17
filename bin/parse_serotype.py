#!/usr/bin/env python3
"""parse_serotype.py <tool> <sample_id> - print the serotype for one typing tool.

Reads the tool's output staged in the current directory as ./<tool>/... and prints
the normalized serotype string, reusing summary_report's per-tool getters.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from summary_report import (
    get_ecoli_serotype, get_klebsiella_serotype, get_legionella_serotype,
    get_salmonella_serotype, get_gas_serotype, get_shigella_serotype,
    get_pneumococcal_serotype, get_acinetobacter_serotype, get_vibrio_serotype,
    get_pseudomonas_serotype, get_listeria_serotype,
)

DISPATCH = {
    'serotypefinder': get_ecoli_serotype,
    'kleborate':      get_klebsiella_serotype,
    'legsta':         get_legionella_serotype,
    'seqsero2':       get_salmonella_serotype,
    'emm_typing':     get_gas_serotype,
    'shigatyper':     get_shigella_serotype,
    'seroba':         get_pneumococcal_serotype,
    'kaptive_ab':     get_acinetobacter_serotype,
    'kaptive_vp':     get_vibrio_serotype,
    'pasty':          get_pseudomonas_serotype,
    'lissero':        get_listeria_serotype,
}


def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <tool> <sample_id>")
    tool, sid = sys.argv[1], sys.argv[2]
    getter = DISPATCH.get(tool)
    value = getter('.', sid) if getter else None
    print(value if value is not None else 'Not detected')


if __name__ == '__main__':
    main()
