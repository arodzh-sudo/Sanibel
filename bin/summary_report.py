#!/usr/bin/env python3
"""summary_report.py - Build Sanibel summary reports from individual tool outputs."""

import argparse
import csv
import glob
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from sanibel_taxonomy import (
    genus_of, detect_contamination, extract_contam_candidates,
    blast16s_top_hit, format_16s, read_skani_rows, SKANI_MIN_ANI, SKANI_MIN_AF,
)


NO_DATA    = 'No data' # tool did not run, or its output was missing or unreadable
NOT_FOUND  = 'Not detected' # tool ran and the locus is absent from the assembly

SP_REVIEW_ANI   = 94.0
QC_MIN_COVERAGE = 40.0
QC_WARN_CONTIGS = 200
QC_FAIL_CONTIGS = 500
QC_MIN_N50      = 15000


def pmga_serotyped(skani_species):
    key = (skani_species or '').strip().replace(' ', '_').lower()
    return key.startswith(('neisseria_meningitidis', 'haemophilus_influenzae'))


def parse_assembly_stats(filepath):
    with open(filepath) as f:
        fields = f.read().strip().split(',')
    return {
        'genus':          fields[0],
        'species':        fields[1],
        'mash_distance':  fields[2],
        'accession':      fields[3],
        'num_contigs':    fields[4],
        'longest_contig': fields[5],
        'n50':            fields[6],
        'l50':            fields[7],
        'total_length':   fields[8],
        'gc_content':     fields[9],
    }


def parse_read_metrics(filepath):
    with open(filepath) as f:
        header = f.readline().strip().split()
        values = f.readline().strip().split()
    rm = dict(zip(header, values))
    return {
        'avg_read_len': rm.get('avgReadLength', NO_DATA),
        'avg_qual':     rm.get('avgQuality',    NO_DATA),
        'num_reads':    rm.get('numReads',      NO_DATA),
        'coverage':     rm.get('coverage',      NO_DATA),
    }


def parse_prokka_txt(filepath):
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split()
            if parts and parts[0] == 'CDS:':
                return parts[1]
    return NO_DATA


def parse_mlst(filepath):
    scheme = st = NO_DATA
    with open(filepath) as f:
        for line in f:
            out = line.strip().split()
            if len(out) >= 3:
                scheme = out[1]
                st     = out[2]
                if scheme in ('-', '') :
                    scheme = NOT_FOUND
                if st in ('-', ''):
                    st = NOT_FOUND
            break
    return {'scheme': scheme, 'st': st}


def parse_kraken_report(filepath):
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6 and parts[3] == 'S':
                return {'percent': parts[0].strip(), 'species': parts[5].strip()}
    return {'percent': NO_DATA, 'species': NO_DATA}


def parse_pmga_prediction(filepath):
    try:
        with open(filepath) as f:
            lines = f.readlines()
        if len(lines) > 1:
            data = lines[1].strip().split('\t')
            return data[2] if len(data) > 2 else ''
    except Exception as e:
        print(f"Warning: Could not parse PMGA file {filepath}: {e}", file=sys.stderr)
    return ''


SKANI_NO_DATA = {'ani': NO_DATA, 'confirmed_species': NO_DATA, 'align_fraction': NO_DATA, 'reference': NO_DATA}


def skani_fields(rows, resolved_species):
    row = next((r for r in rows if r.species == resolved_species), None) if resolved_species else None
    row = row or (rows[0] if rows else None)
    if row is None:
        return {'ani': 'ANI < 80%', 'confirmed_species': 'Inconclusive',
                'align_fraction': NO_DATA, 'reference': NO_DATA}
    return {'ani': f'{row.ani:.3f}', 'confirmed_species': row.species or NO_DATA,
            'align_fraction': row.align_fraction, 'reference': row.reference or NO_DATA}


def _read_text(path):
    try:
        with open(path) as f:
            return f.read().strip()
    except OSError:
        return ''


def parse_blast16s_result(filepath):
    if not os.path.isfile(filepath):
        return {'pident': NO_DATA, 'tophit': NO_DATA}
    top = blast16s_top_hit(filepath)
    if not top:
        return {'pident': NOT_FOUND, 'tophit': NOT_FOUND}
    genera, pident = top
    return {'pident': f"{pident:.3f}", 'tophit': format_16s(genera)}


AMR_TARGETS = [
    ('VIM',    r'blaVIM(-\d+\w*)?'),
    ('KPC',    r'blaKPC(-\d+\w*)?'),
    ('IMP',    r'blaIMP(-\d+\w*)?'),
    ('OXA-48', r'blaOXA-48'),
    ('NDM',    r'blaNDM(-\d+\w*)?'),
]


def amr_target_genes(genes):
    hits = [g for _, pattern in AMR_TARGETS
            for g in genes if re.fullmatch(pattern, g, re.IGNORECASE)]
    return ', '.join(hits) or 'None'


def carbapenemase_family(genes):
    hits = [label for label, pattern in AMR_TARGETS
            if any(re.fullmatch(pattern, g, re.IGNORECASE) for g in genes)]
    return ', '.join(hits) or 'None'


def _amr_field(row, *names):
    for n in names:
        val = row.get(n)
        if val is not None:
            return val.strip()
    return ''


def parse_amrfinder(filepath):
    if not os.path.isfile(filepath):
        return None
    genes, subclasses = set(), set()
    with open(filepath, newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if _amr_field(row, 'Type', 'Element type').upper() != 'AMR':
                continue
            if _amr_field(row, 'Scope').lower() != 'core':
                continue
            symbol = _amr_field(row, 'Element symbol', 'Gene symbol')
            sub    = _amr_field(row, 'Subclass')
            if symbol:
                genes.add(symbol)
            if sub and sub != 'NA':
                subclasses.add(sub)
    return {'genes': sorted(genes), 'subclasses': sorted(subclasses)}


# QC verdicts

def compute_id_qc(sk, min_ani, min_af, review_ani):
    ani = sk.get('ani', NO_DATA)
    af  = sk.get('align_fraction', NO_DATA)
    try:
        ani_v = float(ani)
        af_v  = float(af)
    except (TypeError, ValueError):
        return NO_DATA if ani == NO_DATA else f'NO ID (ANI < {min_ani:g}%)'
    if af_v < min_af:
        return f'NO ID (align fraction < {min_af:g}%)'
    if ani_v >= min_ani:
        return 'PASS'
    if ani_v >= review_ani:
        return 'REVIEW (borderline ANI)'
    return f'NO ID (ANI < {min_ani:g}%)'


def compute_assembly_qc(coverage, num_contigs, n50, min_cov, warn_contigs, fail_contigs, min_n50,
                        contaminated):
    try:
        cov     = float(coverage)
        contigs = int(float(num_contigs))
        n50_v   = int(float(n50))
    except (TypeError, ValueError):
        return 'FAIL: contamination' if contaminated else NO_DATA
    fails, warns = [], []
    if cov < min_cov:
        fails.append(f'coverage <{min_cov:g}x')
    if contigs > fail_contigs:
        fails.append(f'contigs >{fail_contigs}')
    elif contigs >= warn_contigs:
        warns.append(f'contigs >={warn_contigs}')
    if n50_v < min_n50:
        warns.append(f'N50 <{min_n50}')
    if contaminated:
        fails.append('contamination')
    if fails:
        return 'FAIL: ' + '; '.join(fails + warns)
    if warns:
        return 'Warning: ' + '; '.join(warns)
    return 'PASS'


# Per-file parsers - species-specific

def get_ecoli_serotype(sample_id):
    ecoli_json = f'{sample_id}_serotypefinder.json'
    if not os.path.isfile(ecoli_json):
        return None
    try:
        with open(ecoli_json) as f:
            data = json.load(f)
        sf = data.get('serotypefinder', {}).get('results', {})
        o_types = list(dict.fromkeys(
            hit.get('serotype', '') for hit in sf.get('O_type', {}).values() if hit.get('serotype')
        ))
        h_types = list(dict.fromkeys(
            hit.get('serotype', '') for hit in sf.get('H_type', {}).values() if hit.get('serotype')
        ))
        o_str = '/'.join(sorted(set(o_types))) if o_types else 'NT'
        h_str = '/'.join(sorted(set(h_types))) if h_types else 'NT'
        return f"{o_str}:{h_str}"
    except Exception as e:
        print(f"Warning: Could not parse SerotypeFinder JSON for {sample_id}: {e}", file=sys.stderr)
        return NO_DATA


def get_klebsiella_serotype(sample_id):
    matches = sorted(glob.glob(f'{sample_id}_kleborate_*_output.txt'))
    for path in matches:
        try:
            with open(path) as f:
                reader = csv.DictReader(f, delimiter='\t')
                if 'K_locus' not in (reader.fieldnames or []):
                    continue
                for row in reader:
                    k = (row.get('K_locus') or '').strip()
                    o = (row.get('O_locus') or '').strip()
                    if k in ('', '-') and o in ('', '-'):
                        return NOT_FOUND
                    return f"{k or '-'}/{o or '-'}"
        except Exception as e:
            print(f"Warning: Could not parse Kleborate output for {sample_id}: {e}", file=sys.stderr)
    return None


def get_legionella_serotype(sample_id):
    legsta_txt = f'{sample_id}_legsta.txt'
    if not os.path.isfile(legsta_txt):
        return None
    try:
        with open(legsta_txt) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                loci = [
                    row.get('SBT',   '-').strip(), row.get('flaA',  '-').strip(),
                    row.get('pilE',  '-').strip(), row.get('asd',   '-').strip(),
                    row.get('mip',   '-').strip(), row.get('mompS', '-').strip(),
                    row.get('proA',  '-').strip(), row.get('neuA',  '-').strip(),
                ]
                return ','.join(loci)
    except Exception as e:
        print(f"Warning: Could not parse Legsta output for {sample_id}: {e}", file=sys.stderr)
        return None


def get_salmonella_serotype(sample_id):
    salm_tsv = f'{sample_id}_seqsero2.tsv'
    if not os.path.isfile(salm_tsv):
        return None
    try:
        with open(salm_tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                profile  = row.get('Predicted antigenic profile', '').strip()
                serotype = row.get('Predicted serotype', '').strip()
                if profile:
                    return f"{profile}({serotype})" if serotype else profile
                return NOT_FOUND
    except Exception as e:
        print(f"Warning: Could not parse SeqSero2 TSV for {sample_id}: {e}", file=sys.stderr)
        return None


def get_gas_serotype(sample_id):
    gas_txt = f'{sample_id}_emm_typing.txt'
    if not os.path.isfile(gas_txt):
        return None
    try:
        validated = nonvalidated = None
        with open(gas_txt) as f:
            for line in f:
                fields = line.strip().split()
                if len(fields) >= 3:
                    emm_raw = fields[2].replace('.sds', '')
                    if fields[1] == 'EMM_validated' and validated is None:
                        validated = emm_raw
                    elif fields[1] == 'EMM_nonValidated' and nonvalidated is None:
                        nonvalidated = emm_raw
        return validated or nonvalidated or NOT_FOUND
    except Exception as e:
        print(f"Warning: Could not parse emm-typing-tool output for {sample_id}: {e}", file=sys.stderr)
        return None


def get_shigella_serotype(sample_id):
    sh_txt = f'{sample_id}_shigatyper.txt'
    if not os.path.isfile(sh_txt):
        return None
    try:
        found_header = False
        with open(sh_txt) as f:
            for line in f:
                if found_header:
                    fields = line.strip().split('\t')
                    return fields[1].strip() if len(fields) >= 2 else NOT_FOUND
                if line.startswith('sample\t'):
                    found_header = True
        return NOT_FOUND
    except Exception as e:
        print(f"Warning: Could not parse Shigatyper output for {sample_id}: {e}", file=sys.stderr)
        return NO_DATA


def get_pneumococcal_serotype(sample_id):
    pred_csv = f'{sample_id}_seroba.csv'
    if not os.path.isfile(pred_csv):
        return None
    try:
        with open(pred_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                serotype = row.get('Serotype', '').strip()
                return serotype if serotype else NOT_FOUND
    except Exception as e:
        print(f"Warning: Could not parse Seroba pred.csv for {sample_id}: {e}", file=sys.stderr)
        return None


def _parse_kaptive_txt(filepath):
    try:
        with open(filepath) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                locus = row.get('Best match locus', '').strip()
                ktype = row.get('Best match type', '').strip()
                return locus, ktype
    except Exception:
        pass
    return '', ''


def get_acinetobacter_serotype(sample_id):
    k_txt  = f'{sample_id}_ab_k.txt'
    oc_txt = f'{sample_id}_ab_oc.txt'
    if not os.path.isfile(k_txt) and not os.path.isfile(oc_txt):
        return None
    try:
        k_locus,  k_type  = _parse_kaptive_txt(k_txt)  if os.path.isfile(k_txt)  else ('', '')
        oc_locus, oc_type = _parse_kaptive_txt(oc_txt) if os.path.isfile(oc_txt) else ('', '')
        parts = []
        if k_locus or k_type:
            parts.append(f"{k_locus}({k_type})" if k_locus and k_type else k_locus or k_type)
        if oc_locus or oc_type:
            parts.append(f"{oc_locus}({oc_type})" if oc_locus and oc_type else oc_locus or oc_type)
        return '/'.join(parts) if parts else NOT_FOUND
    except Exception as e:
        print(f"Warning: Could not parse Kaptive output for {sample_id}: {e}", file=sys.stderr)
        return None


def get_vibrio_serotype(sample_id):
    k_txt = f'{sample_id}_vp_k.txt'
    o_txt = f'{sample_id}_vp_o.txt'
    if not os.path.isfile(k_txt) and not os.path.isfile(o_txt):
        return None
    try:
        k_locus, k_type = _parse_kaptive_txt(k_txt) if os.path.isfile(k_txt) else ('', '')
        o_locus, o_type = _parse_kaptive_txt(o_txt) if os.path.isfile(o_txt) else ('', '')
        parts = []
        if k_locus or k_type:
            parts.append(f"{k_locus}({k_type})" if k_locus and k_type else k_locus or k_type)
        if o_locus or o_type:
            parts.append(f"{o_locus}({o_type})" if o_locus and o_type else o_locus or o_type)
        return '/'.join(parts) if parts else NOT_FOUND
    except Exception as e:
        print(f"Warning: Could not parse Kaptive VP output for {sample_id}: {e}", file=sys.stderr)
        return None


def get_pseudomonas_serotype(sample_id):
    pasty_tsv = f'{sample_id}.tsv'
    if not os.path.isfile(pasty_tsv):
        return None
    try:
        with open(pasty_tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                stype = row.get('type', '').strip()
                return stype if stype else NOT_FOUND
    except Exception as e:
        print(f"Warning: Could not parse pasty TSV for {sample_id}: {e}", file=sys.stderr)
        return None


def get_listeria_serotype(sample_id):
    liss_txt = f'{sample_id}_lissero.txt'
    if not os.path.isfile(liss_txt):
        return None
    try:
        with open(liss_txt) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                serotype = row.get('SEROTYPE', '').strip()
                return serotype if serotype else NOT_FOUND
    except Exception as e:
        print(f"Warning: Could not parse LisSero output for {sample_id}: {e}", file=sys.stderr)
        return None


# Report headers

HEADER_STANDARD = [
    'sampleID',
    'species_id_qc', 'contamination_flag', 'assembly_qc',
    'skani_species', 'skani_ani', 'skani_align_fraction', 'skani_reference',
    'mash_species', 'mash_reference', 'mash_distance',
    'kraken2_species', 'kraken2_percent',
    'blast_16s_tophit', 'blast_16s_pident',
    'serotype', 'mlst_scheme', 'mlst_st',
    'num_clean_reads', 'avg_read_length', 'avg_read_qual', 'est_coverage',
    'num_contigs', 'longest_contig', 'N50', 'L50', 'total_length', 'gc_content',
    'annotated_cds',
]

HEADER_AMR = ['sampleID', 'carbapenemase_family', 'amr_target',
              'amr_genes', 'amr_subclass']

# MultiQC custom-content emitters

MQC_SPECIES_HEADER = ['Sample', 'skani_species', 'skani_ani', 'skani_align_fraction',
                      'species_id_qc', 'contamination_flag', 'assembly_qc']
MQC_TYPING_HEADER  = ['Sample', 'mlst_scheme', 'mlst_st', 'serotype']

MQC_AMR_HEADER     = ['Sample', 'carbapenemase_family', 'amr_target',
                      'amr_gene_count', 'amr_genes', 'amr_subclass']


def std_cols(mqc_header):
    return [0] + [HEADER_STANDARD.index(c) for c in mqc_header[1:]]


def _mqc_preamble(section_id, section_name, description, pconfig=None, headers=None):
    lines = [
        f"# id: '{section_id}'",
        f"# section_name: '{section_name}'",
        f"# description: '{description}'",
        "# plot_type: 'table'",
    ]
    if pconfig:
        lines.append("# pconfig:")
        for k, v in pconfig.items():
            val = str(v).lower() if isinstance(v, bool) else f"'{v}'"
            lines.append(f"#     {k}: {val}")
    if headers:
        lines.append("# headers:")
        for col, opts in headers.items():
            lines.append(f"#     {col}:")
            for k, v in opts.items():
                val = str(v).lower() if isinstance(v, bool) else f"'{v}'"
                lines.append(f"#         {k}: {val}")
    return lines


def _write_mqc(path, preamble_lines, header, cols, rows):
    with open(path, 'w') as fh:
        for pl in preamble_lines:
            fh.write(pl + '\n')
        fh.write('\t'.join(header) + '\n')
        for row in sorted(rows, key=lambda r: r[0]):
            fh.write('\t'.join(str(row[i]) for i in cols) + '\n')
    print(f"summary_report.py: wrote {path} ({len(rows)} sample(s))")


def emit_sanibel_mqc_tables(rows_std):
    _write_mqc(
        'sanibel_species_mqc.tsv',
        _mqc_preamble(
            'sanibel_species', 'Species ID and QC',
            'skani ANI-confirmed consensus species call with QC verdicts.',
            pconfig={'id': 'sanibel_species_table', 'col1_header': 'Sample',
                     'no_violin': True},
        ),
        MQC_SPECIES_HEADER, std_cols(MQC_SPECIES_HEADER), rows_std,
    )
    _write_mqc(
        'sanibel_typing_mqc.tsv',
        _mqc_preamble(
            'sanibel_typing', 'MLST and Serotyping',
            'MLST scheme and sequence type, and species-specific serotype.',
            pconfig={'id': 'sanibel_typing_table', 'col1_header': 'Sample',
                     'no_violin': True, 'only_defined_headers': False},
            headers={'mlst_st': {'scale': False, 'format': '{}'}},
        ),
        MQC_TYPING_HEADER, std_cols(MQC_TYPING_HEADER), rows_std,
    )


def amr_row(sid, amr):
    if amr is None:
        return [sid, NO_DATA, NO_DATA, NO_DATA, NO_DATA]
    genes = amr['genes']
    if not genes:
        return [sid, 'None', 'None', 'None', 'None']
    return [sid, carbapenemase_family(genes), amr_target_genes(genes),
            ', '.join(genes), ', '.join(amr['subclasses']) or 'None']


def emit_sanibel_amr_mqc_table(amr_by_sample):
    if not amr_by_sample:
        return
    rows = []
    for sid, amr in amr_by_sample.items():
        row = amr_row(sid, amr)
        rows.append(row[:3] + [NO_DATA if amr is None else len(amr['genes'])] + row[3:])
    _write_mqc(
        'sanibel_amr_mqc.tsv',
        _mqc_preamble(
            'sanibel_amr', 'Antimicrobial Resistance',
            'Acquired AMR determinants from AMRFinderPlus. Stress and virulence '
            'elements are excluded.',
            pconfig={'id': 'sanibel_amr_table', 'col1_header': 'Sample',
                     'no_violin': True},
        ),
        MQC_AMR_HEADER, [0, 1, 2, 3, 4, 5], rows,
    )


# Main

def main():
    parser = argparse.ArgumentParser(description='Build Sanibel summary reports.')
    parser.add_argument('--samples', nargs='+', required=True, help='Every submitted sample ID')
    args = parser.parse_args()

    samples = sorted(set(args.samples))

    rows_std = []
    amr_by_sample = {}

    for sid in samples:
        if not os.path.isfile(f'{sid}_assembly_stats.txt'):
            rows_std.append([sid, NO_DATA, NO_DATA, 'FAIL: no assembly']
                            + [NO_DATA] * (len(HEADER_STANDARD) - 4))
            continue

        amr_by_sample[sid] = parse_amrfinder(f'{sid}_amrfinderplus_report.tsv')

        asm  = parse_assembly_stats(f'{sid}_assembly_stats.txt')
        rm   = parse_read_metrics(f'{sid}_readMetrics.txt')
        cds  = parse_prokka_txt(f'{sid}.txt')
        mlst = parse_mlst(f'{sid}.mlst')
        kr   = parse_kraken_report(f'{sid}.report')

        skani_path = f'{sid}_skani.tsv'
        resolved   = _read_text(f'{sid}_species_resolved.txt')
        sk = (skani_fields(read_skani_rows(skani_path), resolved)
              if os.path.isfile(skani_path) else SKANI_NO_DATA)

        skani_ID_val     = sk['confirmed_species']
        skani_ANI_val    = sk['ani']
        skani_align_val  = sk['align_fraction']
        skani_ID_ref_val = sk['reference']

        skani_genus  = genus_of(skani_ID_val) if skani_ID_val not in (NO_DATA, 'NO ID', 'Inconclusive', '', None) else None

        blast16s_path   = f'{sid}_16s_blast.tsv'
        blast16s_result = parse_blast16s_result(blast16s_path)

        if skani_genus and os.path.isfile(blast16s_path):
            contamination_flag = detect_contamination(
                extract_contam_candidates(blast16s_path), skani_genus)
        else:
            contamination_flag = 'Not screened'

        scheme = mlst['scheme']

        if pmga_serotyped(skani_ID_val):
            serotype = parse_pmga_prediction(f'{sid}sta.txt') or NO_DATA
        else:
            serotype = NO_DATA
            complex_getters = [
                lambda: get_ecoli_serotype(sid),
                lambda: get_shigella_serotype(sid),
            ]
            if skani_genus == 'Shigella':
                complex_getters.reverse()
            for getter in complex_getters + [
                lambda: get_klebsiella_serotype(sid),
                lambda: get_legionella_serotype(sid),
                lambda: get_salmonella_serotype(sid),
                lambda: get_gas_serotype(sid),
                lambda: get_pneumococcal_serotype(sid),
                lambda: get_acinetobacter_serotype(sid),
                lambda: get_vibrio_serotype(sid),
                lambda: get_pseudomonas_serotype(sid),
                lambda: get_listeria_serotype(sid),
            ]:
                result = getter()
                if result is not None:
                    serotype = result
                    break

        species_id_qc = compute_id_qc(sk, SKANI_MIN_ANI, SKANI_MIN_AF, SP_REVIEW_ANI)
        in_complex = skani_ID_val == 'Escherichia_coli' or skani_genus == 'Shigella'
        if species_id_qc == 'PASS' and in_complex and not resolved:
            species_id_qc = 'REVIEW (E. coli/Shigella unresolved)'
        contaminated = contamination_flag not in ('None', 'Not screened', NO_DATA)
        assembly_qc = compute_assembly_qc(
            rm['coverage'], asm['num_contigs'], asm['n50'],
            QC_MIN_COVERAGE, QC_WARN_CONTIGS, QC_FAIL_CONTIGS, QC_MIN_N50, contaminated)

        std_row = [
            sid,
            species_id_qc, contamination_flag, assembly_qc,
            skani_ID_val, skani_ANI_val, skani_align_val, skani_ID_ref_val,
            f"{asm['genus']}_{asm['species']}", asm['accession'], asm['mash_distance'],
            kr['species'], kr['percent'],
            blast16s_result['tophit'], blast16s_result['pident'],
            serotype, scheme, mlst['st'],
            rm['num_reads'], rm['avg_read_len'], rm['avg_qual'], rm['coverage'],
            asm['num_contigs'], asm['longest_contig'], asm['n50'], asm['l50'],
            asm['total_length'], asm['gc_content'], cds,
        ]
        rows_std.append(std_row)

    def write_report(path, header, rows):
        with open(path, 'w', encoding='utf-16', newline='\r\n') as fh:
            fh.write('\t'.join(header) + '\n')
            for row in sorted(rows, key=lambda r: r[0]):
                fh.write('\t'.join(str(v).replace(',', ';') for v in row) + '\n')
        print(f"summary_report.py: wrote {path} ({len(rows)} sample(s))")

    rows_amr = [amr_row(sid, amr) for sid, amr in amr_by_sample.items()]

    write_report('sum_report.txt', HEADER_STANDARD, rows_std)
    if rows_amr:
        write_report('amr_report.txt', HEADER_AMR, rows_amr)

    emit_sanibel_mqc_tables(rows_std)
    emit_sanibel_amr_mqc_table(amr_by_sample)


if __name__ == '__main__':
    main()
