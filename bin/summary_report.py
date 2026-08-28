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
    iter_blast16s_rows, BLAST16S_MIN_LENGTH, BLAST16S_MIN_PIDENT,
)


NO_DATA    = 'No data'                       # tool did not run, or its output was missing or unreadable
NOT_FOUND  = 'Not detected'                  # tool ran and the locus is absent from the assembly
UNRESOLVED = 'Present (allele unresolved)'   # locus is in the genome, allele lookup did not resolve
DISRUPTED  = 'Disrupted ORF'                 # gene is present but not intact
NEW_ALLELE = 'New allele'                    # present, no database match

SP_MIN_ANI      = 95.0
SP_MIN_AF       = 50.0
SP_REVIEW_ANI   = 94.0
QC_MIN_COVERAGE = 40.0
QC_WARN_CONTIGS = 200
QC_FAIL_CONTIGS = 500
QC_MIN_N50      = 15000


def normalize_le_value(val):
    if not val:
        return NO_DATA
    if val in ['Not found', 'not found'] or val.startswith('Peptide not found'):
        return NOT_FOUND
    if val.startswith('Allele not identified'):
        return UNRESOLVED
    if val.lower().startswith('incomplete orf'):
        return DISRUPTED
    if val.startswith('New-BLASTonly') or val.startswith('New-PCR'):
        return NEW_ALLELE
    if val.startswith('Error'):
        return NO_DATA
    if val == 'None identified':
        return NOT_FOUND
    return val


def antigen_present(*values):
    """A locus is in the genome unless the lookup said it was absent or said nothing at all."""
    return any(v not in (NO_DATA, NOT_FOUND) for v in values)


def read_bmscan_species(sample_id):
    bmscan_jsons = glob.glob(f'{sample_id}_species_analysis.json')
    if not bmscan_jsons:
        return None
    try:
        with open(bmscan_jsons[0]) as f:
            bmscan = json.load(f)
        for sample_data in bmscan.values():
            if 'mash_results' in sample_data:
                return sample_data['mash_results'].get('species') or None
    except Exception as e:
        print(f"Warning: Could not parse BMScan JSON for {sample_id}: {e}", file=sys.stderr)
    return None


def meningitis_organism(skani_species, bmscan_species):
    """Keyed on species, never the MLST scheme: mlst maps the whole Neisseria genus to 'neisseria'."""
    for label in (skani_species, bmscan_species):
        key = (label or '').strip().replace(' ', '_').lower()
        if key.startswith('neisseria_meningitidis'):
            return 'neisseria'
        if key.startswith('haemophilus_influenzae'):
            return 'hinfluenzae'
    return None


def mutation_drugs(gene):
    """The drugs BMGAP2 itself attaches to each curated mutation, deduplicated."""
    drugs = set()
    for entry in (gene.get('known_mutations') or {}).values():
        drugs.update(d.strip() for d in (entry.get('resistance') or '').split(';') if d.strip())
    return drugs


def read_gene(genes, name):
    """(allele, mutations, phenotype) keyed on BMGAP2's own per-gene status.

    runAST screens a different gene set per species and records a real absence as
    status Absent, so a key missing from the JSON was never screened."""
    gene = find_gene(genes, name)
    if gene is None:
        return NO_DATA, NO_DATA, NO_DATA

    status = gene.get('status')
    if status == 'Absent':
        return NOT_FOUND, NOT_FOUND, NOT_FOUND
    if status != 'Present':
        return gene.get('allele') or NO_DATA, NO_DATA, DISRUPTED

    allele = gene.get('allele') or NO_DATA
    muts = gene.get('known_mutations')
    if muts is None:
        return allele, NO_DATA, NO_DATA
    if not muts:
        return allele, 'None', 'Susceptible'
    return allele, ';'.join(muts), 'Resistant'


def find_gene(gene_dict, gene_name):
    for value in gene_dict.values():
        if value.get('Gene_name') == gene_name:
            return value
    for key, value in gene_dict.items():
        if key == gene_name or f"({gene_name})" in key or key.startswith(gene_name + ' '):
            return value
    return None


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


def lookup_cc(table_path, st, missing_col='NA', default=''):
    try:
        with open(table_path) as tbl:
            for row in tbl:
                cols = row.strip().split('\t')
                if cols[0] == st:
                    return cols[8] if len(cols) >= 9 else missing_col
    except Exception:
        pass
    return default


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


def parse_pmga(filepath):
    result = {'species': '', 'prediction': '', 'serotype_notes': '-'}
    try:
        if os.path.getsize(filepath) > 0:
            with open(filepath) as f:
                lines = f.readlines()
            if len(lines) > 1:
                data = lines[1].strip().split('\t')
                result['species']        = data[1] if len(data) > 1 else ''
                result['prediction']     = data[2] if len(data) > 2 else ''
                result['serotype_notes'] = data[4] if len(data) > 4 else '-'
    except Exception as e:
        print(f"Warning: Could not parse PMGA file {filepath}: {e}", file=sys.stderr)
    return result


def parse_skani(filepath):
    EMPTY  = {'ani': 'ANI < 80%', 'confirmed_species': 'Inconclusive', 'align_fraction': NO_DATA, 'reference': NO_DATA}
    FAILED = {'ani': NO_DATA, 'confirmed_species': NO_DATA, 'align_fraction': NO_DATA, 'reference': NO_DATA}
    try:
        rows = []
        with open(filepath) as f:
            f.readline()
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 5:
                    continue
                try:
                    ani = float(parts[2])
                except ValueError:
                    continue
                rows.append((ani, parts))
        if not rows:
            return EMPTY
        rows.sort(key=lambda x: x[0], reverse=True)
        best_ani, best_parts = rows[0]
        ref_basename = os.path.basename(best_parts[0])
        if ref_basename.endswith('.fna'):
            ref_basename = ref_basename[:-4]
        if '__' in ref_basename:
            species_part, acc_part = ref_basename.split('__', 1)
            confirmed_species = species_part if species_part else NO_DATA
            skani_ref_acc     = acc_part if acc_part else NO_DATA
        else:
            confirmed_species = ref_basename if ref_basename else NO_DATA
            ref_name_col      = best_parts[5].strip() if len(best_parts) > 5 else ''
            skani_ref_acc     = ref_name_col.split()[0] if ref_name_col else NO_DATA
        align_fraction = best_parts[4] if len(best_parts) > 4 else NO_DATA
        return {
            'ani':               f"{best_ani:.3f}",
            'confirmed_species': confirmed_species,
            'align_fraction':    align_fraction,
            'reference':         skani_ref_acc if skani_ref_acc else NO_DATA,
        }
    except Exception as e:
        print(f"Warning: Could not parse skani TSV {filepath}: {e}", file=sys.stderr)
        return FAILED


def _read_text(path):
    try:
        with open(path) as f:
            return f.read().strip()
    except OSError:
        return ''


def apply_resolved_species(sk, skani_path, resolved_path):
    """Re-point the skani fields at the highest-ANI row for the ShigaTyper-resolved species."""
    try:
        with open(resolved_path) as f:
            species = f.read().strip()
    except OSError:
        return sk
    if not species:
        return sk

    best = None
    try:
        with open(skani_path) as f:
            f.readline()
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 5:
                    continue
                name = os.path.basename(parts[0])
                if name.endswith('.fna'):
                    name = name[:-4]
                if name.split('__', 1)[0] != species:
                    continue
                try:
                    ani = float(parts[2])
                except ValueError:
                    continue
                if best is None or ani > best[0]:
                    best = (ani, parts[4], name.split('__', 1)[1] if '__' in name else NO_DATA)
    except OSError:
        return sk
    if best is None:
        return sk

    ani, align_fraction, reference = best
    return {'ani': f'{ani:.3f}', 'confirmed_species': species,
            'align_fraction': align_fraction, 'reference': reference}


def parse_blast16s_result(filepath, anchor_genera=None):
    qualifying = []
    for hit in iter_blast16s_rows(filepath):
        if hit.length < BLAST16S_MIN_LENGTH or hit.pident < BLAST16S_MIN_PIDENT:
            continue
        qualifying.append((hit.pident, hit.length, hit.genus or '', hit.bitscore))
    if not qualifying:
        return {'pident': NO_DATA, 'tophit': NO_DATA}

    qualifying.sort(key=lambda r: (-r[0], abs(r[1] - 1500), -r[3]))

    anchor = {g.lower() for g in (anchor_genera or []) if g and g not in (NO_DATA, 'Unknown')}
    chosen = next((r for r in qualifying if r[2].lower() in anchor), qualifying[0]) if anchor \
             else qualifying[0]

    pident, _len, genus, _bs = chosen
    return {'pident': f"{pident:.3f}", 'tophit': f"{genus} spp." if genus else NO_DATA}


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
        return NO_DATA
    fails, warns = [], []
    if cov < min_cov:
        fails.append(f'coverage <{min_cov:g}x')
    if contigs > fail_contigs:
        fails.append(f'contigs >{fail_contigs}')
    elif contigs >= warn_contigs:
        warns.append(f'contigs >={warn_contigs}')
    if n50_v < min_n50:
        warns.append(f'N50 <{min_n50}')
    if fails:
        return 'FAIL: ' + '; '.join(fails + warns)
    if contaminated:
        return 'FAIL (Contamination)'
    if warns:
        return 'Warning: ' + '; '.join(warns)
    return 'PASS'


# Per-file parsers - species-specific

def get_ecoli_serotype(sample_dir, sample_id):
    ecoli_json = os.path.join(sample_dir, 'serotypefinder', 'data.json')
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


def get_klebsiella_serotype(sample_dir, sample_id):
    matches = sorted(glob.glob(os.path.join(sample_dir, 'kleborate', 'kleborate_out',
                                            '*_output.txt')))
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


def get_legionella_serotype(sample_dir, sample_id):
    legsta_txt = os.path.join(sample_dir, 'legsta', 'legsta_output.txt')
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


def get_salmonella_serotype(sample_dir, sample_id):
    salm_dirs = glob.glob(os.path.join(sample_dir, 'seqsero2', 'SeqSero_result_*'))
    if not salm_dirs:
        return None
    salm_tsv = os.path.join(salm_dirs[0], 'SeqSero_result.tsv')
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


def get_gas_serotype(sample_dir, sample_id):
    gas_txt = os.path.join(sample_dir, 'emm_typing', 'groupAstrep_result.txt')
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


def get_shigella_serotype(sample_dir, sample_id):
    sh_txt = os.path.join(sample_dir, 'shigatyper', 'shigatyper_output.txt')
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


def get_pneumococcal_serotype(sample_dir, sample_id):
    pred_csv = os.path.join(sample_dir, 'seroba', sample_id, 'pred.csv')
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
    """Return (locus, type) from a kaptive v3 TSV output file."""
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


def get_acinetobacter_serotype(sample_dir, sample_id):
    k_txt  = os.path.join(sample_dir, 'kaptive_ab', f'{sample_id}_ab_k.txt')
    oc_txt = os.path.join(sample_dir, 'kaptive_ab', f'{sample_id}_ab_oc.txt')
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


def get_vibrio_serotype(sample_dir, sample_id):
    k_txt = os.path.join(sample_dir, 'kaptive_vp', f'{sample_id}_vp_k.txt')
    o_txt = os.path.join(sample_dir, 'kaptive_vp', f'{sample_id}_vp_o.txt')
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


def get_pseudomonas_serotype(sample_dir, sample_id):
    pasty_tsv = os.path.join(sample_dir, 'pasty', f'{sample_id}.tsv')
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


def get_listeria_serotype(sample_dir, sample_id):
    liss_txt = os.path.join(sample_dir, 'lissero', 'lissero_output.txt')
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


# BMGAP2 data parser

def parse_bmgap2(sample_id, organism, bmscan_species, hinfluenzae_txt=None):
    d = {k: NO_DATA for k in [
        'penA_allele', 'penA_mutations', 'penA_phenotype',
        'gyrA_allele', 'gyrA_mutations', 'gyrA_phenotype',
        'parC_allele', 'parC_phenotype',
        'rpoB_allele', 'rpoB_phenotype',
        'ponA_allele', 'ponA_phenotype',
        'predicted_resistance', 'resistance_markers',
        'bmgap2_species', 'bmgap2_mlst_st', 'bmgap2_mlst_cc',
        'FHbp_variant', 'FHbp_subfamily', 'FHbp_peptide', 'FHbp_Pfizer',
        'NadA_variant', 'NhbA_peptide', 'PorA_type', 'vaccine_antigens_present',
        'folA_allele', 'folA_phenotype',
        'acrR_allele', 'acrR_mutations', 'acrR_phenotype',
        'blaTEM1_status', 'blaROB1_status',
    ]}
    status = {'amr': 'no_output', 'le': 'no_output', 'bmscan': 'no_output'}

    # runAST names its output after the staged PMGA JSON, so the real name is <id>sta-blast_amr_data.json
    amr_jsons = sorted(set(glob.glob(f'{sample_id}sta*_amr_data.json')
                           + glob.glob(f'{sample_id}_amr_data.json')))
    if len(amr_jsons) > 1:
        print(f"Warning: {len(amr_jsons)} AMR JSONs match {sample_id}, using {amr_jsons[0]}",
              file=sys.stderr)
    if amr_jsons:
        try:
            with open(amr_jsons[0]) as f:
                amr = json.load(f)
            status['amr'] = amr.get('bmgap2_status', 'ok')
            genes = amr.get('amr_genes', {})

            d['penA_allele'], d['penA_mutations'], _ = read_gene(genes, 'penA')
            pen_pheno = (amr.get('antimicrobics', {})
                            .get('Penicillins', {}).get('Penicillin', {}))
            if pen_pheno:
                d['penA_phenotype'] = pen_pheno.get('predicted_phenotype', NO_DATA)

            d['gyrA_allele'], d['gyrA_mutations'], d['gyrA_phenotype'] = read_gene(genes, 'gyrA')
            d['parC_allele'], _, d['parC_phenotype'] = read_gene(genes, 'parC')
            d['rpoB_allele'], _, d['rpoB_phenotype'] = read_gene(genes, 'rpoB')
            d['ponA_allele'], _, d['ponA_phenotype'] = read_gene(genes, 'ponA')

            # runAST screens both beta-lactamases for Nm as well as Hi
            for gene, field in (('blaTEM-1', 'blaTEM1_status'), ('blaROB-1', 'blaROB1_status')):
                rec = find_gene(genes, gene)
                d[field] = (rec.get('status') or NO_DATA) if rec else NO_DATA

            if organism == 'hinfluenzae':
                d['penA_allele'], d['penA_mutations'], ftsi_pheno = read_gene(genes, 'ftsI')
                pen_dict = amr.get('antimicrobics', {}).get('Penicillins', {})
                amp = pen_dict.get('Ampicillin') or (next(iter(pen_dict.values())) if pen_dict else None)
                d['penA_phenotype'] = amp.get('predicted_phenotype', NO_DATA) if amp else ftsi_pheno

                d['folA_allele'], _, d['folA_phenotype'] = read_gene(genes, 'folA')
                d['acrR_allele'], d['acrR_mutations'], d['acrR_phenotype'] = read_gene(genes, 'acrR')

            if genes:
                drugs = set()
                for rec in genes.values():
                    drugs |= mutation_drugs(rec)
                d['resistance_markers'] = ';'.join(sorted(drugs)) if drugs else NOT_FOUND

            d['predicted_resistance'] = amr.get('summary', {}).get('predicted_resistance', NO_DATA)

        except Exception as e:
            status['amr'] = 'failed'
            print(f"Warning: Could not parse BMGAP2 AMR JSON for {sample_id}: {e}", file=sys.stderr)

    # LocusExtractor CSV
    le_dirs = glob.glob(f'LE_*_{sample_id}_*')
    if le_dirs:
        le_csv = os.path.join(le_dirs[0], 'Results_text', f'molecular_data_{sample_id}.csv')
        status['le'] = 'ran_no_csv'
        if os.path.isfile(le_csv):
            try:
                with open(le_csv) as f:
                    reader = csv.DictReader(f)
                    all_rows = list(reader)

                prokka_rows = [r for r in all_rows if 'prokka' in r.get('Filename', '')]
                row = prokka_rows[0] if prokka_rows else (all_rows[0] if all_rows else None)
                if row is not None:
                    if organism == 'hinfluenzae':
                        d['bmgap2_mlst_st'] = normalize_le_value(row.get('Hi_MLST_ST', ''))
                        hi_st = d['bmgap2_mlst_st']
                        d['bmgap2_mlst_cc'] = NO_DATA
                        if (hinfluenzae_txt and os.path.isfile(hinfluenzae_txt)
                                and hi_st not in [NO_DATA, NOT_FOUND, 'New', 'NA', '']):
                            d['bmgap2_mlst_cc'] = lookup_cc(
                                hinfluenzae_txt, hi_st, missing_col=NOT_FOUND, default=NO_DATA)
                    else:
                        d['bmgap2_mlst_st'] = normalize_le_value(row.get('Nm_MLST_ST', ''))
                        d['bmgap2_mlst_cc'] = normalize_le_value(row.get('Nm_MLST_cc', ''))

                    d['FHbp_variant']  = normalize_le_value(row.get('FHbp_protein_subvariant_Novartis', ''))
                    d['FHbp_subfamily'] = normalize_le_value(row.get('FHbp_subfamily', ''))
                    d['FHbp_peptide']  = normalize_le_value(row.get('FHbp_protein_subvariant_Oxford', ''))
                    d['FHbp_Pfizer']   = normalize_le_value(row.get('FHbp_protein_subvariant_Pfizer', ''))
                    d['NadA_variant']  = normalize_le_value(row.get('NadA_Protein_subvariant_Novartis', ''))
                    d['NhbA_peptide']  = normalize_le_value(row.get('NhbA_Protein_subvariant_Novartis', ''))
                    d['PorA_type']     = normalize_le_value(row.get('PorA_type', ''))

                    antigen_measured = any(
                        v != NO_DATA for v in (d['FHbp_variant'], d['FHbp_subfamily'],
                                               d['FHbp_peptide'], d['FHbp_Pfizer'],
                                               d['NadA_variant'], d['NhbA_peptide'],
                                               d['PorA_type']))
                    status['le'] = 'ok' if (antigen_measured
                                            or d['bmgap2_mlst_st'] != NO_DATA) else 'no_calls'

                    if organism == 'hinfluenzae':
                        if antigen_measured or d['bmgap2_mlst_st'] != NO_DATA:
                            d['vaccine_antigens_present'] = 'Not applicable'
                    else:
                        # Gene presence only: no peptide identity, so not a coverage prediction.
                        antigens = (('fHbp', (d['FHbp_variant'], d['FHbp_peptide'], d['FHbp_Pfizer'])),
                                    ('NHBA', (d['NhbA_peptide'],)),
                                    ('NadA', (d['NadA_variant'],)),
                                    ('PorA', (d['PorA_type'],)))
                        detected = []
                        for name, vals in antigens:
                            if not antigen_present(*vals):
                                continue
                            if DISRUPTED in vals:
                                detected.append(f'{name}(disrupted)')
                            elif NEW_ALLELE in vals:
                                detected.append(f'{name}(new)')
                            else:
                                detected.append(name)
                        if detected:
                            d['vaccine_antigens_present'] = ';'.join(detected)
                        elif antigen_measured:
                            d['vaccine_antigens_present'] = NOT_FOUND
            except Exception as e:
                status['le'] = 'failed'
                print(f"Warning: Could not parse LocusExtractor CSV for {sample_id}: {e}", file=sys.stderr)

    if bmscan_species:
        status['bmscan'] = 'ok'
        d['bmgap2_species'] = bmscan_species

    d['bmgap2_status'] = ('ok' if all(v == 'ok' for v in status.values())
                          else ';'.join(f'{k}:{v}' for k, v in status.items() if v != 'ok'))

    return d


# Report headers

HEADER_STANDARD = [
    'sampleID',
    'species_id_qc', 'contamination_flag', 'assembly_qc',
    'mash_species', 'mash_reference', 'mash_distance',
    'kraken2_species', 'kraken2_percent',
    'blast_16s_tophit', 'blast_16s_pident',
    'skani_species', 'skani_ani', 'skani_align_fraction', 'skani_reference',
    'serotype', 'mlst_scheme', 'mlst_st',
    'num_clean_reads', 'avg_read_length', 'avg_read_qual', 'est_coverage',
    'num_contigs', 'longest_contig', 'N50', 'L50', 'total_length', 'gc_content',
    'annotated_cds',
]

HEADER_AMR = ['sampleID', 'carbapenemase_family', 'amr_target',
              'amr_genes', 'amr_subclass']

HEADER_NM = [
    'sampleID',
    'pmga_species', 'bmgap2_species', 'bmgap2_mlst_st', 'bmgap2_mlst_cc', 'serotype_notes', 'nm_genogroup',
    'predicted_resistance', 'resistance_markers',
    'penA_allele', 'penA_mutations', 'penA_phenotype',
    'gyrA_allele', 'gyrA_mutations', 'gyrA_phenotype',
    'parC_allele', 'parC_phenotype',
    'rpoB_allele', 'rpoB_phenotype',
    'ponA_allele', 'ponA_phenotype',
    'blaTEM1_status', 'blaROB1_status',
    'FHbp_variant', 'FHbp_subfamily', 'FHbp_peptide', 'FHbp_Pfizer',
    'NadA_variant', 'NhbA_peptide', 'PorA_type', 'vaccine_antigens_present',
    'bmgap2_status',
]

HEADER_HI = [
    'sampleID',
    'pmga_species', 'bmgap2_species', 'bmgap2_mlst_st', 'bmgap2_mlst_cc', 'serotype_notes', 'hi_capsule_genotype',
    'predicted_resistance', 'resistance_markers',
    'ftsI_allele', 'ftsI_mutations', 'ftsI_phenotype',
    'gyrA_allele', 'gyrA_mutations', 'gyrA_phenotype',
    'parC_allele', 'parC_phenotype',
    'rpoB_allele', 'rpoB_phenotype',
    'folA_allele', 'folA_phenotype',
    'acrR_allele', 'acrR_mutations', 'acrR_phenotype',
    'blaTEM1_status', 'blaROB1_status',
    'bmgap2_status',
]

# MultiQC custom-content emitters

MQC_SPECIES_COLS   = [0, 11, 12, 13, 1, 2, 3]
MQC_SPECIES_HEADER = ['Sample', 'skani_species', 'skani_ani', 'skani_align_fraction',
                      'species_id_qc', 'contamination_flag', 'assembly_qc']

MQC_TYPING_COLS    = [0, 16, 17, 15]
MQC_TYPING_HEADER  = ['Sample', 'mlst_scheme', 'mlst_st', 'serotype']

MQC_AMR_HEADER     = ['Sample', 'carbapenemase_family', 'amr_target',
                      'amr_gene_count', 'amr_genes', 'amr_subclass']


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
    if not rows_std:
        return
    _write_mqc(
        'sanibel_species_mqc.tsv',
        _mqc_preamble(
            'sanibel_species', 'Species ID and QC',
            'skani ANI-confirmed consensus species call with QC verdicts.',
            pconfig={'id': 'sanibel_species_table', 'col1_header': 'Sample',
                     'no_violin': True},
        ),
        MQC_SPECIES_HEADER, MQC_SPECIES_COLS, rows_std,
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
        MQC_TYPING_HEADER, MQC_TYPING_COLS, rows_std,
    )


def emit_sanibel_amr_mqc_table(amr_by_sample):
    if not amr_by_sample:
        return
    rows = []
    for sid, amr in amr_by_sample.items():
        if amr is None:
            rows.append([sid, NO_DATA, NO_DATA, NO_DATA, NO_DATA, NO_DATA])
        elif not amr['genes']:
            rows.append([sid, 'None', 'None', 0, 'None', 'None'])
        else:
            rows.append([sid, carbapenemase_family(amr['genes']),
                         amr_target_genes(amr['genes']), len(amr['genes']),
                         ', '.join(amr['genes']),
                         ', '.join(amr['subclasses']) or 'None'])
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
    parser.add_argument('--outdir',          required=True, help='Pipeline output directory')
    parser.add_argument('--hinfluenzae_txt', default=None,  help='H. influenzae MLST CC lookup table')
    args = parser.parse_args()

    outdir          = args.outdir
    hinfluenzae_txt = args.hinfluenzae_txt

    samples = sorted(
        f.replace('_assembly_stats.txt', '')
        for f in glob.glob('*_assembly_stats.txt')
    )

    if not samples:
        print('summary_report.py: no *_assembly_stats.txt files found.', file=sys.stderr)
        sys.exit(1)

    rows_std = []
    rows_nm  = []
    rows_hi  = []
    amr_by_sample = {}

    for sid in samples:
        sample_dir = os.path.join(outdir, sid)

        amr_by_sample[sid] = parse_amrfinder(f'{sid}_amrfinderplus_report.tsv')

        asm  = parse_assembly_stats(f'{sid}_assembly_stats.txt')
        rm   = parse_read_metrics(f'{sid}_readMetrics.txt')
        cds  = parse_prokka_txt(f'{sid}.txt')
        mlst = parse_mlst(f'{sid}.mlst')
        kr   = parse_kraken_report(f'{sid}.report')
        pmga = parse_pmga(f'{sid}sta.txt')

        _sk_empty = {'ani': NO_DATA, 'confirmed_species': NO_DATA, 'align_fraction': NO_DATA, 'reference': NO_DATA}
        skani_path = f'{sid}_skani.tsv'
        sk = parse_skani(skani_path) if os.path.isfile(skani_path) else _sk_empty
        sk = apply_resolved_species(sk, skani_path, f'{sid}_species_resolved.txt')

        skani_ID_val     = sk['confirmed_species']
        skani_ANI_val    = sk['ani']
        skani_align_val  = sk['align_fraction']
        skani_ID_ref_val = sk['reference']

        skani_genus  = genus_of(skani_ID_val) if skani_ID_val not in (NO_DATA, 'NO ID', 'Inconclusive', '', None) else None
        mash_genus   = asm['genus']
        kraken_genus = kr['species'].split()[0] if kr['species'] != NO_DATA else None

        blast16s_path   = f'{sid}_16s_blast.tsv'
        anchor_16s      = [skani_genus] if skani_genus else [mash_genus, kraken_genus]
        blast16s_result = parse_blast16s_result(blast16s_path, anchor_genera=anchor_16s) \
                          if os.path.isfile(blast16s_path) \
                          else {'pident': NO_DATA, 'tophit': NO_DATA}

        if skani_genus and os.path.isfile(blast16s_path):
            contamination_flag = detect_contamination(
                extract_contam_candidates(blast16s_path), skani_genus)
        else:
            contamination_flag = 'Not screened'

        scheme   = mlst['scheme']
        pmga_sp  = pmga['species'] or NO_DATA
        bmscan_species = read_bmscan_species(sid)
        organism       = meningitis_organism(skani_ID_val, bmscan_species)

        if organism:
            serotype = pmga['prediction'] or NO_DATA
        else:
            # Species-specific serotype from published output dirs
            serotype = NO_DATA
            # Both serotypers run for the E. coli/Shigella complex, so the genus picks the winner
            complex_getters = [
                lambda: get_ecoli_serotype(sample_dir, sid),
                lambda: get_shigella_serotype(sample_dir, sid),
            ]
            if skani_genus == 'Shigella':
                complex_getters.reverse()
            for getter in complex_getters + [
                lambda: get_klebsiella_serotype(sample_dir, sid),
                lambda: get_legionella_serotype(sample_dir, sid),
                lambda: get_salmonella_serotype(sample_dir, sid),
                lambda: get_gas_serotype(sample_dir, sid),
                lambda: get_pneumococcal_serotype(sample_dir, sid),
                lambda: get_acinetobacter_serotype(sample_dir, sid),
                lambda: get_vibrio_serotype(sample_dir, sid),
                lambda: get_pseudomonas_serotype(sample_dir, sid),
                lambda: get_listeria_serotype(sample_dir, sid),
            ]:
                result = getter()
                if result is not None:
                    serotype = result
                    break

        species_id_qc = compute_id_qc(sk, SP_MIN_ANI, SP_MIN_AF, SP_REVIEW_ANI)
        # ANI cannot separate E. coli from Shigella, so an unresolved complex call is not confirmed
        in_complex = skani_ID_val == 'Escherichia_coli' or skani_genus == 'Shigella'
        if species_id_qc == 'PASS' and in_complex                 and not _read_text(f'{sid}_species_resolved.txt'):
            species_id_qc = 'REVIEW (E. coli/Shigella unresolved)'
        contaminated = contamination_flag not in ('None', 'Not screened', NO_DATA)
        assembly_qc = compute_assembly_qc(
            rm['coverage'], asm['num_contigs'], asm['n50'],
            QC_MIN_COVERAGE, QC_WARN_CONTIGS, QC_FAIL_CONTIGS, QC_MIN_N50, contaminated)

        std_row = [
            sid,
            species_id_qc, contamination_flag, assembly_qc,
            f"{asm['genus']}_{asm['species']}", asm['accession'], asm['mash_distance'],
            kr['species'], kr['percent'],
            blast16s_result['tophit'], blast16s_result['pident'],
            skani_ID_val, skani_ANI_val, skani_align_val, skani_ID_ref_val,
            serotype, scheme, mlst['st'],
            rm['num_reads'], rm['avg_read_len'], rm['avg_qual'], rm['coverage'],
            asm['num_contigs'], asm['longest_contig'], asm['n50'], asm['l50'],
            asm['total_length'], asm['gc_content'], cds,
        ]
        rows_std.append(std_row)

        if organism == 'neisseria':
            bm = parse_bmgap2(sid, organism, bmscan_species, hinfluenzae_txt)
            nm_row = [
                sid,
                pmga_sp,
                bm['bmgap2_species'], bm['bmgap2_mlst_st'], bm['bmgap2_mlst_cc'],
                pmga['serotype_notes'] or NO_DATA,
                pmga['prediction'] or NO_DATA,
                bm['predicted_resistance'], bm['resistance_markers'],
                bm['penA_allele'], bm['penA_mutations'], bm['penA_phenotype'],
                bm['gyrA_allele'], bm['gyrA_mutations'], bm['gyrA_phenotype'],
                bm['parC_allele'], bm['parC_phenotype'],
                bm['rpoB_allele'], bm['rpoB_phenotype'],
                bm['ponA_allele'], bm['ponA_phenotype'],
                bm['blaTEM1_status'], bm['blaROB1_status'],
                bm['FHbp_variant'], bm['FHbp_subfamily'], bm['FHbp_peptide'], bm['FHbp_Pfizer'],
                bm['NadA_variant'], bm['NhbA_peptide'], bm['PorA_type'], bm['vaccine_antigens_present'],
                bm['bmgap2_status'],
            ]
            rows_nm.append(nm_row)

        elif organism == 'hinfluenzae':
            bm = parse_bmgap2(sid, organism, bmscan_species, hinfluenzae_txt)
            hi_row = [
                sid,
                pmga_sp,
                bm['bmgap2_species'], bm['bmgap2_mlst_st'], bm['bmgap2_mlst_cc'],
                pmga['serotype_notes'] or NO_DATA,
                pmga['prediction'] or NO_DATA,
                bm['predicted_resistance'], bm['resistance_markers'],
                bm['penA_allele'], bm['penA_mutations'], bm['penA_phenotype'],
                bm['gyrA_allele'], bm['gyrA_mutations'], bm['gyrA_phenotype'],
                bm['parC_allele'], bm['parC_phenotype'],
                bm['rpoB_allele'], bm['rpoB_phenotype'],
                bm['folA_allele'], bm['folA_phenotype'],
                bm['acrR_allele'], bm['acrR_mutations'], bm['acrR_phenotype'],
                bm['blaTEM1_status'], bm['blaROB1_status'],
                bm['bmgap2_status'],
            ]
            rows_hi.append(hi_row)

    def write_report(path, header, rows):
        with open(path, 'w', encoding='utf-8-sig') as fh:
            fh.write('\t'.join(header) + '\n')
            for row in sorted(rows, key=lambda r: r[0]):
                fh.write('\t'.join(str(v).replace(',', ';') for v in row) + '\n')
        print(f"summary_report.py: wrote {path} ({len(rows)} sample(s))")

    rows_amr = []
    for sid, amr in amr_by_sample.items():
        if amr is None:
            rows_amr.append([sid, NO_DATA, NO_DATA, NO_DATA, NO_DATA])
        elif not amr['genes']:
            rows_amr.append([sid, 'None', 'None', 'None', 'None'])
        else:
            rows_amr.append([sid, carbapenemase_family(amr['genes']),
                             amr_target_genes(amr['genes']),
                             ', '.join(amr['genes']),
                             ', '.join(amr['subclasses']) or 'None'])

    if rows_std:
        write_report('sum_report.txt',    HEADER_STANDARD, rows_std)
    if rows_nm:
        write_report('nm_sum_report.txt', HEADER_NM,       rows_nm)
    if rows_hi:
        write_report('hi_sum_report.txt', HEADER_HI,       rows_hi)
    if rows_amr:
        write_report('amr_report.txt',    HEADER_AMR,      rows_amr)

    emit_sanibel_mqc_tables(rows_std)
    emit_sanibel_amr_mqc_table(amr_by_sample)

    if not (rows_std or rows_nm or rows_hi):
        print('summary_report.py: no rows generated.', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
