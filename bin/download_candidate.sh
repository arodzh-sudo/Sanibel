#!/usr/bin/env bash
# Download up to N RefSeq genomes for one candidate species into the current dir,
# named <safe_name>__<accession>.fna. Soft-fails to no files if nothing downloads.
set -eu

species="$1"
accession="$2"
safe_name="$3"
n_refs="$4"
ref_dir="$PWD"

extract_all() {
    unz_dir="$1"
    for accdir in "${unz_dir}/ncbi_dataset/data"/*/; do
        [ -d "${accdir}" ] || continue
        acc=$(basename "${accdir}")
        fna=$(find "${accdir}" -name "*.fna" | head -1)
        [ -n "${fna}" ] && cp "${fna}" "${ref_dir}/${safe_name}__${acc}.fna"
    done
}

# RefSeq accessions for this species, plus the sketch representative
acc_list=$(datasets summary genome taxon "${species}" \
    --assembly-source refseq --limit "${n_refs}" \
    --report ids_only --as-json-lines 2>/dev/null \
    | grep -oE 'GC[FA]_[0-9]+\.[0-9]+' || true)
if [ "${accession}" != "NA" ] && [ -n "${accession}" ]; then
    acc_list=$(printf '%s\n%s\n' "${acc_list}" "${accession}")
fi
acc_list=$(printf '%s\n' "${acc_list}" | grep -oE 'GC[FA]_[0-9]+\.[0-9]+' | sort -u || true)

if [ -n "${acc_list}" ]; then
    datasets download genome accession ${acc_list} --include genome --filename dl.zip 2>/dev/null \
        && unzip -o dl.zip -d dl_dir 2>/dev/null \
        && extract_all dl_dir \
        || true
    rm -rf dl.zip dl_dir 2>/dev/null || true
fi

# Fallback: the species reference genome
if ! ls "${ref_dir}/${safe_name}__"*.fna >/dev/null 2>&1; then
    datasets download genome taxon "${species}" --reference \
        --assembly-source refseq --include genome --filename ref.zip 2>/dev/null \
        && unzip -o ref.zip -d ref_dir 2>/dev/null \
        && extract_all ref_dir \
        || true
    rm -rf ref.zip ref_dir 2>/dev/null || true
fi
