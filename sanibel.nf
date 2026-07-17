#!/usr/bin/env nextflow

/*
  Sanibel Pipeline (named after Sanibel Island in southwest Florida)
  Florida's BPHL Nextflow pipeline for Bacterial WGS Analysis
  Authors: Sarah Schemedes, Yibo Dong, Arnold Rodriguez-Hilario, Molly Mitchell
  Email: bphl-sebioinformatics@flhealth.gov 
*/

include { fastqc }                from './modules/fastqc.nf'
include { trimmomatic }           from './modules/trimmomatic.nf'
include { bbtools_adapters }      from './modules/bbtools.nf'
include { bbtools_phix }          from './modules/bbtools.nf'
include { fastqc2 }               from './modules/fastqc2.nf'
include { multiqc }               from './modules/multiqc.nf'
include { multiqc_global }        from './modules/multiqc.nf'
include { mash }                  from './modules/mash.nf'
include { unicycler }             from './modules/unicycler.nf'
include { kraken }                from './modules/kraken.nf'
include { quast }                 from './modules/quast.nf'
include { parse_assembly }        from './modules/parse_assembly.nf'
include { readssum }              from './modules/readssum.nf'
include { prokka }                from './modules/prokka.nf'
include { amrfinder }             from './modules/amrfinder.nf'
include { mlst }                  from './modules/mlst.nf'
include { pmga }                  from './modules/pmga.nf'
include { download_16s_db }       from './modules/blast_16s.nf'
include { blast_16s }             from './modules/blast_16s.nf'
include { download_mlst_tables }  from './modules/download_mlst_tables.nf'
include { aggregate_species_id }  from './modules/aggregate_species_id.nf'
include { build_candidates }      from './modules/build_candidates.nf'
include { download_one_candidate } from './modules/download_one_candidate.nf'
include { collect_references }    from './modules/collect_references.nf'
include { skani }                 from './modules/skani.nf'
include { bmgap2_amr }            from './modules/bmgap2_amr.nf'
include { bmgap2_locusextractor } from './modules/bmgap2_locusextractor.nf'
include { bmgap2_bmscan }         from './modules/bmgap2_bmscan.nf'
include { legsta }                from './modules/legsta.nf'
include { kleborate }             from './modules/kleborate.nf'
include { shigatyper }            from './modules/shigatyper.nf'
include { emm_typing }            from './modules/emm_typing.nf'
include { seqsero2 }              from './modules/seqsero2.nf'
include { serotypefinder }        from './modules/serotypefinder.nf'
include { plasmidfinder }         from './modules/plasmidfinder.nf'
include { seroba }                from './modules/seroba.nf'
include { pasty }                 from './modules/pasty.nf'
include { kaptive as kaptive_ab } from './modules/kaptive.nf'
include { kaptive as kaptive_vp } from './modules/kaptive.nf'
include { lissero }               from './modules/lissero.nf'
include { serotype }              from './modules/serotype.nf'
include { summary_report }        from './modules/summary_report.nf'

def rebind(ch, metaCh) {
    ch.map  { meta, x -> [ meta.id, x ] }
      .join(metaCh)
      .map  { _id, x, emeta -> [ emeta, x ] }
}

workflow {
    log.info """
    Sanibel — Bacterial WGS Analysis Pipeline
    ==========================================================================
    input dir   : ${params.input}
    output dir  : ${params.output}
    bmgap2 db   : ${params.bmgap2_db}
    kraken db   : ${params.kraken_db}
    ==========================================================================
    """

    // FASTQ Input channel
    ch_reads = channel.fromFilePairs(
        ["${params.input}/*_{1,2}.fastq.gz",
         "${params.input}/*_R{1,2}_*.fastq.gz"],
        checkIfExists: false
    )
    .map { id, files ->
        def clean_id = id.replaceAll(/_S\d+_L\d+$/, '')
        def meta = [ id: clean_id, single_end: false ]
        [ meta, files ]
    }

    // MLST CC reference tables (downloaded once, cached via storeDir)
    ch_mlst_tables     = download_mlst_tables()
    ch_neisseria_txt   = ch_mlst_tables.neisseria.first()
    ch_hinfluenzae_txt = ch_mlst_tables.hinfluenzae.first()
    ch_mlst_schemes    = channel.value(file("${projectDir}/assets/mlst_schemes.tsv", checkIfExists: true))

    // Read as a map, not a staged file
    kleborate_presets = file("${projectDir}/assets/kleborate_presets.tsv", checkIfExists: true)
        .readLines()
        .findAll { line -> line.trim() && !line.startsWith('#') }
        .collectEntries { line -> def parts = line.split('\t'); [(parts[0]): parts[1]] }

    // QC & read preprocessing
    ch_fastqc  = fastqc(ch_reads)
    ch_trimmed = trimmomatic(ch_reads)
    ch_adaptertrim = bbtools_adapters(ch_trimmed.reads)
    ch_clean       = bbtools_phix(ch_adaptertrim.reads)
    ch_fastqc2 = fastqc2(ch_clean.reads)

    // Per-sample MultiQC (raw and clean FastQC reports)
    multiqc(
        ch_fastqc.report
            .join(ch_fastqc2.report, by: 0)
            .map { meta, raw_reports, clean_reports ->
                [meta, (raw_reports instanceof List ? raw_reports : [raw_reports]) +
                       (clean_reports instanceof List ? clean_reports : [clean_reports])]
            }
    )

    // Species ID, assembly and read classification
    ch_mash     = mash(ch_clean.reads)
    ch_assembly = unicycler(ch_clean.reads)
    ch_kraken   = kraken(ch_clean.reads)
    ch_quast    = quast(ch_assembly.assembly)

    // Mash top-hit (genus/species/distance/accession)
    ch_tophit = parse_assembly(ch_mash.distances)

    // Meta keyed by sample ID, enriched with the mash species for prokka's fallback
    ch_meta_by_id = ch_tophit
        .splitCsv(header: true, sep: '\t', elem: 1)
        .map { meta, row -> [ meta.id, meta + [ mash_species: row.genus + '_' + row.species ] ] }

    // Rebind clean reads and assembly with enriched meta
    ch_clean_enriched    = rebind(ch_clean.reads,       ch_meta_by_id)
    ch_assembly_enriched = rebind(ch_assembly.assembly, ch_meta_by_id)

    // Read metrics (genome size computed from the assembly inside the process)
    ch_readssum = readssum(
        ch_clean_enriched.map { meta, reads -> [ meta.id, meta, reads ] }
            .join(ch_assembly_enriched.map { meta, asm -> [ meta.id, asm ] })
            .map { _id, meta, reads, asm -> [ meta, reads, asm ] }
    )

    ch_amrfinder = amrfinder(ch_assembly_enriched)

    // Kraken output with enriched meta
    ch_kraken_enriched = rebind(ch_kraken.out, ch_meta_by_id)

    // 16S BLAST (DB downloaded once, cached via storeDir)
    ch_16s_db    = download_16s_db().db
    ch_blast_16s = blast_16s(ch_assembly_enriched, ch_16s_db)

    // Kraken and 16S keyed by sample id (shared by the vote and the candidate pool)
    ch_kraken_by_id = ch_kraken.out.map       { meta, r -> [ meta.id, r ] }
    ch_blast_by_id  = ch_blast_16s.result.map { meta, r -> [ meta.id, r ] }

    // 2-of-3 species vote (Mash + Kraken2 + 16S BLAST)
    ch_aggregate = aggregate_species_id(
        ch_tophit.map { meta, t -> [ meta.id, t ] }
            .join(ch_meta_by_id)
            .join(ch_kraken_by_id)
            .join(ch_blast_by_id)
            .map { _id, tophit, emeta, kreport, blast -> [ emeta, tophit, kreport, blast ] }
    )

    // Build candidate species pool from all tools
    ch_pool = build_candidates(
        ch_mash.distances
            .map  { meta, d -> [ meta.id, d ] }
            .join(ch_kraken_by_id)
            .join(ch_blast_by_id)
            .join(ch_meta_by_id)
            .map  { _id, distances, kreport, blast, emeta -> [ emeta, distances, kreport, blast ] }
    )

    // Download reference genomes: one task per candidate, then collect per sample
    ch_candidates = ch_pool.pool
        .splitCsv(header: true, sep: '\t', elem: 1)
        .map { meta, row -> [ meta, row.species, row.accession ] }
    ch_refs = collect_references(
        download_one_candidate(ch_candidates).fna
            .map    { meta, fnas -> [ meta.id, meta, fnas ] }
            .groupTuple(by: 0)
            .map    { _id, metas, fnas -> [ metas[0], fnas.flatten() ] }
    )

    // Multi-reference ANI confirmation with skani
    ch_skani = skani(
        ch_assembly_enriched
            .map  { meta, asm -> [ meta.id, asm ] }
            .join(ch_refs.references.map { meta, d -> [ meta.id, d ] })
            .join(ch_meta_by_id)
            .map  { _id, asm, refs_dir, emeta -> [ emeta, asm, refs_dir ] }
    )

    // skani-confirmed species drives the species-specific analyses
    ch_meta_typed = ch_meta_by_id
        .join(ch_skani.species.map { meta, f -> [ meta.id, f.text.trim() ] }, remainder: true)
        .map { id, meta, sp ->
            def species = sp ?: 'Unknown'
            def genus   = sp ? sp.tokenize('_')[0] : 'Unknown'
            [ id, meta + [ species: species, genus: genus ] ]
        }

    // Rebind the channels the typing modules use
    ch_assembly_typed = rebind(ch_assembly_enriched, ch_meta_typed)
    ch_clean_typed    = rebind(ch_clean_enriched,    ch_meta_typed)

    // Annotation and MLST
    ch_prokka     = prokka(ch_assembly_typed)
    ch_mlst       = mlst(ch_assembly_typed, ch_mlst_schemes)

    // PMGA + BMGAP2
    ch_pmga = pmga(
        ch_assembly_typed
            .filter { meta, _a -> meta.genus in ['Neisseria', 'Haemophilus'] }
            .join(ch_mlst.out, by: 0)
    )
    ch_bmgap2_amr = bmgap2_amr(ch_mlst.out.join(ch_pmga.out, by: 0))
    ch_bmgap2_le  = bmgap2_locusextractor(ch_bmgap2_amr.out)
    ch_bmgap2_bmscan = bmgap2_bmscan(ch_bmgap2_le.out)

    // Species-specific analyses
    ch_legsta = legsta(ch_assembly_typed.filter { meta, _a -> meta.species == 'Legionella_pneumophila' })
    ch_kleborate = kleborate(
        ch_assembly_typed
            .map    { meta, a  -> [ meta + [kleborate_preset: kleborate_presets[meta.species?.tokenize('_')?.take(2)?.join('_')]], a ] }
            .filter { meta, _a -> meta.kleborate_preset }
    )
    ch_shigatyper     = shigatyper(ch_clean_typed.filter     { meta, _r -> meta.genus   == 'Shigella' })
    ch_emm_typing     = emm_typing(ch_clean_typed.filter     { meta, _r -> meta.species in ['Streptococcus_pyogenes', 'Streptococcus_dysgalactiae'] })
    ch_seqsero2       = seqsero2(ch_clean_typed.filter       { meta, _r -> meta.genus   == 'Salmonella' })
    ch_serotypefinder = serotypefinder(ch_clean_typed.filter { meta, _r -> meta.species == 'Escherichia_coli' })
    ch_plasmidfinder  = plasmidfinder(ch_clean_enriched)
    ch_seroba         = seroba(ch_clean_typed.filter         { meta, _r -> meta.species == 'Streptococcus_pneumoniae' })
    ch_pasty          = pasty(ch_assembly_typed.filter       { meta, _a -> meta.species == 'Pseudomonas_aeruginosa' })
    ch_kaptive_ab     = kaptive_ab(ch_assembly_typed.filter  { meta, _a -> meta.species == 'Acinetobacter_baumannii' }, 'ab')
    ch_kaptive_vp     = kaptive_vp(ch_assembly_typed.filter  { meta, _a -> meta.species == 'Vibrio_parahaemolyticus' }, 'vp')
    ch_lissero        = lissero(ch_assembly_typed.filter     { meta, _a -> meta.species == 'Listeria_monocytogenes' })

    // Normalize each sample's serotype to a single value (one typing tool per sample)
    ch_serotype = serotype(
        ch_legsta.results.map              { meta, d -> [ meta, 'legsta',         d ] }
            .mix(ch_kleborate.results.map      { meta, d -> [ meta, 'kleborate',      d ] },
                 ch_shigatyper.results.map     { meta, d -> [ meta, 'shigatyper',     d ] },
                 ch_emm_typing.results.map     { meta, d -> [ meta, 'emm_typing',     d ] },
                 ch_seqsero2.results.map       { meta, d -> [ meta, 'seqsero2',       d ] },
                 ch_serotypefinder.results.map { meta, d -> [ meta, 'serotypefinder', d ] },
                 ch_seroba.results.map         { meta, d -> [ meta, 'seroba',         d ] },
                 ch_pasty.results.map          { meta, d -> [ meta, 'pasty',          d ] },
                 ch_kaptive_ab.results.map     { meta, d -> [ meta, 'kaptive_ab',     d ] },
                 ch_kaptive_vp.results.map     { meta, d -> [ meta, 'kaptive_vp',     d ] },
                 ch_lissero.results.map        { meta, d -> [ meta, 'lissero',        d ] })
    )

    // Barrier only for the bmgap2 side-channel (still writes into the output dir)
    ch_bmgap2_barrier = ch_bmgap2_bmscan
        .map { _meta, _f -> 1 }
        .collect()
        .map { _ids -> true }

    ch_summary = summary_report(
        ch_bmgap2_barrier,
        ch_quast.report.map         { _meta, q    -> q    }.collect(),
        ch_tophit.map               { _meta, t    -> t    }.collect(),
        ch_readssum.out.map         { _meta, rm   -> rm   }.collect(),
        ch_prokka.cds_txt.map       { _meta, ptxt -> ptxt }.collect(),
        ch_mlst.out.map             { _meta, mlst_file -> mlst_file }.collect(),
        ch_kraken_enriched.map      { _meta, kr   -> kr   }.collect(),
        ch_pmga.out.map             { _meta, pmga_file -> pmga_file }.collect().ifEmpty([]),
        ch_neisseria_txt,
        ch_hinfluenzae_txt,
        ch_aggregate.out.map        { _meta, f -> f }.collect().ifEmpty([]),
        ch_skani.result.map         { _meta, f -> f }.collect().ifEmpty([]),
        ch_blast_16s.result.map     { _meta, f -> f }.collect().ifEmpty([]),
        ch_amrfinder.out.map        { _meta, f -> f }.collect().ifEmpty([]),
        ch_serotype.out.map         { _meta, f -> f }.collect().ifEmpty([])
    )

    // Run-level interactive MultiQC across all samples
    multiqc_global(
        ch_summary.summary,
        channel.value(file("${projectDir}/assets/multiqc_config.yaml",         checkIfExists: true)),
        channel.value(file("${projectDir}/assets/sanibel_pipeline_logo_v2.png", checkIfExists: true)),
        channel.value(file("${projectDir}/assets/sanibel_report.css",           checkIfExists: true)),
        channel.value(file("${projectDir}/nextflow.config",                     checkIfExists: true))
    )
}
