/*
 * Prepares the raw or compressed data holding spectra information for the subsequent database search.
 */

//
// MODULE: Loaded from modules/local/
//

include { MS2RESCORE                                                  } from '../../../modules/local/ms2rescore'
include { OPENMS_PSMFEATUREEXTRACTOR                                  } from '../../../modules/local/openms/psmfeatureextractor'
include {
    OPENMS_PERCOLATORADAPTER ;
    OPENMS_PERCOLATORADAPTER as OPENMS_PERCOLATORADAPTER_GLOBAL
} from '../../../modules/local/openmsthirdparty/percolatoradapter'
include { OPENMS_TEXTEXPORTER as OPENMS_TEXTEXPORTER_GLOBAL           } from '../../../modules/nf-core/openms/textexporter/main'
include { PYOPENMS_FLAGINSOURCE                                       } from '../../../modules/local/pyopenms/flaginsource'
//
// MODULE: Installed directly from nf-core/modules
//

include { OPENMS_IDMERGER as OPENMS_IDMERGER_GLOBAL                   } from '../../../modules/nf-core/openms/idmerger/main'
include { OPENMS_IDSCORESWITCHER                                      } from '../../../modules/nf-core/openms/idscoreswitcher/main.nf'
include {
    OPENMS_IDFILTER as OPENMS_IDFILTER_Q_VALUE ;
    OPENMS_IDFILTER as OPENMS_IDFILTER_Q_VALUE_GLOBAL ;
    OPENMS_IDFILTER as OPENMS_IDFILTER_GLOBAL
} from '../../../modules/nf-core/openms/idfilter/main'

workflow RESCORE {
    take:
    ch_merged_runs
    ch_multiqc_files

    main:
    // Compute features via ms2rescore
    MS2RESCORE(ch_merged_runs)

    if (params.rescoring_engine == 'mokapot') {
        log.warn("The rescoring engine is set to mokapot. This rescoring engine currently only supports psm-level-fdr via ms2rescore.")
        if (params.global_fdr) {
            log.warn("Global FDR is currently not supported by mokapot. The global_fdr parameter will be ignored.")
        }
        // Switch comet e-value to mokapot q-value
        OPENMS_IDSCORESWITCHER(MS2RESCORE.out.idxml)
        ch_rescored_runs = OPENMS_IDSCORESWITCHER.out.idxml

        // Filter by mokapot q-value
        OPENMS_IDFILTER_Q_VALUE(ch_rescored_runs.map { group_meta, idxml -> [group_meta, idxml, []] })
        ch_filter_q_value = OPENMS_IDFILTER_Q_VALUE.out.filtered
    }
    else {
        // Extract PSM features for Percolator
        OPENMS_PSMFEATUREEXTRACTOR(MS2RESCORE.out.idxml.join(MS2RESCORE.out.feature_names))

        // Run Percolator with local FDR
        OPENMS_PERCOLATORADAPTER(OPENMS_PSMFEATUREEXTRACTOR.out.idxml)
        ch_multiqc_files = ch_multiqc_files.mix(OPENMS_PERCOLATORADAPTER.out.feature_weights.map { meta, feature_weights -> feature_weights })
        ch_pout = OPENMS_PERCOLATORADAPTER.out.idxml

        if (params.global_fdr) {
            // Group by search_preset for global FDR. Samples without a preset all share
            // the same params (CLI or defaults), so they correctly group under 'global'.
            OPENMS_IDMERGER_GLOBAL(
                OPENMS_PSMFEATUREEXTRACTOR.out.idxml.map { group_meta, idxml -> [group_meta + [id: group_meta.search_preset ?: 'global'], idxml] }.groupTuple()
            )
            // Run Percolator with global FDR (one per preset group)
            OPENMS_PERCOLATORADAPTER_GLOBAL(OPENMS_IDMERGER_GLOBAL.out.idxml)
            ch_rescored_runs = OPENMS_PERCOLATORADAPTER_GLOBAL.out.idxml
            // Filter by global percolator q-value
            OPENMS_IDFILTER_Q_VALUE_GLOBAL(ch_rescored_runs.map { id, idxml -> [id, idxml, []] })
            // Backfilter: match each local file with its corresponding preset's global FDR file
            OPENMS_IDFILTER_GLOBAL(
                ch_pout.map { group_meta, idxml ->
                    [group_meta.search_preset ?: 'global', group_meta, idxml]
                }.combine(
                    OPENMS_IDFILTER_Q_VALUE_GLOBAL.out.filtered.map { global_meta, idxml -> [global_meta.id, idxml] },
                    by: 0
                ).map { preset, group_meta, local_idxml, global_filtered_idxml ->
                    [group_meta, local_idxml, global_filtered_idxml]
                }
            )
            ch_filter_q_value = OPENMS_IDFILTER_GLOBAL.out.filtered
            // Save globally merged runs in tsv (one per preset group)
            OPENMS_TEXTEXPORTER_GLOBAL(OPENMS_IDFILTER_Q_VALUE_GLOBAL.out.filtered)
        }
        else {
            ch_rescored_runs = ch_pout
            // Filter by percolator q-value
            OPENMS_IDFILTER_Q_VALUE(ch_rescored_runs.map { group_meta, idxml -> [group_meta, idxml, []] })
            ch_filter_q_value = OPENMS_IDFILTER_Q_VALUE.out.filtered
        }
    }

    // Branch the FDR-filtered runs FIRST so empty samples (no peptides passing FDR) are diverted
    // before any flagging, regardless of --flag_in_source.
    ch_filter_q_value
        .map { meta, file -> [[id: meta.id], file] }
        .branch {
            // Empty FDR-filtered idXML (no peptides) is ~120 lines of OpenMS scaffolding.
            non_empty: it[1].countLines() > 130
            empty:     true
        }
        .set { ch_fdr_branched }

    // Optionally flag in-source fragments and re-inject the flagged candidates into the
    // FDR-filtered set so they reach quantification (parents = the FDR-filtered hits, candidates
    // = the 100% FDR rescored runs). Only non-empty samples are flagged; the augmented filtered
    // replaces fdr_filtered and the is_isf-annotated rescored runs replace rescored_runs so the
    // flag survives into quantification.
    if (params.flag_in_source && params.global_fdr) {
        log.warn("'--flag_in_source' is not supported together with '--global_fdr' and will be skipped.")
    }

    if (params.flag_in_source && !params.global_fdr) {
        PYOPENMS_FLAGINSOURCE(
            ch_fdr_branched.non_empty
                .map { meta, file -> [meta.id, meta, file] }
                .join(ch_rescored_runs.map { meta, file -> [meta.id, file] }, by: 0)
                .map { id, meta, filtered, rescored -> [meta, filtered, rescored] }
        )
        // Re-branch the augmented filtered (defensive; appended ISF keep it non-empty).
        PYOPENMS_FLAGINSOURCE.out.idxml
            .map { meta, file -> [[id: meta.id], file] }
            .branch {
                non_empty: it[1].countLines() > 130
                empty:     true
            }
            .set { ch_fdr_flagged }

        ch_fdr_filtered       = ch_fdr_flagged.non_empty
        ch_fdr_filtered_empty = ch_fdr_branched.empty.mix(ch_fdr_flagged.empty)
        ch_rescored_runs      = PYOPENMS_FLAGINSOURCE.out.rescored_idxml
    }
    else {
        ch_fdr_filtered       = ch_fdr_branched.non_empty
        ch_fdr_filtered_empty = ch_fdr_branched.empty
    }

    emit:
    rescored_runs      = ch_rescored_runs.map { meta, file -> [[id: meta.id], file] }
    fdr_filtered       = ch_fdr_filtered
    fdr_filtered_empty = ch_fdr_filtered_empty
    multiqc_files      = ch_multiqc_files
}
