// One table and one verdict across all re-quantified samples.
process SUMMARISE_CONCORDANCE {
    label 'process_low'
    container "${params.r_container}"

    input:
    path jsons
    path lib

    output:
    path "concordance_summary.tsv", emit: tsv
    path "concordance_summary.verdict.json", emit: verdict

    script:
    """
    summarise_concordance.R --lib ${lib} --out-prefix concordance_summary ${jsons}
    """

    stub:
    """
    touch concordance_summary.tsv
    echo '{}' > concordance_summary.verdict.json
    """
}
