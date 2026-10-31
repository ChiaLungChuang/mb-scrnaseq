// One table and one verdict across all re-quantified samples. The verdict
// records whether the criteria were pre-registered and which criteria and
// code produced it.
process SUMMARISE_CONCORDANCE {
    label 'process_low'
    container "${params.r_container}"

    input:
    path jsons
    path lib
    val preregistered
    val criteria_sha256
    val code_sha256

    output:
    path "concordance_summary.tsv", emit: tsv
    path "concordance_summary.verdict.json", emit: verdict

    script:
    """
    summarise_concordance.R --lib ${lib} --out-prefix concordance_summary \\
        --preregistered ${preregistered} --criteria-sha256 ${criteria_sha256} --code-sha256 ${code_sha256} \\
        ${jsons}
    """

    stub:
    """
    touch concordance_summary.tsv
    echo '{}' > concordance_summary.verdict.json
    """
}
