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
    # Threads limited to the task's CPUs (data.table, OpenMP, BLAS); each R script logs them.
    export OMP_NUM_THREADS=${task.cpus} OMP_THREAD_LIMIT=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus} R_DATATABLE_NUM_THREADS=${task.cpus}
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
