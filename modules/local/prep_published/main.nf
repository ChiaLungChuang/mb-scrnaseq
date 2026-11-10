// Stream the published GEO counts CSV into a sparse matrix, once per run.
process PREP_PUBLISHED {
    tag "${counts.name}"
    label 'process_r'
    container "${params.r_container}"

    input:
    path counts
    path lib
    val fingerprint   // SHA-256 of lib/R and bin/*.R: re-run on -resume when the code changes

    output:
    path "published_counts.rds", emit: rds
    path "published_counts.summary.json", emit: summary

    script:
    """
    # code fingerprint: ${fingerprint}
    # Threads limited to the task's CPUs (data.table, OpenMP, BLAS); each R script logs them.
    export OMP_NUM_THREADS=${task.cpus} OMP_THREAD_LIMIT=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus} R_DATATABLE_NUM_THREADS=${task.cpus}
    prep_published.R --lib ${lib} --counts ${counts} \\
        --out-prefix published_counts --chunk-rows ${params.published_chunk_rows}
    """

    stub:
    """
    touch published_counts.rds
    echo '{}' > published_counts.summary.json
    """
}
