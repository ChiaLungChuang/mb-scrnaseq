// Per-sample QC flags and recorded thresholds on STARsolo's called cells.
process SAMPLE_QC {
    tag "${meta.id}"
    label 'process_low'   // 0.8 GB for 5,000 cells, 1.9 GB for 10,000 dense ones (tests/scale/README.md)
    container "${params.r_container}"

    input:
    tuple val(meta), path(filtered)
    path lib
    val fingerprint   // SHA-256 of lib/R and bin/*.R: re-run on -resume when the code changes

    output:
    tuple val(meta), path("${meta.id}.qc_cells.tsv.gz"), emit: cells
    tuple val(meta), path("${meta.id}.qc_thresholds.json"), emit: thresholds
    path "versions.yml", emit: versions

    script:
    """
    # code fingerprint: ${fingerprint}
    # Threads limited to the task's CPUs (data.table, OpenMP, BLAS); each R script logs them.
    export OMP_NUM_THREADS=${task.cpus} OMP_THREAD_LIMIT=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus} R_DATATABLE_NUM_THREADS=${task.cpus}
    sample_qc.R --lib ${lib} --sample ${meta.id} --filtered-dir ${filtered} \\
        --nmads ${params.qc_nmads} --out-prefix ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(Rscript -e 'cat(as.character(getRversion()))')
        scuttle: \$(Rscript -e 'cat(as.character(packageVersion("scuttle")))')
    END_VERSIONS
    """

    stub:
    """
    echo stub | gzip > ${meta.id}.qc_cells.tsv.gz
    echo '{}' > ${meta.id}.qc_thresholds.json
    printf '"%s":\\n    r: stub\\n' "${task.process}" > versions.yml
    """
}
