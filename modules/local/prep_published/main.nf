// Stream the published GEO counts CSV into a sparse matrix, once per run.
process PREP_PUBLISHED {
    tag "${counts.name}"
    label 'process_r_high'
    container "${params.r_container}"

    input:
    path counts
    path lib

    output:
    path "published_counts.rds", emit: rds
    path "published_counts.summary.json", emit: summary

    script:
    """
    prep_published.R --lib ${lib} --counts ${counts} \\
        --out-prefix published_counts --chunk-rows ${params.published_chunk_rows}
    """

    stub:
    """
    touch published_counts.rds
    echo '{}' > published_counts.summary.json
    """
}
