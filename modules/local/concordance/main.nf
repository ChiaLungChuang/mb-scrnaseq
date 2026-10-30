// Compare one sample's re-quantified counts with the published counts and
// judge the result against the pre-registered criteria.
process CONCORDANCE {
    tag "${meta.id}"
    label 'process_r'
    container "${params.r_container}"

    input:
    tuple val(meta), path(raw), path(filtered)
    path published
    path criteria
    path lib

    output:
    tuple val(meta), path("${meta.id}.concordance.json"), emit: json
    tuple val(meta), path("${meta.id}.per_cell.tsv.gz"), emit: per_cell

    script:
    """
    concordance.R --lib ${lib} --sample ${meta.id} --label '${meta.published_label}' \\
        --raw-dir ${raw} --filtered-dir ${filtered} --published ${published} \\
        --criteria ${criteria} --out-prefix ${meta.id}
    """

    stub:
    """
    echo '{}' > ${meta.id}.concordance.json
    echo stub | gzip > ${meta.id}.per_cell.tsv.gz
    """
}
