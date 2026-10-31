// Compare one sample's re-quantified counts with the published counts and
// judge the result against the pre-registered criteria.
//
// `fingerprint` (SHA-256 of lib/R, bin/*.R and the criteria file) is an input
// so that -resume re-runs this step whenever the code or the criteria
// change; Nextflow does not hash the contents of those files itself.
process CONCORDANCE {
    tag "${meta.id}"
    label 'process_r'
    container "${params.r_container}"

    input:
    tuple val(meta), path(raw), path(filtered)
    path published
    path criteria
    path lib
    val fingerprint

    output:
    tuple val(meta), path("${meta.id}.concordance.json"), emit: json
    tuple val(meta), path("${meta.id}.per_cell.tsv.gz"), emit: per_cell

    script:
    """
    # inputs fingerprint: ${fingerprint}
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
