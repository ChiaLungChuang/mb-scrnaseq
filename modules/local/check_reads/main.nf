// Pair barcode and cDNA reads and confirm the declared chemistry from the
// actual read lengths. Fails the sample on any mismatch (see bin/check_reads.sh).
process CHECK_READS {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/10x_bamtofastq:1.4.1--h3ab6199_4'

    input:
    tuple val(meta), path(fastq)
    val fingerprint   // SHA-256 of bin/check_reads.sh: re-run on -resume when the code changes

    output:
    tuple val(meta), path(fastq), path("${meta.id}.reads.tsv"), emit: reads
    path "${meta.id}.check_reads.log", emit: log

    script:
    """
    # code fingerprint: ${fingerprint}
    check_reads.sh ${fastq} ${meta.chemistry} ${meta.id}.reads.tsv | tee ${meta.id}.check_reads.log
    """

    stub:
    """
    printf '%s\\t%s\\n' ${fastq}/x_R1_001.fastq.gz ${fastq}/x_R2_001.fastq.gz > ${meta.id}.reads.tsv
    echo stub > ${meta.id}.check_reads.log
    """
}
