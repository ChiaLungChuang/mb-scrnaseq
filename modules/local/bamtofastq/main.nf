// Recover the original FASTQs from a Cell Ranger BAM (10x bamtofastq).
// The output folder must not exist beforehand; bamtofastq creates it.
process BAMTOFASTQ {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/10x_bamtofastq:1.4.1--h3ab6199_4'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("fastq"), emit: fastq
    path "versions.yml", emit: versions

    script:
    """
    bamtofastq --nthreads ${task.cpus} ${task.ext.args ?: ''} ${bam} fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamtofastq: \$(bamtofastq --help | head -n 1 | sed 's/^bamtofastq v//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p fastq/${meta.id}_0_1_STUB
    touch fastq/${meta.id}_0_1_STUB/bamtofastq_S1_L001_R1_001.fastq.gz
    touch fastq/${meta.id}_0_1_STUB/bamtofastq_S1_L001_R2_001.fastq.gz
    printf '"%s":\\n    bamtofastq: stub\\n' "${task.process}" > versions.yml
    """
}
