// Build the STAR genome index. Human needs ~32 GB RAM; build it once and pass
// it back in with --star_index on later runs.
process STAR_GENOMEGENERATE {
    tag "${fasta.baseName}"
    label 'process_high'
    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_8'

    input:
    path fasta
    path gtf

    output:
    path "star_index", emit: index
    path "versions.yml", emit: versions

    script:
    def ram = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes()}" : ''
    """
    mkdir star_index
    STAR --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir star_index \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang ${params.star_sjdb_overhang} \\
        --genomeSAindexNbases ${params.star_genome_sa_index_nbases} \\
        ${ram}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version)
    END_VERSIONS
    """

    stub:
    """
    mkdir star_index
    touch star_index/SA
    printf '"%s":\\n    star: stub\\n' "${task.process}" > versions.yml
    """
}
