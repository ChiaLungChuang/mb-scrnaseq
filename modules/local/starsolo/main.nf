// Quantify with STARsolo using the options STAR's documentation gives for
// matching Cell Ranger >= 4: CellRanger4 adapter clipping, 1MM_CR UMI
// collapsing, MultiGeneUMI_CR filtering, pseudocount barcode correction and
// EmptyDrops_CR cell calling. No BAM is written (--outSAMtype None).
process STARSOLO {
    tag "${meta.id}"
    label 'process_high'
    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_8'

    input:
    tuple val(meta), path(fastq), path(reads_tsv), path(whitelist)
    path index

    output:
    tuple val(meta), path("${meta.id}.Solo.out/Gene/raw"), path("${meta.id}.Solo.out/Gene/filtered"), emit: counts
    tuple val(meta), path("${meta.id}.Solo.out/Gene/Summary.csv"), emit: summary
    tuple val(meta), path("${meta.id}.Log.final.out"), emit: log
    path "versions.yml", emit: versions

    script:
    def umi_len = meta.chemistry == 'v2' ? 10 : 12
    """
    # Comma-separated lane lists, cDNA first then barcode read, in matching order.
    r1=\$(cut -f1 ${reads_tsv} | paste -sd, -)
    cdna=\$(cut -f2 ${reads_tsv} | paste -sd, -)

    case "${whitelist}" in
        *.gz) gzip -cd ${whitelist} > whitelist.txt ;;
        *)    cp ${whitelist} whitelist.txt ;;
    esac

    STAR --runThreadN ${task.cpus} \\
        --genomeDir ${index} \\
        --readFilesIn \$cdna \$r1 \\
        --readFilesCommand zcat \\
        --soloType CB_UMI_Simple \\
        --soloCBwhitelist whitelist.txt \\
        --soloCBstart 1 --soloCBlen 16 \\
        --soloUMIstart 17 --soloUMIlen ${umi_len} \\
        --soloFeatures Gene \\
        --clipAdapterType CellRanger4 \\
        --outFilterScoreMin 30 \\
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
        --soloUMIfiltering MultiGeneUMI_CR \\
        --soloUMIdedup 1MM_CR \\
        --soloCellFilter ${params.solo_cell_filter} \\
        --outSAMtype None \\
        --outFileNamePrefix ${meta.id}. \\
        ${task.ext.args ?: ''}

    rm -f whitelist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version)
    END_VERSIONS
    """

    stub:
    """
    for d in raw filtered; do
        mkdir -p ${meta.id}.Solo.out/Gene/\$d
        touch ${meta.id}.Solo.out/Gene/\$d/matrix.mtx ${meta.id}.Solo.out/Gene/\$d/barcodes.tsv ${meta.id}.Solo.out/Gene/\$d/features.tsv
    done
    touch ${meta.id}.Solo.out/Gene/Summary.csv ${meta.id}.Log.final.out
    printf '"%s":\\n    star: stub\\n' "${task.process}" > versions.yml
    """
}
