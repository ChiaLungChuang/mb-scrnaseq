#!/usr/bin/env nextflow
/*
 * mb-scrnaseq, phase 1: re-quantify childhood medulloblastoma 10x scRNA-seq
 * (GEO GSE155446) from the deposited Cell Ranger BAMs and test whether the
 * result reproduces the published count matrix.
 *
 *   BAM -> bamtofastq -> read/chemistry check -> STARsolo -> per-sample QC
 *                                                        -> concordance vs published counts
 *
 * See README.md for usage and docs/RUNBOOK.md for running on the St. Jude cluster.
 */

include { BAMTOFASTQ            } from './modules/local/bamtofastq'
include { CHECK_READS           } from './modules/local/check_reads'
include { STAR_GENOMEGENERATE   } from './modules/local/star_genomegenerate'
include { STARSOLO              } from './modules/local/starsolo'
include { SAMPLE_QC             } from './modules/local/sample_qc'
include { PREP_PUBLISHED        } from './modules/local/prep_published'
include { CONCORDANCE           } from './modules/local/concordance'
include { SUMMARISE_CONCORDANCE } from './modules/local/summarise_concordance'

/*
 * Read and validate the samplesheet before anything is submitted, so a typo
 * fails in seconds on the login node rather than hours into a run.
 * Columns: sample, bam, chemistry (v2|v3), published_label (needed only when
 * --published_counts is set).
 */
def readSamplesheet(path, needLabel) {
    def sheet = file(path, checkIfExists: true)
    def lines = sheet.readLines().findAll { line -> line.trim() && !line.startsWith('#') }
    if (lines.size() < 2) {
        error("samplesheet ${path} has no sample rows")
    }
    def header = lines[0].split(',', -1).collect { h -> h.trim() }
    def required = ['sample', 'bam', 'chemistry'] + (needLabel ? ['published_label'] : [])
    def missingCols = required.findAll { c -> !(c in header) }
    if (missingCols) {
        error("samplesheet ${path} is missing column(s): ${missingCols.join(', ')}")
    }
    def seen = [] as Set
    def rows = []
    lines.drop(1).eachWithIndex { line, i ->
        def fields = line.split(',', -1).collect { f -> f.trim() }
        if (fields.size() != header.size()) {
            error("samplesheet line ${i + 2}: expected ${header.size()} fields, found ${fields.size()}")
        }
        def row = [header, fields].transpose().collectEntries { kv -> [(kv[0]): kv[1]] }
        if (!(row.sample ==~ /[A-Za-z0-9._-]+/)) {
            error("samplesheet line ${i + 2}: sample '${row.sample}' may only contain letters, digits, '.', '_' or '-'")
        }
        if (row.sample in seen) {
            error("samplesheet line ${i + 2}: duplicate sample '${row.sample}'")
        }
        seen << row.sample
        if (!(row.chemistry in ['v2', 'v3'])) {
            error("samplesheet line ${i + 2}: chemistry must be v2 or v3, got '${row.chemistry}'")
        }
        if (needLabel && !row.published_label) {
            error("samplesheet line ${i + 2}: published_label is empty; it is required for the concordance check")
        }
        def meta = [id: row.sample, chemistry: row.chemistry, published_label: row.published_label ?: '']
        // Relative BAM paths are read relative to the samplesheet itself.
        def bamPath = (row.bam.startsWith('/') || row.bam.contains('://')) ? row.bam : "${sheet.parent}/${row.bam}"
        rows << [meta, file(bamPath, checkIfExists: true)]
    }
    return rows
}

/*
 * Pre-registration gate. The acceptance criteria must be committed to git,
 * unmodified, and marked `status: locked` before the concordance step may
 * run, so the thresholds demonstrably predate the results. Returns facts
 * recorded in the provenance file.
 */
def checkPreregistration(criteriaPath, allowUnlocked) {
    def f = file(criteriaPath, checkIfExists: true)
    def locked = f.readLines().any { line -> line.replaceAll(/#.*$/, '').trim() == 'status: locked' }
    def dir = f.parent.toString()
    def tracked = runCmd(['git', '-C', dir, 'ls-files', '--error-unmatch', f.name]).code == 0
    def dirty = runCmd(['git', '-C', dir, 'status', '--porcelain', '--', f.name]).out
    def lastCommit = runCmd(['git', '-C', dir, 'log', '-1', '--format=%H %cI', '--', f.name]).out
    def sha = java.security.MessageDigest.getInstance('SHA-256').digest(f.bytes).encodeHex().toString()
    def problems = []
    if (!locked) problems << "it is not marked 'status: locked'"
    if (!tracked) problems << 'it is not tracked by git'
    if (dirty) problems << 'it has uncommitted changes'
    if (problems && !allowUnlocked) {
        error("Refusing to run the concordance check: the acceptance criteria file ${f} is not pre-registered (${problems.join('; ')}).\n" +
              "Review the thresholds, set 'status: locked', commit the file, then re-run. " +
              "(--allow_unlocked_criteria true overrides this for test runs only.)")
    }
    if (problems) {
        log.warn("Acceptance criteria are NOT pre-registered (${problems.join('; ')}); allowed by --allow_unlocked_criteria. Do not report these results.")
    }
    return [path: f.toString(), sha256: sha, locked: locked, tracked: tracked, uncommitted_changes: dirty as boolean,
            last_commit: lastCommit ?: null, preregistered: problems.isEmpty()]
}

workflow {
    main:
    def doConcordance = params.published_counts as boolean
    def samples = readSamplesheet(params.input, doConcordance)
    def chemistries = samples.collect { s -> s[0].chemistry }.unique()
    chemistries.each { c ->
        if (!params["whitelist_${c}"]) {
            error("samples use chemistry ${c} but --whitelist_${c} is not set")
        }
    }
    def prereg = doConcordance ? checkPreregistration(params.criteria, params.allow_unlocked_criteria) : null

    def lib = channel.value(file("${projectDir}/lib/R", checkIfExists: true))

    // Reference index: reuse a pre-built one when given, otherwise build it.
    def index
    if (params.star_index) {
        index = channel.value(file(params.star_index, checkIfExists: true))
    } else {
        STAR_GENOMEGENERATE(file(params.genome_fasta, checkIfExists: true), file(params.genome_gtf, checkIfExists: true))
        index = STAR_GENOMEGENERATE.out.index
    }

    BAMTOFASTQ(channel.fromList(samples))
    CHECK_READS(BAMTOFASTQ.out.fastq)

    def solo_in = CHECK_READS.out.reads.map { meta, fastq, tsv ->
        [meta, fastq, tsv, file(params["whitelist_${meta.chemistry}"], checkIfExists: true)]
    }
    STARSOLO(solo_in, index)

    SAMPLE_QC(STARSOLO.out.counts.map { meta, _raw, filtered -> [meta, filtered] }, lib)

    if (doConcordance) {
        PREP_PUBLISHED(file(params.published_counts, checkIfExists: true), lib)
        CONCORDANCE(STARSOLO.out.counts, PREP_PUBLISHED.out.rds, file(params.criteria), lib)
        SUMMARISE_CONCORDANCE(CONCORDANCE.out.json.map { _meta, json -> json }.collect(), lib)
    }

    BAMTOFASTQ.out.versions
        .mix(STARSOLO.out.versions, SAMPLE_QC.out.versions)
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true, newLine: false)

    workflow.onComplete {
        def info = [
            pipeline: workflow.manifest.name,
            pipeline_version: workflow.manifest.version,
            git_commit: workflow.commitId ?: gitHead(projectDir),
            run_name: workflow.runName,
            session_id: workflow.sessionId.toString(),
            started: workflow.start.toString(),
            completed: workflow.complete.toString(),
            success: workflow.success,
            exit_status: workflow.exitStatus,
            nextflow_version: workflow.nextflow.version.toString(),
            profile: workflow.profile,
            container_engine: workflow.containerEngine,
            command_line: workflow.commandLine,
            params: params.findAll { _k, v -> !(v instanceof Closure) }.collectEntries { k, v -> [(k): v?.toString()] },
            preregistration: prereg,
        ]
        def out = file("${params.outdir}/pipeline_info/run_provenance.json")
        out.parent.mkdirs()
        out.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info))
        log.info("Provenance written to ${out}")
        if (doConcordance && workflow.success) {
            log.info("Concordance verdict: ${params.outdir}/concordance/concordance_summary.verdict.json")
        }
    }
}

// Run a command and capture its exit code and trimmed stdout. A missing
// binary (e.g. no git) is reported as a failed command, not a crash.
def runCmd(cmd) {
    try {
        def p = cmd.execute()
        def out = p.text.trim()
        p.waitFor()
        return [code: p.exitValue(), out: out]
    } catch (java.io.IOException _e) {
        return [code: -1, out: '']
    }
}

// Current commit of a local checkout (workflow.commitId is only set for remote runs).
def gitHead(dir) {
    def r = runCmd(['git', '-C', dir.toString(), 'rev-parse', 'HEAD'])
    return r.code == 0 ? r.out : null
}
