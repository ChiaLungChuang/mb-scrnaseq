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
 * fails in seconds, before any job is submitted, rather than hours into a run.
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
        if (row.published_label && !(row.published_label ==~ /[A-Za-z0-9._-]+/)) {
            error("samplesheet line ${i + 2}: published_label '${row.published_label}' may only contain letters, digits, '.', '_' or '-'")
        }
        def meta = [id: row.sample, chemistry: row.chemistry, published_label: row.published_label ?: '']
        // Relative BAM paths are read relative to the samplesheet itself.
        def bamPath = (row.bam.startsWith('/') || row.bam.contains('://')) ? row.bam : "${sheet.parent}/${row.bam}"
        rows << [meta, file(bamPath, checkIfExists: true)]
    }
    return rows
}

/*
 * Pre-registration gate for the concordance step. On real data it refuses
 * to run unless:
 *   - the criteria file parses as YAML with top-level `status: locked`,
 *   - the file is tracked by git and has no uncommitted changes,
 *   - its last commit is on a remote branch or a tag (i.e. pushed, so the
 *     lock is in the public history, not only in a local commit),
 *   - the pipeline code itself has no uncommitted changes, so the results
 *     can be traced to a commit.
 * This makes the order of events checkable, not tamper-proof: git commit
 * dates are set by the author. Returns the facts recorded in provenance.
 */
def checkPreregistration(criteriaPath, codeDir, allowUnregistered) {
    def f = file(criteriaPath, checkIfExists: true)
    def doc = null
    try {
        doc = new org.yaml.snakeyaml.Yaml().load(f.text)
    } catch (Exception e) {
        error("cannot parse the acceptance criteria ${f} as YAML: ${e.message}")
    }
    def locked = (doc instanceof Map) && doc.status == 'locked'
    def dir = f.parent.toString()
    def tracked = runCmd(['git', '-C', dir, 'ls-files', '--error-unmatch', f.name]).code == 0
    def dirty = runCmd(['git', '-C', dir, 'status', '--porcelain', '--', f.name]).out
    def lastCommit = runCmd(['git', '-C', dir, 'log', '-1', '--format=%H %cI', '--', f.name]).out
    def lastSha = lastCommit ? lastCommit.tokenize(' ')[0] : null
    def pushedTo = []
    if (lastSha) {
        pushedTo += runCmd(['git', '-C', dir, 'branch', '-r', '--contains', lastSha]).out.readLines()*.trim().findAll { b -> b && !b.contains('->') }
        pushedTo += runCmd(['git', '-C', dir, 'tag', '--contains', lastSha]).out.readLines()*.trim().findAll { t -> t }.collect { t -> "tag:${t}" }
    }
    def codeDirty = gitDirty(codeDir)
    def problems = []
    if (!locked) problems << "it does not have a top-level 'status: locked'"
    if (!tracked) problems << 'it is not tracked by git'
    if (dirty) problems << 'it has uncommitted changes'
    if (tracked && !dirty && !pushedTo) problems << 'its last commit is not on any remote branch or tag (push it first)'
    if (codeDirty == null) problems << 'the pipeline directory is not a git checkout, so the code cannot be traced to a commit'
    if (codeDirty) problems << 'the pipeline code has uncommitted changes, so results could not be traced to a commit'
    if (problems && !allowUnregistered) {
        error("Refusing to run the concordance check: the run is not pre-registered (${problems.join('; ')}).\n" +
              "Criteria file: ${f}\n" +
              "Review the thresholds, set 'status: locked', commit and push the file, commit any code changes, then re-run. " +
              "(--allow_unlocked_criteria true overrides this for test runs only; the verdict is then marked preregistered: false.)")
    }
    // With the gate switched off the run is never pre-registered, even when
    // every check above passes (the real criteria may well be locked): a
    // test run must not be mistaken for a reportable one.
    if (allowUnregistered) {
        problems << 'it was run with --allow_unlocked_criteria true (test runs only)'
        log.warn("This run is NOT pre-registered (${problems.join('; ')}). Do not report these results.")
    }
    def sha = java.security.MessageDigest.getInstance('SHA-256').digest(f.bytes).encodeHex().toString()
    return [path: f.toString(), sha256: sha, locked: locked, tracked: tracked, uncommitted_changes: dirty as boolean,
            last_commit: lastCommit ?: null, pushed_to: pushedTo, code_uncommitted_changes: codeDirty,
            problems: problems, preregistered: problems.isEmpty()]
}

/*
 * Launch checks for one server without a scheduler (-profile server,
 * conf/server.config). Refuses to start, before anything runs, when the
 * cores or memory this run may use are missing or more than the machine
 * has, when the memory is below what the largest steps ask for, or when the
 * work directory's file system cannot hold one tumor: the FASTQs of the
 * largest BAM (1.7 times it, the upper end measured in tests/scale), a STAR
 * index built here and its published copy (about 10 times the genome FASTA
 * each), and a margin.
 */
def checkServer(samples, workDir) {
    if (!params.server_cpus || !params.server_memory) {
        error("-profile server needs --server_cpus and --server_memory: the cores and memory this run may use " +
              "on this machine (docs/RUNBOOK.md, \"Route 2: a standalone server\").")
    }
    def cpus = params.server_cpus as int
    def mem = nextflow.util.MemoryUnit.of(params.server_memory as String)
    def minMem = nextflow.util.MemoryUnit.of(params.server_min_memory as String)
    def machineCpus = Runtime.getRuntime().availableProcessors()
    def machineMem = machineMemory()
    def problems = []
    if (cpus < 1 || cpus > machineCpus) {
        problems << "--server_cpus ${cpus} is not between 1 and this machine's ${machineCpus} cores"
    }
    if (mem < minMem) {
        problems << "--server_memory ${mem} is below ${minMem}, what the STAR index build and STARsolo ask for"
    }
    if (machineMem != null && mem > machineMem) {
        problems << "--server_memory ${mem} is more than this machine has (${machineMem})"
    }
    def bamMax = samples.collect { s -> s[1].size() }.max() ?: 0L
    def fastqBytes = (1.7 * bamMax) as long
    def indexBytes = params.star_index ? 0L : 2L * 10L * file(params.genome_fasta).size()
    def marginBytes = nextflow.util.MemoryUnit.of(params.server_disk_margin as String).toBytes()
    def needBytes = fastqBytes + indexBytes + marginBytes
    def freeBytes = usableSpace(workDir)
    if (freeBytes < needBytes) {
        problems << "the work directory's file system has ${sizeText(freeBytes)} free, and one tumor needs about ${sizeText(needBytes)} " +
                    "(FASTQs ${sizeText(fastqBytes)}, STAR index and its copy ${sizeText(indexBytes)}, margin ${sizeText(marginBytes)})"
    }
    if (problems) {
        error("Refusing to start on this server: ${problems.join('; ')}. See docs/RUNBOOK.md, \"Route 2: a standalone server\".")
    }
    log.info("server: using ${cpus} of ${machineCpus} cores and ${mem} of ${machineMem ?: 'unknown'} memory; " +
             "${sizeText(freeBytes)} free for work/, one tumor needs about ${sizeText(needBytes)}")
    return [cpus: cpus, memory: mem.toString(), machine_cpus: machineCpus, machine_memory: machineMem?.toString(),
            work_free: sizeText(freeBytes), one_tumor_needs: sizeText(needBytes),
            needs: [fastq: sizeText(fastqBytes), star_index_and_copy: sizeText(indexBytes), margin: sizeText(marginBytes)]]
}

def sizeText(bytes) {
    return new nextflow.util.MemoryUnit(bytes as long).toString()
}

// Total memory from /proc/meminfo (Linux), or null where there is none.
def machineMemory() {
    def f = new java.io.File('/proc/meminfo')
    if (!f.exists()) {
        return null
    }
    def m = (f.text =~ /MemTotal:\s+(\d+) kB/)
    return m.find() ? new nextflow.util.MemoryUnit((m.group(1) as long) * 1024L) : null
}

// Free space on the file system that holds `path`, or on its nearest
// existing parent (the work directory may not exist yet at launch).
def usableSpace(path) {
    def f = new java.io.File(path.toString()).getAbsoluteFile()
    if (f.exists()) {
        return f.getUsableSpace()
    }
    return f.getParentFile() ? usableSpace(f.getParentFile()) : 0L
}

/*
 * SHA-256 over the files matching `globs` under `dir` (names and contents).
 * Passed as an input to the steps that run that code, so -resume re-runs
 * them after a code change: Nextflow hashes neither bin/ nor the contents
 * of a staged directory such as lib/R. Each step gets a fingerprint of only
 * the code it runs, so an R change does not re-run STARsolo.
 */
def codeFingerprint(dir, globs) {
    def md = java.security.MessageDigest.getInstance('SHA-256')
    def fs = globs.collectMany { g -> files("${dir}/${g}") }.unique().sort { f -> f.toString() }
    fs.each { f ->
        md.update("${dir.relativize(f)}\u0000".getBytes('UTF-8'))
        md.update(f.bytes)
    }
    return md.digest().encodeHex().toString()
}

workflow {
    main:
    def doConcordance = params.published_counts as boolean
    if (!(params.cellranger_compat in ['cr3', 'cr4'])) {
        error("--cellranger_compat must be cr3 or cr4: the Cell Ranger major version that produced the BAMs " +
              "(see the @PG/@CO lines of a BAM header). Got '${params.cellranger_compat}'.")
    }
    def samples = readSamplesheet(params.input, doConcordance)
    if (boolParam('server_check', params.server_check)) {
        writeJson("${params.outdir}/pipeline_info/server_check.json", checkServer(samples, workflow.workDir))
    }
    def chemistries = samples.collect { s -> s[0].chemistry }.unique()
    chemistries.each { c ->
        if (!params["whitelist_${c}"]) {
            error("samples use chemistry ${c} but --whitelist_${c} is not set")
        }
    }
    def codeSha = codeFingerprint(projectDir, ['lib/R/*.R', 'bin/*'])
    def readsCheckSha = codeFingerprint(projectDir, ['bin/check_reads.sh'])
    def rSha = codeFingerprint(projectDir, ['lib/R/*.R', 'bin/*.R'])
    def prereg = null
    if (doConcordance) {
        // Recorded at launch, so the gate's view of the criteria is kept even if the run fails.
        prereg = checkPreregistration(params.criteria, projectDir, boolParam('allow_unlocked_criteria', params.allow_unlocked_criteria))
        writeJson("${params.outdir}/pipeline_info/preregistration_check.json", prereg + [code_sha256: codeSha])
    }

    def lib = channel.value(file("${projectDir}/lib/R", checkIfExists: true))

    // Reference index: reuse a pre-built one when given, otherwise build it.
    def index
    def versions = channel.empty()
    if (params.star_index) {
        index = channel.value(file(params.star_index, checkIfExists: true))
    } else {
        STAR_GENOMEGENERATE(file(params.genome_fasta, checkIfExists: true), file(params.genome_gtf, checkIfExists: true))
        index = STAR_GENOMEGENERATE.out.index
        versions = versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    BAMTOFASTQ(channel.fromList(samples))
    CHECK_READS(BAMTOFASTQ.out.fastq, readsCheckSha)

    def solo_in = CHECK_READS.out.reads.map { meta, fastq, tsv ->
        [meta, fastq, tsv, file(params["whitelist_${meta.chemistry}"], checkIfExists: true)]
    }
    STARSOLO(solo_in, index)

    SAMPLE_QC(STARSOLO.out.counts.map { meta, _raw, filtered -> [meta, filtered] }, lib, rSha)

    if (doConcordance) {
        PREP_PUBLISHED(file(params.published_counts, checkIfExists: true), lib, rSha)
        CONCORDANCE(STARSOLO.out.counts, PREP_PUBLISHED.out.rds, file(params.criteria), lib, "${rSha}:${prereg.sha256}")
        SUMMARISE_CONCORDANCE(CONCORDANCE.out.json.map { _meta, json -> json }.collect(), lib,
                              prereg.preregistered, prereg.sha256, codeSha)
    }

    // One entry per tool: every sample writes an identical versions.yml.
    versions
        .mix(BAMTOFASTQ.out.versions, STARSOLO.out.versions, SAMPLE_QC.out.versions)
        .map { f -> f.text }
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true, newLine: false)

    onComplete:
    def preregFile = file("${params.outdir}/pipeline_info/preregistration_check.json")
    def info = [
        pipeline: workflow.manifest.name,
        pipeline_version: workflow.manifest.version,
        git_commit: workflow.commitId ?: gitHead(projectDir),
        git_uncommitted_changes: gitDirty(projectDir),
        code_sha256: codeFingerprint(projectDir, ['lib/R/*.R', 'bin/*']),
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
        preregistration: preregFile.exists() ? new groovy.json.JsonSlurper().parseText(preregFile.text) : null,
    ]
    writeJson("${params.outdir}/pipeline_info/run_provenance.json", info)
    log.info("Provenance written to ${params.outdir}/pipeline_info/run_provenance.json")
    if (params.published_counts && workflow.success) {
        log.info("Concordance verdict: ${params.outdir}/concordance/concordance_summary.verdict.json")
    }
}

// A boolean parameter given on the command line arrives as the string
// "false", which Groovy treats as true. Accept only true/false.
def boolParam(name, v) {
    if (v instanceof Boolean) {
        return v
    }
    def s = v == null ? 'false' : v.toString().trim().toLowerCase()
    if (s == 'true') {
        return true
    }
    if (s == 'false') {
        return false
    }
    error("--${name} must be true or false, got '${v}'")
}

def writeJson(path, obj) {
    def out = file(path)
    out.parent.mkdirs()
    out.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(obj))
}

// True when the checkout has uncommitted changes (null if git is unavailable).
def gitDirty(dir) {
    def r = runCmd(['git', '-C', dir.toString(), 'status', '--porcelain', '--untracked-files=no'])
    return r.code == 0 ? (r.out as boolean) : null
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
