# Phase 1 results report (lib/R/report.R) on hand-built output directories
# whose answers are known: compute is worked out by hand below.

GB <- 1024^3
H <- 3.6e6   # one hour in ms

# A finished run's outdir: sampleB fails C3, everything else passes.
make_outdir <- function(preregistered = TRUE, trace_rows = NULL, verdict_edit = identity, prov_edit = identity) {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  dir.create(file.path(d, "concordance")); dir.create(file.path(d, "pipeline_info")); dir.create(file.path(d, "qc"))
  crit <- data.frame(
    criterion = c("C1_barcode_recovery", "C3_pseudobulk_pearson"),
    metric = c("barcode_recovery", "pseudobulk_pearson"),
    rule = c(">= 0.99", ">= 0.98"), stringsAsFactors = FALSE)
  tab <- rbind(
    cbind(sample = "sampleA", label = "943", crit, observed = c(0.998, 0.991), pass = c(TRUE, TRUE)),
    cbind(sample = "sampleB", label = "1355", crit, observed = c(0.995, 0.95), pass = c(TRUE, FALSE)))
  tab <- tab[, c("sample", "label", "criterion", "metric", "observed", "rule", "pass")]
  utils::write.table(tab, file.path(d, "concordance", "concordance_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  verdict <- verdict_edit(list(preregistered = preregistered, criteria_sha256 = strrep("c", 64), code_sha256 = strrep("a", 64),
                               passed = FALSE, n_samples = 2L, n_samples_passed = 1L,
                               failing_samples = list("sampleB"), failing_criteria = list("C3_pseudobulk_pearson")))
  write_json(verdict, file.path(d, "concordance", "concordance_summary.verdict.json"))
  prov <- prov_edit(list(
    pipeline_version = "0.1.0", git_commit = strrep("e", 40), git_uncommitted_changes = FALSE,
    code_sha256 = strrep("a", 64), started = "2026-12-01T10:00:00.000Z", completed = "2026-12-01T13:30:00.000Z",
    success = TRUE, nextflow_version = "26.10.4", profile = "stjude", container_engine = "singularity",
    params = list(r_container = "ghcr.io/x/mb-scrnaseq-r:0.1.0@sha256:abc"),
    preregistration = list(path = "/research/lab/mb-scrnaseq/preregistration/acceptance_criteria.yaml",
                           sha256 = strrep("c", 64), locked = preregistered, preregistered = preregistered,
                           last_commit = "f00 2026-11-20", problems = if (preregistered) list() else list("it does not have a top-level 'status: locked'"),
                           pushed_to = list("origin/feature-x", "origin/main", "tag:prereg-phase1"))))
  write_json(prov, file.path(d, "pipeline_info", "run_provenance.json"))
  if (is.null(trace_rows)) {
    trace_rows <- data.frame(
      process = c("STARSOLO", "STARSOLO", "STARSOLO", "CONCORDANCE", "SAMPLE_QC"),
      status = c("COMPLETED", "FAILED", "COMPLETED", "CACHED", "COMPLETED"),
      attempt = c(1, 1, 2, 1, 1), cpus = c(8, 8, 8, 1, 1),
      memory = c(48, 48, 96, 16, 16) * GB,
      realtime = c(2 * H, 0.5 * H, 1 * H, 1 * H, 0.5 * H),
      `%cpu` = c(600, 400, 700, 100, 95),
      peak_rss = c(30, 50, 52, 10, 2) * GB, check.names = FALSE)
  }
  utils::write.table(trace_rows, file.path(d, "pipeline_info", "trace.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write_json(list(nmads = 3, total_umi = list(lower = 501.2), detected_genes = list(lower = 301),
                  pct_mito = list(higher = 12.5), n_cells = 1500, n_flagged = 40),
             file.path(d, "qc", "sampleA.qc_thresholds.json"))
  writeLines(c('"STARSOLO":', '    star: 2.7.11b'), file.path(d, "pipeline_info", "software_versions.yml"))
  d
}

test_that("compute per process is worked out from the raw trace exactly", {
  cp <- compute_by_process(read_trace_raw(file.path(make_outdir(), "pipeline_info", "trace.tsv")))
  s <- cp$by_process[cp$by_process$process == "STARSOLO", ]
  expect_equal(s$tasks, 3)
  expect_equal(s$failed, 1)
  expect_equal(s$retried, 1)
  expect_equal(s$cpu_h_used, 2 * 6 + 0.5 * 4 + 1 * 7)      # 21 h
  expect_equal(s$cpu_h_reserved, (2 + 0.5 + 1) * 8)        # 28 h
  expect_equal(s$realtime_ms_sum, 3.5 * H)
  expect_equal(s$realtime_ms_max, 2 * H)
  expect_equal(s$peak_rss_max, 52 * GB)
  expect_equal(s$memory_max, 96 * GB)
  # The cached CONCORDANCE task ran in an earlier run: not in this run's compute.
  expect_false("CONCORDANCE" %in% cp$by_process$process)
  expect_equal(cp$total$tasks, 4)
  expect_equal(cp$total$cpu_h_used, 21 + 0.5 * 0.95)
  expect_equal(cp$total$cpu_h_reserved, 28 + 0.5)
  expect_equal(cp$cached$tasks, 1)
  expect_equal(cp$cached$cpu_h_used, 1)
})

test_that("the report states the verdict per sample and criterion, the run and the compute", {
  md <- paste(phase1_report(make_outdir()), collapse = "\n")
  expect_match(md, "**FAILED**: 1 of 2 samples passed every criterion.", fixed = TRUE)
  expect_match(md, "Pre-registered: **yes**", fixed = TRUE)
  expect_match(md, "| sampleB | 1355 | PASS 0.9950 | **FAIL** 0.9500 | **FAIL** |", fixed = TRUE)
  expect_match(md, "| sampleA | 943 | PASS 0.9980 | PASS 0.9910 | PASS |", fixed = TRUE)
  expect_match(md, "C3_pseudobulk_pearson (>= 0.98)", fixed = TRUE)
  expect_match(md, "Failing samples: sampleB. Failing criteria: C3_pseudobulk_pearson.", fixed = TRUE)
  expect_match(md, paste0("| code SHA-256 | `", strrep("a", 64), "` |"), fixed = TRUE)
  expect_match(md, "| criteria file | `preregistration/acceptance_criteria.yaml` |", fixed = TRUE)
  expect_match(md, "| criteria pushed to (main, tags) | origin/main, tag:prereg-phase1 |", fixed = TRUE)
  expect_no_match(md, "feature-x", fixed = TRUE)
  expect_match(md, "| wall time | 3 h 30 min |", fixed = TRUE)
  expect_match(md, "This run: 4 tasks (1 failed, 1 retried), 21.475 CPU hours used, 28.500 CPU hours reserved", fixed = TRUE)
  expect_match(md, "Not counted above: 1 cached tasks", fixed = TRUE)
  expect_match(md, "| STARSOLO | 3 | 1 | 1 | 21.000 | 28.000 | 3 h 30 min | 2 h 00 min | 52.0 GB | 96.0 GB |", fixed = TRUE)
  expect_match(md, "| sampleA | 1500 | 40 | 3 | 501.2 | 301.0 | 12.50 |", fixed = TRUE)
  expect_match(md, "| STARSOLO | star | 2.7.11b |", fixed = TRUE)
})

test_that("a run that was not pre-registered is flagged as not reportable, with the reason", {
  md <- paste(phase1_report(make_outdir(preregistered = FALSE)), collapse = "\n")
  expect_match(md, "**Not pre-registered. This verdict is not a reportable result.**", fixed = TRUE)
  expect_match(md, "Reason: it does not have a top-level 'status: locked'.", fixed = TRUE)
  expect_no_match(md, "Pre-registered: **yes**", fixed = TRUE)
})

test_that("outputs that do not belong together are refused", {
  expect_error(phase1_report(make_outdir(prov_edit = function(p) { p$code_sha256 <- strrep("b", 64); p })),
               "code_sha256 differs")
  expect_error(phase1_report(make_outdir(verdict_edit = function(v) { v$criteria_sha256 <- strrep("d", 64); v })),
               "criteria sha256 differs")
  expect_error(phase1_report(make_outdir(verdict_edit = function(v) { v$failing_samples <- list("sampleA"); v })),
               "failing samples")
  expect_error(phase1_report(make_outdir(verdict_edit = function(v) { v$n_samples <- 3L; v })),
               "verdict says 3 samples")
  expect_error(phase1_report(make_outdir(verdict_edit = function(v) { v$passed <- TRUE; v$failing_samples <- list(); v })),
               "disagree|failing samples")
  expect_error(phase1_report(make_outdir(verdict_edit = function(v) { v$preregistered <- FALSE; v })),
               "pre-registration status differs")
})

test_that("a human-formatted or incomplete trace and a missing verdict are refused", {
  human <- data.frame(process = "STARSOLO", status = "COMPLETED", attempt = 1, cpus = 8, memory = "48 GB",
                      realtime = "2h", `%cpu` = "600%", peak_rss = "30 GB", check.names = FALSE)
  expect_error(phase1_report(make_outdir(trace_rows = human)), "trace.raw = true")
  short <- data.frame(process = "STARSOLO", status = "COMPLETED", realtime = 1)
  expect_error(phase1_report(make_outdir(trace_rows = short)), "lacks column")
  d <- make_outdir()
  file.remove(file.path(d, "concordance", "concordance_summary.verdict.json"))
  expect_error(phase1_report(d), "not a finished phase 1 run")
})

test_that("durations and sizes are formatted at their unit boundaries", {
  expect_equal(fmt_duration(c(999, 1000, 59900, 60000, 3599000, 3600000, 5400000)),
               c("999 ms", "1.0 s", "59.9 s", "1 min 00 s", "59 min 59 s", "1 h 00 min", "1 h 30 min"))
  expect_equal(fmt_bytes(c(0, 1023, 1024, 1.5 * 1024^2, 52 * GB, 2 * 1024^4)),
               c("0 B", "1023 B", "1.0 KB", "1.5 MB", "52.0 GB", "2.0 TB"))
})
