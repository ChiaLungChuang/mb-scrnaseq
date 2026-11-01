# Phase 1 results report: turn a finished run's output directory into the
# Markdown that goes to docs/results/phase1/README.md (plan M3 step 5).
#
# Everything is read from the run's own outputs; nothing is computed from
# the published data. The report refuses an output directory whose files
# disagree with each other (verdict vs provenance hashes, verdict vs table),
# so it cannot describe a mixture of runs.

#' Read a Nextflow trace written with `trace.raw = true`.
#'
#' Raw traces give times in milliseconds, memory in bytes and %cpu as a
#' plain number. A human-formatted trace ("1h 2m", "1.2 GB") is refused
#' rather than guessed at.
read_trace_raw <- function(path) {
  tr <- data.table::fread(path, sep = "\t", colClasses = "character", na.strings = c("-", ""))
  need <- c("process", "status", "attempt", "cpus", "memory", "realtime", "%cpu", "peak_rss")
  miss <- setdiff(need, names(tr))
  if (length(miss)) {
    stop("trace ", path, " lacks column(s) ", paste(miss, collapse = ", "),
         "; set them in trace.fields (nextflow.config)")
  }
  num <- function(col) {
    x <- tr[[col]]
    v <- suppressWarnings(as.numeric(x))
    bad <- !is.na(x) & is.na(v)
    if (any(bad)) {
      stop(sprintf("trace column '%s' has non-numeric value '%s': write the trace with trace.raw = true",
                   col, x[bad][1]))
    }
    v
  }
  data.frame(
    process = tr$process, status = tr$status, attempt = num("attempt"), cpus = num("cpus"),
    memory = num("memory"), realtime_ms = num("realtime"), pct_cpu = num("%cpu"),
    peak_rss = num("peak_rss"), stringsAsFactors = FALSE
  )
}

#' Compute per process from a raw trace.
#'
#' CPU hours used = realtime x %cpu / 100; CPU hours reserved = realtime x
#' cpus. Tasks with status CACHED ran in an earlier run (-resume) and are
#' counted separately, not added to this run's compute.
compute_by_process <- function(tr) {
  hours <- function(ms) ms / 3.6e6
  one <- function(d) {
    data.frame(
      tasks = nrow(d),
      failed = sum(d$status %in% c("FAILED", "ABORTED")),
      retried = sum(d$attempt > 1, na.rm = TRUE),
      cpu_h_used = sum(hours(d$realtime_ms) * d$pct_cpu / 100, na.rm = TRUE),
      cpu_h_reserved = sum(hours(d$realtime_ms) * d$cpus, na.rm = TRUE),
      realtime_ms_sum = sum(d$realtime_ms, na.rm = TRUE),
      realtime_ms_max = if (all(is.na(d$realtime_ms))) NA_real_ else max(d$realtime_ms, na.rm = TRUE),
      peak_rss_max = if (all(is.na(d$peak_rss))) NA_real_ else max(d$peak_rss, na.rm = TRUE),
      memory_max = if (all(is.na(d$memory))) NA_real_ else max(d$memory, na.rm = TRUE)
    )
  }
  ran <- tr[tr$status != "CACHED", , drop = FALSE]
  procs <- sort(unique(ran$process))
  by <- do.call(rbind, lapply(procs, function(p) cbind(process = p, one(ran[ran$process == p, , drop = FALSE]))))
  if (is.null(by)) by <- cbind(process = character(), one(ran))[0, ]
  cached <- tr[tr$status == "CACHED", , drop = FALSE]
  list(by_process = by, total = one(ran), cached = one(cached))
}

fmt_bytes <- function(b) {
  vapply(b, function(x) {
    if (is.na(x)) return("NA")
    u <- c("B", "KB", "MB", "GB", "TB")
    i <- if (x < 1) 1 else min(length(u), floor(log(x, 1024)) + 1)
    if (i == 1) sprintf("%d B", as.integer(x)) else sprintf("%.1f %s", x / 1024^(i - 1), u[i])
  }, "")
}

fmt_duration <- function(ms) {
  vapply(ms, function(x) {
    if (is.na(x)) return("NA")
    s <- x / 1000
    if (s < 1) sprintf("%d ms", as.integer(round(x)))
    else if (s < 60) sprintf("%.1f s", s)
    else if (s < 3600) sprintf("%d min %02d s", as.integer(s %/% 60), as.integer(round(s %% 60)))
    else sprintf("%d h %02d min", as.integer(s %/% 3600), as.integer(round((s %% 3600) / 60)))
  }, "")
}

fmt_num <- function(x) ifelse(is.na(x), "NA", formatC(x, digits = 4, format = "g", flag = "#"))
fmt_int <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "d", big.mark = ","))

md_table <- function(df) {
  esc <- function(x) gsub("|", "\\|", as.character(x), fixed = TRUE)
  c(paste0("| ", paste(esc(names(df)), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(esc(r), collapse = " | "), " |")))
}

#' Read everything the report needs and check that it belongs together.
#'
#' `trace` defaults to the run's pipeline_info/trace.tsv; nf-test, for
#' example, writes it elsewhere (-with-trace).
read_phase1_outdir <- function(outdir, trace = NULL) {
  p <- function(...) file.path(outdir, ...)
  need <- c(verdict = p("concordance", "concordance_summary.verdict.json"),
            table = p("concordance", "concordance_summary.tsv"),
            provenance = p("pipeline_info", "run_provenance.json"),
            trace = if (is.null(trace)) p("pipeline_info", "trace.tsv") else trace)
  missing <- need[!file.exists(need)]
  if (length(missing)) stop("not a finished phase 1 run; missing: ", paste(missing, collapse = ", "))

  verdict <- jsonlite::read_json(need[["verdict"]])
  prov <- jsonlite::read_json(need[["provenance"]])
  tab <- data.table::fread(need[["table"]], sep = "\t", colClasses = list(character = c("sample", "label", "criterion", "metric", "rule")))
  tab <- as.data.frame(tab)
  tab$pass <- as.logical(tab$pass)
  if (anyNA(tab$pass)) stop("concordance_summary.tsv has pass values that are not TRUE/FALSE")

  problems <- character()
  if (!identical(verdict$code_sha256, prov$code_sha256)) {
    problems <- c(problems, sprintf("code_sha256 differs: verdict %s, provenance %s", verdict$code_sha256, prov$code_sha256))
  }
  if (!identical(verdict$criteria_sha256, prov$preregistration$sha256)) {
    problems <- c(problems, sprintf("criteria sha256 differs: verdict %s, provenance %s", verdict$criteria_sha256, prov$preregistration$sha256))
  }
  if (!identical(isTRUE(verdict$preregistered), isTRUE(prov$preregistration$preregistered))) {
    problems <- c(problems, "pre-registration status differs between verdict and provenance")
  }
  per_sample <- tapply(tab$pass, tab$sample, all)
  if (!identical(as.integer(verdict$n_samples), length(per_sample))) {
    problems <- c(problems, sprintf("verdict says %s samples, table has %d", verdict$n_samples, length(per_sample)))
  }
  if (!setequal(unlist(verdict$failing_samples), names(per_sample)[!per_sample])) {
    problems <- c(problems, "failing samples in the verdict and the table differ")
  }
  if (!identical(isTRUE(verdict$passed), all(per_sample))) {
    problems <- c(problems, "overall verdict and the table disagree")
  }
  if (length(problems)) stop("inconsistent outputs in ", outdir, ":\n  ", paste(problems, collapse = "\n  "))

  qc_files <- sort(list.files(p("qc"), pattern = "\\.qc_thresholds\\.json$", full.names = TRUE))
  qc <- lapply(qc_files, jsonlite::read_json)
  names(qc) <- sub("\\.qc_thresholds\\.json$", "", basename(qc_files))
  versions <- if (file.exists(p("pipeline_info", "software_versions.yml"))) yaml::read_yaml(p("pipeline_info", "software_versions.yml")) else NULL
  conc_files <- sort(list.files(p("concordance"), pattern = "\\.concordance\\.json$", full.names = TRUE))
  diag <- lapply(conc_files, jsonlite::read_json)

  list(verdict = verdict, provenance = prov, table = tab, trace = read_trace_raw(need[["trace"]]),
       qc = qc, versions = versions, diagnostics = diag)
}

#' The phase 1 report as a character vector of Markdown lines.
phase1_report <- function(outdir, trace = NULL) {
  x <- read_phase1_outdir(outdir, trace)
  v <- x$verdict; prov <- x$provenance; pre <- prov$preregistration; tab <- x$table
  rel <- function(path) sub("^.*/(preregistration/)", "\\1", path)
  yn <- function(b) if (isTRUE(b)) "yes" else "no"
  per_sample <- tapply(tab$pass, tab$sample, all)

  out <- c("# Phase 1 results", "",
           "Generated by `bin/phase1_report.R` from the run's output directory. Do not edit by hand; re-generate it.", "")

  out <- c(out, "## Verdict", "",
           sprintf("**%s**: %d of %d samples passed every criterion.",
                   if (isTRUE(v$passed)) "PASSED" else "FAILED", sum(per_sample), length(per_sample)), "")
  if (isTRUE(v$preregistered)) {
    out <- c(out, "Pre-registered: **yes**. The criteria file was locked, committed, unmodified and pushed before this run, and the code had no uncommitted changes.", "")
  } else {
    why <- unlist(pre$problems)
    out <- c(out, "> **Not pre-registered. This verdict is not a reportable result.**",
             if (length(why)) paste0("> Reason: ", paste(why, collapse = "; "), ".") else NULL, "")
  }
  if (length(unlist(v$failing_criteria))) {
    out <- c(out, sprintf("Failing samples: %s. Failing criteria: %s.",
                          paste(unlist(v$failing_samples), collapse = ", "),
                          paste(unlist(v$failing_criteria), collapse = ", ")), "")
  }

  wall_ms <- tryCatch(1000 * as.numeric(difftime(as.POSIXct(prov$completed, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC"),
                                                  as.POSIXct(prov$started, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC"),
                                                  units = "secs")), error = function(e) NA_real_)
  pushed <- unlist(pre$pushed_to)
  pushed <- pushed[grepl("^tag:|/main$", pushed)]
  run <- data.frame(
    item = c("pipeline version", "git commit", "uncommitted code changes", "code SHA-256",
             "criteria file", "criteria SHA-256", "criteria locked", "criteria last commit",
             "criteria pushed to (main, tags)", "Nextflow", "profile", "container engine", "R image",
             "started (UTC)", "completed (UTC)", "wall time", "run succeeded"),
    value = c(prov$pipeline_version, sprintf("`%s`", prov$git_commit), yn(prov$git_uncommitted_changes),
              sprintf("`%s`", v$code_sha256), sprintf("`%s`", rel(pre$path)), sprintf("`%s`", v$criteria_sha256),
              yn(pre$locked), pre$last_commit %||% "NA",
              if (length(pushed)) paste(pushed, collapse = ", ") else "none",
              prov$nextflow_version, prov$profile, prov$container_engine %||% "none (tools on PATH)",
              sprintf("`%s`", prov$params$r_container %||% "NA"),
              prov$started, prov$completed, fmt_duration(wall_ms), yn(prov$success)),
    stringsAsFactors = FALSE)
  out <- c(out, "## Run", "", md_table(run), "")

  crit <- unique(tab$criterion)
  rules <- vapply(crit, function(c) tab$rule[tab$criterion == c][1], "")
  samples <- sort(unique(tab$sample))
  vt <- data.frame(sample = samples,
                   label = vapply(samples, function(s) tab$label[tab$sample == s][1], ""),
                   stringsAsFactors = FALSE)
  for (i in seq_along(crit)) {
    vt[[sprintf("%s (%s)", crit[i], rules[i])]] <- vapply(samples, function(s) {
      r <- tab[tab$sample == s & tab$criterion == crit[i], ]
      if (!nrow(r)) return("missing")
      sprintf("%s %s", if (r$pass[1]) "PASS" else "**FAIL**", fmt_num(r$observed[1]))
    }, "")
  }
  vt[["all criteria"]] <- ifelse(per_sample[samples], "PASS", "**FAIL**")
  out <- c(out, "## Verdict per sample and criterion", "", md_table(vt), "")

  if (length(x$diagnostics)) {
    keys <- c("gene_match_fraction", "n_published_cells", "n_matched_cells", "n_starsolo_called_cells",
              "published_cells_called_by_starsolo", "per_cell_pearson_na")
    dg <- do.call(rbind, lapply(x$diagnostics, function(d) {
      m <- d$metrics
      data.frame(sample = d$sample, t(vapply(keys, function(k) {
        val <- as.numeric(m[[k]] %||% NA)
        if (startsWith(k, "n_") || k == "per_cell_pearson_na") fmt_int(val) else fmt_num(val)
      }, "")),
                 stringsAsFactors = FALSE, check.names = FALSE)
    }))
    out <- c(out, "## Diagnostics (reported, not judged)", "", md_table(dg), "")
  }

  if (length(x$qc)) {
    qt <- do.call(rbind, lapply(names(x$qc), function(s) {
      q <- x$qc[[s]]
      data.frame(sample = s, cells = q$n_cells, flagged = q$n_flagged, nmads = q$nmads,
                 "low total UMI below" = fmt_num(as.numeric(q$total_umi$lower %||% NA)),
                 "low genes below" = fmt_num(as.numeric(q$detected_genes$lower %||% NA)),
                 "high mito % above" = fmt_num(as.numeric(q$pct_mito$higher %||% NA)),
                 check.names = FALSE, stringsAsFactors = FALSE)
    }))
    out <- c(out, "## QC flags (cells flagged, none removed)", "", md_table(qt), "")
  }

  cp <- compute_by_process(x$trace)
  tot <- cp$total
  out <- c(out, sprintf("## Compute (from `%s`)", if (is.null(trace)) "pipeline_info/trace.tsv" else basename(trace)), "",
           sprintf("This run: %d tasks (%d failed, %d retried), %.3f CPU hours used, %.3f CPU hours reserved; wall time %s.",
                   tot$tasks, tot$failed, tot$retried, tot$cpu_h_used, tot$cpu_h_reserved, fmt_duration(wall_ms)),
           "CPU hours used = task run time x %cpu / 100; reserved = run time x requested CPUs.")
  if (cp$cached$tasks > 0) {
    out <- c(out, sprintf("Not counted above: %d cached tasks (`-resume`) whose compute was spent in an earlier run (%.3f CPU hours used, as recorded).",
                          cp$cached$tasks, cp$cached$cpu_h_used))
  }
  bp <- cp$by_process
  if (nrow(bp)) {
    ct <- data.frame(process = bp$process, tasks = bp$tasks, failed = bp$failed, retried = bp$retried,
                     "CPU h used" = sprintf("%.3f", bp$cpu_h_used), "CPU h reserved" = sprintf("%.3f", bp$cpu_h_reserved),
                     "run time, sum" = fmt_duration(bp$realtime_ms_sum), "longest task" = fmt_duration(bp$realtime_ms_max),
                     "peak memory (max)" = fmt_bytes(bp$peak_rss_max), "memory requested (max)" = fmt_bytes(bp$memory_max),
                     check.names = FALSE, stringsAsFactors = FALSE)
    out <- c(out, "", md_table(ct))
  }
  out <- c(out, "")

  if (length(x$versions)) {
    sv <- do.call(rbind, lapply(names(x$versions), function(proc) {
      tools <- x$versions[[proc]]
      data.frame(process = proc, tool = names(tools), version = vapply(tools, as.character, ""), stringsAsFactors = FALSE)
    }))
    out <- c(out, "## Software versions", "", md_table(sv), "")
  }
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a
