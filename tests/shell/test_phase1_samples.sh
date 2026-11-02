#!/usr/bin/env bash
# The phase 1 sample record (docs/results/phase1/samples.tsv) must agree
# with the barcode-verified run table it was taken from
# (docs/proposals/gse155446_runs.tsv), so a hand edit cannot put a wrong
# label, chemistry or BAM on a run.
# shellcheck disable=SC2016  # the awk programs are single-quoted on purpose
set -euo pipefail
# shellcheck source=tests/shell/lib.sh
source "$(dirname "$0")/lib.sh"
SAMPLES=${SAMPLES:-$REPO/docs/results/phase1/samples.tsv}
RUNS=$REPO/docs/proposals/gse155446_runs.tsv

# check_samples <samples.tsv>: prints one line per problem, nothing if consistent.
check_samples() {
  awk -F'\t' -v RUNS="$RUNS" '
    BEGIN { split("sample=gsm published_label=published_label chemistry=chemistry subgroup=subgroup bam_file=bam_file bam_md5=bam_md5 bam_gb=bam_gb", pairs, " ")
            for (i in pairs) { split(pairs[i], kv, "="); want[kv[1]] = kv[2] } }
    FNR == 1 { for (i = 1; i <= NF; i++) col[FILENAME, $i] = i; next }
    FILENAME == RUNS {
      s = $col[RUNS, "srr"]
      for (k in want) run[s, k] = $col[RUNS, want[k]]
      t = $col[RUNS, "geo_title"]; sub(/-(repeat|recurrence)$/, "", t)   # 934-repeat -> patient 934
      patient[$col[RUNS, "published_label"]] = t
      next
    }
    {
      n++; s = $col[FILENAME, "srr"]; lab = $col[FILENAME, "published_label"]
      if (!((s, "published_label") in run)) { print "unknown run " s; next }
      for (k in want) if ($col[FILENAME, k] != run[s, k]) print s ": " k " is " $col[FILENAME, k] ", the run table says " run[s, k]
      if (seen_label[lab]++) print "label " lab " appears twice"
      p = patient[lab]
      if (seen_patient[p]++) print "patient " p " appears twice (label " lab "): one library per tumor"
      order[$col[FILENAME, "run_order"]] = lab
      sub_seen[$col[FILENAME, "subgroup"]] = 1
    }
    END {
      if (n != 8) print "expected 8 samples, found " n
      for (i = 1; i <= 8; i++) if (!(i in order)) print "run_order " i " missing"
      if (order[1] != "1355") print "run_order 1 is " order[1] ", not the pilot 1355 (the only v2 run)"
      split("WNT|SHH|Group 3|Group 4", g, "|")
      for (k in g) if (!(g[k] in sub_seen)) print "no " g[k] " sample"
    }
  ' "$RUNS" "$1"
}

echo "phase 1 samples"

problems=$(check_samples "$SAMPLES")
if [ -z "$problems" ]; then ok "samples.tsv agrees with the run table: 8 tumors, one library each, pilot 1355 first"
else bad "samples.tsv agrees with the run table: 8 tumors, one library each, pilot 1355 first" "$problems"; fi

# The checker itself must catch what it is for (each case edits a copy).
mutate() { # mutate <name> <awk program applied to a copy of samples.tsv> <expected problem fragment>
  awk -F'\t' -v OFS='\t' "$2" "$SAMPLES" > "$TMP/mutant.tsv"
  local p; p=$(check_samples "$TMP/mutant.tsv")
  if grep -qF -- "$3" <<< "$p"; then ok "$1: caught"; else bad "$1: caught" "expected '$3'; got: ${p:-<nothing>}"; fi
}
mutate "a wrong label on a run (996's BAM labelled 934)" 'NR>1 && $3=="996" {$3="934"} {print}' "published_label is 934"
mutate "a wrong chemistry" 'NR>1 && $3=="1355" {$4="v3"} {print}' "chemistry is v3"
mutate "a wrong BAM checksum" 'NR>1 && $3=="831" {$8="00000000000000000000000000000000"} {print}' "bam_md5 is 0000"
mutate "both libraries of patient 934" 'NR>1 && $3=="1195" {$1="GSM4704529"; $2="SRR12354570"; $3="934"; $7="2_934_possorted_genome_bam.bam"; $8="35f050a7b47c3da7f40726829dd4017b"; $9="10.6"} NR>1 && $3=="996" {$1="GSM4704530"; $2="SRR12354571"; $3="943"; $7="1_943_possorted_genome_bam.bam"; $8="feee5c291094340c49e3d8f633d646f7"; $9="23.1"} {print}' "patient 934 appears twice"
mutate "the pilot not first" 'NR>1 && $6=="1" {$6="9"} NR>1 && $6=="2" {$6="1"} {print}' "not the pilot 1355"
mutate "a sample dropped" 'NR!=9 {print}' "expected 8 samples"

finish
