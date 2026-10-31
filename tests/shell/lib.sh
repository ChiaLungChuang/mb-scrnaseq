# shellcheck shell=bash
# Minimal assertion helpers for the shell tests (sourced, not run).
# Each test file sources this, runs commands with `run`, checks them with
# the expect_* helpers and ends with `finish`, which sets the exit status.

pass=0
fail=0
# shellcheck disable=SC2034  # used by the test files that source this
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

ok()  { pass=$((pass + 1)); printf '  ok    %s\n' "$1"; }
bad() { fail=$((fail + 1)); printf '  FAIL  %s\n' "$1"; [ -z "${2:-}" ] || printf '        %s\n' "$2"; }

# run CMD...: capture combined output in $out and the exit status in $rc.
run() {
  set +e
  out=$("$@" 2>&1)
  rc=$?
  set -e
}

expect_rc()   { if [ "$rc" -eq "$1" ]; then ok "$2"; else bad "$2" "exit $rc, expected $1. Output: $out"; fi; }
expect_out()  { if grep -qF -- "$1" <<< "$out"; then ok "$2"; else bad "$2" "output lacks '$1'. Output: $out"; fi; }
expect_fail() { # expect_fail <message fragment> <name>
  if [ "$rc" -ne 0 ] && grep -qF -- "$1" <<< "$out"; then ok "$2"
  else bad "$2" "exit $rc; expected failure mentioning '$1'. Output: $out"; fi
}
expect_file_eq() { # expect_file_eq <file> <expected content> <name>
  local got
  got=$(cat "$1" 2>/dev/null || echo "<missing $1>")
  if [ "$got" = "$2" ]; then ok "$3"; else bad "$3" "got: $(printf '%q' "$got") expected: $(printf '%q' "$2")"; fi
}

finish() {
  echo "  $pass passed, $fail failed"
  [ "$fail" -eq 0 ]
}
