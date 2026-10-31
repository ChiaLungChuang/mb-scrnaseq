#!/usr/bin/env bash
# Run every tests/shell/test_*.sh; exit non-zero if any fails.
set -uo pipefail
cd "$(dirname "$0")" || exit 1
failed=()
for t in test_*.sh; do
  bash "$t" || failed+=("$t")
done
if [ ${#failed[@]} -gt 0 ]; then
  echo "shell tests FAILED: ${failed[*]}"
  exit 1
fi
echo "all shell tests passed"
