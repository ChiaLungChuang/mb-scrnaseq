#!/usr/bin/env bash
# Tests for scripts/check_versions.sh: manifest.version, the r_container tag
# and digest, and the release tag must agree.
set -euo pipefail
# shellcheck source=tests/shell/lib.sh
source "$(dirname "$0")/lib.sh"
CV="$REPO/scripts/check_versions.sh"

D=sha256:9ce6480a75f1918f0ccd509f12207c642f4ed2bc88898309e187a7cb828367a5
IMG=ghcr.io/chialungchuang/mb-scrnaseq-r
cfg() { # cfg <manifest version> <r_container>: a minimal nextflow.config
  printf "manifest {\n    version         = '%s'\n}\nparams {\n    r_container            = '%s'\n}\n" "$1" "$2" > "$TMP/nextflow.config"
}

echo "check_versions.sh"

run bash "$CV" "$REPO/nextflow.config"
expect_rc 0 "the repository's own nextflow.config passes"

cfg 0.1.0 "$IMG:0.1.0@$D"
run bash "$CV" "$TMP/nextflow.config"
expect_rc 0 "tag + digest, tag equals manifest.version: accepted"
expect_out "(tag 0.1.0, digest $D)" "tag + digest: the tag is the part before '@', not the digest"
expect_out "pinned by digest" "tag + digest: reported as pinned"

cfg 0.1.0 "$IMG:0.1.0@$D"
{ printf "profiles {\n    local_r_image {\n        params.r_container = 'mb-scrnaseq-r:local'\n    }\n}\n"; cat "$TMP/nextflow.config"; } > "$TMP/with_profile.config"
run bash "$CV" "$TMP/with_profile.config"
expect_out "r_container=$IMG:0.1.0@$D " "a profile that overrides r_container is not read as the pipeline's image"

cfg 0.1.0 "$IMG:0.1.0"
run bash "$CV" "$TMP/nextflow.config"
expect_rc 0 "tag without digest (between a version bump and its release): accepted"
expect_out "not pinned by digest" "tag without digest: warned"

cfg 0.1.0 "$IMG:0.2.0@$D"
run bash "$CV" "$TMP/nextflow.config"
expect_fail "manifest.version (0.1.0) and the r_container tag (0.2.0) differ" "tag differs from manifest.version: refused"

cfg 0.1.0 "$IMG:0.1.0@sha256:9ce6480a"
run bash "$CV" "$TMP/nextflow.config"
expect_fail "is not sha256:<64 lowercase hex>" "truncated digest: refused"

cfg 0.1.0 "$IMG@$D"
run bash "$CV" "$TMP/nextflow.config"
expect_fail "has no tag" "digest without a tag: refused (nothing to compare with manifest.version)"

cfg 0.1.0 "localhost:5000/mb-scrnaseq-r:0.1.0"
run bash "$CV" "$TMP/nextflow.config"
expect_rc 0 "registry with a port: the port is not taken as the tag"

cfg 0.1.0 "$IMG:0.1.0@$D"
run bash "$CV" "$TMP/nextflow.config" refs/tags/v0.1.0
expect_rc 0 "release tag v0.1.0 matches: accepted"
run bash "$CV" "$TMP/nextflow.config" refs/tags/v0.2.0
expect_fail "git tag v0.2.0 does not match manifest.version 0.1.0" "release tag differs: refused"
run bash "$CV" "$TMP/nextflow.config" refs/heads/main
expect_rc 0 "a branch ref is not a release: accepted"

cfg 0.2.0 "$IMG:0.2.0"
run bash "$CV" --publish "$TMP/nextflow.config" refs/tags/v0.2.0
expect_rc 0 "publish: unpinned release accepted"
if [ "$(tail -n 1 <<< "$out")" = "$IMG:0.2.0" ]; then ok "publish: prints the image to publish as the last line"
else bad "publish: prints the image to publish as the last line" "last line: $(tail -n 1 <<< "$out")"; fi

cfg 0.1.0 "$IMG:0.1.0@$D"
run bash "$CV" --publish "$TMP/nextflow.config" refs/tags/v0.1.0
expect_fail "already pinned" "publish: refused when the version is already pinned (would move the tag)"

printf "manifest {\n    version = '0.1.0'\n}\n" > "$TMP/nextflow.config"
run bash "$CV" "$TMP/nextflow.config"
expect_fail "no r_container" "no r_container: refused"

finish
