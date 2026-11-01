#!/usr/bin/env bash
# Check that the pipeline version, the R image reference and a release tag
# agree. Run by CI (lint and publish-image jobs).
#
# Usage: scripts/check_versions.sh [--publish] <nextflow.config> [git_ref]
#
# Reads manifest.version and params.r_container from <nextflow.config>.
# r_container must be <image>:<tag> or <image>:<tag>@sha256:<64 hex>, and
# its tag must equal manifest.version. If git_ref is a version tag
# (refs/tags/v*), it must be v<manifest.version>. A missing digest is allowed
# (between a version bump and its release, see docs/RUNBOOK.md section 0) but
# reported.
#
# --publish (release job): also refuses an r_container that is already
# pinned by digest, because publishing again would move the tag away from
# the image the digest pins, and prints the image to publish
# (<image>:<tag>) as the last line of output.
set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "${GITHUB_ACTIONS:+::warning::}WARNING: $*" >&2; }

publish=false
if [ "${1:-}" = --publish ]; then publish=true; shift; fi
config=${1:-}
ref=${2:-}
[ -n "$config" ] || { sed -n '5p' "$0" >&2; exit 2; }
[ -r "$config" ] || die "cannot read $config"

manifest=$(grep -oP "^\s*version\s*=\s*'\K[^']+" "$config" | head -n 1 || true)
image=$(grep -oP "^\s*r_container\s*=\s*'\K[^']+" "$config" | head -n 1 || true)
[ -n "$manifest" ] || die "no manifest version in $config"
[ -n "$image" ] || die "no r_container in $config"

# Split <name>:<tag>[@<digest>]. The tag is after the last ':' of the part
# before '@', and must sit in the last path component (a registry port such
# as localhost:5000/ is not a tag).
name_tag=${image%%@*}
digest=
if [[ "$image" == *@* ]]; then digest=${image#*@}; fi
[[ "${name_tag##*/}" == *:* ]] || die "r_container '$image' has no tag; it must be <image>:<version>[@sha256:<digest>]"
tag=${name_tag##*:}
name=${name_tag%:*}
[[ "$tag" =~ ^[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}$ ]] || die "r_container tag '$tag' is not a valid image tag"
if [ -n "$digest" ] && ! [[ "$digest" =~ ^sha256:[0-9a-f]{64}$ ]]; then
  die "r_container digest '$digest' is not sha256:<64 lowercase hex>"
fi

echo "manifest.version=$manifest r_container=$image (tag $tag${digest:+, digest $digest})"
[ "$tag" = "$manifest" ] || die "manifest.version ($manifest) and the r_container tag ($tag) differ"

if [[ "$ref" == refs/tags/v* ]] && [ "${ref#refs/tags/v}" != "$manifest" ]; then
  die "git tag ${ref#refs/tags/} does not match manifest.version $manifest"
fi

if [ -n "$digest" ]; then
  echo "r_container is pinned by digest"
else
  warn "r_container is not pinned by digest: expected only between a version bump and its release"
fi

if $publish; then
  [ -z "$digest" ] || die "r_container is already pinned to $digest for $tag: publishing again would move tag $tag away from it. Bump the version instead."
  echo "$name:$tag"
fi
