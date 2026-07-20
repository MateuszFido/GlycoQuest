#!/usr/bin/env bash
# Package a release binary into dist/ with the names used in the install docs.
#
# Usage:
#   scripts/ci/package.sh <target> <binary-path>
set -euo pipefail

TARGET=${1:?usage: package.sh <target> <binary-path>}
BIN_PATH=${2:?usage: package.sh <target> <binary-path>}

if [[ ! -f "$BIN_PATH" ]]; then
  echo "error: binary not found: $BIN_PATH" >&2
  exit 1
fi

# Prefer the git tag as-is (e.g. v0.1.0) so Release asset names match the tag.
if [[ -n "${CI_COMMIT_TAG:-}" ]]; then
  VERSION="$CI_COMMIT_TAG"
elif [[ "${GITHUB_REF_TYPE:-}" == "tag" && -n "${GITHUB_REF_NAME:-}" ]]; then
  VERSION="$GITHUB_REF_NAME"
elif [[ -n "${GLYCOQUEST_VERSION:-}" ]]; then
  VERSION="$GLYCOQUEST_VERSION"
else
  VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml | head -1)"
  SHORT_SHA="${CI_COMMIT_SHORT_SHA:-${GITHUB_SHA:-}}"
  SHORT_SHA="${SHORT_SHA:0:7}"
  if [[ -n "$SHORT_SHA" ]]; then
    VERSION="${VERSION}-${SHORT_SHA}"
  fi
fi

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
mkdir -p "$ROOT/dist"
STAGE=$(mktemp -d)
trap 'rm -rf "$STAGE"' EXIT

case "$TARGET" in
  *-pc-windows-*|*-windows-*)
    cp "$BIN_PATH" "$STAGE/glycoquest.exe"
    ARCHIVE="$ROOT/dist/glycoquest-${VERSION}-${TARGET}.zip"
    rm -f "$ARCHIVE"
    (cd "$STAGE" && zip -q "$ARCHIVE" glycoquest.exe)
    ;;
  *)
    cp "$BIN_PATH" "$STAGE/glycoquest"
    chmod +x "$STAGE/glycoquest"
    ARCHIVE="$ROOT/dist/glycoquest-${VERSION}-${TARGET}.tar.gz"
    tar -C "$STAGE" -czf "$ARCHIVE" glycoquest
    ;;
esac

printf '%s\n' "$VERSION" >"$ROOT/dist/VERSION"
ls -la "$ARCHIVE"
echo "Packed $ARCHIVE"
