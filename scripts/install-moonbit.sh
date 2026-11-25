#!/usr/bin/env bash
# Reproducible installation of the newest stable MoonBit toolchain.
# Source: https://www.moonbitlang.com/download
#
# Usage:
#   ./scripts/install-moonbit.sh              # latest stable channel
#   ./scripts/install-moonbit.sh v0.10.19     # a pinned release identifier
#
# The script downloads the official installer, verifies the published SHA-256
# checksum of the binary archive, installs into $MOON_HOME (default ~/.moon),
# and prints `moon version --all`.
set -euo pipefail

VERSION="${1:-latest}"
MOON_HOME="${MOON_HOME:-$HOME/.moon}"
TARGET="$(uname -ms)"

case "$TARGET" in
  'Linux x86_64')   PLATFORM="linux-x86_64" ;;
  'Linux aarch64')  PLATFORM="linux-aarch64" ;;
  'Darwin arm64')   PLATFORM="darwin-aarch64" ;;
  *) echo "Unsupported platform: $TARGET" >&2; exit 1 ;;
esac

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

echo "==> fetching official installer"
curl -fsSL https://cli.moonbitlang.com/install/unix.sh -o "$WORK/install.sh"

echo "==> verifying archive checksum for channel: $VERSION"
curl -fsSL "https://cli.moonbitlang.com/binaries/$VERSION/moonbit-$PLATFORM.tar.gz.sha256" \
  -o "$WORK/moonbit.sha256"
curl -fsSL "https://cli.moonbitlang.com/binaries/$VERSION/moonbit-$PLATFORM.tar.gz" \
  -o "$WORK/moonbit-$PLATFORM.tar.gz"
(cd "$WORK" && sha256sum -c moonbit.sha256)

echo "==> installing to $MOON_HOME"
MOON_HOME="$MOON_HOME" SHELL="${SHELL:-/bin/bash}" bash "$WORK/install.sh" "$VERSION"

echo "==> installed toolchain"
"$MOON_HOME/bin/moon" version --all
"$MOON_HOME/bin/moonc" -v
