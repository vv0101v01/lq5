#!/usr/bin/env bash
# Push this repository to GitHub without ever writing the token to disk.
#
#   GITHUB_PAT=github_pat_xxx ./scripts/push.sh [owner/repo] [branch]
#
# The default remote is vv0101v01/lq5 and the default branch is main.
# The token is only used in the one-shot push URL, so it never lands in
# .git/config, in a tracked file, or in the shell history of this script.
set -euo pipefail
cd "$(dirname "$0")/.."

REPO="${1:-vv0101v01/lq5}"
BRANCH="${2:-main}"
: "${GITHUB_PAT:?set GITHUB_PAT in the environment (never commit it)}"

echo "==> pushing $(git rev-parse --short HEAD) to ${REPO} (${BRANCH})"
git push "https://x-access-token:${GITHUB_PAT}@github.com/${REPO}.git" "HEAD:${BRANCH}"

# leave a token-free remote behind for later fetches
git remote remove origin 2>/dev/null || true
git remote add origin "https://github.com/${REPO}.git"
echo "==> done; 'origin' now points at https://github.com/${REPO}.git"
