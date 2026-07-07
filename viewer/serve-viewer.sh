#!/bin/sh
# Serve the GlycoQuest viewer (required: browsers block fetch() on file:// URLs).
set -euo pipefail
cd "$(dirname "$0")"
PORT="${1:-8080}"
echo "GlycoQuest viewer: http://localhost:${PORT}/"
echo "Press Ctrl+C to stop."
exec python3 -m http.server "$PORT"
