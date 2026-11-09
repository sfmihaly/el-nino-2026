#!/bin/bash
# Run the monitor with the private virtual environment in this folder; log to output/monitor.log
DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$DIR/output"
"$DIR/.venv/bin/python" "$DIR/monitor.py" "$@" >> "$DIR/output/monitor.log" 2>&1
