#!/usr/bin/env bash

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
URL="https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v9/extra/CAGE_peaks_expression/hg38_fair%2Bnew_CAGE_peaks_phase1and2_tpm.osc.txt.gz"
OUT="$HERE/hg38_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz"
EXPECTED_SIZE=837601892

filesize() { stat -f%z "$1" 2>/dev/null || stat -c%s "$1"; }

if [[ -s "$OUT" ]] && [[ "$(filesize "$OUT")" == "$EXPECTED_SIZE" ]] && gunzip -t "$OUT" 2>/dev/null; then
  echo "skip  $(basename "$OUT") (size ok, gzip ok)"
  exit 0
fi

curl --fail -L --connect-timeout 20 --retry 3 --retry-delay 5 -o "$OUT.tmp" "$URL"
actual=$(filesize "$OUT.tmp")
if [[ "$actual" != "$EXPECTED_SIZE" ]]; then
  echo "FAIL  downloaded size $actual != expected $EXPECTED_SIZE"
  rm -f "$OUT.tmp"
  exit 1
fi
if ! gunzip -t "$OUT.tmp" 2>/dev/null; then
  echo "FAIL  bad gzip"
  rm -f "$OUT.tmp"
  exit 1
fi
mv "$OUT.tmp" "$OUT"
echo "OK    $(basename "$OUT")  $actual bytes"
