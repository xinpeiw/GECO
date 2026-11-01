#!/usr/bin/env bash

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
URL="https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
OUT="$HERE/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
EXPECTED_SIZE=322775879

filesize() { wc -c < "$1" | tr -d ' '; }

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
