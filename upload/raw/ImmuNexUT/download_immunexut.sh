#!/usr/bin/env bash

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
URL="https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-398/E-GEAD-398.processed.zip"
ZIP="$HERE/E-GEAD-398.processed.zip"
EXPECTED_MD5="2fc0a61db661d179d1083bb39ec605de"
OUTDIR="$HERE/conditional_eQTL_FDR0.05"
N_EXPECTED=28

md5of() { if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }
n_have() { ls "$OUTDIR"/*_conditional_eQTL_FDR0.05.txt 2>/dev/null | wc -l | tr -d ' '; }

mkdir -p "$OUTDIR"
if [[ "$(n_have)" == "$N_EXPECTED" ]]; then
  echo "skip  $N_EXPECTED conditional eQTL files already present"
  exit 0
fi

if [[ ! -s "$ZIP" ]] || [[ "$(md5of "$ZIP")" != "$EXPECTED_MD5" ]]; then
  curl --fail -L --connect-timeout 20 --retry 3 --retry-delay 5 -o "$ZIP.tmp" "$URL"
  if [[ "$(md5of "$ZIP.tmp")" != "$EXPECTED_MD5" ]]; then
    echo "FAIL  md5 mismatch for $(basename "$ZIP")"
    rm -f "$ZIP.tmp"
    exit 1
  fi
  mv "$ZIP.tmp" "$ZIP"
fi

TMP="$HERE/.extract_tmp"
rm -rf "$TMP"; mkdir -p "$TMP"
unzip -q -o "$ZIP" -d "$TMP"
find "$TMP" -name "*.tar" -exec tar -xf {} -C "$TMP" \;
find "$TMP" -type f -name "*_conditional_eQTL_FDR0.05.txt" -exec mv {} "$OUTDIR"/ \;
rm -rf "$TMP"

if [[ "$(n_have)" != "$N_EXPECTED" ]]; then
  echo "FAIL  expected $N_EXPECTED files in $OUTDIR, found $(n_have)"
  exit 1
fi
echo "OK    $N_EXPECTED conditional eQTL files in $(basename "$OUTDIR")"
