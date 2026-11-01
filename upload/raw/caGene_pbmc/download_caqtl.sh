#!/usr/bin/env bash

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="https://zenodo.org/api/records/14652966/files"

md5of() { if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }

while read -r name md5; do
  out="$HERE/$name"
  if [[ -s "$out" ]] && [[ "$(md5of "$out")" == "$md5" ]]; then
    echo "skip  $name (md5 ok)"
    continue
  fi
  curl --fail -L --connect-timeout 20 --retry 3 --retry-delay 5 -o "$out.tmp" "$BASE/$name/content"
  if [[ "$(md5of "$out.tmp")" != "$md5" ]]; then
    echo "FAIL  md5 mismatch for $name"
    rm -f "$out.tmp"
    exit 1
  fi
  mv "$out.tmp" "$out"
  echo "OK    $name"
done <<'LIST'
full_summstats_B_250kb.txt.gz c332d0e07a125887ecb3aa73d13902e3
full_summstats_CD4.T_250kb.txt.gz f551b6e542cb9944f772727e2b5eaf8d
full_summstats_CD8.T_250kb.txt.gz 475180bb7cfa776e7d9581135e984bfb
full_summstats_DC_250kb.txt.gz 499d3934bcfe127ab5cd8adbd596ad65
full_summstats_other.T_250kb.txt.gz b8c841c67652a01582341ea7da5fa06b
full_summstats_Mono_250kb.txt.gz 3c6211186cd577decafbf7abdb75e3f1
full_summstats_NK_250kb.txt.gz c50f394da7d54b72325d3f32a7f78de5
LIST
