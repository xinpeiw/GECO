#!/usr/bin/env bash

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="$HERE/re2g_files_manifest.tsv"
OUTDIR="$HERE/encode_re2g_multi"

if [[ ! -f "$MANIFEST" ]]; then
  echo "Manifest not found: $MANIFEST"
  exit 1
fi

mkdir -p "$OUTDIR"

_download_one() {
  IFS="|" read -r url out expected_size <<< "$1"
  out_path="$OUTDIR/$out"

  filesize() { stat -f%z "$1" 2>/dev/null || stat -c%s "$1"; }

  if [[ -s "$out_path" ]]; then
    actual=$(filesize "$out_path")
    if [[ "$actual" == "$expected_size" ]] && gunzip -t "$out_path" 2>/dev/null; then
      echo "skip  $out (size=$actual ok, gzip ok)"
      return 0
    else
      echo "stale $out (size=$actual vs expected=$expected_size or bad gzip) - re-downloading"
      rm -f "$out_path"
    fi
  fi

  if ! curl --fail -sSL --connect-timeout 20 --max-time 300 \
      --retry 3 --retry-delay 5 --retry-all-errors \
      -o "$out_path.tmp" "$url"; then
    rm -f "$out_path.tmp"
    echo "FAIL  $out (curl error)"
    return 1
  fi
  actual=$(filesize "$out_path.tmp")
  if [[ "$actual" != "$expected_size" ]]; then
    echo "FAIL  $out (downloaded size $actual != expected $expected_size)"
    rm -f "$out_path.tmp"
    return 1
  fi
  if ! gunzip -t "$out_path.tmp" 2>/dev/null; then
    echo "FAIL  $out (bad gzip)"
    rm -f "$out_path.tmp"
    return 1
  fi
  mv "$out_path.tmp" "$out_path"
  echo "OK    $out  $actual bytes"
}
export -f _download_one
export OUTDIR

awk -F'\t' 'NR>1 {printf "%s|%s_%s_%s.bed.gz|%s\n", $5, $1, $2, $3, $4}' "$MANIFEST" \
  | xargs -P 8 -I{} bash -c '_download_one "$@"' _ {}

echo
echo "Downloaded files: $(ls "$OUTDIR"/*.bed.gz 2>/dev/null | wc -l)"
echo "Total size:      $(du -sh "$OUTDIR" | awk '{print $1}')"
