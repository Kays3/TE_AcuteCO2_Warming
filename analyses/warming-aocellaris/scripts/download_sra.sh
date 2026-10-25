#!/usr/bin/env bash
# Download the 48 poly(A) runs as .sra, validate, checksum and count spots.
# FASTQ extraction is deferred to the processing stage (disk: ~700 GB uncompressed).
# Usage: P=4 download_sra.sh design.tsv runinfo.csv   (run from ~/jim_te_warm, under taskset -c 0-20)
# Fails a run if its spot count differs from runinfo `spots`; records whether the file is .sra or .sralite.
set -uo pipefail
E=~/jim_te_warm/mamba/envs/jim_te_warm/bin
OUT=~/jim_te_warm/sra; LOG=~/jim_te_warm/logs/download; mkdir -p "$OUT" "$LOG"
one() {
  r=$1
  [ -s "$LOG/$r.done" ] && return 0
  /usr/bin/time -v -o "$LOG/$r.prefetch.time" "$E/prefetch" --max-size 100G -O "$OUT" "$r" > "$LOG/$r.prefetch.log" 2>&1 || { echo "$r prefetch_fail" >> "$LOG/failed.txt"; return 1; }
  f=$(ls "$OUT/$r"/*.sra* 2>/dev/null | head -1)
  /usr/bin/time -v -o "$LOG/$r.validate.time" "$E/vdb-validate" "$f" > "$LOG/$r.validate.log" 2>&1 || { echo "$r validate_fail" >> "$LOG/failed.txt"; return 1; }
  sha=$(sha256sum "$f" | cut -d' ' -f1); md5=$(md5sum "$f" | cut -d' ' -f1)
  spots=$("$E/sra-stat" --quick --xml "$f" 2>/dev/null | grep -o 'spot_count="[0-9]*"' | head -1 | grep -o '[0-9]*')
  want=$(awk -F, -v r="$r" 'NR==1{for(i=1;i<=NF;i++)if($i=="spots")c=i} $1==r{print $c}' "$RUNINFO")
  [ "$spots" = "$want" ] || { echo "$r spots_mismatch got=$spots want=$want" >> "$LOG/failed.txt"; return 1; }
  kind=${f##*.}
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$r" "$kind" "$(stat -c %s "$f")" "$sha" "$md5" "$spots" > "$LOG/$r.done"
}
RUNINFO=$(readlink -f "$2")
export -f one; export E OUT LOG RUNINFO
cut -f1 "$1" | tail -n +2 | xargs -P "${P:-2}" -I{} bash -c 'one {}'
cat "$LOG"/*.done | sort > ~/jim_te_warm/sra_checksums.tsv
