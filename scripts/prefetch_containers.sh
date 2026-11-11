#!/usr/bin/env bash
# Optional: pull every container image EPIquence uses into Nextflow's image cache before the first run.
# Files are named exactly as Nextflow names them, so Nextflow finds them and does not pull again.
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: scripts/prefetch_containers.sh [options]

  -p, --profile NAME    all (default) | low_capacity | high_capacity
  -c, --cache-dir DIR   image cache; default: NXF_SINGULARITY_CACHEDIR (or NXF_APPTAINER_CACHEDIR),
                        else cacheDir in the singularity/apptainer block of conf/env.config
  -r, --runtime NAME    apptainer | singularity (default: whichever is installed, apptainer first)
  -j, --jobs N          parallel pulls (default 2)
  -n, --dry-run         list images and target files without pulling
  -h, --help            show this help

Set APPTAINER_CACHEDIR / SINGULARITY_CACHEDIR (layer cache) and APPTAINER_TMPDIR / SINGULARITY_TMPDIR
(build space) to a large disk first; by default they live in your home directory.
EOF
}

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROFILE=all; CACHE=""; RUNTIME=""; JOBS=2; DRY=0
while [ $# -gt 0 ]; do
    case "$1" in
        -p|--profile)   PROFILE="$2"; shift 2 ;;
        -c|--cache-dir) CACHE="$2"; shift 2 ;;
        -r|--runtime)   RUNTIME="$2"; shift 2 ;;
        -j|--jobs)      JOBS="$2"; shift 2 ;;
        -n|--dry-run)   DRY=1; shift ;;
        -h|--help)      usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

case "$PROFILE" in all|low_capacity|high_capacity) ;; *) echo "Invalid --profile: $PROFILE" >&2; exit 2 ;; esac

if [ -z "$RUNTIME" ]; then
    if command -v apptainer >/dev/null 2>&1; then RUNTIME=apptainer
    elif command -v singularity >/dev/null 2>&1; then RUNTIME=singularity
    elif [ "$DRY" -eq 0 ]; then echo "Neither apptainer nor singularity found on PATH" >&2; exit 1
    else RUNTIME=apptainer; fi
fi

# Same order Nextflow uses: environment variable, then cacheDir in the runtime block of conf/env.config.
if [ -z "$CACHE" ]; then
    CACHE="${NXF_SINGULARITY_CACHEDIR:-${NXF_APPTAINER_CACHEDIR:-}}"
fi
if [ -z "$CACHE" ] && [ -f "$REPO/conf/env.config" ]; then
    CACHE="$(awk '
        /^[[:space:]]*(singularity|apptainer)[[:space:]]*\{/ { inblock = 1; next }
        inblock && /^[[:space:]]*\}/ { inblock = 0 }
        inblock && /^[[:space:]]*cacheDir[[:space:]]*=/ {
            if (match($0, /"[^"]+"/) || match($0, /\047[^\047]+\047/)) { print substr($0, RSTART + 1, RLENGTH - 2); exit }
        }' "$REPO/conf/env.config")"
fi
if [ -z "$CACHE" ]; then
    echo "No image cache found: set singularity.cacheDir in conf/env.config, NXF_SINGULARITY_CACHEDIR, or use --cache-dir" >&2
    exit 1
fi

# Images from nextflow.config; each profile skips the images it never uses.
mapfile -t IMAGES < <(grep -oE "docker://[A-Za-z0-9._/:@+-]+" "$REPO/nextflow.config" | sort -u)
case "$PROFILE" in
    low_capacity)  mapfile -t IMAGES < <(printf '%s\n' "${IMAGES[@]}" | grep -vE "ecogenomic/gtdbtk:|staphb/bakta:1\.11\.4$") ;;
    high_capacity) mapfile -t IMAGES < <(printf '%s\n' "${IMAGES[@]}" | grep -vE "staphb/bakta:.*-light$") ;;
esac
[ "${#IMAGES[@]}" -gt 0 ] || { echo "No images found in $REPO/nextflow.config" >&2; exit 1; }

image_file() {   # Nextflow's naming: drop the scheme, replace ':' and '/' with '-', add .img
    local name="${1#*://}"
    name="${name//:/-}"; name="${name//\//-}"
    printf '%s/%s.img\n' "$CACHE" "$name"
}

echo "Runtime:    $RUNTIME"
echo "Image cache: $CACHE"
echo "Profile:    $PROFILE (${#IMAGES[@]} images)"
if [ -z "${APPTAINER_CACHEDIR:-}${SINGULARITY_CACHEDIR:-}" ]; then
    echo "Note: APPTAINER_CACHEDIR / SINGULARITY_CACHEDIR not set; layers are cached under $HOME"
fi

if [ "$DRY" -eq 1 ]; then
    for img in "${IMAGES[@]}"; do
        f="$(image_file "$img")"
        if [ -s "$f" ]; then state=present; else state=missing; fi
        printf '%-8s %s -> %s\n' "$state" "$img" "$f"
    done
    exit 0
fi

mkdir -p "$CACHE"
LOG="$(mktemp)"
trap 'rm -f "$LOG"' EXIT

pull_one() {
    local img="$1" f tmp
    f="$(image_file "$img")"
    if [ -s "$f" ]; then echo "SKIP $img"; return 0; fi
    tmp="$f.pulling.$$.$RANDOM"
    if "$RUNTIME" pull --name "$tmp" "$img" >/dev/null 2>"$tmp.err"; then
        mv -f "$tmp" "$f"; rm -f "$tmp.err"; echo "OK   $img"
    else
        echo "FAIL $img ($(tail -1 "$tmp.err" 2>/dev/null))"; rm -f "$tmp" "$tmp.err"
    fi
}
export -f pull_one image_file
export RUNTIME CACHE

# shellcheck disable=SC2016  # $1 is expanded by the child bash, as intended
printf '%s\n' "${IMAGES[@]}" | xargs -P "$JOBS" -I{} bash -c 'pull_one "$1"' _ {} | tee "$LOG"

ok=$(grep -c '^OK' "$LOG" || true); skip=$(grep -c '^SKIP' "$LOG" || true); fail=$(grep -c '^FAIL' "$LOG" || true)
echo "Done: $ok pulled, $skip already present, $fail failed"
df -h "$CACHE" | tail -1 | awk '{print "Cache disk: " $4 " free of " $2}'
[ "$fail" -eq 0 ]
