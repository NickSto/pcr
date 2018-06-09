#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

TagLenDefault=12
InvariantDefault=5
Usage="Usage: \$ $(basename $0) [-t tag_len] [-i invariant_len] reads_1.fq reads_2.fq > families.tsv
Read raw duplex sequencing reads, extract their barcodes, and group them by barcode.
Input fastq's can be gzipped.
-t: The length of the barcode portion of each read. Default: $TagLenDefault
-i: The length of the invariant (ligation) portion of each read. Default: $InvariantDefault"

function main {

  # Read arguments.
  if [[ "$#" -lt 1 ]] || [[ "$1" == '--help' ]]; then
    fail "$Usage"
  fi
  if [[ "$#" -ge 1 ]] && [[ "$1" == '--version' ]]; then
    version
    return
  fi
  taglen=$TagLenDefault
  invariant=$InvariantDefault
  while getopts ":t:i:vh" opt; do
    case "$opt" in
      t) taglen=$OPTARG;;
      i) invariant=$OPTARG;;
      h) fail "$USAGE";;
      v) version && return;;
    esac
  done
  # Get positionals.
  fastq1="${@:$OPTIND:1}"
  fastq2="${@:$OPTIND+1:1}"

  if ! [[ "$fastq1" ]] || ! [[ "$fastq2" ]]; then
    fail "$Usage
Error: Must provide two input fastq files."
  fi

  script_dir=$(get_script_dir)

  # Are the input files gzipped?
  fq1_is_gzip=$(is_gzip "$fastq1")
  fq2_is_gzip=$(is_gzip "$fastq2")

  # The actual command pipeline that creates the families.
  if [[ "$fq1_is_gzip" ]] && [[ "$fq2_is_gzip" ]]; then
    paste <(gunzip -c "$fastq1") <(gunzip -c "$fastq2") \
      | paste - - - - \
      | awk -f "$script_dir/make-barcodes.awk" -v TAG_LEN=$taglen -v INVARIANT=$invariant \
      | sort
  elif ! [[ "$fq1_is_gzip" ]] && ! [[ "$fq2_is_gzip" ]]; then
    paste "$fastq1" "$fastq2" \
      | paste - - - - \
      | awk -f "$script_dir/make-barcodes.awk" -v TAG_LEN=$taglen -v INVARIANT=$invariant \
      | sort
  else
    fail "Error: Both fastq's must be either gzipped or not. No mixing is allowed."
  fi
}

function is_gzip {
  path="$1"
  ext2="${path:$((${#path}-3))}"
  ext5="${path:$((${#path}-6))}"
  if [[ "$ext2" == '.gz' ]]; then
    echo yes
  elif [[ "$ext2" == '.fq' ]] || [[ "$ext5" == '.fastq' ]]; then
    return
  else
    filetype=$(file -ib "$path")
    if [[ "${filetype:0:16}" == 'application/gzip' ]]; then
      echo yes
    fi
  fi
}

function version {
  script_dir=$(get_script_dir)
  "$script_dir/utillib/version.py" --config-path "$script_dir/VERSION" --repo-dir "$script_dir"
}

function get_script_dir {
  # Find the actual directory this file resides in (resolving links).
  if readlink -f dummy >/dev/null 2>/dev/null; then
    script_path=$(readlink -f "${BASH_SOURCE[0]}")
  else
    # readlink -f doesn't work on BSD systems.
    script_path=$(perl -MCwd -le 'print Cwd::abs_path(shift)' "${BASH_SOURCE[0]}")
  fi
  dirname "$script_path"
}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
