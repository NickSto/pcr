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
-t: The length of the barcode portion of each read. Default: $TagLenDefault
-i: The length of the invariant (ligation) portion of each read. Default: $InvariantDefault"

function main {

  # Read arguments.
  if [[ "$#" -lt 2 ]] || [[ "$1" == '--help' ]]; then
    fail "$Usage"
  fi
  taglen=$TagLenDefault
  invariant=$InvariantDefault
  while getopts ":t:i:h" opt; do
  case "$opt" in
      t) taglen=$OPTARG;;
      i) invariant=$OPTARG;;
      h) fail "$USAGE";;
    esac
  done
  # Get positionals.
  fastq1="${@:$OPTIND:1}"
  fastq2="${@:$OPTIND+1:1}"

  if ! [[ "$fastq1" ]] || ! [[ "$fastq2" ]]; then
    fail "Error: Must provide two input fastq files."
  fi

  script_dir=$(get_script_dir)

  paste "$fastq1" "$fastq2" \
    | paste - - - - \
    | awk -f "$script_dir/make-barcodes.awk" -v TAG_LEN=$taglen -v INVARIANT=$invariant \
    | sort

}

function get_script_dir {
  # Find the actual directory this file resides in (resolving links).
  if readlink -f dummy >/dev/null 2>/dev/null; then
    script_path=$(readlink -f "${BASH_SOURCE[0]}")
  else
    script_path=$(perl -MCwd -le 'print Cwd::abs_path(shift)' "${BASH_SOURCE[0]}")
  fi
  dirname "$script_path"
}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
