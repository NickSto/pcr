#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
# get the name of the test directory
dirname=$(dirname $0)

USAGE="Usage: \$ $(basename $0) [options] [test1 [test2]]"


function main {

  do_all=true
  verbose=true
  # Run the requested tests
  for arg in "$@"; do
    # Check for options
    #TODO: option to keep test data at end instead of removing it.
    if [[ ${arg:0:1} == '-' ]]; then
      case "$arg" in
        -h)
          echo "$USAGE" >&2
          echo "Currently valid tests:" >&2
          list_tests >&2
          exit 1;;
        -q)
          verbose='';;
        -v)
          verbose=true;;
        *)
          echo "Unrecognized option \"$arg\"." >&2;;
      esac
      continue
    fi
    # Execute valid tests (if they're existing functions).
    if [[ $(type -t $arg) == function ]]; then
      do_all=''
      if [[ $verbose ]]; then
        $arg
      else
        $arg 2>/dev/null
      fi
    else
      echo "Unrecognized test \"$arg\"." >&2
      do_all=''
    fi
  done

  # If no tests were specified in arguments, do all tests.
  if [[ $do_all ]]; then
    if [[ $verbose ]]; then
      all
    else
      all 2>/dev/null
    fi
  fi
}

function fail {
  echo "$@" >&2
  exit 1
}

function list_tests {
  while read declare f test; do
    # Filter out functions that aren't tests.
    if echo "$initial_declarations" | grep -qF 'declare -f '"$test"; then
      continue
    else
      echo "$test"
    fi
  done < <(declare -F)
}

# Capture a list of all functions defined before the tests, to tell which are actual functions
# and which are tests.
initial_declarations=$(declare -F)

########## Functional tests ##########

# Do all tests.
function all {
  barcodes
  align
  align_p3
  duplex
  duplex_qual
  stats_diffs
  errstats
}

# make-barcodes.awk
function barcodes {
  echo -e "\tmake-barcodes.awk ::: families.raw_[12].fq"
  paste "$dirname/families.raw_1.fq" "$dirname/families.raw_2.fq" | paste - - - - \
    | awk -f "$dirname/../make-barcodes.awk" -v TAG_LEN=12 -v INVARIANT=5 | sort \
    | diff -s - "$dirname/families.sort.tsv"
}

# align_families.py
function align {
  echo -e "\talign_families.py ::: families.sort.tsv:"
  python "$dirname/../align_families.py" "$dirname/families.sort.tsv" | diff -s - "$dirname/families.msa.tsv"
}

# align_families.py with 3 processes
function align_p3 {
  echo -e "\talign_families.py ::: families.sort.tsv:"
  python "$dirname/../align_families.py" -p 3 "$dirname/families.sort.tsv" | diff -s - "$dirname/families.msa.tsv"
}

# dunovo.py defaults on toy data
function duplex {
  echo -e "\tdunovo.py ::: families.msa.tsv:"
  python "$dirname/../dunovo.py" "$dirname/families.msa.tsv" | diff -s - "$dirname/families.cons.fa"
  python "$dirname/../dunovo.py" --incl-sscs "$dirname/families.msa.tsv" | diff -s - "$dirname/families.cons.incl-sscs.fa"
}

# dunovo.py quality score consideration
function duplex_qual {
  echo -e "\tdunovo.py ::: qual.msa.tsv:"
  python "$dirname/../dunovo.py" --incl-sscs -q 20 "$dirname/qual.msa.tsv" | diff -s - "$dirname/qual.cons.20.fa"
  python "$dirname/../dunovo.py" --incl-sscs -q 10 "$dirname/qual.msa.tsv" | diff -s - "$dirname/qual.cons.10.fa"
}

function duplex_gapqual {
  echo -e "\tdunovo.py ::: gapqual.msa.tsv:"
  python "$dirname/../dunovo.py" --incl-sscs -q 25 "$dirname/gapqual.msa.tsv" | diff -s - "$dirname/gapqual.cons.fa"
}

function stats_diffs {
  echo -e "\tstats.py diffs ::: gaps.msa.tsv:"
  python "$dirname/../utils/stats.py" diffs "$dirname/gaps.msa.tsv" | diff -s - "$dirname/gaps-diffs.out.tsv"
}

function errstats {
  echo -e "\terrstats.py ::: families.msa.tsv:"
  python3 "$dirname/../utils/errstats.py" "$dirname/families.msa.tsv" | diff -s - "$dirname/errstats.out.tsv"
  python3 "$dirname/../utils/errstats.py" -r "$dirname/families.msa.tsv" | diff -s - "$dirname/errstats.-r.out.tsv"
  python3 "$dirname/../utils/errstats.py" -a "$dirname/families.msa.tsv" | diff -s - "$dirname/errstats.-a.out.tsv"
  python3 "$dirname/../utils/errstats.py" -r -a "$dirname/families.msa.tsv" | diff -s - "$dirname/errstats.-r.-a.out.tsv"
}

main "$@"
