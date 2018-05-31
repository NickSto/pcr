#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
# get the name of the test directory
dirname=$(dirname $0)

USAGE="Usage: \$ $(basename $0) [options] [test1 [test2]]"
cmd_prefix="$dirname/../"

function main {

  gave_tests=
  verbose=true
  # Run the requested tests
  for arg in "$@"; do
    # Check for options
    #TODO: option to keep test data at end instead of removing it.
    if [[ ${arg:0:1} == '-' ]]; then
      case "$arg" in
        -h)
          echo "$USAGE" >&2
          echo "Meta tests:" >&2
          list_meta_tests >&2
          echo "Active tests:" >&2
          list_active_tests >&2
          echo "Inactive tests:" >&2
          list_inactive_tests >&2
          exit 1;;
        -p)
          cmd_prefix=;;
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
      gave_tests=true
      if [[ $verbose ]]; then
        $arg
      else
        $arg 2>/dev/null
      fi
    else
      echo "Unrecognized test \"$arg\"." >&2
    fi
  done

  # If no tests were specified in arguments, do all tests.
  if ! [[ $gave_tests ]]; then
    fail "Error: Please specify a valid test to run. Use -h option to list them."
  fi
}

function fail {
  echo "$@" >&2
  exit 1
}

function list_active_tests {
  while read declare f test; do
    if echo "$initial_declarations_plus_meta" | grep -qE "^declare -f $test\$"; then
      # Filter out regular functions and meta tests.
      continue
    elif echo "$test" | grep -qE '^_'; then
      # Filter out functions starting with an underscore.
      continue
    elif ! echo "$all_declarations_minus_inactive" | grep -qE "^declare -f $test\$"; then
      # Filter out inactive tests.
      continue
    else
      echo "  $test"
    fi
  done < <(declare -F)
}

function list_inactive_tests {
  while read declare f test; do
    if echo "$all_declarations_minus_inactive" | grep -qE "^declare -f $test\$"; then
      # Filter out regular functions, meta tests, and active tests.
      continue
    elif echo "$test" | grep -qE '^_'; then
      # Filter out functions starting with an underscore.
      continue
    else
      echo "  $test"
    fi
  done < <(declare -F)
}

function list_meta_tests {
  # Want to list these tests in this order (as long as they exist).
  for test in all active inactive; do
    if declare -F | grep -qE "^declare -f $test\$"; then
      echo "  $test"
    fi
  done
  # Then programmatically list the rest.
  while read declare f test; do
    if echo "$initial_declarations" | grep -qE "^declare -f $test\$"; then
      # Filter out regular functions.
      continue
    elif echo -e "all\nactive\ninactive" | grep -qE "^$test\$"; then
      # Filter out the fixed-order ones we already listed.
      continue
    elif echo "$initial_declarations_plus_meta" | grep -qE "^declare -f $test\$"; then
      # If it matches this list but not the last, it's a meta test.
      echo "  $test"
    fi
  done < <(declare -F)
}

# Capture a list of all functions defined before the tests, to tell which are actual functions
# and which are tests.
initial_declarations=$(declare -F)


########## Meta tests ##########

# Run all tests.
function all {
  active
  inactive
}

function active {
  for test in $(list_active_tests); do
    $test
  done
}

function inactive {
  for test in $(list_inactive_tests); do
    $test
  done
}

# Run the errstats.py-specific tests.
function errstats {
  errstats_simple
  errstats_overlap
}

function varylen {
  varylen_barcodes
  varylen_align
  varylen_consensi
}

# Run the make-consensi.py-specific tests.
function consensi_all {
  declare -a tests
  i=1
  while read declare f test; do
    if echo "$test" | grep -qE '^consensi' && [[ $test != consensi_all ]]; then
      tests[$i]=$test
      i=$((i+1))
    fi
  done < <(declare -F)
  for test in ${tests[@]}; do
    $test
  done
}

# Get the list of functions now that the meta tests have been declared.
initial_declarations_plus_meta=$(declare -F)


########## Functional tests ##########

# make-barcodes.awk
function barcodes {
  echo -e "\t${FUNCNAME[0]}:\tmake-barcodes.awk ::: families.raw_[12].fq"
  paste "$dirname/families.raw_1.fq" "$dirname/families.raw_2.fq" \
    | paste - - - - \
    | awk -f "$dirname/../make-barcodes.awk" -v TAG_LEN=12 -v INVARIANT=5 \
    | sort \
    | diff -s - "$dirname/families.sort.tsv"
}

# align-families.py
function align {
  echo -e "\t${FUNCNAME[0]}:\talign-families.py ::: families.sort.tsv:"
  "${cmd_prefix}align-families.py" --no-check-ids -q "$dirname/families.sort.tsv" \
    | diff -s - "$dirname/families.msa.tsv"
}

# align-families.py with 3 processes
function align_p3 {
  echo -e "\t${FUNCNAME[0]}:\talign-families.py -p 3 ::: families.sort.tsv:"
  "${cmd_prefix}align-families.py" --no-check-ids -q -p 3 "$dirname/families.sort.tsv" \
    | diff -s - "$dirname/families.msa.tsv"
}

# align-families.py smoke test
function align_smoke {
  echo -e "\t${FUNCNAME[0]}:\talign-families.py ::: smoke.families.tsv:"
  "${cmd_prefix}align-families.py" --no-check-ids -q "$dirname/smoke.families.tsv" \
    | diff -s - "$dirname/smoke.families.aligned.tsv"
}

# make-consensi.py defaults on toy data
function consensi {
  _consensi families.msa.tsv families.sscs_1.fa families.sscs_2.fa families.dcs_1.fa \
    families.dcs_2.fa
}

# make-consensi.py with 3 processes
function consensi_p3 {
  _consensi families.msa.tsv families.sscs_1.fa families.sscs_2.fa families.dcs_1.fa \
    families.dcs_2.fa -p 3
}

# make-consensi.py quality score consideration
function consensi_qual {
  _consensi qual.msa.tsv qual.10.sscs_1.fa qual.10.sscs_2.fa empty.txt empty.txt -q 10
  _consensi qual.msa.tsv qual.20.sscs_1.fa qual.20.sscs_2.fa empty.txt empty.txt -q 20
}

function consensi_gapqual {
  _consensi gapqual.msa.tsv gapqual.sscs_1.fa gapqual.sscs_2.fa empty.txt empty.txt -q 25
}

function consensi_consthres {
  _consensi cons.thres.msa.tsv cons.thres.0.5.sscs_1.fa cons.thres.0.5.sscs_2.fa \
          cons.thres.0.5.dcs_1.fa cons.thres.0.5.dcs_2.fa \
          --min-cons-reads 3 --cons-thres 0.5
  _consensi cons.thres.msa.tsv cons.thres.0.7.sscs_1.fa cons.thres.0.7.sscs_2.fa \
          cons.thres.0.7.dcs_1.fa cons.thres.0.7.dcs_2.fa \
          --min-cons-reads 3 --cons-thres 0.7
}

# variable-length reads
# make-barcodes.awk
function varylen_barcodes {
  echo -e "\t${FUNCNAME[0]}:\tmake-barcodes.awk ::: varylen.raw_[12].fq"
  paste "$dirname/varylen.raw_1.fq" "$dirname/varylen.raw_2.fq" \
    | paste - - - - \
    | awk -f "$dirname/../make-barcodes.awk" -v TAG_LEN=12 -v INVARIANT=5 \
    | sort \
    | diff -s - "$dirname/varylen.sort.tsv"
}

# align-families.py
function varylen_align {
  echo -e "\t${FUNCNAME[0]}:\talign-families.py ::: varylen.sort.tsv:"
  "${cmd_prefix}align-families.py" --no-check-ids -q "$dirname/varylen.sort.tsv" \
    | diff -s - "$dirname/varylen.msa.tsv"
}

# make-consensi.py
function varylen_consensi {
  _consensi varylen.msa.tsv varylen.sscs_1.fa varylen.sscs_2.fa varylen.dcs_1.fa varylen.dcs_2.fa
}

# baralign.sh
function baralign {
  echo -e "\t${FUNCNAME[0]}:\tbaralign.sh ::: correct.families.tsv:"
  "${cmd_prefix}baralign.sh" "$dirname/correct.families.tsv" "$dirname/refdir.tmp" 2>/dev/null \
    | _clean_sam | diff -s - "$dirname/correct.sam"
  rm -rf "$dirname/refdir.tmp"
}

# correct.py
function correct {
  echo -e "\t${FUNCNAME[0]}:\tcorrect.py ::: correct.sam"
  "${cmd_prefix}correct.py" --no-check-ids "$dirname/correct.families.tsv" \
      "$dirname/correct.barcodes.fa" "$dirname/correct.sam" \
    | diff -s "$dirname/correct.families.corrected.tsv" -
}

function stats_diffs {
  echo -e "\t${FUNCNAME[0]}:\tstats.py diffs ::: gaps.msa.tsv:"
  "$dirname/../utils/stats.py" diffs "$dirname/gaps.msa.tsv" \
    | diff -s - "$dirname/gaps-diffs.out.tsv"
}


function precheck {
  echo -e "\t${FUNCNAME[0]}:\tprecheck.py ::: families.raw_[12].fq"
  "$dirname/../utils/precheck.py" "$dirname/families.raw_1.fq" "$dirname/families.raw_2.fq" \
    | diff -s - "$dirname/families.precheck.tsv"
}


# All tests below here are considered inactive.
all_declarations_minus_inactive=$(declare -F)

function errstats_simple {
  echo -e "\t${FUNCNAME[0]}:\terrstats.py ::: families.msa.tsv:"
  "$dirname/../utils/errstats.py" "$dirname/families.msa.tsv" | diff -s - "$dirname/errstats.out.tsv"
  "$dirname/../utils/errstats.py" -R "$dirname/families.msa.tsv" | diff -s - "$dirname/errstats.-R.out.tsv"
  "$dirname/../utils/errstats.py" -a "$dirname/families.msa.tsv" | diff -s - "$dirname/errstats.-a.out.tsv"
  "$dirname/../utils/errstats.py" -R -a "$dirname/families.msa.tsv" | diff -s - "$dirname/errstats.-R.-a.out.tsv"
}

function errstats_overlap {
  echo -e "\t${FUNCNAME[0]}:\terrstats.py ::: families.overlap.msa.tsv"
  "$dirname/../utils/errstats.py" --dedup --min-reads 3 --bam "$dirname/families.overlap.sscs.bam" \
    "$dirname/families.overlap.msa.tsv" --overlap-stats "$dirname/overlaps.tmp.tsv" >/dev/null
  diff -s "$dirname/overlaps.tmp.tsv" "$dirname/families.overlap.overlaps.expected.tsv"
  if [[ -f "$dirname/overlaps.tmp.tsv" ]]; then
    rm "$dirname/overlaps.tmp.tsv"
  fi
}

# utility function for all make-consensi.py tests
function _consensi {
  # Read required arguments.
  input=$1
  sscs1=$2
  sscs2=$3
  dcs1=$4
  dcs2=$5
  # Read optional arguments (after the required ones).
  declare -a args
  i=6
  while [[ ${!i} ]]; do
    args[$i]=${!i}
    i=$((i+1))
  done
  echo -e "\t${FUNCNAME[1]}:\tmake-consensi.py ${args[@]} ::: $input:"
  "${cmd_prefix}make-consensi.py" ${args[@]} "$dirname/$input" \
    --sscs1 "$dirname/cons.tmp.sscs_1.fa" --sscs2 "$dirname/cons.tmp.sscs_2.fa" \
    --dcs1  "$dirname/cons.tmp.dcs_1.fa"  --dcs2  "$dirname/cons.tmp.dcs_2.fa"
  diff -s "$dirname/cons.tmp.sscs_1.fa" "$dirname/$sscs1"
  diff -s "$dirname/cons.tmp.sscs_2.fa" "$dirname/$sscs2"
  diff -s "$dirname/cons.tmp.dcs_1.fa"  "$dirname/$dcs1"
  diff -s "$dirname/cons.tmp.dcs_2.fa"  "$dirname/$dcs2"
  for file in cons.tmp.sscs_1.fa cons.tmp.sscs_2.fa cons.tmp.dcs_1.fa cons.tmp.dcs_2.fa; do
    if [[ -f "$dirname/$file" ]]; then
      rm "$dirname/$file"
    fi
  done
}

function _clean_sam {
  # Remove @PG line and XM:i: tags.
  awk -F '\t' -v OFS='\t' '
    $1 !~ /^@PG$/ {
      for (i=1; i<=NF; i++) {
        if (i == 1) {
          printf("%s", $i)
        } else if (i <= 11 || substr($i, 1, 5) != "XM:i:") {
          printf("\t%s", $i)
        }
      }
      printf("\n")
    }'
}

main "$@"
