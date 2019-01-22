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

# Active tests go here.

# All tests below here are considered inactive.
all_declarations_minus_inactive=$(declare -F)

function errstats_simple {
  echo -e "\t${FUNCNAME[0]}:\terrstats.py ::: families.msa.tsv:"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" utils/errstats.py); then return 1; fi
  "${local_prefix}errstats.py" "$dirname/families.msa.tsv" | diff -s - "$dirname/errstats.out.tsv"
  "${local_prefix}errstats.py" -R "$dirname/families.msa.tsv" | diff -s - "$dirname/errstats.-R.out.tsv"
  "${local_prefix}errstats.py" -a "$dirname/families.msa.tsv" | diff -s - "$dirname/errstats.-a.out.tsv"
  "${local_prefix}errstats.py" -R -a "$dirname/families.msa.tsv" | diff -s - "$dirname/errstats.-R.-a.out.tsv"
}

function errstats_overlap {
  echo -e "\t${FUNCNAME[0]}:\terrstats.py ::: families.overlap.msa.tsv"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" utils/errstats.py); then return 1; fi
  "${local_prefix}errstats.py" --dedup --min-reads 3 --bam "$dirname/families.overlap.sscs.bam" \
    "$dirname/families.overlap.msa.tsv" --overlap-stats "$dirname/overlaps.tmp.tsv" >/dev/null
  diff -s "$dirname/overlaps.tmp.tsv" "$dirname/families.overlap.overlaps.expected.tsv"
  if [[ -f "$dirname/overlaps.tmp.tsv" ]]; then
    rm "$dirname/overlaps.tmp.tsv"
  fi
}

function _get_local_prefix {
  local cmd_prefix="$1"
  postfix="$2"
  base=$(basename "$postfix")
  if [[ "$cmd_prefix" ]]; then
    path="${cmd_prefix}$postfix"
    local_prefix=$(dirname "${cmd_prefix}$postfix")/
  else
    # If $cmd_prefix is blank, the user wants to try to execute scripts via the $PATH.
    # Use the basename of the $postfix, removing any directories before the actual script.
    path="$base"
    local_prefix=
  fi
  if which "$path" >/dev/null 2>/dev/null; then
    echo "$local_prefix"
  elif [[ -f "$path" ]]; then
    echo "$local_prefix"
  else
    echo -e "\e[31mError: $base missing!\e[m Searched for: \"$path\"" >&2
    return 1
  fi
}

main "$@"
