#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

Usage="Usage: \$ $(basename $0) [options] families.msa.tsv dcs1.fa dcs2.fa [sscs1.fa sscs2.fa]"

function main {

  # Parse options.
  # Only allow a whitelist of arguments to be passed through to dunovo.py.
  add_next=
  o=0
  opts=()
  p=0
  positionals=()
  for opt in "$@"; do
    add=
    if [[ $add_next ]]; then
      # The previous option took an argument.
      if [[ ${opt:0:1} == '-' ]]; then
        fail "Error: Option \"${opts[$((o-1))]}\" requires an argument (saw \"$opt\" instead)."
      else
        add=true
        add_next=
      fi
    else
      # nargs is how many "args" should be added for a particular option.
      # 1 means only that option itself.
      # 2 means the option and the following one (its argument).
      nargs=0
      case $opt in
        -h) fail "$Usage";;
        --galaxy) nargs=1;;
        --phone-home) nargs=1;;
        --test) nargs=1;;
        -r) nargs=2;;
        -q) nargs=2;;
        -F) nargs=2;;
      esac
      if [[ $nargs == 1 ]]; then
        add=true
      elif [[ $nargs == 2 ]]; then
        add=true
        add_next=true
      elif [[ ${opt:0:1} == '-' ]]; then
        fail "Error: Unrecognized option \"$opt\"."
      else
        # If it's an argument that doesn't start with "-" and it's not the argument of a previous
        # option, assume it's one of the positional arguments for this script.
        if [[ $p -gt 4 ]]; then
          fail "Error: Too many positional arguments. Failed on \"$opt\"."
        else
          positionals[$p]="$opt"
          p=$((p+1))
        fi
      fi
    fi
    if [[ $add ]]; then
      opts[$o]="$opt"
      o=$((o+1))
    fi
  done

  # Get positionals.
  if [[ ${#positionals[@]} -ge 3 ]]; then
    alignments="${positionals[0]}"
    dcs1="${positionals[1]}"
    dcs2="${positionals[2]}"
  else
    fail "Error: Too few positional arguments."
  fi
  if [[ ${#positionals[@]} -ge 5 ]]; then
    sscs1="${positionals[3]}"
    sscs2="${positionals[4]}"
    sscs_args='--sscs-file sscs.fa'
  else
    sscs1=
    sscs2=
    sscs_args=
  fi

  # Find the actual directory this file resides in (resolving links).
  if readlink -f dummy >/dev/null 2>/dev/null; then
    script_path=$(readlink -f "${BASH_SOURCE[0]}")
  else
    script_path=$(perl -MCwd -le 'print Cwd::abs_path(shift)' "${BASH_SOURCE[0]}")
  fi
  script_dir=$(dirname "$script_path")

  # Locate outconv.py.
  # In $script_dir/utils in a normal installation, $script_dir in a Conda installation.
  outconv_script="$script_dir/utils/outconv.py"
  if ! [[ -f "$outconv_script" ]]; then
    outconv_script="$script_dir/outconv.py"
  fi

  python2 "$script_dir/dunovo.py" "${opts[@]}" "$alignments" $sscs_args > duplex.fa
  python2 "$outconv_script" duplex.fa -1 "$dcs1" -2 "$dcs2"

  if [[ $sscs1 ]] && [[ $sscs2 ]]; then
    python2 "$outconv_script" sscs.fa -1 "$sscs1" -2 "$sscs2"
  fi
}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
