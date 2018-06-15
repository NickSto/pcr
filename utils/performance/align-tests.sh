#!/usr/bin/env bash
if [ "x$BASH" = x ] || [ ! "$BASH_VERSINFO" ] || [ "$BASH_VERSINFO" -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

MaxParallelDefault=4
OldPathDefault=~/code/dunovo.old
NewPathDefault=~/code/dunovo.new
Usage="Usage: \$ $(basename "$0") [options] infile"

function main {

  # Get arguments.
  max_parallel="$MaxParallelDefault"
  new_path="$NewPathDefault"
  old_path="$OldPathDefault"
  new_version=
  old_version=
  slurm=
  debug=
  while getopts "m:p:P:v:V:i:sDh" opt; do
    case "$opt" in
      m) max_parallel="$OPTARG";;
      p) new_path="$OPTARG";;
      P) old_path="$OPTARG";;
      v) new_version="$OPTARG";;
      V) old_version="$OPTARG";;
      s) slurm='-s';;
      D) debug=true;;
      h) fail "$Usage";;
      ?) fail "$Usage";;
    esac
  done
  infile="${@:$OPTIND:1}"

  if ! [[ "$infile" ]]; then
    fail "$Usage"
  fi
  if ! [[ -f "$infile" ]]; then
    fail "Error: Could not find input file $infile"
  fi

  # Code paths and versions.

  for path in "$old_path" "$new_path"; do
    if ! [[ -d "$path" ]]; then
      fail "Error: Could not find $path."
    fi
  done

  if ! [[ "$old_version" ]]; then
    old_version=$(get_version "$old_path")
  fi
  if ! [[ "$new_version" ]]; then
    new_version=$(get_version "$new_path")
  fi

  old_commit=$(get_commit "$old_path")
  new_commit=$(get_commit "$new_path")

  echo "Found old version $old_version ($old_commit) and new version $new_version ($new_commit)" >&2

  script_dir=$(get_script_dir)
  measure_cmd="$script_dir/measure-cmd.sh"
  if ! [[ -x "$measure_cmd" ]]; then
    fail "Error: Could not find or execute measure-cmd.sh: $measure_cmd"
  fi

  unfinished=
  stats_files=
  for workers in 32 16 8 4 2 1; do
    for age in old new; do
      if [[ "$age" == old ]]; then
        path="$old_path"
      else
        path="$new_path"
      fi
      for script_name in 'align_families.py' 'align-families.py' ''; do
        if [[ -f "$path/$script_name" ]]; then
          break
        fi
      done
      for algorithm in mafft kalign; do
        if [[ "$age" == old ]]; then
          if [[ "$algorithm" == kalign ]]; then
            continue
          fi
          algo_args=
        elif [[ "$age" == new ]]; then
          algo_args="-a $algorithm"
        fi
        for i in 1 2 3; do
          # Compile results and wait if there are too many concurrent jobs running.
          while [[ $(count_elements "$unfinished") -ge "$max_parallel" ]]; do
            if [[ "$debug" ]]; then
              printf "Too many parallel jobs running (%d). Waiting for them to finish..\n" \
                     $(count_elements "$unfinished") >&2
            fi
            sleep 60
            print_newly_finished "$unfinished"
            unfinished=$(get_unfinished "$stats_files")
          done
          id="align.$age.$algorithm.$workers.$i"
          outfile=$(tempfile --prefix "out." --suffix .tsv)
          stats_file=$(tempfile --prefix "stats." --suffix .tsv)
          job_name="align$workers$age$i$algorithm"
          echo "Running $age script using $algorithm and $workers workers (replicate $i).." >&2
          "$measure_cmd" -i "$id" $slurm -j "$job_name" -o "$outfile" \
            "$path/$script_name" -p "$workers" $algo_args "$infile" > "$stats_file" &
          stats_files="$stats_files $stats_file"
          unfinished="$unfinished $stats_file"
        done
      done
    done
  done

  print_newly_finished "$unfinished"
  for stats_file in $unfinished; do
    rm "$stats_file"
  done
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


function get_version {
  path="$1"
  for script_path in "$path/align_families.py" "$path/align-families.py" ""; do
    if [[ -f "$script_path" ]]; then
      break
    fi
  done
  if [[ "$script_path" ]]; then
    set +e
    "$script_path" --version 2>&1
    set -e
  fi
}


function get_commit {
  path="$1"
  git --work-tree="$path" --git-dir="$path/.git" log -n 1 --format=%h
}


function print_newly_finished {
  unfinished="$1"
  for stats_file in $unfinished; do
    if [[ -s "$stats_file" ]]; then
      cat "$stats_file"
    fi
  done
}


function get_unfinished {
  stats_files="$1"
  for stats_file in $stats_files; do
    if ! [[ -s "$stats_file" ]]; then
      echo "$stats_file"
    fi
  done
}


function count_elements {
  i=0
  for arg in $@; do
    for element in $arg; do
      i=$((i+1))
    done
  done
  echo "$i"
}


function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
