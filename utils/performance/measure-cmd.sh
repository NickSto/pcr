#!/usr/bin/env bash
if [ "x$BASH" = x ] || [ ! "$BASH_VERSINFO" ] || [ "$BASH_VERSINFO" -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

SlurmPrefix='srun -C new --exclusive'
Usage="Usage: \$ $(basename "$0") [options] command [command args]
-i: A string that will be prepended to the output as the first column.
-o: Output file. The command's stdout will be piped here. Default: /dev/null
-s: Run under slurm. Prefixes the command with $SlurmPrefix.
-j: Job name for slurm."

function main {

  # Get arguments.
  id=
  outfile=/dev/null
  job_name=
  slurm=
  while getopts "i:o:j:sh" opt; do
    case "$opt" in
      i) id="$OPTARG";;
      o) outfile="$OPTARG";;
      j) job_name="$OPTARG";;
      s) slurm=true;;
      h) fail "$Usage";;
      ?) fail "$Usage";;
    esac
  done
  command="${@:$OPTIND:1}"
  command_args="${@:$OPTIND+1}"

  if ! [[ "$command" ]]; then
    fail "$Usage"
  fi

  if ! [[ "$outfile" ]]; then
    fail "Error: output filename cannot be empty."
  fi
  if ! [[ -d "$(dirname "$outfile")" ]]; then
    fail "Error: Path invalid: $outfile"
  fi

  if [[ "$job_name" ]]; then
    job_args="-J $job_name"
    time_file=$(tempfile --prefix "time.$job_name." --suffix .tsv)
    mem_file=$(tempfile --prefix "mem.$job_name." --suffix .tsv)
  else
    job_args=
    time_file=$(tempfile --prefix time. --suffix .tsv)
    mem_file=$(tempfile --prefix mem. --suffix .tsv)
  fi

  if [[ "$slurm" ]]; then
    if ! which srun >/dev/null 2>/dev/null; then
      fail "Error: srun command not found. Turn off slurm with -S."
    fi
    slurm_args="$SlurmPrefix $job_args"
  else
    slurm_args=
  fi

  time_cmd=$(which time)
  script_dir=$(get_script_dir)

  # Start the memory monitoring script.
  "$script_dir/mem-mon.sh" -p "$$" > "$mem_file" &

  # Run the actual command.
  $slurm_args "$time_cmd" -f '%e\t%S\t%U\t%M' -o "$time_file" "$command" $command_args \
    > "$outfile" 2>/dev/null

  # Wait for the memory monitoring script to register that the command finished.
  while ! [[ -s "$mem_file" ]]; do
    sleep 5
  done

  # Print the output.
  if [[ "$id" ]]; then
    echo -ne "$id\t"
  fi
  paste "$mem_file" "$time_file"

  # Clean up.
  rm "$mem_file" "$time_file"
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
