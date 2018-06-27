#!/usr/bin/env bash
if [ "x$BASH" = x ] || [ ! "$BASH_VERSINFO" ] || [ "$BASH_VERSINFO" -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

Usage="Usage: \$ $(basename "$0") [options] command [command args]
Measure the resources used by a command and all of its children.
Prints 8 tab-delimited columns:
0. The -i argument, if given.
1. ps %CPU: The maximum percentage of CPU used at any one time.
2. ps %MEM: The maximum percentage of RAM used at any one time.
3. ps  VSZ: The maximum virtual memory size (in KB) used at any one time.
4. ps  RSS: The maximum resident set size (in KB) used at any one time.
5. time -f %e: Total number of elapsed wall clock seconds.
6. time -f %S: Total CPU-seconds used by the kernel on behalf of the command.
7. time -f %U: Total CPU-seconds of user-mode execution time.
8. time -f %M: Maximum resident set size in KB.
-i: A string that will be prepended to the output as the first column.
-o: Output file. The command's stdout will be piped here. Default: /dev/null
-l: Log file. The command's stderr will be piped here. Default: /dev/null
-u: Select processes to monitor for ps stats (the first 4 stat columns) by
    finding all processes matching the current user ($USER). Otherwise, if not
    running via slurm, the processes will be found by their parent pid. If
    running via slurm, the processes will be found by their command line.
-s: Run under slurm. Prefixes the command with \"srun\".
-S: Arguments to give the slurm \"srun\" command.
-D: Turn on debug mode."

function main {

  # Get arguments.
  user=
  id=
  outfile=/dev/null
  logfile=/dev/null
  slurm=
  slurm_args=
  debug=
  while getopts "ui:o:l:sS:Dh" opt; do
    case "$opt" in
      u) user="$USER";;
      i) id="$OPTARG";;
      o) outfile="$OPTARG";;
      l) logfile="$OPTARG";;
      s) slurm=true;;
      S) slurm_args="$OPTARG";;
      D) debug='-D';;
      [h?]) fail "$Usage";;
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

  if [[ "$slurm" ]]; then
    if ! which srun >/dev/null 2>/dev/null; then
      fail "Error: srun command not found. Turn off slurm with -S."
    fi
    node=$(get_free_node)
    while ! [[ "$node" ]]; do
      if [[ "$debug" ]]; then
        echo "No slurm nodes are currently free. Waiting.." >&2
      fi
      sleep 60
      node=$(get_free_node)
    done
    if [[ "$debug" ]]; then
      echo "Launching command on $node" >&2
    fi
    if ! [[ -d ~/tmp ]]; then
      mkdir ~/tmp
    fi
    slurm_cmd="srun -w $node $slurm_args"
    monitor_prefix="ssh $node"
    monitor_selector="$command $command_args"
    time_file=$(tempfile -d ~/tmp --prefix time. --suffix .tsv)
    mem_file=$(tempfile -d ~/tmp --prefix mem. --suffix .tsv)
  else
    slurm_cmd=
    monitor_prefix=
    monitor_selector="-p $$"
    time_file=$(tempfile --prefix time. --suffix .tsv)
    mem_file=$(tempfile --prefix mem. --suffix .tsv)
  fi

  if [[ "$user" ]]; then
    monitor_selector='-u'
  fi

  time_cmd=$(which time)
  script_dir=$(get_script_dir)

  # Start the memory monitoring script.
  $monitor_prefix "$script_dir/mem-mon.sh" $debug $monitor_selector > "$mem_file" &

  # Run the actual command.
  $slurm_cmd "$time_cmd" -f '%e\t%S\t%U\t%M' -o "$time_file" "$command" $command_args 2> "$logfile" \
    | gzip -c - > "$outfile"

  # Wait for the memory monitoring script to register that the command finished.
  while ! [[ -s "$mem_file" ]]; do
    sleep 5
  done

  if ! [[ -s "$time_file" ]]; then
    echo "Error: time stats file $time_file missing or empty." >&2
  fi

  # Print the output.
  if [[ "$id" ]]; then
    echo -ne "$id\t"
  fi
  paste "$mem_file" "$time_file"

  # Clean up.
  rm "$mem_file" "$time_file"
  if [[ -f "$outfile" ]]; then
    rm "$outfile"
  elif [[ "$outfile" != /dev/null ]]; then
    echo "Error: Could not find and delete $outfile" >&2
  fi
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

function get_free_node {
  sinfo -h -p general -t idle -o %n | cut -d . -f 1 | tail -n 1
}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
