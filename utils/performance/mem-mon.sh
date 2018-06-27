#!/usr/bin/env bash
if [ "x$BASH" = x ] || [ ! "$BASH_VERSINFO" ] || [ "$BASH_VERSINFO" -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

STAT_NAMES=(cpu mem vsz rss)
SleepDefault=5
DebugPrintInterval=5 # minutes
Usage="Usage: \$ $(basename "$0") [options] [-u | -k key | -p ppid | command arg1 arg2 etc]
Monitor the resource usage of a set of processes and print the maximums.
You must specify one method of selecting processes to monitor:
-u: Monitor all processes owned by the current user ($USER).
-k: Use this key to search for processes. This should be one of the arguments
    included in the command for each process.
-p: Find all processes descended from the process with this pid
    (children, grandchildren, etc).
Or, give a command as positional arguments to this script. It will monitor any
processes whose command matches those arguments. You can give a partial command,
omitting arguments at the end. If you give X arguments, it will match any
process whose first X arguments match the ones you gave.
Note: This script's process is always excluded.
Options:
-i: An id to prefix the final output line.
    If not given, this column will be omitted.
-s: How often to check the process statistics. Default: $SleepDefault seconds.
-D: Turn on debug printing every $DebugPrintInterval minutes."

function main {

  # Get arguments.
  user=
  key=
  ppid=
  id=
  sleep="$SleepDefault"
  debug=
  while getopts "uk:p:i:s:Dh" opt; do
    case "$opt" in
      u) user="$USER";;
      k) key="$OPTARG";;
      p) ppid="$OPTARG";;
      i) id="$OPTARG";;
      s) sleep="$OPTARG";;
      D) debug=true;;
      [h?]) fail "$Usage";;
    esac
  done
  command_args=${@:$OPTIND}

  if ! [[ "$user" ]] && ! [[ "$key" ]] && ! [[ "$ppid" ]] && ! [[ "$command_args" ]]; then
    fail "Error: You must provide a -k key, -p ppid, or a command."
  fi

  # Initialize the arrays.
  declare -a maxes stats
  for i in 0 1 2 3; do
    maxes[$i]=0
    stats[$i]=0
  done

  # Main loop: monitor the process.
  last_report=$(date +%s)
  stats[0]=1
  first_loop=true
  while [[ "${stats[0]}" -gt 0 ]]; do
    sleep "$sleep"
    if [[ "$user" ]]; then
      pids=$(get_pids_by_user "$user")
    elif [[ "$key" ]]; then
      pids=$(get_pids_by_key "$key")
    elif [[ "$ppid" ]]; then
      pids=$(get_pids_by_ppid "$ppid")
    else
      pids=$(get_pids_by_cmd $command_args)
    fi
    if ! [[ "$pids" ]]; then
      if [[ "$first_loop" ]]; then
        echo "Error: Did not find any processes." >&2
      else
        echo "No more living processes. Ending monitoring." >&2
      fi
      break
    fi
    first_loop=
    read -a stats <<< $(get_stats_by_pids $pids)
    # Check for maximum stats.
    i=0
    while [[ "$i" -lt "${#stats[@]}" ]]; do
      if [[ "${stats[$i]}" -gt "${maxes[$i]}" ]]; then
        maxes[$i]="${stats[$i]}"
      fi
      i=$((i+1))
    done
    # Debug printing.
    now=$(date +%s)
    if [[ "$debug" ]] && [[ $((now-last_report)) -gt $((DebugPrintInterval*60)) ]]; then
      output="$now"
      for stat in "${stats[@]}"; do
        output="$output\t$stat"
      done
      echo -e "$output" >&2
      last_report="$now"
    fi
  done

  # Print the final maximums.
  output=""
  i=0
  while [[ "$i" -lt "${#maxes[@]}" ]]; do
    max="${maxes[$i]}"
    if ! [[ "$max" ]]; then
      fail "Error: No max value detected for stat ${STAT_NAMES[$i]}.
All max values:   ${maxes[@]}
Last stat values: ${stats[@]}"
    fi
    if [[ "$i" -lt 2 ]]; then
      max=$(calc "$max/100")
    fi
    # Prepend with the id, if there is one.
    if [[ "$i" == 0 ]]; then
      if [[ "$id" ]]; then
        output="$id\t$max"
      else
        output="$max"
      fi
    else
      output="$output\t$max"
    fi
    i=$((i+1))
  done
  echo -e "$output"
}


function recursive_pgrep {
  ppid="$1"
  for pid in $(pgrep -P "$ppid"); do
    echo "$pid"
    recursive_pgrep "$pid"
  done
}


function get_pids_by_user {
  user="$1"
  ps aux | awk '$1 == "'"$user"'" && $2 != '"$$"' {print $2}'
}


function get_pids_by_key {
  key="$1"
  ps aux | awk '
    $1 == "'"$USER"'" && $2 != '"$$"' {
      hit = 0
      for (i=11; i<=NF; i++) {
        if ($i == "'"$key"'") {
          print $2
          break
        }
      }
    }'
}


function get_pids_by_ppid {
  ppid="$1"
  pid_conditional=
  for pid in $(recursive_pgrep "$ppid"); do
    if [[ "$pid_conditional" ]]; then
      pid_conditional="$pid_conditional || \$2 == $pid"
    else
      pid_conditional="\$2 == $pid"
    fi
  done
  if ! [[ "$pid_conditional" ]]; then
    echo "no pids found for $ppid!" >&2
    return
  fi
  ps aux | awk '
    $1 == "'"$USER"'" && $2 != '"$$"' && ('"$pid_conditional"') {
      print $2
    }'
}


function get_pids_by_cmd {
  command_args=$@
  i=11
  cmd_conditionals=
  for arg in $command_args; do
    cmd_conditional="\$$i == \"$arg\""
    if [[ "$cmd_conditionals" ]]; then
      cmd_conditionals="$cmd_conditionals && $cmd_conditional"
    else
      cmd_conditionals="$cmd_conditional"
    fi
    i=$((i+1))
  done
  ps aux | awk "\$1 == \"$USER\" && $2 != '"$$"' && $cmd_conditionals {print \$2}"
}


function get_stats_by_pids {
  pids=$@
  if ! [[ "$pids" ]]; then
    echo
    return
  fi
  pid_conditionals=$(get_pid_conditionals $pids)
  ps aux | awk '
    BEGIN {
      OFS = "\t"
    }
    '"$pid_conditionals"' {
      cpu += $3
      mem += $4
      vsz += $5
      rss += $6
    }
    END {
      print int(100*cpu), int(100*mem), vsz, rss
    }'
}


function get_pid_conditionals {
  pids=$@
  pid_conditionals=
  for pid in $pids; do
    pid_conditional="\$2 == $pid"
    if [[ "$pid_conditionals" ]]; then
      pid_conditionals="$pid_conditionals || $pid_conditional"
    else
      pid_conditionals="$pid_conditional"
    fi
  done
  echo "$pid_conditionals"
}


function calc {
  python3 -c "from math import *; print($*)";
}


function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
