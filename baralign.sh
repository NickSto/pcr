#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

Usage="Usage: \$ $(basename $0) [-R] families_file ref_dir out_file
families_file: The families.tsv produced by make-barcodes.awk and sorted.
ref_dir:  The directory to put the reference file (\"barcodes.fa\") and its index
          files in.
out_file: The path to put the output alignment BAM file at.
-R: Don't include reversed barcodes (alpha+beta -> beta+alpha) in the alignment target."

function main {

  # Read in arguments and check them.

  reverse=true
  while getopts ":rh" opt; do
  case "$opt" in
      r) reverse='';;
      h) fail "$Usage";;
    esac
  done
  # Get positional arguments.
  families=${@:$OPTIND:1}
  refdir=${@:$OPTIND+1:1}
  outfile=${@:$OPTIND+2:1}

  if ! [[ -f $families ]]; then
    fail "Error: families_file \"$families\" not found."
  fi
  if ! [[ -d $refdir ]]; then
    echo "Info: ref_dir \"$refdir\" not found. Creating.." >&2
    mkdir $refdir
  fi
  outbase=$(echo $outfile | sed -E 's/\.bam$//')
  if [[ $outbase == $outfile ]]; then
    fail "Error: out_file \"$outfile\" does not end in .bam."
  fi
  if [[ -e $outfile ]] || [[ -e $outbase.sam ]] || [[ -e $outbase.tmp.sam ]]; then
    fail "Error: out_file \"$outfile\" conflicts with existing filename(s)."
  fi

  for cmd in bowtie2 bowtie2-build samtools awk; do
    if ! which $cmd >/dev/null 2>/dev/null; then
      fail "Error: command \"$cmd\" not found."
    fi
  done

  echo "
families: $families
refdir:   $refdir
outfile:  $outfile
outbase:  $outbase" >&2

  # Create FASTA with barcodes as "reads" for alignment.
  awk '$1 != last {
    count++
    print ">" count
    print $1
  }
  {
    last = $1
  }' $families > $refdir/barcodes.fa

  # Create "reference" to align the barcodes to.
  if [[ $reverse ]]; then
    # If we're including reversed barcodes, create a new FASTA which includes reversed barcodes
    # as well as their forward versions.
    awk '
      $1 != last {
        count++
        bar = $1
        print ">" count
        print bar
        print ">" count ":rev"
        print swap_halves(bar)
      }
      {
        last = $1
      }
      function swap_halves(str) {
        half = length(str)/2
        alpha = substr(str, 1, half)
        beta = substr(str, half+1)
        return beta alpha
      }' $families > $refdir/barcodes-ref.fa
  else
    # If we're not including reversed barcodes, the original FASTA is all we need. Just link to it.
    ln -s $refdir/barcodes.fa $refdir/barcodes-ref.fa
  fi

  # Perform alignment.
  bowtie2-build --packed $refdir/barcodes-ref.fa $refdir/barcodes-ref >/dev/null
  bowtie2 -a -x $refdir/barcodes-ref -f -U $refdir/barcodes.fa -S $outbase.sam
  samtools view -Sb $outbase.sam > $outbase.tmp.bam
  samtools sort $outbase.tmp.bam $outbase
  if [[ -s $outfile ]]; then
    samtools index $outbase.bam
    rm $outbase.sam $outbase.tmp.bam
    echo "Success. Output located in \"$outfile\"." >&2
  else
    echo "Warning: No output file \"$outfile\" found." >&2
  fi
}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
