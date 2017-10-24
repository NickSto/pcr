#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

RefdirDefault=refdir

Usage="Usage: \$ $(basename $0) [-R] families.tsv [refdir [outfile.sam|outfile.bam]]
families.tsv: The families file produced by make-barcodes.awk and sorted.
refdir: The directory to put the reference file (\"barcodes.fa\") and its index
        files in. Default: \"$RefdirDefault\".
outfile.bam: Print the output to this path. It will be in SAM format unless the
             path ends in \".bam\". If not given, it will be printed to stdout
             in SAM format.
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

  # Validate arguments.
  if ! [[ $refdir ]]; then
    refdir=$RefdirDefault
  fi
  if ! [[ -d $refdir ]]; then
    echo "Info: ref_dir \"$refdir\" not found. Creating.." >&2
    mkdir $refdir
  fi
  if ! [[ -f $families ]]; then
    fail "Error: families_file \"$families\" not found."
  fi
  # Determine how and where to put the output.
  if [[ ${outfile:${#outfile}-4} == .bam ]]; then
    format=bam
  else
    format=sam
  fi
  outbase=$(echo $outfile | sed -E 's/\.bam$//')
  if [[ $outfile ]]; then
    if [[ -e $outfile ]]; then
      fail "Error: output file \"$outfile\" already exists."
    fi
    if [[ $format == bam ]]; then
      if [[ -e $outbase.sam ]] || [[ -e $outbase.bam.bai ]]; then
        fail "Error: A related filename already exists (.sam/.bam.bai)."
      fi
      bowtie2_outargs="-S $outbase.sam"
    else
      bowtie2_outargs="-S $outfile"
    fi
  else
    bowtie2_outargs=
  fi

  # Check for required commands.
  for cmd in bowtie2 bowtie2-build samtools awk; do
    if ! which $cmd >/dev/null 2>/dev/null; then
      fail "Error: command \"$cmd\" not found."
    fi
  done

  echo "
families: $families
refdir:   $refdir
format:   $format
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
  bowtie2 -a -x $refdir/barcodes-ref -f -U $refdir/barcodes.fa $bowtie2_outargs
  if [[ $format == bam ]]; then
    samtools view -Sb $outbase.sam | samtools sort -o - dummy > $outfile
    if [[ -s $outfile ]]; then
      samtools index $outfile
    fi
  fi
  # Check output.
  if [[ $outfile ]]; then
    if [[ -s $outfile ]]; then
      if [[ $format == bam ]] && [[ -e $outbase.sam ]]; then
        rm $outbase.sam
      fi
      echo "Success. Output located in \"$outfile\"." >&2
    else
      echo "Warning: No output file \"$outfile\" found." >&2
    fi
  fi
}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
