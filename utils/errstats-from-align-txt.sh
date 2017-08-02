#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue
# Get the root directory of this repo.
scripts_dir=$(dirname $(dirname $(readlink -f $0)))
export PYTHONPATH=${PYTHONPATH:=$HOME/bx/code/indels/pybamparser/lib}

Usage="Usage: \$ $(basename $0) testreads.align.txt
Take a human-readable alignment, produce input reads from it using parse-test-align.py, and run it
through the pipeline and errstats.py.
The output will be in a new directory alongside the input file, with the same name but \"txt\"
replaced with \"tmp\" (testreads.align.txt -> testreads.align.tmp/)."

function main {
  if [[ $# -lt 1 ]] || [[ $1 == '-h' ]] || [[ $1 == '--help' ]]; then
    fail "$Usage"
  fi

  infile=$1
  if ! [[ -f $infile ]]; then
    echo "$infile not found." >&2
    return 1
  fi
  dir=$(dirname $1)/$(basename $1 .txt).tmp
  if ! [[ -d $dir ]]; then
    mkdir $dir
  fi
  $scripts_dir/tests/parse-test-align.py $infile -r $dir/ref.fa -1 $dir/reads_1.fq -2 $dir/reads_2.fq
  $scripts_dir/tests/parse-test-align.py $infile -r $dir/ref.fa -b 0 -c 0 \
    -1 $dir/reads.clipped_1.fq -2 $dir/reads.clipped_2.fq
  align $dir/ref.fa $dir/reads.clipped_1.fq
  paste $dir/reads_1.fq $dir/reads_2.fq | paste - - - - | awk -f $scripts_dir/make-barcodes.awk \
    | sort > $dir/families.tsv
  $scripts_dir/align_families.py $dir/families.tsv > $dir/families.msa.tsv
  $scripts_dir/dunovo.py --sscs-file $dir/sscs.fa $dir/families.msa.tsv > /dev/null
  $scripts_dir/utils/outconv.py $dir/sscs.fa -1 $dir/sscs_1.fa -2 $dir/sscs_2.fa
  align $dir/ref.fa $dir/sscs_1.fa
  $scripts_dir/utils/errstats.py -v --all-repeats --dedup --qual-errors -q 25 --min-reads 3 \
    --overlap-stats $dir/overlaps.tsv --dedup-log $dir/dedup.log --bam $dir/sscs.bam \
    $dir/families.msa.tsv > $dir/errstats.dedup.allrepeats.filt3.tsv
}

function align {
  ref=$1
  reads1=$2
  base=$(dirname $reads1)/$(basename $reads1 _1.fq)
  if [[ $base == $reads1 ]]; then
    format=-f
    base=$(dirname $reads1)/$(basename $reads1 _1.fa)
    reads2=${base}_2.fa
    if [[ $base == $reads1 ]]; then
      return 1
    fi
  else
    format=
    reads2=${base}_2.fq
  fi
  bowtie2-build $ref $ref >/dev/null
  bowtie2 -x $ref $format -1 $reads1 -2 $reads2 -S $base.sam 2>/dev/null
  samtools view -Sb $base.sam > $base.tmp.bam
  samtools sort $base.tmp.bam $base
  samtools index $base.bam
  rm $base.sam $base.tmp.bam
}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
