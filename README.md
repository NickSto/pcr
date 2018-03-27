# Du Novo

This is a pipeline for processing duplex sequencing data without the use of a reference genome.

The pipeline was designed for use with the duplex method described in [Kennedy *et al.* 2014](https://dx.doi.org/10.1038/nprot.2014.170), but the assumptions are relatively minimal, so you should be able to apply it to variants of the protocol.

Du Novo 2.0 is released under the GPLv2 license, except for some portions governed by the MIT license. Earlier versions were released under the BSD license. See `LICENSE.txt` for details.


## Running Du Novo from Galaxy

We created a comprehensive [tutorial](https://github.com/galaxyproject/dunovo/wiki) explaining all aspects of interactive use of Du Novo from within [Galaxy](http://usegalaxy.org).


## Running Du Novo on the command line


### Dependencies

#### Required

The pipeline requires a Unix operating system and Python version 2.7. Linux is recommended. OS X and BSD may work, but are untested.

It also requires several standard Unix tools. Version numbers in parentheses are what the software was tested with, but other versions likely work. These must be available on your [`$PATH`](https://en.wikipedia.org/wiki/Search_path):  
 -  the [`gcc`](https://gcc.gnu.org/) command (4.8.4)
 -  the [`make`](https://www.gnu.org/software/make/) command (3.81)
 -  the [`bash`](https://www.gnu.org/software/bash/bash.html) command (4.0)
 -  the [`awk`](https://www.gnu.org/software/gawk/) command (4.0.1)
 -  the [`paste`](https://www.gnu.org/software/coreutils/coreutils.html) command (8.21)
 -  the [`sort`](https://www.gnu.org/software/coreutils/coreutils.html) command (8.21)


#### Optional

To use `align-families.py`'s `-a mafft` option, this must be available on your `$PATH`:  
 - the [`mafft`](http://mafft.cbrc.jp/alignment/software/) command (v7.271 or v7.123b)

To use the barcode error correction script `correct.py`, the following modules must be available from Python and the commands must be on your `$PATH`:
 - the [networkx](https://pypi.python.org/pypi/networkx) Python module (1.9, 1.10, or 1.11)
 - the [`bowtie`](http://bowtie-bio.sourceforge.net/index.shtml) command (1.2.1.1) (nothing below 1.1.2 is confirmed to work)
 - the [`bowtie-build`](http://bowtie-bio.sourceforge.net/index.shtml) command (1.2.1.1) (same)
 - the [`samtools`](http://samtools.sourceforge.net/) command (0.1.18)


### Download

#### Git

    $ git clone --recursive https://github.com/galaxyproject/dunovo.git
    $ cd dunovo
    $ git checkout master
    $ git submodule update --recursive

#### Via the GitHub webpage

Click the [releases](https://github.com/galaxyproject/dunovo/releases) tab at the top of this page, and find the latest release. Download the zip file (the "Source code (zip)" link), as well as `kalign.zip`, `utillib.zip`, `ET.zip`, and `makrutenko.zip`. Extract the first zip file, then unzip the other three into the extracted directory and name those directories `kalign`, `utillib`, `ET`, and `makrutenko`, respectively.

In the end, the organization and names of the three directories should look like this:

    dunovo
    ├─╴kalign
    ├─╴utillib
    ├─╴ET
    └─╴makrutenko

### Installation

    $ cd dunovo
    $ make

The `make` command is needed to compile the C modules and kalign. You need to be in the root source directory (where the file `Makefile` is) before running the command.

### Testing

To check for obvious problems with the installation, you can run:

    $ tests/run.sh active

Successful results should look like

            [script] ::: [input file]:
    Files [test output] and [expected output] are identical

This won't catch every installation problem, but it should check that the basics are working.


### Usage

This example shows how to go from raw duplex sequencing data to the final duplex consensus sequences.

Your raw reads should be in `reads_1.fastq` and `reads_2.fastq`. And the scripts `make-families.sh`, `align-families.py`, `make-consensi.py`, `baralign.sh`, and `correct.py` should be on your `$PATH`.

1. Sort the reads into families based on their barcodes and split the barcodes from the sequence.  
    ```bash
    $ make-families.sh reads_1.fastq reads_2.fastq > families.tsv
    ```

2. (Optional) Correct errors in barcodes.
    ```bash
    $ baralign.sh families.tsv refdir barcodes.sam
    $ correct.py families.tsv refdir/barcodes.fa barcodes.sam | sort > families.corrected.tsv
    ```
    \- If you performed this step, change `families.tsv` below to `families.corrected.tsv`.

3. Do multiple sequence alignments of the read families.  
`$ align-families.py families.tsv > families.msa.tsv`

4. Build duplex consensus sequences from the aligned families.  
`$ make-consensi.py families.msa.tsv -1 duplex_1.fa -2 duplex_2.fa`

See all options for a given command by giving it the `-h` flag.


### Details

#### 1. Sort the reads into families based on their barcodes and split the barcodes from the sequence.  

    $ make-families.sh reads_1.fastq reads_2.fastq > families.tsv

This command will split the 12bp tag off each read, combine the tags from each pair of reads into a combined barcode, and sort them by it. The end result is a file (named `families.tsv` above) listing read pairs, grouped by barcode. See the `make-barcodes.awk` code for the details on the formation of the barcodes and the format.

Note: This step requires your FASTQ files to have exactly 4 lines per read (no multi-line sequences). 5' trimmed sequences of variable length are allowed. Also, in the output, the read sequence does not include the barcode or the 5bp constant sequence after it. You can customize the length of the barcode with the `-t` option or the constant sequence with the `-i` option.

#### 2. (Optional) Correct errors in barcodes.

    $ baralign.sh families.tsv refdir barcodes.sam
    $ correct.py families.tsv refdir/barcodes.fa barcodes.sam | sort > families.corrected.tsv

These commands takes the `families.tsv` file produced in the previous step, "corrects"\* the barcodes in it, and outputs a new version of `families.tsv` with the new barcodes. It does this by aligning all barcodes to themselves, finding pairs of barcodes which differ by only a few edits. Grouping sets of related barcodes gives groups which are likely descended from the same original barcode, differing only because of PCR and/or sequencing errors. By default, only barcodes that differ by 1 edit are allowed. You can allow greater edit distances between barcodes with the `--dist` option to `correct.py`.

\*"corrects" is in scare quotes because the algorithm isn't actually focused on finding the original barcode sequence. Its goal is instead to find barcodes which are all actually descended from the same original barcode, but now have different sequences because of errors. It finds each group of related barcodes and replaces them with a single barcode, so that the following steps identify them as one family.


#### 3. Do multiple sequence alignments of the read families.  

\- If you performed step 3, change `families.tsv` below to `families.corrected.tsv`.

`$ align-families.py families.tsv > families.msa.tsv`

This step aligns each family of reads, but it processes each strand separately. It can be parallelized with the `-p` option.

By default, this uses the [Kalign2](http://msa.sbc.su.se/cgi-bin/msa.cgi) multiple sequence alignment algorithm. Use `-a mafft` to select MAFFT instead. Kalign2 is reccommended, as its results are of similar accuracy and it's 7-8x faster.


#### 4. Build duplex consensus sequences from the aligned families.  

`$ make-consensi.py families.msa.tsv -1 duplex_1.fa -2 duplex_2.fa`

This calls a consensus sequence from the multiple sequence alignments of the previous step. It does this in two steps: First, single-strand consensus sequences (SSCSs) are called from the family alignments, then duplex consensus sequences are called from pairs of SSCSs.

When calling SSCSs, by default 3 reads are required to successfully create a consensus from each strand (change this with `-r`). Quality filtering is done at this step by excluding bases below a quality threshold. By default, no base with a PHRED quality less than 20 will contribute to the consensus (change this with `-q`). If no base passes the threshold or there is no majority base, `N` will be used.

The duplex consensus sequences are created by comparing the two SSCSs. For each base, if they agree, that base will be used. If they disagree, the IUPAC ambiguity code for the two bases will be used. Note that a disagreement between a base and a gap will result in an `N`.

The output of this step is the duplex consensus sequences in FASTA format. To also output all single-strand consensus sequences (including those which didn't produce a duplex consensus), use the `--sscs1` and `--sscs2` options.

The reads will be printed in two files, one per paired-end mate, with this naming format:  
`>{barcode} {# reads in strand 1 family}-{# reads in strand 2 family}`  
e.g.  
`>TTGCGCCAGGGCGAGGAAAATACT 8-13`

## Post-processing

When the consensus-calling process doesn't have enough information to call a base, it inserts an N or another IUPAC ambiguity code. This can happen in several cases, like when the two single-strand consensus sequences disagree, the PHRED quality is low, or the ends of the reads were trimmed.

The script `trimmer.py` in the `makrutenko` directory was written to deal with these bases. It can trim the ends of reads when they contain too many N's or ambiguous bases, or filter out reads with too many of them, or both.

It's a good idea to apply `trimmer.py` to at least remove sequence with a high density of ambiguous bases. This will result from any low-quality region or portion of the consensus where the raw reads were trimmed to different lengths:

`$ trimmer.py --acgt --window 10 --thres 0.3 --min-length 50 duplex_1.fa duplex_2.fa duplex.filt_1.fa duplex.filt_2.fa`

For an explanation of the arguments, see:

### Trimmer usage

This command will trim any read with more than 3 N's in a 10 base window, removing all the sequence after the first N in the offending window:

`$ trimmer.py --filt-bases N --window 10 --thres 0.3 consensus.fa > filtered.fa`

If our reads are 251bp, we can add `--min-length 251` to make it simply remove any of the reads that ever exceed the threshold:

`$ trimmer.py --filt-bases N --window 10 --thres 0.3 --min-length 251 consensus.fa > filtered.fa`

The `--acgt` argument will filter on any non-ACGT base instead of just N's:

`$ trimmer.py --acgt --window 10 --thres 0.3 consensus.fa > filtered.fa`

The script also handles paired-end data, preserving pairs by removing both reads when any one of them needs to be filtered out:

`$ trimmer.py --acgt --window 10 --thres 0.3 cons_1.fa cons_2.fa filtered_1.fa filtered_2.fa`


### Known bugs

Be aware that a [known bug](https://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool/1408476#1408476) in Python when using the [multiprocessing](https://docs.python.org/2/library/multiprocessing.html) module makes it impossible to kill a running process with the Ctrl+C keyboard command. So if you run `align-families.py` or `make-consensi.py` in the foreground, you'll have to exit via Ctrl+Z to stop and background the job, then kill the process (e.g. `$ kill %1`, if it's the only backgrounded job).
