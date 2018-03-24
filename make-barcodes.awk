# The awk code that transforms the one-line fastq record pair into the output that can be sorted
# by barcode.
# Input columns (the 4 FASTQ lines for both reads in a read pair):
#   1: read1 name
#   2: read2 name
#   3: read1 sequence
#   4: read2 sequence
#   5: read1 + line
#   6: read2 + line
#   7: read1 quality
#   8: read2 quality
# Output columns:
#   1: the barcode, put into a canonical form
#   2: the order of the barcode halves ("ab" or "ba")
#   3: read1 name
#   4: sequence of read 1, minus the 12bp barcode and 5bp invariant sequence
#   5: read1 quality scores, minus the same first 17bp
#   6: read2 name
#   7: sequence of read 2, minus the first 17bp
#   8: read2 quality scores, minus the first 17bp
# The canonical form of the barcode is composed of two concatenated tags, one from each read.
# By default, each tag is the first 12bp of the read. The tag from the first read is the "alpha" and
# the tag from the second is the "beta". The barcode is formed by concatenating them in an order
# determined by a string comparison of the two. The lesser tag is first (if they are equal, the
# beta is first, but then you have bigger problems).

BEGIN {
  FS = "\t"
  OFS = "\t"
  # The number of bases from the start of each read that form the two halves of the barcode.
  # (this should be half the size of the full, canonical barcode).
  if (TAG_LEN == "") {
    TAG_LEN = 12
  }
  # The number of bases in the read that are between the barcode and the start of the actual sample
  # sequence (the restriction site in the Loeb 2014 protocol).
  if (INVARIANT == "") {
    INVARIANT = 5
  }
}

# Validate input.
{
  for (i = 1; i <= 8; i++) {
    if (!$i) {
      print "Error on read pair (line) " NR ": Missing or empty columns." > "/dev/stderr"
      exit 1
    }
  }
  if (substr($1, 1, 1) != "@" || substr($2, 1, 1) != "@") {
    print "Error on read pair (line) " NR ": Read name line doesn't begin with '@'." > "/dev/stderr"
    exit 1
  }
  if (substr($5, 1, 1) != "+" || substr($6, 1, 1) != "+") {
    print "Error on read pair (line) " NR ": '+' column does not begin with '+'." > "/dev/stderr"
    exit 1
  }
}

$3 && $4 {
  alpha = substr($3, 1, TAG_LEN)
  beta = substr($4, 1, TAG_LEN)
  if (alpha < beta) {
    barcode = alpha beta
    order = "ab"
  } else {
    barcode = beta alpha
    order = "ba"
  }
  name1 = substr($1, 2)
  name2 = substr($2, 2)
  seq1 = substr($3, TAG_LEN + INVARIANT + 1)
  seq2 = substr($4, TAG_LEN + INVARIANT + 1)
  qual1 = substr($7, TAG_LEN + INVARIANT + 1)
  qual2 = substr($8, TAG_LEN + INVARIANT + 1)
  print barcode, order, name1, seq1, qual1, name2, seq2, qual2
}
