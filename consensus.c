#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

// N.B. This defines the valid bases, but it's also effectively defined in the switches in
// get_votes_simple(), get_votes_qual(), and get_base_prime(), and in the constant IUPAC_BASES.
#define N_BASES 6
const char *BASES = "ACGTN-";
/* A  C   G   T   N   -     A: 2    Compute IUPAC ambiguous base character by representing each base
A  4  6  10  14  22  26     C: 3    with a prime and multiplying. Then use a lookup table (an array
C     9  15  21  33  39     G: 5    where the index is the product of the two primes).
G        25  35  55  65     T: 7
T            49  77  91     N: 11
N               121 143     -: 13    1         2         3         4         5         6         7
-                   169    01234567890123456789012345678901234567890123456789012345678901234567890*/
const char *IUPAC_BASES = "N...A.M..CR...WS.....YN..GN......N.K...N.........T.....N.........N....."
//                                   8         9        10        11        12        13        14
                           "......N.............N.............................N..................."
//                                  15        16        17
                           "..N.........................-";
#define THRES_DEFAULT 0.5
#define MIN_DEFAULT 0
#define WIN_LEN 4
#define GAP_CHAR ' '

int **get_votes_simple(char *align[], int n_seqs, int seq_len);
int **get_votes_qual(char *align[], char *quals[], int n_seqs, int seq_len, char thres);
int init_gap_qual_window(int *window, char *quals, int seq_len);
char get_gap_qual(int *window);
int push_qual(int *window, int win_edge, char *quals, int seq_len);
void print_window(int *window, int win_edge);
int **init_votes(int seq_len);
void free_votes(int *votes[], int seq_len);
void print_votes(char *consensus, int *votes[], int seq_len);
char *rm_gaps(char *consensus, int cons_len);
char *build_consensus(int *votes[], int seq_len, double thres, int min_reads);
char *build_consensus_duplex(int *votes1[], int *votes2[], int seq_len, double thres);
char *build_consensus_duplex_simple(char *cons1, char *cons2, int gapped);
int get_base_prime(char base);
char *get_consensus(char *align[], char *quals[], int n_seqs, int seq_len, double cons_thres,
                    int min_reads, char qual_thres, int gapped);
char *get_consensus_duplex(char *align1[], char *align2[], char *quals1[], char *quals2[],
                           int n_seqs1, int n_seqs2, int seq_len, double cons_thres,
                           int min_reads, char qual_thres, int gapped, char *method);


// Tally the different bases at each position in an alignment.
// Returns an array of arrays: for each position in the alignment, an array of the number of times
// each base occurs at that position. The order of bases is as in the "BASES" constant.
int **get_votes_simple(char *align[], int n_seqs, int seq_len) {
  int **votes = init_votes(seq_len);

  // Tally votes for each base.
  int i, j;
  for (i = 0; i < n_seqs; i++) {
    for (j = 0; j < seq_len; j++) {
      // N.B.: Could write this without hardcoded literals, but it's about 40% slower.
      switch (toupper(align[i][j])) {
        case 'A':
          votes[j][0]++;
          break;
        case 'C':
          votes[j][1]++;
          break;
        case 'G':
          votes[j][2]++;
          break;
        case 'T':
          votes[j][3]++;
          break;
        case 'N':
          votes[j][4]++;
          break;
        case '-':
          votes[j][5]++;
          break;
      }
    }
  }

  return votes;
}


//  Tally votes for each base, ignoring bases with a quality score below "thres".
int **get_votes_qual(char *align[], char *quals[], int n_seqs, int seq_len, char thres) {
  int **votes = init_votes(seq_len);
  int *window = malloc(sizeof(int) * WIN_LEN * 2);
  int win_edge;

  // Tally votes for each base.
  char qual;
  int i, j;
  for (i = 0; i < n_seqs; i++) {
    win_edge = init_gap_qual_window(window, quals[i], seq_len);
    for (j = 0; j < seq_len; j++) {
      // Figure out the quality score of the base (or gap).
      if (align[i][j] == '-') {
        qual = get_gap_qual(window);
      } else {
        win_edge = push_qual(window, win_edge, quals[i], seq_len);
        qual = quals[i][j];
      }
      // Don't count bases whose quality is less than the threshold.
      if (qual < thres) {
        continue;
      }
      // N.B.: Could write this without hardcoded literals, but it's about 40% slower.
      switch (toupper(align[i][j])) {
        case 'A':
          votes[j][0]++;
          break;
        case 'C':
          votes[j][1]++;
          break;
        case 'G':
          votes[j][2]++;
          break;
        case 'T':
          votes[j][3]++;
          break;
        case 'N':
          votes[j][4]++;
          break;
        case '-':
          votes[j][5]++;
          break;
      }
    }
  }

  free(window);
  return votes;
}


/* Tally votes for each base, weighting by the PHRED score of the base.
 * This is based on the theory of PHRED scores representing the literal probability of the base call
 * being erroneous. Thus, if two reads show a C at a position, both with PHRED 20 (1/100 chance of
 * error), then the chances of them both being wrong are 1/100 * 1/100 = 1/10000 (PHRED 40).
 * So, theoretically and intuitively, it makes sense to trust two PHRED 20 C's over one PHRED 30 A.
 * This seems better than arbitrarily deciding not to consider bases below a PHRED score threshold.
 * How to decide when not to call the base? We could just say that if no base's score total is above
 * a certain threshold, we call it an N. Theoretically, this threshold is the confidence we want in
 * our final base calls. This could even replace the arbitrary 3 reads for a consensus threshold.
 */
int **get_votes_weighted(char *align[], char *quals[], int n_seqs, int seq_len) {
  int **votes = init_votes(seq_len);
  int *window = malloc(sizeof(int) * WIN_LEN * 2);
  int win_edge;

  // Tally votes for each base.
  char qual;
  int i, j;
  for (i = 0; i < n_seqs; i++) {
    win_edge = init_gap_qual_window(window, quals[i], seq_len);
    for (j = 0; j < seq_len; j++) {
      // Figure out the quality score of the base (or gap).
      if (align[i][j] == '-') {
        qual = get_gap_qual(window);
      } else {
        win_edge = push_qual(window, win_edge, quals[i], seq_len);
        qual = quals[i][j];
      }
      // N.B.: Could write this without hardcoded literals, but it's about 40% slower.
      switch (toupper(align[i][j])) {
        case 'A':
          votes[j][0] += qual;
          break;
        case 'C':
          votes[j][1] += qual;
          break;
        case 'G':
          votes[j][2] += qual;
          break;
        case 'T':
          votes[j][3] += qual;
          break;
        case 'N':
          votes[j][4] += qual;
          break;
        case '-':
          votes[j][5] += qual;
          break;
      }
    }
  }

  free(window);
  return votes;
}


/* Calculation of gap quality scores:
 * "window" is an array of length 2*WIN_LEN, holding the quality scores of the non-gap bases WIN_LEN
 * from the current one in both directions. When we're at the start or end, fill the slots beyond
 * the edge of the sequence with -1 as a sentinel for "N/A". For example, after the 2nd base in the
 * sequence, if WIN_LEN is 4 and there are no gaps yet, "window" should look something like this:
 *   base coordinates         0  1   2  3  4  5
 *   quality scores           =  9   A  >  E  B
 *   array values     [-1|-1|28|24||32|29|36|33]
 * Usage:
 * Allocate "window" and call init_gap_qual_window() to initialize it by filling its right side with
 * the first WIN_LEN non-gap quality scores in the sequence. It will also return "win_edge", which
 * tracks the coordinate of the rightmost quality score in "window".
 * Then, start traversing the array of quality scores. For each score, if it's a gap character, call
 * get_gap_qual() to get the computed quality score at the gap. Otherwise, call push_qual() to add
 * another quality score to the end of the window. This will keep the window filled with all non-gap
 * quality scores, and each time you call get_gap_qual(), you should be at the gap between the two
 * quality scores at the center of the window.
 */

// This does the initial fill of the window array, adding the first WIN_LEN bases to the right side
// and filling the left side with -1's.
int init_gap_qual_window(int *window, char *quals, int seq_len) {
  // Fill left side with -1's (no quality information).
  int i;
  for (i = 0; i < WIN_LEN; i++) {
    window[i] = -1;
  }
  // Fill right side with first WIN_LEN quality scores. Skip gaps, and if you run out of quality
  // scores (if seq_len < WIN_LEN), fill the rest with -1. Leave win_edge at the last base we added.
  i = WIN_LEN;
  int win_edge = -1;
  int quals_added = 0;
  while (quals_added < WIN_LEN) {
    win_edge++;
    if (win_edge >= seq_len) {
      win_edge = seq_len;
      window[i] = -1;
      i++;
      quals_added++;
    } else if (quals[win_edge] != GAP_CHAR) {
      window[i] = quals[win_edge];
      i++;
      quals_added++;
    }
  }
  return win_edge;
}


// Push the next non-gap quality score onto the right side of the window.
int push_qual(int *window, int win_edge, char *quals, int seq_len) {
  // Find the next quality score that's not a gap.
  char next_qual = GAP_CHAR;
  while (next_qual == GAP_CHAR) {
    win_edge++;
    if (win_edge < seq_len) {
      next_qual = quals[win_edge];
    } else {
      win_edge = seq_len;
      next_qual = -1;
    }
  }
  // Shift all the quality scores left add the new one.
  int last_qual;
  int i;
  for (i = WIN_LEN*2 - 1; i >= 0; i--) {
    last_qual = window[i];
    window[i] = next_qual;
    next_qual = last_qual;
  }
  return win_edge;
}


// Compute the quality of the gap based on a weighted average of the quality scores in the window.
// The scores near the center of the window are weighted higher than the ones further away.
char get_gap_qual(int *window) {
  int score_sum = 0;
  int weight_sum = 0;
  int weight = 1;
  int i;
  for (i = 0; i < WIN_LEN*2; i++) {
    if (window[i] != -1) {
      score_sum += window[i] * weight;
      weight_sum += weight;
    }
    // Increase the weight until we get to the middle of the window (at WIN_LEN), then decrease it.
    if (i < WIN_LEN - 1) {
      weight++;
    } else if (i > WIN_LEN - 1) {
      weight--;
    }
  }
  if (weight_sum > 0) {
    // Divide by the sum of the weights to get the final quality score.
    return (char) (score_sum/weight_sum);
  } else {
    return '\0';
  }
}


void print_window(int *window, int win_edge) {
  printf("[");
  int i;
  for (i = 0; i < WIN_LEN*2; i++) {
    printf("%c", window[i]);
    if (i == WIN_LEN - 1) {
      printf("||");
    } else if (i == WIN_LEN*2 - 1) {
      printf("] %-2d\n", win_edge);
    } else {
      printf("|");
    }
  }
}


int **init_votes(int seq_len) {
  int **votes = malloc(sizeof(int *) * seq_len);
  int i, j;
  for (i = 0; i < seq_len; i++) {
    votes[i] = malloc(sizeof(int) * N_BASES);
    for (j = 0; j < N_BASES; j++) {
      votes[i][j] = 0;
    }
  }
  return votes;
}


void free_votes(int *votes[], int seq_len) {
  int i;
  for (i = 0; i < seq_len; i++) {
    free(votes[i]);
  }
  free(votes);
}


void print_votes(char *consensus, int *votes[], int seq_len) {
  int i, j;
  printf("   ");
  for (j = 0; j < N_BASES; j++) {
    printf(" %c ", BASES[j]);
  }
  printf("\n");
  for (i = 0; i < seq_len; i++) {
    printf("%c: ", consensus[i]);
    for (j = 0; j < N_BASES; j++) {
      if (votes[i][j]) {
        printf("%2d ", votes[i][j]);
      } else {
        printf("   ");
      }
    }
    printf("\n");
  }
}


// Take a consensus sequence which may have gaps ('-' characters) and remove them to produce the
// actual final sequence. "cons_len" should be the length of the original, gapped, sequence.
char *rm_gaps(char *consensus, int cons_len) {
  char *output = malloc(sizeof(char) * cons_len + 1);
  int i;
  int j = 0;
  for (i = 0; i < cons_len; i++) {
    if (consensus[i] != '-') {
      output[j] = consensus[i];
      j++;
    }
  }
  output[j] = '\0';
  return output;
}


char *build_consensus(int *votes[], int seq_len, double thres, int min_reads) {
  char *consensus = malloc(sizeof(char) * seq_len + 1);
  
  int i, j;
  for (i = 0; i < seq_len; i++) {
    int total = 0;
    int max_vote = 0;
    char max_base = 'N';
    for (j = 0; j < N_BASES; j++) {
      total += votes[i][j];
      if (votes[i][j] > max_vote) {
        max_vote = votes[i][j];
        max_base = BASES[j];
      }
      if (total == 0) {
        consensus[i] = 'N';
      } else if ((double)max_vote/total > thres && max_vote >= min_reads) {
        consensus[i] = max_base;
      } else {
        consensus[i] = 'N';
      }
    }
  }

  consensus[seq_len] = '\0';
  return consensus;
}


// Build a consensus sequence from two alignments by weighting each equally and considering only
// the frequency of each base in each alignment.
char *build_consensus_duplex(int *votes1[], int *votes2[], int seq_len, double thres) {
  char *consensus = malloc(sizeof(char) * seq_len + 1);

  int i, j;
  for (i = 0; i < seq_len; i++) {
    // Sum the total votes at this position.
    /*TODO: This does an extra loop through the votes to get the total so it can calculate actual
     *      frequencies in the second pass. Technically, this information could be gathered when
     *      originally tallying the votes in the get_votes functions. Or, the total could be assumed
     *      to be n_seqs if every base always contributes a vote (even when it's not in "ACGTN-").
     */
    int total1 = 0;
    for (j = 0; j < N_BASES; j++) {
      total1 += votes1[i][j];
    }
    int total2 = 0;
    for (j = 0; j < N_BASES; j++) {
      total2 += votes2[i][j];
    }
    double max_freq = 0.0;
    char max_base = 'N';
    for (j = 0; j < N_BASES; j++) {
      // Get the frequency of each base.
      double freq1;
      if (total1 > 0) {
        freq1 = (double)votes1[i][j]/total1;
      }
      double freq2;
      if (total2 > 0) {
        freq2 = (double)votes2[i][j]/total2;
      }
      // frequency of the base = average of frequencies in the two sequences
      double avg_freq;
      if (total1 == 0 && total2 == 0) {
        avg_freq = -1.0;
      } else if (total1 == 0) {
        avg_freq = freq2;
      } else if (total2 == 0) {
        avg_freq = freq1;
      } else {
        avg_freq = (freq1 + freq2) / 2;
      }
      // Track the highest frequency seen.
      if (avg_freq > max_freq) {
        max_freq = avg_freq;
        max_base = BASES[j];
      }
    }
    if (max_freq > thres) {
      consensus[i] = max_base;
    } else {
      consensus[i] = 'N';
    }
  }

  consensus[seq_len] = '\0';
  return consensus;
}


// "cons1" and "cons2" must be null-terminated strings of equal lengths.
char *build_consensus_duplex_simple(char *cons1, char *cons2, int gapped) {
  int seq_len = strlen(cons1);
  char *cons = malloc(sizeof(char) * seq_len + 1);
  int i = 0;
  int base_prime1, base_prime2;
  while (cons1[i] != '\0' && cons2[i] != '\0') {
    base_prime1 = get_base_prime(cons1[i]);
    base_prime2 = get_base_prime(cons2[i]);
    cons[i] = IUPAC_BASES[base_prime1*base_prime2];
    i++;
  }
  cons[seq_len] = '\0';
  if (gapped) {
    return cons;
  } else {
    return rm_gaps(cons, seq_len);
  }
}


int get_base_prime(char base) {
  switch (base) {
    case 'A':
      return 2;
    case 'C':
      return 3;
    case 'G':
      return 5;
    case 'T':
      return 7;
    case 'N':
      return 11;
    case '-':
      return 13;
    default:
      return 0;
  }
}


// Convenience function to create a consensus in one step.
// Give 0 as "quals" to not use quality scores, and -1.0 as "cons_thres" to use the default
// consensus threshold when evaluating base votes.
char *get_consensus(char *align[], char *quals[], int n_seqs, int seq_len, double cons_thres,
                    int min_reads, char qual_thres, int gapped) {
  if (cons_thres == -1.0) {
    cons_thres = THRES_DEFAULT;
  }
  int **votes;
  if (quals == 0) {
    votes = get_votes_simple(align, n_seqs, seq_len);
  } else {
    votes = get_votes_qual(align, quals, n_seqs, seq_len, qual_thres);
  }
  char *consensus_gapped = build_consensus(votes, seq_len, cons_thres, min_reads);
  char *consensus;
  if (gapped) {
    consensus = consensus_gapped;
  } else {
    consensus = rm_gaps(consensus_gapped, seq_len);
  }
  free_votes(votes, seq_len);
  return consensus;
}


char *get_consensus_duplex(char *align1[], char *align2[], char *quals1[], char *quals2[],
                           int n_seqs1, int n_seqs2, int seq_len, double cons_thres,
                           int min_reads, char qual_thres, int gapped, char *method) {
  if (cons_thres == -1.0) {
    cons_thres = THRES_DEFAULT;
  }
  int **votes1;
  int **votes2;
  if (quals1 == 0 || quals2 == 0) {
    votes1 = get_votes_simple(align1, n_seqs1, seq_len);
    votes2 = get_votes_simple(align2, n_seqs2, seq_len);
  } else {
    votes1 = get_votes_qual(align1, quals1, n_seqs1, seq_len, qual_thres);
    votes2 = get_votes_qual(align2, quals2, n_seqs2, seq_len, qual_thres);
  }
  char *consensus_gapped;
  if (!strncmp(method, "freq", 4)) {
    consensus_gapped = build_consensus_duplex(votes1, votes2, seq_len, cons_thres);
  } else if (!strncmp(method, "iupac", 5)) {
    char *cons1 = build_consensus(votes1, seq_len, cons_thres, min_reads);
    char *cons2 = build_consensus(votes2, seq_len, cons_thres, min_reads);
    consensus_gapped = build_consensus_duplex_simple(cons1, cons2, 1);
  } else {
    return "";
  }
  char *consensus;
  if (gapped) {
    consensus = consensus_gapped;
  } else {
    consensus = rm_gaps(consensus_gapped, seq_len);
  }
  free_votes(votes1, seq_len);
  free_votes(votes2, seq_len);
  return consensus;
}


void get_gap_quals(char *quals) {
  int seq_len = strlen(quals);
  int *window = malloc(sizeof(int) * WIN_LEN * 2);
  int win_edge = init_gap_qual_window(window, quals, seq_len);
  print_window(window, win_edge);

  int i;
  char gap_qual;
  for (i = 0; i < seq_len; i++) {
    if (quals[i] == GAP_CHAR) {
      gap_qual = get_gap_qual(window);
      printf("gap %2d: %2d\n", i, gap_qual);
    } else {
      win_edge = push_qual(window, win_edge, quals, seq_len);
      print_window(window, win_edge);
    }
  }
}


int main(int argc, char *argv[]) {
  char **align = malloc(sizeof(char *) * (argc-1));

  int seq_len = INT_MAX;
  int i;
  for (i = 1; i < argc; i++) {
    if (strlen(argv[i]) < seq_len) {
      seq_len = strlen(argv[i]);
    }
    align[i-1] = argv[i];
  }

  if (argc <= 1) {
    return 1;
  }

  // get_gap_quals(align[0]);

  int **votes = get_votes_simple(align, argc-1, seq_len);
  char *consensus = build_consensus(votes, seq_len, THRES_DEFAULT, MIN_DEFAULT);
  print_votes(consensus, votes, seq_len);
  printf("%s\n", consensus);
  free_votes(votes, seq_len);

  return 0;
}
