/*
  The MIT License (MIT)

  Copyright (c) 2016  Eloi Casals Puig, <eloi.casals@cnag.crg.eu>
  Centre Nacional d'Analisi Genomica (CNAG)
  Barcelona, Spain

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

/* Simulator of bisulfite treatment on a methylated DNA sequence. */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <htslib/hfile.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>
#include <cram/sam_header.h>

#include "version.h"
#include "meth_prof.h"

#define SET_HIGH_NIBBLE(a, n) ((((n) & 0xf) << 4) | ((a) & 0xf))
#define SET_LOW_NIBBLE(a, n) (((a) & 0xf0) | ((n) & 0xf))

#define CLEAR_BIT(a, b) ((a)&~(1<<(b)))
#define SET_BIT(a, b) ((a)|(1<<(b)))
#define GET_BIT(a, b) (((a)>>(b))&0x1)

#define NT16_A  1
#define NT16_C  2
#define NT16_M  3
#define NT16_G  4
#define NT16_R  5
#define NT16_S  6
#define NT16_V  7
#define NT16_T  8
#define NT16_W  9
#define NT16_Y 10
#define NT16_H 11
#define NT16_K 12
#define NT16_D 13
#define NT16_B 14
#define NT16_N 15

#define NT16_SYMBOL "=ACMGRSVTWYHKDBN"

static int nt16_complement[] = {
/* = */      0,
  /* NT16_A */ NT16_T,
  /* NT16_C */ NT16_G,
  /* NT16_M */ NT16_K,
  /* NT16_G */ NT16_C,
  /* NT16_R */ NT16_Y,
  /* NT16_S */ NT16_S,
  /* NT16_V */ NT16_B,
  /* NT16_T */ NT16_A,
  /* NT16_W */ NT16_W,
  /* NT16_Y */ NT16_R,
  /* NT16_H */ NT16_D,
  /* NT16_K */ NT16_M,
  /* NT16_D */ NT16_H,
  /* NT16_B */ NT16_V,
  /* NT16_N */ NT16_N
  };


/* Settings to use when processing an aligned read */
typedef struct settings {
  float lambda;      /*  C->T probability */
  float tau;         /* mC->T probability */
  float methylation; /* methylation probability (position independent profile for now) */
  unsigned seed;     /* Seed for the random number generator */
  methylation_profile_t *methylation_profile; /* Handler for methylation probabilities file */
} settings_t;


/* Set a nucleotide in the sequence */
static inline void bam_set_seqi(uint8_t *s, int i, int n) {
  int i_offset = i&1;
  int i_base = i>>1;
  if(i_offset==1) {
    s[i_base] = SET_LOW_NIBBLE(s[i_base], n);
  } else {
    s[i_base] = SET_HIGH_NIBBLE(s[i_base], n);
  }
}

/* Convert a Phred quality score to error probability
   p = 10^(-q/10) 
*/
static inline float phred2prob(float q) {
  return(powf(10.0, (-q/10.0)));
}

/* Weighted random selection.
   Given weighted probabilities for a set of elements, selects one randomly 
   according to its weights.
*/
int weighted_random_selection(int size, int* elem, float* weights) {
  int i = 0;
  double total_weight = 0.0;
  double random_number, rn_start;
    
  /* Normalize random range */
  for(i=0; i<size; i++) {
    total_weight += weights[i];
  }
  //random_number = rn_start = total_weight*(float)rand()/(double)(RAND_MAX+1);
  //random_number = rn_start = (double)rand()*total_weight/((double)RAND_MAX+1.0);
  random_number = rn_start = drand48()*total_weight;

  /* Pick element randomly */
  for(i=0; i<size; i++) {
    if(random_number < weights[i]) {
      return elem[i];
    }
    random_number -= weights[i];
  }
  // return elem[size-1];
  /* Should never get here */
  fprintf(stderr, "Error with random selection of bases.\n#Elements = %d\nWeights are: ", size);
  for(i=0; i<size; i++) {
    fprintf(stderr, "\t%f", weights[i]);
  }
  fprintf(stderr, "\nRandom number: %f\n", rn_start);
  assert(0);
}

#define bam_is_paired(b) (((b)->core.flag&BAM_FPAIRED) != 0)
#define bam_is_proper_pair(b) (((b)->core.flag&BAM_FPROPER_PAIR) != 0)
#define bam_get_target_name(h, b) ((h)->target_name[(b)->core.tid])

/*
Find the position in the reference for a position from the sequence,
parsing the CIGAR string.
*/
int bam_cigar2ref_position(int query_position, int n_cigar, const uint32_t *cigar)
{
    int type_op, len_op;
    int k, l, q;
    l = 0;
    q = 0;
    for (k = 0; k < n_cigar; ++k) {
        type_op = bam_cigar_type(bam_cigar_op(cigar[k]));
        len_op = bam_cigar_oplen(cigar[k]);
        if (type_op&1) {
          /* seq query increment */
          if( q+len_op >= query_position) {
            len_op = query_position-q;
          } 
          q += len_op;
        }
        if (type_op&2) {
          /* ref increment */
          l += len_op;
        }
        if(q == query_position) {
          /* reached position of interest in the query, 
          return position in reference */
          return l;
        }
    }
    return l;
}

/* Apply a statistical model of conversion from cytosines into thymines
   after a bisulfite treatment on the DNA sequence.
   Processes a single read alignment.
   Returns 0 to indicate read should be output 1 otherwise.
*/
static int process_aln(const bam_hdr_t *h, bam1_t *b, settings_t *settings)
{
  //printf("QNAME: %s \tPOS: %d:%d \tLEN: %d \t \n", bam_get_qname(b), b->core.tid, b->core.pos, b->core.l_qseq);
  //printf("QNAME: %s \tPOS: %s:%d \tLEN: %d \t ENDPOS: %d \n", bam_get_qname(b), bam_get_target_name(h,b), b->core.pos, b->core.l_qseq, bam_endpos(b));
  //printf("\tReverse: %d\tPaired in sequencing: %d\t Paired in mapping: %d\n", bam_is_rev(b), bam_is_paired(b), bam_is_proper_pair(b));
  static int warn_nt = 0;
  char* region_name = bam_get_target_name(h,b);
  uint8_t *seq = bam_get_seq(b);
  int32_t len = b->core.l_qseq;
  int i;

  uint8_t *qual = bam_get_qual(b);

  int is_reverse = bam_is_rev(b);

  // for(i=0; i<len; i++) {
  //   printf("%c", NT16_SYMBOL[bam_seqi(seq, i)]);
  // }
  // printf("\n");

  // for(i=0; i<len; i++) {
  //   printf("(%d,'%c') ", qual[i], qual[i]+33);
  //   printf("%0.4f ", phred2prob((float)qual[i]));
  // }
  // printf("\n");

  /* Process each sequence nucleotide */
  for(i=0; i<len; i++) {
    uint8_t nt = bam_seqi(seq, i);

    /* Evaluate probability of observing a certain base after treatment.
       Decide if it was originally a methylated cytosine or 
       a regular cytosine, and if it is converted.
       This considers sequencing errors, methylated position probability,
       and bisulfite conversion efficiency.
    */
    float A2, C2, T2, G2;
    float A3, C3, CM3, T3, G3;
    float A4, C4, T4, G4;
    float A5, C5, T5, G5;
    float seq_error = phred2prob((float)qual[i]);
    float lambda = settings->lambda;
    float tau = settings->tau;
    float methylation = settings->methylation;

    /* Use methylation profile if available */
    if(settings->methylation_profile && !(b->core.flag & BAM_FUNMAP)) {
      int ref_pos = bam_cigar2ref_position(i, b->core.n_cigar, bam_get_cigar(b));
      float meth_prob = methylation_profile_get_position(settings->methylation_profile, region_name, b->core.pos+ref_pos);
      if (meth_prob >= 0) {
        methylation = meth_prob;
      }
    }

    /* Complement the reverse read */
    if(is_reverse) {
      nt = nt16_complement[nt];
    }

    /* Set probability of real base (2) from the observed base (1)
       For Observed A1: A2 = A1*(1-error); C2 = A1*(error/3); T2 = A1*(error/3); G2 = A1*(error/3);
       For Observed C1: A2 = C1*(error/3); C2 = C1*(1-error); T2 = C1*(error/3); G2 = C1*(error/3);
       For Observed T1: A2 = T1*(error/3); C2 = T1*(error/3); T2 = T1*(1-error); G2 = T1*(error/3);
       For Observed G1: A2 = G1*(error/3); C2 = G1*(error/3); T2 = G1*(error/3); G2 = G1*(1-error);
    */
    A2 = C2 = T2 = G2 = (seq_error/3);
    switch(nt) {
    case NT16_A: A2 = (1-seq_error); break;
    case NT16_C: C2 = (1-seq_error); break;
    case NT16_T: T2 = (1-seq_error); break;
    case NT16_G: G2 = (1-seq_error); break;
    default:
      /* Unsupported base, skip without processing */
      if(!GET_BIT(warn_nt, nt)) {
        /* Warn once per non-supported IUPAC code */
        warn_nt = SET_BIT(warn_nt, nt);
        fprintf(stderr, "Bases with symbol \"%c\" are ignored. Processing only A, C, T or G symbols.\n", NT16_SYMBOL[nt]);
      }
      continue;
      break;
    }

    /* Include methylation probability (3).
       CM3 =     prob_meth*C2;
       C3  = (1-prob_meth)*C2;
    */
    A3 = A2;
    CM3 = methylation*C2;
    C3 = (1-methylation)*C2;
    T3 = T2;
    G3 = G2;

    /* Bisulfite treatment conversion probability (4).
       C4 = (1-tau)*CM3 + (1-lambda)*C3;
       T4 =     tau*CM3 +     lambda*C3 + T3;
    */
    A4 = A3;
    C4 = (1-tau)*CM3+(1-lambda)*C3;
    T4 = tau*CM3+lambda*C3+T3;
    G4 = G3;

    /* Sequencing error of the treated sequence.
       Probability of observed base (5) after treatment (4).
       A5 = (1-seq_error)*A4 + (seq_error/3)*C4 + (seq_error/3)*T4 + (seq_error/3)*G4;
       C5 = (seq_error/3)*A4 + (1-seq_error)*C4 + (seq_error/3)*T4 + (seq_error/3)*G4;
       T5 = (seq_error/3)*A4 + (seq_error/3)*C4 + (1-seq_error)*T4 + (seq_error/3)*G4;
       G5 = (seq_error/3)*A4 + (seq_error/3)*C4 + (seq_error/3)*T4 + (1-seq_error)*G4;
    */
    A5 = (1-seq_error)*A4+(seq_error/3)*(C4+T4+G4);
    C5 = (1-seq_error)*C4+(seq_error/3)*(A4+T4+G4);
    T5 = (1-seq_error)*T4+(seq_error/3)*(A4+C4+G4);
    G5 = (1-seq_error)*G4+(seq_error/3)*(A4+C4+T4);

#define BASES_SIZE 4
    static const int BASES[BASES_SIZE] = {NT16_A, NT16_C, NT16_T, NT16_G};
    float weights[BASES_SIZE] = {A5, C5, T5, G5};
    int new_base = weighted_random_selection(BASES_SIZE, (int*)BASES, weights);
    //printf("Pos: %d\tMeth: %f\t", b->core.pos+i, methylation);
    //printf("nt: %c  A: %f  C: %f  T: %f  G: %f  TOTAL: %f  %c-->%c\n", NT16_SYMBOL[nt], A5, C5, T5, G5,  A5+C5+T5+G5, NT16_SYMBOL[nt], NT16_SYMBOL[new_base]);
#undef BASES_SIZE

    if(nt != new_base) {
      /* Complement the reverse read */
      if(is_reverse) {
	       new_base = nt16_complement[new_base];
      }

      /* Update read with new base */
      bam_set_seqi(seq, i, new_base);
    }
  }
  return 0;
}

static inline int check_sam_write1(samFile *fp, const bam_hdr_t *h, const bam1_t *b, const char *fname, int *retp)
{
  int r = sam_write1(fp, h, b);
  if (r >= 0) return r;

  if (fname) fprintf(stderr, "writing to \"%s\" failed.\n", fname);
  else fprintf(stderr, "writing to standard output failed.\n");

  *retp = EXIT_FAILURE;
  return r;
}


#include "docopt.c"

void report_settings(DocoptArgs *args, settings_t *settings) {
  printf("bs_sim version '%s'\n", VERSION);
  printf("Input: %s\n", args->input);
  printf("Output: %s\n", args->output);
  printf("Methylation Profile: %s\n", args->methylation_profile);
  printf("Lambda = %f\t", settings->lambda);
  printf("Tau = %f\t", settings->tau);
  printf("Methylation = %f\t", settings->methylation);
  printf("Random seed = %u\t", settings->seed);
  printf("\n");
}

int main(int argc, char* argv[]) {

  DocoptArgs args = docopt(argc, argv, /* help */ 1, /* version */ VERSION);

  int ret = 0;
  int64_t count = 0;
  char *fn_in = 0;
  char *fn_out = 0;
  samFile *in = 0;
  samFile *out = 0;
  bam_hdr_t *header = NULL;
  char out_mode[5];
  //char *out_format = "";
  htsFormat ga_in = {0};
  htsFormat ga_out = {0};
  int n_threads = 0;

  char *arg_list = stringify_argv(argc, argv);

  if(!args.input || !args.output) {
    fprintf(stderr, "%s\n\n", args.usage_pattern);
    return(EXIT_FAILURE);
  }
  fn_in = args.input;
  fn_out = args.output;

  settings_t settings = {
    .lambda = atof(args.lambda),
    .tau = atof(args.tau),
    .methylation = atof(args.methylation),
    .methylation_profile = NULL,
  };

  /* Intialize random number generator */
  if(args.seed) {
    settings.seed = atoi(args.seed);
    
  } else {
    settings.seed = time(NULL);
  }
  //srand(settings.seed);
  srand48(settings.seed);
  
  /* Report configuration */
  report_settings(&args, &settings);

  strcpy(out_mode, "w");

  /* File format auto-detection first */
  sam_open_mode(out_mode+1, fn_out, NULL);

  /* open file handlers */
  if ((in = sam_open_format(fn_in, "r", &ga_in)) == 0) {
    fprintf(stderr, "Failed to open \"%s\" for reading.\n", fn_in);
    return(EXIT_FAILURE);
  }
  if ((out = sam_open_format(fn_out? fn_out : "-", out_mode, &ga_out)) == 0) {
    fprintf(stderr, "Failed to open \"%s\" for writing.\n", fn_out? fn_out : "standard output");
    return(EXIT_FAILURE);
  }

  /* Check settings */
  if(settings.lambda < 0.0 || settings.lambda > 1.0) {
    fprintf(stderr, "Error: Lambda must be in the interval [0.0, 1.0]\n");
    return(EXIT_FAILURE);
  }
  if(settings.tau < 0.0 || settings.tau > 1.0) {
    fprintf(stderr, "Error: Tau must be in the interval [0.0, 1.0]\n");
    return(EXIT_FAILURE);
  }
  if(settings.methylation < 0.0 || settings.methylation > 1.0) {
    fprintf(stderr, "Error: Methylation must be in the interval [0.0, 1.0]\n");
    return(EXIT_FAILURE);
  }

  /* Open methylation profile (random access) */
  if(args.methylation_profile) {
    methylation_profile_t* meth_prof;
    meth_prof = init_methylation_profile();
    open_methylation_profile(meth_prof, args.methylation_profile);
    settings.methylation_profile = meth_prof;
  }

  /* Get header */
  if ((header = sam_hdr_read(in)) == 0) {
    fprintf(stderr, "Fail to read the header from \"%s\".\n", fn_in);
    return(EXIT_FAILURE);
  }

  /* Add info line into the header */
  SAM_hdr *sh = sam_hdr_parse_(header->text, header->l_text);
  if (sam_hdr_add_PG(sh, "bs_sim",
                     "VN", VERSION,
                     "CL", arg_list,
                     NULL) != 0)
      return -1;

  free(header->text);
  header->text = strdup(sam_hdr_str(sh));
  header->l_text = sam_hdr_length(sh);
  if (!header->text)
      return -1;
  sam_hdr_free(sh);

  /* Write header */
  if (sam_hdr_write(out, header) != 0) {
    fprintf(stderr, "Failed to write the SAM header\n");
    return(EXIT_FAILURE);
  }

  /* Set threads to help compression */
  n_threads = sysconf(_SC_NPROCESSORS_ONLN); /* Get number of cores */
  if (n_threads > 1) {
    if (out) {
      hts_set_threads(out, n_threads);
    }
  }

  /* Parse SAM/BAM file */
  bam1_t *b = bam_init1();
  int r;
  while ((r = sam_read1(in, header, b)) >= 0) { // read one alignment from `in'
    if (!process_aln(header, b, &settings)) {
      if (check_sam_write1(out, header, b, fn_out, &ret) < 0) break;
      count++;
    }
  }
  if (r < -1) {
    fprintf(stderr, "Truncated input file.\n");
    ret = 1;
  }
  bam_destroy1(b);

  if (settings.methylation_profile) destroy_methylation_profile(settings.methylation_profile);
  if (header) bam_hdr_destroy(header);
  if (in) hts_close(in);
  if (out) hts_close(out);

  return ret;
}
