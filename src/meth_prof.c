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

/* Tools for methylation profile handling */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>

#include "meth_prof.h"

methylation_profile_t* init_methylation_profile() {
  /* Initialize a methylation profile from a tabix file */
  methylation_profile_t *methylation_profile = (methylation_profile_t*)calloc(1, sizeof(methylation_profile_t));
  methylation_profile->profile_buffer = (methylation_profile_position_t*)calloc(METHYLATION_PROFILE_BUFFER_SIZE, sizeof(methylation_profile_position_t));
  return methylation_profile;
}

void destroy_methylation_profile(methylation_profile_t* methylation_profile) {
  if(methylation_profile){
    if(methylation_profile->profile_buffer) {
      free(methylation_profile->profile_buffer);
    }
    if(methylation_profile->fp) {
      bgzf_close(methylation_profile->fp);
    }
    if(methylation_profile->tbx) {
      tbx_destroy(methylation_profile->tbx);
    }
    free(methylation_profile);
 }
}

void open_methylation_profile(methylation_profile_t *methylation_profile, char* filename) {
    tbx_t *tbx;
    BGZF *fp;
    if ((tbx = tbx_index_load(filename)) == NULL) {
      fprintf(stderr, "Error loading index file \"%s\"\n", filename);
      exit(EXIT_FAILURE);
    }
    if ((fp = bgzf_open(filename, "r")) == NULL) {
      fprintf(stderr, "Error loading file \"%s\"\n", filename);
      exit(EXIT_FAILURE);
    }
    methylation_profile->tbx = tbx;
    methylation_profile->fp = fp;
}

float methylation_profile_get_position(methylation_profile_t *methylation_profile, char* sequence, long position) {

  methylation_profile_position_t* prof_pos;

  /* If sequence changed, flush cache */
  if(strcmp(sequence, methylation_profile->sequence)) {
    strncpy(methylation_profile->sequence, sequence, 80);
    methylation_profile->last_position = 0;
    memset(methylation_profile->profile_buffer, 0, METHYLATION_PROFILE_BUFFER_SIZE*sizeof(methylation_profile_position_t));
  }

  /* If not cached value yet, read values */
  if(methylation_profile->last_position < position) {
    /* Read a chunk of data */
    kstring_t s;
    s.s = 0; s.l = s.m = 0;
    hts_itr_t *itr;
    char seq[200];
    long int pos;
    float meth_prob;

    /* Seek position */
    sprintf(seq, "%s:%ld", sequence, position);
    if ((itr = tbx_itr_querys(methylation_profile->tbx, seq)) == 0)  return -2;

    /* Read some lines, up to half the buffer */
    while (tbx_bgzf_itr_next(methylation_profile->fp, methylation_profile->tbx, itr, &s) >= 0
      && methylation_profile->last_position < (position + (METHYLATION_PROFILE_BUFFER_SIZE/2))) {
      //puts(s.s);
      sscanf(s.s, "%s\t%ld\t%f", seq, &pos, &meth_prob);
      prof_pos = &(methylation_profile->profile_buffer[pos%METHYLATION_PROFILE_BUFFER_SIZE]);
      prof_pos->position = pos;
      prof_pos->methylation_probability = meth_prob;
      methylation_profile->last_position = pos;
    }
    tbx_itr_destroy(itr);
    free(s.s);
  }

  /* Get cached value of methylation probability */
  prof_pos = &(methylation_profile->profile_buffer[position%METHYLATION_PROFILE_BUFFER_SIZE]);
  if(prof_pos->position == position) {
    return prof_pos->methylation_probability;
  }

  /* Value not available, return error */
  return -1;
}
