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
#ifndef METH_PROF_H
#define METH_PROF_H

#define METHYLATION_PROFILE_BUFFER_SIZE 4096

/* Probability of a methylation position */
typedef struct {
  long int position;
  float methylation_probability;
} methylation_profile_position_t;

/* Object to handle a methylation profile tabix */
typedef struct {
  tbx_t *tbx;
  BGZF *fp;
  char sequence[80];
  long int last_position;
  methylation_profile_position_t* profile_buffer;
} methylation_profile_t;

methylation_profile_t* init_methylation_profile(void);
void destroy_methylation_profile(methylation_profile_t* methylation_profile);
void open_methylation_profile(methylation_profile_t *methylation_profile, char* filename);
float methylation_profile_get_position(methylation_profile_t *methylation_profile, char* sequence, long position);

#endif /* METH_PROF_H */