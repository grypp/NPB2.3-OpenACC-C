/*--------------------------------------------------------------------
  
  NAS Parallel Benchmarks 2.3 OpenACC C versions - EP

  This benchmark is an OpenACC C version of the NPB EP code.
  
  The OpenACC C versions are derived from OpenMP C versions 
  in "NPB 2.3-omp" developed by NAS.

  Permission to use, copy, distribute and modify this software for any
  purpose with or without fee is hereby granted.
  This software is provided "as is" without express or implied warranty.
  
  Information on NAS Parallel Benchmarks 2.3 is available at:
  
           http://www.nas.nasa.gov/NAS/NPB/

--------------------------------------------------------------------*/
/*--------------------------------------------------------------------

  Author: P. O. Frederickson 
          D. H. Bailey
          A. C. Woo

  OpenMP C version: S. Satoh
  OpenACC C version: P. Makpaisit
  
--------------------------------------------------------------------*/

#include "npb-C.h"
#include "npbparams.h"

/* parameters */
#define	MK		16
#define	MM		(M - MK)
#define	NN		(1 << MM)
#define	NK		(1 << MK)
#define	NQ		10
#define EPSILON		1.0e-8
#define	A		1220703125.0
#define	S		271828183.0

#if defined(USE_POW)
#define r23 pow(0.5, 23.0)
#define r46 (r23*r23)
#define t23 pow(2.0, 23.0)
#define t46 (t23*t23)
#else
#define r23 (0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5)
#define r46 (r23*r23)
#define t23 (2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0)
#define t46 (t23*t23)
#endif

inline double RANDLC (double *x, double a) {
    double t1,t2,t3,t4,a1,a2,x1,x2,z;

    t1 = r23 * a;
    a1 = (int)t1;
    a2 = a - t23 * a1;

    t1 = r23 * (*x);
    x1 = (int)t1;
    x2 = (*x) - t23 * x1;
    t1 = a1 * x2 + a2 * x1;
    t2 = (int)(r23 * t1);
    z = t1 - t23 * t2;
    t3 = t23 * z + a2 * x2;
    t4 = (int)(r46 * t3);
    (*x) = t3 - t46 * t4;

    return (r46 * (*x));
}

/*--------------------------------------------------------------------
      program EMBAR
c-------------------------------------------------------------------*/
/*
c   This is the serial version of the APP Benchmark 1,
c   the "embarassingly parallel" benchmark.
c
c   M is the Log_2 of the number of complex pairs of uniform (0, 1) random
c   numbers.  MK is the Log_2 of the size of each batch of uniform random
c   numbers.  MK can be set for convenience on a given system, since it does
c   not affect the results.
*/
int main(int argc, char **argv) {

    double *x, **xx, *q, **qq;

    double Mops, t1, t2, t3, t4, x1, x2, sx, sy, tm, an, tt, gc;
    double dum[3] = { 1.0, 1.0, 1.0 };
    const int TRANSFER_X = 1;
    int np, nn, ierr, node, no_nodes, i, l, k, nit, ierrcode,
    no_large_nodes, np_add, k_offset, j;
    double loc_x,loc_t1,loc_t2,loc_t3,loc_t4;
    double loc_a1,loc_a2,loc_x1,loc_x2,loc_z;
    boolean verified;
    char size[13+1];	/* character*13 */
    
/*     Allocate working memory       */

    x = (double*) malloc(sizeof(double) * 2*NK);
    xx = (double**) malloc(sizeof(double*) * NN);
    xx[0] = (double*) malloc(sizeof(double) * NN * 2*NK);
    for (i = 1; i < NN; i++) xx[i] = xx[i-1] + (2*NK);
    q = (double*) malloc(sizeof(double) * NQ);
    qq = (double**) malloc(sizeof(double*) * NN);
    qq[0] = (double*) malloc(sizeof(double) * NN * NQ);
    for (i = 1; i < NN; i++) qq[i] = qq[i-1] + NQ;

/*
c   Because the size of the problem is too large to store in a 32-bit
c   integer for some classes, we put it into a string (for printing).
c   Have to strip off the decimal point put in there by the floating
c   point print statement (internal file)
*/

    printf("\n\n NAS Parallel Benchmarks 2.3 OpenACC C version"
	   " - EP Benchmark\n");
    sprintf(size, "%12.0f", pow(2.0, M+1));
    for (j = 13; j >= 1; j--) {
	if (size[j] == '.') size[j] = ' ';
    }
    printf(" Number of random numbers generated: %13s\n", size);

    verified = FALSE;

/*
c   Compute the number of "batches" of random number pairs generated 
c   per processor. Adjust if the number of processors does not evenly 
c   divide the total number
*/
    np = NN;

/*
c   Call the random number generator functions and initialize
c   the x-array to reduce the effects of paging on the timings.
c   Also, call all mathematical functions that are used. Make
c   sure these initializations cannot be eliminated as dead code.
*/
#pragma acc data create(qq[0:NN][0:NQ],x[0:2*NK],xx[0:NN][0:2*NK]) \
    copyout(q[0:NQ])
{
    vranlc(0, &(dum[0]), dum[1], &(dum[2]));
    dum[0] = randlc(&(dum[1]), dum[2]);
    for (i = 0; i < 2*NK; i++) x[i] = -1.0e99;
    Mops = log(sqrt(fabs(max(1.0, 1.0))));

    timer_clear(1);
    timer_clear(2);
    timer_clear(3);
    timer_start(1);

    vranlc(0, &t1, A, x);
    #pragma acc update device(x[0:2*NK])

/*   Compute AN = A ^ (2 * NK) (mod 2^46). */

    t1 = A;

    for ( i = 1; i <= MK+1; i++) {
      t2 = randlc(&t1, t1);
    }

    an = t1;
    tt = S;
    gc = 0.0;
    sx = 0.0;
    sy = 0.0;
    
    #pragma acc parallel loop
    for (k = 0; k < np; k++) {
      /* Initialize private q (qq) */
      #pragma acc loop
      for (i = 0; i < NQ; i++)
          qq[k][i] = 0.0;
      /* Initialize private x (xx)  */
      #pragma acc loop
      for (i = 0; i < 2*NK; i++)
          xx[k][i] = x[i];
    }
      
/*
c   Each instance of this loop may be performed independently. We compute
c   the k offsets separately to take into account the fact that some nodes
c   have more numbers to generate than others
*/
    k_offset = -1;

    double t1, t2, t3, t4, x1, x2;
    int kk, i, ik, l;
    double psx, psy;

    #pragma acc parallel loop reduction(+:sx,sy)
    for (k = 1; k <= np; k++) {
      kk = k_offset + k;
      t1 = S;
      t2 = an;

/*      Find starting seed t1 for this kk. */

      #pragma acc loop seq
      for (i = 1; i <= 100; i++) {
          ik = kk / 2;
          if (2 * ik != kk) t3 = RANDLC(&t1, t2);
          if (ik == 0) break;
          t3 = RANDLC(&t2, t2);
          kk = ik;
      }

/*      Compute uniform pseudorandom numbers. */

      loc_t1 = r23 * A;
      loc_a1 = (int)loc_t1;
      loc_a2 = A - t23 * loc_a1;
      loc_x = t1;

      #pragma acc loop seq
      for (i = 1; i <= 2*NK; i++) {
          loc_t1 = r23 * loc_x;
          loc_x1 = (int)loc_t1;
          loc_x2 = loc_x - t23 * loc_x1;
          loc_t1 = loc_a1 * loc_x2 + loc_a2 * loc_x1;
          loc_t2 = (int)(r23 * loc_t1);
          loc_z = loc_t1 - t23 * loc_t2;
          loc_t3 = t23 * loc_z + loc_a2 * loc_x2;
          loc_t4 = (int)(r46 * loc_t3);
          loc_x = loc_t3 - t46 * loc_t4;
          xx[k-1][i-1] = r46 * loc_x;
      }
      t1 = loc_x;

/*
c       Compute Gaussian deviates by acceptance-rejection method and 
c       tally counts in concentric square annuli.  This loop is not 
c       vectorizable.
*/
 
      psx = psy = 0.0;

      #pragma acc loop reduction(+:psx,psy)
      for ( i = 0; i < NK; i++) {
          x1 = 2.0 * xx[k-1][2*i] - 1.0;
          x2 = 2.0 * xx[k-1][2*i+1] - 1.0;
          t1 = pow2(x1) + pow2(x2);
          if (t1 <= 1.0) {
            t2 = sqrt(-2.0 * log(t1) / t1);
            t3 = (x1 * t2);             /* Xi */
            t4 = (x2 * t2);             /* Yi */
            l = max(fabs(t3), fabs(t4));
            qq[k-1][l] += 1.0;                      /* counts */
            psx = psx + t3;  /* sum of Xi */
            psy = psy + t4;               /* sum of Yi */
          }
      }

      sx += psx;
      sy += psy;
      
    }
    
/*      Reduce private qq to q          */
    #pragma acc parallel loop reduction(+:gc)
    for ( i = 0; i < NQ; i++ ) {
      double sumq = 0.0;
      #pragma acc loop reduction(+:sumq)
      for (k = 0; k < np; k++)
          sumq = sumq + qq[k][i];
      q[i] = sumq;
      gc += sumq;
    }

} /* end acc data */

    timer_stop(1);
    tm = timer_read(1);

    nit = 0;
    if (M == 24) {
	if((fabs((sx- (-3.247834652034740e3))/sx) <= EPSILON) &&
	   (fabs((sy- (-6.958407078382297e3))/sy) <= EPSILON)) {
	    verified = TRUE;
	}
    } else if (M == 25) {
	if ((fabs((sx- (-2.863319731645753e3))/sx) <= EPSILON) &&
	    (fabs((sy- (-6.320053679109499e3))/sy) <= EPSILON)) {
	    verified = TRUE;
	}
    } else if (M == 28) {
	if ((fabs((sx- (-4.295875165629892e3))/sx) <= EPSILON) &&
	    (fabs((sy- (-1.580732573678431e4))/sy) <= EPSILON)) {
	    verified = TRUE;
	}
    } else if (M == 30) {
	if ((fabs((sx- (4.033815542441498e4))/sx) <= EPSILON) &&
	    (fabs((sy- (-2.660669192809235e4))/sy) <= EPSILON)) {
	    verified = TRUE;
	}
    } else if (M == 32) {
	if ((fabs((sx- (4.764367927995374e4))/sx) <= EPSILON) &&
	    (fabs((sy- (-8.084072988043731e4))/sy) <= EPSILON)) {
	    verified = TRUE;
	}
    }

    Mops = pow(2.0, M+1)/tm/1000000.0;

    printf("EP Benchmark Results: \n"
	   "CPU Time = %10.4f\n"
	   "N = 2^%5d\n"
	   "No. Gaussian Pairs = %15.0f\n"
	   "Sums = %25.15e %25.15e\n"
	   "Counts:\n",
	   tm, M, gc, sx, sy);
    for (i = 0; i  <= NQ-1; i++) {
	printf("%3d %15.0f\n", i, q[i]);
    }
	  
    c_print_results("EP", CLASS, M+1, 0, 0, nit,
          tm, Mops, "Random numbers generated",
		  verified, NPBVERSION, COMPILETIME,
		  CS1, CS2, CS3, CS4, CS5, CS6, CS7);

    return 0;
}
