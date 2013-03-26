/*--------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  header.h
c
c---------------------------------------------------------------------
c-------------------------------------------------------------------*/
 

/*--------------------------------------------------------------------
c The following include file is generated automatically by the
c "setparams" utility. It defines 
c      maxcells:      the square root of the maximum number of processors
c      problem_size:  12, 64, 102, 162 (for class T, A, B, C)
c      dt_default:    default time step for this problem size if no
c                     config file
c      niter_default: default number of iterations for this problem size
--------------------------------------------------------------------*/

#include "npbparams.h"

#define	AA		0
#define BB		1
#define CC		2
#define	BLOCK_SIZE	5

/* COMMON block: global */
static int grid_points[3];	/* grid_ponts(1:3) */

/* COMMON block: constants */
static double tx1, tx2, tx3, ty1, ty2, ty3, tz1, tz2, tz3;
static double dx1, dx2, dx3, dx4, dx5;
static double dy1, dy2, dy3, dy4, dy5;
static double dz1, dz2, dz3, dz4, dz5;
static double dssp, dt;
static double ce[5][13];	/* ce(5,13) */
static double dxmax, dymax, dzmax;
static double xxcon1, xxcon2, xxcon3, xxcon4, xxcon5;
static double dx1tx1, dx2tx1, dx3tx1, dx4tx1, dx5tx1;
static double yycon1, yycon2, yycon3, yycon4, yycon5;
static double dy1ty1, dy2ty1, dy3ty1, dy4ty1, dy5ty1;
static double zzcon1, zzcon2, zzcon3, zzcon4, zzcon5;
static double dz1tz1, dz2tz1, dz3tz1, dz4tz1, dz5tz1;
static double dnxm1, dnym1, dnzm1, c1c2, c1c5, c3c4, c1345;
static double conz1, c1, c2, c3, c4, c5, c4dssp, c5dssp, dtdssp;
static double dttx1, dttx2, dtty1, dtty2, dttz1, dttz2;
static double c2dttx1, c2dtty1, c2dttz1, comz1, comz4, comz5, comz6;
static double c3c4tx3, c3c4ty3, c3c4tz3, c2iv, con43, con16;

#define	IMAX	PROBLEM_SIZE
#define	JMAX	PROBLEM_SIZE
#define	KMAX	PROBLEM_SIZE

/*
c   to improve cache performance, grid dimensions padded by 1 
c   for even number sizes only.
*/

/* COMMON block: fields */
static double us[IMAX/2*2+1][JMAX/2*2+1][KMAX/2*2+1];
static double vs[IMAX/2*2+1][JMAX/2*2+1][KMAX/2*2+1];
static double ws[IMAX/2*2+1][JMAX/2*2+1][KMAX/2*2+1];
static double qs[IMAX/2*2+1][JMAX/2*2+1][KMAX/2*2+1];
static double rho_i[IMAX/2*2+1][JMAX/2*2+1][KMAX/2*2+1];
static double square[IMAX/2*2+1][JMAX/2*2+1][KMAX/2*2+1];
static double forcing[IMAX/2*2+1][JMAX/2*2+1][KMAX/2*2+1][5+1];
static double u[(IMAX+1)/2*2+1][(JMAX+1)/2*2+1][(KMAX+1)/2*2+1][5];
static double rhs[IMAX/2*2+1][JMAX/2*2+1][KMAX/2*2+1][5];
static double lhs[IMAX/2*2+1][JMAX/2*2+1][KMAX/2*2+1][3][5][5];

/* COMMON block: work_1d */
/*static double cuf[PROBLEM_SIZE];
static double q[PROBLEM_SIZE];
static double ue[PROBLEM_SIZE][5];
static double buf[PROBLEM_SIZE][5];*/

/*
c   to improve cache performance, grid dimensions (first two for these
c   to arrays) padded by 1 for even number sizes only.
*/

/* COMMON block: work_lhs */
static double fjac[IMAX/2*2+1][JMAX/2*2+1][KMAX-1+1][5][5];
/* fjac(5, 5, 0:IMAX/2*2, 0:JMAX/2*2, 0:KMAX-1) */
static double njac[IMAX/2*2+1][JMAX/2*2+1][KMAX-1+1][5][5];
/* njac(5, 5, 0:IMAX/2*2, 0:JMAX/2*2, 0:KMAX-1) */


/* this function returns the exact solution at point xi, eta, zeta */
#define EXACT_SOLUTION(xi,eta,zeta,dtemp) \
{ \
    int m; \
 \
    for (m = 0; m < 5; m++) { \
      (dtemp)[m] =  ce[m][0] + \
        xi*(ce[m][1] + xi*(ce[m][4] + xi*(ce[m][7] \
                      + xi*ce[m][10]))) + \
        eta*(ce[m][2] + eta*(ce[m][5] + eta*(ce[m][8] \
                         + eta*ce[m][11])))+ \
        zeta*(ce[m][3] + zeta*(ce[m][6] + zeta*(ce[m][9] + \
                            zeta*ce[m][12]))); \
    } \
}

/* subtracts bvec=bvec - ablock*avec */
/*--------------------------------------------------------------------
c            rhs(i,ic,jc,kc,ccell) = rhs(i,ic,jc,kc,ccell)
c     $           - lhs[i,1,ablock,ia,ja,ka,acell)*
c-------------------------------------------------------------------*/
#define MATVEC_SUB(ablock,avec,bvec) \
{ \
  int l; \
 \
  for (l = 0; l < 5; l++) { \
    bvec[l] = bvec[l] - ablock[l][0]*avec[0] \
      - ablock[l][1]*avec[1] \
      - ablock[l][2]*avec[2] \
      - ablock[l][3]*avec[3] \
      - ablock[l][4]*avec[4]; \
  } \
}

/* subtracts a(i,j,k) X b(i,j,k) from c(i,j,k) */
#define MATMUL_SUB(ablock,bblock,cblock) \
{ \
  int l; \
 \
  for (l = 0; l < 5; l++) { \
    cblock[0][l] = cblock[0][l] - ablock[0][0]*bblock[0][l] \
      - ablock[0][1]*bblock[1][l] \
      - ablock[0][2]*bblock[2][l] \
      - ablock[0][3]*bblock[3][l] \
      - ablock[0][4]*bblock[4][l]; \
    cblock[1][l] = cblock[1][l] - ablock[1][0]*bblock[0][l] \
      - ablock[1][1]*bblock[1][l] \
      - ablock[1][2]*bblock[2][l] \
      - ablock[1][3]*bblock[3][l] \
      - ablock[1][4]*bblock[4][l]; \
    cblock[2][l] = cblock[2][l] - ablock[2][0]*bblock[0][l] \
      - ablock[2][1]*bblock[1][l] \
      - ablock[2][2]*bblock[2][l] \
      - ablock[2][3]*bblock[3][l] \
      - ablock[2][4]*bblock[4][l]; \
    cblock[3][l] = cblock[3][l] - ablock[3][0]*bblock[0][l] \
      - ablock[3][1]*bblock[1][l] \
      - ablock[3][2]*bblock[2][l] \
      - ablock[3][3]*bblock[3][l] \
      - ablock[3][4]*bblock[4][l]; \
    cblock[4][l] = cblock[4][l] - ablock[4][0]*bblock[0][l] \
      - ablock[4][1]*bblock[1][l] \
      - ablock[4][2]*bblock[2][l] \
      - ablock[4][3]*bblock[3][l] \
      - ablock[4][4]*bblock[4][l]; \
  } \
}

#define BINVRHS(lhs,r) \
{ \
  static double pivot, coeff; \
 \
  pivot = 1.00/lhs[0][0]; \
  lhs[0][1] = lhs[0][1]*pivot; \
  lhs[0][2] = lhs[0][2]*pivot; \
  lhs[0][3] = lhs[0][3]*pivot; \
  lhs[0][4] = lhs[0][4]*pivot; \
  r[0]   = r[0]  *pivot; \
 \
  coeff = lhs[1][0]; \
  lhs[1][1]= lhs[1][1] - coeff*lhs[0][1]; \
  lhs[1][2]= lhs[1][2] - coeff*lhs[0][2]; \
  lhs[1][3]= lhs[1][3] - coeff*lhs[0][3]; \
  lhs[1][4]= lhs[1][4] - coeff*lhs[0][4]; \
  r[1]   = r[1]   - coeff*r[0]; \
 \
  coeff = lhs[2][0]; \
  lhs[2][1]= lhs[2][1] - coeff*lhs[0][1]; \
  lhs[2][2]= lhs[2][2] - coeff*lhs[0][2]; \
  lhs[2][3]= lhs[2][3] - coeff*lhs[0][3]; \
  lhs[2][4]= lhs[2][4] - coeff*lhs[0][4]; \
  r[2]   = r[2]   - coeff*r[0]; \
 \
  coeff = lhs[3][0]; \
  lhs[3][1]= lhs[3][1] - coeff*lhs[0][1]; \
  lhs[3][2]= lhs[3][2] - coeff*lhs[0][2]; \
  lhs[3][3]= lhs[3][3] - coeff*lhs[0][3]; \
  lhs[3][4]= lhs[3][4] - coeff*lhs[0][4]; \
  r[3]   = r[3]   - coeff*r[0]; \
 \
  coeff = lhs[4][0]; \
  lhs[4][1]= lhs[4][1] - coeff*lhs[0][1]; \
  lhs[4][2]= lhs[4][2] - coeff*lhs[0][2]; \
  lhs[4][3]= lhs[4][3] - coeff*lhs[0][3]; \
  lhs[4][4]= lhs[4][4] - coeff*lhs[0][4]; \
  r[4]   = r[4]   - coeff*r[0]; \
 \
  pivot = 1.00/lhs[1][1]; \
  lhs[1][2] = lhs[1][2]*pivot; \
  lhs[1][3] = lhs[1][3]*pivot; \
  lhs[1][4] = lhs[1][4]*pivot; \
  r[1]   = r[1]  *pivot; \
 \
  coeff = lhs[0][1]; \
  lhs[0][2]= lhs[0][2] - coeff*lhs[1][2]; \
  lhs[0][3]= lhs[0][3] - coeff*lhs[1][3]; \
  lhs[0][4]= lhs[0][4] - coeff*lhs[1][4]; \
  r[0]   = r[0]   - coeff*r[1]; \
 \
  coeff = lhs[2][1]; \
  lhs[2][2]= lhs[2][2] - coeff*lhs[1][2]; \
  lhs[2][3]= lhs[2][3] - coeff*lhs[1][3]; \
  lhs[2][4]= lhs[2][4] - coeff*lhs[1][4]; \
  r[2]   = r[2]   - coeff*r[1]; \
 \
  coeff = lhs[3][1]; \
  lhs[3][2]= lhs[3][2] - coeff*lhs[1][2]; \
  lhs[3][3]= lhs[3][3] - coeff*lhs[1][3]; \
  lhs[3][4]= lhs[3][4] - coeff*lhs[1][4]; \
  r[3]   = r[3]   - coeff*r[1]; \
 \
  coeff = lhs[4][1]; \
  lhs[4][2]= lhs[4][2] - coeff*lhs[1][2]; \
  lhs[4][3]= lhs[4][3] - coeff*lhs[1][3]; \
  lhs[4][4]= lhs[4][4] - coeff*lhs[1][4]; \
  r[4]   = r[4]   - coeff*r[1]; \
 \
  pivot = 1.00/lhs[2][2]; \
  lhs[2][3] = lhs[2][3]*pivot; \
  lhs[2][4] = lhs[2][4]*pivot; \
  r[2]   = r[2]  *pivot; \
 \
  coeff = lhs[0][2]; \
  lhs[0][3]= lhs[0][3] - coeff*lhs[2][3]; \
  lhs[0][4]= lhs[0][4] - coeff*lhs[2][4]; \
  r[0]   = r[0]   - coeff*r[2]; \
 \
  coeff = lhs[1][2]; \
  lhs[1][3]= lhs[1][3] - coeff*lhs[2][3]; \
  lhs[1][4]= lhs[1][4] - coeff*lhs[2][4]; \
  r[1]   = r[1]   - coeff*r[2]; \
 \
  coeff = lhs[3][2]; \
  lhs[3][3]= lhs[3][3] - coeff*lhs[2][3]; \
  lhs[3][4]= lhs[3][4] - coeff*lhs[2][4]; \
  r[3]   = r[3]   - coeff*r[2]; \
 \
  coeff = lhs[4][2]; \
  lhs[4][3]= lhs[4][3] - coeff*lhs[2][3]; \
  lhs[4][4]= lhs[4][4] - coeff*lhs[2][4]; \
  r[4]   = r[4]   - coeff*r[2]; \
 \
  pivot = 1.00/lhs[3][3]; \
  lhs[3][4] = lhs[3][4]*pivot; \
  r[3]   = r[3]  *pivot; \
 \
  coeff = lhs[0][3]; \
  lhs[0][4]= lhs[0][4] - coeff*lhs[3][4]; \
  r[0]   = r[0]   - coeff*r[3]; \
 \
  coeff = lhs[1][3]; \
  lhs[1][4]= lhs[1][4] - coeff*lhs[3][4]; \
  r[1]   = r[1]   - coeff*r[3]; \
 \
  coeff = lhs[2][3]; \
  lhs[2][4]= lhs[2][4] - coeff*lhs[3][4]; \
  r[2]   = r[2]   - coeff*r[3]; \
 \
  coeff = lhs[4][3]; \
  lhs[4][4]= lhs[4][4] - coeff*lhs[3][4]; \
  r[4]   = r[4]   - coeff*r[3]; \
 \
  pivot = 1.00/lhs[4][4]; \
  r[4]   = r[4]  *pivot; \
 \
  coeff = lhs[0][4]; \
  r[0]   = r[0]   - coeff*r[4]; \
 \
  coeff = lhs[1][4]; \
  r[1]   = r[1]   - coeff*r[4]; \
 \
  coeff = lhs[2][4]; \
  r[2]   = r[2]   - coeff*r[4]; \
 \
  coeff = lhs[3][4]; \
  r[3]   = r[3]   - coeff*r[4]; \
}

#define BINVCRHS(lhs,c,r) \
{ \
  static double pivot, coeff; \
 \
  pivot = 1.00/lhs[0][0]; \
  lhs[0][1] = lhs[0][1]*pivot; \
  lhs[0][2] = lhs[0][2]*pivot; \
  lhs[0][3] = lhs[0][3]*pivot; \
  lhs[0][4] = lhs[0][4]*pivot; \
  c[0][0] = c[0][0]*pivot; \
  c[0][1] = c[0][1]*pivot; \
  c[0][2] = c[0][2]*pivot; \
  c[0][3] = c[0][3]*pivot; \
  c[0][4] = c[0][4]*pivot; \
  r[0]   = r[0]  *pivot; \
 \
  coeff = lhs[1][0]; \
  lhs[1][1]= lhs[1][1] - coeff*lhs[0][1]; \
  lhs[1][2]= lhs[1][2] - coeff*lhs[0][2]; \
  lhs[1][3]= lhs[1][3] - coeff*lhs[0][3]; \
  lhs[1][4]= lhs[1][4] - coeff*lhs[0][4]; \
  c[1][0] = c[1][0] - coeff*c[0][0]; \
  c[1][1] = c[1][1] - coeff*c[0][1]; \
  c[1][2] = c[1][2] - coeff*c[0][2]; \
  c[1][3] = c[1][3] - coeff*c[0][3]; \
  c[1][4] = c[1][4] - coeff*c[0][4]; \
  r[1]   = r[1]   - coeff*r[0]; \
 \
  coeff = lhs[2][0]; \
  lhs[2][1]= lhs[2][1] - coeff*lhs[0][1]; \
  lhs[2][2]= lhs[2][2] - coeff*lhs[0][2]; \
  lhs[2][3]= lhs[2][3] - coeff*lhs[0][3]; \
  lhs[2][4]= lhs[2][4] - coeff*lhs[0][4]; \
  c[2][0] = c[2][0] - coeff*c[0][0]; \
  c[2][1] = c[2][1] - coeff*c[0][1]; \
  c[2][2] = c[2][2] - coeff*c[0][2]; \
  c[2][3] = c[2][3] - coeff*c[0][3]; \
  c[2][4] = c[2][4] - coeff*c[0][4]; \
  r[2]   = r[2]   - coeff*r[0]; \
 \
  coeff = lhs[3][0]; \
  lhs[3][1]= lhs[3][1] - coeff*lhs[0][1]; \
  lhs[3][2]= lhs[3][2] - coeff*lhs[0][2]; \
  lhs[3][3]= lhs[3][3] - coeff*lhs[0][3]; \
  lhs[3][4]= lhs[3][4] - coeff*lhs[0][4]; \
  c[3][0] = c[3][0] - coeff*c[0][0]; \
  c[3][1] = c[3][1] - coeff*c[0][1]; \
  c[3][2] = c[3][2] - coeff*c[0][2]; \
  c[3][3] = c[3][3] - coeff*c[0][3]; \
  c[3][4] = c[3][4] - coeff*c[0][4]; \
  r[3]   = r[3]   - coeff*r[0]; \
 \
  coeff = lhs[4][0]; \
  lhs[4][1]= lhs[4][1] - coeff*lhs[0][1]; \
  lhs[4][2]= lhs[4][2] - coeff*lhs[0][2]; \
  lhs[4][3]= lhs[4][3] - coeff*lhs[0][3]; \
  lhs[4][4]= lhs[4][4] - coeff*lhs[0][4]; \
  c[4][0] = c[4][0] - coeff*c[0][0]; \
  c[4][1] = c[4][1] - coeff*c[0][1]; \
  c[4][2] = c[4][2] - coeff*c[0][2]; \
  c[4][3] = c[4][3] - coeff*c[0][3]; \
  c[4][4] = c[4][4] - coeff*c[0][4]; \
  r[4]   = r[4]   - coeff*r[0]; \
 \
  pivot = 1.00/lhs[1][1]; \
  lhs[1][2] = lhs[1][2]*pivot; \
  lhs[1][3] = lhs[1][3]*pivot; \
  lhs[1][4] = lhs[1][4]*pivot; \
  c[1][0] = c[1][0]*pivot; \
  c[1][1] = c[1][1]*pivot; \
  c[1][2] = c[1][2]*pivot; \
  c[1][3] = c[1][3]*pivot; \
  c[1][4] = c[1][4]*pivot; \
  r[1]   = r[1]  *pivot; \
 \
  coeff = lhs[0][1]; \
  lhs[0][2]= lhs[0][2] - coeff*lhs[1][2]; \
  lhs[0][3]= lhs[0][3] - coeff*lhs[1][3]; \
  lhs[0][4]= lhs[0][4] - coeff*lhs[1][4]; \
  c[0][0] = c[0][0] - coeff*c[1][0]; \
  c[0][1] = c[0][1] - coeff*c[1][1]; \
  c[0][2] = c[0][2] - coeff*c[1][2]; \
  c[0][3] = c[0][3] - coeff*c[1][3]; \
  c[0][4] = c[0][4] - coeff*c[1][4]; \
  r[0]   = r[0]   - coeff*r[1]; \
 \
  coeff = lhs[2][1]; \
  lhs[2][2]= lhs[2][2] - coeff*lhs[1][2]; \
  lhs[2][3]= lhs[2][3] - coeff*lhs[1][3]; \
  lhs[2][4]= lhs[2][4] - coeff*lhs[1][4]; \
  c[2][0] = c[2][0] - coeff*c[1][0]; \
  c[2][1] = c[2][1] - coeff*c[1][1]; \
  c[2][2] = c[2][2] - coeff*c[1][2]; \
  c[2][3] = c[2][3] - coeff*c[1][3]; \
  c[2][4] = c[2][4] - coeff*c[1][4]; \
  r[2]   = r[2]   - coeff*r[1]; \
 \
  coeff = lhs[3][1]; \
  lhs[3][2]= lhs[3][2] - coeff*lhs[1][2]; \
  lhs[3][3]= lhs[3][3] - coeff*lhs[1][3]; \
  lhs[3][4]= lhs[3][4] - coeff*lhs[1][4]; \
  c[3][0] = c[3][0] - coeff*c[1][0]; \
  c[3][1] = c[3][1] - coeff*c[1][1]; \
  c[3][2] = c[3][2] - coeff*c[1][2]; \
  c[3][3] = c[3][3] - coeff*c[1][3]; \
  c[3][4] = c[3][4] - coeff*c[1][4]; \
  r[3]   = r[3]   - coeff*r[1]; \
 \
  coeff = lhs[4][1]; \
  lhs[4][2]= lhs[4][2] - coeff*lhs[1][2]; \
  lhs[4][3]= lhs[4][3] - coeff*lhs[1][3]; \
  lhs[4][4]= lhs[4][4] - coeff*lhs[1][4]; \
  c[4][0] = c[4][0] - coeff*c[1][0]; \
  c[4][1] = c[4][1] - coeff*c[1][1]; \
  c[4][2] = c[4][2] - coeff*c[1][2]; \
  c[4][3] = c[4][3] - coeff*c[1][3]; \
  c[4][4] = c[4][4] - coeff*c[1][4]; \
  r[4]   = r[4]   - coeff*r[1]; \
 \
  pivot = 1.00/lhs[2][2]; \
  lhs[2][3] = lhs[2][3]*pivot; \
  lhs[2][4] = lhs[2][4]*pivot; \
  c[2][0] = c[2][0]*pivot; \
  c[2][1] = c[2][1]*pivot; \
  c[2][2] = c[2][2]*pivot; \
  c[2][3] = c[2][3]*pivot; \
  c[2][4] = c[2][4]*pivot; \
  r[2]   = r[2]  *pivot; \
 \
  coeff = lhs[0][2]; \
  lhs[0][3]= lhs[0][3] - coeff*lhs[2][3]; \
  lhs[0][4]= lhs[0][4] - coeff*lhs[2][4]; \
  c[0][0] = c[0][0] - coeff*c[2][0]; \
  c[0][1] = c[0][1] - coeff*c[2][1]; \
  c[0][2] = c[0][2] - coeff*c[2][2]; \
  c[0][3] = c[0][3] - coeff*c[2][3]; \
  c[0][4] = c[0][4] - coeff*c[2][4]; \
  r[0]   = r[0]   - coeff*r[2]; \
 \
  coeff = lhs[1][2]; \
  lhs[1][3]= lhs[1][3] - coeff*lhs[2][3]; \
  lhs[1][4]= lhs[1][4] - coeff*lhs[2][4]; \
  c[1][0] = c[1][0] - coeff*c[2][0]; \
  c[1][1] = c[1][1] - coeff*c[2][1]; \
  c[1][2] = c[1][2] - coeff*c[2][2]; \
  c[1][3] = c[1][3] - coeff*c[2][3]; \
  c[1][4] = c[1][4] - coeff*c[2][4]; \
  r[1]   = r[1]   - coeff*r[2]; \
 \
  coeff = lhs[3][2]; \
  lhs[3][3]= lhs[3][3] - coeff*lhs[2][3]; \
  lhs[3][4]= lhs[3][4] - coeff*lhs[2][4]; \
  c[3][0] = c[3][0] - coeff*c[2][0]; \
  c[3][1] = c[3][1] - coeff*c[2][1]; \
  c[3][2] = c[3][2] - coeff*c[2][2]; \
  c[3][3] = c[3][3] - coeff*c[2][3]; \
  c[3][4] = c[3][4] - coeff*c[2][4]; \
  r[3]   = r[3]   - coeff*r[2]; \
 \
  coeff = lhs[4][2]; \
  lhs[4][3]= lhs[4][3] - coeff*lhs[2][3]; \
  lhs[4][4]= lhs[4][4] - coeff*lhs[2][4]; \
  c[4][0] = c[4][0] - coeff*c[2][0]; \
  c[4][1] = c[4][1] - coeff*c[2][1]; \
  c[4][2] = c[4][2] - coeff*c[2][2]; \
  c[4][3] = c[4][3] - coeff*c[2][3]; \
  c[4][4] = c[4][4] - coeff*c[2][4]; \
  r[4]   = r[4]   - coeff*r[2]; \
 \
  pivot = 1.00/lhs[3][3]; \
  lhs[3][4] = lhs[3][4]*pivot; \
  c[3][0] = c[3][0]*pivot; \
  c[3][1] = c[3][1]*pivot; \
  c[3][2] = c[3][2]*pivot; \
  c[3][3] = c[3][3]*pivot; \
  c[3][4] = c[3][4]*pivot; \
  r[3]   = r[3]  *pivot; \
 \
  coeff = lhs[0][3]; \
  lhs[0][4]= lhs[0][4] - coeff*lhs[3][4]; \
  c[0][0] = c[0][0] - coeff*c[3][0]; \
  c[0][1] = c[0][1] - coeff*c[3][1]; \
  c[0][2] = c[0][2] - coeff*c[3][2]; \
  c[0][3] = c[0][3] - coeff*c[3][3]; \
  c[0][4] = c[0][4] - coeff*c[3][4]; \
  r[0]   = r[0]   - coeff*r[3]; \
 \
  coeff = lhs[1][3]; \
  lhs[1][4]= lhs[1][4] - coeff*lhs[3][4]; \
  c[1][0] = c[1][0] - coeff*c[3][0]; \
  c[1][1] = c[1][1] - coeff*c[3][1]; \
  c[1][2] = c[1][2] - coeff*c[3][2]; \
  c[1][3] = c[1][3] - coeff*c[3][3]; \
  c[1][4] = c[1][4] - coeff*c[3][4]; \
  r[1]   = r[1]   - coeff*r[3]; \
 \
  coeff = lhs[2][3]; \
  lhs[2][4]= lhs[2][4] - coeff*lhs[3][4]; \
  c[2][0] = c[2][0] - coeff*c[3][0]; \
  c[2][1] = c[2][1] - coeff*c[3][1]; \
  c[2][2] = c[2][2] - coeff*c[3][2]; \
  c[2][3] = c[2][3] - coeff*c[3][3]; \
  c[2][4] = c[2][4] - coeff*c[3][4]; \
  r[2]   = r[2]   - coeff*r[3]; \
 \
  coeff = lhs[4][3]; \
  lhs[4][4]= lhs[4][4] - coeff*lhs[3][4]; \
  c[4][0] = c[4][0] - coeff*c[3][0]; \
  c[4][1] = c[4][1] - coeff*c[3][1]; \
  c[4][2] = c[4][2] - coeff*c[3][2]; \
  c[4][3] = c[4][3] - coeff*c[3][3]; \
  c[4][4] = c[4][4] - coeff*c[3][4]; \
  r[4]   = r[4]   - coeff*r[3]; \
 \
  pivot = 1.00/lhs[4][4]; \
  c[4][0] = c[4][0]*pivot; \
  c[4][1] = c[4][1]*pivot; \
  c[4][2] = c[4][2]*pivot; \
  c[4][3] = c[4][3]*pivot; \
  c[4][4] = c[4][4]*pivot; \
  r[4]   = r[4]  *pivot; \
 \
  coeff = lhs[0][4]; \
  c[0][0] = c[0][0] - coeff*c[4][0]; \
  c[0][1] = c[0][1] - coeff*c[4][1]; \
  c[0][2] = c[0][2] - coeff*c[4][2]; \
  c[0][3] = c[0][3] - coeff*c[4][3]; \
  c[0][4] = c[0][4] - coeff*c[4][4]; \
  r[0]   = r[0]   - coeff*r[4]; \
 \
  coeff = lhs[1][4]; \
  c[1][0] = c[1][0] - coeff*c[4][0]; \
  c[1][1] = c[1][1] - coeff*c[4][1]; \
  c[1][2] = c[1][2] - coeff*c[4][2]; \
  c[1][3] = c[1][3] - coeff*c[4][3]; \
  c[1][4] = c[1][4] - coeff*c[4][4]; \
  r[1]   = r[1]   - coeff*r[4]; \
 \
  coeff = lhs[2][4]; \
  c[2][0] = c[2][0] - coeff*c[4][0]; \
  c[2][1] = c[2][1] - coeff*c[4][1]; \
  c[2][2] = c[2][2] - coeff*c[4][2]; \
  c[2][3] = c[2][3] - coeff*c[4][3]; \
  c[2][4] = c[2][4] - coeff*c[4][4]; \
  r[2]   = r[2]   - coeff*r[4]; \
 \
  coeff = lhs[3][4]; \
  c[3][0] = c[3][0] - coeff*c[4][0]; \
  c[3][1] = c[3][1] - coeff*c[4][1]; \
  c[3][2] = c[3][2] - coeff*c[4][2]; \
  c[3][3] = c[3][3] - coeff*c[4][3]; \
  c[3][4] = c[3][4] - coeff*c[4][4]; \
  r[3]   = r[3]   - coeff*r[4]; \
}


