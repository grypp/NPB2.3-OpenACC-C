/*--------------------------------------------------------------------
c  Parameter lm (declared and set in "npbparams.h") is the log-base2 of 
c  the edge size max for the partition on a given node, so must be changed 
c  either to save space (if running a small case) or made bigger for larger 
c  cases, for example, 512^3. Thus lm=7 means that the largest dimension 
c  of a partition that can be solved on a node is 2^7 = 128. lm is set 
c  automatically in npbparams.h
c  Parameters ndim1, ndim2, ndim3 are the local problem dimensions. 
c-------------------------------------------------------------------*/

#include "npbparams.h"

/* parameters */
/* actual dimension including ghost cells for communications */
#define	NM	(2+(2<<(LM-1)))
/* size of rhs array */
#define	NV	(2+(2<<(NDIM1-1))*(2+(2<<(NDIM2-1)))*(2+(2<<(NDIM3-1))))
/* size of residual array */
#define	NR	((8*(NV+(NM*NM)+5*NM+7*LM))/7)
/* size of communication buffer */
#define	NM2	(2*NM*NM)
/* maximum number of levels */
#define	MAXLEVEL	11

/*---------------------------------------------------------------------*/
/* common /mg3/ */
static int nx[MAXLEVEL+1], ny[MAXLEVEL+1], nz[MAXLEVEL+1];
/* common /ClassType/ */
static char Class;
/* common /my_debug/ */
static int debug_vec[8];
/* common /fap/ */
/*static int ir[MAXLEVEL], m1[MAXLEVEL], m2[MAXLEVEL], m3[MAXLEVEL];*/
static int m1[MAXLEVEL+1], m2[MAXLEVEL+1], m3[MAXLEVEL+1];
static int lt, lb;
static int m1lt, m2lt, m3lt;

/*c---------------------------------------------------------------------
c  Set at m=1024, can handle cases up to 1024^3 case
c---------------------------------------------------------------------*/
#define	M	1037

/* common /buffer/ */
/*static double buff[4][NM2];*/

#define GET_3D(mat,i) (mat+((i)*m3lt*m2lt*m1lt))
#define ACCESS_3D(mat,i,j,k) (*(mat+((i)*m2lt*m1lt)+((j)*m1lt)+(k)))

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

#define VRANLC(n,x_seed,a) \
{ \
    int i; \
    double x,t1,t2,t3,t4,a1,a2,x1,x2,loc_z; \
 \
    t1 = r23 * a; \
    a1 = (int)t1; \
    a2 = a - t23 * a1; \
    x = x_seed; \
 \
    for (i = 1; i <= n; i++) { \
        t1 = r23 * x; \
        x1 = (int)t1; \
        x2 = x - t23 * x1; \
        t1 = a1 * x2 + a2 * x1; \
        t2 = (int)(r23 * t1); \
        loc_z = t1 - t23 * t2; \
        t3 = t23 * loc_z + a2 * x2; \
        t4 = (int)(r46 * t3); \
        x = t3 - t46 * t4; \
        ACCESS_3D(z,i3,i2,i) = r46 * x; \
    } \
    x_seed = x; \
}
