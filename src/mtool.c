/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "stdmath.h"
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include "globdef.h"
#include "util.h"
#include "mtool.h"
#include "rand.h"
#include "plot.h"

#define  DDC_PI  (3.14159265358979323846)
#define  ONEBYPI ((double) 0.3183098861837906715377675)

#define TRUE  1
#define FALSE 0

#define BITS_PER_WORD   (sizeof(unsigned) * 8)

#define MT_DST_SMP_PER_FRQ	100

#define RND_IA 		16807
#define RND_IM 		2147483647
#define RND_AM  	(1.0/RND_IM)
#define RND_IQ  	127773
#define RND_IR          2836
#define RND_STAB_SZ 	MT_RND_STAB_SZ
#define RND_NDIV       	(1+(RND_IM-1)/RND_STAB_SZ)
#define RND_RNMAX       (1.0-(1.2e-10))
#define RND_MODMUL64( acc, k ) \
{ \
  k = (acc) / RND_IQ; \
  acc = RND_IA*( (acc) - (k)*RND_IQ ) - (k)*RND_IR; \
  if( (acc) < 0 ) acc += RND_IM; \
}  

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define MT_SORT_DST(dst,do_force)\
{\
  if(!(dst)->sorted || do_force)\
  {\
    (dst)->vDist = MtVecHeapSort( (dst)->vDist, (dst)->nSmp );\
    (dst)->sorted = TRUE;\
  }  \
}

/********************************************************************/
/*								    */
/*	typedefs    					            */
/*								    */
/********************************************************************/

typedef struct
{
  int    idx;
  MtReal distance; 
}
MtNLRegNeighbor;

struct MtNLRegAdm
{
  int               n,          /* number of samples */
                    m;          /* dimension of input vectors */
  MtNLRegNeighbor  *neighbors;  /* array of n neighbor points */
  MtVector         *samplesIn;  /* array of m input-sequences */
  MtVector          samplesOut; /* sequence of n output values */
};



/********************************************************************/
/*								    */
/*	global data					            */
/*								    */
/********************************************************************/

static long  rndAcc;
static long  rndStab[ RND_STAB_SZ ];
static long  rndIy = 0; 

/********************************************************************/
/*								    */
/*	forward prototypes                                           */
/*								    */
/********************************************************************/

/* static int tprob(double t, int dof, double *probability);*/
void MtMatLubksb_f(MtMatrix a,int n,int *indx, float *b);

typedef struct
{
  int gw;
  PlScreen ps;
  PlPlot pl;
  PlWindow pw;
  PlLineHd plh;
}
PlGWin;

#if 0
#ifdef MI_WITH_64BIT

#define NR3_ERROR 1

int NR3_SVD_Decomp( 
    		double *A, int n, int m, 
    	      	double *U, double *V,double *S)
{
  TRCERR(("not yet implement for 64 bit\n"));
  return 1;
}
#endif
#endif

#if LMS_TEST
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
PlGWin  *PlOpenGWin( int sx, int sy, char *titel )
{
  int rcThis=0;
  PlGWin *pw=NULL;
  float fx0,fx1,fy0,fy1;

  ALLOCOBJ(pw);

  pw->gw = -1;

  pw->gw = GWopenx(0,sx,sy,19,0,10,titel);
  if(pw->gw==0) TRCERRR(("GWopenx() failed\n"),MI_ERROR);

  GWgetwn( &fx0, &fy0, &fx1, &fy1 );
  PlScreenInit( &pw->ps, 0, 0, fx1-fx0, fy1-fy0 );
  PlWindowInit( &pw->pw, &pw->ps, 0, 0, 100, 100 );
rcCatch:
  return pw;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void PlCloseGWin( PlGWin *pg )
{
  if(pg->gw>0)
  {
    GWclose(pg->gw);
  }
  FREEMEM(pg);
}
#endif

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
double MtDistPoint2LineSegSquared(
    		double px, double py,	// point  
    		double x1, double y1, double x2, double y2 // line
		)
{
  double vx = x1-px, vy = y1-py, ux = x2-x1, uy = y2-y1;
  double length = ux*ux+uy*uy;
  double det = (-vx*ux)+(-vy*uy); //if this is < 0 or > length then its outside the line segment
  if(det<0||length<MT_NUM_EPS)
    return (x1 - px) * (x1 - px) + (y1 - py) * (y1 - py);
  if(det>length)
    return (x2 - px) * (x2 - px) + (y2 - py) * (y2 - py);
  det = ux*vy-uy*vx;
  return (det*det)/length;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
double MtSmoothstep(double a, double b, double x)
{
    if (x < a)
        return 0;
    if (x >= b)
        return 1;
    x = (x - a)/(b - a); /* normalize to [0:1] */
    return (x*x * (3 - 2*x));
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
double MtLinstep(double a, double b, double x)
{
    if (x < a)
        return 0;
    if (x >= b)
        return 1;
    x = (x - a)/(b - a); /* normalize to [0:1] */
    return x;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MtStereographProj( double *v_pq,	// IN/OUT: p,q 
    			  double *v_fg, //  OUT/OUT: f,g 
			  int dir // direction: 1 = p,q -> f,g 
			  	  //           -1 = inverse f,g -> p,q
			  )
{
  double r;
 
  if(dir>0)
  {
    r= sqrt(1+v_pq[0]*v_pq[0]+v_pq[1]*v_pq[1]);
    v_fg[0] = v_pq[0]/(1.+r);
    v_fg[1] = v_pq[1]/(1.+r);
  }
  else
  {
    r = v_fg[0]*v_fg[0]+v_fg[1]*v_fg[1]-1;
    if(fabs(r)<MT_NUM_EPS)
    {
     v_pq[0] = 0;
     v_pq[1] = 0;
    }
    else
    {
     v_pq[0] = -2*v_fg[0]/r;
     v_pq[1] = -2*v_fg[1]/r;
    }
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtDerivPoly2D( 
    	double *P, // IN: polynomial coefficients
	int n, 	// polynomial order
	double x, double y,
	double *PX, double *PY // Derivatives (polynomials of order n-1)
	)
// See also description for MtEvalPoly2D
{
#define MT_EPD2D_TST 0
#define MT_POLY2D_MAXORD 128
  int rcThis=0;
  int n2,i,j,k,k2,kmax,kmax2,k3;
  double tx[MT_POLY2D_MAXORD+1],ty[MT_POLY2D_MAXORD];
  double ax[MT_POLY2D_MAXORD+1],ay[MT_POLY2D_MAXORD],
	 bx,by;

  UT_ASSERT_FATAL(n<=MT_POLY2D_MAXORD);  

  n2=n-1;
  kmax = (n+1)*(n+2)/2-1;
  kmax2 = (n2+1)*(n2+2)/2-1;

  tx[0]=ty[0]=0;
  ax[0]=ay[0]=1;
  k=0,k2=k3=0;
  for(i=0;i<=n;i++)
  {
    if(i>0)
    {
      tx[i]=tx[i-1];
      ty[i]=ty[i-1]+1;
      ax[i]=ax[i-1];
      ay[i]=ay[i-1]*y;      
    }
    for(j=0;j<=i;j++)
    {
      if(j<i) 
      {
	tx[j]++;	
	ax[j] *= x;
      }
      UT_ASSERT_FATAL(k<=kmax);
      if(PX && ((int)(tx[j])>0))
      {
        bx = tx[j]*P[k];
        UT_ASSERT_FATAL(k2<=kmax2);
#if MT_EPD2D_TST  
      printf("%d %g ",(int)tx[j],P[j]);
#endif
        PX[k2] = bx;
	k2++;
      }
      if(PY && ((int)(ty[j])>0))
      {
        by = ty[j]*P[k];
        UT_ASSERT_FATAL(k3<=kmax2);
#if MT_EPD2D_TST  
      printf("%d %g ",(int)ty[j],P[j]);
#endif
        PY[k3] = by;
	k3++;
      }

      k++;
    }
#if MT_EPD2D_TST  
    printf("\n");
#endif
  }
  if(PX)
  {
    UT_ASSERT_FATAL(k2==kmax2+1);
  }
  //rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
double MtEvalPoly2D( double *P, // polynomial coefficients
    		     int no, 	// order
    		     double x, double y,
		     double *B  // OUT: (optional) basis functions at (x,y)
		     		//	(P not used)
		     )
// Z = POLYVAL2(P,X,Y) returns the value of a 2D polynomial P evaluated at (X,Y). P
// is a vector of length (N+1)*(N+2)/2 containing the polynomial coefficients in
// ascending powers:
//
//   P = [p00 p10 p01 p20 p11 p02 p30 p21 p12 p03...]
//
// e.g. For a 3rd order fit, polyval2.m evaluates the matrix equation:
//
//    Z = V*P    or
//
//                   2      2  3  2      2   3
//    Z = [1  x  y  x  xy  y  x  x y  x y   y ]  [p00
//                                                p10
//                                                p01
//                                                p20
//                                                p11
//                                                p02
//                                                p30
//                                                p21
//                                                p12
//                                                p03]
{
#define MT_POLY2D_MAXORD 128
#define MT_EP2D_TST 0
  int i,j,k,kmax;
#if MT_EP2D_TST
  double tx[MT_POLY2D_MAXORD+1],ty[MT_POLY2D_MAXORD];
#endif
  double ax[MT_POLY2D_MAXORD+1],ay[MT_POLY2D_MAXORD],
	 z,b;

  UT_ASSERT_FATAL(no<=MT_POLY2D_MAXORD);  

  kmax = (no+1)*(no+2)/2-1;
#if MT_EP2D_TST  
  tx[0]=ty[0]=0;
#endif
  ax[0]=ay[0]=1;
  z=0;k=0;
  for(i=0;i<=no;i++)
  {
    if(i>0)
    {
#if MT_EP2D_TST  
      tx[i]=tx[i-1];
      ty[i]=ty[i-1]+1;
#endif
      ax[i]=ax[i-1];
      ay[i]=ay[i-1]*y;      
    }
    for(j=0;j<=i;j++)
    {
      if(j<i) 
      {
#if MT_EP2D_TST  
	tx[j]++;	
#endif
	ax[j] *= x;
      }
#if MT_EP2D_TST  
      printf("%d%d (%g,%g)",(int)tx[j],(int)ty[j],ax[j],ay[j]);
#endif
      b = ax[j]*ay[j];
      UT_ASSERT_FATAL(k<=kmax);
      if(P) z += b*P[k];
      if(B) B[k] = b;
      k++;
    }
#if MT_EP2D_TST  
    printf("\n");
#endif
  }

  return z;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MtGrad2NormVec( double *g, double *n )
{
  n[0] = -g[0];
  n[1] = -g[1];
  n[2] = 1;
  MtVecNorm_d(n,3,-1);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MtNormVec2Grad( double *n, double *g )
{
  g[0] = -n[0]/n[2];
  g[1] = -n[1]/n[2];
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtAngular2NormVec( 
    		double tau, // in: tilt (azimuth, angle with x-axis) 
    		double sig, // in: slant (angle with z-axis axis)
    		double *vn  // out: 3D unit normal vector  
		)
{
  double ss=sin(sig);
  vn[0] = ss*cos(tau);
  vn[1] = ss*sin(tau);
  vn[2] = cos(sig);  
  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtNormVec2Angular( double *vn, double *tau, double *sig )
{
  int rcThis=0;
  double tilt,slant,pi=MT_PI,dx=vn[0],dy=vn[1],dz=vn[2];

  //    tilt = atan2(yo,vn[0]+.00000001);

  if(fabs(dx)<MT_NUM_EPS && fabs(dy)<MT_NUM_EPS) 
    dx=MT_NUM_EPS;

  tilt = atan2(dy,dx );

  tilt = tilt * 180.0/pi;
  if (tilt < 0) tilt = tilt + 360;
  tilt = tilt * pi/180.0;

  UT_ASSERT(!(fabs(dz)>1));

  slant = acos(dz);

  *tau = tilt;
  *sig = slant;

rcCatch:
  return rcThis;
}

/********************************************************************/
/*								    */
/*								    */
/********************************************************************/
/********************************************************************/			 
double MtSignum( double x )
/********************************************************************
*********************************************************************/
{
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}

/********************************************************************/			 
MtReal MtSign( MtReal x )
/********************************************************************
*********************************************************************/
{ // DEPRECATED (incorrect for x=0)
  return x >= 0  ? 1 : -1;  
}

/********************************************************************/			 
double MtAbsPow( double x, double a )
/********************************************************************
*********************************************************************/
{
  return MtSign(x)*pow(fabs(x),a);
}

/********************************************************************/			 
float MtAbsPow_f( double x, double a )
/********************************************************************
*********************************************************************/
{
  return MtSign(x)*powf(fabs(x),a);
}

/********************************************************************/
/*								    */
/*								    */
/********************************************************************/

/********************************************************************/			 
MtReal MtVecTrendStrength( MtVector x, int n )
/********************************************************************
  Mann-Kendall Trend-Strength indicator (Tau)
  n   #samples
*********************************************************************/
{
  int     j,k,m;
  MtReal  sum;
  
  sum = 0;
  m = 0;
  for(j=0; j<n-1; j++)
  {
    for(k=j+1;k<=n-1;k++)
    {
      sum += x[k] > x[j] ? 1 : -1;
      m++;
    }
  }
  return sum / (MtReal)m; 
}

/********************************************************************/			 
MtVector MtVecMovTrendStrength( MtVector y, MtVector x, int n, int d )
/********************************************************************
  n   #samples
  d   delay 
  y = x NOT allowed  
*********************************************************************/
{
  int     i,j,k,m,d2;
  MtReal  sum;
  
  for(i=0;i<n;i++)
  {
    d2 = i < d - 1 ? i + 1 : d;
    sum = 0;
    m = 0;
    for(j=i-d2+1; j<i; j++)
    {
      for(k=j+1;k<=i;k++)
      {
        sum += x[k] > x[j] ? 1 : -1;
        m++;
      }
    }
    y[i] = sum / (MtReal)m; 
  }
  return y;
}

/********************************************************************/			 
MtVector MtVecDetrendVola( MtVector y, MtVector x, int n, int d )
/********************************************************************
  n   #samples
  d   delay
  y=x allowed    
*********************************************************************/
{
  int       i,j,k;
  MtVector  v = MtVecOpen(n);
  
  for(i=0;i<n;i++) v[i] = fabs(x[i]);
  MtVecMA_f( v, v, n, d );
  for(i=0;i<n;i++) 
  {
    /**** if v[i] is zero then search next nonzero value */
    for(j=0,k=i; (j<n)&&(v[k]<MT_REAL_EPSILON); j++,k=(k+1)%n );
    y[i] = x[i]/v[k];    
  }  
  MtVecClose( v ); 
  return y;
}

/********************************************************************/			 
double * MtVecMA( double * y, double * x, int n, int d )
/********************************************************************
  y = MovingAverage( x, d )
  n   #samples
  d   delay
  y=x allowed    
*********************************************************************/
{
  int 	     i;
  MtRQueue  *rq;
  double     sum, x0;
  
  d = MIN( n, d );
  rq = MtRQueueOpen( d );
  if( ! rq ) goto fail;
  x0 = x[0];
  for( i = 0; i < d; i++ ) MT_RQUEUE_PUSH( rq, x0 );
  sum = d * x0;
  for( i = 0; i < n; i++ )
  {
    MT_RQUEUE_POP( rq, x0 );
    sum -= x0;
    x0 = x[i];
    MT_RQUEUE_PUSH( rq, x0 );
    sum += x[i];
    y[i] = sum / (double)d;  
  }
  
  MtRQueueClose( rq );
  return y;
  
fail:
  return NULL;  
}

/********************************************************************/			 
MtVector MtVecMA_f( MtVector y, MtVector x, int n, int d )
/********************************************************************
  y = MovingAverage( x, d )
  n   #samples
  d   delay
  y=x allowed    
*********************************************************************/
{
  int 	     i;
  MtRQueue  *rq;
  MtReal     sum, x0;
  
  d = MIN( n, d );
  rq = MtRQueueOpen( d );
  if( ! rq ) goto fail;
  x0 = x[0];
  for( i = 0; i < d; i++ ) MT_RQUEUE_PUSH( rq, x0 );
  sum = d * x0;
  for( i = 0; i < n; i++ )
  {
    MT_RQUEUE_POP( rq, x0 );
    sum -= x0;
    x0 = x[i];
    MT_RQUEUE_PUSH( rq, x0 );
    sum += x[i];
    y[i] = sum / (MtReal)d;  
  }
  
  MtRQueueClose( rq );
  return y;
  
fail:
  return NULL;  
}

/********************************************************************/			 
MtVector MtVecMovMedian( MtVector y, MtVector x, int n, int d )
/********************************************************************
  y = MovingAverage( x, d )
  n   #samples
  d   delay
  y=x allowed    
*********************************************************************/
{
  int 	     i,h;
  
  for( i=0; i<n; i++ )
  {
    if( i < d )
    {
      h = i + 1;
    }
    else
    {
      h = d;
    }
    /* TODO: performance */
    y[i] = MtVecMedian( x + i - h + 1, h, 0 );
  }          
  return y;  
}

/********************************************************************/			 
double * MtVecMovVar( double * y, double * x, int n, int d )
/********************************************************************
  y = MovingVariance( x, d )
  n   #samples
  d   delay
  y=x allowed    
*********************************************************************/
{
  int 	     i;
  double     avg, var;
  MtMAV      mav;            

  mav.max=mav.sum=mav.sumQ=0;
  MtMAVInit( &mav, d, x[0] );
  for(i=0;i<n;i++)
  {
    MtMAVUpdate( &mav, x[i] );
    MtMAVGet( &mav, &avg, &var );
    y[i] = sqrt( fabs(var) );
  }
  return y;  
}

#if 0
/********************************************************************/			 
MtVector MtVecMA( MtVector y, MtVector x, int n, int d )
/********************************************************************
  y = MovingAverage( x, d )
  n   #samples
  d   delay
  y=x NOT allowed    
*********************************************************************/
{
  int 	     i;
  MtReal     sum, x0, xOld;
  
  d = MIN( n, d );
  x0 = x[0];
  sum = d * x0;
  for( i=0; i<n; i++ )
  {
    if( i < d )
      xOld = x0;
    else
      xOld = x[i-d];    
    sum += x[i] - xOld;
    y[i] = sum / (MtReal)d;  
  }  
  return y;  
}
#endif

/********************************************************************/			 
MtVector MtVecMovDiff( MtVector y, MtVector x, int n, int d, 
                       int signum )
/********************************************************************
  y = Moving Differences( x, d )
  n   #samples
  d   delay
  y=x NOT allowed    
*********************************************************************/
{
  int 	     i;
  MtReal     sum;
  
  d = MIN( n, d );
  sum = 0;
  y[0] = 0;
  for( i = 1; i < n; i++ ) 
  {                
    sum += MAX( signum * ( x[i] - x[i-1] ), 0 );
    if( i > d )
      sum -= MAX( signum * ( x[i-d] - x[i-d-1] ), 0 );
    y[i] = sum; 
  }   
  return y;  
}

/********************************************************************/			 
double * MtVecMAE( double * y, double * x, int n, double beta )
/********************************************************************
  y  =  Exponential Moving Average ( x, beta )
  n     #samples
  beta  slowing 
  y=x allowed ?    
*********************************************************************/
{
  int	  i;
  double  mae;

  if(beta<0)
  {
    beta = 1.-2./((double)(-beta)+1);
  }
  
  mae = y[0] = x[0];
  for( i = 1; i < n; i++ )
  {
    mae = beta * mae + ( 1 - beta ) * x[i];
    y[i] = mae;
    /* printf("x[%d]=%f, y[%d]=%f\n",i,x[i],i,y[i]); */
  }
  return y;
}

/********************************************************************/			 
double *MtVecMomentum( double *y, double *x, int n, int d )
/********************************************************************
  y = Momentum( x, d )
  n   #samples
  d   delay
  y=x NOT allowed    
*********************************************************************/
{
  int  i,j;
  
  for( i = 0; i < n; i++ )
  {
    if( i >= d )
      j = i - d;  
    else 
      j = 0;         
    y[i] =  log(x[i]/x[j]);
/*    y[i] =  x[i]/x[j] - 1 ;*/
  }  
  return y;
}

/********************************************************************/			 
MtVector MtVecMomentumDiff( MtVector y, MtVector x, int n, int d )
/********************************************************************
  y = Momentum( x, d )
  n   #samples
  d   delay
  y=x NOT allowed    
*********************************************************************/
{
  int  i,j;
  
  for( i = 0; i < n; i++ )
  {
    if( i >= d )
      j = i - d;  
    else 
      j = 0;         
    y[i] =  x[i] - x[j];
/*    y[i] =  x[i]/x[j] - 1 ;*/
  }  
  return y;
}

/********************************************************************/			 
double *MtVecMomentumLog( double * y, double * x, int n, int d )
/********************************************************************
  y = Momentum( x, d )
  n   #samples
  d   delay
  y=x NOT allowed    
*********************************************************************/
{
  int  i,j;
  
  for( i = 0; i < n; i++ )
  {
    if( i >= d )
      j = i - d;  
    else 
      j = 0;         
    y[i] =  log(x[i]/x[j]);
  }  
  return y;
}

/********************************************************************/			 
MtVector MtVecAdd( MtVector y, MtVector x, int n )
/********************************************************************
  y = y  + x
  n   #samples
  y=x allowed    
*********************************************************************/
{
  int  i;
  
  for( i = 0; i < n; i++ ) y[i] = y[i] + x[i];
  return y;
}

/********************************************************************/			 
MtVector MtVecSub( MtVector y, MtVector x, int n )
/********************************************************************
  y = y - x
  n   #samples
  y=x allowed    
*********************************************************************/
{
  int  i;
  
  for( i = 0; i < n; i++ ) y[i] = y[i] - x[i];
  return y;
}

/********************************************************************/			 
MtVector MtVecConst( MtVector y, int n, MtReal c )
/********************************************************************
  y = c  
  n   #samples
*********************************************************************/
{
  int  i;
  
  for( i = 0; i < n; i++ ) y[i] = c;
  return y;
}

/********************************************************************/			 
double *MtVecScale( double *y, int n, double alpha )
/********************************************************************
  y = alpha * y 
  n   #samples
*********************************************************************/
{
  int  i;
  
  for( i = 0; i < n; i++ ) y[i] = y[i] * alpha;
  return y;
}

/********************************************************************/			 
MtVector MtVecScale_f( MtVector y, int n, MtReal alpha )
/********************************************************************
  y = alpha * y 
  n   #samples
*********************************************************************/
{
  int  i;
  
  for( i = 0; i < n; i++ ) y[i] = y[i] * alpha;
  return y;
}

/********************************************************************/			 
void MtMapNormInterval( MtReal xMin, MtReal xRange, MtReal *t )
/********************************************************************
 * Map interval x=[xMin, xMin+xRange] to y=[0,1]
 * via transformation t:   y = t[0]*x + t[1]
*********************************************************************/
{
  t[0] =  xRange > MT_REAL_EPSILON ? 1./xRange : 0;
  t[1] =  - t[0] * xMin;
}


/********************************************************************/			 
int MtVecSample( MtVector y, MtVector x, int n, int d )
/********************************************************************
  y[i] = x[(i+1)*d-1] 
  n   #samples
  d   spacing
  y = x allowed
  return  actual length of output sample ( = n/d )
*********************************************************************/
{
  int  i,m;
  m = n/d;  
  for( i = 0; i < m; i++ ) y[i] = x[(i+1)*d-1];
  return m;
}

/********************************************************************/			 
MtReal MtVecDistance( MtVector x, MtVector y, int n )
/********************************************************************
  compute Euklidian distance between x and y
  n   #samples
*********************************************************************/
{
  int 	 i;
  MtReal tmp, diff = 0;
  
  for(i=0;i<n;i++)
  {
    tmp = x[i] - y[i];
    diff += tmp * tmp;
  }                                             
  return sqrt( diff );                                                
}

/********************************************************************/			 
double MtVecDistance_d( double *x, double *y, int n )
/********************************************************************
  compute Euklidian distance between x and y
  n   #samples
*********************************************************************/
{
  int 	 i;
  double tmp, diff = 0;
  
  for(i=0;i<n;i++)
  {
    tmp = x[i] - y[i];
    diff += tmp * tmp;
  }                                             
  return sqrt( diff );                                                
}

/********************************************************************/			 
MtReal MtVecDistance2_f( MtVector x, MtVector y, int n )
/********************************************************************
  compute squared Euklidian distance between x and y
  n   #samples
*********************************************************************/
{
  int 	 i;
  MtReal tmp, diff = 0;
  
  for(i=0;i<n;i++)
  {
    tmp = x[i] - y[i];
    diff += tmp * tmp;
  }                                             
  return diff ;                                                
}

/********************************************************************/			 
double MtVecDistance2( double *x, double *y, int n )
/********************************************************************
  compute squared Euklidian distance between x and y
  n   #samples
*********************************************************************/
{
  int 	 i;
  double tmp, diff = 0;
  
  for(i=0;i<n;i++)
  {
    tmp = x[i] - y[i];
    diff += tmp * tmp;
  }                                             
  return diff ;                                                
}

/********************************************************************/			 
MtVector MtVecNorm( MtVector y, int n, MtReal len )
/********************************************************************
 y = normalize(y)
 len (optional) = user supplied vector length 
     if len < 0 => compute length;
*********************************************************************/
{
  int 	 i;
  
  if( len < 0 ) len = MtVecAbs( y, n );
  if( len > MT_REAL_EPSILON )
    for(i=0;i<n;i++) y[i] /= len;
  else MtVecConst( y, n, 0 );
  return y;           
}

/********************************************************************/			 
double *MtVecNorm_d( double *y, int n, double len )
/********************************************************************
 y = normalize(y)
 len (optional) = user supplied vector length 
     if len < 0 => compute length;
*********************************************************************/
{
  int 	 i;

  if( len < 0 ) len = MtVecAbs_d( y, n );
  if( len > MT_REAL_EPSILON )
    for(i=0;i<n;i++) y[i] /= len;
  else for(i=0;i<n;i++) y[i] = 0;
  return y;           
}


/********************************************************************/			 
MtReal MtVecAbs( MtVector x, int n )
/********************************************************************
Absolute value (magnitude) of vector x
*********************************************************************/
{
  int 	 i;
  double diff = 0;
  
  for(i=0;i<n;i++)
  {
    diff += x[i] * x[i];
  }                                             
  return sqrt( diff );                                                
}

/********************************************************************/			 
double MtVecAbs_d( double *x, int n )
/********************************************************************
Absolute value (magnitude) of vector x
*********************************************************************/
{
  int 	 i;
  double diff = 0;
  
  for(i=0;i<n;i++)
  {
    diff += x[i] * x[i];
  }                                             
  return sqrt( diff );                                                
}

/********************************************************************/			 
MtReal MtVecAbs2( MtVector x, int n )
/********************************************************************
Absolute value (magnitude) of vector x
*********************************************************************/
{
  int 	 i;
  MtReal diff = 0;
  
  for(i=0;i<n;i++)
  {
    diff += x[i] * x[i];
  }                                             
  return diff;                                                
}

/********************************************************************/			 
double MtVecAbs2_d( double *x, int n )
/********************************************************************
Absolute value (magnitude) of vector x
*********************************************************************/
{
  int 	 i;
  double diff = 0;
  
  for(i=0;i<n;i++)
  {
    diff += x[i] * x[i];
  }                                             
  return diff;                                                
}

/********************************************************************/			 
MtReal MtVecDotprod( MtVector x, MtVector y, int n )
/********************************************************************
return vector dot product x'y
*********************************************************************/
{
  int 	 i;
  MtReal tmp = 0;
  
  for(i=0;i<n;i++)
    tmp += x[i] * y[i];
  return tmp;                                                
}

/********************************************************************/			 
double MtVecDotprod_d( double *x, double *y, int n )
/********************************************************************
return vector dot product x'y
*********************************************************************/
{
  int 	 i;
  double tmp = 0;
  
  for(i=0;i<n;i++)
    tmp += x[i] * y[i];
  return tmp;                                                
}
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
double *MtVec3CrossProd( double *a, double *b, double *c )
  // a = b x c
{
    a[0] = b[1] * c[2] - b[2] * c[1];
    a[1] = b[2] * c[0] - b[0] * c[2]; 
    a[2] = b[0] * c[1] - b[1] * c[0]; 
  return a;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
float *MtVec3CrossProd_f( float *a, float *b, float *c )
  // a = b x c
{
    a[0] = b[1] * c[2] - b[2] * c[1];
    a[1] = b[2] * c[0] - b[0] * c[2]; 
    a[2] = b[0] * c[1] - b[1] * c[0]; 
  return a;
}

/********************************************************************/			 
MtVector MtVecShuffle_f( MtVector y, int n, int rSeed )
/********************************************************************
  y = random permutation of y 
  n   #samples
  rSeed : random mumber generator starting value ( -1 = no init )
*********************************************************************/
{  
  int 		i,j;
  MtReal        tmp;
    
  if( rSeed >= 0 ) MtRndSeed( rSeed );
  if( n < 2 ) return y;
  for(i=0;i<n;i++)
  {
    j = MtRndGet() * ( n - 1 ); UT_ASSERT0(FALSE);
    if( j <= 0 ) j = 0;
    if( j > n - 1 ) j = n - 1;
    /* printf("%d %d\n",i,j); */
    tmp = y[i];
    y[i] = y[j];
    y[j] = tmp;  
  }  
  return y;
}

/********************************************************************/			 
double *MtVecShuffle( double *y, int n, int rSeed )
/********************************************************************
  y = random permutation of y 
  n   #samples
  rSeed : random mumber generator starting value ( -1 = no init )
*********************************************************************/
{  
  int 		i,j;
  double        tmp,u;
    
  if( rSeed >= 0 ) MtRndSeed( rSeed );
  if( n < 2 ) return y;
  for(i=0;i<n;i++)
  {
    u = MtRndGet();
    j = u*n;
    if(j<0||j>n-1)
    {
      TRCERR((
	    "rand() index out of range u=%g,k=%d,max=%d\n",u,j,n));
      if( j <= 0 ) j = 0;
      if( j > n - 1 ) j = n - 1;
    }
    /* printf("%d %d\n",i,j); */
    tmp = y[i];
    y[i] = y[j];
    y[j] = tmp;  
  }  
  return y;
}

/********************************************************************/			 
void MtIntShuffle( int *v, int n )
/********************************************************************/
{
  int 		i,j,tmp;
    
  if( n < 2 ) return;
  for(i=0;i<n;i++)
  {
    j = n * (float)rand()/(RAND_MAX+1.0);
    if( j <= 0 ) j = 0;
    if( j > n - 1 ) j = n - 1;
    tmp = v[i];
    v[i] = v[j];
    v[j] = tmp;  
  }  
}


/**********************************************************************/
static int MtRealCmpFunc( const void *x, const void *y )
/**********************************************************************
***********************************************************************/
{
if( *((MtReal *)x) >= *((MtReal *)y))
  return  1;
else
  return -1;
}

/**********************************************************************/
static int MtDoubleCmpFunc( const void *x, const void *y )
/**********************************************************************
***********************************************************************/
{
if( *((double *)x) >= *((double *)y))
  return  1;
else
  return -1;
}

/********************************************************************/			 
MtVector MtVecSort( MtVector y, int n )
/********************************************************************
sort vector (ascending order)
  y  vector to be sorted
  n  #samples
*********************************************************************/
{  
  qsort( y, n, sizeof( MtReal ), MtRealCmpFunc );
  return y;
}

/********************************************************************/			 
double *MtVecSort_d( double *y, int n )
/********************************************************************
sort vector (ascending order)
  y  vector to be sorted
  n  #samples
*********************************************************************/
{  
  qsort( y, n, sizeof( double ), MtDoubleCmpFunc );
  return y;
}

/********************************************************************/			 
double MtVecMedian_d( double *x, int n )
/********************************************************************
return median value
*********************************************************************/
{  
  double     median,*y=NULL;
  ALLOCARR(y,n);
  MtVecSort_d( MtVecCopy_d( y,x,n ),n );
  if( n % 2 == 0 )
    median = ( y[n/2] + y[n/2-1] ) / 2;
  else
    median = y[n/2];  
  FREEMEM(y);
  return median;
}

/********************************************************************/			 
MtReal MtVecMedian( MtVector x, int n, int is_sorted )
/********************************************************************
return median value
*********************************************************************/
{  
  MtReal     median;
  MtVector   y=NULL;
  
  if(is_sorted)
  {
    y=x;
  }
  else
  {
    y=MtVecOpen(n);
    MtVecSort( MtVecCopy_f( y,x,n ),n );
  }

  if( n % 2 == 0 )
    median = ( y[n/2] + y[n/2-1] ) / 2;
  else
    median = y[n/2];  
  
  if(!is_sorted) MtVecClose( y );
  return median;
}

#define FSWAP(a,b) temp=(a);(a)=(b);(b)=temp;
/********************************************************************/			 
MtReal MtVecSelect(float arr[], unsigned long n, unsigned long k) 
/********************************************************************
  select k-th smallest element & split/rearrange array accordingly 	   
*********************************************************************/
{
  unsigned long i,ir,j,l,mid;
  float a,temp;
  
  l=1;
  ir=n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
	FSWAP(arr[l],arr[ir])
	}
      return arr[k];
    } else {
      mid=(l+ir) >> 1;
      FSWAP(arr[mid],arr[l+1])
	if (arr[l+1] > arr[ir]) {
	  FSWAP(arr[l+1],arr[ir])
	  }
      if (arr[l] > arr[ir]) {
	FSWAP(arr[l],arr[ir])
	}
      if (arr[l+1] > arr[l]) {
	FSWAP(arr[l+1],arr[l])
	}
      i=l+1;
      j=ir;
      a=arr[l];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	FSWAP(arr[i],arr[j])
	}
      arr[l]=arr[j];
      arr[j]=a;
      if (j >= k) ir=j-1;
      if (j <= k) l=i;
    }
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MtConvToRank(double *w, // IN: sorted array
    // OUT: rank with tie resolution 
    int n
    )
{
  int j=1,ji,jt;
  double t,rank;
  double s=0.0;
  while (j < n) {
    if (w[j] != w[j-1]) {
      w[j-1]=j;
      ++j;
    } else {
      for (jt=j+1;jt<=n && w[jt-1]==w[j-1];jt++);
      rank=0.5*(j+jt-1);
      for (ji=j;ji<=(jt-1);ji++)
	w[ji-1]=rank;
      t=jt-j;
      s += (t*t*t-t);
      j=jt;
    }
  }
  if (j == n) w[n-1]=n;
}

void MtConvToRank_f(float *w, // IN: sorted array
    // OUT: rank with tie resolution 
    int n
    )
{
  int j=1,ji,jt;
  float t,rank;
  float s=0.0;
  while (j < n) {
    if (w[j] != w[j-1]) {
      w[j-1]=j;
      ++j;
    } else {
      for (jt=j+1;jt<=n && w[jt-1]==w[j-1];jt++);
      rank=0.5*(j+jt-1);
      for (ji=j;ji<=(jt-1);ji++)
	w[ji-1]=rank;
      t=jt-j;
      s += (t*t*t-t);
      j=jt;
    }
  }
  if (j == n) w[n-1]=n;
}


/********************************************************************/			 
MtVector MtVecGenMRW( MtVector y, int n, MtReal mu, MtReal sigma, 
                      int rndSeed )
/********************************************************************
  Generate Multiplikative Random Walk 
  y[i+1] = y[i]*r, r = N( mu, sigma ), y0 = 100;  
  n   #samples
*********************************************************************/
{
  int  	  i;
  MtReal  dp;

  if( rndSeed >= 0 ) MtRndSeed( rndSeed );
  y[0] = 100.0;
  for( i = 1; i < n; i++ )
  {
    dp = MtRndNormGet( mu, sigma );    
    if( dp <= -.95 ) dp = -.95;
    y[i] = y[i-1] * ( 1 + dp );
  }
  return y; 
}

/********************************************************************/			 
MtVector MtVecGenChaos( MtVector y, int n, MtReal parm1, 
			int rndSeed )
/********************************************************************
  Generate Chaotic sequence 
  n   #samples
*********************************************************************/
{
  int  	  i;
  MtReal  a;

  if( rndSeed >= 0 ) MtRndSeed( rndSeed );
  a = 4;
  y[0] = 0.3;
  y[1] = 0.1;
  y[2] = 0.7;
  for(i=3;i<n;i++)
  {
    /* y[i] = a*y[i-1]*(1-y[i-1]); */
    
    y[i] = 0.5*(1 + cos(sin(10*y[i-1]) + 4*y[i-3]*cos(5*y[i-2])));
    
    y[i] += 0*MtRndNormGet(0,1);
    if(y[i] < 0) y[i]=0;
    if(y[i] > 1) y[i]=1;
         
    printf("%d %g\n",i,y[i]);          
  }
  return y;
}


/********************************************************************/
/*								    */
/*  RQueue Handling		 				    */
/*								    */
/********************************************************************/

/********************************************************************/
/*								    */
/********************************************************************/			 
MtRQueue *MtRQueueOpen( int max )
{
  MtRQueue *q = MM_malloc( sizeof ( MtRQueue ), MM_DUMMYTAG );
  if( !q ) goto fail;
  q->buffer = MM_malloc( sizeof( MtReal ) * max, MM_DUMMYTAG );
  if( !q->buffer ) goto fail;   
  q->max = max;
  MtRQueueInit( q );
  return q;

fail:   
  MtRQueueClose( q );         
  return(NULL);                                                
}

/********************************************************************/
/*								    */
/********************************************************************/			 
void MtRQueueClose( MtRQueue *q )
{
  if( q )
  {
    if( q->buffer ) MM_free( q->buffer );
    MM_free( q ); 
  } 
}

/********************************************************************/
/*								    */
/********************************************************************/			 

void MtRQueueInit( MtRQueue *q )
{
  q->in = q->out = 0;
  q->len = 0;
}

/********************************************************************/
/*								    */
/*  Vector Handling Functions 					    */
/*								    */
/********************************************************************/

/********************************************************************/
/*								    */
/********************************************************************/			 

MtVector MtVecOpen( int n )
{
  MtVector v = MM_malloc( n * sizeof( MtReal ), MM_DUMMYTAG );
  if( ! v )
  {
    printf("MtVecOpen(): MM_malloc(%ld) failed.\n",
	(long)(n*sizeof(MtReal))); 
    exit(1);
  }
  return v;
}

/********************************************************************/
/*								    */
/*  Matrix Handling Functions 					    */
/*								    */
/********************************************************************/

/********************************************************************/
/*								    */
/********************************************************************/			 
MtVector MtMatOpen( int n, int m )
{
  MtMatrix mat = MM_malloc( n * m * sizeof( MtReal ), MM_DUMMYTAG );
  if( ! mat )
  {
    printf("MtMatOpen(): MM_malloc() failed.\n"); 
  }
  return mat;
}
/********************************************************************/
/*								    */
/********************************************************************/			 
void MtMatClose( MtMatrix mat )
{
  if( mat ) MM_free( mat );
}

/********************************************************************/			 
MtMatrix MtMatMult( MtMatrix c, MtMatrix a, MtMatrix b, 
                    int n, int l, int m )
/********************************************************************
   c(n,m) = a(n,l) * b(l,m)
   -n -> transpose a
   -m -> transpose b 
*********************************************************************/
{
  MtReal  *pa, *pb, *pc, sum;
  int      i, j, k, ia, ib, ta=0, tb=0;

  if( n < 0 )
  {
    n = -n;
    ta = 1;
  } 
  if( m < 0 )
  {
    m = -m;
    tb = 1;
  }  
  pc = c;
  for(i=0; i<n; i++)
  {
    for(j=0; j<m; j++)
    {
      if( ta )
      {
        pa = a + i;
        ia = n; 
      }  
      else
      {
        pa = a + i*l;
        ia = 1;
      }  
      if( tb )
      {
        pb = b + j*l;
        ib = 1; 
      }  
      else
      { 
        pb = b + j;
        ib = m;
      }  
      sum = 0;
      for(k=0; k<l; k++)
      {
         sum += *pa * *pb;
	 pa += ia;
	 pb += ib;
      }
      *pc++ = sum;
    }
  }
  return c;  	
}

/********************************************************************/			 
double *MtMatMult_d( double *c, double *a, double *b, 
                    	int n, int l, int m )
/********************************************************************
   c(n,m) = a(n,l) * b(l,m)
   n : rows of result matrix c
   m : columns of result matrix c
   -n -> transpose a
   -m -> transpose b 
*********************************************************************/
{
  double  *pa, *pb, *pc, sum;
  int      i, j, k, ia, ib, ta=0, tb=0;

  if( n < 0 )
  {
    n = -n;
    ta = 1;
  } 
  if( m < 0 )
  {
    m = -m;
    tb = 1;
  }  
  pc = c;
  for(i=0; i<n; i++)
  {
    for(j=0; j<m; j++)
    {
      if( ta )
      {
        pa = a + i;
        ia = n; 
      }  
      else
      {
        pa = a + i*l;
        ia = 1;
      }  
      if( tb )
      {
        pb = b + j*l;
        ib = 1; 
      }  
      else
      { 
        pb = b + j;
        ib = m;
      }  
      sum = 0;
      for(k=0; k<l; k++)
      {
         sum += *pa * *pb;
	 pa += ia;
	 pb += ib;
      }
      *pc++ = sum;
    }
  }
  return c;  	
}

/********************************************************************/			 
double *MtMatCopy_d( double *a, double *b, int n, int m )
/********************************************************************
   a(n,m) = b(n,m) 
   -n -> transpose b
*********************************************************************/
{
  if(n<0||m<0)
  {
    int i,j;
    if(n<0)n=-n;
    if(m<0)m=-m;
    for(i=0;i<n;i++) 
      for(j=0;j<m;j++) 
        MT_MATEL(a,m,i,j) = MT_MATEL(b,n,j,i);
  }
  else
  {
    MtVecCopy_d(a,b,n*m);
  }
  return a;
}

/********************************************************************/			 
MtMatrix MtMatCholesky_f( MtMatrix a, int n, MtVector p )
/********************************************************************
Cholesky decomposition a = l'l
  a    IN/OUT (n,n) matrix, symmetric, pos-definite
       (only upper triangle of a required on input) 
  n    #rows = #columns
  p    OUT diagonal elements of l
*********************************************************************/
{
  int     i, j, k;
  double  sum; 

  for(i=0;i<n;i++)
  {
    for(j=i;j<n;j++)
    {
      for(sum=MT_MATEL(a,n,i,j),k=i-1;k>=0;k--)
        sum-= MT_MATEL(a,n,i,k) * MT_MATEL(a,n,j,k);
      if(i==j)
      {
        if(sum<1.0E-20)
        {
          printf(
            "ERROR: MtMatCholesky_f(): not positive-definite\n");
        }
        p[i]=sqrt(sum);
      }
      else MT_MATEL(a,n,j,i)=sum/p[i];  
    }    
  }
  return a;
}

/********************************************************************/			 
double *MtMatCholesky( double *a, int n, double *p )
/********************************************************************
Cholesky decomposition a = l'l
  a    IN/OUT (n,n) matrix, symmetric, pos-definite
       (only upper triangle of a required on input) 
  n    #rows = #columns
  p    OUT diagonal elements of l
*********************************************************************/
{
  int     i, j, k;
  double  sum; 

  for(i=0;i<n;i++)
  {
    for(j=i;j<n;j++)
    {
      for(sum=MT_MATEL(a,n,i,j),k=i-1;k>=0;k--)
        sum-= MT_MATEL(a,n,i,k) * MT_MATEL(a,n,j,k);
      if(i==j)
      {
        if(sum<1.0E-20)
        {
          fprintf(stderr,
            "ERROR: MtMatCholesky(): not positive-definite\n");
        }
        p[i]=sqrt(sum);
      }
      else MT_MATEL(a,n,j,i)=sum/p[i];  
    }    
  }
  return a;
}
  
/********************************************************************/			 
MtVector MtMatCholeskySolve( MtMatrix a, int n, MtVector p,
                             MtVector b, MtVector x )
/********************************************************************
Solve a*x=b 
  a    cholesky decomposed matrix returned by MtMatCholesky_f()
       (only lower triangle of a required on input) 
  n    matrix dimension
  p    diagonal of cholesky dec. as returned by MtMatCholesky_f()
  b    right hand side vector
  x    OUT solution vector
*********************************************************************/
{
  int     i, k;
  MtReal  sum;
  
  for(i=0;i<n;i++)
  {
    for(sum=b[i],k=i-1;k>=0;k--) 
      sum -= MT_MATEL(a,n,i,k) * x[k];
    x[i] = sum/p[i];
  }    
  for(i=n-1;i>=0;i--)
  {
    for(sum=x[i],k=i+1;k<n;k++) 
      sum -= MT_MATEL(a,n,k,i) * x[k];
    x[i] = sum/p[i];
  }    
   
  return x;
}

/********************************************************************/			 
MtVector MtMatSolve( MtMatrix a, int n, MtVector b, MtVector x )
/********************************************************************
Solve A*x=b 
  A    matrix
  b    vector
  n    matrix dimension
  x    OUT solution vector
*********************************************************************/
{
  TRCERR(("not implemented")); return NULL;
#if 0  
  int *indx=NULL; 
  float d; 
  MtMatrix a2 = MtMatOpen(n,n);

  MtMatCopy_f( a2, a, n, n );
  ALLOCMEM( indx, 2*sizeof(*indx)*n );

  /* LU decomposition of a */
  MtMatLUdcmp_f(a2,n,indx,&d);
  MtVecCopy_f( x, b, n );
  MtMatLubksb_f(a2,n,indx,x);

  MtMatClose( a2 );
  FREEMEM(indx);
  return x;
#endif
}


/*-----------------------------------------------------------------------*/
/* 									 */
/*-----------------------------------------------------------------------*/
int nrc_gaussj(float *a, int n, float *b, int m) 
{
int *indxc=NULL,*indxr=NULL,*ipiv=NULL; 
int i,icol=0,irow=0,j,k,l,ll; 
float big,dum,pivinv,tmp2; 

ALLOCARR(indxc,n+1); 
ALLOCARR(indxr,n+1); 
ALLOCARR(ipiv,n+1); 

for (j=0;j<n;j++) ipiv[j]=0; 
for (i=0;i<n;i++) 
{  
  big=0.0; 
  for (j=0;j<n;j++) 
    if (ipiv[j] != 1) 
      for (k=0;k<n;k++) 
      { 
	if (ipiv[k] == 0) 
	{ 
	  if (fabs(MT_MATEL(a,n,j,k)) >= big) 
	  { 
	    big=fabs(MT_MATEL(a,n,j,k)); 
	    irow=j; 
	    icol=k; 
	  } 
	} 
      } 
  ++(ipiv[icol]); 
  if (irow != icol) { 
    for (l=0;l<n;l++) SWAP(MT_MATEL(a,n,irow,l),MT_MATEL(a,n,icol,l),tmp2) 
      for (l=0;l<m;l++) SWAP(MT_MATEL(b,m,irow,l),MT_MATEL(b,m,icol,l),tmp2) 
  }
  indxr[i]=irow; 
  indxc[i]=icol; 
  if (MT_MATEL(a,n,icol,icol) == 0.0) 
  {
    TRCERR(("gaussj: Singular Matrix"));
    return MI_EMATH;
  }
  pivinv=1.0/MT_MATEL(a,n,icol,icol); 
  MT_MATEL(a,n,icol,icol)=1.0; 
  for (l=0;l<n;l++) MT_MATEL(a,n,icol,l) *= pivinv; 
  for (l=0;l<m;l++) MT_MATEL(b,m,icol,l) *= pivinv;

  for (ll=0;ll<n;ll++) 
    if (ll != icol) { 
      dum=MT_MATEL(a,n,ll,icol); 
      MT_MATEL(a,n,ll,icol)=0.0; 
      for (l=0;l<n;l++) MT_MATEL(a,n,ll,l) -= MT_MATEL(a,n,icol,l)*dum; 
      for (l=0;l<m;l++) MT_MATEL(b,m,ll,l) -= MT_MATEL(b,m,icol,l)*dum; 
    } 
}

for (l=n-1;l>=0;l--) 
{ 
  if (indxr[l] != indxc[l]) 
    for (k=0;k<n;k++) 
      SWAP(MT_MATEL(a,n,k,indxr[l]),MT_MATEL(a,n,k,indxc[l]),tmp2); 
} 

FREEMEM(ipiv); 
FREEMEM(indxr); 
FREEMEM(indxc); 
return 0;
}


/********************************************************************/			 
int MtMatGaussSolve( MtMatrix a, int n, MtVector b, int m )
/********************************************************************
Solve A*x=b 
  A    IN/OUT matrix, on output: inverse of matrix
  b    IN/OUT right hand side vectors: solution vectors
  n    matrix dimension
  m    number of right hand side vectors
  x    OUT solution vector
*********************************************************************/
{
  return nrc_gaussj( a, n, b, m );
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MtTstMatGaussSolve(int nIter )
{
  int i,j;
  float A[4] = {3,7,1,-2}, b[2] = {-2,3}, x[2], b2[2], A0[4],b0[2],A2[4];

  MtRndSeed( time(NULL));
  for(i=0;i<nIter;i++)
  {
    for(j=0;j<4;j++) A[j] = MtRndNormGet(0,100);
    for(j=0;j<2;j++) b[j] = MtRndNormGet(0,100);

    MtVecCopy_f( b0, b, 2 );
    MtMatCopy_f( A0, A, 2, 2 );

    MtMatGaussSolve( A, 2, b, 1 );
    MtVecCopy_f( x, b, 2 );

    MtMatMult( b2, A0, x, 2, 2, 1 );

    /*** verify with LU */
    MtMatInverseLU( A2, A0, NULL, 2, NULL, NULL, 0 );

    printf("INV(Gauss):\n"); MtMatPrint(A,2,2);
    printf("INV(LU):\n"); MtMatPrint(A2,2,2);

#if 0       
    printf("b=\n"); MtMatPrint( b0, 2, 1 );
    printf("x=\n"); MtMatPrint( x, 2, 1);
    printf("b2=\n"); MtMatPrint( b2, 2, 1 );
    printf("\n");
#endif

    printf("%d> diff = %f\n",i,(float)MtVecDistance(b0,b2,2));
  }
}



/********************************************************************/			 
void MtMatPrint_d( double *a, int n, int m, char *fformat )
/********************************************************************
*********************************************************************/
{
  int i,j;
  char *fformatAct=fformat ? fformat : "%9.5f "; 
  for(i=0;i<n;i++)
  {
    for(j=0;j<m;j++)
    {
      printf(fformatAct,MT_MATEL(a,m,i,j));
    }
  printf("\n");  
  }  
}
/********************************************************************/			 
void MtMatPrint( MtMatrix a, int n, int m )
/********************************************************************
*********************************************************************/
{
  int i,j;
  for(i=0;i<n;i++)
  {
    for(j=0;j<m;j++)
    {
      printf("%6f ",(float)MT_MATEL(a,m,i,j));
    }
  printf("\n");  
  }  
}
/********************************************************************/
/*								    */
/********************************************************************/			 
double MtNormDst( double x )
{
  MtReal y, z, N;
  y = 1.0/(1.0+0.2316419*fabs(x));
  z = 0.3989423*exp(-(x*x)/2.0);
  N = 1.0 - z*y*(0.3193815 + y*(-0.356538
    + y*(1.781478 + y*(-1.821256 + y*1.330274))));
/*  printf("x=%f, y=%f\n",(float)x,(float)N);  */
  return (x > 0) ? N : 1.0 - N;                                     
}

/********************************************************************/
/*								    */
/********************************************************************/			 
double MtNormDstInv( double y )
{
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  InvNorm
//
// Purpose:   Inverse of normal distribution, returning x(y) where
//                 y(x) = Normal(x)
//                           x
//                      = INTEGRAL exp(-u*u/2) du / sqrt(2*pi)
//                         -inf
//                      = ( 1 + erf(x/sqrt(2)) ) / 2
//
//            Rational approximations claimed to be correct to 1.15e-9
//
// History:   JS        13 Sep 2003
//            From Peter J Acklam www//home.online.no/~pjacklam/notes/invnorm
//-----------------------------------------------------------------------------
//
  
    int rcThis=0;
    double  r, s,x=0;
    if( y < 0.02425 )
    {
        s = sqrt(-2.0 * log(y));
        x =  + (2.938163982698783  +
                  (4.374664141464968  -
                  (2.549732539343734  +
                  (2.400758277161838  +
                  (0.3223964580411365 +
                   0.007784894002430293 * s) * s) * s) * s) * s)
                / (1.0                +
                  (3.754408661907416  +
                  (2.445134137142996  +
                  (0.3224671290700398 +
                   0.007784695709041462 * s) * s) * s) * s);
    }
    else if( y > 0.97575 )
    {
        s = sqrt(-2.0 * log(1.0 - y));
        x =  - (2.938163982698783  +
                  (4.374664141464968  -
                  (2.549732539343734  +
                  (2.400758277161838  +
                  (0.3223964580411365 +
                   0.007784894002430293 * s) * s) * s) * s) * s)
                / (1.0                +
                  (3.754408661907416  +
                  (2.445134137142996  +
                  (0.3224671290700398 +
                   0.007784695709041462 * s) * s) * s) * s);
    }
    else
    {
        s = y - 0.5;
        r = s * s;
        x =  s * (2.506628277459239 -
                    (30.66479806614716 -
                    (138.3577518672690 -
                    (275.9285104469687 -
                    (220.9460984245205 -
                     39.69683028665376 * r) * r) * r) * r) * r)
                  / (1.0 -
                    (13.28068155288572 -
                    (66.80131188771972 -
                    (155.6989798598866 -
                    (161.5858368580409 -
                     54.47609879822406 * r) * r) * r) * r) * r);
    }

    UT_ASSERT(!isnan(x));
rcCatch:
    return x;
}



#if 0

        /*  Returns the area under a normal distribution
            from -inf to z standard deviations
         */
        private double pnorm(double z)
        {
            double t, pd;
            t = 1 / (1 + 0.2316419 * z);
            pd = 1 - 0.3989423 *
            System.Math.Exp(-z * z / 2) *
                ((((1.330274429 * t - 1.821255978) * t
                 + 1.781477937) * t - 0.356563782) * t + 0.319381530) * t;
            /* see Gradsteyn & Rhyzik, 26.2.17 p932 */
            return (pd);
        }
#endif

#if 0
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  InvNorm
//
// Purpose:   Inverse of normal distribution, returning x(y) where
//                 y(x) = Normal(x)
//                           x
//                      = INTEGRAL exp(-u*u/2) du / sqrt(2*pi)
//                         -inf
//                      = ( 1 + erf(x/sqrt(2)) ) / 2
//
//            Rational approximations claimed to be correct to 1.15e-9
//
// History:   JS        13 Sep 2003
//            From Peter J Acklam www//home.online.no/~pjacklam/notes/invnorm
//-----------------------------------------------------------------------------
//
double InvNorm(   //   O Inverse normal (= number of standard deviations)
double  y)        // I   Argument       (= cumulative probability)
{
    double  r, s;
    if( y < 0.02425 )
    {
        s = sqrt(-2.0 * log(y));
        return  + (2.938163982698783  +
                  (4.374664141464968  -
                  (2.549732539343734  +
                  (2.400758277161838  +
                  (0.3223964580411365 +
                   0.007784894002430293 * s) * s) * s) * s) * s)
                / (1.0                +
                  (3.754408661907416  +
                  (2.445134137142996  +
                  (0.3224671290700398 +
                   0.007784695709041462 * s) * s) * s) * s);
    }
    else if( y > 0.97575 )
    {
        s = sqrt(-2.0 * log(1.0 - y));
        return  - (2.938163982698783  +
                  (4.374664141464968  -
                  (2.549732539343734  +
                  (2.400758277161838  +
                  (0.3223964580411365 +
                   0.007784894002430293 * s) * s) * s) * s) * s)
                / (1.0                +
                  (3.754408661907416  +
                  (2.445134137142996  +
                  (0.3224671290700398 +
                   0.007784695709041462 * s) * s) * s) * s);
    }
    else
    {
        s = y - 0.5;
        r = s * s;
        return  s * (2.506628277459239 -
                    (30.66479806614716 -
                    (138.3577518672690 -
                    (275.9285104469687 -
                    (220.9460984245205 -
                     39.69683028665376 * r) * r) * r) * r) * r)
                  / (1.0 -
                    (13.28068155288572 -
                    (66.80131188771972 -
                    (155.6989798598866 -
                    (161.5858368580409 -
                     54.47609879822406 * r) * r) * r) * r) * r);
    }
}
#endif
/********************************************************************/
/*								    */
/*  Black/Scholes Option Model                                      */
/*								    */
/********************************************************************/

/********************************************************************/			 
MtReal MtOBSPrice( int type, MtReal s, MtReal t, MtReal p, 
                   MtReal v, MtReal r )
/********************************************************************
  type option type ( MT_OBS_CALL | MT_OBS_PUT )
  s   strike
  t   time remaining ( years )
  p   underlying price
  v   volatility
  r   risk free interest rate
*********************************************************************/
{
  MtReal d1,d2,bs;
  
  if( type == MT_OBS_CALL )
  {
    if( v > MT_REAL_EPSILON )
    {
      d1 = (log(p/s)+(r+v*v/2.0)*t)/(v*sqrt(t));
      d2 = d1-v*sqrt(t);
      bs = p*MtNormDst(d1) - s*exp(-(r*t))*MtNormDst(d2); 
    }
    else
    {
      bs = p - s*exp(-(r*t));
    }
  }  
  else if( type == MT_OBS_PUT )
  {
    if( v > MT_REAL_EPSILON )
    {
      d1 = (log(s/p)-(r+v*v/2.0)*t)/(v*sqrt(t));
      d2 = d1+v*sqrt(t);
      bs = s*exp(-(r*t))*MtNormDst(d2) - p*MtNormDst(d1); 
    }
    else
    {
      bs = s*exp(-(r*t)) - p;
    }
  }
  else
  {
    printf("unexpected if-case.\n");
    bs = -1;
  }  
  return bs;  
}                     

/********************************************************************/
/*								    */
/*  Rnd: pseudo random number generator                             */
/*								    */
/********************************************************************/

#ifdef RND_STDLIB

/********************************************************************/
/*								    */
/********************************************************************/			 

void MtRndSeed( int seed )
{
  srand( seed );
}

/********************************************************************/
/*								    */
/********************************************************************/			 

double MtRndGet( void )
{
  return (double)rand()/((double)RAND_MAX+1);
}  

#else

/********************************************************************/
/*								    */
/********************************************************************/			 

void MtRndSeed( int seed )
{
  int   j;
  long  k;

  if( seed < 0 ) return;
  
  rndAcc = seed ^ 123456789L;
  for( j = RND_STAB_SZ + 7; j >= 0; j--)
  {
    RND_MODMUL64( rndAcc, k );
    if( j < RND_STAB_SZ ) rndStab[j] = rndAcc;  
  }   
  rndIy = rndStab[0];
}

/********************************************************************/
/*								    */
/********************************************************************/			 

double MtRndGet( void )
{
  int    j;		/* linear congruential with shuffling */
  long   k;             /* see "numerical reciepts in C": 'ran1' */
  double  tmp;
  
  RND_MODMUL64( rndAcc, k );
  j = rndIy / RND_NDIV;
  rndIy = rndStab[j];
  rndStab[j] = rndAcc;
  if(( tmp = RND_AM*rndIy ) > RND_RNMAX )
  {
    return RND_RNMAX; 
  }  
  else
  {
    return tmp;  
  }  
}  

#endif

/********************************************************************/
/*								    */
/********************************************************************/			 

int MtRndIntGet( int max )
{
  return ( max + 1 ) * MtRndGet();
}

/**********************************************************************
 *								      *	
 **********************************************************************/
 void MtRndSSeed( MtRndStream *hdl, int seed )
{
  int   j;
  long  k;

  if( seed < 0 ) return;
  
  hdl->rndIy = 0;
  
  hdl->rndAcc = seed ^ 123456789L;
  for( j = RND_STAB_SZ + 7; j >= 0; j--)
  {
    RND_MODMUL64( hdl->rndAcc, k );
    if( j < RND_STAB_SZ ) hdl->rndStab[j] = hdl->rndAcc;  
  }   
  hdl->rndIy = hdl->rndStab[0];
}

/**********************************************************************
 *								      *	
 **********************************************************************/
double MtRndSGet( MtRndStream *hdl )
{
  int    j;		/* linear congruential with shuffling */
  long   k;             /* see "numerical reciepts in C": 'ran1' */
  double  tmp;
  
  RND_MODMUL64( hdl->rndAcc, k );
  j = hdl->rndIy / RND_NDIV;
  hdl->rndIy = hdl->rndStab[j];
  hdl->rndStab[j] = hdl->rndAcc;
  if(( tmp = RND_AM*hdl->rndIy ) > RND_RNMAX )
  {
    return RND_RNMAX; 
  }  
  else
  {
    return tmp;  
  }  
}  

/**********************************************************************
 *								      *	
 **********************************************************************/
int MtRndSIntGet( MtRndStream *hdl, int max )
{
  return MAX( ( max + 1 ) * MtRndSGet( hdl ), max);
}  


/********************************************************************/
/*								    */
/********************************************************************/			 

double MtRndNormGet( double mean, double sigma )
{
  static int    iset = 0;
  static double  gset;
  double         fac, rsq, v1, v2;
  
  if( iset == 0 )
  {
    do
    {
      v1 = 2.0 * MtRndGet() - 1.0;
      v2 = 2.0 * MtRndGet() - 1.0;
      rsq = v1*v1 + v2*v2;
    }
    while(( rsq >= 1.0 )||( rsq == 0.0 ));
    fac = sqrt( -2.0*log( rsq ) / rsq );
    gset = v1 * fac;
    iset = 1;
    return sigma * ( v2 * fac ) + mean;
  }
  else
  {
    iset = 0;
    return sigma * gset + mean;
  }
}

/********************************************************************/
/*								    */
/********************************************************************/			 

double MtRndSNormGet( MtRndStream *hdl, double mean, double sigma )
{
  static int    iset = 0;
  static double  gset;
  double         fac, rsq, v1, v2;
  
  if( iset == 0 )
  {
    do
    {
      v1 = 2.0 * MtRndSGet(hdl) - 1.0;
      v2 = 2.0 * MtRndSGet(hdl) - 1.0;
      rsq = v1*v1 + v2*v2;
    }
    while(( rsq >= 1.0 )||( rsq == 0.0 ));
    fac = sqrt( -2.0*log( rsq ) / rsq );
    gset = v1 * fac;
    iset = 1;
    return sigma * ( v2 * fac ) + mean;
  }
  else
  {
    iset = 0;
    return sigma * gset + mean;
  }
}


/********************************************************************/
/*								    */
/********************************************************************/
void MtDelTab( MtTab *mt )
{
  if(mt)
  {
    FREEMEM(mt->rows);
    FREEMEM(mt->data);
    FREEMEM(mt);
  }
}

/********************************************************************/
/*								    */
/********************************************************************/
MtTab *MtCreateTab( LongInt nRow, LongInt nCol )
{
  LongInt i;
  int	rcThis=0;
  MtTab *mt = NULL;
  float **pRow,*pData;

  ALLOCOBJ(mt);
  mt->rows = NULL;
  mt->data = NULL;
  mt->nRow = mt->nRowMax = nRow;
  mt->nCol = nCol;

  ALLOCARR(mt->rows,nRow);
  if(nCol)
  {
    ALLOCARR(mt->data,nRow*nCol);
    for(i=0,pRow=mt->rows,pData=mt->data;i<nRow;i++)
    {
      *pRow++ = pData;
      pData += nCol;
    }
  }
  else
  {
    // Allow 'zero-column-dim' matrix for storing 
    // row pointers (to foreign data) only
    // (e.g. for representing sub sets of rows)
    for(i=0;i<nRow;i++) mt->rows[i]=NULL;
  }

/*rcCatch:*/
  if(rcThis)
  {
    MtDelTab( mt );
    mt = NULL;
  }
  return mt;
}

/********************************************************************/
/*								    */
/********************************************************************/
void MtDelMat( MtMat *mt )
{
  if(mt)
  {
    FREEMEM(mt->rows);
    FREEMEM(mt->data);
    FREEMEM(mt);
  }
}

/********************************************************************/
/*								    */
/********************************************************************/
MtMat *MtCreateMat( LongInt nRow, LongInt nCol )
{
  int rcThis=0;
  LongInt i;
  MtMat *mt = NULL;
  double **pRow,*pData;

  ALLOCOBJ(mt);
  mt->rows = NULL;
  mt->data = NULL;
  mt->nRow = mt->nRowMax = nRow;
  mt->nCol = nCol;

  ALLOCARR(mt->rows,nRow);
  ALLOCARR(mt->data,nRow*nCol);
  for(i=0,pRow=mt->rows,pData=mt->data;i<nRow;i++)
  {
    *pRow++ = pData;
    pData += nCol;
  }

/*rcCatch:*/
  if(rcThis)
  {
    MtDelMat( mt );
    mt = NULL;
  }
  return mt;
}

/********************************************************************/
/*								    */
/********************************************************************/
MtMat *MtCopyMat( MtMat *dst, MtMat *src )
{
  int rcThis=0;
  UT_ASSERT(src->nCol==dst->nCol&&src->nRow==dst->nRow);
  memcpy(dst->data,src->data,
         sizeof(dst->data[0])*src->nRow*src->nCol);
rcCatch:
  return dst;
}

/********************************************************************/
/*								    */
/********************************************************************/
MtMat *MtSetMat( MtMat *M, double val )
{
  int i;
  for(i=0;i<M->nRow*M->nCol;i++) M->data[i] = val;
  return M;
}

/********************************************************************/
/*								    */
/********************************************************************/			 
MtTab *MtGenRndMultNorm( int m,
    			 MtVector vMean, MtMatrix mSigma,
    		    	 int nSmp,
			 int seed )
/** create nSmp samples of */
{
  int i,k;
  float	*vy,*vx=NULL,*L=NULL,*vd=NULL;
  MtTab *mt=NULL;

  if( seed >= 0 ) MtRndSeed( seed );

  mt = MtCreateTab(nSmp,m);
  ALLOCARR(vx,m); ALLOCARR(vd,m ); ALLOCARR(L,m*m);
  MtMatCopy_f(L,mSigma,m,m);

#if 0
  for(k=0;k<m*m;k++) 
    if(isnan(L[k]))
    {
      TRCERR(("Assertion"));
    }
#endif

  MtMatCholesky_f(L,m,vd);

#if 0
  for(k=0;k<m*m;k++) 
    if(isnan(L[k]))
    {
      TRCERR(("Assertion"));
    }
#endif

  for(i=0;i<nSmp;i++)
  {
    vy = mt->rows[i];
    for(k=0;k<m;k++)
    {
      vx[k] = MtRndNormGet(0,1);
    }
    MtVecAdd( MtMatMult( vy, L, vx, -m, m, 1 ), vMean, m );
  }

  FREEMEM(vx); FREEMEM(L); FREEMEM(vd);
  return mt;
}

/********************************************************************/
/*								    */
/********************************************************************/			 

void MtRnd1Seed( int seed )
{
  srand( seed );
}  
			 

/********************************************************************/
/*								    */
/*  Dst: probability distributions                                  */
/*								    */
/********************************************************************/			 

#if 0
/********************************************************************/			 
MtReal MtDstInverse( MtDst *dst, MtReal u )
/********************************************************************
*********************************************************************/
{
  int       j;
  MtReal    y;

  y = dst->nSmp*u;  
  j = (int)y;
  return dst->vDist[j];  
}
#endif

#if 1
/********************************************************************/			 
MtReal MtDstInverse( MtDst *dst, MtReal u )
/********************************************************************
*********************************************************************/
{
  int       n,j;
  MtReal    x,y;
  MtVector vDist;

  TRCERR(("Deprecated function: MtDstDistr(). "
	  "Use MtDstDistr2() instead\n"));

  n = dst->nSmp; 
  vDist = dst->vDist; 
  y = n*u;
  
  if( y < 0.5 )            /* left tail ? */
  {
    j = (int)( y - 0.5 );
    x = vDist[j] + (vDist[j+1] - vDist[j]) * ( y - j - 0.5 );
  }
  else if( y > n - 0.5 )   /* right tail ? */
  {
    j = n - 2;
    x = vDist[j] + (vDist[j+1] - vDist[j]) * ( y - j - 0.5 );
  }
  else                     /* middle part */
  {  
    j = (int)( y - 0.5 );
    x = vDist[j] + (vDist[j+1] - vDist[j]) * ( y - j - 0.5 );
  }
  return x;   
}
#endif

#if 0
#define BSEARCH_REAL( f, x, i, j ) \
{ \
  int k_; \
  i=0;j=n-1; \
  while( i < j-1 ) \
  { \
    k_ = (i + j) >> 1; \
    if( x < f[k_] ) \
      j = k_; \
    else \
      i = k_; \
  } \
}
#endif

#define BSEARCH_REAL( f, n, x, i, j ) \
{ \
  register int k_; \
  i=0;j=n-1; \
  while( i < j-1 ) \
  { \
    k_ = (i + j) >> 1; \
    if( x - f[k_] > MT_REAL_EPSILON ) \
      i = k_; \
    else \
      j = k_; \
  } \
}

/********************************************************************/			 
MtReal MtDstDistr( MtDst *dst, MtReal x, int use_old )
/********************************************************************
*********************************************************************/
{
  int i,j,n=dst->nSmp;
  MtReal x1,x2,y,*vDist = dst->vDist,dx;

  BSEARCH_REAL( vDist, n, x, i, j );

  if(!use_old)
    TRCERR(("Deprecated function: MtDstDistr()\n"));

#if 0
  while ((j<n-1) && 
      (fabs(vDist[j+1] - vDist[j]) < MT_REAL_EPSILON))
  {
    j++;
  }
#endif

  x1=vDist[i]; x2=vDist[j];
  dx = x2-x1;
  if(dx < MT_REAL_EPSILON ) 
  {
    y = (float)i/(float)n;
  }
  else
  {
    y = (MtReal)i + ((float)j-(float)i)*(x-x1)/dx;
  }

  return y/((MtReal)n);   
}
#if 0
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static MtDst *MtIDstCreateNormalRegul( int nSmp, double mu, double sig,
    			  double domMin, double domMax,
			  double domEps )
{
  int j,i,rcThis=0,nStep,nSmpMax,nSmpAct;
  float *smp=NULL;
  MtDst *dst=NULL;
  double sum,sumOld,f,x,domMin2,domSpan2,domSpan=domMax-domMin,delta,x2;

  USE(rcThis);

  UT_ASSERT(domMax>domMin);

  domSpan2 = 1.5*domSpan;	// Tail heuristic (TODO)
  nStep = nSmpMax = 1.5*nSmp;

  delta = 1.0/(double)nStep;

  ALLOCARR(smp,nSmpMax);

  domMin2 = domMin - (domSpan2-domSpan)/2;
  TRC1(("domMin=%g, domSpan=%g -> domMin2=%g, domSpan2=%g\n",
	domMin,domSpan, domMin2, domSpan2));
  sum=sumOld=0;
  j=0;
  for(i=0;i<nStep;i++)
  {
    x = domMin2 + domSpan2  *  i/((double)nStep-1.0);
    f = MtNormalDens(x,mu,sig);
    sum+=f;
    if(sum-sumOld > delta)
    {
      sumOld = sum;
      x2 = x;
      if(x2<domMin) x2=domMin;
      if(x2>domMax) x2=domMax;
      UT_ASSERT(j<nSmpMax);
      smp[j++] = x2;
      TRC1(("%3d -> %10.4f %10.4f\n",j-1,x,x2));
    }
  }
  nSmpAct = j;
  UT_ASSERT(nSmpAct >= nSmpMax-1);

  dst = MtDstOpen2(smp,nSmpAct,10,domMin,domMax,domEps,MT_SORTED);
  UT_ASSERT(dst);

rcCatch:
  FREEMEM(smp);
  return dst;
}
#endif


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtDelTransTab( MtDstTransTab **p_mtt )
{
  int rcThis=0;
  MtDstTransTab *mtt=*p_mtt;
  if(mtt)
  {
    FREEMEM(mtt->lookup);
    FREEMEM(mtt);
    *p_mtt = NULL;
  }
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
double MtDstLookupTransTab( MtDstTransTab *mtt, double x )
{
  int rcThis=0,i=0, n = mtt->nSteps;

  USE(rcThis);

  i = (int)( (n-1) * (( x - mtt->dstSrc->domMin ) / 
			mtt->srcDomSpan) + 0.5);
  UT_ASSERT_FATAL(i>=0&&i<n);
//rcCatch:
  return mtt->lookup[i];
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/



/********************************************************************/			 
MtReal MtRndGenericGet( MtDst *dst )
/********************************************************************
Return pseudo random number drawn from the generic distribution
described by <dst>.
*********************************************************************/
{
  return 1;
}

#if 0
/********************************************************************/			 
int MtDstTransform( MtDst *dst, MtReal *vy, int n )
/********************************************************************
Generate Transformation function for empirical distribution <dst>
to be used for generating random realizations of this distr.
by the uniform distr. MtRndGet()

  dst        empirical distr.
  vy   (OUT) inverse distribution function [0,1] -> [dst.min,dst.max]
  n          number of function points
*********************************************************************/
{
  int        j,i,ix,ix0,ix1;
  MtVector   vF;
  MtReal     dx,dy,y,sum,x,r,f0,f1,x0,x1;
 
  vF = MtVecOpen( dst->nFreq );
  sum = 0;
  for(i=0;i<dst->nFreq;i++)
  {
    sum += dst->freq[i];
    vF[i] = sum; 
  }
    
  dx = dst->dx;
  dy = sum/(MtReal)n;
  r = dst->nFreq*dx/sum;
  y = 0;
  for(i=0;i<n;i++)
  {
    for(j=0;j<dst->nFreq;j++)
    {
      if( vF[j] > y )
      { 
        ix0 = j-1;
        ix1 = j;       
        f0 = ix0*dx;
        f1 = ix1*dx;
        x0 = ix0 >= 0 ? vF[ix0] : 0;
        x1 = vF[ix1]; 
        x = f0 + (f1-f0)/(x1-x0) * (y-x0); 
        vy[i] = x + dst->xMin;
        break;                                  
      }
    }    
    y += dy;  
  }  
  MtVecClose( vF );
}

/********************************************************************/			 
MtReal MtRndGenericGet( MtReal *vInvDstr, int n )
/********************************************************************
Return pseudo random number drawn from the generic distribution
described by <vInvDstr>

  vInvDstr    invesre of the empiric distribution function as
              generated by MtDstTransform()
  n           number of samples of <vInvDstr>
*********************************************************************/
{
  int    ix,ix0,ix1;
  MtReal y,x,dx,x0,f1,f0;
  
  x = MtRndGet();
  ix = n*x;
  dx = 1/(MtReal)n;
  if(ix<0) ix = 0;
  if(ix>=n) ix = n-1;
  if( ix < n - 1)
  {      
    ix0 = ix;
    ix1 = ix+1;
  }
  else
  {
    ix0 = ix-1;
    ix1 = ix;
  }
  x0 = ix0*dx;
  f0 = vInvDstr[ix0];
  f1 = vInvDstr[ix1];
  y = f0 + ((f1-f0)/dx) * (x-x0);  
  return y;
}
#endif

/********************************************************************/
static MtReal MtDstKernelGauss( MtReal x )
/********************************************************************			 
*********************************************************************/			 
{
  return exp( -x*x / 2 ) /MT_SQRT2PI;
}

/********************************************************************/
MtReal MtKernParzen( MtReal x, MtReal h )
/********************************************************************			 
*********************************************************************/			 
{
  return exp(-x*x / (2*h*h)) / MT_SQRT2PI;
}

/********************************************************************/
double MtGauss( double x )
/********************************************************************			 
*********************************************************************/			 
{
  return exp( (-x*x) / 2. );
}

/********************************************************************/
double MtNormalDens( double x, double mu, double sig )
/********************************************************************			 
*********************************************************************/			 
{
  return (1./(sig*MT_SQRT2PI)) * 
    		exp( (- (x-mu)*(x-mu)) / (2.*sig*sig) );
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtMultStatPrint2( int dim, int nSmp, double *vMean, 
    			double *vSigma, double *vMin, double *vMax,
			double *vCov, double *vCorr,
			char *fformat )
{
  int i,j,k;
  char *fformatAct=fformat ? fformat : "%11.6g "; 
  // TODO: auto detect format string from data
  double r;
  printf("dimensions=%d, samples=%d\n",dim,nSmp);
  if(vMean)
  {
    printf("Mean:  ");
    for(k=0;k<dim;k++) printf(fformatAct,vMean[k]);
    printf("\n");
  }
  if(vSigma)
  {
    printf("Sigma: ");
    for(k=0;k<dim;k++) printf(fformatAct,vSigma[k]);
    printf("\n");
  }
  if(vMin)
  {
    printf("Min:  ");
    for(k=0;k<dim;k++) printf(fformatAct,vMin[k]);
    printf("\n");
  }
  if(vMax)
  {
    printf("Max:  ");
    for(k=0;k<dim;k++) printf(fformatAct,vMax[k]);
    printf("\n");
  }
  if(vCorr)
  {
  printf("Correlation matrix=\n");
  for(i=0;i<dim;i++)
  {
    for(j=0;j<dim;j++)
    {
      r = MT_MATEL(vCorr,dim,i,j);
      printf("%8.5f",r);
    }
    printf("\n");
  }
  printf("\n");
  }
  if(vCov)
  {
  printf("Covariance matrix=\n");
  for(i=0;i<dim;i++)
  {
    for(j=0;j<dim;j++)
    {
      r = MT_MATEL(vCov,dim,i,j);
      printf(fformatAct,r);
    }
    printf("\n");
  }
  printf("\n");
  }
  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtMultStatPrint( int dim, int nSmp, double *vMean, 
    			double *vSigma, double *vCov, double *vCorr,
			char *fformat )
{
  return MtMultStatPrint2( dim, nSmp, vMean, vSigma, NULL, NULL, vCov,
      			   vCorr, fformat );
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtMultStatPrintTwo( MtMultStat *S1, MtMultStat *S2 )
{
  int rcThis=0,i,j,k,dim=S1->dim,ne;
  double f1,f2,esum=0,dr;

  if(S1->dim != S2->dim) 
    TRCERRR(("dimension mismatch\n"),MI_ERROR);

  esum=0;ne=0;
  printf("============ Multivariate Statistics ==================\n");
  printf("dimensions=%d, samples=%d / %d\n", dim,S1->nSmp,S2->nSmp);
  printf("---- mean  \n");
  for(k=0;k<dim;k++) 
  {
    f1 = S1->mean[k]; f2 = S2->mean[k]; 
    dr = MtRelDiff(f1,f2); 
    esum += dr; ne++;
    printf("%3d %11.6f %11.6f  %9.5f%%\n",k,f1,f2,100.*dr);
  }
  printf("\n");
  printf("---- sigma  \n");
  for(k=0;k<dim;k++) 
  {
    f1 = S1->sigma[k]; f2 = S2->sigma[k]; 
    dr = MtRelDiff(f1,f2); 
    esum += dr; ne++;
    printf("%3d %11.6f %11.6f  %9.5f%%\n",k,f1,f2,100.*dr);
  }
  printf("\n");

  printf("---- correlation matrix\n");
  for(i=0;i<dim;i++)
  {
    for(j=0;j<dim;j++)
    {
      f1 = MT_MATEL(S1->corr,dim,i,j);
      f2 = MT_MATEL(S2->corr,dim,i,j);
      dr = MtRelDiff(f1,f2); 
      esum += dr; ne++;
      printf("C(%2d,%2d) %11.6f %11.6f  %9.5f%%\n",
	  i,j,f1,f2,100.*dr);
    }
    printf("\n");
  }
  esum /= (double)ne;
  printf("---- total rel. diff = %12.8f\n",esum);

rcCatch:
  return rcThis;
}

#if 0
/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int MtAddMultStat( 
    	// dst = dst + src
    	MtMultStat *dst, MtMultStat *src )
{
  double n1 = dst->nSmp, n2= src->nSmp, n=n1+n2,w1,w2; 
  int dim=dst->dim,k,i,j;
  double *M1=dst->mean, *M2=src->mean,
	 *C1=dst->cov,  *C2=src->cov;
  dst->nSmp = n;

  w1 = n1/n; w2 = n2/n;

  UT_ASSERT(dst->dim == src->dim);

  // Combine New Covariance
  if(C1&&C2)
  {
    for(i=0;i<dim;i++)
    {
      for(j=i;j<dim;j++)
      {
	MT_MATEL(C1,i,j) = 
	  MT_MATEL(C1,j,i ) = 
	    (MT_MATEL(C1,i,j) + M1[i]*M1[j]) 
      }
    }
  }



  for(i=0;i<dim;i++)
  {
    dst->mean[i] = 
       ((double)n1*dst->mean[i] + (double)n2*src->mean[i])/(double)n;
    xxx
  }

  // Combine Covariances
  if(src->cov)
  {
    for(i=0;i<dim;i++)
    {
      for(j=i;j<dim;j++)
      {
	MT_MATEL(dst->cov
      }
    }
  }
}
#endif
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtMultStatGet( 
    		double **smp,  // IN: array of pointers to sample points
    		float **smp_f, 	  // type float samples (deprecated)
    		int dim,     // IN: sample point dimension, 
		int n,	     // IN: number of samples 
		double *mean,	// OUT: mean vector (dim x 1)
		double *sigma,  // OUT: std. deviations (dim x 1)
		double *cov,	// OUT: covariance matrix ( dim x dim )
		double *corr,	// OUT: correlation matrix (dimxdim)
		int (*pgiCallback)( float fracCompleted,
		  		    void **pState )
			// Optional progress indicator callback
    		)
{
  int i,j,k,rcThis=0,rc;
  double sum,c,s,r,tmp;
  void	*pgiState=NULL;

  // mean vector
  for(k=0;k<dim;k++) mean[k] = 0;
  
  if(smp_f)
    for(i=0;i<n;i++) for(k=0;k<dim;k++) mean[k] += smp_f[i][k];
  else
    for(i=0;i<n;i++) for(k=0;k<dim;k++) mean[k] += smp[i][k];

  for(k=0;k<dim;k++) mean[k] /= (double)n;

  // covariance matrix
  if(cov)
  {
    if(smp_f)
    {
      for(i=0;i<dim;i++)
      {
	if(pgiCallback) 
	{
	  rc = pgiCallback(i/(double)(dim-1),&pgiState);
	  if(rc) raiseRc(rc);
	}
	for(j=i;j<dim;j++)
	{
	  sum=0;
	  for(k=0;k<n;k++)
	  {
	    sum += ( smp_f[k][i] - mean[i] ) * ( smp_f[k][j] - mean[j] ) ;
	  }
	  MT_MATEL(cov,dim,i,j) = 
	    MT_MATEL(cov,dim,j,i) = sum/(double)n;  /* TODO: unbiased */
	}
      }
    }
    else
    {
      for(i=0;i<dim;i++)
      {
	if(pgiCallback) 
	{
	  rc = pgiCallback(i/(double)(dim-1),&pgiState);
	  if(rc) raiseRc(rc);
	}
	for(j=i;j<dim;j++)
	{
	  sum=0;
	  for(k=0;k<n;k++)
	  {
	    sum += ( smp[k][i] - mean[i] ) * ( smp[k][j] - mean[j] ) ;
	  }
	  MT_MATEL(cov,dim,i,j) = 
	    MT_MATEL(cov,dim,j,i) = sum/(double)n;
	}
      }
    }
  }

  if(sigma)
  {
    // standard deviations (diagonal)
    if(cov)
    {
      for(k=0;k<dim;k++)
      {
	c = MT_MATEL(cov,dim,k,k);
	sigma[k] = sqrt(c);
      }
    }
    else
    {
      if(smp_f)
      {
	for(j=0;j<dim;j++)
	{
	  for(sum=0,k=0;k<n;k++)
	  {
	    tmp = smp_f[k][j] - mean[j];
	    sum += tmp*tmp ;
	  }
	  sigma[j]= sqrt(sum/(double)n);
	}      
      }
      else
      {
	for(j=0;j<dim;j++)
	{
	  for(sum=0,k=0;k<n;k++)
	  {
	    tmp = smp[k][j] - mean[j];
	    sum += tmp*tmp ;
	  }
	  sigma[j]= sqrt(sum/(double)n);
	}
      }
    }
  }


  // correlation matrix
  if(corr)
  {
    UT_ASSERT(sigma!=NULL);
    {
      for(i=0;i<dim;i++)
	for(j=i;j<dim;j++)
	{
	  s = sigma[i]*sigma[j];
	  if(s<MT_REAL_EPSILON)
	  {
	    TRC1(("****WARNING: zero standard deviation detected\n"));
	    r = 0.0;
	  }
	  else
	  {
	    r = MT_MATEL(cov,dim,i,j) / s;
	  }
	  MT_MATEL(corr,dim,i,j) = 
	    MT_MATEL(corr,dim,j,i) = r;
	}
    }
  }
rcCatch:
  if(pgiState) FREEMEM(pgiState);
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtCovToCorr(double *corr, int dim, double *cov, double *sigma)
{
  int i,j;
  double s,r;

  for(i=0;i<dim;i++)
    for(j=i;j<dim;j++)
    {
      s = sigma[i]*sigma[j];
      if(s<MT_REAL_EPSILON)
      {
	TRC1(("****WARNING: zero standard deviation detected\n"));
	r = 0.0;
      }
      else
      {
	r = MT_MATEL(cov,dim,i,j) / s;
      }
      MT_MATEL(corr,dim,i,j) = 
	MT_MATEL(corr,dim,j,i) = r;
    }
  return 0;
}
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtCorrToCov(double *cov, int dim, double *corr, double *sigma)
{
  int i,j;

  for(i=0;i<dim;i++)
    for(j=i;j<dim;j++)
    {
      MT_MATEL(cov,dim,i,j) = 
	MT_MATEL(cov,dim,j,i) = 
	  MT_MATEL(corr,dim,j,i) * sigma[i] * sigma[j];
    }
  return 0;
}
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
double MtMultNormDensity(
    		// Evaluate multivariate normal density at point x 
    		MtMultStat *mst, 
		    // Handle to MtMultStat object
		double *x, 
		    // Target point (dim x 1)-vector
		double *mean,
		    // Optional mean vector
		double *sigma,
		    // Optional standard deviations
		unsigned opt 
		    // options:
		    // 	 MT_NONORM | MT_LOGDENS 
		)
{
  int n=mst->dim,i,rcThis=0;
  double g[1],f=0,*z,det,lf,*Q;

  if(!mean) mean=mst->mean;

  if( sigma || (opt & MT_NORMSIG ))
  {
    if(!(opt & MT_NONORM)) 
      TRCERRR(("not implemented\n"),MI_ERROR);
    if(!mst->icorrIsValid)
    {
      if(!mst->icorr)
      {
	ALLOCARR(mst->icorr,n*n);
      }
      if(!mst->corr)
      {
	UT_ASSERT(mst->cov && mst->sigma);
	ALLOCARR(mst->corr,n*n);
	MtCovToCorr(mst->corr,mst->dim,mst->cov,mst->sigma);
      }
      MtMatCopy_fast(mst->icorr,mst->corr,n,n);
      mst->corrDet = MtMatInverse(mst->icorr,n);
      mst->icorrIsValid = TRUE;
    }
    det = mst->corrDet; /* TODO: correct determinant -> covariance matrix */
    Q = mst->icorr;
  }
  else
  {
    if(!mst->icovIsValid)
    {
      if(!mst->icov)
      {
	ALLOCARR(mst->icov,n*n);
      }
      if(!mst->cov)
      {
	UT_ASSERT(mst->corr && mst->sigma);
	ALLOCARR(mst->cov,n*n);
	MtCorrToCov(mst->cov,mst->dim,mst->corr,mst->sigma);
      }
      MtMatCopy_fast(mst->icov,mst->cov,n,n);
      mst->covDet = MtMatInverse(mst->icov,n);
      mst->icovIsValid = TRUE;
    }
    det = mst->covDet;
    Q = mst->icov;
  }

#if 0
  if(!mst->xtmp)
  {
    ALLOCARR(mst->xtmp,n);
    ALLOCARR(mst->ytmp,n);
  }
  z = mst->xtmp;
  double *ytmp = mst->ytmp; 
#else  
  // use local buffer for thread safety !
  double xtmp[1000],
	 ytmp[1000];
  UT_ASSERT_FATAL(n<=ARRAY_LENGTH(xtmp));
  z = xtmp;
#endif

  // Center Mean
  for(i=0;i<n;i++) z[i] = x[i] - mean[i];
  
  // Multiply optional variances
  if(sigma)
  {
    for(i=0;i<n;i++) z[i] /= sigma[i];
  }

  // Compute z' S^(-1) z
  MtMatMult_d(ytmp,Q,z,n,n,1);
  MtMatMult_d(g,z,ytmp,-1,n,1);

  lf = -0.5*g[0];

  if( opt & MT_LOGDENS )
  {
    f =  opt & MT_NONORM ? 	
      		lf : 
      		lf + (-log(det))/2.0-(log(2*MT_PI)*(double)n)/2.0; 
  }
  else
  {
    f = exp(lf);
    f = opt & MT_NONORM ? 
      		f : 
		f / (pow(2*MT_PI,(double)n/2.0)*sqrt(det));
  }
rcCatch:
  return f;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtPrintMultStat( MtMultStat *adm, char *fformat, int vl )
{
  return MtMultStatPrint2( adm->dim, adm->nSmp, adm->mean, 
      			  adm->sigma, adm->min, adm->max,
			  adm->cov, adm->corr,
			  fformat );

}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtCopyMultStat( MtMultStat *dst, MtMultStat *src )
{
  int n=src->dim,rcThis=0; 
  UT_ASSERT(dst->dim==src->dim);
  // whole structure copy
  dst->dim = src->dim;
  dst->nSmp = src->nSmp;
  dst->icovIsValid = src->icovIsValid;
  dst->covDet = src->covDet;
  dst->icorrIsValid = src->icorrIsValid;
  dst->corrDet = src->corrDet;
  // copy arrays
  if(src->mean) MtVecCopy_d(dst->mean,src->mean,n);
  if(src->sigma) MtVecCopy_d(dst->sigma,src->sigma,n);
  if(src->min) MtVecCopy_d(dst->min,src->min,n);
  if(src->max) MtVecCopy_d(dst->max,src->max,n);
  if(src->cov) MtVecCopy_d(dst->cov,src->cov,n*n);
  if(src->corr) MtVecCopy_d(dst->corr,src->corr,n*n);
  if(src->icov && src->icovIsValid) 
    MtVecCopy_d(dst->icov,src->icov,n*n);
  if(src->icorr && src->icorrIsValid) 
    MtVecCopy_d(dst->icorr,src->icorr,n*n);
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtDelMultStat( MtMultStat *mst )
{
  if(mst)
  {
    FREEMEM(mst->ytmp);
    FREEMEM(mst->xtmp);
    FREEMEM(mst->icorr);
    FREEMEM(mst->icov);
    FREEMEM(mst->mean);
    FREEMEM(mst->min); FREEMEM(mst->max);
    FREEMEM(mst->sigma);
    FREEMEM(mst->cov);
    FREEMEM(mst->corr);
    FREEMEM(mst);
  }
  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
MtMultStat *MtCreateMultStat( 
    		double **smp,  // IN: array of pointers to sample points
    		float **smp_f, // for float samples (deprecated) 
    		int dim,     // IN: sample point dimension, 
		int n,	     // IN: number of samples 
		int (*pgiCallback)( float fracCompleted, 
		  		    void **pState ),
			// Optional progress indicator callback
		unsigned opt,
		MtMultStat **p_mst
    		)
{
  int rcThis=0,rc,doSel=opt & MT_SELECT,k,i,
  	doSig=1,doCov=1,doCorr=1,doMinmax=0,doNotCorr=0;
  MtMultStat *mst = p_mst ? *p_mst : NULL;
  double *min=NULL,*max=NULL,x;

  if(doSel)
  {
    doSig = opt & MT_SIG;
    doCov = opt & MT_COV;
    doCorr = opt & MT_CORR;
    doNotCorr = opt & MT_NOCORR;
    doMinmax = opt & MT_MINMAX;
  }

  if(!mst)
  {
    ALLOCOBJ0(mst);
    ALLOCARR(mst->mean,dim);
    if(doSig) ALLOCARR(mst->sigma,dim);
    if(doCov) ALLOCARR(mst->cov,dim*dim);
    if(doCorr) ALLOCARR(mst->corr,dim*dim);
    if(opt & MT_MINMAX )
    {
      ALLOCARR(mst->min,dim); 
      ALLOCARR(mst->max,dim);
    }
  }
  else
  {
    if(mst->dim!=dim)
    {
      TRCERRR(("dimension mismatch in MtMultStat"),MI_ERROR);
    }
  }
  mst->nSmp = n;
  mst->dim = dim;

  if(smp||smp_f)
  {
    mst->icovIsValid =  FALSE;
    mst->icorrIsValid =  FALSE;
    rc = MtMultStatGet(smp,smp_f,dim,n,
		       mst->mean,mst->sigma,mst->cov,
		       doNotCorr ? NULL : mst->corr,
		       pgiCallback );
    if(rc) 
    {
      MtDelMultStat(mst);
      mst = NULL;
      raiseRc(rc);
    }
    if(opt & MT_MINMAX)
    {
      min=mst->min;
      max=mst->max;
      if(smp_f) 
      {
	for(k=0;k<dim;k++) 
	{
	  min[k] = n > 0 ? smp_f[0][k] : FLT_MAX; 
	  max[k] = n > 0 ? smp_f[0][k] : -FLT_MAX;
	}
	for(k=0;k<dim;k++) 
	{
	  for(i=0;i<n;i++) 
	  {
	    x = smp_f[i][k];
	    if(x>max[k]) max[k] = x;
	    if(x<min[k]) min[k] = x;	   
	  }
	}
      }
      else
      {
	for(k=0;k<dim;k++) 
	{
	  min[k] = n > 0 ? smp[0][k] : DBL_MAX; 
	  max[k] = n > 0 ? smp[0][k] : -DBL_MAX;
	}
	for(k=0;k<dim;k++) 
	{
	  for(i=0;i<n;i++) 
	  {
	    x = smp[i][k];
	    if(x>max[k]) max[k] = x;
	    if(x<min[k]) min[k] = x;	   
	  }
	}
      }
    }    
  }

  if(mst && (opt&MT_ICOV))
  {
    // trigger creation of inverse covariance matrix
    UT_ASSERT(mst->mean&&mst->sigma&&mst->cov);
    double f = 
      MtMultNormDensity(mst,mst->mean,NULL,NULL,
		      MT_LOGDENS|MT_NONORM);
    UT_ASSERT(isfinite(f));
    UT_ASSERT(mst->icovIsValid);
  }

rcCatch:

  if(p_mst) 
  {
    *p_mst = mst;
  }
  return mst;
}


/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int MtMultStatGet2(
    		double **smp,  // IN: array of pointers to sample points
    		float **smp_f, 	  // type float samples (deprecated)
		double *wgt,	// IN: Sample weights: 
				   //   must be normalized to 1 ! 
    		int dim,     // IN: sample point dimension, 
		int n,	     // IN: number of samples 
		double *mean,	// OUT: mean vector (dim x 1)
		double *cov,	// OUT: covariance matrix ( dim x dim )
		unsigned opt
		)
{
  int i,j,k,rcThis=0;
  double sum,wsqsum,sumdiv;

  UT_ASSERT(smp_f!=NULL&&smp==NULL);	// TODO: implement for type double

  // mean vector
  for(k=0;k<dim;k++) mean[k] = 0;

  if(smp_f)
    for(i=0;i<n;i++) for(k=0;k<dim;k++) mean[k] += wgt[i] * smp_f[i][k];
  else
    for(i=0;i<n;i++) for(k=0;k<dim;k++) mean[k] += wgt[i] * smp[i][k];

  // covariance matrix
  if(cov)
  {
    for(k=0,wsqsum=0;k<n;k++) wsqsum+=wgt[k]*wgt[k];
    sumdiv = 1.0-wsqsum;

    if(smp_f)
    {
      for(i=0;i<dim;i++)
      {
	for(j=i;j<dim;j++)
	{
	  if((opt & MT_NOCOV)  && i!=j) continue;
	  sum=0;
	  for(k=0;k<n;k++)
	  {
	    sum += wgt[k] * ( smp_f[k][i] - mean[i] ) * 
	      		    ( smp_f[k][j] - mean[j] ) ;
	  }
	  MT_MATEL(cov,dim,i,j) = 
	    MT_MATEL(cov,dim,j,i) = 
	      fabs(sumdiv) > MT_REAL_EPSILON ? sum/sumdiv : 0;
	}
      }
    }
    else
    {
      UT_ASSERT(FALSE); // TODO: type double 
    }
  }

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtGetWeightedMeanCov(
    		double **smp,  // IN: array of pointers to sample points
    		float **smp_f, 	  // type float samples (deprecated)
		double *wgt,	// IN: Sample weights: 
				   //   must be normalized to 1 ! 
    		int dim,     // IN: sample point dimension, 
		int n,	     // IN: number of samples 
		double *mean,	// OUT: mean vector (dim x 1)
		double *cov,	// OUT: covariance matrix ( dim x dim )
		unsigned opt
		)
{
  return MtMultStatGet2( smp, smp_f, wgt, dim, n, mean, cov, opt );
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtSetMultStat(MtMultStat *S,
    		   double *mean, 
		   double *sigma,
		   double *corr,
		   double *cov )
{
  int rcThis=0,n=S->dim;
  S->icovIsValid =  FALSE;
  S->icorrIsValid =  FALSE;
  if(mean) 
  {
    MtVecCopy_d(S->mean,mean,n);
  }
  if(sigma) 
  {
    if(!S->sigma) ALLOCARR(S->sigma,n);
    MtVecCopy_d(S->sigma,sigma,n);
  }
  if(corr) 
  {
    if(!S->corr) ALLOCARR(S->corr,n*n);
    MtVecCopy_d(S->corr,corr,n*n);
  }
  else
  {
    FREEMEM(S->corr);
  }
  if(cov) 
  {
    if(!S->cov) ALLOCARR(S->cov,n*n);
    MtVecCopy_d(S->cov,cov,n*n);
  }
  else
  {
    FREEMEM(S->cov);
  }
  return rcThis;
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
double MtShrinkCovarianceToTarget( 
    		double *Z,  // Sample empirical covariance matrix
		double *T,  // Target covariance
		int dim,    // dimension	
		double nSamples,   // number of samples 
		double alpha  // shrinkage strength factor
			      // (mult. of optimal strength lambda )
					      	
    )
{
  int i,j;
  double lambda1,lambda =  nSamples > 0 ? 
      pow(2,-(double)nSamples / dim) : 1;

  lambda *= alpha;
  if(lambda<0) lambda = 0;
  if(lambda>1) lambda = 1;
  lambda1 = 1-lambda;

  for(i=0;i<dim;i++)
    for(j=0;j<dim;j++)
    {
      MT_MATEL(Z,dim, i,j) = 
	lambda1 * MT_MATEL(Z,dim,i,j) +
	lambda * MT_MATEL(T,dim,i,j);
    }


//rcCatch:
    return lambda;

}
#if 0
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtShrinkCovariance( double **smp,	// Input data matrix
    			float **smp_f,
    			double *wgt,	// Optional data weights
			int	dim,
			int	n,	// Number of data points

			double *pLambda,
			double *mean,	// Out: mean
			double *cov	// Out: (shrinkage) covariance
    			)
{					
  int rcThis=0,rc;
  double lambda = pLambda ? *pLambda : -1;
  double *XS=NULL, *cov0=NULL,s,h1;

  ALLOCARR(C,dim*dim);
  ALLOCARR(XS,n*dim);

  // 
  // Get empirical covariance
  //
  rc = wgt ?  MtMultStatGet2( smp, smp_f, wgt, dim, n, mean, cov0 ) :
	      MtMultStatGet( smp, smp_f, dim, n, mean, NULL, cov0, 
		  	     NULL, NULL );
  if(rc) raiseRc(rc);

  //
  // Center and scale data
  //
  for(j=0;j<dim;j++)
  {
    s = MT_MATEL(cov0,j,j);
    UT_ASSERT(!(s<0));
    s = sqrt(s);
    for(i=0;i<n;i++)
    {
      MT_MATEL(XS,dim,i,j) = (smp[i][j]-mean[j])/s;
    }
  }

  //
  // Bias correction factor
  //
  if(wgt)
  {
    for(i=0,sum=0;i<n;i++) sum+=wgt[i]*wgt[i];
    h1 = 1/(1-sum);
  }
  else
  {
    h1 = (double)n/(double)(n-1);
  }
	
rcCatch:
  FREEMEM(cov0);
  FREEMEM(XS);
  return rcThis;
}
#endif

#if 1
/********************************************************************/
float MtClipInterval( float x, float a, float b, float c,
    		      unsigned opt )
/********************************************************************
 * map interval [0,c) to [0,b] with (smooth) cutoff at x=a 
*********************************************************************/			 
{

  if(opt & MT_LINEAR)
  {
    float y,d,t;

    d = (c-a)/(b-a);
    t = (b-a)/pow(fabs(a-c),d);

    if(x>=c)
    {
      y=b;
    }
    else if(x>a)
    {
      y = a + (x-a)*(b-a)/(c-a);
    }
    else y = x;

    return y;
  }
  else
  {
    double y,A,B,d,e;

    e = c - a;
    d = b - a;
    B = e / d;
    A = d / pow(e,B);

    if(x>=c)
    {
      y=b;
    }
    else if(x>a)
    {
/*      y = b-f(c-x);*/
      y = b-A*pow(c-x,B);
    }
    else y = x;

    return y;
  }

#if 0  
  float y,d,t;
  
  d = (c-a)/(b-a);
  t = (b-a)/pow(fabs(a-c),d);

  if(x>=c)
  {
    y=b;
  }
  else if(x>a)
  {
    y = b-t*pow(fabs(x-c),d);
  }
  else y = x;

  return y;
#endif

#if 0
  float y,t= 1/(b/a-1);
  if(x<a)
  {
    y = x;
  }
  else
  {
    y = a*(1.+1./t-1./(t*pow(x/a,t)));
  }
  return y;
#endif
}
#endif

/********************************************************************/
float MtSmoothClip( float x, float a, float b, float c )
/********************************************************************
 * map interval [0,c) to [0,b] with smooth cutoff at x=a 
*********************************************************************/			 
{
  float y,d,t;
  
  d = (c-a)/(b-a);
  t = (b-a)/pow(fabs(a-c),d);

  if(x>=c)
  {
    y=b;
  }
  else if(x>a)
  {
    y = b-t*pow(fabs(x-c),d);
  }
  else y = x;

  return y;

#if 0
  float y,t= 1/(b/a-1);
  if(x<a)
  {
    y = x;
  }
  else
  {
    y = a*(1.+1./t-1./(t*pow(x/a,t)));
  }
  return y;
#endif
}

/********************************************************************/
MtReal MtDstEstimateDensity( MtDst *dst, MtReal x, float h )
/********************************************************************			 
Estimate probability density of emprical distribution <dst> at <x>
*********************************************************************/			 
{
  int     i,n;
  MtReal  h2,y,sum;

  h2 = h > MT_REAL_EPSILON ? h : dst->stat.var/10;  
/*  h2 = .002;*/
  n=dst->nSmp;
  sum = 0;
  for(i=0;i<n;i++)
  {
    sum += MtDstKernelGauss( ( x - dst->vDist[i] ) / h2 );
  }
  y = sum/((MtReal)n*h);
  return y;
}

/********************************************************************/
MtReal MtVecParzenWinFuncDefault( MtVector x, int n, MtReal h )
/********************************************************************			 
*********************************************************************/			 
{
  MtReal d = MtVecAbs(x,n);
  return exp(-d*d / (2*h*h)) / MT_SQRT2PI;
}

/********************************************************************/
MtReal MtVecDensity( MtVector x, MtVector y, int d, int m,
    			     MtReal h, MtVecParzenWindowFunc winFunc )
/********************************************************************			 
  Estimate probability density at point x for ditribution y_i 
  y: array of samples 
  d: dimension of sample vectors
  m: number of samples
  h: window width
  winFunc: winodw function (NULL: use default)
*********************************************************************/			 
{
  int i;
  MtReal sum;
  MtVector vTmp = MtVecOpen(d);

  if(!winFunc) winFunc = MtVecParzenWinFuncDefault;

  for(i=0,sum=0;i<m;i++)
  {
    /* winFunc( x - y_i, h ) */
    sum += winFunc( MtVecSub( MtVecCopy_f( vTmp, x, d ), y, d ), d, h);
    y+=d;
  }

  MtVecClose(vTmp);
  return sum/(MtReal)m;
}

/********************************************************************/
MtReal MtVecDensityNN( MtVector x, MtVector y, int d, int m,
    			       int k )
/********************************************************************			 
  Estimate probability density at point x for ditribution y_i 
  y: array of samples 
  d: dimension of sample vectors
  m: number of samples
*********************************************************************/			 
{
  int i;
  MtVector vDist = MtVecOpen(m);

  for(i=0;i<m;i++)
  {
    vDist[i] = MtVecDistance(x,y,d); 
    y+=d;
  }

  MtVecSort( vDist, m );

  MtVecClose(vDist);
  return 1/vDist[0];
}



/********************************************************************/
/*								    */
/*  Stat: basic statistics                                          */
/*								    */
/********************************************************************/			 

void MtDstGetStat2_old( MtDst *dst, 
    		    float *skew, float *curt,
		    float *energy, float *entropy )
{
  double skew_d,curt_d,energy_d,entropy_d;
 MtDstGetStat2( dst, NULL,
    		    &skew_d, &curt_d,
		    &energy_d, &entropy_d );
 *skew = skew_d;
 *curt = curt_d;
 *energy = energy_d;
 *entropy = entropy_d;
}


/********************************************************************/
/*								    */
/********************************************************************/			 
void MtDstGetStat2( MtDst *dst, double *median,
    		    double *skew, double *curt,
		    double *energy, double *entropy )
{
  int 	i, n=dst->nSmp,me,rcThis=0;
  float	*x=dst->vDist;
  double f,y,yOld=0,xAvg = dst->stat.avg;
  register double a,b,s2,s3,s4,c,w;
  double nf=1./(double)n,dxe; 

  for(i=0,s2=s3=s4=0;i<n;i++)
  {
    a = b = *x++ - xAvg;
    a *=b;
    s2 += a;	/* += a^2 */
    a *= b;
    s3 += a;	/* += a^3 */
    a*=b;
    s4 += a;	/* += a^4 */
  }
  s2 *= nf;
  s3 *= nf;
  s4 *= nf;

  c = s2*s2*s2;
  b = c > MT_REAL_EPSILON ? sqrt( c ) : 0;
  if(skew) *skew = b > MT_REAL_EPSILON ? s3 / b : 0;
  b = (s2*s2);
  if(curt) *curt = b > MT_REAL_EPSILON ? s4 / b - 3 : 0;

  if(median||energy||entropy)
  {
    if(!dst->sorted)
    {
      MtVecSort( dst->vDist, dst->nSmp );
      dst->sorted = TRUE;
    }
  }

  if(median) 
    *median = MtVecMedian(dst->vDist,dst->nSmp,TRUE);

  if(energy || entropy )
  {
    UT_ASSERT(dst->sorted);
    w = dst->stat.max - dst->stat.min;
    if( w < MT_REAL_EPSILON ) w = MT_REAL_EPSILON;
    dxe = w/15.; /* TODO */
    for(f=dst->stat.min,s2=0,s3=0,me=0,i=0;
	f<dst->stat.max;
	f+=dxe,i++)
    {
      y = MtDstDistr(dst,f,1);
      if(i)
      {
	a = (y - yOld)*dxe/w;
	s2 += a*a;
	s3 += a*log(a+MT_NUM_EPS);
	me++;
      }
      yOld = y;
    }

    if(energy) *energy = s2;
    if(entropy) *entropy = -s3;
  }
rcCatch:
  return;
}

/********************************************************************/
/*								    */
/*  Correlation: linear regression x = alpha1 + alpha2*y            */
/*								    */
/********************************************************************/			 

/********************************************************************/
/*								    */
/********************************************************************/			 

int MtCorrelation_f( MtReal *xSeq, MtReal *ySeq, int n, MtReal *corr, 
                   MtReal *alpha1, MtReal *alpha2 )
{
  int         i;
  MtStatCorr  statCorr;          
  
  MT_STAT_CORR_INIT( &statCorr );
  for( i = 0; i < n; i++ )
  {
    MT_STAT_CORR_UPDATE( &statCorr, xSeq[i], ySeq[i] );
  }    
  MT_STAT_CORR_RESULT( &statCorr );
  *corr = statCorr.corr;
  *alpha1 = statCorr.alpha1; 
  *alpha2 = statCorr.alpha2;   
  return SUCCESS;
}
/********************************************************************/
/*								    */
/********************************************************************/			 
int MtCorrelation( double *xSeq, double *ySeq, int n, double *corr, 
                   double *alpha1, double *alpha2 )
{
  int         i;
  MtStatCorr  statCorr;          
  
  MT_STAT_CORR_INIT( &statCorr );
  for( i = 0; i < n; i++ )
  {
    MT_STAT_CORR_UPDATE( &statCorr, xSeq[i], ySeq[i] );
  }    
  MT_STAT_CORR_RESULT( &statCorr );
  if(corr) *corr = statCorr.corr;
  if(alpha1) *alpha1 = statCorr.alpha1; 
  if(alpha2) *alpha2 = statCorr.alpha2;   
  return SUCCESS;
}

/********************************************************************/			 
MtReal MtCorrConfidence_f( MtReal corr, int n, MtReal alpha )
/********************************************************************
Calc confidence of correlation coefficient estimation.
corr:   empiric correl.
n:      sample size
alpha:  confidence level
return: t>1 => significant, t<1: not significant
*********************************************************************/
{
  int    m;
  MtReal t1,t2;
  
  m = n - 2; 
  t1 = (corr*sqrt(m))/(sqrt(1-corr*corr));
  t2 = MtStudentT(alpha,m);
  /* t2 = 1.96; */
  return fabs(t1)/t2;
}
/********************************************************************/			 
double MtCorrConfidence( double corr, int n, double alpha )
/********************************************************************
Calc confidence of correlation coefficient estimation.
corr:   empiric correl.
n:      sample size
alpha:  confidence level
return: t>1 => significant, t<1: not significant
*********************************************************************/
{
  int    m;
  double t1,t2;
  
  m = n - 2; 
  t1 = (corr*sqrt(m))/(sqrt(1-corr*corr));
  t2 = MtStudentT(alpha,m);
  /* t2 = 1.96; */
  return fabs(t1)/t2;
}

/********************************************************************/			 
int MtMultipleCorrelation( MtVector *pvx, MtVector y, int n, int m, 
	                   MtVector b, MtReal *rss, MtReal *r,
		           MtVector t )
/********************************************************************
Multiple Regression analysis
Solve regression for y = b[0] + b[1..m]'x[1..m]  
  pvx  array of m pointers to the input sequences of regr. variables
  y    n output samples (n-vector)
  n    # of samples 
  m    # of input (regression) variables
  b    OUT estimated regression coefficients (m+1)-vector      
  rss  OUT regression residual
  r    OUT overall fitness (empiric multiple correlation)
  t    OUT t-statistic confidence of coefficients b (m+1)-vector    
*********************************************************************/
{
  MtMatrix  g, x;
  MtVector  c, d, h;
  double    yAvg, ySigma, tmp;
  MtReal   *px;
  int       i,j,p;
  
  p = m + 1;
  x = MtMatOpen( n, p );
  g = MtMatOpen( p, p );
  c = MtVecOpen( p );
  d = MtVecOpen( p );
  h = MtVecOpen( p );
  
  px = x;
  for(i=0;i<n;i++)
  {
    *px++ = 1;
    for(j=0;j<m;j++) *px++ = pvx[j][i]; 
  }
  
  /* MtMatPrint( x, n, p ); printf("\n");
  MtMatPrint( y, n, 1 ); printf("\n"); */
  
  MtMatMult( c, x, y, -p, n, 1 );      /* c = x'y */     
  MtMatMult( g, x, x, -p, n, p );      /* g = x'x */
  
  /* MtMatPrint( c, p, 1 ); printf("\n");
  MtMatPrint( g, p, p ); printf("\n");  */
  
  MtMatCholesky_f( g, p, d );            /* g = l'l (cholesky) */
  MtMatCholeskySolve( g, p, d, c, b);  /* solve g*b = c for b */
  
  for(i=0;i<p;i++)                     /* g = l ( clear upper triangle ) */
  {
    MT_MATEL(g,p,i,i) = d[i];
    for(j=i+1;j<p;j++) MT_MATEL(g,p,i,j) = 0; 
  }
  MtMatMult( h, g, b, -p, p, 1 );                  /* h = l'b */
  *rss = MtVecDotprod(y,y,n)-MtVecDotprod(h,h,p);  /* rss = y'y - h'h */ 
  yAvg = 0; for(i=0;i<n;i++) yAvg += y[i]; yAvg /= n;
  ySigma = 0; 
  for(i=0;i<n;i++) 
  {
    tmp = ( y[i] - yAvg );
    ySigma += tmp*tmp;
  }  
  *r = sqrt( 1 - *rss / ySigma );  
  
  MtVecClose( h );
  MtVecClose( d );
  MtVecClose( c );
  MtMatClose( g );  
  return SUCCESS;
}

/********************************************************************/			 
int MtMultipleCorrelationTest( int n, int m, MtReal noise, 
                               int rseed )
/********************************************************************
*********************************************************************/
{
  int 	      i, j, p = m+1;
  MtReal      sum, r, rss, tmp;
  MtMatrix    x = MtMatOpen( m, n );
  MtVector    y = MtVecOpen( n ),
              b0 = MtVecOpen( p ),
              b = MtVecOpen( p ),
              t = MtVecOpen( p ),
             *pvx;
               
  MtRndSeed(rseed);
  printf("b0 = ");
  for(i=0;i<p;i++) 
  {
    b0[i] = 10*MtRndGet()-5;
    printf("%6g ", b0[i] );
  }  
  printf("\n");
    
  for(i=0;i<n;i++)
  {
    /* printf("%6d: ",i); */
    for(j=0;j<m;j++) 
    {
      tmp = 10*MtRndGet()-5; 
      MT_MATEL(x,n,j,i) = tmp;
      /* printf("%6g ",tmp); */
    }    
    /* printf(" -> "); */
    sum = b0[0];
    for(j=0;j<m;j++) sum += b0[j+1] * MT_MATEL(x,n,j,i);
    y[i] = sum;
    y[i] = y[i] + MtRndNormGet( 0, noise );
    /* printf(" %6g \n",y[i]);  */
  } 
  
  pvx = MM_malloc( m * sizeof( MtVector *), MM_DUMMYTAG);
  for(i=0;i<m;i++) pvx[i] = &MT_MATEL(x,n,i,0);
   
  MtMultipleCorrelation( pvx, y, n, m, b, &rss, &r, t );
  
  printf("b  = ");
  for(i=0;i<p;i++) printf("%6g ",b[i]);
  printf("\n");
  printf("rss=%6g\n",rss);
  printf("r=%6g\n",r);
  
  MM_free(pvx);
  MtVecClose(t);
  MtVecClose(b);
  MtVecClose(b0);
  MtVecClose(y);
  MtMatClose(x);

  return SUCCESS;
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/ 
int MtGetStat2_f( LongInt n,
	       float *data,
	       float validMin, float validMax,
	       MtStat *stat
		 )
{
  int rcThis=0;
  float zmax,zmin;
  double zsum,zsum2,nsmp;

  zmax=-1e20; zmin=1e20; zsum=0; zsum2=0; nsmp=0;
  #pragma omp parallel
  {
    float my_zmax=-1e20, my_zmin=1e20;
    double my_zsum=0, my_zsum2=0, my_nsmp=0;
    LongInt ii;
    #pragma omp for private(ii) 
    for(ii=0;ii<n;ii++)
    {
      float z = data[ii];
#if 0
      if(!isfinite(z))
      {
	TRCERR(("invalid floating point number %g at pos %d\n",z,ii));
      }
#endif
      if(!(z>validMin&&z<validMax)) continue;
      my_nsmp += 1;
      my_zsum += z;
      my_zsum2 += z*(double)z;
      if(z>my_zmax) my_zmax=z;
      if(z<my_zmin) my_zmin=z;
    }
    #pragma omp critical
    {
      nsmp += my_nsmp;
      zmax=MAX(zmax,my_zmax);
      zmin=MIN(zmin,my_zmin);
      zsum+=my_zsum;
      zsum2+=my_zsum2;
    }
  }
  MT_STAT_INIT(stat);
  MT_STAT_UPDATE_N(stat,nsmp,zmin,zmax,zsum,zsum2);
  MT_STAT_RESULT(stat);

  if(!isfinite(zsum))
  {
    TRCERRR(("invalid floating point number %g\n",zsum),
	MI_ENUMERICAL);
  }
  if(!isfinite(zsum2))
  {
    TRCERRR(("invalid floating point number %g\n",zsum2),
	MI_ENUMERICAL);
  }
 
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/ 
int MtGetStat_f( LongInt n,
	       float *data,
	       MtStat *stat
		 )
{
  int rcThis=0;
  float zmax,zmin;
  double zsum,zsum2;

  zmax=-1e20; zmin=1e20; zsum=0; zsum2=0;
  #pragma omp parallel
  {
    float my_zmax=-1e20, my_zmin=1e20;
    double my_zsum=0, my_zsum2=0;
    LongInt ii;
    #pragma omp for private(ii) 
    for(ii=0;ii<n;ii++)
    {
      float z = data[ii];
#if 0
      if(!isfinite(z))
      {
	TRCERR(("invalid floating point number %g at pos %d\n",z,ii));
      }
#endif
      my_zsum += z;
      my_zsum2 += z*(double)z;
      if(z>my_zmax) my_zmax=z;
      if(z<my_zmin) my_zmin=z;
    }
    #pragma omp critical
    {
      zmax=MAX(zmax,my_zmax);
      zmin=MIN(zmin,my_zmin);
      zsum+=my_zsum;
      zsum2+=my_zsum2;
    }
  }
  MT_STAT_INIT(stat);
  MT_STAT_UPDATE_N(stat,n,zmin,zmax,zsum,zsum2);
  MT_STAT_RESULT(stat);

  if(!isfinite(zsum))
  {
    TRCERRR(("invalid floating point number %g\n",zsum),
	MI_ENUMERICAL);
  }
  if(!isfinite(zsum2))
  {
    TRCERRR(("invalid floating point number %g\n",zsum2),
	MI_ENUMERICAL);
  }
 
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/ 
int MtGetStat( LongInt n,
	       double *data,
	       MtStat *stat
		 )
{
  int rcThis=0;
  double zmax,zmin;
  double zsum,zsum2;

  zmax=-1e20; zmin=1e20; zsum=0; zsum2=0;
  #pragma omp parallel
  {
    double my_zmax=-1e20, my_zmin=1e20;
    double my_zsum=0, my_zsum2=0;
    LongInt ii;
    #pragma omp for private(ii) 
    for(ii=0;ii<n;ii++)
    {
      double z = data[ii];
      my_zsum += z;
      my_zsum2 += z*z;
      if(z>my_zmax) my_zmax=z;
      if(z<my_zmin) my_zmin=z;
    }
    #pragma omp critical
    {
      zmax=MAX(zmax,my_zmax);
      zmin=MIN(zmin,my_zmin);
      zsum+=my_zsum;
      zsum2+=my_zsum2;
    }
  }
  MT_STAT_INIT(stat);
  MT_STAT_UPDATE_N(stat,n,zmin,zmax,zsum,zsum2);
  MT_STAT_RESULT(stat);

  if(!isfinite(zsum))
  {
    TRCERRR(("invalid floating point number %g\n",zsum),
	MI_ENUMERICAL);
  }
  if(!isfinite(zsum2))
  {
    TRCERRR(("invalid floating point number %g\n",zsum2),
	MI_ENUMERICAL);
  }
 
rcCatch:
  return rcThis;
}

/********************************************************************/
/*								    */
/********************************************************************/			 

int MtStatGet_f( MtStat *stat, int nSmp, MtReal *samples )
{
  int i; 
  
  MT_STAT_INIT( stat );
  for( i=0; i < nSmp; i++ )
    MT_STAT_UPDATE( stat, samples[i] );
  MT_STAT_RESULT( stat );  
  return SUCCESS;
}

int MtStatGet( MtStat *stat, int nSmp, double *samples )
{
  int i; 
  
  MT_STAT_INIT( stat );
  for( i=0; i < nSmp; i++ )
    MT_STAT_UPDATE( stat, samples[i] );
  MT_STAT_RESULT( stat );  
  return SUCCESS;
}

void MtVecStat_d( MtStat *s, double *x, LongInt n )
{
  LongInt i; 
  
  MT_STAT_INIT( s );
  for( i=0; i < n; i++ ) MT_STAT_UPDATE( s, x[i] );
  MT_STAT_RESULT( s );  
}

void MtVecStat( MtStat *s, float *x, LongInt n )
{
  LongInt i; 
  
  MT_STAT_INIT( s );
  for( i=0; i < n; i++ ) MT_STAT_UPDATE( s, x[i] );
  MT_STAT_RESULT( s );  
}

/*-------------------------------------------------------------------------*/
MtStat *MtMultStatGet1( MtTab *X )
/*-------------------------------------------------------------------------
Get multivariate statistics
IN: Matrix X (nRow=number of samples, nCol=dimension)
---------------------------------------------------------------------------*/
{
  int i,k,dim=X->nCol;
  MtStat *S=NULL;

  ALLOCARR0(S,dim);
  for(k=0;k<dim;k++) MT_STAT_INIT(&S[k]);

  for(i=0;i<X->nRow;i++)
    for(k=0;k<dim;k++)
      MT_STAT_UPDATE(&S[k],X->rows[i][k]);

  for(k=0;k<dim;k++) MT_STAT_RESULT(&S[k]);

  return S;      
}

/********************************************************************/
/*								    */
/********************************************************************/			 
void MtStatPrint( MtStat *stat )
{
  printf("n=%d> mean=%g var=%g min=%g max=%g\n",(int)stat->n,
      (double)stat->avg, (double)stat->var, (double)stat->min,
      (double)stat->max);
}


/********************************************************************/
/*								    */
/********************************************************************/			 

MtReal MtVecSigmaLeft( MtVector x, int n )
{
  int    i,k;
  MtReal sigma2;
  
  MtStat stat;
  MtStatGet_f( &stat, n, x );
  k = 0;
  sigma2 = 0;
  for(i=0;i<n;i++)
  { 
    if( x[i] < stat.avg )
    { 
      sigma2 += (x[i] - stat.avg) * (x[i] - stat.avg);
      k++;
    }  
  }   
  sigma2 = sqrt( sigma2 / k ); 
  return sigma2;
}

/********************************************************************/
/*								    */
/*  Nonparametric nonlinear multiple regression and forecasting     */
/*								    */
/********************************************************************/

/**********************************************************************/
static int MtNLRegNeighborCompareFunc( const void *nb1, 
                                       const void *nb2 )
/**********************************************************************
comparison function for qsort()
***********************************************************************/
{
if(((  MtNLRegNeighbor *)nb1)->distance >= 
    (( MtNLRegNeighbor *)nb2)->distance )
  return  1;
else
  return -1;
}
			 
/********************************************************************/
int MtNLRegOpen( int nSmp, MtNLRegHandle *handle )
/********************************************************************
*********************************************************************/                        
{
  int    rc;
  struct MtNLRegAdm *adm=NULL;
  
  adm = MM_malloc( sizeof (struct MtNLRegAdm), MM_DUMMYTAG );
  if(!adm) 
  {
    rc = MT_ENOMEM;
    goto fail;
  }  
  
  adm->neighbors = NULL;
  adm->samplesIn = NULL;
  adm->samplesOut = NULL;
  
  *handle = adm;
  return MT_OK;
fail:
  MtNLRegClose(adm); 
  return rc;
}

/********************************************************************/
int MtNLRegFit( MtNLRegHandle handle, MtVector *pvx, MtVector y,
                int n, int m )
/********************************************************************
Fit a nonlinear nonparametric model to the function 
  y = F(x_1,..,x_m)
given the n observed samples 
  pvx[0][0]   ... pvx[0][n-1]
  ...             ...
  pvx[m-1][0] ... pvx[m-1][n-1]     
Parameters:    
  pvx  array of m pointers to the input sequences of regr. variables
  y    n output samples (n-vector)
  n    # of samples 
  m    # of input (regression) variables
*********************************************************************/                        
{
  struct MtNLRegAdm *adm = handle;
  
  adm->samplesIn = pvx;
  adm->samplesOut = y;
  adm->n = n;
  adm->m = m; 
  return MT_OK;
}

/********************************************************************/
int MtNLRegForecast( MtNLRegHandle handle, MtVector x, MtReal *y )
/********************************************************************
Forecast the value y = F( x_1,...,x_m ) using the nonlinear model
created by MtNLRegFit()
Parameters:
  handle  MtNLReg object with MtNLRegFit() called previously
  x       input variable
  y       (OUT) estimated (forecasted) value 
*********************************************************************/                        
{
  int                 n,m,nN,i,j;
  struct MtNLRegAdm  *adm = handle;
  MtNLRegNeighbor    *neighbors;
  MtVector            vy,x2,vDistance;
  MtStat              stat,statDistance;

  x2 = MtVecOpen( adm->m );
  vDistance = MtVecOpen( adm->n );
  neighbors = MM_malloc( adm->n * sizeof( MtNLRegNeighbor ), MM_DUMMYTAG );
  if( ! neighbors )
  {
    printf("MtNLRegForecast: MM_malloc() failed.\n");   
  }
  adm->neighbors = neighbors;
  n = adm->n;
  m = adm->m;
  
  /**** compute distance from each sample point */
    
  for(i=0;i<n;i++)
  {
    neighbors[i].idx = i;
    for(j=0;j<m;j++) x2[j] = adm->samplesIn[j][i];
    neighbors[i].distance = MtVecDistance( x, x2, m );
    vDistance[i] = neighbors[i].distance;       
  }
  MtStatGet_f( &statDistance, n, vDistance );
  
  /**** sort all points in increasing distance */
  
  qsort( neighbors, 
         n, 
         sizeof( MtNLRegNeighbor ), 
         MtNLRegNeighborCompareFunc );  

  /**** compute x's y-value according to neighbor's y-values */  

/*  nN = MAX( 6, n/100 );*/
  nN = 20;
  vy = MtVecOpen( nN );
  for(i=0;i<nN;i++) vy[i] = adm->samplesOut[neighbors[i].idx]; 
  MtStatGet_f( &stat, nN, vy );

  *y = stat.avg;

/*  *y = MtVecMedian( vy, nN );*/

  MtVecClose( vy );    
  MtVecClose( x2 );
  MtVecClose( vDistance );
  return MT_OK;
}

/********************************************************************/
int MtNLRegClose( MtNLRegHandle handle )
/********************************************************************
*********************************************************************/                        
{
  struct MtNLRegAdm  *adm = handle;

  if( adm )
  {
    if( adm->neighbors ) MM_free( adm->neighbors );
    MM_free( adm );
  }  
  return MT_OK;
}

/********************************************************************/
/*								    */
/*  MA: simple moving average                                       */
/*								    */
/********************************************************************/			 

/********************************************************************/
/*								    */
/********************************************************************/			 

int MtMAInit( MtMA *ma, int max, double y0 )
{
  int i; 

  if( max > MT_MA_MAX )
  {
    printf(" MtMAInit(): parameter 'max'=%d out of range\n",
         max );
    return 1;
  }
  ma->max = max; 
  ma->sum = 0; 
  ma->beta = 1;
  for( i = 0; i < max; i++ )
  {     
    ma->buffer[ i ] = y0; 
    ma->sum += y0;        
  } 
  ma->val = y0;
  ma->in = ma->out = 0;
  return SUCCESS;
}

/********************************************************************/
/*								    */
/********************************************************************/			 
double MtMAUpdate( MtMA *ma, double y )
{
  ma->sum -= ma->buffer[ ma->out++ ];     /* dequeue oldest entry */
  if( ma->out >= ma->max ) ma->out = 0; 
  ma->buffer[ ma->in++ ] = y;
  if( ma->in >= ma->max ) ma->in = 0;
  ma->sum += y;
  return ( ma->val = ma->sum / ma->max );  
}

/********************************************************************/
/*								    */
/* MAV: moving average + volatility                                 */
/*								    */
/********************************************************************/			 

/********************************************************************/
/*								    */
/********************************************************************/			 

int MtMAVInit( MtMAV *mv, int max, double y0 )
{
  int 	   i; 
  double   tmp;

  if( max > MT_MA_MAX )
  {
    printf(" MtMVInit(): parameter 'max'=%d out of range\n",
         max );
    return 1;
  }
  mv->max = max; 
  mv->sum = 0;
  mv->sumQ = 0; 
  mv->beta = 1;
  tmp = y0*y0;
  for( i = 0; i < max; i++ )
  {     
    mv->buffer[ i ] = y0; 
    mv->sum += y0;
    mv->sumQ += tmp;        
  } 
  mv->in = mv->out = 0;
  return SUCCESS;
}

/********************************************************************/
/*								    */
/********************************************************************/			 
void MtMAVUpdate( MtMAV *mv, double y )
{
  double  tmp;
  
  tmp = mv->buffer[ mv->out++ ];     /* dequeue oldest entry */
  mv->sum -= tmp;
  mv->sumQ -= tmp*tmp; 
  if( mv->out >= mv->max ) mv->out = 0; 
  mv->buffer[ mv->in++ ] = y;
  if( mv->in >= mv->max ) mv->in = 0;
  mv->sum += y;
  mv->sumQ += y*y;
}

/********************************************************************/
/*								    */
/********************************************************************/			 
void MtMAVGet( MtMAV *ma, double *mean, double *variance )
{
  *mean = ma->sum / ma->max;
  *variance = ma->sumQ / ma->max - (*mean)*(*mean);   
}  

/********************************************************************/
/*								    */
/* MAE: exponential moving average                                  */
/*								    */
/********************************************************************/			 

/********************************************************************/
/*								    */
/********************************************************************/			 
int MtInitMAE( MtMA *ma, int max, double beta, double y0 )
{
  int i; 

  if( max > MT_MA_MAX )
  {
    printf(" MtMAInit(): parameter 'max'=%d out of range\n",
         max );
    return 1;
  }
  ma->max = max; 
  ma->sum = 0;
  ma->beta = beta; 
 
  for( i = 0; i < max; i++ )
  {     
    ma->buffer[ i ] = y0; 
    ma->sum += pow( ma->beta, max-i-1 ) * y0;        
  } 
  ma->in = ma->out = 0;
  return SUCCESS;
}

/********************************************************************/
/*								    */
/********************************************************************/			 
double MtUpdateMAE( MtMA *ma, double y )
{
  double yOld, sumBeta;
  
  yOld = ma->buffer[ ma->out++ ];
  if( ma->out >= ma->max ) ma->out = 0; 
  
  ma->sum = y + ma->beta * ma->sum - pow( ma->beta, ma->max ) * yOld;
       
  ma->buffer[ ma->in++ ] = y;
  if( ma->in >= ma->max ) ma->in = 0;
  
  sumBeta = ( pow( ma->beta, ma->max ) - 1 ) / ( ma->beta -1 );
  
  return ma->sum / sumBeta ;  
}

/********************************************************************/
/*								    */
/********************************************************************/			 

int MtIsPow2 ( unsigned x )
{
   unsigned i, y;

   for ( i=1, y=2; i < BITS_PER_WORD; i++, y<<=1 )
   {
      if ( x == y ) return TRUE;
   }

   return FALSE;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtFFT(int dir,int log2n,double *x,double *y)
{
   long nn,i,i1,j,k,i2,l,l1,l2,m=log2n;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   nn = 1;
   for (i=0;i<m;i++)
      nn *= 2;

   /* Do the bit reversal */
   i2 = nn >> 1;
   j = 0;
   for (i=0;i<nn-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0;
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0;
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<nn;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1;
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1)
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<nn;i++) {
         x[i] /= (double)nn;
         y[i] /= (double)nn;
      }
   }

   return(TRUE);
}


/*-------------------------------------------------------------------------
   Perform a 2D FFT inplace given a complex 2D array
   The direction dir, 1 for forward, -1 for reverse
   The size of the array (nx,ny)
   Return false if there are memory problems or
      the dimensions are not powers of 2
*/
int MtFFT2D(double *CR, double *CI, int nx,int ny,int dir)
{
#define FFT2D_MAXDIM 32768

   int i,j,rcThis=0;
   int nx2,ny2,kx,ky;
   double real[FFT2D_MAXDIM],imag[FFT2D_MAXDIM];

  UT_ASSERT(nx<=FFT2D_MAXDIM && ny<=FFT2D_MAXDIM);

  MT_NEAREST_POW2( nx, nx2, kx, 1<<20 );  
  MT_NEAREST_POW2( ny, ny2, ky, 1<<20 );  

  if(nx!=nx2 || ny!=ny2)
  {
    TRCERRR(("MtFFT2D: input dimensions not a power of two (%dx%d)\n",
	  nx,ny),1);
  }

   /* Transform the rows */
   for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
         real[i] = CR[MT_GXY(nx,i,j)];
         imag[i] = CI[MT_GXY(nx,i,j)];
      }
      MtFFT(dir,kx,real,imag);
      for (i=0;i<nx;i++) {
         CR[MT_GXY(nx,i,j)] = real[i];
         CI[MT_GXY(nx,i,j)] = imag[i];
      }
   }

   /* Transform the columns */
   for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
         real[j] = CR[MT_GXY(nx,i,j)];
         imag[j] = CI[MT_GXY(nx,i,j)];
      }
      MtFFT(dir,ky,real,imag);
      for (j=0;j<ny;j++) {
         CR[MT_GXY(nx,i,j)] = real[j];
         CI[MT_GXY(nx,i,j)] = imag[j];
      }
   }
rcCatch:
   return(rcThis);
}

/*-------------------------------------------------------------------------
   Perform a 2D FFT inplace given a complex 2D array
   The direction dir, 1 for forward, -1 for reverse
   The size of the array (nx,ny)
   Return false if there are memory problems or
      the dimensions are not powers of 2
*/
int MtFFT2D_f(float *CR, float *CI, int nx,int ny,int dir)
{
//#define FFT2D_MAXDIM 32768

   int i,j,rcThis=0;
   int nx2,ny2,kx,ky;
   double *real=NULL,*imag=NULL;

//  UT_ASSERT(nx<=FFT2D_MAXDIM && ny<=FFT2D_MAXDIM);

  MT_NEAREST_POW2( nx, nx2, kx, 1<<20 );  
  MT_NEAREST_POW2( ny, ny2, ky, 1<<20 );  

  if(nx!=nx2 || ny!=ny2)
  {
    TRCERRR(("MtFFT2D: input dimensions not a power of two (%dx%d)\n",
	  nx,ny),1);
  }

  ALLOCARR(real,MAX(nx,ny));
  ALLOCARR(imag,MAX(nx,ny));

   /* Transform the rows */
   for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
         real[i] = CR[MT_GXY(nx,i,j)];
         imag[i] = CI[MT_GXY(nx,i,j)];
      }
      MtFFT(dir,kx,real,imag);
      for (i=0;i<nx;i++) {
         CR[MT_GXY(nx,i,j)] = real[i];
         CI[MT_GXY(nx,i,j)] = imag[i];
      }
   }

   /* Transform the columns */
   for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
         real[j] = CR[MT_GXY(nx,i,j)];
         imag[j] = CI[MT_GXY(nx,i,j)];
      }
      MtFFT(dir,ky,real,imag);
      for (j=0;j<ny;j++) {
         CR[MT_GXY(nx,i,j)] = real[j];
         CI[MT_GXY(nx,i,j)] = imag[j];
      }
   }
rcCatch:
   FREEMEM(real);FREEMEM(imag);
   return(rcThis);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtFFT2D_old( double *X,     // input data matrix 
    	     float  *X_f,   //     (optional: single precisison input)
	     double *XI,    //  optional: imaginarey part 
	     float  *XI_f,
    	     int n, int m,  // dimensions (rows,columns) 
	     		    // ATTENTION: Power of two required !
    	     int dir,	    // direction: 1=forward,-1=inverse transform
	     double *FR,    // OUT: real part
	     double *FI,    // OUT: imaginary part
	     double *F	    // OUT: sqrt(FR^2+FI^2)
    	)
{
  int kx,ky,nx,ny,rcThis=0,fftDir=dir,i,j,nx2,ny2,i2,j2;
  double *real=FR,*imag=FI,f,
	 realBuf[FFT2D_MAXDIM],
	 imagBuf[FFT2D_MAXDIM];

  UT_ASSERT(n<=FFT2D_MAXDIM && m<=FFT2D_MAXDIM);

  MT_NEAREST_POW2( m, nx, kx, 1<<20 );  
  MT_NEAREST_POW2( n, ny, ky, 1<<20 );  

  if(m!=nx || n!=ny)
  {
    TRCERRR(("MtFFT2D: input dimensions not a power of two (%dx%d)\n",
	  m,n),1);
  }
  /*** transform rows */
  for (i=0;i<ny;i++) 
  {
    for (j=0;j<nx;j++) 
    {
      MT_MATEL(real,nx,i,j) = X_f ? MT_MATEL(X_f,m,j,i) : 
				    MT_MATEL(X,m,j,i);
      if(XI)
	f = MT_MATEL(XI,m,j,i) ; 
      else if(XI_f)
        f = MT_MATEL(XI_f,m,j,i) ; 
      else f = 0;

      MT_MATEL(imag,nx,i,j) = f;
    }
    MtFFT( fftDir, kx, real+nx*i, imag+nx*i ); 
  }

  /*** transform columns */
  nx2=nx/2; ny2=ny/2;
  for (j=0;j<nx;j++) 
  {
    for (i=0;i<ny;i++) 
    {
      realBuf[i] = MT_MATEL( real,nx,i,j );
      imagBuf[i] = MT_MATEL( imag,nx,i,j );
    }
    MtFFT( fftDir, ky, realBuf, imagBuf ); 
    if(F)
    {
      for(i=0;i<ny;i++)
      {
	if(i<ny2)
	{
	  i2 = i+ny2;
	  j2 = j<nx2 ? j+nx2 : j-nx2; 
	}
	else
	{
	  i2=i-ny2;
	  j2 = j<nx2 ? j+nx2 : j-nx2; 
	}

	f = sqrt( realBuf[i]*realBuf[i]+ imagBuf[i]*imagBuf[i] );
	MT_MATEL(F,m,j2,i2) = f;
      }
    }
  }


rcCatch:
  return rcThis;
}

/********************************************************************/			 
static double MtIStudentProb( double x, double *params )
/********************************************************************
Student's t-probability ( utility func for MtStudentT())
 x: t
 params[0]: degrees of freedom;
 params[1]: alpha;
*********************************************************************/
{
  double prob;
  int    rc;
  
  rc = tprob( x, params[0], &prob );
  if(rc) printf("MtIStudentProb(): tprob() failed");
  prob = 2*(1-prob);
  return ( prob - params[1]);
}

/********************************************************************/			 
MtReal MtStudentT( MtReal alpha, int m )
/********************************************************************
Student's t-distribution
  alpha:  confidence level
  m:	  degrees of freedom      
*********************************************************************/
{
  int 	   rc;
  double   x0, x1, x, y0, y1, epsilon;
  double   parms[2];
  
  parms[0] = m;
  parms[1] = alpha;
  x0 = 0;
  x1 = 10;
  y0 = MtIStudentProb( x0, parms );
  while( (y1 = MtIStudentProb( x1, parms ))*y0 > 0 ) 
    x1 *= 2;  
  epsilon = alpha / 500;  
  rc = 
    MtSolveBin( MtIStudentProb, parms, x0, x1, 0.0001, 1000, &x, 
	NULL, NULL );
  return x;   
}

/********************************************************************/
/*								    */
/********************************************************************/			 
double MtFactorial( double n )
{
  int i; double f=1;
  for(i=1;i<=(int)n;i++) f *= i;
  return f;
}

/********************************************************************/
/*								    */
/********************************************************************/			 
double MtFactorialOdd( double n )
{
  int i; double f=1;
  for(i=1;i<=(int)n;i+=2) f *= i;
  return f;
}
/*************************************************************07.05.1998*/
/* quartic: solution of a biquadratic equation				*/
/* Equations of lesser degree are solved by the appropriate formulas.	*/
/* The solutions are arranged in ascending order.			*/
/*									*/
/* ref.: R. Zurmuehl, Praktische Mathematik fuer Ingenieure und		*/
/* Physiker, Springer, Berlin 1965, p. 60ff				*/

/* Program by Thomas Kraska
   T.Kraska@uni-koeln.de
   http://www.uni-koeln.de/math-nat-fak/phchem/deiters/quartic/quartic.html */

/************************************************************************/
/* a[5]		(i)  vector containing the polynomial coefficients
   x[]		(o)  result vector
   RESULT	(o)  number of valid solutions				*/
/*======================================================================*/
#if 0
#include <math.h>
#include <mathc.h>
#endif

extern double sqrt (double);
extern double pow (double, double);
extern double cos (double);
extern double acos (double);
double fabs (double);

int MtSolveCubic(double *, double *);

#define M_PI 3.14159265358979323846
#define SWAP2(X,Y)	{ h = X; X = Y; Y = h; }

int MtSolveQuartic(double *a, double *x)
{
  int i, nreal;
  double w, b0, b1, b2, c[4], d0, d1, h, t, z, *px;
  double w2;

  if (a[4]==0.0) return(MtSolveCubic(a,x));	/* quartic problem? */

  w = a[3]/(4*a[4]);			/* offset */
  w2 = w * w;
  b2 = -6*(w2) + a[2]/a[4];		/* coeffs. of shifted polynomial */
  b1 = (8*(w2) - 2*a[2]/a[4])*w + a[1]/a[4];
  b0 = ((-3*(w2) + a[2]/a[4])*w - a[1]/a[4])*w + a[0]/a[4];

  c[3] = 1.0;				/* cubic resolvent */
  c[2] = b2;
  c[1] = -4*b0;
  c[0] = (b1*b1) - 4*b0*b2;
  i = MtSolveCubic(c, x);
  z = x[0];				/* only lowermost root needed */

  nreal = 0;
  px = x;
  t = sqrt(0.25*(z*z) - b0);
  for (i=-1; i<=1; i+=2) {
    d0 = -0.5*z + i*t;			/* coeffs. of quadratic factor */
    d1 = (t!=0.0)? -i*0.5*b1/t : i*sqrt(-z - b2);
    h = 0.25*(d1*d1) - d0;
    if (h>=0.0) {
      h = sqrt(h);
      nreal += 2;
      *px++ = -0.5*d1 - h - w;
      *px++ = -0.5*d1 + h - w;
    }
  }

  if (nreal==4) {			/* sort results */
    if (x[2]<x[0]) SWAP2(x[0], x[2]);
    if (x[3]<x[1]) SWAP2(x[1], x[3]);
    if (x[1]<x[0]) SWAP2(x[0], x[1]);
    if (x[3]<x[2]) SWAP2(x[2], x[3]);
    if (x[2]<x[1]) SWAP2(x[1], x[2]);
  }
  return(nreal);
}

/*************************************************************05.05.1998*/
/* cubic: solution of a cubic equation					*/
/* Equations of lesser degree are solved by the appropriate formulas.	*/
/* The solutions are arranged in ascending order.			*/
/************************************************************************/
/* a[4]		(i)  vector containing the polynomial coefficients
   x[]		(o)  result vector
   RESULT	(o)  number of valid solutions				*/
/*======================================================================*/
#if 0
#include <math.h>
#include <mathc.h>
#endif

#define Q1_3		1.0/3.0
#define SWAP2(X,Y)	{ h = X; X = Y; Y = h; }

int MtSolveCubic(double *a, double *x)
{
  int i, nreal;
  double w, p, q, dis, h, phi;

  if (a[3]!=0.0) {			/* cubic problem? */
    w = a[2]/(3*a[3]);
    p = a[1]/(3*a[3])-(w*w);
    p = p*p*p;
    q = -0.5*(2*(w*w*w)-(a[1]*w-a[0])/a[3]);
    dis = (q*q)+p;			/* discriminant */
    if (dis<0.0) {			/* 3 real solutions */
      h = q/sqrt(-p);
      if (h>1.0) h = 1.0;		/* confine the argument of */
      if (h<-1.0) h = -1.0;		/* acos to [-1;+1] */
      phi = acos(h);
      p = 2*pow(-p, 1.0/6.0);
      for (i=0; i<3; i++) x[i] = p*cos((phi+2*i*M_PI)/3.0) - w;
      if (x[1]<x[0]) SWAP2(x[0], x[1]);	/* sort results */
      if (x[2]<x[1]) SWAP2(x[1], x[2]);
      if (x[1]<x[0]) SWAP2(x[0], x[1]);
      nreal = 3;
    }
    else {				/* only one real solution */
      dis = sqrt(dis);
      h = pow(fabs(q+dis), 1.0/3.0);
      p = pow(fabs(q-dis), 1.0/3.0);
      x[0] = ((q+dis>0.0)? h : -h) + ((q-dis>0.0)? p : -p) -  w;
      nreal = 1;
    }
    /* Perform one step of a Newton iteration in order to minimize
       round-off errors */
    for (i=0; i<nreal; i++, x++) {
      if ((h = a[1] + *x * (2 * a[2] + 3 * *x * a[3])) != 0.0)
      *x -= (a[0] + *x * (a[1] + *x * (a[2] + *x * a[3])))/h;
    }
  }

  else if (a[2]!=0.0) {			/* quadratic problem? */
    p = 0.5*a[1]/a[2];
    dis = (p*p) - a[0]/a[2];
    if (dis>=0.0) {			/* two real solutions */
      dis = sqrt(dis);
      x[0] = -p - dis;
      x[1] = -p + dis;
      nreal = 2;
    }
    else				/* no real solution */
      nreal = 0;
  }

  else if (a[1]!=0.0) {			/* linear problem? */ 
    x[0] = -a[0]/a[1];
    nreal = 1;
  }

  else					/* no equation */
    nreal = 0;

  return(nreal);
}

/********************************************************************/
/*								    */
/********************************************************************/			 
double MtSolveDiff( double    (* f)( double x, double *parms ),  
                double    (* df)( double x, double *parms ),
                double   *parms,
                double    x0,
                double    epsilon,
                int       maxIter )
{
  double x,y;
  int    i;
  char   msg[128];

  i = 0;
  x = x0;
  
  while( ( fabs( y = f( x, parms ) ) > epsilon ) && 
         ( i < maxIter ) )
  {       
    /* printf("MtSolve: %d x=%f f(x)=%f\n", i, x, y );  */
    x = x - y / df( x, parms );
    i++;    
  }
  if( i >= maxIter )
  {
    sprintf( msg, "MtSolve(): too many iterations. current error=%f",y);
    printf( "%s",msg );
  }

return x;
}                

/********************************************************************/
/*								    */
/********************************************************************/			 
double MtSolve2( double    (* f)( double x, double *parms ),  
                 double    *parms,
                 double    x0,
                 double    epsilon,
                 int       maxIter )
{
  double x,y,df,step;
  int    i;
  char   msg[128];

  i = 0;
  x = x0;
  
  step = epsilon;
  while( ( fabs( y = f( x, parms ) ) > epsilon ) && 
         ( i < maxIter ) )
  {       
    /* printf("MtSolve: %d x=%f f(x)=%f\n", i, x, y );  */
    df = ( f( x + step, parms ) - y ) / step;
    x = x - y / df;
    i++;    
  }

  if( i >= maxIter )
  {
    sprintf( msg, "MtSolve(): too many iterations. current error=%f",y);
    printf( "%s", msg );
  }

return x;
}                

/********************************************************************/
/*								    */
/********************************************************************/			 
int MtSolve( double    (* f)( double x, double *parms ),  
             double   *parms,
             double    x0,
             double    x1,
             double    xMin,
             double    xMax,
             double    epsilon,
             int       maxIter,
             double   *solution )
{
  double x,y0,y1,y;
  int    i;
  char   msg[128];

  i = 0;
  y0 = f( x0, parms );
  y1 = f( x1, parms );
  
  for(;;)
  {           
    x = x1 - y1 * ( x1 - x0 ) / ( y1 - y0 );
    if( x < xMin ) x = xMin;
    if( x > xMax ) x = xMax;
    y = f( x, parms );
    /* printf("MtSolve: %6d   %6g %6g   %6g %6g   %6g %6g\n", 
      i, x0, y0, x1, y1, x, y );     */
    if(( fabs( y ) < epsilon ) || ( i >= maxIter ))
      break;
    /* if( y*y1 < 0 ) */
    {
      x0 = x1;
      y0 = y1;
    }  
    x1 = x;
    y1 = y;
    i++;        
  }

  *solution = x;
  
  if( i >= maxIter )
  {
    sprintf( msg, "MtSolve(): too many iterations. current error=%f",y);
    printf( "%s", msg );
    return 1;
  }

return SUCCESS;
}                

/********************************************************************/
/*								    */
/********************************************************************/			 
int MtSolveBin( double    (* f)( double x, double *parms ),  
                double   *parms,
                double    x0,
                double    x1,
                double    epsilon,
                int       maxIter,
                double   *solution,
		double   *p_x0,
		double   *p_x1 )
{
  double x=0,y0,y1,y=0;
  int    i;
  //char   msg[128];

  i = 0;
  y0 = f( x0, parms );
  y1 = f( x1, parms );
  
  for( i=0; i<maxIter; i++)
  {
    x = x0 + ( x1 - x0 ) / 2;
    y = f( x, parms );
    if( fabs( y ) < epsilon ) break;    
    if( y*y0 < 0 )
    {
      x1 = x; y1 = y;
    }
    else
    {
      x0 = x; y0 = y;
    }  
  }

  *solution = x;

  if(p_x0) *p_x0 = x0;
  if(p_x1) *p_x1 = x1;

/*  printf(" <%d> ",i);*/
  
/*  if( i >= maxIter )
  {
    sprintf( msg, "MtSolve(): too many iterations. current error=%f",y);
    printf( msg );
    return 1;
  }*/

return SUCCESS;
}                

/********************************************************************/
/*								    */
/********************************************************************/			 
int MtSolveBinFl( float    (* f)( float x, float *parms ),  
                float   *parms,
                float    x0,
                float    x1,
                float    epsilon,
                int       maxIter,
                float   *solution,
		float   *pEps )
{
  float x=0,y0,y1,y=0;
  int    i;

  i = 0;
  y0 = f( x0, parms );
  y1 = f( x1, parms );
  
  for( i=0; i<maxIter; i++)
  {
    x = x0 + ( x1 - x0 ) / 2;
    y = f( x, parms );
    if( fabs( y ) < epsilon ) break;    
    if( y*y0 < 0 )
    {
      x1 = x; y1 = y;
    }
    else
    {
      x0 = x; y0 = y;
    }  
  }

  *solution = x;
  *pEps = fabs(y);

/*  printf(" <%d> ",i);*/
#if 0  
  if( i >= maxIter )
  {
    sprintf( msg, "MtSolve(): too many iterations. current error=%f",y);
    printf( msg );
    return 1;
  }
#endif

return SUCCESS;
}                

/********************************************************************/
/*								    */
/********************************************************************/			 
static double rootN1Func( double x, double *p )
{
  return (( pow( x, p[0] ) - 1 ) / ( x - 1 )) - p[1];
}
/********************************************************************/
/*								    */
/********************************************************************/			 
static double rootN1FuncD( double x, double *p )
{
  return (pow(x,p[0]-1)*((p[0]-1)*x-p[0])+1)/((x-1)*(x-1));
}                   
/********************************************************************/
/*								    */
/********************************************************************/			 
double MtRootN1( int n, double c )
{
/* solve (x^n-1)/(x-1) = c */
double parms[2];

parms[0] = (double)n;
parms[1] = c;

return MtSolveDiff( rootN1Func, rootN1FuncD, parms, 1.1, .001, 1000 );
}
  
/********************************************************************/
/*								    */
/********************************************************************/			 
int tprob(double t, int dof, double *probability)

/* Arguments:
	t : wanna evaluate integral from t to +infinity.
	dof: degrees of freedom of this t distribution.
	probability : computed value might be returned.

	function returns error code 1 if something is wrong,
	0 if all is well.

	This functional interface is not exactly what AS 3 used to be.
*/

				/* Abramowitz & Stegun, of course. */
{
	double d_dof, s, c, f, a, b;
	int fk, ks, im2, ioe, k;

	if (dof < 1) return 1;
	d_dof = (double) dof;	/* d_dof is F of fortran code. */

	a = t / sqrt(d_dof);
	b = d_dof/(d_dof + (t*t));
	im2 = dof - 2;
	ioe = dof % 2;
	s = c = f = 1.0;
	fk = ks = 2 + ioe;
	if (im2 > 2)
		for (k=ks; k<=im2; k+=2) {
			c = c*b*(fk - 1.0)/fk;
			s += c;
			if (s == f) break;	/* == ? */
			f = s;
			fk += 2.0;
		}
	if (ioe != 1) {	/* Label 20 of fortran code. */
		*probability = 0.5 + (0.5*a*sqrt(b)*s);
		return 0;
	}
	else {	/* Label 30 of fortran code. */
		if (dof == 1) s = 0.0;
		*probability = 0.5 + ((a*b*s + atan(a))*ONEBYPI);
		return 0;
	}
}

/********************************************************************/			 
int MtTest( void ) 
/********************************************************************
*********************************************************************/
{
  return 1;
}



/*===================================================

   Routine to perform LU decomposition of a matrix. (float version: deprecated)

=====================================================*/

#define TINY MT_REAL_EPSILON
/********************************************************************/			 
int MtMatLUdcmp_f(MtMatrix a,int n,int *indx, float *d, unsigned opt )
/********************************************************************
*********************************************************************/
{
  MtVector vv = MtVecOpen(n);  
  int		i,imax,j,k,rcThis=0;
  double	big,dum,sum,temp;

	imax = 0;
	*d=1.0;			/* No row interchanges yet. */
/* 
 *	Loop over rows to get the implicit scaling information. 
 */
	for (i=0; i<n; i++) { 
		big=0.0;
		for (j=0; j<n; j++)
			if ((temp=fabs(MT_MATEL(a,n,i,j))) > big) big=temp;
		if (big == 0.0) 
		{
		  if( !(opt & MT_QUIETERR ))
		  {
		    TRCERR(("Singular matrix in routine ludcmp"));
		  }
		  raiseRc( MI_EMATH );
		}
	/* No nonzero largest element. */
		vv[i]=1.0/big;	/* Save the scaling. */
	}
	for (j=0; j<n; j++) {	/* the loop over columns of Crout's method. */
		for (i=0; i<j; i++) {	/* equation (2.3.12) except for i=j */
			sum=MT_MATEL(a,n,i,j);
			for (k=0; k<i; k++) sum -= MT_MATEL(a,n,i,k)*MT_MATEL(a,n,k,j);
				MT_MATEL(a,n,i,j)=sum;
		}
		big=0.0;	/* Initialize search for largest pivot element. */
		for (i=j;i<n;i++) { /* i=j of equation (2.3.12) and i=j+1...N */
							/* of equation (2.3.13). */
			sum=MT_MATEL(a,n,i,j);
			for (k=0; k<j; k++)
				sum -= MT_MATEL(a,n,i,k)*MT_MATEL(a,n,k,j);
			MT_MATEL(a,n,i,j)=sum;
/* 
 * Is the figure of merit for the pivot better than the best so far?
 */
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) { /*	do we need to interchange rows? */
   			for (k=0; k<n; k++) { /* Yes, do so... */
				dum=MT_MATEL(a,n,imax,k);
				MT_MATEL(a,n,imax,k)=MT_MATEL(a,n,j,k);
				MT_MATEL(a,n,j,k)=dum;
			}
			*d = -(*d);		/* ...and change the parity of d. */
			vv[imax]=vv[j];	/* Also interchange the scale factor. */
		}
		indx[j]=imax;
		if (MT_MATEL(a,n,j,j) == 0.0) MT_MATEL(a,n,j,j)=TINY;
/*
 *	If the pivot element is zero the matrix is singular (at least to the 
 *	precision of the algorithm). For some applications on singular 
 *	matrices, it is desirable to substitute TINY for zero.
 */
		if (j < (n-1)) {	
/* 
 *	Now, finally, divide by the pivot element. 
 */
			dum=1.0/(MT_MATEL(a,n,j,j));
			for (i=j+1; i<n; i++) MT_MATEL(a,n,i,j) *= dum;
		}
	}	/* Go back for the next column in the reduction. */
  MtVecClose(vv);
rcCatch:
  return rcThis;
} /* ludcmp() */

#if 0
  int i,imax,j,k;
 float big;
 double dum,sum,temp;
 MtVector vv = MtVecOpen( n );


 *d=1.0;
 for (i=0;i<n;i++) {
  big=0.0;
  for (j=0;j<n;j++)
   if ((temp=fabs(MT_MATEL(a,n,i,j))) > big) big=temp;
  if (big == 0.0) 
  {
    TRCERR(("Singular matrix in routine LUDCMP"));
    exit(1);
  }
  vv[i]=1.0/big;
 }
 for (j=0;j<n;j++) {
  for (i=0;i<j;i++) {
   sum=MT_MATEL(a,n,i,j);
   for (k=0;k<i;k++) sum -= MT_MATEL(a,n,i,k)*MT_MATEL(a,n,k,j);
   MT_MATEL(a,n,i,j)=sum;
  }
  big=0.0;
  for (i=j;i<n;i++) {
   sum=MT_MATEL(a,n,i,j);
   for (k=0;k<j;k++)
    sum -= MT_MATEL(a,n,i,k)*MT_MATEL(a,n,k,j);
   MT_MATEL(a,n,i,j)=sum;
   if ( (dum=vv[i]*fabs(sum)) >= big) {
    big=dum;
    imax=i;
   }
  }
  if (j != imax) {
   for (k=0;k<n;k++) {
    dum=MT_MATEL(a,n,imax,k);
    MT_MATEL(a,n,imax,k)=MT_MATEL(a,n,j,k);
    MT_MATEL(a,n,j,k)=dum;
   }
   *d = -(*d);
   vv[imax]=vv[j];
  }
  indx[j]=imax;
/*  printf("indx[%d]=%d\n",j,imax);*/
  if (MT_MATEL(a,n,j,j) == 0.0) MT_MATEL(a,n,j,j)=TINY;
  if (j != n-1) {
   dum=1.0/( MT_MATEL(a,n,j,j));
   for (i=j+1;i<n;i++)
     MT_MATEL(a,n,i,j) *= dum;
         }
 } 
 MtVecClose( vv );
#endif

#undef TINY

/*===================================================
  =                                                 =
  =  Routine to backsubstitute a matrix solution.   =
  =                                                 =
  ===================================================*/
/********************************************************************/			 
void MtMatLubksb_f(MtMatrix a,int n,int *indx, float *b)
/********************************************************************
*********************************************************************/
{
 int i,ii=-1,ip,j;
 double sum;

 for (i=0;i<n;i++) {
  ip=indx[i];
  sum=b[ip];
  b[ip]=b[i];
  if (ii>=0)
   for (j=ii;j<i;j++) sum -= MT_MATEL(a,n,i,j)*b[j];
  else if (sum) ii=i;
  b[i]=sum;
 }
 for (i=n-1;i>=0;i--) {
  sum=b[i];
  for (j=i+1;j<n;j++) sum -= MT_MATEL(a,n,i,j)*b[j];
  b[i]=sum/MT_MATEL(a,n,i,i);
 }
}


/*==========================================================================
 *
 *
 *
 *
 *
 *
 *
 * ========================================================================*/

 /*===================================================

   Routine to perform LU decomposition of a matrix.

=====================================================*/

#define TINY MT_REAL_EPSILON
/********************************************************************/			 
int MtMatLUdcmp(double *A,int n,int *indx, double *d, unsigned opt )
/********************************************************************
*********************************************************************/
{
  double	*vv=NULL,*a=A;
  int		i,imax,j,k,rcThis=0;
  double	big,dum,sum,temp;

  ALLOCARR(vv,n);

	imax = 0;
	*d=1.0;			/* No row interchanges yet. */
/* 
 *	Loop over rows to get the implicit scaling information. 
 */
	for (i=0; i<n; i++) { 
		big=0.0;
		for (j=0; j<n; j++)
			if ((temp=fabs(MT_MATEL(a,n,i,j))) > big) big=temp;
		if (big == 0.0) 
		{
		  if( !(opt & MT_QUIETERR ))
		  {
		    TRCERR(("Singular matrix in routine ludcmp"));
		  }
		  raiseRc( MI_EMATH );
		}
	/* No nonzero largest element. */
		vv[i]=1.0/big;	/* Save the scaling. */
	}
	for (j=0; j<n; j++) {	/* the loop over columns of Crout's method. */
		for (i=0; i<j; i++) {	/* equation (2.3.12) except for i=j */
			sum=MT_MATEL(a,n,i,j);
			for (k=0; k<i; k++) sum -= MT_MATEL(a,n,i,k)*MT_MATEL(a,n,k,j);
				MT_MATEL(a,n,i,j)=sum;
		}
		big=0.0;	/* Initialize search for largest pivot element. */
		for (i=j;i<n;i++) { /* i=j of equation (2.3.12) and i=j+1...N */
							/* of equation (2.3.13). */
			sum=MT_MATEL(a,n,i,j);
			for (k=0; k<j; k++)
				sum -= MT_MATEL(a,n,i,k)*MT_MATEL(a,n,k,j);
			MT_MATEL(a,n,i,j)=sum;
/* 
 * Is the figure of merit for the pivot better than the best so far?
 */
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) { /*	do we need to interchange rows? */
   			for (k=0; k<n; k++) { /* Yes, do so... */
				dum=MT_MATEL(a,n,imax,k);
				MT_MATEL(a,n,imax,k)=MT_MATEL(a,n,j,k);
				MT_MATEL(a,n,j,k)=dum;
			}
			*d = -(*d);		/* ...and change the parity of d. */
			vv[imax]=vv[j];	/* Also interchange the scale factor. */
		}
		indx[j]=imax;
		if (MT_MATEL(a,n,j,j) == 0.0) MT_MATEL(a,n,j,j)=TINY;
/*
 *	If the pivot element is zero the matrix is singular (at least to the 
 *	precision of the algorithm). For some applications on singular 
 *	matrices, it is desirable to substitute TINY for zero.
 */
		if (j < (n-1)) {	
/* 
 *	Now, finally, divide by the pivot element. 
 */
			dum=1.0/(MT_MATEL(a,n,j,j));
			for (i=j+1; i<n; i++) MT_MATEL(a,n,i,j) *= dum;
		}
	}	/* Go back for the next column in the reduction. */
  FREEMEM(vv);
rcCatch:
  return rcThis;
} /* ludcmp() */

#if 0
  int i,imax,j,k;
 float big;
 double dum,sum,temp;
 MtVector vv = MtVecOpen( n );


 *d=1.0;
 for (i=0;i<n;i++) {
  big=0.0;
  for (j=0;j<n;j++)
   if ((temp=fabs(MT_MATEL(a,n,i,j))) > big) big=temp;
  if (big == 0.0) 
  {
    TRCERR(("Singular matrix in routine LUDCMP"));
    exit(1);
  }
  vv[i]=1.0/big;
 }
 for (j=0;j<n;j++) {
  for (i=0;i<j;i++) {
   sum=MT_MATEL(a,n,i,j);
   for (k=0;k<i;k++) sum -= MT_MATEL(a,n,i,k)*MT_MATEL(a,n,k,j);
   MT_MATEL(a,n,i,j)=sum;
  }
  big=0.0;
  for (i=j;i<n;i++) {
   sum=MT_MATEL(a,n,i,j);
   for (k=0;k<j;k++)
    sum -= MT_MATEL(a,n,i,k)*MT_MATEL(a,n,k,j);
   MT_MATEL(a,n,i,j)=sum;
   if ( (dum=vv[i]*fabs(sum)) >= big) {
    big=dum;
    imax=i;
   }
  }
  if (j != imax) {
   for (k=0;k<n;k++) {
    dum=MT_MATEL(a,n,imax,k);
    MT_MATEL(a,n,imax,k)=MT_MATEL(a,n,j,k);
    MT_MATEL(a,n,j,k)=dum;
   }
   *d = -(*d);
   vv[imax]=vv[j];
  }
  indx[j]=imax;
/*  printf("indx[%d]=%d\n",j,imax);*/
  if (MT_MATEL(a,n,j,j) == 0.0) MT_MATEL(a,n,j,j)=TINY;
  if (j != n-1) {
   dum=1.0/( MT_MATEL(a,n,j,j));
   for (i=j+1;i<n;i++)
     MT_MATEL(a,n,i,j) *= dum;
         }
 } 
 MtVecClose( vv );
#endif

#undef TINY


#define TINY MT_REAL_EPSILON

/********************************************************************/			 
int MtMatLUdcmp_tst(double *A,int n,int *indx, double *d, unsigned opt )
/********************************************************************
*********************************************************************/
{
  double	*vv=NULL,*a=A;
  int		i,imax,j,k,rcThis=0;
  double	big,dum,sum,temp;

  ALLOCARR(vv,n);

	imax = 0;
	*d=1.0;			/* No row interchanges yet. */
/* 
 *	Loop over rows to get the implicit scaling information. 
 */
	for (i=0; i<n; i++) { 
		big=0.0;
		for (j=0; j<n; j++)
			if ((temp=fabs(MT_MATEL(a,n,i,j))) > big) big=temp;
		if (fabs(big) <TINY) 
		{
		  if( !(opt & MT_QUIETERR ))
		  {
		    TRCERR(("Singular matrix in routine ludcmp"));
		  }
		  raiseRc( MI_EMATH );
		}
	/* No nonzero largest element. */
		vv[i]=1.0/big;	/* Save the scaling. */
	}
	for (j=0; j<n; j++) {	/* the loop over columns of Crout's method. */
		for (i=0; i<j; i++) {	/* equation (2.3.12) except for i=j */
			sum=MT_MATEL(a,n,i,j);
			for (k=0; k<i; k++) sum -= MT_MATEL(a,n,i,k)*MT_MATEL(a,n,k,j);
				MT_MATEL(a,n,i,j)=sum;
		}
		big=0.0;	/* Initialize search for largest pivot element. */
		for (i=j;i<n;i++) { /* i=j of equation (2.3.12) and i=j+1...N */
							/* of equation (2.3.13). */
			sum=MT_MATEL(a,n,i,j);
			for (k=0; k<j; k++)
				sum -= MT_MATEL(a,n,i,k)*MT_MATEL(a,n,k,j);
			MT_MATEL(a,n,i,j)=sum;
/* 
 * Is the figure of merit for the pivot better than the best so far?
 */
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) { /*	do we need to interchange rows? */
   			for (k=0; k<n; k++) { /* Yes, do so... */
				dum=MT_MATEL(a,n,imax,k);
				MT_MATEL(a,n,imax,k)=MT_MATEL(a,n,j,k);
				MT_MATEL(a,n,j,k)=dum;
			}
			*d = -(*d);		/* ...and change the parity of d. */
			vv[imax]=vv[j];	/* Also interchange the scale factor. */
		}
		indx[j]=imax;
		if (fabs(MT_MATEL(a,n,j,j)) < TINY) 
		{
		  TRCERR(("singular matrix\n"));
		  MT_MATEL(a,n,j,j)=TINY;
		}
/*
 *	If the pivot element is zero the matrix is singular (at least to the 
 *	precision of the algorithm). For some applications on singular 
 *	matrices, it is desirable to substitute TINY for zero.
 */
		if (j < (n-1)) {	
/* 
 *	Now, finally, divide by the pivot element. 
 */
			dum=1.0/(MT_MATEL(a,n,j,j));
			for (i=j+1; i<n; i++) MT_MATEL(a,n,i,j) *= dum;
		}
	}	/* Go back for the next column in the reduction. */
  FREEMEM(vv);
rcCatch:
  return rcThis;
}


/*===================================================
  =                                                 =
  =  Routine to backsubstitute a matrix solution.   =
  =                                                 =
  ===================================================*/
/********************************************************************/			 
void MtMatLubksb(double *A,int n,int *indx, double *b)
/********************************************************************
*********************************************************************/
{
 int i,ii=-1,ip,j;
 double sum,*a=A;

 for (i=0;i<n;i++) {
  ip=indx[i];
  sum=b[ip];
  b[ip]=b[i];
  if (ii>=0)
   for (j=ii;j<i;j++) sum -= MT_MATEL(a,n,i,j)*b[j];
  else if (sum) ii=i;
  b[i]=sum;
 }
 for (i=n-1;i>=0;i--) {
  sum=b[i];
  for (j=i+1;j<n;j++) sum -= MT_MATEL(a,n,i,j)*b[j];
  b[i]=sum/MT_MATEL(a,n,i,i);
 }
}

/*=========================================================================
 *
 *
 * =======================================================================*/

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
double MtMatInverse(double *matrix, int num_dims)
{ register int i, j, k, n = num_dims;
  register double *dptr1, *dptr2;
  register int *l, *m;                          /* Row and column permutation vectors */
  double det, biga, hold, abs_biga, neg_biga;

  l = (int *) MM_calloc(num_dims, sizeof(int), MM_DUMMYTAG);
  m = (int *) MM_calloc(num_dims, sizeof(int), MM_DUMMYTAG);

  det = 1.0;

  for (k=0; k<num_dims; k++) {
    l[k] = k;
    m[k] = k;
    biga = MT_MATEL(matrix,n, k,k);
    abs_biga = abs(biga);

    /* Find the biggest element in the submatrix */
    for (i=k; i<num_dims; i++) {
      for (j = k, dptr1 = &MT_MATEL(matrix,n, i,j); j<num_dims; j++, dptr1++) {
        if (abs(*dptr1) > abs_biga) {
          biga = *dptr1;
          abs_biga = abs(biga);
          l[k] = i;
          m[k] = j;
        }
      }
    }
    /* Interchange rows */
    i = l[k];
    if (i > k)
      for (j=0, dptr1=&MT_MATEL(matrix,n, i,j), 
	  	dptr2=&MT_MATEL(matrix,n, k,j); j<num_dims; j++) {
        hold = -(*dptr2);
        *dptr2++ = *dptr1;
        *dptr1++ = hold;
      }

    /* Interchange columns */
    j = m[k];
    if (j > k)
      for (i=0; i<num_dims; i++) {
        dptr1 = &MT_MATEL(matrix,n, i,k);
        dptr2 = &MT_MATEL(matrix,n, i,j);
        hold = -(*dptr1);
        *dptr1 = *dptr2;
        *dptr2 = hold;
      }

    /* Divide column by minus pivot (value of pivot element is contained in biga). */
    if (biga == 0.0) return 0.0;

    neg_biga = -biga;

    for (i=0; i<k; i++)
      MT_MATEL(matrix,n, i,k) /= neg_biga;

    for (i=k+1; i<num_dims; i++)
      MT_MATEL(matrix,n, i,k) /= neg_biga;

    /* Reduce matrix */
    for (i=0; i<k; i++) {
      hold = MT_MATEL(matrix,n, i,k);

      for (j=0, dptr1 = &MT_MATEL(matrix,n, i,j), 
	  	dptr2 = &MT_MATEL(matrix,n, k,j); j<k; j++)
        *dptr1++ += hold * *dptr2++;

      for (dptr1++, dptr2++, j = k + 1; j<num_dims; j++)
        *dptr1++ += hold * *dptr2++;
    }

    for (i=k+1; i<num_dims; i++) {
      hold = MT_MATEL(matrix,n, i,k);

      for (j=0, dptr1 = &MT_MATEL(matrix,n, i,j), 
	  	dptr2 = &MT_MATEL(matrix,n, k,j); j<k; j++)
        *dptr1++ += hold * *dptr2++;

      for (dptr1++, dptr2++, j=k+1; j<num_dims; j++)
        *dptr1++ += hold * *dptr2++;
    }

    /* Divide row by pivot */
    for (j=0, dptr1 = &MT_MATEL(matrix,n, k,j); j<k; j++)
      *dptr1++ /= biga;

    for (dptr1++, j=k+1; j<num_dims; j++)
      *dptr1++ /= biga;

    det *= biga;                /* Product of pivots */
    MT_MATEL(matrix,n, k,k) = 1.0 / biga;
  }                             /* K loop */

  /* Final row & column interchanges */
  for (k=num_dims; k-->0; ) {
    i = l[k];
    if (i > k)
      for (j=0; j <num_dims; j++) {
        dptr1 = &MT_MATEL(matrix,n, j,k);
        dptr2 = &MT_MATEL(matrix,n, j,i);
        hold = *dptr1;
        *dptr1 = -(*dptr2);
        *dptr2 = hold;
      }
    j = m[k];
    if (j > k)
      for (i=0, dptr1 = &MT_MATEL(matrix,n, k,i), 
	  	dptr2 = &MT_MATEL(matrix,n, j,i); i<num_dims; i++) {
        hold = *dptr1;
        *dptr1++ = -(*dptr2);
        *dptr2++ = hold;
      }
  }
  MM_free(l);
  MM_free(m);

  return det;                   /* return determinant   */
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtMatInverseLU( MtMatrix b, MtMatrix a, MtVector c, 
    	            int n, int *luIdx, float *d, unsigned opt )
{
  /*
   *  b = a ^(-1)
   */

  int i,j,cAlloc=0,idxAlloc=0;
  float d0;

  if(!d) d=&d0;
  if(!c) 
  {
    cAlloc= 1;
    ALLOCARR(c,n);
  }
  if(!luIdx) 
  {
    idxAlloc= 1;
    ALLOCARR(luIdx,n);
  }

  MtMatLUdcmp_f( a, n, luIdx, d, opt );

  for(j=0;j<n;j++)
  {
    for(i=0;i<n;i++) c[i] = 0.0;
    c[j] = 1.0;
    MtMatLubksb_f( a, n, luIdx, c );
    for(i=0;i<n;i++) MT_MATEL( b, n, i, j ) = c[i];
  }

  if(cAlloc) FREEMEM(c);
  if(idxAlloc) FREEMEM(luIdx);
  return 0;
}



/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtMatInverseLU_d( double *b, double *a, double *c, 
    	            int n, int *luIdx, double *d, unsigned opt )
{
  /*
   *  b = a ^(-1)
   */

  int rcThis=0,rc,i,j,cAlloc=0,idxAlloc=0;
  double d0;

  if(!d) d=&d0;
  if(!c) 
  {
    cAlloc= 1;
    ALLOCARR(c,n);
  }
  if(!luIdx) 
  {
    idxAlloc= 1;
    ALLOCARR(luIdx,n);
  }

  rc = MtMatLUdcmp( a, n, luIdx, d, opt );
  if(rc) raiseRc(rc);

  for(j=0;j<n;j++) *d = (*d) * MT_MATEL(a,n,j,j);

  for(j=0;j<n;j++)
  {
    for(i=0;i<n;i++) c[i] = 0.0;
    c[j] = 1.0;
    MtMatLubksb( a, n, luIdx, c );
    for(i=0;i<n;i++) MT_MATEL( b, n, i, j ) = c[i];
  }

rcCatch:
  if(cAlloc) FREEMEM(c);
  if(idxAlloc) FREEMEM(luIdx);
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtMatInverseLU_d_tst( double *b, double *a, double *c, 
    	            int n, int *luIdx, double *d, unsigned opt )
{
  /*
   *  b = a ^(-1)
   */

  int rcThis=0,rc,i,j,cAlloc=0,idxAlloc=0;
  double d0;

  if(!d) d=&d0;
  if(!c) 
  {
    cAlloc= 1;
    ALLOCARR(c,n);
  }
  if(!luIdx) 
  {
    idxAlloc= 1;
    ALLOCARR(luIdx,n);
  }

  rc = MtMatLUdcmp_tst( a, n, luIdx, d, opt );
  if(rc) raiseRc(rc);

  for(j=0;j<n;j++)
  {
    for(i=0;i<n;i++) c[i] = 0.0;
    c[j] = 1.0;
    MtMatLubksb( a, n, luIdx, c );
    for(i=0;i<n;i++) MT_MATEL( b, n, i, j ) = c[i];
  }

rcCatch:
  if(cAlloc) FREEMEM(c);
  if(idxAlloc) FREEMEM(luIdx);
  return rcThis;
}

#if 0
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void svdcmp(float **a, int m, int n, float *w, float **v)
{
  /* IN: matrix A a[1..m][1..n]
   * OUT: compute SVD  A = U * W * V^T 
   *
   */
  int i,j,k;

  float g,scale,anorm;  
  float *rv1=MM_malloc(sizeof(float)*n)-1;
  g=scale=anorm=0.0;
  for(i=1;i<=n;++i){
    int l=i+1;
    rv1[i]=scale*g;
    g=scale=0.0;
    if(i<=m){
      for(k=i;k<=m;++k) scale+=fabs(a[k][i]);
      if(scale){
        float s=0.0;
        for(k=i;k<=m;++k){
          a[k][i]/=scale;
          s+=a[k][i]*a[k][i];
        }
        float f=a[i][i];
        g=-SIGN(sqrt(s),f);
        float h=f*g-s;
        a[i][i]=f-g;
        for(j=l;j<=n;++j){
          float sum=0.0;
          for(k=i;k<=m;++k) sum+=a[k][i]*a[k][j];
          float fct=sum/h;
          for(k=i;k<=m;++k) a[k][j]+=fct*a[k][i];
        }
        for(k=i;k<=m;++k) a[k][i]*=scale;
      }
    }
    w[i]=scale*g;
    g=scale=0.0;
    if((i<=m)&&(i!=n)){
      for(k=l;k<=n;++k) scale+=fabs(a[i][k]);
      if(scale){
        float s=0.0;
        for(k=l;k<=n;++k){
          a[i][k]/=scale;
          s+=a[i][k]*a[i][k];
        }
        float f=a[i][l];
        g=-SIGN(sqrt(s),f);
        float h=f*g-s;
        a[i][l]=f-g;
        for(k=l;k<=n;++k) rv1[k]=a[i][k]/h;
        for(j=l;j<=m;++j){
          float sum=0.0;
          for(k=l;k<=n;++k) sum+=a[j][k]*a[i][k];
          for(k=l;k<=n;++k) a[j][k]+=sum*rv1[k];
        }
        for(k=l;k<=n;++k) a[i][k]*=scale;
      }
    }
    anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  {
    float f;
    for(i=n,l;i>=1;--i){
      if(i<n){       // this makes f and l not dependent
        if(f){
          for(j=l;j<=n;++j) v[j][i]=(a[i][j]/a[i][l])/f;
            for(j=l;j<=n;++j){
            float sum=0.0;
            for(k=l;k<=n;++k) sum+=a[i][k]*v[k][j];
            for(k=l;k<=n;++k) v[k][j]+=sum*v[k][i];
          }
        }
        for(j=l;j<=n;++j) v[i][j]=v[j][i]=0.0;
      }
      v[i][i]=1.0;
      f=rv1[i];
      l=i;
    }
  }
  for(i=min(m,n);i>=1;--i){
    int l=i+1;
    g=w[i];
    for(j=l;j<=n;++j) a[i][j]=0.0;
    if(g){
      g=1.0/g;
      for(j=l;j<=n;++j){
        float sum=0.0;
        for(k=l;k<=m;++k) sum+=a[k][i]*a[k][j];
        float f=(sum/a[i][i])*g;
        for(k=i;k<=m;++k) a[k][j]+=f*a[k][i];
      }
      for(j=i;j<=m;++j) a[j][i]*=g;
    }else for(j=i;j<=m;++j) a[j][i]=0.0;
    ++a[i][i];
  }
  for(k=n;k>=1;--k){
    for(its=1;its<=30;++its){
      int flag=1,nm,l;
      for(l=k;l>=1;--l){
        nm=l-1;
        if((float)(fabs(rv1[l])+anorm)==anorm){
          flag=0;
          break;
        }
        if((float)(fabs(w[nm])+anorm)==anorm) break;
      }
      if(flag){
        float c=0.0,s=1.0;
        for(i=l;i<=k;++i){
          float f=s*rv1[i];
          rv1[i]=c*rv1[i];
          if((float)(fabs(f)+anorm)==anorm) break;
          g=w[i];
          float h=pythag(f,g);
          w[i]=h;
          h=1.0/h;
          c=g*h;
          s=-f*h;
          for(j=1;j<=m;++j){
            float y=a[j][nm];
            float z=a[j][i];
            a[j][nm]=y*c+z*s;
            a[j][i]=z*c-y*s;
          }
        }
      }
      float z=w[k];
      if(l==k){
        if(z<0.0){
          w[k]=-z;
          for(j=1;j<=n;++j) v[j][k]=-v[j][k];
        }
        break;
      }
      if(its==30) 
      {
	TRCERR(("no convergence in 30 svdcmp iterations"));
	exit(1);
      }
      float x=w[l];
      float y=w[nm=k-1];
      g=rv1[nm];
      float h=rv1[k];
      float f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      float c=1.0,s=1.0;
      for(j=l;j<=nm;++j){
        i=j+1;
        g=rv1[i];
        y=w[i];
        h=s*g;
        g*=c;
        z=pythag(f,h);
        rv1[j]=z;
        c=f/z;
        s=h/z;
        f=x*c+g*s;
        g=g*c-x*s;
        h=y*s;
        y*=c;
        for(jj=1;jj<=n;++jj){
          x=v[jj][j];
          z=v[jj][i];
          v[jj][j]=x*c+z*s;
          v[jj][i]=z*c-x*s;
        }
        z=pythag(f,h);
        w[j]=z;
        if(z){
          z=1.0/z;
          c=f*z;
          s=h*z;
        }
        f=c*g+s*y;
        x=c*y-s*g;
        for(jj=1;jj<=m;++jj){
          y=a[jj][j];
          z=a[jj][i];
          a[jj][j]=y*c+z*s;
          a[jj][i]=z*c-y*s;
        }
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  MM_free (rv1+1);
}

/********************************************************************/			 
void MtMatSVD_Old(MtMatrix a, int n, int m, float *w, float **v )
// TODO: Deprecated: use NR3-Version instead
/********************************************************************
*********************************************************************/
{
  float *a1=NULL,*w1=NULL,*v1=NULL,*src,*dst;
  int	 i,j;

  ALLOCARR(a1,(m+1)*(n+1));
  ALLOCARR(w1,m+1);
  ALLOCARR(v1,(m+1)*(m+1));

  /*** convert to 1..n matrix format (TODO: performance ) */
  for(i=0;i<n;i++)
    MtVecCopy_f( MT_MATEL(a1,m+1,i+1,1), MT_MATEL(a,m,i,0), m );

  svdcmp(a1,n,m,w1,v1);

  /*** convert back to 0..n-1 matrix format */
  for(i=0;i<n;i++)
    MtVecCopy_f( MT_MATEL(a,m,i,0), MT_MATEL(a1,m+1,i+1,1), m );
  for(i=0;i<m;i++)
    MtVecCopy_f( MT_MATEL(v,m,i,0), MT_MATEL(v1,m+1,i+1,1), m );
  MtVecCopy_f( v, v1+1, m );

  FREEMEM(a1); FREEMEM(v1); FREEMEM(w1);
}

#endif


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtMatSimpleSVD(double *U, double *S, double *V,
				int nRow, int nCol)
{
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Original comments from Bryant Marks:
 
   This SVD routine is based on pgs 30-48 of "Compact Numerical Methods
   for Computers" by J.C. Nash (1990), used to compute the pseudoinverse.
   Modifications include:
        Translation from Pascal to ANSI C.
        Array indexing from 0 rather than 1.
        Float replaced by double everywhere.
        Support for the Matrix structure.
        I changed the array indexing so that the matricies (float [][])
           could be replaced be a single list (double *) for more
           efficient communication with Mathematica.

  From: bryant@sioux.stanford.edu (Bryant Marks)
  >
  > A couple of things to note: A needs to have twice as much room
  > allocated for it (2*n + 2*m) since the W in the svd function requires
  > this (part of a rotation algorithm).  After the routine has run W
  > contains two maticies of the decomposition A = USV'.  The first nRow
  > rows contain the product US and the next nCol rows contain V (not V').
  > Z is equal to the vector of the sqares of the diagonal elements of S. */

/* Comments from GWF: The note above is not strictly correct.  To compute
   the SVD of an (m x n) matrix, W must be ((m + n) x n) in size.

   12-26-96: I slightly rewrote things to include V so that the interface
   could at least be close to a lapack-like routine. 

   01-02-97: Rewrote it so that U,S,V was actually returned.  V is now optional. */

#define TOLERANCE 1.0e-12

  int i, j, k, EstColRank, RotCount, SweepCount, slimit, rcThis=0;
  double eps, e2, tol, vt, p, x0, y0, q, r, c0, s0, d1, d2;

  eps = TOLERANCE;
  slimit = nCol / 4;
  if (slimit < 6.0)
    slimit = 6;
  SweepCount = 0;
  e2 = 10.0 * nRow * eps * eps;
  tol = eps * .1;
  EstColRank = nCol;
  if(V)
    for (i = 0; i < nCol; i++)
      for (j = 0; j < nCol; j++) {
	V[nCol * i + j] = 0.0;
	V[nCol * i + i] = 1.0;
      }
  RotCount = EstColRank * (EstColRank - 1) / 2;
  while (RotCount != 0 && SweepCount <= slimit) {
    RotCount = EstColRank * (EstColRank - 1) / 2;
    SweepCount++;
    for (j = 0; j < EstColRank - 1; j++) {
      for (k = j + 1; k < EstColRank; k++) {
	p = q = r = 0.0;
	for (i = 0; i < nRow; i++) {
	  x0 = U[nCol * i + j];
	  y0 = U[nCol * i + k];
	  p += x0 * y0;
	  q += x0 * x0;
	  r += y0 * y0;
	}
	S[j] = q;
	S[k] = r;
	if (q >= r) {
	  if (q <= e2 * S[0] || fabs(p) <= tol * q)
	    RotCount--;
	  else {
	    p /= q;
	    r = 1 - r / q;
	    vt = sqrt(4 * p * p + r * r);
	    c0 = sqrt(fabs(.5 * (1 + r / vt)));
	    s0 = p / (vt * c0);
	    for (i = 0; i < nRow; i++) {
	      d1 = U[nCol * i + j];
	      d2 = U[nCol * i + k];
	      U[nCol * i + j] = d1 * c0 + d2 * s0;
	      U[nCol * i + k] = -d1 * s0 + d2 * c0;
	    }
	    if(V)
	      for (i = 0; i < nCol; i++) {
		d1 = V[nCol * i + j];
		d2 = V[nCol * i + k];
		V[nCol * i + j] = d1 * c0 + d2 * s0;
		V[nCol * i + k] = -d1 * s0 + d2 * c0;
	      }
	  }
	}
	else {
	  p /= r;
	  q = q / r - 1;
	  vt = sqrt(4 * p * p + q * q);
	  s0 = sqrt(fabs(.5 * (1 - q / vt)));
	  if (p < 0)
	    s0 = -s0;
	  c0 = p / (vt * s0);
	  for (i = 0; i < nRow; i++) {
	    d1 = U[nCol * i + j];
	    d2 = U[nCol * i + k];
	    U[nCol * i + j] = d1 * c0 + d2 * s0;
	    U[nCol * i + k] = -d1 * s0 + d2 * c0;
	  }
	  if(V)
	    for (i = 0; i < nCol; i++) {
	      d1 = V[nCol * i + j];
	      d2 = V[nCol * i + k];
	      V[nCol * i + j] = d1 * c0 + d2 * s0;
	      V[nCol * i + k] = -d1 * s0 + d2 * c0;
	    }
	}
      }
    }
    while (EstColRank >= 3 && S[(EstColRank - 1)] <= S[0] * tol + tol * tol)
      EstColRank--;
  }
  for(i = 0; i < nCol; i++)
    S[i] = S[i] > MT_REAL_EPSILON ? sqrt(S[i]) : 0;
  for(i = 0; i < nCol; i++)
    for(j = 0; j < nRow; j++)
      U[nCol * j + i] = 
	S[i] > MT_REAL_EPSILON ? U[nCol * j + i] / S[i] : 0;
//rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtMatSVD(
    double *A, int n, int m, // IN: matrix A (nxm)
    double *U, 	// OUT: matrix U (nxm)
    double *V,	// OUT: matrix V (mxm)
    double *S	// OUT: vector of m singular values
    )
{
#if 1
  UT_ASSERT0(FALSE);
#else
  return NR3_SVD_Decomp(A,n,m,U,V,S);
#endif
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtMatPseudoInverseSVD( double *A, int n, 
    			   double *S, // optional: (nx1) eigenvalues 
			   double *V, // optional: (nxn) eigenvectors (=columns)
			   double eps,  // optional
    			   unsigned opt )
//void pinv(double **A, double **Ainv, unsigned n)
{
#define MACH_EPS 1e-10
  unsigned i, j, k, V_alloc=FALSE,rcThis=0;
  int	rc;
  double *U=NULL, *S2=NULL, smin, smax, sum;

  if(!(eps>0)) eps = MACH_EPS;

  USE(rc);

  opt |= MT_NR3; // always use the NR3 version

  if(opt & MT_CHECK)
  {
    TRCERRR(("option not supported\n"),MI_ERROR);
  }

  ALLOCARR(U,n*n);
  if(!V)
  {
    ALLOCARR(V,n*n);
    V_alloc=TRUE;
  }
 
  ALLOCARR(S2,n);
  
  if(opt & MT_NR3 )
  {
    UT_ASSERT0(FALSE);
//    printf("NR3\n");
//    rc = NR3_SVD_Decomp( A, n, n, U, V, S2 );
    UT_ASSERT(rc==0);
  }
  else
  {
    for(i = 0; i < n; i++)
      for(j = 0; j < n; j++)
	MT_MATEL(U,n,i,j) = MT_MATEL(A,n,i,j);

    MtMatSimpleSVD(U, S2, V, n, n);
  }

#if 0
  {
    int m=n;
    double f,sum,esum=0,err;
    MtStat estat;
    double *A2=NULL;

    ALLOCARR(A2,n*n);
    
    MT_STAT_INIT(&estat);

    for(i=0;i<m;i++)
      for(j=0;j<m;j++)
      {
	f = MT_MATEL(A,m,i,j);
	for(sum=0,k=0;k<m;k++)
	  sum += MT_MATEL(V,m,j,k) *
	    	 MT_MATEL(U,m,i,k) * 
	         S2[k];
	MT_MATEL(A2,n,i,j) = sum;
	err = fabs(f-sum);
	MT_STAT_UPDATE(&estat,err);
	esum += err*err;
      }

    MT_STAT_RESULT(&estat);

    printf(">>>A=\n");MtMatPrint_d(A,n,n,"%6.2f ");
    printf(">>>U=\n");MtMatPrint_d(U,n,n,"%6.2f ");
    printf(">>>S2=\n"); for(i=0;i<n;i++) printf("%6.2f ",S2[i]); printf("\n");
    printf(">>>V=\n");MtMatPrint_d(V,n,n,"%6.2f ");
    printf(">>>A2=\n");MtMatPrint_d(A2,n,n,"%6.2f ");
    FREEMEM(A2);
//    if(opt & RND_VERBOSE)
    {
      TRC1(("SVD (%dx%d) verification:\n",
	    m,m));
      TRC1(("Elem-Err: "));MtStatPrint(&estat);
      TRC1(("Total error: %g\n",
	    sqrt(esum)/(double)(m*m)));
    }

    if(estat.max>0.001)
    {
      TRCERR(("Numerical verification error (e=%g)\n",
	    estat.max));
      if(estat.max>0.01) raiseRc(MT_ERROR);
    }

  }
#endif

  if(S)
  {
    for(i=0;i<n;i++) S[i] = S2[i];
  }

  /* Zero out tiny values */
  smax = 0.0;
  for(i = 0; i < n; i++)
    if(S2[i] > smax)
      smax = S2[i];
  smin = fabs(sqrt(smax) * n * eps);
  for(i = 0; i < n; i++)
    if(fabs(S2[i]) < smin)
      S2[i] = 0.0;
    else
      S2[i] = 1.0 / S2[i];

  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++) {
      sum = 0;
      for(k = 0; k < n; k++)
	sum += MT_MATEL(V,n,i,k) * MT_MATEL(U,n,j,k) * S2[k];
      MT_MATEL(A,n,i,j) = sum;
    }

rcCatch:
  FREEMEM(U);
  if(V_alloc)
  {
    FREEMEM(V);
  }
  FREEMEM(S2);
  return rcThis;
}


/********************************************************************/			 
double MtMatDet(MtMatrix a,int n)
/********************************************************************
*********************************************************************/
{
  int j, *idx=NULL;
  float d;

  ALLOCMEM( idx, 2*sizeof(*idx)*n );
  MtMatLUdcmp_f(a,n,idx,&d,0);
  for(j=0;j<n;j++) d*=MT_MATEL(a,n,j,j);
  FREEMEM( idx );
  return d;
}


#if 0
#include <stdio.h>
#include <SPTK.h>

int toeplitz(t, a, b, n, eps)
double *t, *a, *b, eps;
int n;
{
    register int 	l, k;
    static double 	*c = NULL, *cc;
    static int		size;
    double 	  	rmd, mue, mue2;

    if (c == NULL){
	c = dgetmem(n+n+2);
	cc = c + n;
	size = n;
    }
    if (n > size){
	MM_free(c);
	c = dgetmem(n+n+2);
	cc = c + n;
	size = n;
    }

    if (eps < 0.0) eps = 1.0e-6;

    fillz(c, sizeof(*c), n+1);

    rmd = t[0];
    if (((rmd < 0.0) ? -rmd : rmd) <= eps) return(-1);

    a[0] = b[0] / rmd;

    for (l=1; l<n; l++){
	mue = -t[l];
	for (k=1; k<l; k++)
	    mue -= c[k] * t[l-k];
	mue /= rmd;

	for (k=1; k<l; k++)
	    cc[k] = c[k] + mue * c[l-k];
	cc[l] = mue;

	rmd = (1.0 - mue*mue) * rmd;
	if (((rmd < 0.0) ? -rmd : rmd) <= eps) return(-1);

	for (k=1; k<=l; k++) c[k] = cc[k];

	mue2 = b[l];
	for (k=0; k<=l-1; k++)
	    mue2 += c[l-k] * b[k];
	mue2 /= rmd;

	for (k=0; k<l; k++)
	    a[k] += mue2 * c[l-k];
	a[l] = mue2;
    }
    return(0);
}

#endif

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtLMSRegression( 
    	float  **pSmp,  /* samples for model F(s[0]..s[m-1]) = s[m] */
	int      m,	/* model dimension */
	int      nSmp,	/* # samples */
	float	 alph,	/* step width scaling */
	float    sig,	/* sample variance sigma (-1: not supplied)*/
	int	 maxIter,
	float	 e0,
	float   *w	/* configuration (OUT) */ )
{
/*#define LMS_TEST*/
  int i,k,tMax,is,n,t;
  float *vs,s,sum,sum2,eta,e,y,eLag,beta=0.98;
#ifdef LMS_TEST
  float wOld[1000];
static PlGWin *gwin=NULL;
  PlLineHd plh2;
  if(!gwin) gwin = PlOpenGWin(500,350,"plot");
#endif
  if(sig<0)
  {
    sum  = 0; sum2 =0; n=0;
    for(i=0;i<nSmp;i++)
    {
      vs = pSmp[i];
      for(k=0;k<m;k++)
      {
	s = *vs++;
	sum += s; sum2 += s*s;
	n++;
      }
    }
    sig = sum2 / (float)nSmp;
    if(sig<MT_REAL_EPSILON) sig = 0.001; /* TODO */
  }

  tMax = maxIter;

  /**** step width */
  eta = alph * 2.0 / sig;

  /*TRC1(("sig=%f eta=%f\n",sig,eta));*/

#ifdef LMS_TEST
  PlPlotInit(&gwin->pl,&gwin->pw,0,tMax/2,0,20,0);
  GWselect(gwin->gw); GWclear(-1); GWerase(0,0);
  PlColor(GWkrgb(40,40,40)); 
  PlPlotAxes(&gwin->pl,PL_AXDEFAULT,10,5);
  PlPlotLineInit( &gwin->plh );
  PlPlotLineInit( &plh2 );
#endif

  eLag = 1;

  /**** itarative adaption */
  is=0;
  for(t=0;;t++)
  {
    if(t>=tMax) break;
    vs = pSmp[is++];
    if(is>=nSmp) is=0;

    for(k=0,y=0;k<m;k++) y+=w[k]*vs[k];
    e = vs[m] - y;

#if 1    
    eLag = t>0 ? beta*eLag + (1-beta)*fabs(e) : fabs(e);
#endif

#ifdef LMS_TEST
    MtVecCopy_f(wOld,w,m);
#endif

    for(k=0;k<m;k++) 
    {
      w[k] += eta * e * vs[k];
    }

#ifdef LMS_TEST

    if(t<tMax/2)
    {
      PlColor(GWkrgb(200,200,200));
      PlPlotLine(&gwin->pl,t,fabs(e),&gwin->plh);
      PlColor(GWkrgb(0,200,0));
      PlPlotLine(&gwin->pl,t,eLag,&plh2);
    }
    if(t % 10 == 0 )
    {
/*      printf("%f %f -> %f %f\n",vs[0],vs[1],vs[2],y);*/
      printf("#%d -> dw=%f, w=%f,%f e=%f\n",
	  t,(float)MtVecDistance(w,wOld,m),w[0],w[1],(float)e);
    }

#endif
  }
  return 0;
}


/* 
  LU-decomposition according to Crout's algorithm with pivoting. 
	Description:
   int LU_decompos(double **a,int n,int *indx,int *d,double *vv);
	Parameters:
   a - source matrix (n x n) on input, destination on output;
   n - the matrix size;
   indx - integer array (size n) to remember permutations;
   d - on output, contains +1 or -1 for even or odd permutations number.
   vv - temporary array (size n).
	   Returns: 
   0 - the source matrix is singular (invalid for decomposition),
   1 - if OK.

  Back substitution, using LU decomposed matrix.
	Description:
  void LU_backsub(double **a,int n,int *indx,double *b);
	Parameters:
  a - the matrix decomposed by Crout;
  n - the matrix size;
  indx - permutation order obtained by decomposition algorithm;
  b - the vector (size n) to be substituted on input, the result
      of the substitution on output.
  Note: a and indx are not modified by this routine and could be 
  used in multiple calls.

  Invertation of matrix, using LU decomposed matrix.
	Description:
  void LU_invert(double **a,int n,int *indx,double **inv,double *col);
	Parameters:
  a - the matrix decomposed by Crout;
  n - the matrix size;
  indx - permutation order obtained by decomposition algorithm;
  inv - the destination matrix;
  col - temporary array (size n).
  Note: test for singularity has been already obtained on the 
  matrix decomposition, a and indx are not modified by this routine, 
  the routine uses multiple backsubstitutions (previous algorithm).

  Obtaining the matrix determinant, using LU-decomposed matrix
	Description:
  double LU_determ(double **a,int n,int *indx,int *d);
	Parameters:
  a - the matrix decomposed by Crout;
  n - the matrix size;
  indx - permutation order obtained by decomposition algorithm;
  d - the parity sign (+1 or -1) obtained at decomposition.
	Returns:
  the determinant value. Note: non-zero (the matrix cannot be 
  singular, if decomposed properly); a, indx and d are not modified 
  by this routine.

*/

/* for fabs(); inclusion could be removed if fabs() is implemented inline */ 
/* "small number" to avoid overflow in some cases */
#define TINY2 1.e-30

int LU_decompos(double *a,int n,int *indx,int *d,double *vv) {
 register int i,imax,j,k;
 double big,sum,temp;
 *d=1;
 /* search for the largest element in each row; save the scaling in the 
    temporary array vv and return zero if the matrix is singular */
 for(i=0;i<n;i++) {
  big=0.;
  for(j=0;j<n;j++) if((temp=fabs(MT_MATEL(a,n,i,j)))>big) big=temp;
  if(big==0.) return(0);
  vv[i]=big;
 }
 /* the main loop for the Crout's algorithm */
 for(j=0;j<n;j++) {
  /* this is the part a) of the algorithm except for i==j */
  for(i=0;i<j;i++) {
   sum=MT_MATEL(a,n,i,j);
   for(k=0;k<i;k++) sum-=MT_MATEL(a,n,i,k)*MT_MATEL(a,n,k,j);
   MT_MATEL(a,n,i,j)=sum;
  }
  /* initialize for the search for the largest pivot element */
  big=0.;imax=j;
  /* this is the part a) for i==j and part b) for i>j + pivot search */
  for(i=j;i<n;i++) {
   sum=MT_MATEL(a,n,i,j);
   for(k=0;k<j;k++) sum-=MT_MATEL(a,n,i,k)*MT_MATEL(a,n,k,j);
   MT_MATEL(a,n,i,j)=sum;
   /* is the figure of merit for the pivot better than the best so far? */
   if((temp=vv[i]*fabs(sum))>=big) {big=temp;imax=i;}
  }
  /* interchange rows, if needed, change parity and the scale factor */
  if(imax!=j) {
   for(k=0;k<n;k++) {temp=MT_MATEL(a,n,imax,k);
     		MT_MATEL(a,n,imax,k)=MT_MATEL(a,n,j,k);
		MT_MATEL(a,n,j,k)=temp;}
   *d=-(*d);vv[imax]=vv[j];
  }
  /* store the index */
  indx[j]=imax;
  /* if the pivot element is zero, the matrix is singular but for some 
     applications a tiny number is desirable instead */
  if(MT_MATEL(a,n,j,j)==0.) MT_MATEL(a,n,j,j)=TINY2;
  /* finally, divide by the pivot element */
  if(j<n-1) {
   temp=1./MT_MATEL(a,n,j,j);
   for(i=j+1;i<n;i++) MT_MATEL(a,n,i,j)*=temp;
  }
 }
 return(1);
}

void LU_backsub(double *a,int n,int *indx,double *b) {
 register int i,j,ip,ii=-1;
 double sum;
 /* First step of backsubstitution; the only wrinkle is to unscramble 
    the permutation order. Note: the algorithm is optimized for a 
    possibility of large amount of zeroes in b */
 for(i=0;i<n;i++) {
  ip=indx[i];
  sum=b[ip];b[ip]=b[i];
  if(ii>=0) for(j=ii;j<i;j++) sum-=MT_MATEL(a,n,i,j)*b[j];
  else if(sum) ii=i;	/* a nonzero element encounted */
  b[i]=sum;
 }
 /* the second step */
 for(i=n-1;i>=0;i--) {
  sum=b[i];
  for(j=i+1;j<n;j++) sum-=MT_MATEL(a,n,i,j)*b[j];
  b[i]=sum/MT_MATEL(a,n,i,i);
 }
}

void LU_invert(double *a,int n,int *indx,double *inv,double *col) {
 register int i,j;
 for(j=0;j<n;j++) {
  for(i=0;i<n;i++) col[i]=0.;
  col[j]=1.;
  LU_backsub(a,n,indx,col);
  for(i=0;i<n;i++) MT_MATEL(inv,n,i,j)=col[i];
 }
}

double LU_determ(double *a,int n,int *indx,int *d) {
 register int j;
 double res=(double)(*d);
 for(j=0;j<n;j++) res*=MT_MATEL(a,n,j,j);
 return(res);
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtLMSRegressionTest( void )
{
  int i,nSmp,rcThis=0, m = 3;
  float smp[] = {
                  1,2,2,
		  2,5,1,
		  2,3,2,
		  2,2,2,
		  2,4,1,
		  3,5,3,
		  4,6,2,
		  5,5,3,
		  5,6,4,
		  5,7,3,
		  6,8,4,
		  7,6,2,
		  8,4,4,
		  8,9,3,
		  9,8,4, };
  MtTab *mt=NULL;		  
  float *vw=NULL;

  m = 2;
  nSmp = ARRAY_LENGTH(smp) / (m+1);

  mt = MtCreateTab( nSmp, m+1 );

  for(i=0;i<nSmp;i++) MtVecCopy_f( mt->rows[i], AELEM(smp,m+1,i), m+1 );
  ALLOCARR( vw, m );

  MtLMSRegression( mt->rows, m, nSmp, 0.1, -1, 200, -1, vw ); 

/*rcCatch:*/
  FREEMEM( vw );
  MtDelTab(mt);
  return rcThis;
}

/*=====================================================================
 *
 *  M I S C  
 *
 *  Miscellaneous routines
 *
 * ===================================================================*/

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
double *MtVecFltToDbl( float *x, int n )
{
  int i;
  double *xd=NULL;
  ALLOCARR(xd,n);
  for(i=0;i<n;i++)xd[i]=x[i];
  return xd;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
float *MtVecDblToFlt( double *x, int n )
{
  int i;
  float *xf=NULL;
  ALLOCARR(xf,n);
  for(i=0;i<n;i++)xf[i]=x[i];
  return xf;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
double MtCheckMonotony( double *f,  int	n,
    			       MtStat *ddstat )
{
  int dt,rcThis=0,i,iPeak,iValley,j;
  double dd,peak,valley,ddSum=0,fSum=0,valleySum=0; 
  MtStat ddstat0,*stat=ddstat?ddstat:&ddstat0;

  USE(rcThis);

  peak = f[0]; iPeak = 0; 
  valley = 0; iValley = 0; valleySum=0;

  MT_STAT_INIT(stat);
  for(fSum=0,i=0;i<n;i++)fSum+=f[i];
  
  for(j=0,i=1;i<n;i++)
  {
    dd = f[i] - peak;
    if((dd<0)&&(i<n-1))
    {
      if( dd < valley )
      {
	valley = dd; iValley = i;
      }
      valleySum += f[i];
    }
    else 
    {
      dt = i-iPeak;
      peak = f[i]; iPeak = i;
      if(dt>1)
      {
	MT_STAT_UPDATE(stat,(valleySum/fSum));
	ddSum += valleySum;
	valley = 0;
	valleySum = 0;
      }
    }
  }
  MT_STAT_RESULT(stat);
/*rcCatch:*/
  return ddSum/fSum;
}

/*-------------------------------------------------------------------------*/
/*  MtCompareDistr: compares empirical distributions f1,f2		   */
/*  NOTE: f1, f2 assumed to be sorted arrays				   */		    
/*-------------------------------------------------------------------------*/
int MtCompareDistr( MtVector f1, int n1,   
    	            MtVector f2, int n2,
		    unsigned	opt,
		    double	*difCramer,
		    double	*difKomolSmir )
{
/*#define CVM_ABS*/

  float df,sum1=0,sum2=0,diff=0, df1=1./(float)n1, df2=1./(float)n2,x,xOld,
  	dfMax=0,dx;
  int	i1=0,i2=0,usef1,rcThis=0;

  xOld = 0;
  while(i1<n1||i2<n2)
  {
    if(i1==n1)
      usef1=0;
    else if(i2==n2)
      usef1=1;
    else if( f1[i1] < f2[i2] )
      usef1=1;
    else
      usef1=0;

    if(usef1)
    {
      x = f1[i1++];
      dx = (x-xOld); df = (sum1-sum2);
      sum1+=df1;
    }
    else
    {
      x = f2[i2++];
      dx = (x-xOld); df = (sum1-sum2);
      sum2+=df2;
    }

/*    diff += df*df;*/

#ifdef CVM_ABS    
    diff += fabs(df*dx);
#else
    diff += df*dx * df*dx;
#endif

    if( dx > MT_REAL_EPSILON && fabs(df) > dfMax ) dfMax = fabs(df);
    
/*    printf("x=%f s1=%f s2=%f dx=%f df=%f diff=%f\n",x,sum1,sum2,x-xOld,df,diff);*/

    xOld = x;
  }

/*rcCatch:*/

#ifdef CVM_ABS    
  *difCramer = diff;
#else
  *difCramer = sqrt(diff);
#endif

  *difKomolSmir = dfMax;
  return rcThis;
}


/*-------------------------------------------------------------------------*/
/*  MtCompareDistr: 							   */
/*  compare sum & difference histograms of bivariate distributions         */
/*-------------------------------------------------------------------------*/
int MtCompareDistrBIN256_SumDiff( float *f1, int n1,   
    	            		  float *f2, int n2,
		    		  int		 doSum,
		    		  double	*difCramer,
		    		  double	*difKomolSmir )
{
#define BIN_SIZE  256  
  int m,i,bin1[BIN_SIZE+10], bin2[BIN_SIZE+10];
  float *smp,x,y,df,f,sum1=0,sum2=0,diff=0, dfMax=0;
  int	rcThis=0;

  m = BIN_SIZE;
  memset(bin1,0,ARRAY_LENGTH(bin1)*sizeof(bin1[0]));
  memset(bin2,0,ARRAY_LENGTH(bin2)*sizeof(bin2[0]));

  for(i=0,smp=f1;i<n1;i++) 
  {
    x = *smp++; y = *smp++; 
    f = doSum ? (x + y)/4+128 : fabs(x - y);
    bin1[MAX((int)f,0)]++;  
  }
  for(i=0,smp=f2;i<n2;i++) 
  {
    x = *smp++; y = *smp++; 
    f = doSum ? (x + y)/4+128 : fabs(x - y);
    bin2[MAX((int)f,0)]++;  
  }

  for(i=0;i<m;i++)
  {
    sum1 += bin1[i]/(float)n1; sum2 += bin2[i]/(float)n2;
    df = sum1-sum2;
    diff += df*df;
    if( fabs(df) > dfMax ) dfMax = fabs(df);
  }

  if( difCramer) *difCramer = sqrt(diff);

  if(difKomolSmir) *difKomolSmir = dfMax;
  return rcThis;
}


/*-------------------------------------------------------------------------*/
/*  MtCompareDistr: compares empirical distributions f1,f2		   */
/*-------------------------------------------------------------------------*/
int MtCompareDistrBIN256( MtVector f1, int n1,   
    	            MtVector f2, int n2,
		    unsigned	opt,
		    double	*difCramer,
		    double	*difKomolSmir )
{
/*#define CVM_ABS*/

  int i,bin1[260], bin2[260];
  float df,sum1=0,sum2=0,diff=0, dfMax=0;
  int	rcThis=0, m=ARRAY_LENGTH(bin1);

  m = 256;
  memset(bin1,0,ARRAY_LENGTH(bin1)*sizeof(bin1[0]));
  memset(bin2,0,ARRAY_LENGTH(bin2)*sizeof(bin2[0]));
  for(i=0;i<n1;i++) bin1[MAX((int)f1[i]/2+128,0)]++;  
  for(i=0;i<n2;i++) bin2[MAX((int)f2[i]/2+128,0)]++;  

  for(i=0;i<m;i++)
  {
    sum1 += bin1[i]; sum2 += bin2[i];
    df = sum1-sum2;
    diff += df*df;
    if( fabs(df) > dfMax ) dfMax = fabs(df);
  }

  *difCramer = sqrt(diff);

  *difKomolSmir = dfMax;
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtTstCompareDistr( int nIter )
{
   int	i;
   float mu1,sig1,mu2,sig2;

   for(i=0;i<nIter;i++)
   {
    mu1 = MtRndGet()*100; sig1 = MtRndGet()*50;
    mu2 = MtRndGet()*100; sig2 = MtRndGet()*50;

   }
   return 0;
}


/*=====================================================================
 *
 *  MtHist: simple binned histograms
 *
 * ===================================================================*/

/*-------------------------------------------------------------------------*/
void MtDelHist( MtHist *h )
/*-------------------------------------------------------------------------*/
{
  if(h)
  {
    FREEMEM(h->sm_krn);
    FREEMEM(h->bin0);
    FREEMEM(h);
  }
}

/*-------------------------------------------------------------------------*/
MtHist *MtCreateHist( int nb,		// number of bins
    		      LongInt n_smp,	// expected number of samples
		      int use_zero_eps,	// <>0 : use eps for empty bin
    	      double xmin,
	      double xmax
	      )
/*-------------------------------------------------------------------------*/
{
  size_t binAllocLen;
  int	i;

  MtHist *h=NULL;
  double w;
  ALLOCOBJ0(h);
  h->nb = nb;
  h->xmin = xmin;
  h->xmax = xmax;
  w = xmax-xmin;
  h->n_smp = n_smp;
  if(use_zero_eps) 
    h->eps = 0.1*(double)nb/(double)n_smp;
  else
    h->eps = 0;
  binAllocLen = 2*(nb+2)+1;
  ALLOCARR(h->bin0,binAllocLen);
  h->bin = h->bin0+nb/2;   //offsets for efficient border handling 
  			   // when convolving with smoothing kernels
  for(i=0;i<binAllocLen;i++) 
  {
    h->bin0[i] = 0;
  }
  h->dx = w/(double)nb;
  h->sm_krn=NULL;
  h->sm_sig = -1;
  h->sm_r = 0;

/*rcCatch:*/
  TRC(("MtCreateHist: bins=%d, dx=%g\n",h->nb,(double)h->dx));
  return h;
}

/*-------------------------------------------------------------------------*/
int MtFillHistFun( MtHist *h, double (*func)( double x ))
/*-------------------------------------------------------------------------*/
{
  int i;
  double sum,g,x,x1,x2,xmin=h->xmin,dx=h->dx,f;
  sum=0;
  for(i=0;i<h->nb;i++)
  {
    x1 = xmin+i*dx; x2 = x1+dx;
    x = (x1+x2)/2.;
    f = func(x);
    if(fabs(f)<MT_REAL_EPSILON) f = h->eps;
    g = f*dx;
    h->bin[i] = g;
    sum+=g;    
  }
  for(i=0;i<h->nb;i++) h->bin[i] /= sum;
  return 0;
}

/*-------------------------------------------------------------------------*/
int MtIFillHist( MtHist *h, double *data, float *data_f,
    		LongInt n, double a, int do_norm)
/*-------------------------------------------------------------------------*/
{
  LongInt i,j,m,nb=h->nb;
  int rcThis=0;
  double mu,sig,s,s2,x,*bin=h->bin,eps=h->eps,
	 xmin=h->xmin,xmax=h->xmax,dx=h->dx,sum,f;
  MtStat st;

//  do_norm = 1;
  if(do_norm)
  {
    s=0;s2=0;
    MT_STAT_INIT(&st);
    for(i=0;i<n;i++) 
    {
      x= a * (data ? data[i] : data_f[i]);
      s+=x; s2+=x*x;
    }
    MT_STAT_UPDATE_N(&st,n,0,0,s,s2);
    MT_STAT_RESULT(&st);
    mu=st.avg; sig=st.var;
//    printf("mu=%g sig=%g\n",mu,sig);
  }
  else
  {
    mu=0; sig=1;
  }

  for(i=0;i<nb;i++) bin[i] = eps;
  m=0;
  for(i=0;i<n;i++)
  {
    x = a * (data ? data[i] : data_f[i]);
    if(do_norm) x = (x-mu)/sig;
    if(x<xmin || x>xmax) 
    {
      TRCERR(("sample %g out of bin range (%g,%g)\n",
	    (double)x,(double)h->xmin,(double)h->xmax)); 
      if(x<xmin) x = xmin;
      if(x>xmax) x = xmax-MT_REAL_EPSILON;
    }

    j = (x-xmin)/dx;
    if(j<0||j>nb) 
    {
      TRCERR(("bin %d out of bin range (0,%d)\n",j,h->nb-1));
    }
    else
    {
      bin[j]++;
      m++;
    }

  }
  for(sum=0,i=0;i<nb;i++) 
  {
    f = bin[i] / ((double)m*dx);
    bin[i] = f;
    sum += f;
  }
  if(!(fabs(sum-(double)nb)<0.00001))
  {
    TRCERR(("sum=%g != #bins=%g\n",sum,(double)nb));
    raiseRc(MI_ENUMERICAL);
  }
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
double MtGetPercHist( MtHist *h, double perc)
// Get percentile  
/*-------------------------------------------------------------------------*/  
{
  int i;
  double sum,x=0,xmin=h->xmin,xmax=h->xmax;

  perc *= h->nb;
  for(i=0,sum=0;i<h->nb;i++)
  {
    sum+=h->bin[i];
    if(sum>perc) break;
  }

  x = xmin + i * h->dx;

  if(x<xmin || x>xmax) 
  {
    TRCERR(("sample %g out of bin range (%g,%g)\n",
	  (double)x,(double)h->xmin,(double)h->xmax)); 
    if(x<xmin) x = xmin;
    if(x>xmax) x = xmax-MT_REAL_EPSILON;
  }

//rcCatch:
  return x;
}

/*-------------------------------------------------------------------------*/
int MtSmoothHist( MtHist *h, 	// histogram 
    		  int 	  r,	// smoothing radius (max. neighb. 
		                // bin-offset considered for smoothing )
    		  double  sig   // kernel bandwidth 
		  )
/*-------------------------------------------------------------------------*/
{
  int i,k,nb=h->nb,rcThis=0,lk;
  double *bin=h->bin,*krn,dst,sum,krnsum,*pkrn,*pbin;

  // kernel mask length 
  if(r>nb/2)
  {
    TRCERR(("MtSmoothHist: kernel radius truncated: %d-> %d\n",r,nb/2));
    r = nb/2;
  }

  // make new smoothing kernel mask ?
  if(!(( (r == h->sm_r) && h->sm_krn && 
      (fabs(h->sm_sig - sig)<MT_REAL_EPSILON))))
  {
    if(h->sm_krn) 
    {
      FREEMEM(h->sm_krn);
    }
    ALLOCARR(h->sm_krn,2*r+1);
    krn = h->sm_krn;
    h->sm_r = r;
    h->sm_sig = sig;

    sum=0;
    for(i=0;i<2*r+1;i++)
    {
      dst = fabs(i-r) / r;
      krn[i] = MtKernParzen(dst,sig);
      sum+=krn[i];
      printf("%f ",(float)krn[i]);
    }
    h->sm_krnsum = sum;
  }
  krn = h->sm_krn;
  krnsum = h->sm_krnsum;


  // convolve histogram with smooting kernel
  lk=2*r+1;
  for(i=0;i<nb;i++)
  {
    pkrn=krn; pbin=bin+i-r;
    for(sum=0,k=0;k<lk;k++) sum+=pkrn[k]*pbin[k];
    bin[i] = sum/krnsum;
  }

  return rcThis;
}  

/*-------------------------------------------------------------------------*/
double MtHistJDiv( MtHist *p, MtHist *q )
/*-------------------------------------------------------------------------*/
{
  int i;
  double d=0,f,g;

  for(i=0;i<p->nb;i++)
  {
    f = p->bin[i]; g = q->bin[i];
    if(f<MT_REAL_EPSILON) f = p->eps;
    if(g<MT_REAL_EPSILON) g = q->eps;    
    d += f * log( f/g ) + g * log( g/f );
  }

  return d;
}

/*-------------------------------------------------------------------------*/
double MtHistCMDist( MtHist *p, MtHist *q )
/*-------------------------------------------------------------------------*/
{
  int i;
  double df=0,d,s1=0,s2=0,*pb1=p->bin,*pb2=q->bin;

  for(i=0;i<p->nb;i++)
  {
    s1+=*pb1++; s2+=*pb2++;
    d = s1-s2;
    df+=d*d;
  }
//  printf("s1=%g s2=%g\n",s1,s2);

  return df;
}

/*-------------------------------------------------------------------------*/
double MtHistKSDist( MtHist *p, MtHist *q )
/*-------------------------------------------------------------------------*/
{
  int i;
  double dmax=0,d,s1=0,s2=0,*pb1=p->bin,*pb2=q->bin;

  for(i=0;i<p->nb;i++)
  {
    s1+=*pb1++; s2+=*pb2++;
    d = fabs(s1-s2);
    if(d>dmax) dmax=d;
  }
//  printf("s1=%g s2=%g\n",s1,s2);

  return dmax;
}

/*=====================================================================
 *
 *  Hist3D : 3D binned histograms (e.g. bitmap RGB,LUV data)
 *
 * ===================================================================*/

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int MtDelHist3D( MtHist3D **p_h )
{
  MtHist3D *h=*p_h;
  if(h)
  {
    FREEMEM(h->H);
    FREEMEM(h);
  }
  *p_h = NULL;
  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
MtHist3D *MtCreateHist3D( MtTab *X,   // Samples (3D, X->nCol=3)
    				      // normalized to [0..1]
    			  int nb_1d  // Number of bins in each dim
    			)
{
  int k,i,rcThis=0,d=nb_1d,nb,nb_dd,j,n,jj[3];
  double *H;
  float *p,f;
  MtHist3D *h=NULL;

  USE(rcThis);

  ALLOCOBJ0(h);

  h->n = 0;
  h->nb_d = d = nb_1d;
  h->nb_dd = nb_dd = d*d;
  h->nb = nb = d*d*d;

  ALLOCARR(h->H,nb);
  H = h->H;
  for(i=0;i<nb;i++) H[i]=0;

  n=0;
  for(i=0;i<X->nRow;i++)
  {
    p = X->rows[i];
    for(k=0;k<3;k++) 
    {
      f = p[k];
      UT_ASSERT(!(f<0||f>1));
      j = (int)(d*f); if(j<0)j=0; if(j>=d)j=d-1;
      jj[k] = j;
    }

    j = (&MT_HIST_3D(h,jj[0],jj[1],jj[2])) - H;
//    BMP_RGB_HISTO(h,cr,cg,cb)++;
    UT_ASSERT(j<nb);
    H[j]++;
    n++;
  }
  h->n=n;

  // Normalize
  for(i=0;i<nb;i++) H[i] /= (double)n;

rcCatch:
  return h;
}



/*=====================================================================
 *
 *  MtDistr
 *
 * ===================================================================*/

#define MT_SORT_DISTR( ds ) \
  if(!(ds)->isSorted) \
  { \
    qsort( (ds)->smp, (ds)->nSmp, sizeof( double ), MtDoubleCmpFunc ); \
    (ds)->isSorted = TRUE; \
  }  

/*-------------------------------------------------------------------------*/
void MtDelDistr( MtDistr *ds )
/*-------------------------------------------------------------------------*/
{
  if(ds)
  {
    MtDelHist(ds->hist);
    FREEMEM(ds->smp);
    FREEMEM(ds->dens);
    FREEMEM(ds);
  }
}

/*-------------------------------------------------------------------------*/
MtDistr *MtCreateDistr( double *samples,	// samples
    			int 	n,		// number of samples		
			unsigned options
			)
/*-------------------------------------------------------------------------*/
{
  int i;
  MtDistr *ds=NULL;
  double  x;

  ALLOCOBJ0(ds);
  ALLOCARR(ds->smp,n);

  ds->nSmp = n;
  ds->isSorted = FALSE;
  ds->kbwOpt0 = -1;

  MT_STAT_INIT( &ds->stat );
  for( i=0; i < n; i++ )
  {
    x=samples[i];
    ds->smp[i]=x;
    MT_STAT_UPDATE( &ds->stat, x );
  }
  MT_STAT_RESULT( &ds->stat );  

  return ds;
}

/*-------------------------------------------------------------------------*/
MtHist *MtGetDistrHist( MtDistr *ds,		// distribution 
    		    	int	 numBins 	// number of bins
		    	)
/*-------------------------------------------------------------------------*/
{
  MtHist *h;
  int	rc,rcThis=0;
  
  if( !ds->hist || ds->hist->nb != numBins )
  {
    MtDelHist(ds->hist);
    ds->hist = NULL;
    h = MtCreateHist(numBins,ds->nSmp,FALSE,
		       ds->stat.min,ds->stat.max);
    if(!h) raiseRc(MI_ERROR);
    ds->hist = h;
    rc = MtFillHist( h, ds->smp, ds->nSmp, FALSE );
    if(rc) raiseRc(MI_EMATH);
  }

rcCatch:
return ds->hist;
}

/*-------------------------------------------------------------------------*/
int	MtMakeDistrDensity( 
    		MtDistr *ds,  	// distribution 
		double   xMin,  // minimum value
		double   xMax,  // maximum value
    		int      n,	// number of interpolation points
    		double   bw	// kernel bandwidth times default bw
		)
/*-------------------------------------------------------------------------*/
{
  int	i,j,m,rcThis=0;
  double dx,x,sum,*smp,h,u;

  FREEMEM(ds->dens);
  ALLOCARR(ds->dens,n);
  ds->nDens = n;

  m = ds->nSmp;
  smp = ds->smp;
  dx = (xMax-xMin)/(double)(n-1);

  // Silverman's rule-of-thumb kernel bandwidth  
  h = 1.06 * ds->stat.var / pow(m,1./5);
  h = bw*h;

  for(i=0;i<n;i++)
  {
    x = xMin+i*dx;
    for(sum=0,j=0;j<m;j++) 
    {
      u = (x-smp[j])/h;
      sum += exp(-u*u*0.5);
    }
    ds->dens[i] = sum/((double)m * h * MT_SQRT2PI);
  }
  
/*rcCatch:*/
return rcThis;
}


/*-------------------------------------------------------------------------*/
void MtSortDistr( MtDistr *ds )
/*-------------------------------------------------------------------------*/
{
  MT_SORT_DISTR(ds);
}

/*-------------------------------------------------------------------------*/
double MtDistrInvCDF( MtDistr *ds,	// distribution
    		      double   u,	// prob. value: 0 < u < 1 
		      double   xMin,	// minimum possible sample value
		      double   xMax,	// maximum possible sample value
		      unsigned optIn
		    	)
/*-------------------------------------------------------------------------*/
{
  double x,y,*smp=ds->smp,alphaL,alphaR,beta,x1,x2;
  int j,n=ds->nSmp,nQ;
  unsigned opt;

  opt = optIn ? optIn : MT_LINEAR;

  if(opt & MT_RAW)
  {
    /*** no interpolation: simple resampling */
    j = (int)(u*(n/*+1*/));
    if(j>=n) j=n-1; 
    if(j<0) j=0;
    return smp[j];
  }

  if(opt & MT_KDE)
  {
    double sigKrn,bw,sig;
    /*** no interpolation: simple resampling */
    j = (int)(u*n);
    if(j>=n) j=n-1; 
    if(j<0) j=0;
    if(ds->kbwOpt0 < 0)
    {
      double qrange;
      MT_SORT_DISTR(ds);
      qrange = smp[(int)(0.75*n)]-smp[(int)(0.25*n)];
      sig = ds->stat.var;
      bw = 0.776 * 1.364 * MIN( sig, qrange/1.34 ) *
	   pow( n, -0.2 );
      ds->kbwOpt0 = bw;
      // printf("sig=%g qr=%g bw=%g\n",sig,qrange,bw);
    }
    bw = ds->kbwOpt0; sigKrn = 1.0;
    x = smp[j] + ds->kbwOpt0 * MtRndNormGet(0,sigKrn);
    return x;
  }

  y = u*(n+1);

  beta = 1; /* "fat-tail" factor 1=normal <1:fatter */

  MT_SORT_DISTR(ds);

  j=(int)y;
  if(j<0||j>n)
  {
    TRCERR(("index %d out of range 0..%d y=%g u=%g)\n",j,n,y,u));
    if(j<0) j = 0; 
    if(j>n) j=n;
  }

  /*** left & right tail slopes */
  nQ = MAX(n*0.25,2);
  alphaL = (smp[nQ-1]-smp[0])/(double)(nQ-1);
  alphaR = (smp[n-1]-smp[n-nQ])/(double)(nQ-1);
  
  if(fabs(xMax) < MT_REAL_EPSILON && 
      fabs(xMin) < MT_REAL_EPSILON)
  {
    /*** estimate xMax, u. xMin */
    xMin = smp[0] - alphaL;
    xMax = smp[n-1] + alphaR;
  }

  if( opt & MT_LINEAR )
  {
    /*** linear tail interpolation */
    if(j<1) 
    {
      x1 = xMin; x2 = smp[j];
    }
    else if(j>=n)
    {
      x1 = smp[j-1]; x2 = xMax;
    }
    else 
    {
      x1 = smp[j-1]; x2 = smp[j];
    }
    x = x1 + (x2-x1)*(y-j);
    return x;
  }
  
  if(j<1)
  {
    /*** left tail: exponential interpolation */
    x = smp[0] + alphaL*beta*log(y);
  }
  else if(j>=n)
  {
    /*** right tail */
    x = smp[n-1] - alphaR*beta*log(1-(y-n));
  }
  else
  {
    /*** middle part: linear interpolation */
    x = smp[j-1] + ( smp[j] - smp[j-1] ) * (y - j);
  }

  return x;
}

/*-------------------------------------------------------------------------*/
/*  MtCompareDistr: compares empirical distributions f1,f2		   */
/*  NOTE: f1, f2 assumed to be sorted arrays				   */		    
/*-------------------------------------------------------------------------*/
int MtICompareDistr( double * f1, int n1,   
    	            double * f2, int n2,
		    unsigned	opt,
		    double	*difCramer,
		    double	*difKomolSmir )
{
#define CVM_ABS

  double df,sum1=0,sum2=0,diff=0, df1=1./(double)n1, df2=1./(double)n2,x,xOld,
  	dfMax=0,dx;
  int	i1=0,i2=0,usef1,rcThis=0;

  xOld = 0;
  while(i1<n1||i2<n2)
  {
    if(i1==n1)
      usef1=0;
    else if(i2==n2)
      usef1=1;
    else if( f1[i1] < f2[i2] )
      usef1=1;
    else
      usef1=0;

    if(usef1)
    {
      x = f1[i1++];
      dx = (x-xOld); df = (sum1-sum2);
      sum1+=df1;
    }
    else
    {
      x = f2[i2++];
      dx = (x-xOld); df = (sum1-sum2);
      sum2+=df2;
    }

/*    diff += df*df;*/

#ifdef CVM_ABS    
    diff += fabs(df*dx);
#else
    diff += df*dx * df*dx;
#endif

    if( dx > MT_REAL_EPSILON && fabs(df) > dfMax ) dfMax = fabs(df);
    
/*    printf("x=%f s1=%f s2=%f dx=%f df=%f diff=%f\n",x,sum1,sum2,x-xOld,df,diff);*/

    xOld = x;
  }

/*rcCatch:*/

#ifdef CVM_ABS    
  *difCramer = diff;
#else
  *difCramer = sqrt(diff);
#endif

  *difKomolSmir = dfMax;
#if 0
  *difKomolSmir = sqrt( ((double)(n1)*(double)(n2)) /
      			 ((double)(n1)+(double)(n2)) ) * dfMax;
#endif
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtDistrCompareCVM( 
    	MtDistr *ds1, 		// IN: first distribution
        MtDistr *ds2, 		// IN: second distribution
        unsigned opt,		// IN: options
	double  *pDisCVM,	// OUT: Cramer v. Mises Distance	
	double  *pDisKS		// OUT: Komologrov Smirnov 
	)
{
  int rcThis=0,rc;

  MT_SORT_DISTR(ds1); MT_SORT_DISTR(ds2);
  rc = MtICompareDistr( ds1->smp, ds1->nSmp, 
                   	ds2->smp, ds2->nSmp,
		   	0,
		   	pDisCVM, pDisKS );
/*rcCatch:*/
  return rcThis;
}

#if 0
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtDistrCompareNORM( 
    	MtDistr *dst, 		// IN: empirical distribution
	double   nrmMu,		
	double   nrmSig,
        unsigned opt,		// IN: options
	double  *pDisCVM,	// OUT: Cramer v. Mises Distance	
	double  *pDisKS		// OUT: Komologrov Smirnov 
	)
{
  double df,sum1=0,sum2=0,diff=0, df1, df2,
  	dfMax=0,dx,*f1,f2,x,xmin,xmax;
  int	i1=0,i2=0,usef1,rcThis=0,n1,n2;

  n1 = dst->nSmp;
  df1=1./(double)n1;
  f1 = dst->smp;

  n2 = MAX(n1,10000000);
  df2=1./(double)n2;
  f2 = 0;

  xmin = dst->smp[0];
  xmax = dst->smp[n1-1];
  dx2 = (xmin-xmax)/(double)(n2);

  xOld = 0;
  while(i1<n1||i2<n2)
  {
    if(i1==n1)
      usef1=0;
    else if(i2==n2)
      usef1=1;
    else if( f1[i1] < f2 )
      usef1=1;
    else
      usef1=0;

    if(usef1)
    {
      x = f1[i1++];
      dx = (x-xOld); df = (sum1-sum2);
      sum1+=df1;
    }
    else
    {
      x = f2;
      i2++;
      f2 += MtNormalDens(x,nrmMu,nrmSig);
      dx = (x-xOld); df = (sum1-sum2);
      sum2+=df2;
    }

/*    diff += df*df;*/

#ifdef CVM_ABS    
    diff += fabs(df*dx);
#else
    diff += df*dx * df*dx;
#endif

    if( dx > MT_REAL_EPSILON && fabs(df) > dfMax ) dfMax = fabs(df);
    
/*    printf("x=%f s1=%f s2=%f dx=%f df=%f diff=%f\n",x,sum1,sum2,x-xOld,df,diff);*/

    xOld = x;
  }

/*rcCatch:*/

#ifdef CVM_ABS    
  *difCramer = diff;
#else
  *difCramer = sqrt(diff);
#endif

  *difKomolSmir = dfMax;
#if 0
  *difKomolSmir = sqrt( ((double)(n1)*(double)(n2)) /
      			 ((double)(n1)+(double)(n2)) ) * dfMax;
#endif
  return rcThis;
}
#endif

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtTstDistrInvCDF_old( double *vSmpRef, // reference distribution
    		      int     nSmpRef,
		      double  subSmpFrac,  // fractions sub sample
		      int     rndSeed,
		      unsigned optCDF )
{
  int nSub,i,rcThis=0,m,rowNum,rc,quNum,k,ni,niMax;
  MtDistr *dstRef=NULL,*dst=NULL,*dstSub=NULL;
  double *vSmp=NULL,*vSub=NULL;
  double u,xMax,xMin,q,p,dq=0.1,disKS,disCVM;
  double qu[]={0,0.02,0.25,0.5,0.75,1},qmax,eps;
  
  dstRef = MtCreateDistr( vSmpRef, nSmpRef, 0 );
  if(!dstRef) raiseRc( 1 );

  if(rndSeed >=0) MtRndSeed( rndSeed );

  /*** take a sub sample of the reference distribution */
  MtVecShuffle(vSmpRef,nSmpRef,-1);
  nSub = MAX( nSmpRef * subSmpFrac, 5 );
  printf(
      "===============> ref. samples=%d, sub-samples = %d\n",
      nSmpRef,nSub);
  FREEMEM(vSub); ALLOCARR(vSub,nSub);
  for(i=0;i<nSub;i++)
  {
    vSub[i] = vSmpRef[i];
  }
  dstSub = MtCreateDistr(vSub,nSub,0);
  if(!dstSub) raiseRc(1);

  /*** take Monte-Carlo sample */
  quNum = ARRAY_LENGTH(qu);
  dq = 0.2;
  printf("       "); 
  for(k=0;k<quNum;k++) printf("%8.4f",qu[k]); 
  printf("\n");
  
  printf("       "); 
  for (k=0;k<quNum;k++) printf("--------"); 
  printf("\n");

  printf("       ");
  eps = MT_NUM_EPS;
  qmax = (double)((double)1.0-(double)eps);
  for (k=0;k<quNum;k++) 
  {
    q = qu[k]; 
    p = MtDistrInvCDF(dstRef,MIN(q,qmax),0,0,optCDF);
    printf("%8.4f",p); 
  }
  printf("  %7.3f%7.3f%7.3f%7.3f\n",
      dstRef->stat.min,dstRef->stat.max,dstRef->stat.avg,dstRef->stat.var);
  for (k=0;k<quNum;k++) printf("--------"); 
  printf("\n");

  niMax = 3; ni=0;
  for(m=10, rowNum=0;m<=nSmpRef;m = MIN(10*m,nSmpRef), rowNum++)
  {
    xMin = 0; xMax = 0;
    FREEMEM(vSmp); ALLOCARR(vSmp,m);
    for(i=0;i<m;i++)
    {
      u = MtRndGet();
      vSmp[i] = MtDistrInvCDF(dstSub,u,xMin,xMax,optCDF);
    }
    MtDelDistr(dst); dst=NULL;
    dst = MtCreateDistr(vSmp,m,0);
    if(!dst) raiseRc(1);

    /**** compare quantiles */
    printf("%7d",m);
    for (k=0;k<quNum;k++)
    {
      p = MtDistrInvCDF( dst, MIN(qu[k],qmax), 0,0, optCDF );
      printf("%8.4f", p);
    }
    /**** distribution divergence */
    rc = MtDistrCompareCVM( dst, dstRef, 0, &disCVM, &disKS );
    if(rc) raiseRc(rc);

    printf("  %7.3f%7.3f%7.3f%7.3f",
	dst->stat.min,dst->stat.max,dst->stat.avg,dst->stat.var);

    printf(" -> %7.6f %7.6f",disCVM, disKS );

    printf("\n");
    if(m==nSmpRef)
    {
      if( ++ni > niMax) break;
    }

  }


rcCatch:
  MtDelDistr(dstRef);
  MtDelDistr(dstSub);
  MtDelDistr(dst);
  FREEMEM(vSmp);
  FREEMEM(vSub);
  return rcThis;
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MtTstDistrInvCDF( double *vSmpRef, // reference distribution
    		      int     nSmpRef,
		      double  subSmpFrac,  // fractions sub sample
		      int     nSmpMC,
		      int     nIterMC,	// number of monte carlo rounds
		      int     rndSeed,
		      unsigned optCDF )
{
#define RELERR(x,x0,sd) (fabs((x)-(x0))/(sd))
  int nSub,i,rcThis=0,m,rc,quNum,k,iMC,iMC2,j,sd_n;
  MtDistr *dstRef=NULL,*dst=NULL,*dstSub=NULL;
  double *vSmp=NULL,*vSub=NULL,sdref=0;
  double u,xMax,xMin,q,p,dq=0.1,disKS,disCVM,pRef;
  double qu[]={0,0.02,0.25,0.5,0.75,1},qmax,qmin,eps;

  double s_q[100],sd_q[100],s_q_ref[100],
	 sd_max,sd_min,sd_mu,sd_sig,sd_CVM,sd_KS,
	 s_min,s_max,s_mu,s_sig;
    quNum = ARRAY_LENGTH(qu);

    dq = 0.2;
    eps = MT_NUM_EPS;
    qmax = (double)((double)1.0-(double)eps);
    qmin=eps;
  
  dstRef = MtCreateDistr( vSmpRef, nSmpRef, 0 );
  if(!dstRef) raiseRc( 1 );

  sdref = dstRef->stat.var;

  for (k=0;k<quNum;k++)
  {	
    q = MIN(qu[k],qmax); if(q<eps) q=eps;
    p = MtDistrInvCDF( dstRef, q, 0,0, MT_LINEAR );
    s_q_ref[k] = p;
  }


  if(rndSeed >=0) MtRndSeed( rndSeed );

  nSub = MAX( nSmpRef * subSmpFrac, 5 );
  printf(
      "===============> ref. samples=%d, sub-samples = %d\n",
      nSmpRef,nSub);
  FREEMEM(vSub); ALLOCARR(vSub,nSub);

  sd_max = sd_min = sd_mu = sd_sig = sd_CVM = sd_KS = 
    s_max = s_min = s_mu = s_sig =0;
  sd_n=0;
  for(k=0;k<ARRAY_LENGTH(qu);k++) s_q[k] = sd_q[k] = 0;

  /*** take a sub sample of the reference distribution */
  for(iMC=0;iMC<nIterMC;iMC++)
  {
    //MtVecShuffle(vSmpRef,nSmpRef,-1);
    for(i=0;i<nSub;i++)
    {
      u = MtRndGet();
      j = u*nSmpRef; if(j<0) j=0; if(j>=nSmpRef) j=nSmpRef-1;
      vSub[i] = vSmpRef[j];
    }
    MtDelDistr(dstSub); dstSub=NULL;
    dstSub = MtCreateDistr(vSub,nSub,0);
    if(!dstSub) raiseRc(1);

    /*** take Monte-Carlo sample */

    m = nSmpMC;
    for(iMC2=0;iMC2<nIterMC;iMC2++)
    {
      xMin = 0; xMax = 0;
      FREEMEM(vSmp); ALLOCARR(vSmp,m);
      for(i=0;i<m;i++)
      {
	u = MtRndGet();
	vSmp[i] = MtDistrInvCDF(dstSub,u,xMin,xMax,optCDF);
      }
      MtDelDistr(dst); dst=NULL;
      dst = MtCreateDistr(vSmp,m,0);
      if(!dst) raiseRc(1);

      /**** compare quantiles */
      for (k=0;k<quNum;k++)
      {	
	q = MIN(qu[k],qmax); if(q<eps) q=eps;
	pRef = s_q_ref[k]; 
	p = MtDistrInvCDF( dst, q, 0,0, optCDF );
	s_q[k] += p;
	sd_q[k] += RELERR(p,pRef,sdref);
      }

      /**** distribution divergence */
      rc = MtDistrCompareCVM( dst, dstRef, 0, &disCVM, &disKS );
      sd_CVM += disCVM;
      sd_KS += disKS;

      s_min += dst->stat.min;
      s_max += dst->stat.max;
      s_mu += dst->stat.avg;
      s_sig += dst->stat.var;

      sd_min += RELERR( dst->stat.min, dstRef->stat.min, sdref );
      sd_max += RELERR( dst->stat.max, dstRef->stat.max, sdref );
      sd_mu += RELERR( dst->stat.avg, dstRef->stat.avg, sdref );
      sd_sig += RELERR( dst->stat.var, dstRef->stat.var, sdref );
      sd_n++;
    }
  }

  s_min /= (double)sd_n;
  s_max /= (double)sd_n;
  s_mu /= (double)sd_n;
  s_sig /= (double)sd_n;
  sd_min /= (double)sd_n;
  sd_max /= (double)sd_n;
  sd_mu /= (double)sd_n;
  sd_sig /= (double)sd_n;
  sd_CVM /= (double)sd_n;
  sd_KS /= (double)sd_n;
  for (k=0;k<quNum;k++) 
  {
    sd_q[k] /= (double)sd_n;
    s_q[k] /= (double)sd_n;
  }

    for (k=0;k<quNum;k++) printf("%8.4f",qu[k]);
    printf("\n");

    for (k=0;k<quNum;k++) printf("%8.4f",s_q_ref[k]);
    printf("   %7.3f%7.3f%7.3f%7.3f",
	dstRef->stat.min,dstRef->stat.max,
	dstRef->stat.avg,dstRef->stat.var);
    printf("\n");

    for (k=0;k<quNum;k++) printf("%8.4f",s_q[k]);
    printf("   %7.3f%7.3f%7.3f%7.3f",s_min,s_max,s_mu,s_sig);
    printf("\n");
   
    for (k=0;k<quNum;k++) printf("%8.4f",sd_q[k]);
    printf("   %7.3f%7.3f%7.3f%7.3f",sd_min,sd_max,sd_mu,sd_sig);
    printf(" -> %7.6f %7.6f",sd_CVM, sd_KS );
    printf("\n");


rcCatch:
  MtDelDistr(dstRef);
  MtDelDistr(dstSub);
  MtDelDistr(dst);
  FREEMEM(vSmp);
  FREEMEM(vSub);
  return rcThis;
}

//
double MtRelDiff( double a, double b )
{
  double d;
  a = fabs(a); b=fabs(b);
  if(a<MT_REAL_EPSILON) d = b;
  else if(b<MT_REAL_EPSILON) d = a;
  else d = ((a) > (b) ? ((a)-(b))/(b) :  ((b)-(a))/(a));
  return d;
}

/*=========================================================================
 *
 *
 *   Complex numbers      
 *
 *
 * =======================================================================*/


/*
 * Zero out relatively very small real or imaginary parts of a complex number,
 * because they probably are a result of accumulated floating point inaccuracies.
 *
 * Return true if something was zeroed out.
 */
int cmplx_fixup(cmplx *ap)
{
#define cmplx_epsilon	0.00000000000005
	if (fabs(ap->re * cmplx_epsilon) > fabs(ap->im)) {
		ap->im = 0.0;
		return TRUE;
	}
	if (fabs(ap->im * cmplx_epsilon) > fabs(ap->re)) {
		ap->re = 0.0;
		return TRUE;
	}
	return FALSE;
}

/*
 * Add two complex numbers (a + b) and return result.
 *
 * Complex number subtraction (a - b) is done by
 * cmplx_add(a, cmplx_negate(b)).
 */
cmplx cmplx_add(cmplx a, cmplx b)
{
	a.re += b.re;
	a.im += b.im;
	return(a);
}

// Negate a complex number (-a) and return result.
cmplx cmplx_negate(cmplx a)
{
	a.re = -a.re;
	a.im = -a.im;
	return(a);
}

// Multiply two complex numbers (a * b) and return result.
cmplx cmplx_mult(cmplx a, cmplx b)
{
	cmplx	r;

	r.re = a.re * b.re - a.im * b.im;
	r.im = a.re * b.im + a.im * b.re;
	return(r);
}

// Divide two complex numbers (a / b) and return result.
cmplx cmplx_div(cmplx a, cmplx b)
{
	cmplx	r, num;
	double		denom;

	b.im = -b.im;
	num = cmplx_mult(a, b);
	denom = b.re * b.re + b.im * b.im;
	r.re = num.re / denom;
	r.im = num.im / denom;
	return r;
}

// Take the natural logarithm of a complex number and return result.
cmplx cmplx_log(cmplx a)
{
	cmplx	r;

	r.re = log(a.re * a.re + a.im * a.im) / 2.0;
	r.im = atan2(a.im, a.re);
	return(r);
}

// Raise the natural number to the power of a complex number (e ^ a)
// and return result.
cmplx cmplx_exp(cmplx a)
{
	cmplx	r;
	double		m;

	m = exp(a.re);
	r.re = m * cos(a.im);
	r.im = m * sin(a.im);
	return(r);
}

// Raise complex number "a" to the power of complex number "b" (a ^ b)
// and return result.
cmplx cmplx_pow(cmplx a, cmplx b)
{
	cmplx	r;

	r = cmplx_log(a);
	r = cmplx_mult(r, b);
	r = cmplx_exp(r);
	cmplx_fixup(&r);
	return(r);
}


