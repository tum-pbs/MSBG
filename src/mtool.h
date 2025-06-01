/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#ifndef MTOOL_H
#define MTOOL_H

#include "util.h" 

MSBG_NAMESPACE_BEGIN

typedef	float	  MtReal;

typedef MtReal   *MtVector;

typedef MtReal   *MtMatrix; 

#define MT_REAL_MAX	FLT_MAX;
/*#define MT_REAL_MIN	FLT_MIN;*/
#define MT_REAL_MIN	-1e35   /* TODO */
#define MT_NUM_EPS	(1.0E-12)
#define MT_NUM_EPS_F	(1.0E-12f)
#define MT_DBL_EPS	(1e-15)	
#define MT_REAL_EPSILON   (1.0E-20)
#define MT_EPSILON MT_REAL_EPSILON
#define MT_MA_MAX 8192
#define MT_PI  (3.14159265358979323846)
#define MT_2PI  (6.28318530717958647692)
#define MT_SQRT2PI (2.50662827463100050241)
#define MT_SPHERE_VOLUME(r) ((4.*MT_PI*(r)*(r)*(r))/3.)
#define MT_SQRT2 1.4142135623
#define MT_SQRT3 1.73205080756887
#define MT_SQRT3F 1.73205080756887f
#define MT_SQRT2INV .707106781186
#define MT_RND_STAB_SZ 	32
#define MT_LOG2	.69314718055994530941

#define MT_OBS_CALL 	1
#define MT_OBS_PUT 	2

#define MT_INFINITY	 (1e20)
#define MT_MINUS_INFINITY (-1e20)

#define MT_INIT		(1<<1)
#define MT_NO_COPY	(1<<4)
#define MT_QUIETERR     (1<<5)
#define MT_LINEAR	(1<<6)
#define MT_RAW		(1<<7)
#define MT_KDE		(1<<8)
#define MT_EXP		(1<<9)
#define MT_MINMAX	(1<<10)

#define MT_SELECT	(1<<11)
#define MT_MEAN		(1<<12)
#define MT_SIG		(1<<13)
#define MT_COV		(1<<14)
#define MT_CORR		(1<<15)
#define MT_NONORM	(1<<16)
#define MT_LOGDENS	(1<<17)
#define MT_NORMSIG	(1<<18)
#define MT_NOCORR	(1<<19)
#define MT_CHECK	(1<<20)
#define MT_VERBOSE	(1<<21)
#define MT_NR3		(1<<22)
#define MT_SORTED	(1<<23)
#define MT_ICOV		(1<<24)
#define MT_NOCOV	(1<<25)

#define MT_DC_CVM		1
#define MT_DC_KS		2

/**** error return codes */

#define MT_OK  		0
#define MT_ERROR 	1
#define MT_ENOMEM 	2
#define MT_EDIVZERO 	3
#define MT_SINGULAR	4
  
#ifndef MAX
  #define MAX(a,b)    ((a)>(b)?(a):(b))
#endif
#ifndef MIN
  #define MIN(a,b)    ((a)<(b)?(a):(b))
#endif
#ifndef ABS
  #define ABS(a)      (((a)>0)?(a):(-(a)))
#endif  

#define MT_MAX3(a,b,c) MAX(a,MAX(b,c))

#define MT_SIGNUM(x) (((x) > 0) ? 1 : (((x) < 0) ? -1 : 0))

#define MT_SIGN(x)  ((x) >= 0  ? 1 : -1)

#define MT_SIGNUM_EPS(x,eps) (((x) > (eps)) ? 1 : (((x) < -(eps)) ? -1 : 0))

#define MT_POSPART( a ) ((a) > 0 ? (a) : 0)
#define MT_NEGPART( a ) ((a) < 0 ? -(a) : 0)

#define MT_CLAMP_MIN(x,a) \
{\
  if(unlikely((x)<(a))) x=a;	\
}
#define MT_CLAMP_MAX(x,a) \
{\
  if(unlikely((x)>(a))) x=a;	\
}

#define MT_CLAMP(x,xmin,xmax) \
{\
  MT_CLAMP_MIN(x,(xmin));\
  MT_CLAMP_MAX(x,(xmax));\
}

#define MT_CLIP_REFL_1D(x1,x2, /* IN: min,max */\
    		      xc  /* OUT: coordante to clip */\
    		      )\
{\
  if((xc)<(x1)) xc = (x1) + ((x1)-(xc)); \
  if((xc)>(x2)) xc = (x2) - ((xc)-(x2)); \
  if(xc<x1) xc=x1;\
  if(xc>x2) xc=x2;\
}

// Smooth transition to small value 
//   - delt = thershold for beginning of nonlinearity, anything above delt stays unchanged) 
//   - esp = final small value for x=0
// http://www.iquilezles.org/www/articles/functions/functions.htm
#define MT_SMOOTH_TO_EPS( delt , eps, \
    			  x /* IN/OUT */ )\
{	\
  if( x<(delt) )\
  {\
    double a = 2.0f*(eps) - (delt),\
	   b = 2.0f*(delt) - 3.0f*(eps),\
	   t = (x)/(delt);\
    x = (a*t + b)*t*t + (eps);\
  }\
}

#ifdef __cplusplus

namespace MT
{
  inline double sqr(double x) { return (x * (double)x); };
  inline float sqrf(float x) { return (x * x); };

  inline float absmin2 (float x, float y) {
    if (fabsf(x) < fabsf(y)) return x;
    return y;
  }

  inline double absmin2 (double x, double y) {
    if (fabs(x) < fabs(y)) return x;
    return y;
  }

} // namespace MT

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
inline
void MtTrlinearInterpolationWeights( float u, float v, float w,  // fractions
    				     float *weights8 // OUT: weights (scan-line order)
    )
{
  float *p=weights8,
	u1=1.f-u,
	v1=1.f-v,
	w1=1.f-w;

  p[0] = u1*v1*w1;
  p[1] = u*v1*w1; 
  p[2] = u1*v*w1;  
  p[3] = u*v*w1; 
  p[4] = u1*v1*w;  
  p[5] = u*v1*w;  
  p[6] = u1*v*w;  
  p[7] = u*v*w;
}



template<typename T, typename T_intermediate> inline 
void MtTrilinearInterpolation( 
    	T *g, // The eight interpolation cube corners in lex. order		
	float u, float v, float w, // Fractional coords in x,y, and z direction

	T *outVal,  // OUT: interpolated value
	T *outGrad	// optional OUT: gradient components
    )
{
  /*
  Maxima:
  f(u,v,w) :=
	  ( (   (1-u) * f0 
	  +      u * f1  ) * (1-v)
	  + (   (1-u) * f2 
	  +      u * f3  ) * v ) * (1-w)
	  + ( ( (1-u) * f4 
	  +      u * f5  ) * (1-v)
	  + (   (1-u) * f6 
	  +      u * f7  ) * v ) * w ;

  string(ratsimp(diff(f(u,v,w),u)));
  string(ratsimp(diff(f(u,v,w),v)));
  string(ratsimp(diff(f(u,v,w),w)));

  */

  T_intermediate 
    	 f0=g[0], f1=g[1], f2=g[2], f3=g[3],
  	 f4=g[4], f5=g[5], f6=g[6], f7=g[7];

  if(outVal) *outVal = 
	  ( (   (1-u) * f0 
	  +      u * f1  ) * (1-v)
	  + (   (1-u) * f2 
	  +      u * f3  ) * v ) * (1-w)
	  + ( ( (1-u) * f4 
	  +      u * f5  ) * (1-v)
	  + (   (1-u) * f6 
	  +      u * f7  ) * v ) * w ;

  if(outGrad)
  {
    T_intermediate f01234567 = f7 - f6 - f5 + f4 - 
		       f3 + f2 + f1 - f0,
	   f01234567u = f01234567*u,
	   f01234567uf6420 = f01234567u + f6 - f4 - f2+ f0,
	   f5410 = f5 - f4 - f1 + f0,
	   f3210 = f3 - f2 - f1 + f0;

    outGrad[0] = (f01234567 * v + f5410)*w + f3210*v + f1-f0;
    outGrad[1] = f01234567uf6420 * w + f3210*u + f2-f0;
    outGrad[2] = f01234567uf6420 * v + f5410*u + f4-f0;
  }
}




static inline 
void MtTrilinearInterpolationOld( 
    	float *g, // The eight interpolation cube corners in lex. order		
	float u, float v, float w, // Fractional coords in x,y, and z direction

	float *outVal,  // OUT: interpolated value
	float *outGrad	// optional OUT: gradient components
    )
{
  /*
  Maxima:
  f(u,v,w) :=
	  ( (   (1-u) * f0 
	  +      u * f1  ) * (1-v)
	  + (   (1-u) * f2 
	  +      u * f3  ) * v ) * (1-w)
	  + ( ( (1-u) * f4 
	  +      u * f5  ) * (1-v)
	  + (   (1-u) * f6 
	  +      u * f7  ) * v ) * w ;

  string(ratsimp(diff(f(u,v,w),u)));
  string(ratsimp(diff(f(u,v,w),v)));
  string(ratsimp(diff(f(u,v,w),w)));

  */

  double f0=g[0], f1=g[1], f2=g[2], f3=g[3],
  	 f4=g[4], f5=g[5], f6=g[6], f7=g[7];

  if(outVal) *outVal = 
	  ( (   (1-u) * f0 
	  +      u * f1  ) * (1-v)
	  + (   (1-u) * f2 
	  +      u * f3  ) * v ) * (1-w)
	  + ( ( (1-u) * f4 
	  +      u * f5  ) * (1-v)
	  + (   (1-u) * f6 
	  +      u * f7  ) * v ) * w ;

  if(outGrad)
  {
    double f01234567 = f7 - f6 - f5 + f4 - 
		       f3 + f2 + f1 - f0,
	   f01234567u = f01234567*u,
	   f01234567uf6420 = f01234567u + f6 - f4 - f2+ f0,
	   f5410 = f5 - f4 - f1 + f0,
	   f3210 = f3 - f2 - f1 + f0;

    outGrad[0] = (f01234567 * v + f5410)*w + f3210*v + f1-f0;
    outGrad[1] = f01234567uf6420 * w + f3210*u + f2-f0;
    outGrad[2] = f01234567uf6420 * v + f5410*u + f4-f0;
  }
}

// Quadratic (approximating) uniform B-spline weights 
//
//  Evaluation at position x 
//
//  f = w0 * f[i0] + 
//	w1 * f[i0+1] +
//	w2 * f[i0+2]
//
//  where 
//    i0 = floor(x-0.5-o))
//    t = 1-((x-0.5-o)-floor(x-0.5-o))
//    o = 0 for node centered, 0.5 for cell centered grid sampling 

static inline
void MtQuadBSplineWeights1D( float t, // In: fractional coordinate
    			     float *w, // OUT: weights w[0],...w[3]
			     float *w_  // OUT: derivative weights 
    )    				
{
    // B-spline
    w[0] = t * t / 2.f;
    w[1] = .5f + t * (1.f - t);
    w[2] = (1.f - t) * (1.f - t) / 2.f;

#if 0
    // partial derivative
    w_[0] = t;
    w_[1] = (1.f-2.f*t);
    w_[2] = (t-1.f);
#else
    w_[0] = -t;
    w_[1] = -(1.f-2.f*t);
    w_[2] = -(t-1.f);

#endif
}

static inline
void MtCubicBSplineWeights1D( float t, // In: fractional coordinate
    			     float *w, // OUT: weights w[0],...w[3]
			     float *w_  // OUT: derivative weights 
    )    				
{
  //	cf. Section (4) in 
  //	Ruijters: "Efficient GPU-Based Texture Interpolation using Uniform B-Splines"  
  //
  w[0] = (1.f/6.f)*(1.f-t)*(1.f-t)*(1.f-t);
  w[1] = (2.f/3.f) - 0.5f*t*t*(2.f-t);
  w[2] = (2.f/3.f) - 0.5f*(1.f-t)*(1.f-t)*(1.f+t);
  w[3] = (1.f/6.f)*t*t*t;

  // partial derivative
  w_[0] = -(1.f-t)*(1.f-t)/2.f;
  w_[1] = 0.5f*t*t-1.0f*(2.f-t)*t;
  w_[2] = 1.f*(1.f-t)*(t+1.f)-0.5f*(1.f-t)*(1.f-t);
  w_[3] = t*t/2.f;
}


static inline
void MtCubicBSplineWeightsForSecondDerivs1D( float t, // In: fractional coordinate
    			     float *w, // OUT: weights w[0],...w[3]
			     float *w_,  // OUT: derivative weights 
			     float *w__  // OUT: second derivative weights 
    )    				
{
  w[0] = (1.f/6.f)*(1.f-t)*(1.f-t)*(1.f-t);
  w[1] = (2.f/3.f) - 0.5f*t*t*(2.f-t);
  w[2] = (2.f/3.f) - 0.5f*(1.f-t)*(1.f-t)*(1.f+t);
  w[3] = (1.f/6.f)*t*t*t;

  // 1st derivatives
  w_[0] = -(1.f-t)*(1.f-t)/2.f;
  w_[1] = 0.5f*t*t-1.0f*(2.f-t)*t;
  w_[2] = 1.f*(1.f-t)*(t+1.f)-0.5f*(1.f-t)*(1.f-t);
  w_[3] = t*t/2.f;

  //  2nd derivatives
  w__[0] = 1.f - t;
  w__[1] = 3.f*t - 2.f;
  w__[2] = 1.f - 3.f*t;
  w__[3] = t;
}

inline float MT_SMOOTHSTEPF(float a, float b, float x)
{
    if (x < a)
        return 0.0f;
    if (x >= b)
        return 1.0f;
    x = (x - a)/(b - a); /* normalize to [0:1] */
    return (x*x * (3.0f - 2.0f*x));
}

inline float MT_LINSTEPF(const float a, const float b, float x)
{
  //UT_ASSERT2( b>a );
    if (x < a)
        return 0;
    if (x >= b)
        return 1;
    x = (x - a)/(b - a); /* normalize to [0:1] */
    return x;
}

inline float MT_LINFIT( float val, 
    			const float *r // mapping range r[0]->r[1], r[2]->r[3]
		       )
{
  //UT_ASSERT2(r[2]>r[1]);
  return r[1] + (r[3]-r[1]) * MT_LINSTEPF( r[0], r[2], val );
}

/*  Say you don't want to change a value unless it's too small and screws 
 *  some of your computations up. Then, rather than doing a sharp conditional 
 *  branch, you can blend your value with your threshold, and do it smoothly 
 *  (say, with a cubic polynomial). Set m to be your threshold (anything 
 *  above m stays unchanged), and n the value things will take when your 
 *  value is zero.*/
inline float MtAlmostIdentity(float m, float n, float x )
{
   if( x>m ) return x;
    const float a = 2.0f*n - m;
    const float b = 2.0f*m - 3.0f*n;
    const float t = x*(1.0f/m);
    return (a*t + b)*t*t + n;
}

inline float MtSoftClampToZero( 
    		float m1, // x smaller than m1 => zero
    		float m2, // x larger than m2 => unchanged (identity) 
    		float x ) 
{
  /* maxima:
     f(x,m,n) := if (x > m) then x else \
	( ((2*n-m)*(x/m) + (2*m-3*n))*(x/m)^2 + n);\
     f2(x,m1,m2,n) := if (x < m1) then n else f(x-m1,m2-m1,n-m1)+m1;\
     plot2d( f2(x,0.1,0.3,0), [x,0,1]); 
   */

  return x < m1 ? 0 : (MtAlmostIdentity( m2-m1, 0-m1, x-m1 )+m1); 
}

#endif

#ifdef __cplusplus
extern "C" {
#endif

#define MT_LERP(a,b,t) ( (a) + (t) * ((b) - (a)) )

#define MT_SQR(x) ((x)*(x))
#define MT_POW3(x) ((x)*(x)*(x))
#define MT_POW4(x) ((x)*(x)*(x)*(x))
#define MT_POW6(x) ((x)*(x)*(x)*(x)*(x)*(x))

#define ABSDIFF( a, b ) (a) > (b) ? ( (a) - (b) ) : ( (b) - a )

//
// Matrix elements (row,column) and grid point (x,y) offsets
//

#define MT_MATEL(a,ncol,i,j) (*((a)+(LongInt)(ncol)*(i)+(j)))   
#define MT_MATEL1(a,ncol,i,j) (*((a)+(LongInt)(ncol)*((i)-1)+((j)-1)))   
#define MT_KRNEL(K,dmax,dx,dy) \
      MT_MATEL((K),(2*(dmax)+1),(dmax)+(dx),(dmax)+(dy))

#define MT_GXY( bmp_sx, x, y ) ((x)+((LongInt)(bmp_sx)*(y)))  
#define MT_GXYZ( grid_sx, grid_sxy, x, y, z ) ((x)+(y)*((LongInt)(grid_sx))+(z)*(LongInt)(grid_sxy))

#define MT_GCOORDS( ixy, sx, /*IN*/\
    		    x, y /*OUT*/)\
{ \
  y = (LongInt)(ixy) / (LongInt)(sx); \
  x = (LongInt)(ixy) % (LongInt)(sx); \
}

#define MT_GCOORDS_3D( sx, sxy, \
    		       idx, /*IN*/\
    		       x, y, z /*OUT*/)\
{ \
  z = (LongInt)(idx) / (LongInt)(sxy); \
  LongInt rem_ = (LongInt)(idx) % (LongInt)(sxy); \
  y = rem_ / (LongInt)(sx); \
  x = rem_ % (LongInt)(sx); \
}

/*
 *  Get grid neighborhood with Neumann (i.e. const. derivative) boundary
 *
 *				z7
 *
 *  			    z5	z0  z1
 *  				
 *  				z3
 *
 */

#define MT_GET_GRID_NEIGHBORS( Z,sx,sy, \
  			       x, y, \
    			       z0,z1,z5,z3,z7 /*OUT*/) \
{\
  z0 = Z[MT_GXY(sx,x,y)];\
  z1 = x<sx-1 ? Z[MT_GXY(sx,x+1,y)] : 2*z0 - Z[MT_GXY(sx,x-1,y)]; \
  z5 = x>0 ?    Z[MT_GXY(sx,x-1,y)] : 2*z0 - Z[MT_GXY(sx,x+1,y)]; \
  z3 = y<sy-1 ? Z[MT_GXY(sx,x,y+1)] : 2*z0 - Z[MT_GXY(sx,x,y-1)]; \
  z7 = y>0 ?    Z[MT_GXY(sx,x,y-1)] : 2*z0 - Z[MT_GXY(sx,x,y+1)]; \
}


//
// shift FFT frequency (0,0) to center by quadrant switching
// see matlab's "fftshift"
// x2,y2 = output coordinates
//
#define MT_FFTSHIFT( sx, sy, x, y, x2, y2 ) \
{ \
  if((y)<(sy)/2) \
  { \
    y2 = (y)+(sy)/2; \
    x2 = (x)<(sx)/2 ? (x)+(sx)/2 : (x)-(sx)/2; \
  } \
  else \
  { \
    y2=(y)-(sy)/2; \
    x2 = (x)<(sx)/2 ? (x)+(sx)/2 : (x)-(sx)/2; \
  } \
}

#define MT_FFTFRQ( sx, sy, x, y, fx, fy ) \
{ \
  fx = ((x)-(sx)/2)/(double)(sx); \
  fy = ((y)-(sy)/2)/(double)(sy); \
}

// 
// Compact storage scheme for symmetric matrices: 
// store only upper triangle 
//
#define MT_ISYMATEL_O(i,j)  (((i)+((((j)+1)*((j)))>>1)))
#define MT_SYMATEL_O(i,j) ((j)>=(i) ? MT_ISYMATEL_O(i,j) :\
				 MT_ISYMATEL_O(j,i))
#define MT_SYMATEL(a,i,j) (*((a)+MT_SYMATEL_O(i,j)))

#define MT_ROUND(x_) ((int)((x_)+0.5))

#define MT_VEC3_DOTPROD( u, v ) \
  ((u)[0]*(v)[0]+(u)[1]*(v)[1]+(u)[2]*(v)[2])

#define MT_VEC3_LEN(v) sqrt( ((double)v[0])*v[0]+\
			     ((double)v[1])*v[1]+\
			     ((double)v[2])*v[2] )

#define MT_VEC3_CPY( x_, y_ ) \
   { \
   (x_)[0]=(y_)[0]; (x_)[1]=(y_)[1]; (x_)[2]=(y_)[2]; \
   }
#define MT_VEC3_SET( x_, x1_, x2_, x3_ ) \
   { \
   (x_)[0]=x1_; (x_)[1]=x2_; (x_)[2]=x3_; \
   }
#define MT_VEC4_SET( x_, x1_, x2_, x3_, x4_ ) \
   { \
   (x_)[0]=x1_; (x_)[1]=x2_; (x_)[2]=x3_; (x_)[3]=x4_;\
   }
#define MT_VEC3_GET( x_, x1_, x2_, x3_ ) \
   { \
   x1_=(x_)[0]; x2_=(x_)[1]; x3_=(x_)[2]; \
   }

#define MT_VEC2_ROTATE( vy_, vx_, ss_, cs_ ) \
{ \
  (vy_)[0] = (vx_)[0]*cs_ + (vx_)[1]*ss_; \
  (vy_)[1] = -(vx_)[0]*ss_ + (vx_)[1]*cs_; \
}

#define MT_VEC2_DIST( x1,y1, x2,y2 ) \
  sqrt(((x2)-(x1))*((x2)-(x1))+((y2)-(y1))*((y2)-(y1)))


#define MT_VEC2_ROTATE2( vy0_, vy1_, vx0_, vx1_, ss_, cs_ ) \
{ \
  (vy0_) = (vx0_)*cs_ + (vx1_)*ss_; \
  (vy1_) = -(vx0_)*ss_ + (vx1_)*cs_; \
}

#define MT_VEC_DISTANCE2(x_,y_,n_,dist_) \
{ \
  int 	 i_; \
  MtReal tmp_; \
  for(i_=0,dist_=0;i_<n_;i_++) \
  { \
    tmp_ = x_[i_] - y_[i_]; \
    dist_ += tmp_ * tmp_; \
  } \
}

#define MtCreateMultStatFromTab( mt, opt )\
  MtCreateMultStat(NULL,(mt)->rows,(mt)->nCol,(mt)->nRow,NULL,opt,NULL)


/********************************************************************/
/*								    */
/********************************************************************/			 

typedef struct
{
  long  rndAcc;
  long  rndStab[ MT_RND_STAB_SZ ];
  long  rndIy; 
}
MtRndStream;

/********************************************************************/
/*								    */
/********************************************************************/			 

/* basic Vectors of real numbers */

typedef struct
{
  int 	 n;
  MtReal *elem;
} 
MtSeq;

/********************************************************************/
/*								    */
/********************************************************************/			 
typedef struct
{
  LongInt nRow,nRowMax,
	  nCol;
  float *data;
  float **rows;
}
MtTab;

/********************************************************************/
/*								    */
/********************************************************************/			 
typedef struct
{
  LongInt nRow,nRowMax,
	  nCol;
  double *data;
  double **rows;
}
MtMat;

/********************************************************************/
/*								    */
/********************************************************************/			 

/* basic statistics : mean, variance, min, max, sample size */ 

typedef struct      
{                   
  double n,   nConf,
         avg, avgConf,
         var, varConf,
	 span, min, max,
	 sumSqr ;
}
MtStat;

#define MT_STAT_INIT( s ) \
{ \
  (s)->min = MT_REAL_MAX; \
  (s)->max = MT_REAL_MIN; \
  (s)->avg = (s)->var = (s)->sumSqr = (s)->n = 0; \
  (s)->span=0; (s)->nConf=0; \
  (s)->avgConf = -1; \
  (s)->varConf = -1; \
}

#define MT_STAT_UPDATE( s, sample ) \
{ \
  (s)->n++; \
  if( sample > (s)->max ) (s)->max = (sample); \
  if( sample < (s)->min ) (s)->min = (sample);  \
  (s)->avg += (sample); \
  (s)->var += (sample)*((double)(sample)); \
}

#define MT_STAT_UPDATE_N( s, n_, min_, max_, savg, svar ) \
{ \
  (s)->n += n_; \
  UT_ASSERT2( UT_IMPLIES((n_)>0 && (savg)>0, !((max_)<(min_)) )); \
  if( max_ > (s)->max ) (s)->max = max_; \
  if( min_ < (s)->min ) (s)->min = min_;  \
  (s)->avg += savg; \
  (s)->var += svar; \
}

#define MT_STAT_UPDATE_STAT( s, s2 ) /* s = s+s2 */\
{ \
  MT_STAT_UPDATE_N( s, (s2)->n, (s2)->min, (s2)->max, \
      		       ((s2)->avg*(s2)->n), (s2)->sumSqr ) \
}

#define MT_STAT_RESULT( s ) \
{ \
  (s)->avg = (s)->n > 0  ? ((s)->avg / (s)->n) : 0; \
  (s)->sumSqr = (s)->var; \
  (s)->var =  (s)->n > 0 ? ( ((s)->var / (s)->n) - (s)->avg * (s)->avg ) : 0; \
  if((s)->var < 0) (s)->var = 0; \
  (s)->var = sqrt((s)->var); \
  (s)->span = (s)->max-(s)->min; \
} 

#define MT_STAT_SPRINT(str,st) \
    (sprintf(str,"n=%g, avg=%g, sig=%g, min=%g, max=%g",\
		  	(st)->n,(st)->avg,(st)->var,(st)->min,(st)->max))

#define MT_STAT_SPRINT2(str,st) \
    ((sprintf(str,"n=%g avg=%g, sig=%g, min=%g, max=%g",\
		  	(st)->n, (st)->avg,(st)->var,(st)->min,(st)->max)),str)

#define V2_SQLEN(x,y) ((x)*(x)+(y)*(y))
#define V2_LEN(x,y) sqrt(V2_SQLEN(x,y))

#define MT_SUM2MEANVAR( sum, sum2, n, mean, var ) \
{ \
  if(n>1)\
  {\
    var = ((sum2) - ((double)(sum)/(double)(n))*(sum))/((double)(n)-1);\
  }\
  else var=0;\
  if(var<0) var=0; \
  mean = (sum)/(double)(n); \
}

/********************************************************************/
/*								    */
/********************************************************************/			 

/* correlation coefficient */ 

typedef struct      
{                   
  double n, xVar, yVar, xAvg, yAvg,             
         xxAvg, yyAvg, xyAvg,
         corr, alpha1, alpha2;   /* output */      
}
MtStatCorr;


#define MT_STAT_CORR_INIT( s )	\
{ \
  (s)->xAvg = 0; \
  (s)->yAvg = 0; \
  (s)->xxAvg = 0; \
  (s)->yyAvg = 0; \
  (s)->xyAvg = 0; \
  (s)->n = 0; \
}

#define MT_STAT_CORR_UPDATE( s, x, y ) \
{ \
  (s)->n++; \
  (s)->xAvg += x; \
  (s)->xxAvg += x*x; \
  (s)->yAvg += y; \
  (s)->yyAvg += y*y; \
  (s)->xyAvg += x*y; \
}

#define MT_STAT_CORR_RESULT( s ) \
{ \
  (s)->xAvg /= (s)->n; \
  (s)->yAvg /= (s)->n; \
  (s)->xyAvg /= (s)->n; \
  (s)->xxAvg /= (s)->n; \
  (s)->yyAvg /= (s)->n; \
  (s)->xVar = (s)->xxAvg - (s)->xAvg * (s)->xAvg; \
  (s)->yVar = (s)->yyAvg - (s)->yAvg * (s)->yAvg;  \
  (s)->corr = ( (s)->xyAvg - (s)->xAvg*(s)->yAvg ) / \
    sqrt( (s)->xVar * (s)->yVar );  \
  (s)->alpha1 = ( (s)->xxAvg * (s)->yAvg - (s)->xAvg * (s)->xyAvg ) / \
    (s)->xVar; \
  (s)->alpha2 = ((s)->corr) * sqrt( (s)->yVar / (s)->xVar ); \
} 

/********************************************************************/
/*								    */
/********************************************************************/			 

/* queues of real numbers */

typedef struct
{
  int     len, 
  	  max, 
          in, out;    		/* fifo pointer */
  MtReal *buffer;   
}
MtRQueue;

#define MT_RQUEUE_PUSH( q, x ) \
{ \
  (q)->buffer[(q)->in++] = x; \
  if( (q)->in >= (q)->max ) (q)->in = 0; \
}

#define MT_RQUEUE_POP( q, x ) \
{ \
  x = (q)->buffer[ (q)->out++ ]; \
  if( (q)->out >= (q)->max ) (q)->out = 0; \
}

#define MT_RQUEUE_LEN( q ) ( (q)->len )

#define MT_NEAREST_POW2( n, n2, k, max ) \
{ \
  for(n2=1, k=0; n2<max; n2=n2<<1, k++) \
  { \
    if(n2>=n) break; \
  } \
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
#define MtVecCopy_f( y_, x_, n_ ) \
  memcpy( y_, x_, (n_) * sizeof( MtReal ) )

#define MtVecCopy_d( y_, x_, n_ ) \
  memcpy( y_, x_, (n_) * sizeof( double ) )

#define MtVecDbl2Flt( y_, x_, n_ ) \
{ \
  int k_; \
  for(k_=0;k_<n_;k_++) (y_)[k_] = (x_)[k_]; \
}
  
#define MtMatCopy_f( y_, x_, n_, m_ ) \
  memcpy( y_, x_, (n_)*(m_) * sizeof( MtReal ) )

#define MtMatCopy_fast( y_, x_, n_, m_ ) \
  memcpy( y_, x_, (n_)*(m_) * sizeof( *(x_) ) )

#define MtMatTransp( A, B, n, m )\
  MtMatCopy_d( A, B, -(n), m )

#define    MtVecClose( v ) \
{ \
  if( v ) MM_free(v); \
}  


/********************************************************************/
/*								    */
/********************************************************************/			 

/* probability distributions of sequences */ 

typedef struct
{
  MtStat   stat,
           statFreq;
  int      nFreq, nSmp, sorted;
  unsigned flags;
  MtReal   dx, xMin, xMax,
          *freq,
          *vDist;
  void	  *usrInfo;        
  double   domMin,domMax,domEps;
  float   *auxDstFuncX,
	  *auxDstFuncY;
  int	   auxDstFuncN,
    	   auxDstFuncValid;
}
MtDst;

typedef struct
{
  int nSteps;
  MtDst *dstSrc,
	*dstRef;
  double srcDomSpan;
  double *lookup;
}
MtDstTransTab;

/********************************************************************/
/*								    */
/********************************************************************/			 

/* moving average */

typedef struct
{
  int     max, 
          in, out;    		/* fifo pointer */
  double  val,
          sum,
          beta,
          buffer[ MT_MA_MAX ];   
}
MtMA;

/* moving average + volatility */

typedef struct
{
  int     max, 
          in, out;    		/* fifo pointer */
  double  val,
          sum,
          sumQ,
          beta,
          buffer[ MT_MA_MAX ];
}
MtMAV;

typedef struct MtNLRegAdm *MtNLRegHandle;

typedef MtReal (*MtVecParzenWindowFunc)( MtVector x, int d, MtReal h ); 

/********************************************************************/
/*								    */
/********************************************************************/			 
int MtGetStat2_f( LongInt n,
	       float *data,
	       float validMin, float validMax,
	       MtStat *stat
		 );

int MtGetStat_f( LongInt n,
	       float *data,
	       MtStat *stat
		 );
int MtGetStat( LongInt n,
	       double *data,
	       MtStat *stat
		 );

double MtDistPoint2LineSegSquared(
    		double px, double py,	// point  
    		double x1, double y1, double x2, double y2 // line
		);
double MtSmoothstep(double a, double b, double x);
double MtLinstep(double a, double b, double x);

void MtStereographProj( double *v_pq,	// IN/OUT: p,q 
    			  double *v_fg, //  OUT/OUT: f,g 
			  int dir // direction: 1 = p,q -> f,g 
			  	  //           -1 = inverse f,g -> p,q
			  );

int MtFFT(int dir,int log2n,double *x,double *y);

int MtFFT2D(double *CR, double *CI, int nx,int ny,int dir);
int MtFFT2D_f(float *CR, float *CI, int nx,int ny,int dir);

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
    	);

int    MtIsPow2 ( unsigned x );

double   MtRootN1( int n, double c );

int MtSolveCubic(double *a, double *x);

int MtSolveQuartic(double *a, double *x);

double   MtSolveDiff( double    (* f)( double x, double *parms ),  
                  double    (* df)( double x, double *parms ),
                  double   *parms,
                  double    x0,
                  double    epsilon,
                  int       maxIter );

int MtSolveBin( double    (* f)( double x, double *parms ),  
                double   *parms,
                double    x0,
                double    x1,
                double    epsilon,
                int       maxIter,
                double   *solution,
		double   *p_x0,
		double   *p_x1 );

int MtSolveBinFl( float    (* f)( float x, float *parms ),  
                float   *parms,
                float    x0,
                float    x1,
                float    epsilon,
                int       maxIter,
                float   *solution,
		float   *pEps );

double MtSignum( double x );
MtReal     MtSign( MtReal x );
double MtAbsPow( double x, double a );
float MtAbsPow_f( double x, double a );
double     MtRelDiff( double a, double b );
void 	   MtMapNormInterval( MtReal xMin, MtReal xRange, MtReal *t );
int        MtMAInit( MtMA *ma, int max, double y0 );
double     MtMAUpdate( MtMA *ma, double y );

int        MtMAVInit( MtMAV *mv, int max, double y0 );
void       MtMAVUpdate( MtMAV *ma, double y );
void       MtMAVGet( MtMAV *ma, double *mean, double *variance );

int        MtInitMAE( MtMA *ma, int max, double beta, double y0 );
double     MtUpdateMAE( MtMA *ma, double y );
int        MtStatGet_f( MtStat *stat, int nSmp, MtReal *samples );
int        MtStatGet( MtStat *stat, int nSmp, double *samples );
MtStat 	  *MtMultStatGet1( MtTab *X );
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
    		);
void       MtStatConfidenceGen( MtStat *stat, MtDst *dst, int nIter );
void 	   MtStatPrint( MtStat *stat );
MtDst 	  *MtDstOpen( MtReal *samples, int nSmp, int nFreq, unsigned opt );
MtDst *MtDstOpen2( MtReal *samples, int nSmp, int nFreq, 
    			  double domMin, double domMax,  // domain bounds 
			  double domEps,
    		          unsigned opt );
int        MtDstGet( MtDst *dst, MtReal *samples, int nSmp, MtReal *vFreq, 
                     int nFreq, MtReal *vDist, unsigned opt );
int 	   MtDstGetStat( MtDst *dst, MtReal *samples, int nSmp );
void MtDstGetStat2( MtDst *dst, double *median,
    		    double *skew, double *curt,
		    double *energy, double *entropy );
void MtDstGetStat2_old( MtDst *dst, 
    		    float *skew, float *curt,
		    float *energy, float *entropy );
void 	   MtDstClose( MtDst *dst );
MtReal 	   MtDstMedian( MtDst *dst );
double     MtVecMedian_d( double *x, int n );
MtReal    *MtDstDensity( MtDst *dst, MtReal *vDensY );                     
MtReal 	   MtDstEstimateDensity( MtDst *dst, MtReal x, float h );
MtReal     MtDstDistr( MtDst *dst, MtReal x, int use_old );
MtReal     MtDstInverse( MtDst *dst, MtReal u );
double     MtDstDistr2( MtDst *dst, double x, 
    		    int dir // direction 1=normal, -1=inverse transform
		    );
double     MtDstNormTrans( MtDst *dst, 
    		       double x, 
    		       int dir // direction 1=normal, -1=inverse transform
		      );
double     MtNormDst( double x );
double     MtNormDstInv( double y );
double     MtNormalDens( double x, double mu, double sig );
void       MtRndSeed( int seed );
double     MtRndGet( void );
int 	   MtRndIntGet( int max );
void 	   MtRndSSeed( MtRndStream *hdl, int seed );
double     MtRndSGet( MtRndStream *hdl );
int 	   MtRndSIntGet( MtRndStream *hdl, int max );
double 	   MtRndSNormGet( MtRndStream *hdl, double mean, double sigma );
MtReal     MtStudentT( MtReal alpha, int m );
double 	   MtGauss( double x );
double    *MtMakeGaussKernel( int r, double sig, unsigned opt );
MtReal 	   MtKernParzen( MtReal x, MtReal h );
float 	   MtClipInterval( float x, float a, float b, float c,
    			   unsigned opt );
double 	   MtFactorial( double n );
double 	   MtFactorialOdd( double n );
double     MtRndNormGet( double mean, double sigma );
MtReal     MtRndGenericGet( MtDst *dst );
int 	   MtCorrelation( double *xSeq, double *ySeq, int n, double *corr, 
                   double *alpha1, double *alpha2 );
int        MtCorrelation_f( MtReal *x, MtReal *y, int n, MtReal *corr, 
                          MtReal *alpha1, MtReal *alpha2 );                       
double	   MtCorrConfidence( double corr, int n, double alpha );
MtReal     MtCorrConfidence_f( MtReal corr, int n, MtReal alpha );                          
int        MtMultipleCorrelation( MtVector *pvx, MtVector y, int n, int m, 
	                          MtVector b, MtReal *rss, MtReal *r,
		                  MtVector t );
int	   MtMultipleCorrelationTest( int n, int p, MtReal noise,
                                      int rseed );		                  
MtVector   MtMatOpen( int n, int m );
void       MtMatClose( MtMatrix mat );
MtMatrix   MtMatMult( MtMatrix c, MtMatrix a, MtMatrix b, 
                    int n, int l, int m );
double    *MtMatMult_d( double *c, double *a, double *b, 
                    	int n, int l, int m );
double    *MtMatCopy_d( double *a, double *b, int n, int m );
int 	   MtMatLUdcmp_f(MtMatrix a,int n,int *indx, float *d, unsigned opt );
void 	   MtMatLubksb_f(MtMatrix a,int n,int *indx, float *b);
int 	   MtMatLUdcmp(double *A,int n,int *indx, double *d, unsigned opt );
void 	   MtMatLubksb(double *A,int n,int *indx, double *b);
MtMatrix   MtMatCholesky_f( MtMatrix a, int n, MtVector p );
double    *MtMatCholesky( double *a, int n, double *p );
MtVector   MtMatCholeskySolve( MtMatrix a, int n, MtVector p,
                               MtVector b, MtVector x );
MtVector   MtMatSolve( MtMatrix a, int n, MtVector b, MtVector x );
int MtMatGaussSolve( MtMatrix a, int n, MtVector b, int m );
double 	   MtMatDet(MtMatrix a,int n);
int MtMatSVD(
    double *A, int n, int m, // IN: matrix A (nxm)
    double *U, 	// OUT: matrix U (nxm)
    double *V,	// OUT: matrix V (mxm)
    double *S	// OUT: vector of m singular values
    );
int        MtMatSimpleSVD(double *U, double *S, double *V,
				int nRow, int nCol);
int MtMatPseudoInverseSVD( double *A, int n, 
    			   double *S, // optional: (nx1) eigenvalues 
			   double *V, // optional: (nxn) eigenvectors (=columns)
			   double eps,  // optional
    			   unsigned opt );
int        MtMatInverseLU( MtMatrix b, MtMatrix a, MtVector c, 
    	          	   int n, int *luIdx, float *d, unsigned opt );
int MtMatInverseLU_d( double *b, double *a, double *c, 
    	            int n, int *luIdx, double *d, unsigned opt );
int MtMatInverseLU_d_tst( double *b, double *a, double *c, 
    	            int n, int *luIdx, double *d, unsigned opt );
double 	   MtMatInverse(double *matrix, int num_dims);
void 	   MtMatPrint( MtMatrix a, int n, int m );
void       MtMatPrint_d( double *a, int n, int m, char *fformat );
MtVector   MtVecOpen( int n );
MtVector   MtVecAdd( MtVector y, MtVector x, int n );
MtVector   MtVecSub( MtVector y, MtVector x, int n );
double 	  *MtVecScale( double *y, int n, double alpha );
MtVector   MtVecScale_f( MtVector y, int n, MtReal alpha );
MtReal     MtVecDotprod( MtVector x, MtVector y, int n );
double     MtVecDotprod_d( double *x, double *y, int n );
double    *MtVec3CrossProd( double *a, double *b, double *c );
float 	  *MtVec3CrossProd_f( float *a, float *b, float *c );
MtReal 	   MtVecAbs( MtVector x, int n );
double     MtVecAbs_d( double *x, int n );
MtReal     MtVecAbs2( MtVector x, int n );
double     MtVecAbs2_d( double *x, int n );
MtReal     MtVecDistance( MtVector x, MtVector y, int n );
double     MtVecDistance_d( double *x, double *y, int n );
double     MtVecDistance2( double *x, double *y, int n );
MtReal 	   MtVecDistance2_f( MtVector x, MtVector y, int n );
MtVector   MtVecNorm( MtVector y, int n, MtReal len );
double 	  *MtVecNorm_d( double *y, int n, double len );
MtVector   MtVecConst( MtVector y, int n, MtReal c );
int 	   MtVecSample( MtVector y, MtVector x, int n, int d );
double *MtVecShuffle( double *y, int n, int rSeed );
MtVector   MtVecShuffle_f( MtVector y, int n, int rSeed );
MtVector   MtVecSort( MtVector y, int n );
MtReal     MtVecSelect(float arr[], unsigned long n, unsigned long k); 
void 	   MtVecStat_d( MtStat *s, double *x, LongInt n );
void MtVecStat( MtStat *s, float *x, LongInt n );
double    *MtVecFltToDbl( float *x, int n );
float     *MtVecDblToFlt( double *x, int n );
MtReal     MtVecMedian( MtVector x, int n, int is_sorted );
void MtConvToRank(double *w, // IN: sorted array
    		  int n );
void MtConvToRank_f(float *w, // IN: sorted array
    		  int n );
MtReal     MtVecSigmaLeft( MtVector x, int n );
MtReal     MtVecTrendStrength( MtVector x, int n );	
MtVector   MtVecMovTrendStrength( MtVector y, MtVector x, int n, int d );
MtVector MtVecMA_f( MtVector y, MtVector x, int n, int d );
double * MtVecMA( double * y, double * x, int n, int d );
double * MtVecMAE( double * y, double * x, int n, double beta );
MtVector   MtVecMovMedian( MtVector y, MtVector x, int n, int d );
double *   MtVecMovVar( double * y, double * x, int n, int d );
MtVector   MtVecDetrendVola( MtVector y, MtVector x, int n, int d );
MtVector   MtVecMomentumDiff( MtVector y, MtVector x, int n, int d );
double *MtVecMomentumLog( double * y, double * x, int n, int d );
double *MtVecMomentum( double *y, double *x, int n, int d );
MtVector   MtVecMovDiff( MtVector y, MtVector x, int n, int d, 
                         int signum );
MtVector   MtVecGenMRW( MtVector y, int n, MtReal mu, MtReal sigma, 
                        int rndSeed );
MtVector   MtVecGenChaos( MtVector y, int n, MtReal parm1,
		          int rndSeed );
MtReal 	   MtVecDensity( MtVector x, MtVector y, int d, int m,
    			         MtReal h, MtVecParzenWindowFunc winFunc );
MtReal 	   MtVecDensityNN( MtVector x, MtVector y, int d, int m,
    			       	   int k );

//
// Polynomials
//
double MtEvalPoly2D( double *P, // polynomial coefficients
    		     int no, 	// order
    		     double x, double y,
		     double *B  // OUT: (optional) basis functions at (x,y)
		     		//	(P not used)
		     );
int MtDerivPoly2D( 
    	double *P, // IN: polynomial coefficients
	int n, 	// polynomial order
	double x, double y,
	double *PX, double *PY // Derivatives (polynomials of order n-1)
	);

//
// Trigonometry
//
void MtNormVec2Grad( double *n, double *g );
void MtGrad2NormVec( double *g, double *n );

int MtAngular2NormVec( 
    		double tau, // in: tilt (azimuth, angle with x-axis) 
    		double sig, // in: slant (angle with z-axis)
    		double *vn  // out: 3D unit normal vector  
		);
int MtNormVec2Angular( double *vn, double *tau, double *sig );

#define MT_DEG2RAD( d ) ((((d)*MT_PI)/((double)180.0)))
#define MT_RAD2DEG( r ) ((r)*((double)180.0)/((double)MT_PI))

MtRQueue  *MtRQueueOpen( int max );
void 	   MtRQueueClose( MtRQueue *q );
void 	   MtRQueueInit( MtRQueue *q );
MtReal     MtOBSPrice( int type, MtReal s, MtReal t, MtReal p, 
                       MtReal v, MtReal r );
int        tprob(double t, int dof, double *probability);
void       MtIntShuffle( int *v, int n );
double MtCheckMonotony( double *f,  int	n,
    			       MtStat *ddstat );

int        MtTest( void );

int MtNLRegOpen( int nSmp, MtNLRegHandle *handle );
int MtNLRegFit( MtNLRegHandle handle, MtVector *pvx, MtVector y,
                int n, int m );
int MtNLRegForecast( MtNLRegHandle handle, MtVector x, MtReal *y );
int MtNLRegClose( MtNLRegHandle handle );
void MtDelTab( MtTab *mt );
MtTab *MtCreateTab( LongInt nRow, LongInt nCol );

void MtDelMat( MtMat *mt );
MtMat *MtCreateMat( LongInt nRow, LongInt nCol );
MtMat *MtCopyMat( MtMat *dst, MtMat *src );
MtMat *MtSetMat( MtMat *M, double val );

MtTab *MtGenRndMultNorm( int m,
    			 MtVector vMean, MtMatrix mSigma,
    		    	 int nSmp,
			 int seed );

int MtLMSRegressionTest( void );
int MtLMSRegression( 
    	float  **pSmp,  /* samples for model F(s[0]..s[m-1]) = s[m] */
	int      m,	/* model dimension */
	int      nSmp,	/* # samples */
	float	 alph,	/* step width scaling */
	float    sig,	/* sample variance sigma (-1: not supplied)*/
	int	 maxIter,
	float	 e0,
	float   *w	/* configuration (OUT) */ );

int MtCompareDistr( MtVector f1, int n1,   /* pre-sorted arrays */
    	            MtVector f2, int n2,
		    unsigned	opt,
		    double	*difCramer,
		    double	*difKomolSmir );

int MtCompareDistrBIN256( MtVector f1, int n1,   
    	            MtVector f2, int n2,
		    unsigned	opt,
		    double	*difCramer,
		    double	*difKomolSmir );

int MtCompareDistrBIN256_SumDiff( float *f1, int n1,   
    	            		  float *f2, int n2,
		    		  int		 doSum,
		    		  double	*difCramer,
		    		  double	*difKomolSmir );
#define MT_CD_CRAMER	1
#define MT_CD_KOMSMIR	2	

/*=====================================================================
 *
 *  Distribution transformation (histogram specification) lookup tables
 *
 * ===================================================================*/

MtDstTransTab *MtDstCreateTransTab( 
    			MtDst *dstSrc,  // source distribution
    		        MtDst *dstRef,  // target reference distribution
			int nSteps     // number of lookup table bins
    );

double MtDstLookupTransTab( MtDstTransTab *mtt, double x );

int MtDelTransTab( MtDstTransTab **p_mtt );

/*=====================================================================
 *
 *  MtHist: simple binned histograms
 *
 * ===================================================================*/

typedef struct
{
  int nb;	      // number of bins
  double xmin,xmax,   // support region
	 dx;	      // bin size 
  double *bin,*bin0;   // the bin array
  LongInt n_smp;
  double eps;
  double *sm_krn, sm_krnsum;		// smoothing kernel
  int	  sm_r;			// kernel radius (bins)
  double  sm_sig;		// smoothing kernel sigma	
}
MtHist;

void MtDelHist( MtHist *h );
MtHist *MtCreateHist( int nb,		// number of bins
    		      LongInt n_smp,	// expected number of samples
		      int use_zero_eps,	// <>0 : use eps for empty bin
    	      double xmin,
	      double xmax
	      );
double MtGetPercHist( MtHist *h, double perc);
int MtIFillHist( MtHist *h, double *data, float *data_f,
    		LongInt n, double a, int do_norm);
#define MtFillHist(h,data,n,do_norm) MtIFillHist(h,data,NULL,n,1,do_norm);
#define MtFillHist_f(h,data,n,do_norm) MtIFillHist(h,NULL,data,n,1,do_norm);
#define MtFillHist2_f(h,data,n,a,do_norm) MtIFillHist(h,NULL,data,n,a,do_norm);
double MtHistJDiv( MtHist *p, MtHist *q );
double MtHistCMDist( MtHist *p, MtHist *q );
double MtHistKSDist( MtHist *p, MtHist *q );
int MtFillHistFun( MtHist *h, double (*func)( double x ));
int MtSmoothHist( MtHist *h, 	// histogram 
    		  int 	  r,	// smoothing radius (max. neighb. 
		                // bin-offset considered for smoothing )
    		  double  sig   // kernel bandwidth 
		  );

/*=====================================================================
 *
 *  Hist3D : 3D binned histograms (e.g. bitmap RGB,LUV data)
 *
 * ===================================================================*/

typedef struct
{
  int n,	// number of data points
      nb,	// total number of bins
      nb_d,	// number of bins in each dimension
      nb_dd;	
  double *H;	// histogram data (nb x nb x nb)	
}
MtHist3D;

#define MT_HIST_3D(h,ir,ig,ib)\
  *(h->H + ((h)->nb_dd*((ir)))\
       + ((h)->nb_d*((ig)))\
       + ((ib)))

int MtDelHist3D( MtHist3D **p_h );
MtHist3D *MtCreateHist3D( MtTab *X,   // samples (3D, X->nCol=3)
    			  int nb_1d  // Number of bins in each dim
    			 );


/*=====================================================================
 *
 *  MtDistr: (Empirical) Probability Distribution
 *
 * ===================================================================*/
typedef struct
{
  MtStat   stat;	 	// basic statistics (mean, sigma)
  MtHist  *hist;  		// histogram
  double  *smp;			// data samples (sorted)
  int	   nSmp;		// number of samples
  int	   isSorted;		// samples already sorted ?
  double  *dens;
  int	   nDens;
  double   kbwOpt0;
}
MtDistr;

void MtDelDistr( MtDistr *ds );

MtDistr *MtCreateDistr( double *samples,	// samples
    			int 	n,		// number of samples		
			unsigned options
			);

double MtDistrInvCDF( MtDistr *ds,	// distribution
    		      double   u,	// prob. value: 0 < u < 1 
		      double   xMin,	// minimum possible sample value
		      double   xMax,	// maximum possible sample value
		      unsigned opt
		    	);

MtHist *MtGetDistrHist( MtDistr *ds,		// distribution 
    		    	int	 numBins 	// number of bins
		    	);
void MtSortDistr( MtDistr *ds );

int MtDistrCompareCVM( 
    	MtDistr *ds1, 		// IN: first distribution
        MtDistr *ds2, 		// IN: second distribution
        unsigned opt,		// IN: options
	double  *pDisCVM,	// OUT: Cramer v. Mises Distance	
	double  *pDisKS		// OUT: Komologrov Smirnov 
	);

int MtTstDistrInvCDF( double *vSmpRef, // reference distribution
    		      int     nSmpRef,
		      double  subSmpFrac,  // fractions sub sample
		      int     nSmpMC,
		      int     nIterMC,	// number of monte carlo rounds
		      int     rndSeed,
		      unsigned optCDF );
/*=====================================================================
 *
 *  MtMultStat: Multivariate Statistics
 *
 * ===================================================================*/

typedef struct
{
  /*** ATTENTION: when adding fields: update MtCopyMultStat ! */
  int 	dim,nSmp;
  double *mean,
	 *sigma,
	 *min,*max,
	 *cov,		// Covariance Matrix (dim x dim)
	 *corr;		// Correlation Matrix (dim x dim)
  //
  int 	  icovIsValid,
	  icorrIsValid;	 
  double *icov,		// Inverse Covariance (dim x dim)
	  covDet,	// Determinant of Covariance Matrix
	 *icorr,
	  corrDet;

  double *xtmp,
	 *ytmp;
  /*** ATTENTION: when adding fields: update MtCopyMultStat ! */
}
MtMultStat;

int MtDelMultStat( MtMultStat *mst );

MtMultStat *MtCreateMultStat( 
    		double **smp,  // IN: array of pointers to sample points
    		float **smp_f, // for float samples (deprecated) 
    		int dim,     // IN: sample point dimension, 
		int n,	     // IN: number of samples 
		int (*pgiCallback)( float fracCompleted,
		  		    void **pState ),
			// Optional progress indicator callback
		unsigned options,
		MtMultStat **p_mst
    		);

int MtCopyMultStat( MtMultStat *dst, MtMultStat *src );

int MtSetMultStat(MtMultStat *S,
    		   double *mean, 
		   double *sigma,
		   double *corr,
		   double *cov );

int MtPrintMultStat( MtMultStat *adm, char *fformat, int vl );
int MtMultStatPrint( int dim, int nSmp, double *vMean, 
    			double *vSigma, double *vCov, double *vCorr,
			char *fformat );
int MtMultStatPrintTwo( MtMultStat *S1, MtMultStat *S2 );

int MtCorrToCov(double *cov, int dim, double *corr, double *sigma);
int MtCovToCorr(double *corr, int dim, double *cov, double *sigma);

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
		);

double MtShrinkCovarianceToTarget( 
    		double *Z,  // Sample empirical covariance matrix
		double *T,  // Target covariance
		int dim,    // dimension	
		double nSamples,   // number of samples 
		double alpha  // shrinkage strength factor
			      // (mult. of optimal strength lambda )
					      	
    );

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
		);

/*=========================================================================
 *
 *   Complex numbers      
 *
 * =======================================================================*/

typedef struct cmplx {	/* complex number structure */
	double	re;		/* real part */
	double	im;		/* imaginary part */
} cmplx;

int cmplx_fixup(cmplx *ap);
cmplx cmplx_add(cmplx a, cmplx b);
cmplx cmplx_negate(cmplx a);
cmplx cmplx_mult(cmplx a, cmplx b);
cmplx cmplx_div(cmplx a, cmplx b);
cmplx cmplx_log(cmplx a);
cmplx cmplx_exp(cmplx a);
cmplx cmplx_pow(cmplx a, cmplx b);


#ifdef __cplusplus
}  // extern "C"
#endif

MSBG_NAMESPACE_END

#endif

