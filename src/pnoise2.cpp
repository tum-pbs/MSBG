/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#include <memory>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <assert.h>
#include <xmmintrin.h>
#include <vectorclass/vectorclass.h>
#include <vectorclass/special/vector3d.h>
#include "vectorclass_util.h"
#include "globdef.h"
#include "fastmath.h"
#include "mtool.h"
#include "rand.h"
#include "util.h"
#include "plot.h"
#include "pnoise.h"
#include "bitmap.h"
#include "panel.h"

using namespace MSBG_NAMESPACE;

namespace PNS
{

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
float genFractalScalar3D(
    	int sx, int sy, int sz,  // grid resolution
    	int x, int y, int z, // grid coords
	double seed,
    	double fscale,  // Feature scale in voxels
	double minscale,
	double octaves_leadin,
	double alpha,   // Fractal amplitude persistence
	int do_turb	// Apply fractal 'turbulence' (abs(noise))
	)
{  
  int dim = 3;
  double oct_offs = 0,
	 beta = 1.97;
#if 0
  double amplitude=1,
	 fscale=20,   // feature scale in grid cells
	 minscale=0.1,
	 oct_leadin=1,
	 oct_offs=1,
	 alpha = 1.5,
	 beta = 1.97,
	 do_turb = 1;
#endif

  double octaves = MAX( MAX( logf( fscale / minscale ) / MT_LOG2 , 1 ) + oct_offs , 1),
		 freq = 1./( fscale * powf( 2, octaves_leadin ) );

  //octaves = 2;

  double pos[3] = {(double)x*freq+seed*7.7,
		   (double)y*freq+seed*3.141,
		   (double)z*freq+seed*11.1};
  float f = PNS_Fractal( 
	      pos, dim, alpha, beta, 
	      octaves, octaves_leadin,
	      do_turb, 0.6, 1, 0,
	      0, 
	      PNS_NTYP_FAST, 0, 0);

  return f;
}

/*-------------------------------------------------------------------------------------*/
/* 										       */
/*-------------------------------------------------------------------------------------*/
void genFractalScalarField(  
    	int dim,	// dimension
    	int sx, int sy, int sz,  // grid resolution
	double seed,
    	double fscale,  // Feature scale in voxels
	double minscale,
	double octaves_leadin,
	double alpha,   // Fractal amplitude persistence
	int do_turb,	// Apply fractal 'turbulence' (abs(noise))
	float *field    // OUT: generated field
	)
{
  if(dim==2) 
  {
    sz=1;
  }
  else
  {
    UT_ASSERT0(dim==3);
  }

  #pragma omp parallel for schedule(static,50)
  for(LongInt idx=0;idx<(LongInt)sx*sy*sz; idx++)
  {
    int x,y,z;

    MT_GCOORDS_3D(sx,sx*(LongInt)sy, idx,
		  x,y,z);
  
    double oct_offs = 0,
	   beta = 1.97;
#if 0
    double amplitude=1,
	   fscale=20,   // feature scale in grid cells
	   minscale=0.1,
	   oct_leadin=1,
	   oct_offs=1,
	   alpha = 1.5,
	   beta = 1.97,
	   do_turb = 1;
#endif

    double octaves = MAX( MAX( logf( fscale / minscale ) / MT_LOG2 , 1 ) + oct_offs , 1),
		   freq = 1./( fscale * powf( 2, octaves_leadin ) );

    //octaves = 2;

    double pos[3] = {(double)x*freq+seed*7.7,
      		     (double)y*freq+seed*3.141,
		     (double)z*freq+seed*11.1};
    float f = PNS_Fractal( 
		pos, dim, alpha, beta, 
		octaves, octaves_leadin,
		do_turb, 0.6, 1, 0,
		0, 
		PNS_NTYP_FAST, 0, 0);

    //if(!inDomainRange(x,y,z)||IS_OBSTACLE(_cells[idx])) f =0;

    field[idx]=f;
  }
}  
  
/*=========================================================================
 *
 *
 *   
 * =======================================================================*/

#ifdef USE_RAN1
#if 0  // Choose between Mersenne Twister and SFMT generator
  #include "randomc/mersenne.cpp"
  #define  STOC_BASE CRandomMersenne     // define random number generator base class
#else
  #include "randomc/sfmt.cpp"
  #define  STOC_BASE CRandomSFMT
#endif
#endif

/*=========================================================================
 *
 *
 *   VVN : Vector Valued Noise
 *
 *   Fast SIMD version: fractal_vec3() is three times faster than scalar
 *   version.
 *
 *   During testing, very faint grid-like artefacts occured when summing 
 *   vector components that are not present with the non-SIMD version.
 *   The root cause of the slight artefacts is unknown. (Maybe related to
 *   the fact that in the non-SIMD version not only the integer 
 *   grid-point is shifted but also the fractional coordinate for the 
 *   output vector components)
 *   However: In practice, the output of the fast version is identical and 
 *   no artefacts could be observed in real rendering scenarios 
 *   
 * =======================================================================*/

namespace VVN
{

const int _B_= 0x1000,
	  BM = 0xFFF;

   //
   // Layout of gradient table
   //
   //   node 0: g0x,g1x,g2x,0, g0y,g1y,g2y,0, g0z,g1z,g2z,0, g0w,g1w,g2w,0
   //   node 1: ...
   //

const int GTAB3_ENTRY_LEN=12,
      	  GTAB4_ENTRY_LEN=16;

static float gtab3[(_B_+_B_+2)*GTAB3_ENTRY_LEN+1]
  			__attribute__((aligned(CPU_SIMD_WIDTH))); 

static float gtab4[(_B_+_B_+2)*GTAB4_ENTRY_LEN+1]
  			__attribute__((aligned(CPU_SIMD_WIDTH))); 

// Hash permutation table
static int ptab[_B_ + _B_ + 2];

void normalize3(float *v)
{
   double s = sqrt(v[0] * (double)v[0] + 
                   v[1] * (double)v[1] + 
	           v[2] * (double)v[2]);
   UT_ASSERT_FATAL(s>1e-5);
   v[0] = v[0] / s;
   v[1] = v[1] / s;
   v[2] = v[2] / s;
}

void normalize4(float *v)
{
   double s = sqrt(v[0] * (double)v[0] + 
                   v[1] * (double)v[1] + 
	           v[2] * (double)v[2] +
	           v[3] * (double)v[3]);
   UT_ASSERT_FATAL(s>1e-5);
   v[0] = v[0] / s;
   v[1] = v[1] / s;
   v[2] = v[2] / s;
   v[3] = v[3] / s;
}

// Initialize tables
int init(int rseed)
{
  TRC(("PNS: VVN::init()\n"));
  int rcThis=0,i,j,k;
  RND_Stream rnds;
  pns_rand_init(rseed,&rnds); 
  //if(rc) raiseRc(rc);

  pns_rand_init(rseed,&rnds);

#ifdef USE_RAN1
  STOC_BASE ran1(rseed);               // Make instance of random number generator
#endif

  // Setup gradient table

  // one random stream for vector component
  int rnds_[4]={ 782367, 4711119, 87267, 12345 };

  // Fill table
  for(i=0;i<_B_;i++)
  {
    // Compute 3 gradients per entry (one for each result vector component)
    for(k=0;k<3;k++)
    {
      float g[4],tmp[4];
      for(j=0;j<4;j++) tmp[j] = UtFastRand(&rnds_[k])-0.5;

      // 3d table
      for(j=0;j<4;j++) g[j] = tmp[j];
      normalize3(g);
      for(j=0;j<3;j++)
      {
        gtab3[12*i+4*j+k] = g[j];
      }	

      // 4d table
      for(j=0;j<4;j++) g[j] = tmp[j];
      normalize4(g);
      for(j=0;j<4;j++)
      {
        gtab4[16*i+4*j+k] = g[j];
      }	

    }
  }

   // Setup the permutation table
   for (i = 0 ; i < _B_ ; i++) 
   {
      ptab[i] = i;	
   }

   while (--i) 
   {
      k = ptab[i];
      {
	uint32_t h=i;
	UT_FAST_HASH_32(h);
        j = h % _B_;
      }
      UT_ASSERT_FATAL(j>=0&&j<=_B_-1);
      ptab[i] = ptab[j];
      ptab[j] = k;
   }

   // Fill the second half of the tables (wrap around)
   for (i = 0 ; i < _B_ + 2 ; i++) 
   {
     int j = _B_+i;
     ptab[j] = ptab[i];
     for(k=0;k<GTAB3_ENTRY_LEN;k++) 
       gtab3[j*GTAB3_ENTRY_LEN+k] = gtab3[i*GTAB3_ENTRY_LEN+k];
     for(k=0;k<GTAB4_ENTRY_LEN;k++) 
       gtab4[j*GTAB4_ENTRY_LEN+k] = gtab4[i*GTAB4_ENTRY_LEN+k];
   }

rcCatch:
  return rcThis;
}

// Setup
static inline void setup( double pos, // IN: input position coordinate
        int& b0, int& b1, // OUT: bracketing block indicies
        float& r0, float& r1, // OUT: residual coordinates
	float& s	// OUT: smooth interpolation curve
	)
{
   double p = pos+73856093.314;
   double p0 = floor(p);
   b0 = ((int)p0) & BM;
   b1 = (b0+1) & BM;
   float t = p - p0;     // flot precision for local coordinate
   r0 = t;
   r1 = t - 1.0f;
   s = ( t * t * (3.0f - 2.0f * t) );
}

// Linear interpolation
static Vec4f lerp(float t, const Vec4f& a, const Vec4f& b) 
{ 
  return a + t * (b - a); 
}

// Load random gradients
static inline void grad3( Vec4f &gx, Vec4f& gy, Vec4f& gz, 
    			  int index
			  )
{
  float *q = &gtab3[12*index];
  gx.load_a(q), gy.load_a(q + 4), gz.load_a(q + 8);
}

static inline void grad4( 
    	Vec4f &gx, Vec4f& gy, Vec4f& gz, Vec4f& gw, 
    	int index
	)
{
  float *q = &gtab4[16*index];
  gx.load_a(q), gy.load_a(q+4), gz.load_a(q+8), gw.load_a(q+12);
}


// Dot product
static inline Vec4f dot3( const Vec4f& gx, const Vec4f& gy, const Vec4f& gz,
    		          float& rx, float& ry, float& rz )
{
  return gx*rx + gy*ry + gz*rz;
}

static inline Vec4f dot4( 
    	const Vec4f& gx, const Vec4f& gy, const Vec4f& gz, const Vec4f& gw,
    	const Vec4f& rx, const Vec4f& ry, const Vec4f& rz, const Vec4f& rw )
{
  return gx*rx + gy*ry + gz*rz + gw*rw;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
// Vector valued noise 3d -> 3d
Vec4f noise3_vec3( double pos[3] )
{
  int bx0, bx1, by0, by1, bz0, bz1, b00, b10, b01, b11,
      i,j;
  float rx0,rx1, ry0,ry1, rz0,rz1,
	sx, sy, sz;
  Vec4f gx,gy,gz, u,v, a,b,c,d;

  //IACA_START

  // Setup: bracketing corner indices, 
  // residual coords and smooth curve 
  setup(pos[0], 
      	bx0,bx1, rx0,rx1, sx);

  setup(pos[1], 
        by0,by1, ry0,ry1, sy);

  setup(pos[2], 
        bz0,bz1, rz0,rz1, sz);

  // Corner indices
  i = ptab[ bx0 ];
  j = ptab[ bx1 ];

  b00 = ptab[ i + by0 ];
  b10 = ptab[ j + by0 ];
  b01 = ptab[ i + by1 ];
  b11 = ptab[ j + by1 ];

  // Load gradients and do dot products and perform lerps

  grad3( gx,gy,gz, b00 + bz0 ); u = dot3( gx,gy,gz, rx0, ry0, rz0);
  grad3( gx,gy,gz, b10 + bz0 ); v = dot3( gx,gy,gz, rx1, ry0, rz0);
  a = lerp(sx, u, v);

  grad3( gx,gy,gz,  b01 + bz0 ); u = dot3( gx,gy,gz, rx0,ry1,rz0);
  grad3( gx,gy,gz,  b11 + bz0 ); v = dot3( gx,gy,gz, rx1,ry1,rz0);
  b = lerp(sx, u, v);

  c = lerp(sy, a, b);

  grad3( gx,gy,gz,  b00 + bz1 ); u = dot3( gx,gy,gz, rx0,ry0,rz1);
  grad3( gx,gy,gz,  b10 + bz1 ); v = dot3( gx,gy,gz, rx1,ry0,rz1);
  a = lerp(sx, u, v);

  grad3( gx,gy,gz,  b01 + bz1 ); u = dot3( gx,gy,gz, rx0,ry1,rz1);
  grad3( gx,gy,gz,  b11 + bz1 ); v = dot3( gx,gy,gz, rx1,ry1,rz1);
  b = lerp(sx, u, v);

  d = lerp(sy, a, b);

  //IACA_END
   
  return lerp(sz, c, d);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
// Vector valued noise 4d -> 3d
Vec4f noise4_vec3( double pos[4] )
{
   int bx0, bx1, by0, by1, bz0, bz1, bw0, bw1;
   float rx0, rx1, ry0, ry1, 
	 rz0, rz1, rw0, rw1, 
	  sx, sy, sz, sw;

   Vec4f gx, gy, gz, gw, 
     	 a, b, c, d, e, t, u, v, f;

  // Setup: bracketing corner indices, 
  // residual coords and smooth curve 
  setup(pos[0], 
      	bx0,bx1, rx0,rx1, sx);

  setup(pos[1], 
        by0,by1, ry0,ry1, sy);

  setup(pos[2], 
        bz0,bz1, rz0,rz1, sz);

  setup(pos[3], 
        bw0,bw1, rw0,rw1, sw);

  // Corner indices
   int 
   i = ptab[ bx0 ],
   j = ptab[ bx1 ],
   b00 = ptab[ i + by0 ],
   b10 = ptab[ j + by0 ],
   b01 = ptab[ i + by1 ],
   b11 = ptab[ j + by1 ],
   b000 = ptab[ b00 + bz0 ],
   b100 = ptab[ b10 + bz0 ],
   b010 = ptab[ b01 + bz0 ],
   b110 = ptab[ b11 + bz0 ],
   b001 = ptab[ b00 + bz1 ],
   b101 = ptab[ b10 + bz1 ],
   b011 = ptab[ b01 + bz1 ],
   b111 = ptab[ b11 + bz1 ];

   grad4( gx,gy,gz,gw,  b000 + bw0 ) ; u = dot4( gx,gy,gz,gw,  rx0,ry0,rz0,rw0);   
   grad4( gx,gy,gz,gw,  b100 + bw0 ) ; v = dot4( gx,gy,gz,gw,  rx1,ry0,rz0,rw0);
   a = lerp(sx, u, v);
   grad4( gx,gy,gz,gw,  b010 + bw0 ) ; u = dot4( gx,gy,gz,gw,  rx0,ry1,rz0,rw0);
   grad4( gx,gy,gz,gw,  b110 + bw0 ) ; v = dot4( gx,gy,gz,gw,  rx1,ry1,rz0,rw0);
   b = lerp(sx, u, v);
   c = lerp(sy, a, b);
   grad4( gx,gy,gz,gw,  b001 + bw0 ) ; u = dot4( gx,gy,gz,gw,  rx0,ry0,rz1,rw0);
   grad4( gx,gy,gz,gw,  b101 + bw0 ) ; v = dot4( gx,gy,gz,gw,  rx1,ry0,rz1,rw0);
   a = lerp(sx, u, v);
   grad4( gx,gy,gz,gw,  b011 + bw0 ) ; u = dot4( gx,gy,gz,gw,  rx0,ry1,rz1,rw0);
   grad4( gx,gy,gz,gw,  b111 + bw0 ) ; v = dot4( gx,gy,gz,gw,  rx1,ry1,rz1,rw0);
   b = lerp(sx, u, v);
   d = lerp(sy, a, b);
   e = lerp(sz, c, d);

   grad4( gx,gy,gz,gw,  b000 + bw1 ) ; u = dot4( gx,gy,gz,gw,  rx0,ry0,rz0,rw1);   
   grad4( gx,gy,gz,gw,  b100 + bw1 ) ; v = dot4( gx,gy,gz,gw,  rx1,ry0,rz0,rw1);
   a = lerp(sx, u, v);
   grad4( gx,gy,gz,gw,  b010 + bw1 ) ; u = dot4( gx,gy,gz,gw,  rx0,ry1,rz0,rw1);
   grad4( gx,gy,gz,gw,  b110 + bw1 ) ; v = dot4( gx,gy,gz,gw,  rx1,ry1,rz0,rw1);
   b = lerp(sx, u, v);
   c = lerp(sy, a, b);
   grad4( gx,gy,gz,gw,  b001 + bw1 ) ; u = dot4( gx,gy,gz,gw,  rx0,ry0,rz1,rw1);
   grad4( gx,gy,gz,gw,  b101 + bw1 ) ; v = dot4( gx,gy,gz,gw,  rx1,ry0,rz1,rw1);
   a = lerp(sx, u, v);
   grad4( gx,gy,gz,gw,  b011 + bw1 ) ; u = dot4( gx,gy,gz,gw,  rx0,ry1,rz1,rw1);
   grad4( gx,gy,gz,gw,  b111 + bw1 ) ; v = dot4( gx,gy,gz,gw,  rx1,ry1,rz1,rw1);
   b = lerp(sx, u, v);
   d = lerp(sy, a, b);
   f = lerp(sz, c, d);
   
  return lerp(sw, e, f);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void fractal_vec3( 
    	      int dim,
	      double *pos, 
	      double octaves, 
	      int octaves_leadin,  // see http://www.planetside.co.uk/forums/index.php/topic,17566.msg170871.html#msg170871
              double lacunarity,
	      double alpha,  // roughness
	      int turb_typ,  // 0=fbm, 1=abs, 2=ridged
	      //float ridge_offs, float ridge_mult, float ridge_exp,

	      float *outvec	// 3d vector output
    	   )
{
  int i,k;
  float rem, alpha_inv=1./alpha,scale=1;
  double p[4];

  if(!(dim>=3&&dim<=4))
  {
   UT_ASSERT_NR(FALSE);
   return;
  }

  Vec4f val(0),sum(0);

/*  if(turb_typ==-2)
  {
   turb_typ=2;
   ridge_offs_exp = ridge_exp ? FMA_FastApproxPow(ridge_offs,ridge_exp) : 
				ridge_offs*ridge_offs;
  }*/

  for(k=0;k<dim;k++) p[k]=pos[k];

  octaves += octaves_leadin;

  int n=(int)octaves;

  if((rem = octaves - n) > 1.0E-12) n++;

  for (i=0;i<n;i++) 
  { 
    val = dim==3 ? VVN::noise3_vec3( p ) : 
          dim==4 ? VVN::noise4_vec3( p ) : 
	  0;
    for(k=0;k<dim;k++) p[k] *= lacunarity;

    if(turb_typ==1)
    {
      val = abs(val);
    }
    else if(turb_typ)
    {
      UT_ASSERT_NR(FALSE); // TODO
    }

    if(i>=(int)octaves) val*=rem;
    sum += val * scale;
    if(!(i<octaves_leadin)) scale *= alpha_inv;
  }

  float sum_[4] __attribute__((aligned(32)));
  sum.store_a(sum_);
  for(k=0;k<3;k++) outvec[k] = sum_[k];
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void noise3_vec3_( double pos[3],
    		   float *vec3_out )
{
  Vec4f res = noise3_vec3(pos);
  res.store_partial(3,vec3_out);
}

};  // namespace VVN


/*=========================================================================
 *
 *
 *   SIMD
 *
 *
 * =======================================================================*/

namespace simdtst
{
  typedef __m128i v128i;
  typedef __m128d v128d;
  typedef __m128 v128f;

  enum VectorSelect
  {
      Ax = 0, Ay = 1, Az = 2, Aw = 3,
      Bx = 8, By = 9, Bz = 10, Bw = 11,
  };

  template <VectorSelect S0, VectorSelect S1, VectorSelect S2, VectorSelect S3>
  inline v128f shuffle_ps(v128f x, v128f y)
  {
      return _mm_shuffle_ps(x, y, S0 + S1 * 4 + (S2 - Bx) * 16 + (S3 -Bx) * 64);
  }

  template<VectorSelect S0, VectorSelect S1, VectorSelect S2, VectorSelect S3>
  inline v128f blend_ps(v128f x, v128f y)
  {
      return _mm_blend_ps(x, y, (S0 / Bx) * 1 + (S1 / By) * 2 + (S2 / Bz) * 4 + (S3 / Bw) * 8);
  }

  template<VectorSelect S0, VectorSelect S1, VectorSelect S2>
  inline v128f blend_ps(v128f x, v128f y)
  {
      return _mm_blend_ps(x, y, (S0 / Bx) * 1 + (S1 / By) * 2 + (S2 / Bz) * 4);
  }

  template <VectorSelect S0, VectorSelect S1, VectorSelect S2, VectorSelect S3>
  inline v128i shuffle_epi32(v128i x)
  {
      return _mm_shuffle_epi32(x, S0 + S1 * 4 + S2 * 16 + S3 * 64);
  }


  inline v128f lerp(v128f t, v128f a, v128f b)
  {
      v128f r = _mm_sub_ps(b, a);
      r = _mm_mul_ps(r, t);
      r = _mm_add_ps(r, a);
      return r;
  }
  
  void trcv128f( __m128 f, const char *msg )
  {
    float tmp[4]; _mm_storeu_ps(tmp,f);
    TRCP(("%s: %g %g %g %g\n",msg,tmp[0],tmp[1],tmp[2],tmp[3]));
  }
  void trcv128d( __m128d f, const char *msg )
  {
    double tmp[4]; _mm_storeu_pd(tmp,f);
    TRCP(("%s: %g %g\n",msg,tmp[0],tmp[1]));
  }

  // A series of primes for the hash function
  const int NOISE_X = 1213;
  const int NOISE_Y = 6203;
  const int NOISE_Z = 5237;
  const int NOISE_SEED = 1039;
  const int NOISE_SHIFT = 13;

  const __m128i iONE = _mm_set1_epi32(1);
  const __m128 fONE = _mm_set1_ps(1.0f);
  const __m128 fTWO = _mm_set1_ps(2.0f);
  const __m128 fTHREE = _mm_set1_ps(3.0f);

  const int _B_= 0x1000,
	    BM = 0xFFF;

  float gtab3_v4[(_B_+_B_+2)*4] __attribute__((aligned(CPU_SIMD_WIDTH))); // AOS
  float gtab3_x[_B_+_B_+2] __attribute__((aligned(CPU_SIMD_WIDTH)));  // SOA
  float gtab3_y[_B_+_B_+2] __attribute__((aligned(CPU_SIMD_WIDTH)));  // SOA
  float gtab3_z[_B_+_B_+2] __attribute__((aligned(CPU_SIMD_WIDTH)));  // SOA

  uint32_t ptab[_B_+_B_+2] __attribute__((aligned(CPU_SIMD_WIDTH)));

  void normalize3(float *v)
  {
     double s;
     s = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
     v[0] = v[0] / s;
     v[1] = v[1] / s;
     v[2] = v[2] / s;
  }
   
  void noise_init(int rseed)
  {
    #define PNS_INIT_DO_VISU 0

    #if PNS_INIT_DO_VISU
    BmpBitmap *B=BmpNewBitmap(512,512,BMP_GREY|BMP_CLEAR);
    #endif
    
    RND_Stream rnds;
    pns_rand_init(rseed,&rnds);
    int i,k,j;
    for(i=0;i<_B_;i++)
    {
      volatile int dummy=0;
      for(j=0;j<3;j++) dummy += pns_rand(&rnds);  // force same rand stream as in old version
      float g[4]={0};
      for(int j=0;j<3;j++) 
      {
	g[j] = (double)((pns_rand(&rnds) % (_B_+_B_)) - _B_) / _B_;
	#if PNS_INIT_DO_VISU
	LongInt jj = (i%(B->sx*B->sy/4))*4+j;
	//TRCP(("%d\n",jj));
	B->dataGrey[jj] = g[j];
	#endif
      }

      // populate table for 3D gradients
      normalize3(g);
      for(int j=0;j<4;j++) gtab3_v4[i*4+j] = j<3 ? g[j] : 0;  // AOS
      gtab3_x[i] = g[0];  // SOA
      gtab3_y[i] = g[1];
      gtab3_z[i] = g[2];
      //TRCP(("g3_[%d] = %g %g %g\n",i,gtab3_x[i],gtab3_y[i],gtab3_z[i]));
    }
    

   // Setup the permutation table
   for (i = 0 ; i < _B_ ; i++) {
      ptab[i] = i;		// TODO*: use better PRNG than 'rand()'
   }

   while (--i) {
      k = ptab[i];
      ptab[i] = ptab[j = pns_rand(&rnds) % _B_];
      ptab[j] = k;
   }
#if 0
   for (i = 0 ; i < _B_ ; i++) 
   {
     TRCP(("p_[%d]=%d\n",i,ptab[i]));
   }
#endif

   // Fill the second half of the tables (wrap around)
   for (i = 0 ; i < _B_ + 2 ; i++) 
   {
     int j = _B_+i;
     ptab[j] = ptab[i];
     for(k=0;k<4;k++) gtab3_v4[j*4+k] = gtab3_v4[i*4+k];
     gtab3_x[j] = gtab3_x[i];  
     gtab3_y[j] = gtab3_x[i];  
     gtab3_z[j] = gtab3_x[i];  
   }

    #if PNS_INIT_DO_VISU
    PnlShowBitmap3(MiGlobalPnlTool(),B,1,BMP_GREY,0,NULL,0);
    BmpDeleteBitmap(&B);
    #endif
  }

  float grad_noise3d(simdtst::v128d pos_xy, simdtst::v128d pos_z)
  {
    // Round to the lowest integer boundary
    // DOUBLE PRECISION
    //IACA_START
    v128d pos_xy_0 = _mm_floor_pd(pos_xy);
    v128d pos_z_0 = _mm_floor_pd(pos_z);

    // Convert to integer to get the lower cube corner
    v128i xyi = _mm_cvtpd_epi32(pos_xy_0);
    v128i zi = _mm_cvtpd_epi32(pos_z_0);
    zi = shuffle_epi32<Az, Aw, Ax, Ay>(zi);
    v128i cube_pos0 = _mm_add_epi32(xyi, zi);

    // Add one to get the upper cube corner
    v128i cube_pos1 = _mm_add_epi32(cube_pos0, iONE);

    // Get fractional to lower cube corner
    // DOUBLE PRECISION
    v128d t_xy_d = _mm_sub_pd(pos_xy, pos_xy_0);
    v128d t_z_d = _mm_sub_pd(pos_z, pos_z_0);

    // Convert to higher-throughput float
    v128f t_xy_f = _mm_cvtpd_ps(t_xy_d);
    v128f t_zw_f = _mm_cvtpd_ps(t_z_d);
    v128f t0 = shuffle_ps<Ax, Ay, Bx, By>(t_xy_f, t_zw_f);

    // Get fractional to upper cube corner
    v128f t1 = _mm_sub_ps(t0, fONE);

    int x0 = vget_x(cube_pos0),
        y0 = vget_y(cube_pos0),
        z0 = vget_z(cube_pos0);

    int x1 = vget_x(cube_pos1),
        y1 = vget_y(cube_pos1),
        z1 = vget_z(cube_pos1);

    // Partial hash of extreme cube corner indices
    int ox0 = NOISE_X * x0 + NOISE_SEED;
    int oy0 = NOISE_Y * y0;
    int oz0 = NOISE_Z * z0;
    int ox1 = NOISE_X * x1 + NOISE_SEED;
    int oy1 = NOISE_Y * y1;
    int oz1 = NOISE_Z * z1;

    // Hash all cube corners
    int index0 = ox0 + oy0 + oz0;
    int index1 = ox1 + oy0 + oz0;
    int index2 = ox0 + oy1 + oz0;
    int index3 = ox1 + oy1 + oz0;
    int index4 = ox0 + oy0 + oz1;
    int index5 = ox1 + oy0 + oz1;
    int index6 = ox0 + oy1 + oz1;
    int index7 = ox1 + oy1 + oz1;
    index0 ^= (index0 >> NOISE_SHIFT);
    index1 ^= (index1 >> NOISE_SHIFT);
    index2 ^= (index2 >> NOISE_SHIFT);
    index3 ^= (index3 >> NOISE_SHIFT);
    index4 ^= (index4 >> NOISE_SHIFT);
    index5 ^= (index5 >> NOISE_SHIFT);
    index6 ^= (index6 >> NOISE_SHIFT);
    index7 ^= (index7 >> NOISE_SHIFT);
    index0 &= 0xFF;
    index1 &= 0xFF;
    index2 &= 0xFF;
    index3 &= 0xFF;
    index4 &= 0xFF;
    index5 &= 0xFF;
    index6 &= 0xFF;
    index7 &= 0xFF;

    // Lookup gradients
    v128f grad0 = _mm_load_ps(gtab3_v4+(index0<<2));
    v128f grad1 = _mm_load_ps(gtab3_v4+(index1<<2));
    v128f grad2 = _mm_load_ps(gtab3_v4+(index2<<2));
    v128f grad3 = _mm_load_ps(gtab3_v4+(index3<<2));
    v128f grad4 = _mm_load_ps(gtab3_v4+(index4<<2));
    v128f grad5 = _mm_load_ps(gtab3_v4+(index5<<2));
    v128f grad6 = _mm_load_ps(gtab3_v4+(index6<<2));
    v128f grad7 = _mm_load_ps(gtab3_v4+(index7<<2));

    // Project permuted offsets onto gradient vector
    v128f g0_ = _mm_dp_ps(grad0, (blend_ps<Ax, Ay, Az>(t0, t1)), 0x7F);
    v128f g1_ = _mm_dp_ps(grad1, (blend_ps<Bx, Ay, Az>(t0, t1)), 0x7F);
    v128f g2_ = _mm_dp_ps(grad2, (blend_ps<Ax, By, Az>(t0, t1)), 0x7F);
    v128f g3_ = _mm_dp_ps(grad3, (blend_ps<Bx, By, Az>(t0, t1)), 0x7F);
    v128f g4_ = _mm_dp_ps(grad4, (blend_ps<Ax, Ay, Bz>(t0, t1)), 0x7F);
    v128f g5_ = _mm_dp_ps(grad5, (blend_ps<Bx, Ay, Bz>(t0, t1)), 0x7F);
    v128f g6_ = _mm_dp_ps(grad6, (blend_ps<Ax, By, Bz>(t0, t1)), 0x7F);
    v128f g7_ = _mm_dp_ps(grad7, (blend_ps<Bx, By, Bz>(t0, t1)), 0x7F);

    // mix g0, g2, g4, g6 for lerp
    v128f g02__ = blend_ps<Ax, By, Az, Aw>(g0_, g2_);
    v128f g__46 = blend_ps<Ax, Ay, Az, Bw>(g4_, g6_);
    v128f g0246 = blend_ps<Ax, Ay, Bz, Bw>(g02__, g__46);

    // mix g1, g3, g5, g7 for lerp
    v128f g13__ = blend_ps<Ax, By, Az, Aw>(g1_, g3_);
    v128f g__57 = blend_ps<Ax, Ay, Az, Bw>(g5_, g7_);
    v128f g1357 = blend_ps<Ax, Ay, Bz, Bw>(g13__, g__57);

    // Apply a cubic fade to the near distance parameter for trilinear interpolation
    v128f r = _mm_mul_ps(t0, t0);
    v128f r0 = _mm_mul_ps(fTWO, t0);
    r0 = _mm_sub_ps(fTHREE, r0);
    r = _mm_mul_ps(r, r0);

    // Trilinear interpolation
    v128f rx = shuffle_ps<Ax, Ax, Bx, Bx>(r, r);
    v128f gx0123 = lerp(rx, g0246, g1357);

    v128f ry = shuffle_ps<Ay, Ay, By, By>(r, r);
    v128f gx1032 = shuffle_ps<Ay, Ax, Bw, Bz>(gx0123, gx0123);

    v128f gy0_1_ = lerp(ry, gx0123, gx1032);
    v128f rz = shuffle_ps<Az, Az, Bz, Bz>(r, r);
    v128f gy1_0_ = shuffle_ps<Az, Az, Bx, Bx>(gy0_1_, gy0_1_);
    v128f gz = lerp(rz, gy0_1_, gy1_0_);

    //IACA_END

    return vfget_x(gz);
  }

#if INSTRSET >= 7  // AVX
#ifdef NDEBUG
float noise3_v4( 
    		 //__m128 p
    		 const __m256d& p 
		 )
{
  //IACA_START

  __m256d p0 = _mm256_round_pd(p, 1);
  v128i ip0 = _mm_and_si128( _mm256_cvtpd_epi32(p0), 
      			     _mm_set1_epi32(BM) ),
	ip1 = _mm_and_si128( 
	        _mm_add_epi32( ip0, _mm_set1_epi32(1)),
	        _mm_set1_epi32(BM));
  
  v128f t0 = _mm256_cvtpd_ps( _mm256_sub_pd(p,p0) ),
	t1 = _mm_sub_ps(t0, _mm_set1_ps(1.0f));

  int bx0 = vget_x(ip0),
      by0 = vget_y(ip0),
      bz0 = vget_z(ip0),
      bx1 = vget_x(ip1),
      by1 = vget_y(ip1),
      bz1 = vget_z(ip1);

  int i = ptab[ bx0 ],
      j = ptab[ bx1 ],

      b00 = ptab[ i + by0 ],
      b10 = ptab[ j + by0 ],
      b01 = ptab[ i + by1 ],
      b11 = ptab[ j + by1 ];

  //TRCP(("%d %d %d %d\n",b00,b10,b01,b11));

  // Lookup gradients
  v128f grad0 = _mm_load_ps(gtab3_v4+ ((b00+bz0)<<2));
  v128f grad1 = _mm_load_ps(gtab3_v4+ ((b10+bz0)<<2));
  v128f grad2 = _mm_load_ps(gtab3_v4+ ((b01+bz0)<<2));
  v128f grad3 = _mm_load_ps(gtab3_v4+ ((b11+bz0)<<2));
  v128f grad4 = _mm_load_ps(gtab3_v4+ ((b00+bz1)<<2));
  v128f grad5 = _mm_load_ps(gtab3_v4+ ((b10+bz1)<<2));   
  v128f grad6 = _mm_load_ps(gtab3_v4+ ((b01+bz1)<<2));
  v128f grad7 = _mm_load_ps(gtab3_v4+ ((b11+bz1)<<2));

  // Project permuted offsets onto gradient vector
  v128f g0_ = _mm_dp_ps(grad0, (blend_ps<Ax, Ay, Az>(t0, t1)), 0x7F);
  v128f g1_ = _mm_dp_ps(grad1, (blend_ps<Bx, Ay, Az>(t0, t1)), 0x7F);
  v128f g2_ = _mm_dp_ps(grad2, (blend_ps<Ax, By, Az>(t0, t1)), 0x7F);
  v128f g3_ = _mm_dp_ps(grad3, (blend_ps<Bx, By, Az>(t0, t1)), 0x7F);
  v128f g4_ = _mm_dp_ps(grad4, (blend_ps<Ax, Ay, Bz>(t0, t1)), 0x7F);
  v128f g5_ = _mm_dp_ps(grad5, (blend_ps<Bx, Ay, Bz>(t0, t1)), 0x7F);
  v128f g6_ = _mm_dp_ps(grad6, (blend_ps<Ax, By, Bz>(t0, t1)), 0x7F);
  v128f g7_ = _mm_dp_ps(grad7, (blend_ps<Bx, By, Bz>(t0, t1)), 0x7F);

  // mix g0, g2, g4, g6 for lerp
  v128f g02__ = blend_ps<Ax, By, Az, Aw>(g0_, g2_);
  v128f g__46 = blend_ps<Ax, Ay, Az, Bw>(g4_, g6_);
  v128f g0246 = blend_ps<Ax, Ay, Bz, Bw>(g02__, g__46);
  // mix g1, g3, g5, g7 for lerp
  v128f g13__ = blend_ps<Ax, By, Az, Aw>(g1_, g3_);
  v128f g__57 = blend_ps<Ax, Ay, Az, Bw>(g5_, g7_);
  v128f g1357 = blend_ps<Ax, Ay, Bz, Bw>(g13__, g__57);

  // Apply a cubic fade to the near distance parameter for trilinear interpolation
  v128f r = _mm_mul_ps(t0, t0);
  v128f r0 = _mm_mul_ps(fTWO, t0);
  r0 = _mm_sub_ps(fTHREE, r0);
  r = _mm_mul_ps(r, r0);

  // Trilinear interpolation
  v128f rx = shuffle_ps<Ax, Ax, Bx, Bx>(r, r);
  v128f gx0123 = lerp(rx, g0246, g1357);

  v128f ry = shuffle_ps<Ay, Ay, By, By>(r, r);
  v128f gx1032 = shuffle_ps<Ay, Ax, Bw, Bz>(gx0123, gx0123);

  v128f gy0_1_ = lerp(ry, gx0123, gx1032);
  v128f rz = shuffle_ps<Az, Az, Bz, Bz>(r, r);
  v128f gy1_0_ = shuffle_ps<Az, Az, Bx, Bx>(gy0_1_, gy0_1_);
  v128f gz = lerp(rz, gy0_1_, gy1_0_);

  //TRCP(("%g %g %g -> %d %d %d %d ->%g\n",p[0],p[1],p[2],b00,b01,b10,b11,vfget_x(gz)));
  //IACA_END
  return vfget_x(gz);
}
#endif
#endif

}; // namespace simdtst

float PNS_simdtst_noise3_v4( double *pos )
{
  Vec4d p(pos[0],pos[1],pos[2],0);
//  Vec4f p(pos[0],pos[1],pos[2],0);
#if INSTRSET>=7  // AVX
#ifdef NDEBUG
  float f=simdtst::noise3_v4(p);
  return f;
#else
  UT_ASSERT_FATAL(FALSE);  // no AVX support in debug mode
  return 0;
#endif
#else
  return 0;
#endif
}

float PNS_donw_noise3d( double *p )
{
  __m128d p_xy=_mm_set_pd(p[0],p[1]),
	      p_z=_mm_set_pd(0,p[2]);
  float f=simdtst::grad_noise3d(p_xy,p_z);
  return f;
}

namespace simd
{

static inline Vec4f s_curve(const Vec4f& t) 
{
  return ( t * t * (3.0f - 2.0f * t) );
}

static inline Vec4f lerp(const Vec4f& t, const Vec4f& a, const Vec4f& b)
{
    return a + t * (b - a);
}

static inline Vec4f get_rand( Vec4ui& rnds ) 
{
  // return random gradient component in [-0.5,0.5]
  //  fast PRNG from http://www.iquilezles.org/www/articles/sfrand/sfrand.htm
  rnds *= 16807; 
#if 0
  Vec4f u = _mm_castsi128_ps((rnds)>>9|0x3f800000); 
  u -= 1.5f;
  return u; 
#else
  return (Vec4f)_mm_castsi128_ps((rnds)>>9|0x3f800000)-1.5f;
#endif
}

static inline
Vec4f grad_proj( const Vec4i& index, 
		 const Vec4f& dx, const Vec4f& dy, const Vec4f& dz )
{
#if 1
#if 0
  // seed local generator with index
  Vec4ui rnds = (Vec4ui)(index ^ (index>>13)) & 0xFF;  
#else
//  Vec4ui h_(1356267,34847,25324534,142423);
  Vec4ui h_ = (Vec4ui)index;
  h_ ^= h_ >> 16;	// 'murmur' hash 
  h_ *= 0x85ebca6b;
/*  h_ ^= h_ >> 13;
  h_ *= 0xc2b2ae35;
  h_ ^= h_ >> 16;*/
  Vec4ui rnds = h_;
#endif
  // dot product with random gradient components
  Vec4f 
  sum = ((Vec4f)_mm_castsi128_ps((rnds)>>9|0x3f800000)-1.5f) * dx;
  rnds *= 16807;
  sum += ((Vec4f)_mm_castsi128_ps((rnds)>>9|0x3f800000)-1.5f) * dy;
  rnds *= 16807;
  sum += ((Vec4f)_mm_castsi128_ps((rnds)>>9|0x3f800000)-1.5f) * dz;
  
  return sum;
#else


//  idx = idx ^ (idx >> simdtst::NOISE_HASH_SHIFT);
  Vec4ui idx = (Vec4ui)index;
#if 1
  idx ^= idx >> 16;	// 'murmur' hash 
  idx *= 0x85ebca6b;
  idx ^= idx >> 13;
  idx *= 0xc2b2ae35;
  idx ^= idx >> 16;
#endif
  idx = idx & simdtst::GTAB_MASK;

  __m128 *g = simdtst::gtab3;

  int igx = _mm_extract_epi32(idx,0),
      igy = _mm_extract_epi32(idx,1),
      igz = _mm_extract_epi32(idx,2);
      return g[igx]*dx + g[igy]*dy + g[igz]*dz;

#endif
}

static inline
Vec4f grad_proj2( Vec4i rnds, 
		 const Vec4f& dx, const Vec4f& dy, const Vec4f& dz )
{
#if 1
  rnds ^= (rnds >> 13);
  rnds = rnds & simdtst::BM;
#else
  rnds ^= rnds >> 16;	// 'murmur' hash 
  rnds *= 0x85ebca6b;
  rnds ^= rnds >> 13;
  rnds *= 0xc2b2ae35;
  rnds ^= rnds >> 16;
  rnds = rnds & simdtst::BM;
#endif
  Vec4f 
    gx = vgather4(simdtst::gtab3_x,rnds),
    gy = vgather4(simdtst::gtab3_y,rnds),
    gz = vgather4(simdtst::gtab3_z,rnds);
  return gx*dx+gy*dy+gz*dz;
}

Vec4f noise3d_simd4f( const Vec4f &px, const Vec4f &py, const Vec4f &pz  )
{
  //IACA_START
  Vec4f px0 = floor(px),
	py0 = floor(py),
	pz0 = floor(pz);

  // calculate hash index for lattice neighborhood points
  Vec4i ix0 =  truncate_to_int( px0 ),
        iy0 =  truncate_to_int( py0 ),
        iz0 =  truncate_to_int( pz0 );

  Vec4f dx0 = ( px - px0 ), dx1 = dx0-1.0f,
	dy0 = ( py - py0 ), dy1 = dy0-1.0f,
	dz0 = ( pz - pz0 ), dz1 = dz0-1.0f;

  Vec4f   tx = s_curve(dx0),
	  ty = s_curve(dy0),
  	  tz = s_curve(dz0);

    const int NOISE_X = 1213;
    const int NOISE_Y = 6203;
    const int NOISE_Z = 5237;
    const int NOISE_SEED = 1039;
    Vec4i ox0 = NOISE_X * ix0 + NOISE_SEED,
          oy0 = NOISE_Y * iy0,
          oz0 = NOISE_Z * iz0,
          ox1 = NOISE_X * (ix0+1) + NOISE_SEED,
          oy1 = NOISE_Y * (iy0+1),
          oz1 = NOISE_Z * (iz0+1);
    
  Vec4f g0 = grad_proj2( ox0 + oy0 + oz0, dx0, dy0, dz0  ),
  	g1 = grad_proj2( ox1 + oy0 + oz0, dx1, dy0, dz0  ),
	g2 = grad_proj2( ox0 + oy1 + oz0, dx0, dy1, dz0  ),
	g3 = grad_proj2( ox1 + oy1 + oz0, dx1, dy1, dz0  ),
        g4 = grad_proj2( ox0 + oy0 + oz1, dx0, dy0, dz1  ),
	g5 = grad_proj2( ox1 + oy0 + oz1, dx1, dy0, dz1  ),
	g6 = grad_proj2( ox0 + oy1 + oz1, dx0, dy1, dz1  ),
	g7 = grad_proj2( ox1 + oy1 + oz1, dx1, dy1, dz1  );

  Vec4f t; 

  Vec4f gx0 = lerp(tx, g0, g1),
        gx1 = lerp(tx, g2, g3),
        gx2 = lerp(tx, g4, g5),
        gx3 = lerp(tx, g6, g7);

  Vec4f   gy0 = lerp(ty, gx0, gx1),
          gy1 = lerp(ty, gx2, gx3);

  //IACA_END

  return lerp(tz, gy0, gy1);
}

Vec4f noise3d_simd4( const Vec4d &px, const Vec4d &py, const Vec4d &pz  )
{
  //IACA_START
  Vec4d px0 = floor(px),
	py0 = floor(py),
	pz0 = floor(pz);

  // calculate hash index for lattice neighborhood points
  Vec4i ix0 =  _mm256_cvtpd_epi32( px0 ),
        iy0 =  _mm256_cvtpd_epi32( py0 ),
        iz0 =  _mm256_cvtpd_epi32( pz0 );

  Vec4f dx0 = _mm256_cvtpd_ps( px - px0 ), dx1 = dx0-1.0f,
	dy0 = _mm256_cvtpd_ps( py - py0 ), dy1 = dy0-1.0f,
	dz0 = _mm256_cvtpd_ps( pz - pz0 ), dz1 = dz0-1.0f;

  Vec4f   tx = s_curve(dx0),
	  ty = s_curve(dy0),
  	  tz = s_curve(dz0);

    const int NOISE_X = 1213;
    const int NOISE_Y = 6203;
    const int NOISE_Z = 5237;
    const int NOISE_SEED = 1039;
    Vec4i ox0 = NOISE_X * ix0 + NOISE_SEED,
          oy0 = NOISE_Y * iy0,
          oz0 = NOISE_Z * iz0,
          ox1 = NOISE_X * (ix0+1) + NOISE_SEED,
          oy1 = NOISE_Y * (iy0+1),
          oz1 = NOISE_Z * (iz0+1);
    
  Vec4f g0 = grad_proj2( ox0 + oy0 + oz0, dx0, dy0, dz0  ),
  	g1 = grad_proj2( ox1 + oy0 + oz0, dx1, dy0, dz0  ),
	g2 = grad_proj2( ox0 + oy1 + oz0, dx0, dy1, dz0  ),
	g3 = grad_proj2( ox1 + oy1 + oz0, dx1, dy1, dz0  ),
        g4 = grad_proj2( ox0 + oy0 + oz1, dx0, dy0, dz1  ),
	g5 = grad_proj2( ox1 + oy0 + oz1, dx1, dy0, dz1  ),
	g6 = grad_proj2( ox0 + oy1 + oz1, dx0, dy1, dz1  ),
	g7 = grad_proj2( ox1 + oy1 + oz1, dx1, dy1, dz1  );

  Vec4f t; 

  Vec4f gx0 = lerp(tx, g0, g1),
        gx1 = lerp(tx, g2, g3),
        gx2 = lerp(tx, g4, g5),
        gx3 = lerp(tx, g6, g7);

  Vec4f   gy0 = lerp(ty, gx0, gx1),
          gy1 = lerp(ty, gx2, gx3);

  //IACA_END

  return lerp(tz, gy0, gy1);
}


static inline
Vec4f grad_proj_i( Vec4i idx, 
		 const Vec4f& x, const Vec4f& y, const Vec4f& z )
{
#if 1
  idx ^= idx >> 16;	// 'murmur' hash 
  idx *= 0x85ebca6b;
  idx ^= idx >> 13;
  idx *= 0xc2b2ae35;
  idx ^= idx >> 16;
#endif
  Vec4i h = idx & 15;

  Vec4f u = select( h<8, x, y );
  Vec4f      v = select( h<4, y, select(h==12||h==14,x,z)),
      a = select( (h&1)==0, u, -u ),
      b = select( (h&2)==0, v, -v );

  return a+b;
}

Vec4f noise3d_simd4_i( const Vec4d &px, const Vec4d &py, const Vec4d &pz  )
{
//  IACA_START
  Vec4d px0 = floor(px),
	py0 = floor(py),
	pz0 = floor(pz);

  // calculate hash index for lattice neighborhood points
  Vec4i ix0 =  _mm256_cvtpd_epi32( px0 ),
        iy0 =  _mm256_cvtpd_epi32( py0 ),
        iz0 =  _mm256_cvtpd_epi32( pz0 );

  Vec4f dx0 = _mm256_cvtpd_ps( px - px0 ), dx1 = dx0-1.0f,
	dy0 = _mm256_cvtpd_ps( py - py0 ), dy1 = dy0-1.0f,
	dz0 = _mm256_cvtpd_ps( pz - pz0 ), dz1 = dz0-1.0f;

  Vec4f   tx = s_curve(dx0),
	  ty = s_curve(dy0),
  	  tz = s_curve(dz0);

    const int NOISE_X = 73856093;
    const int NOISE_Y = 19349663;
    const int NOISE_Z = 83492791;
    Vec4i ox0 = NOISE_X * ix0,
          oy0 = NOISE_Y * iy0,
          oz0 = NOISE_Z * iz0,
          ox1 = NOISE_X * (ix0+1),
          oy1 = NOISE_Y * (iy0+1),
          oz1 = NOISE_Z * (iz0+1);
    
  Vec4f g0 = grad_proj_i( ox0 + oy0 + oz0, dx0, dy0, dz0  ),
  	g1 = grad_proj_i( ox1 + oy0 + oz0, dx1, dy0, dz0  ),
	g2 = grad_proj_i( ox0 + oy1 + oz0, dx0, dy1, dz0  ),
	g3 = grad_proj_i( ox1 + oy1 + oz0, dx1, dy1, dz0  ),
        g4 = grad_proj_i( ox0 + oy0 + oz1, dx0, dy0, dz1  ),
	g5 = grad_proj_i( ox1 + oy0 + oz1, dx1, dy0, dz1  ),
	g6 = grad_proj_i( ox0 + oy1 + oz1, dx0, dy1, dz1  ),
	g7 = grad_proj_i( ox1 + oy1 + oz1, dx1, dy1, dz1  );

  Vec4f t; 

  Vec4f gx0 = lerp(tx, g0, g1),
        gx1 = lerp(tx, g2, g3),
        gx2 = lerp(tx, g4, g5),
        gx3 = lerp(tx, g6, g7);

  Vec4f   gy0 = lerp(ty, gx0, gx1),
          gy1 = lerp(ty, gx2, gx3);

//  IACA_END

  return lerp(tz, gy0, gy1);
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
Vec4f fractal_simd4_i( 
    	      int dim, 
	      const Vec4d& px_, const Vec4d& py_, const Vec4d& pz_, 
	      double octaves, 
	      int octaves_leadin,  // see http://www.planetside.co.uk/forums/index.php/topic,17566.msg170871.html#msg170871
              double lacunarity,
	      double alpha,  // roughness
	      int turb_typ  // 0=fbm, 1=abs, 2=ridged
	      //float ridge_offs, float ridge_mult, float ridge_exp,
    	   )
{
  int i;
  float rem, alpha_inv=1./alpha,scale=1;

  Vec4d px=px_,py=py_,pz=pz_;

  Vec4f val(0),sum(0);

/*  if(turb_typ==-2)
  {
   turb_typ=2;
   ridge_offs_exp = ridge_exp ? FMA_FastApproxPow(ridge_offs,ridge_exp) : 
				ridge_offs*ridge_offs;
  }*/

  octaves += octaves_leadin;

  int n=(int)octaves;

  if((rem = octaves - n) > 1.0E-12) n++;

  for (i=0;i<n;i++) 
  { 
    val = dim==3 ? simd::noise3d_simd4_i( px,py,pz ) : 
	  0;

    px *= lacunarity;
    py *= lacunarity;
    pz *= lacunarity;

    if(turb_typ==1)
    {
      val = abs(val);
    }
    else if(turb_typ)
    {
      UT_ASSERT_NR(FALSE); // TODO
    }

    if(i>=(int)octaves) val*=rem;
    sum += val * scale;
    if(!(i<octaves_leadin)) scale *= alpha_inv;
  }

  return sum;
}


static inline void trc_v8ui(const Vec8ui& g, const char *msg)
{
  TRCP(("%s: %d %d %d %d %d %d %d %d\n",msg,
	g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7]));
}
static inline void trc_v8i(const Vec8i &g, const char *msg)
{
  TRCP(("%s: %d %d %d %d %d %d %d %d\n",msg,
	g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7]));
}
static inline void trc_v8f(const Vec8f& g, const char *msg)
{
  TRCP(("%s: %g %g %g %g %g %g %g %g\n",msg,
	g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7]));
}


static 
inline Vec8f grad_proj(const Vec8f* g, const Vec8f& dx, const Vec8f& dy, const Vec8f& dz)
{
  return g[0]*dx+g[1]*dy+g[2]*dz;
}

static inline Vec8f s_curve(const Vec8f& t) 
{
  return ( t * t * (3.0f - 2.0f * t) );
}

static inline Vec8f lerp(const Vec8f& t, const Vec8f& a, const Vec8f& b)
{
    return a + t * (b - a);
}

Vec8f noise3d_simd8( 
    		const Vec4d &px_l /*0123*/, const Vec4d &px_h /*4567*/,  // position x coords
    		const Vec4d &py_l /*0123*/, const Vec4d &py_h /*4567*/,  // position y coords
    		const Vec4d &pz_l /*0123*/, const Vec4d &pz_h /*4567*/  // position z coords
    )
{
  // Get fractional coordinates in double precision
  // then convert deltas to float precision 
  Vec4d dl_d, dh_d;
  
  // x coords
  Vec4d px_l_0 = floor(px_l), 	px_h_0 = floor(px_h);
  dl_d = px_l - px_l_0; 	dh_d = px_h - px_h_0;
  Vec8f dx0( _mm256_cvtpd_ps( dl_d ), _mm256_cvtpd_ps(dh_d) ),
	dx1 = dx0 - 1.0f;

  // y coords
  Vec4d py_l_0 = floor(py_l), 	py_h_0 = floor(py_h); 
  dl_d = py_l - py_l_0; 	dh_d = py_h - py_h_0;
  Vec8f dy0( _mm256_cvtpd_ps( dl_d ), _mm256_cvtpd_ps(dh_d) ),
	dy1 = dy0 - 1.0f;

  // z coords
  Vec4d pz_l_0 = floor(pz_l), 	pz_h_0 = floor(pz_h);
  dl_d = pz_l - pz_l_0; 	dh_d = pz_h - pz_h_0;
  Vec8f dz0( _mm256_cvtpd_ps( dl_d ), _mm256_cvtpd_ps(dh_d) ),
	dz1 = dz0 - 1.0f;

  // calculate hash index for lattice neighborhood points
  Vec8i ix0 = Vec8i( _mm256_cvtpd_epi32( px_l_0 ), _mm256_cvtpd_epi32( px_h_0) ),
        iy0 = Vec8i( _mm256_cvtpd_epi32( py_l_0 ), _mm256_cvtpd_epi32( py_h_0) ),
        iz0 = Vec8i( _mm256_cvtpd_epi32( pz_l_0 ), _mm256_cvtpd_epi32( pz_h_0) );

    const int NOISE_X = 1213;
    const int NOISE_Y = 6203;
    const int NOISE_Z = 5237;
    const int NOISE_SEED = 1039;
    Vec8i ox0 = NOISE_X * ix0 + NOISE_SEED,
          oy0 = NOISE_Y * iy0,
          oz0 = NOISE_Z * iz0,
          ox1 = NOISE_X * (ix0+1) + NOISE_SEED,
          oy1 = NOISE_Y * (iy0+1),
          oz1 = NOISE_Z * (iz0+1);

    return 0;
}

} // namespace simd

void PNS_noise3d_simd4_f( const float *px, const float *py, const float *pz, int n,
    			float *result )
{
#if 0
  UT_ASSERT_FATAL(n<=4);
  Vec4d x,y,z;
  x.load_partial(n,px);
  y.load_partial(n,py);
  z.load_partial(n,pz);
#else
  Vec4f x(px[0]),y(py[0]),z(pz[0]);
#endif
  Vec4f f = simd::noise3d_simd4f( x,y,z );
  f.store(result);
}

void PNS_noise3d_simd8( const double *px, const double *py, const double *pz, int n,
    			float *result )
{
#if 0
  Vec4d xl,xh;
  xl.load_partial(MIN(n,4),px); xh.load_partial(MAX(n-4,0),px+4);
  Vec4d yl,yh;
  yl.load_partial(MIN(n,4),py); yh.load_partial(MAX(n-4,0),py+4);
  Vec4d zl,zh;
  zl.load_partial(MIN(n,4),pz); zh.load_partial(MAX(n-4,0),pz+4);
#else
  Vec4d xl(px[0]),xh(px[0]);
  Vec4d yl(py[0]),yh(py[0]);
  Vec4d zl(pz[0]),zh(pz[0]);

#endif
  Vec8f f = simd::noise3d_simd8( xl, xh,
      		 	   	 yl, yh,
				 zl, zh );
  f.store(result);
}

}; // namespace PNS

int PNS_Init2(int rseed)
{
   TRC(("PNS_Init2(%d)\n",rseed));
   PNS::VVN::init(rseed);
#if 1
   PNS::simdtst::noise_init(rseed);
#endif
   return 0;
}


