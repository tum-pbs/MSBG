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
#include <iacaMarks.h>
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
#include "noise.h"
#include "pnoise.h"
#include "bitmap.h"
#include "panel.h"

namespace CurlNoise
{
#define DO_NORMALIZE_GRAD
  
#define SIMDW 8
#define s_curve2(t) (t * t * t * (t * (t * 6.f - 15.f) + 10.f))

#if 1
#define HASH4D(x,y,z,w) \
    ( ((x)*73856093) ^ \
      ((y)*19349663) ^ \
      ((z)*83492791) ^ \
      ((w)*56598313) ) 
      //h1000 = HASH4D( ix0+1, iy0, iz0, iw0 );
#else
#define HASH4D(x,y,z,w) \
    ( ((x)*93) ^ \
      ((y)*63) ^ \
      ((z)*91) ^ \
      ((w)*13) ) 
#endif

#define RNDF(rnds)\
  (rnds *= 16807,\
   2.0f*((Vec8f)_mm256_castsi256_ps((rnds)>>9|0x3f800000)-1.5f))

#define NORMALIZE(x,y,z,w) \
{ \
  Vec8f len = x*x + y*y + z*z + w*w,\
	s = approx_rsqrt(len+1e-6f);\
  x *= s;\
  y *= s;\
  z *= s;\
  w *= s;\
}
  
#define LERP(t, a, b) mul_add( b-a, t, a )

static inline 
Vec8f dot_rnd_grad( Vec8ui& hash,
    		    const Vec8f& rx, const Vec8f& ry, const Vec8f& rz, const Vec8f& rw )
{
  Vec8f gx=RNDF(hash),
	gy=RNDF(hash), 
	gz=RNDF(hash),
	gw=RNDF(hash);
  #ifdef DO_NORMALIZE_GRAD
  NORMALIZE(gx,gy,gz,gw)
  #else
  UT_ASSERT0(FALSE); // TODO: ensure correct output range 
  #endif
  return rx*gx + ry*gy + rz*gz + rw*gw;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static inline void curlNoisePotential4D( 

   int nPoints,  
   
   const float *__restrict posX,  // 4D input position (assumed zero padded to SIMDW)
   const float *__restrict posY,
   const float *__restrict posZ,
   float posW,  // fixed position in time

   const float *__restrict compOut,  // Component of 3D output vector 0..2 (assumed zero padded)

   float freqX,	// Frequencies per input coordinate
   float freqY,
   float freqZ,
   float freqW,

   float seed,	// Random seed

   float *out   // Output values

   )

{
  for(int iPoint=0; iPoint < nPoints; iPoint += SIMDW)
  {
    Vec8f poff = seed + 47.11*Vec8f().load(compOut+iPoint),
      	  px = poff + freqX * Vec8f().load(posX+iPoint),
	  py = poff + freqY * Vec8f().load(posY+iPoint),
	  pz = poff + freqZ * Vec8f().load(posZ+iPoint),
	  pw = poff + freqW * posW;

    Vec8f px0 = floor(px),
	  py0 = floor(py),
	  pz0 = floor(pz),
	  pw0 = floor(pw);

    Vec8f rx0 = ( px - px0 ), rx1 = rx0-1.0f,
	  ry0 = ( py - py0 ), ry1 = ry0-1.0f,
	  rz0 = ( pz - pz0 ), rz1 = rz0-1.0f,
	  rw0 = ( pw - pw0 ), rw1 = rw0-1.0f;

    Vec8f   sx = s_curve2(rx0),
	    sy = s_curve2(ry0),
	    sz = s_curve2(rz0),
	    sw = s_curve2(rw0);

    Vec8i  ix =  truncate_to_int( px0 ),
	   iy =  truncate_to_int( py0 ),
	   iz =  truncate_to_int( pz0 ),
	   iw =  truncate_to_int( pw0 );

    Vec8ui h;
    Vec8f u,v,
	  a,b,c,d,e,f, res;

    u = dot_rnd_grad( h = HASH4D( ix, 	iy, iz, iw ), 		rx0, ry0, rz0, rw0 );
    v = dot_rnd_grad( h = HASH4D( ix+1, iy, iz, iw ), 		rx1, ry0, rz0, rw0 );
    a = LERP(sx, u, v);

    u = dot_rnd_grad( h = HASH4D( ix, 	iy+1, iz, iw ), 	rx0, ry1, rz0, rw0 );
    v = dot_rnd_grad( h = HASH4D( ix+1, iy+1, iz, iw ), 	rx1, ry1, rz0, rw0 );
    b = LERP(sx, u, v);
    c = LERP(sy, a, b);

    u = dot_rnd_grad( h = HASH4D( ix, 	iy, iz+1, iw ), 	rx0, ry0, rz1, rw0 );
    v = dot_rnd_grad( h = HASH4D( ix+1, iy, iz+1, iw ), 	rx1, ry0, rz1, rw0 );
    a = LERP(sx, u, v);

    u = dot_rnd_grad( h = HASH4D( ix, 	iy+1, iz+1, iw ), 	rx0, ry1, rz1, rw0 );
    v = dot_rnd_grad( h = HASH4D( ix+1, iy+1, iz+1, iw ), 	rx1, ry1, rz1, rw0 );
    b = LERP(sx, u, v);
    d = LERP(sy, a, b);
    e = LERP(sz, c, d);


    u = dot_rnd_grad( h = HASH4D( ix, 	iy, iz, iw+1 ), 	rx0, ry0, rz0, rw1 );
    v = dot_rnd_grad( h = HASH4D( ix+1, iy, iz, iw+1 ), 	rx1, ry0, rz0, rw1 );
    a = LERP(sx, u, v);

    u = dot_rnd_grad( h = HASH4D( ix, 	iy+1, iz, iw+1 ), 	rx0, ry1, rz0, rw1 );
    v = dot_rnd_grad( h = HASH4D( ix+1, iy+1, iz, iw+1 ), 	rx1, ry1, rz0, rw1 );
    b = LERP(sx, u, v);
    c = LERP(sy, a, b);

    u = dot_rnd_grad( h = HASH4D( ix, 	iy, iz+1, iw+1 ), 	rx0, ry0, rz1, rw1 );
    v = dot_rnd_grad( h = HASH4D( ix+1, iy, iz+1, iw+1 ), 	rx1, ry0, rz1, rw1 );
    a = LERP(sx, u, v);

    u = dot_rnd_grad( h = HASH4D( ix, 	iy+1, iz+1, iw+1 ), 	rx0, ry1, rz1, rw1 );
    v = dot_rnd_grad( h = HASH4D( ix+1, iy+1, iz+1, iw+1 ), 	rx1, ry1, rz1, rw1 );
    b = LERP(sx, u, v);
    d = LERP(sy, a, b);
    f = LERP(sz, c, d);

    res = LERP(sw, e, f);

    res.store(out+iPoint);
  }
}

static void curlNoise( 
		float *pos,  	// IN: position (4D)  in grid coords, time_steps
		float freq,   // IN: frequencies (3D x 4D)
		float freq_tm,
		float seed,
		float *velOut   
    	)
{
  float x=pos[0], y=pos[1], z=pos[2], w=pos[3],
#if 0
        h = 0.1f;
#else
//        h = 0.1*1.0/freq;
        h = 0.1*1.0/freq;
#endif

  float posX[32] = 
  {	x,	x, 	x,	x, 	x,	x, 	x+h,	x-h,	x+h,	x-h,	x, 	x };
  float posY[32] = 
  {	y+h,	y-h, 	y,	y, 	y,	y, 	y,	y,	y,	y,	y+h, 	y-h };
  float posZ[32] = 
  {	z,	z, 	z+h,	z-h, 	z+h,	z-h, 	z,	z,	z,	z,	z, 	z };
  const float compOut[32] = 
    {	2.f,	2.f, 	1.f,	1.f, 	0.f,	0.f, 	2.f,	2.f,	1.f,	1.f,	0.f, 	0.f, 0,0,0,0 };
  float f[32];

  //for(int i=12;i<ARRAY_LENGTH(posX);i++) posX[i] = posY[i] = posZ[i] = 0.0f;

  curlNoisePotential4D( 12, posX, posY, posZ, w, 
				   compOut, 
				   freq,freq,freq,freq_tm, 
				   seed,
				   f );

  double scale = 1.0f/(2.0*h*freq);

#if 1
  velOut[0] = ( ((double)f[0] - (double)f[1]) - ( ((double)f[2] - (double)f[3]) )  ) * scale;
  velOut[1] = ( ((double)f[4] - (double)f[5]) - ( ((double)f[6] - (double)f[7]) )  ) * scale;
  velOut[2] = ( ((double)f[8] - (double)f[9]) - ( ((double)f[10] - (double)f[11]) )  ) * scale;
#else
  velOut[0] = f[10];
  velOut[1] = f[9];
  velOut[2] = f[7];
#endif
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void curlNoiseSum( 
		float *pos,  	// IN: position (4D)
		float freq,   // IN: frequency
		float freq_tm,
		float alpha,
		float beta,
		float seed,
		float octaves,     // IN: ocatves 
		int oct_leadin,
		float *ampl,
		float *ampl_df,  // divergence free amplitude factor
		float ampl_hf_gain, // high frequency boost
		float ampl_hf_min_oct,
		float *out	// OUT: result vector (3D)
    	)
{
  float rem,scale=1;
  double sum[3]={0,0,0};
  float val[3];

  octaves += oct_leadin;

  int n=(int)octaves;

  int ampl_hf_i0 = ampl_hf_gain>MT_NUM_EPS ? floor(ampl_hf_min_oct) : INT_MAX;

  if((rem = octaves - n) > 1.0E-12) n++;

  for(int i=0;i<n;i++) 
  { 
#if 1
    curlNoise( pos, freq, freq_tm, seed, val );
#else
    {
      float poff = seed;
      Vec4f pos_ = poff + Vec4f( freq*pos[0], freq*pos[1], freq*pos[2], freq_tm*pos[3] ); 
      Vec4f v = CurlNoise2::curlNoise4D( pos_ );
      val[0] = vfget_x(v);
      val[1] = vfget_y(v);
      val[2] = vfget_z(v);
    }
#endif
    freq *= beta;
    freq_tm *= beta;

    for(int k=0;k<3;k++)
    {
      if(i>=(int)octaves) val[k] *= rem;
      float scaleAct = scale;
      if(i>=ampl_hf_i0)
      {
	scaleAct *= ampl_hf_gain;
	if(i==ampl_hf_i0) 
	{
	  float rem = 1.0f - (ampl_hf_min_oct - ampl_hf_i0);
	  scaleAct *= rem;
	}
      }
      sum[k] += val[k] * scaleAct;
    }
    if(!(i<oct_leadin)) scale /= alpha;
  }

  for(int k=0;k<3;k++)
  {
    out[k] = ampl[k] * sum[k];
  }
}


}  // namespace


