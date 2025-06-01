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
#include "pnoise.h"
#include "plot.h"
#include "noise.h"
#include "bitmap.h"
#include "panel.h"

using namespace MSBG_NAMESPACE;

namespace CurlNoise2
{
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

#define s_curve2(t) (t * t * t * (t * (t * 6.f - 15.f) + 10.f))
#define s_curve2_deriv(t) (30.f * t * t * (t * (t - 2.f) + 1.f)); 
 
//
//	
//
template<int DO_RET_POT_VAL, int S_CURVE_TYP /*1=linear,3=quintic*/, int DO_DEBUG >
Vec4f curlNoise4D_( Vec4f pos, Vec4f *pPotValue )
{
  // 
  // Curl noise based on fast hash & analytical noise derivatives
  //

  if(!(vall(abs(pos)<2000000.f)))  // single precision stability  TODO: UT_ASSERT2
  {
    #if 0
    static int numTraced=0;
    if(numTraced++ <=100)
    #endif
    {
      //xxx ERROR: pos = 43725.2 11656.2 124355 666.999
      TRCERR(("pos = %g %g %g %g\n",pos[0],pos[1],pos[2],pos[3]));
      UT_ASSERT0(FALSE);
    }
  }

  Vec4f posf = floor(pos);
  Vec4i ipos = truncate_to_int(posf);			
  Vec4f r = pos-posf,
	r1 = r-1.0f;

  // Generate a seed hash for each of the 16 cell corners
  Vec8ui hash0,hash1;
  {
    const Vec4i hashmul( 73856093, 19349663, 83492791, 56598313 );
    Vec8i h = Vec8i(ipos,ipos+1) * Vec8i( hashmul, hashmul );

    //	0	1	2	3	4	5	6	7
    // 	x	y	z	w	x1	y1	z1	w1
    //--------------------------------------------------------------
    // hashW0:
    //	x	x1	x	x1	x	x1	x	x1
    //	y	y	y1	y1	y	y	y1	y1
    //	z	z	z	z	z1	z1	z1	z1
    //	w	w	w	w	w	w	w	w
    // hashW1:
    //	x	x1	x	x1	x	x1	x	x1
    //	y	y	y1	y1	y	y	y1	y1
    //	z	z	z	z	z1	z1	z1	z1
    //	w1	w1	w1	w1	w1	w1	w1	w1

    //IACA_START
    Vec8i hx = permute8<0,4,0,4,0,4,0,4>(h),
	  hy = permute8<1,1,5,5,1,1,5,5>(h),
	  hz = permute8<2,2,2,2,6,6,6,6>(h),
	  hw0 = permute8<3,3,3,3,3,3,3,3>(h),
	  hw1 = permute8<7,7,7,7,7,7,7,7>(h);
    //IACA_END

    Vec8i hxyz = hx^hy^hz;
    hash0 = (Vec8ui)(hxyz^hw0); 
    hash1 = (Vec8ui)(hxyz^hw1); 
  }

  // 16 random gradients (corner of 4D cube) for each of the 3 vector value components
  Vec8f gx0[3],gy0[3],gz0[3],gw0[3],
        gx1[3],gy1[3],gz1[3],gw1[3];

  for(int k=0;k<3;k++)
  {
    gx0[k] = RNDF(hash0); gx1[k] = RNDF(hash1);	
    gy0[k] = RNDF(hash0); gy1[k] = RNDF(hash1);
    gz0[k] = RNDF(hash0); gz1[k] = RNDF(hash1);
    gw0[k] = RNDF(hash0); gw1[k] = RNDF(hash1);
    NORMALIZE( gx0[k], gy0[k], gz0[k], gw0[k] );
    NORMALIZE( gx1[k], gy1[k], gz1[k], gw1[k] );
  }

  // Calculate the dot products
  Vec8f r0r1( r, r1 );
  Vec8f rx = permute8<0,4,0,4,0,4,0,4>(r0r1),
	ry = permute8<1,1,5,5,1,1,5,5>(r0r1),
	rz = permute8<2,2,2,2,6,6,6,6>(r0r1),
  	rw0 = permute8<3,3,3,3,3,3,3,3>(r0r1),
	rw1 = permute8<7,7,7,7,7,7,7,7>(r0r1);

  Vec8f dot0[3],dot1[3];

  for(int k=0;k<3;k++)
  {
    dot0[k] = gx0[k]*rx + gy0[k]*ry + gz0[k]*rz + gw0[k] * rw0,
    dot1[k] = gx1[k]*rx + gy1[k]*ry + gz1[k]*rz + gw1[k] * rw1;
  }

  // Interpolation weights
  Vec4f s;

  if constexpr ( S_CURVE_TYP==3 ) 
  {
    s = s_curve2(r);
  }
  else if constexpr( S_CURVE_TYP==2 )
  {
    s = ( r * r * (3.f - 2.f * r) );
  }
  else 
  {
    s = r;
  }
       	
  Vec8f s1s0 = Vec8f(1.0f-s,s);
  Vec8f sx = permute8<0,4,0,4,0,4,0,4>(s1s0),
	sy = permute8<1,1,5,5,1,1,5,5>(s1s0),
	sz = permute8<2,2,2,2,6,6,6,6>(s1s0),
	sw0 = permute8<3,3,3,3,3,3,3,3>(s1s0),
	sw1 = permute8<7,7,7,7,7,7,7,7>(s1s0);    
  /*
  sx0 : 1-s(x0);  sx1 : 1-sx0;
  sy0 : 1-s(y0);  sy1 : 1-sy0;
  sz0 : 1-s(z0);  sz1 : 1-sz0;
  sw0 : 1-s(w0);  sw1 : 1-sw0;

  f_(x0,y0,z0,w0) :=
    sx0 * sy0 * sz0 * sw0 * (gx0000 * x0 + gy0000 * y0 + gz0000 * z0 + gw0000 * w0)		
  + sx1 * sy0 * sz0 * sw0 * (gx1000 * x1 + gy1000 * y0 + gz1000 * z0 + gw1000 * w0)
  + sx0 * sy1 * sz0 * sw0 * (gx0100 * x0 + gy0100 * y1 + gz0100 * z0 + gw0100 * w0) 
  + sx1 * sy1 * sz0 * sw0 * (gx1100 * x1 + gy1100 * y1 + gz1100 * z0 + gw1100 * w0)
  + sx0 * sy0 * sz1 * sw0 * (gx0010 * x0 + gy0010 * y0 + gz0010 * z1 + gw0010 * w0)
  + sx1 * sy0 * sz1 * sw0 * (gx1010 * x1 + gy1010 * y0 + gz1010 * z1 + gw1010 * w0)
  + sx0 * sy1 * sz1 * sw0 * (gx0110 * x0 + gy0110 * y1 + gz0110 * z1 + gw0110 * w0)
  + sx1 * sy1 * sz1 * sw0 * (gx1110 * x1 + gy1110 * y1 + gz1110 * z1 + gw1110 * w0)

  + sx0 * sy0 * sz0 * sw1 * (gx0001 * x0 + gy0001 * y0 + gz0001 * z0 + gw0001 * w1)		
  + sx1 * sy0 * sz0 * sw1 * (gx1001 * x1 + gy1001 * y0 + gz1001 * z0 + gw1001 * w1)
  + sx0 * sy1 * sz0 * sw1 * (gx0101 * x0 + gy0101 * y1 + gz0101 * z0 + gw0101 * w1)
  + sx1 * sy1 * sz0 * sw1 * (gx1101 * x1 + gy1101 * y1 + gz1101 * z0 + gw1101 * w1)
  + sx0 * sy0 * sz1 * sw1 * (gx0011 * x0 + gy0011 * y0 + gz0011 * z1 + gw0011 * w1)
  + sx1 * sy0 * sz1 * sw1 * (gx1011 * x1 + gy1011 * y0 + gz1011 * z1 + gw1011 * w1)
  + sx0 * sy1 * sz1 * sw1 * (gx0111 * x0 + gy0111 * y1 + gz0111 * z1 + gw0111 * w1)
  + sx1 * sy1 * sz1 * sw1 * (gx1111 * x1 + gy1111 * y1 + gz1111 * z1 + gw1111 * w1);

  d/dx0 =  //sx0_ := diff(1-s(x0),x0); sx1_ := diff(s(x),x); 
   sx0_ * sy0 * sz0 * sw0 * (gx0000 * x0 + gy0000 * y0 + gz0000 * z0 + gw0000 * w0) + sx0 * sy0 * sz0 * sw0 * gx0000 +
   sx1_ * sy0 * sz0 * sw0 * (gx1000 * x1 + gy1000 * y0 + gz1000 * z0 + gw1000 * w0) + sx1 * sy0 * sz0 * sw0 * gx1000 +

   sx0_ * sy1 * sz0 * sw0 * (gx0100 * x0 + gy0100 * y1 + gz0100 * z0 + gw0100 * w0) + sx0 * sy1 * sz0 * sw0 * gx0100 +
   ...

  */

  Vec4f ds;

  if constexpr ( S_CURVE_TYP==3 )
  {
    ds = s_curve2_deriv(r);    
  }
  else if constexpr( S_CURVE_TYP==2 )
  {
    ds = 6.f*(1.f-r)*r;
  }
  else
  {
    ds = 0.0f;
  }

  Vec8f ds1s0 = Vec8f(-ds,ds),
   	dsx = permute8<0,4,0,4,0,4,0,4>(ds1s0),
   	dsy = permute8<1,1,5,5,1,1,5,5>(ds1s0),
   	dsz = permute8<2,2,2,2,6,6,6,6>(ds1s0);

  Vec8f sxyzw0 = sx*sy*sz*sw0,
	sxyzw1 = sx*sy*sz*sw1;



  // Note: we rely on compiler for optimizing common sub expressions

  #define df_dx(k) \
    ((double)horizontal_add(dsx*sy*sz*sw0 * dot0[k]) + (double)horizontal_add(sxyzw0 * gx0[k]) +\
     (double)horizontal_add(dsx*sy*sz*sw1 * dot1[k]) + (double)horizontal_add(sxyzw1 * gx1[k]))

  #define df_dy(k) \
    ((double)horizontal_add(sx*dsy*sz*sw0 * dot0[k]) + (double)horizontal_add(sxyzw0 * gy0[k]) +\
     (double)horizontal_add(sx*dsy*sz*sw1 * dot1[k]) + (double)horizontal_add(sxyzw1 * gy1[k]))

  #define df_dz(k) \
    ((double)horizontal_add(sx*sy*dsz*sw0 * dot0[k]) + (double)horizontal_add(sxyzw0 * gz0[k]) +\
     (double)horizontal_add(sx*sy*dsz*sw1 * dot1[k]) + (double)horizontal_add(sxyzw1 * gz1[k]))

  if constexpr( DO_DEBUG )
  {
    return Vec4f( df_dx(0), df_dy(0), df_dz(0), 0.0f );  	
  }

  float vx = df_dy(2) - df_dz(1),
        vy = df_dz(0) - df_dx(2),
        vz = df_dx(1) - df_dy(0);

  #undef df_dx
  #undef df_dy
  #undef df_dz

  if constexpr (DO_RET_POT_VAL )
  {
    #define V(k) (double)horizontal_add(sxyzw0 * dot0[k]) + \
                 (double)horizontal_add(sxyzw1 * dot1[k])

    *pPotValue = Vec4f( V(0),V(1),V(2), 0.0f );
    
    #undef V
  }

  return Vec4f( vx, vy, vz, 0.0f );  	
}



template Vec4f curlNoise4D_<TRUE>( Vec4f pos, Vec4f *pPotValue );
template Vec4f curlNoise4D_<FALSE>( Vec4f pos, Vec4f *pPotValue );
//template Vec4f curlNoise4D_<FALSE,1>( Vec4f pos, Vec4f *pPotValue );
template Vec4f curlNoise4D_<FALSE,2>( Vec4f pos, Vec4f *pPotValue );
template Vec4f curlNoise4D_<TRUE,2>( Vec4f pos, Vec4f *pPotValue );
template Vec4f curlNoise4D_<FALSE,2,TRUE>( Vec4f pos, Vec4f *pPotValue );
template Vec4f curlNoise4D_<FALSE,3,TRUE>( Vec4f pos, Vec4f *pPotValue );
//template Vec4f curlNoise4D_<TRUE,3,TRUE>( Vec4f pos, Vec4f *pPotValue );
template Vec4f curlNoise4D_<true,1>( Vec4f pos, Vec4f *pPotValue );
template Vec4f curlNoise4D_<false,1>( Vec4f pos, Vec4f *pPotValue );

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
Vec4f curlNoise4D( Vec4f pos ) 
{
  return curlNoise4D_<0>(pos,NULL);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
Vec4f curlNoiseSum4D( 
		const Vec4f& pos0,  	 // position (4D)
		float freq_sp,   // frequency (spatial)
		float freq_tm,	 // frequency (time)
		float alpha,
		float beta_sp,
		float beta_tm,
		float seed,
		float octaves,     // IN: ocatves 
		int oct_leadin,
		int lf_min_oct,	// optional cutoff for LF component (out_lf)
		const Vec4f& ampl,
		float ampl_hf_gain, // high frequency boost
		float ampl_hf_min_oct,
		const Vec4f *modfield,  // Modulating scalar field and derivatives f,fx,fy,fz
		unsigned options,
		Vec4f *out_lf  // Optional out: LF filtered (
    	)
{
  //TRCP(("%s: ampl=%g %g %g %g\n",UT_FUNCNAME,ampl[0],ampl[1],ampl[2],ampl[3]));

  bool doScaleAmplHfTimeFreq = options & 2,
       doReturnPotential = options & 4;

  bool doTraceOctaves = options & 8;
  if(doTraceOctaves)
  {
    Vec4ui ipos = _mm_castps_si128(pos0);  // re-interpret as 32 bit integer
    FMA_FastRandSeed(&ipos);
    int32_t rnds = UtHash32Coords4D16(ipos);
    float u = UtFastRand(&rnds);
    doTraceOctaves = u < 1./100000.;
  }

  float rem,scale=1;
  Vec4d sum(0.),
	sum_lf(0.),
	sum_pot(0.);
  Vec4f val(0.f),val_pot(0.0f);
  Vec4d pos = _mm256_cvtps_pd(pos0);

  Vec4d beta( beta_sp, beta_sp, beta_sp, beta_tm ),
	freq( freq_sp, freq_sp, freq_sp, freq_tm );
 
  octaves += oct_leadin;

  if(out_lf)
  {
    //UT_ASSERT0(FALSE); // XXX TODO: sum_pot_lf
  }

  //TRCP(("octaves=%g lf_min_oct=%d\n",octaves,lf_min_oct));

  int n=(int)octaves;

  int ampl_hf_i0 = ampl_hf_gain>MT_NUM_EPS ? floor(ampl_hf_min_oct) + oct_leadin : INT_MAX;

  if((rem = octaves - n) > 1.0E-12) n++;

  if(!out_lf) 
  {
    lf_min_oct = -INT_MAX;
  }
  else lf_min_oct = MIN(lf_min_oct,n-1);


  for(int i=0;i<n;i++) 
  { 
    {
      Vec4f pos_ = _mm256_cvtpd_ps(seed + freq*pos); 

      if(modfield || doReturnPotential)
      {
	val = curlNoise4D_<1>( pos_, &val_pot );
      }
      else
      {
        #if 1
        val = curlNoise4D_<0>( pos_, NULL );
        #else
        val = curlNoise4D_<1>( pos_, &val_pot );
	val = val_pot;
        #endif
      }

      // Check v == curl(pot)
      #if 0
      {
	double h=0.01;
	Vec4f f1,f2;

        curlNoise4D_<TRUE>( pos_-Vec4f(0,h,0,0), &f1 ); 
        curlNoise4D_<TRUE>( pos_+Vec4f(0,h,0,0), &f2 );
	double  fzy = (f2[2] - f1[2])/(2*h); 
        curlNoise4D_<TRUE>( pos_-Vec4f(0,0,h,0), &f1 ); 
        curlNoise4D_<TRUE>( pos_+Vec4f(0,0,h,0), &f2 );
	double  fyz = (f2[1] - f1[1])/(2*h);
	double vx = fzy - fyz;

        curlNoise4D_<TRUE>( pos_-Vec4f(0,0,h,0), &f1 ); 
        curlNoise4D_<TRUE>( pos_+Vec4f(0,0,h,0), &f2 );
	double  fxz = (f2[0] - f1[0])/(2*h); 
        curlNoise4D_<TRUE>( pos_-Vec4f(h,0,0,0), &f1 ); 
        curlNoise4D_<TRUE>( pos_+Vec4f(h,0,0,0), &f2 );
	double  fzx = (f2[2] - f1[2])/(2*h);
	double vy = fxz - fzx;

        curlNoise4D_<TRUE>( pos_-Vec4f(h,0,0,0), &f1 ); 
        curlNoise4D_<TRUE>( pos_+Vec4f(h,0,0,0), &f2 );
	double  fyx = (f2[1] - f1[1])/(2*h); 
        curlNoise4D_<TRUE>( pos_-Vec4f(0,h,0,0), &f1 ); 
        curlNoise4D_<TRUE>( pos_+Vec4f(0,h,0,0), &f2 );
	double  fxy = (f2[0] - f1[0])/(2*h);
	double vz = fyx - fxy;

	TRCP(("%g=%g %g=%g %g=%g\n", vfget_x(v),vx,vfget_y(v),vy, vfget_z(v),vz));        	
      }
      #endif

      // Check div(v) == 0
      #if 0
      {
	float h=0.001*freq;

	Vec4f v[6];
	Vec4f dir[3] = { Vec4f(h,0,0,0), Vec4f(0,h,0,0), Vec4f(0,0,h,0) };

	float velMagAvg=0;
	for(int i=0;i<3;i++)
	{
	  Vec4f v0 = curlNoise4D_<0>(pos_-dir[i]),
	  	v1 = curlNoise4D_<0>(pos_+dir[i]);
	  velMagAvg += sqrtf(v3normSq(0.5f*(v0+v1)));
	  v[2*i+0] = v0;
	  v[2*i+1] = v1;
	}
	velMagAvg *= 1.f/3.f;

	float vxx = ( v[1][0] - (double)v[0][0] ) / (double)(2.f*h),
	      vyy = ( v[3][1] - (double)v[2][1] ) / (double)(2.f*h),
	      vzz = ( v[5][2] - (double)v[4][2] ) / (double)(2.f*h);

	float div = vxx + vyy + vzz;

	//if(div>0.001) 
	{
	  TRCP(("div=%g\n",div));
	}
      }
      #endif
    }

    float gain =1.0f;

    if(i>=(int)octaves) val *= rem;

    float scaleAct=scale;
    if(i>=ampl_hf_i0)
    {
      gain = ampl_hf_gain;
      UT_ASSERT0(!(gain<1.0f));
      if(i==ampl_hf_i0) 
      {
	float rem = 1.0f - (ampl_hf_min_oct - ampl_hf_i0);
	gain = 1 + rem*(gain-1);
        if( doScaleAmplHfTimeFreq )
	{
	  freq *= Vec4d(1,1,1,1./ampl_hf_gain);
	}
      }

      scaleAct *= gain;
    }

#if 0  // TEST
    if( i==n-1 ) scaleAct *= 4; // boost last octave => random
#endif

    Vec4d valAct = _mm256_cvtps_pd( val ) * scaleAct;

    #if 0  // Verify 'alpha' amplitude decay (from energy spectrum)
    {
      static double oct[100],n[100];
      static int iDisp;
      double e = valAct[0]*valAct[0] + valAct[1]*valAct[1] + valAct[2]*valAct[2]; 
      oct[i] += sqrt(e);
      n[i] += 1.0f;
      if((iDisp++)%100000 == 0)
      {
	for(int j=0;j<10;j++)
	{
	  float u = n[j]>0?oct[j]/n[j] : 0,
		uPrev = j>0 ? ( n[j-1]>0?oct[j-1]/n[j-1]:0 ) : 0;
	  TRCP(("i=%d u=%g -> alpha ~ %g (rem=%g)\n",j, u, u>0?uPrev/u:0,rem));
	}
	TRCP(("------\n"));
      }
    }
    #endif

    sum += valAct;

    if(modfield||doReturnPotential) sum_pot += _mm256_cvtps_pd( val_pot ) * scaleAct;

    if(i==lf_min_oct)
    {
      sum_lf = sum;
    }

    if( doTraceOctaves )
    {
      TRCP(("octave %d/%d(%d) -> scale=%g ampl=%g/%g gain=%g\n",
	    i,n,ampl_hf_i0,1.0f/freq[0],scale,scaleAct,gain));
    }

    if(!(i<oct_leadin)) scale /= alpha;

    freq *= beta;
  }  // octaves

  // 
  // Apply spatialy varying modulation field a(x,y,z)
  //
  // curl( a * psi ) = a * curl(psi) + grad(f) X psi        
  //

  if(modfield)
  {
    Vec3d sum_curl_pot = sum,
	  sum_val_pot = sum_pot,
	  modf_grad( vfget_y(*modfield), vfget_z(*modfield), vfget_w(*modfield) );
    double modf = vfget_x(*modfield);
    sum_curl_pot = modf * sum_curl_pot + cross_product( modf_grad, sum_val_pot );
    sum = sum_curl_pot;
  }

  if(out_lf)
  {
    *out_lf = v4f_zero_to<3>(ampl * (Vec4f)_mm256_cvtpd_ps(sum_lf));
  }

  return 
    v4f_zero_to<3>(ampl * 
	(Vec4f)_mm256_cvtpd_ps(doReturnPotential ? sum_pot : sum));
}

} // namespace

#if 0
#if 1
#define _B_ 0x1000
#define BM 0xfff
#else
#define _B_ 0x100
#define BM 0xff
#endif

#ifdef PNS3_DOUBLE_PREC
  #ifdef MI_WITH_64BIT
  #define N 0x80000000
  #else
  #define N 0x10000000
  #endif
#else
  #define N 0x1000
#endif

#ifdef MI_WITH_64BIT
#define CORR_OFFSET_X 12345678.98765	
#define CORR_OFFSET_Y 98765432.12345
#define CORR_OFFSET_Z 31417549.14133
#else
#define CORR_OFFSET_X 123.98765	
#define CORR_OFFSET_Y 987.12345
#define CORR_OFFSET_Z 314.14133
#endif

#define s_curve(t) ( t * t * (3. - 2. * t) )

#define s_curve2(t) (t * t * t * (t * (t * 6 - 15) + 10))
#define s_curve2_d(t) 30 * t * t * (t - 1) * (t - 1)


#define LERP(t, a, b) ( a + t * (b - a) )
#define setup(i,b0,b1,r0,r1)\
        t = vec[i] + N;\
        b0 = ((pns_int_t)t) & BM;\
        b1 = (b0+1) & BM;\
        r0 = t - (pns_int_t)t;\
        r1 = r0 - 1.;
#define at2(rx,ry) ( rx * q[0] + ry * q[1] )
#define at3(rx,ry,rz) ( rx * q[0] + ry * q[1] + rz * q[2] )
#define at4(rx,ry,rz,rw) ( rx * q[0] + ry * q[1] + rz * q[2] + rw * q[3] )

static pns_int_t p[_B_ + _B_ + 2];
static pns_grad_t g4[_B_ + _B_ + 2][4];

template< int deriv1Idx, int deriv2Idx >
float noise4d_derivs(float vec[4], // position
    		     float *deriv1,  // partial derivatives
		     float *deriv2
    		       )
{
   pns_int_t bx0, bx1, by0, by1, bz0, bz1, bw0, bw1;
   pns_grad_t rx0, rx1, ry0, ry1, rz0, rz1, rw0, rw1, 
	      a, b, c, d, e, t, u, v, f;
   pns_int_t i, j;
   pns_grad_t *q;

//   IACA_START
   setup(0, bx0,bx1, rx0,rx1);
   setup(1, by0,by1, ry0,ry1);
   setup(2, bz0,bz1, rz0,rz1);
   setup(3, bw0,bw1, rw0,rw1);

   i = p[ bx0 ];
   j = p[ bx1 ];

   int 
   b00 = p[ i + by0 ],
   b10 = p[ j + by0 ],
   b01 = p[ i + by1 ],
   b11 = p[ j + by1 ],
   b000 = p[ b00 + bz0 ],
   b100 = p[ b10 + bz0 ],
   b010 = p[ b01 + bz0 ],
   b110 = p[ b11 + bz0 ],
   b001 = p[ b00 + bz1 ],
   b101 = p[ b10 + bz1 ],
   b011 = p[ b01 + bz1 ],
   b111 = p[ b11 + bz1 ];



   // 
   // maxima optimize(matrix([diff(f(x,y,z,w),x)
   //

#if 0
pns_grad_t   
h1=x^3,h2=12*x-15,h3=X0-1,h4=gz1000*z0\
-gz0000*z0+gy1000*y0-gy0000*y0-gx0000*x0+gw1000*w0-gw0000*w0+h3*gx1000,h5=x^2,\
h6=x*(6*x-15)+10,h7=y^3,h8=y*(6*y-15)+10,h9=-h1*h2*h4,h10=-3*h5*h6*h4,h11=y0-1\
,h12=gz1100*z0-gz0100*z0+gy1100*h11-gy0100*h11-gx0100*x0+gw1100*w0-gw0100*w0+h\
3*gx1100,h13=3*h5*h6*h12+h1*h2*h12+h10+h9,h14=z^3,h15=z*(6*z-15)+10,h16=z0-1,h\
17=gz1010*h16-gz0010*h16+gy1010*y0-gy0010*y0-gx0010*x0+gw1010*w0-gw0010*w0+h3*\
gx1010,h18=gz1110*h16-gz0110*h16+gy1110*h11-gy0110*h11-gx0110*x0+gw1110*w0-gw0\
110*w0+h3*gx1110,h19=-h7*h8*h13,h20=h19+h10+h9+3*h5*h6*h17+h1*h2*h17+h7*h8*(3*\
h5*h6*h18+h1*h2*h18-3*h5*h6*h17-h1*h2*h17),h21=w0-1,h22=gz1001*z0-gz0001*z0+gy\
1001*y0-gy0001*y0-gx0001*x0+gw1001*h21-gw0001*h21+h3*gx1001,h23=-h1*h2*h22,h24\
=-3*h5*h6*h22,h25=gz1101*z0-gz0101*z0+gy1101*h11-gy0101*h11-gx0101*x0+gw1101*h\
21-gw0101*h21+h3*gx1101,h26=3*h5*h6*h25+h1*h2*h25+h24+h23,h27=gz1011*h16-gz001\
1*h16+gy1011*y0-gy0011*y0-gx0011*x0+gw1011*h21-gw0011*h21+h3*gx1011,h28=gz1111\
*h16-gz0111*h16+gy1111*h11-gy0111*h11-gx0111*x0+gw1111*h21-gw0111*h21+h3*gx111\
1,w^3*(w*(6*w-15)+10)*(h14*h15*((-h7*h8*h26)+h24+h23+3*h5*h6*h27+h1*h2*h27+h7*\
h8*(3*h5*h6*h28+h1*h2*h28-3*h5*h6*h27-h1*h2*h27))-h14*h15*h20+h7*h8*h26+h19+3*\
h5*h6*h22+h1*h2*h22+h10+h9)+h14*h15*h20+h7*h8*h13+3*h5*h6*h4+h1*h2*h4;
#endif

   //xxx SIMD Vec4d
   
   pns_grad_t 
   sx = s_curve2(rx0),
   sy = s_curve2(ry0),
   sz = s_curve2(rz0),
   sw = s_curve2(rw0);

   pns_grad_t 
   sxd = s_curve2_d(rx0),
   syd = s_curve2_d(ry0),
   szd = s_curve2_d(rz0),
   swd = s_curve2_d(rw0);
   
   q = g4[ b000 + bw0 ] ; u = at4(rx0,ry0,rz0,rw0);   
   q = g4[ b100 + bw0 ] ; v = at4(rx1,ry0,rz0,rw0);
   a = LERP(sx, u, v);
   q = g4[ b010 + bw0 ] ; u = at4(rx0,ry1,rz0,rw0);
   q = g4[ b110 + bw0 ] ; v = at4(rx1,ry1,rz0,rw0);
   b = LERP(sx, u, v);
   c = LERP(sy, a, b);
   q = g4[ b001 + bw0 ] ; u = at4(rx0,ry0,rz1,rw0);
   q = g4[ b101 + bw0 ] ; v = at4(rx1,ry0,rz1,rw0);
   a = LERP(sx, u, v);
   q = g4[ b011 + bw0 ] ; u = at4(rx0,ry1,rz1,rw0);
   q = g4[ b111 + bw0 ] ; v = at4(rx1,ry1,rz1,rw0);
   b = LERP(sx, u, v);
   d = LERP(sy, a, b);
   e = LERP(sz, c, d);

   q = g4[ b000 + bw1 ] ; u = at4(rx0,ry0,rz0,rw1);   
   q = g4[ b100 + bw1 ] ; v = at4(rx1,ry0,rz0,rw1);
   a = LERP(sx, u, v);
   q = g4[ b010 + bw1 ] ; u = at4(rx0,ry1,rz0,rw1);
   q = g4[ b110 + bw1 ] ; v = at4(rx1,ry1,rz0,rw1);
   b = LERP(sx, u, v);
   c = LERP(sy, a, b);
   q = g4[ b001 + bw1 ] ; u = at4(rx0,ry0,rz1,rw1);
   q = g4[ b101 + bw1 ] ; v = at4(rx1,ry0,rz1,rw1);
   a = LERP(sx, u, v);
   q = g4[ b011 + bw1 ] ; u = at4(rx0,ry1,rz1,rw1);
   q = g4[ b111 + bw1 ] ; v = at4(rx1,ry1,rz1,rw1);
   b = LERP(sx, u, v);
   d = LERP(sy, a, b);
   f = LERP(sz, c, d);

   //IACA_END
   return LERP(sw, e, f);
}
#endif

#if 0
static void normalize4(pns_grad_t v[4])
{
   double s;

   s = sqrt((double)v[0] * (double)v[0] + 
            (double)v[1] * (double)v[1] + 
	    (double)v[2] * (double)v[2] + 
	    (double)v[3] * (double)v[3]);
   v[0] = v[0] / s;
   v[1] = v[1] / s;
   v[2] = v[2] / s;
   v[3] = v[3] / s;
}
#endif

int PNS_Init3(int rseed)
{
#if 0
  UtTimer tm;

  TIMER_START(&tm);

   pns_int_t i, j, k;
   RND_Stream rnds;
 
//   TRCP(("RAND_MAX = %ld\n",(long)RAND_MAX));

   pns_rand_init(rseed,&rnds); 

   for (i = 0 ; i < _B_ ; i++) {
      p[i] = i;		// TODO*: use better PRNG than 'rand()'
   }

   while (--i) {
      k = p[i];
      p[i] = p[j = pns_rand(&rnds) % _B_];
      p[j] = k;
   }

   // do the 4D gradients separately for compatibility of g3 with old seed 
   for (i = 0 ; i < _B_ ; i++) 
   {  
     for (j = 0 ; j < 4 ; j++)
	 g4[i][j] = (double)((pns_rand(&rnds) % (_B_ + _B_)) - _B_) / _B_;
      normalize4(g4[i]);
    }

   // fill wrap around
   for (i = 0 ; i < _B_ + 2 ; i++) {
      p[_B_ + i] = p[i];
      for (j = 0 ; j < 4 ; j++) g4[_B_ + i][j] = g4[i][j];
   }

   PNS_Init2(rseed);

   TIMER_STOP(&tm);

  TRC(("PNS_Init-> CPU = %.0f seconds\n", 
	  (double)(TIMER_DIFF_MS(&tm)/1000.0)));
#endif
   return 0;
}


