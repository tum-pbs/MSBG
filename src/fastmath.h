#ifndef FASTMATH_H
#define FASTMATH_H

#ifdef __cplusplus
#include <float.h>
#include <vectorclass/vectorclass.h>
#include "vectorclass_util.h"
#include "util.h"
#endif

/*=========================================================================
 *
 *
 *
 *
 * =======================================================================*/

#ifdef __cplusplus

inline bool FMA_isfinite( float x )
{
  return std::isfinite(x);
}

inline bool FMA_isfinite( const Vec4f& x )
{
  return visfinite(x);
}

inline int FMA_FAST_ROUND_LOG2(double x) 
{
  union
  {
      double f;
      uint64_t i;
  } u;
  u.f = x*1.4142135623;
  u.i = ((u.i >> 52) & 0x7ff) - 1023;
  return u.i;
}

inline int FMA_FAST_INT_LOG2(double x) 
{
  union
  {
      double f;
      uint64_t i;
  } u;
  u.f = x;
  u.i = ((u.i >> 52) & 0x7ff) - 1023;
  return u.i;
}

inline Vec4f FMA_FastRand( Vec4ui* rnds )
{
  // http://www.iquilezles.org/www/articles/sfrand/sfrand.htm
#if 0
  (*rnds) *= 16807; // http://www.iquilezles.org/www/articles/sfrand/sfrand.htm
  Vec4f u = _mm_castsi128_ps((*rnds)>>9|0x3f800000); 
  return u-1.0f;
#else
  Vec4f u = _mm_castsi128_ps((*rnds)>>9|0x3f800000); 
  (*rnds) *= 16807; 
  return u-1.0f;
#endif
}

inline Vec8f FMA_FastRand( Vec8ui* rnds )
{
  Vec8f u = _mm256_castsi256_ps((*rnds)>>9|0x3f800000); 
  (*rnds) *= 16807; 
  return u-1.0f;
}

inline void FMA_FastRandSeed( Vec4ui* rnds )
{
  //UT_ASSERT2( vall((*rnds)<100000));
  Vec4ui h_(135623,34847,25324534,142423);
  h_ += *rnds;
  h_ ^= h_ >> 16;	// 'murmur' hash 
  h_ *= 0x85ebca6b;
  h_ ^= h_ >> 13;
  h_ *= 0xc2b2ae35;
  h_ ^= h_ >> 16;
  h_ = select( h_==0, // ensure not 0 (random integers)
        Vec4ui(659809188,687332320,534213693,107384577),
	h_ );
  *rnds = h_; 
  (void)FMA_FastRand( rnds );
}

inline void FMA_FastRandSeed( Vec8ui* rnds )
{
  Vec8ui h_(135623,34847,25324534,142423,963511,337616,977104,238094);
  h_ += *rnds;
  h_ ^= h_ >> 16;	// 'murmur' hash 
  h_ *= 0x85ebca6b;
  h_ ^= h_ >> 13;
  h_ *= 0xc2b2ae35;
  h_ ^= h_ >> 16;
  h_ = select( h_==0, // ensure not 0 (random integers)
        Vec8ui(659809188,687332320,534213693,107384577,26243257,6254649,9134374,1149582),
	h_ );
  *rnds = h_; 
  (void)FMA_FastRand( rnds );
}

inline float FMA_FastRand( uint32_t* rnds )
{
  // http://www.iquilezles.org/www/articles/sfrand/sfrand.htm
  union
  {
      float fres;
      uint32_t ires;
  } u;

  rnds[0] *= 16807;

  u.ires = ((((unsigned int)rnds[0])>>9 ) | 0x3f800000);
  return u.fres - 1.0f;
}

#if 0
inline float FMA_FastGammaRand( int a  /* assuming a <= 5 */, float b, uint32_t* rnds )
{
  // cf. gsl_ran_gamma_int 
  UT_ASSERT2( a>0 && a<=5 );
  float prod = 1;

  for (int i = 0; i < a; i++)
  {
    float u = FMA_FastRand (rnds); 
    prod *= u;
    //TRCP(("%d/%d>prod=%g\n",i,a,prod));
  }

  prod = vmax( std::numeric_limits<float>::min(), prod );

  float u = -logf (prod) * b;

  UT_ASSERT2( FMA_isfinite(u) );
  return u;
}
#endif

inline float FMA_FastGamma5Rand( uint32_t* rnds )
{
  // Gamma distribution with mean 1 and shape parameter a=5 (scale b=1/a for mu=1)
  // cf. gsl_ran_gamma_int    
  float prod = 1.0f;

  prod *= FMA_FastRand (rnds); 
  prod *= FMA_FastRand (rnds); 
  prod *= FMA_FastRand (rnds); 
  prod *= FMA_FastRand (rnds); 
  prod *= FMA_FastRand (rnds); 

  prod = vmax( std::numeric_limits<float>::min(), prod );
  float u = -logf (prod) * 1.f/5.f;  // a=5, b=1/a 

  UT_ASSERT2( FMA_isfinite(u) );
  return u;
}


template<typename T, typename T_seed>
inline T FMA_FastApproxGaussRand( T_seed* rnds )
{
  // http://c-faq.com/lib/gaussian.html
  T s1=0.0f,s2=0.0f,s3=0.0f,s4=0.0f;  // mult. accumulators
  s1 += FMA_FastRand(rnds);
  s2 += FMA_FastRand(rnds);
  s3 += FMA_FastRand(rnds);
  s4 += FMA_FastRand(rnds);
  s1 += FMA_FastRand(rnds);
  s2 += FMA_FastRand(rnds);
  s3 += FMA_FastRand(rnds);
  s4 += FMA_FastRand(rnds);
  T s = ((s1+s2+s3+s4)-8.f/2.f)*1.224744f;  /* s = (SUM-N/2)/sqrt(N/12) */ 
  UT_ASSERT2( FMA_isfinite( s ) );
  return s;
}

inline Vec4f FMA_FastRandVec3D( Vec4ui* rnds )
{
  // Uniformly distributed directions, uniformly distributed length in [0,1]
  //https://stackoverflow.com/questions/9750908/how-to-generate-a-unit-vector-pointing-in-a-random-direction-with-isotropic-dist
  Vec4f u = FMA_FastApproxGaussRand<Vec4f,Vec4ui>( rnds );
  u = v4f_zero_to<3>(u);
  float mag = sqrtf(v3normSq(u));
  mag = vmax(mag,1e-7f);
  float scale = vfget_x(FMA_FastRand(rnds));
  return u * scale * 1.0f/mag;
}

inline void FMA_FastRandSeed( uint32_t* rnds )
{
  uint32_t h_ = *rnds;
  h_ ^= h_ >> 16;  // 'murmur' hash 
  h_ *= 0x85ebca6b;
  h_ ^= h_ >> 13;
  h_ *= 0xc2b2ae35;
  h_ ^= h_ >> 16;
  if( unlikely(h_==0) ) h_=0xdeadbeef; // ensure not 0
}

inline float FMA_FAST_APPROX_LOG2(float x) 
{
  // Maximum error ~= 7%, average ~=0.2% for x in range 1e-3 to 1e6
  // http://stackoverflow.com/questions/9411823/fast-log2float-x-implementation-c
  union { float f; uint32_t i; } u;
  u.f = x;
  const int  log_2 = ((u.i >>  23) & 255) - 128;
  u.i &= ~(255 << 23);
  u.i += 127 << 23;
  float f = u.f;
  #if 1
  f = ((-1.0f/3.0f) * f + 2.0f) * f - 2.0f/3.0f;  
  #endif
  return f+log_2;
}


// Fast approximate square root: ~11 bit accuracy , ~1.7x faster than mm_sqrt_ss
inline void FMA_FastSqrtf( float *__restrict pOut, 
    			         float *__restrict pIn )
{
   __m128 in = _mm_load_ss( pIn );
   // 1/recip_sqrt: this correctly handles the special case 0  
   // see also: https://forums.developer.nvidia.com/t/implementation-of-sqrt/1792/3
   _mm_store_ss( pOut,  _mm_rcp_ps ( _mm_rsqrt_ss( in ) ) );
}

#if 0  // simple sqrtf (translates to one vsrtss assembly instruction) is faster and more accurate
template<int i_accuracy> float FMA_FAST_APPROX_SQRTF_Q( float x )
{
    union // get bits for floating value
    {
      float x;
      int i;
    } u;
    u.x = x+FLT_MIN;
    float xhalf = 0.5f*x;
    u.i = 0x5f3759df - (u.i >> 1);  // gives initial guess y0
    float y = u.x;
    for(int i=0;i<i_accuracy;i++)
    {
      y = y*(1.5f - xhalf*y*y);
    }
    return y*x;
}
#endif


// Ankerl/Kmett: https://github.com/ekmett/approximate/blob/master/cbits/fast.c
// max. error about 5% (?)
inline float FMA_FAST_APPROX_POW(float a, float b) 
{
  if(a<1e-20f) if(a>=0.0f) return fabs(b)<1e-20f ? 1.0f : 0.0f;
  
  union { float f; int x; } u = { a };

  int flipped = 0;
  if (b < 0) 
  {
    flipped = 1;
    b = -b;
  }

  int e = (int) b;
  float b0 = b-e;

  float r = 1.0f;
  while (e) 
  {
    if (e & 1) r *= a;
    a *= a;
    e >>= 1;
  }

  if(b0)
  {
    u.x = (int)((b0) * (u.x - 1065353216) + 1065353216);
    r *= u.f;
  }

  return flipped ? 1.0f/r : r;
}

namespace FMA
{

inline Vec4uq morton64_split3( const Vec4i& x_ )
{
  UT_ASSERT2( vall(x_ < 0x1fffff) );  // coordinates up to 2097151
  //IACA_START
  Vec4uq x =  _mm256_cvtepi32_epi64(x_);
  x = (x | x << 32) & 0x1f00000000ffff; // shift left 32 bits, OR with self, and 00011111000000000000000000000000000000001111111111111111
  x = (x | x << 16) & 0x1f0000ff0000ff; // shift left 32 bits, OR with self, and 00011111000000000000000011111111000000000000000011111111
  x = (x | x << 8) & 0x100f00f00f00f00f; // shift left 32 bits, OR with self, and 0001000000001111000000001111000000001111000000001111000000000000
  x = (x | x << 4) & 0x10c30c30c30c30c3; // shift left 32 bits, OR with self, and 0001000011000011000011000011000011000011000011000011000100000000
  x = (x | x << 2) & 0x1249249249249249;
  //IACA_END
  return x;
}

inline uint64_t morton64_encode( const Vec4i& pos_ )
{
//  IACA_START
  Vec4uq p = morton64_split3( pos_ );    
  uint64_t p_[4] __attribute__((aligned(32)));
  p.store_a(p_);
  uint64_t key = p_[0] | (p_[1]<<1) | (p_[2]<<2); 
//  IACA_END
  return key;
}

inline Vec4i morton64_decode( uint64_t m )
{
/* Morton encoding in binary (components 21-bit: 0..2097151)
                0zyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyxzyx */
#define BITMASK_0000000001000001000001000001000001000001000001000001000001000001 UINT64_C(18300341342965825)
#define BITMASK_0000001000001000001000001000001000001000001000001000001000001000 UINT64_C(146402730743726600)
#define BITMASK_0001000000000000000000000000000000000000000000000000000000000000 UINT64_C(1152921504606846976)
/*              0000000ccc0000cc0000cc0000cc0000cc0000cc0000cc0000cc0000cc0000cc */
#define BITMASK_0000000000000011000000000011000000000011000000000011000000000011 UINT64_C(844631138906115)
#define BITMASK_0000000111000000000011000000000011000000000011000000000011000000 UINT64_C(126113986927919296)
/*              00000000000ccccc00000000cccc00000000cccc00000000cccc00000000cccc */
#define BITMASK_0000000000000000000000000000000000001111000000000000000000001111 UINT64_C(251658255)
#define BITMASK_0000000000000000000000001111000000000000000000001111000000000000 UINT64_C(1030792212480)
#define BITMASK_0000000000011111000000000000000000000000000000000000000000000000 UINT64_C(8725724278030336)
/*              000000000000000000000000000ccccccccccccc0000000000000000cccccccc */
#define BITMASK_0000000000000000000000000000000000000000000000000000000011111111 UINT64_C(255)
#define BITMASK_0000000000000000000000000001111111111111000000000000000000000000 UINT64_C(137422176256)
/*                                                         ccccccccccccccccccccc */

  const uint64_t  
    mask0 = BITMASK_0000000001000001000001000001000001000001000001000001000001000001,
    mask1 = BITMASK_0000001000001000001000001000001000001000001000001000001000001000,
    mask2 = BITMASK_0001000000000000000000000000000000000000000000000000000000000000,
    mask3 = BITMASK_0000000000000011000000000011000000000011000000000011000000000011,
    mask4 = BITMASK_0000000111000000000011000000000011000000000011000000000011000000,
    mask5 = BITMASK_0000000000000000000000000000000000001111000000000000000000001111,
    mask6 = BITMASK_0000000000000000000000001111000000000000000000001111000000000000,
    mask7 = BITMASK_0000000000011111000000000000000000000000000000000000000000000000,
    mask8 = BITMASK_0000000000000000000000000000000000000000000000000000000011111111,
    mask9 = BITMASK_0000000000000000000000000001111111111111000000000000000000000000;

  Vec4uq x( m, m>>1, m>>2, 0 );

  x = (x & mask0) | ((x & mask1) >> 2) | ((x & mask2) >> 4);
  x = (x & mask3) | ((x & mask4) >> 4);
  x = (x & mask5) | ((x & mask6) >> 8) | ((x & mask7) >> 16);
  x = (x & mask8) | ((x & mask9) >> 16);

  uint64_t p_[4] __attribute__((aligned(32)));
  x.store_a(p_);

  return Vec4i(p_[0],p_[1],p_[2],0);
}



#if 0
static const uint32_t morton3dLUT[256] =
{
    0x00000000,
    0x00000001, 0x00000008, 0x00000009, 0x00000040, 0x00000041, 0x00000048, 0x00000049, 0x00000200,
    0x00000201, 0x00000208, 0x00000209, 0x00000240, 0x00000241, 0x00000248, 0x00000249, 0x00001000,
    0x00001001, 0x00001008, 0x00001009, 0x00001040, 0x00001041, 0x00001048, 0x00001049, 0x00001200,
    0x00001201, 0x00001208, 0x00001209, 0x00001240, 0x00001241, 0x00001248, 0x00001249, 0x00008000,
    0x00008001, 0x00008008, 0x00008009, 0x00008040, 0x00008041, 0x00008048, 0x00008049, 0x00008200,
    0x00008201, 0x00008208, 0x00008209, 0x00008240, 0x00008241, 0x00008248, 0x00008249, 0x00009000,
    0x00009001, 0x00009008, 0x00009009, 0x00009040, 0x00009041, 0x00009048, 0x00009049, 0x00009200,
    0x00009201, 0x00009208, 0x00009209, 0x00009240, 0x00009241, 0x00009248, 0x00009249, 0x00040000,
    0x00040001, 0x00040008, 0x00040009, 0x00040040, 0x00040041, 0x00040048, 0x00040049, 0x00040200,
    0x00040201, 0x00040208, 0x00040209, 0x00040240, 0x00040241, 0x00040248, 0x00040249, 0x00041000,
    0x00041001, 0x00041008, 0x00041009, 0x00041040, 0x00041041, 0x00041048, 0x00041049, 0x00041200,
    0x00041201, 0x00041208, 0x00041209, 0x00041240, 0x00041241, 0x00041248, 0x00041249, 0x00048000,
    0x00048001, 0x00048008, 0x00048009, 0x00048040, 0x00048041, 0x00048048, 0x00048049, 0x00048200,
    0x00048201, 0x00048208, 0x00048209, 0x00048240, 0x00048241, 0x00048248, 0x00048249, 0x00049000,
    0x00049001, 0x00049008, 0x00049009, 0x00049040, 0x00049041, 0x00049048, 0x00049049, 0x00049200,
    0x00049201, 0x00049208, 0x00049209, 0x00049240, 0x00049241, 0x00049248, 0x00049249, 0x00200000,
    0x00200001, 0x00200008, 0x00200009, 0x00200040, 0x00200041, 0x00200048, 0x00200049, 0x00200200,
    0x00200201, 0x00200208, 0x00200209, 0x00200240, 0x00200241, 0x00200248, 0x00200249, 0x00201000,
    0x00201001, 0x00201008, 0x00201009, 0x00201040, 0x00201041, 0x00201048, 0x00201049, 0x00201200,
    0x00201201, 0x00201208, 0x00201209, 0x00201240, 0x00201241, 0x00201248, 0x00201249, 0x00208000,
    0x00208001, 0x00208008, 0x00208009, 0x00208040, 0x00208041, 0x00208048, 0x00208049, 0x00208200,
    0x00208201, 0x00208208, 0x00208209, 0x00208240, 0x00208241, 0x00208248, 0x00208249, 0x00209000,
    0x00209001, 0x00209008, 0x00209009, 0x00209040, 0x00209041, 0x00209048, 0x00209049, 0x00209200,
    0x00209201, 0x00209208, 0x00209209, 0x00209240, 0x00209241, 0x00209248, 0x00209249, 0x00240000,
    0x00240001, 0x00240008, 0x00240009, 0x00240040, 0x00240041, 0x00240048, 0x00240049, 0x00240200,
    0x00240201, 0x00240208, 0x00240209, 0x00240240, 0x00240241, 0x00240248, 0x00240249, 0x00241000,
    0x00241001, 0x00241008, 0x00241009, 0x00241040, 0x00241041, 0x00241048, 0x00241049, 0x00241200,
    0x00241201, 0x00241208, 0x00241209, 0x00241240, 0x00241241, 0x00241248, 0x00241249, 0x00248000,
    0x00248001, 0x00248008, 0x00248009, 0x00248040, 0x00248041, 0x00248048, 0x00248049, 0x00248200,
    0x00248201, 0x00248208, 0x00248209, 0x00248240, 0x00248241, 0x00248248, 0x00248249, 0x00249000,
    0x00249001, 0x00249008, 0x00249009, 0x00249040, 0x00249041, 0x00249048, 0x00249049, 0x00249200,
    0x00249201, 0x00249208, 0x00249209, 0x00249240, 0x00249241, 0x00249248, 0x00249249
 };

inline uint64_t morton64_encode_lut( const Vec4i& pos0 )
{
  //IACA_START
  Vec4ui pos = (__m128i)pos0;
  Vec4ui  k1 = ( pos >> 16 ) & 0xFF,
	  k2 = ( pos >> 8 ) & 0xFF,
	  k3 = pos & 0xFF;

  uint32_t k1_[4] __attribute__((aligned(32))),
	   k2_[4] __attribute__((aligned(32))),
           k3_[4] __attribute__((aligned(32)));

  k1.store_a(k1_);
  k2.store_a(k2_);
  k3.store_a(k3_);


  uint32_t l1 = morton3dLUT[k1_[2]] << 2 |
      		morton3dLUT[k1_[1]] << 1 |
                morton3dLUT[k1_[0]],

           l2 = morton3dLUT[k2_[2]] << 2 |
      		morton3dLUT[k2_[1]] << 1 |
                morton3dLUT[k2_[0]],

           l3 = morton3dLUT[k3_[2]] << 2 |
      		morton3dLUT[k3_[1]] << 1 |
                morton3dLUT[k3_[0]];

  uint64_t key = l1;		

  key = key << 24 | l2;
  key = key << 24 | l3;
  //IACA_END
  return key;
}
#endif

#if 0
inline uint64_t morton64_encode_pdep (const Vec4i& pos)
{
#ifdef  __BMI2__
https://gitlab.mpcdf.mpg.de/mtr/ducc/-/blob/718634a34fa1d3ef515c1bc4bbe20edecb94201b/mr_util/space_filling.h
UT_ASSERT0(FALSE);  // extremely slow (at least on AMD architecture)
return   _pdep_u64(vget_x(pos),0x1249249249249249u)
        |_pdep_u64(vget_y(pos),0x2492492492492492u)
        |_pdep_u64(vget_z(pos),0x4924924924924924u);
#else
return morton64_encode(pos);
#endif
}
#endif

inline uint64_t morton64_encode_ref( const Vec4i& pos_ )
{
  uint64_t 
    x = vget_x(pos_),
    y = vget_y(pos_),
    z = vget_z(pos_),
    answer = 0;
  for (uint64_t i = 0; i < (sizeof(uint64_t)* CHAR_BIT)/3; ++i) 
  {
    answer |= ((x & ((uint64_t)1 << i)) << 2*i) | 
              ((y & ((uint64_t)1 << i)) << (2*i + 1)) | 
	      ((z & ((uint64_t)1 << i)) << (2*i + 2));
  }
  return answer;
}

inline Vec4i morton32_split3( Vec4i x )
{
  //x &= 1023;
  UT_ASSERT2( vall(x < 1024) );
  x = (x | (x << 16)) & 0x030000FF;
  x = (x | (x <<  8)) & 0x0300F00F;
  x = (x | (x <<  4)) & 0x030C30C3;
  x = (x | (x <<  2)) & 0x09249249;
  return x;
}

inline uint32_t morton32_encode( Vec4i pos )
{
  pos = morton32_split3( pos );    
  return vget_x(pos) | (vget_y(pos)<<1) | (vget_z(pos)<<2);
}
  
  const float NUM_EPS_F = 1.0E-12f;


} // namespace FMA


//
// Explicitely vectorized functions 
//

inline 
Vec4f fast_isqrtf_general_Vec4f(Vec4f x, const uint32_t ISQRT_ITERATIONS) 
{
  const Vec4f x2 = x * 0.5f;	
  Vec4f y  = x;
  Vec4i i = _mm_castps_si128(y);
  i  = 0x5f3759df - ( i >> 1 );	
  y =  _mm_castsi128_ps(i);
  for (uint32_t j=0;j<ISQRT_ITERATIONS;++j)
	  y  *= ( 1.5f - ( x2 * y * y ) );
  return y;
}

inline 
Vec4f fast_isqrtf_approx_Vec4f(Vec4f x) 
{
  const Vec4f x2 = x * 0.5f,
	      threehalf(1.5f);
  Vec4f y  = x;
  Vec4i i = _mm_castps_si128(y);
  i  = 0x5f3759df - ( i >> 1 );	
  y =  _mm_castsi128_ps(i);
  y  *= ( threehalf - ( x2 * y * y ) );
  y  *= ( threehalf - ( x2 * y * y ) ); //second iteration
  return y;
}

inline Vec4f FMA_FAST_ISQRT(Vec4f x) {return fast_isqrtf_general_Vec4f(x,2);}

inline Vec4f FMA_FAST_APPROX_ISQRT(Vec4f x) {return fast_isqrtf_approx_Vec4f(x);}

inline Vec4f FMA_FAST_APPROX_SQRT(Vec4f x) 
  {return x*fast_isqrtf_approx_Vec4f(x+FLT_MIN);}


inline 
Vec4f FMA_FAST_EXP( Vec4f initial_x )
{
  const float MAXLOGF = 88.72283905206835f;
  const float MINLOGF = -88.f;

  const float C1F =   0.693359375f;
  const float C2F =  -2.12194440e-4f;

  const float PX1expf = 1.9875691500E-4f;
  const float PX2expf =1.3981999507E-3f;
  const float PX3expf =8.3334519073E-3f;
  const float PX4expf =4.1665795894E-2f;
  const float PX5expf =1.6666665459E-1f;
  const float PX6expf =5.0000001201E-1f;

  const float LOG2EF = 1.44269504088896341f;

  Vec4f x = initial_x;
  Vec4f z = floor( LOG2EF * x +0.5f ); /* floor() truncates toward -infinity. */
  x -= z * C1F;
  x -= z * C2F;
  const Vec4i n = truncate_to_int(z);
  const Vec4f x2 = x * x;
  z = x*PX1expf;
  z += PX2expf;
  z *= x;
  z += PX2expf;
  z *= x;
  z += PX3expf;
  z *= x;
  z += PX4expf;
  z *= x;
  z += PX5expf;
  z *= x;
  z += PX6expf;
  z *= x2;
  z += x + 1.0f;

  /* multiply by power of 2 */
  z *=  _mm_castsi128_ps((n+0x7f)<<23);

  z = select( initial_x < MINLOGF, 0.0f, z ); 
  z = select( initial_x > MAXLOGF, std::numeric_limits<float>::infinity(), z ); 
  return z;
}


inline 
__attribute__((always_inline))  // mingw-gcc bug: wrong stack alignment for AVX  
Vec8f FMA_FAST_EXP( Vec8f initial_x )
{
  const float MAXLOGF = 88.72283905206835f;
  const float MINLOGF = -88.f;

  const float C1F =   0.693359375f;
  const float C2F =  -2.12194440e-4f;

  const float PX1expf = 1.9875691500E-4f;
  const float PX2expf =1.3981999507E-3f;
  const float PX3expf =8.3334519073E-3f;
  const float PX4expf =4.1665795894E-2f;
  const float PX5expf =1.6666665459E-1f;
  const float PX6expf =5.0000001201E-1f;

  const float LOG2EF = 1.44269504088896341f;

  Vec8f x = initial_x;
  Vec8f z = floor( LOG2EF * x +0.5f ); /* floor() truncates toward -infinity. */
  x -= z * C1F;
  x -= z * C2F;
  const Vec8i n = truncate_to_int(z);
  const Vec8f x2 = x * x;
  z = x*PX1expf;
  z += PX2expf;
  z *= x;
  z += PX2expf;
  z *= x;
  z += PX3expf;
  z *= x;
  z += PX4expf;
  z *= x;
  z += PX5expf;
  z *= x;
  z += PX6expf;
  z *= x2;
  z += x + 1.0f;

  /* multiply by power of 2 */
  z *=  _mm256_castsi256_ps((n+0x7f)<<23);

  z = select( initial_x < MINLOGF, 0.0f, z ); 
  z = select( initial_x > MAXLOGF, std::numeric_limits<float>::infinity(), z ); 
  return z;
}

inline Vec4f FMA_FAST_LOG( Vec4f x ) 
{
  const float MAXNUMF = 3.4028234663852885981170418348451692544e38f;
  const float LOGF_UPPER_LIMIT = MAXNUMF;
  const float LOGF_LOWER_LIMIT = 0;
  const Vec4f original_x = x;
  Vec4ui n = _mm_castps_si128(x);  // re-interpret as 32 bit integer
  Vec4i e = (n>>23)-127;
  Vec4f fe =  to_float(e);
  {
    Vec4ui t1 = _mm_castps_si128(x);   // reinterpret as 32-bit integer
    Vec4ui t2 = Vec4ui((t1 & 0x807fffff) | 0x3f000000); // set exponent to 0 + bias
    x = _mm_castsi128_ps(t2);
  }
  const float SQRTHF = 0.707106781186547524f;
  fe = select( x > SQRTHF, fe+1.f, fe );
  x = select( x > SQRTHF, x, x+x );
  x -= 1.0f;
  const Vec4f x2 = x*x;
  Vec4f res;
  {
    const float 
     PX1logf = 7.0376836292E-2f,
     PX2logf = -1.1514610310E-1f,
     PX3logf = 1.1676998740E-1f,
     PX4logf = -1.2420140846E-1f,
     PX5logf = 1.4249322787E-1f,
     PX6logf = -1.6668057665E-1f,
     PX7logf = 2.0000714765E-1f,
     PX8logf = -2.4999993993E-1f,
     PX9logf = 3.3333331174E-1f;

    Vec4f y = x*PX1logf;
    y += PX2logf;
    y *= x;
    y += PX3logf;
    y *= x;
    y += PX4logf;
    y *= x;
    y += PX5logf;
    y *= x;
    y += PX6logf;
    y *= x;
    y += PX7logf;
    y *= x;
    y += PX8logf;
    y *= x;
    y += PX9logf;
    res = y;
  }
  res *= x2*x;
  res += -2.12194440e-4f * fe;
  res +=  -0.5f * x2;
  res= x + res;
  res += 0.693359375f * fe;
  res = select( original_x > LOGF_UPPER_LIMIT, std::numeric_limits<float>::infinity(), res );
  res = select( original_x < LOGF_LOWER_LIMIT, -std::numeric_limits<float>::quiet_NaN(), res );
  return res;
}

inline Vec4f FMA_FAST_POW( Vec4f x, Vec4f y ) 
{
  return FMA_FAST_EXP( y*FMA_FAST_LOG(x));
}

#else

#define FMA_FAST_APPROX_LOG2(x) FMA_FastApproxLog2(x)
#define FMA_FAST_SQRT(x) FMA_FastSqrt(x)
#define FMA_FAST_APPROX_SQRT(x) FMA_FastApproxSqrt(x)

#endif

#ifdef __cplusplus
extern "C" {
#endif

float FMA_FastApproxLog2( float x );
float FMA_FastApproxSqrt( float x );
double FMA_FastSqrt( double x );
double FMA_FastInvSqrt( double x );
float FMA_FastPow( float a, float b );
float FMA_FastApproxPow( float a, float b );
float FMA_FastLog( float x );
float FMA_FastExp( float x );
void FMA_FastSinCos( float x, float *s, float *c );
float FMA_FastAtan2( float x, float y );
int FMA_FastRoundLog2( double x ); 
int FMA_FastIntLog2( double x ); 

#ifdef __cplusplus
}
#endif

#endif /* FASTMATH_H */

