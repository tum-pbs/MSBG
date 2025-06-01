/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/
#ifndef VECTORCLASS_UTIL_H
#define VECTORCLASS_UTIL_H  

#include "vectorclass_util2.h"

/////////////////////////////////////////////////////
//  Misc. utility
/////////////////////////////////////////////////////

#if 0
static inline int horizontal_add(const Vec4i& x) {
    // Calculates the sum of SSE Register - https://stackoverflow.com/a/35270026/195787
    __m128i shufReg, sumsReg;

    shufReg = _mm_castps_si128(_mm_movehdup_ps(_mm_castsi128_ps(x)));        // Broadcast elements 3,1 to 2,0
    sumsReg = _mm_add_epi32(x, shufReg);    
    shufReg = _mm_castps_si128( _mm_movehl_ps ( _mm_castsi128_ps(shufReg), _mm_castsi128_ps(sumsReg))); // High Half -> Low Half
    sumsReg = _mm_add_epi32(sumsReg, shufReg);
    return  _mm_cvtsi128_si32(sumsReg); // Result in the lower part of the SSE Register
}
#endif

#if 0
static inline int horizontal_logical_or(const Vec4i& x) {
    __m128i shufReg, sumsReg;

    // Calculates the sum of SSE Register - https://stackoverflow.com/a/35270026/195787
    shufReg = _mm_castps_si128(_mm_movehdup_ps(_mm_castsi128_ps(x)));        // Broadcast elements 3,1 to 2,0
    sumsReg = _mm_or_si128(x, shufReg);    
    shufReg = _mm_castps_si128( _mm_movehl_ps ( _mm_castsi128_ps(shufReg), _mm_castsi128_ps(sumsReg))); // High Half -> Low Half
    sumsReg = _mm_or_si128(sumsReg, shufReg);
    return  _mm_cvtsi128_si32(sumsReg); // Result in the lower part of the SSE Register
}
#endif

//--------------------------------------------

#if 0
inline char *vprint( char *chp, const Vec2fb v ) 
{ 
  sprintf( chp, "%d,%d", v[0],v[1] );
  return chp;
}
inline char *vprint( char *chp, const Vec2f v ) 
{ 
  sprintf( chp, "%g,%g", v[0],v[1] );
  return chp;
}
#endif

inline char *vprint( char *chp, const Vec4fb v ) 
{ 
  sprintf( chp, "%d,%d,%d,%d", v[0],v[1],v[2],v[3] );
  return chp;
}
inline char *vprint( char *chp, const Vec4f v ) 
{ 
  sprintf( chp, "%g,%g,%g,%g", v[0],v[1],v[2],v[3] );
  return chp;
}
inline char *vprint( char *chp, const Vec4d v ) 
{ 
  sprintf( chp, "%g,%g,%g,%g", v[0],v[1],v[2],v[3] );
  return chp;
}

inline char *vprint( char *chp, const Vec8fb v ) 
{ 
  sprintf( chp, "%d,%d,%d,%d,%d,%d,%d,%d", v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7] );
  return chp;
}
inline char *vprint( char *chp, const Vec8f v ) 
{ 
  sprintf( chp, "%g,%g,%g,%g,%g,%g,%g,%g", v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7] );
  return chp;
}

inline char *vprint( char *chp, const Vec8ui v ) 
{ 
  sprintf( chp, "%d,%d,%d,%d,%d,%d,%d,%d", v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7] );
  return chp;
}

inline Vec4i vzero() { return _mm_setzero_si128(); }

inline Vec4f vfzero() { return _mm_setzero_ps(); }


inline Vec4f vlinstep(Vec4f a, Vec4f b, Vec4f x) 
{
  x = (x - a) * 1.0f/(b - a); 
  return max(0.0f,min(x,1.0f));
}

inline Vec8f vlinstep(Vec8f a, Vec8f b, Vec8f x) 
{
  x = (x - a) * 1.0f/(b - a); 
  return max(0.0f,min(x,1.0f));
}

inline Vec4f vsmoothstep(Vec4f a, Vec4f b, Vec4f x) 
{
  Vec4f y = (x - a) * 1.0f/(b - a); 
  y = (y*y * (3.0f - 2.0f*y));
  y = select( x<a, 
      	0.0f,
        select( x>=b, 1.0f,
		      y ));
  return y;
}
    
inline Vec4i square( const Vec4i& vec ) { return  vec*vec; }


#if 0
inline Vec4i vshift(Vec4i x, Vec4i s )
{
  return _mm_sllv_epi32( x, s );
}
#endif


enum VectorSelect
{
    Ax = 0, Ay = 1, Az = 2, Aw = 3,
    Bx = 8, By = 9, Bz = 10, Bw = 11,
};
template <VectorSelect S0, VectorSelect S1, VectorSelect S2, VectorSelect S3>
inline __m128 vshuffle_ps(__m128 x, __m128 y)
{
    return _mm_shuffle_ps(x, y, S0 + S1 * 4 + (S2 - Bx) * 16 + (S3 -Bx) * 64);
}

inline Vec4f vlerp_ps(const Vec4f& t, const Vec4f& a, const Vec4f& b)
{
  return mul_add( b-a, t, a );
}

template <VectorSelect S0, VectorSelect S1, VectorSelect S2, VectorSelect S3>
inline Vec8f vshuffle_ps(Vec8f x, Vec8f y)
{
    #if INSTRSET >= 7		// AVX
    return _mm256_shuffle_ps(x, y, S0 + S1 * 4 + (S2 - Bx) * 16 + (S3 -Bx) * 64);
    #else
    return Vec8f( vshuffle_ps<S0,S1,S2,S3>(x.get_low(),y.get_low()),
		  vshuffle_ps<S0,S1,S2,S3>(x.get_high(),y.get_high()));
    #endif
}

inline Vec8f vlerp_ps(Vec8f t, Vec8f a, Vec8f b)
{
  return mul_add( b-a, t, a );
}

/////////////////////////////////////////////////////
//  fractional part
/////////////////////////////////////////////////////
inline Vec4f vfrac( const Vec4f& x ) 
{
  return x - floor(x);
}


/////////////////////////////////////////////////////
//  VPOSPART / VNEGPART
/////////////////////////////////////////////////////

inline Vec8f vpospart( const Vec8f& x ) 
{ 
  return max(x,
        #if INSTRSET >= 7		// AVX
      	_mm256_setzero_ps()
	#else
	0.f
	#endif
	); 
} 

inline Vec8f vnegpart( const Vec8f& x ) 
{ 
  return min(x,
        #if INSTRSET >= 7		// AVX
      	_mm256_setzero_ps()
	#else
	0.f
	#endif
      );
} 

/////////////////////////////////////////////////////
//  branchless min/max/clamp-to-zero
/////////////////////////////////////////////////////
#if 0
inline float vclampz( float a )
{
  return _mm_cvtss_f32( max( Vec4f(a), 
			_mm_setzero_ps()) );
}
#else
inline float vclampz(float a )
{
  float result = a;
  asm ("maxss %1, %0" : "+x" (result) : "x" (0.0f));
  return result;
}
#endif

inline __m128 vclampz(const Vec4f &a )
{
  return max(a,_mm_setzero_ps());
}

inline __m128 vflipsign(const Vec4f &a )
{
  // https://stackoverflow.com/questions/3361132/flipping-sign-on-packed-sse-floats
  return _mm_xor_ps(a, _mm_set1_ps(-0.f));
}

inline Vec4f vapprox_sqrt( const Vec4f &a ) 
{
  return _mm_rcp_ps ( _mm_rsqrt_ps( a ) );
}
inline Vec8f vapprox_sqrt( const Vec8f &a ) 
{
#if INSTRSET >= 8  // AVX2
  return _mm256_rcp_ps ( _mm256_rsqrt_ps( a ) );
#else
  return sqrt(a);
#endif
}

inline float vmin(float a, float b)
{
  float result = a;
  asm ("minss %1, %0" : "+x" (result) : "x" (b));
  return result;
}

inline float vmax(float a, float b)
{
  float result = a;
  asm ("maxss %1, %0" : "+x" (result) : "x" (b));
  return result;
}

/////////////////////////////////////////////////////
//  V I S F I N I T E
/////////////////////////////////////////////////////
// https://stackoverflow.com/questions/30674291/how-to-check-inf-for-avx-intrinsic-m256
inline bool visfinite( Vec4f a )
{
  __m128 self_sub = _mm_sub_ps( a, a );
  return  !_mm_movemask_epi8(_mm_castps_si128(self_sub));
}

inline bool visfinite( Vec8f a )
{
  #if INSTRSET >= 7		// AVX
  __m256 self_sub = _mm256_sub_ps( a, a );
  return  !_mm256_movemask_epi8(_mm256_castps_si256(self_sub));
  #else
  return visfinite( a.get_low() ) &&  
         visfinite( a.get_high() );
  #endif
}

inline bool visfinite( Vec2d a )
{
  __m128d self_sub = _mm_sub_pd( a, a );
  return  !_mm_movemask_epi8(_mm_castpd_si128(self_sub));
}

inline bool visfinite( Vec4d a )
{
  #if INSTRSET >= 7		// AVX
  __m256d self_sub = _mm256_sub_pd( a, a );
  return  !_mm256_movemask_epi8(_mm256_castpd_si256(self_sub));
  #else
  return visfinite( a.get_low() ) &&  
         visfinite( a.get_high() );
  #endif
}

/////////////////////////////////////////////////////
//  V A L L / VANY
/////////////////////////////////////////////////////

inline int vany3( const __m128& mask )
{
  return _mm_movemask_ps(mask) & 7; 
}
inline int vany3( const __m128i& mask )
{
  return _mm_movemask_ps(_mm_castsi128_ps(mask)) & 7; 
}

inline int vany( const __m128& mask )
{
  return _mm_movemask_ps(mask); 
}

inline int vany( const __m128i& mask )
{
  return _mm_movemask_ps(_mm_castsi128_ps(mask)); 
}

inline __attribute__((always_inline)) int vall( const __m128& mask )
{
  return _mm_movemask_ps(mask)==0xF; 
}

inline __attribute__((always_inline)) int vall( const __m128i& mask )
{
  return _mm_movemask_ps(_mm_castsi128_ps(mask))==0xF; 
}

inline __attribute__((always_inline)) int vall3( const __m128i& mask )
{
  return (_mm_movemask_ps(_mm_castsi128_ps(mask))&7) == 7;
}

#if INSTRSET >= 7		// AVX
inline __attribute__((always_inline)) int vall( const __m256d& mask )
{
  return _mm256_movemask_pd(mask)==0xF; 
}
#else
inline int vall( const Vec4db& m )
{
  return _mm_movemask_pd(m.get_low())==0x3 &&
	 _mm_movemask_pd(m.get_high())==0x3;
}
#endif

#if INSTRSET >= 7		// AVX
inline int vany( const __m256d& mask )
{
  return _mm256_movemask_pd(mask); 
}

inline int vany( const __m256& mask )
{
  return _mm256_movemask_ps(mask); 
}

inline int vall( const __m256& mask )
{
  return _mm256_movemask_ps(mask)==0xFF; 
}

#else
inline int vany( const Vec4db& m )
{
  return _mm_movemask_pd(m.get_low()) ||
	 _mm_movemask_pd(m.get_high());
}

inline int vany( const Vec8fb& m )
{
  return _mm_movemask_ps(m.get_low()) ||
	 _mm_movemask_ps(m.get_high());
}

inline int vany( const Vec8i &m )
{
  return vany( m.get_low() ) || vany( m.get_high() ); 
}

/*inline int vall( const Vec8i &m )
{
  return vall( m.get_low() ) && vall( m.get_high() ); 
}*/

inline int vall( const Vec8f &m )
{
  return vall( m.get_low() ) && vall( m.get_high() ); 
}

#endif


// AVX




/////////////////////////////////////////////////////
//  V G E T
/////////////////////////////////////////////////////


#if 0
inline int32_t vget_x(const __m128i& vec){return _mm_extract_epi32(vec,0);}
#else
inline int32_t vget_x(const __m128i& vec){return _mm_cvtsi128_si32(vec);}
#endif

inline int32_t vget_y(const __m128i& vec){return _mm_extract_epi32(vec,1);}
inline int32_t vget_z(const __m128i& vec){return _mm_extract_epi32(vec,2);}
inline int32_t vget_w(const __m128i& vec){return _mm_extract_epi32(vec,3);}

inline Vec4i vset_x_(const __m128i& vec, int x) { return _mm_insert_epi32(vec,x,0); } 
#define vset_x(vec,x) (vec = vset_x_(vec,x) )
inline Vec4i vset_w_(const __m128i& vec, int x) { return _mm_insert_epi32(vec,x,3); } 
#define vset_w(vec,x) (vec = vset_w_(vec,x) )
inline Vec4ui vsetui_w_(const __m128i& vec, int x) { return _mm_insert_epi32(vec,x,3); } 


inline float vfget_x(const Vec4f& vec){return  _mm_cvtss_f32( vec );}
inline float vfget_y(const Vec4f& vec) {
        Vec4f t = permute4f<1,-256,-256,-256>(vec);
        return vfget_x(t); }
inline float vfget_z(const Vec4f& vec) {
        Vec4f t = permute4f<2,-256,-256,-256>(vec);
        return vfget_x(t); }
inline float vfget_w(const Vec4f& vec) {
        Vec4f t = permute4f<3,-256,-256,-256>(vec);
        return vfget_x(t); }

inline double vdget_x(const Vec4d& vec){return  _mm_cvtsd_f64(vec.get_low());}
inline double vdget_y(const Vec4d& vec) {
        Vec4d t = permute4d<1,-256,-256,-256>(vec);
        return vdget_x(t); }
inline double vdget_z(const Vec4d& vec) {
        Vec4d t = permute4d<2,-256,-256,-256>(vec);
        return vdget_x(t); }
inline double vdget_w(const Vec4d& vec) {
        Vec4d t = permute4d<3,-256,-256,-256>(vec);
        return vdget_x(t); }


inline double vfget_x(const Vec4d& vec) { return vdget_x(vec); }
inline double vfget_y(const Vec4d& vec) { return vdget_y(vec); }
inline double vfget_z(const Vec4d& vec) { return vdget_z(vec); }
inline double vfget_w(const Vec4d& vec) { return vdget_w(vec); }


template <int pos>
void vinsert(Vec8f& v, float x) 
{
  #if INSTRSET >= 7		// AVX
  Vec8f xv =  _mm256_broadcast_ss(&x);
  v = _mm256_blend_ps(v, xv, 1<<pos);
  #else
  v.insert(pos,x);
  #endif
}

inline Vec8f vuint2float( Vec8ui x )
{
  return Vec8f( _mm_cvtepi32_ps(x.get_low()), _mm_cvtepi32_ps(x.get_high()));
}

inline Vec4f vuint2float( Vec4ui x )
{
  return _mm_cvtepi32_ps(x);
}

inline Vec4f vint2float( Vec4i x )
{
  return _mm_cvtepi32_ps(x);
}



/////////////////////////////////////////////////////
//  convert from/to 8/16 bit memory
/////////////////////////////////////////////////////

inline void vload_from_uint16_a( Vec8ui &z, const uint16_t *p )
{
  Vec8us y = Vec8us().load_a( p );
  Vec4ui  y4l = _mm_unpacklo_epi16(y, _mm_set1_epi16(0)),
	  y4h = _mm_unpackhi_epi16(y, _mm_set1_epi16(0));
  z = Vec8ui(y4l,y4h);
}

inline void vload_from_uint16_a( Vec4ui &z, const uint16_t *p )
{
  Vec8us y(p[0],p[1],p[2],p[3],0,0,0,0);
  z = _mm_unpacklo_epi16(y, _mm_set1_epi16(0));
}

inline void vload_from_uint16_a( Vec8f &w, const uint16_t *p )
{
#if 1
  Vec8us y = Vec8us().load_a( p );
  Vec4ui  y4l = _mm_unpacklo_epi16(y, _mm_set1_epi16(0)),
	  y4h = _mm_unpackhi_epi16(y, _mm_set1_epi16(0));

  #if 1
  w = Vec8f( _mm_cvtepi32_ps(y4l),
	      _mm_cvtepi32_ps(y4h) );

  #else
  #if INSTRSET >= 7		// AVX
  Vec8ui z(y4l,y4h);
  w = _mm256_cvtepi32_ps(z);
  #else
  w = Vec8f( _mm_cvtepi32_ps(y4l),
	        _mm_cvtepi32_ps(y4h) );
  #endif  
  #endif
#else
  w = Vec8f( p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]);
#endif
}

inline void vload_from_uint8_a(Vec8ui &w, uint8_t *source)
{
#if INSTRSET >= 8  // AVX2
    #if 0
    // Load 8 bytes into a 128-bit XMM register using _mm_loadl_epi64
    __m128 temp = _mm_loadl_pi(_mm_setzero_ps(),(const __m64*)source);
    // Zero-extend to 256-bits in a YMM register
    __m256i temp256 = _mm256_cvtepu8_epi32(_mm_castps_si128(temp));
    // Convert unsigned 32-bit integers to single-precision floats
    w = temp256;
     #else
    // Load 8 bytes into a 128-bit XMM register using _mm_loadl_epi64
    __m128i temp = _mm_loadl_epi64((const __m128i*)source);
    // Zero-extend to 256-bits in a YMM register
    __m256i temp256 = _mm256_cvtepu8_epi32(temp);
    // Convert unsigned 32-bit integers to single-precision floats
    w = temp256;
    #endif
#else
    w = Vec8ui(source[0],source[1],source[2],source[3],
	      source[4],source[5],source[6],source[7] );
#endif
}

inline void vload_from_uint8_a(Vec8f &w, uint8_t *source)
{
#if INSTRSET >= 8  // AVX2
    // Load 8 bytes into a 128-bit XMM register using _mm_loadl_epi64
    __m128i temp = _mm_loadl_epi64((const __m128i*)source);
    // Zero-extend to 256-bits in a YMM register
    __m256i temp256 = _mm256_cvtepu8_epi32(temp);
    // Convert unsigned 32-bit integers to single-precision floats
    w = _mm256_cvtepi32_ps(temp256);
#else
    w = Vec8f(source[0],source[1],source[2],source[3],
	      source[4],source[5],source[6],source[7] );
#endif
}

inline void vload_from_uint16_a( Vec4f &z, const uint16_t *p )
{
  Vec8us y(p[0],p[1],p[2],p[3],0,0,0,0);
  Vec4ui w = _mm_unpacklo_epi16(y, _mm_set1_epi16(0));
  z = _mm_cvtepi32_ps(w);
}
inline void vload_from_uint8_a( Vec4f &z, const uint8_t *p )
{
  z = Vec4f( p[0], p[1], p[2], p[3] );
}


inline void vload8(Vec8ui &w, const uint8_t *source)
{
  vload_from_uint8_a( w, (uint8_t*)(source) );
}

inline void vload8(Vec8ui &w, const uint32_t *source)
{
  w.load(source);
}



inline void vstore_to_uint16_a( const Vec8f &x, uint16_t *p )
{
  Vec8i xi = truncate_to_int(round(x));
  Vec8us xus = compress( (Vec4ui)xi.get_low(), (Vec4ui)xi.get_high() );
  xus.store_a((void*)p);
}

inline void vstore_to_uint16_a( const Vec4f &x, uint16_t *p )
{
  Vec4i xi = truncate_to_int(round(x));
  p[0] = vget_x(xi);
  p[1] = vget_y(xi);
  p[2] = vget_z(xi);
  p[3] = vget_w(xi);
}

inline void vstore_to_uint16_a( const Vec8ui &xi, uint16_t *p )
{
  Vec8us xus = compress( (Vec4ui)xi.get_low(), (Vec4ui)xi.get_high() );
  xus.store_a((void*)p);
}

inline void vstore_to_uint16_a( const Vec4ui &xi, uint16_t *p )
{
  p[0] = vget_x(xi);
  p[1] = vget_y(xi);
  p[2] = vget_z(xi);
  p[3] = vget_w(xi);
}

#if 0
inline void vstore_to_uint16_nt( const Vec8f &x, uint16_t *p )
{
  Vec8i xi = truncate_to_int(round(x));
  Vec8us xus = compress( (Vec4ui)xi.get_low(), (Vec4ui)xi.get_high() );
  xus.store_a((void*)p);
}
#endif

inline void vstore_to_uint16_nt( const Vec4f &x, uint16_t *p )
{
  vstore_to_uint16_a(x,p);
}
inline void vstore_to_uint16_nt( const Vec8ui &xi, uint16_t *p )
{
  Vec8us xus = compress( (Vec4ui)xi.get_low(), (Vec4ui)xi.get_high() );
  _mm_stream_si128((__m128i*) p, xus);
}
inline void vstore_to_uint16_nt( const Vec4ui &xi, uint16_t *p )
{
  vstore_to_uint16_a(xi,p);
}

inline void vstore_to_uint8_a( const Vec8ui &xi, uint8_t *p )
{
#if INSTRSET >= 8  // AVX2

#if 1
  Vec8ui shuffle_control = _mm256_set_epi8( 
				    -1,-1,-1,-1, -1,-1,-1,-1,
				    -1,-1,-1,-1, 12,8,4,0,

				    -1,-1,-1,-1, -1,-1,-1,-1,
				    -1,-1,-1,-1, 12,8,4,0 );    

  Vec8ui tmp = _mm256_shuffle_epi8(xi, shuffle_control);
  Vec4ui tmp_low = tmp.get_low();
  Vec4ui tmp_high = tmp.get_high();
  Vec4ui packed = blend4i<0,4,-1,-1>( tmp_low, tmp_high );       
  _mm_storel_pi( (__m64*)p, _mm_castsi128_ps(packed));
#else

//  { uint8_t buf[32]; xi.store((int32_t*)buf); printf("xi= "); for(int i=0;i<32;i++) printf("%d ",buf[i]); printf("\n"); }

  Vec8ui tmp = _mm256_shuffle_epi8(xi, _mm256_set_epi8( 
							-1,-1,-1,-1,
							-1,-1,-1,-1,
							-1,-1,-1,-1,
							12,8,4,0,

							-1,-1,-1,-1,
							-1,-1,-1,-1,
							-1,-1,-1,-1,
							12,8,4,0

							));
  //{ uint8_t buf[32]; tmp.store((int32_t*)buf); printf("tmp= "); for(int i=0;i<32;i++) printf("%d ",buf[i]); printf("\n"); } 

  _mm_maskmove_si64( _mm_movepi64_pi64(tmp.get_low()), 
      		     _mm_set_pi16(0x0000, 0x0000, 0xFFFF, 0xFFFF), 
		    (char*)p );

  _mm_maskmove_si64( _mm_movepi64_pi64(tmp.get_high()), 
      		     _mm_set_pi16(0x0000, 0x0000, 0xFFFF, 0xFFFF), 
		    ((char*)p)+4 );
#endif
#else
  int x[8] __attribute__((aligned(64)));
  xi.store_a(x);
  uint8_t *s=(uint8_t*)x;
  p[0] = s[0];
  p[1] = s[4];
  p[2] = s[8];
  p[3] = s[12];
  p[4] = s[16];
  p[5] = s[20];
  p[6] = s[24];
  p[7] = s[28];
#endif
}

inline void vstore_to_uint8_( const Vec8ui &xi, uint8_t *p )
{
#if INSTRSET >= 8  // AVX2
  Vec8ui shuffle_control = _mm256_set_epi8( 
				    -1,-1,-1,-1, -1,-1,-1,-1,
				    -1,-1,-1,-1, 12,8,4,0,

				    -1,-1,-1,-1, -1,-1,-1,-1,
				    -1,-1,-1,-1, 12,8,4,0 );    

  Vec8ui tmp = _mm256_shuffle_epi8(xi, shuffle_control);
  Vec4ui tmp_low = tmp.get_low();
  Vec4ui tmp_high = tmp.get_high();
  Vec4ui packed = blend4i<0,4,-1,-1>( tmp_low, tmp_high );       
  _mm_storel_pi( (__m64*)p, _mm_castsi128_ps(packed));
#else
  int x[8] __attribute__((aligned(64)));
  xi.store_a(x);
  uint8_t *s=(uint8_t*)x;
  p[0] = s[0];
  p[1] = s[4];
  p[2] = s[8];
  p[3] = s[12];
  p[4] = s[16];
  p[5] = s[20];
  p[6] = s[24];
  p[7] = s[28];
#endif
}

inline void vstore8(Vec8ui &w, uint32_t *dest)
{
  w.store(dest);
}

/////////////////////////////////////////////////////
//  horizontal sum,min,max
/////////////////////////////////////////////////////

static inline double horizontal_max( const Vec4d& a )
{
  return std::max(std::max(vdget_x(a),vdget_y(a)),
	          std::max(vdget_z(a),vdget_w(a)));
}

static inline float horizontal_max(__m128 x) {
    // Calculates the sum of SSE Register - http://stackoverflow.com/a/35270026/195787
    __m128 shufReg, sumsReg;
    shufReg = _mm_movehdup_ps(x);        // Broadcast elements 3,1 to 2,0
    sumsReg = _mm_max_ps(x, shufReg);
    shufReg = _mm_movehl_ps(shufReg, sumsReg); // High Half -> Low Half
    sumsReg = _mm_max_ss(sumsReg, shufReg);
    return  _mm_cvtss_f32(sumsReg); // Result in the lower part of the SSE Register
}

inline float horizontal_max( const Vec8f& x )
{
  return _mm_cvtss_f32( 
      		_mm_max_ss( Vec4f(horizontal_max(x.get_low())),
      	             	    Vec4f(horizontal_max(x.get_high()))));
}

static inline float horizontal_min(__m128 x) {
    // Calculates the sum of SSE Register - http://stackoverflow.com/a/35270026/195787
    __m128 shufReg, sumsReg;
    shufReg = _mm_movehdup_ps(x);        // Broadcast elements 3,1 to 2,0
    sumsReg = _mm_min_ps(x, shufReg);
    shufReg = _mm_movehl_ps(shufReg, sumsReg); // High Half -> Low Half
    sumsReg = _mm_min_ss(sumsReg, shufReg);
    return  _mm_cvtss_f32(sumsReg); // Result in the lower part of the SSE Register
}


// Faster Horizontal add: Calculates the sum of all vector elements.
// https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-float-vector-sum-on-x86
static inline float horizontal_sum (Vec4f const & a) {
    __m128 t1 = _mm_movehl_ps(a,a);
    __m128 t2 = _mm_add_ps(a,t1);
    __m128 t3 = _mm_shuffle_ps(t2,t2,1);
    __m128 t4 = _mm_add_ss(t2,t3);
    return _mm_cvtss_f32(t4);
}

#if 0  // TODO: test
inline float horizontal_sum(const Vec4i &mABCD)
{ // http://microperf.blogspot.de/2016/12/the-horizontal-sum.html
    __m128i mCDCD = _mm_castps_si128(_mm_movehl_ps(_mm_castsi128_ps(mABCD), 
	  					   _mm_castsi128_ps(mABCD)));
    __m128i mApCBpD = _mm_add_epi32(mABCD, mCDCD);
    __m128i mBpD = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(mApCBpD), 
	  				           _mm_castsi128_ps(mApCBpD), 0x55));
    __m128i mApBpCpD = _mm_add_epi32(mApCBpD, mBpD);
    return vget_x(mApBpCpD);
}
#endif

static inline int vhadd3( Vec4i const & a )
{
  return vget_x(a)+vget_y(a)+vget_z(a); 
}

#if 0
/////////////////////////////////////////////////////
//  S C A T T E R
/////////////////////////////////////////////////////

// Store 4 integers from SSE vector using offsets from another vector
inline void vscatter( float* rdi, __m128i idx, __m128 data )
{
    rdi[ (uint32_t)_mm_cvtsi128_si32( idx ) ] = _mm_cvtss_f32( data );
    rdi[ (uint32_t)_mm_extract_epi32( idx, 1 ) ] = _mm_extract_ps( data, 1 );
    rdi[ (uint32_t)_mm_extract_epi32( idx, 2 ) ] = _mm_extract_ps( data, 2 );
    rdi[ (uint32_t)_mm_extract_epi32( idx, 3 ) ] = _mm_extract_ps( data, 3 );
}

// Store 8 integers from AVX vector using offsets from another vector
inline void vscatter( float* rdi, __m256i idx, __m256 data )
{
    vscatter( rdi, _mm256_castsi256_si128( idx ), _mm256_castps256_ps128( data ) );
    vscatter( rdi, _mm256_extracti128_si256( idx, 1 ), _mm256_extractf128_ps( data, 1 ) );
}
#endif

/////////////////////////////////////////////////////
//  G A T H E R
/////////////////////////////////////////////////////

inline __m128 vgather4(const float *base, const __m128i &indices)
{
#if INSTRSET >= 8  // AVX2
  return _mm_i32gather_ps( base, indices, 4 );
#else
__m128 res = _mm_load_ss(base+_mm_cvtsi128_si32(indices));

res = _mm_insert_ps(res,_mm_load_ss(base+_mm_extract_epi32(indices,1)),_MM_MK_INSERTPS_NDX(0,1,0));

res = _mm_insert_ps(res,_mm_load_ss(base+_mm_extract_epi32(indices,2)),_MM_MK_INSERTPS_NDX(0,2,0));

res = _mm_insert_ps(res,_mm_load_ss(base+_mm_extract_epi32(indices,3)),_MM_MK_INSERTPS_NDX(0,3,0));

return res;
#endif
}

inline __m128 vgather4(const float **bases, const int64_t idx, const int n)
{
__m128 res;
if(n>0)
  res = _mm_load_ss(bases[0]+idx);
else
  res = _mm_setzero_ps();
if(n>1) res = _mm_insert_ps(res,_mm_load_ss(bases[1]+idx),_MM_MK_INSERTPS_NDX(0,1,0));
if(n>2) res = _mm_insert_ps(res,_mm_load_ss(bases[2]+idx),_MM_MK_INSERTPS_NDX(0,2,0));
if(n>3) res = _mm_insert_ps(res,_mm_load_ss(bases[3]+idx),_MM_MK_INSERTPS_NDX(0,3,0));
return res;
}

// masked gather
inline __m128 vgather4(const float *base, const __m128i& indices, const __m128& mask)

{
  const int m = _mm_movemask_ps(mask); 
  __m128 res;

  if(m&1) 
    res = _mm_load_ss(base+_mm_cvtsi128_si32(indices));
  else 
    res = _mm_setzero_ps();

  if(m&2) res = _mm_insert_ps(res,_mm_load_ss(base+_mm_extract_epi32(indices,1)),_MM_MK_INSERTPS_NDX(0,1,0));

  if(m&4) res = _mm_insert_ps(res,_mm_load_ss(base+_mm_extract_epi32(indices,2)),_MM_MK_INSERTPS_NDX(0,2,0));

  if(m&8) res = _mm_insert_ps(res,_mm_load_ss(base+_mm_extract_epi32(indices,3)),_MM_MK_INSERTPS_NDX(0,3,0));

  return res;

}

inline __m128 vgather4(const float **bases, const Vec4i& indices, const int n)
{
__m128 res;
if(n>0)
  res = _mm_load_ss(bases[0]+vget_x(indices));
else
  res = _mm_setzero_ps();
if(n>1) res = _mm_insert_ps(res,_mm_load_ss(bases[1]+vget_y(indices)),_MM_MK_INSERTPS_NDX(0,1,0));
if(n>2) res = _mm_insert_ps(res,_mm_load_ss(bases[2]+vget_z(indices)),_MM_MK_INSERTPS_NDX(0,2,0));
if(n>3) res = _mm_insert_ps(res,_mm_load_ss(bases[3]+vget_w(indices)),_MM_MK_INSERTPS_NDX(0,3,0));
return res;
}


#if INSTRSET >= 7		// AVX

inline __m256 vgather8(const float *base, const Vec8i &indices)

{
#if INSTRSET >= 8 && defined(VECTORI256_H) // AVX2
    return _mm256_i32gather_ps(base, indices, 4);
#else
const __m128 low = vgather4(base,indices.get_low()),
	     high = vgather4(base,indices.get_high());

return _mm256_insertf128_ps(_mm256_castps128_ps256(low),high,1);
#endif
}


inline __m256 vgather8(const float *base, const Vec8i &indices, const Vec8fb& mask)

{

const __m128 low = vgather4(base,indices.get_low(),
    				 _mm256_castps256_ps128(mask)),

             high = vgather4(base,indices.get_high(),
		 		  _mm256_extractf128_ps(mask,1));

return _mm256_insertf128_ps(_mm256_castps128_ps256(low),high,1);

}

//
inline __m256 vgather8(const float **bases, const int64_t idx, const int n)

{
const __m128 low = vgather4(bases,idx,n),
             high = n>4 ? vgather4(bases+4,idx,n-4) : _mm_setzero_ps();
return _mm256_insertf128_ps(_mm256_castps128_ps256(low),high,1);
}

//
inline __m256 vgather8(const float **bases, const Vec8i& indices, const int n)

{
const __m128 low = vgather4(bases,indices.get_low(),n),
             high = n>4 ? vgather4(bases+4,indices.get_high(),n-4) : _mm_setzero_ps();
return _mm256_insertf128_ps(_mm256_castps128_ps256(low),high,1);
}

#else				// SSE

inline Vec8f vgather8(const float *base, const Vec8i &indices)

{

const __m128 low = vgather4(base,indices.get_low()),
	     high = vgather4(base,indices.get_high());

Vec8f res(low,high);

return res;

}

inline Vec8f vgather8(const float **bases, const int64_t idx, const int n)

{
const __m128 low = vgather4(bases,idx,n),
             high = vgather4(bases+4,idx,n>=4?n-4:0);
Vec8f res(low,high);
return res;
}

inline Vec8f vgather8(const float **bases, const Vec8i& indices, const int n)

{
const __m128 low = vgather4(bases,indices.get_low(),n),
             high = vgather4(bases+4,indices.get_high(),n>=4?n-4:0);
Vec8f res(low,high);
return res;
}

#endif  

/////////////////////////////////////////////////////
//  S E T _ E L E M
/////////////////////////////////////////////////////
template<int pos>
static inline Vec4f v4f_set_elem( const Vec4f& v, float f ) 
{
  return _mm_insert_ps( v, _mm_set_ss( f ),  _MM_MK_INSERTPS_NDX(0,pos,0)   );
}

/////////////////////////////////////////////////////
//  L O A D _ A T
/////////////////////////////////////////////////////
template<int pos>
static inline Vec4f v4f_load_to( float const * p, const Vec4f& v ) 
{
  return _mm_insert_ps( v, _mm_load_ss( p ),  _MM_MK_INSERTPS_NDX(0,pos,0)   );
}


template<int pos>
static inline Vec4f v4f_zero_to( const Vec4f& v ) 
{
  return _mm_insert_ps( v, _mm_setzero_ps(),  _MM_MK_INSERTPS_NDX(0,pos,0)   );
}


// v4_zero_to
template<int pos>
static inline Vec4f v4_zero_to( const Vec4f& v ) 
{
  return _mm_insert_ps( v, _mm_setzero_ps(),  _MM_MK_INSERTPS_NDX(0,pos,0)   );
}

template<int pos>
static inline Vec4d v4_zero_to( const Vec4d& v ) 
{
  Vec4d w = v;
  return w.insert( pos, 0. );
}

static inline void vzero_w( Vec4f& v ) 
{
  v = v4f_zero_to<3>(v);
}



// v4_load_to
template<int pos>
static inline Vec4f v4_load_to( float const * p, const Vec4f& v ) 
{
  return _mm_insert_ps( v, _mm_load_ss( p ),  _MM_MK_INSERTPS_NDX(0,pos,0)   );
}

template<int pos>
static inline Vec4d v4_load_to( double const * p, const Vec4d& v ) 
{
  Vec4d w = v;
  return w.insert( pos, *p );
}

template<int pos>
static inline Vec4d v4_load_to( float const * p, const Vec4d& v ) 
{
  Vec4d w = v;
  return w.insert( pos, *p );
}


inline float v3normSq( const Vec4f& vec ) 
{ 
  Vec4f v = v4f_zero_to<3>(vec);
  v *= v;
  return vfget_x(v)+vfget_y(v)+vfget_z(v);
}


template<int pos>
static inline Vec8f vload_to( Vec8f& x, float const *p )
{
#if INSTRSET>=7  // AVX
  x = _mm256_blend_ps (x,  _mm256_broadcast_ss(p), (1<<pos));
#else
  if (pos < 4) 
    return Vec8f( v4f_load_to<pos>( p, x.get_low() ), x.get_high() );
  else 
    return Vec8f( x.get_low(), v4f_load_to<pos-4>( p, x.get_high() ) );
#endif
  return x;
}


#if INSTRSET >= 7		// AVX

#if 0

inline void vload_to(float * p, Vec8f& ymm, uint32_t mask /* (mask=1<<index) */) 
{
  ymm = _mm256_blend_ps (ymm,  _mm256_broadcast_ss(p), mask);
}
#else




#define vload_to(p, ymm, mask) \
{ \
  ymm = _mm256_blend_ps (ymm,  _mm256_broadcast_ss(p), mask);\
}

#endif

#else				// SSE

inline void vload_to(float * p, Vec8f& ymm, uint32_t mask /* (mask=1<<index) */) 
{
  int i=0;
  while(mask>>=1)  i++;
  ymm.insert(i,*p);
}
#endif  

/////////////////////////////////////////////////////
//  S T R E A M I N G    L O A D
/////////////////////////////////////////////////////

#if INSTRSET >= 8  // AVX2

inline Vec8f vstream_load(const float * p) 
{
  return _mm256_castsi256_ps (_mm256_stream_load_si256((const __m256i *)p));
}

#else

inline Vec8f vstream_load(const float * p) 
{
  Vec8f v;
  v.load_a(p);
  return v;
}

#endif

/////////////////////////////////////////////////////
//  S T R E A M I N G    S T O R E
/////////////////////////////////////////////////////

#if 0 // performance ?
// stream single float
inline void vstream( float *p, float f )
{
  __m128 v = _mm_set_ss(f);
  _mm_maskmoveu_si128(reinterpret_cast<__m128i>(v),
      		      _mm_set_epi32(0,0,0,-1), 
		      (char*)(p));
}

inline void vstream( uint16_t *p, uint16_t val )
{
  __m128i v = _mm_set1_epi16(val);
  _mm_maskmoveu_si128(reinterpret_cast<__m128i>(v),
      		      _mm_set_epi32(0,0,0,0x0000FFFF), 
		      (char*)(p));
}
#endif

inline void vstream(float * p, Vec4f ymm) 
{
  _mm_stream_ps(p, ymm); 
}

#if INSTRSET >= 7		// AVX

inline void vstream(float * p, Vec8f ymm) 
{
  _mm256_stream_ps(p, ymm); 
}

inline void vstream(double * p, Vec4d ymm) 
{
  _mm256_stream_pd(p, ymm); 
}

#else				// SSE

inline void vstream(float * p, Vec8f ymm) 
{
  ymm.store(p); 
}

inline void vstream(double * p, Vec4d ymm) 
{
  ymm.store(p); 
}

#endif  

/////////////////////////////////////////////////////
//  non-AVX emulations
/////////////////////////////////////////////////////

#if INSTRSET < 7		// non-AVX


#undef _mm256_castsi256_ps
static inline Vec8f _mm256_castsi256_ps( const Vec8ui& x )
{
  return Vec8f( _mm_castsi128_ps(x.get_low()), _mm_castsi128_ps(x.get_high()) ); 
}

#undef _mm256_cvtepi32_epi64
static inline Vec4uq _mm256_cvtepi32_epi64( const Vec4i& x )
{
  int32_t x_[4];
  x.store(x_);
  return Vec4uq(x_[0],x_[1],x_[2],x_[3]);
}

#define _mm_maskstore_ps( addr, msk, a )\
\
_mm_maskmoveu_si128( _mm_castps_si128(a), msk, \
		       (char*)(addr));  
#undef _mm256_extractf128_ps
#define _mm256_extractf128_ps( v8, m ) ((m)?(v8).get_high():(v8).get_low())

#undef _mm256_cvtpd_ps
static inline __m128 _mm256_cvtpd_ps (Vec4d a)
{
  return blend4f<0,1,4,5>(_mm_cvtpd_ps(a.get_low()),
      			  _mm_cvtpd_ps(a.get_high()));
}

#undef _mm256_cvtps_pd
static inline Vec4d _mm256_cvtps_pd (Vec4f a)
{
  return Vec4d( extend_low(a), extend_high(a) );
}


#if 0
#undef _mm256_cvtpd_epi32
static inline __m128i _mm256_cvtpd_epi32 (Vec4d a)
{
  return blend4i<0,1,4,5>(_mm_cvtpd_epi32(a.get_low()),
      			  _mm_cvtpd_epi32(a.get_high()));
}
#else
static inline __m128i _mm256_cvtpd_epi32 (Vec4d a)
{
  printf("fatal error: %s:%d\n",__FILE__,__LINE__);
  exit(0);
}
#endif

#undef _mm256_insertf128_ps
static inline Vec8f _mm256_insertf128_ps( Vec8f a, Vec4f b, int m )
{
  return m ? Vec8f( a.get_low(), b ) : Vec8f( b, a.get_high());
}

#define _mm256_castps128_ps256(x) Vec8f( x, 0 )

#endif  // non-AVX

static inline Vec4d vcastps_pd( const Vec4f& x )
{
  return _mm256_cvtps_pd( x );
}

// vmask4f
static inline Vec4f vmask4( const Vec4f& x, const Vec4ui mask )
{
  return _mm_and_ps( x ,_mm_castsi128_ps(mask) );
}

static inline Vec4d vmask4( const Vec4d& x, const Vec4ui mask )
{
  #if INSTRSET >= 7		// AVX
  return x & (Vec4d)_mm256_cvtps_pd( _mm_castsi128_ps(mask) );
  #else
  Vec4f mf = _mm_castsi128_ps(mask);
  Vec4d md = _mm256_cvtps_pd(mf);
  return x &  md;  
  #endif
}

static inline void vmaskstore( float *dest, const Vec4f& x, const Vec4ui mask )
{
  _mm_maskstore_ps( dest, mask, x );
}

static inline void vmaskstore( float *dest, const Vec8f& x, const Vec8ui mask )
{
  #if INSTRSET >= 7		// AVX
  _mm256_maskstore_ps( dest, mask, x );
  #else
  float x_[8];
  x.store(x_);
  uint32_t m_[8];
  mask.store(m_);
  if(m_[0]) dest[0] = x_[0];
  if(m_[1]) dest[1] = x_[1];
  if(m_[2]) dest[2] = x_[2];
  if(m_[3]) dest[3] = x_[3];
  if(m_[4]) dest[4] = x_[4];
  if(m_[5]) dest[5] = x_[5];
  if(m_[6]) dest[6] = x_[6];
  if(m_[7]) dest[7] = x_[7];
  #endif
}


// vmaskstore4f
static inline void vmaskstore4( float *dest, const Vec4f& x, const Vec4ui mask )
{
  _mm_maskstore_ps( dest, mask, x );
}

// vmaskstore4f
static inline void vmaskstore4( double *dest, const Vec4d& x, const Vec4ui mask )
{
  #if INSTRSET >= 7		// AVX
  Vec4uq maskuq = _mm256_cvtepi32_epi64(mask);    
  _mm256_maskstore_pd( dest, maskuq, x );
  #else
  double x_[4];
  x.store(x_);
  uint32_t m_[4];
  mask.store(m_);
  if(m_[0]) dest[0] = x_[0];
  if(m_[1]) dest[1] = x_[1];
  if(m_[2]) dest[2] = x_[2];
  if(m_[3]) dest[3] = x_[3];
  #endif
}


/////////////////////////////////////////////////////
//  M I S C
/////////////////////////////////////////////////////

inline Vec4f vconv_pd_ps( const Vec4d& a )
{
  #if INSTRSET < 7		// non-AVX
  return blend4f<0,1,4,5>(_mm_cvtpd_ps(a.get_low()),
			    _mm_cvtpd_ps(a.get_high()));
  #else 
  return _mm256_cvtpd_ps( a );
  #endif
}

inline Vec4d vconv_ps_pd( const Vec4f& a )
{
  #if INSTRSET < 7		// non-AVX
  return Vec4d( extend_low(a), extend_high(a) );
  #else 
  return _mm256_cvtps_pd( a );
  #endif
}

/*static inline Vec8i truncate_to_uint(Vec8f const & a) {
    return _mm256_cvttps_epu32(a);
}*/

// Load Vec4d from single prec. float ptr

static inline Vec4d loadf( float const * p) 
{
#if INSTRSET < 7		// non-AVX
  Vec4f a = _mm_loadu_ps(p);
  return Vec4d( extend_low(a), extend_high(a) );  
#else
  return _mm256_cvtps_pd(_mm_loadu_ps(p));
#endif
}
static inline Vec4d loadf_a( float const * p) 
{
#if INSTRSET < 7		// non-AVX
  Vec4f a = _mm_load_ps(p);
  return Vec4d( extend_low(a), extend_high(a) );  
#else
  return _mm256_cvtps_pd(_mm_load_ps(p));
#endif
}

static inline Vec4d loadf_a( double const * p) 
{
  return Vec4d().load_a(p);
}


static inline void vstreamf_a( float * p, const Vec4d& v ) 
{
  _mm_stream_ps(p,_mm256_cvtpd_ps(v));
}
static inline void vstreamf_a( double * p, const Vec4d& v ) 
{
  v.store_a(p);
}

static inline void vstoref_a( float * p, const Vec4d& v ) 
{
  _mm_store_ps(p,_mm256_cvtpd_ps(v));
}
static inline void vstoref_a( double * p, const Vec4d& v ) 
{
  v.store_a(p);
}

static inline char *vecsprintf( char *buf, const Vec4d& v )
{
  sprintf(buf,"%g %g %g %g",v[0],v[1],v[2],v[3]);
  return buf;
}
static inline char *vecsprintf( char *buf, const Vec4f& v )
{
  sprintf(buf,"%g %g %g %g",v[0],v[1],v[2],v[3]);
  return buf;
}
static inline char *vecsprintf( char *buf, const Vec4ui& v )
{
  sprintf(buf,"%u %u %u %u",v[0],v[1],v[2],v[3]);
  return buf;
}

#endif  // VECTORCLASS_UTIL_H

