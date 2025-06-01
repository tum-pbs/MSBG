/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#ifndef UTIL_H
#define UTIL_H

#ifndef MIMP_ON_WINDOWS
#include <errno.h>
#include <pthread.h>
#endif


#include <time.h> 
#include "math.h"
#ifdef __cplusplus
#include <vector>
#include <vectorclass/vectorclass.h>
#include "vectorclass_util.h"
#endif

MSBG_NAMESPACE_BEGIN

#define UT_MAXPATHLEN 1024

#define UT_CLEAR	(1<<1)
#define UT_ALL		(1<<2)
#define UT_NOSCALEXY	(1<<3)
#define UT_WAIT		(1<<4)
#define UT_VERBOSE	(1<<5)
#define UT_NOERRTRC	(1<<6)

#define UT_FUNCNAME	__FUNCTION__	/* gcc */
#define UT_CLASSNAME  typeid(*this).name()
#define UT_TEMPDIR	"c:/tmp"

#define UT_NEXTTOKEN(tok) ((tok=strtok(NULL," \t")))

#ifdef MIMP_ON_WINDOWS
#define UtStrCaseCmp( s1_, s2_ ) stricmp( s1_,s2_)
#else
#define UtStrCaseCmp( s1_, s2_ ) strcasecmp( s1_,s2_)
#endif

#define UT_SWAP( x, y, tmp ) \
{ \
  tmp = x; \
  x = y; \
  y = tmp; \
}

#define UT_SWAP_( type_, x, y) \
{ \
  type_ tmp__ = x; \
  x = y; \
  y = tmp__; \
}

#define UT_SWAP_PTR( p1, p2 ) \
{ \
  void *tmp__ = p1; \
  p1 = p2; \
  p2 = tmp__; \
}

#define UT_SWAP_PTR_( type_, p1, p2 ) \
{ \
  type_ *tmp__ = p1; \
  p1 = p2; \
  p2 = tmp__; \
}

#define UT_SWAP_DBL( p1, p2 ) \
{ \
  double tmp__ = p1; \
  p1 = p2; \
  p2 = tmp__; \
}

#define UT_NUM_EPS	(1.0E-12)

#define UT_ISPOW2(n_) (!((n_) & ((n_) - 1)))

#define UT_NEXT_POW2( n, n2, k, max ) \
{ \
  for(n2=1, k=0; n2<max; n2=n2<<1, k++) \
  { \
    if(n2>=n) break; \
  } \
}

#define UT_GXYZ( grid_sx, grid_sxy, x, y, z ) ((x)+(y)*((LongInt)(grid_sx))+(z)*(LongInt)(grid_sxy))

#ifdef __cplusplus

namespace UT
{

#define UT_VECTOR_CONTAINS( vec, a ) \
 ( std::find((vec).begin(), (vec).end(), (a)) != (vec).end() )
    
#define UT_VECTOR_APPEND( a, b ) \
  (a).insert((a).end(), (b).begin(), (b).end());

/*-----------------------------------------------------------------------------------*/
/* 										     */
/*-----------------------------------------------------------------------------------*/
template <typename... T>
std::vector<T...> ReservedVector(size_t n=256) {
    std::vector<T...> vec;
    vec.reserve(n);
    return vec;
}

#ifdef MIMP_ON_WINDOWS
/*-----------------------------------------------------------------------*/ 
/* 									 */
/*-----------------------------------------------------------------------*/
class SparsePagedArray
{
  public:

  SparsePagedArray( const char *name, 
		    size_t nMaxElem, size_t elemSz, size_t pageSzLog2 );
  ~SparsePagedArray();

  void *allocElem( LongInt idxElem ); 

  bool isPageCommited( uint32_t iPage ) const
  {
    uint32_t idx = (iPage)>>5;
    return (_bitmap)[idx] & (((uint32_t)1)<<((iPage)&31));
  }

  void setPageCommited( uint32_t iPage ) 
  {
    uint32_t idx = (iPage)>>5;
    _bitmap[idx] |= (((uint32_t)1)<<((iPage)&31));
  }

  void lock( void );

  void unlock( void );

  private:

  char _name[80];
  void *_baseAllocated,
       *_base;
  uint32_t *_bitmap;  // bitmap
  size_t _bitmapLen;
  size_t _nMaxElem,
	 _elemSz,
	 _pageSzLog2,
	 _pageSz,
	 _nMaxPages,
	 _nReservedPages,
	 _nCommitedPages,
	 _nSentinelPages;
  long int _lock;
  LongInt _nLockReq,
          _nLockCol = 0;
};
#endif

/*-----------------------------------------------------------------------*/ 
/* 									 */
/*-----------------------------------------------------------------------*/
class LargeContBuf
{
  public:
    LargeContBuf( size_t initialSz, size_t maxSz, int useLargePages=0 );
    ~LargeContBuf();

    size_t extend( size_t newSz=0, int noErrorTrace=0 );

    void *base( void ) const { return _base; }
    size_t getSize( void ) const { return _actSz; }
    size_t getMaxSize( void ) const { return _reservedPages*_pageSz; }
    size_t getCommitedSize( void ) const { return _commitedPages*_pageSz; }
    static size_t getTotalCommitedSz( void ) { return _totalCommitedSz; }
    //size_t getMaxCommitableExtension( void ) ;

  private:
    void *_base;
    size_t _initialSz,
	   _actSz,
	   _pageSz,
      	   _reservedPages,
	   _commitedPages;

    static size_t _totalCommitedSz;
};

/*-----------------------------------------------------------------------*/ 
/* 									 */
/*-----------------------------------------------------------------------*/
template <typename T, unsigned B>
inline T signextend(const T x)
{
  struct {T x:B;} s;
  return s.x = x;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
inline uint32_t UtFloatToIntBits( float f )
{
  union 
  {
    float fl;
    uint32_t ui;
  }
  u;
  u.fl = f;
  return u.ui;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
inline float UtIntBitsToFloat( uint32_t ui )
{
  union 
  {
    float f_;
    uint32_t ui_;
  }
  u;

  u.ui_ = ui;
  return u.f_;
}


/*-----------------------------------------------------------------------*/ 
/* 									 */
/*-----------------------------------------------------------------------*/
class StrTok
{
  private:
    char *str_,
	 *strOrig_;
    std::vector<char *> tokens_;
  public:
    // return number of tokens
    size_t n() { return tokens_.size(); }
    // return token
    char *get(size_t i) { return i<tokens_.size()?tokens_[i]:NULL; }
    //
    double getFloat(size_t i, double deflt=-1e20 ); 
    // constructor
    StrTok(const char*str, const char *sep);
    ~StrTok();
};

inline void sort3val(float &a, float &b, float &c)
{
  float temp;
  if(a>b){ temp=a; a=b; b=temp; }
  if(a>c){ temp=a; a=c; c=temp; }
  if(b>c){ temp=b; b=c; c=temp; }
}


template <typename T>
void sort3(T& a, T& b, T& c)
{
    if (a < b) {
        if (b < c) {
            return;
        } else if (a < c) {
            std::swap(b, c);
        } else {
            T tmp = std::move(a);
            a = std::move(c);
            c = std::move(b);
            b = std::move(tmp);
        }
    } else {
        if (a < c) {
            std::swap(a, b);
        } else if (c < b) {
            std::swap(a, c);
        } else {
            T tmp = std::move(a);
            a = std::move(b);
            b = std::move(c);
            c = std::move(tmp);
        }
    }
}

#if 0
inline void sort3valFast(float &a, float &b, float &c)
{
  UT_ASSERT0(FALSE); // inaccurate when numbers differ by large amount (e.g. SDF_INFINITY)
  // https://www.geeksforgeeks.org/sort-3-integers-without-using-condition-using-max-function/
  float a0=a,b0=b,c0=c;
  a = std::min(std::min(a0,b0),c0);  // min3
  c = std::max(std::max(a0,b0),c0);  // max3
  b = a0+b0+c0 - (a+c);    // med3
}
#endif


}  // namespace UT

#ifdef __cplusplus
extern "C" {
#endif

#endif
typedef struct __attribute__((aligned(CPU_NO_FALSE_SHARING_SIZE)))
{
  size_t bytes,
	 flops;
}
UtPerfStats;

typedef struct 
{ 
  struct timeval start, stop; 
  float diff;
  
  uint64_t PrevTotal, 
           PrevIdle,
           PrevUser;
  float cpuLoadPercent;
}
UtTimer;

typedef struct 
{ 
  double hrStart,hrStop,hrDiffSec;
}
UtHTimer;

#define UT_RECT_IN_RANGE(x,y,xmin,ymin,xmax,ymax) \
 ((x)>=(xmin)&&(x)<=(xmax)&&(y)>=(ymin)&&(y)<=(ymax))

/*-------------------------------------------------------------------------*/
/* simple inter-thread mutexes						   */
/*-------------------------------------------------------------------------*/
//typedef CRITICAL_SECTION UtMutex;

#ifdef MIMP_ON_LINUX

typedef pthread_mutex_t UtMutex;

#define UtMutexInit( m_ ) pthread_mutex_init( m_, NULL )

#define UtMutexLock(m_) pthread_mutex_lock( m_ )

#define UtMutexUnlock(m_) pthread_mutex_unlock( m_ )

#define UtMutexDelete( m_ ) pthread_mutex_destroy( m_ )

#else

#define UtMutex CRITICAL_SECTION

#define UtMutexInit( m_ ) \
{ \
  InitializeCriticalSection(m_); \
}

#define UtMutexLock(m_) \
{ \
  EnterCriticalSection(m_); \
}

#define UtMutexUnlock(m_) \
{ \
  LeaveCriticalSection(m_); \
}

#define UtMutexDelete( m_ ) \
{ \
  DeleteCriticalSection(m_); \
}

#endif

/*-------------------------------------------------------------------------*/
/* spinlocks */
/*-------------------------------------------------------------------------*/
void *UtCreateSpinlock(  );
void UtAcquireSpinlock( void *s );
void UtReleaseSpinlock( void *s );


/*-------------------------------------------------------------------------*/
/* millisecond timer (elapsed time)					   */
/*-------------------------------------------------------------------------*/

//#define UT_TIMER_WITH_CPULOAD

#ifdef MI_ON_NATWIN
int gettimeofday(struct timeval *tv, struct timezone *tz);
#endif

#ifdef UT_TIMER_WITH_CPULOAD

#define TIMER_START(tm_) \
{ \
  if( gettimeofday( &((tm_)->start), NULL ) !=0 ) \
  { \
    TRCERR(("ERROR: gettimeofday() failed. errno=%d\n",errno)); \
  } \
  (tm_)->PrevTotal = 0;\
  (tm_)->PrevIdle = 0;\
  (tm_)->PrevUser = 0;\
  UtGetAverageCPULoad( NULL, tm_ );\
}

#define TIMER_STOP(tm_) \
{ \
  (tm_)->cpuLoadPercent = UtGetAverageCPULoad( NULL, tm_ );\
  \
  if( gettimeofday( &((tm_)->stop), NULL ) != 0) \
  { \
    TRCERR(("ERROR: gettimeofday() failed. errno=%d\n",errno)); \
  } \
  else \
  { \
    (tm_)->diff = TIMER_TV2MS((tm_)->stop) - \
    		      TIMER_TV2MS((tm_)->start); \
  } \
  TRC(("%s: Avg. CPU load was %g%%\n", UT_FUNCNAME, (tm_)->cpuLoadPercent ));\
} \

#else

#define TIMER_START(tm_) \
{ \
  if( gettimeofday( &((tm_)->start), NULL ) !=0 ) \
  { \
    TRCERR(("ERROR: gettimeofday() failed. errno=%d\n",errno)); \
  } \
}

#define TIMER_STOP(tm_) \
{ \
  if( gettimeofday( &((tm_)->stop), NULL ) != 0) \
  { \
    TRCERR(("ERROR: gettimeofday() failed. errno=%d\n",errno)); \
  } \
  else \
  { \
    (tm_)->diff = TIMER_TV2MS((tm_)->stop) - \
    		      TIMER_TV2MS((tm_)->start); \
  } \
}

#endif

#define GTIMER_START() \
{ \
  UT_ASSERT0(!UtGlobTimerActive);\
  TIMER_START(&UtGlobTimer);\
  UtGlobTimerActive=1;\
}

#define GTIMER_STOP(msg) \
{ \
  UT_ASSERT0(UtGlobTimerActive);\
  TIMER_STOP(&UtGlobTimer);\
  UtGlobTimerActive=0; \
  TRC(("GTIMER [%s:%d] %s => CPU %.2f sec\n",\
	__FILE__,__LINE__,msg?msg:"",\
	(double)TIMER_DIFF_MS(&UtGlobTimer)/1000.));\
}


#define TIMER_STOP_INT(tm_) TIMER_STOP(tm_);

#define TIMER_DIFF_MS(tm_) ((tm_)->diff)

#define HTIMER_START(tm_) UtHTimerStart(tm_)

#define HTIMER_STOP(tm_) UtHTimerStop(tm_)

#define HTIMER_DIFF_US(tm_) ((tm_)->hrDiffSec*1000000.0)

#if 0
#define CTIMER_START(tm_) \
{ \
  (tm_)->tsc = ReadTSC();\
}

#define CTIMER_STOP(tm_) \
{ \
  (tm_)->tsc_diff = ReadTSC() - (tm_)->tsc;\
}

#define CTIMER_DIFF_CLK(tm_) ((tm_)->tsc_diff)
#endif

/*-------------------------------------------------------------------------*/
/* assertion/check macros						   */
/*-------------------------------------------------------------------------*/

// 
// UT_ASSERT_LEVEL_1 : Default
//
// -DUT_ASSERT_LEVEL_0 : Switch off all checks for max. performance
//
// -DUT_ASSERT_LEVEL_2 : Checks with little/moderate performance penalty (mk_assert2)
//
// -DUT_ASSERT_LEVEL_3 : Extended checks that may have larger performance penalty
//

#ifndef UT_ASSERT_LEVEL_0
  #define UT_ASSERT_LEVEL_1
#endif

//
// UT_ASSERT (level 1)
//
#ifdef UT_ASSERT_LEVEL_1

#define UT_ASSERT(exp)  if (likely(exp)) ; \
        else TRCERRR(("Assertion failed.\n"),MI_EASSERT); 

#define UT_ASSERT_FATAL(exp)  if (likely(exp)) ; \
  else \
{ \
  TRCERR(("Assertion failed (fatal) -> Exit\n")); \
  exit(1);\
}

#define UT_ASSERT_NR(exp)  if (likely(exp)) ; \
        else TRCERR(("Assertion failed.\n")); 

#define UT_ASSERT0(exp) UT_ASSERT_FATAL(exp)

#else

#define UT_ASSERT0(exp)  
#define UT_ASSERT(exp)  
#define UT_ASSERT_NR(exp)  
#define UT_ASSERT_FATAL(exp) 

#endif 

//
// UT_ASSERT2 (level 2; extended checks)
//
#ifdef UT_ASSERT_LEVEL_3
#ifndef UT_ASSERT_LEVEL_2
#define UT_ASSERT_LEVEL_2
#endif
#endif

#ifdef UT_ASSERT_LEVEL_2

#if 0
#define UT_ASSERT2(exp)  if (likely(exp)) ; \
        else TRCERRR(("Assertion failed.\n"),MI_EASSERT); 
#endif

#define UT_ASSERT2(exp)  if (likely(exp)) ; \
  else \
{ \
  TRCERR(("Assertion failed (fatal) -> Exit\n")); \
  exit(1);\
}

#define UT_ASSERT2_NF(exp)  if (likely(exp)) ; \
  else \
{ \
  TRCERR(("Assertion failed \n")); \
}

#define UT_ASSERT2_FATAL(exp)  if (likely(exp)) ; \
  else \
{ \
  TRCERR(("Assertion failed (fatal) -> Exit\n")); \
  exit(1);\
}

#else

#define UT_ASSERT2(exp)  
#define UT_ASSERT2_NF(exp)  
#define UT_ASSERT2_FATAL(exp)  

#endif

#define UT_THROW_ASSERTION() TrcThrowError(__FILE__,__LINE__)

// 
// UT_CHECK: for checking user input etc.
//
#define UT_CHECK(exp)  if (likely(exp)) ; \
        else TRCERRR(("Assertion failed.\n"),MI_EASSERT); 

#define UT_CHECK_WARN(exp)  if (likely(exp)) ; \
        else TRCWARN(("Assertion failed.\n")); 

#define UT_CHECK_THROW(exp)  if (likely(exp)) ; \
        else TRCERRT(("Assertion failed.\n"),MI_EASSERT); 

#define UT_CHECK_FATAL(exp)  if (likely(exp)) ; \
  else \
{ \
  TRCERR(("Assertion failed (fatal) -> Exit\n")); \
  exit(1);\
}

#define UT_IMPLIES(A,B) (!(A) || (B))
#define UT_EQUIVALENT(A,B) (!!(A) == !!(B))

#ifndef UINT32_MAX
#define UINT32_MAX  (0xffffffff)
#endif

#ifndef UINT64_MAX
#define UINT64_MAX (18446744073709551615ULL)
#endif

/*-------------------------------------------------------------------------*/
/* Progress Indicator 							   */
/*-------------------------------------------------------------------------*/
typedef struct
{
  LongInt iOld;
  int	trig,doTime;
  float percCompletedOld,percStep,actSpeed;
  void *userInfo;
  UtTimer tm;
  int	tmStepSec;
  time_t tmStampSec;
}
UtProgrInd;

#define UT_PROGRIND_INIT( pgi ) \
{ \
  (pgi)->trig = 0; \
  (pgi)->iOld = -1;    \
  (pgi)->percCompletedOld = -1;    \
  (pgi)->userInfo = NULL; \
  (pgi)->doTime = 0; \
  (pgi)->percStep = 0.05; \
  (pgi)->tmStampSec = (time_t)(-1L); \
  (pgi)->actSpeed = 0; \
  (pgi)->tmStepSec = -1; \
  TIMER_START(&(pgi)->tm); \
}

#define UT_PROGRIND_INIT2( pgi, percStep_, tmStepSec_ ) \
{ \
  UT_PROGRIND_INIT(pgi);\
  (pgi)->doTime = 0; \
  (pgi)->percStep = percStep_; \
  (pgi)->tmStepSec = tmStepSec_; \
  (pgi)->tmStampSec = time(NULL); \
}

#define UT_PROGRIND_CONF( pgi, percStep_, doTime_ ) \
{ \
  (pgi)->doTime = doTime_; \
  (pgi)->percStep = percStep_; \
  if((pgi)->doTime ) \
  { \
    TIMER_START(&(pgi)->tm); \
  } \
}

#define UT_PROGRIND( pgi, i_, iFrom_, iTo_, msg_, callback ) \
{ \
  if( (i_) > (pgi)->iOld ) \
  { \
    (pgi)->iOld = i_; \
    (pgi)->trig = FALSE; \
      double percCompleted = \
    			    (double)((i_) - (iFrom_)) / \
      			    (double)((iTo_) - (iFrom_)); \
      if( percCompleted <= 0.0001 || \
	  percCompleted >= 0.9999 || \
	  percCompleted >= \
	    (pgi)->percCompletedOld + (pgi)->percStep ) \
      { \
	UtProgrIndDisplay( pgi, percCompleted, (char*)(msg_), \
	    		  callback ? callback : GWXsetmsg ); \
        (pgi)->percCompletedOld = percCompleted; \
      } \
  } \
}

#define UT_PROGRIND_LAPSED( pgi )  ((pgi)->trig)

#define UT_USER_BREAK(msg_,rc_) \
{ \
  int rcGw,gwich,gwncnt,gwiflg; \
    rcGw = GWkybrd(&gwich, &gwncnt, &gwiflg, 0); \
    if(rcGw==220 /*GWK_ROOF*/) \
    { \
      if( GWmsgboxl("sure to abort operation ?") == 1) \
      { \
	TRC1(("*** WARNING: operation aborted by user !\n")); \
	  raiseRc( rc_ ); \
      } \
    } \
}

#define UT_USER_BREAK2(msg_,lab_) \
{ \
  int rcGw,gwich,gwncnt,gwiflg; \
    rcGw = GWkybrd(&gwich, &gwncnt, &gwiflg, 0); \
    if(rcGw==GWK_ROOF) \
    { \
      if( GWmsgboxl("sure to abort operation ?") == 1) \
      { \
	TRC1(("*** WARNING: operation aborted by user !\n")); \
	  goto lab_; \
      } \
    } \
}

#define UT_USER_BREAK3(rc_) \
{ \
  int rcGw,gwich,gwncnt,gwiflg; \
    rc_=0; \
  rcGw = GWkybrd(&gwich, &gwncnt, &gwiflg, 0); \
    if(rcGw==220 /*GWK_ROOF*/) \
    { \
      if( GWmsgboxl((char*)"sure to abort operation ?") == 1) \
      { \
	TRC1(("*** WARNING: operation aborted by user !\n")); \
        rc_ = MI_ECANCEL; \
      } \
    } \
}

#define UT_BSEARCH_ARR( f, elsz, n, x, i, j ) \
{ \
  int k_; \
  i=0;j=n-1; \
  while( i < j-1 ) \
  { \
    k_ = (i + j) >> 1; \
    if( x < f[k_*(elsz)] ) \
      j = k_; \
    else \
      i = k_; \
  } \
}

#define UtLineFromStr( str_, pp_ ) \
  strtok_r(str_,"\n",pp_)

#define UT_GETINT16( buf ) \
  (((unsigned long)(*(((unsigned char*)(buf))))     << 8 ) + \
   ((unsigned long)(*((((unsigned char*)(buf)))+1))))

#define UT_GETINT32( buf ) \
  (((unsigned long)(*(((unsigned char*)(buf))))     << 24 ) + \
   ((unsigned long)(*((((unsigned char*)(buf)))+1)) << 16 ) + \
   ((unsigned long)(*((((unsigned char*)(buf)))+2)) << 8  ) + \
   ((unsigned long)(*((((unsigned char*)(buf)))+3))       ))

#define UT_PUTINT32( buf, x ) \
{ \
  *((unsigned char*)(buf))++ =  (unsigned long)(x) >> 24 ; \
  *((unsigned char*)(buf))++ =  (unsigned long)(x) >> 16 ; \
  *((unsigned char*)(buf))++ = (unsigned long)(x) >> 8 ; \
  *((unsigned char*)(buf))++ = (unsigned long)x; \
}

#define UT_FAST_HASH_32( h_ ) /* IN: uint32_t*/ \
{\
    /* 'murmur' hash */\
    h_ ^= h_ >> 16;\
    h_ *= 0x85ebca6b;\
    h_ ^= h_ >> 13;\
    h_ *= 0xc2b2ae35;\
    h_ ^= h_ >> 16;\
}

#define UT_FAST_HASH_32F( i, /* IN: uint32_t*/ \
    		         out     /* OUT: float in [0,1] */\
		        )\
{\
    union { float fres; unsigned int ires; } u_;\
    uint32_t h_=i; \
    UT_FAST_HASH_32(h_);\
    /* http://www.iquilezles.org/www/articles/sfrand/sfrand.htm */\
    u_.ires = ((((uint32_t)(h_))>>9 ) | 0x3f800000);\
    out = u_.fres - 1.0f;\
}

#define UT_FAST_HASH_64( k_ )\
{ \
    k_ ^= k_ >> 33;\
    k_ *= 0xff51afd7ed558ccdULL;\
    k_ ^= k_ >> 33;\
    k_ *= 0xc4ceb9fe1a85ec53ULL;\
    k_ ^= k_ >> 33;\
}

#ifdef __cplusplus

inline uint64_t UtFastHash64(uint64_t k)
{
  UT_FAST_HASH_64( k );
  return k;
}
inline uint32_t UtFastHash32(uint32_t k)
{
  UT_FAST_HASH_32( k );
  return k;
}

inline uint32_t UtHash64To32(uint64_t key)
{
  key = (~key) + (key << 18); // key = (key << 18) - key - 1;
  key = key ^ (key >> 31);
  key = key * 21; // key = (key + (key << 2)) + (key << 4);
  key = key ^ (key >> 11);
  key = key + (key << 6);
  key = key ^ (key >> 22);
  return (uint32_t) key;
}

inline float UtHash64ToFloat( uint64_t k )
{
  uint32_t kl = UtHash64To32(k);
  float f;
  UT_FAST_HASH_32F( kl, f ); 
  return f;
}

inline uint32_t UtHash32Coords4D8( uint32_t seed,
    				    uint32_t x,  // 4D coords (mod 256)
    			            uint32_t y,
				    uint32_t z,
				    uint32_t w )
{
  uint32_t h = seed ^
	       ( ((((uint32_t)(z))&255) <<24) | 
		 ((((uint32_t)(y))&255) <<16) |
		 ((((uint32_t)(x))&255) <<8) |
		 ((((uint32_t)(w))&255) <<0) );       
  UT_FAST_HASH_32( h );
  return h;
}

#if 0
inline uint32_t UtHash32Coords4D16( Vec4ui pos4d )
{
  Vec4uq x =  _mm256_cvtepi32_epi64(pos4d & 65535);
  //x = _mm256_srlv_epi64( x, Vec4uq( 0,16,32,48 ) );
  uint64_t x_[4] __attribute__((aligned(32)));
  x.store_a(x_);
  //uint64_t h = (x_[0]*73856093) ^ (x_[1]*19349663) ^ (x_[2]*83492791) ^ (x_[3]*56598313);
  uint64_t h = x_[0] ^ (x_[1]<<16) ^ (x_[2]<<32) ^ (x_[3]<<48);
  uint32_t h32 = UtHash64To32(h);
  return h32;
}
#else

inline uint32_t UtHash32Coords4D16( Vec4i pos4d )
{
  static const Vec4i hashmul( 73856093, 19349663, 83492791, 56598313 );
  pos4d *= hashmul;
  uint32_t h = vget_x(pos4d) ^ vget_y(pos4d) ^ vget_z(pos4d) ^ vget_w(pos4d);
  return h;
}
#endif

#endif


/*-------------------------------------------------------------------------*/
/* fast fixed size dynamic memory management				   */
/*-------------------------------------------------------------------------*/
typedef struct
{
  void 	   **stack;
  size_t    sp;
  void      *data;
  size_t    nMaxElem;
}
UtFSM;

#define UT_FALLOC( fsm ) \
  ( (fsm)->sp ? ( (fsm)->stack[--(fsm)->sp] ) : NULL)

#define UT_FFREE( fsm, ptr ) \
{ \
  (fsm)->stack[(fsm)->sp++] = ptr; \
}

void UtDeleteFSM( UtFSM *fsm );
UtFSM *UtCreateFSM( size_t maxBlocks, size_t blockSize );

#define UtGetMaxElemFSM( fsm ) (fsm)->nMaxElem;

/*-----------------------------------------------------------------------*/
/* 									 */
/*-----------------------------------------------------------------------*/
extern int UtCheckLevel;
extern int UtGlobDbgCheckCount;
extern UtTimer UtGlobTimer;
extern int UtGlobTimerActive;
extern int UtGlobDbgErrCode;
typedef struct
{
  int dummy;
}
UtGlobDbgErrInfoS;
extern UtGlobDbgErrInfoS UtGlobDbgErrInfo;

/*-----------------------------------------------------------------------*/
/* 									 */
/*-----------------------------------------------------------------------*/
void UtSetTempDir( char *tmpDir );

char *UtGetTempDir( void ); 

float UtGetAverageCPULoad( double *deltaMilliSec, UtTimer *tm );

uint64_t UtGetCurrentThreadId(void);

uint32_t UtAtomicCompareExchange( uint32_t volatile *dest, 
    			          uint32_t exchange, uint32_t comparand,
				  unsigned opt );

uint32_t UtAtomicExchange( uint32_t volatile *dest, 
    			          uint32_t exchange,
				  unsigned opt );

void *UtAllocLargePageMem( size_t size );
void UtFreeLargePageMem( void *address );
int UtFileLastModified( char *path, 		// IN
    			time_t *pTimeModified   // OUT 
			);
int UtGetFileSize( char *fpath, size_t *p_size );

int UtReadWriteBinData2( int do_write, 
    			     const void *data_ptr, size_t data_sz, 
    			     char *fname,
			     int doUnbuffered );

int UtReadWriteBinData( int do_write, 
    			     const void *data_ptr, size_t data_sz, 
    			     char *fname, int fd,
			     unsigned opt );
int UtSleep( unsigned millisec );

static inline float UtFastRand( 
    	int *seed // use as per-thread state (initialize e.g. with seed=HASH(thread_id)+1)
		  // XXX: seed must not be zero !
    )
{
  // http://www.iquilezles.org/www/articles/sfrand/sfrand.htm
    union
    {
        float fres;
        uint32_t ires;
    } u;

    seed[0] *= 16807;

    u.ires = ((((unsigned int)seed[0])>>9 ) | 0x3f800000);
    return u.fres - 1.0f;
}


// 
// Normal random numbers generator - Marsaglia algorithm.
static inline float UtFastRandGauss( int *seed, // random seed
    				     int *state_i, // state: must be initial zero
				     float *state_g,
    				     float mean, float sigma )
{
  float         fac, rsq, v1, v2;
  
  if( *state_i == 0 )
  {
    do
    {
      v1 = 2.0f * UtFastRand(seed) - 1.0f;
      v2 = 2.0f * UtFastRand(seed) - 1.0f;
      rsq = v1*v1 + v2*v2;
    }
    while(( rsq >= 1.0f )||( rsq == 0.0f ));
    fac = sqrtf( -2.0f*logf( rsq ) / rsq );
    *state_g = v1 * fac;
    *state_i = 1;
    return sigma * ( v2 * fac ) + mean;
  }
  else
  {
    *state_i = 0;
    return sigma * (*state_g) + mean;
  }  
}

int UtMemDisclaim( void *addr, size_t size,
    		   int doDecommit,
		   int isPageAligned
		   );
int UtMemProtect( void *addr, size_t size,
    		  int doProtectWrite );
void *UtAllocMemNuma( size_t sizeBytes );
void UtFreeMemNuma( void **ptr );
size_t UtLargeContbufTotalUsedSize( void );
char *UtStrtok(char *s, const char *delim, char **save_ptr);

#define UT_INT2STR( x, base, outbuf )itoa( x, outbuf, base )
char *UtThousandSep( int64_t x, char *outbuf );
int UtGetKeyPressed( int *keyOut, int *codeOut );
int UtParseFloat( char *tok, double *pdbl );
int UtParseInt( char *tok, unsigned opt, LongInt *plong );
int UtParseFloatArray( double *f, int *pLen, char *tok_in, char *delim );
size_t UtGetSystemPageSize( void );
void UtGetSystemMemSize( size_t *p_phys_total, 
    			 size_t *p_phys_avail,
			 size_t *p_commit_avail,  // max. per process 
			 size_t *p_virtual_avail );
int UtOpenLogFile(char *fname, int do_copy_old, int do_append);
int UtOpenLogFile2(char *fname);
int UtCloseLogFile();
int UtCloseLogFile2();

int UtProgrIndDisplay( UtProgrInd *pgi, float percCompleted, char *chbuf0, 
    		     int (*callback)( char *msg, int l) );

int UtSimpleRegexMatch(char	*string,
			char	*pattern,
			int          accept_sub_string,
			int          case_insensitive);

char *MiHexdump2Str( unsigned char *memaddr, size_t memlen, 
		      char *buf, int buflen );
size_t MiStripTrailingCRLF( char *str, size_t slen );
size_t UtStripTrailingWS( char *str, size_t slen ); 
int UtWriteCRLF( FILE *fp );
int UtReadCRLF( FILE *fp );
int UtGetTextLine( FILE *fi, char *line, size_t maxLineLen, unsigned opt,
    		   char **pline, int *plen );
char *UtGetNextLineFromString( char **textBegin, char *textEnd );
char *UtSetFilenameSuffix( char *fnameCompl, 
    			    char *fname, 
			    char *suffix );
int UtCheckFileExists( char *path );
char *UtFileSuffix( char *fname );
char *MiTranspathCyg2Win( char *fname );
char *UtPathCyg2Win( char *pin );
char *UtPathWin2Cyg( char *pin, int maxlen );
int UtShellExec( char *command, char *params, unsigned options,
    		 int *pRetCode );
int UtGetCurDir( char *pbuf, size_t *plen );
int UtSetCurDir( char *pbuf );
int UtDiffFiles( char *fn1, char *fn2 );
int UtDeleteFile( char *fname );
int UtCopyFile( char *src, char *dst, int doFailIfExists );
char *UtReadFileToString( char *path );
int UtCreatePath( char *pathIn );
int UtIsAbsFilePath( char *path );
int UtParseFilePath( char *fullPath, 
    		     char **pDir, char **pFile, char **pSuffix );
int UtStripLastDirSep( char *path );
int UtHTimerStart(UtHTimer *tm);
int UtHTimerStop(UtHTimer *tm);

#ifdef __cplusplus
}
#endif

MSBG_NAMESPACE_END

#endif /* UTIL_H */

