#ifndef GLOBDEF_H
#define GLOBDEF_H

#define MSBG_NAMESPACE MSBG
#ifdef __cplusplus
#define MSBG_NAMESPACE_BEGIN namespace MSBG_NAMESPACE { 
#define MSBG_NAMESPACE_END } 
#else
#define MSBG_NAMESPACE_BEGIN
#define MSBG_NAMESPACE_END
#endif

//#define MSBG_RENDERDENS_BSX_LOG2 5
#define MSBG_RENDERDENS_BSX_LOG2 4
#define MSBG_RENDERDENS_BSX (1<<MSBG_RENDERDENS_BSX_LOG2)

#define DO_USE_GLOB_HEAP 	0
//#define DO_USE_GLOB_HEAP 	1

#if 0
//
// set to undefined to avoid accidental use
// of libc malloc/free
// if needed (e.g. for GNU readline lib) those should 
// be explicitely called
// via MM_libc_malloc/MM_libc_free 
//
#ifndef MIMP_NOREDEFINEMALLOC
#ifdef MIMP_ON_WINDOWS
#define strdup xxxstrdup  
#define malloc xxmalloc   
#define calloc xxcalloc   
#define realloc xxrealloc  
#define free(p) xxfree(p)
#endif
#endif
#else
#define calloc xxcalloc   
#define realloc xxrealloc  
#endif

//#ifdef MI_WITH_ICC
//#include <mathimf.h>
//#endif

#include "windows_linux.h"

#ifdef MI_WITH_MSVC
#include "stdint_msvc.h"
#else
#include <stdint.h>
#endif

#ifdef __cplusplus
#include <cstring>
#else
#include <string.h>
#endif

#include <limits.h>

#ifndef MI_ON_NATWIN
#include <sys/time.h> 
#else
#include <time.h> 
#endif

#include <ctype.h>

#ifdef __cplusplus
extern "C" {
#endif

int getch(void);

#undef CRA_TEST

#define MI_MAX_THREADS		128 
#define MI_MAX_FNAME_LEN	512


#define ASM_MARK_( MARK_ID )						\
__asm__ __volatile__ (									\
					  "\n\t  movq $"#MARK_ID", %%rax"	\
					  : : : "memory" );
#define ASM_MARK_START {ASM_MARK_(1229782938247303441)}
#define ASM_MARK_END {ASM_MARK_(2459565876494606882)}

/* usage: 
 * >mk asm_mark
 *or:
 * > mk foo.o 
 * objdump -d --prefix-addresses -M intel mimp.exe | sed -n '/1111111111111111/,/2222222222222222/p'	
 *
 */

#ifdef MI_WITH_64BIT
typedef int64_t LongInt;
#else
typedef int LongInt;
#endif

#if 1
#ifdef MI_WITH_64BIT
#define strtok_r( _s, _sep, _lasts ) \
        ( *(_lasts) = strtok( (_s), (_sep) ) )
#endif
#endif

#ifndef TRUE
  #define TRUE 1
#endif

#ifndef FALSE
  #define FALSE 0
#endif

#ifndef MAX
  #define MAX(a,b)    ((a)>(b)?(a):(b))
#endif
#ifndef MIN
  #define MIN(a,b)    ((a)<(b)?(a):(b))
#endif


/*-------------------------------------------------------------------------*/
/* global error codes							   */
/*-------------------------------------------------------------------------*/
#define	 MI_OK	0
#define  MI_ERROR	1
#define	 MI_ENOMEM	2
#define	 MI_ESYNTAX	3
#define  MI_ENOTIMPL	4
#define  MI_EINVARG	5
#define MI_ENOTFOUND	6
#define MI_EVERI	7
#define MI_EDONE	8
#define MI_ECANCEL	9
#define MI_EEXISTS	10
#define MI_EOVERFLOW    11
#define MI_EAGAIN	12
#define MI_ECONTINUE	13
#define MI_ENOMATCH	14
#define MI_EMATH	15
#define MI_EQUIT	16
#define MI_EIO		17
#define MI_EOF		18
#define MI_EASSERT	19
#define MI_EINTERNAL	20
#define MI_EOUTOFRANGE	21  
#define MI_ENUMERICAL	22
#define MI_EOUTOFGAMUT  23
#define MI_EFATAL	24
#define MI_ENOTAVAIL	25
#define MI_ERAYMISS	26
#define MI_EDIFFER	27
//#define MI_ECRITICAL	28
#define MI_EINVPIX	29
#define MI_ETIMEOUT	30
#define MI_ELOCKED	31

#define SUCCESS 0

#define ONE_KB 1024
#define ONE_MB (1024*ONE_KB)
#define ONE_GB (((LongInt)1024)*ONE_MB)

#define SWAP_PTR( p1, p2 ) \
{ \
  void *tmp__ = p1; \
  p1 = p2; \
  p2 = tmp__; \
}

#define SWAP_PTR_( type_, p1, p2 ) \
{ \
  type_ *tmp__ = p1; \
  p1 = p2; \
  p2 = tmp__; \
}

/*--------------------------------------------------------------------*/
/* fixed bit width integral types		                      */
/*								      */	
/*    CRA_UINT32 						      */
/*								      */
/*--------------------------------------------------------------------*/
#define CRA_UINT32_MAX 0xFFFFFFFFu

#if UINT_MAX == 0xFFFFFFFFu
  typedef unsigned int CRA_UINT32;                   
#elif ULONG_MAX == 0xFFFFFFFFu
  typedef unsigned long CRA_UINT32;                
#else
  #error "We need an unsigned int type with exactly 4 bytes"
#endif

#include "log.h"

/*--------------------------------------------------------------------*/
/* utility macros						      */	
/*--------------------------------------------------------------------*/
#define ARRAY_LENGTH(Array) (int)((sizeof(Array) / sizeof(Array[0])))

#define BYTE_OFFSET(ptr, offset)  \
  ((void *)((unsigned char *)(void *)(ptr) + (offset)))
#define BYTE_OFFSET_NEG(ptr, offset)  \
  ((void *)((unsigned char *)(void *)(ptr) - (offset)))  

#define ALIGN(i,a) (((i)%(a))?((i)+((a)-((i)%(a)))):(i))

#define ALIGN2(i,a) (((i)&((a)-1))?((i)+((a)-((i)&((a)-1)))):(i))

#define IS_ALIGNED(i,a) ((uint64_t)(i) & ((a)-1) == 0)

#define raiseRc(rc) { rcThis = rc; goto rcCatch; }

#define FLAG_SET( word, flag ) ((word) |= (flag))
#define FLAG_RESET( word, flag ) ((word) &= (~(flag)))
#define FLAG_COPY( word, flag, trueFalse ) \
{\
  FLAG_RESET( word, (flag) );\
  if((trueFalse)) \
  {\
    FLAG_SET(word, (flag) );\
  }\
}
#define FLAG_IS_SET( word, flag ) ((word) & (flag))
#define FLAG_ONOFF( word, onoff_, flag ) \
{ \
  if( onoff_ ) \
    FLAG_SET( word, flag ); \
  else \
    FLAG_RESET( word, flag ); \
}

/*--------------------------------------------------------------------*/
/* GrWin shortcuts						      */
/*--------------------------------------------------------------------*/

// 
// redefine GW calls to generic 'GWX' implementation
//

#include "gwx.h"

#ifndef MIMP_NOREDEFINE_GWX

#define GWsavebmp( NB, FN,  l) \
    	GWXsavebmp( NB, FN,  l)

#define GWline2( X,  Y)\
        GWXline2( X,  Y)

#define GWmove2( X,  Y) \
        GWXmove2( X,  Y)

#define GWcapvec( X1,  Y1, X2, Y2, TEXT,  l) \
        GWXcapvec( X1,  Y1, X2, Y2, TEXT,  l)

#define GWcaplin(X1, Y1, X2, Y2, TEXT,  l)\
        GWXcaplin(X1, Y1, X2, Y2, TEXT,  l)

#define GWclose( NW) \
        GWXclose( NW)

#define GWclear(K) \
        GWXclear(K)

#if 0
#define GWrefresh() \
        GWXrefresh()
#endif

#define GWpolylin( POINTS, N ) \
        GWXpolylin( POINTS, N )

#define GWpolygon(POINTS,  N,  MF) \
        GWXpolygon(POINTS,  N,  MF)

#define GWquitx(MQ) \
        GWXquitx(MQ)

#define GWsetogn( ISRC,  IDST) \
        GWXsetogn( ISRC,  IDST)

#define GWellipse( X1,  Y1,  X2,  Y2) \
        GWXellipse( X1,  Y1,  X2,  Y2)

#define GWsrect( X1,  Y1,  X2,  Y2,  K) \
        GWXsrect( X1,  Y1,  X2,  Y2,  K)

#define GWsetpxl( X,  Y,  K) \
        GWXsetpxl( X,  Y,  K)

#define GWline( X1,  Y1,  X2,  Y2) \
        GWXline( X1,  Y1,  X2,  Y2)

#define GWrect( X1,  Y1,  X2,  Y2) \
        GWXrect( X1,  Y1,  X2,  Y2)

#define GWmakebmp( NB,  IW,  IH,  IBC, IBITS) \
        GWXmakebmp( NB,  IW,  IH,  IBC, IBITS)

#define GWerase( N,  LRF) \
        GWXerase( N,  LRF)

#define GWclipimg( X1,  Y1,  X2,  Y2) \
        GWXclipimg( X1,  Y1,  X2,  Y2)

#define GWshowwn( NW,  IS) \
        GWXshowwn( NW,  IS)

#define GWinput(TITLE,  l1, TXT,  l2)\
        GWXinput(TITLE,  l1, TXT,  l2)

#define GWarrange( M) \
        GWXarrange( M)

#define GWgetpen(IPC, IPS, IPW, MX) \
        GWXgetpen(IPC, IPS, IPW, MX)

#define GWcolor( K,  IDST) \
        GWXcolor( K,  IDST)

#define GWopenx( NW,  IW,  IH,  IFC,  IBC,  M, FN,  l) \
        GWXopenx( NW,  IW,  IH,  IFC,  IBC,  M, FN,  l)

#define GWcappnt(X, Y, TEXT, l) \
    	GWXcappnt(X, Y, TEXT, l)

#define GWcaprect(X1, Y1, X2, Y2, TEXT, l) \
        GWXcaprect(X1, Y1, X2, Y2, TEXT, l)

#define GWdelbmp( NM) \
        GWXdelbmp( NM)

#define GWevent(X, Y) \
        GWXevent(X, Y)

#define GWfiledlg(FN,  l) \
        GWXfiledlg(FN,  l)

#define GWgetwn(X1, Y1, X2, Y2) \
        GWXgetwn(X1, Y1, X2, Y2);

#define GWkybrd(ICH, NCNT, IFLG,  M)\
        GWXkybrd(ICH, NCNT, IFLG,  M)

#define GWputbmp( NB,  X,  Y,  IBK) \
        GWXputbmp( NB,  X,  Y,  IBK)

#define GWputtxt( X,  Y, TXT,  l) \
        GWXputtxt( X,  Y, TXT,  l)

#define GWmsgbox(TXT,l) \
        GWXmsgbox(TXT,l)

#define GWselect(NW) \
        GWXselect(NW)
  
#define GWsetbmp( NB,  W,  H,  MX,  ITR,  IOF) \
        GWXsetbmp( NB,  W,  H,  MX,  ITR,  IOF)

#define GWsetpen( IPC,  IPS,  IPW,  MX) \
        GWXsetpen( IPC,  IPS,  IPW,  MX)

#define GWsettxt( H,  A,  IO,  K,  KB,  FACE,  l) \
        GWXsettxt( H,  A,  IO,  K,  KB,  FACE,  l)

#define GWkrgb( IR,  IG,  IB) \
    	GWXkrgb( IR,  IG,  IB)

#define GWinitx( IRB,  IX,  IY,  IW,  IH,  MA,  MM, MZ,  ND) \
    	GWXinitx( IRB,  IX,  IY,  IW,  IH,  MA,  MM, MZ,  ND);

#define GWsetmsg(t,l) \
	GWXsetmsg(t,l)

#endif


// 'l' variants
#define GWopenxl(a_, b_, c_, d_, e_, f_, g_ ) \
        GWopenx(a_, b_, c_, d_, e_, f_, g_, strlen(g_) ) 

#define GWsettxtl( H,  A,  IO,  K,  KB,  FACE ) \
	GWsettxt( H,  A,  IO,  K,  KB,  (char*)(FACE),  strlen(FACE) )

#define GWgettxtl(W, H, X, Y, TXT ) \
	GWgettxt(W, H, X, Y, TXT, strlen(TXT))

#define GWputtxtl( X,  Y,  TXT ) \
  	GWputtxt( X,  Y,  TXT, strlen(TXT))

#define GWsetmsgl(TXT) \
	GWsetmsg(TXT, strlen(TXT) )

#define GWmsgboxl(TXT)\
	GWmsgbox(TXT, strlen(TXT))

#define GWcappntl( X,  Y,  TEXT )\
	GWcappnt( X,  Y,  TEXT, strlen(TEXT) )

#define GWcaprectl( X1,  Y1,  X2,  Y2,  TEXT )\
 	GWcaprect( X1,  Y1,  X2,  Y2,  TEXT, strlen(TEXT) )

#define GWpausel(TXT)\
  	GWpause(TXT,strlen(TXT))

/*--------------------------------------------------------------------*/
/*  Compiler runtime optimizations		                      */
/*--------------------------------------------------------------------*/

/* likely() and unlikely():                                                
   Use them to indicate the likelyhood of the truthfulness of a condition. 
   This serves two purposes - newer compilers (e.g. gcc) will be able to   
   optimize for branch predication, which could yield siginficant          
   performance gains in frequently executed sections of the code, and the  
   other reason to use them is for documentation  */

# if defined(__GNUC__) && (__GNUC__ >= 3)
#  define likely(x)   __builtin_expect(!!(x), 1)
#  define unlikely(x) __builtin_expect(!!(x), 0)
# else
#  define likely(x)   (x)
#  define unlikely(x) (x)
# endif

/*--------------------------------------------------------------------*/
/* memory allocation macros					      */	
/*--------------------------------------------------------------------*/
#include "mm.h"

#define AELEM( pa_, elemsz_, i_ ) ((pa_)+(elemsz_)*((LongInt)(i_)))
    
#define SWAP( x, y, tmp ) \
{ \
  tmp = x; \
  x = y; \
  y = tmp; \
}

/*--------------------------------------------------------------------*/
/* string handling utilities					      */	
/*--------------------------------------------------------------------*/

inline void STRNZCPY( char *dst, char *src, int n )
{
  #if defined(__GNUC__) && !defined (__clang__)
  #pragma GCC diagnostic push 
  #pragma GCC diagnostic ignored "-Wstringop-truncation" 
  #endif

  strncpy( dst, src, n ); 
  dst[n] = '\0'; 

  #if defined(__GNUC__) && !defined (__clang__)
  #pragma GCC diagnostic pop
  #endif
}

#define STRCPY_CHK( dest, maxlen, src, srclen ) \
{ \
  if( (srclen) > (maxlen) ) \
  { \
    TRCERR(("ERROR: string overflow %d > %d '%.*s'\n",\
    	(int)(srclen),(int)(maxlen),MIN(srclen,50),src)); \
    raiseRc( CRA_EOVERFLOW ); \
  } \
  else \
  { \
    memcpy( dest , src, srclen ); \
    (dest)[srclen] = '\0'; \
  }  \
}

#define STRCCAT( dest, destlen_, maxlen, src, srclen_ ) \
{ \
  size_t destlen = destlen_ >= 0 ? destlen_ : strlen( dest ), \
	 srclen = srclen_ >= 0 ? srclen_ : strlen( src ); \
 \
  if( (destlen) + (srclen) > (maxlen) ) \
  { \
    TRCERR(("ERROR: string overflow %d > %d '%.*s'\n",\
    	(int)((srclen)+(destlen)),(int)(maxlen),MIN(srclen,50),src)); \
    raiseRc( MI_ERROR ); \
  } \
  else \
  { \
    memcpy( (dest) + (destlen), src, srclen ); \
    destlen += (srclen); \
    dest[destlen] = '\0'; \
  } \
}

#ifdef MI_WITH_ICC
#define USE(x) \
{ \
/* #pragma warning ( disable: 592 ) */\
  x=0;\
}
#else
#define USE(x) (void)x
#endif

//#define USE(param)           ((param)=(param))  
//#define USE(param)           ((param)=0)  

#define TIMER_TV2MS( tv ) \
  ((double)((tv).tv_sec) * 1000.0 + ((double)((tv).tv_usec)/1000.0))

#define CLOCK_START() \
{ \
  globClockStart = clock(); \
    if( globClockStart==(clock_t)(-1)) TRCERR(("clock() failed (-1)")); \
}

#define CLOCK_STOP() \
{ \
  globClockStop = clock(); \
  if( globClockStop==(clock_t)(-1)) \
  { \
    TRCERR(("clock() failed (-1)")); \
    globClockDiffMs = 0; \
  } \
  else \
  { \
    globClockDiffMs = \
     1000.0*(((float)globClockStop - (float)globClockStart) / \
   	  (float)CLOCKS_PER_SEC); \
    if( globClockDiffMs < 0 ) \
    { \
      TRCERR(("negative timediff: %f\n",(float)globClockDiffMs)); \
      globClockDiffMs = 0; \
    } \
  } \
}

#define CLOCK_DIFF_MS() globClockDiffMs

#define CLOCK_STOP_TRC(tl_,msg) \
{ \
  CLOCK_STOP(); \
  TRCL(tl_,("CLOCK: '%s'  >>> CPU = %.2f ms\n", \
	msg, (float)CLOCK_DIFF_MS())); \
}

extern clock_t globClockStart,
	       globClockStop;

extern float   globClockDiffMs;

extern int MiGlobalProfileVersion;

extern int	      vlevel;

#ifndef MI_WITH_64BIT
extern int errno;
#endif

#ifdef __cplusplus
}
#endif


#endif
