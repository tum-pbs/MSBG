/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#ifndef MM_H
#define MM_H

#ifndef MIMP_ON_WINDOWS
#define MEMORY_ALLOCATION_ALIGNMENT 16
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define MM_MEMALIGN	32

#define CPU_SIMD_WIDTH			32		// AVX
#define CPU_CACHELINE_SZ 		64
#define CPU_NO_FALSE_SHARING_SIZE	64
  
#define CPU_BYPASS_CACHE_SIZE		100000000	

#define MM_SYSTEM_PAGE_SIZE	4096

#define MM_MK_UID(uid)		MM_MakeUserId(__FILE__,__LINE__,uid)
#define MM_MK_TMP_UID()		MM_MakeUserId(__FILE__,__LINE__,\
    					      mm_glob_heap_dummy_tag)

extern char mm_glob_heap_dummy_tag[];

#define MM_DUMMYTAG		mm_glob_heap_dummy_tag

#define DO_OVERRIDE_CPP_NEW	0

#define MM_CHK_START_DIFF	(1<<2)
#define MM_CHK_STOP_DIFF	(1<<3)

//
// (re)-define active allocation functions
//

#if DO_USE_GLOB_HEAP

#define MM_malloc(size,tag) 	   MM_GlobHeapAlloc(size,tag)
#define MM_calloc(n,size,tag)      MM_GlobHeapCAlloc(n,size,tag)
#define MM_realloc(ptr,size,tag)   MM_GlobHeapReAlloc(ptr,size,tag)
#define MM_strdup(str) 	   	   MM_GlobHeapStrdup(str,MM_DUMMYTAG)
#define MM_free(ptr) 	           MM_GlobHeapFree(ptr,MM_DUMMYTAG) 

#else

#if 1

#define MM_malloc(size,tag) 	   MM_libc_aligned_malloc(size)
#define MM_calloc(n,size,tag)      MM_libc_aligned_calloc(n,size)
#define MM_realloc(ptr,size,tag)   xxrealloc
#define MM_strdup(str) 	           MM_libc_aligned_strdup(str)
#define MM_free(ptr) 		   MM_libc_aligned_free(ptr)

#else
#define MM_malloc(size,tag) 	   malloc(size.MM_MEMALIGN)
#define MM_calloc(n,size,tag)      calloc(n,size)
#define MM_realloc(ptr,size,tag)   realloc(ptr,size)
#define MM_strdup(str) 	           strdup(str)
#define MM_free(ptr) 		   free(ptr)
#endif

#endif

#define MM_LARGE_PAGES		(1<<0)
#define MM_DYNAMIC_EXTEND	(1<<1)
#define MM_DIAGMODE		(1<<2)
#define MM_NONVERBOSE		(1<<3)
#define MM_INIT_RANDOM		(1<<4)

#define CLEAROBJ(objptr) \
{ \
  memset((void*)objptr,0,sizeof(*objptr));\
}

#define TRCALLOC(size_)\
{ \
  if((size_)>=ONE_GB) \
  { \
    TRC(("ALLOCMEM %.1f MB  (%s:%d)\n",\
	(size_)/(1024*1024.),__FILE__,__LINE__));\
  } \
}

#ifndef __cplusplus

#define ALLOCMEM( newptr_, size_) ALLOCMEM_TAG_(newptr_,void,size_,NULL)
#define ALLOCMEM_TAG( newptr_, size_,tag_) ALLOCMEM_TAG_(newptr_,void,size_,tag_)
#define ALLOCARR( objPtr_, nObj_) ALLOCARR_( objPtr_, void, nObj_) 
#define ALLOCARR0( objPtr_, nObj_) ALLOCARR0_TAG_( objPtr_, void, nObj_, NULL) 
#define ALLOCARR0_TAG( objPtr_, nObj_,tag_) ALLOCARR0_TAG_( objPtr_, void, nObj_, tag_) 
#define ALLOCOBJ0(objPtr_) ALLOCOBJ0_TAG_(objPtr_,void,NULL) 
#define ALLOCOBJ0_TAG(objPtr_,tag_) ALLOCOBJ0_TAG_(objPtr_,void,tag_) 
#define RALLOCMEM( newptr_, ptr_, size_) \
{ \
  char tag[64];\
  MM_MK_UID(tag);\
  newptr_ = MM_realloc( ptr_, size_, tag); \
  if(!(newptr_)) \
  { \
    TRCERR(("out of virtual memory: MM_realloc(%lu) failed -> exit.\n", \
	  (unsigned long)(size_))); \
    exit(1); \
  } \
}
#define RALLOCARR( objPtr_, nObj_ ) \
  RALLOCMEM( objPtr_, objPtr_, (nObj_)*sizeof(*(objPtr_)))

#define ALLOCOBJ( objPtr_ ) \
  ALLOCARR( objPtr_, 1 )

#define ALLOCOBJS( objPtr_, nObj_ ) \
  ALLOCARR( objPtr_, nObj_ )

#endif /* !__cplusplus */

#define ALLOCMEM_TAG_( newptr_, type_, size_, tag_) \
{ \
  TRCALLOC(size_); \
  char tag[64];\
  if(!tag_)\
  {\
    MM_MK_UID(tag);\
  }\
  newptr_ = (type_*)MM_malloc(size_,tag_?tag_:tag); \
  if(!(newptr_)) \
  { \
    TRCERR(("out of virtual memory: MM_malloc(%.0f) failed -> exit.\n", \
	  (double)(size_))); \
    exit(1); \
  } \
}
#define ALLOCMEM_(newptr_,type_,size_) \
  ALLOCMEM_TAG_(newptr_,type_,size_,NULL)

#define ALLOCMEM2_TAG_( newptr_, type_, size_, noexit, tag_ ) \
{ \
  TRCALLOC(size_); \
  char tag[64];\
  if(!tag_)\
  {\
  MM_MK_UID(tag);\
  }\
  newptr_ = (type_*)MM_malloc(size_,tag_?tag_:tag); \
  if(!(newptr_)) \
  { \
    TRCERR(("out of virtual memory: MM_malloc(%.0f) failed -> exit.\n", \
	  (double)(size_))); \
    exit(1); \
  } \
}

#define ALLOCMEM2_(newptr_,type_,size_) \
  ALLOCMEM2_TAG_(newptr_,type_,size_,NULL)

#define ALLOCARR0_TAG_( objPtr_, type_, nObj_, tag_) \
{ \
  char tag[64];\
  if(!tag_)\
  {\
    MM_MK_UID(tag);\
  }\
  TRCALLOC(sizeof(*(objPtr_))*(size_t)(nObj_)); \
  objPtr_ = (type_*)MM_calloc( nObj_, sizeof(*(objPtr_)),tag_?tag_:tag); \
  if(!(objPtr_)) \
  { \
    TRCERR(( \
	  "out of virtual memory: MM_calloc(%ld,%lu) failed -> exit.\n", \
	  (unsigned long)nObj_,(unsigned long)(sizeof(*(objPtr_))))); \
    exit(1); \
  } \
}
#define ALLOCARR0_( objPtr_, type_, nObj_) \
  	ALLOCARR0_TAG_( objPtr_, type_, nObj_, NULL)


#define ALLOCARR_( objPtr_, type_, nObj_ ) \
  ALLOCMEM_( objPtr_, type_, (nObj_)*sizeof(*(objPtr_)))

#define ALLOCARR_ALIGNED_( objPtr_, type_, nObj_ ) \
  ALLOCMEM_ALIGNED_( objPtr_, type_, (nObj_)*sizeof(*(objPtr_)), CPU_CACHELINE_SZ )

#define ALLOCARR_ALIGNED0_( objPtr_, type_, nObj_ ) \
{ \
  ALLOCARR_ALIGNED_( objPtr_, type_, nObj_ )\
  memset((void*)objPtr_, 0, (nObj_)*sizeof(*(objPtr_)) );\
}

#define ALLOCOBJ0_TAG_( objPtr_, type_, tag_ ) \
  ALLOCARR0_TAG_( objPtr_, type_, 1, tag_ )

#define ALLOCOBJ0_( objPtr_, type_ ) \
  ALLOCARR0_( objPtr_, type_, 1 )

#define FREEMEM( objptr_ ) \
{ \
  if( objptr_ ) \
  { \
    MM_free( objptr_); \
    objptr_ = NULL; \
  } \
}

// 
// Aligned memory management
//
#define ALLOCMEM_ALIGNED_( newptr_, type_, size_, align_ )\
{\
  uint8_t *p0_,*p_;\
  UT_ASSERT_FATAL(((align_)>=16)&&((align_)<=255) && \
                 !(((align_) & ((align_) - 1)))); /* check power of two*/\
  ALLOCMEM_( p0_, uint8_t, ((size_) + (align_)));\
  p_ = (uint8_t*)(((uintptr_t)(p0_) + (align_)) & (~(uintptr_t)((align_) - 1)));\
  *(p_-1) = p_-p0_;\
  newptr_ = (type_*)p_;\
}

#define ALLOCMEM_ALIGNED0_( newptr_, type_, size_, align_ )\
{\
  ALLOCMEM_ALIGNED_( newptr_, type_, size_, align_ );\
  memset((void*)newptr_,0,size_);\
}

#define ALLOCMEM_ALIGNED2_( newptr_, type_, size_, align_, opt_noexit, tag_)\
{\
  uint8_t *p0_,*p_;\
  UT_ASSERT_FATAL(((align_)>=16)&&((align_)<=255) && \
                 !(((align_) & ((align_) - 1)))); /* check power of two*/\
  ALLOCMEM2_TAG_( p0_, uint8_t, ((size_) + (align_)), opt_noexit, tag_ );\
  p_ = (uint8_t*)(((uintptr_t)(p0_) + (align_)) & (~(uintptr_t)((align_) - 1)));\
  *(p_-1) = p_-p0_;\
  newptr_ = (type_*)p_;\
}

#define FREEMEM_ALIGNED( p_ ) \
{\
  if(p_)\
  {\
    void *p0_ = (uint8_t*)(p_) - *(((uint8_t *)p_)-1); \
    FREEMEM(p0_);\
    p_=NULL;\
  }\
}

#define ALLOCMEM_LARGE_ALIGNED_( newptr_, type_, size_, align_, opt_noexit, tag_)\
{\
  uint8_t *p0_,*p_;\
  UT_ASSERT_FATAL(((align_)>=16)&&((align_)<=4096) && \
                 !(((align_) & ((align_) - 1)))); /* check power of two*/\
  ALLOCMEM2_TAG_( p0_, uint8_t, ((size_) + (align_)), opt_noexit, tag_ );\
  p_ = (uint8_t*)(((uintptr_t)(p0_) + (align_)) & (~(uintptr_t)((align_) - 1)));\
  *(((uint16_t*)p_)-1) = p_-p0_;\
  newptr_ = (type_*)p_;\
  UT_ASSERT0(\
      ALIGN( ((uint64_t)(newptr_)), ((uint64_t)(align_))) == (uint64_t)(newptr_));\
}
#define FREEMEM_LARGE_ALIGNED( p_ ) \
{\
  if(p_)\
  {\
    void *p0_ = (uint8_t*)(p_) - *(((uint16_t *)p_)-1); \
    FREEMEM(p0_);\
    p_=NULL;\
  }\
}

int MM_GlobHeapInit( size_t total_sz, unsigned options );
void *MM_GlobHeapAlloc( size_t sz, const char *caller );
void *MM_GlobHeapCAlloc( size_t n_obj, size_t sz_obj, const char *caller );
char *MM_GlobHeapStrdup (const char *s, const char *caller); 
void *MM_GlobHeapReAlloc( void *oldptr, size_t sz, const char *caller );
int MM_GlobHeapFree( void *ptr, const char *caller );
void *MM_GlobHeapThreadAlloc( int tid /*thread id*/, 
    			      size_t sz, const char *caller );
void MM_GlobHeapThreadFree( int tid /*thread id*/, void *ptr );
int MM_GlobHeapTest( size_t total_sz );
int MM_GlobHeapClose( void );
int MM_GlobHeapCheck( unsigned opt );
int MM_GlobHeapLeakDetectorCheckPoint( 
    		unsigned opt   // MM_CHK_START_DIFF | MM_CHK_STOP_DIFF
		);
int MM_GlobHeapResetRecentPeak( void );
int MM_GlobHeapInfo( size_t *p_total_used_sz, 
    		       size_t *p_peak_used_sz,
		       size_t *p_recent_peak_used_sz, 
		       size_t *p_total_heap_sz,
		       size_t *p_largest_free_sz		       
		       );
char *MM_MakeUserId( const char *file, int line,
    		       char *uid // OUT
		      );
char *MM_Str2UserId( const char *str,
    		       char *uid // OUT
		      );
void *MM_libc_malloc(size_t s);
void  MM_libc_free(void *p);
char *MM_libc_strdup(const char *s);
void *MM_libc_aligned_malloc(size_t s);
void *MM_libc_aligned_calloc(size_t num, size_t size);
void  MM_libc_aligned_free(void *p);
char *MM_libc_aligned_strdup(const char *s);
void *MM_libc_mm_malloc(size_t s, size_t align);
void MM_libc_mm_free(void *p);
void MM_Lock(void);
void MM_Unlock(void);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
#if DO_OVERRIDE_CPP_NEW
#include <cstddef>
void* operator new(size_t size, const char* file, int line);
void* operator new[](size_t size, const char* file, int line);
//void operator delete(void* pointer, const char* file, int line) throw();
//void operator delete[](void* pointer, const char* file, int line) throw();

#ifndef MIMP_NOREDEFINEMALLOC
#define MM_NEW new(__FILE__, __LINE__)
#define new MM_NEW
#endif

#endif
#endif


#endif /* MM_H */

