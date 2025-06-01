/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <inttypes.h>
#include <unistd.h>
#ifdef MIMP_ON_LINUX
#include <mm_malloc.h>
#endif

//#define MIMP_NOREDEFINEMALLOC
#include "globdef.h"

#include "util.h"

char	mm_glob_heap_dummy_tag[100];


// thread-local storage for on-the-fly tag creation
char    mm_glob_heap_temp_tag[100];  
#pragma omp threadprivate( mm_glob_heap_temp_tag )

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void *MM_libc_malloc(size_t s)
{
  return malloc(s);
}

void MM_libc_free(void *p)
{
  free(p);
}

char *MM_libc_strdup(const char *s)
{
  return strdup(s);
}

void *MM_libc_mm_malloc(size_t s, size_t align)
{
  return _mm_malloc(s,align);
}

void MM_libc_mm_free(void *p)
{
  _mm_free(p);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void *MM_libc_aligned_malloc(size_t s)
{
  return _mm_malloc(s,MM_MEMALIGN);
}

void *MM_libc_aligned_calloc(size_t num, size_t size) 
{
  size_t total_size = num * size;
  void *ptr = MM_libc_aligned_malloc(total_size);
  if (ptr == NULL) return NULL;
  memset(ptr, 0, total_size); 
  return ptr;
}

void MM_libc_aligned_free(void *p)
{
  _mm_free(p);
}

char *MM_libc_aligned_strdup(const char *s)
{
  if (s == NULL) return NULL;
  size_t szNeeded = strlen(s) + 1;
  char *copy = MM_libc_aligned_malloc(szNeeded);
  if (copy == NULL) return NULL;
  memcpy(copy, s, szNeeded);
  return copy;
}


#define MM_UID_SZ	   13

#define MM_MAKE_UID( file, line, buf )\
{\
  int i=MM_UID_SZ;\
  ((uint8_t *)(uid))[--i]=((unsigned)(line))&255;\
  ((uint8_t *)(uid))[--i]=((unsigned)(line))>>8;\
  int l = strlen(file);\
  while(l>0&&i>0) uid[--i] = file[--l];\
  while(i>0) uid[--i] = 0;\
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
char *MM_MakeUserId( const char *file, int line,
    		       char *uid // OUT
		      )
{
  MM_MAKE_UID(file,line,uid);
  return uid;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
char *MM_Str2UserId( const char *str,
    		       char *uid // OUT
		      )
{
  int len=strlen(str);
  for(int i=0;i<MM_UID_SZ;i++) uid[i] = i<len?str[i]:'_';
  uid[MM_UID_SZ] = 0;
  return uid;
}

