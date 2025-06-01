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
#include <unistd.h>

#define MIMP_NOREDEFINEMALLOC
#include "globdef.h"

#if DO_OVERRIDE_CPP_NEW

void* operator new(size_t size, const char* file, int line)
{
  TRCALLOC(size); 
  char tag[64];
  MM_MakeUserId(file,line,tag);  
  void *newptr_ = MM_malloc(size,tag); 
  if(!(newptr_)) 
  { 
    TRCERR(("out of virtual memory: MM_malloc(%.0f) failed -> exit.\n", \
	  (double)(size)));
    //throw std::bad_alloc();
    exit(1); 
  } 
  return newptr_;
}

void* operator new[](size_t size, const char* file, int line)
{
  TRCALLOC(size); 
  char tag[64];
  MM_MakeUserId(file,line,tag);  
  void *newptr_ = MM_malloc(size,tag); 
  if(!(newptr_)) 
  { 
    TRCERR(("out of virtual memory: MM_malloc(%.0f) failed -> exit.\n", \
	  (double)(size))); 
    //throw std::bad_alloc();
    exit(1); 
  } 
  return newptr_;
}

void* operator new(size_t size) throw(std::bad_alloc)
{
    return operator new(size, (char*)"N/A", 0);
}

void* operator new[](size_t size) throw(std::bad_alloc)
{
    return operator new[](size, (char*)"N/A", 0);
}

void* operator new(size_t size, const std::nothrow_t&) throw()
{
    return operator new(size, (char*)"N/A", 0);
}

void* operator new[](size_t size, const std::nothrow_t&) throw() 
{
    return operator new[](size, (char*)"N/A", 0);
}

void operator delete(void* pMemory) throw()
{
    MM_free(pMemory);
}

void operator delete[](void* pMemory) throw()
{
    MM_free(pMemory);
}

void operator delete(void* pMemory, const std::nothrow_t&) throw()
{
    MM_free(pMemory);
}

void operator delete[](void* pMemory, const std::nothrow_t&) throw() 
{
    MM_free(pMemory);
}
#endif
