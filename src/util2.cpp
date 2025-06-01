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
#ifdef MIMP_ON_LINUX
#include <sys/mman.h> 
#endif
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <errno.h>
#include "globdef.h"
#include "mtool.h"
#include "util.h"

MSBG_NAMESPACE_BEGIN

namespace UT
{

/*=========================================================================
 *
 *
 * 	class UT:SparsePagedArray
 *
 *
 * =======================================================================*/

#ifdef MIMP_ON_WINDOWS

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
SparsePagedArray::SparsePagedArray( const char *name, 
    				    size_t nMaxElem, size_t elemSz, 
    				    size_t pageSzLog2 ) :
  _base(NULL),
  _bitmap(NULL),
  _nMaxElem(nMaxElem),
  _elemSz(elemSz),
  _pageSzLog2(pageSzLog2),
  _pageSz(1<<pageSzLog2),
  _nReservedPages(0),
  _nCommitedPages(0),
  _nSentinelPages(0),
  _lock(0),
  _nLockReq(0),
  _nLockCol(0)
{
  TRC(("%s %.0f \n",UT_FUNCNAME,(double)nMaxElem));

  UT_ASSERT0(_nSentinelPages==0);

  UtTimer tm;
  TIMER_START(&tm);

  UT_ASSERT0(strlen(name)<ARRAY_LENGTH(_name)-1);
  strcpy(_name,name);

  SYSTEM_INFO sSysInfo;
  GetSystemInfo(&sSysInfo);
  int sysPageSz = sSysInfo.dwPageSize;
  TRC(("System page size = %d\n",(int)sysPageSz));

  UT_ASSERT0( (size_t)sysPageSz == _pageSz );

  size_t maxSz = _elemSz * nMaxElem;
  UT_ASSERT0(maxSz>0);

  _nReservedPages = ( maxSz + (size_t)_pageSz - 1 ) / (size_t)_pageSz;

  UT_ASSERT0(_pageSz*_nReservedPages>=maxSz);

  _nReservedPages += 2*_nSentinelPages; // Add buffer page to safely access beyond array
  
  size_t sizeNeeded = _nReservedPages*_pageSz;

  TRC(("%s: Reserving %.0f bytes of virtual address space\n",
	UT_FUNCNAME,
	(double)sizeNeeded));

  _baseAllocated = VirtualAlloc( NULL, sizeNeeded,
      			MEM_RESERVE 
			//| MEM_TOP_DOWN			
			,
			PAGE_NOACCESS );

  TRC(("%s: base address = %p\n",UT_FUNCNAME,_baseAllocated));

  if(!_baseAllocated)
  {
    TRCERR(("VirtualAlloc(%.0f) failed. errno=%d\n",
      (double)_nReservedPages*_pageSz,
      (int)GetLastError()));
  }

  _base = BYTE_OFFSET( _baseAllocated, _nSentinelPages*_pageSz );

  UT_ASSERT0(sizeof(*_bitmap)==1<<5);
  _bitmapLen = (_nReservedPages + sizeof(uint32_t)-1 )/ sizeof(uint32_t);
  ALLOCARR_ALIGNED_(_bitmap, uint32_t, _bitmapLen );
  memset(_bitmap,0,_bitmapLen*sizeof(uint32_t));

  TIMER_STOP(&tm);
  TRC(("%s: CPU %.2f sec, %.0f MB/s)\n",UT_FUNCNAME,
      (double)TIMER_DIFF_MS(&tm)/1000.,
      (((double)_nReservedPages*_pageSz)/(double)ONE_MB) / 
	(double)(TIMER_DIFF_MS(&tm)/1000.0)));
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
SparsePagedArray::~SparsePagedArray( void ) 
{
  TRC(("%s %p\n",UT_FUNCNAME,_baseAllocated));
  UT_ASSERT0(_baseAllocated);
  BOOL ok = VirtualFree(_base, 0, MEM_RELEASE);
  if(!ok)
  {
    TRCERR(("VirtualFree(%p) failed. errno=%d\n",_base,
	  (int)GetLastError()));
  }
  FREEMEM_ALIGNED(_bitmap);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void SparsePagedArray::lock( void )
{
  LONG volatile *pLock = (LONG volatile *)&_lock;
  _nLockReq++;
  while(!((InterlockedCompareExchange( 
	    pLock, 1, 0 ) == 0)))
  {
    _nLockCol++;
    Sleep(0);
  }      	
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void SparsePagedArray::unlock( void ) 
{
  InterlockedExchange((volatile LONG*)(&_lock),0);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void *SparsePagedArray::allocElem( LongInt idxElem ) 
{
  UT_ASSERT2(idxElem>=0&&idxElem<(LongInt)_nMaxElem);
  size_t offs = idxElem*_elemSz,
	 iPage = offs>>_pageSzLog2;

  #ifdef UT_ASSERT_LEVEL_2
  {
    int idx = (iPage)>>5;
    UT_ASSERT0((idx>=0&&idx<(int)_bitmapLen));
  }
  #endif

  int doCommitPage = FALSE;

  lock();

  if( ! isPageCommited( iPage ) )
  {
    doCommitPage = TRUE;
    setPageCommited( iPage );
  }

  unlock();

  if(doCommitPage)
  {
    void *startPointer  = BYTE_OFFSET(_base,iPage*_pageSz);
    size_t sizeNeeded = _pageSz;
    LPVOID lpvResult = VirtualAlloc(
			  startPointer,		     
			  sizeNeeded,        
			  MEM_COMMIT,        
			  PAGE_READWRITE);    
    if (lpvResult == NULL )
    {
      char chbuf[1000];
      sprintf(chbuf,"SPA '%s': VirtualAlloc(%.0f,MEM_COMMIT) failed. errno=%d",
	      _name,
	      (double)sizeNeeded,
	      (int)GetLastError());
      TRCERR(("%s\n",chbuf));
      return 0;
    }
    else
    {
      _nCommitedPages++;
    }
  }

  return BYTE_OFFSET(_base,offs);
}

#endif  // MIMP_ON_WINDOWS

/*=========================================================================
 *
 *
 * 	class UT::LargeContBuf
 *
 *
 * =======================================================================*/

const int sentinelPage = 1;

size_t LargeContBuf::_totalCommitedSz=0;
  
  // See https://msdn.microsoft.com/de-de/library/windows/desktop/aa366803(v=vs.85).aspx
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
LargeContBuf::LargeContBuf( size_t initialSz, size_t maxSz, 
    			    int useLargePages ) :
  _base(NULL),
  _initialSz(initialSz),
  _actSz(0),
  _pageSz(0),
  _reservedPages(0),
  _commitedPages(0)
{
  TRC(("%s %.0f %.0f\n",UT_FUNCNAME,(double)initialSz,(double)maxSz));
  UT_ASSERT0(!useLargePages);

  UtTimer tm;
  TIMER_START(&tm);

  #ifdef MIMP_ON_WINDOWS

  SYSTEM_INFO sSysInfo;
  GetSystemInfo(&sSysInfo);
  _pageSz = sSysInfo.dwPageSize;
  
  #else

  _pageSz = sysconf(_SC_PAGESIZE);

  #endif

  TRC(("System page size = %d\n",(int)_pageSz));

  _reservedPages = ( maxSz + _pageSz - 1 ) / _pageSz;
  UT_ASSERT0(_pageSz*_reservedPages>=maxSz);
  _reservedPages += sentinelPage; // Add buffer page to safely access beyond array end
  
  size_t sizeNeeded = _reservedPages*_pageSz;

  TRC(("Reserving %.0f bytes of virtual address space\n",
	(double)sizeNeeded));

  #ifdef MIMP_ON_WINDOWS

  _base = VirtualAlloc( NULL, sizeNeeded,
      			MEM_RESERVE 
			//| MEM_TOP_DOWN			
			,
			PAGE_NOACCESS );
  if(!_base)
  {
    TRCERR(("VirtualAlloc(%.0f) failed. errno=%d\n",
      (double)_reservedPages*_pageSz,
      (int)GetLastError()));
  }

  #else

  _base = mmap(nullptr, sizeNeeded, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  if( _base == MAP_FAILED )
  {
    _base = nullptr;
  }

  if(!_base)
  {
    TRCERR(("mmap(%.0f) failed. errno=%d\n",
      (double)_reservedPages*_pageSz,
      (int)errno));
  }

  #endif

  TRC(("%s: base address = %p\n",UT_FUNCNAME,_base));

  if( _base )
  {
    if(initialSz) extend( initialSz );
  }

  TIMER_STOP(&tm);
  TRC(("%s: CPU %.2f sec, %.0f MB/s)\n",UT_FUNCNAME,
      (double)TIMER_DIFF_MS(&tm)/1000.,
      (((double)_reservedPages*_pageSz)/(double)ONE_MB) / 
	(double)(TIMER_DIFF_MS(&tm)/1000.0)));
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
LargeContBuf::~LargeContBuf( void ) 
{
  TRC(("%s %p\n",UT_FUNCNAME,_base));
  UT_ASSERT0(_base);

  int rcSys=0;

  #ifdef MIMP_ON_WINDOWS

  BOOL ok = VirtualFree(_base, 0, MEM_RELEASE);
  if(!ok)
  {
    rcSys = MI_ERROR;
    TRCERR(("VirtualFree(%p) failed. errno=%d\n",_base,
	  (int)GetLastError()));
  }

  #else

  int rc = munmap(_base, _reservedPages * _pageSz);
  if(rc)
  {
    rcSys = MI_ERROR;
    TRCERR(("munmap(%p) failed. errno=%d\n",_base,(int)errno));
  }
  
  #endif

  if(!rcSys)
  {
    UT_ASSERT0(_totalCommitedSz >= _commitedPages*_pageSz);
    _totalCommitedSz -= _commitedPages*_pageSz;
  }

}

#if 0
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
size_t LargeContBuf::getMaxCommitableExtension( void ) 
{
  size_t physmem_sz, physmem_avail_sz,commit_avail,virt_avail;
  UtGetSystemMemSize(&physmem_sz,&physmem_avail_sz,&commit_avail,&virt_avail);
  UT_ASSERT0(physmem_sz>0);
  return commit_avail;
}
#endif

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
size_t LargeContBuf::extend( size_t newSz, int noErrorTrace ) 
{
  {
    char chbuf1[256],chbuf2[256];
    TRC(("LargeContBuf::extend %s -> %s\n",
	  UtThousandSep(_actSz,chbuf1), UtThousandSep(newSz,chbuf2)));
  }

  UT_ASSERT0(_base);
  UT_ASSERT0(newSz>=_actSz);
  UtTimer tm;

  TIMER_START(&tm);

  // Trace info about currently available memory
  UtGetSystemMemSize(NULL,NULL,NULL,NULL);
  size_t commitedPagesOld = _commitedPages;

  // calculate actual number of pages needed
  size_t newPages = ( newSz + _pageSz - 1 ) / _pageSz;
  newPages += sentinelPage;  // Buffer page for safe access beyond array end
  UT_ASSERT0(newPages>=_commitedPages);

  if(newPages>_commitedPages)
  {
    UT_ASSERT0(_pageSz*(newPages-sentinelPage)>=newSz);
    if(newPages>_reservedPages) 
    {
      newPages =_reservedPages;
      newSz = (newPages-sentinelPage)*_pageSz;
    }

    size_t sizeNeeded = (newPages-_commitedPages)*_pageSz;
    void *startPointer = BYTE_OFFSET(_base, _commitedPages*_pageSz);
 
    TRC(("%s: ptr=%p size=%.0f\n",UT_FUNCNAME,
	  startPointer,(double)sizeNeeded));

    int rcSys;

    #ifdef MIMP_ON_WINDOWS

    LPVOID lpvResult = VirtualAlloc(
			  startPointer,		     
			  sizeNeeded,        
			  MEM_COMMIT,        
			  PAGE_READWRITE);   

    rcSys = (lpvResult == NULL ) ? MI_ERROR : 0;

    #else

    rcSys = mprotect(startPointer, sizeNeeded, PROT_READ | PROT_WRITE);

    rcSys = rcSys ? MI_ERROR : 0;

    #endif

    if( rcSys )
    {
      char chbuf[1000];
      sprintf(chbuf,"VirtualAlloc(%.0f,MEM_COMMIT) failed. errno=%d",
	      (double)sizeNeeded,
	      (int)GetLastError());
      if(noErrorTrace)
      {
	TRC(("%s\n",chbuf));
      }
      else
      {
	TRCERR(("%s\n",chbuf));
      }
      return 0;
    }
    else
    {
      _commitedPages = newPages;
    }
  }
  _actSz = newSz;
  UT_ASSERT0(_actSz<=_commitedPages*_pageSz && 
             _actSz<=_reservedPages*_pageSz);
  _totalCommitedSz += (newPages-commitedPagesOld)*_pageSz;

  TIMER_STOP(&tm);

  if(TIMER_DIFF_MS(&tm)>20)
  {
    TRC(("%s: CPU %.2f sec, %.0f MB/s)\n",UT_FUNCNAME,
	(double)TIMER_DIFF_MS(&tm)/1000.,
	(((double)newPages*_pageSz)/(double)ONE_MB) / 
	  (double)(TIMER_DIFF_MS(&tm)/1000.0)));
  }

  return _commitedPages*_pageSz;
}

/*=========================================================================
 *
 *
 * 	class UT::StrTok
 *
 *
 * =======================================================================*/
  
// constructor
StrTok::StrTok(const char*str, const char *sep)
{
  //TRCP(("%s: alloc(str_)\n",UT_FUNCNAME));
  str_ = strOrig_ = NULL;
  if(str)
  {
    TRC3(("UT::StrTok str='%s' sep='%s'\n",str,sep));
    strOrig_ = MM_strdup(str);
    str_ = MM_strdup(str);
    UT_ASSERT_FATAL(str_);
    for(char *sp,*tok=UtStrtok(str_,sep,&sp); 
	tok; tok=UtStrtok(NULL,sep,&sp))
    {
      tokens_.push_back(tok);
    }
  }
}

double StrTok::getFloat(size_t i, double deflt ) 
{
  double val=deflt;
  char *t=get(i);
  if(t) 
  {
    int rc = UtParseFloat(t,&val);	
    if(rc) 
    {
      val=-1e20;
      TRCERR(("StrTok::getFloat('%s',%d): UtParseFloat('%s') failed.\n",strOrig_,i,t));
      throw MI_ERROR;
    }
  }
  else
  {
    if(!(deflt>-1e20))
    {
      TRCERR(("StrTok::getFloat('%s',%d) = NULL\n",strOrig_,i));
      throw MI_ERROR;
    }       
  }
  return val;
}      

// Destructor
StrTok::~StrTok()
{
  //TRCP(("%s: free(str_)\n",UT_FUNCNAME));
  MM_free(str_);
  MM_free(strOrig_);
}

}  // namespace UT


size_t UtLargeContbufTotalUsedSize( void )
{
  return UT::LargeContBuf::getTotalCommitedSz();
}

MSBG_NAMESPACE_END

