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
#include <atomic>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <vectorclass/vectorclass.h>
#include "globdef.h"
#include "mtool.h"
#include "fastmath.h"
#include "util.h"
#include "thread.h"
#include "blockpool.h"

#if 0
// 
// Must declare InterlockedPushListSList manually here
// See: http://stackoverflow.com/questions/11980721/interlockedpushlistslist-is-missing-where-is-it
//
extern "C" 
{
  PSLIST_ENTRY  FASTCALL InterlockedPushListSList(
      PSLIST_HEADER ListHead,
      PSLIST_ENTRY List,
      PSLIST_ENTRY ListEnd,
	 ULONG Count
  );
}
#endif

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BlockPool::extend_( int iExtend, void **pBlockOut, int nBlocksIn, 
    		       bool doPushToFreelist,
		       int *pIdxExtOut )  
{
  int rcThis=0;
  void *chunkOfBlocks=NULL,
       *chunkOfBlocks0=NULL;
  size_t szNeeded=0,
	 szChunkPad=0;
  if(pBlockOut) *pBlockOut=NULL;
  char chbuf[256];

  #ifdef BLOCKPOOL_FAST_MONOTONIC
  doPushToFreelist = false;
  #else
  UT_ASSERT0(FALSE);
  #endif

  int nBlocks = nBlocksIn ? nBlocksIn : _nBlocksExtend,
      trcLevel = 3; //pBlockOut ? 2 : 3;
      //trcLevel = 2; //pBlockOut ? 2 : 3;

  nBlocks = MIN( nBlocks, _nBlocksMax-_nBlocksTotal );
  UT_ASSERT(nBlocks>0);

  if(_maxTotalSize>0)
  {
    int nBlocksQuota = ceil( _maxTotalSize / (double)_blockSize );

    nBlocks = MIN( nBlocks, nBlocksQuota-_nBlocksTotal );
    if(nBlocks<=0)
    {
      _outOfTotalMem = TRUE;
      TRCERRR(("BlockPool '%s' reached maximum allowed size of %.0f MB\n",
	    _name,(double)_maxTotalSize/ONE_MB),MI_ENOMEM);
    }
  }

  // Pad for safe SIMD loads beyond & before block boundaries
  szChunkPad = MAX(_blockAlign,(CPU_SIMD_WIDTH));
  szChunkPad = ALIGN(szChunkPad,_blockAlign);

  szNeeded = _blockSize*nBlocks + 2*szChunkPad; 

  _szChunkPad = szChunkPad;

  UT_ASSERT(_blockSize==ALIGN(_blockSize,_blockAlign));
  UT_ASSERT(_blockAlign>=CPU_CACHELINE_SZ);

  UT_ASSERT0( iExtend < ARRAY_LENGTH(_extends));

#if 1
  TRCL(trcLevel,("SBG::BlockPool::extend() '%s' size=%.0f+%.0f MB\n",
	_name,
	(double)szNeeded/(double)ONE_MB,
	(double)_totalSize/(double)ONE_MB ));
#else
  TRC(("SBG::BlockPool::extend() '%s' size=%.0f+%.0f MB (%p) blocks=%d/%d\n",
	_name,
	(double)szNeeded/(double)ONE_MB,
	(double)_totalSize/(double)ONE_MB,
	memResource, 
	getNumBlocksAllocated(),_nBlocksMax
	));
#endif

  // 
  // Allocate new chunk of blocks
  //
  {
    char chbuf2[256], *puid=MM_Str2UserId(chbuf,chbuf2);
    if(!chunkOfBlocks0)
    {
      //if(!UtGlobDbgCheckCount==4711)  // XXX
      {
        ALLOCMEM_LARGE_ALIGNED_(chunkOfBlocks0, void, szNeeded, _blockAlign, 1, puid );
      }
    }		       
  }
  if(!chunkOfBlocks0)
  {
    _outOfTotalMem = TRUE;
    TRCERRR(("SBG::BlockPool '%s' out of total memory\n",_name),MI_ENOMEM);
  }
  chunkOfBlocks = BYTE_OFFSET(chunkOfBlocks0,szChunkPad);
  UT_ASSERT((uint64_t)chunkOfBlocks==ALIGN(((uint64_t)((char*)chunkOfBlocks)),
				            ((uint64_t)_blockAlign)));
  memset(chunkOfBlocks0,0,szChunkPad);
  memset(BYTE_OFFSET(chunkOfBlocks0,szNeeded-szChunkPad),0,szChunkPad);

  UT_ASSERT(_extends[iExtend]==NULL);
  
  _extends[iExtend] = chunkOfBlocks0;

  ExtendInfo ei;
  ei.extendSize = szNeeded;
  ei.extendNumBlocks = nBlocks;

  _extendsInfo.resize(iExtend+1);
  _extendsInfo[iExtend] = ei;

  if(pIdxExtOut) *pIdxExtOut = iExtend;

  _totalSize += szNeeded;
  _nBlocksTotal += nBlocks;
  UT_ASSERT(_nBlocksTotal<=_nBlocksMax);

  #if 0
  //
  // Push new blocks to the pool's freelist
  //
  if(doPushToFreelist)
  {
    // assemble local temp. list first
    SLIST_HEADER tmpList
		 __attribute__((aligned(SBG_BLOCK_ALIGN)));
    
    SLIST_ENTRY *firstEntry=NULL,
		*lastEntry=NULL;

    UT_ASSERT(2*sizeof(SLIST_ENTRY)<_blockSize);

    int iBlock0 = 0;

    // Reserve one block for the calling thread itself
    if(pBlockOut)
    {
      *pBlockOut = BYTE_OFFSET(chunkOfBlocks, _blockSize*iBlock0);
      iBlock0++;
      nBlocks--;
    }
    
    // Build local temp. list of new blocks
    for(int i=0;i<nBlocks;i++) 
    {
      // push in reverse order (to be later pop'ed in original order)
      SLIST_ENTRY *entry = 
	(SLIST_ENTRY*)BYTE_OFFSET(chunkOfBlocks, 
				  _blockSize*(nBlocks-1-i+iBlock0) 
				  + sizeof(SLIST_ENTRY)  );

      InterlockedPushEntrySList( &tmpList, entry );

      if(i==0) lastEntry = entry;
      if(i==nBlocks-1) firstEntry = entry;
      TRC3((" push block %p to temp.\n",i,entry)); 
    }
    TRC3(("tmp_list: first=%p last=%p\n",firstEntry,lastEntry));

    // Insert into shared freelist at once
    if(nBlocks)
    {
      TRCL(trcLevel,(" '%s': Push temp list of %d blocks to freelist\n",
	    _name, nBlocks));
      InterlockedPushListSList( &_freelist, 
 	                         firstEntry, lastEntry, nBlocks );
    }

#if 0
    for(int i=0;i<nBlocks;i++)
    {
      void *block=alloc(i);
      TRC(("alloc -> %p\n",block));
    }
#endif

    TRCL(trcLevel,("  CPU ( BlockPool::extend '%s' )\n",_name ));
  }
  #endif

  // 
  // Apply extension size policy 
  //
#if 0
  _nBlocksExtend = MIN( 1.414*_nBlocksExtend,
	  	     (int)((512*(size_t)ONE_MB)/(size_t)_blockSize));
#endif

rcCatch:
  if(rcThis)
  {
    FREEMEM_LARGE_ALIGNED(chunkOfBlocks0);
  }
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void *BlockPool::allocBlock( int bid )
{
  void *block = NULL;

  #ifdef BLOCKPOOL_FAST_MONOTONIC
  
  {
    int iBlock = _bp_next_free++;
    UT_ASSERT0(iBlock<_nBlocksMax);
    uint32_t iSeg = iBlock >> _bp_blocks_per_seg_log2,
	     iBlockInSeg = iBlock & (_bp_blocks_per_seg-1);
    UT_ASSERT2( iSeg>=0 );

    if( ! (iSeg < ARRAY_LENGTH(_extends)))
    {
      TRCERR(("BlockPool '%s': Number of extends %d exceeds maximum %d.\n",
	  _name,iSeg+1,ARRAY_LENGTH(_extends)));
      UT_ASSERT_FATAL(FALSE);  
    }

    UT_ASSERT2( iBlockInSeg>=0 && iBlockInSeg < _bp_blocks_per_seg );
    std::atomic<void*> seg = _extends[iSeg].load(std::memory_order_acquire);
    if(!seg)
    {
      int tid = UtGetCurrentThreadId(),
	  taCount;
      USE(tid);

      UtMutexLock(&_lockExtend);
      taCount = _lockExtendTaCnt++;

      // re-check if not already extended by other threads
      seg = _extends[iSeg].load(std::memory_order_acquire);
      if(!seg) 
      {
	extend_(iSeg,NULL,0,false,NULL);
	seg = _extends[iSeg].load(std::memory_order_acquire); 
	_extends[iSeg].store(seg,std::memory_order_release);	
      }

      UtMutexUnlock(&_lockExtend);
    }

    UT_ASSERT0( seg != nullptr );
    size_t offs = _blockSize * iBlockInSeg + _szChunkPad;
    UT_ASSERT2( offs + _blockSize - 1 <  _blockSize*_nBlocksExtend + _szChunkPad );    
    block = BYTE_OFFSET( seg, offs );
  }

  #else

  UT_ASSERT0(FALSE);

  #if 0
  block = (void*)InterlockedPopEntrySList( &_freelist );
  if( block ) block = BYTE_OFFSET_NEG( block, sizeof(SLIST_ENTRY));

  if(!block)
  {
    int tid = GetCurrentThreadId(),
	taCount;

    TRC3(("BlockPool '%s': %p Thread %d entering "
	  " critical section (lockExtend)\n",
	  _name,this,tid));

    UtMutexLock(&_lockExtend);
    taCount = _lockExtendTaCnt++;

    // re-check if not already extended by other threads
    block = (void*)InterlockedPopEntrySList( &_freelist );
    if( block ) block = BYTE_OFFSET_NEG( block, sizeof(SLIST_ENTRY));
    if(!block && ! _outOfTotalMem)
    {
      TRC3(("BlockPool '%s': %p Thread %d performing free list extension (cnt=%d)\n",
	    _name,this,tid,taCount));

      extend( &block );
    }

    UtMutexUnlock(&_lockExtend);

    TRC3(("BlockPool '%s': %p Thread %d leaved "
	  " critical section (lockExtend, count=%d)\n",
	  _name,this,tid,taCount));

    /*if(!block)
    {
      block = (void*)InterlockedPopEntrySList( &_freelist );
    }*/
    TRC3(("BlockPool '%s': %p Thread %d alloc block %d -> %p\n",
	  _name,this,tid,bid,block));
  }
  else
  {
    TRC3(("BlockPool '%s': %p alloc block %d -> %p\n",
	  _name,this,bid,block));
  }
  #endif

  #endif

  if(!block)
  {
    TRCERR(("out-of-memory for block pool '%s'\n",_name));
    UT_ASSERT_FATAL(FALSE);  // TODO: handle gracefully in upper layers
  }

  auto b = (Block<uint8_t> *)block;
  b->uheader.flags = 0;
  b->uheader.eyecatcher = BLOCK_EYECATCHER;
  b->uheader.parent = NULL;

  return block;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void BlockPool::freeBlock_( void *block, int bid )
{
  UT_ASSERT0(FALSE); // deprecated
  #if 0
  TRC3(("BlockPool '%s', %p free block %d -> %p\n", _name,this,bid,block));
  UT_ASSERT_NR((uint64_t)block==ALIGN((uint64_t)block,CPU_CACHELINE_SZ));
  if(block) 
  {
    invalidateEyecatcher((Block<uint8_t>*)block);
    block = BYTE_OFFSET(block,sizeof(SLIST_ENTRY));
  }
  else
  {
    UT_ASSERT0(FALSE);
  }
  InterlockedPushEntrySList( &_freelist, (SLIST_ENTRY*)block );
  #endif
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BlockPool *BlockPool::create( const char *name,
    			      int blockSize,
    			      int nBlocksMax,
    			      int nBlocksInitial,
			      int maxSizeMB,
			      unsigned options
			       )  
{
  int rcThis=0;
  BlockPool *p=NULL;
  char chbuf[1024];

  #if 0
  TRCP(("sizeof(SLIST_ENTRY) = %d\n",sizeof(SLIST_ENTRY)));
  TRCP(("sizeof(SLIST_HEADER) = %d\n",sizeof(SLIST_HEADER)));
  TRCP(("MEMORY_ALLOCATION_ALIGNMENT=%d\n",MEMORY_ALLOCATION_ALIGNMENT));
  #endif

  #ifndef BLOCKPOOL_FAST_MONOTONIC
  UT_ASSERT0(FALSE); // Impl. outdated
  #endif

  bool isProtectable = options & OPT_PROTECTABLE;

  UT_ASSERT0( offsetof(Block<float>,uheader)==16 );
  UT_ASSERT0( ALIGN(offsetof(Block<float>,_data),CPU_SIMD_WIDTH) == 
      	      offsetof(Block<float>,_data) );

  UT_ASSERT0(strlen(name)<ARRAY_LENGTH(_name)-1);

  size_t blockAlign = SBG_BLOCK_ALIGN;

  if(isProtectable)
  {
    blockAlign = ALIGN( blockAlign, MM_SYSTEM_PAGE_SIZE );    
  }

  UT_ASSERT(CPU_CACHELINE_SZ>=MEMORY_ALLOCATION_ALIGNMENT);
  UT_ASSERT(blockAlign>=16);


  // This is not intended for very small blocks 
  UT_ASSERT(blockSize>=CPU_CACHELINE_SZ);  
  // Align 
  blockSize = ALIGN(blockSize,blockAlign);

  // Allocate memory for object adm (zero initialized)
  {
    sprintf(chbuf,"BPC#%s",name);
    char chbuf2[256];
    ALLOCMEM_ALIGNED2_(p, BlockPool, sizeof(*p), CPU_CACHELINE_SZ,1,
      		     MM_Str2UserId(chbuf,chbuf2));
  }
  if(!p) 
  {
    TRCERRR(("BlockPool::create() '%s' out of memory\n",name),MI_ENOMEM);
  }

  #ifdef BLOCKPOOL_FAST_MONOTONIC
  p = new ((void*)p) BlockPool();  // placement new + value initialization

  #ifdef UT_ASSERT_LEVEL_2
  for(size_t i=0;i<sizeof(*p);i++)
  {
    UT_ASSERT0(((uint8_t*)p)[i]==0);
  }
  #endif

  #else
  memset(p,0,sizeof(*p));
  p = new ((void*)p) BlockPool();  // placement new
  #endif


  UT_ASSERT_FATAL(p);   // must be able to call destructor now 

  // Initialize members
  strcpy(p->_name,name);
  
  p->_isProtectable = isProtectable;
  p->_blockAlign = blockAlign;
  p->_blockSize = blockSize;
  p->_nBlocksMax = nBlocksMax;
  p->_maxTotalSize = maxSizeMB*(size_t)ONE_MB;
  
  #ifndef BLOCKPOOL_FAST_MONOTONIC
  InitializeSListHead( &p->_freelist );
  #else
  p->_freelist = {};
  #endif
  
  UT_ASSERT((uint64_t)&p->_freelist == 
      ALIGN((uint64_t)&p->_freelist, CPU_CACHELINE_SZ));

  UtMutexInit(&p->_lockExtend);

  /*TRC(("SBG::BlockPool::create(): Initially allocate %d blocks @ %.0f bytes\n",
	nBlocksInitial, (double)p->_blockSize));*/


  #ifndef BLOCKPOOL_FAST_MONOTONIC
  p->_nBlocksExtend = nBlocksInitial;
  if(nBlocksInitial)
  {
    p->extend( NULL );
  }
  #endif

  UT_ASSERT0( UT_IMPLIES( options & OPT_FIXED_EXTEND, nBlocksInitial==nBlocksMax ));

  if( options & OPT_FIXED_EXTEND )
  {
    p->_nBlocksExtend = nBlocksInitial;
  }
  else
  {
    // Use constant byte size for all block pool extends to prevent global heap fragmentation
    size_t extendSize_MB = 64;    
    p->_nBlocksExtend = ((extendSize_MB*(size_t)ONE_MB)/(size_t)p->_blockSize);
  }

  if(!(options & OPT_FIXED_EXTEND))
  { // avoid too much overhead for coarse grids (small blocks)
    p->_nBlocksExtend = MIN( p->_nBlocksExtend, 16384 );
  }
  
  #ifdef BLOCKPOOL_FAST_MONOTONIC 
  {
    // Round to (upper) power of 2
    uint32_t n = p->_nBlocksExtend,
	     lg2,n2;

    UT_NEXT_POW2(n,n2,lg2,1<<30);
    UT_ASSERT0( n2>=n && ((((uint32_t)(1))<<lg2) == n2));
    p->_nBlocksExtend = n2;
    p->_bp_blocks_per_seg_log2 = lg2;
    p->_bp_blocks_per_seg = 1<<lg2;
    p->_bp_max_blocks = p->_nBlocksMax;
    UT_ASSERT0((int)p->_bp_blocks_per_seg == p->_nBlocksExtend);
  }
  #endif

  TRC3(("%s: _nBlocksExtend=%d\n",p->_name, p->_nBlocksExtend));

rcCatch:
  if(rcThis)
  {
    destroy(&p);
  }
  return p;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void BlockPool::destroy( BlockPool** p_, int disclaimMem )
{
  BlockPool *p=p_?*p_:NULL;
  if(p)
  {    
   
    int nExtends = p->getMaxActExtensions();
    UT_ASSERT0(nExtends>=0&& nExtends<=ARRAY_LENGTH(p->_extends));

    #ifdef BLOCKPOOL_DETAIL_TRACE   // Trace effective mem usage
    {
      int nBlocks = p->getNumBlocksAllocated();
      if(nExtends>0 && nBlocks>16)
      {
	size_t szTotalAllocated = nExtends * p->_extendsInfo[0].extendSize,
	       szUsed = nBlocks * p->_blockSize;
        TRC1(("BlockPool::destroy '%s' eff. mem usage was %.2f%% (waste=%.2f MB)\n",
	      p->_name, 
	      (szUsed/(double)szTotalAllocated)*100.,
	      (szTotalAllocated-szUsed)/((double)ONE_MB)));
      }
    }
    #endif

    #ifdef UT_ASSERT_LEVEL_2
    nExtends = ARRAY_LENGTH(p->_extends);
    #endif

    for(int i=0;i<nExtends;i++)
    {
      UT_ASSERT2( UT_IMPLIES( i>=nExtends, p->_extends[i]==NULL ) );
      if(!p->_extends[i]) continue;

      UT_ASSERT0( (size_t)i < p->_extendsInfo.size() );

      #ifdef UT_ASSERT_LEVEL_2
      // Invalidate all eyecatchers
      void *chunkOfBlocks = BYTE_OFFSET(p->_extends[i],p->_szChunkPad);	
      int nBlocks = p->_extendsInfo[i].extendNumBlocks;
      for(int i=0;i<nBlocks;i++)
      {
	Block<uint8_t> *block = 
	  (Block<uint8_t>*)BYTE_OFFSET(chunkOfBlocks,i*p->_blockSize);
	p->invalidateEyecatcher(block);
      }
      #endif

      if(disclaimMem && p->_extends[i])
      {
	int doReleasePhysical = disclaimMem & 2,
	    isAligned = FALSE;
	UtMemDisclaim(p->_extends[i],p->_extendsInfo[i].extendSize,
	    	      doReleasePhysical,isAligned);
      }

      {
	void *ptr = p->_extends[i];
        FREEMEM_LARGE_ALIGNED(ptr);
	p->_extends[i] = NULL;
      }
    }
    UtMutexDelete(&p->_lockExtend);
    p->~BlockPool();
    FREEMEM_ALIGNED(p);
  }
}

