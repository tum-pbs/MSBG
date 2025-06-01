#ifndef BLOCKPOOL_H
#define BLOCKPOOL_H

#include "globdef.h"
#include "thread.h"

#define BLOCKPOOL_FAST_MONOTONIC
//#define BLOCKPOOL_DETAIL_TRACE

/*=========================================================================
 *
 *  Block Pools for lightweight / low contention fixed size 
 *  memory management	
 *
 * =======================================================================*/

#define SBG_BLOCK_ALIGN	MAX( CPU_CACHELINE_SZ, \
    			     MEMORY_ALLOCATION_ALIGNMENT )

#define SBG_BLOCK_DATA_ALIGN CPU_SIMD_WIDTH

#if 0
#define BLOCKPOOL_MAX_EXTENDS 1600
#else
#define BLOCKPOOL_MAX_EXTENDS (55000/MIN(MSBG_RENDERDENS_BSX,16))
#endif

/*------------------------------------------------------------------------*/
/* 	 								  */
/*------------------------------------------------------------------------*/
class BlockPool
{
  public:
  
  // Create 
  static 
  BlockPool *create( 
	// Name tag for monitoring/tracing convenience
        const char *name,
        // Block size in bytes (32 bit)
      	int blockSizeBytes,  
	// Max. number of virtual blocks
	int nBlocksMax,  
	// Initially allocated  number of blocks
	int nBlocksInitial,  
	// Maximum size in MB to which the pool is allowed to grow (0=unlimited)
	int maxSizeMB,
	// Options (OPT_PROTECTABLE)
	unsigned options=0
	);    

  enum
  {
    OPT_PROTECTABLE = (1<<0),
    OPT_FIXED_EXTEND = (1<<1)
  };

  // Destroy 
  static void destroy(BlockPool** bp, int disclaimMem=0);

  // Extend by allocating a new large chunk of blocks
  int extend_( int iExt, void **pBlockOut, int nBlocksIn=0, bool doPushTofreelist=true,
      	       int *pIdxExtOut=NULL );

  // Allocate block from pool
  void *allocBlock(  int bid // block id for tracing only
      		);

  // Free block
  void freeBlock_(  void *block, // Block to be freed
      	      int bid=0 // Block id for tracing only
      		);

  size_t getNumBlocksAllocated(void) const
  {
    return _bp_next_free.load(std::memory_order_acquire);
  }

  size_t getMaxActExtensions(void) const
  {
    return (getNumBlocksAllocated() >> _bp_blocks_per_seg_log2) + 1;
  }
      
  size_t getBlockSize(void) const { return _blockSize; }  
  
  size_t getMemUsage(void) const { return _totalSize; }

  size_t getBlocksAvailable(void) const { return _nBlocksTotal; }

  int isOutOfMem(void) const { return _outOfTotalMem; }

  template<typename Data_T>
  static int offsetToBlockData( ) 
  { 
    return offsetof( Block<Data_T>, _data );
  }

  enum
  {
    BLOCK_EYECATCHER = 0x0123
  };

  template<typename Data_T> 
  struct Block
  {
    // reserve space for list node 
    // (e.g. Win32 SLIST_ENTRY for lockfree stack
    // (Maybe not necessary because it's superimposed 
    // over unused blocks only ))
    //
    SLIST_ENTRY zeropad1;  // for zero padding (SIMD access accross block boundaries)

    // Header valid for used blocks only (superimposed with free list node)  
    struct
    {
      uint16_t eyecatcher;
      uint16_t flags;
      uint32_t pad;
      void *parent;
    }
    uheader;

    static_assert( sizeof( uheader ) == 16, "" );
    static_assert( sizeof( uheader) == sizeof(SLIST_ENTRY), "");

    // 
    // Voxel Data
    //
    uint8_t zeropad2[8];  // safe SIMD stencil access of one beyond block border
    Data_T _data[0] __attribute__((aligned(SBG_BLOCK_DATA_ALIGN)));;
  };

  template<typename Data_T>
  static
  void checkEyecatcher( Block<Data_T> *block )
  {
    #ifdef UT_ASSERT_LEVEL_2
    if( block->uheader.eyecatcher != BLOCK_EYECATCHER )
    {
      char chbuf1[1000],chbuf2[1000],*p;
      int len;

      p= (char*)&block->zeropad1; len = sizeof(block->zeropad1);
      for(int i=0;i<len;i++) chbuf1[i] = isprint(p[i]) ? p[i] : '.'; 
      chbuf1[len] = 0;

      p= (char*)&block->uheader; len = sizeof(block->uheader);
      for(int i=0;i<len;i++) chbuf2[i] = isprint(p[i]) ? p[i] : '.'; 
      chbuf2[len] = 0;

      TRCERR(("%s: invalid block eyecatcher: '%s' '%s' \n",UT_FUNCNAME,chbuf1,chbuf2));
      UT_ASSERT0(FALSE);
    }
    #endif
  }

  private:

  template<typename Data_T>
  void invalidateEyecatcher( Block<Data_T> *block )
  {
    memcpy( &block->uheader, _name, MIN(sizeof(block->uheader),sizeof(_name)) );
    memcpy( &block->zeropad1, _name, MIN(sizeof(block->zeropad1),sizeof(_name)) );
  }

  char	    _name[80];
   
  size_t  _blockAlign,
          _blockSize,
	  _totalSize,
	  _maxTotalSize;

  int	  _nBlocksTotal,
	  _nBlocksMax,
	  _nBlocksExtend,
	  _outOfTotalMem,
	  _isProtectable;

  static const int MAX_EXTENDS=BLOCKPOOL_MAX_EXTENDS;
  
  #ifdef BLOCKPOOL_FAST_MONOTONIC
  std::atomic<int32_t> _bp_next_free;
  std::atomic<void*> _extends[MAX_EXTENDS];
  size_t _bp_blocks_per_seg_log2,
         _bp_blocks_per_seg,
	 _bp_max_blocks;
  #else
  void	*_extends[MAX_EXTENDS];
  #endif

  typedef struct
  {
    int extendNumBlocks;
    size_t extendSize;
  }
  ExtendInfo;

  std::vector<ExtendInfo> _extendsInfo;

  size_t _szChunkPad;

  SLIST_HEADER _freelist 
    	       __attribute__((aligned(SBG_BLOCK_ALIGN)));

  UtMutex  _lockExtend;
  int	   _lockExtendTaCnt;
};

#endif // BLOCKPOOL_H

