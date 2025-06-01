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
#include <errno.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <tbb/tbb.h>
#include <vectorclass/vectorclass.h>
#include "globdef.h"
#include "mtool.h"
#include "fastmath.h"
#include "util.h"
#include "blockpool.h"
#include "grid.h"
#include "sbg.h"


#define CLIP_BOX(x1,y1,z1, x2,y2,z2, /* IN: region */\
    		      x,y,z  /* OUT: coordantes to clip */\
    		     )\
{\
  if(unlikely((x)>(x2))) x = x2; \
  if(unlikely((x)<(x1))) x = x1; \
  if(unlikely((y)>(y2))) y = y2; \
  if(unlikely((y)<(y1))) y = y1; \
  if(unlikely((z)>(z2))) z = z2; \
  if(unlikely((z)<(z1))) z = z1; \
}


namespace SBG 
{
/*=========================================================================
 *
 * 	
 * 	S B G   
 *
 * 	Sparse Block Grids
 *
 *
 * =======================================================================*/

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template<typename Data_T>
SparseGrid<Data_T> *SparseGrid<Data_T>::create( 
      const char *name,
      int sx, int sy, int sz, 
      int blockSize, 
      int maxSizeMB,
      int initialSizeMB,
      unsigned options,
      size_t entrySize
      )  
{
  int rcThis=0;
  SparseGrid *sg=NULL;
  int simdAlign = 0;
  char chbuf[1024];

  UT_ASSERT0( 
      UT_IMPLIES( !(options & OPT_DENSE_GRID),
      		  (options & (OPT_NO_BORDPAD)) 
		  /*&& (options & (OPT_BLK_ALIGNED))*/));  

  if(options & OPT_DENSE_GRID )
  {
    simdAlign = 8;
/*    simdAlign = MAX(4,CPU_SIMD_WIDTH/(MAX(sizeof(Data_T),sizeof(float))));
    UT_ASSERT0(simdAlign==8||simdAlign==4);*/
  }

  // we assume that bit shift is handled 'arithemtically' for signed
  // integers, which should be true for most compilers 
  UT_ASSERT0((-2)>>1==-1);

  UT_ASSERT0(SDF_INFINITY<1e7);  // must be representable with 24 bit floats

  if(options & OPT_DENSE_GRID)
  {
    //UT_ASSERT(options & OPT_NO_BORDPAD);
    blockSize = 1;
  }
  else
  {
    UT_ASSERT(blockSize >= 2);  
  }

  // Check input
  UT_ASSERT(UT_ISPOW2(blockSize));
  UT_ASSERT(sx>0&&sy>0&&sz>0);
  // restrict min. blocks size to facilitate 4-point stencils
  
  // Allocate memory for object adm (zero initialized)
  ALLOCOBJ0_TAG_(sg,SparseGrid,MM_Str2UserId(name,chbuf));
  if(!sg) 
  {
    TRCERRR(("SparseGrid::create() '%s' out of memory\n",name),MI_ENOMEM);
  }

    sg = new ((void*)sg) SparseGrid();  // placement new

  UT_ASSERT_FATAL(sg);   // must be able to call destructor now 

  // Initialize
  UT_ASSERT(strlen(name)<ARRAY_LENGTH(_name)-1);
  strcpy(sg->_name,name);

  UtMutexInit(&sg->_lock);

  sg->_optionsOrig = options;
  
  sg->_isDenseGrid = options & OPT_DENSE_GRID;

  sg->_doProtectConstBlocks = options & OPT_PROTECT_CONST_BLOCKS;

  sg->_maxBordOff = (options & OPT_DENSE_GRID) ? 2 :
    		    (options & OPT_NO_BORDPAD) ? 0 : 
		    3;

  sg->_bsx = blockSize;
  sg->_bsx2 = sg->_bsx*sg->_bsx;

  sg->_sx = sg->_gridRes.x = sx;
  sg->_sy = sg->_gridRes.y = sy;
  sg->_sz = sg->_gridRes.z = sz;

  sg->_sxyzMax = MAX(MAX(sx,sy),sz);
  sg->_sxyzMin = MIN(MIN(sx,sy),sz);

  sg->_nbDist = -1.f;

  sg->_v4iDomMax = Vec4i(sx-1,sy-1,sz-1,0x7fffffffL);
  sg->_v4iDomMin = Vec4i(0,0,0,INT32_MIN);

  if( !(options & OPT_NO_BORDPAD) && !(options&OPT_DENSE_GRID) )
  {
    UT_ASSERT_FATAL(sg->_maxBordOff==3); // XXX file layout VCL_CompressGrid relies on this
  }

  {
    int blkoff = options & OPT_BLK_ALIGNED ? 0 : blockSize;
    sg->_nx = sg->_gridResBlk[0] = 
      ((sx+2*sg->_maxBordOff+1)+blkoff-1)/blockSize;
    sg->_ny = sg->_gridResBlk[1] = 
      ((sy+2*sg->_maxBordOff+1)+blkoff-1)/blockSize;
    sg->_nz = sg->_gridResBlk[2] = 
      ((sz+2*sg->_maxBordOff+1)+blkoff-1)/blockSize;
  }

  sg->_nxy = (LongInt)sg->_nx*(LongInt)sg->_ny;

  sg->_v4iBlkCoordsMax = Vec4i(sg->_nx-1,sg->_ny-1,sg->_nz-1,0x7fffffffL);

  // Calculate actual (physical) grid dimensions
  if(options & OPT_DENSE_GRID)
  {
    sg->_sxAct = 
      simdAlign * (((sx+2*sg->_maxBordOff)+simdAlign-1) / simdAlign);
    sg->_syAct = sg->_sy + 2*sg->_maxBordOff;
    sg->_szAct = sg->_sz + 2*sg->_maxBordOff;
  }
  else
  {
    sg->_sxAct = sg->_nx * sg->_bsx;
    sg->_syAct = sg->_ny * sg->_bsx;
    sg->_szAct = sg->_nz * sg->_bsx;
  }

  sg->_nhx = sg->_nx+2;
  sg->_nhy = sg->_ny+2;
  sg->_nhz = sg->_nz+2;
  sg->_nhxy = sg->_nhx * sg->_nhy;

  sg->_sxyAct = sg->_sxAct*(LongInt)sg->_syAct;
  
  //sg->_nTotalVoxels = (LongInt)sg->_sx*(LongInt)sg->_sy*(LongInt)sg->_sz;
  sg->_nTotActVoxels = (LongInt)sg->_sxAct*(LongInt)sg->_syAct*(LongInt)sg->_szAct;
  sg->_nTotVirtVoxels = (LongInt)sg->_sx*(LongInt)sg->_sy*(LongInt)sg->_sz;

  {
    LongInt tmp = (LongInt)sg->_nx * (LongInt)sg->_ny * (LongInt)sg->_nz;
    UT_ASSERT0(tmp<2000000000);
    sg->_nBlocks = tmp;
  }
    
  sg->_bsxLog2 = (int)(log(sg->_bsx)/log(2)+0.5);
  sg->_bsx2Log2 = 2*sg->_bsxLog2;
  UT_ASSERT(sg->_bsx==(1<<sg->_bsxLog2));
  sg->_bsxMask = sg->_bsx-1;
  
  sg->_v4iStridesLerp0246 = Vec4i( 0,  
                                    sg->_bsx,  
      				    sg->_bsx2, 
				    (sg->_bsx2 + sg->_bsx));

  sg->_v4iStridesLerp1357 = Vec4i( (0+1),  
                                    (sg->_bsx+1),  
      				    (sg->_bsx2+1), 
				    (sg->_bsx2 + sg->_bsx+1));

  sg->_v4iStrides3Lerp0246 = Vec4i( 3*0,  
                                    3*sg->_bsx,  
      				    3*sg->_bsx2, 
				    3*(sg->_bsx2 + sg->_bsx));

  sg->_v4iStrides3Lerp1357 = Vec4i( 3*(0+1),  
                                    3*(sg->_bsx+1),  
      				    3*(sg->_bsx2+1), 
				    3*(sg->_bsx2 + sg->_bsx+1)),

  sg->_v8iStridesLerp02461357 = Vec8i( 
      				    0,  		//0
                                    sg->_bsx,  	//2
      				    sg->_bsx2, 	//4
				    (sg->_bsx2 + sg->_bsx), //6
				    (0+1),  //1
                                    (sg->_bsx+1), //3  
      				    (sg->_bsx2+1), //5
				    (sg->_bsx2 + sg->_bsx+1) //7
				    ); 
  
  sg->_v8iStrides3Lerp02461357 = Vec8i( 
      				    3*0,  		//0
                                    3*sg->_bsx,  	//2
      				    3*sg->_bsx2, 	//4
				    3*(sg->_bsx2 + sg->_bsx), //6
				    3*(0+1),  //1
                                    3*(sg->_bsx+1), //3  
      				    3*(sg->_bsx2+1), //5
				    3*(sg->_bsx2 + sg->_bsx+1) //7
				    ); 

  sg->_v8iStrides3Lerp01234567 = Vec8i( 
      				    3*0,  		//0
				    3*(0+1),  //1
                                    3*sg->_bsx,  	//2
                                    3*(sg->_bsx+1), //3  
      				    3*sg->_bsx2, 	//4
      				    3*(sg->_bsx2+1), //5
				    3*(sg->_bsx2 + sg->_bsx), //6
				    3*(sg->_bsx2 + sg->_bsx+1) //7
				    ); 



  sg->_v4iGridStridesBlk= Vec4i(1,sg->_nx,sg->_nxy,0);

  {
    // float domain boundaries such that (int)x<sx
    Vec3Float min0(0.0f),
	      max0((float)sg->_sx, (float)sg->_sy, (float)sg->_sz ),
	      min,max;	  

    #ifdef PARTICLE_QUANT_POS

    float beps = 0.005f; // avoid underflow for 24 bit quantized values 
    min = min0+beps;
    max = max0-beps;

    auto align24bitQuant = [&]( float x ) -> float
    { // XXX: Keep consistent with code in PARTICLES::setPosition()/getPosition
      double scale = ((1<<(24))-1) * (1.0/(double)(MSBG_MAX_RESOLUTION));
      return ((int32_t)(floor( scale*x + 0.5 ))) / scale;      
    };

    for(int k=0;k<3;k++)
    {
      min[k] = align24bitQuant( min[k] );
      max[k] = align24bitQuant( max[k] );
    }

    #else

    float beps=1e-4;

    for(int iTry=0;;iTry++)
    {
      UT_ASSERT0(iTry<10 && beps>0.0f && beps<0.05f);
      min = min0+beps;
      max = max0-beps;
      TRC3(("beps=%.12f => %.12f %.12f %.12f - %.12f %.12f %.12f\n",
	    beps,min.x,min.y,min.z, max.x,max.y,max.z));
      if(( min.x>0 && min.y>0 && min.z>0 &&
		max.x < sg->_sx &&
		max.y < sg->_sy &&
		max.z < sg->_sz )) break;      
      beps *= 2;
    }

    #endif

    //UT_ASSERT0( beps==0.01f ); // avoid underflow for 24 bit quantized values

    // Make sure no negative coordinates as we use (int)() instead of floor()
    UT_ASSERT0( min.x>0.0f && min.y>0.0f && min.z>0.0f );   
    UT_ASSERT0( (int)max.x<sg->_sx && (int)max.y<sg->_sy && (int)max.z<sg->_sz );   

    // Save domain boundaries as SIMD vector

    sg->_v4fDomMin = Vec4f(min.x,min.y,min.z,-1e12f);
    sg->_v4fDomMax = Vec4f(max.x,max.y,max.z,1e12f);

    TRC3(("_v4fDomMin = %g,%g,%g, _v4fDomMax = %g,%g,%g\n",
	sg->_v4fDomMin[0],sg->_v4fDomMin[1],sg->_v4fDomMin[2],
	sg->_v4fDomMax[0],sg->_v4fDomMax[1],sg->_v4fDomMax[2] ));

    // Boundaries for linear interpolation stencil (cell-centered)
    sg->_v4fDomMin2 = Vec4f(0.5f+2*beps, 0.5f+2*beps, 0.5f+2*beps, -1e12f);
    sg->_v4fDomMax2 = Vec4f( sg->_sx-0.5f-2*beps, 
			     sg->_sy-0.5f-2*beps, 
			     sg->_sz-0.5f-2*beps, 1e12f);

    {
      Vec4i ip = truncate_to_int(floor(sg->_v4fDomMin2-0.5f));
      UT_ASSERT0( vget_x(ip)>=0 && vget_y(ip)>=0 && vget_z(ip)>=0 );
      ip = truncate_to_int(floor(sg->_v4fDomMax2-0.5f));
      UT_ASSERT0( vget_x(ip)<sg->_sx-1 && vget_y(ip)<sg->_sy-1 && 
	  	  vget_z(ip)<sg->_sz-1 );
    }
  }

  {
    LongInt tmp = sg->_bsx*(LongInt)sg->_bsx*(LongInt)sg->_bsx;
    UT_ASSERT(tmp<1500000000);
    sg->_nVoxelsInBlock = tmp;
  }
  sg->_nVoxelsInBlockAct = sg->_nVoxelsInBlock;
  if( options & OPT_BLOCK_HALO )
  {
    sg->_nVoxelsInBlockAct = (sg->_bsx+2)*(sg->_bsx+2)*(sg->_bsx+2); 
  }

  sg->_emptyValue = zeroValue<Data_T>();  // XXX: do not change, assume zero for now
  sg->_invalidValue = zeroValue<Data_T>();
  sg->_fullValue = fullValue<Data_T>();

  // Set default domain boundary conditions
  {
    DomainBC dbc;
    for(int i=0;i<3;i++) 
      for(int j=0;j<2;j++) dbc.bcType[i][j] = OPT_IPBC_NEUMAN;
    sg->setDomainBC( &dbc );
  }

  for(int i=0;i<3;i++) for(int j=0;j<2;j++) 
  {
    UT_ASSERT0(sg->_domainBC.bcType[i][j]==OPT_IPBC_NEUMAN);
  }

  sg->_valueScale = 1.0f;
  sg->_valueOffs = 0.0f;

  UT_ASSERT0( entrySize < 2*ONE_GB );
  sg->_entrySize = entrySize ? entrySize : (unsigned)sizeof(Data_T);
  sg->_hasCompTimeEntrySz = entrySize > 0 ? FALSE : TRUE;

  sg->_blockSizeBytes = sg->_entrySize*sg->_nVoxelsInBlockAct +
    			sizeof(Block<Data_T>);

  sg->_blockSizeBytes = ALIGN( sg->_blockSizeBytes,
			   CPU_CACHELINE_SZ );

  // Precalculate 2x2x2 nieghborhood offsets
  //	-> x
  //			
  //	0  --- 	1			
  //	|	|	4  ---	5
  //    |	|       |       |
  //	2  ---	3	|       |
  //			6 ---	7	
  //	| y			   z			    
  //    v
  //
  {
    int yStrid = options & OPT_DENSE_GRID ? sg->_sxAct : sg->_bsx,
        zStrid = options & OPT_DENSE_GRID ? sg->_sxyAct : sg->_bsx2;
    sg->_neighVoxOffs[0] = 0;
    sg->_neighVoxOffs[1] = 1;
    sg->_neighVoxOffs[2] = yStrid;
    sg->_neighVoxOffs[3] = yStrid + 1;
    sg->_neighVoxOffs[4] = zStrid;
    sg->_neighVoxOffs[5] = zStrid + 1;
    sg->_neighVoxOffs[6] = zStrid + yStrid;
    sg->_neighVoxOffs[7] = zStrid + yStrid + 1;
  }

  if(maxSizeMB<0) // default
  {
    maxSizeMB = 0;  
  }
  sg->_maxSizeMB = maxSizeMB;

  if(initialSizeMB<0)
  {
    initialSizeMB = 0;
  }
  sg->_initialSizeMB = initialSizeMB;

  TRC(("SBG::SparseGrid::create '%s' %dx%dx%d grid of %dx%dx%d blocks @" 
	" %d%sx%d%sx%d%s voxels @ %d bytes. Initial/max size = %.0f/%.0f MB, opt=%d\n",
	sg->_name,
	sg->_sx, sg->_sy, sg->_sz, 
	sg->_nx, sg->_ny, sg->_nz, 
	sg->_bsx, (options & OPT_BLOCK_HALO) ? "+2":"",
	sg->_bsx, (options & OPT_BLOCK_HALO) ? "+2":"",
	sg->_bsx, (options & OPT_BLOCK_HALO) ? "+2":"",
	(int)sg->_entrySize,
	(double)sg->_initialSizeMB, (double)sg->_maxSizeMB,
	options )); 

  if(options & OPT_DENSE_GRID)
  {
    TRC(("Dense grid padded to maxBordOffs%d, simdwidth=%d -> %dx%dx%d\n",
	  sg->_maxBordOff,simdAlign,sg->_sxAct,sg->_syAct,sg->_szAct));

  }

  {
    if(!(options & OPT_DENSE_GRID ))
    {
      if(!(options & OPT_NO_DATA))
      {

      if(options & OPT_TRANSPARENCY_BUF)
      {
	#if TLB_USE_SKIPMAP
	size_t szNeeded = sg->_nBlocks*(LongInt)sizeof(sg->_skipmap[0]);
	TRC3(("Allocating skipmap for %.0f blocks (%.0f MB)\n",
	      (double)sg->_nBlocks,(double)szNeeded/(double)ONE_MB));
	void *buf=NULL;
	ALLOCMEM_ALIGNED2_(buf, void, szNeeded, CPU_CACHELINE_SZ,1,
	    		   MM_Str2UserId(sg->_name,chbuf));
	if(!buf)
	{
	  TRCERRR(("Out of memory allocating block index\n"),MI_ENOMEM);
	}
	sg->_skipmap = (SkipMapEntry*)buf;
	#endif
      }

      {
	Block<Data_T> block;
	UT_ASSERT(sizeof(block.uheader)>=sizeof(SLIST_ENTRY));
	UT_ASSERT((uint8_t*)(&block._data[0]) - 
	    	  (uint8_t*)(&block) >= CPU_SIMD_WIDTH );
      }

      //
      // Create block pool
      //
      
      
      #if 0
      {
	Block<Data_T> B;
	TRCP(("%s: name='%s' dim=%dx%dx%d bsx=%d\n",UT_FUNCNAME,sg->_name,sg->_sx,sg->_sy,sg->_sz,sg->_bsx));
	TRCP(("%s: sizeof(SLIST_ENTRY)=%d\n",UT_FUNCNAME,sizeof(SLIST_ENTRY)));
	TRCP(("%s: sizeof(Block<%d>)=%d\n",UT_FUNCNAME,sizeof(Data_T),sizeof(Block<Data_T>)));
	TRCP(("BLOCK %p %d\n",&B, (char*)&B - (char*)&B ));
	TRCP(("zeropad1 %p %d\n",&B, (char*)&B.zeropad1 - (char*)&B ));
	TRCP(("uheader %p %d\n",&B, (char*)&B.uheader - (char*)&B ));
	TRCP(("zeropad2 %p %d\n",&B, (char*)&B.zeropad2 - (char*)&B ));
	TRCP(("data %p %d\n",&B, (char*)&B._data - (char*)&B ));
	TRCP(("_blockSizeBytes = %d\n",sg->_blockSizeBytes ));
	TRCP(("\n"));
      }
      #endif

      sg->_blockPool = BlockPool::create( 
			    sg->_name,
			    sg->_blockSizeBytes,
			    sg->_nBlocks,
			    initialSizeMB ?
			      MAX( (initialSizeMB*(double)ONE_MB)/
				    (double)sg->_blockSizeBytes, 4) : 0,
			    sg->_maxSizeMB
			     );

      if(!sg->_blockPool) raiseRc(MI_ENOMEM);
           
      }

      // Special block pool for constant value blocks      
      {
        snprintf(chbuf,ARRAY_LENGTH(chbuf),"%s_CONST",sg->_name);
	int nConstBlocks = 3+2;
	
	sg->_blockPoolConst = BlockPool::create( 
	    			chbuf,sg->_blockSizeBytes,nConstBlocks,nConstBlocks,0,
			        BlockPool::OPT_FIXED_EXTEND | 
				  ( sg->_doProtectConstBlocks ? 
				    BlockPool::OPT_PROTECTABLE : 0 )								
				);

	if(!sg->_blockPoolConst) raiseRc(MI_ENOMEM);

        sg->_blockProtectSize = ((BlockPool*)sg->_blockPoolConst)->getBlockSize();

	Block<Data_T>*block = sg->allocBlock(0,0,TRUE/*isConstBlock*/); 
	UT_ASSERT0(block->uheader.flags | SBLK_CONSTANT);
        block->uheader.flags |= SBLK_FULL;
	for(int i=0;i<sg->_nVoxelsInBlockAct;i++) 
	  *sg->getValueInBlockPtr0(block,i)=sg->_fullValue;
	sg->_fullBlock = block;
	if(sg->_doProtectConstBlocks) sg->protectBlock(block);

	block = sg->allocBlock(0,0,TRUE/*isConstBlock*/); 
	UT_ASSERT0(block->uheader.flags | SBLK_CONSTANT);
        block->uheader.flags |= SBLK_EMPTY;
	for(int i=0;i<sg->_nVoxelsInBlockAct;i++) 
	  *sg->getValueInBlockPtr0(block,i)=sg->_emptyValue;
	sg->_emptyBlock = block;
	if(sg->_doProtectConstBlocks) sg->protectBlock(block);

	block = sg->allocBlock(0,0,TRUE/*isConstBlock*/); 
	UT_ASSERT0(block->uheader.flags | SBLK_CONSTANT);
        block->uheader.flags |= SBLK_INVALID;
	for(int i=0;i<sg->_nVoxelsInBlockAct;i++) 
	  *sg->getValueInBlockPtr0(block,i)=sg->_invalidValue;
	sg->_invalidBlock = block;
	if(sg->_doProtectConstBlocks) sg->protectBlock(block);
      }

      if(!(options & OPT_NO_DATA))
      {
        sg->initializeBlockmap( );
      }
    }
  }

  if(options & OPT_PREP_DATA_ACCESS )
  {
    sg->prepareDataAccess();
  }


rcCatch:
  if(rcThis)
  {
    destroy(sg);
  }

  return sg;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
SparseGrid<float>* createFromDenseGrid( 
      const char *name,
      int sx, int sy, int sz,       
      int blockSize,
      float *data, 
      float zeroThreshold,
      int options
      )  
{
  TRC(("this method is for testing only\n"));
  int rcThis=0;
  UtTimer tm;
#define GRID_INDEX( grid_sx, grid_sxy, x, y, z ) \
  ((x)+(y)*((LongInt)(grid_sx))+(z)*(LongInt)(grid_sxy))

  LongInt nVoxelsSet=0,nVoxelsTotal=0;

  SparseGrid<float>*sg = SparseGrid<float>::create( name, sx,sy,sz, blockSize, 0, -1, options );
  if(!sg) raiseRc(MI_ERROR);

  TIMER_START(&tm);
  #pragma omp parallel
  {
    #pragma omp for collapse(3) schedule(dynamic) \
			     reduction(+:nVoxelsTotal,nVoxelsSet)
    for(int bk=0;bk<sg->nbz();bk++)
      for(int bj=0;bj<sg->nby();bj++)
	for(int bi=0;bi<sg->nbx();bi++)
	{
	  // Loop over voxels in block
	  double f_max=0;
	  for(int i_pass=0;i_pass<2;i_pass++)	    
	  {
	    for(int vk=0;vk<sg->bsx();vk++)
	      for(int vj=0;vj<sg->bsx();vj++)
		for(int vi=0;vi<sg->bsx();vi++)
	      {
		int rc;
		int ix,iy,iz;
		sg->getGridCoords( bi,bj,bk, vi,vj,vk,
				   ix,iy,iz );

		if(!sg->inRange(ix,iy,iz))  continue; // skip padding
		LongInt ii= GRID_INDEX(sx,sx*sy, ix,iy,iz); 
		float f = data[ii];

		if(i_pass==0)
		{
		  f_max = MAX(f_max,fabs(f));		
		  nVoxelsTotal++;		    
		}
		else
		{
		  rc = sg->setValue(ix,iy,iz, f);
		  UT_ASSERT0(rc==0);
		  nVoxelsSet++;
		}
	      }
	    if(i_pass==0 && f_max<zeroThreshold) 
	    {
	      break;  // disacard near-zero block
	    }
	  }
	}

  } // omp parallel
  TIMER_STOP(&tm);
  TRC(("SBG::createFromDensGrid:  CPU = %.4f seconds (%.0f voxels/sec)\n", 
	  (double)(TIMER_DIFF_MS(&tm)/1000.0), 
	  (sg->sx()*sg->sy()*sg->sz())/(double)(TIMER_DIFF_MS(&tm)/1000.0)));
  TRC(("SBG: nVoxelsSet=%.0f/%.0f, mem usage=%.0f MB\n",
	(double)nVoxelsSet,(double)nVoxelsTotal,
	(double)sg->getMemUsage()/(double)ONE_MB));
  UT_ASSERT(nVoxelsTotal== (sg->sx()*sg->sy()*sg->sz()));  
rcCatch:
  return sg;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
void SparseGrid<Data_T>::destroy( SparseGrid*& sg )
{
  if(sg)
  {    
    TRC3(("SBG::SparseGrid::destroy '%s' %p\n",sg->_name,sg));

    sg->checkConstBlocks(1);

    //if(sg->_doProtectConstBlocks)
    {
      Block<Data_T> *blocks[] = { sg->_fullBlock, sg->_emptyBlock, sg->_invalidBlock };
      for(int i=0;i<ARRAY_LENGTH(blocks);i++) if(blocks[i]) sg->unprotectBlock(blocks[i]);
    }

    BlockPool::destroy( &sg->_blockPoolConst );
    BlockPool::destroy( &sg->_blockPool );
    FREEMEM_ALIGNED(sg->_denseData);
    FREEMEM_ALIGNED(sg->_blockBitmapExtPtr);  
    FREEMEM_ALIGNED(sg->_blockmap);  
    FREEMEM_ALIGNED(sg->_skipmap);
    sg->~SparseGrid();
    UtMutexDelete(&sg->_lock);
    FREEMEM(sg);
  }
}

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
void SparseGrid<Data_T>::allocBlockmap( bool doInit )
{
  char chbuf[1024];
  if(!_blockmap)
  {
    size_t szNeeded = _nBlocks*(LongInt)sizeof(void*);
    TRC3(("Allocating index for %.0f blocks (%.0f MB)\n",
	  (double)_nBlocks,(double)szNeeded/(double)ONE_MB));
    void *buf=NULL;
    ALLOCMEM_ALIGNED2_(buf, void, szNeeded, CPU_CACHELINE_SZ,1,
		       MM_Str2UserId(_name,chbuf));
    if(!buf)
    {
      TRCERR(("Out of memory allocating block index\n"));
      UT_ASSERT0(FALSE);
    }
    _blockmap = (Block<Data_T>**)buf;
    if(doInit) initializeBlockmap();
    _dataGenCount++;
  }
}

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
void SparseGrid<Data_T>::prepareDataAccess( unsigned accType )
{  
  TRC3(("%s %d\n",UT_FUNCNAME,accType));   

  /*if(hasData())
  {
    TRCWARN(("%s: '%s' data already in use\n",UT_FUNCNAME,_name));
  }*/

  if(_isDenseGrid)
  {
    getDenseData();
  }
  else
  {
    allocBlockmap();
    if(!_blockPool->getBlocksAvailable()) _blockPool->extend_(0,NULL);
  }
}

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
void SparseGrid<Data_T>::reset( unsigned opt )
{  
  if(_optionsOrig & OPT_NO_DATA) return;

  TRC3(("SBG::%s '%s' %d\n",UT_FUNCNAME,_name,opt));
 
  if( opt & (4|8) )
  {
    prepareDataAccess( (opt & 4 ? ACC_READ:0) | (opt & 8 ? ACC_WRITE:0) );
    return;
  }

  if(!hasData()) return;

  int disclaimMem = opt & 3;
  BlockPool::destroy( &_blockPool, disclaimMem );
  
  FREEMEM_ALIGNED(_denseData);

  if(!_isDenseGrid)
  {
    checkConstBlocks(0);

    FREEMEM_ALIGNED(_blockmap);  

    _blockPool = BlockPool::create( 
			  _name,
			  _blockSizeBytes,
			  _nBlocks,
			  _initialSizeMB ?
			    MAX( (_initialSizeMB*(double)ONE_MB)/
				  (double)_blockSizeBytes, 4) : 0,
			   _maxSizeMB,
			  0
			);
  }
  else
  {
    UT_ASSERT0(_skipmap==NULL);
    UT_ASSERT0(_blockmap==NULL);
    UT_ASSERT0(_blockPool==NULL);
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
Block<Data_T>* SparseGrid<Data_T>::allocBlock( int bid, int doZeroData, 
    					       int isConstBlock )
{  
  Block<Data_T>*block = NULL;
  void *blockPool = isConstBlock ? _blockPoolConst : _blockPool;

  UT_ASSERT_FATAL(bid>=-1&&bid<_nBlocks);
  
  if( ! (((BlockPool*) blockPool)->isOutOfMem()))
  {
    block = (Block<Data_T>*)((BlockPool*) blockPool)->allocBlock(bid);

    if(block)
    {
      UT_ASSERT2( (uint64_t)&block->_data == ALIGN((uint64_t)&block->_data,
	    					   CPU_SIMD_WIDTH)); 
      UT_ASSERT2( (uint64_t)&block->_data == ALIGN((uint64_t)&block->_data,
	    					   CPU_CACHELINE_SZ)); 
      setParentPtr(block,this);
      block->uheader.flags = 0;

      // Allow for safe SIMD reads of one element past block boundariew
      memset(&block->zeropad1,0,sizeof(block->zeropad1));
      memset(&block->zeropad2,0,sizeof(block->zeropad2));

      initializeBlock(block, doZeroData);

      if(!isConstBlock)
      {
	if(bid>=0)
	{
	  #if 0
	  UT_ASSERT0(_blockmap[bid]==NULL);
	  #else
	  if(isValueBlock(bid) /*&& !isForeignBlock(bid)*/)
	  {
	    TRCERR(("name='%s' bsx=%d bid=%d %d %p\n",_name,_bsx,bid,
		  doZeroData,_blockmap[bid]));
	  }
	  #endif
	  _blockmap[bid] = block;
	}
      }
      else
      {
	UT_ASSERT0( block->uheader.eyecatcher == BlockPool::BLOCK_EYECATCHER );
	block->uheader.flags |= SBLK_CONSTANT;
      }
    }
    else
    {
      TRCERR(("SBG::allocBlock '%s' failed for block index %d\n",_name,bid));
    }
  }
  
rcCatch:
  return block;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
int SparseGrid<Data_T>::freeBlock_( int bid )
{
  int rcThis=0;
  UT_ASSERT0(bid>=0&&bid<_nBlocks);
  Block<Data_T> *block = _blockmap[bid];
  if(block && !isConstBlock(block))
  {
    TRCERR(("freeBlock_() is deprecated\n"));
    UT_ASSERT0(!(block->uheader.flags & SBLK_FIXED));
    UT_ASSERT0( getParentPtr(block) == this );
    //((BlockPool*) _blockPool)->freeBlock_(block,bid);
  }
  _blockmap[bid] = NULL;

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
void SparseGrid<Data_T>::protectBlock( Block<Data_T>* block )
{
  UT_ASSERT0(_doProtectConstBlocks);
  FLAG_SET(block->uheader.flags,SBLK_PROTECTED);
  UtMemProtect(block, _blockProtectSize, TRUE);
}

template <typename Data_T>
void SparseGrid<Data_T>::unprotectBlock( Block<Data_T>* block )
{
  if( block->uheader.flags & SBLK_PROTECTED )
  {
    UtMemProtect(block, _blockProtectSize, FALSE);
    FLAG_RESET(block->uheader.flags,SBLK_PROTECTED);
  }
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
size_t SparseGrid<Data_T>::getMemUsage( void ) const
{
return  (!_isDenseGrid) ?  
	    (_blockPool ? ((BlockPool*) _blockPool)->getMemUsage() : 0) +
	    (_blockPoolConst ? ((BlockPool*) _blockPoolConst)->getMemUsage() : 0) +
	      (_blockmap ? sizeof(*_blockmap)*_nBlocks : 0) :
	  _denseData ? 
	    sizeof(_denseData[0])*_sxAct*_syAct*(size_t)_szAct : 
	    0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
int SparseGrid<Data_T>::isOutOfMem( void ) const
{
  return  ((BlockPool*) _blockPool)->isOutOfMem();
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
double SparseGrid<Data_T>::getFillRate( void ) const
{
  LongInt nUsed=0;
  #pragma omp parallel for reduction(+:nUsed)
  for(int i=0;i<_nBlocks;i++) if(isValueBlock(i)) nUsed++;
  return (double)nUsed/(double)_nBlocks;
}


#if 0 // TODO: testing
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
int SparseGrid<Data_T>::getValueStats(double zeroThreshold, // IN 
					  double *pN,
					  double *pMin,
					  double *pMax,
					  double *pAvg,
					  double *pSig ) const
{
  using namespace std;
  LongInt nErr=0;
  MtStat stat;
  MT_STAT_INIT(&stat);
  #pragma omp parallel
  {
    MtStat stat_;
    MT_STAT_INIT(&stat_); 
    LongInt nErr_=0;

    #pragma omp for collapse(3) schedule(dynamic) 
    for(int bz=0;bz<_nz;bz++)  // Loop over blocks
    for(int by=0;by<_ny;by++)
    for(int bx=0;bx<_nx;bx++)
    {
      int bid = getBlockIndex(bx,by,bz);
      Data_T *data = getBlockDataPtr(bid);
      if(!data) continue;  // Skip empty blocks
      for(int vz=0;vz<_bsx;vz++)
      for(int vy=0;vy<_bsx;vy++)
      for(int vx=0;vx<_bsx;vx++)
      {
	int ix,iy,iz;
	getGridCoords( bx,by,bz, vx,vy,vz,
		       ix,iy,iz );
	if(!inRange(ix,iy,iz))  continue; // skip padding
	int vid = getVoxelInBlockIndex(vx,vy,vz);
	double z = data[vid];
	if(!(std::isfinite(z)))
	{
	  TRCERR(("invalid number %g\n",z));
	  nErr_++;
	}
	if(!(zeroThreshold<0) && z>zeroThreshold)	
	{
	  MT_STAT_UPDATE(&stat_,z);
	}
      }
    }
    #pragma omp critical
    {
      MT_STAT_RESULT(&stat_); 
      MT_STAT_UPDATE_STAT(&stat, &stat_);
      nErr += nErr_;
    }
  } // omp parallel

  if(nErr)
  {
    TRCERR(("Detected %.0f invalid numbers.\n",nErr));
  }

  MT_STAT_RESULT(&stat);
  if(pN) *pN = stat.n;
  if(pMin) *pMin = stat.min;
  if(pMax) *pMax = stat.max;
  if(pAvg) *pAvg = stat.avg;
  if(pSig) *pSig = stat.var;

  return nErr ? MI_ERROR : 0;
}
#endif

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
int SparseGrid<Data_T>::setBoundary(int typ)
{
  int rcThis=0;

  UT_ASSERT(FALSE); //not yet impl.
  // loop over blocks first
  #pragma omp parallel for
  for(int bk=0;bk<_nz;bk++) 
    for(int bj=0;bj<_ny;bj++)
      for(int bi=0;bi<_nx;bi++)
      {

      }

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
void SparseGrid<Data_T>::toDenseGrid(Data_T *dataOut)
{
  TRC(("%s '%s' -> %p\n",UT_FUNCNAME,_name,dataOut));
  InterlockedIterator iit(_nBlocks);
  #pragma omp parallel 
  {    
    for(int bid; (bid = iit.getNextBlock()) < iit.nBlocks(); )
    {
      int bx,by,bz,bsx=_bsx;
      getBlockCoordsById(bid, bx,by,bz);
      Data_T *data = getBlockDataPtr(bid);
      if(!data) data = &_emptyBlock->_data[0];
      for(int vid=0, z=bz*bsx; z<(bz+1)*bsx; z++)
      for(int y=by*bsx; y<(by+1)*bsx; y++)
      for(int x=bx*bsx; x<(bx+1)*bsx; x++, vid++)
      {
	dataOut[getVirtGridIndex(x,y,z)] = data[vid];
      }
    }
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
void SparseGrid<Data_T>::fromDenseGrid(Data_T *dataSrc )
{
  TRC(("%s '%s' -> %p\n",UT_FUNCNAME,_name,dataSrc));
  InterlockedIterator iit(_nBlocks);
  #pragma omp parallel 
  {    
    for(int bid; (bid = iit.getNextBlock()) < iit.nBlocks(); )
    {
      int bx,by,bz,bsx=_bsx;
      getBlockCoordsById(bid, bx,by,bz);
      Data_T *data = getBlockDataPtr(bid,TRUE);
      for(int vid=0, z=bz*bsx; z<(bz+1)*bsx; z++)
      for(int y=by*bsx; y<(by+1)*bsx; y++)
      for(int x=bx*bsx; x<(bx+1)*bsx; x++, vid++)
      {
	data[vid] = dataSrc[getVirtGridIndex(x,y,z)];
      }
    }
  }
}


/*=========================================================================
 *
 *
 *  Utility
 *
 *
 * =======================================================================*/

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template<typename Data_T>  Data_T zeroValue( void ) 
{
    return  static_cast<Data_T>(0.0);
}
template<> float zeroValue( void ) { return 0.0f; }
template<> uint16_t zeroValue( void ) { return 0; }
template<> uint8_t zeroValue( void ) { return 0; }
template<> Vec3Float zeroValue( void ) { return 0.0f; }

template<typename Data_T>  Data_T fullValue( void ) 
{
    return  static_cast<Data_T>(1.0);
}
template<> float fullValue( void ) { return 1.0f; }
template<> Vec3Float fullValue( void ) { return 1.0f; }


template<> float diffValues( const float& val1, 
    			     const float& val2 )
{
  return fabs(val1-val2);
}

template<> float diffValues( const Vec3Float& val1, 
    			     const Vec3Float& val2 ) 
{
  return normNoSqrt(val1-val2);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static 
int intersectRayCell(
		float *pos,	// IN: Ray origin
		float *dir,	// IN: Ray direction
		// IN: cell position (lower corner) and edge length
		float *cell, float cellSz,  
		// OUT: Ray entry/exit times
		float& tmin, float& tmax  
		)
{
   float tymin, tymax, tzmin, tzmax,
   	  eps = 1e-7,   
          dx_inv = fabs(dir[0])>eps ? 1/dir[0] : 1/eps,
          dy_inv = fabs(dir[1])>eps ? 1/dir[1] : 1/eps,
          dz_inv = fabs(dir[2])>eps ? 1/dir[2] : 1/eps;
    
   if (dx_inv >= 0) 
   {
      tmin = (cell[0] -pos[0]) * dx_inv;
      tmax = (cell[0]+cellSz -pos[0]) * dx_inv;
   }
   else 
   {
      tmin = (cell[0]+cellSz -pos[0]) * dx_inv;
      tmax = (cell[0] -pos[0]) * dx_inv;
   }
   
   if (dy_inv >= 0) 
   {
      tymin = (cell[1] -pos[1]) * dy_inv;
      tymax = (cell[1]+cellSz -pos[1]) * dy_inv;
   }
   else 
   {
      tymin = (cell[1]+cellSz -pos[1]) * dy_inv;
      tymax = (cell[1] -pos[1]) * dy_inv;
   }
   
   if ((tmin > tymax) || (tymin > tmax)) return MI_ENOTFOUND;
   if (tymin > tmin) tmin = tymin;
   if (tymax < tmax) tmax = tymax;
   
   if (dz_inv >= 0) 
   {
      tzmin = (cell[2] - pos[2]) * dz_inv;
      tzmax = (cell[2]+cellSz - pos[2]) * dz_inv;
   }
   else 
   {
      tzmin = (cell[2]+cellSz - pos[2]) * dz_inv;
      tzmax = (cell[2] - pos[2]) * dz_inv;
   }
   if ((tmin > tzmax) || (tzmin > tmax)) return MI_ENOTFOUND;
   if (tzmin > tmin) tmin = tzmin;
   if (tzmax < tmax) tmax = tzmax;
    
   return MI_OK;
}

/*=========================================================================
 *
 *
 *  Sparse Transparency Buffer (x-axis point towards light-dir)
 *
 *
 * =======================================================================*/

template< typename Data_T >
float SparseGrid<Data_T>::lookupTransparencyBuf( 
    				            float *pos // Grid position
				           ) 
{
  UT_ASSERT0(_blockmap);
  #ifdef TLB_USE_SKIPMAP
  UT_ASSERT0(_skipmap);
  #endif

  // Neuman BC on whole grid -> clip 
  float beps=1e-2;
  CLIP_BOX( .5f+beps, .5f+beps, .5f+beps, 
      	   _sx-.5f-beps, _sy-.5f-beps, _sz-.5f-beps,
	   pos[0],pos[1],pos[2]);

  int x=floorf(pos[0]),
      y=floorf(pos[1]),
      z=floorf(pos[2]);

  //applyBorderOffset( x, y, z );

  int bx = x>>_bsxLog2, bx0=bx,
      by = y>>_bsxLog2,
      bz = z>>_bsxLog2;

  int bid = getBlockIndex(bx,by,bz);
  UT_ASSERT2(bid>=0&&bid<_nBlocks);


  // skip to last non-constant predecessor
  float fullVal = getFullValue(),
        f_skip = 0;

  #ifdef TLB_USE_SKIPMAP
  
  int dvx = (x & (_bsx-1))+1;
  SkipMapEntry *sme = &_skipmap[bid];
  bx = bx0 + sme->nSkip;
  f_skip = dvx * fullVal + 
    	   _bsx * fullVal * sme->nFull;
  if(bx<0) return f_skip; //0.0f;

  #else
  int dvx = (x & (_bsx-1))+1;
  for( ; bx>=0 && !isValueBlock(bid); bx--, bid-=1 )
  {
    if(isFullBlock(bid))
    {
      f_skip += fullVal * dvx;      
    }
    dvx = _bsx; 
  }
  if(bx<0) return f_skip; //0.0f;
  #endif

  UT_ASSERT2(bx>=0&&bid>=0);
  UT_ASSERT2(_maxBordOff==0.0f);

  if(bx!=bx0)
  {
    // project to face of non-empty predecessor block
    pos[0] =  (bx+1)*_bsx -0.5f - beps;
    UT_ASSERT2(!(pos[0]<0));
  }

  // Interpolate with neuman-BC on empty neighbor blocks

  {
    float buf_x[2],buf_y[2],buf_z[2];

    float  x = pos[0]-.5f, // po=0.5 -> interpolate at pixel centers
	   y = pos[1]-.5f, 
	   z = pos[2]-.5f;

    int    ix = floorf(x), 
	   iy = floorf(y), 
	   iz = floorf(z);

    int ix0 = ix,
	iy0 = iy,
	iz0 = iz;
  
    //applyBorderOffset( ix0, iy0, iz0 );

    float u = x - ix,
	  v = y - iy,
	  w = z - iz;

    // Locate the block that contains the origin of the stencil window
    int bx0, by0, bz0, bid;
    int vx0, vy0, vz0;
    getBlockCoords(ix0, iy0, iz0, 
		   bx0, by0, bz0);
    bid = getBlockIndex(bx0, by0, bz0);
    getVoxelInBlockCoords(ix0, iy0, iz0, 
			  vx0, vy0, vz0 );

    // Check for fast case: whole stencil window within block
    if( ! (((vx0+1)>>_bsxLog2) |
	   ((vy0+1)>>_bsxLog2) |
	   ((vz0+1)>>_bsxLog2)))
    {
      Block<Data_T>*block = getBlock(bid);
      UT_ASSERT0(block);
      for(int k=0;k<2;k++)
      {
	for(int j=0;j<2;j++)
	{
	  Data_T *data = block->_data + 
		    getVoxelInBlockIndex(vx0+0,vy0+j,vz0+k);
	  buf_x[0] = data[0];
	  buf_x[1] = data[1];
	  buf_y[j] = Interpolate1D<float,float,IP_LINEAR>(buf_x,u);
	}
	buf_z[k] = Interpolate1D<float,float,IP_LINEAR>(buf_y,v);	       
      }
      float f_out = Interpolate1D<float,float,IP_LINEAR>(buf_z,w);      	
      f_out += f_skip;
      return f_out;
    }
    else
    {
      // Not within single block:
      // Take into account the neighboring blocks 
      // (A total maximum of 2x2x2 blocks)
      float f[8],fsum=0;
      int l=0,nNeigh=0;
      for(int k=0;k<2;k++)
      {
	int ibz = ( bz0 + ((vz0+k) >> _bsxLog2) ) *_nxy,
	    ivz = (((vz0+k) & _bsxMask) << _bsx2Log2);
	for(int j=0;j<2;j++)
	{
	  int iby = ibz + (by0 + ((vy0+j) >> _bsxLog2)) *_nx,
	      ivy = ivz + (((vy0+j) & _bsxMask) <<_bsxLog2);
	  for(int i=0;i<2;i++,l++)
	  {
	    int ib = iby + bx0 + ((vx0+i) >> _bsxLog2),
		iv = ivy + ((vx0+i) & _bsxMask);
	    UT_ASSERT2(ib>=0&&ib<_nBlocks);
	    Block<Data_T>* block = _blockmap[ib]; 	  	    
	    if(block)
	    {
	      f[l] = block->_data[ iv ];
	      fsum += f[l];
	      nNeigh++;
	    }
	    else f[l] = -1e20;
	  }
	}
      }

      UT_ASSERT2(nNeigh>0);
      if(nNeigh<8)
      {
	// Set missing neighbors to avg. of all known ones
	float favg = fsum/(float)nNeigh;
	for(int l=0;l<8;l++) if (f[l]<-1e19) f[l] = favg;
      }

       float f_out = 
	  ( (   (1-u) * f[0] 
	  +      u * f[1]  ) * (1-v)
	  + (   (1-u) * f[2] 
	  +      u * f[3]  ) * v ) * (1-w)
	  + ( ( (1-u) * f[4] 
	  +      u * f[5]  ) * (1-v)
	  + (   (1-u) * f[6] 
	  +      u * f[7]  ) * v ) * w ;

       f_out += f_skip;

      return f_out;
    }    
  }
}

/*=========================================================================
 *
 *
 *  Filtering
 *
 *
 * =======================================================================*/
template<typename Data_T>
int SparseGrid<Data_T>::filter( int typ, 
    				float krn_r,  // Kernel radius
 				float krn_s,  // Sharpening strength
				unsigned options,
	      			SparseGrid<Data_T> **pDstFieldOut )
{
  int rcThis=0;
  
  TRC(("%s %d %g %g\n",UT_FUNCNAME,typ,krn_r,krn_s));

  UT_ASSERT0(typ==1); // gaussian blur
  UT_ASSERT0(_maxBordOff==0);
  
  UtTimer tm;
  TIMER_START(&tm);

  SparseGrid<Data_T> *dstField = NULL,  
		     *srcField = this;

  dstField = srcField->clone( "FILTERTMP", FALSE /* do not copy data */ );

  //
  // Prepare kernel coefficients
  //
  double coeffs[7*7*7];

  int krnRadius = 1,
      krnWidth = krnRadius*2+1,
      doBlur = fabs(krn_s) < MT_NUM_EPS;

  USE(doBlur);

  UT_ASSERT0(krnWidth==3);

  krn_s *= 2;
  krn_r *= 0.75;
  double radius=krn_r,
	 sharpness=krn_s;
  double sigma = sharpness > 0 ? radius*2.0 *sharpness / 3. : 
				 radius*2.0;
  int i_pass,width=krnWidth;
  double sum=0;
  for(i_pass=0;i_pass<2;i_pass++)
  {
    for ( int w = -width/2,i=0; w <= width/2; ++w)
    for ( int v = -width/2; v <= width/2; ++v)
    for ( int u = -width/2; u <= width/2; ++u,++i)
    {
      double a = exp(-(u*u+v*v+w*w)/(2.*sigma*sigma));
      if(i_pass==0)
      {
	sum += a;
      }
      else
      {
	if(sharpness>0)
	{
	  a = u==0&&v==0&&w==0 ? -(a/sum)*sharpness+sharpness+1: 
				 -(a/sum)*sharpness;
	}
	else a = a/sum;
	coeffs[i] = a;
	TRC(("%7.4f %s",a,v==width/2&&u==width/2?"\n\n":u==width/2?"\n":""));
      }
    }
    UT_ASSERT_FATAL(sum>0);
  }

  {
    double sum=0;
    for(int i=0;i<krnWidth*krnWidth*krnWidth;i++) sum += coeffs[i];
    UT_ASSERT0(fabs(sum-1)<MT_NUM_EPS);
  }

  //
  // Loop over all blocks and apply filter
  //
  //
  int bsx = _bsx,
      bsxLog2 = _bsxLog2,
      bsx2Log2 = _bsx2Log2,
      sx = _sx,
      sy = _sy,
      sz = _sz,
      nx = nbx(),
      nxy = nbxy(),
      bsxMask = bsx-1;

  #pragma omp parallel for schedule( static, 10 )
  for(int bid=0; bid < _nBlocks; bid++)
  {
    int bx,by,bz;
    getBlockCoordsById(bid, bx, by, bz);

    if(srcField->isValueBlock(bid)) continue;

    Data_T *dataSrc = srcField->getBlockDataPtr(bid,0,0);
    if(!dataSrc) continue;

    Data_T *dataDst = dstField->getBlockDataPtr(bid,1,0);

    for(int vid=0,vz=0;vz<_bsx;vz++)
    for(int vy=0;vy<_bsx;vy++)
    for(int vx=0;vx<_bsx;vx++,vid++)
    {
      float val = 0.0f;
	 //   f0 = dataSrc[vid];

      int x,y,z;
      getGridCoords( bx, by, bz, vx, vy, vz,
		     x, y, z );	  

      UT_ASSERT2(krnWidth==3);

      int ix0=x-1,
	  iy0=y-1,
	  iz0=z-1;

      if( (options & OPT_ZERO_BORDER) && 
	 (x==0||y==0||z==0||
	  x==_sx-1||y==_sy-1||z==_sz-1) )
      {
        dataDst[vid] = 0.0f;   	
	continue;
      }

      for(int k=0,iKrn=0;k<3;k++)
      {
	int z2 = iz0+k;
	MT_CLAMP( z2, 0, sz-1 );
	int bz2 = z2 >> bsxLog2,
	    vz2 = z2 & bsxMask,
	    bidZ = bz2*nxy,
	    vidZ = (vz2<<bsx2Log2);

	for(int j=0;j<3;j++)
	{
	  int y2 = iy0+j;
	  MT_CLAMP( y2, 0, sy-1 );
	  int by2 = y2 >> bsxLog2,
	      vy2 = y2 & bsxMask,
	      bidY = bidZ+by2*nx,
	      vidY = vidZ + (vy2<<bsxLog2); 

	  for(int i=0;i<3;i++,iKrn++)
	  {
	    int x2 = ix0+i;
	    MT_CLAMP( x2, 0, sx-1 );
	    int bx2 = x2 >> bsxLog2,
		vx2 = x2 & bsxMask;
	    int bid2 = bidY + bx2,
		vid2 = vidY + vx2;

	    double w = coeffs[iKrn];

	    UT_ASSERT2(bid2>=0&&bid2<_nBlocks);
	    UT_ASSERT2(vid2>=0&&vid2<nVoxelsInBlock());

	    SBG::Block<Data_T>* block = getBlock(bid2); 	  

	    float f = block ? block->_data[vid2] : 0.0;

	    val += f*w;
	  } // kernel x
	} // kernel y
      }  // kernel z

      dataDst[vid] = val;   	
    }  // Voxels
  } // Blocks

  TIMER_STOP(&tm);
  TRC(("CPU(%s) =%.2f sec, %.0f virt. voxels/sec/iter\n",UT_FUNCNAME,
	(double)TIMER_DIFF_MS(&tm)/1000.,
      (double)_sx*_sy*_sz/(double)(TIMER_DIFF_MS(&tm)/1000.0)));


  // Swap data pointers to make the final result field the new current field
  if(rcThis==0)
  {
    if(!pDstFieldOut)
    {
      swapData( dstField );
    }
  }

rcCatch:

  // Cleanup intermediate and old fields
  if(pDstFieldOut) 
  {
    if( rcThis )
    {
      destroy( dstField );
      *pDstFieldOut = NULL;
    }
    else *pDstFieldOut = dstField;
  }
  else
  {
    destroy( dstField );
  }

  return rcThis;
}
template<>
int SparseGrid<Vec3Float>::filter( int typ, 
    				float krn_r,  // Kernel radius
 				float krn_s,  // Sharpening strength
				unsigned options,
	      			SparseGrid<Vec3Float> **pDstFieldOut )
{
  return MI_ENOTIMPL;
}
template<>
int SparseGrid<Vec3Uint16>::filter( int typ, 
    				float krn_r,  // Kernel radius
 				float krn_s,  // Sharpening strength
				unsigned options,
	      			SparseGrid<Vec3Uint16> **pDstFieldOut )
{
  return MI_ENOTIMPL;
}

/*=========================================================================
 *
 *
 *
 *
 * =======================================================================*/

template<> float SparseGrid<Vec3Float>::lookupTransparencyBuf( float *pos )
{
  return 0;
}
template<> float SparseGrid<Vec3Uint16>::lookupTransparencyBuf( float *pos )
{
  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template class SparseGrid<float>;
template class SparseGrid<double>;
template class SparseGrid<uint8_t>;
template class SparseGrid<uint16_t>;
//template class SparseGrid<int16_t>;
template class SparseGrid<Vec3Float>;
template class SparseGrid<Vec3Uint16>;

} // namespace SBG

