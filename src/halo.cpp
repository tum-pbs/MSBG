/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#define __STDC_LIMIT_MACROS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <vectorclass/vectorclass.h>
#include "vectorclass_util.h"
#include "globdef.h"
#include "mtool.h"
#include "util.h"
#include "plot.h"
#include "bitmap.h"
#include "panel.h"
#include "sbg.h"
#include "msbg.h"
#include "halo.h"



namespace MSBG 
{

using namespace SBG;


/*=========================================================================
 *
 *
 *  Halo Block Set
 *
 *
 * =======================================================================*/

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
HaloBlockSet::HaloBlockSet( MultiresSparseGrid *msbg,
		  	    int dmax, int nMaxBlocks,
			    unsigned options
			    ) :
  _msbg(msbg), 
  _dmax(dmax), 
  _nMaxBlocks(nMaxBlocks),
  _nVecCmp( options & OPT_HALO_VEC3 ? 3 : 1),
  _haloBlocks(NULL),
  _dataValid(NULL)
{
  SparseGrid<Vec3Float> *sg0=_msbg->sg0();
  size_t 
      bsx = sg0->bsx(),
      bsxAct = bsx+2*_dmax,
      bsyAct = bsx+2*_dmax,
      bszAct = bsx+2*_dmax,
      bufLen = bsxAct*bsyAct*bszAct;

  bufLen += 16;   // pad for SIMD beyond-border safety

  ALLOCARR0_(_haloBlocks,float*,nMaxBlocks*_nVecCmp);
  ALLOCARR0_(_dataValid,int,nMaxBlocks*_nVecCmp);

  _bufLenAllocated = bufLen;
  for(int i=0;i<_nMaxBlocks;i++)
  {
    /*if(i==0)
    {
      TRC(("%s: Allocating %.0f bytes for %dx%dx%d halo block\n",
	    UT_FUNCNAME,(double)(bufLen*sizeof(float)),
	    bsxAct,bsyAct,bszAct));
    }*/

    for(int j=0;j<_nVecCmp;j++)
    {
      float *ptr=NULL;
      ALLOCARR_ALIGNED_( ptr, float, bufLen );
      memset(ptr,0,bufLen*sizeof(*ptr));
      _haloBlocks[haloBlockIdx(i,j)] = ptr;
    }    
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
HaloBlockSet::~HaloBlockSet( void ) 
{
  if(_haloBlocks)
  {
    for(int i=0;i<_nMaxBlocks*_nVecCmp;i++)
    {
      FREEMEM_ALIGNED(_haloBlocks[i]);
    }
    FREEMEM(_haloBlocks);
  }
  FREEMEM(_dataValid);
}

/*-------------------------------------------------------------------------*/ 
/* 									   */
/*-------------------------------------------------------------------------*/
float **HaloBlockSet::fillHaloBlockFast16Uint8(  SparseGrid<uint8_t> *sg,
    		        int bid,  
    			int tid,   // halo block nr. (thread id)
			unsigned options
    			)
{
  TRCWARN(("experimental\n"));  // only ~10% improvement for lapl-smoothing 
  
  int sx = sg->sx(),
      sy = sg->sy(),
      sz = sg->sz();

  constexpr int bsxLog2 = 4,
	    	bsx2Log2 = 2*bsxLog2,
		bsx = 1<<bsxLog2,
		bsxMask = bsx-1,
		hsx = bsx+2,
		hsx2 = hsx*hsx;
  int nbx = sg->nbx(),
      nbxy = sg->nbxy();

  int hix0 = haloBlockIdx(tid,0);
  float **dataDst0 = &_haloBlocks[hix0],
	 *dataDst = dataDst0[0];

  int *dataValid0 = &_dataValid[hix0];
  dataValid0[0]=0;

  // Loop over 3*3*3 block neighborhood
  int bx0,by0,bz0;
  sg->getBlockCoordsById(bid, bx0,by0,bz0);

  Vec4i pos0 = Vec4i(bx0,by0,bz0,0)*bsx,
	pos1 = pos0-_dmax,
	pos2 = pos0+bsx+_dmax-1;

  int x0=vget_x(pos0), x1=vget_x(pos1), x2=vget_x(pos2),
                       y1=vget_y(pos1), y2=vget_y(pos2),
                       z1=vget_z(pos1), z2=vget_z(pos2);

  uint8_t *dataBlock0 = sg->getBlockDataPtr0(bid);

  for(int z=z1; z<=z2; z++)
  {
    int inRangeZ = z>=0&&z<sz,
        ibz = (z >> bsxLog2 ) * nbxy,
	ivz = (z & bsxMask) << bsx2Log2,
	ihz = (z-z1) * hsx2;

    for(int y=y1; y<=y2; y++)
    {		
      int inRangeZY = inRangeZ && y>=0 && y<sy,
	  ibzy = ibz + (y >> bsxLog2) * nbx,
	  ivzy = ivz + ((y & bsxMask) << bsxLog2),
	  ihzy = ihz + (y-y1)*hsx;		

      // 
      // First line segment
      //
      for(int x=x1; x<x0; x++)
      {
	int inRange = inRangeZY && x>=0 && x<sx,
	    ib = ibzy + ( x >> bsxLog2 ),
	    iv = ivzy + ( x & bsxMask ),
	    ih = ihzy + x-x1;
	float val; 
	if( inRange && sg->isValueBlock(ib) )
	{
	  val = sg->getBlockDataPtr0(ib)[iv];
	}
	else if( options & OPT_BC_NEUMANN )
	{
          Vec4i pos = min(max( Vec4i(x,y,z,0), pos0 ), pos0+bsx-1 );
          int vid = sg->getVoxelInBlockIndex(sg->getVoxelInBlockCoords(pos));
          UT_ASSERT2(vid>=0&&vid<sg->nVoxelsInBlock());
          val = dataBlock0[ vid ];
	}
	else val = 0; 

	RENDER_DENS_TO_FLOAT_UI8(val);
	dataDst[ih] = val;
      }

      // 
      // Middle line segment 
      //
      {
        int x=x0;
	int inRange = inRangeZY && x>=0 && x<sx,
	    ib = ibzy + ( x >> bsxLog2 ),
	    iv = ivzy,
	    ih = ihzy + x0-x1;

	if( inRange && sg->isValueBlock(ib) )
	{
	  transferInteriorLine<uint8_t>( iv, ih, bsx, hsx, sg, 
					sg->getBlockDataPtr0(ib),
					dataDst, NULL, NULL );
	}
	else if( options & OPT_BC_NEUMANN )
	{
          Vec4i pos = min(max( Vec4i(x,y,z,0), pos0 ), pos0+bsx-1 );
          int vid = sg->getVoxelInBlockIndex(sg->getVoxelInBlockCoords(pos));
	  transferInteriorLine<uint8_t>( vid, ih, bsx, hsx, sg, 
					sg->getBlockDataPtr0(bid),
					dataDst, NULL, NULL );
	}
	else
	{
	  Vec8f zero(0.0f);
	  for(int i=0;i<bsx;i+=8) zero.store(dataDst+ih+i);
	}
      }

      // 
      // Last line segment
      //
      for(int x=x0+bsx; x<=x2; x++)
      {
	int inRange = inRangeZY && x>=0 && x<sx,
	    ib = ibzy + ( x >> bsxLog2 ),
	    iv = ivzy + ( x & bsxMask ),
	    ih = ihzy + x-x1;
	float val; 
	if( inRange && sg->isValueBlock(ib) )
	{
	  val = sg->getBlockDataPtr0(ib)[iv];
	}
	else if( options & OPT_BC_NEUMANN )
	{
          Vec4i pos = min(max( Vec4i(x,y,z,0), pos0 ), pos0+bsx-1 );
          int vid = sg->getVoxelInBlockIndex(sg->getVoxelInBlockCoords(pos));
          UT_ASSERT2(vid>=0&&vid<sg->nVoxelsInBlock());
          val = dataBlock0[ vid ];
	}
	else val = 0; 

        RENDER_DENS_TO_FLOAT_UI8(val);
	dataDst[ih] = val;
      }
    }
  }

  dataValid0[0]=TRUE;

  return _haloBlocks+hix0;
}


/*-------------------------------------------------------------------------*/ 
/* 									   */
/*-------------------------------------------------------------------------*/
template< typename Data_T, int nVecCmp_T, int do1stOrderOnly_T > 
float **HaloBlockSet::fillHaloBlock_(  int iChan,
    		        int bid,   // MSBG block ID
			int levelMg,
    			int tid,   // halo block nr. (thread id)
			unsigned options
    			)
{
//  #define TEST_NO_BORDERS

  UT_ASSERT2( nVecCmp_T==1 || nVecCmp_T==3);
  UT_ASSERT2(bid>=0&&bid<_msbg->nBlocks());
  UT_ASSERT2(levelMg>=0&&levelMg<_msbg->getNumLevels());
  BlockInfo *bi=_msbg->getBlockInfo(bid,levelMg);  
  SparseGrid<Data_T> *sg = getChannel<Data_T>(_msbg,iChan,bi->level,levelMg),
                     *sgLo = bi->level<_msbg->getNumLevels()-1 ? 
		       getChannel<Data_T>(_msbg,iChan,bi->level+1,levelMg) : NULL;

  //SparseGrid<CellFlags> *sgFlagsLo = _msbg->getFlagsChannel(bi->level+1,levelMg);
  //
  if(options & OPT_INTERPOLATE)
  {
    UT_ASSERT0(_dmax<=1);  
  }

  int bsx=sg->bsx();
  int hsx=bsx+2*_dmax,
      hsx2=hsx*hsx;
  USE(hsx);

  Data_T sgEmptyVal = sg->getEmptyValue();

  int hix0 = haloBlockIdx(tid,0);
  float **dataDst0 = &_haloBlocks[hix0],
	 *dataDstX = dataDst0[0],
	 *dataDstY = nVecCmp_T==3 ? dataDst0[1] : NULL,
	 *dataDstZ = nVecCmp_T==3 ? dataDst0[2] : NULL;

  int *dataValid0 = &_dataValid[hix0];
  for(int k=0;k<nVecCmp_T;k++) dataValid0[k]=0;

  int sx = sg->sx(),
      sy = sg->sy(),
      sz = sg->sz(),    
      bsxLog2 = sg->bsxLog2(),
      bsx2Log2 = sg->bsx2Log2(),
      bsxMask = sg->bsx()-1,
      nbx = sg->nbx(),
      nbxy = sg->nbxy();

  // Loop over 3*3*3 block neighborhood
  int bx0,by0,bz0;
  sg->getBlockCoordsById(bid, bx0,by0,bz0);

  Vec4i pos0 = Vec4i(bx0,by0,bz0,0)*bsx,
	pos1 = pos0-_dmax,
	pos2 = pos0+bsx+_dmax-1;

  int x0=vget_x(pos0), x1=vget_x(pos1), x2=vget_x(pos2),
      y0=vget_y(pos0), y1=vget_y(pos1), y2=vget_y(pos2),
      z0=vget_z(pos0), z1=vget_z(pos1), z2=vget_z(pos2);
  USE(y0);
  USE(z0);
    
  Data_T *dataBlock0 = sg->getBlockDataPtr0(bid);
  UT_ASSERT2(dataBlock0);

  if( do1stOrderOnly_T )
  {
    UT_ASSERT0( _dmax==1 );
  }

  for(int z=z1; z<=z2; z++)
  {
    int inRangeZ = z>=0&&z<sz,
        ibz = (z >> bsxLog2 ) * nbxy,
	ivz = (z & bsxMask) << bsx2Log2,
	ihz = (z-z1) * hsx2;

    for(int y=y1; y<=y2; y++)
    {		
      int inRangeZY = inRangeZ && y>=0 && y<sy,
	  ibzy = ibz + (y >> bsxLog2) * nbx,
	  ivzy = ivz + ((y & bsxMask) << bsxLog2),
	  ihzy = ihz + (y-y1)*hsx;		

      bool doHorizontalHalo = true;
      if constexpr ( do1stOrderOnly_T )
      {
        doHorizontalHalo = z>=z0 && z<z0+bsx && y>=y0 && y<y0+bsx;
      }

      // 
      // First line segment
      //
      
      #ifndef TEST_NO_BORDERS      
      if( doHorizontalHalo )
      {
	for(int x=x1; x<x0; x++)
	{
	  int inRange = inRangeZY && x>=0 && x<sx,
	      ib = ibzy + ( x >> bsxLog2 ),
	      iv = ivzy + ( x & bsxMask ),
	      ih = ihzy + x-x1;

	  Data_T val;
	  bool isValid=FALSE;

	  /*if(x==-1&&y==23&&z==-1)
	  {
	    TRCERR(("DBG\n"));
	  }*/


	  if(inRange)
	  {
	    UT_ASSERT2(ib>=0&&ib<sg->nBlocks());	
	    SBG::Block<Data_T> *block = sg->getBlock(ib); 	  	    
	    UT_ASSERT2(iv>=0&&iv<sg->nVoxelsInBlock());
	    if( block )
	    {
	      val = block->_data[iv];
	      isValid = TRUE;
	    }
	  }
	  if(!isValid) 
	    val = getOutOfBlockValue<Data_T>(iChan,inRange,Vec4i(x,y,z,0),levelMg,
			  sg,sgLo,pos0,bsx,dataBlock0,sgEmptyVal,options);

	  UT_ASSERT2(ih>=0&&ih<hsx*hsx*hsx);

	  transferValue<Data_T>( val, ih, dataDstX, dataDstY, dataDstZ ); 
	}
      }
      #endif 

      // 
      // Middle line segment 
      //
      int inRange = inRangeZY;
      SBG::Block<Data_T> *block=NULL;
      if(inRange)
      {
        int ib = ibzy + ( x0 >> bsxLog2 );
        UT_ASSERT2(ib>=0&&ib<sg->nBlocks());	
        block = sg->getBlock(ib); 
      }

      if(!block)
      {
	#ifndef TEST_NO_BORDERS
	for(int x=x0; x<x0+bsx; x++) 
	{
	  Data_T val = getOutOfBlockValue<Data_T>( 
	      		iChan, inRange, Vec4i(x,y,z,0), levelMg,
			sg,sgLo,pos0,bsx,dataBlock0,sgEmptyVal,options);
	  transferValue<Data_T>( val, ihzy + x-x1, dataDstX, dataDstY, dataDstZ ); 
	}
	#endif
      }
      else 
      {
	int iv = ivzy,
	    ih = ihzy + x0-x1;
	transferInteriorLine<Data_T>( iv, ih, bsx, hsx, sg, 
				      block->_data,
				      dataDstX, dataDstY, dataDstZ );
      }
      
      // 
      // Last line segment
      //
      #ifndef TEST_NO_BORDERS
      if( doHorizontalHalo )
      {
	for(int x=x0+bsx; x<=x2; x++)
	{
	  int inRange = inRangeZY && x>=0 && x<sx,
	      ib = ibzy + ( x >> bsxLog2 ),
	      iv = ivzy + ( x & bsxMask ),
	      ih = ihzy + x-x1;

	  Data_T val;
	  int isValid = FALSE;

	  if(inRange)
	  {
	    UT_ASSERT2(ib>=0&&ib<sg->nBlocks());	
	    SBG::Block<Data_T> *block = sg->getBlock(ib); 	  	    
	    UT_ASSERT2(iv>=0&&iv<sg->nVoxelsInBlock());
	    if(block)
	    {
	      val = block->_data[iv];
	      isValid = TRUE;
	    }
	  }
	  if(!isValid) 
	    val = getOutOfBlockValue<Data_T>( iChan, inRange, Vec4i(x,y,z,0), levelMg,
			  sg,sgLo,pos0,bsx,dataBlock0,sgEmptyVal,options);

	  UT_ASSERT2(ih>=0&&ih<hsx*hsx*hsx);

	  transferValue<Data_T>( val, ih, dataDstX, dataDstY, dataDstZ ); 
	}
      }
      #endif
    }
  }

  for(int k=0;k<nVecCmp_T;k++) dataValid0[k]=TRUE;

  return _haloBlocks+hix0;
}

template 
float **HaloBlockSet::fillHaloBlock_<float,1>(  int iChan,
    		        int bid,   // MSBG block ID
			int levelMg,
    			int tid,   // halo block nr. (thread id)
			unsigned options
    			);
template 
float **HaloBlockSet::fillHaloBlock_<Vec3Float,3>(  int iChan,
    		        int bid,   // MSBG block ID
			int levelMg,
    			int tid,   // halo block nr. (thread id)
			unsigned options
    			);

template 
float **HaloBlockSet::fillHaloBlock_<uint8_t,1>(  int iChan,
    		        int bid,   // MSBG block ID
			int levelMg,
    			int tid,   // halo block nr. (thread id)
			unsigned options
    			);

template 
float **HaloBlockSet::fillHaloBlock_<uint8_t,1,1>(  int iChan,
    		        int bid,   // MSBG block ID
			int levelMg,
    			int tid,   // halo block nr. (thread id)
			unsigned options
    			);

template 
float **HaloBlockSet::fillHaloBlock_<uint16_t,1>(  int iChan,
    		        int bid,   // MSBG block ID
			int levelMg,
    			int tid,   // halo block nr. (thread id)
			unsigned options
    			);
template 
float **HaloBlockSet::fillHaloBlock_<float,1,1>(  int iChan,
    		        int bid,   // MSBG block ID
			int levelMg,
    			int tid,   // halo block nr. (thread id)
			unsigned options
    			);
template 
float **HaloBlockSet::fillHaloBlock_<uint16_t,1,1>(  int iChan,
    		        int bid,   // MSBG block ID
			int levelMg,
    			int tid,   // halo block nr. (thread id)
			unsigned options
    			);
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
float getOutOfBlockValue3(
        MSBG::MultiresSparseGrid *msg,
    	int inRange,
    	Vec4i pos,    // Coords to project
    	const SparseGrid<float> *sg,
    	const SparseGrid<float> *sgLo,
    	const SparseGrid<CellFlags> *sgFlagsLo
	)
{
  if( inRange )
  {
    pos >>= 1;
    int bid = sgLo->getBlockIndex(sgLo->getBlockCoords(pos)),
	vid = sgLo->getVoxelInBlockIndex(sgLo->getVoxelInBlockCoords(pos));
    CellFlags flags = sgFlagsLo->getValueInBlock(bid,vid);
    if(flags & CELL_VOID) 
      return SPEC_FLOAT_DIRICHLET;
    else if(flags & CELL_SOLID) 
      return SPEC_FLOAT_NEUMANN;
    else
      return sgLo->getValueInBlock(bid,vid);
  }
  else 
  {
    return sg->getOutOfDomValue(pos,msg->_ghostFloatSpecVals3x3x3);      		     
  }
}

/*-------------------------------------------------------------------------*/ 
/* 									   */
/*-------------------------------------------------------------------------*/
float *HaloBlockSet::fillHaloBlock3(  int iChan,
    		        int bid, 
			int bx0, int by0, int bz0,   // MSBG block ID
			int levelMg,
    			int tid   // halo block nr. (thread id)
    			)
{
  UT_ASSERT2(_dmax==1);
  UT_ASSERT2(bid>=0&&bid<_msbg->nBlocks());
  UT_ASSERT2(levelMg>=0&&levelMg<_msbg->getNumLevels());
  BlockInfo *bi=_msbg->getBlockInfo(bid,levelMg);  
  int level = bi->level;
  SparseGrid<float> *sg = getChannel<float>(_msbg,iChan,level,levelMg),
                     *sgLo = bi->level<_msbg->getNumLevels()-1 ? 
		       getChannel<float>(_msbg,iChan,level+1,levelMg) : NULL;
  SparseGrid<CellFlags> 
    *sgFlags = _msbg->getFlagsChannel(CH_CELL_FLAGS,level,levelMg),
    *sgFlagsLo = _msbg->getFlagsChannel(CH_CELL_FLAGS,level+1,levelMg);

  float *dataDst = getHaloBlockPtr(tid);

  int bsx=sg->bsx();
  int hsx=bsx+2*_dmax,
      hsx2=hsx*hsx;
  USE(hsx);

  int sx = sg->sx(),
      sy = sg->sy(),
      sz = sg->sz(),    
      bsxLog2 = sg->bsxLog2(),
      bsx2Log2 = sg->bsx2Log2(),
      bsxMask = sg->bsx()-1,
      nbx = sg->nbx(),
      nbxy = sg->nbxy();

  // Loop over 3*3*3 block neighborhood

  Vec4i pos0 = Vec4i(bx0,by0,bz0,0)*bsx,
	pos1 = pos0-_dmax,
	pos2 = pos0+bsx+_dmax-1;

  int x0=vget_x(pos0), x1=vget_x(pos1), x2=vget_x(pos2),
      y0=vget_y(pos0), y1=vget_y(pos1), y2=vget_y(pos2),
      z0=vget_z(pos0), z1=vget_z(pos1), z2=vget_z(pos2);
  USE(y0);
  USE(z0);  

  for(int z=z1; z<=z2; z++)
  {
    int inRangeZ = z>=0&&z<sz,
        ibz = (z >> bsxLog2 ) * nbxy,
	ivz = (z & bsxMask) << bsx2Log2,
	ihz = (z-z1) * hsx2;

    for(int y=y1; y<=y2; y++)
    {		
      int inRangeZY = inRangeZ && y>=0 && y<sy,
	  ibzy = ibz + (y >> bsxLog2) * nbx,
	  ivzy = ivz + ((y & bsxMask) << bsxLog2),
	  ihzy = ihz + (y-y1)*hsx;		

      // 
      // First line segment
      //
      for(int x=x1; x<x0; x++)
      {
        int inRange = inRangeZY && x>=0 && x<sx,
	    ib = ibzy + ( x >> bsxLog2 ),
	    iv = ivzy + ( x & bsxMask ),
	    ih = ihzy + x-x1;
	UT_ASSERT2(ih>=0&&ih<hsx*hsx*hsx);

	if(inRange)
	{
	  UT_ASSERT2(ib>=0&&ib<sg->nBlocks());	
	  UT_ASSERT2(iv>=0&&iv<sg->nVoxelsInBlock());
	  if( _msbg->getBlockInfo(ib,levelMg)->level==levelMg )
	  {
	    CellFlags flags = sgFlags->getValueInBlock0(ib,iv);
	    dataDst[ih] = flags & CELL_VOID ? SPEC_FLOAT_DIRICHLET :
	          	  flags & CELL_SOLID ? SPEC_FLOAT_NEUMANN :
		          sg->getValueInBlock0(ib,iv);
	    continue;
	  }
	}
	dataDst[ih] = getOutOfBlockValue3(_msbg,inRange,Vec4i(x,y,z,0),
			                  sg,sgLo,sgFlagsLo);
      }
      
      // 
      // Middle line segment 
      //

      int inRange = inRangeZY;
      SBG::Block<float> *block=NULL;
      CellFlags *dataFlags=NULL;
      if(inRange)
      {
        int ib = ibzy + ( x0 >> bsxLog2 );

	/*if(levelMg==1 && y==8&&z==0)
        { 
	  TRCERR(("DBG \n"));
        }*/
	
        UT_ASSERT2(ib>=0&&ib<sg->nBlocks());	
	if( _msbg->getBlockInfo(ib,levelMg)->level==levelMg )
	{
	  block = sg->getBlock(ib); 	  
	  dataFlags = sgFlags->getBlockDataPtr0(ib);
	}
      }

      if(!block)
      {
	for(int x=x0; x<x0+bsx; x++) 
	  dataDst[ihzy + x-x1] = 
	    getOutOfBlockValue3( _msbg,inRange, Vec4i(x,y,z,0), sg,sgLo,sgFlagsLo );
      }
      else 
      {
	int iv = ivzy,
	    ih = ihzy + x0-x1;
	float *data = block->_data;


	for(int j=0; j<bsx; j+=4/*SIMDW*/)
	{
	  int iSrc = iv+j,
	      iDst = ih+j;
	  UT_ASSERT2(iSrc>=0&&iSrc<sg->nVoxelsInBlock());	    
	  UT_ASSERT2(iDst>=0&&iDst<hsx*hsx*hsx);
	  Vec4f val = Vec4f().load_a( data+iSrc );
	  setBcSpecValsSIMD<Vec4f>( dataFlags + iSrc, 
				    CELL_SOLID, CELL_VOID,
				    val );
	  val.store( dataDst+iDst );
	}
      }
      
      // 
      // Last line segment
      //
      for(int x=x0+bsx; x<=x2; x++)
      {
        int inRange = inRangeZY && x>=0 && x<sx,
	    ib = ibzy + ( x >> bsxLog2 ),
	    iv = ivzy + ( x & bsxMask ),
	    ih = ihzy + x-x1;
	UT_ASSERT2(ih>=0&&ih<hsx*hsx*hsx);

	if(inRange)
	{
	  UT_ASSERT2(ib>=0&&ib<sg->nBlocks());	
	  UT_ASSERT2(iv>=0&&iv<sg->nVoxelsInBlock());
	  if( _msbg->getBlockInfo(ib,levelMg)->level==levelMg )
	  {
	    CellFlags flags = sgFlags->getValueInBlock0(ib,iv);
	    dataDst[ih] = flags & CELL_VOID ? SPEC_FLOAT_DIRICHLET :
	          	  flags & CELL_SOLID ? SPEC_FLOAT_NEUMANN :
		          sg->getValueInBlock0(ib,iv);
	    continue;
	  }
	}
	dataDst[ih] = getOutOfBlockValue3(_msbg,inRange,Vec4i(x,y,z,0),
			                  sg,sgLo,sgFlagsLo);
      }
    }
  }

  return dataDst;
}

} // namespace MSBG

