/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#ifndef HALO_H
#define HALO_H

#include "msbg.h"


namespace MSBG
{
using namespace SBG;

/*=========================================================================
 *
 *
 *  B L O C K  P R O C E S S O R
 *
 *
 * =======================================================================*/

enum
{
  BPROC_TEST = (1<<1)
};


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/

template< int bsxLog2 >
class BlockProcessor
{
  public:

  void init( SBG::SparseGrid<float> *sg, int bid )
  {
    const Vec4i v4iStencil7[7] = 
	    { Vec4i( 0, 0, 0, 0),
	      Vec4i(-1, 0, 0, 0), Vec4i(1, 0, 0, 0),
	      Vec4i( 0,-1, 0, 0), Vec4i(0, 1, 0, 0),
	      Vec4i( 0, 0,-1, 0), Vec4i(0, 0, 1, 0) };

    Vec4i bpos = sg->getBlockCoordsById(bid);

    for(int i=0;i<7;i++)	
    {
      Vec4i bdir = v4iStencil7[i];
      Vec4i bpos2 = bpos+bdir;
      float *data=NULL;
      if(!sg->blockCoordsInRange(bpos2))
      {
	data = sg->getEmptyBlock()->_data;
      }
      else
      {
        int bid2 = sg->getBlockIndex(bpos2);
	data = sg->getBlockDataPtr(bid2);
      }
      bpos2 = Vec4i(1,1,1,0) + bdir;
      int idx = MT_GXYZ(3,3*3,vget_x(bpos2),vget_y(bpos2),vget_z(bpos2));
      UT_ASSERT0(idx>=0 && idx<ARRAY_LENGTH(_blocksNeigh));
      _blocksNeigh[idx] = data;
    }
  }

  template< int opmode,
	    int options,
	    int SIMDW,
	    typename VecNf,
	    int segT /*(&1=first,&2=last SIMD segment in block*/
		  >
  inline
#ifdef FORCE_INLINE_KERNELS
  __attribute__((always_inline))
#endif
  void procLineSegmentSIMD(
		  int dby0, int dby1, int dbz0, int dbz1,
		  int vx, int vy, int vz,
		  float dt
		  )
  {
    constexpr int bsx = 1<<bsxLog2,
		  bsxMsk = bsx-1;

    constexpr int bid = MT_GXYZ(3,3*3,1,1,1);

  //  char chbuf[128];
    constexpr int DBX=1; //DBY=3,DBZ=3*3;

    VecNf F0,F1,F2,F3,F4,F5,F6,
	  T;

    int vid = VID(vx,vy,vz);

    // Load center cell
    F0.load_a( _blocksNeigh[bid]+vid );

    // 
    // Load Left/Right
    //
    if constexpr (segT==(1|2))  // Segment at left & right border
    {
      // Left
      T.load( _blocksNeigh[bid]+vid-1 ); 
      F1 =  vload_to<0>( T, _blocksNeigh[bid-DBX]+vid+(bsx-1) );
      // Right
      T.load(_blocksNeigh[bid]+vid+1);
      F2 = vload_to<SIMDW-1>( T, _blocksNeigh[bid+DBX]+vid+SIMDW-bsx );
    }
    else if constexpr (segT==1) // Segment at left border
    {
      // Left
      T.load( _blocksNeigh[bid]+vid-1 );
      F1 = vload_to<0>( T, _blocksNeigh[bid-DBX]+vid+(bsx-1));
      // Right
      F2.load( _blocksNeigh[bid]+vid+1 );
    }
    else if constexpr (segT==2)	// Segment at right border
    {
      // Left
      F1.load( _blocksNeigh[bid]+vid-1 );
      // Right
      T.load(_blocksNeigh[bid]+vid+1);
      F2 =  vload_to<SIMDW-1>( T, _blocksNeigh[bid+DBX]+vid+SIMDW-bsx);
    }
    else if(segT==0)  // Segment touches no border
    {
      // Left
      F1.load( _blocksNeigh[bid]+vid-1 );
      // Right
      F2.load( _blocksNeigh[bid]+vid+1 );
    }

    /* 
	  4 6
	1 0 2
	5 3
    */
    int vid3 = VID(vx, (vy-1)&bsxMsk, vz),
	vid4 = VID(vx, (vy+1)&bsxMsk, vz),

	vid5 = VID(vx ,vy, (vz-1)&bsxMsk),
	vid6 = VID(vx ,vy, (vz+1)&bsxMsk);

    // Down
    F3.load_a( _blocksNeigh[bid + dby0] + vid3 );
    // Up
    F4.load_a( _blocksNeigh[bid + dby1] + vid4 );

    // Front
    F5.load_a( _blocksNeigh[bid + dbz0] + vid5 );
    // Back
    F6.load_a( _blocksNeigh[bid + dbz1] + vid6 );

    //
    //  Apply kernel 
    //
    VecNf Y;

    if constexpr (opmode==BPROC_TEST)
    {
      VecNf S = (1.0f/6.0f) * (F1 + F2 +F3 + F4 + F5 + F6),
	    D = dt*(S-F0);
      Y = F0 + D;
    }

    //
    // Store result (optionally using mask & non-temp. hint)
    //
    if constexpr( options & MSBG::OPT_STREAM_STORE )
    {
      vstream(_dataDst+vid,Y);
    }
    else
    {
      Y.store_a(_dataDst+vid);
    }
  }

  //template<int opmode, int options > 
  void processBlock( float dt, float *dataDst )
  {
    constexpr int bsx = 1<<bsxLog2,
	      byStrid = 3,
	      bzStrid = 3*3;

    _dataDst = dataDst;

    for(int vz=0;vz<bsx;vz++)
    {
      const int dbz0 = ((vz-1) >> bsxLog2)*bzStrid,
		dbz1 = ((vz+1) >> bsxLog2)*bzStrid;

      for(int vy=0;vy<bsx;vy++)
      {
	const int dby0 = ((vy-1) >> bsxLog2)*byStrid,
		  dby1 = ((vy+1) >> bsxLog2)*byStrid;

	if constexpr (bsx==16)
	{
	  procLineSegmentSIMD<BPROC_TEST, 0/*OPT_STREAM_STORE*/, 8, Vec8f, 1>(
	      dby0, dby1, dbz0, dbz1, 0, vy, vz, dt );
	  procLineSegmentSIMD<BPROC_TEST, 0/*OPT_STREAM_STORE*/, 8, Vec8f, 2>(
	      dby0, dby1, dbz0, dbz1, 8, vy, vz, dt );
	}
      }
    }
  }

  private:

  inline int VID(int vx, int vy, int vz) 
  {
    return ((vx) + ((vy)<<bsxLog2) + ((vz)<<(2*bsxLog2)));
  }
  
  const float * __restrict _blocksNeigh[3*3*3];
  float *_dataDst;
};
}

/*=========================================================================
 *
 *
 *  Halo Block Set
 *
 *
 * =======================================================================*/

namespace MSBG
{
using namespace SBG;


//
// Halo blocks
//

enum
{
  OPT_HALO_VEC3 = (1<<0),
};


class HaloBlockSet
{    
  public: 
    
  HaloBlockSet( MultiresSparseGrid *msbg,
		int dmax, int nMaxBlocks,
		unsigned options=0 );

  /*-----------------------------------------------------------------------------------*/
  /* 										       */
  /*-----------------------------------------------------------------------------------*/

  template< typename Data_T > inline
  Data_T getOutOfBlockValue(
	  int iChan,
	  int inRange,
	  Vec4i pos,    // Coords to project
	  int levelMg,
	  const SparseGrid<Data_T> *sg,
	  const SparseGrid<Data_T> *sgLo,
	  const Vec4i& block0Pos, int bsx, 
	  Data_T *block0Data, 
	  Data_T emptyVal,
	  unsigned opt
	  )
  {
    if(opt & OPT_BC_NEUMANN)
    {
      pos = min(max( pos, block0Pos ), block0Pos+bsx-1 );
      int vid = sg->getVoxelInBlockIndex(sg->getVoxelInBlockCoords(pos));
      UT_ASSERT2(vid>=0&&vid<sg->nVoxelsInBlock());
      Data_T val = block0Data[ vid ];
      return val;
    }
    else if((opt & OPT_BC_COARSE_LEVEL))
    {
      int x = vget_x(pos),
	  y = vget_y(pos),
	  z = vget_z(pos);

      if(!inRange)
      {
	sg->clipGridCoords(x,y,z);
	Data_T *pVal = sg->getValuePtr( x,y,z );
	if( pVal ) return *pVal;
      }

      if( opt & OPT_INTERPOLATE )
      {
	Data_T val;
	Vec4f posLo( (x+.5f)/2, (y+.5f)/2, (z+.5f)/2, 0.f );
	#if 1
	val = sgLo->interpolateGhost(posLo,0);
	#else
	val = sgLo->template interpolateGhost(posLo,0);
	#endif
	return val;
      }
      else
      { 
	#if 1
	UT_ASSERT2(sgLo);
	#else
	if(!sgLo)
	{
	  Vec4i bpos = sg->getBlockCoords(Vec4i(x,y,z,0)),
		bpos0 = sg->getBlockCoords(block0Pos);

	  Vec4i blockPos = bpos*sg->bsx();

	  TRCERR(("sg->bsx=%d, %d,%d,%d -> %d,%d,%d, inRange=%d, bpos0=%d,%d,%d bflags0=%d/l=%d, bpos=%d,%d,%d bflags=%d/l=%d, \n",
		sg->bsx(), pos[0],pos[1],pos[2], x,y,z,inRange,
		
		block0Pos[0],block0Pos[1],block0Pos[2],
		sg->blockCoordsInRange(bpos0) ?
		  _msbg->getBlockInfo(sg->getBlockIndex(bpos0))->flags : 
		  -1,
		sg->blockCoordsInRange(bpos0) ?
		  _msbg->getBlockInfo(sg->getBlockIndex(bpos0))->level : 
		  -1,

		blockPos[0],blockPos[1],blockPos[2],
		sg->blockCoordsInRange(bpos) ?
		  _msbg->getBlockInfo(sg->getBlockIndex(bpos))->flags : 
		  -1,
		sg->blockCoordsInRange(bpos) ?
		  _msbg->getBlockInfo(sg->getBlockIndex(bpos))->level : 
		  -1));
	}
	#endif
	int xl = x/2,
	    yl = y/2,
	    zl = z/2;
	sgLo->clipGridCoords(xl,yl,zl);
	Data_T *pVal = sgLo->getValuePtr( xl,yl,zl );
	UT_ASSERT2(pVal);
	return *pVal;
      }

    }
    return emptyVal;
  }

  template<typename Data_T> 
  void transferValue( const Data_T& val, 
		      int idxDst,
		      float *dataDstX, float *dataDstY, float *dataDstZ );

  template<typename Data_T> 
  void transferInteriorLine( int iv, int ih, int bsx, int hsx, SparseGrid<Data_T> *sg,  
      			     Data_T *blockData,
		             float *dataDstX, float *dataDstY, float *dataDstZ );
  
  /*-----------------------------------------------------------------------------------*/
  /* 										       */
  /*-----------------------------------------------------------------------------------*/
  int haloBlockIdx( int tid, int iVecCmp=0 ) const 
  { 
    int idx = tid*_nVecCmp + iVecCmp;
    UT_ASSERT0(idx>=0 && idx<_nMaxBlocks*_nVecCmp);
    return idx;
  }

  float *getHaloBlockPtr( int tid, int iVecCmp=0 ) const 
  { 
    return _haloBlocks[haloBlockIdx(tid,iVecCmp)]; 
  }


  float **fillHaloBlockFast16Uint8(  SparseGrid<uint8_t> *sg,
    		        int bid,  
    			int tid,   // halo block nr. (thread id)
			unsigned options
    			);

  template<typename Data_T, int nVecCmp_T=1, int do1stOrderOnly_T=0 > 
  float **fillHaloBlock_(  int iChan,
			  int bid,   // MSBG block ID
			  int levelMg,
			  int tid,   // halo block nr. (thread id)
			  unsigned options=0
			  );
  float *fillHaloBlock3(  int iChan,
		      int bid, 
		      int bx0, int by0, int bz0,   // MSBG block ID
		      int levelMg,
		      int tid   // halo block nr. (thread id)
		      );

  size_t bufLenAllocated( void ) const { return _bufLenAllocated; }

  ~HaloBlockSet();

  private:
     //CellFlags *_hblockFlags[256];
     MultiresSparseGrid *_msbg;
     int _dmax,
	 _nMaxBlocks,
	 _nVecCmp,
         _maxHBlockSize,
         _bufLenAllocated;
     float **_haloBlocks;
     int *_dataValid;

    /*float *_tmpBlkDataOutOfDom=NULL;
    CellFlags *_tmpBlkDataFlagsOutOfDom=NULL;*/
};

template<> inline
void HaloBlockSet::transferValue( const Vec3Float& val, 
		    int idxDst,
		    float *dataDstX, float *dataDstY, float *dataDstZ )
{
  dataDstX[idxDst] = val.x;
  dataDstY[idxDst] = val.y;
  dataDstZ[idxDst] = val.z;    
}

template<> inline
void HaloBlockSet::transferValue( const float& val, 
		    int idxDst,
		    float *dataDstX, float *dataDstY, float *dataDstZ )
{
  dataDstX[idxDst] = val;
}

template<> inline
void HaloBlockSet::transferValue( const uint16_t& val, 
		    int idxDst,
		    float *dataDstX, float *dataDstY, float *dataDstZ )
{
  dataDstX[idxDst] = renderDensToFloat_<uint16_t>(val);
}

template<> inline
void HaloBlockSet::transferValue( const uint8_t& val, 
		    int idxDst,
		    float *dataDstX, float *dataDstY, float *dataDstZ )
{
  dataDstX[idxDst] = renderDensToFloat_<uint8_t>(val);
}

template<> inline
void HaloBlockSet::transferInteriorLine( int iv, int ih, int bsx, int hsx, 
    			   SparseGrid<float> *sg,       
			   float *blockData,
			   float *dataDstX, float *dataDstY, float *dataDstZ )
{
  if(bsx>2)
  {
    for(int j=0; j<bsx; j+=4/*SIMDW*/)
    {
      Vec4f val;
      int iSrc = iv+j,
	  iDst = ih+j;
      UT_ASSERT2(iSrc>=0&&iSrc<sg->nVoxelsInBlock());	    
      val.load_a( blockData+iSrc );
      UT_ASSERT2(iDst>=0&&iDst<hsx*hsx*hsx);
      val.store( dataDstX+iDst );
    }
  }
  else
  {
    int iSrc = iv,
	iDst = ih;
    dataDstX[iDst] = blockData[iSrc];
    dataDstX[iDst+1] = blockData[iSrc+1];
  }
}

template<> inline
void HaloBlockSet::transferInteriorLine( int iv, int ih, int bsx, int hsx, 
    			   SparseGrid<uint16_t> *sg,       
			   uint16_t *blockData,
			   float *dataDstX, float *dataDstY, float *dataDstZ )
{
  for(int j=0; j<bsx; j+=4/*SIMDW*/)
  {
    Vec4f val;
    int iSrc = iv+j,
	iDst = ih+j;
    UT_ASSERT2(iSrc>=0&&iSrc<sg->nVoxelsInBlock());	        
    vload_from_uint16_a( val, blockData+iSrc );
    RENDER_DENS_TO_FLOAT_UI16(val);
    UT_ASSERT2(iDst>=0&&iDst<hsx*hsx*hsx);
    val.store( dataDstX+iDst );
  }
}

template<> inline
void HaloBlockSet::transferInteriorLine( int iv, int ih, int bsx, int hsx, 
    			   SparseGrid<uint8_t> *sg,       
			   uint8_t *blockData,
			   float *dataDstX, float *dataDstY, float *dataDstZ )
{
  #if INSTRSET >= 8 // AVX 2
  if(bsx>=16)
  {
    // Need at least 8 bytes to fully utilize memory bus transfers 
    for(int j=0; j<bsx; j+=16)
    {
      Vec8f val;
      int iSrc = iv+j,
	  iDst = ih+j;
      UT_ASSERT2(iSrc>=0&&iSrc<sg->nVoxelsInBlock());	    
      UT_ASSERT2(iDst>=0&&iDst<hsx*hsx*hsx);

      __m128i temp = _mm_load_si128((__m128i*)(blockData+iSrc));
      __m128i temp2 = temp;
      __m256i temp256 = _mm256_cvtepu8_epi32(temp2);
      val = _mm256_cvtepi32_ps(temp256);
      RENDER_DENS_TO_FLOAT_UI8(val);
      UT_ASSERT2(iDst>=0&&iDst<hsx*hsx*hsx);
      val.store( dataDstX+iDst );

      temp2 = permute4<2,3,-1,-1>((Vec4ui)temp);   // upper 64bit to lower
      temp256 = _mm256_cvtepu8_epi32(temp2);
      val = _mm256_cvtepi32_ps(temp256);
      RENDER_DENS_TO_FLOAT_UI8(val);
      UT_ASSERT2(iDst>=0&&iDst<hsx*hsx*hsx);
      val.store( dataDstX+iDst+8 );

      #if 0
      {
	for(int i=0;i<16;i++)
	{
	  uint8_t 
	    x = blockData[iSrc+i],
	    y = renderDensFromFloat_<uint8_t>(dataDstX[iDst+i]);
	  if(x!=y)
	  {
	    TRCERR(("%d %d\n",x,y));
	  }
	}
      }
      #endif
    }
  }
  else
  #endif
  {
    for(int j=0; j<bsx; j+=4/*SIMDW*/)
    {
      Vec4f val;
      int iSrc = iv+j,
	  iDst = ih+j;
      UT_ASSERT2(iSrc>=0&&iSrc<sg->nVoxelsInBlock());	    
      vload_from_uint8_a( val, blockData+iSrc );
      RENDER_DENS_TO_FLOAT_UI8(val);
      UT_ASSERT2(iDst>=0&&iDst<hsx*hsx*hsx);
      val.store( dataDstX+iDst );
    }
  }
}

template<> inline
void HaloBlockSet::transferInteriorLine( int iv, int ih, int bsx, int hsx, 	
    			   SparseGrid<Vec3Float> *sg, 
			   Vec3Float *blockData,
			   float *dataDstX, float *dataDstY, float *dataDstZ )
{
  for(int j=0; j<bsx; j++)
  {
    Vec3Float val;
    int iSrc = iv+j,
	iDst = ih+j;

    UT_ASSERT2(iSrc>=0&&iSrc<sg->nVoxelsInBlock());
    
    val = blockData[iSrc];
    UT_ASSERT2(iDst>=0&&iDst<hsx*hsx*hsx);

    transferValue<Vec3Float>( val, iDst, dataDstX, dataDstY, dataDstZ ); 
  }
}

} // namespace MSBG


#endif // HALO_H
