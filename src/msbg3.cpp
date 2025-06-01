/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#include "common_headers.h"

#include "omp.h"

//#define DBG_RELAX

//#define MSBG_AR_SYMMETRIC

//#define DISCARD_BLOCK_SMALL_CHANGE

namespace MSBG 
{


namespace TST
{


} // namespace TST

static PSFloat *tmpBlkDataXBuf[MI_MAX_THREADS][3] = {0};
static float *tmpBlkDataWBuf[MI_MAX_THREADS][6] = {0};
static uint8_t *tmpBlkDataFBuf[MI_MAX_THREADS] = {0};
static int tmpBlkErrorStatus[MI_MAX_THREADS] = {0};
static int tmpBlkMaxLen,
	   tmpBlkMaxThreads;

/*=========================================================================
 *
 *  
 *  
 *  
 *
 * =======================================================================*/

#define SUB_CONTRIB( dx, dy, dz, sumVal, sumCoeff0 ) \
{ \
    int sbidx = sbix<dx,dy,dz>(bsxHi,bsx2Hi); \
    vid2 = vid2bas + sbidx;  \
    a = h * DPMR_GRAD_COEFF_CF;\
    if( dataFA2 ) a *= dataFA2[vid2bas0 + sbidx]; \
    sumVal += a * data2[vid2];\
    sumCoeff0 += a;\
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template< int bsx, 
  	  int dir_T, // normal direction (x=0,y=1,z=2)
  	  int face_T >  // face (0,1)
static inline
void procNeighMixedRes2( 
    		   const int&level, const int&level2,
    		   const double& h,
    		   CellFlags *dataFlg2, PSFloat *data2, float *dataFA2,
		   int vx, int vy, int vz, int vid,
    		   double& sumVal, double& sumCoeff0 )
{
  int vid2,vid2bas,vid2bas0;
  PSFloat a;

  constexpr int bsx2 = bsx*bsx;
  constexpr int bsxHi = 2*bsx;
  constexpr int bsx2Hi = bsxHi*bsxHi;
  
  if(level2<level)  // coarse - fine
  {
    switch(dir_T)
    {
      case 0:  // left/right 
        vid2bas0 = (8*bsx2)*vz + (4*bsx)*vy;
	vid2bas = face_T==0 ? (8*bsx2)*vz + (4*bsx)*vy + (2*bsx-1) : vid2bas0;
	SUB_CONTRIB(0,0,0, sumVal, sumCoeff0 );
	SUB_CONTRIB(0,1,0, sumVal, sumCoeff0 );
	SUB_CONTRIB(0,0,1, sumVal, sumCoeff0 );
	SUB_CONTRIB(0,1,1, sumVal, sumCoeff0 );
	break;

      case 1: // down/up
	vid2bas0 = (8*bsx2)*vz+2*vx;
	vid2bas = face_T==0 ? (8*bsx2)*vz + 2*vx + (2*bsx*(2*bsx-1)) : vid2bas0;
	SUB_CONTRIB(0,0,0, sumVal, sumCoeff0 );
	SUB_CONTRIB(1,0,0, sumVal, sumCoeff0 );
	SUB_CONTRIB(0,0,1, sumVal, sumCoeff0 );
	SUB_CONTRIB(1,0,1, sumVal, sumCoeff0 );
	break;

      case 2: // front/back
	vid2bas0 = (4*bsx)*vy + 2*vx;
	vid2bas = face_T==0 ? (4*bsx)*vy + 2*vx + (4*bsx2*(2*bsx-1)) : vid2bas0;
	SUB_CONTRIB(0,0,0, sumVal, sumCoeff0 );
	SUB_CONTRIB(1,0,0, sumVal, sumCoeff0 );
	SUB_CONTRIB(0,1,0, sumVal, sumCoeff0 );
	SUB_CONTRIB(1,1,0, sumVal, sumCoeff0 );
	break;
    }
  }
  else  // fine - coarse
  {
    switch(dir_T)
    {
      case 0: vid2 = face_T==0 ? 
		((bsx2)/4) * (vz/2) + (bsx/2) * (vy/2) + bsx/2-1 : 
		((bsx2)/4) * (vz/2) + bsx/2 * (vy/2);
		break;

      case 1: vid2 = face_T==0 ?
		(bsx2/4)*(vz/2) + (vx/2) + ((bsx/2-1)*bsx/2) :
		(bsx2/4)*(vz/2) + (vx/2);
		break;

      case 2: vid2 = face_T==0 ? 
		(bsx/2)*(vy/2) + (vx/2) + ((bsx/2-1)*bsx2/4) :
		(bsx/2)*(vy/2) + (vx/2);
		break;		  		
    
    }
    
    int vidFA = vid;
    if(face_T==1)
    {
      switch(dir_T)
      {
	case 0: vidFA = MT_GXYZ(bsx,bsx2, 0,vy,vz); break;
	case 1: vidFA = MT_GXYZ(bsx,bsx2, vx,0,vz); break;
	case 2: vidFA = MT_GXYZ(bsx,bsx2, vx,vy,0); break;
      }
    }

    a = h * 2 * DPMR_GRAD_COEFF_CF;
    if(dataFA2) a *= dataFA2[vidFA];

    sumVal += a * data2[vid2];
    sumCoeff0 += a;

#if 0
#if 0
      UT_ASSERT2(std::isfinite(sumVal));
#else
      if(!std::isfinite(sumVal)||fabs(sumVal)>1e9)
      {
	TRCP(("%d/%d: bsx=%d %d/%d, a=%g sumVal=%g, %d,%d,%d, vid=%d, vid2=%d, data2[vid2]=%g\n", 
	      dir_T,face_T,
	      bsx,level,level2,a,sumVal,vx,vy,vz,vid,vid2,data2[vid2]));
	UT_ASSERT0(FALSE);
      }
#endif
#endif
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static inline 
void ImultiplyLaplacianMatrixProcessSliceDenseLevel(
    		        int iIter, int pass_red_black,
    			int z0,  // z-slice to process			
			int doJacobi,
			int doReverseOrder,
    			int doBoundZoneOnly,
			int doCalcResidual,
			int sx, int sy, int sz,  // grid dimensions
			int yStrid, int zStrid,
			PSFloat h, PSFloat omega,
    			SparseGrid<CellFlags>*sg,
			CellFlags *dataFlags,
			PSFloat *dataX,
			PSFloat *dataB,
			PSFloat *dataD,
			float **dataFaceArea,
			PSFloat *dataDstY,
			bool doBypassCache,
			bool usePhysicalAirVelocity,
			UtPerfStats *pfs=NULL
			)
{
  UT_ASSERT0(dataD);
  const	int SIMDW=4;
  int z = z0;

  const Vec4psf HROmega( omega * 1.0f/h ),
	        oneMinusOmega( 1.0f - omega );

  #ifdef MSBG_PERF_STATS
  UtPerfStats pfs0={};
  if(!pfs) pfs=&pfs0;
  #endif

#ifdef MG_BLOCK_GAUSS_SEIDEL_DENSELEVELS
    int x1,x2, y1,y2, dx,dy;

    if( !!doReverseOrder == !!(iIter&1) )
    {
      x1=ALIGN(sx,SIMDW)-SIMDW; x2=-SIMDW; dx=-SIMDW;
      y1=sy-1; y2=-1; dy=-1;
    }
    else
    {
      x1=0; x2=ALIGN(sx,SIMDW); dx=SIMDW;
      y1=0; y2=sy; dy=1;
    }

    for(int y=y1; y!=y2; y+=dy)
    {
      int vid =  sg->getGridIndex(x1,y,z);
      for(int x=x1; x!=x2; x+=dx, vid+=dx)
      {

#else

  for(int y=0;y<sy;y++)
  {
    int vid = sg->getGridIndex(0,y,z);
    for(int x=0; x<sx; x+=SIMDW,vid+=SIMDW)
    {	

#endif

      Vec4psf F0,F1,F2,F3,F4,F5,F6,F7,
	      Y,D;
      Vec4ui M(-1),  // mask non-fluid cells (set zero)
	     MS(-1),
	     Mflg;

      // Set mask (zero for non-fluid zells)
      {
	Vec8us M8us;
	// Load cell flags (8*16 bit) and convert to 4*32 bit
	M8us.load( dataFlags+vid );
	
	Mflg = _mm_unpacklo_epi16(M8us, _mm_set1_epi16(0));

	M = _mm_cmpeq_epi32( Mflg & (CELL_SOLID|CELL_VOID), 
			    _mm_set1_epi32(0)) ;    

        #ifdef MSBG_PERF_STATS
	pfs->bytes += SIMDW*2;
	pfs->flops += 1*SIMDW;
	#endif

	// skip non-fluid or fixed cells 
	if(_mm_movemask_ps(_mm_castsi128_ps(M))==0x0)
	{
	  //UT_ASSERT2( vall(Y.load(dataDstY+vid) == 0.0f)); 
	  continue;  
	}
	
	#if 0
	// do not touch fixed or non-boundzone cells
	{
	  MS =_mm_cmpeq_epi32( Mflg & CELL_FIXED, _mm_setzero_si128());
	  if( doBoundZoneOnly )
	  {
	    MS &= _mm_cmpeq_epi32( Mflg & (CELL_BOUNDARY_ZONE), 
				  _mm_set1_epi32(CELL_BOUNDARY_ZONE));
	  }

	  if( _mm_movemask_ps(_mm_castsi128_ps(MS)) == 0x0) continue;
	}

        #else
	if(doBoundZoneOnly)
	{
	  // Set storage mask to not touch cells outside bundary zone or fixed cells
	  MS = _mm_cmpeq_epi32( Mflg & (CELL_BOUNDARY_ZONE), 
			    _mm_set1_epi32(CELL_BOUNDARY_ZONE));
	  if( _mm_movemask_ps(_mm_castsi128_ps(MS))==0x0 ||
	      _mm_movemask_ps(_mm_castsi128_ps(M))==0x0 )
	    continue;
	}
	#endif

      }
      D.load(dataD+vid);

      // Load cell neighbors
      PSFloat *p=dataX+vid;
      F0.load( p );
      F1.load( p - 1 );
      F2.load( p + 1 );
      F3.load( p - yStrid );
      F4.load( p + yStrid );	 
      F5.load( p - zStrid );
      F6.load( p + zStrid );

      if( usePhysicalAirVelocity || 
	  vany(    // any non-fully weighted neighbor ?
	    (Mflg & (CELL_SOLID|CELL_PARTIAL_SOLID)) != (Vec4ui)_mm_set1_epi32( 0 )))
      {
	Vec4psf W1 = vecf2vecpsf(Vec4f().load( dataFaceArea[0] + vid )),   // U left
	        W2 = vecf2vecpsf(Vec4f().load( dataFaceArea[0] + vid + 1)), // U right 
	        W3 = vecf2vecpsf(Vec4f().load( dataFaceArea[1] + vid )),   // V down
	        W4 = vecf2vecpsf(Vec4f().load( dataFaceArea[1] + vid + yStrid)), // V up 
	        W5 = vecf2vecpsf(Vec4f().load( dataFaceArea[2] + vid )),   // W front
	        W6 = vecf2vecpsf(Vec4f().load( dataFaceArea[2] + vid + zStrid)); // W back 

        #if 0
	#ifdef UT_ASSERT_LEVEL_2
	//#ifdef UT_ASSERT_LEVEL_3
	// Check diagonally dominant
	if(iIter==0)
	{
	  //UT_ASSERT0(!vany( W1+W2+W3+W4+W5+W6 > + 1e-3f + 1.0f/D ));
	  UT_ASSERT0(!vany( W1+W2+W3+W4+W5+W6 > 1.0f/D ));
	}
	#endif
        #endif
	
	F1 *= W1;
	F2 *= W2;
	F3 *= W3;
	F4 *= W4;
	F5 *= W5;
	F6 *= W6;

	#ifdef MSBG_PERF_STATS
	pfs->bytes += 6*SIMDW*4;
	pfs->flops += 6*SIMDW;
	#endif
      }

      // Calculate result
      Vec4psf H(h),
	      S = ( F1+F2+F3+F4+F5+F6 ) * H,
	      B;

      #ifdef MSBG_PERF_STATS
      pfs->bytes += 8*SIMDW*4;  // F0,D
      pfs->flops += (1+6)*SIMDW;
      #endif

      if(doJacobi)
      {
	B.load(dataB+vid);

	Y = (B+S)*D*HROmega  + oneMinusOmega * F0;
	Y = vmask4(Y,M);

	#ifdef MSBG_PERF_STATS
	pfs->bytes += 1*SIMDW*4;  // B
	pfs->flops += 5*SIMDW;
	#endif
      }
      else
      {
	Vec4psf N = 1.0f/D;
	Y = H*N*F0 - S;
	Y = vmask4(Y,M);

	#ifdef MSBG_PERF_STATS
	pfs->flops += 4*SIMDW;
	#endif

	if(doCalcResidual)
	{
	  B.load(dataB+vid);
	  Y = B - Y;
	  #ifdef MSBG_PERF_STATS
	  pfs->bytes += SIMDW*4;  // B
	  pfs->flops += 1*SIMDW;
	  #endif
	}
      }
      // 
      // Store result	
      //
      Y.cutoff(sx-x); // make sure to leave buffer zone zero     

      UT_ASSERT2(visfinite(Y));
      
      #ifndef MG_BLOCK_GAUSS_SEIDEL_DENSELEVELS
      if(doBypassCache)
      {
	UT_ASSERT0(FALSE); // verify correct masking 
	#if 0
	_mm_maskmoveu_si128( _mm_castps_si128(Y), MS,  
			     (char*)(dataDstY+vid));  
	#endif
      }
      else
      #endif
      {
        vmaskstore4( dataDstY+vid, Y, MS );
      }
      
      #ifdef MSBG_PERF_STATS
      pfs->bytes += SIMDW*4;
      #endif
    }  // x
  }  // y      
}

/*=========================================================================
 *
 *
 *  Block relaxation with halo-cells and cache-local sub-iterations
 *
 *
 * =======================================================================*/

//#define ADAPTIVE_BLOCK_LOCAL_RELAX

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template< int opmode, // LAPL_RELAX, LAPL_MULTIPLY
	  int bsx, 
    	  int SIMDW, 
	  typename VecNpfs, typename VecNf, typename VecNui >
void MultiresSparseGrid::IprocessBlockLaplacian(
		    int tid,	// Thread number (-1=init buffers)
		    unsigned options,
    		    int bid, int bx, int by, int bz,
		    BlockInfo *bi,
		    int levelMg,
		    int chanX,
		    int chanB,
		    int chanDstY, 
		    int nIter,
		    Vec4d *v4dAlpha,
		    PSFloat *pMaxScaledResid
		    )
{
  //#define TEST_COPY_ONLY

  #undef SKIP_NONFLUID_SIMD_STATS
  //#define SKIP_NONFLUID_SIMD_STATS

  #ifdef SKIP_NONFLUID_SIMD_STATS
  MT::Stats<float>  skipstats;
  #endif

  PSFloat residThetaAct = pMaxScaledResid ? *pMaxScaledResid : -1;
  USE(residThetaAct);
  const bool blockHasFaceCoeffs = (bi->flags & BLK_CUTS_SOLID);

  int doBoundZoneOnly = options & OPT_BOUNDZONE_ONLY,
      blockHasSkippableSIMD = 0;
  USE(doBoundZoneOnly);
  //bool blockHasFullDiag = bi->flags & BLK_FULL_DIAG;

  PSFloat *residDataX = NULL;
  if( (opmode & LAPL_RELAX) && (options & OPT_CALC_PER_BLOCK_RESID ))
  {
    residDataX = tmpBlkDataXBuf[tid][2];
  }
    
  UT_ASSERT2(levelMg<_nLevels);  
  UT_ASSERT2(bid<_nBlocks);
  int level = bi->level;

  PSFloat h = (1<<level)*_dx0;

  SparseGrid<CellFlags> *sg = getFlagsChannel(CH_CELL_FLAGS,level,levelMg),
			*sgHi = getFlagsChannel(CH_CELL_FLAGS,level-1,levelMg),
			*sgLo = getFlagsChannel(CH_CELL_FLAGS,level+1,levelMg);

  UT_ASSERT0( bsx==sg->bsx());

  const int bsx2 = bsx*bsx;
  
  UT_ASSERT2(sg->blockCoordsInRange(bx,by,bz));

  SparseGrid<PSFloat> 
    *sgX = getPSFloatChannel(chanX, level, levelMg),
    *sgXHi = getPSFloatChannel(chanX, level-1, levelMg),
    *sgXLo = getPSFloatChannel(chanX, level+1, levelMg),
    *sgB = getPSFloatChannel(chanB, level, levelMg),
    *sgDiag = getPSFloatChannel(CH_DIAGONAL,level,levelMg),
    *sgDstY = getPSFloatChannel(chanDstY, level, levelMg );

  PSFloat *dataDstY = sgDstY->getBlockDataPtr(bid,1,0),
	  *dataB = sgB ? sgB->getBlockDataPtr(bid,0,0) : NULL,
	  *dataDiag = sgDiag ? sgDiag->getBlockDataPtr(bid,0,0) : NULL;

  SparseGrid<float> *sgFaceArea[3];

  for(int k=0;k<3;k++) sgFaceArea[k] = getFaceAreaChannel(k,level,levelMg);

  float 
    *dataFaceAreaU = sgFaceArea[0] ? sgFaceArea[0]->getBlockDataPtr(bid) : NULL,
    *dataFaceAreaV = sgFaceArea[1] ? sgFaceArea[1]->getBlockDataPtr(bid) : NULL,
    *dataFaceAreaW = sgFaceArea[2] ? sgFaceArea[2]->getBlockDataPtr(bid) : NULL;

  UT_ASSERT0(dataDiag);

  CellFlags *dataFlags = sg->getBlockDataPtr(bid); 

  PSFloat *dataX = sgX->getBlockDataPtr(bid); 

  UT_ASSERT2(!(bi->flags & BLK_NO_FLUID ));
  UT_ASSERT2(!(bi->flags & BLK_FIXED ));

  //
  // Prepare halo block with ghost cell layers
  //
  UT_ASSERT0(tid<MI_MAX_THREADS);

  PSFloat *readDataX = tmpBlkDataXBuf[tid][0], 
	  *writeDataX = readDataX;

  if( opmode & LAPL_RELAX )
  {
    #ifdef RELAX_BLOCK_GAUSS_SEIDEL
    writeDataX = readDataX;
    #else
    UT_ASSERT0(FALSE);  // no longer up to date
    writeDataX = tmpBlkDataXBuf[tid][1];
    #endif
  }

  float	**tmpDataW = tmpBlkDataWBuf[tid];
  uint8_t *tmpDataF = tmpBlkDataFBuf[tid];

  const int nVoxelsInBlock = bsx*bsx*bsx;
  const int hsx = bsx+2,
      hsxAct = ALIGN( bsx+2*SIMDW, SIMDW ),
      hsx2Act = hsx*hsxAct,
      hsxActLen = hsxAct*hsx*hsx;

  USE(hsxActLen);
  USE(nVoxelsInBlock);


#ifdef TEST_COPY_ONLY
    if(opmode & LAPL_RELAX)
    {
      int iSimd=0;
      for(int vz=0,vid=0;vz<bsx;vz++)
      for(int vy=0;vy<bsx;vy++)
      {
	int vid2 = MT_GXYZ(hsxAct,hsx2Act, SIMDW,1+vy,1+vz);
        for(int vx=0; vx<bsx; vx+=SIMDW,vid+=SIMDW,vid2+=SIMDW,iSimd++)
	{
	    VecNf X;
	    //X.load_a( readDataX+vid2 );
	    X.load_a( dataB+vid );
	    X.store_a( dataDstY+vid);
	}
      }
      return;
    }
#endif


  VecNpfs X,H(h);
  int idxSrc,idxDst;
  PSFloat *dataSrc;
  VecNf W;

  // Copy interior data
  for(int vz=0, iSimd=0; vz<bsx; vz++)
  for(int vy=0; vy<bsx; vy++)
  {
    idxSrc = MT_GXYZ(bsx,bsx2, 0, vy,vz),
    idxDst = MT_GXYZ(hsxAct,hsx2Act, SIMDW,1+vy,1+vz);
    
    UT_ASSERT2(idxSrc>=0&&idxSrc+bsx-1<nVoxelsInBlock);
    UT_ASSERT2(idxDst>=0&&idxDst+bsx-1<hsxActLen);

    for(int vx=0; vx<bsx; vx+=SIMDW, iSimd++)
    {
      VecNui flags;
      vload_from_uint16_a( flags, dataFlags+idxSrc+vx );

      int doUseFaceCoeffs = horizontal_or( (flags & (CELL_SOLID|CELL_PARTIAL_SOLID)) != 0 );
      
      // No fixed cells allowed in non-fixed blocks
      UT_ASSERT2( !horizontal_or( (flags & (CELL_FIXED)) != 0 ));

      uint8_t tmpFlags = !!doUseFaceCoeffs;

      int doSkip = 0;

      //if constexpr ( opmode & LAPL_RELAX )
      {
	doSkip = ! horizontal_or( (flags & (CELL_SOLID|CELL_VOID)) == 0 );
	tmpFlags |= ((!!doSkip)<<1);
	blockHasSkippableSIMD |= doSkip;
      }

      tmpDataF[ iSimd ] = tmpFlags;
 
      #ifdef SKIP_NONFLUID_SIMD_STATS
      skipstats.update( doSkip ? 1 : 0 );
      #endif

      if(doSkip || !dataX ) 
      {
	X = 0.0f;
	if( residDataX ) X.store_a( residDataX+idxSrc+vx);
      }
      else
      {
        X.load_a(dataX+idxSrc+vx);
	if( residDataX ) X.store_a( residDataX+idxSrc+vx);
        X *= H;
      }

      X.store_a(readDataX+idxDst+vx);
    }

#ifdef TEST_COPY_ONLY
    if(opmode & LAPL_RELAX)
    {
      iSimd=0;
      for(int vz=0,vid=0;vz<bsx;vz++)
      for(int vy=0;vy<bsx;vy++)
      {
	int vid2 = MT_GXYZ(hsxAct,hsx2Act, SIMDW,1+vy,1+vz);
        for(int vx=0; vx<bsx; vx+=SIMDW,vid+=SIMDW,vid2+=SIMDW,iSimd++)
	{
	    VecNf X;
	    //X.load_a( readDataX+vid2 );
	    X.load_a( dataB+vid );
	    X.store_a( dataDstY+vid);
	}
      }
      return;
    }
#endif

    if(blockHasFaceCoeffs)
    {
      for(int vx=0; vx<bsx; vx+=SIMDW)
      {
	int idx = idxSrc+vx;
	W.load_a( dataFaceAreaU + idx );
	W.store_a( tmpDataW[0] + idx );
	W.load( dataFaceAreaU + idx + 1 );  // fetching one over block boundary assumed safe
	W.store_a( tmpDataW[1] + idx  );

	W.load_a( dataFaceAreaV + idx );
	W.store_a( tmpDataW[2] + idx );
	if( vy < bsx-1 )  // right hand block border is filled later
	{
	  W.load_a( dataFaceAreaV + idx + bsx);
	  W.store_a( tmpDataW[3] + idx );
	}

	W.load_a( dataFaceAreaW + idx );
	W.store_a( tmpDataW[4] + idx );
	if( vz < bsx-1 )
	{
	  W.load_a( dataFaceAreaW + idx + bsx2);
	  W.store_a( tmpDataW[5] + idx );
	}
      }
    }
  }
      
  // 
  // Copy the border slabs
  //
  for(int iDir=0; iDir<3; iDir++)
  {
    const int offsets[3][2][3] = { { {-1,0,0}, {1,0,0} },
				 {   {0,-1,0}, {0,1,0} },
				 {   {0,0,-1}, {0,0,1} } };
    for(int iFace=0; iFace<2; iFace++)
    {
      // 
      // Gather info about the the neighbor block
      //
      const int *offs=&offsets[iDir][iFace][0];
      int bx2=bx+offs[0],
	  by2=by+offs[1],
	  bz2=bz+offs[2];

      PSFloat *data2=NULL;
      float *dataFaceArea2=NULL;
      CellFlags *dataFlg2=NULL;
      BlockInfo *bi2 = NULL;
      int level2=level;

      if(sg->blockCoordsInRange(bx2,by2,bz2))
      {
	int bid2 = sg->getBlockIndex(bx2,by2,bz2);
	bi2 = getBlockInfo(bid2,levelMg);
	level2 = bi2->level;

	if(level2==level)  // FF
	{
	  // FINE-FINE
	  data2 = sgX->getBlockDataPtr(bid2, 0,0);
	  dataFaceArea2 = sgFaceArea[iDir]->getBlockDataPtr(bid2,0,0);
	}
	else if(level2>level)  // FC
	{
	  data2 = sgXLo->getBlockDataPtr(bid2, 0,0); 
	  dataFlg2 = sgLo->getBlockDataPtr(bid2, 0,0);
	  if( blockHasFaceCoeffs ) 
	  {
	    auto sg = getFaceAreaChannel(iDir,level,levelMg);

	    dataFaceArea2 = iFace ? sg->getBlockDataPtr(bid2) :
	    			    sg->getBlockDataPtr(bid);
	  }
	}
	else   // CF
	{
	  data2 = sgXHi->getBlockDataPtr(bid2, 0,0); 
	  dataFlg2 = sgHi->getBlockDataPtr(bid2, 0,0);
	  if( blockHasFaceCoeffs ) 
	  {
	    auto sg = getFaceAreaChannel(iDir,level2,levelMg);		

	    dataFaceArea2 = iFace ? sg->getBlockDataPtr(bid2) :
	    			    sg->getBlockDataPtr(bid);
	  }
	}

	UT_ASSERT2(data2);
      }
      else 
      {
	data2 = _sparseGrids[level].zeroBlockPSFloat;
	UT_ASSERT2(_sparseGrids[level].zeroBlockPSFloatSz == 
		    bsx*bsx*bsx*sizeof(PSFloat));
	int domBc = domainBc(iDir,iFace);
	dataFlg2 = getGhostFlagsBlock(level,domBc);

	if( blockHasFaceCoeffs )
	  dataFaceArea2 = iFace ?
	      domBc != DBC_OPEN ? _sparseGrids[level].zeroBlock :				  
			          _sparseGrids[level].constFaceCoeffAirBlock :
	      			  sgFaceArea[iDir]->getBlockDataPtr(bid,0,0);
	level2=level;
      }


      UT_ASSERT2(ABS(level2-level)<=1);

      // 
      // Copy data from neighbor to temp. halo-block
      //

      if(level2==level)  // ======== FINE-FINE ============================
      {
	if(iDir==0) // X - Z,Y
	{
	  if(iFace==0)
	  {
	    for(int vz=0; vz<bsx; vz++)
	    for(int vy=0; vy<bsx; vy++)
	    {
	      idxSrc = MT_GXYZ(bsx,bsx2, bsx-1, vy,vz),
	      idxDst = MT_GXYZ(hsxAct,hsx2Act, SIMDW-1, 1+vy,1+vz);

	      UT_ASSERT2(idxSrc>=0&&idxSrc<nVoxelsInBlock);
	      UT_ASSERT2(idxDst>=0&&idxDst+SIMDW-1<hsxActLen);

	      PSFloat val = data2[idxSrc] * h;
	      readDataX[idxDst] = val;
	      #ifndef RELAX_BLOCK_GAUSS_SEIDEL
	      if(opmode & LAPL_RELAX ) writeDataX[idxDst] = val;
	      #endif
	    }
	  }
	  else
	  {
	    for(int vz=0; vz<bsx; vz++)
	    for(int vy=0; vy<bsx; vy++)
	    {
	      idxSrc = MT_GXYZ(bsx,bsx2, 0, vy,vz),
	      idxDst = MT_GXYZ(hsxAct,hsx2Act, hsxAct-SIMDW, 1+vy,1+vz);

	      UT_ASSERT2(idxSrc>=0&&idxSrc<nVoxelsInBlock);
	      UT_ASSERT2(idxDst>=0&&idxDst+SIMDW-1<hsxActLen);

	      PSFloat val = data2[idxSrc] * h;
	      readDataX[idxDst] = val;
	      #ifndef RELAX_BLOCK_GAUSS_SEIDEL
	      if(opmode & LAPL_RELAX ) writeDataX[idxDst] = val;
	      #endif

	      if(blockHasFaceCoeffs) 
		*(tmpDataW[1] + MT_GXYZ(bsx,bsx2, bsx-1, vy,vz)) = 
			dataFaceArea2[idxSrc]; 
	    }
	  }
	}
	else if(iDir==1)  // Y - Z,X
	{
	  if(iFace==0)
	  {
	    for(int vz=0; vz<bsx; vz++)
	    {
	      idxSrc = MT_GXYZ(bsx,bsx2, 0,bsx-1,vz);
	      idxDst = MT_GXYZ(hsxAct,hsx2Act, SIMDW,0,1+vz);
	      dataSrc = data2 + idxSrc;

	      UT_ASSERT2(idxSrc>=0&&idxSrc+bsx-1<nVoxelsInBlock);
	      UT_ASSERT2(idxDst>=0&&idxDst+bsx-1<hsxActLen);

	      for(int vx=0; vx<bsx; vx+=SIMDW)
	      {
		X.load_a(dataSrc+vx);
		X *= H;
		X.store_a(readDataX+idxDst+vx);
	        #ifndef RELAX_BLOCK_GAUSS_SEIDEL
		if(opmode & LAPL_RELAX ) X.store_a(writeDataX+idxDst+vx);
		#endif
	      }
	    }
	  }
	  else
	  {
	    for(int vz=0; vz<bsx; vz++)
	    {
	      idxSrc = MT_GXYZ(bsx,bsx2, 0,0,vz);
	      idxDst = MT_GXYZ(hsxAct,hsx2Act, SIMDW,hsx-1,1+vz);
	      dataSrc = data2 + idxSrc;

	      UT_ASSERT2(idxSrc>=0&&idxSrc+bsx-1<nVoxelsInBlock);
	      UT_ASSERT2(idxDst>=0&&idxDst+bsx-1<hsxActLen);

	      for(int vx=0; vx<bsx; vx+=SIMDW)
	      {
		X.load_a(dataSrc+vx);
		X *= H;
		X.store_a(readDataX+idxDst+vx);
	        #ifndef RELAX_BLOCK_GAUSS_SEIDEL
		if(opmode & LAPL_RELAX ) X.store_a(writeDataX+idxDst+vx);
		#endif
		if(blockHasFaceCoeffs) 
		{
		  W.load_a( dataFaceArea2+idxSrc+vx );
		  W.store_a( tmpDataW[3]+MT_GXYZ(bsx,bsx2, vx, bsx-1,vz) );
		}
	      }
	    }
	  }
	}
	else if(iDir==2) // Z - Y,X
	{
	  if(iFace==0)
	  {
	    for(int vy=0; vy<bsx; vy++)
	    {
	      idxSrc = MT_GXYZ(bsx,bsx2,0,vy,bsx-1);
	      idxDst = MT_GXYZ(hsxAct,hsx2Act, SIMDW,1+vy,0);
	      dataSrc = data2 + idxSrc;

	      UT_ASSERT2(idxSrc>=0&&idxSrc+bsx-1<nVoxelsInBlock);
	      UT_ASSERT2(idxDst>=0&&idxDst+bsx-1<hsxActLen);

	      for(int vx=0; vx<bsx; vx+=SIMDW)
	      {
		X.load_a(dataSrc+vx);
		X *= H;
		X.store_a(readDataX+idxDst+vx);
	        #ifndef RELAX_BLOCK_GAUSS_SEIDEL
		if(opmode & LAPL_RELAX ) X.store_a(writeDataX+idxDst+vx);
		#endif
	      }
	    }
	  }
	  else
	  {
	    for(int vy=0; vy<bsx; vy++)
	    {
	      idxSrc = MT_GXYZ(bsx,bsx2,0,vy,0);
	      idxDst = MT_GXYZ(hsxAct,hsx2Act, SIMDW,1+vy,hsx-1);
	      dataSrc = data2 + idxSrc;

	      UT_ASSERT2(idxSrc>=0&&idxSrc+bsx-1<nVoxelsInBlock);
	      UT_ASSERT2(idxDst>=0&&idxDst+bsx-1<hsxActLen);

	      for(int vx=0; vx<bsx; vx+=SIMDW)
	      {
		X.load_a(dataSrc+vx);
		X *= H;
		X.store_a(readDataX+idxDst+vx);
	        #ifndef RELAX_BLOCK_GAUSS_SEIDEL
		if(opmode & LAPL_RELAX ) X.store_a(writeDataX+idxDst+vx);
		#endif
		if(blockHasFaceCoeffs) 
		{
		  W.load_a(dataFaceArea2+idxSrc+vx);
		  W.store_a( tmpDataW[5] + MT_GXYZ(bsx,bsx2, vx, vy,bsx-1) );
		}
	      }
	    }
	  }
	}  
      }
      else  // ======== Mixed Resolution =======================
      {
	int vx,vy,vz,vid;
	double sumVal,sumCoeff0;

	if(iDir==0)	
	{
	  if(iFace==0)
	  {
	    for(vz=0; vz<bsx; vz++) 
	    for(vy=0; vy<bsx; vy++) 
	    {
	      vx = 0; 
	      vid = MT_GXYZ(bsx,bsx2,vx,vy,vz);
	      sumVal=sumCoeff0=0;
	      procNeighMixedRes2<bsx,0,0>( level,level2, 1, 
		  		       dataFlg2, data2, dataFaceArea2,
				       vx, vy, vz, vid,
				       sumVal, sumCoeff0 );
	      PSFloat val = sumVal * h;
	      idxDst = MT_GXYZ(hsxAct,hsx2Act, SIMDW-1, 1+vy,1+vz);
	      readDataX[idxDst] = val;
	      #ifndef RELAX_BLOCK_GAUSS_SEIDEL
	      if(opmode & LAPL_RELAX ) writeDataX[idxDst] = val;
	      #endif	     
	      // values from neighbor block already pre-multiplied by procNeighMixedRes2()
	      if(blockHasFaceCoeffs) *(tmpDataW[0] + MT_GXYZ(bsx,bsx2, 0,vy,vz)) = 1.0f;
	    }
	  }
	  else
	  {
	    for(vz=0; vz<bsx; vz++) 
	    for(vy=0; vy<bsx; vy++) 
	    {
	      vx = bsx-1; 
	      vid = MT_GXYZ(bsx,bsx2,vx,vy,vz);
	      sumVal=sumCoeff0=0;	      
	      procNeighMixedRes2<bsx,0,1>( level,level2, 1, 
		  		       dataFlg2, data2, dataFaceArea2,
				       vx, vy, vz, vid,
				       sumVal, sumCoeff0 );
	      PSFloat val = sumVal * h;
	      idxDst = MT_GXYZ(hsxAct,hsx2Act, hsxAct-SIMDW, 1+vy,1+vz);
	      readDataX[idxDst] = val;
	      #ifndef RELAX_BLOCK_GAUSS_SEIDEL
	      if(opmode & LAPL_RELAX ) writeDataX[idxDst] = val;
	      #endif	     
	      if(blockHasFaceCoeffs) *(tmpDataW[1] + MT_GXYZ(bsx,bsx2, bsx-1,vy,vz)) = 1.0f;
	    }
	  }
	}
	else if(iDir==1)
	{
	  if(iFace==0)
	  {
	    for(vz=0; vz<bsx; vz++) 
	    for(vx=0; vx<bsx; vx++) 
	    {
	      vy = 0; 
	      vid = MT_GXYZ(bsx,bsx2,vx,vy,vz);
	      sumVal=sumCoeff0=0;
	      procNeighMixedRes2<bsx,1,0>( level,level2, 1, 
		  		       dataFlg2, data2, dataFaceArea2,
				       vx, vy, vz, vid,
				       sumVal, sumCoeff0 );
	      idxDst = MT_GXYZ(hsxAct,hsx2Act, SIMDW+vx,0,1+vz);
	      PSFloat val = sumVal * h;
	      readDataX[idxDst] = val;
	      #ifndef RELAX_BLOCK_GAUSS_SEIDEL
	      if(opmode & LAPL_RELAX ) writeDataX[idxDst] = val;
	      #endif	     
	      if(blockHasFaceCoeffs) *(tmpDataW[2]+MT_GXYZ(bsx,bsx2, vx,0,vz)) = 1.0f;
	    }
	  }
	  else
	  {
	    for(vz=0; vz<bsx; vz++) 
	    for(vx=0; vx<bsx; vx++) 
	    {
	      vy = bsx-1; 
	      vid = MT_GXYZ(bsx,bsx2,vx,vy,vz);
	      sumVal=sumCoeff0=0;
	      procNeighMixedRes2<bsx,1,1>( level,level2, 1, 
		  		       dataFlg2, data2, dataFaceArea2,
				       vx, vy, vz, vid,
				       sumVal, sumCoeff0 );
	      idxDst = MT_GXYZ(hsxAct,hsx2Act, SIMDW+vx,hsx-1,1+vz);
	      PSFloat val = sumVal * h;
	      readDataX[idxDst] = val;
	      #ifndef RELAX_BLOCK_GAUSS_SEIDEL
	      if(opmode & LAPL_RELAX ) writeDataX[idxDst] = val;
	      #endif	     
	      if(blockHasFaceCoeffs) *(tmpDataW[3]+MT_GXYZ(bsx,bsx2, vx,bsx-1,vz)) = 1.0f;
	    }
	  }
	}
	else if(iDir==2)
	{
	  if(iFace==0)
	  {
	    for(vy=0; vy<bsx; vy++) 
	    for(vx=0; vx<bsx; vx++) 
	    {
	      vz = 0; 
	      vid = MT_GXYZ(bsx,bsx2,vx,vy,vz);
	      sumVal=sumCoeff0=0;
	      procNeighMixedRes2<bsx,2,0>( level,level2, 1, 
		  		       dataFlg2, data2, dataFaceArea2,
				       vx, vy, vz, vid,
				       sumVal, sumCoeff0 );
	      idxDst = MT_GXYZ(hsxAct,hsx2Act, SIMDW+vx,1+vy,0);
	      PSFloat val = sumVal * h;
	      readDataX[idxDst] = val;
	      #ifndef RELAX_BLOCK_GAUSS_SEIDEL
	      if(opmode & LAPL_RELAX ) writeDataX[idxDst] = val;
	      #endif	     
	      if(blockHasFaceCoeffs) *(tmpDataW[4] + MT_GXYZ(bsx,bsx2, vx,vy,0)) = 1.0f;
	    }
	  }
	  else
	  {
	    for(vy=0; vy<bsx; vy++) 
	    for(vx=0; vx<bsx; vx++) 
	    {
	      vz = bsx-1; 
	      vid = MT_GXYZ(bsx,bsx2,vx,vy,vz);
	      sumVal=sumCoeff0=0;
	      procNeighMixedRes2<bsx,2,1>( level,level2, 1, 
		  		       dataFlg2, data2, dataFaceArea2,
				       vx, vy, vz, vid,
				       sumVal, sumCoeff0 );
	      idxDst = MT_GXYZ(hsxAct,hsx2Act, SIMDW+vx,1+vy,hsx-1);
	      PSFloat val = sumVal * h;
	      readDataX[idxDst] = val;
	      #ifndef RELAX_BLOCK_GAUSS_SEIDEL
	      if(opmode & LAPL_RELAX ) writeDataX[idxDst] = val;
	      #endif	     
	      if(blockHasFaceCoeffs) *(tmpDataW[5] + MT_GXYZ(bsx,bsx2, vx,vy,bsx-1)) = 1.0f;
	    }
	  }
	}
      }
    }
  }

  /*if(levelMg==1 && level==1 && bid==144 )
  {
    TRCERR(("DBG\n"));
  }*/

#if 0
    for(int vz=0;vz<bsx;vz++)
    for(int vy=0;vy<bsx;vy++)
    for(int vx=0;vx<bsx;vx++)
    {
      Vec4i vpos(vx,vy,vz,0);
      for(int iFace=0;iFace<=6;iFace++)
      {
	Vec4i vpos2 = vpos + _v4iStencil7[iFace];
	int vid2 = MT_GXYZ(hsxAct,hsx2Act, SIMDW+vpos2[0],1+vpos2[1],1+vpos2[2]);
	float *p = readDataX+vid2;
	if(!( *p==0.0f ))
	//if(!( fabsf(*p)<1e-12 ))
	{
	  Vec4i bpos=sgX->getBlockCoordsById(bid);
	  TRCERR(("L=%d, l=%d, bid=%d bpos=%d,%d,%d, vpos=%d,%d,%d, val=%g\n",levelMg,level,
	      bid,bpos[0],bpos[1],bpos[2],vpos2[0],vpos2[1],vpos2[2],*p));
	  Vec4i pos = sgX->getGridCoords(bid,sgX->getVoxelInBlockIndex(vx,vy,vz));
	  traceCellNeighborhoodInfo(levelMg,level,pos[0],pos[1],pos[2]);
	}
      }
    }
#endif

  #ifdef SKIP_NONFLUID_SIMD_STATS
  skipstats.result();
  if(bi->level==0)
  {
    TRCP(("skipstats (bid=%d, SIMDW=%d) : n=%d, avg=%g\n",bid,SIMDW,(int)skipstats._n,skipstats._avg));  
  }
  #endif

  #ifdef SKIP_NONFLUID_SIMD_STATS
  skipstats.reset();
  #endif

  //
  //
  // opmode LAPL_RELAX
  //
  //

  //
  // Loop over destination block
  //      
  //UT_ASSERT2(!(nIter&1));


  #ifdef ADAPTIVE_BLOCK_LOCAL_RELAX
  float maxResid=0,maxResidPrev=-1,maxResid0=-1;
  VecNpfs maxResid_(0);
  #endif

  #if defined( RELAX_BLOCK_GAUSS_SEIDEL ) && !defined( RELAX_BLOCK_GAUSS_SEIDEL_RB)
  int vx1,vx2, vy1,vy2, vz1,vz2, dx,dy,dz,iSimd1;

  #endif
  
  for(int iIter=0; iIter < ((opmode & LAPL_RELAX) ? nIter : 1); iIter++)
  {
    #ifdef ADAPTIVE_BLOCK_LOCAL_RELAX
    maxResid_=0;
    #endif

    int iSimd;

    const VecNpfs HR(1.f/h);

    const VecNpfs omega(_mgSmOmega),
	        oneMinusOmega = 1.0f - omega;

    #if defined( RELAX_BLOCK_GAUSS_SEIDEL_RB )
    for( int pass_red_black=0; pass_red_black < (opmode & LAPL_RELAX ? 2:1); pass_red_black++)
    {

    iSimd=0;
    for(int vz=0,vid=0;vz<bsx;vz++)
    for(int vy=0;vy<bsx;vy++)
    {
      int vid2 = MT_GXYZ(hsxAct,hsx2Act, SIMDW,1+vy,1+vz);
      float rbMask_[17] = {0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0};

      bool rb = (vy+vz) & 1; 
      if( pass_red_black ) rb = !rb;
      if( options & OPT_REVERSE_ORDER ) rb = !rb;


      VecNf rbMask = VecNf().load( rbMask_ + (rb?1:0) );
      VecNf vzero(0.0f);

      for(int vx=0; vx<bsx; vx+=SIMDW,vid+=SIMDW,vid2+=SIMDW,iSimd++)
      {	
    

    #elif defined( RELAX_BLOCK_GAUSS_SEIDEL )
    {

    if( (opmode & LAPL_RELAX) && (!!(options & OPT_REVERSE_ORDER) == !!(iIter&1) ) )
    {
      vx1=bsx-SIMDW; vx2=-SIMDW; dx=-SIMDW;
      vy1=bsx-1; vy2=-1; dy=-1;
      vz1=bsx-1, vz2=-1; dz=-1;
      iSimd1 = bsx*bsx*bsx/SIMDW-1;
    }
    else
    {
      vx1=0; vx2=bsx; dx=SIMDW;
      vy1=0; vy2=bsx; dy=1;
      vz1=0, vz2=bsx; dz=1;
      iSimd1 = 0;
    }

    iSimd = iSimd1;
    for(int vz=vz1; vz!=vz2; vz+=dz)
    for(int vy=vy1; vy!=vy2; vy+=dy)
    {
      int vid =  MT_GXYZ(bsx,bsx2, vx1,vy,vz),
	  vid2 = MT_GXYZ(hsxAct,hsx2Act, SIMDW+vx1,1+vy,1+vz);
      for(int vx=vx1; vx!=vx2; vx+=dx,vid+=dx,vid2+=dx,iSimd+=dy)
      {
	/*if(tid==0 && bsx==8)
	{
	  TRCP(("RO=%d  %d %d %d  iSimd=%d vid=%d vid2=%d\n",options & OPT_REVERSE_ORDER,vx,vy,vz,iSimd,vid,vid2));
	}*/
    #else
    {
    iSimd=0;
    for(int vz=0,vid=0;vz<bsx;vz++)
    for(int vy=0;vy<bsx;vy++)
    {
      int vid2 = MT_GXYZ(hsxAct,hsx2Act, SIMDW,1+vy,1+vz);

      for(int vx=0; vx<bsx; vx+=SIMDW,vid+=SIMDW,vid2+=SIMDW,iSimd++)
      {
    #endif
	
	//for(int ii=0;ii<2;ii++)
	{

	//IACA_START

	VecNpfs F0,F1,F2,F3,F4,F5,F6,F7,
	        W1,W2,W3,W4,W5,W6,W7,
	        Y,D,S;

	#ifdef DBG_RELAX
	W1=W2=W3=W4=W5=W6=W7 = 0;
	S=0;
        #endif


	if(( blockHasSkippableSIMD && (tmpDataF[iSimd] & (1<<1) ))) 
	{
	  #ifdef SKIP_NONFLUID_SIMD_STATS
	  skipstats.update(1);
	  #endif
	  Y = 0;
	  continue;
	}
	else
	{
	  PSFloat *p = readDataX+vid2;
	  F0.load_a( p );

	  #ifdef SKIP_NONFLUID_SIMD_STATS
	  skipstats.update(0);
	  #endif
	  D.load_a( dataDiag+vid );

	  VecNpfs B;
	  F1.load( p - 1 );
	  F2.load( p + 1 );
	  F3.load_a( p - hsxAct );
	  F4.load_a( p + hsxAct );	  
	  F5.load_a( p - hsx2Act );
	  F6.load_a( p + hsx2Act );
	  
	  if( blockHasFaceCoeffs && tmpDataF[iSimd] )
	  {
	    W1 = vecf2vecpsf( VecNf().load_a( tmpDataW[0]+vid ) ),
	    W2 = vecf2vecpsf( VecNf().load_a( tmpDataW[1]+vid ) ),
	    W3 = vecf2vecpsf( VecNf().load_a( tmpDataW[2]+vid ) ),
	    W4 = vecf2vecpsf( VecNf().load_a( tmpDataW[3]+vid ) ),
	    W5 = vecf2vecpsf( VecNf().load_a( tmpDataW[4]+vid ) ),
	    W6 = vecf2vecpsf( VecNf().load_a( tmpDataW[5]+vid ) );
	    
            #if 0
	    #ifdef DBG_RELAX
	    //#ifdef UT_ASSERT_LEVEL_2
	    //#ifdef UT_ASSERT_LEVEL_3
	    // Check diagonally dominant
	    if( iIter==0
		//&& !tmpBlkErrorStatus[tid]
		)
	    {
              #if 0
	      UT_ASSERT0(!vany( W1+W2+W3+W4+W5+W6 > 1e-3f + 1.0f/D ));
	      #else
	      VecNpfs WSUM = W1+W2+W3+W4+W5+W6;
	      for(int k=0;k<SIMDW;k++)
	      {
		float e = WSUM[k] - 1.0f/D[k];		
		if(e>0.0f )
		{
		  TRCERR(("k=%d/%d: %.10g > %.10g e=%.10g %.10g\n",
			k,SIMDW,
			WSUM[k],1.0f/D[k],e,e/WSUM[k]));
		  Vec4i pos = sg->getGridCoords(bid,vid+k);
		  traceCellNeighborhoodInfo(levelMg,level,pos[0],pos[1],pos[2]);
		}
	      }

	      #endif
	      //tmpBlkErrorStatus[tid] = 1;
	    }
	    #endif
	    #endif

	    S = 0.f;
	    S = mul_add(F1,W1, S);
	    S = mul_add(F2,W2, S);
	    S = mul_add(F3,W3, S);
	    S = mul_add(F4,W4, S);
	    S = mul_add(F5,W5, S);
	    S = mul_add(F6,W6, S);
	  }
	  else 
	  {
	    S = F1 + F2 + F3 + F4 + F5 + F6;
	  }

	  if constexpr (opmode & LAPL_MULTIPLY )
	  {

	    Y = select( D!=0.0f, (1.f/D) * F0 - S, 
				 0.0f );
	    if constexpr( opmode & LAPL_CALC_RESIDUAL )
	    {
	      B.load_a(dataB+vid);
	      Y = B - Y;
	    }
	    
	    if constexpr( opmode & LAPL_WITH_DOTPROD )
	    {
	      VecNpfs F0_ = F0 * HR;
	      if constexpr (SIMDW==2)
	      {
		#ifdef SOLVE_PRESSURE_DOUBLE
		UT_ASSERT0(FALSE);
		#else
		Vec4f Y4 = Y, F0_4 = F0_;
		Y4.cutoff(2); F0_4.cutoff(2);
                *v4dAlpha += vecpsf2vec4d( Y4 ) * vecpsf2vec4d( F0_4 );
		#endif
	      }
	      else if constexpr (SIMDW==4)
	      {
                *v4dAlpha += vecpsf2vec4d( Y ) * vecpsf2vec4d( F0_ );
	      }
	      else
	      {
		UT_ASSERT2(SIMDW==8)
                *v4dAlpha +=
		     vecpsf2vec4d( Y.get_low() ) * vecpsf2vec4d( F0_.get_low() ) +
		     vecpsf2vec4d( Y.get_high() ) * vecpsf2vec4d( F0_.get_high() );
	      }
	    }
	    UT_ASSERT2(visfinite(Y));
	    vstream(dataDstY+vid, Y);
	    continue;
	  }
	  else
	  {
	    B.load_a(dataB+vid);
	    Y = omega*(B+S)*D  + oneMinusOmega * F0;

	    #ifdef ADAPTIVE_BLOCK_LOCAL_RELAX
	    //UT_ASSERT0( !vany(D==0.0f) );
	    //maxResid_ = max( maxResid_, abs( (Y-F0)*select( D!=0.0f, 1.0f/D, 0.0f ) ) );
	    maxResid_ = max( maxResid_, abs( (Y-F0) ) );
	    //UT_ASSERT0(visfinite(maxResid_));
	    #endif
	  }
	     
	  UT_ASSERT2(visfinite(Y));

	  //IACA_END	
	}

        #ifdef RELAX_BLOCK_GAUSS_SEIDEL_RB
	if constexpr( opmode & LAPL_RELAX )
	{
	  Y = select( rbMask > vzero, Y, F0 );
	  Y.store_a( writeDataX+vid2 );

	  #if 1
	  if(levelMg==0 && iIter<=1 && vy<=3 && vz<=3)
	  {
	    char chbufM[1024]; vprint( chbufM, rbMask );
	    char chbufY[1024]; vprint( chbufY, Y );
	    char chbufF0[1024]; vprint( chbufF0, F0 );
	    TRCP(("L=%d i=%d/%d rb=%d vx=%d,vy=%d,vz=%d rbMask=%s F0=%s Y=%s\n",
		  levelMg,iIter,nIter,pass_red_black,vx,vy,vz,chbufM,chbufF0,chbufY));
	  }
	  #endif
	}
	else
	#endif
	{

	// Store result	(last iteration writes to final dest.)
	if(likely(iIter<nIter-1))
	{
	  Y.store_a( writeDataX+vid2 );
	}
	else
	{
	  Y *= HR;

	  UT_ASSERT2(visfinite(Y));

	  #ifdef DBG_RELAX
	  if( 
	      !tmpBlkErrorStatus[tid] && 
	      ( !visfinite(Y) || vany(abs(Y)>1e7) ) )
	  {
	    {
	      char chbuf[100];
	      TRCERR(("Y=%s\n",vprint(chbuf,Y)));
	    }
	    int k0 = -1;
	    for(int k=0;k<SIMDW;k++)
	    {
	      if(abs(Y[k])>1e7 || !isfinite(Y[k]) ) 
	      {
		k0=k;
	        break;
	      }
	    }
	    UT_ASSERT0(k0>=0);

	    Vec4i pos = sg->getGridCoords(bid,vid) + Vec4i(k0,0,0,0);
	    traceCellNeighborhoodInfo(levelMg,level,pos[0],pos[1],pos[2]);

	    VecNpfs B = VecNpfs().load(dataB+vid);
	    VecNpfs WSUM = W1+W2+W3+W4+W5+W6;
	    PSFloat *dataX0 = residDataX; 
            VecNpfs X0 = VecNpfs().load_a(dataX0+vid);

	    #if 0
	    char chbuf[100];
	    TRCERR(("L=%d l=%d bid=%d vid=%d Y=%s kSIMD=%d\n",levelMg,level,bid,vid,
		  vprint(chbuf,Y),kSIMD));
	    { char chbuf[100]; TRCERR(("B=%s\n",chbuf,vprint(chbuf,B))); }
	    { char chbuf[100]; TRCERR(("WSUM=%s\n",chbuf,vprint(chbuf,WSUM))); }
	    { char chbuf[100]; TRCERR(("W1=%s\n",chbuf,vprint(chbuf,W1))); }
	    { char chbuf[100]; TRCERR(("W2=%s\n",chbuf,vprint(chbuf,W2))); }
	    { char chbuf[100]; TRCERR(("W3=%s\n",chbuf,vprint(chbuf,W3))); }
	    { char chbuf[100]; TRCERR(("W4=%s\n",chbuf,vprint(chbuf,W4))); }
	    { char chbuf[100]; TRCERR(("W5=%s\n",chbuf,vprint(chbuf,W5))); }
	    { char chbuf[100]; TRCERR(("W6=%s\n",chbuf,vprint(chbuf,W6))); }
	    { char chbuf[100]; TRCERR(("W7=%s\n",chbuf,vprint(chbuf,W7))); }
	    #endif

	    TRCP(("dir=%d L=%d l=%d bid=%d iBlkIter=%d vid=%d op=%d opt=%d\n",
		  options & OPT_REVERSE_ORDER,levelMg,level,bid,iIter,vid,opmode,options));

            TRCP(("Y=%g X0=%g D=%g B=%g S=%g flags=%d\n", Y[k0],X0[k0],D[k0],B[k0],S[k0],dataFlags[vid+k0]));

	    TRCP(("X = %g, %g %g, %g %g, %g %g\n",F0[k0],F1[k0],F2[k0],F3[k0],F4[k0],F5[k0],F6[k0]));

	    for(int k=0;k<SIMDW;k++)
	    {
	      float e = WSUM[k] - 1.0f/D[k];		
	      if(e>0.0f )
	      {
		TRCERR(("k=%d/%d %g > %g e=%g %g\n",WSUM[k],1.0f/D[k],e,e/WSUM[k],k,SIMDW));
	        traceCellNeighborhoodInfo(levelMg,level,pos[0]+k,pos[1],pos[2]);
	      }
	    }
	    tmpBlkErrorStatus[tid] = 1;
	    //UT_ASSERT0(FALSE);
	  }
	  #endif

	  Y.store_a( dataDstY+vid );

	  #if 0
	  for(int l=0;l<4;l++)
	  {
	    UT_ASSERT0( std::isfinite(dataDstY[vid+l]) && fabs(dataDstY[vid+l])<1e9   );
	  }
	  #endif
	}
	}

	//IACA_END
	} // SIMD-line
      }
    } // Cells

    } // red-black

    #ifdef RELAX_BLOCK_GAUSS_SEIDEL_RB
    for(int vz=0; vz<bsx; vz++)
    for(int vy=0; vy<bsx; vy++)
    for(int vx=0; vx<bsx; vx+=SIMDW)
    {
      int vid = MT_GXYZ(bsx,bsx2, vx, vy,vz),
	  vid2 = MT_GXYZ(hsxAct,hsx2Act, SIMDW+vx,1+vy,1+vz);
      VecNf Y = VecNf().load( writeDataX + vid2 );
      Y *= HR;
      Y.store_a( dataDstY+vid );
    }
    #endif

    #ifdef ADAPTIVE_BLOCK_LOCAL_RELAX    
    if( (opmode & LAPL_RELAX) && _psConf.ps_mg_ar_testval>0)
    {
      maxResid = horizontal_max(maxResid_);
      #if 0
      if(maxResid>residThetaAct)
      {
        TRCP(("%d maxResid=%g %g\n",iIter,maxResid,residThetaAct));
      }
      #endif

      if((iIter>0) && iIter<nIter-2)
      {
	#if 1
	if( maxResid <  _psConf.ps_mg_ar_testval*maxResid0 )
	{
	  //TRCP(("%d/%d -> break\n",iIter,nIter));
	  iIter = nIter-2;
	}
	#else
	if(!(maxResidPrev<0))
	{
	  float residReductionRate = maxResid / maxResidPrev;
	  //TRCP(("%d -> %g\n",iIter,residReductionRate));
	  if(residReductionRate > _psConf.ps_mg_ar_testval /*.85*/) 
	  {
	    iIter = nIter-2;
	  }
	}
	#endif
      }
      maxResidPrev = maxResid;
      if(iIter==0) maxResid0 = maxResid;
    }
    #endif

    #ifndef RELAX_BLOCK_GAUSS_SEIDEL
    UT_SWAP_PTR_(float, readDataX, writeDataX );
    #endif
  } // block local iterations

  if( (opmode & LAPL_RELAX) && (options & OPT_CALC_PER_BLOCK_RESID))
  {
    VecNpfs M(0.0f);

    PSFloat *dataX0 = residDataX; 
    for(int vz=0; vz<bsx; vz++)
    for(int vy=0; vy<bsx; vy++)
    for(int vx=0; vx<bsx; vx+=SIMDW)
    {
      int vid = MT_GXYZ(bsx,bsx2, vx, vy,vz);
          //hid = MT_GXYZ(hsxAct,hsx2Act, 1+vx+1,1+vy,1+vz);
      VecNpfs X0 = VecNpfs().load_a(dataX0+vid);
              X = VecNpfs().load(dataDstY+vid);

      #ifdef DISCARD_BLOCK_SMALL_CHANGE
      X0.store_a(&tmpBlkDataXBuf[tid][1][vid]);
      #endif

      X = abs(X-X0);      

      #if 1
      VecNpfs D = VecNpfs().load_a(dataDiag+vid);
      X = select( D>1e-5f, X/D, 0.0f );
      #endif

      UT_ASSERT2(visfinite(X));
      
      M = max( M, X );

      #if 0
      #ifdef DBG_RELAX
      if( !visfinite(X) || vany(abs(X)>1e7) )
      {
	char chbuf1[100],chbuf2[100];
	TRCERR(("L=%d l=%d bid=%d vid=%d D=%s X=%s\n",levelMg,level,bid,vid,
	      vprint(chbuf1,D),vprint(chbuf2,X)));
      }
      #endif
      #endif

      X.store_a(dataX0+vid);
    }

    // TODO: subtract from smoothed residual to get HF-residual
    #if 0

    float maxResid = 0;
    for(int vz=0; vz<bsx; vz+=2)
    for(int vy=0; vy<bsx; vy+=2)
    for(int vx=0; vx<bsx; vx+=2)
    {
      float avg=0;
      for(int dz=0;dz<2;dz++)
      for(int dy=0;dy<2;dy++)
      for(int dx=0;dx<2;dx++)
      {
	int vid2 = MT_GXYZ(bsx,bsx2, vx+dx,vy+dy,vz+dz);
	avg += residDataX[vid2];	
      }
      avg *= 1.0f/8.0f;
      for(int dz=0;dz<2;dz++)
      for(int dy=0;dy<2;dy++)
      for(int dx=0;dx<2;dx++)
      {
	int vid2 = MT_GXYZ(bsx,bsx2, vx+dx,vy+dy,vz+dz);
	residDataX[vid2] = fabsf(residDataX[vid2]-avg);
	maxResid = std::max(maxResid,residDataX[vid2]);
      }
    }

    if(pMaxScaledResid) *pMaxScaledResid = maxResid;

    #else

    // Determine maximum
    float maxResid = horizontal_max(M);
    //_mgPerBlockResid[bid] = maxResid;
    if(pMaxScaledResid) *pMaxScaledResid = maxResid;
    #endif

  }

  #if 0
  #ifdef SKIP_NONFLUID_SIMD_STATS
  skipstats.result();
  if(bi->level==0)
  {
    TRCP(("skipstats_ (bid=%d) : n=%d, avg=%g\n\n",bid,(int)skipstats._n,skipstats._avg));  
  }
  #endif
  #endif

}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template< int opmode /* LAPL_RELAX, LAPL_MULTIPLY */ >
void MultiresSparseGrid::processBlockLaplacian(
		    int tid,	// Thread number (-1=init buffers)
		    unsigned options,
		    int bid, int bx, int by, int bz,
		    BlockInfo *bi,
		    int levelMg,
		    int chanX,
		    int chanB,
		    int chanDstY, 
		    int nIter,
		    Vec4d *v4dAlpha,
		    PSFloat *pMaxScaledResid
		  )
{
  #ifdef SKIP_NONFLUID_SIMD_STATS
  MT::Stats<> skipstats;
  #endif

  // 
  // Initialize per-thread buffers
  //
  
  if(tid<0)
  {
    const int SIMDWMAX = 8;   
    int bsxAct = ALIGN( _sg0->bsx()+2*SIMDWMAX, SIMDWMAX ),
	bsyAct = _sg0->bsx()+2,
	bszAct = _sg0->bsx()+2,
        bufLen = bsxAct*bsyAct*bszAct;

    for(int i=0;i<ARRAY_LENGTH(tmpBlkErrorStatus);i++) tmpBlkErrorStatus[i]=0;

    if( !tmpBlkDataXBuf[0][0] || 
	bufLen != tmpBlkMaxLen || 
	nMaxThreads != tmpBlkMaxThreads  )
    {
      tmpBlkMaxLen = bufLen;
      tmpBlkMaxThreads = nMaxThreads;
      for(int i=0;i<nMaxThreads;i++)
      {
	UT_ASSERT0(i<MI_MAX_THREADS);
	// 
	// Allocate temp. block to be used with ghost cells
	// TODO: NUMA-AWARE
	//
	{
	  if(i==0)
	  {
	    TRC(("%s: Allocating memory for %dx%dx%d temp. halo block\n",
		  UT_FUNCNAME,
		  bsxAct,bsyAct,bszAct));
	  }
	  for(int k=0;k<3;k++)
	  {
	    PSFloat *ptr=NULL;
	    size_t szNeeded = tmpBlkMaxLen*sizeof(*ptr);
	    FREEMEM_ALIGNED( tmpBlkDataXBuf[i][k] );
	    ALLOCARR_ALIGNED_( ptr, PSFloat, tmpBlkMaxLen );
	    memset(ptr,0,szNeeded);
	    tmpBlkDataXBuf[i][k] = ptr;
	  }
	  for(int k=0;k<6;k++)
	  {
	    float *ptr=NULL;
	    size_t szNeeded = tmpBlkMaxLen*sizeof(*ptr);
	    FREEMEM_ALIGNED( tmpBlkDataWBuf[i][k] );
	    ALLOCARR_ALIGNED_( ptr, float, tmpBlkMaxLen );
	    memset(ptr,0,szNeeded);
	    tmpBlkDataWBuf[i][k] = ptr;
	  }

	  size_t szNeeded = tmpBlkMaxLen/SIMDWMAX+1;
	  FREEMEM_ALIGNED( tmpBlkDataFBuf[i] );
	  ALLOCARR_ALIGNED_( tmpBlkDataFBuf[i], uint8_t, szNeeded );
	  memset(tmpBlkDataFBuf[i],0,szNeeded);
	}
      }
    }
    return;
  }

  //if constexpr (!(opmode & LAPL_INIT_THREAD_LOCALS))
  {

  const int bsx = getFlagsChannel(CH_CELL_FLAGS,bi->level,levelMg)->bsx();

  switch(bsx)
  {
    #ifdef UT_ASSERT_LEVEL_2
    case 32: 
          UT_ASSERT0(FALSE);
	  break;   	
    #endif

    case 16: IprocessBlockLaplacian
	     #ifdef SOLVE_PRESSURE_DOUBLE
	     <opmode, 16, 4, Vec4d, Vec4f, Vec4ui >
             #else
	     <opmode, 16, 8, Vec8f, Vec8f, Vec8ui >
	     #endif
	     ( tid, options, 
	       bid, bx,by,bz, bi, levelMg, chanX, chanB, chanDstY, nIter, v4dAlpha,
	       pMaxScaledResid );  
	  break;   

    case 8: IprocessBlockLaplacian
	    #ifdef SOLVE_PRESSURE_DOUBLE
	    < opmode, 8, 4, Vec4d, Vec4f, Vec4ui >
	    #else
	    < opmode, 8, 8, Vec8f, Vec8f, Vec8ui >
	    #endif
	    ( tid, options,
		 bid, bx,by,bz, bi, levelMg, chanX, chanB, chanDstY, nIter, v4dAlpha,  
		 pMaxScaledResid );  
	  break;   

    case 4: IprocessBlockLaplacian
	    #ifdef SOLVE_PRESSURE_DOUBLE
	    < opmode, 4, 4, Vec4d, Vec4f, Vec4ui >
	    #else
	    < opmode, 4, 4, Vec4f, Vec4f, Vec4ui >
	    #endif
	      ( tid, options,
		bid, bx,by,bz, bi, levelMg, chanX, chanB, chanDstY, nIter, v4dAlpha,
		pMaxScaledResid );  
	  break;
#if 1

    case 2: IprocessBlockLaplacian
	    #ifdef SOLVE_PRESSURE_DOUBLE
	    < opmode, 2, 2, Vec4d, Vec4f, Vec4ui >
	    #else
	    < opmode, 2, 2, Vec2f, Vec2f, Vec2ui >
	    #endif
	      ( tid, options,
		bid, bx,by,bz, bi, levelMg, chanX, chanB, chanDstY, nIter, v4dAlpha,
		pMaxScaledResid );  
	  break;
#endif

    default:
	  UT_ASSERT0(FALSE);
	  break;
  }

  }
}

/*-----------------------------------------------------------------------*/
/* 		 							 */
/*-----------------------------------------------------------------------*/
void MultiresSparseGrid::multiplyLaplacianMatrixOpt(  
			  unsigned options, 
			  int levelMg,
			  int chanX,
			  int chanB,
			  float kap,	// optional diffusion coeff
			  int chanDstY,
			  long double* pDotprod,  // OUT
			  std::vector<int> *blockList
			  )
{
  int doBoundZoneOnly = options & OPT_BOUNDZONE_ONLY,
      doCalcResidual = options & OPT_CALC_RESIDUAL;

  long double alpha=0;

  //
  // special case: dense (coarse) grid implementation
  //
  if(levelMg>=_nLevels)  
  {
    UT_ASSERT2(nBlocks(levelMg)==1);
    SparseGrid<CellFlags> *sg = getFlagsChannel(CH_CELL_FLAGS,levelMg,levelMg);
    SparseGrid<PSFloat> *sgX = getPSFloatChannel(chanX,levelMg,levelMg),
		        *sgB = getPSFloatChannel(chanB,levelMg,levelMg),
		        *sgD = getPSFloatChannel(CH_DIAGONAL,levelMg,levelMg),
		        *sgDstY = getPSFloatChannel(chanDstY,levelMg,levelMg);
    double h = (1<<levelMg)*_dx0;

    CellFlags *dataFlags = sg->getDataPtrGen_(0,0,0);
    PSFloat *dataX = sgX->getDataPtrGen_(0,0,0),
	    *dataB = sgB ? sgB->getDataPtrGen_(0,0,0) : NULL,
	    *dataD = sgD ? sgD->getDataPtrGen_(0,0,0) : NULL,
	    *dataDstY = sgDstY->getDataPtrGen_(0,1,0,1);

    float *dataFaceArea[3];
    for(int k=0;k<3;k++) 
      dataFaceArea[k] = getFaceAreaChannel(k,levelMg,levelMg)->getDataPtrGen_();

    int sx=sg->sx(), sy=sg->sy(), sz=sg->sz();
    int yStrid = sg->yStride(), 
	zStrid = sg->zStride();
      
    UT_ASSERT0((LongInt)sx*(LongInt)sy*(LongInt)sz<2000000000);

    const int SIMDW=4;
    
    UT_ASSERT0(!pDotprod);

    UT_ASSERT0((sg->sxAct() & (SIMDW-1))==0);

    bool doBypassCache = false; // levelMg < _mgLevels-3;  TODO

    // Loop over all cells in block
    ThrRunParallel<int>( sz, nullptr,
      [&](int &tls, int tid, int z) 
      {
	ImultiplyLaplacianMatrixProcessSliceDenseLevel( 0, 0,
	    z,
	    false, false, doBoundZoneOnly, doCalcResidual,
	    sx,sy,sz, yStrid,zStrid,
	    h, _mgSmOmega,
	    sg, 
	    dataFlags, dataX, dataB, dataD,
	    dataFaceArea,
	    dataDstY,
	    doBypassCache,
	    false
	    );
      }
      //,nullptr ,0,/*doSerial=*/true
    );

    return;
  } // dense grid implementation


  // Initialze per-thread data
  processBlockLaplacian<LAPL_MULTIPLY|LAPL_INIT_THREAD_LOCALS>
    (-1,0,0,0,0,0,NULL,0,0,0,0,0,NULL);  

  // 
  // Main loop over all blocks
  //
  typedef struct
  {
    Vec4d alpha;
  }
  ThreadLocals;

  int nActBlocks = blockList ? blockList->size() : _nBlocks;

  ThrRunParallel<ThreadLocals>( nActBlocks, 
    [&]( ThreadLocals &tls, int tid )  
    {
      tls.alpha = 0;
    },

    [&](ThreadLocals &tls, int tid, int bid) 
    {
      if(blockList) bid = (*blockList)[bid];
      UT_ASSERT2(bid>=0&&bid<_nBlocks);

      BlockInfo *bi = getBlockInfo(bid,levelMg);

      if((doBoundZoneOnly) && !(bi->flags & BLK_BOUNDARY_ZONE)) 
      {
	return;
      }

      if (bi->flags & (BLK_NO_FLUID | BLK_FIXED))
      {
        #ifdef UT_ASSERT2
	if( ( options & OPT_CALC_RESIDUAL ) && (bi->flags & BLK_FIXED ) )
	{
	  auto sg = getPSFloatChannel( chanDstY, bi->level, levelMg ); 
	  PSFloat *dataY = sg->getBlockDataPtr(bid,bi->level);
	  for(int i=0;i<sg->nVoxelsInBlock();i++) 
	  {
	    UT_ASSERT0(dataY[i]==0.0f);
	  }
	}
	#endif
	return; 
      }

      int bx,by,bz;
      _sg0->getBlockCoordsById( bid, bx, by, bz );

      if( options & OPT_CALC_RESIDUAL )
      {
        UT_ASSERT2(!pDotprod);
        processBlockLaplacian< LAPL_MULTIPLY | LAPL_CALC_RESIDUAL >( 
		  tid, options, bid, bx, by, bz, bi, levelMg, 
		  chanX, chanB, chanDstY, 
		  _mgSmBlockIter, NULL );
      }
      else
      {
        processBlockLaplacian< LAPL_MULTIPLY | LAPL_WITH_DOTPROD >( 
		  tid, options, bid, bx, by, bz, bi, levelMg, 
		  chanX, chanB, chanDstY, 
		  _mgSmBlockIter, &tls.alpha );

      }
    },

    [&]( ThreadLocals &tls, int tid ) // reduce
    {
      alpha += horizontal_add( tls.alpha );
    }
    
    //,0,/*doSerial=*/true 
  );

  if(pDotprod) *pDotprod = alpha;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::ImultiplyLaplacianMatrixDenseLevel(  
			  unsigned options, 
			  int levelMg,
			  int chanX0,
			  int chanB,
			  int chanTmp,
			  int nJacobiIter,
			  int chanDstY0
			  )

{  
  UT_ASSERT0(levelMg >= _nLevels );  

  int doJacobi = options & OPT_JACOBI,
      doCalcResidual = options & OPT_CALC_RESIDUAL,
      doBoundZoneOnly = options & OPT_BOUNDZONE_ONLY,
      doReverseOrder = options & OPT_REVERSE_ORDER;

  bool doBypassCache = false; //levelMg < _mgLevels-3;   TODO
    
  if(!doJacobi)
  {
    nJacobiIter = 1;
    UT_ASSERT0(chanDstY0!=CH_NULL);
  }
  else
  {
    UT_ASSERT0(chanDstY0==CH_NULL);
  }

  volatile LONG chanX = chanX0,
                chanDstY = chanDstY0;

  if(doJacobi)
  {
    #ifdef MG_BLOCK_GAUSS_SEIDEL_DENSELEVELS
    chanX = chanDstY = chanX0;
    #else
    chanX = chanTmp;
    chanDstY = chanX0;
    if(doBoundZoneOnly)
    {
      copyChannel( chanX0, chanTmp, TRUE, -1, levelMg);
    }
    #endif
  }

  UT_ASSERT2(nBlocks(levelMg)==1);
  int level = levelMg;
  SparseGrid<CellFlags> *sg = getFlagsChannel(CH_CELL_FLAGS,levelMg,levelMg);
  SparseGrid<PSFloat> *sgB = getPSFloatChannel(chanB,levelMg,levelMg),
                      *sgD = getPSFloatChannel(CH_DIAGONAL,levelMg,levelMg);
  CellFlags *dataFlags = sg->getDataPtrGen_(0,0,0);
  PSFloat *dataB = sgB ? sgB->getDataPtrGen_(0,0,0) : NULL,
          *dataD = sgD ? sgD->getDataPtrGen_(0,0,0) : NULL;

  int sx=sg->sx(), sy=sg->sy(), sz=sg->sz();
  int yStrid = sg->yStride(), 
      zStrid = sg->zStride();


  //TRCP(("nMaxThreads=%d\n",nMaxThreads));

 // 
 // Determine 'optimal' number of threads. 
 // For 'small' problem sizes it is beneficial to use fewer threads
 // to minimize fine grained synchronization overhead (barrier)
 // TODO: make dependence on problem size less crude

#if 0
  int nThreads = nMaxThreads;
#else
  int nThreads = 
    sx*sy*sz >= 128*128*128 ? nMaxThreads/2 :   
			      nMaxThreads/4 ;

  if(nThreads<4) nThreads = MIN(nMaxThreads,4);

  nThreads = MIN(nThreads,sg->sz());

  if(sx*sy*sz<1500) nThreads=1;   // Assuming single-core performance ~= 300M/sec 
#endif


  if(_psTraceDetail)
  {
    TRC(("%s: %d threads used @ L%d (%dx%dx%d)\n",UT_FUNCNAME,nThreads, levelMg,
	  sx,sy,sz));
  }

#if 0
  /////////////
xxxx   int nThreads = MAX(nMaxThreads/4,4); // TODO: optimze dependent on problem size
   - use more threads for "large" dens grids  e.g. 1024^3@gsx=16 leads to
     first dense level with 256^3 resolution !
////////////
#endif


  ThrSpinBarrier iterationBarrier(nThreads);
  volatile LONG iTask;

  double h = (1<<level)*_dx0,
	 h2r = 1./(h*h);

  UT_ASSERT0((LongInt)sx*(LongInt)sy*(LongInt)sz<2000000000);
  USE(h2r);

  const int SIMDW=4;
  
  UT_ASSERT0((sg->sxAct() & (SIMDW-1))==0);

  PSFloat omegaAct = _mgSmOmega;

  #ifdef MSBG_PERF_STATS
  UtPerfStats perfstats[THR_MAX_THREADS] = {};
  UtTimer tm;
  TIMER_START(&tm);
  #endif

  //TRCP(("%s level=%d jacobi=%d doBoundZoneOnly=%d\n",UT_FUNCNAME,level,doJacobi,doBoundZoneOnly));

  #ifdef NDEBUG  // gdb mingw64 hangs here when runnung multiple threads
  #pragma omp parallel num_threads(nThreads)
  #endif
  {
    int tid = omp_get_thread_num();
    USE(tid);

    //UT_ASSERT0( omp_get_num_threads() == nThreads );      
    //
    
    for(int iIter=0;iIter<nJacobiIter;iIter++)
    {
      int pass_red_black=0;
      #ifdef MG_BLOCK_GAUSS_SEIDEL_DENSELEVELS
      for( pass_red_black=0; pass_red_black<2; pass_red_black++)
      #endif
      {

      #if 1
      #ifdef NDEBUG
      if( iterationBarrier.enter( tid, nThreads ) == 0 ) 
      #endif
      {
	#ifdef MG_SCHEDULED_JACOBI_DENSELEVELS
	UT_ASSERT0( (nJacobiIter & 1) == 0);
	
	omegaAct = levelMg == _mgLevels-1 ? 
	  		_mgSmOmega : 
	  		(iIter & 1) == !!(options & OPT_REVERSE_ORDER) ? 
			  _mgSmOmegaSched1 : 
			  _mgSmOmegaSched2;

	/*TRCP(("iIter=%d/%d, reverseOrde=%d omegaAct=%g\n",iIter,nJacobiIter,
	      !!(options & OPT_REVERSE_ORDER), omegaAct));*/
        #endif

	#ifndef MG_BLOCK_GAUSS_SEIDEL_DENSELEVELS
        if(doJacobi) std::swap( chanX, chanDstY );      
	#endif


	iTask = 0;
	#ifdef NDEBUG
	iterationBarrier.release();
	#endif
      }
      #else
      #pragma omp single
      {
        if(doJacobi) std::swap( chanX, chanDstY );      
	iTask = 0;
      }
      #endif
      // 
      // Continue here after thread barrier
      //      
      SparseGrid<PSFloat> 
	*sgX = getPSFloatChannel(chanX,levelMg,levelMg),
	*sgDstY = getPSFloatChannel(chanDstY,levelMg,levelMg);

      PSFloat 
	*dataX = sgX->getDataPtrGen_(0,0,0),
	*dataDstY = sgDstY->getDataPtrGen_(0,1,0);

      float *dataFaceArea[3];

      for( int k=0;k<3; k++ ) 
	dataFaceArea[k] = getFaceAreaChannel(k,levelMg,levelMg)->getDataPtrGen_(0);

      int z;
      while( (z = InterlockedIncrement((LONG volatile *)&iTask)-1) < sz)	  
      {
	#ifdef MG_BLOCK_GAUSS_SEIDEL_DENSELEVELS
	{
	  bool doSkip = (z & 1) == pass_red_black;
	  if( options & OPT_REVERSE_ORDER) doSkip = !doSkip;
	  if(doSkip) continue;
	}
	#endif

	ImultiplyLaplacianMatrixProcessSliceDenseLevel( iIter, pass_red_black,
	    z, 
	    doJacobi, doReverseOrder, doBoundZoneOnly, doCalcResidual,
	    sx,sy,sz, yStrid,zStrid,
	    h, omegaAct,
	    sg, 
	    dataFlags, dataX, dataB, dataD,
	    dataFaceArea,
	    dataDstY,
	    doBypassCache,
	    false
	    #ifdef MSBG_PERF_STATS
	    ,&perfstats[tid]
	    #endif
	    );
	    					      
      }  // z
      } // red-black
    }  // jacobi iterations
  }  // OMP parallel

  #ifdef MSBG_PERF_STATS
  TIMER_STOP(&tm);
  UtPerfStats pfs={};
  for(int i=0;i<THR_MAX_THREADS;i++) 
  {
    pfs.bytes += perfstats[i].bytes;
    pfs.flops += perfstats[i].flops;
  }

  double t_sec = ((double)TIMER_DIFF_MS(&tm)/1000.);
  TRCP(("%s: CPU=%.2f iter=%d res=%dx%dx%d footprint=%.0f bytes"
	" GB/s=%.2f  GFLOPs/sec=%.2f\n",UT_FUNCNAME,
	t_sec,nJacobiIter,sx,sy,sz, (double)pfs.bytes/(double)nJacobiIter,
	(pfs.bytes / (double)ONE_GB) / t_sec,
	(pfs.flops / 1000000000.) / t_sec));
  #endif

  if(doJacobi)
  {
    if(chanDstY != chanX0)
    {
      copyChannel( chanDstY, chanX0, TRUE,-1,levelMg);
    }
  }

}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
LongInt MultiresSparseGrid::relaxBlockList(  
    			  std::vector<int> *blocks,
			  int levelMg,
			  int iIter,
			  int chanX,  // source channel
			  int chanB,  // right hand side
			  unsigned options,
			  int chanDstY // destination channel				
			  )
{
  // Initialize per-thread buffers
  processBlockLaplacian<LAPL_RELAX>(-1,0,0,0,0,0,NULL,0,0,0,0,0,NULL); 
  LongInt nTotalCells = 0;

  #ifdef RELAX_BLOCKS_RED_BLACK
  for( int pass_red_black=0; pass_red_black<2; pass_red_black++)
  #endif
  {
    typedef struct
    {
      LongInt nTotalCells;
    } 
    ThreadLocals;

    ThrRunParallel<ThreadLocals> ( blocks->size(), 
      [&]( ThreadLocals &tls, int tid )  
      {
	tls.nTotalCells = 0;
      },

      [&](ThreadLocals &tls, int tid, int ibid) 
      {
	int bid = (*blocks)[ibid];
	BlockInfo *bi = getBlockInfo(bid,levelMg);
	UT_ASSERT2(!( bi->flags & BLK_NO_FLUID ));

	int bx,by,bz;
	_sg0->getBlockCoordsById( bid, bx, by, bz );

	// TODO: prepare blocklist accordingly
	if( bi->flags & BLK_FIXED ) return;

	#ifdef RELAX_BLOCKS_RED_BLACK
	{
	  bool doSkip = ((bx+by+bz) & 1) == pass_red_black;
	  if( options & OPT_REVERSE_ORDER) doSkip = !doSkip;
	  if(doSkip) return;
	}
	#endif

	processBlockLaplacian<LAPL_RELAX>( 
	    	    tid, options, bid, bx, by, bz, bi, levelMg, 
		    chanX, chanB, chanDstY, 
		    _mgSmBlockIter, NULL );

        tls.nTotalCells += getFlagsChannel0(bi->level,levelMg)->nVoxelsInBlock();
      },

      [&]( ThreadLocals &tls, int tid )  // Reduce locals
      {
	nTotalCells += tls.nTotalCells;
      }

    //,0,/*doSerial=*/true
    
    ); // ThrRunParallel
  }

  return nTotalCells;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::relax(   int boundaryZoneOnly, 
    				  int reverseOrder,
				   int levelMg,     					    
				   int chanX,
				   int chanB,
				   int chanTmp,
				   int nIter,
				   unsigned options
					     )
{
  //nIter = 200;  
  //if(levelMg==3&&boundaryZoneOnly) nIter=100000;
#if 0
  if(levelMg==2 /*&& boundaryZoneOnly*/)  // XXX ZZZZ
  {    
    nIter=100;
    TRCP(("*** Testing relax() L=%d, bs=%d -> nIter=%d\n",
	  levelMg,boundaryZoneOnly,nIter));
  }
#endif
//  UT_ASSERT0((nIter&1)==0); // Even number avoids extra copy 

  if(!nIter) return;
//  if(levelMg>=_nLevels && boundaryZoneOnly) nIter*=1.5;

//  nIter = ALIGN(nIter,2);   // Even number of iterations avoids final mem-copy


 //if(levelMg>=_nLevels) nIter *= 100;


  UtTimer tm;
  int doTimer=FALSE;
  if(_psTraceDetail)
  {
    TRC(("%s %d %d %d %d\n",
	UT_FUNCNAME,levelMg,chanX,chanB,nIter));
     doTimer = levelMg<=TIMER_MAX_LEVEL;
  }

  #ifdef CHECK_NONFLUID_ZERO_PRE
  checkChannelNonFluidZero(chanX,levelMg);
  #endif

  LongInt nTotalCells=0;

  if(doTimer)
  {
    TIMER_START(&tm);
  }
  
#ifdef FINE_GRAIN_THREAD_BARRIER
  if(levelMg>=_nLevels)
  {
    ImultiplyLaplacianMatrixDenseLevel( 
		OPT_JACOBI 
		  | (boundaryZoneOnly ? OPT_BOUNDZONE_ONLY:0)
		  | (reverseOrder ? OPT_REVERSE_ORDER:0 ),
		levelMg, 
		chanX, chanB, chanTmp, 
		
		//100*nIter, 
		nIter, 

		CH_NULL );
  }
  else
#endif
  {
    std::vector<int> *blocksRelax;

    if( boundaryZoneOnly )
    {
      blocksRelax = &_blocksBoundaryRelax[levelMg];
    }
    else
    {
      blocksRelax = &_blocksRelax[levelMg];
    }

    #ifdef RELAX_BLOCKS_RED_BLACK
    int chRead = chanX,
        chWrite = chanX;
    #else
    // Copy the fixed boundary blocks of the relaxation zone
    UT_ASSERT0(FALSE);
    #endif

    unsigned opt = OPT_JACOBI; 
    if(boundaryZoneOnly) opt |= OPT_BOUNDZONE_ONLY;
    if(reverseOrder) opt  |= OPT_REVERSE_ORDER;


    {
      for(int iIter=0; iIter < nIter; iIter++)   
      {
	std::swap( chWrite, chRead );      

	nTotalCells += 
	  relaxBlockList( blocksRelax, levelMg, iIter, chRead, chanB, opt,
		   chWrite ); 
      }

      if(chWrite != chanX)
      {
	UT_ASSERT0(FALSE);  // Avoid memcopy for performance
	//copyChannel( chWrite, chanX, TRUE,-1,levelMg);
      }
    }
  }

  if(doTimer)
  {
    TIMER_STOP(&tm);
    TRC(("relax (L%d, bs=%d, nIter=%d) CPU=%.2f sec, %.0f cells, %.0f cells/sec/it\n",
	  levelMg,
	boundaryZoneOnly,nIter,
	  (double)TIMER_DIFF_MS(&tm)/1000.,(double)nTotalCells,
	(double)nTotalCells/(double)(TIMER_DIFF_MS(&tm)/1000.0)));

  }

  #ifdef CHECK_NONFLUID_ZERO_POST
  checkChannelNonFluidZero(chanX,levelMg);
  #endif
}

} // namespace MSBG

