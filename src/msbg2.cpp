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

extern "C" { 
  extern int MiUserBreak( void ); 
}

namespace MSBG 
{

#if 1
const float FACE_AREA_CLAMP_ONE = 0.9995f,
//	    FACE_AREA_CLAMP_ZERO = 0.01f;
//	    	   
            FACE_AREA_CLAMP_ZERO = 0.1f;
//	    FACE_AREA_CLAMP_ZERO = 0.05f;
#ifdef USE_CELL_VOLUME_FRACTIONS
const float CELL_VOLUME_CLAMP_ZERO = 0.01f;
#endif
#else
const float FACE_AREA_CLAMP_ONE = 0.9995f,
	    FACE_AREA_CLAMP_ZERO = 0.001f;

#endif

typedef struct
{
  int x,y,z,
      level;
}
CellLocationGen;


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::sortBlockListMorton( std::vector<int>& blockList )
{
  tbb::parallel_sort( blockList.begin(),blockList.end(), 
		      [this](const int & bid1, const int & bid2) -> bool
    {
      Vec4i bpos1 = _sg0->getBlockCoordsById(bid1),
            bpos2 = _sg0->getBlockCoordsById(bid2);
      uint64_t m1 = FMA::morton64_encode( bpos1 ),
               m2 = FMA::morton64_encode( bpos2 );
      return m1 < m2;
    } 
  );
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
size_t MultiresSparseGrid::showMemoryUsage( int trcLevel )
{

  size_t sum=0;

  for(int l=0;l<MSBG_MAX_LEVELS;l++)
  {      
     if(_sparseGrids[l].distFineCoarse) 
       sum += _sparseGrids[l].distFineCoarse->showMemUsage(trcLevel);

     if(_sparseGrids[l].genUint16_2) 
       sum += _sparseGrids[l].genUint16_2->showMemUsage(trcLevel);

     if(_sparseGrids[l].genUint8) 
       sum += _sparseGrids[l].genUint8->showMemUsage(trcLevel);
     if(_sparseGrids[l].genUint8_2) 
       sum += _sparseGrids[l].genUint8_2->showMemUsage(trcLevel);
     if(_sparseGrids[l].uint8Tmp) 
       sum += _sparseGrids[l].uint8Tmp->showMemUsage(trcLevel);
     if(_sparseGrids[l].uint8Tmp2) 
       sum += _sparseGrids[l].uint8Tmp2->showMemUsage(trcLevel);

     if(_sparseGrids[l].floatTmp3)
	 sum += _sparseGrids[l].floatTmp3->showMemUsage(trcLevel);

     if(_sparseGrids[l].divergenceAdj)
	 sum += _sparseGrids[l].divergenceAdj->showMemUsage(trcLevel);

    for(int j=0;j<ARRAY_LENGTH(_sparseGrids[l].cellFlags);j++)      
    { 
      if(_sparseGrids[l].genUint16[j]) 
       sum += _sparseGrids[l].genUint16[j]->showMemUsage(trcLevel);

      if(_sparseGrids[l].float8[j]) 
       sum += _sparseGrids[l].float8[j]->showMemUsage(trcLevel);

      if(_sparseGrids[l].cgQ[j] )
       sum += _sparseGrids[l].cgQ[j]->showMemUsage(trcLevel);

      if( _sparseGrids[l].cgP[j])
        sum += _sparseGrids[l].cgP[j]->showMemUsage(trcLevel);

      if(_sparseGrids[l].pressureOld[j] )
       sum += _sparseGrids[l].pressureOld[j]->showMemUsage(trcLevel);

      if( _sparseGrids[l].cellFlags[j])
        sum += _sparseGrids[l].cellFlags[j]->showMemUsage(trcLevel);

      if( _sparseGrids[l].cellFlagsTmp[j])
       sum += _sparseGrids[l].cellFlagsTmp[j]->showMemUsage(trcLevel);

      if( _sparseGrids[l].massDensity[j])
       sum += _sparseGrids[l].massDensity[j]->showMemUsage(trcLevel);

      if( _sparseGrids[l].float4[j])
       sum += _sparseGrids[l].float4[j]->showMemUsage(trcLevel);
    }

    for(int j=0;j<ARRAY_LENGTH(_sparseGrids[l].pressure);j++)      
      if( _sparseGrids[l].pressure[j] ) 
        sum += _sparseGrids[l].pressure[j]->showMemUsage(trcLevel);

    for(int j=0;j<ARRAY_LENGTH(_sparseGrids[l].float2);j++)      
      if( _sparseGrids[l].float2[j] )
        sum += _sparseGrids[l].float2[j]->showMemUsage(trcLevel);

    for(int j=0;j<ARRAY_LENGTH(_sparseGrids[l].float3);j++)      
      if( _sparseGrids[l].float3[j] )
        sum += _sparseGrids[l].float3[j]->showMemUsage(trcLevel);

    for(int j=0;j<ARRAY_LENGTH(_sparseGrids[l].diagonal);j++)      
      if( _sparseGrids[l].diagonal[j] )
        sum += _sparseGrids[l].diagonal[j]->showMemUsage(trcLevel);

    for(int j=0;j<ARRAY_LENGTH(_sparseGrids[l].float7);j++)      
      if( _sparseGrids[l].float7[j] )
        sum += _sparseGrids[l].float7[j]->showMemUsage(trcLevel);

    for(int j=0;j<ARRAY_LENGTH(_sparseGrids[l].floatTmpPS);j++)      
      if( _sparseGrids[l].floatTmpPS[j] )
        sum += _sparseGrids[l].floatTmpPS[j]->showMemUsage(trcLevel);

    for(int j=0;j<ARRAY_LENGTH(_sparseGrids[l].divergence);j++)
      if( _sparseGrids[l].divergence[j] )
        sum += _sparseGrids[l].divergence[j]->showMemUsage(trcLevel);

    for(int j=0;j<ARRAY_LENGTH(_sparseGrids[l].float5);j++)
      if( _sparseGrids[l].float5[j] )
        sum += _sparseGrids[l].float5[j]->showMemUsage(trcLevel);

    if(_sparseGrids[l].densityDiff )
      sum += _sparseGrids[l].densityDiff->showMemUsage(trcLevel);

    if(_sparseGrids[l].sootDiff )
      sum += _sparseGrids[l].sootDiff->showMemUsage(trcLevel);

    if(_sparseGrids[l].heatDiff )
      sum += _sparseGrids[l].heatDiff->showMemUsage(trcLevel);

    for(int j=0;j<ARRAY_LENGTH(_sparseGrids[l].float1);j++)
      if( _sparseGrids[l].float1[j] )
        sum += _sparseGrids[l].float1[j]->showMemUsage(trcLevel);

    for(int j=0;j<ARRAY_LENGTH(_sparseGrids[l].float6);j++)
      if( _sparseGrids[l].float6[j] )
        sum += _sparseGrids[l].float6[j]->showMemUsage(trcLevel);

    for(int k=0;k<3;k++)
    {
      for(int j=0;j<ARRAY_LENGTH(_sparseGrids[l].faceArea[k]);j++)
	if( _sparseGrids[l].faceArea[k][j] )
	  sum += _sparseGrids[l].faceArea[k][j]->showMemUsage(trcLevel);
    }
    if(_sparseGrids[l].heat )
      sum += _sparseGrids[l].heat->showMemUsage(trcLevel);

    if(_sparseGrids[l].curvature )
      sum += _sparseGrids[l].curvature->showMemUsage(trcLevel);

    for(int j=0;j<ARRAY_LENGTH(_sparseGrids[l].vec3_4);j++)
      if( _sparseGrids[l].vec3_4[j] )
	sum += _sparseGrids[l].vec3_4[j]->showMemUsage(trcLevel);

    for(int j=0;j<ARRAY_LENGTH(_sparseGrids[l].vec3_2);j++)
      if( _sparseGrids[l].vec3_2[j] )
	sum += _sparseGrids[l].vec3_2[j]->showMemUsage(trcLevel);

    if( _sparseGrids[l].velocityAvg )
      sum += _sparseGrids[l].velocityAvg->showMemUsage(trcLevel);

    if( _sparseGrids[l].velocityAirDiff )
      sum += _sparseGrids[l].velocityAirDiff->showMemUsage(trcLevel);

    for(int j=0;j<ARRAY_LENGTH(_sparseGrids[l].vec3_1);j++)
      if( _sparseGrids[l].vec3_1[j] )
        sum += _sparseGrids[l].vec3_1[j]->showMemUsage(trcLevel);

    for(int j=0;j<ARRAY_LENGTH(_sparseGrids[l].vec3_3);j++)
      if( _sparseGrids[l].vec3_3[j] )
         sum += _sparseGrids[l].vec3_3[j]->showMemUsage(trcLevel);

    for(int j=0;j<ARRAY_LENGTH(_sparseGrids[l].faceDensity);j++)
      if( _sparseGrids[l].faceDensity[j] )
         sum += _sparseGrids[l].faceDensity[j]->showMemUsage(trcLevel);
  }

  TRC(("%s: MSBG '%s' -> total = %.1f G\n",UT_FUNCNAME,_name,
	(double)sum/(double)ONE_GB));

  return sum;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template< int bsx, int SIMDW >
inline 
__attribute__((always_inline))

void prepareHaloBlock( 	
    	int bx0, int by0, int bz0,
        SparseGrid<float> *sg,  // Source data  
	float **haloBlockmap,  // Block map with halo
	float *dataDst )  // Destination halo block: (bsx+2*SIMDW)*bsx^2
{
  const int bsx2=bsx*bsx;
  const int hsx = bsx+2,
      hsxAct = ALIGN( bsx+2*SIMDW, SIMDW ),
      hsx2Act = hsx*hsxAct;

  int nbx = sg->nbx()+2,
      nby = sg->nby()+2,
      nbz = sg->nbz()+2,
      nbxy = nbx*nby;
  USE(nbz);

  bx0 += 1;
  by0 += 1;
  bz0 += 1;

  int x0 = bx0*bsx-SIMDW,
      y0 = by0*bsx-1,
      z0 = bz0*bsx-1;


  for(int hz=0; hz<hsx; hz++)
  {
    int z = hz+z0;
    int ibz = (z/bsx) * nbxy,
	ivz = (z & (bsx-1)) * bsx2;

    for(int hy=0; hy<hsx; hy++)
    {
      int y = hy+y0;
      int iby = ibz + (y/bsx) * nbx,
	  ivy = ivz + (y & (bsx-1)) * bsx;

      float *lineDst = dataDst + MT_GXYZ(hsxAct,hsx2Act, 0,hy,hz);

      for(int hx=0; hx<hsxAct; hx+=SIMDW)
      {
	int x = hx+x0;
	int ib = iby + (x/bsx),
	    iv = ivy + (x & (bsx-1));

	UT_ASSERT2(ib>=0&&ib<nbx*nby*nbz);
	UT_ASSERT2(iv>=0&&iv<bsx*bsx*bsx);

	float *dataSrc = haloBlockmap[ib];
	Vec4f X;
	X.load_a( dataSrc+ iv );
	X.store_a( lineDst + hx );
      }
    }
  }
}


} // namespace MSBG

