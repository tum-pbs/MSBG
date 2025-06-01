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

#include "pnoise.h"
#include <vectorclass/special/vector3d.h>

#define EIK_DBG_X 80
#define EIK_DBG_Y 96
#define EIK_DBG_Z 147

#define UT_ASSERT2_(x) 

namespace MSBG 
{


void **MultiresSparseGrid::getChannelAddr( int chanId, int level, int levelMg )
{
  UT_ASSERT2(chanId>=0&&chanId<MSBG_MAX_CHANNELS);

  if(!(level>=0&&level<_mgLevels)) return NULL;
  const Level *lev=&_sparseGrids[level]; 
  int levelMgAct = levelMg<_nLevels?levelMg:0;

  #define RET_CHAN_ADDR( T, A )\
  {\
    union \
    {\
      SparseGrid<T> *const*sg;\
      void **addr;\
    }\
    sg_addr;\
    sg_addr.sg = A;\
    return sg_addr.addr;\
  }

  switch(chanId)
  {
    case CH_NULL: return (void**)&lev->auxNullPtr;
    case CH_FLOAT_1: RET_CHAN_ADDR( float, &lev->float1[levelMgAct] );	
    case CH_DENSITY_DIFF: RET_CHAN_ADDR(  float, levelMgAct?NULL:&lev->densityDiff );			    
    case CH_FLOAT_6: RET_CHAN_ADDR( float, &lev->float6[levelMgAct] );	
    case CH_SOOT_DIFF: RET_CHAN_ADDR(  float, levelMgAct?NULL:&lev->sootDiff );			    
    case CH_FLOAT_8: RET_CHAN_ADDR(  float, &lev->float8[levelMgAct] ); 
    case CH_MASS_DENSITY: RET_CHAN_ADDR(  MassDensity, &lev->massDensity[levelMgAct] ); 
    case CH_DIVERGENCE: RET_CHAN_ADDR(  PSFloat, &lev->divergence[levelMgAct] ); 
    case CH_PRESSURE: RET_CHAN_ADDR(  PSFloat, &lev->pressure[levelMgAct] ); 
    case CH_PRESSURE_OLD: RET_CHAN_ADDR(  PSFloat, &lev->pressureOld[levelMgAct] ); 
    case CH_FLOAT_2: RET_CHAN_ADDR(  float, &lev->float2[levelMgAct] ); 
    #ifdef SOLVE_PRESSURE_DOUBLEsoo
    case CH_FLOAT_TMP_PS: RET_CHAN_ADDR(  PSFloat, &lev->floatTmpPS[levelMgAct] ); 
    #endif
    case CH_FLOAT_3: RET_CHAN_ADDR(  float, &lev->float3[levelMgAct] ); 
    case CH_FLOAT_TMP_3: RET_CHAN_ADDR(  float, levelMgAct?NULL:&lev->floatTmp3 ); 
    case CH_DIVERGENCE_ADJ: RET_CHAN_ADDR( float, levelMgAct?NULL:&lev->divergenceAdj ); 
    case CH_DIAGONAL: RET_CHAN_ADDR(  PSFloat, &lev->diagonal[levelMgAct] ); 
    case CH_FLOAT_7: RET_CHAN_ADDR(  float, &lev->float7[levelMgAct] ); 
    case CH_CG_P: RET_CHAN_ADDR(  PSFloat, &lev->cgP[levelMgAct] ); 
    case CH_CG_Q: RET_CHAN_ADDR(  PSFloat, &lev->cgQ[levelMgAct] ); 
    case CH_FLOAT_5: RET_CHAN_ADDR( float, &lev->float5[levelMgAct] ); 
    case CH_FLOAT_4: RET_CHAN_ADDR( float, &lev->float4[levelMgAct] ); 
    case CH_CURVATURE: RET_CHAN_ADDR( float, levelMgAct?NULL:&lev->curvature ); 
    case CH_HEAT: RET_CHAN_ADDR( float, levelMgAct?NULL:&lev->heat ); 
    case CH_HEAT_DIFF: RET_CHAN_ADDR(  float, levelMgAct?NULL:&lev->heatDiff );			    
    case CH_CELL_FLAGS: RET_CHAN_ADDR( CellFlags, &lev->cellFlags[levelMgAct] ); 
    case CH_CELL_FLAGS_TMP: RET_CHAN_ADDR( CellFlags, &lev->cellFlagsTmp[levelMgAct] ); 
    case CH_DIST_FINECOARSE: RET_CHAN_ADDR(  uint16_t, levelMgAct?NULL:&lev->distFineCoarse ); 
    case CH_UINT16_1: RET_CHAN_ADDR( uint16_t, &lev->genUint16[levelMgAct] ); 
    case CH_UINT16_2: RET_CHAN_ADDR( uint16_t, levelMgAct?NULL:&lev->genUint16_2 ); 
    case CH_FACE_DENSITY: RET_CHAN_ADDR( FaceDensity, &lev->faceDensity[levelMgAct] ); 
    case CH_VEC3_1: RET_CHAN_ADDR( Vec3Float, &lev->vec3_1[levelMgAct] ); 
    case CH_VELOCITY_AIR: RET_CHAN_ADDR( Vec3Float, &lev->velocityAir[levelMgAct] ); 
    case CH_VEC3_4: RET_CHAN_ADDR( Vec3Float, &lev->vec3_4[levelMgAct] ); 
    case CH_VEC3_3: RET_CHAN_ADDR( Vec3Float, &lev->vec3_3[levelMgAct] ); 
    case CH_VELOCITY_AVG: RET_CHAN_ADDR( Vec3Float, levelMgAct?NULL:&lev->velocityAvg ); 
    case CH_VELOCITY_AIR_DIFF: RET_CHAN_ADDR( Vec3Float, levelMgAct?NULL:&lev->velocityAirDiff ); 
    case CH_VEC3_2: RET_CHAN_ADDR( Vec3Float, &lev->vec3_2[levelMgAct] ); 

    case CH_UINT8_1: RET_CHAN_ADDR( uint8_t, levelMgAct?NULL:&lev->genUint8 ); 
    case CH_UINT8_2: RET_CHAN_ADDR( uint8_t, levelMgAct?NULL:&lev->genUint8_2 ); 
    case CH_UINT8_TMP: RET_CHAN_ADDR( uint8_t, levelMgAct?NULL:&lev->uint8Tmp ); 
    case CH_UINT8_TMP_2: RET_CHAN_ADDR( uint8_t, levelMgAct?NULL:&lev->uint8Tmp2 ); 

    default: return NULL; 
	     break;
  }
}


void MultiresSparseGrid::protectChannel( int iChan, bool protRead, bool protWrite )
{
  TRC3(("%s: %d %d %d\n",UT_FUNCNAME,iChan,protRead,protWrite));
  if(iChan)
  {
    FLAG_COPY(_protectedRead,(((int64_t)1)<<iChan),protRead);
    FLAG_COPY(_protectedWrite,(((int64_t)1)<<iChan),protWrite);      
  }
}


void MultiresSparseGrid::resetChannel( int ch, int ignoreNotFound, int levelMg,
		   int opt /* &1=reset, 
				&2=release physically,
				&4=prepareDataAccess-read 
				&8=prepareDataAccess-write */ )
{
  int doPrepareDataAccess = opt & (4|8);
  if(ch==-1)
  {
    if(_options & OPT_SINGLE_CHANNEL_FLOAT) 
    {
      ch = CH_GEN_FLOAT;
    }
    else if(_options & OPT_SINGLE_CHANNEL_VEC) 
    {
      ch = CH_GEN_VEC3_FLOAT;
    }
    else
    {
      UT_ASSERT0(FALSE);
    }
  }
  if(ch==CH_NULL) return;

  #if 0
  TRC(("%s %d %d %d %s\n",UT_FUNCNAME,(int)ch,levelMg,opt, 
	opt==0||(opt&1) ? "doReset":""));
  #endif

  if(!doPrepareDataAccess)
  {
    FLAG_RESET(_validForInterpolation,(((int64_t)1)<<ch));
    FLAG_RESET(_protectedRead,(((int64_t)1)<<ch));      
    FLAG_RESET(_protectedWrite,(((int64_t)1)<<ch));      
  }


  for(int l=0; l<_mgLevels; l++) 
  {
    if( isDenseLevel(l) && levelMg>=0 && levelMg!=l ) continue;

    int nLevAct = isDenseLevel(l) ? 1 : _nLevels,
	levMgAct = isDenseLevel(l) ? -1 : levelMg;

    SparseGrid<float> *sg=NULL;
    SparseGrid<PSFloat> *sgPsf=NULL;

    if(ch==CH_UINT8_1)
    {
      if(_sparseGrids[l].genUint8) _sparseGrids[l].genUint8->reset(opt);
    }
    else if(ch==CH_UINT8_2)
    {
      if(_sparseGrids[l].genUint8_2) _sparseGrids[l].genUint8_2->reset(opt);
    }
    else if(ch==CH_UINT8_TMP)
    {
      if(_sparseGrids[l].uint8Tmp) _sparseGrids[l].uint8Tmp->reset(opt);
    }
    else if(ch==CH_UINT8_TMP_2)
    {
      if(_sparseGrids[l].uint8Tmp2) _sparseGrids[l].uint8Tmp2->reset(opt);
    }
    else if(ch==CH_UINT16_1)
    {
      for(int lMg=0;lMg<nLevAct;lMg++)
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if(_sparseGrids[l].genUint16[lMg])
	  _sparseGrids[l].genUint16[lMg]->reset(opt);
      }
    }
    else if(ch==CH_UINT16_2)
    {
      if(_sparseGrids[l].genUint16_2)
	_sparseGrids[l].genUint16_2->reset(opt);
    }
    else if(ch==CH_DIST_FINECOARSE)
    {
      if(_sparseGrids[l].distFineCoarse)
	_sparseGrids[l].distFineCoarse->reset(opt);
    }
    else if(ch==CH_FLOAT_8)
    {
      for(int lMg=0;lMg<nLevAct;lMg++)
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if(_sparseGrids[l].float8[lMg])
	  _sparseGrids[l].float8[lMg]->reset(opt);
      }
    }
    else if(ch==CH_MASS_DENSITY)
    {
      for(int lMg=0;lMg<nLevAct;lMg++)
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if(_sparseGrids[l].massDensity[lMg])
	  _sparseGrids[l].massDensity[lMg]->reset(opt);
      }
    }
    else if(ch==CH_CELL_FLAGS)
    {
      for(int lMg=0;lMg<nLevAct;lMg++)
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if(_sparseGrids[l].cellFlags[lMg])
	  _sparseGrids[l].cellFlags[lMg]->reset(opt);
      }
    }
    else if(ch==CH_CELL_FLAGS_TMP)
    {
      for(int lMg=0;lMg<nLevAct;lMg++)
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if(_sparseGrids[l].cellFlagsTmp[lMg])
	  _sparseGrids[l].cellFlagsTmp[lMg]->reset(opt);
      }
    }
    else if(ch==CH_VEC3_1)
    {
      for(int lMg=0;lMg<nLevAct;lMg++)
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if(_sparseGrids[l].vec3_1[lMg])
	   _sparseGrids[l].vec3_1[lMg]->reset(opt);
      }
    }
    else if(ch==CH_VEC3_4)
    {
      for(int lMg=0;lMg<nLevAct;lMg++)
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if(_sparseGrids[l].vec3_4[lMg])
	  _sparseGrids[l].vec3_4[lMg]->reset(opt);
      }
    }
    else if(ch==CH_VEC3_2)
    {
      for(int lMg=0;lMg<nLevAct;lMg++)
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if(_sparseGrids[l].vec3_2[lMg])
	  _sparseGrids[l].vec3_2[lMg]->reset(opt);
      }
    }
    else if(ch==CH_VELOCITY_AIR)
    {
      for(int lMg=0;lMg<nLevAct;lMg++)
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if(_sparseGrids[l].velocityAir[lMg])
	  _sparseGrids[l].velocityAir[lMg]->reset(opt);
      }
    }
    else if(ch==CH_VELOCITY_AVG)
    {
      if(_sparseGrids[l].velocityAvg)
	_sparseGrids[l].velocityAvg->reset(opt);
    }
    else if(ch==CH_VELOCITY_AIR_DIFF)
    {
      if(_sparseGrids[l].velocityAirDiff)
	_sparseGrids[l].velocityAirDiff->reset(opt);
    }
    else if(ch==CH_VEC3_3)
    {
      for(int lMg=0;lMg<nLevAct;lMg++)
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if(_sparseGrids[l].vec3_3[lMg])
	  _sparseGrids[l].vec3_3[lMg]->reset(opt);
      }
    }
    else if(ch==CH_FACE_DENSITY)
    {
      for(int lMg=0;lMg<nLevAct;lMg++)
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if(_sparseGrids[l].faceDensity[lMg])
	  _sparseGrids[l].faceDensity[lMg]->reset(opt);
      }
    }
#if 0
    else
    {
      for(int lMg=0;lMg<nLevAct;lMg++)
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	SparseGrid<float> *sg = getFloatChannel(ch, l, lMg);
	if(sg) sg->reset(opt);
      }
    }
#else
    else if(ch==CH_FLOAT_1)
    {
      for(int lMg=0;lMg<nLevAct;lMg++) 
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if((sg = _sparseGrids[l].float1[lMg] )!=NULL) sg->reset(opt);
      }
    }
    else if(ch==CH_FLOAT_6)
    {
      for(int lMg=0;lMg<nLevAct;lMg++) 
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if((sg = _sparseGrids[l].float6[lMg] )!=NULL) sg->reset(opt);
      }
    }
    else if(ch==CH_FACE_AREA)
    {
      for(int k=0;k<3;k++)
      {
	for(int lMg=0;lMg<nLevAct;lMg++) 
	{
	  if(levMgAct>=0 && lMg!=levMgAct) continue;
	  if((sg = _sparseGrids[l].faceArea[k][lMg] )!=NULL) sg->reset(opt);
	}
      }
    }
    else if(ch==CH_FACE_COEFF)
    {
      for(int k=0;k<3;k++)
      {
	for(int lMg=0;lMg<nLevAct;lMg++) 
	{
	  if(levMgAct>=0 && lMg!=levMgAct) continue;
	  if((sg = _sparseGrids[l].faceCoeff[k][lMg] )!=NULL) sg->reset(opt);
	}
      }
    }
    else if(ch==CH_DENSITY_DIFF)
    {
      if((sg = _sparseGrids[l].densityDiff )!=NULL) sg->reset(opt);
    }
    else if(ch==CH_SOOT_DIFF)
    {
      if((sg = _sparseGrids[l].sootDiff )!=NULL) sg->reset(opt);
    }
    else if(ch==CH_HEAT_DIFF)
    {
      if((sg = _sparseGrids[l].heatDiff )!=NULL) sg->reset(opt);
    }
    else if(ch==CH_DIVERGENCE)
    {
      for(int lMg=0;lMg<nLevAct;lMg++) 
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if((sgPsf=_sparseGrids[l].divergence[lMg])!=NULL) sgPsf->reset(opt);
      }
    }
    else if(ch==CH_PRESSURE)
    {
      for(int lMg=0;lMg<nLevAct;lMg++) 
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if((sgPsf=_sparseGrids[l].pressure[lMg])!=NULL) sgPsf->reset(opt);
      }
    }
    else if(ch==CH_FLOAT_2)
    {
      for(int lMg=0;lMg<nLevAct;lMg++) 
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if((sg=_sparseGrids[l].float2[lMg])!=NULL) sg->reset(opt);
      }
    }
    else if(ch==CH_FLOAT_TMP_PS)
    {
      for(int lMg=0;lMg<nLevAct;lMg++) 
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	auto sg = _sparseGrids[l].floatTmpPS[lMg]; if(sg) sg->reset(opt);
      }
    }
    else if(ch==CH_FLOAT_3)
    {
      for(int lMg=0;lMg<nLevAct;lMg++) 
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if((sg=_sparseGrids[l].float3[lMg])!=NULL) sg->reset(opt);
      }
    }
    else if(ch==CH_FLOAT_TMP_3)
    {
      if((sg = _sparseGrids[l].floatTmp3 )!=NULL) sg->reset(opt);
    }
    else if(ch==CH_DIVERGENCE_ADJ)
    {
      if((sg = _sparseGrids[l].divergenceAdj )!=NULL) sg->reset(opt);
    }
    else if(ch==CH_FLOAT_5)
    {
      for(int lMg=0;lMg<nLevAct;lMg++) 
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if((sg=_sparseGrids[l].float5[lMg])!=NULL) sg->reset(opt);
      }
    }
    else if(ch==CH_DIAGONAL)
    {
      for(int lMg=0;lMg<nLevAct;lMg++) 
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if((sgPsf=_sparseGrids[l].diagonal[lMg])!=NULL) sgPsf->reset(opt);
      }
    }
    else if(ch==CH_FLOAT_7)
    {
      for(int lMg=0;lMg<nLevAct;lMg++) 
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if((sg=_sparseGrids[l].float7[lMg])!=NULL) sg->reset(opt);
      }
    }
    else if(ch==CH_CG_P)
    {
      for(int lMg=0;lMg<nLevAct;lMg++) 
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if((sgPsf=_sparseGrids[l].cgP[lMg])!=NULL) sgPsf->reset(opt);
      }
    }
    else if(ch==CH_CG_Q)
    {
      for(int lMg=0;lMg<nLevAct;lMg++) 
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if((sgPsf=_sparseGrids[l].cgQ[lMg])!=NULL) sgPsf->reset(opt);
      }
    }
    else if(ch==CH_PRESSURE_OLD)
    {
      for(int lMg=0;lMg<nLevAct;lMg++) 
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if((sgPsf=_sparseGrids[l].pressureOld[lMg])!=NULL) sgPsf->reset(opt);
      }
    }
    else if(ch==CH_FLOAT_4)
    {
      for(int lMg=0;lMg<nLevAct;lMg++) 
      {
	if(levMgAct>=0 && lMg!=levMgAct) continue;
	if((sg=_sparseGrids[l].float4[lMg])!=NULL) sg->reset(opt);
      }
    }
    else if(ch==CH_CURVATURE)
    {
      if((sg = _sparseGrids[l].curvature )!=NULL) sg->reset(opt);
    }
    else if(ch==CH_HEAT)
    {
      if((sg = _sparseGrids[l].heat )!=NULL) sg->reset(opt);
    }
#endif
    else if(!ignoreNotFound)
    {
      UT_ASSERT0(FALSE);
    }
  }
}


template< typename Data_T >
void MultiresSparseGrid::setChannelZero(  int iChan, int levelMg )
{
  prepareDataAccess( iChan );
  if(levelMg>=_nLevels)  // dense (coarse) grid implementation
  {
    UT_ASSERT0(FALSE);
  }
  else
  {
    ThrRunParallel<int>( _nBlocks, nullptr, [&](int &tls, int tid, int bid) 
    {
      BlockInfo *bi=getBlockInfo(bid,levelMg);

      SparseGrid<Data_T> *sg =  getChannel<Data_T>(this,iChan,bi->level,levelMg);
      if(!(bi->flags & BLK_EXISTS) 
	  //|| (bi->flags & BLK_NO_FLUID )
	    ) 
      {
	sg->setEmptyBlock(bid);
	return;
      }
#if 1
      Data_T *data = sg->getBlockDataPtr( bid, 1, 0 );
      SBG::setBlockDataZero( data, sg->nVoxelsInBlock());
#else
      Data_T *data = sg->getBlockDataPtr( bid, 1, 0 );
      //SBG::setBlockDataZero( data, sg->nVoxelsInBlock());

#endif
      }
    ); 
  }
}
template
void MultiresSparseGrid::setChannelZero<float>(  int iChan, int levelMg );
#ifdef SOLVE_PRESSURE_DOUBLE
template
void MultiresSparseGrid::setChannelZero<PSFloat>(  int iChan, int levelMg );
#endif
template
void MultiresSparseGrid::setChannelZero<Vec3Float>(  int iChan, int levelMg );


void MultiresSparseGrid::traceCellNeighborhoodInfo( 
	      int levelMg, int level, int x, int y, int z,  // Cell location
	      unsigned opt, 
	      CellNeighborhoodInfo<0> *cniIn
	      )
{
  char chbuf[1024];

  CellNeighborhoodInfo<0> *cni = cniIn;

  if(!cni)
  {
    cni = new CellNeighborhoodInfo<0>(this,levelMg,level,x,y,z,
		  CH_PRESSURE,
		  _liqGpActive == GP_TYPE_STD ? CH_FLOAT_1 : CH_NULL,
		  CH_VEC3_1 ); 
  }
  else
  {
    levelMg = cni->_levelMg;
    level = cni->_level;
    x = cni->_x;
    y = cni->_y;
    z = cni->_z;
  }

  CellNeighborhoodInfo<0>::Neighbor *cn0 = &cni->neighbors[0];

  #pragma omp critical
  {
  
  CellFlags flags0 = cn0->flags[0];

  if(!opt&1)
  {
    _lastErrNum++;
    _lastErrX = x;
    _lastErrY = y;
    _lastErrZ = z;
    _lastErrLevel = level;
    _lastErrLevelMG = levelMg;
  }

  TRC(("%s: levelMg=%d level=%d ipos=%d,%d,%d => flags=%s, dens=%g gen=%g\n", 
	UT_FUNCNAME,
	levelMg, level, x, y, z, 
	UT_INT2STR(flags0,2,chbuf), cn0->phi[0], cn0->genFloat[0] ));

  Vec4i pos1(x,y,z,0);

  for(int iDir=0;iDir<3;iDir++)
  {
    for(int iSide=0;iSide<2;iSide++)
    {
      int idxNeigh = idxStencil7(iDir,iSide);
      CellNeighborhoodInfo<0>::Neighbor *cn = &cni->neighbors[idxNeigh];

      Vec4i dpos = _v4iStencil7[idxNeigh],
	    pos2 = pos1+dpos;
      int x2=vget_x(pos2),
	  y2=vget_y(pos2),
	  z2=vget_z(pos2);

      TransResType transRes = cn->transRes;

      int nNeigh2 = transRes == COARSE_FINE ? 4 : 1,
	  level2 = transRes == COARSE_FINE ? level-1 :
		   transRes == FINE_COARSE ? level+1 :
		   level;

      TRC(("%d/%d: %d:%d,%d,%d => f=%s a=%.9g t=%g c=%.9g d=%g v=%g/%g\n", 
	    iDir,iSide,level2, x2, y2, z2, 
	    UT_INT2STR(cn->flags[0],2,chbuf), 
	    cn->faceArea[0], cn->theta[0], cn->faceCoeff[0], 
	    cn->phi[0], cn->faceVel[0], 
	    cn->faceVel[0]
	    ));

      if(nNeigh2>1)
      {
	for(int i=0;i<nNeigh2;i++)
	{
	  TRC(("  #%d => f=%s a=%.9g t=%g c=%.9g d=%g v=%g/%g\n", 
		i, UT_INT2STR(cn->flags[i],2,chbuf), 
		cn->faceArea[i], cn->theta[i], cn->faceCoeff[i],
		cn->phi[i], cn->faceVel[0], 
		cn->faceVel[0] 
		));
	}
      }
    }  // iSide 
  }  // iDir
  TRC(("\n"));

  } // omp critical

  if(!cniIn) delete(cni);
}

void MultiresSparseGrid::initStep( void )
{
  TRC(("%s\n",UT_FUNCNAME));
  _doDumpSimulationState = 0;
  _mgSmOmega = _mgSmOmegaDefault;
  _chDivAdjValid = 0;
  {
    for(int l=0;l<ARRAY_LENGTH(_blocksNonFrozen);l++) _blocksNonFrozen[l].clear();
    _coarseDistInterfaceLevel = -1;
  }

}

} // namespace MSBG

