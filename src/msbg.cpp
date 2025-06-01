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

#ifdef UT_ASSERT_LEVEL_2
#define PROTECT_CONSTBLOCKS
#endif

//#define TEST_WITH_FLUID_3D
#ifdef TEST_WITH_FLUID_3D
#include "wturb/FLUID_3d.h"
#endif

extern "C" { 
  extern int MiUserBreak( void ); 
}

#define GET_IP_CHANNEL( chanId, level, \
    		     sgVel, sgFloat )\
{\
  UT_ASSERT2(level>=0&&level<_nLevels);\
  UT_ASSERT2( ! ( (_protectedRead & (((int64_t)1)<<(chanId)) ) ||\
		  (_protectedWrite & (((int64_t)1)<<(chanId)) )));\
  \
  if(chanId==CH_VEC3_1)\
  {  /*Most frequent case first*/\
    sgVel = _sparseGrids[level].vec3_1[0];\
    sgFloat = NULL;\
  }\
  else\
  {\
    Level *lev=&_sparseGrids[level]; \
    sgVel = NULL;\
    sgFloat = NULL;\
\
    switch(chanId)\
    {\
      case CH_VEC3_2: sgVel = lev->vec3_2[0]; break;\
      case CH_VELOCITY_AIR_DIFF: sgVel = lev->velocityAirDiff; break;\
      case CH_VELOCITY_AIR: sgVel = lev->velocityAir[0]; break;\
      case CH_VEC3_3: sgVel = lev->vec3_3[0]; break;\
      case CH_VEC3_4: sgVel = lev->vec3_4[0]; break;\
      case CH_VELOCITY_AVG: sgVel = lev->velocityAvg; break;\
      case CH_FACE_DENSITY: sgVel = NULL /*lev->faceDensity[0]*/; break;\
      case CH_FLOAT_1: sgFloat = lev->float1[0]; break;\
      case CH_FLOAT_6: sgFloat = lev->float6[0]; break;\
      case CH_FACE_AREA: sgFloat = lev->faceArea[0][0]; break;\
      case CH_FACE_COEFF: sgFloat = lev->faceCoeff[0][0]; break;\
      case CH_DIVERGENCE: sgFloat = (SparseGrid<float>*)lev->divergence[0]; break;\
      case CH_PRESSURE: sgFloat = (SparseGrid<float>*)lev->pressure[0]; break;\
      case CH_PRESSURE_OLD: sgFloat = (SparseGrid<float>*)lev->pressureOld[0]; break;\
      case CH_CG_P: sgFloat = (SparseGrid<float>*)lev->cgP[0]; break;\
      case CH_CG_Q: sgFloat = (SparseGrid<float>*)lev->cgQ[0]; break;\
      case CH_FLOAT_2: sgFloat = lev->float2[0]; break;\
      case CH_FLOAT_3: sgFloat = lev->float3[0]; break;\
      case CH_FLOAT_TMP_3: sgFloat = lev->floatTmp3; break;\
      case CH_DIVERGENCE_ADJ: sgFloat = lev->divergenceAdj; break;\
      case CH_DIAGONAL: sgFloat = (SparseGrid<float>*)lev->diagonal[0]; break;\
      case CH_FLOAT_7: sgFloat = lev->float7[0]; break;\
      case CH_FLOAT_5: sgFloat = lev->float5[0]; break;\
      case CH_DENSITY_DIFF: sgFloat = lev->densityDiff; break;\
      case CH_SOOT_DIFF: sgFloat = lev->sootDiff; break;\
      case CH_HEAT_DIFF: sgFloat = lev->heatDiff; break;\
      case CH_FLOAT_8: sgFloat = lev->float8[0]; break;\
      case CH_MASS_DENSITY: sgFloat = lev->massDensity[0]; /*UT_ASSERT0(FALSE);*/ break;\
      case CH_FLOAT_4: sgFloat = lev->float4[0]; break;\
      case CH_CURVATURE: sgFloat = lev->curvature; break;\
      case CH_HEAT: sgFloat = lev->heat; break;\
      case CH_DIST_FINECOARSE: break;\
      case CH_UINT16_1: break;\
      case CH_UINT16_2: break;\
      default:\
	UT_ASSERT2(FALSE);\
	break;\
    }\
  }\
}

#define ZEROBLK_PADSZ (2*CPU_SIMD_WIDTH)

#define SG_DIM_EQUAL(sg1,sg2)\
  ( ((sg1)->sx()==(sg2)->sx()) && \
    ((sg1)->sy()==(sg2)->sy()) && \
    ((sg1)->sz()==(sg2)->sz()) && \
    ((sg1)->sxAct()==(sg2)->sxAct()) &&\
    ((sg1)->syAct()==(sg2)->syAct()) &&\
    ((sg1)->szAct()==(sg2)->szAct()) )

namespace MSBG 
{

using namespace SBG;

/*=========================================================================
 *
 * 	
 * 	M S B G   
 *
 * 	Multiresolution Sparse Block Grids
 *
 *
 * =======================================================================*/
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
MultiresSparseGrid
	*MultiresSparseGrid::create( 
	        const char *name,
		// Grid resolution at finest level 
		int sx0, int sy0, int sz0, 
		// Block size at fienst level (power-of-two)
		int blockSize0, 
		 // Width of transition zones between resolutin levels
		 // in units of number of finer-resolution cells (-1=auto)
		int dTransRes, 
		int maxSizeMB,
		int initialSizeMB,
		unsigned options,
		WTURB_MiscConf *miscConf,
		WTURB_PSolverConf *psConf,
		int nResLevels
		)  
{
  int rcThis=0;
  MultiresSparseGrid *msg=NULL;
  unsigned sgOptions0=0;
  if(nResLevels<0) nResLevels = 3;

#if 0
  {
    void testSetBcSpecValsSIMD();

    testSetBcSpecValsSIMD();
  }
#endif

#if 0
  for( int i=2; i>=-2; i--)
  {
    float phi = i,
	  a =  1.0f - areaFractionInside( phi, phi, phi, phi );
    TRCP(("%d => phi=%g a = %g\n",i,phi,a));
  }
    exit(1);
#endif

  {
    //
    // we assume that bit shift is handled 'arithemtically' for signed
    // integers, which should be true for most compilers 
    int a=-10;
    a >>= 1;  // 
    UT_ASSERT_FATAL(a==-5);  
  }

  UT_ASSERT0( MSBG_MAX_CHANNELS <= 64 );

  // Ensure enough increments per cell at 24 bit fixed point quantization
  UT_ASSERT0(MSBG_MAX_RESOLUTION<=65536);

  UT_ASSERT0(INVALID_VELOCITY < SBG::PCELL_INF_QUANTITY);

  UT_ASSERT0( // Constants must be consistent with config file
      MSBG::LS_TYPE_PARTICLE_SDF==2 &&  
      MSBG::LS_TYPE_PARTICLE_DENSITY==1);

  USE(axesPermutation);

  USE(maxSizeMB);  // TODO
  USE(initialSizeMB); 
  
  //UT_ASSERT_FATAL(MIN(MIN(sx0,sy0),sz0)>=48);

  if(blockSize0<0)
  {
    blockSize0 = MIN(MIN(sx0,sy0),sz0) >= 400 || 
      		 MAX(MAX(sx0,sy0),sz0) >= 1000 ? 32 : 
		 			 	 16;
  }

  int optSingleChannel = options & (OPT_SINGLE_CHANNEL_VEC | 
      				    OPT_SINGLE_CHANNEL_FLOAT
				    );

  int optSingleLevel = options & OPT_SINGLE_LEVEL;

#if 0
  dTransRes = dTransRes < 0 ? blockSize0/3 : 
    			      dTransRes;
#else
  dTransRes = 6;
//  dTransRes = 5;
#endif

  int blockSizeMin = blockSize0 >> (nResLevels-1);

  UT_ASSERT0( blockSize0 >= 2 );

  TRCP(("%s: name='%s', res0=%dx%dx%d bsx0=%d bsxMin=%d resLevels=%d "
	"dTransRes=%d staggeredVel=%d opt=%d\n",
	UT_FUNCNAME,name,sx0,sy0,sz0,blockSize0,blockSizeMin,nResLevels,dTransRes,
	1,options));

  UT_ASSERT0( nResLevels >=2 && nResLevels <= MSBG_MAXRESLEVELS );

  // Check input
  UT_ASSERT(sx0>0&&sx0<=std::numeric_limits<uint16_t>::max());
  UT_ASSERT(sy0>0&&sy0<=std::numeric_limits<uint16_t>::max());
  UT_ASSERT(sz0>0&&sz0<=std::numeric_limits<uint16_t>::max());
  UT_ASSERT(UT_ISPOW2(blockSize0));
  //UT_ASSERT(blockSize0==16||blockSize0==32);
  UT_ASSERT_FATAL(blockSize0<=32);  // uint16_t for voxel-in-block ccords

  if(!(sx0%blockSize0==0 && sy0%blockSize0==0 && sz0%blockSize0==0))
  {
    TRCERRR(("MultiresSparseGrid dimensions %dx%dx%d not multiple of blocksize %d\n",
	  sx0,sy0,sz0,blockSize0),MI_EINVARG);
  }

  // Allocate memory for object adm (zero initialized)
  ALLOCOBJ0_(msg,MultiresSparseGrid);
  if(!msg) 
  {
    TRCERRR(("MultiresSparseGrid: out of memory\n"),MI_ENOMEM);
  }

  msg = new ((void*)msg) MultiresSparseGrid();  // placement new

  UT_ASSERT_FATAL(msg);   // must be able to call destructor now 

  msg->_chAdvectorVelocity = MSBG::CH_NULL;
  
  if(miscConf)
  {
    msg->_miscConf = *miscConf;
  }

  if(psConf)
  {
    msg->_psConf = *psConf;
  }

  {
    msg->_ghost_pressure_min_theta = 0.005f;
    msg->_gp_theta_out_of_dom_liq = 1.0f;
    msg->_gp_theta_out_of_dom_air = 1.0f;
  }

  UT_ASSERT(strlen(name)<ARRAY_LENGTH(_name)-1);
  strcpy(msg->_name,name);
  msg->_nLevels = nResLevels;  
  TRCP(("number of resolution levels = %d\n",msg->_nLevels));
  UT_ASSERT0(msg->_nLevels<=MSBG_MAXRESLEVELS);
  UT_ASSERT0(msg->_nLevels <= MSBG_MAX_LEVELS );

//  UT_ASSERT(dTransRes < (blockSize0>>(msg->_nLevels-2)));

  msg->_options = options;
  msg->_totalSteps = 0;
  msg->_validForInterpolation = 0;
  msg->_protectedRead = 0;
  msg->_protectedWrite = 0;

  MT_STAT_INIT( &msg->_velMagStatLiq );
  MT_STAT_INIT( &msg->_velMagStatLiqPrev );
  MT_STAT_INIT( &msg->_velMagStatLiqFrame );
  MT_STAT_INIT( &msg->_velMagStatLiqFramePrev );

  msg->_mgUseFaceArea = 0;

  for(int i=0;i<3;i++)
  {
    for(int j=0;j<2;j++)
    {
      msg->_domainBc[i][j] = DBC_INVALID;
      msg->_domainBcOpt[i][j] = 0;
    }
  }

  msg->_dTransRes = dTransRes;

  // 
  // Misc. global auxillary data
  //
  {
    //double t0=0,t1=0,t2=0;

    for(int iWgtTyp=0;iWgtTyp<4;iWgtTyp++)
    {
      int iVelComp = iWgtTyp;
      TRC3(("iWgtTyp=%d ======================= \n",iWgtTyp));
      float *weights=NULL;
      Vec3Float goff(0.5f);
      Vec3Int dMax(2);
      if(iVelComp<3)
      {
	// MAC staggered
	goff[iVelComp] = 0.0f;
	weights = msg->_krnDownSample4x4x4[1+iVelComp];
	//dMax[iVelComp] = 1;
      }
      else
      {
	// Cell centered
	weights = msg->_krnDownSample4x4x4[0];
      }
      
      double wsum=0;
      for(int iPass=0;iPass<2;iPass++)
      {
	for(int dz=-1,ik=0;dz<=dMax.z;dz++)
	for(int dy=-1;dy<=dMax.y;dy++)
	for(int dx=-1;dx<=dMax.x;dx++,ik++)
	{
	  if(iPass==0)
	  {
	    double w = (1-0.5*fabs((goff.x)*2 - (dx+goff.x))) *
		       (1-0.5*fabs((goff.y)*2 - (dy+goff.y))) *
		       (1-0.5*fabs((goff.z)*2 - (dz+goff.z)));
	    weights[ik] = w;
	    wsum += w;
	  }
	  else 
	  {
	    weights[ik] /= wsum;
	    TRC3(("%5.4f %s%s%s",weights[ik],
		  dx==dMax.x?"\n":"",
		  dx==dMax.x&&dy==dMax.y?"\n":"",
		  dx==dMax.x&&dy==dMax.y&&dz==dMax.z?"\n":""));		
	  }
	}
      }
    }

    for(int io=0;io<4;io++)
    {    
      QuadSplineWeights *q=NULL;
      double t0=0,t1=0,t2=0;
      switch(io)
      {	
	case 0: q = &msg->_quadSplineWgtCellCenter; 
                t0=0.5,t1=0.5,t2=0.5; break;
	case 1: q = &msg->_quadSplineWgtStaggered_[io-1]; 
                t0=1.,t1=0.5,t2=0.5; break;
	case 2: q = &msg->_quadSplineWgtStaggered_[io-1];
                t0=0.5,t1=1.,t2=0.5; break;
	case 3: q = &msg->_quadSplineWgtStaggered_[io-1];
                t0=0.5,t1=0.5,t2=1.; break;
      }

      //t0=0.5,t1=0.5,t2=0.5; 

      TRC3(("%d -> %g %g %g\n",io,t0,t1,t2));

      #define GET_WEIGHTS( t, w, w_ )\
      {\
        MtQuadBSplineWeights1D( t,  w, w_ );\
	if(t==1.) /* staggered uses only w0 & w1 for data points 1 & 2*/\
	{\
	  UT_ASSERT0(fabsf(w[2])<MT_NUM_EPS);\
	  UT_ASSERT0(fabsf(w_[2])<MT_NUM_EPS);\
	  /* shift to right*/\
	  w[2]=w[1]; w[1]=w[0]; w[0]=0;\
	  w_[2]=w_[1]; w_[1]=w_[0]; w_[0]=0;\
	}\
      }

      GET_WEIGHTS( t0,  q->w0, q->w0_  );
      GET_WEIGHTS( t1,  q->w1, q->w1_  );
      GET_WEIGHTS( t2,  q->w2, q->w2_  );
      #undef GET_WEIGHTS

      for (int idx=0, oz = 0; oz <3; oz++)
      for (int oy = 0; oy <3; oy++)
      for (int ox = 0; ox <3; ox++, idx++)
      {      
	q->w[idx] = q->w0[ox] * (double)q->w1[oy] *  (double)q->w2[oz],	     
	q->wDerivs[idx][0] = - q->w0_[ox] * (double)q->w1[oy] *  (double)q->w2[oz],	      
	q->wDerivs[idx][1] = - q->w0[ox] *  (double)q->w1_[oy] * (double)q->w2[oz],
	q->wDerivs[idx][2] = - q->w0[ox] *  (double)q->w1[oy] *  (double)q->w2_[oz];

	if(io-1==0)  q->wDerivs[idx][0] *= -1;
	if(io-1==1)  q->wDerivs[idx][1] *= -1;
	if(io-1==2)  q->wDerivs[idx][2] *= -1;

	TRC3(("%3.3f %3.3f %3.3f %3.3f  %s%s",
	      q->w[idx], q->wDerivs[idx][0], q->wDerivs[idx][1], q->wDerivs[idx][2],
	      ox==2?"\n":"",ox==2&&oy==2?"\n":""))
      }
    }
  }

  // Staggered velocity grid offsets
  msg->_goffVel[0][0]=0.5; msg->_goffVel[0][1]=0.5; msg->_goffVel[0][2]=0.5;
  msg->_goffVel[1][0]=0.5; msg->_goffVel[1][1]=0.5; msg->_goffVel[1][2]=0.5;
  msg->_goffVel[2][0]=0.5; msg->_goffVel[2][1]=0.5; msg->_goffVel[2][2]=0.5;
  msg->_goffVel[0][0] = msg->_goffVel[1][1] = msg->_goffVel[2][2] = 0;
  TRC3(("Velocity grid offsets:\n"));
  for(int i=0;i<3;i++) 
    TRC3(("%g %g %g \n",msg->_goffVel[i][0], msg->_goffVel[i][1], msg->_goffVel[i][2]));

  strcpy(msg->_visOutDir,"");
  msg->_visSaveAsBitmap = 0;

  msg->_visSlicePos = msg->_visSlicePosAndRoi;
  msg->_visSliceRoiPos = msg->_visSlicePosAndRoi+3;
  msg->_visSliceRoiSize = msg->_visSlicePosAndRoi+6;

  msg->_visSlicePos[0] = 0.52;
  msg->_visSlicePos[1] = 0.52;
  msg->_visSlicePos[2] = 0.52;
  for(int k=0;k<3;k++)
  {
    msg->_visSliceRoiPos[k] = 0;
    msg->_visSliceRoiSize[k] = 1;
  }

  for(int k=0;k<ARRAY_LENGTH(msg->_obsVelocity);k++) msg->_obsVelocity[k]=0;

  // 
  // Initialize sparse grids for each level
  // 
  sgOptions0 = 0;
  #ifdef UT_ASSERT_LEVEL_2
  sgOptions0 |= OPT_PROTECT_CONST_BLOCKS;
  TRCP(("MSBG: Activate write protection of SBG const blocks\n"));
  #endif

  {
    int sx = sx0,
	sy = sy0,
	sz = sz0,
	blockSize = blockSize0;

    // restrict actually occuring block sizes 
    UT_ASSERT(blockSize>=blockSizeMin && blockSize<=32);   
    
    int sgOptions = sgOptions0 |
      		    OPT_NO_BORDPAD |
      		    OPT_BLK_ALIGNED,
        sgMaxSizeMB = 0,
	sgInitialSizeMB = -2;
    	
    for(int l=0;l<msg->_nLevels;l++)
    {
      char chbuf[100];

      UT_ASSERT(sx>=4&&sy>=4&&sz>=4);
      UT_ASSERT(blockSize>=blockSizeMin);

      if(!( optSingleLevel && l>0 ))
      {

      TRC(("Creating sparse grid at level %d : %dx%dx%d,bsx=%d\n",
	    l,sx,sy,sz,blockSize));

      MultiresSparseGrid::Level *sparseGrids = &msg->_sparseGrids[l];

      if( !(options & OPT_NO_CHAN_FCDIST) )
      {
	sprintf(chbuf,"MSG:FCDIST:%d",l);
	sparseGrids->distFineCoarse = 
	  SparseGrid<uint16_t>::create( chbuf, sx, sy, sz, blockSize,
					   sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
      }

      sprintf(chbuf,"MSG:DENSD:%d",l);
      sparseGrids->densityDiff = 
	SparseGrid<float>::create( chbuf, sx, sy, sz, blockSize,
					 sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

      //if(msg->USE_PHYSICAL_AIR_VELOCITY())
      {
	sprintf(chbuf,"MSG:SOOTD:%d",l);
	sparseGrids->sootDiff = 
	  SparseGrid<float>::create( chbuf, sx, sy, sz, blockSize,
					   sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
      }

      if(msg->_miscConf.heat_on)
      {
	sprintf(chbuf,"MSG:HEATD:%d",l);
	sparseGrids->heatDiff = 
	  SparseGrid<float>::create( chbuf, sx, sy, sz, blockSize,
					   sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
      }

      //if(! optSingleLevel )
      {
	sprintf(chbuf,"MSG:UIN16_2:%d",l);
	sparseGrids->genUint16_2 = 
	  SparseGrid<uint16_t>::create( chbuf, sx, sy, sz, blockSize,
					   sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
      }

      #ifdef RSURF_8_BIT
      if(optSingleLevel)
      {
	sprintf(chbuf,"MSG:UINT8:%d",l);
	sparseGrids->genUint8 = 
	  SparseGrid<uint8_t>::create( chbuf, sx, sy, sz, blockSize,
				       sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	sprintf(chbuf,"MSG:UINT8_2:%d",l);
	sparseGrids->genUint8_2 = 
	  SparseGrid<uint8_t>::create( chbuf, sx, sy, sz, blockSize,
				       sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	sprintf(chbuf,"MSG:UINT8T:%d",l);
	sparseGrids->uint8Tmp = 
	  SparseGrid<uint8_t>::create( chbuf, sx, sy, sz, blockSize,
				       sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	sprintf(chbuf,"MSG:UINT8T2:%d",l);
	sparseGrids->uint8Tmp2 = 
	  SparseGrid<uint8_t>::create( chbuf, sx, sy, sz, blockSize,
				       sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
      }
      #endif

      if(! optSingleLevel )
      {
	sprintf(chbuf,"MSG:FTMP3:%d",l);
	sparseGrids->floatTmp3 = 
	  SparseGrid<float>::create( chbuf, sx, sy, sz, blockSize,
					   sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	{
	  sprintf(chbuf,"MSG:VELAD:%d",l);
	  sparseGrids->velocityAirDiff = 
	    SparseGrid<Vec3Float>::create( chbuf, sx, sy, sz, blockSize,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	}
      }

      for(int lMg=0; lMg<=l; lMg++)
      {
	sprintf(chbuf,"MSG:VEC1:%d:%d",l,lMg);
	sparseGrids->vec3_1[lMg] = 
	  SparseGrid<Vec3Float>::create( chbuf, sx, sy, sz, blockSize,
					sgMaxSizeMB, sgInitialSizeMB, 					
					sgOptions | (optSingleLevel ? OPT_NO_DATA : 0)
					);     

	sprintf(chbuf,"MSG:P:%d:%d",l,lMg);
	sparseGrids->pressure[lMg] = 
	  SparseGrid<PSFloat>::create( chbuf, sx, sy, sz, blockSize,
				       sgMaxSizeMB, sgInitialSizeMB, sgOptions );     

	sprintf(chbuf,"MSG:SOBJID:%d:%d",l,lMg);
	sparseGrids->genUint16[lMg] = 
	  SparseGrid<uint16_t>::create( chbuf, sx, sy, sz, blockSize,
					   sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	sprintf(chbuf,"MSG:FLG:%d:%d",l,lMg);
	sparseGrids->cellFlags[lMg] = 
	  SparseGrid<CellFlags>::create( chbuf, sx, sy, sz, blockSize,
					   sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	sparseGrids->cellFlags[lMg]->setEmptyValue( CELL_EMPTY_VAL );

	sprintf(chbuf,"MSG:FLGTMP:%d:%d",l,lMg);
	sparseGrids->cellFlagsTmp[lMg] = 
	  SparseGrid<CellFlags>::create( chbuf, sx, sy, sz, blockSize,
					   sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	sparseGrids->cellFlagsTmp[lMg]->setEmptyValue( CELL_EMPTY_VAL );

	sprintf(chbuf,"MSG:FLT2:%d:%d",l,lMg);
	sparseGrids->float2[lMg] = 
	  SparseGrid<float>::create( chbuf, sx, sy, sz, blockSize,
					   sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	sprintf(chbuf,"MSG:FLT3:%d:%d",l,lMg);
	sparseGrids->float3[lMg] = 
	  SparseGrid<float>::create( chbuf, sx, sy, sz, blockSize,
					   sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	if(!optSingleLevel)
	{
	  sprintf(chbuf,"MSG:FLT8:%d:%d",l,lMg);
	  sparseGrids->float8[lMg] = 
	    SparseGrid<float>::create( chbuf, sx, sy, sz, blockSize,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	}

	if( ! (options & OPT_SINGLE_CHANNEL_VEC ))
	{
	  sprintf(chbuf,"MSG:DENS:%d:%d",l,lMg);
	  sparseGrids->float1[lMg] = 
	    SparseGrid<float>::create( chbuf, sx, sy, sz, blockSize,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	}
      }

      if( !optSingleChannel )
      {
	sprintf(chbuf,"MSG:VELA:%d",l);
	sparseGrids->velocityAvg = 
	  SparseGrid<Vec3Float>::create( chbuf, sx, sy, sz, blockSize,
					   sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	sprintf(chbuf,"MSG:TURBEN:%d",l);
	sparseGrids->curvature = 
	  SparseGrid<float>::create( chbuf, sx, sy, sz, blockSize,
					   sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	sprintf(chbuf,"MSG:HEAT:%d",l);
	sparseGrids->heat = 
	  SparseGrid<float>::create( chbuf, sx, sy, sz, blockSize,
					   sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	for(int lMg=0; lMg<=l; lMg++)
	{
	  sprintf(chbuf,"MSG:VEC4:%d:%d",l,lMg);
	  sparseGrids->vec3_4[lMg] = 
	    SparseGrid<Vec3Float>::create( chbuf, sx, sy, sz, blockSize,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	  sprintf(chbuf,"MSG:VEC2:%d",l);
	  sparseGrids->vec3_2[lMg] = 
	    SparseGrid<Vec3Float>::create( chbuf, sx, sy, sz, blockSize,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	  sprintf(chbuf,"MSG:VEC3:%d:%d",l,lMg);
	  sparseGrids->vec3_3[lMg] = 
	    SparseGrid<Vec3Float>::create( chbuf, sx, sy, sz, blockSize,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	  sprintf(chbuf,"MSG:FDENS:%d:%d",l,lMg);
	  sparseGrids->faceDensity[lMg] = 
	    SparseGrid<FaceDensity>::create( chbuf, sx, sy, sz, blockSize,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	  for(int iDir=0;iDir<3;iDir++)
	  {
	    sprintf(chbuf,"MSG:FAREA_%d:%d:%d",iDir,l,lMg);
	    sparseGrids->faceArea[iDir][lMg] = 
	      SparseGrid<float>::create( chbuf, sx, sy, sz, blockSize,
					       sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	    sparseGrids->faceArea[iDir][lMg]->setEmptyValue( 0.0f );
	    sparseGrids->faceArea[iDir][lMg]->setFullValue( 1.0f );
	  }

	  sprintf(chbuf,"MSG:FLT6:%d:%d",l,lMg);
	  sparseGrids->float6[lMg] = 
	    SparseGrid<float>::create( chbuf, sx, sy, sz, blockSize,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	  sprintf(chbuf,"MSG:DIV:%d:%d",l,lMg);
	  sparseGrids->divergence[lMg] = 
	    SparseGrid<PSFloat>::create( chbuf, sx, sy, sz, blockSize,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	  sprintf(chbuf,"MSG:MASSDENS:%d:%d",l,lMg);
	  sparseGrids->massDensity[lMg] = 
	    SparseGrid<MassDensity>::create( chbuf, sx, sy, sz, blockSize,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	  #if 0 
	  sprintf(chbuf,"MSG:CELLV:%d:%d",l,lMg);
	  sparseGrids->cellVolume[lMg] = 
	    SparseGrid<float>::create( chbuf, sx, sy, sz, blockSize,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	  sparseGrids->cellVolume[lMg]->setEmptyValue( 0.0f );
	  sparseGrids->cellVolume[lMg]->setFullValue( 1.0f );
          #endif

	  sprintf(chbuf,"MSG:FLT4:%d:%d",l,lMg);
	  sparseGrids->float4[lMg] = 
	    SparseGrid<float>::create( chbuf, sx, sy, sz, blockSize,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	  sprintf(chbuf,"MSG:DIAG:%d:%d",l,lMg);
	  sparseGrids->diagonal[lMg] = 
	    SparseGrid<PSFloat>::create( chbuf, sx, sy, sz, blockSize,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	  #ifdef SOLVE_PRESSURE_DOUBLE
	  sprintf(chbuf,"MSG:FTMPPS:%d:%d",l,lMg);
	  sparseGrids->floatTmpPS[lMg] = 
	    SparseGrid<PSFloat>::create( chbuf, sx, sy, sz, blockSize,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	  #endif


	  sprintf(chbuf,"MSG:CGP:%d:%d",l,lMg);
	  sparseGrids->cgP[lMg] = 
	    SparseGrid<PSFloat>::create( chbuf, sx, sy, sz, blockSize,
					 sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	  sprintf(chbuf,"MSG:CGQ:%d:%d",l,lMg);
	  sparseGrids->cgQ[lMg] = 
	    SparseGrid<PSFloat>::create( chbuf, sx, sy, sz, blockSize,
					 sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	  sprintf(chbuf,"MSG:POLD:%d:%d",l,lMg);
	  sparseGrids->pressureOld[lMg] = 
	    SparseGrid<PSFloat>::create( chbuf, sx, sy, sz, blockSize,
					 sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	}

	sprintf(chbuf,"MSG:DIVADJ:%d",l);
	sparseGrids->divergenceAdj = 
	  SparseGrid<float>::create( chbuf, sx, sy, sz, blockSize,
					   sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

      }

      }

      msg->_sg0 = msg->_sparseGrids[0].vec3_1[0];;

      msg->_domMin = msg->_sg0->v4fDomMin();
      msg->_domMax = msg->_sg0->v4fDomMax();

      if(!optSingleLevel)
      {
	msg->_scaleToNormCoords = 1./msg->_sg0->sxyzMax();
        UT_ASSERT0( msg->_sg0->nbx() * msg->_sg0->nby() * (LongInt)msg->_sg0->nbz() 
	    <= 2000000000 );
      }
      else msg->_scaleToNormCoords = 1./msg->getFlagsChannel0(0)->sxyzMax();


      if(!( optSingleLevel && l>0 ))
      {
	SparseGrid<CellFlags> *sg=msg->_sparseGrids[l].cellFlags[0];
	
	// Prepare zero data ghost blocks for safe access beyond domain border

	{
	  {
	    PSFloat *ptr;
	    ALLOCARR_ALIGNED_(ptr,PSFloat,sg->nVoxelsInBlock()+2*ZEROBLK_PADSZ);
	    memset(ptr,0,sizeof(*ptr)*(sg->nVoxelsInBlock()+2*ZEROBLK_PADSZ));
	    msg->_sparseGrids[l].zeroBlockPSFloat = ptr+ZEROBLK_PADSZ;
	    msg->_sparseGrids[l].zeroBlockPSFloatSz = sizeof(*ptr)*(sg->nVoxelsInBlock());
	  }


	  float *ptr;
	  ALLOCARR_ALIGNED_(ptr,float,sg->nVoxelsInBlock()+2*ZEROBLK_PADSZ);
	  memset(ptr,0,sizeof(*ptr)*(sg->nVoxelsInBlock()+2*ZEROBLK_PADSZ));
	  msg->_sparseGrids[l].zeroBlock = ptr+ZEROBLK_PADSZ;
	  msg->_sparseGrids[l].zeroBlockSz = sizeof(*ptr)*(sg->nVoxelsInBlock());


	  ALLOCARR_ALIGNED_(ptr,float,sg->nVoxelsInBlock()+2*ZEROBLK_PADSZ);
	  msg->_sparseGrids[l].constOneBlock = ptr+ZEROBLK_PADSZ;
	  for(int i=0;i<sg->nVoxelsInBlock();i++) 
	    msg->_sparseGrids[l].constOneBlock[i] = 1.0f;

	  ALLOCARR_ALIGNED_(ptr,float,sg->nVoxelsInBlock()+2*ZEROBLK_PADSZ);
	  msg->_sparseGrids[l].constFaceCoeffAirBlock = ptr+ZEROBLK_PADSZ;
	  for(int i=0;i<sg->nVoxelsInBlock();i++) 
	    msg->_sparseGrids[l].constFaceCoeffAirBlock[i] = 1;

	  ALLOCARR_ALIGNED_(ptr,float,sg->nVoxelsInBlock()+2*ZEROBLK_PADSZ);
	  msg->_sparseGrids[l].sdfInfiniteBlock = ptr+ZEROBLK_PADSZ;
	  for(int i=0;i<sg->nVoxelsInBlock();i++) 
	    msg->_sparseGrids[l].sdfInfiniteBlock[i] = SBG::SDF_INFINITY;
	}

	{
	  uint16_t *ptr;
	  ALLOCARR_ALIGNED_(ptr,uint16_t,sg->nVoxelsInBlock()+2*ZEROBLK_PADSZ);
	  msg->_sparseGrids[l].sdfInfiniteBlockUint16 = ptr+ZEROBLK_PADSZ;
	  for(int i=0;i<sg->nVoxelsInBlock();i++) 
	    msg->_sparseGrids[l].sdfInfiniteBlockUint16[i] = SBG::SDF_INFINITY_UI16;
	}

	{
	  Vec3Float *ptr;
	  ALLOCARR_ALIGNED_(ptr,Vec3Float,sg->nVoxelsInBlock()+2*ZEROBLK_PADSZ);
	  memset((void*)ptr,0,sizeof(*ptr)*(sg->nVoxelsInBlock()+2*ZEROBLK_PADSZ));
	  msg->_sparseGrids[l].zeroBlockVec3 = ptr+ZEROBLK_PADSZ;
	}

	{
	  CellFlags *ptr;
	  ALLOCARR_ALIGNED_(ptr,CellFlags,sg->nVoxelsInBlock()+2*ZEROBLK_PADSZ);
	  memset((void*)ptr,0,sizeof(*ptr)*(sg->nVoxelsInBlock()+2*ZEROBLK_PADSZ));
	  msg->_sparseGrids[l].ghostFlagsBlockNeumann = ptr+ZEROBLK_PADSZ;
	  for(int i=0;i<sg->nVoxelsInBlock();i++) 
	    msg->_sparseGrids[l].ghostFlagsBlockNeumann[i] = 
	      CELL_SOLID | CELL_BLK_BORDER | CELL_OUT_OF_DOMAIN;
	}

	{
	  CellFlags *ptr;
	  ALLOCARR_ALIGNED_(ptr,CellFlags,sg->nVoxelsInBlock()+2*ZEROBLK_PADSZ);
	  memset(ptr,0,sizeof(*ptr)*(sg->nVoxelsInBlock()+2*ZEROBLK_PADSZ));
	  msg->_sparseGrids[l].ghostFlagsBlockDirichlet = ptr+ZEROBLK_PADSZ;
	  for(int i=0;i<sg->nVoxelsInBlock();i++) 
	    msg->_sparseGrids[l].ghostFlagsBlockDirichlet[i] = 
	      CELL_VOID | CELL_BLK_BORDER | CELL_OUT_OF_DOMAIN;
	}

	if(l>0)
	{
          SparseGrid<CellFlags> *sg0 = msg->_sparseGrids[0].cellFlags[0];
	  UT_ASSERT( sg->nbx()==sg0->nbx() &&
		     sg->nby()==sg0->nby() &&
		     sg->nbz()==sg0->nbz() );
	}
      }

      blockSize /= 2;
      sx = sx/2;
      sy = sy/2;
      sz = sz/2;
    }

    // Pressure solver / Multigrid 
    msg->_psTraceDetail = 0;  
    msg->_mgOn = TRUE;
    msg->_mgRestrictTyp = 2;
    msg->_mgSmType = 1;  //1=jacobi,2=GS
    msg->_mgSmOmega = msg->_mgSmOmegaDefault = 6.f/7.f;  // Optimal for 3D 
    msg->_mgSmOmegaSched1 =  1.7319f;  // Yang et. al 2017. "Efficient relaxed-Jacobi smoothers for multigrid ..."
    msg->_mgSmOmegaSched2 =  0.5695f;

#if 0
    msg->_mgSmOmegaSchedM3[0] =  2.2473; 
    msg->_mgSmOmegaSchedM3[1] =  0.8571; 
    msg->_mgSmOmegaSchedM3[2] =  0.5296; 
#endif

//    msg->_mgBoundZoneWidth = 4;
    msg->_mgBoundZoneWidth = psConf ? psConf->ps_mg_bz_width : 2;

//    msg->_mgVcycPerCgIter = 2;
    msg->_mgVcycPerCgIter = 1;

    msg->_mgSmIter = psConf ? psConf->ps_mg_sm_iter : 2;
    msg->_mgSmBoundIter = psConf ? psConf->ps_mg_sm_bound_iter : 2;           
    msg->_mgSmBlockIter = psConf ? psConf->ps_mg_sm_block_iter : 5;           

    //if(msg->_mgSmIter!=2) TRCWARN(("_mgSmIter=%d\n",msg->_mgSmIter));
    //if(msg->_mgSmBlockIter!=5) TRCWARN(("_mgSmBlockIter=%d\n",msg->_mgSmBlockIter));      

    msg->_mgSxyzMin = psConf ? psConf->ps_mg_sxyzmin : 6; //6;
    msg->_mgCoarsestMaxIter = psConf ? psConf->ps_mg_coarsest_sm_iter : 4000;

    //
    // Prepare coarse level MG hierarchy of dense grids
    //


    //    msg->_mgSxyzMin = 8;
//    msg->_mgSxyzMin = 4;
    msg->_mgLevels = msg->_nLevels;

    if(!optSingleLevel)
    {
      int l = msg->_nLevels;
      int sgOptions = sgOptions0 | 
		      OPT_DENSE_GRID;

      blockSize = 1;

      while( (MIN(MIN(sx,sy),sz) >= msg->_mgSxyzMin ))
      {
	UT_ASSERT0(l<MSBG_MAX_LEVELS);
	UT_ASSERT0( 
	    (LongInt)(sx+10)*(LongInt)(sy+10)*(LongInt)(sz+10) < 2000000000 );

	TRC(("Creating dense grid at MG level %d : %dx%dx%d\n",
	      l,sx,sy,sz));
	char chbuf[100];
	MultiresSparseGrid::Level *lev = &msg->_sparseGrids[l];


	sprintf(chbuf,"MSG:SOBJID:%d",l);
	lev->genUint16[0] = 
	  SparseGrid<uint16_t>::create( chbuf, sx, sy, sz, 1,
					   sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	sprintf(chbuf,"MSG:FLG:%d",l);
	lev->cellFlags[0] = 
	  SparseGrid<CellFlags>::create( chbuf, sx, sy, sz, 1,
					 sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	lev->cellFlags[0]->setEmptyValue( CELL_EMPTY_VAL );

	SparseGrid<CellFlags> *sgFlags = lev->cellFlags[0];
	USE(sgFlags);

	sprintf(chbuf,"MSG:FLGTMP:%d",l);
	lev->cellFlagsTmp[0] = 
	  SparseGrid<CellFlags>::create( chbuf, sx, sy, sz, 1,
					 sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	lev->cellFlagsTmp[0]->setEmptyValue( CELL_EMPTY_VAL );

	sprintf(chbuf,"MSG:FLT2:%d",l);
	lev->float2[0] = 
	  SparseGrid<float>::create( chbuf, sx, sy, sz, 1,
					   sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	sprintf(chbuf,"MSG:FLT3:%d",l);
	lev->float3[0] = 
	  SparseGrid<float>::create( chbuf, sx, sy, sz, 1,
					   sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	if( ! (options & OPT_SINGLE_CHANNEL_VEC ))
	{
	  sprintf(chbuf,"MSG:DENS:%d",l);
	  lev->float1[0] = 
	    SparseGrid<float>::create( chbuf, sx, sy, sz, 1,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	}

	if( !optSingleChannel )
	{  
	  sprintf(chbuf,"MSG:FDENS:%d",l);
	  lev->faceDensity[0] = 
	    SparseGrid<FaceDensity>::create( chbuf, sx, sy, sz, 1,
					   sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	  sprintf(chbuf,"MSG:DIV:%d",l);
	  lev->divergence[0] = 
	    SparseGrid<PSFloat>::create( chbuf, sx, sy, sz, 1,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	  {
	    sprintf(chbuf,"MSG:FLT6:%d",l);
	    lev->float6[0] = 
	      SparseGrid<float>::create( chbuf, sx, sy, sz, 1,
					       sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	  }

	  for(int iDir=0;iDir<3;iDir++)
	  {
	    sprintf(chbuf,"MSG:FAREA_%d:%d",iDir,l);
	    lev->faceArea[iDir][0] = 
	      SparseGrid<float>::create( chbuf, sx, sy, sz, 1,
					       sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	  }

	  {
	    for(int iDir=0;iDir<3;iDir++)
	    {
	      sprintf(chbuf,"MSG:FCOEFF_%d:%d",iDir,l);
	      lev->faceCoeff[iDir][0] = 
		SparseGrid<float>::create( chbuf, sx, sy, sz, 1,
						 sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	    }
	  }

	  sprintf(chbuf,"MSG:P:%d",l);
	  lev->pressure[0] = 
	    SparseGrid<PSFloat>::create( chbuf, sx, sy, sz, 1,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	  sprintf(chbuf,"MSG:FLT4:%d",l);
	  lev->float4[0] = 
	    SparseGrid<float>::create( chbuf, sx, sy, sz, 1,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	  sprintf(chbuf,"MSG:MASSDENS:%d",l);
	  lev->massDensity[0] = 
	    SparseGrid<MassDensity>::create( chbuf, sx, sy, sz, 1,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	  sprintf(chbuf,"MSG:DIAG:%d",l);
	  lev->diagonal[0] = 
	    SparseGrid<PSFloat>::create( chbuf, sx, sy, sz, 1,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	  #ifdef SOLVE_PRESSURE_DOUBLE
	  sprintf(chbuf,"MSG:FTMPPS:%d",l);
	  lev->floatTmpPS[0] = 
	    SparseGrid<PSFloat>::create( chbuf, sx, sy, sz, 1,
					     sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
          #endif

	  sprintf(chbuf,"MSG:CGP:%d",l);
	  lev->cgP[0] = 
	    SparseGrid<PSFloat>::create( chbuf, sx, sy, sz, 1,
					 sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	  sprintf(chbuf,"MSG:CGQ:%d",l);
	  lev->cgQ[0] = 
	    SparseGrid<PSFloat>::create( chbuf, sx, sy, sz, 1,
					 sgMaxSizeMB, sgInitialSizeMB, sgOptions );      

	  sprintf(chbuf,"MSG:POLD:%d",l);
	  lev->pressureOld[0] = 
	    SparseGrid<PSFloat>::create( chbuf, sx, sy, sz, 1,
					 sgMaxSizeMB, sgInitialSizeMB, sgOptions );      
	} // !optSingleChannel

	sx = (sx+1)/2;
	sy = (sy+1)/2;
	sz = (sz+1)/2;

	l++;
      }
      msg->_mgLevels = l;
    }

    TRC(("Number of multigrid levels = %d\n",msg->_mgLevels));
    UT_ASSERT0(msg->_mgLevels>=msg->_nLevels);

    if(msg->_mgSmBoundIter<0)
    {
      msg->_mgSmBoundIter = MAX((int)(1.1*MT_SQR(msg->_mgLevels)+0.5),15);
      TRCP(("_mgSmBoundIter = %d\n",msg->_mgSmBoundIter));
    }

    SparseGrid<CellFlags> *sg0 = msg->_sparseGrids[0].cellFlags[0];

    msg->_nBlocks = sg0->nBlocks();
    
    UT_ASSERT0(msg->_nBlocks<2000000000);  // use of 32-bit types ( e.g. LONG for InterlockedIncrement etc.) 

    msg->_v4iDomMax = Vec4i( sg0->sx()-1,
			     sg0->sy()-1,
			     sg0->sz()-1, 
			     0x7fffffffL // INT32_MAX 
			     );
    msg->_bsx0 = sg0->bsx();
    msg->_bsx0Log2 = sg0->bsxLog2();
    msg->_blkStridY = sg0->nbx();
    msg->_blkStridZ = sg0->nbxy();

    ALLOCARR_ALIGNED_(msg->_blockTimeLevels,int8_t,sg0->nBlocks());
    #if 0
    ALLOCARR_ALIGNED_(msg->_mgPerBlockResid,float,sg0->nBlocks());
    #endif

    for(int lMg=0;lMg<msg->_nLevels;lMg++)
    {
      ALLOCARR_ALIGNED_(msg->_blockmap[lMg],BlockInfo,sg0->nBlocks());
      for(LongInt i=0;i<sg0->nBlocks();i++) 
      {
	BlockInfo *bi=&msg->_blockmap[lMg][i];
	memset(bi,0,sizeof(*bi));
	bi->level = msg->_nLevels-1;
	if(lMg==0) 
	{
	  msg->_blockTimeLevels[i] = 0;
	}
      }
    }

    msg->_blockListTmpMaxLen = sg0->nBlocks();
    ALLOCARR_ALIGNED_(msg->_blockListTmp,uint32_t,
		      msg->_blockListTmpMaxLen);

    msg->_dx0 = 1;
    
    msg->_h0 = 1./sg0->sxyzMax();

    TRC(("grid spacing at finest resolution _dx0 = %g / %g\n",
	  msg->_dx0,msg->_h0));

  }

  UT_ASSERT0(msg->_mgLevels<=MSBG_MAX_LEVELS);
  for(int level=0; level<msg->_mgLevels; level++)
  {
    msg->_nActCells[level] = 
    msg->_nFluidCells[level] = 
    msg->_nSolidCells[level] =
    msg->_nVoidCells[level] = -1;
    Level *lev = &msg->_sparseGrids[level];
    lev->auxInvScale = 1./((double)(1<<level));

    // Register channel pointers for fast access
    for(int iChan=0;iChan<MSBG_MAX_CHANNELS;iChan++)
    {
      UT_ASSERT0( iChan < ARRAY_LENGTH(lev->channelPointers));
      void **pChan = msg->getChannelAddr( iChan,level, 0 );
      lev->channelPointers[iChan] = pChan;
    }
  }

  msg->setDomainBoundarySpec_( 
	DBC_OPEN, DBC_OPEN, // left,right (X)
	DBC_OPEN, DBC_OPEN, // bottom,top (Y)
	DBC_OPEN, DBC_OPEN // front,back (Z)
	);

#if 0
  {
    void *buffers[MI_MAX_THREADS]={0};
    #pragma omp parallel
    {
      int tid=omp_get_thread_num();
      void *p=UtAllocMemNuma( 100000 );
      TRCP(("tid=%d -> %p\n",tid,p));
      buffers[tid]=p;
    }
    for(int i=0;i<MI_MAX_THREADS;i++) UtFreeMemNuma(&buffers[i]);
  }
#endif

  TRC(("sizeof(_blocksAdaptiveRelax) = %d %d\n",
	sizeof(_blocksAdaptiveRelax),sizeof(std::vector<int>)));

rcCatch:
  if(rcThis)
  {
    destroy(msg);
  }
  return msg;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::checkChannelNonFluidZero( int chan, int levelMg )
{
  //TRC(("%s\n",UT_FUNCNAME));
  if(levelMg>=_nLevels)  // special case: dense (coarse) grid
  {
    SparseGrid<PSFloat> *sg = getPSFloatChannel(chan,levelMg,levelMg);
    PSFloat *data = sg->getDataPtrGen_();
    SparseGrid<CellFlags> *sgFlags = getFlagsChannel(CH_CELL_FLAGS,levelMg,levelMg);      
    CellFlags *dataFlags = sgFlags->getDataPtrGen_();

    ThrRunParallel<int> ( sg->sz(), nullptr,   
      [&](int &tls, int tid, int z)
      {
        for(int y=0;y<sg->sy();y++)
        for(int x=0;x<sg->sx();x++)
        {	
	  int i = sg->getGridIndex(x,y,z);
          UT_ASSERT0(
	      UT_IMPLIES( !CELL_IS_FLUID_(dataFlags[i]) || (dataFlags[i]&CELL_FIXED),
			  data[i]==(PSFloat)(0)));
        } 
      });
  }
  else
  {
    ThrRunParallel<int> ( _nBlocks, nullptr,   
      [&](int &tls, int tid, int bid)
      {
	BlockInfo *bi=getBlockInfo(bid,levelMg);
	SparseGrid<CellFlags> *sgFlags=getFlagsChannel(CH_CELL_FLAGS,bi->level,levelMg);      
	CellFlags *dataFlags=sgFlags->getBlockDataPtr(bid);
	
	PSFloat *data = 
	  getPSFloatChannel( chan, bi->level, levelMg )->getBlockDataPtr( bid, 0, 0);
			

	if( data && dataFlags )
	{
	  for(int i=0;i<sgFlags->nVoxelsInBlock();i++)
	  {
	    if(!(UT_IMPLIES(
		    !CELL_IS_FLUID_(dataFlags[i]) || (dataFlags[i]&CELL_FIXED),
		    data[i]==(PSFloat)(0))))
	    {
	      Vec4i pos = sgFlags->getGridCoords(bid,i);
	      traceCellNeighborhoodInfo(levelMg-1,bi->level,pos[0],pos[1],pos[2]);
	      traceCellNeighborhoodInfo(levelMg,bi->level,pos[0],pos[1],pos[2]);
	      UT_ASSERT0(FALSE);
	    }
	  }
	}
	else
	{
	  UT_ASSERT0( (bi->flags & BLK_NO_FLUID) || !(bi->flags & BLK_EXISTS));
	}
      }
    #if 0
    ,nullptr
    ,0,/*doSerial=*/true 
    #endif

    );
  }
}

/*-------------------------------------------------------------------------------------*/
/* 										       */
/*-------------------------------------------------------------------------------------*/
template<typename Data_T>
void MultiresSparseGrid::setChannelConstValue( int chan, int constval, Data_T val, 
    					       int levelMg )
{
  if(chan==CH_NULL) return;
  UT_ASSERT2(chan>0&&chan<MSBG_MAX_CHANNELS);
  UT_ASSERT2(levelMg==-1 || (levelMg>=0&&levelMg<_nLevels));

    int lMg1 = levelMg==-1 ? 0 : levelMg,
	lMg2 = levelMg==-1 ? _nLevels-1 : levelMg;

    for(int levelMg=lMg1;levelMg<=lMg2;levelMg++)
      for(int level=levelMg; level<_nLevels; level++)
      {
	auto sg = getChannel<Data_T>(this,chan,level,levelMg);
	if(sg) sg->setConstValue( constval, val );
      }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
// Enforce resolution transitions not larger than one level
void MultiresSparseGrid::regularizeRefinementMap( int *blockmapIn,
    					      int *blockmapOut, 
					      unsigned *blockFlags,
					      int rbzWidthBlk
					      )
{
  SparseGrid<Vec3Float> *sg0 = _sg0;
  int8_t *blockmapTmp=NULL,*blockmapTmp2=NULL;
  unsigned *blockFlagsTmp=NULL;  
  
  int nLevels = _nLevels;
  UtTimer tm;
  LongInt nActVoxels = 0, nBlocks0=0, nBlocks1=0, nBlocks2=0, 
	  		  nBlocks3=0, nBlocks4=0;

  if(!blockmapOut) blockmapOut=blockmapIn;

  TRC(("%s %p %p %p %d\n",UT_FUNCNAME,blockmapIn, blockmapOut,blockFlags,rbzWidthBlk));

  UT_ASSERT0( UT_EQUIVALENT( blockFlags, strcmp(_name,"FLUIDSIM")==0 ));
      
  TIMER_START(&tm);

  UT_ASSERT0(nLevels>=2 && nLevels<=5);

  if(blockFlags)
  {
    ALLOCARR_(blockFlagsTmp,unsigned,sg0->nBlocks());
  }
  ALLOCARR_(blockmapTmp,int8_t,sg0->nBlocks());
  ALLOCARR_(blockmapTmp2,int8_t,sg0->nBlocks());

  ThrRunParallel<int>(sg0->nBlocks(), nullptr,
    [&](int &tls, int tid, int i) 
    {
      blockmapTmp[i] = blockmapIn[i];
      if(blockFlags) blockFlagsTmp[i] = blockFlags[i];
    } );

  ThrRunParallel<int>(sg0->nBlocks(), nullptr,
    [&](int &tls, int tid, int bid) 
    {
      int bx,by,bz;
      sg0->getBlockCoordsById(bid,bx,by,bz);
      
      if(blockFlags)
      {
	FLAG_RESET( blockFlags[bid], BLK_TMP_MARK );
	#ifdef USE_BLK_NEAR_FLUID
	FLAG_RESET( blockFlags[bid], BLK_EX_NEAR_FLUID );
	#endif
      }
      int level = blockmapTmp[bid];
      int level0 = 0;
      if(level>level0)
      {
	#if 0
	if(!( blockFlags && !(blockFlagsTmp[bid]&BLK_NO_FLUID)))
	{
	  int level2=nLevels-1,
	      nActiveNeigh=0;
	  for(int dz=-1;dz<=1;dz++)
	  for(int dy=-1;dy<=1;dy++)
	  for(int dx=-1;dx<=1;dx++)
	  {
	    if(dx==0&&dy==0&&dz==0) continue;
	    int bx2=bx+dx,
		by2=by+dy,
		bz2=bz+dz;
	    if(!sg0->blockCoordsInRange(bx2,by2,bz2)) continue;
	    LongInt bid2 = sg0->getBlockIndex(bx2,by2,bz2);
	    UT_ASSERT2(bid2>=0&&bid2<sg0->nBlocks());
	    level2 = MIN(level2, blockmapTmp[bid2]);	  
	    nActiveNeigh += (blockFlags && !(blockFlagsTmp[bid2] & BLK_NO_FLUID)) ? 1 : 0;
	  }
	  if(level2 == level0) 
	  {
	    level=level0;
	    if(blockFlags) blockFlags[bid] |= BLK_TMP_MARK;
	  }
	  if(nActiveNeigh && blockFlags) blockFlags[bid] |= BLK_TMP_MARK;
	}
	#endif
      }
      blockmapTmp2[bid] = level;
    } );

    LongInt nBlocksPresent=0;
    #pragma omp parallel for if (sg0->nBlocks()>10000) \
    			     reduction(+:nBlocksPresent)
    for(LongInt i=0;i<sg0->nBlocks();i++) 
    {
      blockmapTmp[i] = -1;
      if( blockFlags )
      {
	FLAG_RESET( blockFlags[i], BLK_EXISTS );
	if( (!(blockFlags[i] & BLK_NO_FLUID) || (blockFlags[i] & BLK_TMP_MARK)))
	{
	  blockFlags[i] |= BLK_EXISTS;
	  nBlocksPresent++;
	}
      }
    }
    TRC(("%s '%s': %.1f%% blocks physically present in grid\n",UT_FUNCNAME,_name,
	  100.*(double)nBlocksPresent/(double)sg0->nBlocks()));

  for(int ll=1; ll<=nLevels-2; ll++)
  {
    int isLastPass = (ll==nLevels-2);

    typedef struct
    {
      LongInt nActVoxels,nBlocks0,nBlocks1,nBlocks2,nBlocks3,nBlocks4;
    } 
    ThreadLocals;

    ThrRunParallel<ThreadLocals>(sg0->nBlocks(), nullptr,
      [&](ThreadLocals &tls, int tid, int bid) 
      {
	int bx,by,bz;
	sg0->getBlockCoordsById(bid,bx,by,bz);
	int level = blockmapTmp2[bid];
	if(level>=ll)
	{
	  // level2 = level of finest resolution neighbor
	  int level2=nLevels-1;
	  for(int dz=-1;dz<=1;dz++)
	  for(int dy=-1;dy<=1;dy++)
	  for(int dx=-1;dx<=1;dx++)
	  {
	    if(dx==0&&dy==0&&dz==0) continue;
	    int bx2=bx+dx,
		by2=by+dy,
		bz2=bz+dz;
	    if(!sg0->blockCoordsInRange(bx2,by2,bz2)) continue;
	    LongInt bid2 = sg0->getBlockIndex(bx2,by2,bz2);
	    UT_ASSERT2(bid2>=0&&bid2<sg0->nBlocks());
	    level2 = MIN(level2, blockmapTmp2[bid2]);	  

	    #ifdef USE_BLK_NEAR_FLUID
	    if(isLastPass && blockFlags && !(blockFlagsTmp[bid2] & BLK_NO_FLUID))
	    {
	      blockFlags[bid] |= BLK_EX_NEAR_FLUID;
	    }
            #endif
	  }
	  UT_ASSERT0(ll>0&&level2>=ll-1);
	  if(level2 == ll-1) 
	  {
	    level=ll;
	  }
	}
	blockmapOut[bid] = level;
	//msg->setBlockLevel(bid,level);

	if(isLastPass)
	{
	  int m = (sg0->bsx())>>level;
	  tls.nActVoxels += m*m*m;
	  switch(level)
	  {
	    case 0: tls.nBlocks0++; break;
	    case 1: tls.nBlocks1++; break;
	    case 2: tls.nBlocks2++; break;
	    case 3: tls.nBlocks3++; break;
	    case 4: tls.nBlocks4++; break;
	    default:
	      UT_ASSERT0(FALSE);
	  }
	}
      },

      [&]( ThreadLocals &tls, int tid ) // reduce
      {
        nActVoxels += tls.nActVoxels;
	nBlocks0 += tls.nBlocks0;
	nBlocks1 += tls.nBlocks1;
	nBlocks2 += tls.nBlocks2;
	nBlocks3 += tls.nBlocks3;
	nBlocks4 += tls.nBlocks4;
      } );
    
    if(!isLastPass)
    {
      for(int i=0;i<sg0->nBlocks();i++) blockmapTmp2[i] = blockmapOut[i];
    }
  }

  LongInt nTotalVoxels0 = (LongInt)sg0->nBlocks() * 
			  (LongInt)sg0->bsx()*sg0->bsx()*sg0->bsx();

  UT_ASSERT0(nActVoxels<=nTotalVoxels0);
  if(nLevels>2)
  {
    UT_ASSERT0(nBlocks0+nBlocks1+nBlocks2+nBlocks3==sg0->nBlocks());
  }

#ifdef UT_ASSERT_LEVEL_2
  #pragma omp parallel for
  for(LongInt i=0;i<sg0->nBlocks();i++) 
  {
    int level = blockmapOut[i]; 
    UT_ASSERT2(level>=0 && level<=nLevels-1);
  }
#endif

  TIMER_STOP(&tm);
  TRC(("CPU (%s) = %.2f sec, %.0f blocks/sec)\n",UT_FUNCNAME,
      (double)TIMER_DIFF_MS(&tm)/1000.,
      (sg0->nBlocks())/
      (double)(TIMER_DIFF_MS(&tm)/1000.0)));

  TRC(("MSBG: voxel usage = %g%%, blocks/level = %.0f %.0f %.0f %.0f %.0f of %.0f\n",
      100.*((double)nActVoxels/((double)nTotalVoxels0)),
      (double)nBlocks0,(double)nBlocks1,(double)nBlocks2,
      (double)nBlocks3,(double)nBlocks4,
      (double)sg0->nBlocks()));

  FREEMEM(blockFlagsTmp);
  FREEMEM(blockmapTmp);
  FREEMEM(blockmapTmp2);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::setRefinementMap( 
    					   int *blockLevels, 
					   unsigned *blockFlags,
					   int levelMgMax,
					   MSBG::MultiresSparseGrid *msbgObstacles,
					   bool doInitCellFlags,
					   bool doResetChannels
					   )
{
//#define SRFM_TRACE_DETAIL

  TRC(("%s %d%d \n",UT_FUNCNAME,_options&OPT_SINGLE_LEVEL, doInitCellFlags));

  UtTimer tm;

  TIMER_START(&tm);

  bool blockLevelsAllocated = FALSE;
  if(!blockLevels)
  {
    ALLOCARR0_(blockLevels,int,nBlocks());
    blockLevelsAllocated = TRUE;
  }

  int doVisualize = !!(strcmp(_name,"FLUIDSIM")==0);  // avoid overwrite static pnl !

  bool doBlockLists = !(_options & OPT_SINGLE_LEVEL);

  if(_options & OPT_SINGLE_LEVEL)
  {
    levelMgMax = 0;
  }
  else if(levelMgMax==-1) levelMgMax=_mgLevels-1;

  UT_ASSERT0(levelMgMax>=0&&levelMgMax<_mgLevels);

#ifdef UT_ASSERT_LEVEL_2
    // Check refinement map
    #pragma omp parallel for collapse(2) 
    for(int bz=0;bz<_sg0->nbz();bz++) 
    for(int by=0;by<_sg0->nby();by++)
    for(int bx=0;bx<_sg0->nbx();bx++)
    {
      LongInt bid = _sg0->getBlockIndex(bx,by,bz);
      int level = blockLevels[bid];
      UT_ASSERT0(level>=0&&level<_nLevels);
      for(int dz=-1;dz<=1;dz++)
      for(int dy=-1;dy<=1;dy++)
      for(int dx=-1;dx<=1;dx++)
      {
        if(dx==0&&dy==0&&dz==0) continue;
	int bx2=bx+dx,
	    by2=by+dy,
	    bz2=bz+dz;
	if(!_sg0->blockCoordsInRange(bx2,by2,bz2)) continue;
	int bid2 = _sg0->getBlockIndex(bx2,by2,bz2);
	UT_ASSERT2(bid2>=0&&bid2<_sg0->nBlocks());
        int level2 = blockLevels[bid2];
	UT_ASSERT0( ABS(level-level2) <= 1);
	/*if(bid==33) 
	{
	  TRCP(("(%d,%d,%d) + %2d %2d %2d -> %d/%d\n",bx,by,bz,dx,dy,dz,bid2,level2));
	}*/
      }
    }
#endif

  if(doResetChannels)
  {
    for(int chan=1;chan<MSBG_MAX_CHANNELS;chan++)
    {

      /*if( chan != CH_CELL_FLAGS &&
	  chan != CH_DIST_FINECOARSE )*/
      #ifdef SRFM_TRACE_DETAIL
      TRC1(("reset channel %d \n",chan));
      #endif

      {
	FLAG_RESET(_validForInterpolation,(((int64_t)1)<<chan));
	for(int levelMg=0; levelMg<=levelMgMax; levelMg++)
	  resetChannel(chan,TRUE,levelMg); 
      }
    }
  }
  _validForInterpolation = 0;
  _protectedWrite = 0;
  _protectedRead = 0;


  for(int levelMg=0; levelMg<=levelMgMax; levelMg++)
  {
    LongInt nActCells=0,nActCellsBlkFluid=0;

    int isDenseLevel = levelMg >= getNumLevels();

    resetChannel(CH_CELL_FLAGS,0,levelMg);

    auto doAllocFlagsData = [&](int level, SparseGrid<CellFlags> *sgFlags) -> bool 
    {
      if( _options & OPT_SINGLE_LEVEL )
      {
	if( !((!isDenseLevel && level==0 && sgFlags )) ) return false;
      }
      return true;
    };
    
    if( doInitCellFlags && ( ! (_options & OPT_SINGLE_LEVEL) ))
    {
      prepareDataAccess( CH_CELL_FLAGS, -1, levelMg );

      if(levelMg==0)
      {
	// Allocate constant flags blocks for non-BLK_EXISTS blocks
	for(int l=0;l<_nLevels;l++)
	{
          if(!_sparseGrids[l].constBlockFlagsNonexSolid)
	  {
	    SparseGrid<CellFlags> *sg = getFlagsChannel(CH_CELL_FLAGS,l,0);
	    Block<CellFlags> *block = sg->allocBlock(0,0,/*isConstBlock=*/true);
	    CellFlags *data = sg->getBlockDataPtr( block );
	    for(int i=0;i<sg->nVoxelsInBlock();i++) 
	    {
	      data[i] = CELL_SOLID | CELL_OUT_OF_DOMAIN;
	    }
	    _sparseGrids[l].constBlockFlagsNonexSolid = block;
	    #ifdef PROTECT_CONSTBLOCKS
	    sg->protectBlock(block);
	    #endif

	    block = sg->allocBlock(0,0,/*isConstBlock=*/true);
	    data = sg->getBlockDataPtr( block );
	    for(int i=0;i<sg->nVoxelsInBlock();i++) 
	    {
	      data[i] = CELL_VOID | CELL_OUT_OF_DOMAIN;
	    }
	    _sparseGrids[l].constBlockFlagsNonexVoid = block;
	    #ifdef PROTECT_CONSTBLOCKS
	    sg->protectBlock(block);
	    #endif
	  }
	}
      }
    }

    if(!(_options & OPT_SINGLE_LEVEL) && 
       !isDenseLevel && levelMg==0)
    {
      for(int l=0;l<_nLevels;l++) 
      {
	SparseGrid<uint16_t> *sg = getDistFineCoarseChannel( l );
	if(sg) sg->prepareDataAccess();
      }
    }

    if(!isDenseLevel)
    {
      if(levelMg==0) 
      {
        _blocksExisting.clear();
	for(int l=0;l<_nLevels;l++) _blocksFineCoarsePerLevel[l].clear();
      }
    }

    int bxmax = _sg0->nbx()-1,
	bymax = _sg0->nby()-1,
	bzmax = _sg0->nbz()-1;

    int nbx = _sg0->nbx(),
	//nby = _sg0->nby(),
	nbxy = _sg0->nbxy();

    {

    typedef struct
    {
      LongInt nActCells,
	      nActCellsBlkFluid;
      std::vector<int> 
	blocksExisting,
	blocksFineCoarsePerLevel[MSBG_MAXRESLEVELS];
    } 
    ThreadLocals;

    /*---------------------------------------------------------------------*/
    /* Parallel loop over blocks					   */
    /*---------------------------------------------------------------------*/
    ThrRunParallel<ThreadLocals> ( 
      
      nBlocks(), // Number of tasks

      [&]( ThreadLocals &tls, int tid )  // Initialize locals
      {
        if(levelMg==0)
	{
          tls.blocksExisting.reserve(256);
          for(int i=0;i<_nLevels;i++) tls.blocksFineCoarsePerLevel[i].reserve(256);
	}
        //TRCP(("L=%d tid=%d nActCells(init) = %g\n",levelMg,tid,(double)tls.nActCells));
      },

      // Run parallel
      [&](ThreadLocals &tls, int tid, int bid)
      {
	UT_ASSERT0(blockLevels[bid]>=0&&blockLevels[bid]<_nLevels);

	BlockInfo *bi=NULL,
		  *biObst=NULL;
	int level;
	int bx=0,by=0,bz=0;
	int blockCrossSolid = 0,
	    blockInSolid = 0;
	// 
	// Initialize block flags
	//
	if( isDenseLevel )
	{
	  level = levelMg;
	}
	else
	{
	  bi=getBlockInfo(bid,levelMg);
	  bi->level =  level = MAX( blockLevels[bid], levelMg );
	  bi->flags = blockFlags ? blockFlags[bid] : 0;

	  if( msbgObstacles )
	  {
	    biObst = msbgObstacles->getBlockInfo(bid);
	    UT_ASSERT0(!((biObst->flags & BLK_ONLY_FLUID) &&
			 (biObst->flags & BLK_NO_FLUID)));
	    blockInSolid = biObst->flags & BLK_ONLY_FLUID;
	    blockCrossSolid = !(biObst->flags & (BLK_ONLY_FLUID|BLK_NO_FLUID));
	  }

	  if(levelMg>0) 
	  {
	    if( _options & OPT_SINGLE_LEVEL) 
	    {
	      bi->flags |= BLK_EXISTS;
	    }
	    else
	    {
	      // Block contains fluid or solid mass (within safety distance)	      
	      if(( bi->flags & BLK_EXISTS ) || blockCrossSolid )
	      {
		FLAG_SET(bi->flags, BLK_EXISTS);
	      }
	      else
	      {
		FLAG_RESET(bi->flags,BLK_EXISTS);
	      }
	    }
	  }

	  if(levelMg==0 && (bi->flags & BLK_EXISTS))
	  {
	    tls.blocksExisting.push_back(bid);
	  }

	  _sg0->getBlockCoordsById(bid,bx,by,bz); 
	  FLAG_COPY( bi->flags, BLK_DOM_BORDER,
		     bx==0 || bx==_sg0->nbx()-1 ||
		     by==0 || by==_sg0->nby()-1 ||
		     bz==0 || bz==_sg0->nbz()-1 );
	  #if 0
	  if(levelMg==1&&level==1 && bx==0 && by==0 && bz==0)
	  {
	    TRCERR(("DBG\n"));
	  }
	  #endif
	  

	  int level2Max=0,
	      level2Min=_nLevels-1;
	  for(int dz=-1;dz<=1;dz++)
	  for(int dy=-1;dy<=1;dy++)
	  for(int dx=-1;dx<=1;dx++)
	  {
	    if(dx==0&&dy==0&&dz==0) continue;
	    int bx2=bx+dx,
		by2=by+dy,
		bz2=bz+dz;
	    if(!_sg0->blockCoordsInRange(bx2,by2,bz2)) continue;
	    LongInt bid2 = _sg0->getBlockIndex(bx2,by2,bz2);
	    UT_ASSERT2(bid2>=0&&bid2<_sg0->nBlocks());
	    int level2 = MAX( blockLevels[bid2], levelMg );	    	    
	    level2Max = MAX(level2Max, level2);	  
	    level2Min = MIN(level2Min, level2);	  
	  }
	  
	  UT_ASSERT2( (level2Max - level <= 1) &&
		      (level - level2Min <= 1) );
	      	  
	  if(level2Max > level) 
	  {
	    bi->flags |= BLK_FINE_COARSE;
	    if(!isDenseLevel && doBlockLists ) 
	    {
	      if(levelMg==0) tls.blocksFineCoarsePerLevel[level].push_back(bid);
	      /*TRCP(("capcity = %d, size = %d\n",
		    blocksResBorder_.capacity(),
		    blocksResBorder_.size()));*/
	    }
	  }
	  if(level2Min < level) 
	  {
	    bi->flags |= BLK_COARSE_FINE;
	  }

	  FLAG_RESET(bi->flags,BLK_UNIFORM_EMPTY);
	  // Empty block that crosses no obstacle and no resolution border 
	  if( ! (bi->flags & BLK_EXISTS) && 
	      ! BLK_IS_RES_BORDER( bi->flags ) &&	  
	      UT_IMPLIES( biObst, (biObst->flags & (BLK_ONLY_FLUID|BLK_NO_FLUID)) ))
	  {
	    bi->flags |= BLK_UNIFORM_EMPTY;
	  }
	}

	//
	// Initialize cell flags & distFineCoarse
	// 
	SparseGrid<CellFlags> *sgFlags = 
	  getFlagsChannel(CH_CELL_FLAGS,level,levelMg);

	if(doInitCellFlags)
	{
	  if( !doAllocFlagsData( level, sgFlags ) ) return;

	  if(!isDenseLevel && (bi->flags & BLK_EXISTS)) 
	    tls.nActCells += sgFlags->nVoxelsInBlock();

          #ifdef CONST_FLAGS_UNIFORM_BLOCKS
	  if(!isDenseLevel && strcmp(_name,"FLUIDSIM")==0)
	  {
	    if( bi->flags & BLK_UNIFORM_EMPTY )
	    {
	      sgFlags->setBlock(bid, biObst->flags & BLK_ONLY_FLUID ? 
		    _sparseGrids[level].constBlockFlagsNonexSolid :
		    _sparseGrids[level].constBlockFlagsNonexVoid);
	      return;
	    }
	  }
	  #endif

	  #if 0 
	  if(!isDenseLevel && strcmp(_name,"FLUIDSIM")==0 && levelMg<1)
	  {
	    if(!( bi->flags & BLK_EXISTS) && !blockCrossSolid )
	    {
	      // TODO: set const-block with all CELL_OPEN | CELL_OUT_OF_DOMAIN | CELL_NON_EXISTENT
	      // with CELL_NON_EXISTENT = CELL_COARSE_FINE | CELL_FINE_COARSE
	      UT_ASSERT0(!sgFlags->getBlockDataPtr(bid));
	      return;
	    }
	  }
	  #endif
	  
	  CellFlags *dataFlags = sgFlags->getDataPtrGen_(bid,1,0);

	  if(levelMg==0 && bi && !( bi->flags & BLK_NO_FLUID )) 
	    tls.nActCellsBlkFluid += sgFlags->nVoxelsInBlock();

	  if(isDenseLevel)
	  {
	    SBG_FOR_EACH_VOXEL_IN_BLOCK_GEN( sgFlags, bid, x, y, z, idx )
	    {
	      CellFlags flags=0;
	      dataFlags[idx] = flags;
	    }	
	  }
	  else
	  {
	    int bsx = sgFlags->bsx(),
		level2,bid2;

	    for(int vz=0,vid=0;vz<bsx;vz++) 
	    {
	      CellFlags flagsZ=0;	
	      int bstridZ=0;
	      if(vz==0)
	      {
		flagsZ |= CELL_BLK_BORDER;
		if(bz!=0) bstridZ = -nbxy;
	      }
	      else if(vz==bsx-1)
	      {
		flagsZ |= CELL_BLK_BORDER;
		if(bz!=bzmax) bstridZ = nbxy;
	      }

	      for(int vy=0;vy<bsx;vy++) 
	      {
		CellFlags flagsY = 0;
		int bstridY = 0;
		if(vy==0)
		{
		  flagsY |= CELL_BLK_BORDER;
		  if(by!=0) bstridY = -nbx;
		}
		else if(vy==bsx-1)
		{
		  flagsY |= CELL_BLK_BORDER;
		  if(by!=bymax) bstridY = nbx;
		}

		for(int vx=0;vx<bsx;vx++,vid++) 
		{
		  CellFlags flagsX = 0;
		  int bstridX = 0;
		  if(vx==0)
		  {
		    flagsX |= CELL_BLK_BORDER;
		    if(bx!=0) bstridX = -1;
		  }
		  else if(vx==bsx-1)
		  {
		    flagsX |= CELL_BLK_BORDER;
		    if(bx!=bxmax) bstridX = 1;
		  }

		  CellFlags flags = flagsX | flagsY | flagsZ;

		  if(!(_options & OPT_SINGLE_LEVEL) && 
		      BLK_IS_RES_BORDER(bi->flags))
		  {
		    #define UPD_TRANS_RES_FLAGS( bstrid )\
		    if(bstrid)\
		    { \
		      bid2 = bid + bstrid;\
		      UT_ASSERT2(bid2>=0&&bid2<_nBlocks);\
		      level2 = MAX( blockLevels[bid2], levelMg );	 \
		      if(level2 > level) flags |= CELL_FINE_COARSE; \
		      if(level2 < level) flags |= CELL_COARSE_FINE; \
		    }

		    UPD_TRANS_RES_FLAGS(bstridX);
		    UPD_TRANS_RES_FLAGS(bstridY);
		    UPD_TRANS_RES_FLAGS(bstridZ);

		    UPD_TRANS_RES_FLAGS(bstridX + bstridY)
		    UPD_TRANS_RES_FLAGS(bstridX + bstridZ)
		    UPD_TRANS_RES_FLAGS(bstridY + bstridZ)

		    UPD_TRANS_RES_FLAGS(bstridX + bstridY + bstridZ);
		  }

		  if(!(bi->flags & BLK_EXISTS)) 
		  {
		    flags |= ( blockInSolid ? CELL_SOLID : CELL_VOID );
		  }
		  dataFlags[vid] = flags;
		}
	      }
	    }
	  } // Sparse level 
	} // doInitCellFlags

	//
	// Create CH_DIST_FINECOARSE for smooth level transitions 
	//	
	SparseGrid<uint16_t> *sgDistFineCoarse = getDistFineCoarseChannel( level );
	if( !(_options & OPT_SINGLE_LEVEL) && 
	    !isDenseLevel && levelMg==0 && sgDistFineCoarse )
	{
	  /*if(!(bi->flags & BLK_EXISTS))
	  {
	    sgDistFineCoarse->setEmptyBlock(bid);
	  }
	  else*/ if( bi->flags & BLK_FINE_COARSE )
	  {
	    uint16_t *dataDistFineCoarse = sgDistFineCoarse->getBlockDataPtr(bid,1,0);

	    int bsx=sgFlags->bsx();
	    int nbx=_sg0->nbx(), 
		nby=_sg0->nby(), 
		nbz=_sg0->nbz();

	    Vec4f blkActPos[27];
	    int blkActNum=0;

	    for(int dz=-1;dz<=1;dz++)   
	    {
	      int bz2=bz+dz;
	      if(bz2<0||bz2>nbz-1) continue;
	      int isCenter = dz==0;
	      int bid2z = bz2*nbx*nby;
	      for(int dy=-1;dy<=1;dy++)   
	      {
		int by2=by+dy;
		if(by2<0||by2>nby-1) continue;
		isCenter &= dy==0;
		int bid2y = bid2z+by2*nbx;
		for(int dx=-1;dx<=1;dx++)
		{
		  int bx2=bx+dx;
		  if(bx2<0||bx2>nbx-1) continue;
		  isCenter &= dx==0;
		  if(isCenter) continue;
		  UT_ASSERT2((_sg0->blockCoordsInRange(bx2,by2,bz2)));
		  int bid2 = bid2y+bx2,
		      level2 =blockLevels[bid2];
		  if(level2 <= level) continue;
		  UT_ASSERT0(blkActNum<ARRAY_LENGTH(blkActPos));
		  blkActPos[blkActNum++] = Vec4f(bx2,by2,bz2,0.0f);
		}
	      }
	    }

	    for(int vid=0, z1=bz*bsx; z1<(bz+1)*bsx; z1++)
	    for(int y1=by*bsx; y1<(by+1)*bsx; y1++)
	    for(int x1=bx*bsx; x1<(bx+1)*bsx; x1++, vid++)
	    {
	      Vec4f pos1(x1+0.5f,y1+0.5f,z1+0.5f,0.0f);
	      float dist=-1;
	      float dist2=1e20;
	      for(int i=0;i<blkActNum;i++)   
	      {
	        Vec4f bpos2 = blkActPos[i];
		dist = distToBoxSq( bpos2*bsx, (bpos2+1.f)*bsx, pos1 );					 
		dist2 = MIN(dist2,dist);
	      }
	      UT_ASSERT2(!(dist<0));
	      dist = sqrtf(dist2);
	      if(dist<0 || dist>_dTransRes)
	      {
		dist = _dTransRes;
	      }
	      UT_ASSERT2(dist<65535/1024);

	      uint16_t uiDist = dist*1024;
	      dataDistFineCoarse[vid] = uiDist;
	    }

	  }
	}
      },

      [&]( ThreadLocals &tls, int tid )  // Reduce locals
      {
	if(!isDenseLevel)
	{
	  if(levelMg==0)
	  {
	    UT_VECTOR_APPEND( _blocksExisting, tls.blocksExisting );
	    for(int l=0;l<_nLevels;l++) 
	      UT_VECTOR_APPEND( _blocksFineCoarsePerLevel[l], 
		                tls.blocksFineCoarsePerLevel[l] );
	  }
	}
	nActCells += tls.nActCells;
	nActCellsBlkFluid += tls.nActCellsBlkFluid;
      } );

    } 
    
    _nActCells[levelMg] = isDenseLevel ? 
      getFlagsChannel(CH_CELL_FLAGS,levelMg,levelMg)->nTotVirtVoxels() :
      nActCells;

    #ifdef SRFM_TRACE_DETAIL
    TRC(("_nActCells[%d] = %g\n",levelMg,(double)_nActCells[levelMg]));
    #endif
    //
    // Make sure blocks are defined to support
    // interpolation stencil for the distance field
    // 
    if(  !(_options & OPT_SINGLE_LEVEL) && 
	!isDenseLevel && levelMg==0 && getDistFineCoarseChannel( ) )
    {
      for(int level=0; level<_nLevels; level++)
      {
	SparseGrid<uint16_t> *sgDistFineCoarse = getDistFineCoarseChannel( level );
	sgDistFineCoarse->setFullValue(_dTransRes*1024);

	ThrRunParallel<int>( nBlocks(), nullptr,
	  [&](int &tls, int tid, int bid) 
	  {
	    BlockInfo *bi=getBlockInfo0(bid);
	    int level1 = bi->level;
	    if(level1>level||level1<level-1) return;	  
	    if(sgDistFineCoarse->isValueBlock(bid)) return;	  
	    int bx,by,bz;
	    _sg0->getBlockCoordsById(bid, bx,by,bz);
	    int found=FALSE;
	    for(int dz=-1;dz<=1;dz++)   
	    for(int dy=-1;dy<=1;dy++)   
	    for(int dx=-1;dx<=1;dx++)
	    {
	      if(dx==0&&dy==0&&dz==0) continue;
	      int bx2=bx+dx,
		  by2=by+dy,
		  bz2=bz+dz;
	      if(!_sg0->blockCoordsInRange(bx2,by2,bz2)) continue;
	      int bid2 = _sg0->getBlockIndex(bx2, by2, bz2);
	      if( getBlockLevel(bid2) == level &&
		 (_blockmap[0][bid2].flags & BLK_FINE_COARSE))
	      {
		found=TRUE;
	      }
	    }
	    if(found) sgDistFineCoarse->setFullBlock(bid);
	  }
	);
      }
    }

    if(levelMg == 0 ) 
    {

      #ifdef MSBG_SORT_BLOCKLISTS
      {
	UtTimer tm;
	TIMER_START(&tm);

	sortBlockListMorton( _blocksExisting );

	TIMER_STOP(&tm);
	int nBlocks = _blocksExisting.size();
	TRC(("L%d sort block list (_blocksExisting)  CPU=%.2f sec, %.0f blocks, %.0f blocks/sec\n",
	      levelMg, (double)TIMER_DIFF_MS(&tm)/1000.,
	      (double)nBlocks, nBlocks/(double)(TIMER_DIFF_MS(&tm)/1000.0)));
      }
      #endif

      _nActCellsBlkFluid = nActCellsBlkFluid;

      char chbuf[256];
      sprintf(chbuf,"MSBG usage efficiency ('%s') = %g%% (%g%%)\n", _name,
	    100.*(double)_nActCells[0]/(double)_sg0->nTotVirtVoxels(),
	    100.*(double)_nActCellsBlkFluid/(double)_sg0->nTotVirtVoxels());
      if(_options & (OPT_SINGLE_CHANNEL_FLOAT | OPT_SINGLE_CHANNEL_VEC ))
      {
        TRC(("%s",chbuf));
      }
      else
      {
        TRCP(("%s",chbuf));
      }
    }

    // Visualization

    if( doVisualize && (levelMg < _nLevels && (_visTopology & 8) ))
    {
      //
      // Refinement level map
      //
      {
	static int pnl[MSBG_MAX_LEVELS]={-2};
	if(pnl[0]==-2) for(int i=0;i<ARRAY_LENGTH(pnl);i++) pnl[i]=-1;
	char chbuf[200];
	sprintf(chbuf,"Refinement level after setRefinementMap L=%d",levelMg);
	visualizeSlices( MiGetAuxPanel2(&pnl[levelMg],chbuf),
	  CH_CELL_FLAGS,IP_NEAREST,NULL,0,levelMg,0,0,NULL,
	  
	  [&]( float f, MtStat *fstat,
	      int flags, int blockFlags, int level, const Vec3Int& ipos,
	      int &colorOut ) -> void
	  {
	    Vec3Int c(254*(((_nLevels-1)-level)/(float)(_nLevels-1)));
	    colorOut = BMP_MKRGB( c.x, c.y, c.z );
	  } );
      }

      //
      // Cell flags
      //
      {
	static int pnl[MSBG_MAX_LEVELS]={-2};
	if(pnl[0]==-2) for(int i=0;i<ARRAY_LENGTH(pnl);i++) pnl[i]=-1;
	char chbuf[200];
	sprintf(chbuf,"Topology (cell flags) after setRefinementMap L=%d",levelMg);
	visualizeSlices( MiGetAuxPanel2(&pnl[levelMg],chbuf),
	  CH_CELL_FLAGS,IP_NEAREST,NULL,0,levelMg,0,0,NULL,
	  
	  [&]( float f, MtStat *fstat,
	      int flags, int blockFlags, int level, const Vec3Int& ipos,
	      int &colorOut ) -> void
	  {
	    Vec3Int c(0);

	    if( flags & CELL_FINE_COARSE ) c.x = 255;
	    if( flags & CELL_COARSE_FINE ) c.y = 255;
	    if( flags & CELL_BLK_BORDER ) c.z = 255;

	    if( !(flags & (CELL_FINE_COARSE|CELL_COARSE_FINE|CELL_BLK_BORDER) ))
	    {
	      c.x=c.y=c.z = 255;
	    }

	    int cl = 255 - level * (255.0/_nLevels);

	    for(int k=0;k<3;k++) c[k] = (float)c[k] * cl/255.;

	    colorOut = BMP_MKRGB( c.x, c.y, c.z );
	  } );
      }

      {
	static int pnl[MSBG_MAX_LEVELS]={-2};
	if(pnl[0]==-2) for(int i=0;i<ARRAY_LENGTH(pnl);i++) pnl[i]=-1;
	char chbuf[200];
	sprintf(chbuf,"Topology cell flags(II) after setRefinementMap L=%d",levelMg);
	visualizeSlices( MiGetAuxPanel2(&pnl[levelMg],chbuf),
	  CH_CELL_FLAGS,IP_NEAREST,NULL,VIS_OPT_TEST,levelMg);
      }

      //
      // Block flags
      //      
      {
	static int pnl[MSBG_MAX_LEVELS]={-2};
	if(pnl[0]==-2) for(int i=0;i<ARRAY_LENGTH(pnl);i++) pnl[i]=-1;
	char chbuf[200];
	sprintf(chbuf,"Bblocks after setRefinementMap L=%d  (R=uniform-empty,G=res-border,B=existing)",levelMg);

	visualizeBlockFlags( MiGetAuxPanel2(&pnl[levelMg], chbuf),	
	    BLK_UNIFORM_EMPTY, BLK_COARSE_FINE|BLK_FINE_COARSE, BLK_EXISTS, levelMg );
      }
    }
  } // levelMg

  TIMER_STOP(&tm);
  TRC(("CPU (%s '%s') = %.2f sec, %.0f cells/sec)\n",UT_FUNCNAME,getName(),
      (double)TIMER_DIFF_MS(&tm)/1000.,
      _nActCells[0]/
      (double)(TIMER_DIFF_MS(&tm)/1000.0)));


  if(blockLevelsAllocated)
  {
    FREEMEM(blockLevels);
  }

  resetChannel(CH_VEC3_3);
  resetChannel(CH_FLOAT_2);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void 
MultiresSparseGrid::destroy( MultiresSparseGrid*& msg ) 
{
  if(msg)
  {
    TRC(("%s %p\n",UT_FUNCNAME,(void*)msg));

    #ifdef PROTECT_CONSTBLOCKS
    for(int l=0;l<msg->_nLevels;l++)
    {
      {
	SparseGrid<CellFlags> *sg = msg->getFlagsChannel(CH_CELL_FLAGS,l,0);	 
	Block<CellFlags> *block;
	block = msg->_sparseGrids[l].constBlockFlagsNonexSolid;
	if(block) sg->unprotectBlock(block);
	block = msg->_sparseGrids[l].constBlockFlagsNonexVoid;
	if(block) sg->unprotectBlock(block);
      }
      {
	SparseGrid<PSFloat> *sg = msg->getPSFloatChannel(CH_DIAGONAL,l,0);	 
	Block<PSFloat> *block;
	block = msg->_sparseGrids[l].constBlockDiagFullAir;
	if(block) sg->unprotectBlock(block);
	block = msg->_sparseGrids[l].constBlockDiagFullLiq;
	if(block) sg->unprotectBlock(block);	
      }
      {
	SparseGrid<float> *sg = msg->getFaceCoeffChannel(0,0,0);
	if(sg)
	{
	  Block<float> *block;
	  block = msg->_sparseGrids[l].constBlockFaceCoeffFullAir;
	  if(block) sg->unprotectBlock(block);
	  block = msg->_sparseGrids[l].constBlockFaceCoeffFullLiq;
	  if(block) sg->unprotectBlock(block);	
	}
      }
    }
    #endif

    FREEMEM_ALIGNED(msg->_blockListTmp);

    FREEMEM_ALIGNED(msg->_blockTimeLevels);
    FREEMEM_ALIGNED(msg->_mgPerBlockResid);

    for(int lMg=0;lMg<msg->_nLevels;lMg++)
    {
      FREEMEM_ALIGNED(msg->_blockmap[lMg]);
    }

    for(int l=0;l<MSBG_MAX_LEVELS;l++)
    {      
      // 
      // XXX: when adding new channels keep in sync with showMemoryUsage()
      //
      SparseGrid<uint16_t>::destroy( msg->_sparseGrids[l].distFineCoarse );
      SparseGrid<uint16_t>::destroy( msg->_sparseGrids[l].genUint16_2 );

      SparseGrid<uint8_t>::destroy( msg->_sparseGrids[l].genUint8 );
      SparseGrid<uint8_t>::destroy( msg->_sparseGrids[l].genUint8_2 );
      SparseGrid<uint8_t>::destroy( msg->_sparseGrids[l].uint8Tmp );
      SparseGrid<uint8_t>::destroy( msg->_sparseGrids[l].uint8Tmp2 );

      SparseGrid<float>::destroy( msg->_sparseGrids[l].floatTmp3 );
      SparseGrid<float>::destroy( msg->_sparseGrids[l].divergenceAdj );

      for(int j=0;j<ARRAY_LENGTH(msg->_sparseGrids[l].cellFlags);j++)      
      {
        SparseGrid<uint16_t>::destroy( msg->_sparseGrids[l].genUint16[j] );
	SparseGrid<float>::destroy( msg->_sparseGrids[l].float5[j] );
        SparseGrid<float>::destroy( msg->_sparseGrids[l].float8[j] );
	SparseGrid<PSFloat>::destroy( msg->_sparseGrids[l].cgQ[j] );
	SparseGrid<PSFloat>::destroy( msg->_sparseGrids[l].cgP[j] );
	SparseGrid<PSFloat>::destroy( msg->_sparseGrids[l].pressureOld[j] );
        SparseGrid<CellFlags>::destroy( msg->_sparseGrids[l].cellFlags[j] );
        SparseGrid<CellFlags>::destroy( msg->_sparseGrids[l].cellFlagsTmp[j] );
      }

      for(int j=0;j<ARRAY_LENGTH(msg->_sparseGrids[l].pressure);j++)      
        SparseGrid<PSFloat>::destroy( msg->_sparseGrids[l].pressure[j] );

      for(int j=0;j<ARRAY_LENGTH(msg->_sparseGrids[l].massDensity);j++)      
        SparseGrid<MassDensity>::destroy( msg->_sparseGrids[l].massDensity[j] );

      for(int j=0;j<ARRAY_LENGTH(msg->_sparseGrids[l].float4);j++)      
        SparseGrid<float>::destroy( msg->_sparseGrids[l].float4[j] );

      for(int j=0;j<ARRAY_LENGTH(msg->_sparseGrids[l].float2);j++)      
      {
	SparseGrid<float>::destroy( msg->_sparseGrids[l].float3[j] );
        SparseGrid<float>::destroy( msg->_sparseGrids[l].float2[j] );
      }

      for(int j=0;j<ARRAY_LENGTH(msg->_sparseGrids[l].diagonal);j++)
        SparseGrid<PSFloat>::destroy( msg->_sparseGrids[l].diagonal[j] );

      for(int j=0;j<ARRAY_LENGTH(msg->_sparseGrids[l].float7);j++)
        SparseGrid<float>::destroy( msg->_sparseGrids[l].float7[j] );

      for(int j=0;j<ARRAY_LENGTH(msg->_sparseGrids[l].floatTmpPS);j++)
        SparseGrid<PSFloat>::destroy( msg->_sparseGrids[l].floatTmpPS[j] );

      for(int j=0;j<ARRAY_LENGTH(msg->_sparseGrids[l].float1);j++)
        SparseGrid<float>::destroy( msg->_sparseGrids[l].float1[j] );

      for(int j=0;j<ARRAY_LENGTH(msg->_sparseGrids[l].float6);j++)
        SparseGrid<float>::destroy( msg->_sparseGrids[l].float6[j] );

      for(int iDir=0;iDir<3;iDir++)
      {
	for(int j=0;j<ARRAY_LENGTH(msg->_sparseGrids[l].faceArea[iDir]);j++)
	  SparseGrid<float>::destroy( msg->_sparseGrids[l].faceArea[iDir][j] );
      }

      for(int iDir=0;iDir<3;iDir++)
      {
	for(int j=0;j<ARRAY_LENGTH(msg->_sparseGrids[l].faceCoeff[iDir]);j++)
	  SparseGrid<float>::destroy( msg->_sparseGrids[l].faceCoeff[iDir][j] );
      }

      for(int j=0;j<ARRAY_LENGTH(msg->_sparseGrids[l].divergence);j++)
        SparseGrid<PSFloat>::destroy( msg->_sparseGrids[l].divergence[j] );

      SparseGrid<float>::destroy( msg->_sparseGrids[l].densityDiff );
      SparseGrid<float>::destroy( msg->_sparseGrids[l].sootDiff );
      SparseGrid<float>::destroy( msg->_sparseGrids[l].heatDiff );
      SparseGrid<float>::destroy( msg->_sparseGrids[l].curvature );
      SparseGrid<float>::destroy( msg->_sparseGrids[l].heat );
      
      for(int j=0;j<ARRAY_LENGTH(msg->_sparseGrids[l].vec3_4);j++)
        SparseGrid<Vec3Float>::destroy( msg->_sparseGrids[l].vec3_4[j] );

      for(int j=0;j<ARRAY_LENGTH(msg->_sparseGrids[l].vec3_2);j++)
        SparseGrid<Vec3Float>::destroy( msg->_sparseGrids[l].vec3_2[j] );

      for(int j=0;j<ARRAY_LENGTH(msg->_sparseGrids[l].velocityAir);j++)
        SparseGrid<Vec3Float>::destroy( msg->_sparseGrids[l].velocityAir[j] );

      for(int j=0;j<ARRAY_LENGTH(msg->_sparseGrids[l].faceDensity);j++)
        SparseGrid<FaceDensity>::destroy( msg->_sparseGrids[l].faceDensity[j] );

      SparseGrid<Vec3Float>::destroy( msg->_sparseGrids[l].velocityAvg );
      SparseGrid<Vec3Float>::destroy( msg->_sparseGrids[l].velocityAirDiff );

      for(int j=0;j<ARRAY_LENGTH(msg->_sparseGrids[l].vec3_1);j++)
      {
        SparseGrid<Vec3Float>::destroy( msg->_sparseGrids[l].vec3_1[j] );
        SparseGrid<Vec3Float>::destroy( msg->_sparseGrids[l].vec3_3[j] );
      }

      // 
      // XXX: when adding new channels keep in sync with showMemoryUsage()
      //
      {
	{
	  {
	    PSFloat *ptr=msg->_sparseGrids[l].zeroBlockPSFloat;
	    if(ptr)
	    {
	      ptr -= ZEROBLK_PADSZ;
	      FREEMEM_ALIGNED( ptr );
	      msg->_sparseGrids[l].zeroBlockPSFloat = NULL;
	    }
	  }

	  float *ptr=msg->_sparseGrids[l].zeroBlock;
	  if(ptr)
	  {
	    ptr -= ZEROBLK_PADSZ;
	    FREEMEM_ALIGNED( ptr );
	    msg->_sparseGrids[l].zeroBlock = NULL;
	  }

	  ptr=msg->_sparseGrids[l].constFaceCoeffAirBlock;
	  if(ptr)
	  {
	    ptr -= ZEROBLK_PADSZ;
	    FREEMEM_ALIGNED( ptr );
	    msg->_sparseGrids[l].constFaceCoeffAirBlock = NULL;
	  }

	  ptr=msg->_sparseGrids[l].constOneBlock;
	  if(ptr)
	  {
	    ptr -= ZEROBLK_PADSZ;
	    FREEMEM_ALIGNED( ptr );
	    msg->_sparseGrids[l].constOneBlock = NULL;
	  }

	  ptr=msg->_sparseGrids[l].sdfInfiniteBlock;
	  if(ptr)
	  {
	    ptr -= ZEROBLK_PADSZ;
	    FREEMEM_ALIGNED( ptr );
	    msg->_sparseGrids[l].sdfInfiniteBlock = NULL;
	  }
	}

	{
	  uint16_t *ptr=msg->_sparseGrids[l].sdfInfiniteBlockUint16;
	  if(ptr)
	  {
	    ptr -= ZEROBLK_PADSZ;
	    FREEMEM_ALIGNED( ptr );
	    msg->_sparseGrids[l].sdfInfiniteBlockUint16 = NULL;
	  }
	}

	{
	  Vec3Float *ptr=msg->_sparseGrids[l].zeroBlockVec3;
	  if(ptr)
	  {
	    ptr -= ZEROBLK_PADSZ;
	    FREEMEM_ALIGNED( ptr );
	    msg->_sparseGrids[l].zeroBlockVec3 = NULL;
	  }
	}

	{
	  CellFlags *ptr=msg->_sparseGrids[l].ghostFlagsBlockNeumann;
	  if(ptr)
	  {
	    ptr -= ZEROBLK_PADSZ;
	    FREEMEM_ALIGNED( ptr );
	    msg->_sparseGrids[l].ghostFlagsBlockNeumann = NULL;
	  }
	}
	{
	  CellFlags *ptr=msg->_sparseGrids[l].ghostFlagsBlockDirichlet;
	  if(ptr)
	  {
	    ptr -= ZEROBLK_PADSZ;
	    FREEMEM_ALIGNED( ptr );
	    msg->_sparseGrids[l].ghostFlagsBlockDirichlet = NULL;
	  }
	}
      }
    }

    SparseGrid<float>::destroy( msg->_sgCoarseDistToRefineRoi );

    msg->~MultiresSparseGrid();
    FREEMEM( msg );
  }
}

template<>
SparseGrid<uint8_t> *getChannel<uint8_t>( MSBG::MultiresSparseGrid *msbg,
				              int chanId, int level, int levelMg )
{
  return msbg->getUint8Channel(chanId,level,levelMg);
}

template<>
SparseGrid<CellFlags> *getChannel<CellFlags>( MSBG::MultiresSparseGrid *msbg,
				              int chanId, int level, int levelMg )
{
  return msbg->getUint16Channel(chanId,level,levelMg);
}

template<>
SparseGrid<float> *getChannel<float>( MSBG::MultiresSparseGrid *msbg,
				      int chanId, int level, int levelMg )
{
  return msbg->getFloatChannel(chanId,level,levelMg);
}

template<>
SparseGrid<double> *getChannel<double>( MSBG::MultiresSparseGrid *msbg,
				      int chanId, int level, int levelMg )
{
  return (SparseGrid<double> *)msbg->getGenChannel(chanId,level,levelMg);
}

template<>
SparseGrid<Vec3Float> *getChannel<Vec3Float>( MSBG::MultiresSparseGrid *msbg,
					      int chanId, int level, int levelMg )
{
  return msbg->getVecChannel(chanId,level,levelMg);
}

template<>
SparseGrid<Vec3Uint16> *getChannel<Vec3Uint16>( MSBG::MultiresSparseGrid *msbg,
					        int chanId, int level, int levelMg )
{
  UT_ASSERT2(chanId == CH_FACE_DENSITY);
  return (SparseGrid<Vec3Uint16> *)msbg->getGenChannel(chanId,level,levelMg);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MultiresSparseGrid::setDomainBoundarySpec_( 
    			int dbLeft,  int dbRight, // left,right (X)
    			int dbBottom, int dbTop, // bottom,top (Y)
			int dbFront, int dbBack,  // front,back (Z)
			int *dbOpt // int[3][2]
			)
{
  int rcThis=0;

  _domainBc[0][0] = dbLeft;    _domainBc[0][1] = dbRight;
  _domainBc[1][0] = dbBottom;  _domainBc[1][1] = dbTop;
  _domainBc[2][0] = dbFront;   _domainBc[2][1] = dbBack;

  for(int i=0;i<3;i++)
  {
    for(int j=0;j<2;j++)
    {
      _domainBcOpt[i][j] = dbOpt ? dbOpt[i*2+j] : 0;
      TRC(("db[%d][%d] = %d/%d\n",i,j,_domainBc[i][j],_domainBcOpt[i][j]));
    }
  }

  for(int i=0;i<3;i++)
  {
    for(int j=0;j<2;j++)
    {
      if(_domainBc[i][j] == DBC_PERIODIC)
      {
	UT_ASSERT0(_domainBc[i][1-j] == DBC_PERIODIC);
      }
    }
  }

  for(int i=0;i<3;i++) 
  {
    for(int j=0;j<2;j++)
    {
      int dbc = _domainBc[i][j];
      CellFlags c;
      switch(dbc)
      {
	case DBC_OPEN: c =  CELL_VOID;
	  break;
	case DBC_INFLOW: c =   CELL_VOID;
	  break;
	case DBC_OUTFLOW: c =  CELL_VOID;
	  break;
	case DBC_SOLID: c =  CELL_SOLID;
	  break;
	case DBC_PERIODIC: c = 0;
	  UT_ASSERT0(FALSE);  // TODO		   
	  break;
	default:
	  TRCERRR(("invalid domain boundary type at [%d,%d] = %d\n",i,j,dbc),
	      MI_EINVARG);
	  break;
      }

      //UT_ASSERT0(!dbSolidFreeSlip); // TODO

      _domainBcFlags[i][j] = c;
    }
  }

  for(int z=-1;z<=1;z++)
  for(int y=-1;y<=1;y++)
  for(int x=-1;x<=1;x++)
  {
    int iDir=-1,iFace=-1;
    if(z<0) { iDir=2; iFace=0; } 
    if(z>0) { iDir=2; iFace=1; }

    if(y<0) { iDir=1; iFace=0; } 
    if(y>0) { iDir=1; iFace=1; } 

    if(x<0) { iDir=0; iFace=0; } 
    if(x>0) { iDir=0; iFace=1; } 	

    float val=0;

    if(iDir>=0)
    {
      val = _domainBc[iDir][iFace] == DBC_OPEN ? SPEC_FLOAT_DIRICHLET :
	                                         SPEC_FLOAT_NEUMANN;
    }

    _ghostFloatSpecVals3x3x3[MT_GXYZ(3,3*3,x+1,y+1,z+1)] = val;

    //TRCP(("%g %s%s%s", val, x==1?"\n":"", x==1&&y==1?"\n":"",x==1&&y==1&&z==1?"\n":""));
  }

  /*TRC(("%s: X:%d-%d, Y:%d-%d, Z:%d-%d \n",UT_FUNCNAME,
	_domainBcFlags[0][0],_domainBcFlags[0][1],
	_domainBcFlags[1][0],_domainBcFlags[1][1],
	_domainBcFlags[2][0],_domainBcFlags[2][1] ));*/

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
// Fast linear interpolation of scalar field
float MultiresSparseGrid::interpolateScalarFast2( 
    	 int chanId,		
    	 Vec4f pos0, 		// Position on finest resolution
	 unsigned options,	// SBG::OPT_IPCORNER, SBG::OPT_IPGRIDALGN
	 Vec4f *pOutGrad
    	)
{
  #ifdef UT_ASSERT_LEVEL_2
  if(!( isValidForInterpolation( chanId ) ))
  {
    TRCERR(("chanId=%d\n",chanId));
    UT_ASSERT0(FALSE);
  }
  #endif

  UT_ASSERT2(!( options & OPT_IPCORNER ));

  float out=0;
  pos0 = v4f_zero_to<3>(pos0);
  Vec4f pos = pos0;
  // Clip coords to global domain grid
  pos = max(pos,_sg0->v4fDomMin());
  pos = min(pos,_sg0->v4fDomMax());

  Vec4i ipos = truncate_to_int(pos);

  #if 0
  if(!_sg0->inRange(ipos))
  {
    TRCERR(("ipos=%d,%d,%d out of range\n",ipos[0],ipos[1],ipos[2]));
  }
  #endif

  int ix0=vget_x(ipos),
      iy0=vget_y(ipos),
      iz0=vget_z(ipos);
  UT_ASSERT2(_sg0->inRange(ix0,iy0,iz0));
  int bx,by,bz;
  _sg0->getBlockCoords(ix0,iy0,iz0, bx,by,bz);
  int bid = _sg0->getBlockIndex(bx,by,bz);

  BlockInfo *bi = getBlockInfo0(bid);
  int level = bi->level;
  UT_ASSERT2(level>=0&&level<_nLevels);

  Level *lev = &_sparseGrids[level];
  SparseGrid<float> *sg = getFloatChannel(chanId,level);
  SparseGrid<CellFlags> *sgFlags = lev->cellFlags[0];

  #if 0
  if(!sg || !sgFlags) 
  {
    TRCERR(("%p %p level=%d, ipos=%d,%d,%d\n",sg,sgFlags,level,ipos[0],ipos[1],ipos[2]));
  }
  #endif

  pos = pos0;  // pass un-clipped coords to SBG::interpolate
  	       // to enable domain-BC handling there

  // Convert coordinates to current level grid
  pos *= lev->auxInvScale;  // 1.0f/(1<<level)

  ipos >>= level;	
  Vec4i ivpos = ipos & (sgFlags->bsx()-1);
  int ivx=vget_x(ivpos),
      ivy=vget_y(ivpos),
      ivz=vget_z(ivpos),    
      vid = sgFlags->getVoxelInBlockIndex(ivx, ivy, ivz);

  #if 0
  if(!sgFlags->getBlock(bid)) 
  {
    TRCERR(("__ level=%d, ipos=%d,%d,%d\n",level,ipos[0],ipos[1],ipos[2]));
  }   
  #endif

  CellFlags cellFlags = sgFlags->getValueInBlock0(bid,vid);

  //IACA_END

  unsigned ipOpt=OPT_IPCORNER;
  if(!(bi->flags & BLK_DOM_BORDER))      		
  {
    ipOpt |= OPT_IP_NO_BORDER_CHECK;	  
  }

  UT_ASSERT2(options & (OPT_FAST_DISCONTINUOUS_FINECOARSE ))

  SparseGrid<float> *sgAct = sg;
  Vec4f posAct = pos;

  if(cellFlags & CELL_FINE_COARSE)
  {
    UT_ASSERT2(level<_nLevels-1);
    SparseGrid<float> *sgLo = getFloatChannel(chanId,level+1);
    UT_ASSERT2((!sgLo)||(sgLo&&sgLo->getBlock(bid)));

    #if 0
    if(!((!sgLo)||(sgLo&&sgLo->getBlock(bid)))) 
    {
      TRCERR(("sgLo=%p level=%d, ipos=%d,%d,%d\n",sgLo,level,ipos[0],ipos[1],ipos[2]));
    }   
    #endif

    sgAct = sgLo;
    posAct = pos*0.5f;
  }

  if(!(options & OPT_IPCORNER) ) posAct -= Vec4f(.5f,.5f,.5f,.0f);
  ipOpt |= OPT_IPCORNER;
  
  if(pOutGrad)
  {
    float grad_[4];
    sgAct->interpolateLinearWithGradient<1,1>(posAct,ipOpt,
					      &out,grad_);
    grad_[3] = 0.0f;
    (*pOutGrad).load(grad_);
  }
  else
  {
    out = sgAct->interpolateFloatFast(posAct,ipOpt);
  }

  return out;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template<int ipolType_T> 
void MultiresSparseGrid::interpolate( 
    	 int chanId,
    	 float *pos0, 		// Position on finest resolution
	 unsigned options,	// SBG::OPT_IPCORNER, SBG::OPT_IPGRIDALGN
	 Vec3Float& dataOut,      // OUT
	 bool doOnlyLevel0
    	)
{
  SparseGrid<Vec3Float> *sg0 = _sg0;

  #ifdef UT_ASSERT_LEVEL_2
  if(ipolType_T!=IP_NEAREST && ipolType_T!=IP_LINEAR && chanId!=CH_DIST_FINECOARSE)
  {
    if(!(_validForInterpolation & (((int64_t)1)<<chanId)))
    {
      TRCERR(("ipolType_T=%d, chanId=%d\n",ipolType_T,chanId));
      UT_ASSERT0( FALSE );
    }
  }
  #endif

  if( options & OPT_IPEXTRAPOL )
  {
    float beps=1e-2;
    G3D_CLIP_CONST(beps, beps, beps, 
		   sg0->sx()-beps, sg0->sy()-beps, sg0->sz()-beps,
		   pos0[0],pos0[1],pos0[2]);
  }

  /*if(chanId!=CH_VEC3_1)
  {
    UT_ASSERT0(ipolType_T==IP_NEAREST);
  }*/

  float pos[4];
  for(int k=0;k<3;k++) pos[k] = pos0[k];

  // 
  // Determine block id and level
  //
  int bx, by, bz;
  LongInt bid;
  {
    int ix0 = floorf(pos[0]),
	iy0 = floorf(pos[1]),
    	iz0 = floorf(pos[2]);
   
    MT_CLAMP(ix0,0,sg0->sx()-1);
    MT_CLAMP(iy0,0,sg0->sy()-1);
    MT_CLAMP(iz0,0,sg0->sz()-1);
   

    sg0->getBlockCoords( ix0, iy0, iz0, 
			 bx, by, bz);
    UT_ASSERT2(sg0->blockCoordsInRange(bx,by,bz));


    bid = sg0->getBlockIndex(bx, by, bz);
  }


  int level=0;
  if(!doOnlyLevel0)
  {
    level = _blockmap[0][bid].level;
  }
  else
  {
    UT_ASSERT0(ipolType_T==IP_NEAREST);
  }


  UT_ASSERT2(level>=0&&level<_nLevels);
  Level *lev= &_sparseGrids[level];

  SparseGrid<Vec3Float> *sgVel;
  SparseGrid<float> *sgFloat;

  GET_IP_CHANNEL( chanId, level,
		  sgVel, sgFloat );

  int isStaggered = TRUE && sgVel && 
                     chanId != CH_VEC3_4 ;
  //
  // Convert coordinates to current level grid
  //
  for(int k=0;k<3;k++) pos[k] *= 1.0f/(1<<level);

  //UT_ASSERT0(!(options & OPT_IPCORNER));

  if(ipolType_T==IP_NEAREST)
  {
    if(sgVel) dataOut = sgVel->getValue(floorf(pos[0]),floorf(pos[1]),floorf(pos[2]));
    if(sgFloat) dataOut[0] = sgFloat->getValue(floorf(pos[0]),floorf(pos[1]),floorf(pos[2]));
    return;
  }

  //UT_ASSERT0( UT_IMPLIES(sgVel, !isStaggered) );

  unsigned ipOpt = 
    options & ( OPT_IPCORNER | OPT_IPBC_NEUMAN | OPT_IPBC_DIRICHLET | 
		OPT_IPBC_DIRICHLET_OLD | OPT_IP_FADE_WEIGHTS | 
      		OPT_IP_VCOMP_X | OPT_IP_VCOMP_Y | OPT_IP_VCOMP_Z );
  ipOpt |= OPT_IPCHKBLK;

  if(!((_blockmap[0][bid].flags & BLK_DOM_BORDER))) 
    ipOpt |= OPT_IP_NO_BORDER_CHECK;
    
  int stenWidth= ipolType_T==IP_LINEAR || ipolType_T==IP_NEAREST ? 2 : 4;
  UT_ASSERT0(stenWidth>0);
  int dmax = stenWidth-1;

  /*if( fabs(pos[0]-95.5)<MT_NUM_EPS &&
      fabs(pos[1]-63.5)<MT_NUM_EPS &&
      fabs(pos[2]-95.5)<MT_NUM_EPS )
  {
    TRCP(("DBG\n"));
  }*/

  //if(( lev->distFineCoarse->getBlock(bx,by,bz) ) && ipolType_T!=IP_NEAREST)
  if  (((_blockmap[0][bid].flags & BLK_FINE_COARSE) && ipolType_T!=IP_NEAREST) 
    || chanId==CH_DIST_FINECOARSE)
    
  {
    /*if(level==1 && (int)(pos[0])==33 && (int)(pos[1])==15 && (int)(pos[2])==0)
    {
      TRCP(("DBG\n"));
    }*/

//    double dist = lev->distFineCoarse->interpolateToFloat<ipolType_T>(pos,0) / 1024.0;
    double dist = 
      lev->distFineCoarse->interpolateToFloat<IP_LINEAR>(
	  					pos,
						OPT_IPBC_NEUMAN
						) / 1024.0;
  double 
     alphaFine = MtLinstep( (stenWidth/2-(isStaggered?.0f:.5f))*MT_SQRT3*1.01, 
			   0.99*(_dTransRes+2)/MT_SQRT3,
			  // 0.99*(_dTransRes+2-isStaggered?1.f:0.f)/MT_SQRT3,
			   dist);    
			

    if(chanId==CH_DIST_FINECOARSE)
    {
      dataOut[0]=dist;
      //dataOut[0]=alphaFine;
      return;
    }

    dataOut = 0;
    if(alphaFine>0.00001)
    {
      int ipOpt2 = ipOpt; 
      if(level==0 && (options & OPT_IPGRIDALGN)) ipOpt2 |= OPT_IPGRIDALGN;

      if(sgVel) dataOut = sgVel->interpolateVec3Float<ipolType_T>(dmax,pos,-1,ipOpt2);

      if(sgFloat) 
	dataOut[0] = sgFloat->interpolate<ipolType_T>(pos,ipOpt2);
    }
    if(alphaFine<.9999)
    {
      UT_ASSERT0(level<_nLevels-1);

      SparseGrid<Vec3Float> *sgVelLo;
      SparseGrid<float> *sgFloatLo;
      GET_IP_CHANNEL( chanId, level+1,
		   sgVelLo, sgFloatLo );

      UT_ASSERT2((!sgVelLo)||(sgVelLo&&sgVelLo->getBlock(bx,by,bz)));

      float posLo[3] = { pos[0]*0.5f, pos[1]*0.5f, pos[2]*0.5f };

      if(sgVelLo)
      {
        Vec3Float velLo = sgVelLo->interpolateVec3Float<ipolType_T>(dmax,posLo,-1,ipOpt);	  		
        dataOut = dataOut*alphaFine + velLo*(1-alphaFine);
      }
      if(sgFloatLo)
      {
        float dataLo = sgFloatLo->interpolate<ipolType_T>(posLo,ipOpt);
        dataOut[0] = dataOut[0]*alphaFine + dataLo*(1-alphaFine);
      }      
    }
    return;
  }
  else
  {
    if(chanId==CH_DIST_FINECOARSE)
    {
      dataOut[0]=0;//_dTransRes;
      return;
    }

    int ipOpt2 = ipOpt; 
    if(level==0 && (options & OPT_IPGRIDALGN)) ipOpt2 |= OPT_IPGRIDALGN;
    if(sgVel) dataOut = sgVel->interpolateVec3Float<ipolType_T>(dmax,pos,-1,ipOpt2);
    if(sgFloat) dataOut[0] = sgFloat->interpolate<ipolType_T>(pos,ipOpt2);
    return;
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::downsampleGhostBlockFloat( 
    			  float **dataHaloBlocks,
    			  int bid, int level, int levelMg, // hi-res source level
    			  int iChanDataFloat )
{
  float *krnDown = _krnDownSample4x4x4[0/*cell centerd*/];
  float *dataHaloBlock = dataHaloBlocks[0];

  int bx,by,bz;
  _sg0->getBlockCoordsById(bid, bx,by,bz);

  SparseGrid<float> *sgFloatLo = getFloatChannel(iChanDataFloat,level+1,levelMg);
  float *dataFloatLo = sgFloatLo->getBlockDataPtr(bid,1);

  int bsxHi = _sg0->bsx() >> level,
      bsxLo = bsxHi/2;
  
  int hsx = bsxHi+2,
      hsx2 = hsx*hsx;

  UT_ASSERT0(dataFloatLo);

  for(int vz=0;vz<bsxLo;vz++)
  for(int vy=0;vy<bsxLo;vy++)
  for(int vx=0;vx<bsxLo;vx++)
  {
    int vid = sgFloatLo->getVoxelInBlockIndex(vx,vy,vz);

    Vec4f S0(0.0f),S1(0.0f),S2(0.0f),S3(0.0f), // mult. accumulators
	  W,X;
    int vxh = 2*vx-1 + 1,
	vyh = 2*vy-1 + 1,
	vzh = 2*vz-1 + 1;

    int hid = MT_GXYZ(hsx,hsx2, vxh,vyh,vzh);
    for(int k=0,iKrn=0; k<4; k++, iKrn+=16, hid += hsx2 )
    {
      float *pW=krnDown+iKrn,
	    *pX=dataHaloBlock+hid;
      S0 += Vec4f().load(pW+0) * Vec4f().load(pX+0);
      S1 += Vec4f().load(pW+4) * Vec4f().load(pX+hsx);
      S2 += Vec4f().load(pW+8) * Vec4f().load(pX+2*hsx);
      S3 += Vec4f().load(pW+12) * Vec4f().load(pX+3*hsx);

#if 0
    if( UtGlobDbgCheckCount==DBG_CELL_DBGCOUNT &&
	bid == DBG_CELL_BID &&
	vx==DBG_CELL_VX && vy==DBG_CELL_VY && vz==DBG_CELL_VZ)
    {
      Vec4f X0=Vec4f().load(pX+0),
	    X1=Vec4f().load(pX+hsx),
	    X2=Vec4f().load(pX+2*hsx),
	    X3=Vec4f().load(pX+3*hsx);

      Vec4f W0=Vec4f().load(pW+0),
	    W1=Vec4f().load(pW+4),
	    W2=Vec4f().load(pW+8),
	    W3=Vec4f().load(pW+12);

      TRCP(("%g:%g %g:%g %g:%g %g:%g\n",W0[0],X0[0],W0[1],X0[1],W0[2],X0[2],W0[3],X0[3]));
      TRCP(("%g:%g %g:%g %g:%g %g:%g\n",W1[0],X1[0],W1[1],X1[1],W1[2],X1[2],W1[3],X1[3]));
      TRCP(("%g:%g %g:%g %g:%g %g:%g\n",W2[0],X2[0],W2[1],X2[1],W2[2],X2[2],W2[3],X2[3]));
      TRCP(("%g:%g %g:%g %g:%g %g:%g\n",W3[0],X3[0],W3[1],X3[1],W3[2],X3[2],W3[3],X3[3]));
    }
#endif

    }

    float val = horizontal_add(S0+S1+S2+S3);
    UT_ASSERT2( std::isfinite(val));
    dataFloatLo[vid] = val;
  }
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::downsampleGhostBlockSOAtoVec3Float( 
    			  float **dataHaloBlocks,
    			  int bid, int level, int levelMg, // hi-res source level
    			  int iChanDataVec )
{
  SparseGrid<Vec3Float> *sgVelLo=NULL;
  Vec3Float *dataVelLo=NULL;
  int bx,by,bz;
  _sg0->getBlockCoordsById(bid, bx,by,bz);
  
  sgVelLo = getVecChannel(iChanDataVec,level+1,levelMg);      
  dataVelLo = sgVelLo->getBlockDataPtr(bid,1);

  int bsxHi = _sg0->bsx() >> level,
      bsxLo = bsxHi/2;

  for(int iVelCmp=0; iVelCmp<3; iVelCmp++)
  {
    float *krnDown = _krnDownSample4x4x4[1+iVelCmp];
    int hsx = bsxHi+2,
	hsx2 = hsx*hsx;

    for(int vz=0;vz<bsxLo;vz++)
    for(int vy=0;vy<bsxLo;vy++)
    for(int vx=0;vx<bsxLo;vx++)
    {
      int vid = sgVelLo->getVoxelInBlockIndex(vx,vy,vz);

      Vec4f S0(0.0f),S1(0.0f),S2(0.0f),S3(0.0f), // mult. accumulators
	    W,X;
      int vxh = 2*vx-1 + 1,
	  vyh = 2*vy-1 + 1,
	  vzh = 2*vz-1 + 1;

      int hid = MT_GXYZ(hsx,hsx2, vxh,vyh,vzh);
      for(int k=0,iKrn=0; k<4; k++, iKrn+=16, hid += hsx2 )
      {
	float *pW=krnDown+iKrn,
	      *pX=dataHaloBlocks[iVelCmp]+hid;
	S0 += Vec4f().load(pW+0) * Vec4f().load(pX+0);
	S1 += Vec4f().load(pW+4) * Vec4f().load(pX+hsx);
	S2 += Vec4f().load(pW+8) * Vec4f().load(pX+2*hsx);
	S3 += Vec4f().load(pW+12) * Vec4f().load(pX+3*hsx);
      }
      float val = horizontal_add(S0+S1+S2+S3);
      UT_ASSERT2( std::isfinite(val));
      dataVelLo[vid][iVelCmp] = val;
    }
  } // iVelCmp
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::resetEMPTYFlags( std::vector<int> *blockList )
{
  TRC(("%s\n",UT_FUNCNAME));
  ThrRunParallel<int> ( blockList ? blockList->size() : _nBlocks, nullptr,   
    [&](int &tls, int tid, int ibid)
    {
      int bid = blockList ? (*blockList)[ibid] : ibid;
      BlockInfo *bi=getBlockInfo(bid,0);
      FLAG_RESET( bi->flags, BLK_HAS_UNDEFINED | BLK_HAS_DEFINED );
      //if(!(bi->flags & BLK_FIXED)) return;
      if(!(bi->flags & BLK_EXISTS)) return;
      SparseGrid<CellFlags> *sg = getFlagsChannel0( bi->level );
      UT_ASSERT2( sg->isValueBlock(bid ));
      CellFlags *dataFlags = sg->getBlockDataPtr(bid);
      for( int i=0;i<sg->nVoxelsInBlock();i++) 
      {
        FLAG_RESET( dataFlags[i], CELL_EMPTY_C | 
	    			  CELL_EMPTY_U | CELL_EMPTY_V | CELL_EMPTY_W );
      }          
    } );
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::prepare( int iChanData,
				   int setSideObstacles,
				   int doGhostCellsOnly
				)
{
  int rcThis=0;
  USE(rcThis);
  int levelMg = 0;

  doGhostCellsOnly = TRUE;

  int iChanDataVec=0,
      iChanDataFloat=0;

  if(iChanData)
  {
    if(isVecChannel(iChanData))
      iChanDataVec = iChanData;
    else
      iChanDataFloat = iChanData;
  }

  TRC3(("%s '%s' %d chan=%d (%d %d)\n",UT_FUNCNAME,_name,
	doGhostCellsOnly,
	iChanData,iChanDataFloat,iChanDataVec));

  UT_ASSERT0(!( iChanDataFloat && iChanDataVec ));

  UtTimer tm;
  TIMER_START(&tm);

  //
  // Prepare ghost blocks
  //    
  std::unique_ptr<HaloBlockSet> haloBlocks( nullptr );
  haloBlocks.reset( 
      new HaloBlockSet( this, 1 /* halo size*/, nMaxThreads,
			       iChanDataVec ? OPT_HALO_VEC3 : 0) );

  // Loop over each block belonging to current level  

  for(int level0=0; level0<_nLevels; level0++)
  {
    std::vector<int> *blocks = &_blocksFineCoarsePerLevel[level0];    
    ThrRunParallel<int> ( blocks->size(), nullptr, 
      [&](int &tls, int tid, int bidx)
      {
	UT_ASSERT2(bidx>=0&&bidx<(int)blocks->size());
	int bid = (*blocks)[bidx];
	BlockInfo *bi = getBlockInfo(bid,levelMg);

	int level = bi->level;

	UT_ASSERT2(level == level0);
	UT_ASSERT2(bi->flags & BLK_FINE_COARSE);

	UT_ASSERT0( level<_nLevels-1);
	if( iChanDataFloat )
	{
	  if(!(bi->flags & BLK_EXISTS)) 
	  {
	    getFloatChannel(iChanDataFloat,level+1)->setEmptyBlock(bid);
	    return;
	  }

	  float **dataHaloBlock = 
	    haloBlocks->fillHaloBlock_<float>( iChanDataFloat, bid, 0, tid,
				       OPT_BC_COARSE_LEVEL );
	  downsampleGhostBlockFloat( dataHaloBlock, bid, level, levelMg,
				     iChanDataFloat );
	}
	else if( iChanDataVec )
	{
	  if(!(bi->flags & BLK_EXISTS)) 
	  {
	    getVecChannel(iChanDataVec,level+1)->setEmptyBlock(bid);
	    return;
	  }

	  float **dataHaloBlocks =
	    haloBlocks->fillHaloBlock_<Vec3Float,3>( iChanDataVec, bid, 0, tid,
							OPT_BC_COARSE_LEVEL );
	  downsampleGhostBlockSOAtoVec3Float( dataHaloBlocks, bid, level, levelMg,
					      iChanDataVec );
	}
	else
	{
	  UT_ASSERT0(FALSE);
	}
      } );
  }
 
  // Cleanup temp. ressources
  resetChannel( CH_CELL_FLAGS_TMP );

  _validForInterpolation |= (((int64_t)1) << iChanData);
  
  // 
  TIMER_STOP(&tm);
  TRC(("%s %d: CPU = %.2f sec, %.0f act. voxels/sec)\n",
    UT_FUNCNAME,iChanData,
    (double)TIMER_DIFF_MS(&tm)/1000.,
    (double)_nActCells[0]/(double)(TIMER_DIFF_MS(&tm)/1000.0)));
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::prepareGhostBlocks( int iChanData )
{
  return prepare( iChanData, 1, TRUE );
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static void copyBlockFlagsFromToFloatChan( 
    			MultiresSparseGrid *msg,
    			int level, int levelMg,
			int bid,
    			int chSrc, int chDst,
    			int doParallel )
{
  UT_ASSERT2(sizeof(CellFlags)==sizeof(uint16_t));
  if(msg->isFlagsChannel(chSrc))
  {
    SparseGrid<CellFlags> *sgSrc = msg->getFlagsChannel(chSrc, level, levelMg );
    SparseGrid<float> *sgDst = msg->getFloatChannel(chDst, level, levelMg ); 
    CellFlags *dataSrc = sgSrc->getDataPtrGen_(bid,0,0,1);
    float *dataDst = sgDst->getDataPtrGen_(bid,1,0,1);
    LongInt n = sgDst->isDenseGrid() ? sgDst->nTotActVoxels() :
      		sgDst->nVoxelsInBlock();
    if(doParallel)
    {
      #pragma omp parallel for 
      for(LongInt i=0;i<n;i++) dataDst[i] = dataSrc[i];
    }
    else for(LongInt i=0;i<n;i++) dataDst[i] = dataSrc[i];
  }
  else
  {
    SparseGrid<float> *sgSrc = msg->getFloatChannel(chSrc, level, levelMg );
    SparseGrid<CellFlags> *sgDst = msg->getFlagsChannel(chDst, level, levelMg ); 
    float *dataSrc = sgSrc->getDataPtrGen_(bid,0,0,1);
    CellFlags *dataDst = sgDst->getDataPtrGen_(bid,1,0,1);
    LongInt n = sgDst->isDenseGrid() ? 
      		  sgDst->nTotActVoxels() :
      		  sgDst->nVoxelsInBlock();
    if(doParallel)
    {
      #pragma omp parallel for
      for(LongInt i=0;i<n;i++) dataDst[i] = dataSrc[i];    
    }
    else for(LongInt i=0;i<n;i++) dataDst[i] = dataSrc[i];
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::copyChannelFromMSBG( 
    				      const MSBG::MultiresSparseGrid *msgSrc,
    				      int chSrc, int chDst )
{
  TRC(("%s: %d %d\n",UT_FUNCNAME,(int)chSrc,(int)chDst));
  UT_ASSERT0( dimEqual(msgSrc));

  UtTimer tm;
  TIMER_START(&tm);

  resetChannel(chDst);
  prepareDataAccess(chDst);

  int isVecChan = isVecChannel(chSrc),
      isPSFloatChan = isPSFloatChannel(chSrc);

  ThrRunParallel<int> ( nBlocks(), nullptr,
    [&](int &tls, int tid, int bid)
    { 
      BlockInfo *bi = getBlockInfo(bid);
      if(isVecChan)
      {
	SparseGrid<Vec3Float> *sgSrc = msgSrc->getVecChannel(chSrc, bi->level),
			      *sgDst = getVecChannel(chDst, bi->level ); 
	UT_ASSERT0( SBG::dimEqual<Vec3Float>( sgSrc, sgDst ) );

	Vec3Float *dataDst = sgDst->getBlockDataPtr(bid,TRUE,TRUE),
		   *dataSrc = sgSrc->getBlockDataPtr(bid);
	for(LongInt i=0; i<sgDst->nVoxelsInBlock(); i++) 
	{
	  dataDst[i] = dataSrc[i];
	}
      }
      else if(isPSFloatChan)
      {
	SparseGrid<PSFloat> *sgSrc = msgSrc->getPSFloatChannel(chSrc, bi->level),
			    *sgDst = getPSFloatChannel(chDst, bi->level ); 
	UT_ASSERT0( SBG::dimEqual<PSFloat>( sgSrc, sgDst ) );

	PSFloat *dataDst = sgDst->getBlockDataPtr(bid,TRUE,TRUE),
		   *dataSrc = sgSrc->getBlockDataPtr(bid);
	for(LongInt i=0; i<sgDst->nVoxelsInBlock(); i++) 
	{
	  dataDst[i] = dataSrc[i];
	}
      }
      else
      {
	UT_ASSERT0(FALSE); // TODO
      }
    } );
  invalidateForInterpolation(chDst);

  TIMER_STOP(&tm);
  TRC(("%s: CPU = %.2f sec, %.0f cells/sec)\n",UT_FUNCNAME,
	(double)TIMER_DIFF_MS(&tm)/1000.,
	(_nActCells[0])/(double)(TIMER_DIFF_MS(&tm)/1000.0)));

}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::copyChannel( int chSrc, int chDst, 
        			      int includeDomBorder,
    				      int levelMax,
    				      int levelMg,
				      unsigned options,
				      std::vector<int> *blockList
				      )
{
  TRC3(("%s %d %d\n",UT_FUNCNAME,(int)chSrc,(int)chDst));
  SparseGrid<Vec3Float> *sg0 = _sparseGrids[0].vec3_1[0];
  UT_ASSERT0(includeDomBorder);

  prepareDataAccess( chDst, -1, levelMg );

  if(levelMg>=_nLevels)  // special case: dense (coarse) grid
  {
    UT_ASSERT0(nBlocks(levelMg)==1);
    if(isFlagsChannel(chSrc)||isFlagsChannel(chDst))
    {
      if(!isFlagsChannel(chSrc)||!isFlagsChannel(chDst))
      {
	copyBlockFlagsFromToFloatChan(this,levelMg,levelMg,-1,chSrc,chDst,TRUE);
      }
      else
      {
	SparseGrid<CellFlags> *sgSrc = getFlagsChannel(chSrc, levelMg, levelMg ),
			      *sgDst = getFlagsChannel(chDst, levelMg, levelMg ); 
	CellFlags *dataDst = sgDst->getDataPtrGen_(0,1,0,1),
		  *dataSrc = sgSrc->getDataPtrGen_(0,0,0,1);
	LongInt n = sgDst->nTotActVoxels();
	#pragma omp parallel for if( n > N_CELLS_PARALLEL)
	for(LongInt i=0;i<n;i++) dataDst[i] = dataSrc[i];
      }
    }
    else if(isVecChannel(chSrc))
    {
      SparseGrid<Vec3Float> *sgDst = getVecChannel(chDst,levelMg,levelMg),
			*sgSrc = getVecChannel(chSrc,levelMg,levelMg);
      Vec3Float *dataDst = sgDst->getDataPtrGen_(0,1,0,1),
	    *dataSrc = sgSrc->getDataPtrGen_(0,0,0,1);
      LongInt n = sgDst->nTotActVoxels();
      #pragma omp parallel for if( n > N_CELLS_PARALLEL)
      for(LongInt i=0;i<n;i++) dataDst[i] = dataSrc[i];
    }
    else if( isPSFloatChannel(chSrc ) || isPSFloatChannel(chDst) )
    {
      #ifdef SOLVE_PRESSURE_DOUBLE
      UT_ASSERT2( isPSFloatChannel(chSrc) && isPSFloatChannel( chDst )) 
      #endif
      SparseGrid<PSFloat> *sgDst = getPSFloatChannel(chDst,levelMg,levelMg),
			*sgSrc = getPSFloatChannel(chSrc,levelMg,levelMg);
      PSFloat *dataDst = sgDst->getDataPtrGen_(0,1,0,1),
	      *dataSrc = sgSrc->getDataPtrGen_(0,0,0,1);
      LongInt n = sgDst->nTotActVoxels();
      ThrRunParallel<int>(n, nullptr, 
        [&](int &tls, int tid, size_t i)
	{
	  dataDst[i] = dataSrc[i];
	} );
    }
    else
    {
      SparseGrid<float> *sgDst = getFloatChannel(chDst,levelMg,levelMg),
			*sgSrc = getFloatChannel(chSrc,levelMg,levelMg);
      float *dataDst = sgDst->getDataPtrGen_(0,1,0,1),
	    *dataSrc = sgSrc->getDataPtrGen_(0,0,0,1);
      LongInt n = sgDst->nTotActVoxels();
      ThrRunParallel<int>(n, nullptr, 
        [&](int &tls, int tid, size_t i)
	{
	  dataDst[i] = dataSrc[i];
	} );
    }
    return;
  } // dense grid

  // 
  // Sparse grid
  //
  int nActBlocks = blockList ? blockList->size() : sg0->nBlocks();

  ThrRunParallel<int>( nActBlocks,
    nullptr, // initialize
    // Run parallel
    [&](int &tls, int tid, int ibid) 
    {
      int bid = blockList ? (*blockList)[ibid] : ibid;
      BlockInfo *bi = getBlockInfo(bid,levelMg);

      if((options & OPT_EXISTING_BLOCKS_ONLY) && 
	  !(bi->flags & BLK_EXISTS) ) return;

      int level = MAX( _blockmap[0][bid].level, levelMg);
      UT_ASSERT0(level==bi->level);
      
      if(isFlagsChannel(chSrc)||isFlagsChannel(chDst))
      {
	UT_ASSERT2(includeDomBorder);
	if(!isFlagsChannel(chSrc)||!isFlagsChannel(chDst))
	{
	  copyBlockFlagsFromToFloatChan(this,level,levelMg,bid,chSrc,chDst,TRUE);
	}
	else
	{
	  SparseGrid<CellFlags> *sgSrc = getFlagsChannel(chSrc, level, levelMg ),
	    *sgDst = getFlagsChannel(chDst, level, levelMg );

	  if(!(sgSrc&&sgDst)) return;

	  CellFlags *dataDst = sgDst->getBlockDataPtr(bid,TRUE,0),
		    *dataSrc = sgSrc->getBlockDataPtr(bid);

	  for(int i=0; i<sgDst->nVoxelsInBlock(); i++) dataDst[i] = dataSrc[i];
	}
      }
      else if(isVecChannel(chSrc))
      {
	UT_ASSERT0(isVecChannel(chDst));
	SparseGrid<Vec3Float> *sgSrc = getVecChannel(chSrc, level, levelMg),
		     	      *sgDst = getVecChannel(chDst, level, levelMg );
	if(sgSrc->isEmptyBlock(bid)) 
	{
	  sgDst->setEmptyBlock(bid); 
	  return;
	}
	else if(sgSrc->isFullBlock(bid)) 
	{
	  sgDst->setFullBlock(bid);
	  return;
	}
	else if(sgSrc->isInvalidBlock(bid)) 
	{
	  sgDst->setInvalidBlock(bid);
	  return;
	}

	Vec3Float *dataDst = sgDst->getBlockDataPtr(bid,TRUE,0),
	           *dataSrc = sgSrc->getBlockDataPtr(bid);
	for(LongInt i=0; i<sgDst->nVoxelsInBlock(); i++) 
	{
	  dataDst[i] = dataSrc[i];
	}
      }
      else if(isPSFloatChannel(chSrc))
      {
	SparseGrid<PSFloat> *sgSrc = getPSFloatChannel(chSrc, level, levelMg ),
			    *sgDst = getPSFloatChannel(chDst, level, levelMg ); 

	if(!(sgSrc&&sgDst)) return;

	if(sgSrc->isEmptyBlock(bid)) 
	{
	  sgDst->setEmptyBlock(bid); 
	  return;
	}
	else if(sgSrc->isFullBlock(bid)) 
	{
	  sgDst->setFullBlock(bid);
	  return;
	}
	else if(sgSrc->isInvalidBlock(bid)) 
	{
	  sgDst->setInvalidBlock(bid);
	  return;
	}
	
	PSFloat *dataSrc = sgSrc->getBlockDataPtr(bid);
	  
	if(!dataSrc) return;

	PSFloat *dataDst = sgDst->getBlockDataPtr(bid,TRUE,0);

	for(int i=0; i<sgDst->nVoxelsInBlock(); i++) dataDst[i] = dataSrc[i];
      }
      else
      {
	SparseGrid<float> *sgSrc = getFloatChannel(chSrc, level, levelMg ),
			  *sgDst = getFloatChannel(chDst, level, levelMg ); 

	if(!(sgSrc&&sgDst)) return;

	if(sgSrc->isEmptyBlock(bid)) 
	{
	  sgDst->setEmptyBlock(bid); 
	  return;
	}
	else if(sgSrc->isFullBlock(bid)) 
	{
	  sgDst->setFullBlock(bid);
	  return;
	}
	else if(sgSrc->isInvalidBlock(bid)) 
	{
	  sgDst->setInvalidBlock(bid);
	  return;
	}
	
	float *dataSrc = sgSrc->getBlockDataPtr(bid);
	  
	if(!dataSrc) return;

	float *dataDst = sgDst->getBlockDataPtr(bid,TRUE,0);

	if( _nActCells[levelMg] > 500000 )
	{
	  // Bypass cache for large data channels
	  UT_ASSERT2((sgDst->nVoxelsInBlock()&3)==0); 
	  for(int i=0;i<sgDst->nVoxelsInBlock();i+=4)
	  {
	    Vec4f X;
	    X.load_a(dataSrc+i);
	    //X.store(dataDstY+i);
	    _mm_stream_ps( dataDst+i, X );
	  }
	}
	else
	{
	  for(int i=0; i<sgDst->nVoxelsInBlock(); i++) dataDst[i] = dataSrc[i];
	}		
      }
    });  // ThrRunParallel

  invalidateForInterpolation(chDst);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::scalarChannelFromToVecComp( 
    			int dir,  // direction 1=vec2scalar, -1=scalar2vec
			int iVecComp,
    			int chScalar, int chVec
    )

{
  int levelMg = 0;
  TRC3(("%s dir=%d chScalar=%d chVec=%d\n",UT_FUNCNAME,dir,chScalar,chVec));

  UT_ASSERT0(levelMg==0);

  // 
  // Sparse grid
  //
  LONG iTask=0;
  #pragma omp parallel 
  {
    for(;;)
    {
      uint32_t bid = InterlockedIncrement((LONG volatile *)&iTask) - 1;
      if(bid>=_sg0->nBlocks()) break;
      int level = getBlockLevel(bid);
      
      SparseGrid<CellFlags> *sg = getFlagsChannel(CH_CELL_FLAGS,level);

      SparseGrid<float> *sgFloat = getFloatChannel(chScalar,level);
      float *dataFloat = sgFloat->getBlockDataPtr(bid,dir==1?1:0,0);

      SparseGrid<Vec3Float> *sgVec = getVecChannel(chVec,level);
      Vec3Float *dataVec = sgVec->getBlockDataPtr(bid,dir==-1?1:0,0);

      int n = sg->nVoxelsInBlock();
      if(dir==1) // vec2scalar
      {
	for(int i=0; i<n; i++) dataFloat[i] = dataVec[i][iVecComp];
      }
      else if(dir==-1)  // scalar2vec
      {
	for(int i=0; i<n; i++) dataVec[i][iVecComp] = dataFloat[i];
      }
      else
      {
	UT_ASSERT0(FALSE);
      }
    }
  }
  invalidateForInterpolation( dir==1?chScalar:chVec );
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::fillTestPatternChannel( 
    				     int typ,
				     int chanDst, 
				     int levelMg,
				     double rndSeed,
				     unsigned options
				     )
{
  // 
  TRC3(("%s %d %d %d\n",UT_FUNCNAME,typ,(int)chanDst,levelMg));

  if(isDenseLevel(levelMg))
  {
    if(isVecChannel(chanDst)) 
      getVecChannel(chanDst, levelMg,levelMg)->getDenseData();
    else
      getFloatChannel(chanDst, levelMg,levelMg)->getDenseData();
  }

  BlockGridInfo bgi( this, levelMg );

  // Loop over each block of levelMg
  #pragma omp parallel for schedule( static, 1 )
  for(int bid=0; bid<bgi._nBlocks; bid++)
  {
    int level = getBlockLevel(bid, levelMg );  
      
    SparseGrid<CellFlags> 
      *sgFlags = getFlagsChannel(CH_CELL_FLAGS,level,levelMg);

    SparseGrid<float> *sg = isVecChannel(chanDst) ? NULL : 
      				getFloatChannel(chanDst, level,levelMg);
    SparseGrid<Vec3Float> *sgVec = isVecChannel(chanDst) ? 
      				getVecChannel(chanDst, level,levelMg) : NULL;

    CellFlags *dataFlags = sgFlags->getDataPtrGen_(bid,0,0);
    USE(dataFlags);

    float *data = sg ? sg->getDataPtrGen_(bid,1,0):NULL;
    Vec3Float *dataVec  = sgVec ? sgVec->getDataPtrGen_(bid,1,0):NULL;

    // Loop over each voxel in block
    BlockAccessor bac(&bgi,sgFlags,bid,level);
    for(int vz=0;vz<bac._bsz;vz++)
    for(int vy=0;vy<bac._bsy;vy++)
    for(int vx=0;vx<bac._bsx;vx++)
    {
      int vid = bac.getCellIndex( vx, vy, vz );
      int ix,iy,iz;
      bac.getGridCoords( vx, vy, vz,
			 ix, iy, iz );	  
      Vec3Float pos( ix, iy, iz );
      pos = (pos+0.5f)*(1<<level);

      //CellFlags flags = dataFlags[vid];
      float f=0;
      Vec3Float fv(0.0f);
      if(typ==0||typ==1)
      {
        LongInt idx = level<_nLevels ? bid*1000+vid : vid;
	f = (float)idx;
      }
      else if(typ==2)
      {
        uint32_t idx = level<_nLevels ? bid*1000+vid : vid;
	UT_FAST_HASH_32F(idx, f);
      }
      else if(typ==3) // Fractal
      {
	if((options&OPT_ALL_CELLS) || CELL_IS_FLUID_(dataFlags[vid]))
	{
	  //UT_ASSERT0(sgVec);

	  double 
		 oct_offs = 2,
		 //oct_offs = -2,	       	       
		 beta = 1.97,
		 fscale = _sg0->sxyzMax()/2.,   // feature scale in grid cells
		 minscale =0.25,
		 oct_leadin = 1,
		 alpha = 1.5,
		 do_turb = 0;

	  double octaves = MAX( MAX( logf( fscale / minscale ) / MT_LOG2 , 1 ) + oct_offs , 1),
		       freq = 1./( fscale * powf( 2, oct_leadin ) );

	  //octaves = 2;


	  for(int k=0;k<3;k++)
	  {
	    pos.x = (ix+_goffVel[k][0])*(1<<level);
	    pos.y = (iy+_goffVel[k][1])*(1<<level);
	    pos.z = (iz+_goffVel[k][2])*(1<<level);

	    double seed = k*5.123;
	    double p[3] = {(double)pos.x*freq+seed*7.7,
			     (double)pos.y*freq+seed*3.141,
			     (double)pos.z*freq+seed*11.1};

	    float f = PNS_Fractal( 
			p, 3, alpha, beta, 
			octaves, oct_leadin,
			do_turb, 0.6, 1, 0,
			0, 
			PNS_NTYP_FAST, rndSeed, 0);
	    fv[k] = f;
	    UT_ASSERT0(std::isfinite(fv[k]));
	  }
	  f = fv[0];
	}
	else
	{
	  for(int k=0;k<3;k++) fv[k] = _obsVelocity[k];
	}
      }
      else if(typ==4)
      {
	for(int k=0;k<3;k++)
	{
	  {
	    pos.x = (ix+_goffVel[k][0])*(1<<level);
	    pos.y = (iy+_goffVel[k][1])*(1<<level);
	    pos.z = (iz+_goffVel[k][2])*(1<<level);
	  }

	  pos /= _sg0->sxyzMin();
	  double r = norm(pos);
	  fv[k]  = 0.5*sin(r*20);
	}
      }
      else if(typ==5)
      {
        pos /= _sg0->sxyzMax();
	f = norm( pos - Vec3Float(0.3,0.3,0.5) ) - 0.2;
      }

      if(sgVec)
	dataVec[vid] = fv;
      else
        data[vid] = f;           
    }
  }
  invalidateForInterpolation(chanDst);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
double MultiresSparseGrid::compareChannel( 
		       int chanA,  // first channel to compare (reference)
		       int chanB,  // second channel to compare
		       int levelMg,
		       int doVisu,
		       int doOnlyLevel0,
		       int chVisDiff,
		       double errThr  // (absolute) error threshold
				     )
{
  TRC(("%s %d %d %d\n",UT_FUNCNAME,(int)chanA,(int)chanB,levelMg));
  double diffMax = -1e20;
  Vec3Int iposMax(-1);
  int levelMax=-1;

  if(doVisu)
  {
    if(chVisDiff)
    {
      resetChannel(chVisDiff,false,levelMg);
      prepareDataAccess( chVisDiff, -1, levelMg );
    }

    #if 0
    static int pnl1=-1,pnl2=-1;
    visualizeSlices( MiGetAuxPanel2(&pnl1,"compareChannels_A"), chanA, 
      SBG::IP_NEAREST, NULL, VIS_OPT_TEST );
    visualizeSlices( MiGetAuxPanel2(&pnl2,"compareChannels_B"), chanB, 
      SBG::IP_NEAREST, NULL, VIS_OPT_TEST );
    #endif
  }

  BlockGridInfo bgi( this, levelMg );
  #pragma omp parallel
  {
    double diffMax_ = -1e20;
    Vec3Int iposMax_(-1);
    int levelMax_=-1;

    #pragma omp  for schedule( static, 1 ) 
    for(int bid=0; bid<bgi._nBlocks; bid++)
    {
      int level = getBlockLevel(bid, levelMg );  

      if(doOnlyLevel0 && level>0) continue;

      CellFlags *dataUintA=NULL,*dataUintB=NULL;
      float *dataA=NULL,*dataB=NULL;
      Vec3Float *dataVecA=NULL,
		*dataVecB=NULL;

      SparseGrid<CellFlags> *sgFlags = 
	getFlagsChannel(CH_CELL_FLAGS,level,levelMg);

      if(isVecChannel(chanA))
      {
	dataVecA = getVecChannel(chanA,level,levelMg)->getDataPtrGen_(bid,0,0);
	dataVecB = getVecChannel(chanB,level,levelMg)->getDataPtrGen_(bid,0,0);
      }
      else if(isFlagsChannel(chanA))
      {
	dataUintA = getFlagsChannel(chanA,level,levelMg)->getDataPtrGen_(bid,0,0);
	dataUintB = getFlagsChannel(chanB,level,levelMg)->getDataPtrGen_(bid,0,0);
      }
      else
      {
	dataA = getFloatChannel(chanA,level,levelMg)->getDataPtrGen_(bid,0,0);
	dataB = getFloatChannel(chanB,level,levelMg)->getDataPtrGen_(bid,0,0);
      }

      SparseGrid<float> 
	*sgVisDiff = chVisDiff ? getFloatChannel(chVisDiff,level,levelMg) : NULL;
      float *dataVisDiff = sgVisDiff ? sgVisDiff->getDataPtrGen_(bid,1) : NULL;

      // Loop over each voxel in block
      BlockAccessor bac(&bgi,sgFlags,bid,level);
      for(int vz=0;vz<bac._bsz;vz++)
      for(int vy=0;vy<bac._bsy;vy++)
      for(int vx=0;vx<bac._bsx;vx++)
      {
	int vid = bac.getCellIndex( vx, vy, vz );
	double diff,valA=0,valB=0;

	if(dataVecA)
	{
	  Vec3Float d = dataVecA[vid]-dataVecB[vid];
	  diff = fabs(d.x) + fabs(d.y) + fabs(d.z);
	  valA = norm(dataVecA[vid]);
	}
	else
	{
	  if(dataA)
	  {
	    valA = dataA[vid];
	    valB = dataB[vid];
	  }
	  else
	  {
	    valA = dataUintA[vid];
	    valB = dataUintB[vid];
	  }

	  diff = fabs(valA-valB);
	}
 
	if(diff > diffMax_)
	{
	  diffMax_ = diff;
	  int x,y,z; bac.getGridCoords( vx,vy,vz, x,y,z );
	  iposMax_ = Vec3Int(x,y,z);
	  levelMax_ = level;
	}

	if(dataVisDiff) dataVisDiff[vid] = diff;

	if(diff>errThr)
//	if(diff>1e1)
//	if(diff>1e-4)
//	if(diff>1e-7)
	{
	  double relErr = fabs(valA)>MT_NUM_EPS ? diff/fabs(valA) : 0;
	   /* int x,y,z;
	    bac.getGridCoords( vx, vy, vz,
			       x, y, z );*/
	  char chbuf[100];

	  if(dataVecA)
	    sprintf(chbuf,"%g,%g,%g <-> %g,%g,%g",
		dataVecA[vid].x, dataVecA[vid].y, dataVecA[vid].z,
		dataVecB[vid].x, dataVecB[vid].y, dataVecB[vid].z);
	  else
	    sprintf(chbuf,"%g <-> %g",valA,valB);

	  #if 0
	  if((relErr>0.05 && fabs(valA)>1e-4) || 
	      ((fabs(valA)<MT_NUM_EPS) != (fabs(valB)<MT_NUM_EPS)))
	  #endif
	  {
	    int x,y,z;
	    bac.getGridCoords( vx, vy, vz,
			       x, y, z );	  
	    TRCERR(("L%d/%d %d:%d vpos=%d,%d,%d pos=%d,%d,%d => %s e=%g relErr=%g,"
		  "dbgCnt=%d   errThr=%g\n",
		    levelMg,level,bid,vid, vx,vy,vz, x,y,z, chbuf ,diff, relErr,
		    UtGlobDbgCheckCount,errThr));
	  }
	}
      }
    }
    #pragma omp critical
    {
      if(diffMax_ > diffMax )
      {
	diffMax = diffMax_;
	iposMax = iposMax_;
	levelMax = levelMax_;
      }
    }
  }

  if(doVisu && chVisDiff )
  {
    Vec3Float pos0 = Vec3Float(iposMax.x,iposMax.y,iposMax.z);
    pos0 = (pos0+0.5f)*(1<<levelMax);
    Vec3Int ipos0(pos0.x, pos0.y, pos0.z);
    _sg0->clipGridCoords(ipos0.x,ipos0.y,ipos0.z);
    double slicePos[100] = 
      { (double)(-ipos0.x), (double)(-ipos0.y), (double)(-ipos0.z) };

    static int pnl=-1;
    visualizeSlices( MiGetAuxPanel2(&pnl,"compareChannel: Abs. difference"),
		chVisDiff, IP_NEAREST, slicePos
		);
    resetChannel(chVisDiff);
  }

  return diffMax;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::sqrtChannel( int iChanDst )
{
  TRC3(("%s\n",UT_FUNCNAME));

  ThrRunParallel<int>( nBlocks(),nullptr, 
    [&](int &tls, int tid, int bid) 
    {
      BlockInfo *bi = getBlockInfo(bid);
      SparseGrid<float> *sg = getFloatChannel(iChanDst,bi->level);
      float *data = sg->getBlockDataPtr(bid,1,1);
      for(int i=0;i<sg->nVoxelsInBlock();i++) data[i] = sqrtf(MAX(data[i],0));
    } );
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template< typename Data_T >
void MultiresSparseGrid::setChannelConstBlock( SBG::ConstVal constVal, 
    					  int iChanDstY, int levelMg )
{
  TRC3(("%s %d %d %d\n",
	UT_FUNCNAME,(int)constVal,(int)iChanDstY,levelMg));

  UT_ASSERT0( levelMg < _nLevels );
  UT_ASSERT0( constVal == SBG::CONSTVAL_FULL || 
              constVal == SBG::CONSTVAL_EMPTY );

  ThrRunParallel<int>( _nBlocks, nullptr, 
    [&](int &tls, int tid, int bid) 
    {
      BlockInfo *bi=getBlockInfo(bid,levelMg);
      SparseGrid<Data_T> *sg = getChannel<Data_T>(this,iChanDstY,bi->level,levelMg);
      if(constVal==SBG::CONSTVAL_FULL)
	sg->setFullBlock(bid);
      else
	sg->setEmptyBlock(bid);
    } );

  invalidateForInterpolation(iChanDstY);
}
template 
void MultiresSparseGrid::setChannelConstBlock<float>( SBG::ConstVal constVal, 
    					  int iChanDstY, int levelMg );
template 
void MultiresSparseGrid::setChannelConstBlock<uint16_t>( SBG::ConstVal constVal, 
    					  int iChanDstY, int levelMg );
template 
void MultiresSparseGrid::setChannelConstBlock<uint8_t>( SBG::ConstVal constVal, 
    					  int iChanDstY, int levelMg );

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::setChannel( double a, 
				     int iChanDstY, 
				     int levelMax,
				     int levelMg,
				     CellFlags cellMask,
				     std::vector<int> *blockList,
				     unsigned skipBlockFlags
				     )
{
  // Y = a
  TRC3(("%s %g %d %d\n",UT_FUNCNAME,a,(int)iChanDstY,levelMg));
  UT_ASSERT0(levelMax==-1);

  prepareDataAccess( iChanDstY );

  int isVecChan = isVecChannel(iChanDstY);//==CH_VEC3_1||iChanDstY==CH_FORCE;

  size_t numActCells = nActCells(levelMg);

  if(levelMg>=_nLevels)  // special case: dense (coarse) grid
  {
    UT_ASSERT0(nBlocks(levelMg)==1);

    CellFlags *dataFlags = cellMask ? 
      getFlagsChannel(CH_CELL_FLAGS,levelMg,levelMg)->getDataPtrGen_() : 
      NULL;

    if(isVecChan)
    {
      SparseGrid<Vec3Float> *sgDst = getVecChannel(iChanDstY,levelMg,levelMg);
      Vec3Float *dataDst = sgDst->getDataPtrGen_(0,1,0,1);
      LongInt n = sgDst->nTotActVoxels();
      if(cellMask)
      {
        #pragma omp parallel for  if( n > N_CELLS_PARALLEL)
        for(LongInt i=0;i<n;i++) if(dataFlags[i]&cellMask) dataDst[i] = a;
      }
      else
      {
        #pragma omp parallel for  if( n > N_CELLS_PARALLEL)
        for(LongInt i=0;i<n;i++) dataDst[i] = a;
      }
    }
    else if (isPSFloatChannel(iChanDstY))
    {
      UT_ASSERT2(!cellMask);
      SparseGrid<PSFloat> *sgDst = getPSFloatChannel(iChanDstY,levelMg,levelMg);
      PSFloat *dataDst = sgDst->getDataPtrGen_(0,1,0,1);
      LongInt n = sgDst->nTotActVoxels();
      #pragma omp parallel for  if( n > N_CELLS_PARALLEL)
      for(LongInt i=0;i<n;i++) dataDst[i] = a;
    }
    else
    {
      SparseGrid<float> *sgDst = getFloatChannel(iChanDstY,levelMg,levelMg);
      float *dataDst = sgDst->getDataPtrGen_(0,1,0,1);
      LongInt n = sgDst->nTotActVoxels();
      if(cellMask)
      {
        #pragma omp parallel for  if( n > N_CELLS_PARALLEL)
        for(LongInt i=0;i<n;i++) if(dataFlags[i]&cellMask) dataDst[i] = a;
      }
      else
      {
        #pragma omp parallel for  if( n > N_CELLS_PARALLEL)
        for(LongInt i=0;i<n;i++) dataDst[i] = a;
      }
    }
    return;
  } // dense grid

  int nActBlocks = blockList ? blockList->size() : _nBlocks;

  tbb::parallel_for( tbb::blocked_range<size_t>(0,nActBlocks), 
  [=](const tbb::blocked_range<size_t>& tbb_range) 
  {
    for(size_t ibid=tbb_range.begin(); ibid != tbb_range.end(); ++ibid)
    {

    int bid = blockList ? (*blockList)[ibid] : ibid;
    UT_ASSERT2(bid>=0&&bid<_nBlocks);

    int level = MAX(_blockmap[0][bid].level,levelMg);

    BlockInfo *bi = getBlockInfo(bid,levelMg);
    if( bi->flags & skipBlockFlags ) continue;

    CellFlags *dataFlags = cellMask ? 
      getFlagsChannel(CH_CELL_FLAGS,level,levelMg)->getBlockDataPtr(bid) : 
      NULL;

    if(isVecChan)
    {
      SparseGrid<Vec3Float> *sgDstY = getVecChannel(iChanDstY, level, levelMg ); 
      Vec3Float *dataDstY = sgDstY->getBlockDataPtr( bid, TRUE, FALSE );

      Vec3Float va(a);
      if(cellMask)
      {
	for(int i=0; i<sgDstY->nVoxelsInBlock(); i++) 
	  if(dataFlags[i]&cellMask) dataDstY[i] = va;
      }
      else
      {
	for(int i=0; i<sgDstY->nVoxelsInBlock(); i++) 
	  dataDstY[i] = va;
      }
    }
    else if( isPSFloatChannel( iChanDstY ))
    {
      UT_ASSERT2(!cellMask);
      SparseGrid<PSFloat> *sgDstY = getPSFloatChannel(iChanDstY, level, levelMg ); 
      PSFloat *dataDstY = sgDstY->getBlockDataPtr( bid, TRUE, FALSE );
      Vec4psf A4(a);	
      for(int i=0; i<sgDstY->nVoxelsInBlock(); i+=4) A4.store_a(dataDstY+i);
    }
    else
    {
      SparseGrid<float> *sgDstY = getFloatChannel(iChanDstY, level, levelMg ); 
      float *dataDstY = sgDstY->getBlockDataPtr( bid, TRUE, FALSE );

      float af = a;
      if(cellMask)
      {
	for(int i=0; i<sgDstY->nVoxelsInBlock(); i++) 
	  if(dataFlags[i]&cellMask) dataDstY[i] = af;
      }
      else
      {
	Vec8f A(af);
	UT_ASSERT0( (sgDstY->nVoxelsInBlock() & 7) == 0 );
	
	size_t n=sgDstY->nVoxelsInBlock();

	if( numActCells * sizeof(*dataDstY) > CPU_BYPASS_CACHE_SIZE )
	{
	  for(size_t i=0; i<n; i+=8) vstream(dataDstY+i,A);
	}
	else
	{
	  for(size_t i=0; i<n; i+=8) A.store_a(dataDstY+i);
	}
      }
    }
  }

  _mm_mfence();  // non-temp load/store

  } ); // parallel

  invalidateForInterpolation(iChanDstY);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::scaleChannel(  int iChanX, double a, double b,
    				        int iChanDstY,
					std::vector<int> *blockList
					)
{
  int levelMg = 0;
  // Y = a*X+b
  TRC3(("%s %d %d %g %g\n",UT_FUNCNAME,iChanX,iChanDstY,a,b));

  ThrRunParallel<int>( blockList ? blockList->size() : _nBlocks, nullptr,   
    [&](int &tls, int tid, int ibid)
    {
      int bid = blockList ? (*blockList)[ibid] : ibid;
      int level = getBlockLevel(bid,levelMg); 
      SparseGrid<CellFlags> *sg = getFlagsChannel(CH_CELL_FLAGS,level); 
      float *dataX = getFloatBlockDataPtr(iChanX,bid,level,levelMg,0,0),
	    *dataDstY = getFloatBlockDataPtr(iChanDstY,bid,level,levelMg,1,1);	   
      for(int i=0; i<sg->nVoxelsInBlock(); i++) dataDstY[i] = a*dataX[i] +b;
    } );

  invalidateForInterpolation(iChanDstY);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::multChannel(  int iChanX1, int iChanX2,
    				       int iChanDstY )
{
  // Y = X1 * X2;
  TRC3(("%s %d %d %d\n",UT_FUNCNAME,iChanX1,iChanX2,iChanDstY));

  SparseGrid<Vec3Float> *sg0 = _sparseGrids[0].vec3_1[0];
  LONG iTask=0;
  #pragma omp parallel 
  {
    for(;;)
    {
      uint32_t bid = InterlockedIncrement((LONG volatile *)&iTask) - 1;
      if(bid>=sg0->nBlocks()) break;

      int level = _blockmap[0][bid].level;
      CellFlags *cellFlags = 
	_sparseGrids[level].cellFlags[0]->getBlockDataPtr(bid);
      UT_ASSERT0(cellFlags);

      {
	SparseGrid<float> *sgX1 = getFloatChannel(iChanX1, level),
			  *sgX2 = getFloatChannel(iChanX2, level),
			  *sgDstY = getFloatChannel(iChanDstY, level ); 
	float *dataX1 = sgX1->getBlockDataPtr(bid),
	      *dataX2 = sgX2->getBlockDataPtr(bid),
	      *dataDstY = sgDstY->getBlockDataPtr(bid,TRUE,TRUE);

	for(LongInt i=0; i<sgDstY->nVoxelsInBlock(); i++) 
	  dataDstY[i] = (double)dataX1[i] * (double)dataX2[i];
      }
    }
  }
  invalidateForInterpolation(iChanDstY);
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::scaleByCellVolume( int chan )
{
#pragma omp parallel for schedule( static, 10 )
  for(int bid=0; bid < _nBlocks; bid++)
  {
    int level = _blockmap[0][bid].level;
    SparseGrid<float> *sg = getFloatChannel(chan,level);
    float *data = sg->getBlockDataPtr(bid,0,0);
    double h = (1<<level)*_dx0,
	   scale = 1./(h*h*h);
    for(int i=0;i<sg->nVoxelsInBlock();i++)
    {
	data[i] *= scale;
    }
  }
  invalidateForInterpolation(chan);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template< typename Float_T >
double MultiresSparseGrid::maxAbsChannel_(  int iChan, 
    				           int levelMax, 
    					   int levelMg,
					   int doScaleByCellVolume,
					   Vec3Int *outPosMax,
					   int *outLevelMax,
					   std::vector<int> *blockList 
					   )
{
  // 
  TRC3(("%s %d %d\n",UT_FUNCNAME,(int)iChan,levelMg));
  Vec3Int iposTotalMax(-1);
  int levelTotalMax=-1;
  double max=-1e20;

  if(levelMg>=_nLevels)  // special case: dense (coarse) grid
  {
    UT_ASSERT0(nBlocks(levelMg)==1);
    int level = levelMg;
    SparseGrid<Float_T> *sg = getChannel<Float_T>(this, iChan,levelMg,levelMg);
    Float_T *data = sg->getDataPtrGen_(0,0,0,1);
    LongInt n = sg->nTotActVoxels();
    USE(n);
    double h = (1<<levelMg)*_dx0,
	   scale = doScaleByCellVolume ? 1./(h*h*h) : 1.;
    #pragma omp parallel if( n>N_CELLS_PARALLEL )
    {
      double max_=-1e20;
      Vec3Int iposTotalMax_(-1);
      int levelTotalMax_=-1;

      #pragma omp for 
      for(int z=0;z<sg->sz();z++)
      for(int y=0;y<sg->sy();y++)
      for(int x=0;x<sg->sx();x++)
      {
	LongInt i = sg->getGridIndex(x,y,z);
	double f = fabs( scale*data[i] );
	if( f > max_ )
	{
	  max_ = f;
	  iposTotalMax_ = Vec3Int(x,y,z);
	  levelTotalMax_ = level;
	}
      }
      #pragma omp critical
      {
	if(max_>max)
	{
	  max = max_;
	  iposTotalMax = iposTotalMax_;
	  levelTotalMax = levelTotalMax_;
	}
      }
    }
  }
  else
  {
    int nActBlocks = blockList ? blockList->size() : _nBlocks;

    LONG iTask=0;
    #pragma omp parallel 
    {
      double max_=-1e20;
      Vec3Int iposTotalMax_(-1);
      int levelTotalMax_=-1;

      for(;;)
      {
	uint32_t ibid = InterlockedIncrement((LONG volatile *)&iTask) - 1;
	if((int)ibid>=nActBlocks) break;

	int bid = blockList ? (*blockList)[ibid] : ibid;

	int level = MAX( _blockmap[0][bid].level,levelMg );

	if(levelMax<_nLevels&&levelMax>=0)
	{
	  UT_ASSERT0(levelMg==0);
	  if(level>levelMax) continue;
	  level = levelMax;
	}

	double h = (1<<level)*_dx0,
	       scale = doScaleByCellVolume ? 1./(h*h*h) : 1.;

	SparseGrid<Float_T> *sg = getChannel<Float_T>(this, iChan, level, levelMg);
	Float_T *data = sg->getBlockDataPtr(bid);
	/*CellFlags *cellFlags = 
	  getFlagsChannel(CH_CELL_FLAGS,level,levelMg)->getBlockDataPtr(bid);*/

	SBG_FOR_EACH_VOXEL_IN_BLOCK( sg, bid, x, y, z, i )
	{
	  Float_T f = std::abs((Float_T)(scale*data[i]));
	  if( f > max_ )
	  {
	    max_ = f;
	    iposTotalMax_ = Vec3Int(x,y,z);
	    levelTotalMax_ = level;
	  }
	}
      }
      #pragma omp critical
      {
	if(max_>max)
	{
	  max = max_;
	  iposTotalMax = iposTotalMax_;
	  levelTotalMax = levelTotalMax_;
	}
      }
    }
  }
  if( outPosMax ) *outPosMax = iposTotalMax;
  if( outLevelMax ) *outLevelMax = levelTotalMax;
  return max;
}
template
double MultiresSparseGrid::maxAbsChannel_<float>(  int iChan, 
    				           int levelMax, 
    					   int levelMg,
					   int doScaleByCellVolume,
					   Vec3Int *outPosMax,
					   int *outLevelMax,
					   std::vector<int> *blockList 
					   );
template
double MultiresSparseGrid::maxAbsChannel_<double>(  int iChan, 
    				           int levelMax, 
    					   int levelMg,
					   int doScaleByCellVolume,
					   Vec3Int *outPosMax,
					   int *outLevelMax,
					   std::vector<int> *blockList 
					   );


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::absChannel(  int iChan )
{
  // 
  TRC3(("%s %d\n",UT_FUNCNAME,(int)iChan));

  SparseGrid<Vec3Float> *sg0 = _sparseGrids[0].vec3_1[0];

  LONG iTask=0;
  #pragma omp parallel 
  {
    for(;;)
    {
      uint32_t bid = InterlockedIncrement((LONG volatile *)&iTask) - 1;
      if(bid>=sg0->nBlocks()) break;
      int level = _blockmap[0][bid].level;
      SparseGrid<float> *sg = getFloatChannel(iChan, level);
      if(!sg->isValueBlock(bid)) continue;
      float *data = sg->getBlockDataPtr(bid);
      for(LongInt i=0; i<sg->nVoxelsInBlock(); i++) data[i] = fabsf(data[i]);
    }
  }
  invalidateForInterpolation(iChan);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::statChannel(  int iChan, 
    				       int doCountFluidCellsOnly,
				       double zeroEps,
				       MtStat *stat, // OUT
				       std::vector<int> *blockList
					)
{
  // 
  TRC3(("%s %d\n",UT_FUNCNAME,(int)iChan));
  int levelMg=0; 
  UT_ASSERT0(levelMg==0);

  int nActBlocks = blockList ? blockList->size() : _nBlocks;

  MT_STAT_INIT(stat);

  LONG iTask=0;
  #pragma omp parallel 
  {
    MtStat stat_;
    MT_STAT_INIT(&stat_);
    for(;;)
    {
      uint32_t ibid = InterlockedIncrement((LONG volatile *)&iTask) - 1;
      if((int)ibid>=nActBlocks) break;

      int bid = blockList ? (*blockList)[ibid] : ibid;

      int level = MAX( _blockmap[0][bid].level,levelMg );

      SparseGrid<CellFlags> *sgFlags = 
	getFlagsChannel(CH_CELL_FLAGS,level,levelMg);

      float *data=NULL;
      Vec3Float *dataVec=NULL;

      if(isVecChannel(iChan))
        dataVec = getVecChannel(iChan, level, levelMg)->getBlockDataPtr(bid);
      else
	data = getFloatChannel(iChan, level, levelMg)->getBlockDataPtr(bid);

      CellFlags *cellFlags = 
	sgFlags->getBlockDataPtr(bid);

      for(LongInt i=0; i<sgFlags->nVoxelsInBlock(); i++) 
      {
	if( !doCountFluidCellsOnly ||  CELL_IS_FLUID_( cellFlags[i] ))
	{
	  if(dataVec)
	  {
	    for(int k=0;k<3;k++)
	    {
	      if(!(zeroEps>0) || fabs(dataVec[i][k])>zeroEps)
	      {
		MT_STAT_UPDATE(&stat_,dataVec[i][k]);
	      }
	    }
	  }
	  else
	  {
	    if(!(zeroEps>0) || fabs(data[i])>zeroEps)
	    {
	      MT_STAT_UPDATE(&stat_,data[i]);
	    }
	  }
	}
      }
    }
    #pragma omp critical
    {
      MT_STAT_RESULT(&stat_); 
      MT_STAT_UPDATE_STAT(stat, &stat_);
    }
  }
  MT_STAT_RESULT(stat);
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::normalizeChannel(  int chan )
{
  TRC3(("%s %d\n",UT_FUNCNAME,(int)chan));
  MtStat stat;  
  statChannel( chan, FALSE, -1, &stat );
  double scale,offs;
  if(stat.span>MT_NUM_EPS)
  {
    scale = 1./stat.span;
    offs = -stat.min/stat.span;
  }
  else scale=offs=0;

  scaleChannel(chan,scale,offs,
	       chan);

  invalidateForInterpolation(chan);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::checkFloatChannel(  int iChan, int levelMg=0 )
{
  // 
  TRC3(("%s %d %d\n",UT_FUNCNAME,(int)iChan,levelMg));

  BlockGridInfo bgi( this, levelMg );

  // Loop over each block of levelMg
  #pragma omp parallel for schedule( static, 1 )
  for(int bid=0; bid<bgi._nBlocks; bid++)
  {
    int level = getBlockLevel(bid, levelMg );  
      
    SparseGrid<CellFlags> 
      *sgFlags = getFlagsChannel(CH_CELL_FLAGS,level,levelMg);
    SparseGrid<float> 
      *sg = getFloatChannel(iChan, level,levelMg);

    CellFlags *dataFlags = sgFlags->getDataPtrGen_(bid,0,0);
    float *data = sg->getDataPtrGen_(bid,0,0);

    BlockAccessor bac(&bgi,sgFlags,bid,level);

    // Loop over each voxel in block
    for(int vz=0;vz<bac._bsz;vz++)
    for(int vy=0;vy<bac._bsy;vy++)
    for(int vx=0;vx<bac._bsx;vx++)
    {
      int vid = bac.getCellIndex( vx, vy, vz );
      CellFlags flags = dataFlags[vid];
      float val = data[vid];
      UT_ASSERT0(std::isfinite(val));
      if(CELL_IS_FLUID_(flags))
      {
	// Interior cell
      }
      else
      {
	// Assume exactly zero for all non-interior cells !
	#if 1
	UT_ASSERT0(val==0.0f);
        #else
	if(!(val==0.0f))
	{
	  int x,y,z;
	  bac.getGridCoords( vx, vy, vz,
			     x, y, z );	  
	  TRCERR(("l=%d/%d %d,%d,%d flags=%d -> %g\n",level,levelMg,
		x,y,z,flags,val));
	}
	#endif
      }

    }
  }
}

// Fast versions

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template< typename Float_T >
long double MultiresSparseGrid::dotProdChannel_( int iChan1,
    				           int iChan2,
					   int levelMg,
					   std::vector<int> *blockList,
					   int *pErrCode )
{
  TRC3(("%s %d %d\n",UT_FUNCNAME,(int)iChan1,(int)iChan2));

  if(levelMg>=_nLevels)  // special case: dense (coarse) grid
  {
    UT_ASSERT0(nBlocks(levelMg)==1);
    SparseGrid<Float_T> *sg1 = getChannel<Float_T>(this,iChan1,levelMg,levelMg),
		        *sg2 = getChannel<Float_T>(this,iChan2,levelMg,levelMg);
    Float_T *data1 = sg1->getDataPtrGen_(),
	    *data2 = sg2->getDataPtrGen_();
    double sum=0.0f;

    #pragma omp parallel for
    for(int z=0;z<sg1->sz();z++)
    for(int y=0;y<sg1->sy();y++)
    for(int x=0;x<sg1->sx();x++)
    {	  
      LongInt i = sg1->getGridIndex(x,y,z);
      sum += data1[i]*data2[i];
      UT_ASSERT0(FALSE);  // check if we are are called here (unsafe reduction)
    }

    return sum;
  } // dense grid

  long double sum=0;

  int nActBlocks = blockList ? blockList->size() : _nBlocks;

  typedef struct
  {
    long double mySum;
  } 
  ThreadLocals;

  ThrThreadLocals<ThreadLocals> threadLocals( 
  [&]( ThreadLocals *tls, int tid )  // initialize
  {
    tls->mySum = 0;
  } );

  tbb::parallel_for( tbb::blocked_range<size_t>(0,nActBlocks), 
  [&](const tbb::blocked_range<size_t>& tbb_range) 
  {
    int tid = tbb::this_task_arena::current_thread_index();
    auto tls = threadLocals[tid];

    for(size_t ibid=tbb_range.begin(); ibid != tbb_range.end(); ++ibid)
    {
      int bid = blockList ? (*blockList)[ibid] : ibid;
      BlockInfo *bi=getBlockInfo(bid,levelMg);

      int level = bi->level;
      SparseGrid<Float_T> *sg1 = getChannel<Float_T>(this,iChan1, level, levelMg),
			  *sg2 = getChannel<Float_T>(this,iChan2, level, levelMg);
      Float_T *data1= sg1->getBlockDataPtr(bid),
	      *data2= sg2->getBlockDataPtr(bid);
      /*CellFlags *cellFlags = 
	_sparseGrids[level].cellFlags->getBlockDataPtr(bid);*/

      double mySum2=0;

      Vec4d S(0.0);

      #if 0
      UT_ASSERT0((sg1->nVoxelsInBlock() & 15) ==0);
      for(LongInt i=0; i<sg1->nVoxelsInBlock(); i+=8) 
      {
        Vec8f X1,X2;
	X1.load_a(data1+i);
	X2.load_a(data2+i);
	Vec4d A1,A2;

	A1 = extend_low(X1);
	A2 = extend_low(X2);
	S += A1*A2;

	A1 = extend_high(X1);
	A2 = extend_high(X2);
	S += A1*A2;
      }
      #else

      for(LongInt i=0; i<sg1->nVoxelsInBlock(); i+=4) 
      {
	Vec4d X1 = loadf_a(data1+i),
	      X2 = loadf_a(data2+i);
	S += X1*X2;
      }

      #endif

      mySum2 = horizontal_add(S);

      tls->mySum += mySum2;
    }
  } );
  
  // Reduce thread locals
  for(int tid=0; tid < threadLocals.size(); tid++)
  {
    auto tls = threadLocals[tid];
    sum += tls->mySum;
  }    

  int errCode = MI_OK;
  if( !(std::isfinite(sum)) )
  {
    TRCERR(("%s: invalid sum %Lf\n",UT_FUNCNAME,sum));
    errCode = MI_ERROR;
    _doDumpSimulationState = 1;
  }
  if(pErrCode) *pErrCode = errCode;

  //UT_ASSERT0(std::isfinite(sum));
  return sum;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::AXPBYChannelCombined(  double a, double b,    						
    						int iChanQ,
						int iChanX,
						int iChanP,
						std::vector<int> *blockList
					)
{
  UT_ASSERT2( isPSFloatChannel(iChanX) && isPSFloatChannel(iChanQ));

  #define AXPBYChannelCombined_SKIP_NONFLUID
  //#define AXPBYChannelCombined_DO_TIMER

  #ifdef AXPBYChannelCombined_DO_TIMER
  UtTimer tm;
  TIMER_START(&tm);
  #endif
  // X = X+aP
  // P = Q+bP
  TRC3(("%s %g %g %d %d %d\n",UT_FUNCNAME,a,b,iChanQ,iChanX,iChanP));

  int levelMg=0;

  SparseGrid<Vec3Float> *sg0 = _sparseGrids[0].vec3_1[0];
  int nActBlocks = blockList ? blockList->size() : _nBlocks;

  tbb::parallel_for( tbb::blocked_range<size_t>(0,nActBlocks), 
  [&](const tbb::blocked_range<size_t>& tbb_range) 
  {
    for(size_t ibid=tbb_range.begin(); ibid != tbb_range.end(); ++ibid)
    {
      int bid = blockList ? (*blockList)[ibid] : ibid;
      UT_ASSERT2(bid>=0&&bid<_nBlocks);
      if(bid>=sg0->nBlocks()) break;
      int level = MAX(_blockmap[0][bid].level,levelMg);
    
      SparseGrid<PSFloat> *sgQ = getPSFloatChannel(iChanQ,level),
			  *sgX = getPSFloatChannel(iChanX,level),
			  *sgP = getPSFloatChannel(iChanP,level);

      PSFloat *dataQ = sgQ->getBlockDataPtr(bid,0,0),
	      *dataX = sgX->getBlockDataPtr(bid,1,0),
	      *dataP = sgP->getBlockDataPtr(bid,1,0);

     #ifdef AXPBYChannelCombined_SKIP_NONFLUID
     BlockInfo *bi = getBlockInfo(bid,0);
     CellFlags *dataFlags = ( bi->flags & BLK_ONLY_FLUID ) ? NULL :
       			     getFlagsChannel(CH_CELL_FLAGS,level)->getBlockDataPtr(bid);	    
     #endif

      UT_ASSERT0(dataX&&dataQ&&dataP);
      UT_ASSERT0((sgX->nVoxelsInBlock() & 7) ==0);

      const Vec4ui MS(-1);
      Vec4d A(a),B(b);
      for(int i=0; i<sgX->nVoxelsInBlock(); i+=8) 
      {

	#ifdef AXPBYChannelCombined_SKIP_NONFLUID
	if(dataFlags)
	{
	  Vec8ui flags;
	  vload_from_uint16_a( flags, dataFlags+i );
	  bool doSkip = ! horizontal_or((flags & (CELL_SOLID|CELL_VOID)) == 0 );
	  if(doSkip)
	  {
	    #ifdef UT_ASSERT_LEVEL_2
	    #ifndef SOLVE_PRESSURE_DOUBLE
	    UT_ASSERT2( vall( (Vec8f)(Vec8f().load(dataX+i)==0.0f) ) );
	    UT_ASSERT2( vall( (Vec8f)(Vec8f().load(dataP+i)==0.0f) ) );
	    #endif
	    #endif
	    continue;
	  }
	}
	#endif

	Vec4d P1=loadf_a(dataP+i) ,
	      P2=loadf_a(dataP+i+4) ,

	      Q1=loadf_a(dataQ+i),
	      Q2=loadf_a(dataQ+i+4),

	      X1=loadf_a(dataX+i),
	      X2=loadf_a(dataX+i+4);

	      X1 = X1 + A*P1;
	      X2 = X2 + A*P2;

	      P1 = Q1+B*P1;
	      P2 = Q2+B*P2;

	 #ifdef SOLVE_PRESSURE_DOUBLE
	 X1.store_a(dataX+i);
	 X2.store_a(dataX+i+4);
	 P1.store_a(dataP+i);
	 P2.store_a(dataP+i+4);
         #else
	 vstreamf_a(dataX+i,X1);
	 vstreamf_a(dataX+i+4,X2);
	 vstreamf_a(dataP+i,P1);
	 vstreamf_a(dataP+i+4,P2);
	 #endif
      }
    }
    _mm_mfence();  // non-temp load/store    
  } ); // parallel
  
  invalidateForInterpolation(iChanX);
  invalidateForInterpolation(iChanP);

  #ifdef AXPBYChannelCombined_DO_TIMER
    TIMER_STOP(&tm);
    TRC(("%s CPU=%.2f sec, %.0f cells/sec/iter\n",UT_FUNCNAME,
	  (double)TIMER_DIFF_MS(&tm)/1000.,
	(double)_nActCells[0]/(double)(TIMER_DIFF_MS(&tm)/1000.0)));

  #endif
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::AXPBYChannel(  double a, int iChanX,	
    					double b,  int iChanDstY,
					int levelMax,  
					int levelMg,
					std::vector<int> *blockList
					)
{
  // Y = a*X + b*Y
  TRC3(("%s %g %d %g %d\n",UT_FUNCNAME,a,(int)iChanX,b,(int)iChanDstY));

  UT_ASSERT2( isPSFloatChannel(iChanX) && isPSFloatChannel(iChanDstY));
      
  if(levelMg>=_nLevels)  // special case: dense (coarse) grid
  {
    UT_ASSERT0(nBlocks(levelMg)==1);
    SparseGrid<PSFloat> *sgDst = getPSFloatChannel(iChanDstY,levelMg,levelMg),
		        *sgSrc = getPSFloatChannel(iChanX,levelMg,levelMg);
    PSFloat *dataDst = sgDst->getDataPtrGen_(0,1,0,1),
	    *dataSrc = sgSrc->getDataPtrGen_(0,0,0);
    LongInt n = sgDst->nTotActVoxels();
    #pragma omp parallel for if( n > N_CELLS_PARALLEL)
    for(LongInt i=0;i<n;i++) dataDst[i] = a*dataSrc[i] + b*dataDst[i];
    return;
  } // dense grid

  int nActBlocks = blockList ? blockList->size() : _nBlocks;

  tbb::parallel_for( tbb::blocked_range<size_t>(0,nActBlocks), 
  [&](const tbb::blocked_range<size_t>& tbb_range) 
  {
    for(size_t ibid=tbb_range.begin(); ibid != tbb_range.end(); ++ibid)
    {

    int bid = blockList ? (*blockList)[ibid] : ibid;
    UT_ASSERT2(bid>=0&&bid<_nBlocks);

    int level = MAX(_blockmap[0][bid].level,levelMg);
    
    if(levelMax<_nLevels&&levelMax>=0)
    {
      UT_ASSERT0(levelMg==0);
      if(level>levelMax) continue;
      level = levelMax;
    }

    CellFlags *cellFlags = 
      getFlagsChannel( CH_CELL_FLAGS, level, levelMg )->getBlockDataPtr(bid);
    UT_ASSERT0(cellFlags);

    if(iChanX==CH_VEC3_1||iChanDstY==CH_VEC3_1)
    {
      UT_ASSERT0(levelMg==0);
      SparseGrid<Vec3Float> *sgX = getVecChannel(iChanX, level),
			    *sgDstY = getVecChannel(iChanDstY, level ); 
      Vec3Float *dataDstY = sgDstY->getBlockDataPtr(bid,TRUE,TRUE),
		 *dataX = sgX->getBlockDataPtr(bid);
      for(LongInt i=0; i<sgDstY->nVoxelsInBlock(); i++) 
      {
	dataDstY[i][0] = (double)(dataX[i][0])*a + (double)(dataDstY[i])[0]*b;
	dataDstY[i][1] = (double)(dataX[i][1])*a + (double)(dataDstY[i])[1]*b;
	dataDstY[i][2] = (double)(dataX[i][2])*a + (double)(dataDstY[i])[2]*b;
      }
    }
    else
    {
      SparseGrid<PSFloat> *sgX = getPSFloatChannel(iChanX, level, levelMg),
	                  *sgDstY = getPSFloatChannel(iChanDstY, level, levelMg ); 
      PSFloat *dataX = sgX->getBlockDataPtr(bid),
	      *dataDstY = sgDstY->getBlockDataPtr(bid,TRUE,TRUE);

      UT_ASSERT0(dataX);
      //UT_ASSERT0((sgDstY->nVoxelsInBlock() & 15) ==0);

      const Vec4ui MS(-1);
      Vec4d A(a),B(b);

      for(int i=0; i<sgDstY->nVoxelsInBlock(); i+=8) 
      {
	Vec4d X1=loadf_a(dataX+i),
	      Y1=loadf_a(dataDstY+i),
	      X2=loadf_a(dataX+i+4),
	      Y2=loadf_a(dataDstY+i+4);
	Y1 = A*X1+B*Y1,
	Y2 = A*X2+B*Y2;

	#ifdef SOLVE_PRESSURE_DOUBLE
	Y1.store_a( dataDstY+i );
	Y2.store_a( dataDstY+i+4 );
        #else
	vstreamf_a(dataDstY+i,Y1);
	vstreamf_a(dataDstY+i+4,Y2);
	#endif
      }
    }
  }
   _mm_mfence();  // non-temp load/store    
  });

  invalidateForInterpolation(iChanDstY);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MultiresSparseGrid::visGetCellColor( 
    				double g, double g0,
    			       uint16_t cellFlags, uint16_t blockFlags,
			       int level, int level0,
			       unsigned options
			       )
{
  using namespace MSBG;

  /*if( options & VIS_OPT_NULLSPACE_BLOCKS )
  {
    UT_ASSERT0(blockId>=0&&blockId<_nBlocks);
    NullSpaceInfo *nsi=_blockNullSpace[bid];
    UT_ASSERT0(nsi->id<64)
    return nsi ? BMP_MKRGB( (((1+nsi->id)&3)>>0) *100, 
			    (((1+nsi->id)&12)>>2) *100, 
			    (((1+nsi->id)&48)>>4) *100 ) : 
      		 BMP_MKRGB(0,0,0);
  }*/
  if( options & VIS_OPT_COLOR_RGB )
  {
    uint32_t color = UT::UtFloatToIntBits( g0 );
    return color;
  }

  if( !(g0<PCELL_INF_QUANTITY)) return VIS_COL_INVALID_QUANTITY;

  if( options & VIS_OPT_COLOR_PAL )
  {
    int c = round(g0);

    if( c==0 )
    {
      if( cellFlags & CELL_SOLID ) return BMP_MKRGB(200,200,200);
      if( cellFlags & CELL_AIR ) return BMP_MKRGB(10,10,0);
      if( cellFlags & CELL_VOID ) return BMP_MKRGB(0,0,0);
      return BMP_MKRGB(128,128,128);
    }
    else
    {
      switch(c)
      {
	case 1 : return BMP_MKRGB(0xFF,0x20,0x20);
	case 2 : return BMP_MKRGB(0x70,0xFF,0x70);
	case 3 : return BMP_MKRGB(0x60,0x70,0xFF);
	case 4 : return BMP_MKRGB(0xFF,0x0,0xFF);
	case 5 : return BMP_MKRGB(0xFF,0xFF,0x00);
	case 6 : return BMP_MKRGB(0x40,0xB0,0x40);
	default : return BMP_MKRGB(0xFF,0,0);
      }
    }
  }

  UT_ASSERT0(!(g<-.0001||g>1.0001));
  MT_CLAMP(g,0,1);

  int doVisTmpMark = options & VIS_OPT_FLAGTMPMARK;

  float hue = 0, sat=-0.5;

#if 1
  if( cellFlags & CELL_SOLID )
  {
    hue = 0.2;
  }
  else if( cellFlags & CELL_AIR )
  {
    hue = 0.3;
  }
  else if( cellFlags & CELL_VOID )
  {
    hue = 0.4;
  }
  /*else if( cellFlags & CELL_DOM_BORDER )
  {
    hue = 0.1;
    sat = -0.5;
  }*/
  else if( doVisTmpMark && (cellFlags & CELL_TMP_MARK_2))
  {
    hue = -0.4;
  }
  else
#endif
  {
    hue = 0;
    sat = -1;
  }
  
#if 0
  if(blockFlags & BLK_NO_FLUID)
  {
    hue = -.3;
    sat = -0.5;
    g = 0.5;
  }
  if(blockFlags & BLK_ONLY_FLUID)
  {
    hue = .3;
    sat = -0.5;
    g=0.5;
  }
#endif

  if( level>=0 )
  {
    hue = fmod( hue+0.8*(level/(double)getNumLevels()), 1 );
    sat = level==0 && CELL_IS_FLUID_(cellFlags) ? -1:-0.7;
  }

#if 0
  if((blockFlags & BLK_TMP_MARK) || (blockFlags & BLK_INACTIVE) )
  {
    hue = -.2;
    sat = -0.5;
    g = MAX(g,0.1);
  }
#endif

  float vHCL[3] = { 
#if 0
    			(float) (level==0 ? g : (0.1+0.8999*g)),   // brightness [0,1]
#else
    			(float) (level==0 ? g : (0.1+0.8*g)),   // brightness [0,1]
#endif
		    hue, 	// Hue [-1,1]
		    sat },	// Saturation [-1,1]
	vRGB[3];
	CoConvHcl2Rgb( vRGB, vHCL );

   int	color;
   
   if(options & VIS_OPT_SIGN_NZ) 
   {
     color = g0>0 ? 
       BMP_MKRGB( 0, vRGB[1]*255, options & VIS_OPT_RESLEV ? level0*80 : 0 ) :
       BMP_MKRGB( vRGB[0]*255, 0, options & VIS_OPT_RESLEV ? level0*80 : 0 );
   }
   else color = BMP_MKRGB( vRGB[0]*255, vRGB[1]*255, vRGB[2]*255 );
  
   return color;
}

#if 0
static int visGetCellColor2( double g, 
    			       uint16_t cellFlags, uint16_t blockFlags )
{
  using namespace MSBG;
  UT_ASSERT0(!(g<-.0001||g>1.0001));
  MT_CLAMP(g,0,1);

  float hue = 0, sat=-0.5;

  if( cellFlags & CELL_SOLID )
  {
    hue = 0.2;
  }
  else if( cellFlags & CELL_DOM_BORDER )
  {
    hue = 0.6;
  }
  else
  {
    hue = 0;
    sat = -1;
  }

  float vHCL[3] = { (float)g,   // brightness [0,1]
		    hue, 	// Hue [-1,1]
		    sat },	// Saturation [-1,1]
	vRGB[3];
	CoConvHcl2Rgb( vRGB, vHCL );
   int	color = BMP_MKRGB( vRGB[0]*255, vRGB[1]*255, vRGB[2]*255 );
   return color;
}
#endif

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MultiresSparseGrid::getSlices2D( 
	int chanId,
	int ipolType,
	double *slicePos0,
	double roiSize,    // if >0 then use cubic region-of-interest
	double roiZoom,
	unsigned opt,
	int levelMg,
	BmpBitmap **pOutXZ,   // result in dataFloat[0..2] channel
	BmpBitmap **pOutXY,
	BmpBitmap **pOutZY
    ) 
{
  int rcThis=0,doUseOldRoi = opt & OPT_USE_OLD_BEHAVIOR;

  UT_ASSERT0(levelMg==0);
  int isVecChan = isVecChannel(chanId);
  int nChan = isVecChannel(chanId) ? 3 : 1;
  double *slicePos = slicePos0 ? slicePos0 : _visSlicePos;

  SparseGrid<Vec3Float> *sg0 = _sparseGrids[0].vec3_1[0];

  G3D_Grid G0 = 
    G3D_CreateGrid( sg0->sx(), sg0->sy(), sg0->sz(),0);

  int sxRoi = roiSize>0 ? roiSize*G0->smin : 0;

  TRC(("%s slice=%g,%g,%g, roi_sz=%g\n",UT_FUNCNAME,
	slicePos[0],slicePos[1],slicePos[2],roiSize));

  Vec3Int iSlicePos( (int)(sg0->sx()*(double)slicePos[0]),
      		     (int)(sg0->sy()*(double)slicePos[1]),
      		     (int)(sg0->sz()*(double)slicePos[2]) );

  BmpBitmap *BXZ=NULL,*BXY=NULL,*BZY=NULL;
  BXZ = BmpNewBitmap( sxRoi ? sxRoi:sg0->sx(), sxRoi ? sxRoi:sg0->sz(), 0);
  BXY = BmpNewBitmap( sxRoi ? sxRoi:sg0->sx(), sxRoi ? sxRoi:sg0->sy(), 0);
  BZY = BmpNewBitmap( sxRoi ? sxRoi:sg0->sz(), sxRoi ? sxRoi:sg0->sy(), 0);

  BmpRectangle roiXZ,roiXY,roiZY;
  if(doUseOldRoi)
  {
    BMP_SET_RECT( &roiXZ, (G0->sx-sxRoi)/2, (G0->sz-sxRoi)/2, sxRoi, sxRoi );
    BMP_SET_RECT( &roiXY, (G0->sx-sxRoi)/2, (G0->sy-sxRoi)/2, sxRoi, sxRoi );
    BMP_SET_RECT( &roiZY, (G0->sz-sxRoi)/2, (G0->sy-sxRoi)/2, sxRoi, sxRoi );
  }
  else
  {
    BMP_SET_RECT( &roiXZ, iSlicePos.x-sxRoi/2, iSlicePos.z-sxRoi/2, sxRoi, sxRoi );
    BMP_SET_RECT( &roiXY, iSlicePos.x-sxRoi/2, iSlicePos.y-sxRoi/2, sxRoi, sxRoi );
    BMP_SET_RECT( &roiZY, iSlicePos.z-sxRoi/2, iSlicePos.y-sxRoi/2, sxRoi, sxRoi );
  }

  for(int k=0;k<nChan;k++)
  {
    BmpGetFloatChannel(BXZ,k,BMP_CLEAR);
    BmpGetFloatChannel(BXY,k,BMP_CLEAR);
    BmpGetFloatChannel(BZY,k,BMP_CLEAR);
  }

  #pragma omp parallel for schedule( dynamic )
  for(int iz=0;iz<sg0->sz();iz++) 
  for(int iy=0;iy<sg0->sy();iy++)
  for(int ix=0;ix<sg0->sx();ix++)
  {
    int yIsSlice = (iy==iSlicePos.y),
	zIsSlice = (iz==iSlicePos.z),
	xIsSlice = (ix==iSlicePos.x);
    if (xIsSlice || yIsSlice || zIsSlice )
    {
      int bx,by,bz;
      sg0->getBlockCoords(ix,iy,iz, bx,by,bz);
      UT_ASSERT0(sg0->blockCoordsInRange(bx,by,bz));
      LongInt bid = sg0->getBlockIndex(bx,by,bz);
      UT_ASSERT2(bid>=0&&bid<sg0->nBlocks());
      BlockInfo *bi=getBlockInfo(bid,levelMg);
      int levelBlk = bi->level,
	  level = levelBlk;

      if(levelMg>0)
      {
	level = MAX(level,levelMg);
      }

      int ix1=ix>>level,
	  iy1=iy>>level,
	  iz1=iz>>level;

      SparseGrid<CellFlags> 
	*sgFlags = getFlagsChannel(CH_CELL_FLAGS,level,levelMg);

      CellFlags  cellFlags = sgFlags->getValueGen_(ix1,iy1,iz1);
      USE(cellFlags);

      Vec3Float data;
      float pos[3] = { (float)ix+0.5f, (float)iy+0.5f, (float)iz+0.5f };

      int ipolOpt = OPT_IPCHKBLK|
		    OPT_IPEXTRAPOL|
		    OPT_IPBC_NEUMAN;
      switch(ipolType)
      {
	case IP_NEAREST:
	  if(levelMg>0)
	  {
	    if(isVecChan)
	      data = getVecChannel(chanId,level,levelMg)->getValueGen_(ix1,iy1,iz1);
	    else
	      data[0] = getFloatChannel(chanId,level,levelMg)->getValueGen_(ix1,iy1,iz1);
	  }
	  else 
	  {
	    interpolate<IP_NEAREST>( chanId, pos, 0, data ); 
	  }
	  break;
	case IP_LINEAR: interpolate<IP_LINEAR>( chanId, pos, ipolOpt, data ); break;
	case IP_CUBIC_MONO_2: interpolate<IP_CUBIC_MONO_2>( chanId, pos, ipolOpt, data ); break;
	case IP_CUBIC: interpolate<IP_CUBIC>( chanId, pos, ipolOpt, data ); break;
	case IP_WENO4: interpolate( IP_WENO4, chanId, pos, ipolOpt, data ); break;
	default:
	  UT_ASSERT0(FALSE);
	  break;
      }
      if(yIsSlice)
      {
	BmpBitmap *B=BXZ;
	BmpRectangle *roi=&roiXZ;
	int ix2 = ix - roi->x0,
	    iz2 = iz - roi->y0;
	if(BMP_IN_RANGEB(ix2,iz2,B))
	{
	  LongInt i=BMP_XYOFF(B,ix2,iz2);
	  for(int k=0;k<nChan;k++)  B->dataFloat[k][i] = data[k];
	}
      }
      if(zIsSlice)
      {
	BmpBitmap *B=BXY;
	BmpRectangle *roi=&roiXY;
	int ix2 = ix - roi->x0,
	    iy2 = iy - roi->y0;
	if(BMP_IN_RANGEB(ix2,iy2,B))
	{
	  LongInt i=BMP_XYOFF(B,ix2,iy2);
	  for(int k=0;k<nChan;k++) B->dataFloat[k][i] = data[k];
	}
      }
      if(xIsSlice)
      {
	BmpBitmap *B=BZY;
	BmpRectangle *roi=&roiZY;
	int iz2 = iz - roi->x0,
	    iy2 = iy - roi->y0;
	if(BMP_IN_RANGEB(iz2,iy2,B))
	{
	  LongInt i=BMP_XYOFF(B,iz2,iy2);
	  for(int k=0;k<nChan;k++) B->dataFloat[k][i] = data[k];
	}
      }

    }  // isSlice
  } // full resolution cells


rcCatch:
  
  G3D_DeleteGrid(&G0);
  
  if(rcThis)
  {
    BmpDeleteBitmap(&BXZ);
    BmpDeleteBitmap(&BZY);
    BmpDeleteBitmap(&BXY);
  }

  if(pOutXZ) *pOutXZ = BXZ;
  if(pOutXY) *pOutXY = BXY;
  if(pOutZY) *pOutZY = BZY;

  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::getSumSlices( 
    	        int chan,
    		BmpBitmap *BZY, BmpBitmap *BXZ, BmpBitmap *BXY )
{

  TRC(("%s \n",UT_FUNCNAME));

  UtTimer tm;

  TIMER_START(&tm);

  int chanSum = CH_FLOAT_2;

  BmpBitmap *bitmaps[3] = { BZY, BXZ, BXY };
  
  for(int iNormDir=0;iNormDir<3;iNormDir++) // Loop over normal direction
  {
    BmpBitmap *B=bitmaps[iNormDir];

    memset( B->dataFloat[0], 0, (B->sx*B->sy)*sizeof(float));

    setChannel(0.0f,chanSum);

    const int *perm = axesPermutation[iNormDir];

    UT_ASSERT0(iNormDir==perm[0]);

    // Sum over blocks along normal dir
    //
    // i = normal direction
    // k,j = tangential directions
    //
    #pragma omp parallel for collapse(2)
    for(int bk=0; bk<_sg0->gridResBlk()[perm[2]]; bk++)
    for(int bj=0; bj<_sg0->gridResBlk()[perm[1]]; bj++)
    {
      int biLastLevel0 = -1;
      // Sum blocks along normal dir
      for(int bi=0; bi<_sg0->gridResBlk()[perm[0]]; bi++)
      {
	int bx=bi, by=bj, bz=bk;
	permuteCoords(perm,bx,by,bz);
	int bid = _sg0->getBlockIndex(bx,by,bz),
	    level = getBlockInfo(bid)->level;
	UT_ASSERT2(level>=0&&level<_nLevels);

	// only consider finest resolution
	if(level>0) 
	{
	  continue;  
	}

	SparseGrid<CellFlags> *sg=getFlagsChannel(CH_CELL_FLAGS,level);
	float *data = getFloatBlockDataPtr(chan,bid,level,0,0,0),
	      *dataSum = getFloatBlockDataPtr(chanSum,bid,level,0,1,1);

	// Get initial block-slice from previous level0-block
	if(biLastLevel0>=0)
	{
	  // Initialize first slice in block from predecessor (level-0) block
	  int bx2=biLastLevel0, by2=bj, bz2=bk;
	  permuteCoords(perm,bx2,by2,bz2);
	  int bid2 = _sg0->getBlockIndex(bx2,by2,bz2),
	      level2 = getBlockInfo(bid2)->level;	
	  UT_ASSERT2(level2==0);
	  SparseGrid<float> *sg2=getFloatChannel(chanSum,level2);
	  float *dataSum2 = sg2->getBlockDataPtr(bid2);

	  for(int vk=0; vk<sg->bsx(); vk++)
	  for(int vj=0; vj<sg->bsx(); vj++)
	  {	
	    int vx=0, vy=vj, vz=vk;
	    permuteCoords(perm,vx,vy,vz);	
	    int vid = sg->getVoxelInBlockIndex(vx,vy,vz);

	    vx=sg2->bsx()-1; vy=vj; vz=vk;
	    permuteCoords(perm,vx,vy,vz);	
	    int vid2 = sg2->getVoxelInBlockIndex(vx,vy,vz);

	    dataSum[vid] = dataSum2[vid2];
	  }
	}
	biLastLevel0 = bi;

	int bsx = sg->bsx();
	for(int vk=0; vk<bsx; vk++)
	for(int vj=0; vj<bsx; vj++)
	{	
	  double sum=0;
	  for(int vi=0; vi<bsx; vi++)
	  {
	    int vx=vi, vy=vj, vz=vk;
	    permuteCoords(perm,vx,vy,vz);	
	    int vid = sg->getVoxelInBlockIndex(vx,vy,vz);
	    sum = vi==0 ? dataSum[vid] : sum + data[vid];
	    dataSum[vid] = sum;
	  }
	}
      } // normal blocks (sum)
    
      if(biLastLevel0>=0)
      {
	int bx=biLastLevel0, by=bj, bz=bk;
	permuteCoords(perm,bx,by,bz);
	int bid = _sg0->getBlockIndex(bx,by,bz),
	    level = getBlockInfo(bid)->level;	
	UT_ASSERT2(level==0);
	SparseGrid<float> *sg=getFloatChannel(chanSum,level);
	float *dataSum = sg->getBlockDataPtr(bid);
	int bsx = sg->bsx();

	for(int vk=0; vk<bsx; vk++)
	for(int vj=0; vj<bsx; vj++)
	{	
	  int vx=bsx-1, vy=vj, vz=vk;
	  permuteCoords(perm,vx,vy,vz);	
	  int vid = sg->getVoxelInBlockIndex(vx,vy,vz);
	  double sum = dataSum[vid];

	  int y = (bj*bsx+vj)<<level,
	      z = (bk*bsx+vk)<<level;
	  if(iNormDir==0) std::swap(y,z);
	  UT_ASSERT0(BMP_IN_RANGEB(y,z,B));
	  B->dataFloat[0][BMP_XYOFF(B,y,z)] = sum;		
	}
      }
    } // tangential blocks 
  
  
      /*{
      static int pnl[3]={-1,-1,-1};
      char chbuf[200];
      sprintf(chbuf,"sum_slices_%d",iNormDir);
      double slicePos[3]={0.5,.5,.5};
      visualizeSlices( MiGetAuxPanel2(&pnl[iNormDir],chbuf),
		  chanSum,
		  IP_NEAREST,
		  slicePos, VIS_OPT_TEST
		  );
    }*/

  
  
  } // Normal direction 

  resetChannel(chanSum);

  TIMER_STOP(&tm);
  TRC(("%s: CPU %.2f sec, %.0f voxels/sec)\n",UT_FUNCNAME,
      (double)TIMER_DIFF_MS(&tm)/1000.,
      _nActCells[0]/(double)(TIMER_DIFF_MS(&tm)/1000.0)));
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::visualizeDenseField( 
    	PnlPanel *pnl, 
	float **velField,  // dense grid data 
	uint8_t *obstacles,
	double *slicePos,
	int doVisDistField
    ) 
{
  int nChan=3;

  TRC(("%s slice=%g,%g,%g\n",UT_FUNCNAME,
	slicePos[0],slicePos[1],slicePos[2]));

  SparseGrid<Vec3Float> *sg0 = _sparseGrids[0].vec3_1[0];

  G3D_Grid G0 = 
    G3D_CreateGrid( sg0->sx(), sg0->sy(), sg0->sz(),0);

  BmpBitmap *BXZ=NULL,*BXY=NULL,*BZY=NULL;
  BXZ = BmpNewBitmap(sg0->sx(),sg0->sz(),BMP_XYZ|BMP_GREY|BMP_CLEAR);
  BXY = BmpNewBitmap(sg0->sx(),sg0->sy(),BMP_XYZ|BMP_GREY|BMP_CLEAR);
  BZY = BmpNewBitmap(sg0->sz(),sg0->sy(),BMP_XYZ|BMP_GREY|BMP_CLEAR);

  for(int k=0;k<2;k++)
  {
    BmpGetUshortChannel(BXZ,k,BMP_CLEAR);
    BmpGetUshortChannel(BXY,k,BMP_CLEAR);
    BmpGetUshortChannel(BZY,k,BMP_CLEAR);
  }

  for(int k=3;k<3+2;k++)
  {
    BmpGetFloatChannel(BXZ,k,BMP_CLEAR);
    BmpGetFloatChannel(BXY,k,BMP_CLEAR);
    BmpGetFloatChannel(BZY,k,BMP_CLEAR);
  }

  #pragma omp parallel for schedule( dynamic )
  for(int iz=0;iz<sg0->sz();iz++) 
  for(int iy=0;iy<sg0->sy();iy++)
  for(int ix=0;ix<sg0->sx();ix++)
  {
    int bx,by,bz;
    sg0->getBlockCoords(ix,iy,iz, bx,by,bz);
    UT_ASSERT0(sg0->blockCoordsInRange(bx,by,bz));
    LongInt bid = sg0->getBlockIndex(bx,by,bz);
    UT_ASSERT2(bid>=0&&bid<sg0->nBlocks());
    int level = _blockmap[0][bid].level;
    UT_ASSERT0(!(level<0||level>_nLevels-1));	      
    Level *lev = &_sparseGrids[level];
    int ix1=ix>>level,
	  iy1=iy>>level,
	  iz1=iz>>level;
    CellFlags  cellFlags = lev->cellFlags[0]->getValue(ix1,iy1,iz1);
    uint16_t blockFlags = _blockmap[0][bid].flags;

    LongInt gidx = G3D_INDEX(G0,ix,iy,iz);

    level = (level+1)*(obstacles[gidx]?-1:1);

    int ix0 = (ix==(int)(sg0->sx()*(double)slicePos[0])),
	iy0 = (iy==(int)(sg0->sy()*(double)slicePos[1])),
	iz0 = (iz==(int)(sg0->sz()*(double)slicePos[2]));

    int yIsSlice = iy0,
	zIsSlice = iz0,
	xIsSlice = ix0;

    float dist=-1;
    if(xIsSlice||yIsSlice||zIsSlice)
    {
      if(doVisDistField)
      {
        float pos[3]={ix+0.5f,iy+0.5f,iz+0.5f};
        dist = interpolateScalar<IP_LINEAR>(CH_DIST_FINECOARSE,pos,0) / 1024.0;
      }
    }

    if(yIsSlice)
    {
      LongInt i=BMP_XYOFF(BXZ,ix,iz);
      for(int k=0;k<3;k++) 
	BXZ->dataFloat[k][i] = velField[k][gidx];
      BXZ->dataGrey[i] = level; 
      BXZ->dataUshort[0][i] = cellFlags;
      BXZ->dataUshort[1][i] = blockFlags;
      if(doVisDistField) BXZ->dataFloat[2][i] = dist;
    }
    if(zIsSlice)
    {
      LongInt i=BMP_XYOFF(BXY,ix,iy);
      for(int k=0;k<3;k++) 
	BXY->dataFloat[k][i] = velField[k][gidx];
      BXY->dataGrey[i] = level;
      BXY->dataUshort[0][i] = cellFlags;
      BXY->dataUshort[1][i] = blockFlags;
      if(doVisDistField) BXY->dataFloat[2][i] = dist;
    }
    if(xIsSlice)
    {
      LongInt i=BMP_XYOFF(BZY,iz,iy);
      for(int k=0;k<3;k++) 
	BZY->dataFloat[k][i] = velField[k][gidx];
      BZY->dataGrey[i] = level;
      BZY->dataUshort[0][i] = cellFlags;
      BZY->dataUshort[1][i] = blockFlags;
      if(doVisDistField) BZY->dataFloat[2][i] = dist;
    }
  }

  MtStat gstat_[4*4];
  // Color according to refinement level
  for(int k=0;k<3;k++)  // Slices
  {
    BmpBitmap *bitmaps[3] = { BZY, BXZ, BXY };   
    BmpBitmap *B=bitmaps[k];
    for(int j=0;j<nChan;j++)  // grey channels
    {
      MtStat *gstat = gstat_+j*4+1+k; 
      int igchan = j;
      MtGetStat_f(B->sx*B->sy,B->dataFloat[igchan],gstat);
      if(gstat->span<MT_NUM_EPS) gstat->span=1;

      #pragma omp parallel for
      for(LongInt i=0;i<B->sx*B->sy;i++)
      {
	int level = B->dataGrey[i], obst=FALSE;	
	unsigned cellFlags = B->dataUshort[0][i],
		 blockFlags = B->dataUshort[1][i];
	USE(cellFlags);
	USE(blockFlags);

	if(level<0)
	{
	  level = -level;
	  obst = TRUE;
	}
	level -= 1;
	UT_ASSERT0(level>=0&&level<_nLevels);

	obst = cellFlags & CELL_SOLID;

	int cr=0,cg=0,cb=0;
	if(level==0)
	{
	  cr=255; cg=200; cb=50;
	  if(obst) { cr=220; cg=220; cb=220; }
	}
	else if(level==1)
	{
	  cr=80; cg=255; cb=150; 
	  if(obst) { cr=150; cg=150; cb=150; }
	}
	else if(level==2)
	{
	  cr=170; cg=170; cb=255;
//	  cr=100; cg=50; cb=255;
	  if(obst) { cr=80; cg=80; cb=80; }
	}
	double g = B->dataFloat[igchan][i];
	g = (g-gstat->min)/gstat->span;
	UT_ASSERT0(!(g<-0.001||g>1.001));
	MT_CLAMP(g,0,1);
	if(j==1) g = sqrtf(g);
	g = 0.2+0.8*g;
	MT_CLAMP(g,0,1);

	((uint32_t*)B->dataFloat[igchan])[i] = 

#if 1
	  // 
	  // Color level
	  //
	  BMP_MKRGB(  (int)(g*cr), (int)(g*cg), (int)(g*cb) );
#elif 0
	  //
	  // Color cell flags
	  //	 
	#if 1
	  (cellFlags & (CELL_SOLID|
	              (CELL_COARSE_FINE|CELL_FINE_COARSE)|
		      CELL_DOM_BORDER|
		      CELL_BLK_BORDER)) == 0 ?
	    BMP_MKRGB(255*g,255*g,255*g) :
	    BMP_MKRGB( (int)((cellFlags&CELL_SOLID)?240:0)*g,
		       (int)(((cellFlags&(CELL_COARSE_FINE|CELL_FINE_COARSE))||
			      (cellFlags&CELL_DOM_BORDER))?240:0)*g,
		       (int)((cellFlags&CELL_BLK_BORDER)?240:0)*g	     	      
	      );
	#else
	  (cellFlags & (CELL_SOLID|
	              CELL_FINE_COARSE|
		      CELL_COARSE_FINE)) == 0 ?
	    BMP_MKRGB(255*g,255*g,255*g) :
	    BMP_MKRGB( (int)((cellFlags&CELL_SOLID)?240:0)*g,
		       (int)((cellFlags&CELL_FINE_COARSE)?240:0)*g,
		       (int)((cellFlags&CELL_COARSE_FINE)?240:0)*g	     	      
	      );
	#endif
#elif 0
	  //
	  // Color obstacles and dom-boundaries
	  //	 
	  (cellFlags & (CELL_SOLID|
		        CELL_DOM_BORDER)) == 0 ?
	    BMP_MKRGB(255*g,255*g,255*g) :
	    BMP_MKRGB( (int)((cellFlags&CELL_SOLID)?240:0)*g,
		       (int)((cellFlags&CELL_DOM_BORDER)?240:0)*g,
		       50*g	     	      
	      );
#elif 1
	  //
	  // Color block flags
	  //	 
	  (blockFlags & (BLK_RES_BORDER)) == 0 ?
	    BMP_MKRGB(255*g,255*g,255*g) :
	    BMP_MKRGB( (int)((blockFlags&BLK_RES_BORDER)?240:0)*g,
		       (int)((level==1)?240:0)*g,
		       (int)((level==2)?240:0)*g);
#elif 0
	  //
	  // Black/White
	  //
	  BMP_MKRGB(  (int)((g-0.15)*255*1.1),  
	      	      (int)((g-0.15)*255*1.1), 
		      (int)((g-0.15)*255*1.1) );
#endif
      }
    }
  }
  {
    BmpBitmap *B=NULL;

    PnlVisualize3DFieldSlices( pnl, NULL, 
			   BZY->dataFloat, BXZ->dataFloat, BXY->dataFloat,
			   nChan, sg0->sx(), sg0->sy(), sg0->sz(), 
			   slicePos[0],slicePos[1],slicePos[2], gstat_, 
			   PNL_RGB | PNL_PROVIDE_SLABSTATS, 
			   &B 
			   );
#if 0
    char path[UT_MAXPATHLEN+1];
    sprintf(path, "visu.png");
    TRC(("saving image: '%s'\n",path));
    int rc = BmpSaveBitmapPNG( B, path, NULL, BMP_PNG);
    if(rc) 
    {
      TRCERR(("BmpSaveBitmapPNG('%s') failed %d",path,rc));
    }
#endif
  }

  G3D_DeleteGrid(&G0);
  BmpDeleteBitmap(&BXZ);
  BmpDeleteBitmap(&BZY);
  BmpDeleteBitmap(&BXY);
}

/*=========================================================================
 *
 * 	
 * 	Testing
 *
 *
 * =======================================================================*/

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::tstCreateRefinementMap( int refinementTyp,
    					    int *blockmap,
					    double x0, double y0, double z0
					    )
{
  SparseGrid<Vec3Float> *sg0 = _sparseGrids[0].vec3_1[0];

  int *blockmapTmp=NULL;
  ALLOCARR_(blockmapTmp,int32_t,sg0->nBlocks());

  #pragma omp parallel for
  for(LongInt i=0;i<sg0->nBlocks();i++)
  {
    blockmap[i] = INVALID_LEVEL;
  }
 
  //
  // Create block refinement map
  //
  #pragma omp parallel for collapse(3)
  for(int bz=0;bz<sg0->nbz();bz++) 
  for(int by=0;by<sg0->nby();by++)
  for(int bx=0;bx<sg0->nbx();bx++)
  {
    LongInt bid = sg0->getBlockIndex(bx,by,bz);
    double bxyzMin = MIN(sg0->nbx(),MIN(sg0->nby(),sg0->nbz()));
    double x = (bx+0.5)/(double)bxyzMin,
	   y = (by+0.5)/(double)bxyzMin,
	   z = (bz+0.5)/(double)bxyzMin;


    UT_ASSERT0(bid>=0&&bid<sg0->nBlocks());
    
    if(refinementTyp==1)
    {
    double r = sqrt(MT_SQR(x-x0)+MT_SQR(y-y0)+MT_SQR(z-z0)); 
      blockmap[bid] = blockmapTmp[bid] = (r < 1./2.5) || 
      			(bx==sg0->nbx()/2&&bz==sg0->nbz()/2) ? 0 : 
      			_nLevels-1;   
    }
    else if(refinementTyp==2)  // MSBG usage ~ 60%
    {
      Vec3Float center(0.2,0.2,.2);
      double r = sqrt(MT_SQR(x-center.x)+MT_SQR(y-center.y)+MT_SQR(z-center.z)); 

      blockmap[bid] = blockmapTmp[bid] = 
		//(r < .5) ? 0 : 
		(r < .8) ? 0 : 
			  _nLevels-1;   
    }
    else if(refinementTyp==3)  // MSBG USAGE ~ 10%
    {
      Vec3Float center(0.2,0.2,.2);
      double r = sqrt(MT_SQR(x-center.x)+MT_SQR(y-center.y)+MT_SQR(z-center.z)); 

      blockmap[bid] = blockmapTmp[bid] = 
		//(r < .5) ? 0 : 
		(r < .4) ? 0 : 
			  _nLevels-1;  
    }
    else if(refinementTyp==0)
    {
      blockmap[bid]=blockmapTmp[bid]=0;   // for checking completely the non-adaptive case
    }
    else
    {
      UT_ASSERT0(FALSE);
    }

    // Add isolated refined blocks
    int bx0=(sg0->nbx()-3), 
	by0=(sg0->nby()-3), 
	bz0=sg0->nbz()/2;

    if( ABS(bx-bx0)<=0 && ABS(by-by0)<=0 && ABS(bz-bz0) <=0 )  
    {
      blockmapTmp[bid] = blockmap[bid] = 0;
    }
    // Refined block is diagonally adjacent
    if( ABS(bx-(bx0-1))<=0 && ABS(by-(by0-1))<=0 && ABS(bz-(bz0-1))<=0 )  
    {
      blockmapTmp[bid] = blockmap[bid] = 0;
    }

  }

#if 0
  // 
  // Post process (regularize) the refinement map
  //
  {
    UT_ASSERT0(_nLevels==3);  // for now only works with three levels

    LongInt nActVoxels = 0, nBlocks0=0, nBlocks1=0, nBlocks2=0;
    #pragma omp parallel for collapse(2) \
			    reduction(+:nActVoxels,nBlocks0,nBlocks1,nBlocks2)
    for(int bz=0;bz<sg0->nbz();bz++) 
    for(int by=0;by<sg0->nby();by++)
    for(int bx=0;bx<sg0->nbx();bx++)
    {
	LongInt bid = sg0->getBlockIndex(bx,by,bz);
	int level = blockmapTmp[bid];
	UT_ASSERT0(level>=0&&level<_nLevels);
#if 1
	if(level>1)
	{
	  int level2=_nLevels-1;
	  for(int dz=-1;dz<=1;dz++)
	  for(int dy=-1;dy<=1;dy++)
	  for(int dx=-1;dx<=1;dx++)
	  {
	    if(dx==0&&dy==0&&dz==0) continue;
	    int bx2=bx+dx,
		by2=by+dy,
		bz2=bz+dz;
	    if(!sg0->blockCoordsInRange(bx2,by2,bz2)) continue;
	    LongInt bid2 = sg0->getBlockIndex(bx2,by2,bz2);
	    UT_ASSERT2(bid2>=0&&bid2<sg0->nBlocks());
	    level2 = MIN(level2, blockmapTmp[bid2]);	  
	  }
	  if(level2 == 0) level=1;
	}
#endif
	blockmap[bid] = level;

	{
	  int m = (sg0->bsx())>>level;
	  nActVoxels += m*m*m;
	  switch(level)
	  {
	    case 0: nBlocks0++; break;
	    case 1: nBlocks1++; break;
	    case 2: nBlocks2++; break;
	    default:
	      UT_ASSERT0(FALSE);
	  }
	}
      }

      LongInt nTotalVoxels0 = sg0->nBlocks() * 
			      (LongInt)sg0->bsx()*sg0->bsx()*sg0->bsx();

      UT_ASSERT0(nActVoxels<=nTotalVoxels0);
      UT_ASSERT0(nBlocks0+nBlocks1+nBlocks2==sg0->nBlocks());
	
#ifdef UT_ASSERT_LEVEL_2
    #pragma omp parallel for
    for(LongInt i=0;i<sg0->nBlocks();i++) 
    {
      int level = blockmap[i];
      UT_ASSERT2(level>=0 && level<=_nLevels-1);
    }
#endif
  }
#endif
  regularizeRefinementMap( blockmapTmp, 
      			    blockmap );

  FREEMEM(blockmapTmp);
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int MultiresSparseGrid::testBlockIterator( )
{
  int chanSrc = CH_FLOAT_1,
      chanDst = CH_FLOAT_2;

  UtTimer tm;

  fillTestPatternChannel( 3, chanSrc, 0, 0, OPT_ALL_CELLS );
  {
    static int pnl=-1;
    visualizeSlices( MiGetAuxPanel2(&pnl,"Test scalar field"),
      chanSrc, SBG::IP_NEAREST, NULL, VIS_OPT_TEST, 0 );
  }

    TIMER_START(&tm);

    int nIter=200;
  for(int i=0;i<nIter;i++)
  {
    if(i>0) std::swap( chanSrc,chanDst);

    processFloatChannel( chanSrc, chanDst, 0, OPT_NEIGHBORHOOD,
      [&]( BlockIterator& bit, const float& valSrc, float& valDst,
	  LoopReduction* reduction )
      {
        valDst = valSrc*0.01;      
     } );  
  }

    TIMER_STOP(&tm);
    TRCP(("testBlockIterator: CPU %.2f sec, %.0f voxels/sec)\n",
	(double)TIMER_DIFF_MS(&tm)/1000.,
	(nIter*(double)getNumActCells())/(double)(TIMER_DIFF_MS(&tm)/1000.0)));

  {
    static int pnl=-1;
    visualizeSlices( MiGetAuxPanel2(&pnl,"Test processed scalar field"),
      chanDst, SBG::IP_NEAREST, NULL, VIS_OPT_TEST, 0 );
  }

  return MI_OK;
}

/*=========================================================================
 *
 *
 *  Misc Utility
 *
 *
 * =======================================================================*/

static void cycle_array(float* arr, int size) {
   float t = arr[0];
   for(int i = 0; i < size-1; ++i)
      arr[i] = arr[i+1];
   arr[size-1] = t;
}

//Given two signed distance values (line endpoints), determine what fraction of a connecting segment is "inside"
float lineFractionInside(float phi_left, float phi_right) {
   if(phi_left < 0 && phi_right < 0)
      return 1;
   if (phi_left < 0 && phi_right >= 0)
      return phi_left / (phi_left - phi_right);
   if(phi_left >= 0 && phi_right < 0)
      return phi_right / (phi_right - phi_left);
   else
      return 0;
}

//Given four signed distance values (square corners), determine what fraction of the square is "inside"
float areaFractionInside(float phi_bl, float phi_br, float phi_tl, float phi_tr) 
{
   
   int inside_count = (phi_bl<0?1:0) + (phi_tl<0?1:0) + (phi_br<0?1:0) + (phi_tr<0?1:0);
   float list[] = { phi_bl, phi_br, phi_tr, phi_tl };

   if(inside_count == 4)
      return 1;
   else if (inside_count == 3) {
      //rotate until the positive value is in the first position
      while(list[0] < 0) {
         cycle_array(list,4);
      }

      //Work out the area of the exterior triangle
      float side0 = 1-lineFractionInside(list[0], list[3]);
      float side1 = 1-lineFractionInside(list[0], list[1]);
      return 1 - 0.5f*side0*side1;
   }
   else if(inside_count == 2) {
      
      //rotate until a negative value is in the first position, and the next negative is in either slot 1 or 2.
      while(list[0] >= 0 || !(list[1] < 0 || list[2] < 0)) {
         cycle_array(list,4);
      } 
      
      if(list[1] < 0) { //the matching signs are adjacent
         float side_left = lineFractionInside(list[0], list[3]);
         float side_right = lineFractionInside(list[1], list[2]);
         return  0.5f*(side_left + side_right);
      }
      else { //matching signs are diagonally opposite
         //determine the centre point's sign to disambiguate this case
         float middle_point = 0.25f*(list[0] + list[1] + list[2] + list[3]);
         if(middle_point < 0) {
            float area = 0;

            //first triangle (top left)
            float side1 = 1-lineFractionInside(list[0], list[3]);
            float side3 = 1-lineFractionInside(list[2], list[3]);

            area += 0.5f*side1*side3;

            //second triangle (top right)
            float side2 = 1-lineFractionInside(list[2], list[1]);
            float side0 = 1-lineFractionInside(list[0], list[1]);
            area += 0.5f*side0*side2;
            
            return 1-area;
         }
         else {
            float area = 0;

            //first triangle (bottom left)
            float side0 = lineFractionInside(list[0], list[1]);
            float side1 = lineFractionInside(list[0], list[3]);
            area += 0.5f*side0*side1;

            //second triangle (top right)
            float side2 = lineFractionInside(list[2], list[1]);
            float side3 = lineFractionInside(list[2], list[3]);
            area += 0.5f*side2*side3;
            return area;
         }
         
      }
   }
   else if(inside_count == 1) {
      //rotate until the negative value is in the first position
      while(list[0] >= 0) {
         cycle_array(list,4);
      }

      //Work out the area of the interior triangle, and subtract from 1.
      float side0 = lineFractionInside(list[0], list[3]);
      float side1 = lineFractionInside(list[0], list[1]);
      return 0.5f*side0*side1;
   }
   else
      return 0;

}

#if 0
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static inline void setBcSpecVal( SparseGrid<CellFlags> *sgFlags, 
    				 int ib, int iv,
    			         float& val )
{
  UT_ASSERT2(ib>=0&&ib<sgFlags->nBlocks());	
  UT_ASSERT2(iv>=0&&iv<sgFlags->nVoxelsInBlock());	

  SBG::Block<CellFlags> *block = sgFlags->getBlock(ib); 	  	    
  UT_ASSERT2(block);
  CellFlags flags = block->_data[iv];
  if( flags & (CELL_OUT_OF_DOMAIN | CELL_SOLID) )
  {
    val = SPEC_FLOAT_NEUMANN;
  }
  else if( flags & CELL_OPEN )
  {
    val = SPEC_FLOAT_DIRICHLET;
  }
}
#endif

void testSetBcSpecValsSIMD()
{
    uint16_t flags[] = { CELL_SOLID, 0, CELL_AIR, CELL_BLK_BORDER,
      			 CELL_OUT_OF_DOMAIN, 0, CELL_BLK_BORDER, CELL_AIR, 
			 CELL_AIR, 
			 4711,4711,4711,4711,4711,4711,4711,4711,
			 4711,4711,4711,4711,4711,4711,4711,4711,
			 4711,4711,4711,4711,4711,4711,4711,4711,
			 4711,4711,4711,4711,4711,4711,4711,4711
    };

    //
    // Vec4f
    //

    Vec4f X,X0;
    
    X = Vec4f(1,2,3,4); X0=X;
    setBcSpecValsSIMD<Vec4f>( flags, CELL_SOLID|CELL_OUT_OF_DOMAIN,
				     CELL_AIR, 
			      X );		
    TRCP(("%g %g %g %g -> %g %g %g %g\n",
	  	X0[0],X0[1],X0[2],X0[3],
	  	X[0],X[1],X[2],X[3] ));

    X = Vec4f(5,6,7,8); X0=X;
    setBcSpecValsSIMD<Vec4f>( flags+4, CELL_SOLID|CELL_OUT_OF_DOMAIN,
				     CELL_AIR, 
			      X );		
    TRCP(("%g %g %g %g -> %g %g %g %g\n",
	  	X0[0],X0[1],X0[2],X0[3],
	  	X[0],X[1],X[2],X[3] ));

    X = Vec4f(9,-1,-1,-1); X0=X;
    setBcSpecValsSIMD<Vec4f>( flags+8, CELL_SOLID|CELL_OUT_OF_DOMAIN,
				     CELL_AIR, 
			      X );		
    TRCP(("%g %g %g %g -> %g %g %g %g\n",
	  	X0[0],X0[1],X0[2],X0[3],
	  	X[0],X[1],X[2],X[3] ));


    //
    // Vec4f
    //
    {

    Vec8f X,X0;
    
    X = Vec8f(1,2,3,4,5,6,7,8); X0=X;
    setBcSpecValsSIMD<Vec8f>( flags, CELL_SOLID|CELL_OUT_OF_DOMAIN,
				     CELL_AIR, 
			      X );		
    TRCP(("%g %g %g %g %g %g %g %g -> %g %g %g %g %g %g %g %g\n",
	  	X0[0],X0[1],X0[2],X0[3],X0[4],X0[5],X0[6],X0[7],
	  	X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7] ));

    X = Vec8f(9,-2,-3,-4,-5,-6,-7,-8); X0=X;
    setBcSpecValsSIMD<Vec8f>( flags+8, CELL_SOLID|CELL_OUT_OF_DOMAIN,
				     CELL_AIR, 
			      X );		
    TRCP(("%g %g %g %g %g %g %g %g -> %g %g %g %g %g %g %g %g\n",
	  	X0[0],X0[1],X0[2],X0[3],X0[4],X0[5],X0[6],X0[7],
	  	X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7] ));

    }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/

// 
// Force template instancing 
//
//
template
void MultiresSparseGrid::setChannelConstValue<float>( 
    int chan, int constval, float val, int levelMg );

template
void MultiresSparseGrid::setChannelConstValue<double>( 
    int chan, int constval, double val, int levelMg );

template
void MultiresSparseGrid::setChannelConstValue<Vec3Float>( 
    int chan, int constval, Vec3Float val, int levelMg );

template
void MultiresSparseGrid::setChannelConstValue<Vec3Uint16>( 
    int chan, int constval, Vec3Uint16 val, int levelMg );

template
void MultiresSparseGrid::setChannelConstValue<uint16_t>( 
    int chan, int constval, uint16_t val, int levelMg );

template
void MultiresSparseGrid::setChannelConstValue<uint8_t>( 
    int chan, int constval, uint8_t val, int levelMg );

} // namespace MSBG

