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

namespace MSBG 
{
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
inline int isInRange( float *pos_, const Vec4f& min, Vec4f& max )
{
  Vec4f pos;
  pos.load_partial(3,pos_);
  return !( _mm_movemask_ps( pos < min || pos > max )  );
}

//#define AO_UPPER_HEMISPHERE_ONLY 
//#define AO_DO_JITTER_1
#define AO_DO_JITTER_2
//#define AO_MULTISCALE


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template <typename T> void clampZero( T& val );
template<> void clampZero(float& f)
{
  if(f<MT_NUM_EPS) f=0.0f;
}
template<> void clampZero(Vec3Float& v)
{
  UT_ASSERT0(FALSE);
}

template <typename T> void clampOne( T& val );
template<> void clampOne(float& f)
{
  if(f>1.0f) f=1.0f;
}
template<> void clampOne(Vec3Float& v)
{
  UT_ASSERT0(FALSE);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template<typename Data_T>
void MultiresSparseGrid::applyChannelPdeFast( 
    						   int chPhi, int chTmp, int chTmp2,
						   std::vector<int> *blockList,
    						   int laplTyp,
						   int doConstrZeroOne,
    						   float densThr,
						   float densThrTol1,
					           float densThrTol2,
    						   int numIter, 
						   float gridSpacing,
						   float timeStep,
						   float velMagThr,
						   MultiresSparseGrid *msgVel,
						   float lambda,
						   float lambdaSDF,
						   float taubinStrength,
						   float redist,
						   float anisDiffAlpha,
						   int veldep_on,
						   float *veldep_fit,
						   int visOpt
						   )
{
  //#define USE_RB_SEP_LISTS
  #ifdef USE_RB_SEP_LISTS
  std::vector<int> blockListsPerCol[2];
  for(int i=0; i<(int)blockList->size();i++)
  {
    int bid = (*blockList)[i];
    Vec4i bpos = _sg0->getBlockCoordsById(bid);
    int bx=vget_x(bpos),
	by=vget_y(bpos),
	bz=vget_z(bpos);
    int icol = ((bx+by+bz) & 1);
    UT_ASSERT0(icol>=0&&icol<2);
    blockListsPerCol[icol].push_back(bid);
  }

  #endif

  bool doUseRedBlack = false,
       doUseColor8 = false;
  if( laplTyp < 0)
  {
    laplTyp *= -1;
    if(laplTyp & OPT_8_COLOR_SCHEME)
    {
      FLAG_RESET( laplTyp, OPT_8_COLOR_SCHEME );
      doUseColor8 = true;
    }
    else
    {
      doUseRedBlack=true;
    }
  }


  TRCP(("doUseRedBlack=%d doUseColor8=%d typ=%d\n",
	doUseRedBlack,doUseColor8,laplTyp));

  bool doTest = false;
  if( laplTyp & 1024 )
  {
    FLAG_RESET( laplTyp, 1024 );
    doTest = true;
  }

  int chPhiTmp = chTmp,
      chPhi0 = CH_NULL, //CH_FLOAT_3;
      chPhiPrevRK2 = CH_NULL,
      chMask = CH_NULL;

  float nbDist = 3;
  float maxCFL = 0.5f;

  int doPreProcessSDF=FALSE,
      doPostProcessSDF=FALSE;

  if(laplTyp==11)
  {
    doPreProcessSDF = doPostProcessSDF = TRUE;
  }

  float epsComprRange = 0.0f;
  if(numIter<0)
  {
    numIter = -numIter;
    epsComprRange = 0.01;
  }

  int nMaxIter0 = numIter, 
      nMaxIter = nMaxIter0;
  if((nMaxIter % 2) == 1) nMaxIter++;

  TRC(("%s: %d %d:%d %d %g:%g:%g %d/%d %g %g %g %g %g %g %g vd=%d\n",
	UT_FUNCNAME,chPhi,doUseRedBlack,doUseColor8,laplTyp,densThr,densThrTol1,densThrTol2,
	nMaxIter0,nMaxIter,
	gridSpacing,timeStep,lambda,lambdaSDF,taubinStrength,anisDiffAlpha,redist,
	veldep_on));

  //float hInv = 1.0f/gridSpacing;    

  int redistFreq = redist>MT_NUM_EPS ? ceil(redist) : 0;

  lambda = 0.0f;

  int doSoftConstraints = 0;  
  if( lambda>MT_NUM_EPS )
  {
    doSoftConstraints = 1;
  }
  else if(lambda<-MT_NUM_EPS)
  {
    doSoftConstraints = 2;
    lambda *= -1;
  }

  int do2ndOrderTimeStepping = FALSE;
  if( timeStep < 0 )
  {
    UT_ASSERT0(FALSE);
    UT_ASSERT2(!redistFreq);
    do2ndOrderTimeStepping = TRUE;
    timeStep *= -1;
    chPhiPrevRK2 = CH_FLOAT_TMP_3;
    resetChannel(chPhiPrevRK2);
  }


  int doCorrSDF = fabsf(lambdaSDF)>MT_NUM_EPS;

  int doLaplTaubin = FALSE;
  if(fabsf(taubinStrength)>MT_REAL_EPSILON)
  {
    UT_ASSERT2(!redistFreq);
    doLaplTaubin = TRUE;
  }

  int doHardConstraints = !doSoftConstraints && 
    	(densThr>MT_NUM_EPS || densThrTol1>MT_NUM_EPS);

  if(doHardConstraints||doSoftConstraints)
  {
    if(veldep_on)
    {
      UT_ASSERT_NR(FALSE);
      return;
    }
    chPhi0 = chTmp2;
  }
  else if(veldep_on)
  {
    chMask = chTmp2;
  }


  LongInt nActBytesMoved=0;

  UtTimer tm;
  TIMER_START(&tm);

  //prepare(chPhi);

  int doVisualize = visOpt&1;

  #ifdef DBG_MEM_USAGE
  {
    TRC(("MSBG memory applyLaplacianSmoothing\n"));
    showMemoryUsage(2);
  }
  #endif

  //
  // Main loop over iterations
  //
  
  float dt = timeStep; //0.01;

  SparseGrid<Data_T> *sgPhi = getChannel<Data_T>(this,chPhi);
  UT_ASSERT0( sgPhi->hasData());

  if(epsComprRange)
  {
    UT_ASSERT0(FALSE);
    sgPhi->setEmptyValue( renderDensFromFloat_<Data_T>( epsComprRange ) );
    sgPhi->setFullValue( renderDensFromFloat_<Data_T>( 1.0f - epsComprRange ) );

    ThrRunParallel<int> ( sgPhi->nBlocks(), nullptr,        
      [&](int &tls, int tid, int bid)
      {        
	if(sgPhi->isValueBlock(bid))
	{
	  Data_T *data = sgPhi->getBlockDataPtr(bid,1);
	  for(int i=0; i<sgPhi->nVoxelsInBlock();i++) 
	  {
	    data[i] =  renderDensFromFloat_<Data_T>( 
	      epsComprRange + 1.f/(1.f+2.f*epsComprRange) * 
	      	renderDensToFloat_<Data_T>(data[i]));
	  }
	}
      } );
  }

  if(chPhi0) 
  {
    getChannel<Data_T>(this,chPhi0,0)->prepareDataAccess();
    getChannel<Data_T>(this,chPhi,0)
      ->copyDataTo(getChannel<Data_T>(this,chPhi0,0));
  }

  int chSrc,chDst;

  if(doUseRedBlack||doUseColor8)
  {
    chSrc = chPhi;
    chDst = chSrc;
  }
  else
  {
    chDst = chPhiTmp;
    chSrc = chPhi;
  }

  SparseGrid<Data_T> *sgSrc = getChannel<Data_T>(this,chSrc,0),
    		     *sgDst = getChannel<Data_T>(this,chDst,0);

  if(!(doUseRedBlack||doUseColor8))
  {
    // Prepare temp. destination
    sgDst->reset();
    sgDst->prepareDataAccess();
    sgSrc->copyConstValues( sgDst );
    ThrRunParallel<int> ( sgDst->nBlocks(), nullptr,        
      [&](int &tls, int tid, int bid)
      {
	if(!sgSrc->isValueBlock(bid)) sgSrc->copyConstBlock( sgDst, bid );
      } );
  }


  std::swap( chSrc, chDst );

  LongInt nActCells=0;

  //
  // Create halo blocks for each thread
  //
  const int haloSize = laplTyp==2 ? 2 :
    		       1;
  HaloBlockSet haloBlocks( this, haloSize, nMaxThreads );

  #ifdef DBG_MEM_USAGE
  {
    TRC(("MSBG memory applyLaplacianSmoothing before main loop\n"));
    showMemoryUsage(2);
  }
  #endif

  if(laplTyp==10)
  {
    UT_ASSERT0(numIter==1);
    UT_ASSERT0(!doHardConstraints);
    UT_ASSERT0(!doSoftConstraints);
  }
  else
  {
    UT_ASSERT0(gridSpacing==1.0f);
  }


  ////////////////////////////////////////////////////////////////////////////
  //
  //    Main iterations
  //
  ////////////////////////////////////////////////////////////////////////////
  float maxChange=-1e20, maxChangeLast=-1e20;
  float dtau = maxCFL * timeStep;
  int nMaxIterAuto = ceil(nbDist/dtau);

  if(nMaxIter<0) 
  {
    UT_ASSERT0(laplTyp==11);
    nMaxIter=nMaxIterAuto;
  }

  TRC(("nbDist=%g, maxCFL=%g, dtau=%g, nMaxIter=%d\n",
	nbDist,maxCFL,dtau,nMaxIter));

  for(int iIter=0;; iIter++)
  {
    #if 1
    UtTimer tm2;
    TIMER_START(&tm2);
    #endif

    float maxChange=-1e20;
    #ifdef DBG_MEM_USAGE
    TRCP(("applyLaplacianSmoothing iter %d/%d\n",iIter,nMaxIter));
    #endif

    if(doVisualize && ((iIter%50==0) || iIter==0 || iIter==nMaxIter))
    {
      static int pnl=-1;
      visualizeSlices( MiGetAuxPanel2(&pnl,"Smoothed phi"), CH_NULL,
	    IP_NEAREST,NULL, 0/*VIS_OPT_TEST*/, 0,0,0, 
	    getFloatChannel(chDst,0)
	    );
    }  

    if(iIter==nMaxIter) break;

    if(!(maxChange<0) && fabsf( maxChange - maxChangeLast )< 0.01f/nbDist ) break;
    maxChangeLast = maxChange;

    std::swap( chSrc, chDst );


    if(redistFreq && (iIter % redistFreq ==0)) 
    {
      UT_ASSERT0(FALSE);
    }

    LongInt nActBytesMovedIter = 0;

    for( int pass_red_black=0; 
	 pass_red_black < (doUseColor8 ? 8 : 
	   		   doUseRedBlack ? 2: 
			   1); 
	 pass_red_black++)
    {
      std::vector<int> 
	*blockListAct =
	#ifdef USE_RB_SEP_LISTS
	doUseRedBlack ? &blockListsPerCol[pass_red_black] :
	#endif
					blockList;


    SparseGrid<Data_T>
      *sgSrc = getChannel<Data_T>(this,chSrc,0),
      *sgDst = getChannel<Data_T>(this,chDst,0),
      *sgPhi0 = chPhi0 ? getChannel<Data_T>(this,chPhi0,0) : NULL;

    nActCells=0;

    typedef struct
    {
      LongInt nActCells;
      size_t nActBytesMoved;
      float maxChange;
      Vec8ui rnds;
    } 
    ThreadLocals;

    ThrRunParallel<ThreadLocals>
    ( 
      blockListAct ? blockListAct->size() : nBlocks(),
      
      // initialize
      [&]( ThreadLocals &tls, int tid )  // Initialize
      {
        tls.nActCells = 0;
	tls.nActBytesMoved = 0;
	tls.maxChange = -1e20;
        tls.rnds = ((((uint32_t)(tid))<<8) ^ 
	    (2*(uint32_t)_totalSteps+pass_red_black));
        FMA_FastRandSeed(&tls.rnds);
      },

      // 
      // Run parallel
      //
      [&](ThreadLocals &tls, int tid, size_t iTask) 
      {
        int bid;

	if( !blockListAct )
	{
          bid = iTask;
	  if(!sgSrc->isValueBlock(bid)) return;
	}
	else
	{
          bid = (*blockListAct)[iTask];
	  UT_ASSERT0(sgSrc->isValueBlock(bid));
	}

	#ifndef USE_RB_SEP_LISTS
	if( doUseRedBlack )
	{
	  int bx,by,bz;
	  sgSrc->getBlockCoordsById(bid, bx, by, bz );
	  if( ((bx+by+bz) & 1) == pass_red_black ) return;
	}
        #endif

	if( doUseColor8 )
	{
	  UT_ASSERT0(!doUseRedBlack);
	  int bx,by,bz;
	  sgSrc->getBlockCoordsById(bid, bx, by, bz );
	  if(getBlockColor8(bx,by,bz) != pass_red_black) return;
	}

	const int BSX=sgSrc->bsx();
	UT_ASSERT0(sgSrc->bsx()==BSX);

	BlockInfo *bi=getBlockInfo(bid);
	int level=bi->level;
	
	if(level!=0) return;

	UT_ASSERT2(getChannel<Data_T>(this,chSrc,level)==sgSrc);

	Data_T 
	  *dataDst = sgDst->getBlockDataPtr(bid,1,0),
	  *dataPhi0 = sgPhi0 ? sgPhi0->getBlockDataPtr(bid) : NULL;

	// odd number of iterations: last round = copy back 
	if(!(doUseRedBlack||doUseColor8) && ( (iIter==nMaxIter-1) && ((nMaxIter0 % 2)==1)))
	{
	  Data_T *dataSrc = sgSrc->getBlockDataPtr(bid);
	  for(int i=0;i<sgDst->nVoxelsInBlock();i++) dataDst[i] = dataSrc[i];
	  return;
	}

	tls.nActCells += sgDst->nVoxelsInBlock();

	//#define BSX_ MSBG_RENDERDENS_BSX
	int bsx = sgSrc->bsx();
	const int BSX_ = bsx;
	UT_ASSERT0(bsx==BSX_);
        constexpr int SIMDW=8;
	UT_ASSERT0( (bsx&(SIMDW-1))==0);
	//UT_ASSERT0(bsx==16);

	int hsx = bsx+2*haloSize,
	    hsx2 = hsx*hsx;

	// Fill thread local Halo block (assumed in CPU cache)
	float *dataHaloBlock;

	{
	  dataHaloBlock = doTest ?
	    haloBlocks.fillHaloBlock_<Data_T,1,1>( chSrc, bid, 0, 
							   tid, OPT_BC_NEUMANN )[0] : 
	    haloBlocks.fillHaloBlock_<Data_T>( chSrc, bid, 0, 
							   tid, OPT_BC_NEUMANN )[0];
	}

	tls.nActBytesMoved += hsx*hsx*hsx*sizeof(Data_T) * 2;

	// Loop over each voxel in block


	#ifdef F
	#error "F already defined"
	#endif

	#ifdef TEST_ISPC_LAPSM
	uint16_t *dataDstVeri=NULL;
	if constexpr (sizeof(Data_T)==2 )
	{
	  UT_ASSERT0(!doHardConstraints);

	  uint16_t *dataDstAct = dataDst;

	  if(laplTyp==4)
	  {
	    ispc::ispc_meancurv_smooth_halo_block( dt, 
						    bsx, haloSize, dataHaloBlock,
						    dataDstAct );
	    if(!dataDstVeri) return;
	  }
	}

	#endif // TEST_ISPC_LAPSM

	/*-----------------------------------------------------------------*/
	/* 								   */
	/*-----------------------------------------------------------------*/

	for(int vz=0,vid=0;vz<bsx;vz++) 
	for(int vy=0;vy<bsx;vy++) 
	for(int vx=0, 
	    hid=MT_GXYZ(hsx,hsx2, haloSize, haloSize+vy, haloSize+vz);
	    vx<bsx; vx+=SIMDW,vid+=SIMDW,hid+=SIMDW) 
	{
	  const float * __restrict ph=dataHaloBlock+hid;
	  Vec8f phi = Vec8f().load(ph),
		phiPrev = phi,
		phi0;

	  if( doHardConstraints||doSoftConstraints )
	  {
	    phi0 = renderDensToFloat_simd8<Data_T>(dataPhi0+vid);
	  }
	  else phi0 = 0.0f;

	  if(laplTyp==4)
	  {
	    //
	    // Mean Curvature (old)
	    //
	    #undef HSX_
	    #define HSX_ (BSX_+1*2)
	    #undef F_
	    #define F_(dx,dy,dz) Vec8f().load((ph+(dx)+(dy)*(HSX_)+(dz)*(HSX_)*(HSX_)))
	    UT_ASSERT2(BSX_==bsx&&HSX_==hsx);


	    Vec8f f0 = F_(0,0,0),
		  f1 = F_(-1,0,0), f2 = F_(1,0,0),
		  f3 = F_(0,-1,0), f4 = F_(0,1,0),
		  f5 = F_(0,0,-1), f6 = F_(0,0,1);

	    Vec8f  fx = 0.5f*(f2-f1),
		   fy = 0.5f*(f4-f3),
		   fz =	0.5f*(f6-f5);

	    Vec8f fxx = f1 + f2 - 2*f0,
		  fyy = f3 + f4 - 2*f0,
		  fzz = f5 + f6 - 2*f0;

	    Vec8f 
	      fy_right = (F_(1,1,0) - F_(1,-1,0) ) * 0.5f,
	      fy_left = (F_(-1,1,0) - F_(-1,-1,0) ) * 0.5f,		  
	      fz_right = (F_(1,0,1) - F_(1,0,-1) ) * 0.5f,
	      fz_left = (F_(-1,0,1) - F_(-1,0,-1) ) * 0.5f,	    
	      fz_up = (F_(0,1,1) - F_(0,1,-1) ) * 0.5f,
	      fz_down = (F_(0,-1,1) - F_(0,-1,-1) ) * 0.5f,		 	  	    
	      fyx = (fy_right - fy_left ) * 0.5f,
	      fzx = (fz_right - fz_left ) * 0.5f,
	      fzy = (fz_up - fz_down ) * 0.5f;

  	    Vec8f H =  ( fy*fy+fz*fz )*fxx + 
	      ( fx*fx+fz*fz )*fyy + 
	      ( fx*fx+fy*fy )*fzz 
	      - 2.0f * (fx*fy*fyx + 
		        fx*fz*fzx + 
		        fy*fz*fzy );

	    Vec8f gradMagSq = fx*fx + fy*fy + fz*fz;

	    H = select(gradMagSq>1e-7f,H/gradMagSq,0.0f);

	    if( doCorrSDF )
	    {
	      UT_ASSERT0(FALSE);
	    }


	    if(doSoftConstraints==1)
	    {
	      H += lambda * sqrt( gradMagSq ) * (phi0-densThr);
	    }
	    else if(doSoftConstraints==2)
	    {
	      H += lambda * sqrt( gradMagSq ) * 
				  	MT_POW3((phi0-densThr)/densThrTol1);
	    }

	    Vec8f D;


	    D = dt * H;

	    D = min(max(D,-0.1f),0.1f); // simple limiter for more stability

	    if(doLaplTaubin)
	    {
	      phi = iIter&1 ? phi + D : phi - taubinStrength*D;
	    }
	    else
	    {
	      phi += D;
	    }
	  }
	  else if(laplTyp==2)
	  {
	    //
	    // Bi-laplacian smoothing
	    //
	    #undef HSX_
	    #define HSX_ (BSX_+2*2)
	    #undef F_
	    #define F_(dx,dy,dz) Vec8f().load((ph+(dx)+(dy)*(HSX_)+(dz)*(HSX_)*(HSX_)))
	    UT_ASSERT2(BSX_==bsx&&HSX_==hsx);

	    
	    Vec8f F0 = F_(0,0,0);
	    Vec8f H = 42.f * F_(0,0,0) 
	      
		    - 12.f * ( F_(1,0,0) + F_(0,1,0) + F_(0,0,1) + 
			     F_(-1,0,0) + F_(0,-1,0) + F_(0,0,-1) )
		    
		    + 1.f * ( F_(2,0,0) + F_(0,2,0) + F_(0,0,2) + 
			     F_(-2,0,0) + F_(0,-2,0) + F_(0,0,-2) )

		    + 2.f * ( F_(1,0,1) + F_(0,1,1) + F_(-1,0,1) +
			    F_(0,-1,1) + F_(1,0,-1) + F_(0,1,-1) + 
			    F_(-1,0,-1) + F_(0,-1,-1) + F_(1,1,0) + 
			    F_(-1,1,0) + F_(-1,-1,0) + F_(1,-1,0) );	

	    Vec8f  GX = 0.5f*(F_(1,0,0)-F_(-1,0,0)),
		   GY = 0.5f*(F_(0,1,0)-F_(0,-1,0)),
		   GZ =	0.5f*(F_(0,0,1)-F_(0,0,-1)),
		   G = sqrt( GX*GX + GY*GY + GZ*GZ );


	    Vec8f D = doSoftConstraints == 1 ? 
			      - dt * ( H*G + lambda * G * (phi0-densThr) ) :
			doSoftConstraints == 2 ? 
			      - dt * ( H*G + lambda * G * 
				  	MT_POW3((phi0-densThr)/densThrTol1) ) :
			      - dt * H * G;

	    D = min(max(D,-0.1f),0.1f); // simple limiter for more stability

	    phi = F0 + D;
	  }
	  else if(laplTyp==1)
	  {
	    // 
	    // Regular laplacian smothing
	    //
	    Vec8f F0=phi,F1,F2,F3,F4,F5,F6;
	   
	    F1.load( ph-1 );
	    F2.load( ph+1 );
	    F3.load( ph-18 );
	    F4.load( ph+18 );
	    F5.load( ph-18*18 );
	    F6.load( ph+18*18 );

	    Vec8f L = F1 + F2 +F3 + F4 + F5 + F6 -6.0f*F0,
		  GX = 0.5f*(F2-F1),
		  GY = 0.5f*(F4-F3),
		  GZ = 0.5f*(F6-F5),
		  G = sqrt( GX*GX + GY*GY + GZ*GZ );

	    Vec8f D = doSoftConstraints ? 
	      			dt * ( L*G + lambda * G * (phi0-densThr) ) :	      			
	      			dt * L*G;

	    if(doLaplTaubin)
	    {
	      phi = iIter&1 ? phi + D : phi - taubinStrength*D;
	    }
	    else
	    {
	      phi += D;
	    }

	  }
	  else if(laplTyp==3)
	  {
	    // 
	    // Simple diffusion
	    //
	    Vec8f F0=phi,F1,F2,F3,F4,F5,F6;	  
	    F1.load( ph-1 );
	    F2.load( ph+1 );
	    F3.load( ph-18 );
	    F4.load( ph+18 );
	    F5.load( ph-18*18 );
	    F6.load( ph+18*18 );

	    Vec8f S = (1.0f/6.0f) * (F1 + F2 +F3 + F4 + F5 + F6),
		  D = dt*(S-phi);
	    phi += D;

	    if(doSoftConstraints)
	    {
	      D = lambda*(phi0-phi);
	      phi += D;
	    }	     
	  }
	  else if(laplTyp==14 || laplTyp==15 )  // Jacobi iteration version
	  {
	    // 
	    // Anisotropic diffusion (Jacobi version)
	    // cf. Perona, Malik "Scale-Space and Edge Detection Using Anisotropic Diffusion"
	    //
	    //
	    #undef HSX_
	    #define HSX_ (BSX_+1*2)
	    #undef F_
	    #define F_(dx,dy,dz) Vec8f().load((ph+(dx)+(dy)*(HSX_)+(dz)*(HSX_)*(HSX_)))
	    UT_ASSERT2(BSX_==bsx&&HSX_==hsx);
            
	    Vec8f f0=F_(0,0,0),
	      	  f1=F_(-1,0,0), f2=F_(1,0,0), 
	          f3=F_(0,-1,0), f4=F_(0,1,0), 
	          f5=F_(0,0,-1), f6=F_(0,0,1);

	    #if 0
	    auto G = [&]( const Vec8f& f ) 
	    { 
	      return anisDiffAlpha < 0 ? 
			FMA_FAST_EXP( -MT_SQR(abs(f)*(1.0f/-anisDiffAlpha)) ) : 
			1.0f/(1.0f+MT_SQR(abs(f)/anisDiffAlpha));
	    };  
	    #else

	    #undef G_
	    #define G_( f ) \
	      ( anisDiffAlpha < 0 ? \
			FMA_FAST_EXP( -MT_SQR(abs(f)*(1.0f/-anisDiffAlpha)) ) : \
			1.0f/(1.0f+MT_SQR(abs(f)/anisDiffAlpha)) )

	   #endif

	    Vec8f g1 = f1-f0, g2 = f2-f0,
		  g3 = f3-f0, g4 = f4-f0,
		  g5 = f5-f0, g6 = f6-f0,
	          G1=G_(g1), G2=G_(g2), G3=G_(g3), G4=G_(g4), G5=G_(g5), G6=G_(g6);

	    if( laplTyp == 15 ) 
	    {
	      // Jacobi formulation 
	      // cf. Basran 2018 "Numerical Assessment of Anisotropic Diffusion Equation..."
	      Vec8f 
		sum = G1*f1 + G2*f2 + G3*f3 + G4*f4 + G5*f5 + G6*f6,
		wsum = G1 + G2 + G3 + G4 + G5 + G6;
	      phi = ( phi + dt * sum ) * 1.0f/(1.0f + dt * wsum);
	    }
	    else
	    {
	      // Original Perona-Malik formulation
	      Vec8f sum = G1*g1 + G2*g2 + G3*g3 + G4*g4 + G5*g5 + G6*g6;
	      phi += dt * sum;	      	    	   	    
	    }
	  }
	  else if(laplTyp==17 )  // Simple Diffusion (Jacobi)
	  {
	    Vec8f F0=phi,F1,F2,F3,F4,F5,F6;	  
	    F1.load( ph-1 );
	    F2.load( ph+1 );
	    F3.load( ph-18 );
	    F4.load( ph+18 );
	    F5.load( ph-18*18 );
	    F6.load( ph+18*18 );
		  	      
	    phi = ( phi + dt * (F1+F2+F3+F4+F5+F6) ) * 1.0f/(1.0f + dt * 6.0f );	    
	  }

	  int doHardConstraintsAct = doHardConstraints;
	  if( do2ndOrderTimeStepping )
	  {
	    UT_ASSERT0(FALSE);
	  }

	  if(doHardConstraintsAct)
	  {
	    Vec8f phiMin = max(0.0f,phi0-densThrTol1),
		  phiMax = min(1.0f,phi0+densThrTol2);	    
	    phi = max( min( phi, phiMax ), phiMin );
	  }
	  else if(doConstrZeroOne)
	  {
	    if(epsComprRange)
	      phi = min(max(phi,epsComprRange),1.0f-epsComprRange);
	    else
	      phi = min(max(phi,0.0f),1.0f);
	  }

	  UT_ASSERT2(visfinite(phi));

          #if 1

	  if constexpr (sizeof(Data_T)==1 )
	  {
	    // Store with stochastic rounding
	    renderDensFromFloat_storeSimd8_SR<Data_T>( phi,dataDst+vid, tls.rnds );
	  }
	  else
	  {
	    renderDensFromFloat_storeSimd8<Data_T>( phi,dataDst+vid );
	  }

          #else

          #if 0 
	  vstream( dataDst+vid, phi );
          #else
	  renderDensFromFloat_storeSimd8<Data_T>( phi,dataDst+vid );
          #endif

	  #endif

	} // cells      

      
	if constexpr (sizeof(Data_T)!=1 )
	{
	  _mm_mfence();  // non-temp load/store
	}

	tls.nActBytesMoved += bsx*bsx*bsx*sizeof(Data_T);

	#ifdef TEST_ISPC_LAPSM
	if constexpr (sizeof(Data_T)==2 )
	{
	  if( dataDstVeri && laplTyp==4 )
	  {
	    for(int vz=0;vz<bsx;vz++)
	    for(int vy=0;vy<bsx;vy++)
	    for(int vx=0;vx<bsx;vx++)
	    {
	      int i = MT_GXYZ(bsx,bsx*bsx,vx,vy,vz);

	      float fVeri =  renderDensToFloat_<Data_T>(dataDstVeri[i]),
		    f = renderDensToFloat_<Data_T>(dataDst[i]);
	      float e = fabsf( fVeri - f );
	      if(e>2e-5)
	      {
		TRCERR(("%d,%d,%d -> %d %d -> e=%g\n",
		      vx,vy,vz,dataDstVeri[i], dataDst[i], e));
	      }
	    }
	    FREEMEM(dataDstVeri);
	  }
	}
	#endif

      }, 

      // 
      // Reduce
      //
      [&]( ThreadLocals &tls, int tid )
      {
	nActCells += tls.nActCells;
	nActBytesMovedIter += tls.nActBytesMoved;
        maxChange = std::max(maxChange,tls.maxChange);
      }

    ); // ThrRunParallel

    } // pass_red_black;

    nActBytesMoved += nActBytesMovedIter;
    /*MT_STAT_RESULT(&dphiStat);
    char chbuf[100];
    TRCP(("dphiStat: %s\n",MT_STAT_SPRINT2(chbuf,&dphiStat)));*/

    #if 0
    TIMER_STOP(&tm2);
    TRCP(("%s: Iter=%d CPU=%.2f sec, %.0f act. voxels/sec (%.1f GB/s)\n",
	  UT_FUNCNAME,iIter,(double)TIMER_DIFF_MS(&tm2)/1000.,
	(double)nActCells/(double)(TIMER_DIFF_MS(&tm2)/1000.0),
	(nActBytesMovedIter/(double)(TIMER_DIFF_MS(&tm2)/1000.0))/((double)ONE_GB) ));
    #endif

  } // Iterations

  UT_ASSERT0(chDst == chPhi );
  TIMER_STOP(&tm);

  TRC(("maxChange=%g\n",maxChange));

  TRC(("%s: CPU=%.2f sec, %.0f act. voxels/sec/iter (%.1f GB/s)\n",
	UT_FUNCNAME,(double)TIMER_DIFF_MS(&tm)/1000.,
      (double)nMaxIter*(double)nActCells/(double)(TIMER_DIFF_MS(&tm)/1000.0),
      (nActBytesMoved/(double)(TIMER_DIFF_MS(&tm)/1000.0))/((double)ONE_GB) ));

  if(epsComprRange)
  {
    sgPhi->setEmptyValue( renderDensFromFloat_<Data_T>( 0.0f ) );
    sgPhi->setFullValue( renderDensFromFloat_<Data_T>( 1.0f ) );

    ThrRunParallel<int> ( sgPhi->nBlocks(), nullptr,        
      [&](int &tls, int tid, int bid)
      {        
	if(sgPhi->isValueBlock(bid))
	{
	  Data_T *data = sgPhi->getBlockDataPtr(bid,1);	  
	  for(int i=0; i<sgPhi->nVoxelsInBlock();i++) 
	  {
	    float f = renderDensToFloat_<Data_T>(data[i]);
	    f = std::max(std::min(f,1.0f-epsComprRange),epsComprRange);
	    f = (f - epsComprRange) * (1.0f+2*epsComprRange);
	    f = std::max(std::min(f,1.0f),0.0f);
	    data[i] = renderDensFromFloat_<Data_T>(f);
	  }
	}
      } );
  }

  invalidateForInterpolation(chPhi);

  resetChannel( chMask );
  resetChannel( chPhiTmp );
  resetChannel( chPhi0 );
  resetChannel( chPhiPrevRK2 );

}

//
// Explicit instantiations 
//
template
void MultiresSparseGrid::applyChannelPdeFast<float>( 
    						   int chPhi, int chTmp, int chTmp2,
						   std::vector<int> *blockList,
    						   int laplTyp,
						   int doConstrZeroOne,
    						   float densThr,
						   float densThrTol1,
					           float densThrTol2,
    						   int numIter, 
						   float gridSpacing,
						   float timeStep,
						   float velMagThr,
						   MultiresSparseGrid *msgVel,
						   float lambda,
						   float lambdaSDF,
						   float taubinStrength,
						   float redist,
						   float anisDiffAlpha,
						   int	veldep_on,
						   float *veldep_fit,
						   int visOpt
						   );
template
void MultiresSparseGrid::applyChannelPdeFast<uint16_t>( 
    						   int chPhi, int chTmp, int chTmp2,
						   std::vector<int> *blockList,
    						   int laplTyp,
						   int doConstrZeroOne,
    						   float densThr,
						   float densThrTol1,
					           float densThrTol2,
    						   int numIter, 
						   float gridSpacing,
						   float timeStep,
						   float velMagThr,
						   MultiresSparseGrid *msgVel,
						   float lambda,
						   float lambdaSDF,
						   float taubinStrength,
						   float redist,
						   float anisDiffAlpha,
						   int	veldep_on,
						   float *veldep_fit,
						   int visOpt
						   );
template
void MultiresSparseGrid::applyChannelPdeFast<uint8_t>( 
    						   int chPhi, int chTmp, int chTmp2,
						   std::vector<int> *blockList,
    						   int laplTyp,
						   int doConstrZeroOne,
    						   float densThr,
						   float densThrTol1,
					           float densThrTol2,
    						   int numIter, 
						   float gridSpacing,
						   float timeStep,
						   float velMagThr,
						   MultiresSparseGrid *msgVel,
						   float lambda,
						   float lambdaSDF,
						   float taubinStrength,
						   float redist,
						   float anisDiffAlpha,
						   int	veldep_on,
						   float *veldep_fit,
						   int visOpt
						   );

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::downsampleFloatChannelNew( int levelMg,  // target level
					     int chSrc,
					     int chDst,
					     unsigned options
					     )
{
  int doTimer=FALSE,doSimpleAverage = options & OPT_SIMPLE_AVERAGE;
  UtTimer tm;

  UT_ASSERT0(!(options & OPT_BC_LINEAR_EXTRAPOL));
  UT_ASSERT0( options & OPT_ALL_CELLS );
  UT_ASSERT0( UT_IMPLIES( doSimpleAverage, levelMg<_nLevels )); 

  UT_ASSERT0(levelMg<_nLevels);  // incorrect for halo blocks do not match blocks here
  
  if(_psTraceDetail)
  {
    TRC(("%s l=%d %d %d\n",UT_FUNCNAME,levelMg,chSrc,chDst));
    doTimer=levelMg<=TIMER_MAX_LEVEL+1;
  }

  if(doTimer)
  {
    TIMER_START(&tm);
  }

  UT_ASSERT0(levelMg>0);

  prepareDataAccess( chDst, -1, levelMg );

  std::unique_ptr<HaloBlockSet> haloBlocks( nullptr );
  if((levelMg <= _nLevels - 1) && !doSimpleAverage)
  {
    // Create halo blocks for each thread
    haloBlocks.reset( 
	new HaloBlockSet( this, 1 /* halo size*/, nMaxThreads ) );
  }

  // downsampling kernel 
  // (Use even sized kernel for cell-centered downsampling
  float krnDown[4*4*4];
  int krnDownWidth;

  if(doSimpleAverage)
  {
    krnDownWidth=2;
    for(int i=0;i<8; i++) krnDown[i]=1.0f/8.0f;
  }
  else
  {
    G3D_GetGausslikeKernel4x4x4(krnDown);
    krnDownWidth = 4;
  }

  #ifdef USE_EXT_BLOCK_REFS
  syncChannelsDataGenCount<float>( chSrc, levelMg-1, // For checking ext. block refs.
      			           chDst, levelMg );
  #endif

  /*for upsampling: syncChannelsDataGenCount( chSrc, levelMg+1, // For checking ext. block refs.
      			    chDst, levelMg );*/
  
  //
  // Loop over all blocks of target MG level
  //
  BlockGridInfo bgi( this, levelMg );

  std::vector<int> *blockList = NULL;

  typedef struct
  {
    float *dataTmp;
  } 
  ThreadLocals;

  ThrRunParallel<ThreadLocals>( 
    blockList ? blockList->size() : bgi._nBlocks, // number of tasks

    [&]( ThreadLocals &tls, int tid )  // Initialize locals
    {
      ALLOCARR_(tls.dataTmp,float,_sg0->nVoxelsInBlock());
    },

    // Run parallel
    [&](ThreadLocals &tls, int tid, int ibid) 
    {
      int bid = blockList ? (*blockList)[ibid] : ibid;
      int level,level0;
      getBlockLevel( bid, levelMg,
		     level0, level );

      // Get sparse grid block data pointers
      SparseGrid<CellFlags> 
	*sgFlags = getFlagsChannel(CH_CELL_FLAGS,level,levelMg),
	*sgFlagsHi = NULL; // Source channel at next finer MG level
      UT_ASSERT2( sgFlags );
      SparseGrid<float>
	*sgDataDst = getChannel<float>(this,chDst,level,levelMg),
	*sgDataSrcHi = NULL,
	*sgDataSrc = NULL;
      UT_ASSERT2( sgDataDst );


      if(levelMg<_nLevels)
      {
	BlockInfo *bi = getBlockInfo(bid,levelMg);
	if((options & OPT_EXISTING_BLOCKS_ONLY) && 
	    !(bi->flags & BLK_EXISTS)) 
	{
	  BlockInfo *bi2 = getBlockInfo(bid,levelMg-1);
          sgDataSrcHi = getChannel<float>(this,chSrc,bi2->level,levelMg-1);
	  UT_ASSERT0(!sgDataSrcHi->isValueBlock(bid));
	  sgDataSrcHi->copyConstBlock( sgDataDst, bid );
	  return;
	}
      }

      /*CellFlags *dataFlags = sgFlags->getDataPtrGen_(bid,0,0),
		*dataFlagsHi=NULL;*/

      float *dataDst = NULL,
	    *dataSrcHi = NULL;

    if(level0>=levelMg)
    {
      // No corresponding finer-res block -> simple copy from levelMg-1      
      sgDataSrcHi = getChannel<float>(this,chSrc,level,levelMg-1);      

#ifndef USE_EXT_BLOCK_REFS
      if( levelMg < getNumMgLevels() && 
	  !sgDataSrcHi->isValueBlock(bid) )
      {
        // Special handling for constant blocks of sparse grids 
	if( sgDataSrcHi->isFullBlock(bid) ) sgDataDst->setFullBlock(bid);
	if( sgDataSrcHi->isEmptyBlock(bid) ) sgDataDst->setEmptyBlock(bid);
	UT_ASSERT2( sgDataSrcHi->isFullBlock(bid) || sgDataSrcHi->isEmptyBlock(bid) );
      }
      else
#else
	if(levelMg < getNumMgLevels())
	{
	  sgDataDst->setExternalBlockDataPtr(bid, sgDataSrcHi->getBlock(bid));
	}
	else

#endif
      {
	dataDst = sgDataDst->getDataPtrGen_(bid,1,0);
	dataSrcHi = sgDataSrcHi->getDataPtrGen_(bid,0,0);
	UT_ASSERT2(!sgFlags->isDenseGrid());
	UT_ASSERT2(sgFlags->bsx()>=4);
	for(int vid=0;vid<sgFlags->nVoxelsInBlock();vid+=4) 
	{
	  Vec4f X=Vec4f().load_a(dataSrcHi+vid);
	  _mm_stream_ps( dataDst+vid, X );	  
	}
      }
    }
    else
    {
      // Downsample/restrict from finer resolution parents
      sgFlagsHi = getFlagsChannel(CH_CELL_FLAGS,level-1,levelMg-1);
      sgDataSrcHi = getChannel<float>(this,chSrc,level-1,levelMg-1);
      sgDataSrc = getChannel<float>(this,chSrc,level,levelMg-1);

      dataSrcHi = sgDataSrcHi->getDataPtrGen_(bid,0,0);
      //dataFlags = sgFlags->getDataPtrGen_(bid,0,0);
      //dataFlagsHi = sgFlagsHi->getDataPtrGen_(bid,0,0);

      dataDst = sgDataDst->getDataPtrGen_(bid,1,0);

      if(doSimpleAverage)
      {
	UT_ASSERT0(levelMg<_nLevels);
	int bsx = sgDataDst->bsx(),	  
	    bsx2 = bsx*2,
	    bsxbsx2 = bsx2*bsx2,
	    iv0=0,
	    iv1=1,
	    iv2=bsx2+0,
	    iv3=bsx2+1,
	    iv4=bsxbsx2+0,
	    iv5=bsxbsx2+1,
	    iv6=bsxbsx2+bsx2+0,
	    iv7=bsxbsx2+bsx2+1;

	for(int vz=0;vz<bsx;vz++) 
	for(int vy=0;vy<bsx;vy++) 
	for(int vx=0;vx<bsx;vx++) 
	{
	  int vid = sgDataDst->getVoxelInBlockIndex(vx,vy,vz),
	      vidHi = sgDataSrcHi->getVoxelInBlockIndex(2*vx,2*vy,2*vz);

	  float *p = dataSrcHi+vidHi;

	  float sum = 
	    p[iv0] + p[iv1] + p[iv2] + p[iv3] + 
	    p[iv4] + p[iv5] + p[iv6] + p[iv7];
	  dataDst[vid] = sum * 1.f/8.f;
	}		  	
      }
      else
      {
	
	BlockAccessor bac(&bgi,sgFlags,bid,level);

	//TRCP(("%s: L%d -> haloBlocks=%p\n",UT_FUNCNAME,levelMg,haloBlocks.get()));

	float *dataHaloBlock = NULL;
	if(haloBlocks)
	{
	  haloBlocks->fillHaloBlock_<float>( chSrc, bid, levelMg-1, 
					       tid, 
					       OPT_BC_COARSE_LEVEL );
	  dataHaloBlock = haloBlocks->getHaloBlockPtr(tid);
	}

	int dmax = krnDownWidth;

	// Loop over each voxel in block
	for(int vz=0;vz<bac._bsz;vz++)
	for(int vy=0;vy<bac._bsy;vy++)
	for(int vx=0;vx<bac._bsx;vx++)
	{
	  LongInt vid = bac.getCellIndex( vx, vy, vz );
	  float val = 0.0f;
	  int x,y,z;
	  bac.getGridCoords( vx, vy, vz,
			     x, y, z );	  
	  int x0 = 2*x,
	      y0 = 2*y,
	      z0 = 2*z,

	      ix0=x0-(dmax/2-1),
	      iy0=y0-(dmax/2-1),
	      iz0=z0-(dmax/2-1);
	      
	  if(level-1>=_nLevels) // Dense grid source level
	  {
	    for(int dz=0,ik=0;dz<dmax;dz++)
	    {
	      int z2 = iz0+dz;
	      for(int dy=0;dy<dmax;dy++)
	      {
		int y2 = iy0+dy;
		for(int dx=0;dx<dmax;dx++,ik++)
		{
		  int x2 = ix0+dx;		  
		  sgDataSrcHi->clipGridCoords(x2,y2,z2);
		  float val2 = dataSrcHi[sgDataSrcHi->getGridIndex(x2,y2,z2)];
		  float w = krnDown[ik];
		  val += val2 * w;
		}
	      }
	    }
	  }  
	  else  // Sparse grid source level
	  {
	    int bsx = sgFlagsHi->bsx(),
		hsx = bsx+2,
		hsx2 = hsx*hsx;
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
	    }
	    val = horizontal_add(S0+S1+S2+S3);
	    UT_ASSERT2(std::isfinite(val));
	  }
	  dataDst[vid] = val;   	
	} // cells      
      }
    }
  },

  [&]( ThreadLocals &tls, int tid )  // Reduce locals
  {
    FREEMEM(tls.dataTmp);
  } );

  //invalidateForInterpolation(chDst,levelMg);

#if 1
  #ifdef COARSE_MG_LEVELS_NONEXISTING_BLOCKS
  if(levelMg <= _nLevels - 1 )
  {
    UtTimer tm;
    TIMER_START(&tm);

    ThrRunParallel<int>( nBlocks(), nullptr, 
      [&](int &tls, int tid, int bid) 
      {
	BlockInfo *bi=getBlockInfo(bid,levelMg);
        SparseGrid<float> *sg = getFloatChannel(chDst,bi->level,levelMg);
	if(!sg->isEmptyBlock(bid)) return;
	int bx, by, bz;
	_sg0->getBlockCoordsById(bid,bx,by,bz);
	Vec4i bpos(bx,by,bz,0);
	bool hasNonEmptyNeighbor=false;
	for(int iNeigh=1;iNeigh<=6;iNeigh++)
	{
	  Vec4i bpos2 = bpos+_v4iStencil7[iNeigh];
	  if(!_sg0->blockCoordsInRange(bpos2)) continue;
	  int bid2=_sg0->getBlockIndex(bpos2);
	  BlockInfo *bi2=getBlockInfo(bid2,levelMg);
          SparseGrid<float> *sg2 = getFloatChannel(chDst,bi2->level,levelMg);
	  if(!sg2->isEmptyBlock(bid2))
	  {
	    hasNonEmptyNeighbor = true;
	    break;
	  }
	}
	if(!hasNonEmptyNeighbor) 
	{
	  FLAG_RESET( bi->flags, BLK_EXISTS );
	}
      });

    TIMER_STOP(&tm);
    TRCP(("CPU(determine nonexisting blocks) =%.2f sec, L=%d, %.0f blocks/sec\n",
	  (double)TIMER_DIFF_MS(&tm)/1000.,levelMg,
	(double)nBlocks()/(double)(TIMER_DIFF_MS(&tm)/1000.0)));
  }
  #endif
#endif

  if(doTimer)
  {
    TIMER_STOP(&tm);
    TRC(("CPU(downsampleVel L%d) =%.2f sec, %.0f act. voxels/sec/iter\n",levelMg,
	  (double)TIMER_DIFF_MS(&tm)/1000.,
	(double)_nActCells[levelMg]/(double)(TIMER_DIFF_MS(&tm)/1000.0)));

  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template< typename Data_T, typename DataIntern_T >
void MultiresSparseGrid::downsampleFaceDensity( int levelMg,  // target level
					     int chSrc,
					     int chDst,
					     unsigned options
					     )
{
  int doTimer=FALSE;
  UtTimer tm;
  
  UT_ASSERT0( options & OPT_SIMPLE_AVERAGE );


  UT_ASSERT0(!( options & ( OPT_BC_LINEAR_EXTRAPOL)));

  UT_ASSERT0( options & OPT_ALL_CELLS );
  UT_ASSERT0( isVecChannel( chSrc ) && isVecChannel( chDst ));
  UT_ASSERT0(levelMg>0);
  
  if(_psTraceDetail)
  {
    TRC(("%s l=%d %d %d\n",UT_FUNCNAME,levelMg,chSrc,chDst));
    doTimer=levelMg<=TIMER_MAX_LEVEL+1;
  }

  if(doTimer)
  {
    TIMER_START(&tm);
  }

  prepareDataAccess( chDst, -1, levelMg );

  bool isDenseLevel = levelMg >= getNumLevels();

  ThrRunParallel<int> ( nBlocks(), nullptr,
    [&](int &tls, int tid, int bid)
    {
      int level0,level;
      getBlockLevel(bid,levelMg,level0,level);
      BlockInfo *bi = getBlockInfo(bid,levelMg);

      if(!isDenseLevel)
      {
        if(!(bi->flags & BLK_EXISTS))
	{
	  getChannel<Data_T>(this,chDst,level,levelMg)->setEmptyBlock(bid);
	  return;
	}
      }

      bool hasFineBlock = isDenseLevel || level0<levelMg;

      SparseGrid<CellFlags> *sgFlags = getFlagsChannel(CH_CELL_FLAGS,level,levelMg);
      SparseGrid<Data_T> 
        *sgDst = getChannel<Data_T>(this,chDst,level,levelMg),
	*sgSrcHi = hasFineBlock ? 
	  getChannel<Data_T>(this,chSrc,level-1,levelMg-1):
	  getChannel<Data_T>(this,chSrc,level,levelMg-1);
      Data_T *dataDst = sgDst->getDataPtrGen_( bid, 1, 0 ),
	     *dataSrcHi = sgSrcHi->getDataPtrGen_( bid );

      if(!hasFineBlock)
      {
	// No corresponding finer-res block -> simple copy from levelMg-1  	
	for(int i=0;i<sgDst->nVoxelsInBlock();i++) dataDst[i] = dataSrcHi[i];
	return;
      }

      SBG_FOR_EACH_VOXEL_IN_BLOCK_GEN( sgFlags, bid, x, y, z, idx )
      {
	DataIntern_T f=0.0f;

	for(int iDir=0;iDir<3;iDir++)
	{
	  // Loop over 4 hi-res cells
	  for(int i=0;i<2;i++)  
	  for(int j=0;j<2;j++)
	  {
	    int d[3] = {0,i,j};	
	    int xh = x*2+d[axesPermutationInv[iDir][0]],
		yh = y*2+d[axesPermutationInv[iDir][1]],
		zh = z*2+d[axesPermutationInv[iDir][2]];
	    UT_ASSERT2(
		UT_IMPLIES(!isDenseLevel,sgSrcHi->inRange(xh,yh,zh)));
	    DataIntern_T fHi = sgSrcHi->getValueGen_(xh,yh,zh);
	    /*UT_ASSERT2( std::isfinite(f) && 
			!(f<0.0f||f>1.0f) );*/
	    f[iDir] += fHi[iDir];
	  }
	  f[iDir] *= 1.0f/4.0f;
	}
	dataDst[idx] = f;
      }
    } ); // parallel blocks

  if(doTimer)
  {
    TIMER_STOP(&tm);
    TRC(("CPU(downsampleVel L%d) =%.2f sec, %.0f act. voxels/sec/iter\n",levelMg,
	  (double)TIMER_DIFF_MS(&tm)/1000.,
	(double)_nActCells[levelMg]/(double)(TIMER_DIFF_MS(&tm)/1000.0)));

  }
}

template<>
void MultiresSparseGrid::downsampleFaceDensity<float,float>( int levelMg,  // target level
					     int chSrc,
					     int chDst,
					     unsigned options
					     )
{
  UT_ASSERT0(FALSE);  // not implemented
}
template<>
void MultiresSparseGrid::downsampleFaceDensity<uint16_t,float>( int levelMg,  // target level
					     int chSrc,
					     int chDst,
					     unsigned options
					     )
{
  UT_ASSERT0(FALSE);  // not implemented
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template<typename data_T, typename dataIntern_T >
void MultiresSparseGrid::downsampleChannel( int levelMg,  // target level
					     int chSrc,
					     int chDst,
					     unsigned options
					     )
{

#if 0   // 

  if(options & OPT_SIMPLE_AVERAGE)
  {
    downsampleChannel2x2x2<data_T>(levelMg,chSrc,chDst,options);
    return;

    /*int chTmp = isVecChannel(chDst) ? CH_VELOCITY_TMP :
      				      CH_PRESSURE;
    copyChannel(chDst,chTmp,1,-1,levelMg);
    downsampleChannel2x2x2<data_T>(levelMg,chSrc,chDst,options);
    double diffMax = compareChannel(chDst,chTmp,levelMg);
    TRCP(("diffMax=%g\n",diffMax));*/
  }

#endif

  int doTimer=FALSE,doSimpleAverage = options & OPT_SIMPLE_AVERAGE;
  UtTimer tm;
  
  if(_psTraceDetail)
  {
    TRC(("%s l=%d %d %d\n",UT_FUNCNAME,levelMg,chSrc,chDst));
    doTimer=levelMg<=TIMER_MAX_LEVEL+1;
  }

  if(doTimer)
  {
    TIMER_START(&tm);
  }

  UT_ASSERT0(levelMg>0);

  bool doUseFlagsChannel = ! (options & OPT_ALL_CELLS );
  
  prepareDataAccess( chDst, -1, levelMg );

  // downsampling kernel 
  // (Use even sized kernel for cell-centered downsampling
  float krnDown[4*4*4];
  int krnDownWidth;

  if(doSimpleAverage)
  {
    krnDownWidth=2;
    for(int i=0;i<8; i++) krnDown[i]=1.0f/8.0f;
  }
  else
  {
    G3D_GetGausslikeKernel4x4x4(krnDown);
    krnDownWidth = 4;
  }

  int doBcLinearExtrapol = options & OPT_BC_LINEAR_EXTRAPOL;
  if(doBcLinearExtrapol)
  {
    UT_ASSERT0(krnDownWidth==3);
  }

  #ifdef USE_EXT_BLOCK_REFS
  syncChannelsDataGenCount<data_T>( chSrc, levelMg-1, // For checking ext. block refs.
      			            chDst, levelMg );
  #endif

  //
  // Loop over all blocks of target MG level
  //
  BlockGridInfo bgi( this, levelMg );

  ThrRunParallel<int>( bgi._nBlocks, nullptr,
    [&](int &tls, int tid, int bid) 
    {
      int level,level0;
      getBlockLevel( bid, levelMg,
		     level0, level );

      // Get sparse grid block data pointers
      SparseGrid<CellFlags> 
	*sgFlags = getFlagsChannel(CH_CELL_FLAGS,level,levelMg),
	*sgFlagsHi = NULL; // Source channel at next finer MG level
      UT_ASSERT2( sgFlags );
      SparseGrid<data_T>
	*sgDataDst = getChannel<data_T>(this,chDst,level,levelMg),
	*sgDataSrcHi = NULL,
	*sgDataSrc = NULL;
      UT_ASSERT2( sgDataDst );

      CellFlags *dataFlags = doUseFlagsChannel ? sgFlags->getDataPtrGen_(bid,0,0) : NULL,
		*dataFlagsHi=NULL;

      data_T *dataSrcHi = NULL;

      if(level0>=levelMg)
      {
	// No corresponding finer-res block -> simple copy from levelMg-1  
	
	//#define USE_EXT_BLOCK_REFS_OLD

	sgDataSrcHi = getChannel<data_T>(this,chSrc,level,levelMg-1);      

	if( levelMg < getNumMgLevels() && 
	    !sgDataSrcHi->isValueBlock(bid) )
	{
	  // Special handling for constant blocks of sparse grids 
	  sgDataSrcHi->copyConstBlock( sgDataDst, bid );
	}
	else
	{
	  dataSrcHi = sgDataSrcHi->getDataPtrGen_(bid,0,0);
	  
	  sgFlagsHi = getFlagsChannel(CH_CELL_FLAGS,level,levelMg-1);
	  dataFlagsHi = doUseFlagsChannel ? sgFlagsHi->getDataPtrGen_(bid,0,0) : NULL;
	  
	  //UT_ASSERT2(sgFlagsHi && dataFlagsHi && dataSrcHi);
	  UT_ASSERT2(sgFlags->bsx()==sgFlagsHi->bsx());
	  UT_ASSERT2(!sgFlags->isDenseGrid());


	  if(!(options&OPT_ALL_CELLS))
	  {
	    data_T *dataDst = sgDataDst->getDataPtrGen_(bid,1,0);
	    for(int vid=0;vid<sgFlags->nVoxelsInBlock();vid++)
	      dataDst[vid] = 
		(!(dataFlags&&(dataFlags[vid] & (CELL_SOLID|CELL_VOID)))) ? 
					      dataSrcHi[vid] : 0.0f;
	  }
	  else 
	  {
	    #ifdef USE_EXT_BLOCK_REFS
	    UT_ASSERT0(sgDataDst->bsx()==sgDataSrcHi->bsx());
	    sgDataDst->setExternalBlockDataPtr(bid, sgDataSrcHi->getBlock(bid));
	    #else
	    data_T *dataDst = sgDataDst->getDataPtrGen_(bid,1,0);
	    for(int vid=0;vid<sgFlags->nVoxelsInBlock();vid++) 
	      dataDst[vid]=dataSrcHi[vid];
	    #endif
	  }
	}
      }
      else
      {
	// Downsample/restrict from finer resolution parents
	sgFlagsHi = getFlagsChannel(CH_CELL_FLAGS,level-1,levelMg-1);
	sgDataSrcHi = getChannel<data_T>(this,chSrc,level-1,levelMg-1);
	sgDataSrc = getChannel<data_T>(this,chSrc,level,levelMg-1);

	dataSrcHi = sgDataSrcHi->getDataPtrGen_(bid,0,0);
	dataFlags = doUseFlagsChannel ? sgFlags->getDataPtrGen_(bid,0,0) : NULL;
	dataFlagsHi = doUseFlagsChannel ? sgFlagsHi->getDataPtrGen_(bid,0,0) : NULL;

	SparseGrid<data_T> *sg=sgDataSrcHi;
	int bsx = sg->bsx(),
	    bsxLog2 = sg->bsxLog2(),
	    bsx2Log2 = sg->bsx2Log2(),
	    sx = sg->sx(),
	    sy = sg->sy(),
	    sz = sg->sz(),
	    nx = sg->nbx(),
	    ny = sg->nby(),
	    nz = sg->nbz(),
	    nxy = sg->nbxy(),
	    bsxMask = bsx-1;
	USE(ny);USE(nz);

	data_T *dataDst = sgDataDst->getDataPtrGen_(bid,1,0);

	// Loop over each voxel in block
	BlockAccessor bac(&bgi,sgFlags,bid,level);
	for(int vz=0;vz<bac._bsz;vz++)
	for(int vy=0;vy<bac._bsy;vy++)
	for(int vx=0;vx<bac._bsx;vx++)
	{
	  LongInt vid = bac.getCellIndex( vx, vy, vz );
	  dataIntern_T val = 0.0f;
	  //if(CELL_IS_FLUID(dataFlags[vid]))
	  if(  (options&OPT_ALL_CELLS) || 
	       !(dataFlags&&(dataFlags[vid] & CELL_SOLID)) )
	  {
	    int x,y,z;
	    bac.getGridCoords( vx, vy, vz,
			       x, y, z );	  

	    #if 0
	    if(levelMg==DBG_CELL_LEVEL_MG && 
	       x==DBG_CELL_X && y==DBG_CELL_Y && z==DBG_CELL_Z)
	    {
	      TRCERR(("DBG\n"));
	    }
	    #endif

	    int x0 = 2*x,
		y0 = 2*y,
		z0 = 2*z;

	    int dmax = krnDownWidth,ix0,iy0,iz0;

	    if(dmax==3)
	    {
	      ix0=x0-1;
	      iy0=y0-1;
	      iz0=z0-1;
	    }
	    else
	    {
	      ix0=x0-(dmax/2-1);
	      iy0=y0-(dmax/2-1);
	      iz0=z0-(dmax/2-1);
	    }
		

	    if(level-1>=_nLevels) // Dense grid source level
	    {
	      // Loop over high-res cells
	      for(int dz=0,ik=0;dz<dmax;dz++)
	      {
		int z2 = iz0+dz;
		for(int dy=0;dy<dmax;dy++)
		{
		  int y2 = iy0+dy;
		  for(int dx=0;dx<dmax;dx++,ik++)
		  {
		    int x2 = ix0+dx;
		    dataIntern_T val2;
		    

#if 0
		    sgDataSrcHi->clipGridCoords(x2,y2,z2);
		    val2 = dataSrcHi[sgDataSrcHi->getGridIndex(x2,y2,z2)];
#else
		    int isClipped;
		    sgDataSrcHi->clipGridCoords(x2,y2,z2,isClipped);
		    if(isClipped && doBcLinearExtrapol )
		    {
		      int x1=x2,
			  y1=y2,
			  z1=z2;
		      int x2=ix0+dx,
			  y2=iy0+dy,
			  z2=iz0+dz;
		      int x3 = x2<0 ? 1 : x2 > sx-1 ? sx-2 : x2,
			  y3 = y2<0 ? 1 : y2 > sy-1 ? sy-2 : y2,
			  z3 = z2<0 ? 1 : z2 > sz-1 ? sz-2 : z2;
		      UT_ASSERT2( sgDataSrcHi->inRange(x3,y3,z3));
		      dataIntern_T 
			     val1 = dataSrcHi[sgDataSrcHi->getGridIndex(x1,y1,z1)],
			     val3 = dataSrcHi[sgDataSrcHi->getGridIndex(x3,y3,z3)];
		      val2 = val1*2 - val3;
		    }
		    else
		    {
		      val2 = dataSrcHi[sgDataSrcHi->getGridIndex(x2,y2,z2)];
		    }
#endif

		    float w = krnDown[ik];
		    val += val2 * w;
		  }
		}
	      }
	    }  // Dense grid source level
	    else
	    {
#if 0
	      if(blockNoBorderFast && ( (ix0 & bsxMask) <= bsx-4))
	      {
		// TODO
		// Fast SIMD case
	      }
	      else
#endif
	      {
		for(int k=0,iKrn=0;k<dmax;k++)
		{
		  int z2 = iz0+k, clampZ=FALSE;
		  if(unlikely(z2<0))
		  {
		    z2=0; clampZ=TRUE;
		  }
		  else if(unlikely(z2>sz-1))
		  {
		    z2=sz-1; clampZ=TRUE;
		  }
		  int bz2 = z2 >> bsxLog2,
		      vz2 = z2 & bsxMask,
		      bidZ = bz2*nxy,
		      vidZ = (vz2<<bsx2Log2);

		  for(int j=0;j<dmax;j++)
		  {
		    int y2 = iy0+j, clampYZ=clampZ;
		    if(unlikely(y2<0))
		    {
		      y2=0; clampYZ=TRUE;
		    }
		    else if(unlikely(y2>sy-1))
		    {
		      y2=sy-1; clampYZ=TRUE;
		    }
		    int by2 = y2 >> bsxLog2,
			vy2 = y2 & bsxMask,
			bidY = bidZ+by2*nx,
			vidY = vidZ + (vy2<<bsxLog2); 

		    for(int i=0;i<dmax;i++,iKrn++)
		    {
		      int x2 = ix0+i,clampXYZ = clampYZ;
		      if(unlikely(x2<0))
		      {
			x2=0; clampXYZ=TRUE;
		      }
		      else if(unlikely(x2>sx-1))
		      {
			x2=sx-1; clampXYZ=TRUE;
		      }

		      double w = krnDown[iKrn];

		      dataIntern_T f;

		      if(clampXYZ && doBcLinearExtrapol)
		      {
			int x1=x2, y1=y2, z1=z2;
			int x2=ix0+i,
			    y2=iy0+j,
			    z2=iz0+k;
			int x3 = x2<0 ? 1 : x2 > sx-1 ? sx-2 : x2,
			    y3 = y2<0 ? 1 : y2 > sy-1 ? sy-2 : y2,
			    z3 = z2<0 ? 1 : z2 > sz-1 ? sz-2 : z2;
			UT_ASSERT2( sgDataSrcHi->inRange(x3,y3,z3));
			data_T *f1 = sgDataSrcHi->getValuePtr(x1,y1,z1),
			       *f3 = sgDataSrcHi->getValuePtr(x3,y3,z3);
			// fallback: take coarser level source values 
			if(!f1)  f1 = sgDataSrc->getValuePtr(x1/2,y1/2,z1/2);
			if(!f3)  f3 = sgDataSrc->getValuePtr(x3/2,y3/2,z3/2);
			UT_ASSERT2(f1&&f3)
			f = (dataIntern_T)(*f1) * 2 - (dataIntern_T)(*f3);
		      }
		      else
		      {		    		   
			int bx2 = x2 >> bsxLog2,
			    vx2 = x2 & bsxMask;
			int bid2 = bidY + bx2,
			    vid2 = vidY + vx2;


			UT_ASSERT2(bid2>=0&&bid2<_nBlocks);
			UT_ASSERT2(vid2>=0&&vid2<sg->nVoxelsInBlock());
#if 0
			SBG::Block<CellFlags> 
			  *blockFlags2 = sgFlagsHi->getBlock(bid2); 	  
			CellFlags flags2 = blockFlags2 ? 
				    blockFlags2->_data[vid2] :
				    getCellFlags( x2<<level,y2<<level,z2<<level );
#endif
			SBG::Block<data_T>* block = sg->getBlock(bid2); 	  
			if(block)
			{
			  f = block->_data[ vid2 ]; 
			}
			else
			{
			  int x2l=x2/2, y2l=y2/2, z2l=z2/2;
			  
			  UT_ASSERT2(sgDataSrc->getBlock(bid2));
			  UT_ASSERT2(sgDataSrc->inRange(x2l,y2l,z2l));
			  
			  f = sgDataSrc->getValue(x2l,y2l,z2l);
			}
		      }

		      val += f*w;
		    } // kernel x
		  } // kernel y
		} // kernel z

	      }  // non-SIMD case
	    
	    } // Sparse grid source level

	  } // Fluid cell

	  dataDst[vid] = val;   	
	} // cells      
      }
    } ); // blocks

  if(doTimer)
  {
    TIMER_STOP(&tm);
    TRC(("CPU(downsampleVel L%d) =%.2f sec, %.0f act. voxels/sec/iter\n",levelMg,
	  (double)TIMER_DIFF_MS(&tm)/1000.,
	(double)_nActCells[levelMg]/(double)(TIMER_DIFF_MS(&tm)/1000.0)));

  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::downsampleVelocity( int levelMg,  // target level
					     int chSrc,
					     int chDst )
{
  int doTimer=FALSE;
  UtTimer tm;
  
  if(_psTraceDetail)
  {
    TRC(("%s l=%d %d %d\n",UT_FUNCNAME,levelMg,chSrc,chDst));
    doTimer=levelMg<=TIMER_MAX_LEVEL+1;
  }

  if(doTimer)
  {
    TIMER_START(&tm);
  }

  UT_ASSERT0(levelMg>0);

  UT_ASSERT0(isVecChannel(chSrc)&&isVecChannel(chDst));

  if(isDenseLevel(levelMg))
  {
    getVecChannel(chDst,levelMg,levelMg)->getDenseData();
  }

  // downsampling kernel 
  // (Use even sized kernel for cell-centered downsampling
  float krnDown[4*4*4];
  int krnDownWidth=4;
  G3D_GetGausslikeKernel4x4x4(krnDown);
  krnDownWidth = 4;

  //
  // Loop over all blocks of target MG level
  //
  BlockGridInfo bgi( this, levelMg );

  #pragma omp parallel for schedule( static, 1 )
  for(int bid=0; bid<bgi._nBlocks; bid++)
  {
    int level,level0;
    getBlockLevel( bid, levelMg,
		   level0, level );

    // Get sparse grid block data pointers
    SparseGrid<CellFlags> 
      *sgFlags = getFlagsChannel(CH_CELL_FLAGS,level,levelMg),
      *sgFlagsHi = NULL; // Source channel at next finer MG level
    UT_ASSERT2( sgFlags );
    SparseGrid<Vec3Float>
      *sgDataDst = getVecChannel(chDst,level,levelMg),
      *sgDataSrcHi = NULL,
      *sgDataSrc = NULL;
    UT_ASSERT2( sgDataDst );

    CellFlags *dataFlags = sgFlags->getDataPtrGen_(bid,0,0),
	      *dataFlagsHi=NULL;

    Vec3Float 
      	  *dataDst = sgDataDst->getDataPtrGen_(bid,1,0),
	  *dataSrcHi = NULL;

    UT_ASSERT2(dataDst);

    if(level0>=levelMg)
    {
      // No corresponding finer-res block -> simple copy from levelMg-1      
      sgDataSrcHi = getVecChannel(chSrc,level,levelMg-1);      
      dataSrcHi = sgDataSrcHi->getDataPtrGen_(bid,0,0);
      
      sgFlagsHi = getFlagsChannel(CH_CELL_FLAGS,level,levelMg-1);
      dataFlagsHi = sgFlagsHi->getDataPtrGen_(bid,0,0);
      
      UT_ASSERT2(sgFlagsHi && dataFlagsHi && dataSrcHi);
      UT_ASSERT2(sgFlags->bsx()==sgFlagsHi->bsx());
      UT_ASSERT2(!sgFlags->isDenseGrid());

      for(int vid=0;vid<sgFlags->nVoxelsInBlock();vid++)
      {
	//UT_ASSERT2(dataFlagsHi[vid]==dataFlags[vid]);
	dataDst[vid] = (!(dataFlags[vid] & (CELL_SOLID|CELL_VOID))) ? 
	  				dataSrcHi[vid] : 0.0f;
      }
    }
    else
    {
      // Downsample/restrict from finer resolution parents
      sgFlagsHi = getFlagsChannel(CH_CELL_FLAGS,level-1,levelMg-1);
      sgDataSrcHi = getVecChannel(chSrc,level-1,levelMg-1);
      sgDataSrc = getVecChannel(chSrc,level,levelMg-1);

      dataSrcHi = sgDataSrcHi->getDataPtrGen_(bid,0,0);
      dataFlags = sgFlags->getDataPtrGen_(bid,0,0);
      dataFlagsHi = sgFlagsHi->getDataPtrGen_(bid,0,0);

      SparseGrid<Vec3Float> *sg=sgDataSrcHi;
      int bsx = sg->bsx(),
	  bsxLog2 = sg->bsxLog2(),
	  bsx2Log2 = sg->bsx2Log2(),
	  sx = sg->sx(),
	  sy = sg->sy(),
	  sz = sg->sz(),
	  nx = sg->nbx(),
	  ny = sg->nby(),
	  nz = sg->nbz(),
	  nxy = sg->nbxy(),
	  bsxMask = bsx-1;
      USE(ny);USE(nz);

      int blockNoBorderFast = FALSE;

      if(level-1<MSBG_MAX_LEVELS && bsx>=8)
      {
        BlockInfo *bi = getBlockInfo(bid,levelMg-1);
	blockNoBorderFast = ( !(bi->flags & (BLK_DOM_BORDER | 
			      	      BLK_COARSE_FINE | BLK_FINE_COARSE)));
      }

      // Loop over each voxel in block
      BlockAccessor bac(&bgi,sgFlags,bid,level);
      for(int vz=0;vz<bac._bsz;vz++)
      for(int vy=0;vy<bac._bsy;vy++)
      for(int vx=0;vx<bac._bsx;vx++)
      {
	LongInt vid = bac.getCellIndex( vx, vy, vz );
	Vec3Float val(0.0f);
	//if(CELL_IS_FLUID(dataFlags[vid]))
	if(!(dataFlags[vid] & CELL_SOLID))
	{
	  int x,y,z;
	  bac.getGridCoords( vx, vy, vz,
			     x, y, z );	  

	  UT_ASSERT2(krnDownWidth==4);

	  int ix0=2*x-1,
	      iy0=2*y-1,
	      iz0=2*z-1;

	  if(level-1>=_nLevels) // Dense grid source level
	  {
	    // Loop over high-res cells
	    for(int dz=0,ik=0;dz<4;dz++)
	    {
	      int z2 = iz0+dz;
	      for(int dy=0;dy<4;dy++)
	      {
		int y2 = iy0+dy;
		for(int dx=0;dx<4;dx++,ik++)
		{
		  int x2 = ix0+dx;
		  Vec3Float val2;
		  sgDataSrcHi->clipGridCoords(x2,y2,z2);
		  val2 = dataSrcHi[sgDataSrcHi->getGridIndex(x2,y2,z2)];
		  float w = krnDown[ik];
		  val += val2 * w;
		}
	      }
	    }
	  }  // Dense grid source level
	  else
	  {
#if 0
	    if(blockNoBorderFast && ( (ix0 & bsxMask) <= bsx-4))
	    {
	      // TODO
	      // Fast SIMD case
	    }
	    else
#endif
	    {
	      for(int k=0,iKrn=0;k<4;k++)
	      {
		int z2 = iz0+k;
		MT_CLAMP( z2, 0, sz-1 );
		int bz2 = z2 >> bsxLog2,
		    vz2 = z2 & bsxMask,
		    bidZ = bz2*nxy,
		    vidZ = (vz2<<bsx2Log2);

		for(int j=0;j<4;j++)
		{
		  int y2 = iy0+j;
		  MT_CLAMP( y2, 0, sy-1 );
		  int by2 = y2 >> bsxLog2,
		      vy2 = y2 & bsxMask,
		      bidY = bidZ+by2*nx,
		      vidY = vidZ + (vy2<<bsxLog2); 

		  for(int i=0;i<4;i++,iKrn++)
		  {
		    int x2 = ix0+i;
		    MT_CLAMP( x2, 0, sx-1 );
		    int bx2 = x2 >> bsxLog2,
			vx2 = x2 & bsxMask;
		    int bid2 = bidY + bx2,
			vid2 = vidY + vx2;

		    double w = krnDown[iKrn];

		    Vec3Float f;

		    UT_ASSERT2(bid2>=0&&bid2<_nBlocks);
		    UT_ASSERT2(vid2>=0&&vid2<sg->nVoxelsInBlock());

		    SBG::Block<Vec3Float>* block = sg->getBlock(bid2); 	  
		    if(block)
		    {
		      f = block->_data[ vid2 ]; 
		    }
		    else
		    {
		      int x2l=x2/2, y2l=y2/2, z2l=z2/2;
		      
		      UT_ASSERT2(sgDataSrc->getBlock(bid2));
		      UT_ASSERT2(sgDataSrc->inRange(x2l,y2l,z2l));
		      
		      f = sgDataSrc->getValue(x2l,y2l,z2l);
		    }		    
		    val += f*w;
		  } // kernel x
		} // kernel y
	      } // kernel z

	    }  // non-SIMD case
	  
	  } // Sparse grid source level

	} // Fluid cell

	dataDst[vid] = val;   	
      } // cells      
    }
  } // blocks

  if(levelMg<=0) invalidateForInterpolation(chDst);

  if(doTimer)
  {
    TIMER_STOP(&tm);
    TRC(("CPU(downsampleVel L%d) =%.2f sec, %.0f act. voxels/sec/iter\n",levelMg,
	  (double)TIMER_DIFF_MS(&tm)/1000.,
	(double)_nActCells[levelMg]/(double)(TIMER_DIFF_MS(&tm)/1000.0)));

  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::buildGaussianPyramidNew( 
    		int chan,
    		int levelMgMax,
		unsigned options )
{ 
  if(levelMgMax<0) levelMgMax=_mgLevels-1;
  TRC(("%s %d %d\n",UT_FUNCNAME,chan,levelMgMax));

  UT_ASSERT0(levelMgMax<_nLevels);

  UtTimer tm;
  UT_ASSERT0(levelMgMax<_mgLevels);

  int levelMg;

  TIMER_START(&tm);

  for( levelMg=0; levelMg < levelMgMax;  levelMg++ )
  {
//#define VERIFY_DOWNSAMPLE_CHANNEL
#ifdef VERIFY_DOWNSAMPLE_CHANNEL
    copyChannel( chan, CH_FLOAT_4, TRUE,-1,levelMg);
#endif

    downsampleFloatChannelNew(levelMg+1,chan,chan,options|OPT_ALL_CELLS);
    
    if(options & OPT_VISUALIZE)
    {
      static int pnlNr[]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
      UT_ASSERT0(levelMg<ARRAY_LENGTH(pnlNr));
      char chbuf[200];
      sprintf(chbuf,"Gaussian pyramid L=%d",levelMg+1);
      visualizeSlices( 	MiGetAuxPanel2(&pnlNr[levelMg],chbuf),
		  chan,
		  IP_NEAREST,
		  _visSlicePos, VIS_OPT_TEST,
		  levelMg+1
		  );
    }

    #if 0
    if((options & OPT_FREE_INTERMEDIATE)  && levelMg>0 && levelMg <levelMgMax)
    {
      resetChannel(chan,0,levelMg);
    }
    #endif

#ifdef VERIFY_DOWNSAMPLE_CHANNEL
    copyChannel( chan, CH_FLOAT_3,TRUE,-1,levelMg+1); // save result in tmp2
    copyChannel( CH_FLOAT_4, chan, TRUE,-1,levelMg); // restore source
    downsampleChannel<float,float>( levelMg+1, chan, chan, OPT_ALL_CELLS );
    double diffMax = compareChannel(chan,CH_FLOAT_3,levelMg+1);
    TRCP(("%s: diffMax = %g\n",UT_FUNCNAME,diffMax));
#endif


  }
  TIMER_STOP(&tm);
 
  TRC(("%s CPU=%.2f sec, %.0f act. voxels/sec\n",UT_FUNCNAME,
      (double)TIMER_DIFF_MS(&tm)/1000.,
      (double)_nActCells[0]/(double)(TIMER_DIFF_MS(&tm)/1000.0)));
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template<typename Data_T, typename DataIntern_T>
void MultiresSparseGrid::buildGaussianPyramid( 
    		int chan,
    		int levelMgMax,
		unsigned options )
{ 
  if(levelMgMax<0) levelMgMax=_mgLevels-1;
  TRC(("%s %d %d\n",UT_FUNCNAME,chan,levelMgMax));

  UtTimer tm;
  UT_ASSERT0(levelMgMax<_mgLevels);

  int levelMg;

  TIMER_START(&tm);

  for( levelMg=0; levelMg <= levelMgMax;  levelMg++ )
  {
    if(levelMg<levelMgMax)
    {
      {
	#ifdef DOWNSAMPLE_FACEDENS_FACE_ORIENTED
	if( chan == CH_FACE_DENSITY )
	{
          downsampleFaceDensity<Data_T, DataIntern_T>( levelMg+1, chan, chan, options );
	}
	else
	#endif
	{
	  downsampleChannel<Data_T, DataIntern_T>( levelMg + 1, chan, chan, options );
	}
      }
    }

    if(options & OPT_VISUALIZE)
    {
      static int pnlNr[]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
      UT_ASSERT0(levelMg<ARRAY_LENGTH(pnlNr));
      char chbuf[200];
      sprintf(chbuf,"Gaussian pyramid L=%d",levelMg);
      visualizeSlices( 	MiGetAuxPanel2(&pnlNr[levelMg],chbuf),
		  chan,
		  IP_NEAREST,
		  _visSlicePos, 0*VIS_OPT_TEST,
		  levelMg
		  );
    }

    if(levelMg<levelMgMax)
    {
      #if 0
      if((options & OPT_FREE_INTERMEDIATE)  && levelMg>0 && levelMg <levelMgMax)
      {
	resetChannel(chan,0,levelMg);
      }
      #endif
    }

  }
  TIMER_STOP(&tm);
 
  TRC(("%s CPU=%.2f sec, %.0f act. voxels/sec\n",UT_FUNCNAME,
      (double)TIMER_DIFF_MS(&tm)/1000.,
      (double)_nActCells[0]/(double)(TIMER_DIFF_MS(&tm)/1000.0)));
}

template void MultiresSparseGrid::buildGaussianPyramid<float,float>( 
    		int chan,int levelMgMax,unsigned options );

template void MultiresSparseGrid::buildGaussianPyramid<Vec3Float,Vec3Float>( 
    		int chan,int levelMgMax,unsigned options );

template void MultiresSparseGrid::buildGaussianPyramid<Vec3Uint16,Vec3Float>( 
    		int chan,int levelMgMax,unsigned options );

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::freeGaussianPyramid( int chan )
{
  TRC(("%s\n",UT_FUNCNAME));
  int levelMgMax=_mgLevels-1;

  for( int levelMg = levelMgMax-1; levelMg >= 1; levelMg-- )
  {
    resetChannel(chan,FALSE,levelMg);
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::interpolateWithDerivs( int ipolType,
    						int chan, 
    					        Vec4f pos0,
						int levelMg,
    					        unsigned opt,
					        float *outVal,
					        Vec4f *outGrad
					   	)
{
  UT_ASSERT2(_validForInterpolation & (((int64_t)1)<<chan));

  int stenWidth = ipolType == IP_BSPLINE_QUAD ? 3 :
    		  ipolType == IP_BSPLINE_CUBIC ? 4 :
		  -1;
  UT_ASSERT0(stenWidth>0);

  pos0 = v4f_zero_to<3>(pos0);
  Vec4f pos = pos0;

  // 
  // Clip coords to global domain grid
  // and lookup the containing block 
  //
  pos = max(pos,_sg0->v4fDomMin());
  pos = min(pos,_sg0->v4fDomMax());

  Vec4i ipos = truncate_to_int(pos);

  int ix0=vget_x(ipos),
      iy0=vget_y(ipos),
      iz0=vget_z(ipos);
  UT_ASSERT2(_sg0->inRange(ix0,iy0,iz0));
  int bx,by,bz;
  _sg0->getBlockCoords(ix0,iy0,iz0, bx,by,bz);
  int bid = _sg0->getBlockIndex(bx,by,bz);

  BlockInfo *bi = getBlockInfo(bid,levelMg);
  int level = bi->level;
  UT_ASSERT2(level>=0&&level<_nLevels);

  Level *lev = &_sparseGrids[level];
  SparseGrid<float> *sg= getFloatChannel(chan,level,levelMg); 

  pos = pos0;  // pass un-clipped coords to SBG::interpolate
  	       // to enable domain-BC handling there

  // Convert coordinates to current level grid
  pos *= lev->auxInvScale;  // 1.0f/(1<<level)

  //
  // Check for resolution transition and blend two interpolants if necessary
  //
  unsigned ipOpt=0;
  if(!(bi->flags & BLK_DOM_BORDER))
  {
    ipOpt |= OPT_IP_NO_BORDER_CHECK;
  }

  if( bi->flags & BLK_FINE_COARSE )
  {
    UT_ASSERT2( ! ( (_protectedRead & (1<<((int64_t)CH_DIST_FINECOARSE)) ) ||
	            (_protectedWrite & (1<<((int64_t)CH_DIST_FINECOARSE)) )));
    
    float dist = 
      lev->distFineCoarse->interpolateShortToFloatFast(
					pos,
					OPT_IPBC_NEUMAN |
					  ((bi->flags & BLK_DOM_BORDER) ? 0 : 
						OPT_IP_NO_BORDER_CHECK)
					) / 1024.0;
    float 
      alphaFine = MtLinstep( (stenWidth/2-(.5f))*MT_SQRT3*1.01, 
			     0.99*(_dTransRes+2)/MT_SQRT3,
			    // 0.99*(_dTransRes+2-isStaggered?1.f:0.f)/MT_SQRT3,
			     dist);    
    if(alphaFine>0.00001f)
    {
      sg->interpolateWithDerivs<0>(ipolType,pos,ipOpt, 
	  			outVal,outGrad);
      if(outGrad) *outGrad *= lev->auxInvScale;
    }
    if(alphaFine<.9999f)
    {
      UT_ASSERT2(level<_nLevels-1);
      SparseGrid<float> *sgLo = getFloatChannel(chan,level+1);	
      UT_ASSERT2((!sgLo)||(sgLo&&sgLo->getBlock(bid)));
      if(sgLo)
      {
        Vec4f posLo = 0.5f*pos;
	Vec4f outGradLo;
	float outValLo;
        sgLo->interpolateWithDerivs<0>(ipolType, posLo, ipOpt, 
	    			  outVal ? &outValLo:NULL, 
				  outGrad ? &outGradLo:NULL);
	if(outGrad) outGradLo *= lev->auxInvScale*0.5;
	if(outVal) *outVal = (*outVal)*alphaFine + outValLo*(1-alphaFine);
  	if(outGrad) 
	  *outGrad = (*outGrad)*alphaFine + outGradLo*(1-alphaFine);
      }
    }
  }
  else
  {
    sg->interpolateWithDerivs<0>(ipolType,pos,ipOpt, 
			      outVal,outGrad);
    if(outGrad) *outGrad *= lev->auxInvScale;
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::getVelocityFieldDerivs( BlockIterator& bit,
	    unsigned options, 
	    double *dvx,  // OUT
	    double *dvy,
	    double *dvz
	)
{
  if(TRUE)
  {
    UT_ASSERT0(FALSE);
#if 0
    double h = (1<<bit.level)*_h0;
    int x,y,z;
    bit.sg->getGridCoords( bit.bx,bit.by,bit.bz, bit.x,bit.y,bit.z,
			   x,y,z );      
    Vec3Float p0((x+0.5f) * (1<<bit.level),
		 (y+0.5f) * (1<<bit.level),
		 (z+0.5f) * (1<<bit.level));

    int ipMeth = IP_LINEAR;
    Vec3Float vNeigh[6];

#if 0
    float offs[6][3] = { -0.5f, 0.0f, 0.0f,
      		         +0.5f, 0.0f, 0.0f,

			 0.0f, -0.5f, 0.0f,
      		         0.0f, +0.5f, 0.0f,

			 0.0f, 0.0f, -0.5f,
      		         0.0f, 0.0f, +0.5f };
    double scale = 1./(h);
#else
    float offs[6][3] = { -1.0f, 0.0f, 0.0f,
      		         +1.0f, 0.0f, 0.0f,

			 0.0f, -1.0f, 0.0f,
      		         0.0f, +1.0f, 0.0f,

			 0.0f, 0.0f, -1.0f,
      		         0.0f, 0.0f, +1.0f };
    double scale = 1./(2*h);
#endif
    for(int i=0;i<6;i++)
      interpolate( ipMeth, bit.chanData, p0+offs[i], 0, vNeigh[i] );

    for(int i=0;i<3;i++)
    {
      dvx[i] = ( vNeigh[2*i+1].x - (double)vNeigh[2*i].x ) * scale;
      dvy[i] = ( vNeigh[2*i+1].y - (double)vNeigh[2*i].y ) * scale;
      dvz[i] = ( vNeigh[2*i+1].z - (double)vNeigh[2*i].z ) * scale;
    }    
#endif
  }
  else
  {
    // Get cell neighborhod 
    Vec3Float velNeigh[7][4];

    getCellNeighborhood<8>( bit, NULL,  NULL, velNeigh );

    if(bit.cfnFlags[0].flags[0] & CELL_COARSE_FINE)
    {
      for(int i=2; i<=6; i+=2 ) 
      {
	if(bit.cfnFlags[i].transResType==COARSE_FINE)
	{
	  // Average over 4 finer-res faces 
	  Vec3Float sum(0.0f);
	  for(int j=0;j<4;j++) sum += velNeigh[i][j];
	  velNeigh[i][0] = sum*(1./4.);
	}
	else
	{
	  UT_ASSERT2(bit.cfnFlags[i].transResType==FINE_FINE);
	}
      }
    }

    double h = (1<<bit.level)*_h0,
	 scale = 1./(2*h);

    dvx[0] = ( velNeigh[2][0].x - (double)velNeigh[0][0].x ) * scale*2.;  // staggered
    dvx[1] = ( velNeigh[4][0].x - (double)velNeigh[3][0].x ) * scale; // non-staggered
    dvx[2] = ( velNeigh[6][0].x - (double)velNeigh[5][0].x ) * scale;

    dvy[0] = ( velNeigh[2][0].y - (double)velNeigh[1][0].y ) * scale;
    dvy[1] = ( velNeigh[4][0].y - (double)velNeigh[0][0].y ) * scale*2.;
    dvy[2] = ( velNeigh[6][0].y - (double)velNeigh[5][0].y ) * scale;

    dvz[0] = ( velNeigh[2][0].z - (double)velNeigh[1][0].z ) * scale;
    dvz[1] = ( velNeigh[4][0].z - (double)velNeigh[3][0].z ) * scale;
    dvz[2] = ( velNeigh[6][0].z - (double)velNeigh[0][0].z ) * scale*2.;
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
Vec3Float MultiresSparseGrid::getScalarFieldDerivs( BlockIterator& bit, 
    						    unsigned opt )
{
  // Get cell neighborhod
  
  float valNeigh[7][4]={{0}};  
  getCellNeighborhood<2>( bit, NULL,  valNeigh, NULL );

  if(bit.cfnFlags[0].flags[0] & CELL_COARSE_FINE)
  {
    for(int i=1;i<=6;i++) 
    {
      if(bit.cfnFlags[i].transResType==COARSE_FINE)
      {
	// Average over 4 finer-res faces 
	double sum=0;
        for(int j=0;j<4;j++) sum += valNeigh[i][j];
	valNeigh[i][0] = sum*(1./4.);
      }
      else
      {
	UT_ASSERT2(bit.cfnFlags[i].transResType==FINE_FINE);
      }
    }
  }

  double h = (1<<bit.level)*_h0,
	 scale = 1./(2*h);

  double vn[7];

  for(int i=1;i<=6;i++)
  {
    if(CELL_IS_FLUID_((bit.cfnFlags[i].flags[0])))
    {
      vn[i] = valNeigh[i][0];
    }
    else
    {
      if( opt & OPT_BC_NEUMANN )
      {
	int i2 = i + 2*(i&1)-1;
	UT_ASSERT2(i2>=1&&i2<=6);
	vn[i] = 2*valNeigh[0][0] - valNeigh[i2][0];
      }
    }

  }

  return Vec3Float( vn[2] - (double)vn[1], 
      		    vn[4] - (double)vn[3],	
		    vn[6] - (double)vn[5] ) * scale;
}

/*=========================================================================
 *
 *  
 *  Boundary Block Smoothing
 *  
 *
 * =======================================================================*/

#define BSX BSB_SX
#define BSY BSB_SY
#define BSZ BSB_SZ
  

  
/*-------------------------------------------------------------------------*/ 
/* 									   */
/*-------------------------------------------------------------------------*/
static inline 
void smoothBlockLineSegSIMD4f( 
    		       uint32_t i,  // index into (padded) buffer
		       float h,
    		       float *bufDiag,  
    		       float *bufB, 	
		       float *bufSrc, 	
    		       float *bufDst 	 
		       )
{
  Vec4f F0,F1,F2,F3,F4,F5,F6,
  	B,D,Y;

  float *p=bufSrc+i;    
  F0.load ( p );   // center
  F1.load ( p-1 );   // left 
  F2.load ( p+1 );   // right 
  F3.load ( p-(BSX+2) );   // down
  F4.load ( p+(BSX+2) );   // up
  F5.load ( p-(BSX+2)*(BSY+2) );   // front
  F6.load ( p+(BSX+2)*(BSY+2) );   // back

  B.load( bufB+i );  // RHS
  D.load( bufDiag+i);  // diagonal factor

  Vec4f S = ( F1+F2+F3+F4+F5+F6 ) * h;
  Y = /*(2.f/3.f) *     Assume D is already pre-multiplied with 2/3*/ 
      (B+S)*D  + (1.f/3.f) * F0;

  Y.store( bufDst+i );
}
  
/*-------------------------------------------------------------------------*/ 
/* 									   */
/*-------------------------------------------------------------------------*/
static inline 
void smoothBlockSIMD4f( 
    		       float h,
    		       float *bufDiag,  // unpadded buffer
    		       float *bufB, 	// unpadded buffer
		       float *bufSrc, 	// padded buffer (BSX+2)x(BSY+2)x(BSZ+2)
    		       float *bufDst 	// padded buffer 
		       )
{
  uint32_t i=G3D_INDEX0(BSX+2,(BSX+2)*(BSY+2),1,1,1); // starting index into padded buffer
  UT_ASSERT2(BSY==4);
  for(int dz=0;dz<BSZ;dz++)
  {
    smoothBlockLineSegSIMD4f( i,h,bufDiag, bufB, bufSrc, bufDst );
    i+=BSX+2;
    smoothBlockLineSegSIMD4f( i,h,bufDiag, bufB, bufSrc, bufDst );
    i+=BSX+2;
    smoothBlockLineSegSIMD4f( i,h,bufDiag, bufB, bufSrc, bufDst );
    i+=BSX+2;
    smoothBlockLineSegSIMD4f( i,h,bufDiag, bufB, bufSrc, bufDst );
    i+=BSX+2;

    i+=2*(BSX+2); // Skip padding lines
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/

template void MultiresSparseGrid::downsampleChannel<float,float>( 
    int levelMg, int chSrc, int chDst, unsigned options );
template void MultiresSparseGrid::downsampleChannel<Vec3Float,Vec3Float>( 
    int levelMg, int chSrc, int chDst, unsigned options );
template void MultiresSparseGrid::downsampleChannel<Vec3Uint16,Vec3Float>( 
    int levelMg, int chSrc, int chDst, unsigned options );


} // namespace MSBG

