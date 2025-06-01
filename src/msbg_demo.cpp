/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/
#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <tuple> 
#include <iostream>
#include <fstream>
#include <tbb/tbb.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <errno.h>
#include <unistd.h>

#include "msbg.h"

#include "render.h"

#define VIS_OUTPUT_DIR "out_msbg_demo"

typedef uint32_t ParticleIdx;

typedef struct
{
  Vec3Float pos;
  ParticleIdx idxNext;
}
Particle;

static uint64_t nMaxParticles=0;
static Particle *particles=NULL;

inline static ParticleIdx getParticleIdx( int nInst, int iInst, int i )
{
  return nInst * (ParticleIdx)(iInst) + (ParticleIdx)(i) 
    	 + 1  /* slot 0 reserved for NULL-idx */ ; 
}

inline static Particle *getParticle( ParticleIdx ixp )
{
  UT_ASSERT2(ixp>0 && ixp<nMaxParticles+1);
  return &particles[ixp];
}

inline static Vec4f getPosition( Particle *p )
{
  return Vec4f(p->pos[0],p->pos[1],p->pos[2],0.0f);
}

inline static void pushToListMonotonic( 
		      Particle *particle,
		      ParticleIdx idxParticle,
		      ParticleIdx *listRoot )
{
  ParticleIdx idx = 
    InterlockedExchange((LONG volatile *)(listRoot), idxParticle );
  particle->idxNext = idx;
}


static 
std::tuple< std::vector<Vec3Float>*,
  	    Vec3Float,
	    Vec3Float >
readVerticesFromPLY( const char *filename )
{
  std::ifstream file(filename);
  std::string line;
  std::vector<Vec3Float>* vertices = new std::vector<Vec3Float>(); 
  Vec3Float bbMin(1e20),
	    bbMax(-1e20);


  bool parsingHeader = true;
  int vertexCount = 0;
  int verticesParsed = 0;

    // Read through the PLY file
  while (std::getline(file, line)) {
      std::istringstream iss(line);

      // Parse the header
      if (parsingHeader) {
	  if (line.substr(0, 10) == "element ve") {
	      iss >> line >> line >> vertexCount;  // Extract vertex count
	      TRCP(("vertex count = %d\n",vertexCount));
	  } else if (line == "end_header") {
	      parsingHeader = false;  // End of the header
	  }
      } else {
	  // Parse vertices
	  if (verticesParsed < vertexCount) {
	      Vec3Float vertex;
	      iss >> vertex.x >> vertex.y >> vertex.z;
	      vertices->push_back(vertex);

	      // Update bounding box
	      for(int k=0;k<3;k++)
	      {
	        bbMin[k] = std::min(bbMin[k],vertex[k]);
	        bbMax[k] = std::max(bbMax[k],vertex[k]);
	      }
	      verticesParsed++;
	  }
      }
    }

  return std::make_tuple( vertices, bbMin, bbMax );
}


float camPos[3] = {.5f,.8f,.7f},
      camLookAt[3] = {.5,.5,.5},
      camLight[3] = {.5,5,5},
      camZoom =1.0f;
int camRes[2] = {3840,2160};

static void render_scene( int testCase, MSBG::MultiresSparseGrid *msbg, int chan, 
     			  PnlPanel *pnl )
{
  SBG::SparseGrid<uint16_t> *sg = msbg->getUint16Channel(chan,0);

  UtTimer tm;
  int sx=camRes[0], sy=camRes[1];

  TRCP(("Rendering image %dx%d.\n",sx,sy));

  RDR::RaymarchRenderer renderer(sx,sy, sg);
#if 1

  #if 1
  renderer.setCamera( {camPos[0],camPos[1],camPos[2]}, 
      		      {camLookAt[0],camLookAt[1],camLookAt[2]}, camZoom); 


//  renderer.setCamera({.5,.8,.7}, {.5,.5,.5}, 1); 
  renderer.setSunLight({camLight[0],camLight[1],camLight[2]});
  renderer.setSurfaceColor({0.8f,0.6f,0.4f});
  #else

  renderer.setCamera({.5,.8,.8}, {.5,.5,.5}, 1.f);
//  renderer.setCamera({.5,1,1}, {.5,.5,.5}, 1.5f);
  

  renderer.setSunLight({.5,5,5});
  renderer.setSurfaceColor({0.8f,0.6f,0.4f});
  #endif

#else
  renderer.setCamera({.5,.5,.65}, {.5,.5,.5}, 0.8f);
//  renderer.setCamera({0,.5,1}, {0.4,0.5,0.5}, 2.5f);
  renderer.setSunLight({.5,5,5});
  renderer.setSurfaceColor({0.8f,0.6f,0.4f});
#endif


  TIMER_START(&tm);

  BmpBitmap *B = renderer.render();

  TIMER_STOP(&tm);
  TRCP(("CPU (rendering) %.2f sec, %.0f pixels/sec)\n",
    (double)TIMER_DIFF_MS(&tm)/1000.,
    ((double)sy*sy)/(double)(TIMER_DIFF_MS(&tm)/1000.0)));	

  PnlShowBitmap2( pnl, B );  

  pnl->totalTime = 0;
  msbg->saveVisualizationBitmap( B, pnl->title );

  TRCP(("Output images saved to '%s/'.\n",VIS_OUTPUT_DIR));
  
  BmpDeleteBitmap(&B);
}

/////////////////////////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////////////////////////
int msbg_test_multires(int bsx0, int sx, int sy, int sz)
{
  using namespace MSBG;
  using namespace SBG;
  using MSBG::MultiresSparseGrid;

  /////////////////////////////////////////////////////////////////////////////////
  //
  // 	Create MSBG grid
  //
  /////////////////////////////////////////////////////////////////////////////////

  MSBG::MultiresSparseGrid *msbg = 
    MSBG::MultiresSparseGrid::create( "TEST", sx,sy,sz, bsx0,
				      -1, 0,-1, 
				      MSBG::OPT_SINGLE_LEVEL |
					MSBG::OPT_SINGLE_CHANNEL_FLOAT );

  msbg->setVisOutDir( 2, (char*)VIS_OUTPUT_DIR );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // 		Free resources 
  //
  /////////////////////////////////////////////////////////////////////////////////

  MSBG::MultiresSparseGrid::destroy(msbg);

  return 0;
}

/////////////////////////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////////////////////////
int msbg_test_sparse(int testCase, const char *basePointsFile, 
    		     int bsx0, int sx, int sy, int sz)
{
  using namespace MSBG;
  using namespace SBG;
  using MSBG::MultiresSparseGrid;

  float rParticle = 2,
	nbDist = 2;

  UtTimer tm;

  /////////////////////////////////////////////////////////////////////////////////
  //
  // 	Load base points
  //
  /////////////////////////////////////////////////////////////////////////////////

  auto [ basePoints, basePointsBMin, basePointsBMax ] = 
    readVerticesFromPLY( basePointsFile );

  Vec3Float basePointsSpan;
  float basePointsSpanMax = -1e20;
  for(int k=0;k<3;k++) 
  {
    basePointsSpan[k] = basePointsBMax[k] - basePointsBMin[k];
    UT_ASSERT0(std::isfinite(basePointsSpan[k] && 
	       basePointsSpan[k] >= 0.f));
    basePointsSpanMax = std::max(basePointsSpanMax,basePointsSpan[k]);
  }

  TRCP(("read %d points from '%s'. bbox = %g,%g,%g, spanMax=%g\n",
	basePoints->size(),basePointsFile,
	basePointsSpan.x,basePointsSpan.y,basePointsSpan.z,
	basePointsSpanMax));

  tbb::parallel_sort( basePoints->begin(),basePoints->end(), 
		      [&](const Vector3Dim<float>& p1_, 
			 	   const Vector3Dim<float>& p2_ ) -> bool
    {
      Vec3Float p1 = (p1_-basePointsBMin)/basePointsSpan,
                p2 = (p2_-basePointsBMin)/basePointsSpan;
      Vec4i ipos1 = Vec4i(p1.x,p1.y,p1.z,0),
            ipos2 = Vec4i(p2.x,p2.y,p2.z,0);
      uint64_t m1 = FMA::morton64_encode( ipos1 ),
               m2 = FMA::morton64_encode( ipos2 );
      return m1 < m2;
    } );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // 	Create MSBG grid
  //
  /////////////////////////////////////////////////////////////////////////////////

  MSBG::MultiresSparseGrid *msbg = 
    MSBG::MultiresSparseGrid::create( "TEST", sx,sy,sz, bsx0,
				      -1, 0,-1, 
				      MSBG::OPT_SINGLE_LEVEL |
					MSBG::OPT_SINGLE_CHANNEL_FLOAT );

  msbg->setVisOutDir( 2, (char*)VIS_OUTPUT_DIR );

  LONG *blockActive=NULL;
  ParticleIdx *particlesPerBlock = NULL;

  ALLOCARR0_(blockActive,LONG,msbg->nBlocks());
  ALLOCARR0_(particlesPerBlock, ParticleIdx, msbg->nBlocks());

  int chan = CH_UINT16_1;
  SparseGrid<uint16_t> *sg0 = msbg->getUint16Channel( chan, 0 );


  uint64_t nBasePoints = basePoints->size(),
	   nInst = nBasePoints;

  nMaxParticles = nInst * (uint64_t)nBasePoints;

  UT_ASSERT0( nMaxParticles+1 < 4*((uint64_t)(ONE_GB)) );  // 32 bit particle index

  TRCP(("nBasePoints=%.0f  nInst=%.0f nMaxParticles=%.0f\n",
	(double)nBasePoints,(double)nInst,(double)nMaxParticles));

  /////////////////////////////////////////////////////////////////////////////////
  //
  // 	Compute particle positions and determine active blocks 
  //
  /////////////////////////////////////////////////////////////////////////////////

  float rScan = rParticle + nbDist;

  TRCP(("rParticle=%g, rScan=%g\n",rParticle,rScan));

  ALLOCARR_(particles, Particle, nMaxParticles + 1/* slot 0 reserved for NULL-idx */);

  uint64_t nActParticles = 0;

  TIMER_START(&tm);
  {
    int bxMax = sg0->nbx()-1,
	byMax=  sg0->nby()-1,
	bzMax = sg0->nbz()-1;
    const float scale2DestBlockGrid = 1.0f/(float)sg0->bsx(),
		rScanBsx = rScan/(float)sg0->bsx();
    Vec4i bposMax(bxMax,byMax,bzMax,INT32_MAX);

    using ThreadLocals = struct
    {
      size_t nActParticles;
    };

    ThrRunParallel<ThreadLocals>( nInst, nullptr,
      [&](ThreadLocals &tls, int tid, size_t iInst) 
      {
        Vec4f baseMin(basePointsBMin.x,basePointsBMin.y,basePointsBMin.z,0);
	float //domMax = sg0->sxyzMax(),
	      baseScale = 1.0f/basePointsSpanMax,
	      scale;
	Vec4f pos0;
	scale = sg0->sxyzMin() * (testCase==2 ? 0.005 : 0.01); //0.1;
	Vec3Float &pos_ = (*basePoints)[iInst];
	Vec4f pos = sg0->sxyzMax() * baseScale * (Vec4f(pos_.x,pos_.y,pos_.z,0.f) - baseMin); // normalized
//	pos0 = 0.1+0.8*pos;
	pos0 =  sg0->sxyzMax()*.2 + 0.6*pos;

        for( int i=0; i<(int)nBasePoints; i++)
	{
	  Vec3Float &pos_ = (*basePoints)[i];
	  Vec4f pos = baseScale * (Vec4f(pos_.x,pos_.y,pos_.z,0.f) - baseMin); // normalized
	  pos = pos0 + scale * pos;
	  if(!msbg->isInDomainRange(pos)) continue;
	  Vec4i ipos = truncate_to_int(pos);
	  if( !sg0->inRange(ipos) ) continue;

	  ParticleIdx ixp = getParticleIdx( nInst, iInst, i );
	  Particle *p = getParticle(ixp);

	  p->pos = Vec3Float(pos);

	  // 
	  // Mark all blocks touched by the particle
	  //
	  Vec4f bpos = pos * scale2DestBlockGrid;	
	  Vec4i bpos1 = max(truncate_to_int(bpos-rScanBsx),0),
		bpos2 = min(truncate_to_int(bpos+rScanBsx),bposMax);
	  int bx1=vget_x(bpos1), by1=vget_y(bpos1), bz1=vget_z(bpos1);
	  int bx2=vget_x(bpos2), by2=vget_y(bpos2), bz2=vget_z(bpos2);
	  for( int bz=bz1; bz<=bz2; bz++)
	  for( int by=by1; by<=by2; by++)
	  for( int bx=bx1; bx<=bx2; bx++)
	  {
	    int bid = sg0->getBlockIndex(bx,by,bz);
	    UT_ASSERT2(bid>=0&&bid<sg0->nBlocks());
	    InterlockedIncrement((volatile LONG*)(&blockActive[bid]));
	  }

	  //
	  // Update particles-per-block list
	  //
	  int bid = sg0->getBlockIndex(sg0->getBlockCoords(ipos));	  
	  pushToListMonotonic(p, ixp, &particlesPerBlock[bid]);
	  tls.nActParticles++;
	}
      },

      [&]( ThreadLocals &tls, int tid )  // Reduce thread locals
      {
	nActParticles += tls.nActParticles;
      }
    );
  }

  TIMER_STOP(&tm);
  TRCP(("CPU (determine active blocks) %.2f sec, %.0f points/sec)\n",
  (double)TIMER_DIFF_MS(&tm)/1000.,
  ((double)nActParticles)/(double)(TIMER_DIFF_MS(&tm)/1000.0)));	

  TRCP(("nActParticles = %.0f\n",(double)nActParticles));

  /////////////////////////////////////////////////////////////////////////////////
  //
  // 	Prepare MSBG block refinement map
  //
  /////////////////////////////////////////////////////////////////////////////////

  LongInt nActiveBlocks=0;
  {

    TIMER_START(&tm);
    
    std::unique_ptr<int[]> blockLevels( new int[msbg->nBlocks()] );
    for(int i=0;i<msbg->nBlocks();i++) 
    {
      int level=msbg->getNumLevels()-1;
      if( blockActive[i] )
      {
	level=0;
	nActiveBlocks++;
      }
      blockLevels[i] = level;
    }
    msbg->regularizeRefinementMap(blockLevels.get());
    msbg->setRefinementMap( blockLevels.get(),NULL,-1,NULL,false );

    TIMER_STOP(&tm);
    TRCP(("CPU (Set refinement map) %.2f sec, %.0f points/sec)\n",
    (double)TIMER_DIFF_MS(&tm)/1000.,
    ((double)nActParticles)/(double)(TIMER_DIFF_MS(&tm)/1000.0)));	
  }

  TRCP(("%s: nActiveBlocks=%d/%d %.2f%%\n",UT_FUNCNAME,
	nActiveBlocks,msbg->nBlocks(),
	100.*nActiveBlocks/(double)msbg->nBlocks()
	));

  FREEMEM( blockActive );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // 	Allocate and activate data channel blocks
  //
  /////////////////////////////////////////////////////////////////////////////////

  TRCP(("Allocating and initializing %d data blocks \n",nActiveBlocks));
  
  std::vector<int> activeBlocks;

  sg0->reset();
  sg0->prepareDataAccess( chan );
  sg0->setEmptyValue( renderDensFromFloat_<RenderDensity>(0.0f) );
  sg0->setFullValue( renderDensFromFloat_<RenderDensity>(1.0f) );
  activeBlocks.clear();

  TIMER_START(&tm);

  {

  using ThreadLocals = struct
  {
    std::vector<int> activeBlocks;
  };

  ThrRunParallel<ThreadLocals>( sg0->nBlocks(), 

    [&]( ThreadLocals &tls, int tid )  // Init thread locals
    {
      tls.activeBlocks.reserve(256);
    },

    [&](ThreadLocals &tls, int tid, int bid) 
    {
      BlockInfo *bi = msbg->getBlockInfo(bid);
      FLAG_RESET(bi->flags,BLK_EXISTS);
      if(bi->level>0) 
      {
        sg0->setEmptyBlock(bid);
        return;
      }

      tls.activeBlocks.push_back(bid);

      FLAG_SET(bi->flags,BLK_EXISTS);
      RenderDensity *data = sg0->getBlockDataPtr(bid,1,0);       
      for(int i=0;i<sg0->nVoxelsInBlock();i++) 
	data[i] = renderDensFromFloat_<RenderDensity>(1.0f) ;      
    },

    [&]( ThreadLocals &tls, int tid )  // Reduce thread locals
    {
      UT_VECTOR_APPEND( activeBlocks, tls.activeBlocks );
    },
    
    0,/*doSerial=*/true // Execute serially to reduce contention caused by memory page initialization  
	
  );

  }

  TIMER_STOP(&tm);
  TRCP(("CPU (allocate & initialize %d active blocks) %.2f sec, %.0f points/sec)\n",
  (int)activeBlocks.size(),(double)TIMER_DIFF_MS(&tm)/1000.,
  ((double)nActParticles)/(double)(TIMER_DIFF_MS(&tm)/1000.0)));	

  UT_ASSERT0( (size_t)nActiveBlocks == activeBlocks.size());

  /////////////////////////////////////////////////////////////////////////////////
  //
  // 	Splat the particles to the MSBG grid
  // 	
  /////////////////////////////////////////////////////////////////////////////////

  TRCP(("Splatting particles.\n"));

  TIMER_START(&tm);

  //float rParticleInvSq = 1.0f/(rParticle*rParticle);
  const float distSqMax = MT::sqf(rParticle + nbDist),
	      distSqMaxInv = 1.0f/distSqMax;

  std::vector<int> activeBlocksPerCol[8];

  for(int i=0; i<(int)activeBlocks.size();i++)
  {
    int bid = activeBlocks[i];
    Vec4i bpos = sg0->getBlockCoordsById(bid);
    int bx=vget_x(bpos),
	by=vget_y(bpos),
	bz=vget_z(bpos);
    int icol = getBlockColor8(bx,by,bz);
    UT_ASSERT0(icol>=0&&icol<8);
    activeBlocksPerCol[icol].push_back(bid);
  }

  for(int icol = 0; icol<8; icol++)
  {
    std::vector<int> *blockList = &activeBlocksPerCol[icol];
    ThrRunParallel<int>( blockList->size(), nullptr,
      [&](int &tls, int tid, int ibid) 
      {
	int bid = (*blockList)[ibid];
	UT_ASSERT0(sg0->isValueBlock(bid));
	
	Vec4i bpos = sg0->getBlockCoordsById(bid);
	int bx=vget_x(bpos),
	    by=vget_y(bpos),
	    bz=vget_z(bpos);

        UT_ASSERT0(getBlockColor8(bx,by,bz) == icol );

	Particle *p=NULL;
	for( ParticleIdx ixp = particlesPerBlock[bid]; 
	     ixp; 
	     ixp = p->idxNext )
	{
	  //TRCP(("bid=%d,ixp=%d\n",bid,ixp));
	  p = getParticle(ixp);

	  Vec4f pos0 = getPosition(p);

	  Vec4i ipos1 = max(truncate_to_int(ceil( pos0-rScan-0.5f )),0),
		ipos2 = min(truncate_to_int(floor( pos0+rScan-0.5f )),
			    sg0->v4iDomMax());

		int ix1=vget_x(ipos1), ix2=vget_x(ipos2),
	      iy1=vget_y(ipos1), iy2=vget_y(ipos2),
	      iz1=vget_z(ipos1), iz2=vget_z(ipos2);
	  Vec4f pshift = pos0-0.5f;	
	  float x0=vfget_x(pshift),
		y0=vfget_y(pshift),
		z0=vfget_z(pshift);

	  const int bsxLog2 = sg0->bsxLog2(),
	       bsx2Log2 = sg0->bsx2Log2(),
	       bsxMask = sg0->bsx()-1,
	       nx = sg0->nbx(),
	       nxy = sg0->nbxy();

	  for(int iz=iz1; iz<=iz2; iz++)
	  {
	    float dz = iz - z0;

	    int ibz = ( iz >> bsxLog2 ) * nxy,
		ivz = (iz & bsxMask) << bsx2Log2;

	    for(int iy=iy1; iy<=iy2; iy++)
	    {
	      float dy = iy - y0,
		    distSqZY = dz*dz + dy*dy;	    

	      int ibzy = ibz + (iy >> bsxLog2) * nx,
		  ivzy = ivz + ((iy & bsxMask) << bsxLog2);

	      for(int ix=ix1; ix<=ix2; ix++)
	      {
		float dx = ix - x0,
		      distSq = ( distSqZY + dx*dx );

		if(distSq > distSqMax) continue;  

		int ib = ibzy + ( ix >> bsxLog2 ),
		    iv = ivzy + ( ix & bsxMask );

		UT_ASSERT2(ib>=0&&ib<sg0->nBlocks());
		UT_ASSERT2(iv>=0&&iv<sg0->nVoxelsInBlock());

		SBG::Block<RenderDensity>* block = sg0->getBlock(ib); 	  	    
		
		if(!sg0->isValueBlock(block))
		{
		  goto loop_exit;
		}
		#if 0
		uint16_t val = 
		  renderDensFromFloat_<RenderDensity>(
		      std::max(distSqMax-distSq,0.0f)/distSqMax);
		block->_data[iv] = std::max(block->_data[iv],val); 
		#else
		uint16_t val = 
		  renderDensFromFloat_<RenderDensity>(distSq*distSqMaxInv);
		block->_data[iv] = std::min(block->_data[iv] , val);
		#endif
	      }
	    }
	  }
	  loop_exit:;
	}
      } 
    );
  }

  TIMER_STOP(&tm);
  TRCP(("CPU (splatting particles) %.2f sec, %.0f points/sec)\n",
  (double)TIMER_DIFF_MS(&tm)/1000.,
  ((double)nActParticles)/(double)(TIMER_DIFF_MS(&tm)/1000.0)));	

  FREEMEM( particlesPerBlock );
  FREEMEM( particles );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // 	Finalize distances and map to density [0,1] with iso-surfcae=0.5
  //
  /////////////////////////////////////////////////////////////////////////////////
  TRCP(("Finalizing splatting pass.\n"));

  TIMER_START(&tm);

  ThrRunParallel<int>( activeBlocks.size(), nullptr,
    [&](int &tls, int tid, int ibid) 
    {
      int bid = activeBlocks[ibid];
      UT_ASSERT0( sg0->isValueBlock(bid) );
      //BlockInfo *bi = msgDst->getBlockInfo(bid);
      //
      RenderDensity *data = sg0->getBlockDataPtr(bid);

      for(int vid=0;vid<sg0->nVoxelsInBlock();vid++)
      {
	float distSq = renderDensToFloat_(data[vid]);
	uint16_t uiVal = 0;
	if( distSq < 0.999f )
	{
	  distSq *= distSqMax;
	  float dist = sqrtf(distSq) - rParticle;
	  //f = std::max(std::min(dist,nbDist),-nbDist);
	  float f = 1.0f-MT_LINSTEPF(-rParticle,(rParticle+nbDist),dist);
	  //f = std::min(std::max(f,0.f),1.f);
	  uiVal = renderDensFromFloat_<RenderDensity>(f);
	}
	data[vid] = uiVal;
      }
    }
  );

  TIMER_STOP(&tm);

  TRCP(("CPU (finalizing splatting) %.2f sec\n",(double)TIMER_DIFF_MS(&tm)/1000.));

  /////////////////////////////////////////////////////////////////////////////////
  //
  //	Visualize intermediate result
  //
  /////////////////////////////////////////////////////////////////////////////////
  #if 0
  {
    static int pnl=-1;
    msbg->visualizeSlices( MiGetAuxPanel2(&pnl,"After splatting"),chan);
  }
  #endif

  #if 0
  {
    static int pnl = -1;
    render_scene( testCase, msbg, chan, MiGetAuxPanel2(&pnl,"scene_after_particle_splatting") ); 
  }
  #endif

  /////////////////////////////////////////////////////////////////////////////////
  //
  // 	Apply Mean Curvature Smoothing PDE
  //
  /////////////////////////////////////////////////////////////////////////////////
  int nIterPde = 5;

  TRCP(("Applying mean curvature smoothing PDE (%d iterations)\n",nIterPde));

  TIMER_START(&tm);

  msbg->applyChannelPdeFast<RenderDensity>( 
		       chan, CH_NULL, CH_NULL,
		       &activeBlocks,
		       -(PDE_MEAN_CURVATURE + OPT_8_COLOR_SCHEME ),
		      TRUE,
		      0, 0, 0, 
		      nIterPde, // number of iterations
		      1.0f,		     
		      testCase==2 ? 0.1 : 0.05, // time step		   
		      0, NULL, 0,0,0,0, 0, 0, 0,
		      0  // vis opt*/
		      );

  TIMER_STOP(&tm);

  TRCP(("CPU (Mean Curvature PDE) %.2f seconds.\n",(double)TIMER_DIFF_MS(&tm)/1000.));

  /////////////////////////////////////////////////////////////////////////////////
  //
  //	Visualize 2D slices and render 3D scene
  //
  /////////////////////////////////////////////////////////////////////////////////

  #if 0
  {
    static int pnl=-1;
    msbg->visualizeSlices( MiGetAuxPanel2(&pnl,
	  "surface after applyChannelPdeFast"),chan);
  }
  #endif

  TRCP(("Rendering final scene\n"));

  {
    static int pnl = -1;
    render_scene( testCase, msbg, chan, 
		  MiGetAuxPanel2(&pnl,"scene_after_PDE_smoothing") ); 
  }

  /////////////////////////////////////////////////////////////////////////////////
  //
  // 		Free resources 
  //
  /////////////////////////////////////////////////////////////////////////////////

  MSBG::MultiresSparseGrid::destroy(msbg);
  FREEMEM( particlesPerBlock );
  FREEMEM( blockActive );
  FREEMEM( particles );

  return 0;
}


int test_MSBG_0(int opcode, char *fpath, int bsx0, int sx, int sy, int sz)
{
#if 1
  return 1;
#else
  using namespace SBG;
  //int rcThis=0;
  //char *tok;

  //bool doWrite = opcode & 1;

//  bool doLevelSet = true;
  bool doLevelSet = false;

  //#define DENSITY_8_BIT

  std::vector<SparseGrid<float>*> sbgGrids;
  std::vector<char *> sbgNames;
  std::vector<unsigned> sbgOptions;

  int bsx = bsx0;
  MSBG::MultiresSparseGrid *msbg = 
    MSBG::MultiresSparseGrid::create( "TEST", sx,sy,sz, bsx,
				      -1, 0,-1, 
				      MSBG::OPT_SINGLE_LEVEL |
					MSBG::OPT_SINGLE_CHANNEL_FLOAT );
  UT_ASSERT0(msbg);

  int *blockRefinementMap = NULL;
  ALLOCARR_(blockRefinementMap,int,msbg->nBlocks());

 {

  #ifdef DENSITY_8_BIT
  typedef uint8_t DataType;
  int chan = MSBG::CH_GEN_UINT8;
  #else
  typedef float DataType;
  int chan = MSBG::CH_FLOAT_1;
  #endif
  msbg->resetChannel( chan );
  msbg->prepareDataAccess( chan );
  
  #ifdef DENSITY_8_BIT
  SparseGrid<uint8_t> *sbg = msbg->getUint8Channel(chan,0);
  #else
  SparseGrid<float> *sbg = msbg->getFloatChannel(chan,0);
  #endif

  sbg->setEmptyValue(0.0f);
  sbg->setFullValue(1.0f);

  #if 1
  {
    //
    // Fill with sphere surface (Enright)
    //

#if 1
    int fbmBaseShape = 1;
    Vec4f center_(0.5,0.5,.5,0.0f),
	  radius_(0.3,.3,.3,0);
    float resolution = msbg->sg0()->sxyzMax();
    float fbmSeed = 4711,
          fbmDepth_ = .1,
	  fbmFscale_ = 0.05,  
	  fbmMinscale = 0.5, // in fractions of voxel
	  fbmAmpl = .2, //0, //1, //1,
	  fbmAlpha = 1.5,
	  fbmOctLeadin = 1,
	  fbmThr = 0;
#endif
    USE(fbmSeed); USE(fbmFscale); USE(fbmMinscale); USE(fbmAlpha);
    USE(fbmOctLeadin); 
    

#if 0
    Vec4f sphereCenter(0.5,0.5,.5,0.0f);
    float sphereRadius = 0.4,
	  resolution = msbg->sg0()->sxyzMax();
    float fbmSeed = 4711,
	  fbmDepth = 4,
	  fbmFscale = 0.4*resolution,
	  fbmMinscale = 0.5,
	  fbmAmpl = 0, //1, //1,
	  fbmAlpha = 2,
	  fbmOctLeadin = 1,
	  fbmThr = 0;
#endif

#if 0
    Vec4f sphereCenter(0.5,0.5,.5,0.0f);
    float sphereRadius = 0.4,
	  resolution = msbg->sg0()->sxyzMax();
    float fbmSeed = 4711,
	  fbmDepth = 0.1*resolution,
	  fbmFscale = 0.1*resolution,
	  fbmMinscale = 0.5,
	  fbmAmpl = 1, //1, //1,
	  fbmAlpha = 1.5,
	  fbmOctLeadin = 1,
	  fbmThr = 1.2;
#endif

#if 0
    Vec4f sphereCenter(0.5,0.5,.5,0.0f);
    float sphereRadius = 0.2,
	  resolution = msbg->sg0()->sxyzMax();
    float fbmSeed = 4711,
	  fbmDepth = 0.4*resolution,
	  fbmFscale = 0.1*resolution,
	  fbmMinscale = 0.5,
	  fbmAmpl = 0, //1, //1,
	  fbmAlpha = 1.5,
	  fbmOctLeadin = 1,
	  fbmThr = .8;
#endif

#if 0
    Vec4f sphereCenter(0.5,0.5,.5,0.0f);
    float sphereRadius = 0.4,
	  resolution = msbg->sg0()->sxyzMax();
    float fbmSeed = 4711,
	  fbmDepth = 10, //0.2*resolution,
	  fbmFscale = 0.1*resolution,
	  fbmMinscale = 0.5,
	  fbmAmpl = 0, //1, //1,
	  fbmAlpha = 1.5,
	  fbmOctLeadin = 1,
	  fbmThr = .8;
#endif
    sphereCenter *= resolution;
    sphereRadius *= resolution;

    float sphereRadius0 = sphereRadius - fbmDepth,
	  sphereRadius1 = sphereRadius + fbmDepth;

    int sx = sbg->sx(),
	sy = sbg->sx(),
	sz = sbg->sx();
    USE(sx);USE(sy);USE(sz);

    UtTimer tm;

    TIMER_START(&tm);

    LongInt nActVoxels = 0;

    using ThreadLocals = struct
    {
      DataType *dataTmp;
      size_t nActVoxels;
    };

    ThrRunParallel<ThreadLocals>( sbg->nBlocks(), 
      [&]( ThreadLocals &tls, int tid )  // Initialize locals
      {
	ALLOCARR_(tls.dataTmp,DataType,sbg->nVoxelsInBlock());
	tls.nActVoxels = 0;
      },

      [&](ThreadLocals &tls, int tid, int bid) 
      {
        float phiMin=1e20,
	      phiMax=-1e20;
	// 
	// Perlin 'hypertexture'
	//	
	Vec4f center = resolution * center_,
	      radius = resolution * radius_;
	float depth = resolution * depth_;      

	auto baseShape = [&]( const Vec4f& pos ) -> float
	{
	  switch( fbmBaseShape )
	  {
	    case 1:  // Sphere
	    {
	      const float rSq = v3normSq( pos ),
			  r0Sq = sqf(radius[0]),
			  r1Sq = sqf(radius[0]-depth);

	      float f = rSq <= r1Sq ? 1.0f :
			rSq >= r0Sq ? 0.0f :
			    ( r0Sq - rSq ) * 1.0f/(r0Sq - r1Sq);
	      return f;
	    }

	    case 4:  // Ellipsoid (https://iquilezles.org/articles/distfunctions/)
	    {
	      const float k0 = sqrtf(v3normSq( pos * 1.0f/radius )),
		          k1 = sqrtf(v3normSq( pos * 1.0f/(radius*radius) ));
	      const float dist = k0*(k0-1.0f)/k1;
	      float f = 1.0f - MT_LINSTEPF( -depth, 0, dist ); 
	      return f;
	    }

	    case 7:  // Plane
	    {
	      Vec4f orig = center,
		    normal = radius * 1.0f/sqrtf(v3normSq(radius)),
		    tmp = normal * (pos-orig);

	      float dist = vfget_x(tmp)+vfget_y(tmp)+vfget_z(tmp);

	      float f = 1.0f-MT_LINSTEPF( 0, depth, dist ); 
	      return f;
	    }

	    default:
	    UT_ASSERT0(FALSE);
	    return 0;
	    break;
	  }
	};

	unsigned bflags = 0;

	Vec4i bpos = sbg->getBlockCoordsById(bid);
	Vec4f boxCenter = ((vint2float(bpos) + 0.5f) * (float)bsx);
	vzero_w(boxCenter);
	float dist = baseShape( boxCenter );

	blockRefinementMap[bid] = msbg->getNumLevels()-1;

	if( dist > (float)MT_SQRT3*bsx*0.5f )
	{
	  sbg->setEmptyBlock(bid);
	  return;
	}
	if( dist < (float)MT_SQRT3*bsx*0.5f )
	{
	  sbg->setFullBlock(bid);
	  return;
	}




	#if 1
	Vec4i bpos = sbg->getBlockCoordsById(bid);
	Vec4f boxCenter = ((vint2float(bpos) + 0.5f) * (float)bsx);
	vzero_w(boxCenter);
	float dist = sqrtf(v3normSq( boxCenter - sphereCenter ));

	blockRefinementMap[bid] = msbg->getNumLevels()-1;

	if( doLevelSet )
	{
	  if( dist > sphereRadius1 + (float)MT_SQRT3*bsx*0.5f )
	  {
	    sbg->setEmptyBlock(bid);
	    return;
	  }
	  if( dist < sphereRadius0 - (float)MT_SQRT3*bsx*0.5f )
	  {
	    sbg->setFullBlock(bid);
	    return;
	  }
	}
	else
	{
	  if( dist > sphereRadius1 + (float)MT_SQRT3*bsx*0.5f )
	  {
	    sbg->setEmptyBlock(bid);
	    return;
	  }
	  if( dist < sphereRadius0 - (float)MT_SQRT3*bsx*0.5f )
	  {
	    sbg->setFullBlock(bid);
	    return;
	  }
	}
	#endif

	{
	  float phiMin=1e20,
	        phiMax=-1e20;

	  SBG_FOR_EACH_VOXEL_IN_BLOCK(sbg,bid,x,y,z,vid)
	  {
	    Vec4f pos = Vec4f( x+0.5f, y+0.5f, z+0.5f, 0.0f);
	    float dist = sqrtf( v3normSq( pos - sphereCenter )) - sphereRadius;

	    float f = doLevelSet ? 	    
	      1.0f - MT_LINSTEPF( -fbmDepth, fbmDepth, dist ) : 
	      1.0f - MT_LINSTEPF( -fbmDepth, fbmDepth, dist );
	    
	    UT_ASSERT2(!( f<0.0f || f>1.0f || !std::isfinite(f)));
#if 1
	    float u = PNS::genFractalScalar3D( sx, sy, sz,
				pos[0]+sx,pos[1]+sy,pos[2]+sz,
				fbmSeed, fbmFscale, fbmMinscale, fbmOctLeadin, fbmAlpha, 0
				//713647, fscale, 2, 0, 2, 0 
				);
	    //TRCP(("f=%g u=%g\n",f,u));
#else
	    float u = 0;
#endif
	    f += fbmAmpl * u;
	    float phi = std::min(std::max(0.0f,f - fbmThr),1.0f);

	    //phi = 0.5; 
	    //phi = powf(x/(double)sbg->sx(),2.);

	    #ifdef DENSITY_8_BIT
	    tls.dataTmp[vid] = renderDensFromFloat_<DataType,false,false>(phi);
	    /*if(phi>0) 
	    {
	      TRCP(("%g -> %d\n",phi,tls.dataTmp[vid]));
	    }*/
	    #else
	    tls.dataTmp[vid] = phi;
	    #endif

	    phiMin = std::min(phiMin,phi);
	    phiMax = std::max(phiMax,phi);
	  }

	  if( phiMax < MT_NUM_EPS )
	  {
	    sbg->setEmptyBlock(bid);
	  }
	  else
	  {
	    blockRefinementMap[bid] = 0;
	    DataType *data = sbg->getBlockDataPtr(bid,1,0);
	    tls.nActVoxels += sbg->nVoxelsInBlock();
	    for(int i=0;i<sbg->nVoxelsInBlock();i++) data[i] = tls.dataTmp[i];
	  }
	}
      },

      [&]( ThreadLocals &tls, int tid )  // Reduce locals
      {
	FREEMEM(tls.dataTmp);
	nActVoxels += tls.nActVoxels;	
      }
    );

    msbg->regularizeRefinementMap( blockRefinementMap );  
    msbg->setRefinementMap( blockRefinementMap, NULL, -1, NULL,			
			       false, false );
    msbg->_blocksValue[0].clear();

    TIMER_STOP(&tm);
    TRCP(("Create: CPU=%.2f sec, %.0f voxels/sec. (%.0f/%.1f%%) voxels\n",
	  (double)TIMER_DIFF_MS(&tm)/1000.,
	(double)nActVoxels/(double)(TIMER_DIFF_MS(&tm)/1000.0),
	(double)nActVoxels,
	100.*nActVoxels/(double)sbg->nTotVirtVoxels()      
	));
  
    TRCP(("CPU %.2f sec  \n",UT_FUNCNAME,(double)TIMER_DIFF_MS(&tm)/1000.));

    sbg->showMemUsage(1);

    #if 1
    {
      static int pnl=-1;
      msbg->visualizeSlices( MiGetAuxPanel2(&pnl,"TEST"), chan );	
    }
    #endif
    #if 0
    {
      static int pnl=-1;
      msbg->visualizeSlices( MiGetAuxPanel2(&pnl,"TEST"),
				 0,IP_NEAREST,NULL, 0,0,0,0,
				 sbg );	
    }
    #endif
  }
  #endif	


  sbgGrids.push_back( sbg );
  sbgNames.push_back( doLevelSet ? (char*)"LevelSet" : (char*)"density" );
  sbgOptions.push_back( doLevelSet ? 1 : 0 );
  
  int rc = SBG::writeAsVDB(fpath,256,sbgGrids,sbgNames,sbgOptions);
  TRCP(("SBG::writeAsVDB('%s') -> %d\n",fpath,rc));

  }

rcCatch:
  MSBG::MultiresSparseGrid::destroy(msbg);
  FREEMEM(blockRefinementMap);
  return 0;
#endif
}

