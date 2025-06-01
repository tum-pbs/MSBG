/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#ifndef SBG_H
#define SBG_H

#include <functional>
#include <limits>
#include "vectorclass_util.h"
#include "VEC3.h"
//#include "msbgcell.h"
//#include "thread.h"
#include "blockpool.h"
#include "mtool2.h"
#include "fastmath.h"
#include "grid.h"

//#define TLB_USE_SKIPMAP

/*=========================================================================
 *
 *
 *  Sparse Block Grids
 *
 *
 * =======================================================================*/
/*=========================================================================
 *
 *
 *  Convenience macros
 *
 * =======================================================================*/

#define FOR_EACH_VOXEL_IN_HALO_BLOCK_SIMD( SIMDW, bsx, hsx, haloSize )\
for(int vz=0,vid=0;vz<bsx;vz++) \
for(int vy=0;vy<bsx;vy++) \
for(int vx=0, \
    hid=MT_GXYZ(hsx,hsx*hsx, haloSize, haloSize+vy, haloSize+vz);\
    vx<bsx; vx+=SIMDW,vid+=SIMDW,hid+=SIMDW) 

// A simple convenience macro for looping over all cells of a block 
// by the cells' (x,y,z) grid coordinates and the voxel-in-block index 'vid'
#define SBG_FOR_EACH_VOXEL_IN_BLOCK( sg_, bid_, x,y,z, vid )\
      int bx,by,bz, bsx=(sg_)->bsx(); \
      (sg_)->getBlockCoordsById(bid_,bx,by,bz); \
      UT_ASSERT2(sg_->blockCoordsInRange(bx,by,bz));\
      for(int vid=0,\
	      z=bz*bsx; z<(bz+1)*bsx; z++)\
      for(int y=by*bsx; y<(by+1)*bsx; y++)\
      for(int x=bx*bsx; x<(bx+1)*bsx; x++, vid++)


#define SBG_FOR_EACH_VOXEL_IN_BLOCK_GEN_2( sg_, bid_, xStrid, x,y,z, vid )\
      int x1,y1,z1,x2,y2,z2; \
      if( sg_->isDenseGrid())\
      { \
	x1=0; x2 = sg_->sx()-1;\
	y1=0; y2 = sg_->sy()-1;\
	z1=bid; z2=bid; \
      }\
      else\
      {\
	int bx,by,bz, bsx=(sg_)->bsx(); \
	(sg_)->getBlockCoordsById(bid_,bx,by,bz); \
	UT_ASSERT2(sg_->blockCoordsInRange(bx,by,bz));\
	x1=bx*bsx; x2=(bx+1)*bsx-1;\
	y1=by*bsx; y2=(by+1)*bsx-1;\
	z1=bz*bsx; z2=(bz+1)*bsx-1;\
      } \
      for(int z=z1; z<=z2; z++) \
      for(int y=y1; y<=y2; y++) \
      for(LongInt x=x1, \
	  vid = sg_->isDenseGrid() ? sg_->getGridIndex(x1,y,z) : \
	  		sg_->getVoxelInBlockIndex(0,y-y1,z-z1) ; \
	  x<=x2; x+=xStrid, vid+=xStrid )

#define SBG_FOR_EACH_VOXEL_IN_BLOCK_GEN( sg_, bid_, x,y,z, vid )\
  SBG_FOR_EACH_VOXEL_IN_BLOCK_GEN_2( sg_, bid_, 1, x,y,z, vid )

#ifdef __cplusplus

#define SBG_MAX_CHANNELS 32

#define SBG_MAX_BORDOFF	3

//#define SBG_IAC_STATS

namespace SBG
{
  class IpolAccessCache; 
  class IpolAccessCache2; 

  const float PCELL_INF_QUANTITY = (1e25f);
  
  //typedef Vector3Dim<float> Vec3Float; 
    
  typedef enum
  {
    CONSTVAL_NONE = 0,  // do not change
    CONSTVAL_EMPTY = 1,
    CONSTVAL_FULL = 2,
    CONSTVAL_INVALID = 3
  }
  ConstVal;

  enum 
  {
    IP_NEAREST = 1,	// XXX Do not change values
    IP_LINEAR = 2,	// as they are literally used in profile params
    IP_CUBIC = 3,
    IP_CUBIC_MONO = 4,
    IP_CUBIC_MONO_2 = 5,
    IP_BSPLINE_QUAD = 6,
    IP_BSPLINE_CUBIC = 7,
    IP_WENO4 = 8,
    IP_LINEAR_FADE = 9,
    IP_LINEAR_DIVERGENCE_FREE = 10
  };

  typedef enum
  {
    BC_DIRICHLET = 0,
    BC_NEUMANN = 1
  }
  BoundCondType;

  enum
  {
    ACC_READ = (1<<0),
    ACC_WRITE = (1<<1)
  };

  enum
  {
    OPT_IPCORNER = (1<<0),
    OPT_NO_DATA	= (1<<1),
    OPT_NO_BORDPAD = (1<<2),
    OPT_BLK_ALIGNED = (1<<3),
    OPT_IPGRIDALGN = (1<<4),
    OPT_IPCHKBLK = (1<<5),
    OPT_DENSE_GRID = (1<<6),
    OPT_IPEXTRAPOL = (1<<7),
    OPT_IPBC_NEUMAN = (1<<8),
    OPT_TRANSPARENCY_BUF = (1<<9),
    OPT_IP_VCOMP_X = (1<<10),
    OPT_IP_VCOMP_Y = (1<<11),
    OPT_IP_VCOMP_Z = (1<<12),
    OPT_IPBC_DIRICHLET = (1<<13),
    OPT_IPBC_DIRICHLET_OLD = (1<<14),
    OPT_IP_NO_BORDER_CHECK = (1<<15),
    OPT_TEST = (1<<16),
    OPT_IP_FADE_WEIGHTS = (1<<17),
    OPT_USE_OLD_BEHAVIOR = (1<<19),
    OPT_IPBC_PERIODIC = (1<<20),
    OPT_IPBC_FREESLIP = (1<<21),
    OPT_IPBC_GENERIC = (1<<22),
    OPT_ZERO_BORDER = (1<<23),
    OPT_SYNTHETIC_VELOCITY = (1<<24),
    OPT_FAST_DISCONTINUOUS_FINECOARSE = (1<<25),
    OPT_BLOCK_HALO = (1<<26),
    OPT_PROTECT_CONST_BLOCKS = (1<<27),
    OPT_IPDIVFREE = (1<<28),
    OPT_IP_VEC_GEN = (1<<29),
    OPT_PREP_DATA_ACCESS = (1<<30)
  };

  enum
  {
    SBLK_CONSTANT = (1<<0),
    SBLK_FULL = (1<<1),
    SBLK_EMPTY = (1<<2),
    SBLK_INVALID = (1<<3),
    SBLK_EXT_REFERENCE = (1<<4),
    SBLK_PROTECTED = (1<<5),
    SBLK_FIXED = (1<<6)
  };

  static int idxStencil7( int iDir, int iSide ) 
  {
    UT_ASSERT2(iDir>=0&&iDir<3&&iSide>=0&&iSide<2);
    return 1+iDir*2+iSide; 
  }

  inline int SDF_SIGN(float x) 
  {
    if (x>0.0f) return 1;
	        else return -1;
  }

  // For signed distance functions on grids 
  // (to be representable with 24 bit floats)
  const float SDF_INFINITY = 1e6f;

  //
  //
  //


  const uint16_t SDF_INFINITY_UI16 = ((1<<15)-1);
  const uint8_t SDF_INFINITY_UI8 = ((1<<7)-1);

  #define SDF_FP16_MANTBITS 15
  #define SDF_FP8_MANTBITS 7

  #define SDF_FP16_ORIGSIGN 15
  #define SDF_FP8_ORIGSIGN 7

  #define SDF_FP16_FACTOR_FLT2INT( narrowBandDist ) \
    ((1.0f/(narrowBandDist)) * ((1<<SDF_FP16_MANTBITS)-1))

  #define SDF_FP8_FACTOR_FLT2INT( narrowBandDist ) \
    ((1.0f/(narrowBandDist)) * ((1<<SDF_FP8_MANTBITS)-1))

  // 
  // renderDensToFloat 
  //
  
  #define RENDER_DENS_TO_FLOAT_UI16(x)  \
  {\
   (x) *= (1.0f/((1<<16)-1));\
  }
  #define RENDER_DENS_TO_FLOAT_UI8(x)  \
  {\
   (x) *= (1.0f/((1<<8)-1));\
  }

  template<typename T, bool doSqrtCompr=false> float renderDensToFloat_( T fi )
  {
    UT_ASSERT2( fi>=0 && fi<=std::numeric_limits<T>::max());
    float f = fi * (1.f/std::numeric_limits<T>::max());
    if constexpr ( doSqrtCompr )
    {
      f = f*f;
    }
    UT_ASSERT2(std::isfinite(f) && !((f)<0.0f||(f)>1.0f));
    return f;
  }

  template<typename T> Vec8f renderDensToFloat_simd8( T *pfi_simd8 );

  // uint8_t  
  template<> inline Vec8f renderDensToFloat_simd8( uint8_t *pfui_simd8 )
  {
    Vec8f val;
    vload_from_uint8_a( val, pfui_simd8 );
    RENDER_DENS_TO_FLOAT_UI8( val );
    return val;
  }

  // uint16_t  
  template<> inline Vec8f renderDensToFloat_simd8( uint16_t *pfui_simd8 )
  {
    Vec8f val;
    vload_from_uint16_a( val, pfui_simd8 );
    RENDER_DENS_TO_FLOAT_UI16( val );
    return val;
  }

  // float
  //template<> inline float renderDensToFloat_( float f ) { return f; }
  template<> inline Vec8f renderDensToFloat_simd8( float *pf_simd8 )
  {
    return Vec8f().load(pf_simd8);
  }


  // 
  // renderDensFromFloat 
  //

  template< typename T, bool doSR > 
  T renderDensRoundFromFloat( float f, 
      				float randValSR /*for doSR=true (stochastic rounding)*/ )
  {
    int fi;
    if constexpr ( doSR )  // Stochastic rounding
    { 
      float f_ = f,
	    fFloor_ = floor( f_ ),
	    fFrac_ = f_ - fFloor_;

      fi = randValSR < fFrac_ ? ((int)fFloor_) + 1 :
				((int)fFloor_);
    }
    else
    {
       fi = roundf( f );
    }

    UT_ASSERT2( fi>=0 && fi<=std::numeric_limits<T>::max());

    return fi;
  }

  template<typename T, bool doSqrtCompr=false, bool doSR=false > 
  T renderDensFromFloat_( float f, float randValSR=0 /* for stochastic rounding */ )
  {
    UT_ASSERT2(std::isfinite(f) && !((f)<0.0f||(f)>1.0f));
    if constexpr ( doSqrtCompr && sizeof(T)==1 )
    {
      float f_;
      FMA_FastSqrtf(&f_/*out*/,&f);
      f = f_;
    }
    return renderDensRoundFromFloat<T,doSR>(
	  	f * std::numeric_limits<T>::max(), randValSR);
  }

  template<typename T> void renderDensFromFloat_storeSimd8( Vec8f vf, T *pfi_simd8 );

  template<typename T> void renderDensFromFloat_storeSimd8_SR( Vec8f vf, T *pfi_simd8, Vec8ui &rnds );

  // uint8_t
  template<> inline void renderDensFromFloat_storeSimd8( Vec8f vf, uint8_t *pfui_simd8 ) 
  {
    Vec8ui tmp = (Vec8ui)round_to_int( vf * ((1<<8)-1)); 
    vstore_to_uint8_a( tmp, (uint8_t*)pfui_simd8 );
  }
  template<> inline void renderDensFromFloat_storeSimd8_SR( Vec8f vf, 
  /* with stochastic rounding */    		uint8_t *pfui_simd8, Vec8ui &rnds ) 
  {
    Vec8f vf_ = 255*vf,
	  vfFloor_ = floor( vf_ ),
	  vfFrac_ = vf_ - vfFloor_;

    Vec8f u = FMA_FastRand( &rnds );

    vf_ = select( u < vfFrac_, vfFloor_ + 1.0f,
			       vfFloor_ );

    Vec8ui vi = (Vec8ui)truncate_to_int( vf_ );
    vstore_to_uint8_a( vi, pfui_simd8 );
  }

  // uint16_t
  template<> inline void renderDensFromFloat_storeSimd8( Vec8f vf, uint16_t *pfui_simd8 ) 
  {
    #if 1
    Vec8ui tmp = (Vec8ui)round_to_int( vf * ((1<<16)-1)); // faster but depends on MXCSR 
    #else
    Vec8ui tmp = (Vec8ui)truncate_to_int(0.5f+( vf * ((1<<16)-1)));
    #endif
    vstore_to_uint16_a( tmp, (uint16_t*)pfui_simd8 );
  }
  template<> inline void renderDensFromFloat_storeSimd8_SR( Vec8f vf, 
      uint16_t *pfui_simd8, Vec8ui &rnds ) {}

  // float
  //template<> inline float renderDensFromFloat_( float f ) { return f; }
  template<> inline void renderDensFromFloat_storeSimd8( Vec8f vf, float *pf_simd8 ) 
  {
#if 0
    vf.store_a(pf_simd8);
#else
    vstream((float*)pf_simd8,vf);
#endif
  }
  template<> inline void renderDensFromFloat_storeSimd8_SR( Vec8f vf, 
      float *pfui_simd8, Vec8ui &rnds ) {}

  //
  //
  //

  typedef struct
  {
    int bcType[3][2];   // X:left-right, Y:bottom-top, Z:front-back
    // Supported types: OPT_IPBC_DIRICHLET, 
    // 			OPT_IPBC_NEUMAN,
    // 			OPT_IPBC_PERIODIC,
    // 			OPT_IPBC_FREESLIP (for vector valued grids)  
  }
  DomainBC;

  inline 
  void setDomainBC( DomainBC *dbc,
      	       int bcLeft, int bcRight, 
	       int bcBottom, int bcTop,
	       int bcFront, int bcBack )
  {
    dbc->bcType[0][0] = bcLeft; dbc->bcType[0][1] = bcRight;
    dbc->bcType[1][0] = bcBottom; dbc->bcType[1][1] = bcTop;
    dbc->bcType[2][0] = bcFront; dbc->bcType[2][1] = bcBack;
  }

/*=========================================================================
 *
 *  Utility
 *  
 * =======================================================================*/
inline void fadeIpolWeights( float& u, float& v, float& w)
{
  //http://www.iquilezles.org/www/articles/texture/texture.htm
  #if 1
  u = u*u*u*(u*(u*6.0f-15.0f)+10.0f);
  v = v*v*v*(v*(v*6.0f-15.0f)+10.0f);
  w = w*w*w*(w*(w*6.0f-15.0f)+10.0f);
  #else
  u =  (u*u * (3.0f - 2.0f*u));
  v =  (v*v * (3.0f - 2.0f*v));
  w =  (w*w * (3.0f - 2.0f*w));
  #endif
  MT_CLAMP(u,0,1);
  MT_CLAMP(v,0,1);
  MT_CLAMP(w,0,1);  
}

static inline void setValueComponent( Vec3Float& val, int iCmp, float a )
{
  val[iCmp] = a;
}

static inline void setValueComponent( Vec3Uint16& val, int iCmp, float a )
{
  val[iCmp] = a;
}

static inline void setValueComponent( float& val, int iCmp, float a )
{
  val = a;
}

static inline void setValueComponent( double& val, int iCmp, float a )
{
  val = a;
}

static inline void setValueComponent( uint16_t& val, int iCmp, float a )
{
  val = a;
}


static inline void setValueComponent( int16_t& val, int iCmp, float a )
{
  val = a;
}

static inline void setValueComponent( uint8_t& val, int iCmp, float a )
{
  val = a;
}


/*=========================================================================
 *
 *  Block
 *  
 * =======================================================================*/
template<typename Data_T> using Block = BlockPool::Block<Data_T>;
template<typename Data_T>  Data_T zeroValue( void ); 
template<typename Data_T>  Data_T fullValue( void ); 
template<typename T> float diffValues( const T& val1, const T& val2 );

/*=========================================================================
 *
 *  SparseGrid
 *  
 * =======================================================================*/


template <typename Data_T>
class SparseGrid
{
  public:

  // Create 
  static SparseGrid<Data_T>* create(
	  // Name tag for monitoring/tracing convenience
          const char *name,
	  // Size of grid in each dimension in voxels
	  int sx, int sy, int sz,
	  // Block size in voxels (power-of-two)
	  int  blockSize = 16,
	  // Maximum size to which the grid is allowed to grow (0=unlimited)
	  int  maxSizeMB = 0,
	  // Initially allocated memory in MB (default: ~2% of virtual size)
	  int  initialSizeMB = -1,
	  // Misc. options: OPT_NO_DATA 
	  unsigned options = 0,
	  // optional: specify entry size in bytes
	  size_t entrySize = 0
	  );

  // Clone
  SparseGrid<Data_T>* clone( const char *name=NULL, 
      			     int doCopyData=FALSE,
			     SparseGrid<Data_T>*sgDest=NULL )
  {
    TRC(("%s '%s' -> '%s'\n",UT_FUNCNAME,_name,name));
    UT_ASSERT0(!_isDenseGrid); 	// !OPT_DENSE_GRID

    SparseGrid<Data_T>* sg=NULL;

    if(!sgDest)
    {
      sg = create( name?name:_name, 
		   _sx, _sy, _sz, 
		   _bsx, _maxSizeMB, _initialSizeMB,
		   _optionsOrig
		   );
    }
    else
    {
      sg = sgDest;
      sg->reset();
    }

    sg->setFullValue( _fullValue );
    sg->setEmptyValue( _emptyValue );
    sg->setInvalidValue( _invalidValue );

    if(doCopyData)
    {      
      sg->prepareDataAccess();
    ThrRunParallel<int> ( _nBlocks, nullptr,        
      [&](int &tls, int tid, int bid)
      {
	Data_T *blockData = getBlockDataPtr(bid);
	if(!blockData) return;  // skip empty blocks	
	if(isFullBlock(bid))
	{
	  sg->setFullBlock(bid);
	  return;
	}
	if(isEmptyBlock(bid))
	{
	  sg->setEmptyBlock(bid);
	  return;
	}

	Data_T *blockDataDst = sg->getBlockDataPtr(bid,1,0);
	for(int vid=0;vid<_nVoxelsInBlock;vid++)
	{
	  blockDataDst[vid] = blockData[vid];
	}
      } );
    }
    return sg;
  }

  // Copy Data
  void copyDataTo( SparseGrid<Data_T>* sgDst )
  {
    TRC(("%s '%s' -> '%s'\n",UT_FUNCNAME,getName(),sgDst->getName()));
    UT_ASSERT0(!(isDenseGrid() || sgDst->isDenseGrid())); 
    UT_ASSERT0( gridDimEqual(this,sgDst) );
    UT_ASSERT0( sgDst->bsx() == _bsx );
    clone( "N/A", 1, sgDst );
  }

  // 
  void setNullBlocksToEmptyVal( void )
  {
    TRC(("%s '%s'\n",UT_FUNCNAME,_name));
    #if 1
    ThrRunParallel<int> ( _nBlocks, nullptr,        
      [&](int &tls, int tid, int bid)
      {
        if(!_blockmap[bid]) setEmptyBlock(bid);
      } );
    #else
    #pragma omp parallel for if(_nBlocks>5000)
    for(int bid=0;bid<_nBlocks;bid++)
      if(!_blockmap[bid]) setEmptyBlock(bid);
    #endif

  }

  void setBlocksToEmptyVal( void )
  {
    ThrRunParallel<int> ( _nBlocks, nullptr,        
      [&](int &tls, int tid, int bid)
      {
        setEmptyBlock(bid);
      } );
  }

  static void setParentPtr( void *block, void *ptr )
  {
    Block<Data_T> *block_ = (Block<Data_T>*)block;
    block_->uheader.parent = ptr;
  }

  static void *getParentPtr( void *block )
  {
    Block<Data_T> *block_ = (Block<Data_T>*)block;
    return block_->uheader.parent;
  }

  //
  int writeAsVDBGrid( int fd /* file descriptor */ );

  // 
  int loadDataFromOpenVDB( char *fname );

  // 
  int saveToOpenVDB( char *fname );

  // Destroy 
  static void destroy(SparseGrid<Data_T>*& sg);

  // Check if voxel coordinates in grid range 
  int inRange(const Vec4f& pos) const
  {
    return !( _mm_movemask_ps( pos < _v4fDomMin || pos > _v4fDomMax )  );
  }

  int inRange(const Vec4i& ipos) const
  {
//    return vall( ipos >= _mm_setzero_si128() && ipos <= _v4iDomMax );
    return vall( ipos >= vzero() && ipos <= _v4iDomMax );
  }


  // Check if voxel coordinates in grid range 
  int inRange(int i, int j,int k) const
  {
    return (i>=0&&i<_sx) && (j>=0&&j<_sy) && (k>=0&&k<_sz);
  }

  // Check if voxel coordinates in grid range 
  int inRange(int bi, int bj, int bk,  // IN: block coords
      		     int vi, int vj, int vk  // IN: voxel-in-block coords
			) const
  {
    int i,j,k;
    getGridCoords( bi,bj,bk, vi,vj,vk,
		   i,j,k );
    return inRange(i,j,k);
  }

  // Check block coordinates 
  int blockCoordsInRange(int bi, int bj, int bk) const
  {
    return (bk>=0&&bk<_nz) && (bj>=0&&bj<_ny) && (bi>=0&&bi<_nx);
  }

  int blockCoordsInRange(const Vec4i& bpos) const
  {
    return vall( bpos >= vzero() && bpos <= _v4iBlkCoordsMax );
  }

  // Check voxel-in-block coords 
  int voxelInBlockRange(int vi, int vj, int vk) const
  {
    return (vi>=0&&vi<_bsx) && (vj>=0&&vj<_bsx) && (vk>=0&&vk<_bsx);
  }

  //
  bool blockIntersectsBox( int bx, int by, int bz,
      			   const Vec4f& boxMin, // XXX: assumed boxMin[3] = -1e20 
			   const Vec4f& boxMax  // boxMax[3] = 1e20
			   ) const
  {
    Vec4f bpos0 = Vec4f(bx,by,bz,0)*bsx(),
	  bpos1 = bpos0 + bsx() + Vec4f(0.f,0.f,0.f,1e20f); 

    return vall( bpos1 > boxMin && bpos0 < boxMax ); 
  }

  // Get grid coords from block and voxel-in-block coords
  void getGridCoords(int bi, int bj, int bk,  // IN: block coords
      		     int vi, int vj, int vk,  // IN: voxel-in-block coords
      		     int &i, int &j, int &k  // OUT: grid coords
			) const
  {
    //UT_ASSERT2( blockCoordsInRange(bi,bj,bk));
    UT_ASSERT2( voxelInBlockRange(vi,vj,vk));
    i = bi*_bsx+vi/*-_maxBordOff*/;
    j = bj*_bsx+vj/*-_maxBordOff*/;
    k = bk*_bsx+vk/*-_maxBordOff*/;
    //UT_ASSERT2(i<_nx*_bsx && j<_ny*_bsx && k<_nz*_bsx);
  }

  void getGridCoordsGen(int bi, int bj, int bk,  // IN: block coords
      		        int vi, int vj, int vk,  // IN: voxel-in-block coords
      		        int &i, int &j, int &k  // OUT: grid coords
			) const
  {
    if(_isDenseGrid) 
    {
      UT_ASSERT2(_maxBordOff==0);
      i = vi; j = vj; k = vk;
    }
    else getGridCoords( bi,bj,bk, vi, vj, vk, i, j, k);
  }
    

  Vec4i getGridCoords( int bx, int by, int bz, 
      		       int vx, int vy, int vz )
  {
    int x,y,z;
    getGridCoords( bx, by, bz, vx, vy, vz,
		   x, y, z );
    return Vec4i(x,y,z,0);
  }

  Vec4i getGridCoords( int bid, int vid )
  {
    int bx,by,bz,vx,vy,vz,x,y,z;
    getBlockCoordsById(bid, bx,by,bz);
    getVoxelInBlockCoordsById(vid, vx,vy,vz);
    getGridCoords( bx, by, bz, vx, vy, vz,
		   x, y, z );
    return Vec4i(x,y,z,0);
  }


  // Clip coords to grid boundary
  void clipBox( 
      int x1, int y1, int z1,
      int x2, int y2, int z2,      
      int &x, int &y, int &z   // OUT
      ) const
  {
    if(unlikely((x)>(x2))) x = x2; 
    if(unlikely((x)<(x1))) x = x1; 
    if(unlikely((y)>(y2))) y = y2; 
    if(unlikely((y)<(y1))) y = y1; 
    if(unlikely((z)>(z2))) z = z2;
    if(unlikely((z)<(z1))) z = z1; 
  }

  void clipBlockCoords( Vec4i& bpos ) const
  {
    bpos = max( bpos, _mm_setzero_si128() );
    bpos = min( bpos, _v4iBlkCoordsMax );
  }

  void clipGridCoords( Vec4i& ip ) const
  {
    ip = max( ip, _mm_setzero_si128() );
    ip = min( ip, _v4iDomMax );
  }

  void clipGridCoords( int &x, int &y, int &z ) const
  {
    clipBox(0,0,0, _sx-1,_sy-1,_sz-1,
	    x,y,z);
  }

  int clipGridCoords( int &x, int &y, int &z,
      		      int &clipped ) const
  {
    clipped=0;
    if(unlikely((x)>=(_sx))) { x = _sx-1; clipped=1; } 
    if(unlikely((x)<(0))) { x = 0; clipped=1; } 
    if(unlikely((y)>=(_sy))) { y = _sy-1; clipped=1; } 
    if(unlikely((y)<(0))) { y = 0; clipped=1; } 
    if(unlikely((z)>=(_sz))) { z = _sz-1; clipped=1; }
    if(unlikely((z)<(0))) { z = 0; clipped=1; } 
    return clipped;
  }
 
  int clipGridCoords( int &x, int &y, int &z,
      		      int &clipped,
		      int &iDir, int &iFace		      
		      ) const
  {
    clipped=0; iDir=-1,iFace=-1;

    int clippedZ = 0;
    if(unlikely(z<0)) {    z=0;     clippedZ=1; iDir=2; iFace=0; } 
    if(unlikely(z>=_sz)) { z=_sz-1; clippedZ=1; iDir=2; iFace=1; }

    int clippedY = clippedZ;
    if(unlikely(y<0)) {     y=0;     clippedY=1; iDir=1; iFace=0; } 
    if(unlikely(y>=_sy)) {  y=_sy-1; clippedY=1; iDir=1; iFace=1; } 

    clipped = clippedY;
    if(unlikely(x<0)) {    x = 0;     clipped=1; iDir=0; iFace=0; } 
    if(unlikely(x>=_sx)) { x = _sx-1; clipped=1; iDir=0; iFace=1; } 

    return clipped;
  }


  // Apply domain boundary conditions 
  Data_T applyDomainBC( int ix, int iy, int iz,  // Original coords
			int ix2, int iy2, int iz2,   // Clipped coords
			int iDir, int iFace ) const
  {
    Data_T val;
    int bcTyp = _domainBC.bcType[iDir][iFace];
    switch( bcTyp )
    {
      case OPT_IPBC_DIRICHLET: 
	val = _emptyValue; 
	break;
      case OPT_IPBC_NEUMAN: 
	val = getValue(ix2,iy2,iz2); 
	break;
      case OPT_IPBC_PERIODIC: 
	ix2 = ix<0 ? _sx - ((_sx-ix)%_sx) : ix%_sx;		        
	iy2 = iy<0 ? _sy - ((_sy-iy)%_sy) : iy%_sy;		        
	iz2 = iz<0 ? _sz - ((_sz-iz)%_sz) : iz%_sz;		        
	val = getValue(ix2,iy2,iz2);
	break;
      case OPT_IPBC_FREESLIP:
	val = getValue(ix2,iy2,iz2); 
	setValueComponent(val,iDir,0);
	break;
      default: UT_ASSERT0(FALSE);
    }
    return val;
  }

  Data_T getValueWithBC( int x, int y, int z, int& isClipped ) const
  {
    Data_T val;
    int iDir,iFace,
	x2=x,y2=y,z2=z;
    clipGridCoords( x2, y2, z2, 
	  	    isClipped, iDir, iFace );
    if(unlikely(isClipped))
    {
      val = applyDomainBC( x,y,z, x2,y2,z2, iDir, iFace );
    }
    else val = getValue(x,y,z);
    return val;
  }

  // Get block pointer by block coords
  Block<Data_T>* getBlock( int bi, int bj, int bk ) const
  {
    return blockCoordsInRange(bi,bj,bk) ? 
      		_blockmap[ getBlockIndex(bi,bj,bk) ] : 
		NULL;
  }

  // Get Access to block data pointer
  Data_T *getDataPtrGen_( int bid=0, 
      		      int doAlloc=FALSE, int doZeroDataOnInit=TRUE, 
		      int doAllocDenseData=FALSE )
  { 
    if( _isDenseGrid )
    {
      if(doAlloc && doZeroDataOnInit)
      {
	UT_ASSERT0(FALSE);  // multithreading
      }
      if(doAllocDenseData) getDenseData();
      UT_ASSERT2(_denseData);
      return _denseData;
    }
    else return getBlockDataPtr(bid,doAlloc,doZeroDataOnInit);
  }

  Data_T getValueDense_( int x, int y, int z ) const
  {
    UT_ASSERT2(inRange(x,y,z));
    applyBorderOffset(x,y,z);
    return _denseData[ G3D_INDEX0(_sxAct,_sxyAct, x,y,z) ];
  }

#if 0
  void setValueGen_( int x, int y, int z, Data_T val )
  {
    if(_isDenseGrid)
    {
      applyBorderOffset(x,y,z);
      UT_ASSERT2( G3D_IN_BOX(0,0,0,_sxAct-1,_syAct-1,_szAct-1, x,y,z) );
      _denseData[ G3D_INDEX0(_sxAct,_sxyAct, x,y,z) ] = val;
    }
    else
    {
      setValue(x,y,z,val);
    }
  }
#endif

  Data_T getValueGen_( int x, int y, int z )
  {
    if(_isDenseGrid)
    {
      applyBorderOffset(x,y,z);
      return G3D_IN_BOX(0,0,0,_sxAct-1,_syAct-1,_szAct-1, x,y,z) ? 
		_denseData[ G3D_INDEX0(_sxAct,_sxyAct, x,y,z) ] :
		_emptyValue;
    }
    else
    {
      return getValue(x,y,z);
    }
  }

  Data_T getValueGen_( Vec4i ipos ) 
  {
    return getValueGen_( vget_x(ipos), vget_y(ipos), vget_z(ipos) );
  }

  void checkExtBlockPtr( Block<Data_T> *block, int bid ) const
  {
    if(block && !isConstBlock(block))
    {
      SparseGrid<Data_T> *sgParent = (SparseGrid<Data_T>*)getParentPtr(block);
      UT_ASSERT0( sgParent->_dataGenCount == _dataGenCount );
      /*UT_ASSERT0( sgParent->_blockmap );
      UT_ASSERT0( sgParent->_blockmap[bid] ==  block );*/
    }
  }

  //
  Data_T *getBlockDataPtr( int bid, 
      			   int doAlloc=FALSE, 

			   int doZeroDataOnInit=FALSE
			   //int doZeroDataOnInit=TRUE
			   
			   ) 
  {
#if 0
    Block<Data_T>* block_ = (Block<Data_T>*)_blockPool->allocBlock(bid);
    return &block_->_data[0];
#endif
    
    /*int doReadWrite = doAlloc & READWRITE; xxx 
      *** PRIO: UtMemProtect() for constant blocks (level 0)*/

    #ifdef UT_ASSERT_LEVEL_2
    checkConstBlocks(0);
    #endif

    Block<Data_T>* block = getBlock(bid);
    
    #ifdef UT_ASSERT_LEVEL_2
    if(block)
    {
      BlockPool::checkEyecatcher(block); 
    }
    #endif

    if( doAlloc && 
	( (!block)  /* "alloc-on-write" TODO: copy-on-write => alloc and copy content from old block*/
	  || (block && isReadOnlyBlock( block ))))
    {
      if(doAlloc)
      {
	if(block)
	{
	  // Do not allow writes into foreign blocks          
	  //UT_ASSERT0( getParentPtr((void*)block) == this );
	  //? Do not allow writes into const blocks         
          #if 0
	  if(isConstBlock(block)&&!isEmptyBlock(block))
	  {
	    TRCWARN(("%s '%s' doAlloc=1 overwriting existing ?\n",
		  UT_FUNCNAME,_name));
	  }
	  #endif
	  
	  //UT_ASSERT0( !isConstBlock(block) );
	}
      }

      /*if(doZeroDataOnInit)
      {
	TRCWARN(("getBlockDataPtr() doZeroDataOnInit=TRUE\n"));
      }*/

      block = allocBlock(bid, doZeroDataOnInit );
    }

    #ifdef UT_ASSERT_LEVEL_2
    checkExtBlockPtr(block,bid);
    #endif

    return block ? &block->_data[0] : NULL;
  }

  // Get Access to block data pointer
  Data_T *getBlockDataPtr0( int bid ) const
  { 
    UT_ASSERT2(bid>=0&&bid<_nBlocks);
    #ifdef UT_ASSERT_LEVEL_2
    Block<Data_T>* block = getBlock(bid);
    checkExtBlockPtr(block,bid);
    #endif
    return &_blockmap[bid]->_data[0];
  }

  Data_T *getBlockDataPtr( Block<Data_T> *block ) const
  { 
    return &block->_data[0];
  }

  Data_T* getBlockDataPtrAtomicZeroInit( int bid, volatile LONG *blockLock,
				   LongInt *nLockCol )
  { 
    // Lock
    while(!((InterlockedCompareExchange((volatile LONG*)(blockLock), 1, 0 ) == 0)))
    {
      if(nLockCol) *nLockCol += 1;
      Sleep(0);	  
    }
    Data_T* block = getBlockDataPtr( bid, TRUE, TRUE );
    InterlockedExchange((volatile LONG*)(blockLock), 0 ); // unlock
    return block;
  }

  void setExternalBlockDataPtr( int bid, Block<Data_T> *block ) 
  { 
#if 0  // fast monotonic blockpool allocation does not support "freeBlock"
    Block<Data_T> *blockPrev = getBlock(bid);
    if( blockPrev && getParentPtr(blockPrev) == this )
    {
      freeBlock_(bid);     
    }
#endif
    _blockmap[bid] = block;
    block->uheader.flags |= SBLK_EXT_REFERENCE;
  }
  
  //
  //
  //
  #if 0
  void prepareDomBorderGhostCells( 
         std::function< BoundCondType( int iDir, int iSide ) > getBC 
		         // rturns BC_DIRICHLET or BC_NEUMANN	 
	 )
  {
    if(!_isDenseGrid) return;

    for(int iDir=0;iDir<3;iDir++)
    {
      for(int iSide=0;iSide<2;iSide++)
      {
        BoundCondType bcType = getBC(iDir,iSide);
      }
    }
  }
  #endif

  // Misc. getters/setters
  LongInt sx() const { return _sx; }
  LongInt sy() const { return _sy; }
  LongInt sz() const { return _sz; }
  LongInt sxAct() const { return _sxAct; }
  LongInt syAct() const { return _syAct; }
  LongInt szAct() const { return _szAct; }
  LongInt bsxLog2() const { return _bsxLog2; }
  LongInt bsx2Log2() const { return _bsx2Log2; }
  LongInt sxyzMax() const { return _sxyzMax; }
  LongInt sxyzMin() const { return _sxyzMin; }
  const Vec4i &v4iGridStridesBlk() const { return _v4iGridStridesBlk; }

  LongInt nTotActVoxels() const { return _nTotActVoxels; }
  LongInt nTotVirtVoxels() const { return _nTotVirtVoxels; }
  LongInt nBlocks() const { return _nBlocks; }
  LongInt nVoxelsInBlock() const { return _nVoxelsInBlock; }
  int isDenseGrid() const { return _isDenseGrid; }

  void setNbDist( float nb ) { _nbDist = nb; }
  float getNbDist( ) const { return _nbDist; }
  
  // Misc.
  void lock(void) { UtMutexLock(&_lock); }
  void unlock(void) { UtMutexUnlock(&_lock); }

  Data_T *getDenseData() 
  {
    if(!_isDenseGrid) return NULL;
    lock();
    if(!_denseData)
    {
      size_t szNeeded = _nTotActVoxels*sizeof(Data_T);
      szNeeded += CPU_SIMD_WIDTH;  // for safety
      ALLOCMEM_ALIGNED_(_denseData,Data_T,szNeeded,CPU_CACHELINE_SZ);
      TRC3(("Allocated dense grid data '%s' %.0f\n",_name,(double)szNeeded));
      memset((void*)_denseData,0,szNeeded);
    }
    unlock();
    return _denseData;
  }
  
  int initialSizeMB() { return _initialSizeMB; }
  int maxSizeMB() { return _maxSizeMB; }
  unsigned optionsOrig() { return _optionsOrig; }

  int nbx() const { return _nx; }
  int nby() const { return _ny; }
  int nbz() const { return _nz; }
  int nbxy() const { return _nxy; }

  size_t getEntrySize() const { return _entrySize; }
  const int *getNeighVoxOffs() const { return _neighVoxOffs; }

  Vec3Int gridRes() { return _gridRes; }
  Vec3Int gridResBlk() { return _gridResBlk; }

  int xStride() const { return 1; };
  int yStride() const { return _isDenseGrid ? _sxAct : _bsx; };
  int zStride() const { return _isDenseGrid ? _sxyAct: _bsx2; };

  Vec4i v4iDomMax() const { return _v4iDomMax; };

  Vec4f v4fDomMin() const { return _v4fDomMin; };
  Vec4f v4fDomMax() const { return _v4fDomMax; };

  Vec4f v4fDomMin2() const { return _v4fDomMin2; };
  Vec4f v4fDomMax2() const { return _v4fDomMax2; };

  int bsx() const { return _bsx; }

  const char *getName() const { return _name; }

  void setDomainBC( DomainBC *dbc )
  {
    TRC3(("%s '%s' => %d:%d %d:%d %d:%d\n",UT_FUNCNAME,_name,
      dbc->bcType[0][0],dbc->bcType[0][1], 
      dbc->bcType[1][0],dbc->bcType[1][1], 
      dbc->bcType[2][0],dbc->bcType[2][1]));

    for(int i=0;i<3;i++) 
      for(int j=0;j<2;j++) _domainBC.bcType[i][j] = dbc->bcType[i][j];
  }

  int maxBordOffs() const { return _maxBordOff; }

#if 0
  void setWorldSpaceInfo( double *pos, double dx )
  {
    for(int k=0;k<3;k++) _wsInfo.pos[k] = pos[k];
    _wsInfo.dx = dx;
    _wsInfo.dx_r = 1./dx;
  }

  void worldToGrid( double *w, // IN: world coords
      		    double *g  // OUT: grid coords
		    )
  {
      g[0] = ((w)[0]-(_wsInfo.pos)[0])*(_wsInfo.dx_r);
      g[2] = ((w)[1]-(_wsInfo.pos)[1])*(_wsInfo.dx_r);  /* SWAP Y/Z axis !*/
      g[1] = ((w)[2]-(_wsInfo.pos)[2])*(_wsInfo.dx_r);
  }

  void gridToWorld( double *g, // IN: world coords
      		    double *w  // OUT: grid coords
		    )
  {\
      w[0] = (g)[0]*(_wsInfo.dx)+(_wsInfo.pos)[0]; \
      w[1] = (g)[2]*(_wsInfo.dx)+(_wsInfo.pos)[1];  /* SWAP Y/Z axis !*/\
      w[2] = (g)[1]*(_wsInfo.dx)+(_wsInfo.pos)[2]; \
  }
#endif

  // Get value by grid coordinates
  Data_T getValue( int i, int j, int k ) const;
  Data_T getValue0( int i, int j, int k ) const;

  Data_T getValue( Vec3Int gridCoords ) const
  {
    return getValue(gridCoords.x, gridCoords.y, gridCoords.z);
  }

  Data_T getValue( Vec4i ipos ) const
  {
    return getValue(vget_x(ipos),vget_y(ipos),vget_z(ipos));
  }

  Data_T getValue( Vec4f pos ) const
  {
    Vec4i ipos = truncate_to_int(floor(pos));
    ipos = max(ipos,_mm_setzero_si128());
    ipos = min(ipos,_v4iDomMax);
    return getValue(vget_x(ipos),vget_y(ipos),vget_z(ipos));
  }

  Block<Data_T> *getBlock( Vec4f pos ) const
  {
    Vec4i ipos = truncate_to_int(floor(pos));
    ipos = max(ipos,_mm_setzero_si128());
    ipos = min(ipos,_v4iDomMax);
    int bi, bj, bk; 
    getBlockCoords(vget_x(ipos), vget_y(ipos), vget_z(ipos), 
	           bi, bj, bk);  
    int id = getBlockIndex(bi, bj, bk);
    Block<Data_T>* block = _blockmap[id];
    return block;
  }

  Data_T* getValuePtr( int i, int j, int k ) const;

  Data_T* getValuePtrGen( int x, int y, int z ) const
  {
    if(_isDenseGrid)
    {
      applyBorderOffset(x,y,z);
      return G3D_IN_BOX(0,0,0,_sxAct-1,_syAct-1,_szAct-1, x,y,z) ? 
		_denseData + G3D_INDEX0(_sxAct,_sxyAct, x,y,z)  :
		NULL;
    }
    else return getValuePtr(x,y,z);  
  }

  Data_T *getValueInBlockPtr0(Block<Data_T>* block, int vid) const
  { 
    return (Data_T*)((uint8_t*)(&block->_data[0])+vid*_entrySize);
  }

  Data_T *getValueInBlockPtr0(int bid, int vid) const
  { 
    return getValueInBlockPtr0(getBlock(bid),vid);  
  }

  // Get value in block
  Data_T getValueInBlock(int bid, int vid) const
  { 
    Block<Data_T>* block = getBlock(bid);
    return block ? *getValueInBlockPtr0(block,vid) : _emptyValue;
  }

  Data_T getValueInBlockGen(int bid, int vid) const
  { 
    if(_isDenseGrid)
    {
      UT_ASSERT2( vid>=0 && vid<_nTotActVoxels );
      return _denseData[ vid ];
    }
    else return getValueInBlock( bid, vid );
  }

  Data_T getValueInBlock0(int bid, int vid) const
  { 
    Block<Data_T>* block = getBlock(bid);
    return *getValueInBlockPtr0(block,vid);
  }

  Data_T *getValueInBlockPtr(int bid, int vid) const
  { 
    Block<Data_T>* block = getBlock(bid);
    UT_ASSERT2(block);
    UT_ASSERT2(vid>=0&&vid<_nVoxelsInBlock);
#if 0
    return &block->_data[vid];
#else
    return getValueInBlockPtr0(block,vid);
#endif
  }

  void *getValueInBlockPtr( void *blockPtr, int vid ) const
  {
    UT_ASSERT2(vid>=0&&vid<_nVoxelsInBlock);
    return (void*)((uint8_t*)(blockPtr)+vid*_entrySize);
  }

  // Get value in block
  Data_T getValueInBlock(int bid, int vx, int vy, int vz) const
  { 
    Block<Data_T>* block = getBlock(bid);
    return block ? *getValueInBlockPtr0(block,getVoxelInBlockIndex(vx,vy,vz)) : 
      		   _emptyValue;
  }

  // Set value by grid coords (returns error if out-of-memory)
  int setValue( int i, int j, int k,
    		Data_T value, 
		int doZeroDataOnBlockInit = TRUE );

  // Set value by block and voxel id
  int setValueInBlock( int bid, int vid, 
    		       Data_T value,
		       int doZeroDataOnBlockInit = TRUE ); 

  // Apply filter callback function to whole grid
  template<typename func_T>
  void applyFunc( func_T func )
  {
    UT_ASSERT0(_maxBordOff==0);
    // Special case: transform constant-value block
    {
      Data_T fullValue = _fullValue;	  
      UT_ASSERT2( memcmp( (char *)&_fullValue, 
			  (char*)getValueInBlockPtr0(_fullBlock,0),
			  sizeof(_fullValue) ) == 0);
      fullValue = func( fullValue, Vec4i(-1) );
      setFullValue(fullValue);
    }
    {
      Data_T emptyValue = _emptyValue;
      UT_ASSERT2( memcmp( (char *)&_emptyValue, 
			  (char*)getValueInBlockPtr0(_emptyBlock,0),
			  sizeof(_emptyValue) ) == 0);
      emptyValue = func( emptyValue, Vec4i(-1) );
      setEmptyValue(emptyValue);
    }
    // Transform regular blocks
    ThrRunParallel<int>( _nBlocks, nullptr,
      [&](int &tls, int tid, int bid) 
      {
	Data_T *data = getBlockDataPtr(bid);
	if(!data || isConstBlock(getBlock(bid)) ) return;
	SBG_FOR_EACH_VOXEL_IN_BLOCK( this, bid, x, y, z, vid )
	{
	  data[vid] = func( data[vid], Vec4i(x,y,z,0) );
	}
      } 
    );
  }

  // Set Boundary
  int setBoundary( int typ );
  
  // Transfer from/to other grid formats
  void toDenseGrid( Data_T *data );
  void fromDenseGrid( Data_T *data );

  template < typename DataRet_T,
	     typename DataInt_T, // type of intermediate results
	     int ipType_T, 
	     int doRecordMinMax_T > 
  DataRet_T interpolateGen0( float *p0, 
			   unsigned options,
			   float *pMin,
			   float *pMax				
			   ) const;

  // With min,max over interpolation neighborhood
  template <int ipType_T, int doRecordMinMax_T > 
  Data_T interpolateGen( float *p0, unsigned options,
      		      // Optional record min/max of 
		      // the interpolation neighborhood
      		      float *pMin, float *pMax
      			) const
  {
    return interpolateGen0<Data_T,Data_T,ipType_T,doRecordMinMax_T>( 
		p0, options, pMin, pMax );
  }
 
  // Regular Interpolation 

  Data_T interpolateGhost( const Vec4f& pos, unsigned options ) const;

  template <int ipType_T > 
  Data_T interpolate( float *p0, unsigned options ) const
  {
    return interpolateGen<ipType_T,0>( p0, options, NULL, NULL );
  }

  template <int ipType_T > 
  Data_T interpolate( const Vec4f& pos, unsigned options ) const
  {
    float pos_[4] __attribute__((aligned(16)));;
    pos.store_a(pos_);
    return interpolateGen<ipType_T,0>( pos_, options, NULL, NULL );
  }

  template <int ipType_T > 
  float interpolateFloat( const Vec4f& pos, unsigned options ) const
  {
    float val = interpolateFloatFast(pos,options);
    #if 0
    #ifdef UT_ASSERT_LEVEL_2
    float pos_[4] __attribute__((aligned(16)));;
    pos.store_a(pos_);
    float val_= interpolateGen0<float,float,ipType_T,0>( 
		pos_, options, NULL, NULL );

    float err = MT::relativeError( val, val_ );
    if(err>1e-2)
    {
      TRCERR(("err=%g %g %g\n",err,val_,val));
      // for debugging
      /*volatile float val2 = interpolateFloatFast(pos,options),
	      val2_= interpolateGen0<float,float,ipType_T,0>( 
	                 pos_, options, NULL, NULL );
      USE(val); USE(val2_);*/
    }
    #endif
    #endif
    return val;
  }

  template <int ipType_T > 
  float interpolateToFloat( float *p0, unsigned options ) const
  {
    return interpolateGen0<float, float, ipType_T, 0>( 
			p0, options, NULL, NULL );
  }

  template <int doRetValue_T, 
	    int doRetGradient_T > 
  void interpolateLinearWithGradient( 
			      const Vec4f& p0, unsigned options,
			      float *pOutVal,
			      float *pOutGrad ) const;

  template < typename DataRet_T,
	     typename DataInt_T,
	     int dummy_T > // type of intermediate results
  void interpolateWithDerivs0( 
      	int ipType,
	const Vec4f& p0,
	unsigned options,
	float *outVal,
	Vec4f *outGrad
	) const;

  template < int dummy_T >
  void interpolateWithDerivs( 
        int ipType,
	const Vec4f& p0,
	unsigned options,
	float *outVal,
	Vec4f *outGrad
	) 
  {
    interpolateWithDerivs0<float,float,dummy_T>( 
      ipType, p0, options, outVal, outGrad );
  }

  template <int options=0 >
  void interpolateWithSecondDerivs( 
    const Vec4f& p0,
    unsigned opt,
    float *out_f, // OUT: function value
    float *out_grad, // OUT: gradient [fx fy fz]
    float *out_hess // OUT: 2nd derivs [fxx fyy fzz fxy fxz fyz]
    ) const;

  // Sparse Transparency buffer  (x-axis point towards light-dir)
  float lookupTransparencyBuf( float *pos );

  template <int ipType_T> 
  Vec3Float interpolateVec3Float(     		    
	      int dmax,   // stenWidth-1
	      float *pos0,	// Position
	      int iVelComp,	// -1= cell-centerd, >=0: staggered/MAC
	      int options
	      ) const;

  Vec4f interpolateVec3FloatFast( const Vec4f& pos, int opt ) const;

  float interpolateFloatFast( const Vec4f& pos, int opt=0, 
      			      IpolAccessCache *iac=NULL ) const;
  float interpolateUint16ToFloatFast( const Vec4f& pos, int opt=0 ) const;
  float interpolateUint8ToFloatFast( const Vec4f& pos, int opt=0 ) const;


  float interpolateShortToFloatFast( 
      			const Vec4f& pos, int opt ) const;


  template <int ipType_T> 
  Data_T interpolateRef( float *p0, unsigned options ) const;

  float interpolateDenseScalarFast( const Vec4f& p0 ) const;

  float interpolateDenseScalar( float *p0, unsigned options ) const;

  template< int doHalfNeigh = 0 >
  int getSevenPointStencil( 
      		Vec4i pos,  // target point
		Data_T *neighbors    // OUT: pointer to neighbors
      )
  {
    Vec4i bpos = pos >> _bsxLog2;
    int bid = getBlockIndex( vget_x(bpos), vget_y(bpos), vget_z(bpos));
    Data_T *block = getBlockDataPtr0(bid);
    Vec4i vpos = pos & _bsxMask;
    int vid = getVoxelInBlockIndex(vpos);
    UT_ASSERT2(vid>=0&&vid<nVoxelsInBlock());
    Data_T 
      *bp = block+vid,
      *out = neighbors,
      val0 = *bp;
    out[0] = val0;
    if( vall( vpos>vzero() && vpos<_bsxMask ))
    {
      // Fast case: all within block
      out[1] = *(bp-1); 
      out[3] = *(bp-_bsx);  
      out[5] = *(bp-_bsx2);  
      if constexpr(! doHalfNeigh )
      {
	out[2] = *(bp+1); 
	out[4] = *(bp+_bsx);
	out[6] = *(bp+_bsx2);
      }
    }
    else
    {
      int x=vget_x(pos), 
	  y=vget_y(pos), 
	  z=vget_z(pos);

      Vec4i pos1 = pos-1;
      clipGridCoords( pos1 ),
      out[1] = getValue0( vget_x(pos1),y,z );
      out[3] = getValue0( x,vget_y(pos1),z );
      out[5] = getValue0( x,y,vget_z(pos1) );

      if constexpr( !doHalfNeigh )
      {
	Vec4i pos2 = pos+1;
	clipGridCoords( pos2 );
	out[2] = getValue0( vget_x(pos2),y,z );
	out[4] = getValue0( x,vget_y(pos2),z );
	out[6] = getValue0( x,y,vget_z(pos2) );
      }
    }
    return MI_OK;
  }


  template<int width_T, int doAllocateBlocks_T, int isFastCase_T, int writeAccess_T=0 > 
  int getNeighborhood( 
      		int ix0, int iy0, int iz0,  // target point
		Data_T **neighbors    // OUT: pointer to neighbors
      )
  {
    if(doAllocateBlocks_T)
    {
      UT_ASSERT0(FALSE);
    }

    int dmax=width_T-1;
    UT_ASSERT2(_maxBordOff==0);

    if(width_T>2)
    {
      ix0 -= dmax/2;  
      iy0 -= dmax/2;  
      iz0 -= dmax/2;        
    }

    // Locate the block that contains the origin of the stencil window
    int bx0, by0, bz0;
    getBlockCoords(ix0, iy0, iz0, 
		   bx0, by0, bz0);

    if( isFastCase_T )
    {
      // Fast case: 2x2x2 neighborhood entirely within single block
      int vx0, vy0, vz0;
      getVoxelInBlockCoords(ix0, iy0, iz0, 
			    vx0, vy0, vz0 );
      int bid = getBlockIndex(bx0, by0, bz0);
      Data_T *block = getBlockDataPtr(bid,FALSE,FALSE);
      UT_ASSERT2(block);
      if(width_T!=2)
      {
	UT_ASSERT0(FALSE);
      }

      if constexpr  ( writeAccess_T ) 
      {
	if( isReadOnlyBlock(bid) )
	{
	  for(int i=0;i<8;i++) neighbors[0] = NULL;
	  return 0;
	}
      }

      int i = getVoxelInBlockIndex(vx0,vy0,vz0);
      UT_ASSERT2(i<_nVoxelsInBlock-1-_bsx-_bsx2);

      unsigned sz=_entrySize;
      uint8_t *p = (uint8_t*)BYTE_OFFSET(block,i*sz);

      neighbors[0] = (Data_T*)( p + 0);		//	x	y	z
      neighbors[1] = (Data_T*)( p + sz*(1));	//	x+1	y	z
      neighbors[2] = (Data_T*)( p + sz*(_bsx));	//	x	y+1	z
      neighbors[3] = (Data_T*)( p + sz*(_bsx+1));	//	x+1	y+1	z
      neighbors[4] = (Data_T*)( p + sz*(_bsx2));		//	x	y	z+1
      neighbors[5] = (Data_T*)( p + sz*(_bsx2+1));		//	x+1	y	z+1
      neighbors[6] = (Data_T*)( p + sz*(_bsx2+_bsx));	//	x	y+1	z+1
      neighbors[7] = (Data_T*)( p + sz*(_bsx2+_bsx+1));	//	x+1	y+1	z+1
      
      return 8;
    }
    else
    {
      // Generic case
      int iNeigh=0,
          nActNeigh=0;
      for(int k=0;k<width_T;k++)  // TODO unroll loop
      {
	int z = iz0+k,
	    validZ = (z>=0&&z<_sz);
	for(int j=0;j<width_T;j++)
	{
	  int y = iy0+j,
	      validY = validZ && (y>=0&&y<_sy);
	  for(int i=0;i<width_T;i++)
	  {
	    int x = ix0+i,
		validX = validY && (x>=0&&x<_sx);	  
	    if(validX)
	    {
              int ib = ( z >> _bsxLog2 ) *_nxy +
		       ( y >> _bsxLog2 ) *_nx +
		       ( x >> _bsxLog2 );

	      int iv = ((z & _bsxMask ) << _bsx2Log2 )  + 
		       ((y & _bsxMask ) << _bsxLog2 ) + 
			(x & _bsxMask ) ;

	      Data_T* block = NULL;
	      UT_ASSERT2(ib>=0&&ib<_nBlocks);
	      UT_ASSERT2(iv>=0&&iv<_nVoxelsInBlock);
	      block = getBlockDataPtr( ib, FALSE, FALSE );
	      if(block)
	      {
		Data_T *ptr = (Data_T*)BYTE_OFFSET(block,_entrySize*iv);
		if constexpr ( writeAccess_T ) 
		{
		  if( isReadOnlyBlock(ib) )
		  {
	            neighbors[iNeigh] = NULL;
		  }
		  else
		  {
	            neighbors[iNeigh] = ptr;
	            nActNeigh++;
		  }
		}
		else
		{
	          neighbors[iNeigh] = ptr;
	          nActNeigh++;
		}
	      }
	      else neighbors[iNeigh] = NULL;
	    }
	    else neighbors[iNeigh] = NULL;
	    iNeigh++;	   
	  }
	}
      }
      return nActNeigh;
    }
  }

  void setSkipMap( int bid, int bx, int dbx, int dbFull ) 
  {
    UT_ASSERT2(dbx>=SHRT_MIN && dbx<=SHRT_MAX);
    UT_ASSERT2(dbFull>=SHRT_MIN && dbFull<=SHRT_MAX);
    UT_ASSERT2(bid>=0&&bid<_nBlocks);
    UT_ASSERT2(bx+dbx>=-1&&bx+dbx<_nx);
    UT_ASSERT2(_skipmap[bid].nSkip==0);
    _skipmap[bid].nSkip = dbx;
    _skipmap[bid].nFull = dbFull;
  }

  //Data_T lookupTransparencyBufTotal( float *pos );

  //
  // Filtering
  //
  int filter( int typ, 
	      float krn_r,  // Kernel radius
	      float krn_s,  // Sharpening strength
	      unsigned options = 0,
	      SparseGrid<Data_T> **pDstFieldOut = NULL );

  //
  // Misc.
  //

#if 0
  // Return statistics of (nonzero) grid values
  int getValueStats(double zeroThreshold, // IN 
      		    double *pN,
      		    double *pMin,
		    double *pMax,
		    double *pAvg,
		    double *pSig
      ) const;
#endif

  // Return the ratio of non-empty blocks to the total number of blocks	      
  double getFillRate(void) const;

  // Return total used memory 
  size_t getMemUsage(void) const;

  size_t showMemUsage(int trcLevel) const 
  { 
    size_t szUsed = getMemUsage(),
	   szElem = sizeof(_emptyValue);

    if(trcLevel>=0)
    {
      char chbuf[1024];
      sprintf(chbuf, "SBG: '%s' (virt. %dx%dx%d@%d=%.0f) -> %.0f MB\n",
	      _name, (int)_sx,(int)_sy,(int)_sz, (int)szElem, 
	      szElem*_sx*_sy*_sz/((double)ONE_MB),
	      szUsed/((double)ONE_MB));

      switch(trcLevel)
      {
	case 1: { TRCP(("%s",chbuf)); } break;
	case 2: { TRC(("%s",chbuf)); } break;
	case 3: { TRC3(("%s",chbuf)); } break;
      }
    }
    return szUsed;
  }

  // Return true if out-of-memory
  int isOutOfMem( void ) const;

  template< int bsxLog2>
  int getBlockIndexFast( const Vec4f& pos,  /* IN */
      			 Vec4i& ipos, Vec4i& bpos /* OUT (auxillary) */ )
  {
    ipos = truncate_to_int(pos);
    ipos = max(ipos,vzero());
    ipos = min(ipos,_v4iDomMax);
    bpos = ipos >> bsxLog2;
    int bid = horizontal_add( bpos * _v4iGridStridesBlk );
    UT_ASSERT2(bid>=0&&bid<_nBlocks);
    return bid;
  }

  int getBlockIndex( Vec4f posf )
  {
    Vec4i ipos = truncate_to_int(posf);
    UT_ASSERT2( inRange(ipos) );
    Vec4i bpos = getBlockCoords( ipos );
    return getBlockIndex( vget_x(bpos), vget_y(bpos), vget_z(bpos));
  }

  int getBlockIndexTest( Vec4i ipos )
  {
    UT_ASSERT2( inRange(ipos) );
    Vec4i bpos = ipos >> _bsxLog2;
    bpos *= _v4iGridStridesBlk;
    int bid = horizontal_add(bpos);
    return bid;
  }

  // Calculate linear block index 
  int getBlockIndex(int bi, int bj, int bk) const
  {
    UT_ASSERT2(blockCoordsInRange(bi,bj,bk));
    int bid = bk*_nxy + bj*_nx + bi;
    return  bid;
  }

  int getBlockIndex(const Vec4i& ipos) const
  {
    return getBlockIndex(vget_x(ipos),vget_y(ipos),vget_z(ipos));
  }

  void getBlockCoordsFromGridCoords(int i, int j, int k,
      				   int& bi, int& bj, int& bk) const
  {
    //applyBorderOffset(i,j,k);
    getBlockCoords( i, j, k, 
	            bi, bj, bk);
  }

  // Get block coords 
  void getBlockCoords( 
      int i, int j, int k, // IN: grid coords with border offset applied
      int &bi, int &bj, int &bk  // OUT:
      ) const
  {
    bi = i >> _bsxLog2;
    bj = j >> _bsxLog2;
    bk = k >> _bsxLog2;
  }

  Vec4i getBlockCoords( const Vec4i& ipos
      ) const
  {
    return ipos >> _bsxLog2;
  }

  // Get block coords from block id (slow!)
  void getBlockCoordsById( 
      int bid, 
      int &bi, int &bj, int &bk  // OUT:
      ) const
  {
    UT_ASSERT2(bid>=0&&bid<_nBlocks);
    bk = bid / (_nxy);
    int rem = bid % (_nxy);
    bj = rem / _nx;
    bi = rem % _nx;
  }

  Vec4i getBlockCoordsById( int bid )
  {
    UT_ASSERT2(bid>=0&&bid<_nBlocks);
    int bi,bj,bk;
    getBlockCoordsById(bid, bi,bj,bk);
    return Vec4i( bi,bj,bk,0 );
  }

  float getOutOfDomValue( Vec4i pos, float *outOfDomData3x3x3 ) const
  {
    Vec4i dp = pos;
    dp = select( dp < vzero() , Vec4i(-1), dp );
    dp = select( dp >= vzero() && dp <= _v4iDomMax, vzero(), dp );
    dp = select( dp > _v4iDomMax, Vec4i(1), dp);
    dp += 1;
    dp *= Vec4i(1,3,3*3,0);
    int idx = horizontal_add(dp);
    UT_ASSERT2(idx>=0&&idx<3*3*3);
    return outOfDomData3x3x3[idx];
  }

  // Calculate linear block index 
  int getVoxelInBlockIndex(int vi, int vj, int vk) const
  {
    UT_ASSERT2(voxelInBlockRange(vi,vj,vk));
    return (vk<<_bsx2Log2)  + (vj<<_bsxLog2) + vi;
  }

  int getVoxelInBlockIndex(Vec4i ivp ) const
  {
    int vi=vget_x(ivp),
    	vj=vget_y(ivp),
	vk=vget_z(ivp);
    UT_ASSERT2(voxelInBlockRange(vi,vj,vk));
    return (vk<<_bsx2Log2)  + (vj<<_bsxLog2) + vi;
  }
  
  // Get dense grid index
  template< int doBoundCheck=true >
  LongInt getGridIndex(int i, int j, int k) const
  {
    UT_ASSERT2(_isDenseGrid);
    UT_ASSERT2(UT_IMPLIES(doBoundCheck,inRange(i,j,k)));
    applyBorderOffset(i,j,k);
    return (i)+(j)*((LongInt)(_sxAct))+(k)*(LongInt)(_sxyAct);
  }

  LongInt getVirtGridIndex(int i, int j, int k) const
  {
    UT_ASSERT2(inRange(i,j,k));
    return (i)+(j)*((LongInt)(_sx))+(k)*(LongInt)(_sx)*_sy;
  }

  LongInt getVirtGridIndex( const Vec4i& ipos ) const
  {
    return getVirtGridIndex( 
		vget_x(ipos), vget_y(ipos), vget_z(ipos));
  }

  // Get voxel-in-block coords from physical grid coords
  Vec4i getVoxelInBlockCoords( const Vec4i& ipos ) const
  {
    return ipos & _bsxMask;
  }

  void getVoxelInBlockCoords(
	      int i, int j, int k, // IN: grid coords with border offset applied
	      int &vi, int &vj, int &vk) const
  {
    //UT_ASSERT2(i>=0 && j>=0 && k>=0);
    vi = i & _bsxMask;
    vj = j & _bsxMask;
    vk = k & _bsxMask;
  }

  /*void getVoxelInBlockCoords(
      	      const Vec4i& ipos,
	      int &vi, int &vj, int &vk) const
  {
    //UT_ASSERT2(i>=0 && j>=0 && k>=0);
    Vec4i ipv = ipos & _bsxMask;
    vi = vget_x(ipv);
    vj = vget_y(ipv);
    vk = vget_z(ipv);
  }*/

  void getVoxelInBlockCoordsById(
	      int vid, // IN: voxel-in-block-index
	      int &vi, int &vj, int &vk) const
  {
    UT_ASSERT2(vid>=0&&vid<_nVoxelsInBlock);
    vk = vid >>  _bsx2Log2; 
    int rem = vid & (_bsx2-1);
    vj = rem >> _bsxLog2; 
    vi = rem & _bsxMask;
  }

  int distToBorder( int x, int y, int z )
  {
    UT_ASSERT2(inRange(x,y,z));
    int dx0=x, dx1=_sx-1-x,
	dy0=y, dy1=_sy-1-y,
	dz0=z, dz1=_sz-1-z,
	dxDom = MIN( dx0, dx1 ),
	dyDom = MIN( dy0, dy1 ),
	dzDom = MIN( dz0, dz1 );
    return MIN(MIN(dxDom,dyDom),dzDom);
  }

  // Get block pointer by block id
  Block<Data_T>* getBlock( int bid ) const
  {
    UT_ASSERT2(bid>=0&&bid<_nBlocks);
    return _blockmap[bid];
  }

  // Allocate and initialize new block from pool 
  Block<Data_T>* allocBlock( int bid, int doZeroData=0, 
      			     int isConstBlock=0 );

  Block<Data_T>* allocAndZeroInitBlock( int bid, volatile LONG *blockLock,
					uint64_t *nLockCol );

  size_t getBlockSizeBytes( void ) const 
  {
    return _blockPool->getBlockSize();
  }

  #if 0
  int extendDataBlocks( int nBlocks, bool doPushTofreelist )
  {
    void *dummy=NULL;
    int idxExt=-1;
    int rc = _blockPool->extend( &dummy, nBlocks, doPushTofreelist,
				 &idxExt );
    UT_ASSERT0(!rc);
    return rc ? -1 : idxExt;
  }
  #endif

  // Free block
  int freeBlock_( int bid );

  void reset( unsigned options=0 );

  void allocBlockmap( bool doInit=TRUE );
  void prepareDataAccess( unsigned accType = ACC_READ | ACC_WRITE );
  bool hasData( void ) const 
  { 
    return _isDenseGrid ? _denseData != NULL : _blockmap!=NULL; 
  }

  void protectBlock( Block<Data_T>* block );

  void unprotectBlock( Block<Data_T>* block );

  bool isForeignBlock( Block<Data_T>* block )
  {
    return getParentPtr((void*)block) != this;
  }

  bool isForeignBlock( int bid )
  {
    UT_ASSERT2(bid>=0&&bid<_nBlocks);
    return isForeignBlock(_blockmap[bid]);
  }

  bool isReadOnlyBlock( Block<Data_T>* block )
  {
    return (!isValueBlock( block )) || isForeignBlock( block );
  }

  bool isReadOnlyBlock( int bid )
  {
    UT_ASSERT2(bid>=0&&bid<_nBlocks);
    Block<Data_T>* block=_blockmap[bid];
    return isReadOnlyBlock(block);
  }

  bool isValueBlock( int bid ) const
  {
    UT_ASSERT2(bid>=0&&bid<_nBlocks);
    Block<Data_T>* block=_blockmap[bid];
    return isValueBlock(block);
  }

  bool isValueBlock( Block<Data_T>* block ) const
  {
    return !(block==NULL || isConstBlock(block));
  }

  bool isConstBlock( Block<Data_T>* block ) const
  {
    return block==NULL || (block->uheader.flags & SBLK_CONSTANT);
  }

  bool isEmptyBlock( int bid )
  {
    UT_ASSERT2(bid>=0&&bid<_nBlocks);
    Block<Data_T>* block=_blockmap[bid];
    return isEmptyBlock(block);
  }

  bool isEmptyBlock( Block<Data_T>* block )
  {
    return block==NULL || (block->uheader.flags & SBLK_EMPTY);
  }

  bool isFullBlock( int bid )
  {
    UT_ASSERT2(bid>=0&&bid<_nBlocks);
    Block<Data_T>* block=_blockmap[bid];
    return isFullBlock(block);
  }

  bool isFullBlock( Block<Data_T>* block )
  {
    return block && (block->uheader.flags & SBLK_FULL);
  }


  bool isInvalidBlock( int bid )
  {
    UT_ASSERT2(bid>=0&&bid<_nBlocks);
    Block<Data_T>* block=_blockmap[bid];
    return isInvalidBlock(block);
  }

  bool isInvalidBlock( Block<Data_T>* block )
  {
    return block==NULL || (block->uheader.flags & SBLK_INVALID);
  }

  void setBlock( int bid, Block<Data_T>* block )
  {
    UT_ASSERT2(block);
    UT_ASSERT2(!isValueBlock(bid));
    _blockmap[bid] = block;
  }

  // Set constant value block
  void setFullBlock( int bid )
  {
    UT_ASSERT2(_fullBlock);
    UT_ASSERT2(!isValueBlock(bid));
    _blockmap[bid] = _fullBlock;
  }

  void setEmptyBlock( int bid )
  {
    UT_ASSERT2(_emptyBlock);
    UT_ASSERT2(!isValueBlock(bid));
    _blockmap[bid] = _emptyBlock;
  }

  void setInvalidBlock( int bid )
  {
    UT_ASSERT2(_invalidBlock);
    UT_ASSERT2(!isValueBlock(bid));
    _blockmap[bid] = _invalidBlock;
  }

  void copyConstValues( SparseGrid<Data_T>*sgDest )
  {  
    sgDest->setFullValue( _fullValue );
    sgDest->setEmptyValue( _emptyValue );
    sgDest->setInvalidValue( _invalidValue );
  }

  void copyConstBlock( SparseGrid<Data_T>*sgDest, int bid )
  {
    if( isFullBlock(bid) ) 
    {
      UT_ASSERT2( compareValuesMem(&sgDest->_fullBlock->_data[0],
	    			   &_fullBlock->_data[0])==0 );
      sgDest->setFullBlock(bid);
    }
    else if(isEmptyBlock(bid)) 
    {
      UT_ASSERT2( compareValuesMem(&sgDest->_emptyBlock->_data[0],
	    			   &_emptyBlock->_data[0])==0 );
      sgDest->setEmptyBlock(bid);
    }
    else if(isInvalidBlock(bid)) 
    {
      UT_ASSERT2( compareValuesMem(&sgDest->_invalidBlock->_data[0],
	    			   &_invalidBlock->_data[0])==0 );
      sgDest->setInvalidBlock(bid);
    }
    else
    {
      UT_ASSERT0(FALSE);
    }
  }

  //
  Data_T getFullValue( void ) const
  { 
    #if 0
    UT_ASSERT2( memcmp( (char *)&_fullValue, 
	  		(char*)getValueInBlockPtr0(_fullBlock,0),
			sizeof(_fullValue) ) == 0);
    #endif
    return _fullValue; 
  }

  //
  Data_T getEmptyValue( void ) const
  { 
    #if 0
    UT_ASSERT2( memcmp( (char *)&_emptyValue, 
	  		(char*)getValueInBlockPtr0(_emptyBlock,0),
			sizeof(_emptyValue) ) == 0);
   #endif
    return _emptyValue; 
  }

  Block<Data_T> *getEmptyBlock( void ) const
  { 
    return _emptyBlock; 
  }

  int compareValuesMem(  Data_T *pval1, Data_T *pval2 ) const
  {   
    if(sizeof(Data_T) > sizeof(Vec3Float)) return 0; // cant compare structs bc padding
    return memcmp( (char *)pval1, (char*)pval2, sizeof(*pval1) );
  }

  void checkConstBlocks(bool doFullCheck)
  {
    Block<Data_T> *blocks[] = { _fullBlock, _emptyBlock, _invalidBlock };
    Data_T *vals[] = { &_fullValue, &_emptyValue, &_invalidValue };
    for(int ib=0; ib<ARRAY_LENGTH(blocks); ib++)
    {
      Block<Data_T> *block = blocks[ib];
      Data_T *val = vals[ib];
      if(block)
      {
        for(int i=0;i<(doFullCheck ? _nVoxelsInBlock : 1);i++)
	{
	  if(!( compareValuesMem(getValueInBlockPtr0(block,i),val)==0 ))
	  {
	    TRCERR(("%s: name='%s' ib=%d i=%d\n",
		  UT_FUNCNAME,_name,ib,i)); 
	  }
	}
      }
    }
  }

  //
  void setConstBlockValue( Block<Data_T> *block, Data_T *pValue, Data_T val ) 
  { 
    Data_T oldVal = *pValue;
    if(block && (compareValuesMem( pValue, &val) !=0) )
    {
      *pValue = val;
      UT_ASSERT0( compareValuesMem( getValueInBlockPtr0(block,0), &oldVal)==0 );
      unprotectBlock(block);
      for(int i=0;i<_nVoxelsInBlock;i++) *getValueInBlockPtr0(block,i)=val;
      if(_doProtectConstBlocks) protectBlock(block);
    }
  }

  //
  void setFullValue( Data_T val ) 
  { 
    setConstBlockValue( _fullBlock, &_fullValue, val );
  }

  //
  void setEmptyValue( Data_T val ) 
  { 
    setConstBlockValue( _emptyBlock, &_emptyValue, val );
  }

  //
  void setInvalidValue( Data_T val ) 
  { 
    setConstBlockValue( _invalidBlock, &_invalidValue, val );
  }

  //
  void setConstValue( int constval, Data_T val ) 
  { 
    switch(constval)
    {
      case CONSTVAL_EMPTY: setEmptyValue( val ); break;
      case CONSTVAL_FULL: setFullValue( val ); break;
      case CONSTVAL_INVALID: setInvalidValue( val ); break;
      default: UT_ASSERT0(FALSE); break;
    }
  }

  //
  void setValueScaleAndOffs( float scale, float offs )
  {
    _valueScale = scale;
    _valueOffs = offs;
  }

  uint64_t getDataGenCount( void ) { return _dataGenCount; }
  void setDataGenCount( uint64_t count ) { _dataGenCount = count; }

  //
  // Memory / Data related
  //
  void swapData( SparseGrid<Data_T>* sg )
  {
    TRC(("SBG::swapData() '%s' <-> '%s'\n",_name,sg->_name));
    UT_SWAP_PTR_(Block<Data_T>*, _blockmap, sg->_blockmap);
    UT_SWAP_PTR_(SkipMapEntry, _skipmap, sg->_skipmap);
    UT_SWAP_PTR_(BlockPool, _blockPool, sg->_blockPool);
    UT_SWAP_PTR_(BlockPool, _blockPoolConst, sg->_blockPoolConst);
    UT_SWAP_(Data_T, _emptyValue, sg->_emptyValue);
    UT_SWAP_(Data_T, _fullValue, sg->_fullValue);
    UT_SWAP_(Data_T, _invalidValue, sg->_invalidValue);
    UT_SWAP_(float,_valueScale,sg->_valueScale);
    UT_SWAP_(float,_valueOffs,sg->_valueOffs);
    UT_SWAP_(uint64_t,_dataGenCount,sg->_dataGenCount);
  }
    
  private:

  void initializeBlockmap( ConstVal constBlock = CONSTVAL_NONE )
  {
    if(_blockmap||_skipmap)
    {
      Block<Data_T> *initBlock = NULL;
      switch(constBlock)
      {
	case CONSTVAL_EMPTY: initBlock = _emptyBlock; break;
	case CONSTVAL_FULL: initBlock = _fullBlock; break;
	case CONSTVAL_INVALID: initBlock = _invalidBlock; break;
	default: break;
      }
      UT_ASSERT0( UT_IMPLIES(constBlock != CONSTVAL_NONE, initBlock  ));
      UT_ASSERT0(_emptyBlock);

      TRC3(("%s: name='%s' initBlock = %p\n",UT_FUNCNAME,_name,initBlock));

      ThrRunParallel<int>( _nBlocks, nullptr, 
      [&](int &tls, int tid, int bid) 
      {
	  if(_blockmap) _blockmap[bid] = initBlock;
	  if(_skipmap) _skipmap[bid]={};
      } );
    }
  }

  void bmapResetMap( uint32_t *bmap, int bid )
  {
    if(_blockBitmapExtPtr) memset(bmap, 0, _blockBitmapLen*4);
  }

  void bmapSet( uint32_t *bmap, int bid ) 
  {
    uint32_t iWord = bid>>5,
	     iBit = bid&31;
    UT_ASSERT2(iWord<_blockBitmapLen);
    bmap[iWord] |= (((uint32_t)1)<<(iBit));
  }

  void bmapReset( uint32_t *bmap, int bid ) 
  {
    uint32_t iWord = bid>>5,
	     iBit = bid&31;
    UT_ASSERT2(iWord<_blockBitmapLen);
    bmap[iWord] &= ~(((uint32_t)1)<<(iBit));
  }

  bool bmapIsSet( uint32_t *bmap, int bid ) const
  {
    uint32_t iWord = bid>>5,
	     iBit = bid&31;
    UT_ASSERT2(iWord<_blockBitmapLen);
    return bmap[iWord] & (((uint32_t)1)<<(iBit));
  }

  //
  // Block address calculations
  //

  // Apply border offset to convert logical to 'physical' grid coords
  void applyBorderOffset( int& i, int& j, int& k ) const
  {
    i += _maxBordOff;
    j += _maxBordOff;
    k += _maxBordOff;
  }
    
  void initializeBlock( Block<Data_T>* block, int doZeroData )
  {
    // Data
    if(doZeroData)
    {
      for(int i=0;i<_nVoxelsInBlock;i++) 
	*getValueInBlockPtr0(block,i)=_emptyValue;
    }
  }
  
  // 
  //  Member variables
  //
  char	    _name[80];
  
  int      _optionsOrig;
  int      _isDenseGrid,
	   _hasCompTimeEntrySz;
  int	   _maxBordOff; // Max. border offset used for padding

  int	   _bsx,	// block size in voxels (at max. resolution level) power-of-two
	   _nx,_ny,_nz, // number of blocks in each dimension
	   _nhx, _nhy, _nhz, _nhxy,  // for halo blockmap
	   _sx,_sy,_sz, // logical grid size in voxels
	   _sxyzMax,
	   _sxyzMin;

  int	   _sxAct,_syAct,_szAct;  // actual grid resolution incl. border/padding
  LongInt  _sxyAct;
  LongInt  _nxy;

  float    _nbDist;  /* narow band distance for distance fields (-1=N/A) */

  Vec3Int  _gridRes,
	   _gridResBlk;

  Vec4i    _v4iStridesLerp0246,
           _v4iStridesLerp1357;

  Vec4i    _v4iStrides3Lerp0246,
           _v4iStrides3Lerp1357;
  Vec8i    _v8iStridesLerp02461357,
           _v8iStrides3Lerp02461357,
           _v8iStrides3Lerp01234567;


  Vec4i	   _v4iGridStridesBlk,
	   _v4iBlkCoordsMax;
  Vec4f    _v4fDomMin,
    	   _v4fDomMax;
  Vec4f    _v4fDomMin2,
    	   _v4fDomMax2;
  Vec4i    _v4iDomMin,
    	   _v4iDomMax;

  DomainBC _domainBC;

  int	   _bsx2,    // bsx^2 = stride in z-direction
    	   _bsxLog2,
	   _bsx2Log2,
	   _bsxMask;

  int	   _neighVoxOffs[8];  

  LongInt  _nTotActVoxels,
	   _nTotVirtVoxels;
  int      _nVoxelsInBlock,
	   _nVoxelsInBlockAct,
	   _nBlocks;

  unsigned _entrySize;

  size_t   _blockSizeBytes,
	   _blockProtectSize;

  LongInt  _initialSizeMB,
	   _maxSizeMB;

  int	   _isLowOnMemory,
	   _doProtectConstBlocks;
#if 0
  struct			// Mapping grid from/to world space coords
  {
    double pos[4],
	   dx, dx_r;
  }
  _wsInfo;
#endif

  // 
  // Memory/Data Management  (TODO: multiple channels)
  // XXX: Keep consistent with swapData() method !
  //
  float    _valueScale,  // Optional linear mapping for
	   _valueOffs;  // transforming uint16_t to float
  Data_T   _emptyValue,
	   _fullValue,
	   _invalidValue;
  
  Data_T   *_denseData;
  UtMutex  _lock;
  
  Block<Data_T> **_blockmap;

  uint64_t _dataGenCount;

#if 0
  Block<Data_T> **_blockmapDB[6];  // beyond domain bonudary 
#endif

  uint32_t *_blockBitmapExtPtr,
	   _blockBitmapLen;
  
  //uint8_t *_blockFlags;
    
  Block<Data_T> *_fullBlock,
    	        *_emptyBlock,
		*_invalidBlock;

  typedef struct
  {
    int16_t nSkip, 
	    nFull;
  }
  SkipMapEntry;

  SkipMapEntry *_skipmap;
   
  BlockPool *_blockPool,
            *_blockPoolConst;
};

//
// Utility
//
//

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
class IpolAccessCache
{
  public:

  IpolAccessCache(void) { reset(); };

  void reset(int tid=0) { 
    _iposl = Vec4i(INT32_MIN, /* x */
		  INT32_MIN, /* y */
		  INT32_MIN, /* z */
		  INT32_MIN  /* level */
		  ); 
    _nHit=0; 
    _nAcc=0;
    _tid = tid;
  }

  Vec4f __attribute__((aligned(CPU_CACHELINE_SZ)))
    _g0246,
    _g1357;

  Vec4i _iposl;

  LongInt _nHit, _nAcc;

  int _tid;
};

#ifdef SBG_IAC_MULT_TEST
class IpolAccessCache2
{
  public:

  IpolAccessCache2(void) { reset(); };

  void reset(int tid=0) { 
    for(int i=0;i<ARRAY_LENGTH(_entries);i++)
    {
      _entries[i].iposl = Vec4i(INT32_MIN, /* x */
		    INT32_MIN, /* y */
		    INT32_MIN, /* z */
		    INT32_MIN  /* level */
		    ); 
    }
    _iEntryLastUsed=0;
    _nHit=0; 
    _nAcc=0;
    _tid = tid;
  }

  struct
  {
    Vec8f __attribute__((aligned(CPU_SIMD_WIDTH)))
      a0246b0246,
      a1357b1357;
    Vec4i iposl;
  }
  _entries[2];

  int _iEntryLastUsed;

  LongInt _nHit, _nAcc;

  int _tid;
};

#else
class IpolAccessCache2
{
  public:

  IpolAccessCache2(void) { reset(); };

  void reset(int tid=0) { 
    _iposl = Vec4i(INT32_MIN, /* x */
		  INT32_MIN, /* y */
		  INT32_MIN, /* z */
		  INT32_MIN  /* level */
		  ); 
    _nHit=0; 
    _nAcc=0;
    _tid = tid;
  }

  Vec8f __attribute__((aligned(CPU_CACHELINE_SZ)))
    _a0246b0246,
    _a1357b1357;

  Vec4i _iposl;

  LongInt _nHit, _nAcc;

  int _tid;
};
#endif

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/

#if 0
class AccessCache
{
  public:
  
  AccessCache( void ) :
    bpos(Vec4i(-1,-1,-1,0)),            
    bid(-1),
    nAcc(0),
    nHit(0)
    { };

  Vec4i bpos;
  int bid;
  int level;
  float auxInvScale;
  int bsx,
      bsxMask;
  Vec4i bsxStrides;

  void *dataFlags;
  Vec3Float *dataVec,
	    *dataVecLo;
  float *dataFloat,
	*dataFloatLo;

  unsigned ipOpt;

  int	   oo0,oo1,oo2,oo3,
	   oo4,oo5,oo6,oo7;

  LongInt nAcc,
	  nHit;
};
#endif
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
class InterlockedIterator
{
  volatile LONG _iTask;
  int _nBlocks;

  public:
  InterlockedIterator( int nBlocks ) : 
    _iTask(-1), _nBlocks(nBlocks)  { }

  int nBlocks() const { return _nBlocks; }

  int getNextBlock() 
  { 
    return InterlockedIncrement((LONG volatile *)&_iTask); 
  }

  void reset(void)
  {
    _iTask = -1;
  }
};

/*-------------------------------------------------------------------------*/
/* getValue() 								   */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
Data_T* SparseGrid<Data_T>::getValuePtr( int i, int j, int k ) const
{
  UT_ASSERT2( inRange(i,j,k) );
  
  // First step ist to apply the border offset. From now on we are only
  // working with 'internal' coords on the larger, padded grid
  //applyBorderOffset( i, j, k );
  
  // Find block 
  int bi, bj, bk;
  getBlockCoords(i, j, k, bi, bj, bk);
  
  LongInt id = getBlockIndex(bi, bj, bk);

  Block<Data_T>* block = _blockmap[id];
  if(block)
  {
    // Find voxel in block
    int vi, vj, vk;
    getVoxelInBlockCoords(i, j, k, 
			 vi, vj, vk);
    int vidx = getVoxelInBlockIndex(vi,vj,vk);
    UT_ASSERT2(vidx>=0&&vidx<_nVoxelsInBlock);
    return getValueInBlockPtr0(block,vidx);
  }
  else
  {
    return NULL;
  }
}

/*-------------------------------------------------------------------------*/
/* getValue() 								   */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
Data_T SparseGrid<Data_T>::getValue( int i, int j, int k ) const
{
  UT_ASSERT2( inRange(i,j,k) );
  
  // First step ist to apply the border offset. From now on we are only
  // working with 'internal' coords on the larger, padded grid
  //applyBorderOffset( i, j, k );
  
  // Find block 
  int bi, bj, bk;
  getBlockCoords(i, j, k, bi, bj, bk);
  
  LongInt id = getBlockIndex(bi, bj, bk);

  Block<Data_T>* block = _blockmap[id];
  if(block)
  {
    // Find voxel in block
    int vi, vj, vk;
    getVoxelInBlockCoords(i, j, k, 
			 vi, vj, vk);
    int vidx = getVoxelInBlockIndex(vi,vj,vk);
    UT_ASSERT2(vidx>=0&&vidx<_nVoxelsInBlock);
    return *getValueInBlockPtr0(block,vidx);
  }
  else
  {
    return _emptyValue;
  }
}

/*-------------------------------------------------------------------------*/
/* getValue0() 								   */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
Data_T SparseGrid<Data_T>::getValue0( int i, int j, int k ) const
{
  UT_ASSERT2( inRange(i,j,k) );
  
  //applyBorderOffset( i, j, k );
  UT_ASSERT2(_maxBordOff==0);
  
  // Find block 
  int bi, bj, bk; 
  getBlockCoords(i, j, k, bi, bj, bk);  
  int id = getBlockIndex(bi, bj, bk);
  Block<Data_T>* block = _blockmap[id];
  UT_ASSERT2(block);
  // Find voxel in block
  int vi, vj, vk;
  getVoxelInBlockCoords(i, j, k, 
		       vi, vj, vk);
  int vidx = getVoxelInBlockIndex(vi,vj,vk);
  UT_ASSERT2(vidx>=0&&vidx<_nVoxelsInBlock);
  return *getValueInBlockPtr0(block,vidx);
}

/*-------------------------------------------------------------------------*/
/* setValue() 								   */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
int SparseGrid<Data_T>::setValueInBlock( int bid, int vid, 
    				   	 Data_T value,
					 int doZeroDataOnInit ) 
{
  int rcThis=0;

  UT_ASSERT2(bid>=0&&bid<_nBlocks);
  UT_ASSERT2(vid>=0&&vid<_nVoxelsInBlock);

  // Obtain access to the block pointer
  Block<Data_T>* block = _blockmap[bid];
  if(block)
  {
    UT_ASSERT0( !isReadOnlyBlock(block));
  }
  if(!block)
  {
    // Allocate missing block
    block = allocBlock(bid,doZeroDataOnInit);
    if(!block) return MI_ENOMEM;
  }

  // set value
  *getValueInBlockPtr0(block,vid) = value;
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* setValue() 								   */
/*-------------------------------------------------------------------------*/
template <typename Data_T>
int SparseGrid<Data_T>::setValue( int i, int j, int k,
    				   Data_T value,
				   int doZeroDataOnInit ) 
{
  UT_ASSERT2( inRange(i,j,k) );

  // First step ist to apply the border offset. From now on we are only
  // working with 'internal' coords on the larger, padded grid
  //applyBorderOffset( i, j, k );

  // Find block 
  int bi, bj, bk;
  getBlockCoords(i, j, k, bi, bj, bk);
  
  int bid = getBlockIndex(bi, bj, bk);

  // Find voxel in block
  int vi, vj, vk;
  getVoxelInBlockCoords(i, j, k, 
		        vi, vj, vk);
  int vid = getVoxelInBlockIndex(vi,vj,vk);

  return setValueInBlock( bid, vid, value, doZeroDataOnInit );
}

/*=========================================================================
 *
 *
 *  Interpolation methods   
 *
 *
 * =======================================================================*/

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static inline double MING( float& x, double& y )
{
  return MIN(x,y);
}
static inline double MAXG( float& x, double& y )
{
  return MAX(x,y);
}

static inline float MING( float& x, float& y )
{
  return MIN(x,y);
}
static inline float MAXG( float& x, float& y )
{
  return MAX(x,y);
}

static inline float MING( float& x, Vec3Float& y )
{
  return 0;
}
static inline float MAXG( float& x, Vec3Float& y )
{
  return 0;
}
static inline float MING( float& x, uint16_t& y )
{
  return 0;
}
static inline float MAXG( float& x, uint16_t& y )
{
  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static inline void makeMonotone( double& d1, double& d2 )
{
  if(!(d1*d2>0.0f))  //  if(MtSign(d1)!=MtSign(d2))
  {
    d1=d2=0.0f;
  }
}

static inline void makeMonotone( float& d1, float& d2 )
{
  if(!(d1*d2>0.0f))  //  if(MtSign(d1)!=MtSign(d2))
  {
    d1=d2=0.0f;
  }
}

static inline void makeMonotone( Vec3Float& d1, Vec3Float& d2 )
{
  for(int k=0;k<3;k++)
  {
    if(!(d1[k]*d2[k]>0.0f))  //  if(MtSign(d1)!=MtSign(d2))
    {
      d1[k]=d2[k]=0.0f;
    }
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static inline void makeMonotone2( double d0, double d1, double d2,
    				  double& m0, double& m1 )
{
  UT_ASSERT0(FALSE);
}

static inline void makeMonotone2( float d0, float d1, float d2,
    				  float& m0, float& m1 )
{
  // Fritsch-Carlsen
  m0 = ( d1*d0 > (1.0E-12f) ) ? 6.0f / ( 3.0f/d0 + 3.0f/d1 ) :
    			   	0;
  m1 = ( d2*d1 > (1.0E-12f) ) ? 6.0f / ( 3.0f/d1 + 3.0f/d2 ) :
    			   	0;
}

static inline void makeMonotone2( Vec3Float& d0, Vec3Float& d1, Vec3Float& d2,
    				  Vec3Float& m0, Vec3Float& m1 )
{
  for(int k=0;k<3;k++)
  {
    makeMonotone2(d0[k],d1[k],d2[k],m0[k],m1[k]);
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static inline void ipol1DClamp( float& p1, float& p2, float& val )
{
  if(p1<p2)
  {
    if(val<p1) val=p1;
    if(val>p2) val=p2;
  }
  else
  {
    if(val<p2) val=p2;
    if(val>p1) val=p1;
  }
}

static inline void ipol1DClamp( Vec3Float& p1, Vec3Float& p2, Vec3Float &val )
{
  for(int k=0;k<3;k++)
  {
    if(p1[k]<p2[k])
    {
      if(val[k]<p1[k]) val=p1[k];
      if(val[k]>p2[k]) val=p2[k];
    }
    else
    {
      if(val[k]<p2[k]) val=p2[k];
      if(val[k]>p1[k]) val=p1[k];
    }
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static inline
void wenoParabola( double f1, double f2, double f3, float x, // IN
    		     double *w, double *g // OUT 
		     )
{
  UT_ASSERT0(FALSE);
}

static inline
void wenoParabola( float f1, float f2, float f3, float x, // IN
    		     float *w, float *g // OUT 
		     )
{
  float d = (f3-f1)*0.5f,	// 1st derivative at 0
	 dd = (f1-2.f*f2+f3),	// 2nd derivative
	 S = d*(d+dd) + (4.f/3.f)*dd*dd;  // Smoothness over x=[0,1]
   *w = (2.f-x)/MT_SQR(1e-6f+S);	// Relative weight
   *g = f2 + x*(d+0.5f*x*dd);	// Value of parabola at x
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static inline
void wenoParabola( Vec3Float f1, Vec3Float f2, Vec3Float f3, float x, // IN
    		     Vec3Float *w, Vec3Float *g // OUT 
		     )
{
  Vec3Float d = (f3-f1)*0.5f,	// 1st derivative at 0
	    dd = (f1-f2*2.0f+f3),	// 2nd derivative
	    S = d*(d+dd) + dd*dd*(4.f/3.f);  // Smoothness over x=[0,1]
   *w = Vec3Float(2.0f-x)/((S+1e-6f)*(S+1e-6f));	// Relative weight
   *g = f2 + (d+0.5f*dd*x)*x;	// Value of parabola at x
}

static inline
float Interpolate1D_MINMOD(float *p, float t) 
{
  float a = p[1]-p[0],
	b = p[2]-p[1];
  float sgn_a = copysignf(1.0f,a),
	abs_a = fabsf(a), 
	sgn_b = copysignf(1.0f,b),
	abs_b = fabsf(b);
  return 0.5f * (sgn_a + sgn_b) * std::min(abs_a, abs_b);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template<typename Data_T, typename DataInt_T, int ipType_T>
DataInt_T Interpolate1D(DataInt_T *p, float t) 
{
  if(ipType_T==IP_LINEAR)
  {
    return p[1]*t+p[0]*(1.0f-t);
  }
  else if(ipType_T==IP_CUBIC)
  {
   return  p[1] + (p[2] - p[0] + (p[0]*2.0f - p[1]*5.0f + 
      p[2]*4.0f - p[3] + ((p[1] - p[2])*3.0f + p[3] - p[0])*t )*t  )* 0.5f * t;
  }
  else if(ipType_T==IP_CUBIC_MONO)
  {
    DataInt_T  
      d1 = (p[2]-p[0])*0.5f,
      d2 = (p[3]-p[1])*0.5f,
      d0 = p[2]-p[1];

    makeMonotone( d1, d2 );
    
    DataInt_T 
      a0 = p[1],
      a1 = d1,
      a2 = d0*3.0f-d1*2.0f-d2,
      a3 = d1+d2-d0*2.0f,

      g = a3*t*t*t+a2*t*t+a1*t+a0;
    return g;    
  }
  else if(ipType_T==IP_CUBIC_MONO_2)
  {
    DataInt_T  
        d0 = p[1]-p[0],
	d1 = p[2]-p[1],
	d2 = p[3]-p[2];

    DataInt_T m0,m1;

    makeMonotone2(d0,d1,d2, m0, m1);

    DataInt_T 
      s = d1,
      a = s*3.0f - m0*2.0f - m1,
      b = m0 + m1 - s*2.0f,

      g = p[1] + m0*t + a*t*t + b*t*t*t;
    
    return g;
  }
  else if(ipType_T==IP_WENO4)
  {
    DataInt_T w1,p1;
    wenoParabola( p[0], p[1], p[2], t,
		  &w1, &p1 );
    DataInt_T w2,p2;
    wenoParabola( p[3], p[2], p[1], 1-t,
		  &w2, &p2 );
    return (w1*p1 + w2*p2) / (w1+w2);  
  }
  else
  {
    UT_ASSERT0(FALSE);
  }
}

/*-------------------------------------------------------------------------*/
/*									   */
/*-------------------------------------------------------------------------*/
template < typename Data_T > // main data type of the interpolant
Data_T SparseGrid<Data_T>::interpolateGhost( 
    				const Vec4f& p0,
				unsigned options
				) const
{
  Data_T buf_x[4],buf_y[4],buf_z[4];
  Data_T f_out = zeroValue<Data_T>();
  float po=options&OPT_IPCORNER?0.0f:0.5f;

  float
	 x = p0[0]-po, 	// po=0.5 -> assume cell centered grid
	 y = p0[1]-po, 
	 z = p0[2]-po;

  int    ix = floorf(x), 
	 iy = floorf(y), 
	 iz = floorf(z);

  int ix0 = ix,
      iy0 = iy,
      iz0 = iz;

  constexpr int dmax=1;

  float u = x - ix,
	v = y - iy,
	w = z - iz;

  ix0 -= dmax/2;  
  iy0 -= dmax/2;  
  iz0 -= dmax/2;  

  int bx0, by0, bz0, bid;
  int vx0, vy0, vz0;

  getBlockCoords(ix0, iy0, iz0, 
		 bx0, by0, bz0);
  bid = getBlockIndex(bx0, by0, bz0);
  getVoxelInBlockCoords(ix0, iy0, iz0, 
			vx0, vy0, vz0 );

  Data_T val0=_emptyValue; 
  bool haveVal0=FALSE;

  //Vec4i block0Pos( bx0*_bsx, by0*_bsx, bz0*_bsx, 0 );

  for(int k=0;k<2;k++)
  {
    int ibz = ( bz0 + ((vz0+k) >> _bsxLog2) ) *_nxy,
	ivz = (((vz0+k) & _bsxMask) << _bsx2Log2),
	iz = iz0 + k;
    for(int j=0;j<2;j++)
    {
      int iby = ibz + (by0 + ((vy0+j) >> _bsxLog2)) *_nx,
	  ivy = ivz + (((vy0+j) & _bsxMask) <<_bsxLog2),
	  iy = iy0+j;
      for(int i=0;i<2;i++)
      {
	int ib = iby + bx0 + ((vx0+i) >> _bsxLog2),
	    iv = ivy + ((vx0+i) & _bsxMask),
	    ix = ix0+i;

	Data_T val;

	Vec4i ipos(ix,iy,iz,0);

	if(!inRange(ipos))
	{
	  clipGridCoords(ipos);
	  val = getValue(ipos);
	}
	else if(!_blockmap[ib])
	{
	  if(!haveVal0)
	  {
	    val0 = getValue(p0);
	    haveVal0=TRUE;
	  }
	  val = val0;
	}
	else
	{
	  UT_ASSERT2(ib>=0&&ib<_nBlocks);	
	  Block<Data_T>* block = _blockmap[ib]; 	  
	  UT_ASSERT2(iv>=0&&iv<_nVoxelsInBlock);
	  val = block->_data[ iv ]; 
	}
	buf_x[i] = val;
      }
      buf_y[j] = buf_x[1]*u + buf_x[0]*(1.0f-u);
    }
    buf_z[k] = buf_y[1]*v + buf_y[0]*(1.0f-v);  
  }
  f_out = buf_z[1]*w + buf_z[0]*(1.0f-w);

  return f_out;
}

/*-------------------------------------------------------------------------*/
/*									   */
/*-------------------------------------------------------------------------*/
template < typename Data_T > // main data type of the interpolant
template < typename DataRet_T,
	   typename DataInt_T, // type of intermediate results
           int ipType_T, 
	   int doRecordMinMax_T > 
DataRet_T SparseGrid<Data_T>::interpolateGen0( 
    				float *p0, 
				unsigned options,
				float *pMin,
				float *pMax				
				) const
{
  //FLAG_RESET(options,OPT_IPCHKBLK);
  //
  UT_ASSERT2(_hasCompTimeEntrySz);

  DataInt_T buf_x[4],buf_y[4],buf_z[4];
  DataRet_T f_out = zeroValue<DataRet_T>();

  if(doRecordMinMax_T)
  {
    *pMin = 0.0f;
    *pMax = 0.0f;
  }

  float po=options&OPT_IPCORNER?0.0f:0.5f;

  //if(ipType_T==IP_NEAREST) po-=1;

  float
	 x = p0[0]-po, 	// po=0.5 -> interpolate at pixel centers
	 y = p0[1]-po, 
	 z = p0[2]-po;


  int    ix = floorf(x), 
	 iy = floorf(y), 
	 iz = floorf(z);

  int ix0 = ix,
      iy0 = iy,
      iz0 = iz;

  // First step ist to apply the border offset. From now on we are only
  // working with 'internal' coords on the larger, padded grid
  
  // dmax = 'stencil window size' - 1 
  int dmax = ipType_T==IP_NEAREST ? 0 :
    	     ipType_T==IP_LINEAR ? 1 :
	     ( ipType_T==IP_CUBIC || ipType_T==IP_CUBIC_MONO || 
	       			     ipType_T==IP_CUBIC_MONO_2 ||
				     ipType_T==IP_WENO4
				     ) ? 3 : 0;      
  float u = x - ix,
	v = y - iy,
	w = z - iz;

  if( options & OPT_IP_FADE_WEIGHTS )
  {
    fadeIpolWeights(u,v,w);
  }

  ix0 -= dmax/2;  
  iy0 -= dmax/2;  
  iz0 -= dmax/2;  

  // 
  // Stencil window out of total range 
  //
  if( ! ( ix0>=0 && iy0>=0 && iz0>=0 && 
	ix0+dmax<_sx && iy0+dmax<_sy && iz0+dmax<_sz )) 
  {
    // 
    // Generic version with domain boundary condition handling
    //
    if(options & OPT_IPBC_DIRICHLET_OLD) return _emptyValue;
 
    UT_ASSERT0(!doRecordMinMax_T);

    int iDir=-1,iFace=-1;
    for(int k=0;k<=dmax;k++)
    {
      int iz = iz0+k,
	  iz2 = iz,
	  clippedZ = 0;
      if(iz2<0) {    iz2=0; 	clippedZ=1; iDir=2; iFace=0; } 
      if(iz2>=_sz) { iz2=_sz-1; clippedZ=1; iDir=2; iFace=1; }

      for(int j=0;j<=dmax;j++)
      {
	int iy = iy0+j,
	    iy2 = iy,
	    clippedY = clippedZ;
	if(iy2<0) { 	iy2=0;     clippedY=1; iDir=1; iFace=0; } 
	if(iy2>=_sy) {  iy2=_sy-1; clippedY=1; iDir=1; iFace=1; } 

	for(int i=0;i<=dmax;i++)
	{
	  int ix = ix0+i,
	      ix2 = ix,
	      clipped = clippedY;
	  if(ix2<0) {    ix2 = 0;     clipped=1; iDir=0; iFace=0; } 
	  if(ix2>=_sx) { ix2 = _sx-1; clipped=1; iDir=0; iFace=1; } 
	  
	  if(options & OPT_IPBC_GENERIC )
	  {
	    Data_T val;
	    if( clipped )
	    {
	      val = applyDomainBC( ix, iy, iz, 
		  	           ix2, iy2, iz2, iDir, iFace );
	    }
	    else 
	    {
	      //UT_ASSERT2(ix2==ix0+i && iy2==iy0+j && iz2==iz0+k);
	      val = getValue(ix2,iy2,iz2);
	    }
	    buf_x[i] = val;
	  }
	  else
	  {
	    const Data_T *valuePtr = 
	     (options & OPT_IPBC_DIRICHLET) && clipped ? &_emptyValue :
		       getValuePtr(ix2,iy2,iz2);
	    if(!valuePtr)
	    {
	      if(options & OPT_IPCHKBLK)
	      {
		TRCERR(("name = '%s', pos=%g,%g,%g, opt=%d dmax=%d",_name,
		      p0[0],p0[1],p0[2],options,dmax));
		UT_ASSERT0(FALSE);
	      }
	      else valuePtr = &_emptyValue;
	    }
	    buf_x[i] = *valuePtr;
	  }
	}
	buf_y[j] = Interpolate1D<Data_T,DataInt_T,ipType_T>(buf_x,u);
      }
      buf_z[k] = Interpolate1D<Data_T,DataInt_T,ipType_T>(buf_y,v);	       
    }
    f_out = Interpolate1D<Data_T,DataInt_T,ipType_T>(buf_z,w);	       
    return f_out;
  }

  //applyBorderOffset( ix0, iy0, iz0 );

  int bx0, by0, bz0, bid;
  int vx0, vy0, vz0;

  // fast special case: exactly grid aligned 
  if(options & OPT_IPGRIDALGN)
  {
    ix0 += dmax/2;  
    iy0 += dmax/2;  
    iz0 += dmax/2;  
    getBlockCoords(ix0, iy0, iz0, 
		   bx0, by0, bz0);
    bid = getBlockIndex(bx0, by0, bz0);
    getVoxelInBlockCoords(ix0, iy0, iz0, 
			vx0, vy0, vz0 );

    Block<Data_T> *block = getBlock(bid);    
    #ifdef UT_ASSERT_LEVEL_2
    if(!( fabs(u)<MT_NUM_EPS && fabs(v)<MT_NUM_EPS && fabs(w)<MT_NUM_EPS))
    {
      TRCERR(("%g %g %g -> %g %g %g\n",p0[0],p0[1],p0[2],u,v,w));
    }
    if(options&OPT_IPCHKBLK)
    {
      UT_ASSERT2(block);
    }
    #endif
    return block ? block->_data[ getVoxelInBlockIndex(vx0,vy0,vz0) ] : 
		   _emptyValue;
  }

  // Locate the block that contains the origin of the stencil window

  getBlockCoords(ix0, iy0, iz0, 
		 bx0, by0, bz0);
  bid = getBlockIndex(bx0, by0, bz0);
  getVoxelInBlockCoords(ix0, iy0, iz0, 
			vx0, vy0, vz0 );

  if(ipType_T==IP_NEAREST)
  {
    UT_ASSERT0(FALSE); // TODO: correct ?
    Block<Data_T> *block = getBlock(bid);    
    if(doRecordMinMax_T)
    {
      UT_ASSERT0(FALSE);
    }

    #ifdef UT_ASSERT_LEVEL_2
    if(options&OPT_IPCHKBLK)
    {
      UT_ASSERT2(block);
    }
    #endif

    return block ? block->_data[ getVoxelInBlockIndex(vx0,vy0,vz0) ] : 
		   _emptyValue;
  }
  
  // Check for fast case: whole stencil window within block
  if( ! (((vx0+dmax)>>_bsxLog2) |
         ((vy0+dmax)>>_bsxLog2) |
         ((vz0+dmax)>>_bsxLog2)))
  {
    Block<Data_T>*block = getBlock(bid);
    if(!block) 
    {
      if(options&OPT_IPCHKBLK)
      {
	UT_ASSERT0(FALSE);
      }
      // fastest possible case 
      return _emptyValue;  
    }

    // interpolate within this block
    // TODO: SIMD as in G3D_InterpolateCubicFast

    if(doRecordMinMax_T)
    {	    
      *pMin = 1e20f; *pMax = -1e20f;
    }

    if(ipType_T==IP_LINEAR)
    {
      for(int k=0;k<2;k++)
      {
	for(int j=0;j<2;j++)
	{
	  Data_T *data = block->_data + 
		  getVoxelInBlockIndex(vx0+0,vy0+j,vz0+k);
	  buf_x[0] = data[0];
	  buf_x[1] = data[1];
	  if(doRecordMinMax_T)
	  {	    
	    *pMin=MING(*pMin, buf_x[0]); *pMax=MAXG(*pMax, buf_x[0]); 	    	    
	    *pMin=MING(*pMin, buf_x[1]); *pMax=MAXG(*pMax, buf_x[1]); 	    	    
	  }	  
	  buf_y[j] = Interpolate1D<Data_T,DataInt_T,ipType_T>(buf_x,u);
	}
	buf_z[k] = Interpolate1D<Data_T,DataInt_T,ipType_T>(buf_y,v);	       
      }
      f_out = Interpolate1D<Data_T,DataInt_T,ipType_T>(buf_z,w);
    }
    else if(ipType_T==IP_CUBIC || 
	    ipType_T==IP_CUBIC_MONO||
	    ipType_T==IP_CUBIC_MONO_2 ||
	    ipType_T==IP_WENO4) 
    {
      for(int k=0;k<4;k++)
      {
	for(int j=0;j<4;j++)
	{
	  Data_T *data = block->_data + 
		  getVoxelInBlockIndex(vx0+0,vy0+j,vz0+k);
	  buf_x[0] = data[0];
	  buf_x[1] = data[1];
	  buf_x[2] = data[2];
	  buf_x[3] = data[3];
	  if(doRecordMinMax_T)
	  {	    
	    *pMin=MING(*pMin, buf_x[0]); *pMax=MAXG(*pMax, buf_x[0]); 	    	    
	    *pMin=MING(*pMin, buf_x[1]); *pMax=MAXG(*pMax, buf_x[1]); 	    	    
	    *pMin=MING(*pMin, buf_x[2]); *pMax=MAXG(*pMax, buf_x[2]); 	    	    
	    *pMin=MING(*pMin, buf_x[3]); *pMax=MAXG(*pMax, buf_x[3]); 	    	    
	  }	  
	  buf_y[j] = Interpolate1D<Data_T,DataInt_T,ipType_T>(buf_x,u);
	}
	buf_z[k] = Interpolate1D<Data_T,DataInt_T,ipType_T>(buf_y,v);	       
      }
      f_out = Interpolate1D<Data_T,DataInt_T,ipType_T>(buf_z,w);
    }
  }
  else  // not within single block
  {
    // 
    // Take into account the neighboring blocks 
    // (A total maximum of 2x2x2 blocks)
    //    
    if(doRecordMinMax_T)
    {	    
      *pMin = 1e20f; *pMax = -1e20f;
    }

    if(ipType_T==IP_LINEAR)
    {
      for(int k=0;k<2;k++)
      {
	int ibz = ( bz0 + ((vz0+k) >> _bsxLog2) ) *_nxy,
	    ivz = (((vz0+k) & _bsxMask) << _bsx2Log2);
	for(int j=0;j<2;j++)
	{
	  int iby = ibz + (by0 + ((vy0+j) >> _bsxLog2)) *_nx,
	      ivy = ivz + (((vy0+j) & _bsxMask) <<_bsxLog2);
	  for(int i=0;i<2;i++)
	  {
	    int ib = iby + bx0 + ((vx0+i) >> _bsxLog2),
		iv = ivy + ((vx0+i) & _bsxMask);
	    UT_ASSERT2(ib>=0&&ib<_nBlocks);
	    
	    Block<Data_T>* block = _blockmap[ib]; 	  
	    UT_ASSERT2(iv>=0&&iv<_nVoxelsInBlock);
	    #ifdef UT_ASSERT_LEVEL_2
	    if(options&OPT_IPCHKBLK)
	    {
	      UT_ASSERT2(block);
	    }
	    #endif

	    buf_x[i] = block ? block->_data[ iv ] : 
			       _emptyValue;
	    if(doRecordMinMax_T)
	    {	    
	      *pMin=MING(*pMin, buf_x[i]); *pMax=MAXG(*pMax, buf_x[i]); 	    	    
	    }	  
	  }
	  buf_y[j] = Interpolate1D<Data_T,DataInt_T,ipType_T>(buf_x,u);
	}
	buf_z[k] = Interpolate1D<Data_T,DataInt_T,ipType_T>(buf_y,v);	       
      }
      f_out = Interpolate1D<Data_T,DataInt_T,ipType_T>(buf_z,w);	       

    }
    else if(ipType_T==IP_CUBIC || 
	    ipType_T==IP_CUBIC_MONO||
	    ipType_T==IP_CUBIC_MONO_2 ||
	    ipType_T==IP_WENO4) 
    {
      for(int k=0;k<4;k++)
      {
	int ibz = ( bz0 + ((vz0+k) >> _bsxLog2) ) *_nxy,
	    ivz = (((vz0+k) & _bsxMask) << _bsx2Log2);
	for(int j=0;j<4;j++)
	{
	  int iby = ibz + (by0 + ((vy0+j) >> _bsxLog2)) *_nx,
	      ivy = ivz + (((vy0+j) & _bsxMask) <<_bsxLog2);
	  for(int i=0;i<4;i++)
	  {
	    int ib = iby + bx0 + ((vx0+i) >> _bsxLog2),
		iv = ivy + ((vx0+i) & _bsxMask);
	    UT_ASSERT2(ib>=0&&ib<_nBlocks);
	    
	    Block<Data_T>* block = _blockmap[ib]; 	  
	    UT_ASSERT2(iv>=0&&iv<_nVoxelsInBlock);
	    #ifdef UT_ASSERT_LEVEL_2
	    if(options&OPT_IPCHKBLK)
	    {
#if 1
	      if(!block)
	      {
		TRCERR(("%g %g %g -> %d %d %d, b0=%d,%d,%d, v0=%d,%d,%d\n",
		      p0[0],p0[1],p0[2],ix0,iy0,iz0,bx0,by0,bz0,vx0,vy0,vz0));
	      }
#else
	      UT_ASSERT2(block);
#endif
	    }
	    #endif

	    buf_x[i] = block ? block->_data[ iv ] : 
			       _emptyValue;

	    if(doRecordMinMax_T)
	    {	    
	      *pMin=MING(*pMin, buf_x[i]); *pMax=MAXG(*pMax, buf_x[i]); 	    	    
	    }	  
	  }
	  buf_y[j] = Interpolate1D<Data_T,DataInt_T,ipType_T>(buf_x,u);
	}
	buf_z[k] = Interpolate1D<Data_T,DataInt_T,ipType_T>(buf_y,v);	       
      }
      f_out = Interpolate1D<Data_T,DataInt_T,ipType_T>(buf_z,w);	       
    }
  }

  return f_out;
}

/*-------------------------------------------------------------------------*/
/*									   */
/*-------------------------------------------------------------------------*/
template<>
template <int doRetValue_T, 
	  int doRetGradient_T>
void SparseGrid<float>::interpolateLinearWithGradient( 
    			    const Vec4f& p0, unsigned options,
			    float *pOutVal,
			    float *pOutGrad ) const
{
  UT_ASSERT2( options & OPT_IPCORNER );
  Vec4f pos = min(max(p0,_v4fDomMin),_v4fDomMax-1.0f),
	posf = floor(pos);

  Vec4i ipos = truncate_to_int(posf);

  UT_ASSERT2( inRange(ipos) );

  Vec4i bpos = getBlockCoords( ipos ),
	vpos = getVoxelInBlockCoords( ipos );
  int bid = getBlockIndex( vget_x(bpos), vget_y(bpos), vget_z(bpos));
  Vec4f	uvw = pos-posf;
  int vx0=vget_x(vpos),
      vy0=vget_y(vpos),
      vz0=vget_z(vpos);

  int vxMax=_bsx-1,
      vid = getVoxelInBlockIndex(vx0,vy0,vz0);

  float C[8];

  UT_ASSERT2( vget_w(vpos)==0 ); 
  if( !vany( vpos==vxMax ) ) 
  {
    // Fast case: all within single block
    float* p = getBlockDataPtr0(bid)+vid;
    C[0]= p[0];
    C[1]= p[1];
    C[2]= p[_bsx];
    C[3]= p[_bsx + 1];
    C[4]= p[_bsx2];
    C[5]= p[_bsx2 + 1];
    C[6]= p[_bsx2 + _bsx];
    C[7]= p[_bsx2 + _bsx + 1];
  }
  else
  {
    int dx,dbx,
	dy,dby,
	dz,dbz;

    if(vx0<vxMax)
    {
      dx = 1; dbx = 0;
    }
    else
    {
      dx = -vxMax; dbx = 1;
    }

    if(vy0<vxMax)
    {
      dy = _bsx; dby = 0;
    }
    else
    {  
      dy = -vxMax*_bsx; dby = _nx;
    }

    if(vz0<vxMax)
    {
      dz = _bsx2; dbz = 0;
    }
    else
    {
      dz = -vxMax*_bsx2; dbz = _nxy;
    }

    C[0]=getBlockDataPtr0(	bid+0)[		vid+0]; 
    C[1]=getBlockDataPtr0(	bid+dbx)[	vid+dx]; 
    C[2]=getBlockDataPtr0(	bid+dby+0)[	vid+dy+0]; 
    C[3]=getBlockDataPtr0(	bid+dby+dbx)[	vid+dy+dx]; 
    C[4]=getBlockDataPtr0(	bid+dbz+0)[	vid+dz+0]; 
    C[5]=getBlockDataPtr0(	bid+dbz+dbx)[	vid+dz+dx]; 
    C[6]=getBlockDataPtr0(	bid+dbz+dby+0)[	vid+dz+dy+0]; 
    C[7]=getBlockDataPtr0(	bid+dbz+dby+dbx)[vid+dz+dy+dx]; 
  }

  float U=vfget_x(uvw), V=vfget_y(uvw), W=vfget_z(uvw);

  MtTrilinearInterpolationOld( C, U, V, W,
      			    doRetValue_T ? pOutVal : NULL,
			    doRetGradient_T ? pOutGrad : NULL );
}

/*-------------------------------------------------------------------------*/
/*									   */
/*-------------------------------------------------------------------------*/
template< typename Data_T >
template< int options >
void SparseGrid<Data_T>::interpolateWithSecondDerivs( 
    const Vec4f& p0,
    unsigned opt,
    float *out_f, // OUT: function value
    float *out_grad, // OUT: gradient [fx fy fz]
    float *out_hess // OUT: 2nd derivs [fxx fyy fzz fxy fxz fyz]
    ) const
{
  //UT_ASSERT0(FALSE); // TODO: Check correctness
  
  // based on cubic spline in 4x4x4 window
  constexpr int WINSZ=4;
  Data_T win[WINSZ*WINSZ*WINSZ];

  // dmax = 'stencil window size' - 1 
  Vec4f pos,t;
  Vec4i ipos;

  const int dmax = 3;
  pos = p0-0.5f;   
  Vec4f posf = floor(pos);
  t = pos-posf;
  ipos = truncate_to_int(posf) - dmax/2;

  int ix0 = vget_x(ipos),
      iy0 = vget_y(ipos),
      iz0 = vget_z(ipos);
  
  // Fetch local neighborhood window

  if( !( opt & OPT_IP_NO_BORDER_CHECK) && 
      (! ( ix0>=0 && iy0>=0 && iz0>=0 && 
	ix0+dmax<_sx && iy0+dmax<_sy && iz0+dmax<_sz ))) 
  {
    // Border cases
    int iDir=-1,iFace=-1;
    for(int k=0, iWin=0; k<=dmax; k++ )
    {
      int iz = iz0+k,
	  iz2 = iz,
	  clippedZ = 0;
      if(iz2<0) {    iz2=0; 	clippedZ=1; iDir=2; iFace=0; } 
      if(iz2>=_sz) { iz2=_sz-1; clippedZ=1; iDir=2; iFace=1; }

      for(int j=0;j<=dmax;j++)
      {
	int iy = iy0+j,
	    iy2 = iy,
	    clippedY = clippedZ;
	if(iy2<0) { 	iy2=0;     clippedY=1; iDir=1; iFace=0; } 
	if(iy2>=_sy) {  iy2=_sy-1; clippedY=1; iDir=1; iFace=1; } 

	for(int i=0;i<=dmax;i++, iWin++)
	{
	  int ix = ix0+i,
	      ix2 = ix,
	      clipped = clippedY;
	  if(ix2<0) {    ix2 = 0;     clipped=1; iDir=0; iFace=0; } 
	  if(ix2>=_sx) { ix2 = _sx-1; clipped=1; iDir=0; iFace=1; } 
	  
	  if(opt & OPT_IPBC_GENERIC )
	  {
	    Data_T val;
	    if( clipped )
	    {
	      val = applyDomainBC( ix, iy, iz, 
		  	           ix2, iy2, iz2, iDir, iFace );
	    }
	    else 
	    {
	      //UT_ASSERT2(ix2==ix0+i && iy2==iy0+j && iz2==iz0+k);
	      val = getValue(ix2,iy2,iz2);
	    }
	    win[iWin] = val;
	  }
	  else
	  {
	    const Data_T *valuePtr = 
	     (opt & OPT_IPBC_DIRICHLET) && clipped ? &_emptyValue :
		       getValuePtr(ix2,iy2,iz2);
	    if(!valuePtr)
	    {
	      if(opt & OPT_IPCHKBLK)
	      {
		UT_ASSERT0(FALSE);
	      }
	      else valuePtr = &_emptyValue;
	    }
	    win[iWin] = *valuePtr;
	  }
	}
      }
    }
  }
  else
  {
    // 
    // Non domain-border Case
    //
    int bx0, by0, bz0, bid;
    int vx0, vy0, vz0;
    getBlockCoords(ix0, iy0, iz0, 
		   bx0, by0, bz0);
    bid = getBlockIndex(bx0, by0, bz0);
    getVoxelInBlockCoords(ix0, iy0, iz0, 
			  vx0, vy0, vz0 );

    UT_ASSERT2(_maxBordOff==0);

    if( ! (((vx0+dmax)>>_bsxLog2) |
	   ((vy0+dmax)>>_bsxLog2) |
	   ((vz0+dmax)>>_bsxLog2)))
    {
      // 
      // Fast case: whole stencil window within block
      //
      Block<Data_T>*block = getBlock(bid);
      if(!block) 
      {
	if(opt&OPT_IPCHKBLK)
	{
	  UT_ASSERT0(FALSE);
	}
	for(int i=0; i<(dmax+1)*(dmax+1)*(dmax+1); i++) win[i] = _emptyValue;
      }
      else
      {
	for(int k=0, iWin=0; k<=dmax;k++)
	for(int j=0;j<=dmax;j++)
	{
	  int idxLine = getVoxelInBlockIndex(vx0+0,vy0+j,vz0+k);
	  for(int i=0;i<=dmax;i++) win[iWin++] = block->_data[idxLine+i];
	}
      }
    }
    else
    {
      // Not within single block: Take into account the neighboring blocks 
      // (A total maximum of 2x2x2 blocks)
      for(int k=0, iWin=0; k<=dmax; k++ )
      {
	int ibz = ( bz0 + ((vz0+k) >> _bsxLog2) ) *_nxy,
	    ivz = (((vz0+k) & _bsxMask) << _bsx2Log2);
	for(int j=0;j<=dmax;j++)
	{
	  int iby = ibz + (by0 + ((vy0+j) >> _bsxLog2)) *_nx,
	      ivy = ivz + (((vy0+j) & _bsxMask) <<_bsxLog2);
	  for(int i=0;i<=dmax;i++)
	  {
	    int ib = iby + bx0 + ((vx0+i) >> _bsxLog2),
		iv = ivy + ((vx0+i) & _bsxMask);
	    UT_ASSERT2(ib>=0&&ib<_nBlocks);
	      
	    Block<Data_T>* block = _blockmap[ib]; 	  
	    UT_ASSERT2(iv>=0&&iv<_nVoxelsInBlock);
	    #ifdef UT_ASSERT_LEVEL_2
	    if(opt&OPT_IPCHKBLK)
	    {
	      UT_ASSERT2(block);
	    }
	    #endif

	    {
	      float f;
	      if(block)
	      {
		f = block->_data[ iv ];
	      }
	      else
	      {
		#if 0
		// Neumann
		int x=ix0+i,
		    y=iy0+j,
		    z=iz0+k;
    		clipBox( bx0<<_bsxLog2, by0<<_bsxLog2, bz0<<_bsxLog2, 
		    	 _bsx-1, _bsx-1, _bsx-1,
			 x, y, z );
	        Block<Data_T>* block = _blockmap[bid]; 	  
		int iv = getVoxelInBlockIndex( 
		    		getVoxelInBlockCoords( Vec4i(x,y,z,0) ) );
		f = block->_data[iv];
		#else
		// Dirichlet
		f = _emptyValue;  
		#endif
	      }
	      win[iWin++] = f;	     
	    }
	  }
	}
      }
    }
  }

#if 0
  // Smooth the window
  {
    Data_T win2[WINSZ*WINSZ*WINSZ];
    for(int i=0;i<ARRAY_LENGTH(win2);i++) win2[i] = win[i];

    float filtRadius = 2.0f,
	  filtRadiusInvSq = 1.0f/(filtRadius*filtRadius);

    for (int z = 0; z <WINSZ; z++)
    for (int y = 0; y <WINSZ; y++)
    for (int x = 0; x <WINSZ; x++)
    {
      float sum=0,wsum=0;
      for(int dx=-1;dx<=1;dx++)
      for(int dy=-1;dy<=1;dy++)
      for(int dz=-1;dz<=1;dz++)
      {
	int x2=x+dx,
	    y2=y+dy,
	    z2=z+dz;

        #if 1
	if( !G3D_IN_BOX(0,0,0, WINSZ-1,WINSZ-1,WINSZ-1, x2,y2,z2) ) continue;
        #else
	clipBox( 0, 0, 0, WINSZ-1, WINSZ-1, WINSZ-1,
		 x2, y2, z2 );
        #endif

	int idx2 = MT_GXYZ(WINSZ,WINSZ*WINSZ, x2, y2, z2 );
	UT_ASSERT2(idx2>=0&&idx2<WINSZ*WINSZ*WINSZ);
	float f = win2[idx2];     
	float w = MAX(0,1.0f-(dx*dx+dy*dy+dz*dz)*filtRadiusInvSq);
	w = w*w*w;
	//TRCP(("%d,%d,%d -> %g\n",dx,dy,dz,w));
	wsum += w;
	sum += w*f;
      }
      int idx = MT_GXYZ(WINSZ,WINSZ*WINSZ, x, y, z );
      UT_ASSERT2(idx>=0&&idx<WINSZ*WINSZ*WINSZ);
      float f = sum * (1.0f/wsum);
      win[idx] = f;
    }  
  }
#endif


  //
  // Perform calculation on fetched neighborhood window
  //
  float wx[4],wy[4],wz[4],
	wx_[4],wy_[4],wz_[4],
	wx__[4],wy__[4],wz__[4];

  MtCubicBSplineWeightsForSecondDerivs1D( vfget_x(t), wx, wx_, wx__); 
  MtCubicBSplineWeightsForSecondDerivs1D( vfget_y(t), wy, wy_, wy__); 
  MtCubicBSplineWeightsForSecondDerivs1D( vfget_z(t), wz, wz_, wz__); 

  Vec4d grad(0.f);
  Vec8f hess(0.0f);

  UT_ASSERT2(sizeof(Data_T)>=4); // TODO: support 16-bit 

  for (int dz = 0, iWin=0; dz <=dmax; dz++)
  {
    Vec4d wgz_( wz[dz], wz[dz], wz[dz], wz_[dz] );
    Vec8f wgz__( wz[dz], wz[dz], wz__[dz], wz_[dz], wz_[dz], wz[dz], 0.f, 0.f ); 

    for (int dy = 0; dy <=dmax; dy++)
    {
      Vec4d wgy_ = wgz_ * Vec4d( wy[dy], wy[dy], wy_[dy], wy[dy] );
      Vec8f wgy__ = wgz__ * Vec8f( wy[dy], wy__[dy], wy[dy], wy_[dy], wy[dy], wy_[dy], 0.f, 0.f );

      for (int dx = 0; dx <=dmax; dx++)
      {      
	Vec4d w_ = wgy_ * Vec4d( wx[dx], wx_[dx],  wx[dx],  wx[dx] );
	Vec8f w__ = wgy__ * Vec8f( wx__[dx], wx[dx], wx[dx], wx[dx], wx_[dx], wx_[dx], 0.f, 0.f );

	float f = /*renderDensToFloat_*/ ( win[iWin++] );

	grad += w_*f;
	hess += w__*f;
      }
    }
  }

  {
    double tmp[8];
    grad.store(tmp);

    *out_f = tmp[0];

    out_grad[0] = tmp[1];
    out_grad[1] = tmp[2];
    out_grad[2] = tmp[3];
  }

  {
    float tmp[8];
    hess.store(tmp);

    out_hess[0] = tmp[0];
    out_hess[1] = tmp[1];
    out_hess[2] = tmp[2];

    out_hess[3] = tmp[5];
    out_hess[4] = tmp[4];
    out_hess[5] = tmp[3];
  }
}


/*-------------------------------------------------------------------------*/
/*									   */
/*-------------------------------------------------------------------------*/
template< typename Data_T >
template < typename DataRet_T,
	   typename DataInt_T,
	   int dummy_T >
void SparseGrid<Data_T>::interpolateWithDerivs0( 
    				int ipType,
    				const Vec4f& p0,
				unsigned opt,
				float *outVal,
				Vec4f *outGrad
				) const
{
  DataInt_T win[4*4*4];
  DataRet_T outVal0;
  Vec4f outGrad0;
  if(!outVal) outVal = &outVal0;
  if(!outGrad) outGrad =&outGrad0;

  // dmax = 'stencil window size' - 1 
  int dmax;
  Vec4f pos,t;
  Vec4i ipos;

  if( ipType == IP_BSPLINE_QUAD )
  {
    dmax = 2;
    pos = p0;
    Vec4f posf = floor(pos);
    t = ceil(p0)-p0; 
    ipos = truncate_to_int(posf) - dmax/2;
  }
  else if( ipType == IP_BSPLINE_CUBIC )
  {
    dmax = 3;
    pos = p0-0.5f;
    Vec4f posf = floor(pos);
    t = pos-posf;
    ipos = truncate_to_int(posf) - dmax/2;
  }
  else
  {
    UT_ASSERT0(FALSE);
  }


  int ix0 = vget_x(ipos),
      iy0 = vget_y(ipos),
      iz0 = vget_z(ipos);
  
  // Fetch local neighborhood window

  if( !( opt & OPT_IP_NO_BORDER_CHECK) && 
      (! ( ix0>=0 && iy0>=0 && iz0>=0 && 
	ix0+dmax<_sx && iy0+dmax<_sy && iz0+dmax<_sz ))) 
  {
    // Border cases
    int iDir=-1,iFace=-1;
    for(int k=0, iWin=0; k<=dmax; k++ )
    {
      int iz = iz0+k,
	  iz2 = iz,
	  clippedZ = 0;
      if(iz2<0) {    iz2=0; 	clippedZ=1; iDir=2; iFace=0; } 
      if(iz2>=_sz) { iz2=_sz-1; clippedZ=1; iDir=2; iFace=1; }

      for(int j=0;j<=dmax;j++)
      {
	int iy = iy0+j,
	    iy2 = iy,
	    clippedY = clippedZ;
	if(iy2<0) { 	iy2=0;     clippedY=1; iDir=1; iFace=0; } 
	if(iy2>=_sy) {  iy2=_sy-1; clippedY=1; iDir=1; iFace=1; } 

	for(int i=0;i<=dmax;i++, iWin++)
	{
	  int ix = ix0+i,
	      ix2 = ix,
	      clipped = clippedY;
	  if(ix2<0) {    ix2 = 0;     clipped=1; iDir=0; iFace=0; } 
	  if(ix2>=_sx) { ix2 = _sx-1; clipped=1; iDir=0; iFace=1; } 
	  
	  if(opt & OPT_IPBC_GENERIC )
	  {
	    Data_T val;
	    if( clipped )
	    {
	      val = applyDomainBC( ix, iy, iz, 
		  	           ix2, iy2, iz2, iDir, iFace );
	    }
	    else 
	    {
	      //UT_ASSERT2(ix2==ix0+i && iy2==iy0+j && iz2==iz0+k);
	      val = getValue(ix2,iy2,iz2);
	    }
	    win[iWin] = val;
	  }
	  else
	  {
	    const Data_T *valuePtr = 
	     (opt & OPT_IPBC_DIRICHLET) && clipped ? &_emptyValue :
		       getValuePtr(ix2,iy2,iz2);
	    if(!valuePtr)
	    {
	      if(opt & OPT_IPCHKBLK)
	      {
		UT_ASSERT0(FALSE);
	      }
	      else valuePtr = &_emptyValue;
	    }
	    win[iWin] = *valuePtr;
	  }
	}
      }
    }
  }
  else
  {
    // 
    // Non domain-border Case
    //
    int bx0, by0, bz0, bid;
    int vx0, vy0, vz0;
    getBlockCoords(ix0, iy0, iz0, 
		   bx0, by0, bz0);
    bid = getBlockIndex(bx0, by0, bz0);
    getVoxelInBlockCoords(ix0, iy0, iz0, 
			  vx0, vy0, vz0 );

    UT_ASSERT2(_maxBordOff==0);

    if( ! (((vx0+dmax)>>_bsxLog2) |
	   ((vy0+dmax)>>_bsxLog2) |
	   ((vz0+dmax)>>_bsxLog2)))
    {
      // 
      // Fast case: whole stencil window within block
      //
      Block<Data_T>*block = getBlock(bid);
      if(!block) 
      {
	if(opt&OPT_IPCHKBLK)
	{
	  UT_ASSERT0(FALSE);
	}
	for(int i=0; i<(dmax+1)*(dmax+1)*(dmax+1); i++) win[i] = _emptyValue;
      }
      else
      {
	for(int k=0, iWin=0; k<=dmax;k++)
	for(int j=0;j<=dmax;j++)
	{
	  int idxLine = getVoxelInBlockIndex(vx0+0,vy0+j,vz0+k);
	  for(int i=0;i<=dmax;i++) win[iWin++] = block->_data[idxLine+i];
	}
      }
    }
    else
    {
      // Not within single block: Take into account the neighboring blocks 
      // (A total maximum of 2x2x2 blocks)
      for(int k=0, iWin=0; k<=dmax; k++ )
      {
	int ibz = ( bz0 + ((vz0+k) >> _bsxLog2) ) *_nxy,
	    ivz = (((vz0+k) & _bsxMask) << _bsx2Log2);
	for(int j=0;j<=dmax;j++)
	{
	  int iby = ibz + (by0 + ((vy0+j) >> _bsxLog2)) *_nx,
	      ivy = ivz + (((vy0+j) & _bsxMask) <<_bsxLog2);
	  for(int i=0;i<=dmax;i++)
	  {
	    int ib = iby + bx0 + ((vx0+i) >> _bsxLog2),
		iv = ivy + ((vx0+i) & _bsxMask);
	    UT_ASSERT2(ib>=0&&ib<_nBlocks);
	      
	    Block<Data_T>* block = _blockmap[ib]; 	  
	    UT_ASSERT2(iv>=0&&iv<_nVoxelsInBlock);
	    #ifdef UT_ASSERT_LEVEL_2
	    if(opt&OPT_IPCHKBLK)
	    {
	      UT_ASSERT2(block);
	    }
	    #endif
	    win[iWin++] = block ? block->_data[ iv ] : 
				 _emptyValue;
	  }
	}
      }
    }
  }

  //
  // Perform calculation on fetched neighborhood window
  //
  float wx[4],wy[4],wz[4],
	wx_[4],wy_[4],wz_[4];

  if(ipType == IP_BSPLINE_QUAD)
  {
    MtQuadBSplineWeights1D( vfget_x(t), wx, wx_); 
    MtQuadBSplineWeights1D( vfget_y(t), wy, wy_); 
    MtQuadBSplineWeights1D( vfget_z(t), wz, wz_); 
    wx[3]=wx_[3]=wy[3]=wy_[3]=wz[3]=wz_[3]=0.0f;
  }
  else if(ipType == IP_BSPLINE_CUBIC )
  {
    MtCubicBSplineWeights1D( vfget_x(t), wx, wx_); 
    MtCubicBSplineWeights1D( vfget_y(t), wy, wy_); 
    MtCubicBSplineWeights1D( vfget_z(t), wz, wz_); 
  }
  else
  {
    UT_ASSERT0(FALSE);
  }

  if(outVal) *outVal=0;
  if(outGrad) *outGrad=0;

  for (int dz = 0, iWin=0; dz <=dmax; dz++)
  {
    float wgz = wz[dz];
    Vec4f wgz_( wz[dz], wz[dz], wz_[dz], 0.0f );

    for (int dy = 0; dy <=dmax; dy++)
    {
      float wgy = wgz * wy[dy];
      Vec4f wgy_ = wgz_ * Vec4f( wy[dy], wy_[dy], wy[dy], 0.0f );

      for (int dx = 0; dx <=dmax; dx++)
      {      
	float w = wgy * wx[dx];
	Vec4f w_ = wgy_ * Vec4f( wx_[dx], wx[dx], wx[dx], 0.0f );

	float f = win[iWin++];
	*outVal += w * f;
	*outGrad += w_*f;
      }
    }
  }

}

/*-------------------------------------------------------------------------*/
/*									   */
/*-------------------------------------------------------------------------*/
template<>
template <int ipType_T> 
Vec3Float SparseGrid<Vec3Float>::interpolateVec3Float(     		    
	      int dmax,   // stenWidth-1
	      float *p0,	// Position
	      int iVelComp,	// -1= cell-centerd, >=0: staggered/MAC
	      int options
	      ) const
{
  if( ipType_T!=IP_LINEAR || !(options & OPT_IP_NO_BORDER_CHECK) )
    return interpolate<ipType_T>(p0,options);
    
  Vec3Float out(0.0f);
	
  int bx0, by0, bz0, bid;
  int vx0, vy0, vz0;

  // fast special case: exactly grid aligned 
  if(options & OPT_IPGRIDALGN)
  {
    int ix0 = floorf(p0[0]),
	iy0 = floorf(p0[1]),
	iz0 = floorf(p0[2]);
    // Border-check already done by caller (OPT_IP_NO_BORDER_CHECK)
    UT_ASSERT2(( ix0>=0&&iy0>=0&&iz0>=0 && ix0<_sx&&iy0<_sy&&iz0<_sz )); 
    getBlockCoords(ix0, iy0, iz0, 
		   bx0, by0, bz0);
    bid = getBlockIndex(bx0, by0, bz0);
    getVoxelInBlockCoords(ix0, iy0, iz0, 
			vx0, vy0, vz0 );
    Block<Vec3Float> *block = getBlock(bid);    
    #ifdef UT_ASSERT_LEVEL_2
    {
      float u = fabs(p0[0]-.5f-floorf(p0[0]-.5f)),
	    v = fabs(p0[1]-.5f-floorf(p0[1]-.5f)),
	    w = fabs(p0[2]-.5f-floorf(p0[2]-.5f));
      if(!( u<MT_NUM_EPS && v<MT_NUM_EPS && w<MT_NUM_EPS))
      {
	TRCERR(("%g %g %g -> %g %g %g\n",p0[0],p0[1],p0[2],u,v,w));
	UT_ASSERT0(FALSE);
      }
      if(options&OPT_IPCHKBLK)
      {
	UT_ASSERT2(block);
      }
    }
    #endif
    return block ? block->_data[ getVoxelInBlockIndex(vx0,vy0,vz0) ] : 
		   _emptyValue;
  }

  float  x = p0[0]-0.5f, 	// po=0.5 -> interpolate at pixel centers
	 y = p0[1]-0.5f, 
	 z = p0[2]-0.5f,
  	
  	 xf = floorf(x),
	 yf = floorf(y),
	 zf = floorf(z);

  float u = x - xf,
	v = y - yf,
	w = z - zf;

  if( options & OPT_IP_FADE_WEIGHTS )
  {
    fadeIpolWeights(u,v,w);
  }

  UT_ASSERT2(_maxBordOff==0);

  int ix0 = (int)(xf) - dmax/2,
      iy0 = (int)(yf) - dmax/2,  
      iz0 = (int)(zf) - dmax/2;  

  // Border-check already done by caller (OPT_IP_NO_BORDER_CHECK)
  UT_ASSERT2(( ix0>=0 && iy0>=0 && iz0>=0 && 
	ix0+dmax<_sx && iy0+dmax<_sy && iz0+dmax<_sz )); 

  getBlockCoords(ix0, iy0, iz0, 
		 bx0, by0, bz0);
  bid = getBlockIndex(bx0, by0, bz0);
  getVoxelInBlockCoords(ix0, iy0, iz0, 
			vx0, vy0, vz0 );

  if(ipType_T==IP_LINEAR)
  {
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
	  //

    int vxMax=_bsx-1,
	vid = getVoxelInBlockIndex(vx0,vy0,vz0);

    Vec4f C0,C1,C2,C3,C4,C5,C6,C7;

    if(vx0<vxMax && vy0<vxMax && vz0<vxMax) 
    {
      // Fast case: all within single block
      Vec3Float* p = getBlockDataPtr0(bid)+vid;
      C0.load( p[0] );
      C1.load( p[1]);
      C2.load( p[_bsx] );
      C3.load( p[_bsx + 1]);
      C4.load( p[_bsx2] );
      C5.load( p[_bsx2 + 1]);
      C6.load( p[_bsx2 + _bsx] );
      C7.load( p[_bsx2 + _bsx + 1]);
    }
    else
    {
      int dx,dbx,
	  dy,dby,
	  dz,dbz;

      if(vx0<vxMax)
      {
	dx = 1; dbx = 0;
      }
      else
      {
	dx = -vxMax; dbx = 1;
      }

      if(vy0<vxMax)
      {
	dy = _bsx; dby = 0;
      }
      else
      {
	dy = -vxMax*_bsx; dby = _nx;
      }

      if(vz0<vxMax)
      {
	dz = _bsx2; dbz = 0;
      }
      else
      {
	dz = -vxMax*_bsx2; dbz = _nxy;
      }

      C0.load(getBlockDataPtr0(	bid+0)[		vid+0]); 
      C1.load(getBlockDataPtr0(	bid+dbx)[	vid+dx]); 
      C2.load(getBlockDataPtr0(	bid+dby+0)[	vid+dy+0]); 
      C3.load(getBlockDataPtr0(	bid+dby+dbx)[	vid+dy+dx]); 

      C4.load(getBlockDataPtr0(	bid+dbz+0)[	vid+dz+0]); 
      C5.load(getBlockDataPtr0(	bid+dbz+dbx)[	vid+dz+dx]); 
      C6.load(getBlockDataPtr0(	bid+dbz+dby+0)[	vid+dz+dy+0]); 
      C7.load(getBlockDataPtr0(	bid+dbz+dby+dbx)[vid+dz+dy+dx]); 
    }

    Vec4f U(u),V(v),W(w),
      	  U1 = 1.0f-U,
	  V1 = 1.0f-V,
	  W1 = 1.0f-W;

    Vec4f out_ = 
	( (   U1 * C0 
	+      U * C1  ) * V1
	+ (   U1 * C2 
	+      U * C3  ) * V ) * W1
	+ ( ( U1 * C4 
	+      U * C5  ) * V1
	+ (   U1 * C6 
	+      U * C7  ) * V ) * W ;

    out.x = vfget_x(out_);
    out.y = vfget_y(out_);
    out.z = vfget_z(out_);
  }
  else
  {
    UT_ASSERT0(FALSE);
  }

  return out;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template<> inline
float SparseGrid<uint16_t>::interpolateShortToFloatFast( 
    				const Vec4f& p0, int opt )
const
{
  Vec4f pos = p0-.5f;
  pos = v4f_zero_to<3>(pos);
  Vec4f posf = floor(pos);

  Vec4i ipos = truncate_to_int(posf);
    int ix0=vget_x(ipos),
	iy0=vget_y(ipos),
	iz0=vget_z(ipos);

  if(!(opt & OPT_IP_NO_BORDER_CHECK))
  {
    const int dmax=1;
    if( ! ( ix0>=0 && iy0>=0 && iz0>=0 && 
	    ix0+dmax<_sx && iy0+dmax<_sy && iz0+dmax<_sz ))
    {
      // Border cases handled in generic interpolation routine
      float p0_[4]; p0.store(p0_);
      return interpolateToFloat<IP_LINEAR>(p0_,opt);
    }
  }

  Vec4i bpos = getBlockCoords( ipos ),
	vpos = getVoxelInBlockCoords( ipos );

  int bid = getBlockIndex( vget_x(bpos), vget_y(bpos), vget_z(bpos));

  Vec4f	uvw = pos-posf;

  int vx0=vget_x(vpos),
      vy0=vget_y(vpos),
      vz0=vget_z(vpos);

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
	//

  int vxMax=_bsx-1,
      vid = getVoxelInBlockIndex(vx0,vy0,vz0);
  float C0,C1,C2,C3,C4,C5,C6,C7;

  UT_ASSERT2( vget_w(vpos)==0 ); 
  if( !vany( vpos==vxMax ) ) 
  {
    Block<uint16_t> *b;
    if( (b = getBlock(bid))!=NULL )
    {    
      uint16_t *p=b->_data+vid;
      C0 = p[0];
      C1 = p[1];
      C2 = p[_bsx];
      C3 = p[_bsx + 1];
      C4 = p[_bsx2];
      C5 = p[_bsx2 + 1];
      C6 = p[_bsx2 + _bsx];
      C7 = p[_bsx2 + _bsx + 1];
    }
    else
    { 
      return _emptyValue;
    }
  }
  else
  {
    int dx,dbx,
	dy,dby,
	dz,dbz;

    if(vx0<vxMax)
    {
      dx = 1; dbx = 0;
    }
    else
    {
      dx = -vxMax; dbx = 1;
    }

    if(vy0<vxMax)
    {
      dy = _bsx; dby = 0;
    }
    else
    {
      dy = -vxMax*_bsx; dby = _nx;
    }

    if(vz0<vxMax)
    {
      dz = _bsx2; dbz = 0;
    }
    else
    {
      dz = -vxMax*_bsx2; dbz = _nxy;
    }

    Block<uint16_t> *b;

    b =getBlock( bid+0 ); 
    C0 = b ? b->_data[vid+0 ] : _emptyValue;

    b =getBlock( bid+dbx ); 
    C1 = b ? b->_data[ vid+dx ] : _emptyValue;
    
    b =getBlock( bid+dby+0); 
    C2 = b ? b->_data[	vid+dy+0 ] : _emptyValue;

    b =getBlock( bid+dby+dbx); 
    C3 = b ? b->_data[vid+dy+dx ] : _emptyValue;

    b =getBlock( bid+dbz+0 ); 
    C4 = b ? b->_data[ vid+dz+0 ] : _emptyValue;

    b =getBlock( bid+dbz+dbx); 
    C5 = b ? b->_data[vid+dz+dx] : _emptyValue;

    b =getBlock( bid+dbz+dby+0 ); 
    C6 = b ? b->_data[vid+dz+dy+0] : _emptyValue;

    b =getBlock( bid+dbz+dby+dbx ); 
    C7 = b ? b->_data[vid+dz+dy+dx] : _emptyValue;
  }

  float U = vfget_x(uvw), V = vfget_y(uvw), W = vfget_z(uvw),
	U1 = 1.0f-U,
	V1 = 1.0f-V,
	W1 = 1.0f-W;

  float out = 
      ( (   U1 * C0 
      +      U * C1  ) * V1
      + (   U1 * C2 
      +      U * C3  ) * V ) * W1
      + ( ( U1 * C4 
      +      U * C5  ) * V1
      + (   U1 * C6 
      +      U * C7  ) * V ) * W ;

  return out;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template<> inline 
Vec4f SparseGrid<Vec3Float>::interpolateVec3FloatFast( 
    				const Vec4f& p0, int opt )
const
{
  Vec4f pos = p0-0.5f,
	posf = floor(pos);

  Vec4i ipos = truncate_to_int(posf);
    int ix0=vget_x(ipos),
	iy0=vget_y(ipos),
	iz0=vget_z(ipos);

  if(!(opt & OPT_IP_NO_BORDER_CHECK))
  {
    const int dmax=1;
    if( ! ( ix0>=0 && iy0>=0 && iz0>=0 && 
	    ix0+dmax<_sx && iy0+dmax<_sy && iz0+dmax<_sz ))
    {
      // Border cases handled in generic interpolation routine
      float p0_[4]; p0.store(p0_);
      Vec3Float f = interpolate<IP_LINEAR>(p0_,opt);
      return Vec4f( f.x, f.y, f.z, 0.0f );
    }
  }

  Vec4i bpos = getBlockCoords( ipos ),
	vpos = getVoxelInBlockCoords( ipos );

  int bid = getBlockIndex( vget_x(bpos), vget_y(bpos), vget_z(bpos));

  Vec4f	uvw = pos-posf;

  int vx0=vget_x(vpos),
      vy0=vget_y(vpos),
      vz0=vget_z(vpos);


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
	//

  int vxMax=_bsx-1,
      vid = getVoxelInBlockIndex(vx0,vy0,vz0);

  Vec4f C0,C1,C2,C3,C4,C5,C6,C7;


  UT_ASSERT2( vget_w(vpos)==0 ); 
  if( !vany( vpos==vxMax ) ) 
  {
    // Fast case: all within single block
    Vec3Float* p = getBlockDataPtr0(bid)+vid;
    C0.load( p[0] );
    C1.load( p[1] );
    C2.load( p[_bsx] );
    C3.load( p[_bsx + 1]);
    C4.load( p[_bsx2] );
    C5.load( p[_bsx2 + 1]);
    C6.load( p[_bsx2 + _bsx] );
    C7.load( p[_bsx2 + _bsx + 1]);
  }
  else
  {
    int dx,dbx,
	dy,dby,
	dz,dbz;

    if(vx0<vxMax)
    {
      dx = 1; dbx = 0;
    }
    else
    {
      dx = -vxMax; dbx = 1;
    }

    if(vy0<vxMax)
    {
      dy = _bsx; dby = 0;
    }
    else
    {
      dy = -vxMax*_bsx; dby = _nx;
    }

    if(vz0<vxMax)
    {
      dz = _bsx2; dbz = 0;
    }
    else
    {
      dz = -vxMax*_bsx2; dbz = _nxy;
    }

    C0.load(getBlockDataPtr0(	bid+0)[		vid+0]); 
    C1.load(getBlockDataPtr0(	bid+dbx)[	vid+dx]); 
    C2.load(getBlockDataPtr0(	bid+dby+0)[	vid+dy+0]); 
    C3.load(getBlockDataPtr0(	bid+dby+dbx)[	vid+dy+dx]); 
    C4.load(getBlockDataPtr0(	bid+dbz+0)[	vid+dz+0]); 
    C5.load(getBlockDataPtr0(	bid+dbz+dbx)[	vid+dz+dx]); 
    C6.load(getBlockDataPtr0(	bid+dbz+dby+0)[	vid+dz+dy+0]); 
    C7.load(getBlockDataPtr0(	bid+dbz+dby+dbx)[vid+dz+dy+dx]); 
  }

  Vec4f U(vfget_x(uvw)), V(vfget_y(uvw)), W(vfget_z(uvw)),
	U1 = 1.0f-U,
	V1 = 1.0f-V,
	W1 = 1.0f-W;

  Vec4f out = 
      ( (   U1 * C0 
      +      U * C1  ) * V1
      + (   U1 * C2 
      +      U * C3  ) * V ) * W1
      + ( ( U1 * C4 
      +      U * C5  ) * V1
      + (   U1 * C6 
      +      U * C7  ) * V ) * W ;

  return v4f_zero_to<3>(out);;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template<> inline
float SparseGrid<float>::interpolateFloatFast( 
    				const Vec4f& p0, int opt,
				IpolAccessCache *iac )
const
{
  constexpr int doUseIAC = true;
  const int bsxLog2 = _bsxLog2;
  //UT_ASSERT2( _bsxLog2 == bsxLog2 );
  const int bsx = 1<<bsxLog2;
  const int bsx2 = bsx*bsx;
  const int vxMax = bsx-1;

  Vec4f pos = p0;
  if(!( opt & OPT_IPCORNER )) pos -= Vec4f(.5f,.5f,.5f,0.f);

  Vec4f posf = floor(pos);

  Vec4i ipos = truncate_to_int(posf);
  
  if(!(opt & OPT_IP_NO_BORDER_CHECK))
  {
    if( vany3( ipos < vzero() || ipos >= _v4iDomMax )) 
    {
      vset_w(ipos,0);
      float p0_[4]; p0.store(p0_);
      return interpolate<IP_LINEAR>(p0_,opt);      
    }
  }

  Vec4f	uvw = pos-posf;

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
	//
  Vec4f g0246,
	g1357;

  if constexpr (doUseIAC)
  {
    if(iac)
    {
      #ifdef SBG_IAC_STATS
      iac->_nAcc++;
      #endif      
      if( vall3( iac->_iposl == ipos ) )
      {
	#ifdef SBG_IAC_STATS
	iac->_nHit++;
	#endif
	g0246 = iac->_g0246;
	g1357 = iac->_g1357;
	goto cache_hit;
      }
    }
  }

  {

  Vec4i bpos = ipos >> bsxLog2;
  int bid = horizontal_add( bpos * _v4iGridStridesBlk );
  Vec4i vpos = ipos & Vec4i(bsx-1,bsx-1,bsx-1,0);    
  int vid = (vget_z(vpos)<<(2*bsxLog2))  | (vget_y(vpos)<<bsxLog2) | vget_x(vpos);
  UT_ASSERT2(bid>=0&&bid<_nBlocks);

  UT_ASSERT2( vget_w(vpos)==0 ); 
  if( !vany( vpos==vxMax ) ) 
  {
    // Fast case: all within single block
    float* p = getBlockDataPtr0(bid)+vid;
    Vec8f tmp = vgather8( p, _v8iStridesLerp02461357 );
    g0246 =  tmp.get_low();
    g1357 =  tmp.get_high();
  }
  else
  {
    int dx,dbx,
	dy,dby,
	dz,dbz;

    Vec4i db = select( vpos < vxMax, 
	               0, Vec4i(1,_nx,_nxy,0) );
    Vec4i d = select( vpos < vxMax, 
	              Vec4i(1,bsx,bsx2,0), -vxMax*Vec4i(1,bsx,bsx2,0) );

     dx = vget_x(d);
     dy = vget_y(d);
     dz = vget_z(d);
     dbx = vget_x(db);
     dby = vget_y(db);
     dbz = vget_z(db);

#if 0
     auto getValueInBlock_dbg = [&]( int bid, int vid ) -> float
     {
      if(!( bid>=0 && bid< _nBlocks && vid>=0 && vid<_bsx*_bsx*_bsx && _blockmap[bid] ))
      {
	TRCERR(("'%s': ipos=%d,%d,%d vpos=%d,%d,%d, dim=%dx%dx%d) bdim=%dx%dx%d _blockmap[bid]=%p\n",
	      _name,
	     ipos[0],ipos[1],ipos[2], 
	     vpos[0],vpos[1],vpos[2], 
	     _sx,_sy,_sz, _nx,_ny,_nz,
	     bid>=0 && bid< _nBlocks && vid>=0 && vid<_bsx*_bsx*_bsx ? _blockmap[bid] : NULL));
      }
      Block<float>* block = getBlock(bid);
      return *getValueInBlockPtr0(block,vid);
    };

    g0246 = Vec4f(
      getValueInBlock_dbg(	bid+0,		vid+0),
      getValueInBlock_dbg(	bid+dby+0,	vid+dy+0), 
      getValueInBlock_dbg(	bid+dbz+0,	vid+dz+0),
      getValueInBlock_dbg(	bid+dbz+dby+0,	vid+dz+dy+0)); 

    g1357 = Vec4f(
      getValueInBlock_dbg(	bid+dbx,	vid+dx), 
      getValueInBlock_dbg(	bid+dby+dbx,	vid+dy+dx), 
      getValueInBlock_dbg(	bid+dbz+dbx,	vid+dz+dx), 
      getValueInBlock_dbg(	bid+dbz+dby+dbx, vid+dz+dy+dx) );
#else
    g0246 = Vec4f(
      getBlockDataPtr0(	bid+0)[		vid+0],
      getBlockDataPtr0(	bid+dby+0)[	vid+dy+0], 
      getBlockDataPtr0(	bid+dbz+0)[	vid+dz+0],
      getBlockDataPtr0(	bid+dbz+dby+0)[	vid+dz+dy+0]); 

    g1357 = Vec4f(
      getBlockDataPtr0(	bid+dbx)[	vid+dx], 
      getBlockDataPtr0(	bid+dby+dbx)[	vid+dy+dx], 
      getBlockDataPtr0(	bid+dbz+dbx)[	vid+dz+dx], 
      getBlockDataPtr0(	bid+dbz+dby+dbx)[vid+dz+dy+dx] );

#endif
  }

  if constexpr ( doUseIAC )
  {
    if(iac)
    {
      iac->_iposl = ipos;
      iac->_g0246 = g0246;
      iac->_g1357 = g1357;
    }
  }

  }

cache_hit:

    Vec4f r = uvw;

    Vec4f rx = vshuffle_ps<Ax, Ax, Bx, Bx>(r, r); 
    Vec4f gx0123 = vlerp_ps(rx, g0246, g1357);       

    Vec4f ry = vshuffle_ps<Ay, Ay, By, By>(r, r);
    Vec4f gx1032 = vshuffle_ps<Ay, Ax, Bw, Bz>(gx0123, gx0123);

    Vec4f gy0_1_ = vlerp_ps(ry, gx0123, gx1032);
    Vec4f rz = vshuffle_ps<Az, Az, Bz, Bz>(r, r);
    Vec4f gy1_0_ = vshuffle_ps<Az, Az, Bx, Bx>(gy0_1_, gy0_1_);
    Vec4f gz = vlerp_ps(rz, gy0_1_, gy1_0_);

    float out = vfget_x(gz);

  return out;
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template<> inline
float SparseGrid<uint8_t>::interpolateUint8ToFloatFast( 
    				const Vec4f& p0, int opt )
const
{
  Vec4f pos = p0;

  if(!( opt & OPT_IPCORNER )) pos -= Vec4f(.5f,.5f,.5f,0.f);
  
  Vec4f posf = floor(pos);

  Vec4i ipos = truncate_to_int(posf);

  if(!(opt & OPT_IP_NO_BORDER_CHECK))
  {
    if( vany( ipos < vzero() || ipos >= _v4iDomMax )) 
    {
      float p0_[4]; p0.store(p0_);
      return interpolateToFloat<IP_LINEAR>(p0_,opt);
    }
  }

  Vec4i bpos = getBlockCoords( ipos ),
	vpos = getVoxelInBlockCoords( ipos );

  int bid = getBlockIndex( vget_x(bpos), vget_y(bpos), vget_z(bpos));

  UT_ASSERT2(bid>=0&&bid<_nBlocks);

  Vec4f	uvw = pos-posf;

  int vx0=vget_x(vpos),
      vy0=vget_y(vpos),
      vz0=vget_z(vpos);


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
	//

  int vxMax=_bsx-1,
      vid = getVoxelInBlockIndex(vx0,vy0,vz0);

  Vec4f g0246,
	g1357;

  UT_ASSERT2( vget_w(vpos)==0 ); 
  if( !vany( vpos==vxMax ) ) 
  {
    // Fast case: all within single block
    uint8_t* p = getBlockDataPtr0(bid)+vid;
    int s[4];
    _v4iStridesLerp0246.store(s);
    g0246 = Vec4f( *(p+s[0]), *(p+s[1]), *(p+s[2]), *(p+s[3]) );
    _v4iStridesLerp1357.store(s);
    g1357 = Vec4f( *(p+s[0]), *(p+s[1]), *(p+s[2]), *(p+s[3]) );
  }
  else
  {
    int dx,dbx,
	dy,dby,
	dz,dbz;

    Vec4i db = select( vpos < vxMax, 
	               0, Vec4i(1,_nx,_nxy,0) );
    Vec4i d = select( vpos < vxMax, 
	              Vec4i(1,_bsx,_bsx2,0), -vxMax*Vec4i(1,_bsx,_bsx2,0) );

     dx = vget_x(d);
     dy = vget_y(d);
     dz = vget_z(d);
     dbx = vget_x(db);
     dby = vget_y(db);
     dbz = vget_z(db);

    g0246 = Vec4f(
      getBlockDataPtr0(	bid+0)[		vid+0],
      getBlockDataPtr0(	bid+dby+0)[	vid+dy+0], 
      getBlockDataPtr0(	bid+dbz+0)[	vid+dz+0],
      getBlockDataPtr0(	bid+dbz+dby+0)[	vid+dz+dy+0]); 

    g1357 = Vec4f(
      getBlockDataPtr0(	bid+dbx)[	vid+dx], 
      getBlockDataPtr0(	bid+dby+dbx)[	vid+dy+dx], 
      getBlockDataPtr0(	bid+dbz+dbx)[	vid+dz+dx], 
      getBlockDataPtr0(	bid+dbz+dby+dbx)[vid+dz+dy+dx] );

  }

    Vec4f r = uvw;

    Vec4f rx = vshuffle_ps<Ax, Ax, Bx, Bx>(r, r); 
    Vec4f gx0123 = vlerp_ps(rx, g0246, g1357);       

    Vec4f ry = vshuffle_ps<Ay, Ay, By, By>(r, r);
    Vec4f gx1032 = vshuffle_ps<Ay, Ax, Bw, Bz>(gx0123, gx0123);

    Vec4f gy0_1_ = vlerp_ps(ry, gx0123, gx1032);
    Vec4f rz = vshuffle_ps<Az, Az, Bz, Bz>(r, r);
    Vec4f gy1_0_ = vshuffle_ps<Az, Az, Bx, Bx>(gy0_1_, gy0_1_);
    Vec4f gz = vlerp_ps(rz, gy0_1_, gy1_0_);

    float out = vfget_x(gz);

  return out;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template<> inline
float SparseGrid<uint16_t>::interpolateUint16ToFloatFast( 
    				const Vec4f& p0, int opt )
const
{
  Vec4f pos = p0;

  if(!( opt & OPT_IPCORNER )) pos -= Vec4f(.5f,.5f,.5f,0.f);
  
  Vec4f posf = floor(pos);

  Vec4i ipos = truncate_to_int(posf);

  if(!(opt & OPT_IP_NO_BORDER_CHECK))
  {
    if( vany( ipos < vzero() || ipos >= _v4iDomMax )) 
    {
      float p0_[4]; p0.store(p0_);
      return interpolateToFloat<IP_LINEAR>(p0_,opt);
    }
  }

  Vec4i bpos = getBlockCoords( ipos ),
	vpos = getVoxelInBlockCoords( ipos );

  int bid = getBlockIndex( vget_x(bpos), vget_y(bpos), vget_z(bpos));

  UT_ASSERT2(bid>=0&&bid<_nBlocks);

  Vec4f	uvw = pos-posf;

  int vx0=vget_x(vpos),
      vy0=vget_y(vpos),
      vz0=vget_z(vpos);


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
	//

  int vxMax=_bsx-1,
      vid = getVoxelInBlockIndex(vx0,vy0,vz0);

  Vec4f g0246,
	g1357;

  UT_ASSERT2( vget_w(vpos)==0 ); 
  if( !vany( vpos==vxMax ) ) 
  {
    // Fast case: all within single block
    uint16_t* p = getBlockDataPtr0(bid)+vid;
    int s[4];
    _v4iStridesLerp0246.store(s);
    g0246 = Vec4f( *(p+s[0]), *(p+s[1]), *(p+s[2]), *(p+s[3]) );
    _v4iStridesLerp1357.store(s);
    g1357 = Vec4f( *(p+s[0]), *(p+s[1]), *(p+s[2]), *(p+s[3]) );
  }
  else
  {
    int dx,dbx,
	dy,dby,
	dz,dbz;

    Vec4i db = select( vpos < vxMax, 
	               0, Vec4i(1,_nx,_nxy,0) );
    Vec4i d = select( vpos < vxMax, 
	              Vec4i(1,_bsx,_bsx2,0), -vxMax*Vec4i(1,_bsx,_bsx2,0) );

     dx = vget_x(d);
     dy = vget_y(d);
     dz = vget_z(d);
     dbx = vget_x(db);
     dby = vget_y(db);
     dbz = vget_z(db);

    g0246 = Vec4f(
      getBlockDataPtr0(	bid+0)[		vid+0],
      getBlockDataPtr0(	bid+dby+0)[	vid+dy+0], 
      getBlockDataPtr0(	bid+dbz+0)[	vid+dz+0],
      getBlockDataPtr0(	bid+dbz+dby+0)[	vid+dz+dy+0]); 

    g1357 = Vec4f(
      getBlockDataPtr0(	bid+dbx)[	vid+dx], 
      getBlockDataPtr0(	bid+dby+dbx)[	vid+dy+dx], 
      getBlockDataPtr0(	bid+dbz+dbx)[	vid+dz+dx], 
      getBlockDataPtr0(	bid+dbz+dby+dbx)[vid+dz+dy+dx] );

  }

    Vec4f r = uvw;

    Vec4f rx = vshuffle_ps<Ax, Ax, Bx, Bx>(r, r); 
    Vec4f gx0123 = vlerp_ps(rx, g0246, g1357);       

    Vec4f ry = vshuffle_ps<Ay, Ay, By, By>(r, r);
    Vec4f gx1032 = vshuffle_ps<Ay, Ax, Bw, Bz>(gx0123, gx0123);

    Vec4f gy0_1_ = vlerp_ps(ry, gx0123, gx1032);
    Vec4f rz = vshuffle_ps<Az, Az, Bz, Bz>(r, r);
    Vec4f gy1_0_ = vshuffle_ps<Az, Az, Bx, Bx>(gy0_1_, gy0_1_);
    Vec4f gz = vlerp_ps(rz, gy0_1_, gy1_0_);

    float out = vfget_x(gz);

  return out;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template<> inline
float SparseGrid<float>::interpolateDenseScalarFast( const Vec4f& p0 ) const
{    			 
  Vec4f pos = p0,
	posf = floor(pos);				
  Vec4i ipos = truncate_to_int(posf);			
  Vec4f uvw = pos-posf;					

  // 
  // compute weights
  //
  //IACA_START
  Vec8f   U1(vfget_x(uvw)),U0=1.0f-U1,
  	  U = blend8f<0,8, 1,9, 2,10, 3,11>(U0, U1);	

  Vec8f   V1(vfget_y(uvw)),V0=1.0f-V1,			
	  V = blend8f<0,1, 8,9, 2,3, 10,11>(V0, V1);	

  Vec8f   W1(vfget_z(uvw)),W0=1.0f-W1,			
	  W = blend8f<0,1,2,3, 8,9,10,11>(W0, W1);	
  //IACA_END

  // 
  // gather quantities (SOA)
  //
  LongInt idx = getGridIndex(vget_x(ipos),vget_y(ipos),vget_z(ipos));

  float *pf = _denseData+idx;

  Vec4f Q_0123,
	Q_4567;
	
  Q_0123 = _mm_loadh_pi( _mm_loadl_pi( Q_0123, (__m64 const*)(pf+0) ), 
      					       (__m64 const*)(pf+_sxAct) );

  Q_4567 = _mm_loadh_pi( _mm_loadl_pi( Q_4567, (__m64 const*)(pf+0+_sxyAct) ), 
      					       (__m64 const*)(pf+_sxAct+_sxyAct) );

  Vec8f Q(Q_0123,Q_4567);

  float res = horizontal_add(Q*U*V*W);

  return res;
}

/*-------------------------------------------------------------------------*/
/* Dense grid interpolation						   */
/*-------------------------------------------------------------------------*/
template<> inline
float SparseGrid<float>::interpolateDenseScalar( 
    				float *p0, 
				unsigned options ) const
{
  int i,j,k;
  float buf_x[4],buf_y[4],buf_z[4],
	 f_out=0.0f;

  //UT_ASSERT0(options==0);

  float  x = p0[0],
	 y = p0[1], 
	 z = p0[2];

  if(!(options & OPT_IPCORNER))
  {
    x -= 0.5f;
    y -= 0.5f;
    z -= 0.5f;
  }

  int    ix = floorf(x), 
	 iy = floorf(y), 
	 iz = floorf(z);
  {
    float u = x - ix,
	   v = y - iy,
	   w = z - iz;

    for(k=0;k<2;k++)
    {
      for(j=0;j<2;j++)
      {
	for(i=0;i<2;i++)
	{
	  int ix2 = ix+i,
	      iy2 = iy+j,
	      iz2 = iz+k;
	  clipGridCoords(ix2,iy2,iz2);
	  buf_x[i] = getValueDense_(ix2,iy2,iz2); 
	}
	buf_y[j] = Interpolate1D<float,float,IP_LINEAR>(buf_x,u);
      }
      buf_z[k] = Interpolate1D<float,float,IP_LINEAR>(buf_y,v);	       
    }
    f_out = Interpolate1D<float,float,IP_LINEAR>(buf_z,w);	       
  }

rcCatch:
  return f_out;
}

/*-------------------------------------------------------------------------*/
/* interpolate	(generic, slow reference implementation)		   */
/*-------------------------------------------------------------------------*/
template <typename Data_T> 
template < int ipType_T > 
Data_T SparseGrid<Data_T>::interpolateRef( 
    				float *p0, 
				unsigned options ) const
{
  UT_ASSERT2(_hasCompTimeEntrySz);

  int i,j,k;
  Data_T buf_x[4],buf_y[4],buf_z[4],
	 f_out=zeroValue<Data_T>();
  float po=options&OPT_IPCORNER?0.0f:0.5f,
	 x = p0[0]-po, 	// po=0.5 -> interpolate at pixel centers
	 y = p0[1]-po, 
	 z = p0[2]-po;

  int    ix = floorf(x), 
	 iy = floorf(y), 
	 iz = floorf(z);

  {
    float u = x - ix,
	   v = y - iy,
	   w = z - iz;

    if(ipType_T==IP_LINEAR)
    {
      for(k=0;k<2;k++)
      {
	for(j=0;j<2;j++)
	{
	  for(i=0;i<2;i++)
	  {
	    int ix2 = ix+i,
		iy2 = iy+j,
		iz2 = iz+k;
	    clipGridCoords(ix2,iy2,iz2);
	    buf_x[i] = getValue(ix2,iy2,iz2); 
	  }
	  buf_y[j] = Interpolate1D<Data_T,Data_T,ipType_T>(buf_x,u);
	}
	buf_z[k] = Interpolate1D<Data_T,Data_T,ipType_T>(buf_y,v);	       
      }
      f_out = Interpolate1D<Data_T,Data_T,ipType_T>(buf_z,w);	       
    }
    else if(ipType_T==IP_CUBIC || 
	    ipType_T==IP_CUBIC_MONO||
	    ipType_T==IP_CUBIC_MONO_2 ) 
    {
      for(k=0;k<4;k++)
      {
	for(j=0;j<4;j++)
	{
	  for(i=0;i<4;i++)
	  {
	    int ix2 = ix+i-1,
		iy2 = iy+j-1,
		iz2 = iz+k-1;
	    clipGridCoords(ix2,iy2,iz2);
	    buf_x[i] = getValue(ix2,iy2,iz2); 
	  }
	  buf_y[j] = Interpolate1D<Data_T,Data_T,ipType_T>(buf_x,u);
	}
	buf_z[k] = Interpolate1D<Data_T,Data_T,ipType_T>(buf_y,v);	       
      }
      f_out = Interpolate1D<Data_T,Data_T,ipType_T>(buf_z,w);
    }
  }

rcCatch:
  return f_out;
}

/*=========================================================================
 *
 *
 *  I/O
 *
 *
 * =======================================================================*/

template<typename Data_T>
int writeAsVDB( char *fpath, double domainSize,
    		std::vector<SparseGrid<Data_T>*> sbgGrids,
    		std::vector<char*> sbgNames,
		std::vector<unsigned> sbgOptions
		);

/*=========================================================================
 *
 *
 *  Utility
 *
 *
 * =======================================================================*/

inline void setBlockDataZero( Vec3Float *data, int nVoxelsInBlock )
{
  for(int i=0; i<nVoxelsInBlock; i++) data[i] = 0.0f;
}

inline void setBlockDataZero( float *data, int nVoxelsInBlock )
{
  Vec8f Z(0);
  for(int i=0; i<nVoxelsInBlock; i+=8) Z.store_a(data+i);  
}

inline void setBlockDataZero( double *data, int nVoxelsInBlock )
{
  Vec4d Z(0);
  for(int i=0; i<nVoxelsInBlock; i+=4) Z.store_a(data+i);  
}


inline int clipBox( int x1, int y1, int z1,
	     int x2, int y2, int z2,      
	     int &x, int &y, int &z,
	     int &iDir, int &iFace ) 
{
  int clipped=0; 
  iDir=-1,iFace=-1;

  if(unlikely(z<z1)) { z=z1; clipped=1; iDir=2; iFace=0; } 
  if(unlikely(z>z2)) { z=z2; clipped=1; iDir=2; iFace=1; }

  if(unlikely(y<y1)) { y=y1; clipped=1; iDir=1; iFace=0; } 
  if(unlikely(y>y2)) { y=y2; clipped=1; iDir=1; iFace=1; } 

  if(unlikely(x<x1)) { x=x1; clipped=1; iDir=0; iFace=0; } 
  if(unlikely(x>x2)) { x=x2; clipped=1; iDir=0; iFace=1; } 

  return clipped;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
template<typename Data_T>
int gridDimEqual( SparseGrid<Data_T>* sg1,
		         SparseGrid<Data_T>* sg2 )
{
  return sg1->sx()==sg2->sx() && 
	 sg1->sy()==sg2->sy() && 
	 sg1->sz()==sg2->sz();
}

template<typename Data_T>
int dimEqual( SparseGrid<Data_T>* sg1,
		         SparseGrid<Data_T>* sg2 )
{
  return sg1->sx()==sg2->sx() && 
	 sg1->sy()==sg2->sy() && 
	 sg1->sz()==sg2->sz() &&

	 sg1->nbx()==sg2->nbx() && 
	 sg1->nby()==sg2->nby() && 
	 sg1->nbz()==sg2->nbz() &&

	 sg1->bsx() == sg2->bsx();
}

SparseGrid<float>* createFromDenseGrid( 
      const char *name,
      int sx, int sy, int sz, 
      int blockSize,
      float *data, 
      float zeroThreshold,
      int options = 0
      );  

} // namespace SBG

#endif  // __cplusplus

#endif // SBG_H
