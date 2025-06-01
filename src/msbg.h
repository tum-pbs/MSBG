/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#ifndef MSBG_H
#define MSBG_H

/*=========================================================================
 *
 * 	
 * 	M S B G   
 *
 * 	Multiresolution Sparse Block Grids
 *
 *
 * =======================================================================*/

#include <cmath>
//#include <functional>
#include <vectorclass/vectorclass.h>


#include "globdef.h"
#include "fastmath.h"
#include "wturbglob.h"
#include "mtool.h"
#include "util.h"
#include "plot.h"
#include "pnoise.h"
#include "bitmap.h"
#include "panel.h"
#include "msbgcell.h"
#include "thread.h"
#include "sbg.h"

/*=========================================================================
 *
 * 	
 *
 *
 * =======================================================================*/

#define MSBG_MAX_LEVELS	12  // Max. including Multigrid levels

#define MSBG_TIME_LEVELS 3

#define MSBG_MAX_CHANNELS 50  //must be less than 64  

#define MSBG_AR_MAX_ITER 100 // max. number of adaptive pressure solver relaxation sweeps

#define N_BLK_PARALLEL 8 // parallelize loops if more than 
			  // N_BLK_PARALLEL blocks
#define N_CELLS_PARALLEL 8000
#define MIN_BLOCKS_PARALLEL 5000

#define VOL_FRAC_MIN	1e-5  

#define MSBG_MAX_ASM_STEPS 200 // Max. adaptive smoothing steps

#define MSBG_MAX_CSDF_DIST_DOM 0.1f

#define MSBG_AIRDRAG_VCFL 0.5f

#ifndef NDEBUG
#undef FORCE_INLINE_KERNELS
#endif

//
// Approximate grid spacing used for pressure gradient discretization
// at mixed-resolution cell interfaces as of [Losasso2004] Section 4.2
//
//#define DPMR_GRAD_SPACING 	1.	  // for h_grad = dx
//#define DPMR_GRAD_SPACING .82915619758884996227  // Eucl. distance


#define DPMR_GRAD_SPACING  (3./4.)  // for h_grad = 0.75*dx	


#define DPMR_GRAD_COEFF_CF (1./(4*DPMR_GRAD_SPACING))  // (1/3)

#define DPMR_GRAD_COEFF_FC (2*DPMR_GRAD_COEFF_CF)  // (2/3)

#define BSB_SX	4    // Boundary smoothing blocks dimensions
#define BSB_SY	4
#define BSB_SZ	4

#define GP_THETA_TRANS_RES 1.0f

#ifdef __cplusplus

namespace MSBG
{
using namespace SBG;

  enum
  {
    OPT_LEVEL0_ONLY = (1<<0),
    OPT_EXISTING_BLOCKS_ONLY = (1<<1),
    OPT_VISUALIZE = (1<<2),
    OPT_JACOBI = (1<<3),
    OPT_BOUNDZONE_ONLY = (1<<4),
    OPT_CUTS_SOLID = (1<<5),
    OPT_STREAM_STORE = (1<<6),
    OPT_INTERPOLATE = (1<<7),
    OPT_NO_CHAN_FCDIST = (1<<8),    
    OPT_BC_DIRICHLET = (1<<9),
    OPT_BC_NEUMANN = (1<<10),
    OPT_VERI_RESULT = (1<<11),
    OPT_CALC_RESIDUAL = (1<<12),
    OPT_CALC_PER_BLOCK_RESID = (1<<13),
    OPT_SINGLE_LEVEL = (1<<14),
    OPT_TEST = (1<<15),
    OPT_CLAMP_ZERO = (1<<16),
    OPT_MAC_GRID = (1<<17),
    OPT_REVERSE_ORDER = (1<<18),
    OPT_SINGLE_CHANNEL_FLOAT = (1<<19),
    OPT_SINGLE_CHANNEL_VEC = (1<<20),
    OPT_ALL_CELLS = (1<<21),
    OPT_READ_ONLY = (1<<22),
    OPT_WRITE_ALL = (1<<23),
    OPT_NO_NEIGHBORHOOD = (1<<24),
    OPT_NEIGHBORHOOD = (1<<25),
    OPT_BC_COARSE_LEVEL = (1<<26),
    OPT_SIMPLE_AVERAGE = (1<<27),
    OPT_UNIFORM = (1<<28),
    OPT_8_COLOR_SCHEME = (1<<29),
    OPT_BC_LINEAR_EXTRAPOL = (1<<30)
  };

  enum
  {
    PDE_BIHARMONIC = 2,
    PDE_ANISOTROPIC_DIFFUSION = 3,
    PDE_MEAN_CURVATURE = 4
  };

  enum
  {
    LAPL_MULTIPLY = (1<<1),
    LAPL_RELAX = (1<<2),
    LAPL_CALC_RESIDUAL = (1<<4),
    LAPL_WITH_DOTPROD = (1<<5),
    LAPL_INIT_THREAD_LOCALS = (1<<6)
  };

  enum
  {
    CHK_MULTIGRID = (1<<1)
  };

  enum  // for usage with MultiresSparseGrid::solveEikonalFIM
  {
    EIK_EXTRAPOLATE_FLOAT = (1<<0),  
    EIK_EXTRAPOLATE_VEC3 = (1<<1),
    EIK_CALC_SDF = (1<<2),
    EIK_DO_ZERO_UNDEFINED = (1<<4),
    EIK_PROVIDE_INITIAL_SDF = (1<<5),
    EIK_LEVEL0_ONLY = (1<<6),
    EIK_OPT_16_BIT = (1<<7),
    EIK_CLAMP_TO_NBDIST = (1<<8),
    EIK_OPT_8_BIT = (1<<9),     
    EIK_OPT_L1_DIST = (1<<10),
    EIK_OPT_RED_BLACK = (1<<11),
    EIK_OPT_INTERNAL_RED_BLACK = (1<<12)
  };

  enum
  {
    SDF_OPT_FLUID = (1<<1),
    SDF_OPT_SOLID = (1<<2)
  };

  enum 
  {
    PROJECT_SOLID = (1<<0),
    PROJECT_FLUID = (1<<1)
  };

  // 
  // Domain boundary types
  //
  enum
  {			
    DBC_INVALID = -1,// XXX: do not change !
    DBC_OPEN = 0,	
    DBC_SOLID = 1,	
    DBC_INFLOW = 2,
    DBC_OUTFLOW = 3,
    DBC_PERIODIC = 4,

    DBC_OUT_OF_DOMAIN = 1<<30
  };

  enum
  {			
    DBC_OPT_KILL_AT_SOLID = (1<<1),
    DBC_OPT_KEEP_AIR_PART_AT_OPEN = (1<<2)
  };

  //
  // Operation mode For setVelocityBoundaryConditions()
  //
  enum
  {
    BC_FOR_INTERPOLATION = 1,
    BC_FOR_PROJECTION = 2,
    BC_FOR_GRID2PARTICLES = 3,
    BC_FOR_VEL_DIFFUSION = 4
  };
  
  enum
  {
    GRAD_CENTRALDIFF = 1,
    GRAD_SOEBEL_3 = 2,
    GRAD_QUAD_BSPLINE = 3,
    GRAD_CUBIC_BSPLINE = 4,
    GRAD_LINEAR = 5
  };

  //
  // Time stepping schemes
  //
  enum
  {
    TMS_STANDARD = 1,   // XXX: do not change
    TMS_SCHECHTER_BDF2 = 8,  // Schechter/bridson "Evolving Sub-Grid 
    			     // Turbulence for Smoke Animation"
    TMS_REFLECTION = 9,
    TMS_REFLECTION_2 = 10,
    TMS_REFLECTION_PARTICLE_BASED = 11
  };

  //
  // Free surface type for iquid simulation
  //
  enum
  {
    LS_TYPE_PARTICLE_DENSITY = 1,   // XXX: do not change
    LS_TYPE_PARTICLE_SDF = 2 
  };

  // Ghost pressure for second order accurate free BC
  enum
  {
    GP_TYPE_NONE = 0,
    GP_TYPE_STD = 1, // XXX: do not change
    GP_TYPE_FACE_DENSITY = 2
  };

  enum
  {
    REFINE_OBST_ALWAYS = 0, // XXX: do not change
    REFINE_OBST_AT_LIQUID_SURFACE_ONLY = 1
  };

  // 
  enum
  {
    LIQ_TYPE_SOFT_LEVELSET = 1,   // XXX: do not change

    LIQ_TYPE_PARTICLE_POS = 3  // Surface based on particle 
      				     // positions only
  };

  // 
  // Data channel ids
  //
  enum  
  {
    CH_NULL = 0, // XXX do not change CH_NULL==0

    // 
    CH_FLOAT_1 = 1,
    CH_DIVERGENCE = 2,
    CH_PRESSURE = 3,
    CH_CG_P = 4,
    CH_CG_Q = 5,
    CH_FLOAT_2 = 6,
    CH_FLOAT_3 = 7,
    CH_FLOAT_TMP_3 = 8,
    CH_DENSITY_DIFF = 9,
    CH_FLOAT_4 = 10,
    CH_CURVATURE = 11,
    CH_HEAT = 12,
    CH_DIAGONAL = 13, // Diagonal coeff for resolution-transition cells
    CH_FLOAT_8 = 14, // Signed distance to liquid surface
    CH_DIVERGENCE_ADJ = 15,
    CH_FACE_AREA = 16,
    CH_FLOAT_5 = 17,
    CH_MASS_DENSITY = 18,
    CH_UINT16_2 = 19,

    // XXX: ATTENTION when adding new channels: update getChannelAddr
    
    // Vector channels
    CH_VEC3_1 = 20,
    CH_VEC3_2 = 21,
    CH_VEC3_3 = 22,
    CH_VEC3_4 = 23,
    CH_VELOCITY_AVG = 24,
    CH_FACE_DENSITY = 25,

    // XXX: ATTENTION when adding new channels: update getChannelAddr

    // Other channels
    CH_CELL_FLAGS = 26,
    CH_CELL_FLAGS_TMP = 27,
    CH_DIST_FINECOARSE = 28,
    CH_PARTICLE_CELLS = 29,
    CH_UINT16_1 = 30,

    CH_VELOCITY_AIR = 31,
    CH_VELOCITY_AIR_DIFF = 32,
    CH_FACE_COEFF = 33,

    CH_FLOAT_7 = 34,

    #ifdef SOLVE_PRESSURE_DOUBLE
    CH_FLOAT_TMP_PS = 35,
    #endif

    CH_FLOAT_6 = 36,
    CH_SOOT_DIFF = 37,

    CH_UINT8_1 = 38,
    CH_UINT8_2 = 39,
    CH_UINT8_TMP = 40,
    CH_UINT8_TMP_2 = 41,

    CH_PRESSURE_OLD = 42,

    CH_HEAT_DIFF = 43

    // XXX: ATTENTION when adding new channels: update getChannelAddr
    //

  };

  #ifndef SOLVE_PRESSURE_DOUBLE
  enum
  {
    CH_FLOAT_TMP_PS = CH_FLOAT_2,
  };
  #endif

  enum // Generic channel-IDs for OPT_SINGLE_CHANNEL
  {
    CH_GEN_FLOAT = CH_FLOAT_1,
    CH_GEN_FLOAT_2 = CH_DENSITY_DIFF,
    CH_GEN_VEC3_FLOAT = CH_VEC3_1,
  };
  #ifdef RSURF_8_BIT

  typedef uint8_t RenderDensity;
  enum 
  {
    CH_RENDER_DENSITY = CH_UINT8_1,
    CH_RENDER_WHITEWATER = CH_UINT8_2,
    CH_RENDER_TMP = CH_UINT8_TMP,
    CH_RENDER_TMP_2 = CH_UINT8_TMP_2
  };

  #elif defined ( RSURF_16_BIT )

  typedef uint16_t RenderDensity;
  enum 
  {
    CH_RENDER_DENSITY = CH_CELL_FLAGS,
    CH_RENDER_WHITEWATER = CH_CELL_FLAGS_TMP,
    CH_RENDER_TMP = CH_UINT16_1,
    CH_RENDER_TMP_2 = CH_DIST_FINECOARSE
  };

  #else

  typedef float RenderDensity;
  enum 
  {
    CH_RENDER_DENSITY = CH_GEN_FLOAT,
    CH_RENDER_WHITEWATER = CH_GEN_FLOAT_2,
    CH_RENDER_TMP = CH_FLOAT_2,
    CH_RENDER_TMP_2 = CH_FLOAT_3
  };

  #endif

  template< typename DataSDF_T > DataSDF_T EIK_SDF_INFINITY( void );
  template<> inline float EIK_SDF_INFINITY<float>( void ) 
  { 
    return SDF_INFINITY; 
  }
  template<> inline uint16_t EIK_SDF_INFINITY<uint16_t>( void ) 
  { 
    return SDF_INFINITY_UI16; 
  }
  
  //
  // FaceDensity
  //

  #ifdef FACE_DENSITY_16_BIT

  typedef Vec3Uint16 FaceDensity;
  
  #define FACE_DENSITY_TO_FLOAT( fui ) ((fui)*(1.0f/UINT16_MAX))

  #if 0
  inline Vec3Uint16 FACE_DENSITY_FROM_FLOAT( Vec3Float& f ) 
  {
    Vec4f vf = Vec4f().load_partial(3,f.getAddress());
    vf = min(max(vf,0.f),1.f);
    Vec4ui tmp = (Vec4ui)round_to_int( vf * (UINT16_MAX)); 
    return Vec3Uint16(vget_x(tmp),vget_y(tmp),vget_z(tmp));
  }
  #endif

  inline uint16_t FACE_DENSITY_FROM_FLOAT( float f ) 
  {
    UT_ASSERT2( std::isfinite(f) && !(f<0.0f||f>1.0f));
    uint32_t tmp = roundf(std::max(0.f,std::min(1.f,f))*UINT16_MAX);
    return (uint16_t)tmp;
  }

  #else
  
  typedef Vec3Float FaceDensity;
  #define FACE_DENSITY_TO_FLOAT( f ) (f)
  #define FACE_DENSITY_FROM_FLOAT( f ) (f)
  
#endif

  //
  // MassDensity
  //
  
  #ifdef MASS_DENSITY_16_BIT

  typedef uint16_t MassDensity;
  
  #define MASS_DENSITY_FP16_FACTOR (2048.0f)
  #define MASS_DENSITY_TO_FLOAT( fui ) ((fui)*(1.0f/MASS_DENSITY_FP16_FACTOR))
  #define MASS_DENSITY_FROM_FLOAT( f ) \
    ((uint16_t)(std::max(0.f, std::min( (float)(UINT16_MAX), \
	    roundf((f)*MASS_DENSITY_FP16_FACTOR)))))

  #else

  typedef float MassDensity;

  #define MASS_DENSITY_TO_FLOAT( f ) (f)
  #define MASS_DENSITY_FROM_FLOAT( f ) (f)

  #endif

  //
  // RenderDensity
  //
  

  // 
  // interpolateRenderDensity
  //
  float interpolateRenderDensity( SparseGrid<RenderDensity> *sgDens, 
      				  const Vec4f& pos, int opt=0 );

  inline float interpolateRenderDensity( SparseGrid<float> *sgDens, 
      				  const Vec4f& pos, int opt ) 
  {
    return sgDens->interpolateFloatFast( pos, opt );
  }

  inline float interpolateRenderDensity( SparseGrid<uint16_t> *sgDens, 
      				  const Vec4f& pos, int opt ) 
  {
    float f = sgDens->interpolateUint16ToFloatFast( pos, opt );
    RENDER_DENS_TO_FLOAT_UI16(f);
    UT_ASSERT2(std::isfinite(f) && !( f<0.0f || f>1.0f));
    return f;
  }

  inline float interpolateRenderDensity( SparseGrid<uint8_t> *sgDens, 
      				  const Vec4f& pos, int opt ) 
  {
    float f = sgDens->interpolateUint8ToFloatFast( pos, opt );
    RENDER_DENS_TO_FLOAT_UI8(f);
    UT_ASSERT2(std::isfinite(f) && !( f<0.0f || f>1.01f));
    return f;
  }

  // 
  // interpolateRenderDensityWithDerivs
  //
  void interpolateRenderDensityWithDerivs( 
      	SparseGrid<RenderDensity> *sgDens, int ipType_T, const Vec4f& p0, unsigned opt,
	float *outVal, Vec4f *outGrad ); 

  inline void interpolateRenderDensityWithDerivs( 
      	SparseGrid<float> *sgDens, int ipType_T, const Vec4f& p0, unsigned opt,
	float *outVal, Vec4f *outGrad )
  {
    sgDens->interpolateWithDerivs<IP_LINEAR>(ipType_T, p0, opt, outVal, outGrad); 
  }

  inline void interpolateRenderDensityWithDerivs( 
      	SparseGrid<uint16_t> *sgDens, int ipType_T, const Vec4f& p0, unsigned opt,
	float *outVal, Vec4f *outGrad )
  {
    sgDens->interpolateWithDerivs<IP_LINEAR>(ipType_T, p0, opt, outVal, outGrad); 
    RENDER_DENS_TO_FLOAT_UI16(*outVal);
    RENDER_DENS_TO_FLOAT_UI16(*outGrad);    
  }

  // 
  // interpolateRenderDensityWithSecondDerivs
  //
  void interpolateRenderDensityWithSecondDerivs( 
      	SparseGrid<RenderDensity> *sgDens, const Vec4f& p0, unsigned opt,
	float *outVal, float *outGrad, float *outHess ); 

  inline void interpolateRenderDensityWithSecondDerivs( 
      	SparseGrid<float> *sgDens, const Vec4f& p0, unsigned opt,
	float *outVal, float *outGrad, float *outHess )
  {
    sgDens->interpolateWithSecondDerivs(p0, opt, outVal, outGrad, outHess); 
  }

  inline void interpolateRenderDensityWithSecondDerivs( 
      	SparseGrid<uint16_t> *sgDens, const Vec4f& p0, unsigned opt,
	float *outVal, float *outGrad, float *outHess ) 
  {
    sgDens->interpolateWithSecondDerivs(p0, opt, outVal, outGrad, outHess);
    RENDER_DENS_TO_FLOAT_UI16(outVal[0]);
    for(int k=0;k<3;k++)
    {
      RENDER_DENS_TO_FLOAT_UI16(outGrad[k]);    
    }
    for(int k=0;k<6;k++)
    {
      RENDER_DENS_TO_FLOAT_UI16(outHess[k]);    
    }
  }


  //

  const float INVALID_VELOCITY = (1e20f);

  // Block-flags
  
  typedef uint16_t BlockFlags;

  enum  
  {
    BLK_HAS_AIR = (1<<0),
    BLK_DOM_BORDER = (1<<1),
    BLK_COARSE_FINE = (1<<2),
    BLK_FINE_COARSE = (1<<3),
    BLK_EMPTY = (1<<4),
    BLK_NO_FLUID = (1<<5),
    BLK_ONLY_FLUID = (1<<6),
    BLK_BOUNDARY_ZONE = (1<<7),
    BLK_HAS_UNDEFINED = (1<<8), 
    BLK_SLOW = (1<<9),
    BLK_TMP_MARK = (1<<10),
    BLK_FIXED = (1<<11),
    BLK_EXISTS = (1<<12),  // touched by fluid inclusive (advection) buffer zone 
    BLK_UNIFORM_EMPTY = (1<<13), // either uniform SOLID or OPEN and no resolution trans  
    BLK_HAS_DEFINED = (1<<14), // for Eikonal solver has 'inside' cell(s)
    BLK_CUTS_SOLID = (1<<15)
  };

  const int BLK_EX_SLOW_FAST = BLK_BOUNDARY_ZONE;
  const int BLK_EX_NEAR_FLUID = BLK_SLOW;

  static inline bool BLK_IS_RES_BORDER( BlockFlags blk ) 
  { 
    return blk & (BLK_COARSE_FINE | BLK_FINE_COARSE); 
  }
  
  typedef enum
  {
    INVALID = -1,
    FINE_FINE = 0,
    COARSE_FINE = 1,
    FINE_COARSE = 2
  }
  TransResType;

  // 
  // Cell flags
  //
  typedef uint16_t CellFlags;

  enum	
  {			
    CELL_SOLID  =	(1<<0),  // cell flags. XXX do not change !
    CELL_FIXED =    	 (1<<1),  // Dirichlet condition  
    CELL_AIR 	 =	(1<<2),
    CELL_AIR_BORDER =  (1<<3),
    CELL_EMPTY_U =  	(1<<4),
    CELL_EMPTY_V = 	(1<<5),
    CELL_PARTIAL_SOLID = (1<<6),
    CELL_OUT_OF_DOMAIN = (1<<7),
    CELL_BLK_BORDER =	(1<<8),  // Cell is part of block border
    CELL_COARSE_FINE =	(1<<9),  // Cell is part of coarse-fine block interafce
    CELL_FINE_COARSE =	(1<<10),  // Cell is part of fine-coarse block interface
    CELL_EMPTY_C = 	(1<<11),   
    CELL_VOID =         (1<<12),
    CELL_BOUNDARY_ZONE = (1<<13),
    CELL_TMP_MARK_2 = (1<<14), 
    CELL_EMPTY_W = (1<<15)
  };


  enum
  {
    CELL_EMPTY_VAL = CELL_VOID | 
      		     CELL_EMPTY_U | CELL_EMPTY_V | CELL_EMPTY_W | CELL_EMPTY_C
  };
  
  static inline 
    int CELL_IS_LIQUID( CellFlags cell ) 
    { 
      return !(cell & (CELL_SOLID|CELL_VOID|CELL_AIR)); 
    }

  static inline // XXX do not change
    int CELL_IS_FLUID_( CellFlags cell ) { return !(cell & (CELL_SOLID|CELL_VOID)); }

  static inline // XXX do not change
    int CELL_IS_AIR( CellFlags cell ) { return (cell & CELL_AIR); }

  typedef float float4[4];
  typedef Vec3Float Vec3Float4[4];

  static inline int getBlockColor8(int bx, int by, int bz )
  {
    // int icol2 = (( bx & 1 )<<2) + ((by & 1)<<1) + ((bz & 1));  // ??? Fang et. al
    return (( bx & 1 )) + ((by & 1)<<1) + ((bz & 1)<<2);
  }


  // 2x2 sub-block index
  template <int x, int y, int z> int sbix(int bsx, int bsx2)
  {
    return x+y*bsx+z*bsx2;
  }

  typedef struct
  {
    int16_t x,y,z;
  }
  CoordsInt16;

  typedef struct 
  {
    float  w0[3], w0_[3],   // 1D weights
	   w1[3], w1_[3],
	   w2[3], w2_[3];
	  
    float wDerivs[3*3*3][3],     // 3D weights
          w[3*3*3];     // 3D weights
  }
  QuadSplineWeights;

  // Block info
  typedef struct // __attribute__ aligned(CPU_CACHELINE_SZ)
  {
    uint16_t level,
             flags;
  }
  BlockInfo;

  typedef struct
  {
    CellFlags flags;
    int level;
    int bid,vid;
  }
  CellInfo;

  typedef struct
  {
    float *data;
    float *dataFaceArea[3];
    Vec3Float *dataVec;
    CellFlags *dataFlg;
    uint16_t flags,
             level;
  }
  NeighBlkInfo2;

  typedef struct
  {
    int	bsx,
	  bsx2,
	  bsxHi,
	  bsx2Hi,
	  bsxMsk,
	  bsxbsxMsk,
	  bsx2bsxMsk,
	  level;
    float h,hr;
    double h2r,
	   h2r8,
	   h2r4r;
  }
  BlkAuxInfo;


  class MultiresSparseGrid;

  template<typename T>
  SparseGrid<T> *getChannel( MultiresSparseGrid *msbg,
			     int chanId, int level=0, int levelMg=0 );

  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
  class LoopReduction
  {
    public:

    LoopReduction()
    {
      this->init();
    }

    void init( void )
    {
      _valid = 0;
      MT_STAT_INIT(&_stat);
      _longIntSum=0; 
      _longDoubleSum=0; 
    }

    void reduce( LoopReduction& reductionPriv )
    {
      MT_STAT_RESULT(&reductionPriv._stat);
      MT_STAT_UPDATE_STAT(&_stat, &reductionPriv._stat);
      _longIntSum += reductionPriv._longIntSum;
      _longDoubleSum += reductionPriv._longDoubleSum;
    }

    void updateStat( float f )
    {
      MT_STAT_UPDATE( &_stat, f );
    }

    void updateSumLongInt( LongInt val ) { _longIntSum+=val; }

    void finalize()
    {
      MT_STAT_RESULT(&_stat);
      _valid = 1;
    }

    const MtStat& getStat() const { UT_ASSERT2(_valid); return _stat; }
    LongInt getSumLongInt() const { UT_ASSERT2(_valid); return _longIntSum; }

    private:

    int _valid;
    MtStat _stat;
    LongInt _longIntSum;
    long double _longDoubleSum;
  };

  // Visualization
  enum
  {
    VIS_OPT_TEST = (1<<0),    // 
    VIS_OPT_OBST = (1<<1),    // Visualize obstacles
    VIS_OPT_RESLEV = (1<<2),   // Visualize resolution levels 
    VIS_OPT_DISTFIELD = (1<<3),
    VIS_OPT_NONRM = (1<<4),  // do not normalize
    VIS_OPT_BLOCKS = (1<<5),
    VIS_OPT_BOUND_BLOCKS = (1<<6),
    //VIS_OPT_USE_ROI_ = (1<<7),
    VIS_OPT_THRESHCOL = (1<<8),
    VIS_OPT_SUM_SLICES = (1<<9),
    VIS_OPT_DOMAIN_BC = (1<<10),
    VIS_OPT_AMBIENTOCCLUSION = (1<<11),
    VIS_OPT_ZOOM = (1<<12),
    VIS_OPT_GRADIENT = (1<<13),
    VIS_OPT_FLAGTMPMARK = (1<<14),
    VIS_OPT_BLKBORDER = (1<<15),
    VIS_OPT_COLOR_PAL = (1<<16),
    VIS_OPT_TIMELEVEL = (1<<17),
    VIS_OPT_COLOR_RGB = (1<<18),
    VIS_OPT_RENDERDENS = (1<<19),
    VIS_OPT_SIGN_NZ = (1<<20),
    VIS_OPT_DIR_X = (1<<21),
    VIS_OPT_DIR_Y = (1<<22),
    VIS_OPT_DIR_Z = (1<<23),
    VIS_OPT_GHOST_FINE = (1<<24),
    VIS_OPT_ABSMAX_POS = (1<<25),
    VIS_OPT_LEVEL0_ONLY = (1<<26),
    VIS_OPT_DISPLAY_GAMMA = (1<<27),
    VIS_OPT_AIR_VELOCITY = (1<<28),
    VIS_OPT_COLOR_PHASES = (1<<29)
  };

#define VIS_COL_INVALID_QUANTITY BMP_MKRGB(70,0,0)

  typedef struct 
  {
    int transResType;	// Resolution transition type
    CellFlags flags[4];  
  }
  CellFaceNeighbors;

  // Permute from generic ( dir_normal, dir_tan1, dir_tan2 ) 
  // to physical ( x, y, z ) grid coordinates
  static const int
     axesPermutation[3][3] = { {0,1,2}, {1,0,2}, {2,0,1} },
     axesPermutationInv[3][3] = { {0,1,2}, {1,0,2}, {1,2,0} };

/*-----------------------------------------------------------------------*/
/* 									 */
/*-----------------------------------------------------------------------*/
class MultiresSparseGrid
{
  public:
  class BlockIterator;
//  class EikonalSolver;

  static
  MultiresSparseGrid *create( 
                  const char *name,
		  // Grid resolution at finest level 
		  int sx0, int sy0, int sz0, 
		  // Block size at fienst level (power-of-two)
		  int blockSize0 = -1, 
		   // Width of transition zones between resolutin levels
		   // in units of number of finer-resolution cells 
		  int dTransRes = -1, 
		  int maxSizeMB = 0,
		  int initialSizeMB = -1,
		  unsigned options = 0,
		  WTURB_MiscConf *miscConf=NULL,
		  WTURB_PSolverConf *psConf=NULL,
		  int nLevels = -1
		  );  
  static
  void destroy( MultiresSparseGrid*& msbg ); 

  void initStep( void );

  template<typename Data_T>
  void setChannelConstValue( int chan, int constval, Data_T val, int levelMg=0 );

  template<typename Data_T>
  void setChannelMemResource( int chan, void *memResource, int levelMg=0 );

  void downsampleGhostBlockFloat( 
    			  float **dataHaloBlocks,
    			  int bid, int level, int levelMg, // hi-res source level
    			  int iChanDataVec );

  void downsampleGhostBlockSOAtoVec3Float( 
    			  float **dataHaloBlocks,
    			  int bid, int level, int levelMg, // hi-res source level
    			  int iChanDataVec );

  void prepare(
      int iChanData,	    // If data channel(s) changed
      int setSideObstacles = 1,
      int doGhostCellsOnly = 1
      );

  void prepareGhostBlocks( int iChanData );

  int dimEqual( const MultiresSparseGrid *msg2 ) const
  {
    return 
        msg2->sg0()->sx() == _sg0->sx() &&
        msg2->sg0()->sy() == _sg0->sy() &&
	msg2->sg0()->sz() == _sg0->sz()  && 
        msg2->sg0()->nbx() == _sg0->nbx() &&
        msg2->sg0()->nby() == _sg0->nby() &&
	msg2->sg0()->nbz() == _sg0->nbz() &&
	msg2->nBlocks() == nBlocks();
  }

  inline float sharpen_profile( float a, float f ) const
  {    
    return f<0.5f ? (powf(2.f*f,a))/2.f : ((1.f-powf(2.f*(1.f-f),a))/2.f+0.5f);
  };

  void regularizeRefinementMap( int *blockmapIn,
    			         int *blockmapOut=NULL, 
      				 unsigned *blockFlags=NULL,
				 int rbzWidthBlk=0
				);

  void setRefinementMap( int *blockLevels, unsigned *blockFlags=NULL, 
			 int levelMgMax=-1,
			 MSBG::MultiresSparseGrid *msbgObst=NULL,
			 bool doInitCellFlags = true,
			 bool doResetChannels = true
			 );

  int setDomainBoundarySpec_( 
	int dbLeft,  int dbRight, // left,right (X)
	int dbBottom, int dbTop, // bottom,top (Y)
	int dbFront, int dbBack,  // front,back (Z)
	int *dbOpt=NULL
      );

  /*-----------------------------------------------------------------------*/
  /* 								           */
  /*-----------------------------------------------------------------------*/
    
  // Compute Coordinates of cell center at level 0
  void globCellCenterCoords( int ix, int iy, int iz, int level,// IN 
      			     float *posOut
      				) 
  {
    UT_ASSERT2(level>=0&&level<_nLevels);
    posOut[0] = (ix+0.5f)*(1<<level);
    posOut[1] = (iy+0.5f)*(1<<level);
    posOut[2] = (iz+0.5f)*(1<<level);    
  }

  Vec3Float globCellCenterCoords( 
      		int ix, int iy, int iz, int level=0 ) 
  {
    return Vec3Float( (ix+0.5f)*(1<<level),
		      (iy+0.5f)*(1<<level),
		      (iz+0.5f)*(1<<level) );
  }

  Vec3Float globCellCenterCoords( const Vec3Int& gc, int level=0 )
  {
    return globCellCenterCoords( gc.x, gc.y, gc.z, level );
  }
  
  BlockInfo *getBlockInfo( int bid, int levelMg=0 ) const 
  {
    if(levelMg>=_nLevels) return NULL;
    UT_ASSERT2(bid>=0&&bid<_nBlocks);
    return &_blockmap[levelMg][bid];
  }

  BlockInfo *getBlockInfo0( int bid ) const 
  {
    UT_ASSERT2(bid>=0&&bid<_nBlocks);
    return &_blockmap[0][bid];
  }

  BlockInfo *_blockmap[ MSBG_MAXRESLEVELS ];
  
  //

  int8_t   *_blockTimeLevels;
  float *_mgPerBlockResid;
  //
  // Liquid simulation settings
  //
  int _liqActive,
      _liqLsTyp,
      _liqLsExtNbSDF;
  int _liqGpActive;
  double _liqDensThr,
	 _liqMass2DensClip;
  double _liqSurfParticleRadius;

  // 
  // BlockIterator: fast iteration over all cells in a block
  //
  // 
  // Obtain cell neighbors (over 6 faces)
  //
 
  class BlockIterator
  {
    public:
      
    enum{ NBX=3,
	  NBY=3,
	  NBZ=3 };

    uint32_t x,y,z,
	vid,
	bx,by,bz,
	bid0;

    uint32_t bsx,
	bsxMsk,bsx2,
	bsxbsxMsk,bsx2bsxMsk;
    
    int dx0,dx1,
	dBx0,dBx1;

    int dy0,dy1,
	dBy0,dBy1;

    int dz0,dz1,
	dBz0,dBz1;

    MultiresSparseGrid *msg;
    SparseGrid<CellFlags> *sg;
    int chanData;
    BlkAuxInfo bai;
    NeighBlkInfo2 neighBlks[3*3*3];
    int bid,
	levelMg,
	level;

    float *data0;
    Vec3Float *dataVec0;
    CellFlags *dataFlg0;

    int isNull;

    TransResType transResType;
    CellFlags cfnFlagsFF[7];
    CellFaceNeighbors cfnFlags[7];

    float cfnFloatValsFF[7];
    float4 cfnFloatVals[7+1];

    Vec3Float4 cfnVecVals[7+1];
    int cfnBlockInfo[7];  // Relative block ids

    // Constructor
    BlockIterator( MultiresSparseGrid *msg, 
		   int chanFlags, int chanX,   // Flags & Data channel id
		   int bid, int bx_=-1, int by_=-1, int bz_=-1,
		   int levelMg=0,  //  multigrid level
		   unsigned options=0
		   ) 
    {
      const int offs[6][3] = { {1,0,0}, {-1,0,0},
			       {0,1,0}, {0,-1,0},
			       {0,0,1}, {0,0,-1} };

/*      isNull = msg==NULL;
      if(isNull) return;*/
      
      this->bid = bid;
      this->levelMg = levelMg;
      UT_ASSERT2(bid>=0&&bid<msg->_nBlocks);
      UT_ASSERT2(levelMg>=0&&levelMg<msg->getNumLevels());

      BlockInfo *bi = msg->getBlockInfo(bid,levelMg);
      level = bi->level;

      SparseGrid<CellFlags> *sgHi=NULL, *sgLo=NULL;
      SparseGrid<float> *sgX=NULL,*sgXHi=NULL,*sgXLo=NULL;
      SparseGrid<Vec3Float> *sgVel=NULL,*sgVelHi=NULL,*sgVelLo=NULL;

      sg = msg->getFlagsChannel(chanFlags,level,levelMg);
      if(!( options & OPT_LEVEL0_ONLY ))
      {
	sgHi = msg->getFlagsChannel(chanFlags,level-1,levelMg),
	sgLo = msg->getFlagsChannel(chanFlags,level+1,levelMg);
      }

      if( (chanData = chanX) != CH_NULL )
      {
	if(msg->isVecChannel(chanX))
	{
	  sgVel = msg->getVecChannel(chanX, level, levelMg);
	  if(!( options & OPT_LEVEL0_ONLY ))
	  {
	    sgVelHi = msg->getVecChannel(chanX, level-1, levelMg),
	    sgVelLo = msg->getVecChannel(chanX, level+1, levelMg);
	    UT_ASSERT2(sgVel && (level<=levelMg||sgVelHi) && 
			      (level>=msg->getNumLevels()-1||sgVelLo));
	  }
	}
	else if(chanX)
	{
	  sgX = msg->getFloatChannel(chanX, level, levelMg);
	  if(!( options & OPT_LEVEL0_ONLY ))
	  {
	    sgXHi = msg->getFloatChannel(chanX, level-1, levelMg),
	    sgXLo = msg->getFloatChannel(chanX, level+1, levelMg);
	    UT_ASSERT2(sgX && (level<=levelMg||sgXHi) && 
			      (level>=msg->getNumLevels()-1||sgXLo));
	  }
	}
      }

      bsx = sg->bsx();
      bsxMsk=bsx-1;
      bsx2=bsx*bsx;
      bsxbsxMsk = bsx*bsxMsk;
      bsx2bsxMsk = bsx2*bsxMsk;

      bai.bsx = bsx;
      bai.bsx2 = bsx2;
      bai.bsxHi = 2*bsx;
      bai.bsx2Hi = bai.bsxHi*bai.bsxHi;
      bai.bsxMsk = bsx-1;
      bai.bsxbsxMsk = bsx*bai.bsxMsk;
      bai.bsx2bsxMsk = bai.bsx2*bai.bsxMsk;
      bai.level = level;
      double h = (1<<level)*msg->_dx0; 
      bai.h = h;
      bai.hr = 1./h;
      bai.h2r = 1./(h*h);
      bai.h2r8 = 8*bai.h2r;
      bai.h2r4r = bai.h2r/4;

      if(bx_<0)
      {
	sg->getBlockCoordsById(bid, bx_, by_, bz_);	
      }
      bx = bx_;
      by = by_;
      bz = bz_;

      if(!(options&OPT_NO_NEIGHBORHOOD))
      {
	for(int t=0;t<6;t++)
	{
	  int dx=offs[t][0],
	      dy=offs[t][1],
	      dz=offs[t][2],
	      bx2=bx+dx,
	      by2=by+dy,
	      bz2=bz+dz;
	  float *data2=NULL;
	  Vec3Float *dataVec2=NULL;
	  BlockInfo *bi2 = NULL;
	  CellFlags *dataFlg2=NULL;
	  int level2=level;
	  int flags2=0,
	      bid2=-1;
	  if(sg->blockCoordsInRange(bx2,by2,bz2))
	  {
	    bid2 = sg->getBlockIndex(bx2,by2,bz2);
	    bi2 = msg->getBlockInfo(bid2,levelMg);
	    level2 = bi2->level;
	    flags2 = bi2->flags;
	    if(level2==level)
	    {
	      data2 = sgX ? sgX->getBlockDataPtr(bid2, 0,0) : NULL;
	      dataVec2 = sgVel ? sgVel->getBlockDataPtr(bid2, 0,0) : NULL; 
	      dataFlg2 = sg->getBlockDataPtr(bid2, 0,0);
	    }
	    else if(level2>level)
	    {
	      data2 = sgXLo ? sgXLo->getBlockDataPtr(bid2, 0,0) : 
		              msg->_sparseGrids[level].zeroBlock; 
	      dataVec2 = sgVelLo ? sgVelLo->getBlockDataPtr(bid2, 0,0) : NULL; 
	      dataFlg2 = sgLo ? sgLo->getBlockDataPtr(bid2, 0,0) : NULL;
	    }
	    else
	    {
	      data2 = sgXHi ? sgXHi->getBlockDataPtr(bid2, 0,0) : 
			      msg->_sparseGrids[level].zeroBlock; 
	      dataVec2 = sgVelHi ? sgVelHi->getBlockDataPtr(bid2, 0,0) : NULL; 
	      dataFlg2 = sgHi ? sgHi->getBlockDataPtr(bid2, 0,0) : NULL;
	    }
	    if(!dataFlg2) dataFlg2 = msg->getGhostFlagsBlock(level, DBC_OPEN ); 
	    UT_ASSERT2(dataFlg2);
	    UT_ASSERT2(!sgX || data2);
	    UT_ASSERT2(!sgVel || dataVec2);
	  }
	  else 
	  {
	    int iDir = t/2,
		iSide = 1-(t&1),
	        domBc = msg->domainBc(iDir,iSide);

	    data2=msg->_sparseGrids[level].zeroBlock;
	    dataVec2=msg->_sparseGrids[level].zeroBlockVec3;
	    UT_ASSERT2(msg->_sparseGrids[level].zeroBlockSz == 
			bsx*bsx*bsx*sizeof(float));	    
	    dataFlg2 = msg->getGhostFlagsBlock(level,domBc);
	    flags2 = 0;//BLK_IS_GHOST;
	    level2=level;
	  }
	  int nbiIdx = MT_GXYZ(3,3*3,1+dx,1+dy,1+dz);
	  NeighBlkInfo2 *nbi=&neighBlks[nbiIdx];
	  nbi->data = data2;
	  nbi->dataVec = dataVec2;
	  nbi->dataFlg = dataFlg2;
	  nbi->flags = flags2;
	  nbi->level = level2;
	}
      }
      int nbiIdx = MT_GXYZ(3,3*3,1+0,1+0,1+0);
      NeighBlkInfo2 *nbi0=&neighBlks[nbiIdx];
      nbi0->data = sgX ? sgX->getBlockDataPtr(bid) : NULL;
      nbi0->dataVec = sgVel ? sgVel->getBlockDataPtr(bid) : NULL;
      nbi0->dataFlg = sg->hasData() ? sg->getBlockDataPtr(bid) : NULL; 
      nbi0->flags = bi->flags;
      nbi0->level = level;

      //UT_ASSERT2(nbi0->dataFlg);
      UT_ASSERT2(!sgX || nbi0->data);
      UT_ASSERT2(!sgVel || nbi0->dataVec);

      data0 = nbi0->data;
      dataVec0 = nbi0->dataVec;
      dataFlg0 = nbi0->dataFlg;
    
      bid0 = MT_GXYZ(3,3*3,1,1,1);


      //
      // Reset iterator state
      //
      reset();
    }

    // 
    Vec3Int getGridCoords( void ) const
    {
      int ix,iy,iz;
      sg->getGridCoords( bx, by, bz, x, y, z,
			     ix, iy, iz );      
      return Vec3Int(ix,iy,iz);
    }

    //
    CellFlags cellFlags( void ) const
    {
      return neighBlks[bid0].dataFlg[vid];
    }

    void setCellFlags( CellFlags f )
    {
      neighBlks[bid0].dataFlg[vid] = f;
    }

    // 
    inline void reset()
    {
      xReset();
      yReset();
      zReset();
      vid=0;
    }

    //
    inline int valid(void) { return  z<bsx; } 

    //
    inline void next(void)
    {

     vid++;

     if((x<bsxMsk))
     {
       if(x==0)
       {
	 x++; dBx0=0; dx0=-1; return; 
       }
       else if(x<bsx-2)
       {
	 x++; return;
       }
       else 
       {
	 x++; dBx1=1; dx1=-bsxMsk; return;
       }
     }
     else
     {
       xReset();
     }

     if(y<bsxMsk)
     {
       if(y==0)
       {
	 y++; dBy0=0; dy0=-bsx; return; 
       }
       else if(y<bsx-2)
       {
	 y++; return;
       }
       else 
       {
	 y++; dBy1=NBX; dy1=-bsxbsxMsk; return;
       }
     }
     else
     {
       yReset();
     }

     if(z<bsxMsk)
     {
       if(z==0)
       {
	 z++; dBz0=0; dz0=-bsx2; return; 
       }
       else if(z<bsx-2)
       {
	 z++; return;
       }
       else 
       {
	 z++; dBz1=NBX*NBY; dz1=-bsx2bsxMsk; return;
       }
     }
     else
     {
       z++; return;
     }
    }

    // Random access (slower)
    void set( int x0, int y0, int z0 )
    {
      x = x0;
      y = y0;
      z = z0;
      vid = z0*bsx2  + y0*bsx + x0;

      dx0=-1;
      dx1=1;
      dBx0=0;
      dBx1=0;
      if(x==0) 
      {
	dx0=bsxMsk;
	dBx0=-1;
      }
      else if(x==bsx-1)
      {
	dx1=-bsxMsk;
	dBx1=1;
      }

      dy0=-bsx;
      dy1=bsx;
      dBy0=0;
      dBy1=0;
      if(y==0) 
      {
	dy0=bsxbsxMsk;
	dBy0=-NBX;
      }
      else if(y==bsx-1)
      {
	dy1=-bsxbsxMsk;
	dBy1=NBX;
      }

      dz0=-bsx2;
      dz1=bsx2;
      dBz0=0;
      dBz1=0;
      if(z==0) 
      {
	dz0=bsx2bsxMsk;
	dBz0=-NBX*NBY;
      }
      else if(z==bsx-1)
      {
	dz1=-bsx2bsxMsk;
	dBz1=NBX*NBY;
      }
    }
    private:

    inline void xReset( )
    {
      x=0; 
      dBx0=-1; dBx1=0;
      dx0=bsxMsk; dx1=1;
    }

    inline void yReset( )
    {
      y=0; 
      dBy0=-NBX; dBy1=0;
      dy0=bsxbsxMsk; dy1=bsx;
    }

    inline void zReset( )
    {
      z=0; 
      dBz0=-NBX*NBY; dBz1=0;
      dz0=bsx2bsxMsk; dz1=bsx2;
    } 

  };  // BlockIterator

  /*-------------------------------------------------------------------------*/
  /* 					   				     */
  /*-------------------------------------------------------------------------*/
  template< int doRetVals_T,
    	    int dir_T, // normal direction (x=0,y=1,z=2)
	    int face_T >  // face (0,1)
  static inline
#ifdef FORCE_INLINE_KERNELS
  __attribute__((always_inline))
#endif
  void getNeighMixedRes( BlockIterator& bit, 
      		     const NeighBlkInfo2 *neighBlks,  int bid0,
		     int strid, int stridBlk,
		     const BlkAuxInfo *bai,
		     const int& vx, const int& vy, const int& vz, 
		     const int &vid,
		     CellFaceNeighbors *cfn, // OUT
		     float *floatVals,
		     Vec3Float *vec3Vals
		     )
 {
    #define GET_SUB_CONTRIB( k, dx, dy, dz ) \
    { \
	vid2 = vid2bas + sbix<dx,dy,dz>(bsxHi,bsx2Hi);  \
	cfn->transResType = COARSE_FINE; \
	cfn->flags[k] = dataFlg2[vid2];\
	if(doRetVals_T&2) floatVals[k] = data2[vid2];\
	if(doRetVals_T&8) vec3Vals[k] = nbi->dataVec[vid2];\
    }

    int bid2 = bid0 + stridBlk;
    const  NeighBlkInfo2 *nbi = &neighBlks[bid2];

    float *data2;
    int vid2,vid2bas;

    int	bsx=bai->bsx,
	  bsx2=bai->bsx2,
	  bsxHi=bai->bsxHi,
	  bsx2Hi=bai->bsx2Hi,
	  level=bai->level;

    //if(!(nbi->flags & BLK_IS_GHOST))  // within domain ?
    {
      data2 = nbi->data;
      CellFlags *dataFlg2=nbi->dataFlg;
      if(nbi->level<level)  // coarse - fine
      {
	switch(dir_T)
	{
	  case 0:  // left/right 
	    vid2bas = face_T==0 ? (8*bsx2)*vz + (4*bsx)*vy + (2*bsx-1) :
				  (8*bsx2)*vz + (4*bsx)*vy;
	    GET_SUB_CONTRIB(0, 0,0,0 );
	    GET_SUB_CONTRIB(1, 0,1,0 );
	    GET_SUB_CONTRIB(2, 0,0,1 );
	    GET_SUB_CONTRIB(3, 0,1,1 );
	    break;

	  case 1: // down/up
	    vid2bas = face_T==0 ? (8*bsx2)*vz + 2*vx + (2*bsx*(2*bsx-1)) :
				  (8*bsx2)*vz+2*vx;
	    GET_SUB_CONTRIB(0, 0,0,0 );
	    GET_SUB_CONTRIB(1, 1,0,0 );
	    GET_SUB_CONTRIB(2, 0,0,1 );
	    GET_SUB_CONTRIB(3, 1,0,1 );
	    break;

	  case 2: // front/back
	    vid2bas = face_T==0 ? (4*bsx)*vy + 2*vx + (4*bsx2*(2*bsx-1)) :
				  (4*bsx)*vy + 2*vx;
	    GET_SUB_CONTRIB(0, 0,0,0 );
	    GET_SUB_CONTRIB(1, 1,0,0 );
	    GET_SUB_CONTRIB(2, 0,1,0 );
	    GET_SUB_CONTRIB(3, 1,1,0 );
	    break;
	}
      }
      else if(nbi->level>level ) // fine - coarse
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

	cfn->transResType = FINE_COARSE;
	cfn->flags[0] = dataFlg2[vid2];
	if(doRetVals_T&2) floatVals[0] = data2[vid2];
	if(doRetVals_T&8) vec3Vals[0] = nbi->dataVec[vid2];
      }
      else  // fine-fine
      {
	vid2 = vid+strid;
	cfn->transResType = FINE_FINE;
	cfn->flags[0] = dataFlg2[vid2];
	if(doRetVals_T&2) floatVals[0] = data2[vid2];
	if(doRetVals_T&8) vec3Vals[0] = nbi->dataVec[vid2];
      }
    }  
  }
  enum
  {
    // XXX do not change !
    CN_FLAGS = 1,  
    CN_FLOAT_CHAN = 2, 
    CN_VEC_CHAN = 8,
    CN_HALF_NEIGH = 32
    // XXX do not change !
  };

  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
  template<     	
    	int doRetVals_T  // &CN_FLAGS =cell-flags, 
			 // &2 = float-chan1, 
			 // &4(CN_FLOAT_CHAN)=float-chan2, 
			 // &8(CN_VEC_CHAN=Vec3-chan 
			 // &16=CellAccessors
	>
  static inline 
//#ifdef FORCE_INLINE_KERNELS
  __attribute__((always_inline))
//#endif
  void 
    getCellNeighborhood( 
	  BlockIterator& bit, 
	  CellFaceNeighbors *cfn=NULL, // out: 0=center, 
	  			  //      1=left,2=right, 
				  //      3=down,4=up, 5=front,6=back
	  float (*floatVals)[4] = NULL,		// 6x4 result values  
	  Vec3Float (*vec3Vals)[4] = NULL
	)
  {
    if(!cfn) cfn=bit.cfnFlags;
    if(!floatVals) floatVals=bit.cfnFloatVals;
    if(!vec3Vals) vec3Vals=bit.cfnVecVals;
      //IACA_START
    const NeighBlkInfo2* neighBlks=bit.neighBlks;
    int bid0 = bit.bid0,
	vid = bit.vid;
    // Center value
    const NeighBlkInfo2 *nbi0 = &neighBlks[bid0];
    CellFlags flags = nbi0->dataFlg[vid];
    cfn[0].transResType = FINE_FINE;
    cfn[0].flags[0] = flags;
    if(doRetVals_T & 2 ) floatVals[0][0] = nbi0->data[vid];
    if(doRetVals_T & 8) vec3Vals[0][0] = nbi0->dataVec[vid];

    //IACA_END
    // 6 face neighbors
    if(!(flags & CELL_BLK_BORDER))  
    {
      //IACA_START
      // Fastest case (neighborhood does not cross any block boundary)      
      uint32_t  i1=vid-1, i2=vid+1,
		i3=vid-bit.bsx, i4=vid+bit.bsx,
		i5=vid-bit.bsx2, i6=vid+bit.bsx2;

      CellFlags *dataFlg = nbi0->dataFlg;

      cfn[1].transResType = FINE_FINE;
      cfn[1].flags[0] = dataFlg[i1];
      cfn[2].transResType = FINE_FINE;
      cfn[2].flags[0] = dataFlg[i2];

      cfn[3].transResType = FINE_FINE;
      cfn[3].flags[0] = dataFlg[i3];
      cfn[4].transResType = FINE_FINE;
      cfn[4].flags[0] = dataFlg[i4];

      cfn[5].transResType = FINE_FINE;
      cfn[5].flags[0] = dataFlg[i5];
      cfn[6].transResType = FINE_FINE;
      cfn[6].flags[0] = dataFlg[i6];

      if(doRetVals_T & 2 )
      {
	float *data = nbi0->data;
        floatVals[1][0] = data[i1];
        floatVals[2][0] = data[i2];
        floatVals[3][0] = data[i3];
        floatVals[4][0] = data[i4];
        floatVals[5][0] = data[i5];
        floatVals[6][0] = data[i6];
      }
      if(doRetVals_T & 8)
      {
	Vec3Float *data = nbi0->dataVec;
        vec3Vals[1][0] = data[i1];
        vec3Vals[2][0] = data[i2];
        vec3Vals[3][0] = data[i3];
        vec3Vals[4][0] = data[i4];
        vec3Vals[5][0] = data[i5];
        vec3Vals[6][0] = data[i6];
      }
    }
    else if(!(flags & (CELL_COARSE_FINE|CELL_FINE_COARSE)))
    {
      //IACA_START
      // Medium fast case: FINE-FINE block boundary (no resolution boundary)
      int       i1 = vid + bit.dx0,  i2 = vid + bit.dx1,
		i3 = vid + bit.dy0,  i4 = vid + bit.dy1,
		i5 = vid + bit.dz0,  i6 = vid + bit.dz1;
      int       ib1 = bid0 + bit.dBx0,  ib2 = bid0 + bit.dBx1,
		ib3 = bid0 + bit.dBy0,  ib4 = bid0 + bit.dBy1,
		ib5 = bid0 + bit.dBz0,  ib6 = bid0 + bit.dBz1;		    

      cfn[1].transResType = FINE_FINE;
      cfn[1].flags[0] = neighBlks[ib1].dataFlg[i1];
      cfn[2].transResType = FINE_FINE;
      cfn[2].flags[0] = neighBlks[ib2].dataFlg[i2];

      cfn[3].transResType = FINE_FINE;
      cfn[3].flags[0] = neighBlks[ib3].dataFlg[i3];
      cfn[4].transResType = FINE_FINE;
      cfn[4].flags[0] = neighBlks[ib4].dataFlg[i4];

      cfn[5].transResType = FINE_FINE;
      cfn[5].flags[0] = neighBlks[ib5].dataFlg[i5];
      cfn[6].transResType = FINE_FINE;
      cfn[6].flags[0] = neighBlks[ib6].dataFlg[i6];

      if(doRetVals_T & 2)
      {
        floatVals[1][0] = neighBlks[ib1].data[i1];
	floatVals[2][0] = neighBlks[ib2].data[i2];
	floatVals[3][0] = neighBlks[ib3].data[i3]; 
	floatVals[4][0] = neighBlks[ib4].data[i4]; 
	floatVals[5][0] = neighBlks[ib5].data[i5]; 
	floatVals[6][0] = neighBlks[ib6].data[i6];
      }
      if(doRetVals_T & 8)
      {
        vec3Vals[1][0] = neighBlks[ib1].dataVec[i1];
	vec3Vals[2][0] = neighBlks[ib2].dataVec[i2];
	vec3Vals[3][0] = neighBlks[ib3].dataVec[i3]; 
	vec3Vals[4][0] = neighBlks[ib4].dataVec[i4]; 
	vec3Vals[5][0] = neighBlks[ib5].dataVec[i5]; 
	vec3Vals[6][0] = neighBlks[ib6].dataVec[i6];
      }
      //IACA_END
    }
    else
    {
      //IACA_START
      // Mixed resolution case
      
      // Left-Right
      getNeighMixedRes<doRetVals_T,0,0>( bit, neighBlks, bid0, bit.dx0, bit.dBx0, 
	  		   &bit.bai, bit.x, bit.y, bit.z, vid, 
			   &cfn[1], 
			   (doRetVals_T&2) ? floatVals[1] : NULL, 
			   (doRetVals_T&8) ? vec3Vals[1] : NULL 
			   );

      getNeighMixedRes<doRetVals_T,0,1>( bit, neighBlks, bid0, bit.dx1, bit.dBx1, 
	  		   &bit.bai, bit.x, bit.y, bit.z, vid, 
			   &cfn[2], 
			   (doRetVals_T&2) ? floatVals[2] : NULL, 
			   (doRetVals_T&8) ? vec3Vals[2] : NULL 
			   );
      // Down-Up
      getNeighMixedRes<doRetVals_T,1,0>( bit, neighBlks, bid0, bit.dy0, bit.dBy0, 
	  		   &bit.bai, bit.x, bit.y, bit.z, vid, 
			   &cfn[3], 
			   (doRetVals_T&2) ? floatVals[3] : NULL, 
			   (doRetVals_T&8) ? vec3Vals[3] : NULL 
			   );

      getNeighMixedRes<doRetVals_T,1,1>( bit, neighBlks, bid0, bit.dy1, bit.dBy1, 
	  		   &bit.bai, bit.x, bit.y, bit.z, vid, 
			   &cfn[4], 
			   (doRetVals_T&2) ? floatVals[4] : NULL, 
			   (doRetVals_T&8) ? vec3Vals[4] : NULL 
			   );
      // Front-Back
      getNeighMixedRes<doRetVals_T,2,0>( bit, neighBlks, bid0, bit.dz0, bit.dBz0, 
	  		   &bit.bai, bit.x, bit.y, bit.z, vid, 
			   &cfn[5], 
			   (doRetVals_T&2) ? floatVals[5] : NULL, 
			   (doRetVals_T&8) ? vec3Vals[5] : NULL 
			   );

      getNeighMixedRes<doRetVals_T,2,1>( bit, neighBlks, bid0, bit.dz1, bit.dBz1, 
	  		   &bit.bai, bit.x, bit.y, bit.z, vid, 
			   &cfn[6], 
			   (doRetVals_T&2) ? floatVals[6] : NULL, 
			   (doRetVals_T&8) ? vec3Vals[6] : NULL 
			   );
      //IACA_END
    }
 }

  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
  template<typename func_T>
  void processFloatChannel( int chanSrc, int chanDst, int levelMg, 
      			    unsigned options, func_T func,
			    LoopReduction *reduction=NULL )
  {
    TRC(("%s %d %d %d %d\n",UT_FUNCNAME,chanSrc,chanDst,levelMg,options));
    int doWriteAll = options&OPT_WRITE_ALL;
    int withNeighborhood = options&OPT_NEIGHBORHOOD;
    if(reduction) reduction->init();

    typedef struct
    {
      LoopReduction reductionPriv,
		    *reductionPrivAct;
    } 
    ThreadLocals;

    ThrRunParallel<ThreadLocals> ( 

      nBlocks(),

      [&]( ThreadLocals &tls, int tid )  // Initialize locals
      {
        tls.reductionPrivAct = reduction ? &tls.reductionPriv : NULL;
      },
    
      // Run parallel
      [&](ThreadLocals &tls, int tid, int bid)
      {
	int isSparseGrid = levelMg < _nLevels;
	int level = isSparseGrid ? 
	  		_blockmap[levelMg][bid].level : levelMg;	  	  
	SparseGrid<float> *sgSrc=getFloatChannel(chanSrc,level,levelMg),
	  		  *sgDst=getFloatChannel(chanDst,level,levelMg);
	if(!(sgSrc&&sgDst)) return;
	if(sgSrc->isDenseGrid())
	{
	  float *dataDst = sgDst->getDataPtrGen_(bid,1,
					  doWriteAll||(!isSparseGrid)?0:1),
		*dataSrc = sgSrc->getDataPtrGen_(bid);

	  MultiresSparseGrid::BlockIterator bitDummy(this,CH_CELL_FLAGS,CH_NULL,0,
	      					     -1,-1,-1,0,OPT_NO_NEIGHBORHOOD); 
	  SBG_FOR_EACH_VOXEL_IN_BLOCK_GEN( sgSrc, bid, x, y, z, idx )
	  {
	    func( bitDummy, dataSrc[idx], dataDst[idx], 
		  tls.reductionPrivAct );
	  }
	}
	else
	{
	  BlockInfo *bi=getBlockInfo(bid,levelMg);
	  if((options & OPT_LEVEL0_ONLY)) 
	  {
	    if( bi->level!=0  ||
	       !sgSrc->isValueBlock(bid) ||
	       !sgDst->isValueBlock(bid) ) return;
	  }

	  float *dataDst = sgDst->getDataPtrGen_(bid,1,
					  doWriteAll||(!isSparseGrid)?0:1),
		*dataSrc = sgSrc->getDataPtrGen_(bid);

	  UT_ASSERT2( sgDst->isValueBlock(bid));
	  for( MultiresSparseGrid::BlockIterator 
		    bit( this, CH_CELL_FLAGS, withNeighborhood ? chanSrc : CH_NULL, 
			 bid, -1,-1,-1, levelMg,
			 withNeighborhood ? 0 : OPT_NO_NEIGHBORHOOD );
	       bit.valid(); bit.next() )
	  {
	    func( bit, dataSrc[bit.vid], dataDst[bit.vid], 
		  tls.reductionPrivAct );
	  }
	}
      },

      [&]( ThreadLocals &tls, int tid )  // Reduce locals
      {
	if(reduction)
	{
	  reduction->reduce( tls.reductionPriv );
	}
      }
    );

    if(reduction) reduction->finalize();
    if(!(options & OPT_READ_ONLY)) invalidateForInterpolation(chanDst);
  }

  // 
  // BlockGridInfo for (multithreaded) iteration over blocks
  //
  class BlockGridInfo
  {
    public:
    int _nbx,_nby,_nbz,   // Number of blocks in each direction
	_bsx0,_bsy0,_bsz0;   // Number of cells in un-clipped block         
    int _nbxy,_nBlocks;
    MultiresSparseGrid *_msg;

    // Create block grid info for given (multigrid) level
    BlockGridInfo( MultiresSparseGrid *msg, int levelMg )
    {
      _msg = msg;
      if(levelMg<msg->getNumLevels())
      {
	_nbx = msg->_sg0->nbx();
	_nby = msg->_sg0->nby();
	_nbz = msg->_sg0->nbz();
	_bsx0= _bsy0= _bsz0= msg->_sg0->bsx();
      }
      else
      {
	_bsx0= _bsy0= _bsz0=16;
	SparseGrid<CellFlags> *sg = msg->getFlagsChannel(CH_CELL_FLAGS,levelMg);
	_nbx = MAX((sg->sx()+_bsx0-1)/_bsx0,1); 
	_nby = MAX((sg->sy()+_bsy0-1)/_bsy0,1); 
	_nbz = MAX((sg->sz()+_bsz0-1)/_bsz0,1); 
      }
      _nbxy = _nbx*_nby;
      _nBlocks = _nbx*_nby*_nbz;
      UT_ASSERT0((LongInt)_nbx*(LongInt)_nby*(LongInt)_nbz<2000000000);
      UT_ASSERT0((LongInt)_bsx0*(LongInt)_bsy0*(LongInt)_bsz0<2000000000);
      TRC3(("%s: level=%d %dx%dx%d @ %dx%dx%d\n",UT_FUNCNAME,
	    levelMg,_nbx,_nby,_nbz,_bsx0,_bsy0,_bsz0));
    }

    // Get block access info for a given block id
    void getBlockAccessInfo( 
	  int level, int bid, // IN: block grid info and block id
	  SparseGrid<CellFlags> *sg,
	  int &bx, int &by, int &bz,
	  int &bsx, int& bsy, int& bsz
	     // OUT: block coords: bx,by,bz
	     //      block resolution: bsx, bsy, bsz
	  ) const
    {
      int sx = sg->sx(),
	  sy = sg->sy(),
	  sz = sg->sz();

      // Block coords
      MT_GCOORDS_3D( _nbx,_nbxy, bid,
		     bx,by,bz);

      UT_ASSERT2(G3D_IN_BOX(0,0,0,_nbx-1,_nby-1,_nbz-1,
			    bx,by,bz));
      
      #ifdef UT_ASSERT_LEVEL_2
      if(level<_msg->getNumLevels())
      {
	int bx_,by_,bz_;
	_msg->_sg0->getBlockCoordsById(bid, bx_,by_,bz_);
	UT_ASSERT0(bx_==bx && by_==by && bz_==bz);
      }
      #endif

      // Block resolution
      if(level<_msg->getNumLevels())
      {
	bsx=bsy=bsz=sg->bsx();
      }
      else
      {
	bsx = _bsx0;
	bsy = _bsy0;
	bsz = _bsz0;
	bsx = MIN((bx+1)*bsx-1, sx-1) - bx*bsx + 1;
	bsy = MIN((by+1)*bsy-1, sy-1) - by*bsy + 1;
	bsz = MIN((bz+1)*bsz-1, sz-1) - bz*bsz + 1;
      }
      UT_ASSERT2(bsx>0 && bsy>0 && bsz>0 && 
		 bsx<=sx && bsy<=sy && bsz<=sz);

      TRC3(("%s: level=%d bid=%d %d,%d,%d @ %dx%dx%d\n",UT_FUNCNAME,
	    level,bid, bx,by,bz, bsx,bsy,bsz));
    }
  };

  //
  // Get cell level access to an individual 
  // block (to be used as thread private iterator)
  //
  class BlockAccessor
  {
    public:
    BlockAccessor( BlockGridInfo* bgi, SparseGrid<CellFlags> *sg, 
	           int bid, int level)
    {
      _sg = sg;
      _bsx0 = bgi->_bsx0;
      _bsy0 = bgi->_bsy0;
      _bsz0 = bgi->_bsz0;
      bgi->getBlockAccessInfo(level,bid,sg,_bx,_by,_bz,_bsx,_bsy,_bsz);	
    }

    int getCellIndex( int vx, int vy, int vz ) const
    {
      return _sg->isDenseGrid() ? 
	     _sg->getGridIndex( _bx*_bsx0+vx, 
		                _by*_bsy0+vy, 
				_bz*_bsz0+vz) : 
	     _sg->getVoxelInBlockIndex(vx,vy,vz);
    }

    void getGridCoords( int vx, int vy, int vz,
			int& x, int& y, int& z ) const
    {
      if(_sg->isDenseGrid())
      {
	x = _bx*_bsx0+vx; 
	y = _by*_bsy0+vy;
	z = _bz*_bsz0+vz;
      }
      else _sg->getGridCoords( _bx,_by,_bz, vx, vy, vz, x, y, z);
    }

    SparseGrid<CellFlags> *_sg;
    int _bx,_by,_bz,
	_bsx0,_bsy0,_bsz0,
	_bsx,_bsy,_bsz;	
  };

  class CellAccessor
  {
    public:
    
    //CellAccessor( void ) { };

    int level,
    	bid,	// Block index
	vid,	// Voxel-in-block index
	x,y,z,
	vx,vy,vz;
  };


  CellAccessor  
    getCellAccessor0( Vec4i ipos // IN: global grid pos
		)
  {
    constexpr int bsx0Log2 = MSBG_BSXLOG2_0;
    UT_ASSERT2(bsx0Log2==_bsx0Log2);
    Vec4i bpos = ipos >> bsx0Log2; 
    int bid = _sg0->getBlockIndex(vget_x(bpos),vget_y(bpos),vget_z(bpos));
    int level = getBlockInfo0(bid)->level;
    ipos >>= level;

    int bsxLog2 = bsx0Log2-level, 
	bsx = (1<<bsxLog2);

    CellAccessor ca;

    ca.level = level;
    ca.bid = bid;
    ca.x = vget_x(ipos);
    ca.y = vget_y(ipos);
    ca.z = vget_z(ipos);
    ca.vx = ca.x & (bsx-1);
    ca.vy = ca.y & (bsx-1);
    ca.vz = ca.z & (bsx-1);
    ca.vid = ca.vx | (ca.vy << bsxLog2) | (ca.vz << (2*bsxLog2));

    UT_ASSERT2(ca.vid>=0&&ca.vid<bsx*bsx*bsx);

    return ca;
  }

  CellAccessor  
    getCellAccessor( Vec4f pos // IN: global grid pos
		) 
  {
    Vec4i ipos = truncate_to_int(floor(pos));
    ipos = max(ipos,_mm_setzero_si128());
    ipos = min(ipos,_v4iDomMax);

    return getCellAccessor0(ipos);
  }

  void getCellAccessor( float x0, float y0, float z0, // IN: global grid pos
      		        CellAccessor *ca // OUT
			)
  {
    SparseGrid<Vec3Float> *sg0 = _sparseGrids[0].vec3_1[0];
    int x = floorf(x0),
	y = floorf(y0),
	z = floorf(z0);
    sg0->clipGridCoords(x,y,z);
    int bx,by,bz;
    sg0->getBlockCoords( x, y, z, 
			 bx, by, bz);
    int bid = sg0->getBlockIndex(bx,by,bz);
    UT_ASSERT2(bid>=0&&bid<sg0->nBlocks());
    int level = getBlockLevel(bid);
    UT_ASSERT2(level>=0&&level<_nLevels);
    SparseGrid<Vec3Float> *sg = _sparseGrids[level].vec3_1[0];
    ca->x = x >> level;
    ca->y = y >> level;
    ca->z = z >> level;
    ca->level = level;
    ca->bid = bid;
    sg->getVoxelInBlockCoords(ca->x,ca->y,ca->z, ca->vx,ca->vy,ca->vz);
    ca->vid = sg->getVoxelInBlockIndex(ca->vx,ca->vy,ca->vz);
  }

  //
  // Return globally (i.e. over a given multigrid level) unique 
  // cell identifier as the index of that cell at the finest resolution
  // (for coarse cells this is simply the index of
  // the lower left corner level 0 cell)
  //
  LongInt virtualCellIndex( int x, int y, int z, int level, int levelMg=0 )
  {
    if(levelMg<_nLevels)
    {
      UT_ASSERT2(level<_nLevels);
      SparseGrid<CellFlags> *sg0 = _sparseGrids[0].cellFlags[0];
      x <<= level;
      y <<= level;
      z <<=level;
      UT_ASSERT2(sg0->inRange(x,y,z));
      return sg0->getVirtGridIndex( x,y,z );
    }
    else 
    {
      UT_ASSERT2(level>=_nLevels);
      return _sparseGrids[level].cellFlags[0]->getVirtGridIndex(x,y,z);
    }
  }

  LongInt nVirtualCells( int levelMg )
  {
    UT_ASSERT0(levelMg>=0&&levelMg<_mgLevels);
    return levelMg < _nLevels ? _sg0->nTotVirtVoxels() :
      	_sparseGrids[levelMg].cellFlags[0]->nTotVirtVoxels();
  }

  void interpolateWithDerivs( int ipolType,
			      int chan, 
			      Vec4f pos0, 
			      int levelMg,
			      unsigned opt,
			      float *outVal,
			      Vec4f *outGrad );

  template<int ipolType_T> 
  void interpolate( 
	   int chanId,
	   float *pos0, 		// Position on finest resolution
	   unsigned options,	// SBG::OPT_IPCORNER, SBG::OPT_IPGRIDALGN
	   Vec3Float& dataOut,      // OUT
	   bool doOnlyLevel0=false
	  );

  void interpolate( 
      	   int ipTyp,
	   int chanId,
	   float *pos0, 		// Position on finest resolution
	   unsigned options,	// SBG::OPT_IPCORNER, SBG::OPT_IPGRIDALGN
	   Vec3Float& dataOut      // OUT
	  )
  {
    if(ipTyp==IP_LINEAR_FADE)
    {
      ipTyp = IP_LINEAR;
      options |= OPT_IP_FADE_WEIGHTS;
    }
    if(ipTyp==IP_LINEAR)
      interpolate<IP_LINEAR>(chanId,pos0,options,dataOut);
    else if(ipTyp==IP_CUBIC_MONO_2)
      interpolate<IP_CUBIC_MONO_2>(chanId,pos0,options,dataOut);
    else if(ipTyp==IP_CUBIC)
      interpolate<IP_CUBIC>(chanId,pos0,options,dataOut);
    else if(ipTyp==IP_WENO4)
      interpolate<IP_WENO4>(chanId,pos0,options,dataOut);
    else
    {
      UT_ASSERT0(FALSE);
    }
  }


  float interpolateScalarFast2( int chanId, Vec4f pos, unsigned opt,
      				Vec4f *pOutGrad=NULL );

  template<int ipolType_T> 
  double interpolateScalar( 
	   int chanId,
	   float *pos0, 		// Position on finest resolution
	   unsigned options	// SBG::OPT_IPCORNER, SBG::OPT_IPGRIDALGN
	  )
  {
    Vec3Float tmp;
    interpolate<ipolType_T>(chanId,pos0,options,tmp);
    return tmp[0];
  }

  double interpolateScalar( int ipTyp,
	   int chanId,
	   float *pos0, 		// Position on finest resolution
	   unsigned options	// SBG::OPT_IPCORNER, SBG::OPT_IPGRIDALGN
	  )
  {
    Vec3Float tmp;
    if(ipTyp==IP_LINEAR)
      interpolate<IP_LINEAR>(chanId,pos0,options,tmp);
    else if(ipTyp==IP_CUBIC_MONO_2)
      interpolate<IP_CUBIC_MONO_2>(chanId,pos0,options,tmp);
    else if(ipTyp==IP_CUBIC)
      interpolate<IP_CUBIC>(chanId,pos0,options,tmp);
    else if(ipTyp==IP_WENO4)
      interpolate<IP_WENO4>(chanId,pos0,options,tmp);
    return tmp[0];
  }

  bool isValidForInterpolation( int iChan )
  {
    return _validForInterpolation & (((int64_t)1)<<iChan);
  }

  void invalidateForInterpolation( int iChan )
  {
    FLAG_RESET(_validForInterpolation,(1<<(uint64_t)iChan));      
  }

  void setChannelModified( int iChan )
  {
    invalidateForInterpolation( iChan );
  }
    

  void computeDistanceFieldSimpleUint8( int ch /* channel */, 
					int chTmp,
					std::vector<int> *blockList,
					uint8_t uiPhiThr,
					int nIter );

  void computeDistanceFieldUint8( int ch /* target channel  in/out */, 
    						    int chTmp, int chTmp2, 
    					            std::vector<int> *blockList,
						    float phiThr,  // phi < piThr => fixed value
    					            uint8_t distMax_ );


  template< typename Data_T >
  void applyChannelPdeFast( 
      				int chPhi, int chTmp, int chTmp2,
				std::vector<int> *blockList,
				int laplTyp, 
      				int doConstrZeroOne,
      				float densThr, 
      				float densThrTol1,
      				float densThrTol2,
      				int numIter,      				
      				float gridSpacing=1,
				float timeStep=0.025, 
				float velMagThr=0,
				MultiresSparseGrid *msbgVel=NULL,
				float lambda=0,
				float lambdaSDF=0,
				float taubin=0,
				float redist=0,
				float anisDiffAlpha=0,
				int veldep_on=0,
				float *veldep_fit=NULL,
				int visOpt=0
				);

  void applyChannelPdeCnf( int chan, WTURB_LapSmConf *cnf )
  {
    applyChannelPdeFast<float>( chan, CH_FLOAT_2, CH_FLOAT_3,
		      NULL,
		      cnf->typ, 
		      TRUE,
		      cnf->thr, 
		      cnf->thr_tol1,
		      cnf->thr_tol2,
		      cnf->iter, 
		      1.0f,
		      cnf->dt,
		      cnf->velmag_thr,
		      this,
		      cnf->lambda,
		      cnf->lambda_sdf,
		      cnf->taubin,
		      cnf->redist,
		      0,
		      0
		      );  
  }


  void downsampleFloatChannelNew( int levelMg,  // target level
			       int chSrc,
			       int chDst,
			       unsigned options );

  template<typename data_T, typename dataIntern_T >
  void downsampleFaceDensity( int levelMg,  // target level
			       int chSrc,
			       int chDst,
			       unsigned options );

  template<typename data_T, typename dataIntern_T >
  void downsampleChannel( int levelMg,  // target level
			  int chSrc,
			  int chDst,
			  unsigned options = 0			  
			  );

  void downsampleVelocity( int levelMg,  // target level
			   int chSrcHi,
			   int chDst );
 
  void buildGaussianPyramidNew( 
    		int chan,
    		int levelMgMax=-1,
		unsigned options = 0 );			  

  template<typename Data_T, typename DataInt_T >
  void buildGaussianPyramid( 
    		int chan,
    		int levelMgMax=-1,
		unsigned options = 0 );			  

  void freeGaussianPyramid( int chan );
 
  int getCoarseSDFLevel( void ) const
  {
    return getNumLevels()-1;
  }
  void getVelocityFieldSplineDerivs( 
      		BlockIterator& bit,
		unsigned options,
		// OUT
		double *dvx,  
		double *dvy,
		double *dvz );

  void getVelocityFieldDerivs( BlockIterator& bit,
      			unsigned options,
			// OUT
			double *dvx,  
			double *dvy,
			double *dvz );

  Vec3Float getScalarFieldDerivs( BlockIterator& bit, unsigned opt=0 );

  Vec3Float getScalarFieldSplineDerivs( BlockIterator& bit,
			unsigned opt=0 );

  Vec3Float computeCellFaceAreaFractionsGhost(
					int levelMg, int level, int level0,
					const Vec3Int& ipos,
					SparseGrid<float> **sgFaceAreaHi,
					SparseGrid<float> *sgObstHi,
					int isGhostBlock,
					int doVoxelized = FALSE );

  float getFaceAreaRightDomBorder( int domBc ) const 
  { 
    return domBc != DBC_OPEN ? 0.0f : 1.0f;
  }

  float getFaceCoeffRightDomBorder( CellFlags flags, int domBc ) const 
  { 
    UT_ASSERT0(FALSE);
    return 0;
  }

  float getVelCompDomBorder( 
      SparseGrid<Vec3Float> *sgVel, int x, int y, int z, // center coords
      int iDir, int iSide, int domBc ) const 
  { 
    float v = sgVel->getValueGen_(x,y,z)[iDir];
    return  iSide ? domBc == DBC_OPEN ? v :
					0 :		
		    v;
  }
  
  /*-----------------------------------------------------------------------------------*/
  /* 										       */
  /*-----------------------------------------------------------------------------------*/
  // 
  // Slow cell neighborhood accessor for reference & diagnostics
  //
  //			y(iDir=1)
  //					z(iDir=2)
  //
  //			4	6
  //			
  //
  //		1	0	2	x(iDir=0)
  //
  //
  //		5	3
  //
  //			
  enum
  {
    CNI_OPT_CELLPOS = 1<<0,  
    CNI_OPT_NO_FACEAREA = 1<<1,
    CNI_OPT_HALF_NEIGHBORHOOD = 1<<2,
    CNI_OPT_FINE_FINE_ONLY = 1<<3
  };

  template < int options_T  > 
  class CellNeighborhoodInfo
  {    
    public:

    typedef struct
    {
      TransResType transRes;
      CellFlags flags[4];  
      float faceCoeff[4];   // Face coefficient ( to be scaled by 
      			    //   h * ( 1 | DPMR_GRAD_COEFF_CF | DPMR_GRAD_COEFF_FC )	
      				
      float faceArea[4];  // Face area fractions wrt. to solid obstacles
      float faceVel[4];	// Face velocity
      float faceVelAir[4];	// Face velocity in air phase
      float theta[4];	// Ghost-pressure theta
      float theta2;	// 
      float phi[4];	// Liquid surface level set
      PSFloat pressure[4];
      float genFloat[4];
      struct
      {
	int x,y,z,level;
      }
      cell[4];
    }
    Neighbor;	
    
    Neighbor neighbors[7];  //0=center 

    int _levelMg, _level, _x, _y, _z;
    
    CellNeighborhoodInfo( MultiresSparseGrid *msg, 
		  int levelMg, int level, int x, int y, int z,  // Cell location
		  int chPressure = CH_NULL, 
		  int chPhiLiquid = CH_NULL,
		  int chVelocity = CH_NULL,
		  int chFlags = CH_NULL,
		  int chGenFloat = CH_NULL
		  )
    {
      UT_ASSERT0(FALSE);
    }  
  }; // CellNeighborhoodInfo

  void traceCellNeighborhoodInfo( 
		int levelMg, int level, int x, int y, int z,  // Cell location
		unsigned opt=0, 
		CellNeighborhoodInfo<0> *cniIn=NULL
		);

  /*-----------------------------------------------------------------------------------*/
  /* 										       */
  /*-----------------------------------------------------------------------------------*/
  template<int transRes_T> 
  float getFaceAreaFractionGen( 
      		SparseGrid<float> *sgFaceArea,
      		int levelMg, int level, int x, int y, int z,  // Cell location
		int iDir, int iSide,  int x2, int y2, int z2, // Neighbor direction
		float *aOut=NULL // output (array of 4 sub faces [0..3] 
  			    // in case of coarse-fine transitions )
		)
  {  
    float a=0.0f;
    if(levelMg>=_nLevels || transRes_T == FINE_FINE )
    {
      if( iSide )
      {
	a = sgFaceArea->inRange(x2,y2,z2) ?
		    sgFaceArea->getValueGen_(x2,y2,z2) :
		    getFaceAreaRightDomBorder( domainBc(iDir,iSide));
      }
      else a = sgFaceArea->getValueGen_(x,y,z);				
      if( aOut ) aOut[0] = a;
    }
    else if( transRes_T == FINE_COARSE )
    {
      // Use fine-res 'ghost' blocks overlaying the low-res block
      a = iSide ? sgFaceArea->getValue(x2,y2,z2) :
			sgFaceArea->getValue(x,y,z);
      if(aOut) aOut[0] = a;
    }
    else if( transRes_T == COARSE_FINE )
    {
      SparseGrid<float> *sg = getFaceAreaChannel(iDir,level-1,levelMg);
      Vec3Int ip2h0 = iSide ? Vec3Int(x2,y2,z2)*2 : Vec3Int(x,y,z)*2;
      a = 0.0f;
      for(int i=0;i<2;i++)  
      for(int j=0;j<2;j++)
      {
	int idx=MT_GXY(2,i,j);
	UT_ASSERT2(idx>=0&&idx<4);
	int d[3] = {0,i,j};	
	int x2h = ip2h0.x+d[axesPermutationInv[iDir][0]],
	    y2h = ip2h0.y+d[axesPermutationInv[iDir][1]],
	    z2h = ip2h0.z+d[axesPermutationInv[iDir][2]];
	float a2 = sg->getValue(x2h,y2h,z2h);
	if(aOut) aOut[idx] = a2;
	a += a2;
      }
      a *= .25f;
    }
    else
    {
      a = -1;
      UT_ASSERT2(FALSE);
    }
    return a;
  }

  /*-----------------------------------------------------------------------------------*/
  /* 										       */
  /*-----------------------------------------------------------------------------------*/
  float getFaceAreaGen( 
		SparseGrid<CellFlags> *sgFlags,
		SparseGrid<float> **sgFaceArea,
      		int levelMg, int level, int x, int y, int z,  // Cell location
		int iDir, int iSide,  int x2, int y2, int z2 // Neighbor direction
      		)
  {
    switch( getTransResType( level, x2, y2, z2, levelMg ))
    {
      case FINE_FINE: 
	return getFaceAreaFractionGen<FINE_FINE>( sgFaceArea[iDir],levelMg,level,x,y,z,
					    iDir,iSide,x2,y2,z2 );
      case COARSE_FINE:
	return getFaceAreaFractionGen<COARSE_FINE>( sgFaceArea[iDir],levelMg,level,x,y,z,
					      iDir,iSide,x2,y2,z2 ); 		
      case FINE_COARSE:  
	return getFaceAreaFractionGen<FINE_COARSE>( sgFaceArea[iDir],levelMg,level,x,y,z,
					      iDir,iSide,x2,y2,z2 );
      default:
	UT_ASSERT2(FALSE);
	return 0;
	break;
    }
  }    

  // Misc. getter / setters 
  bool isRelaxationBlock( int levelMg, const BlockInfo *bi ) const
  {
    return !(bi->flags & (BLK_NO_FLUID|BLK_FIXED)) && 
      	    ( bi->level<=levelMg || 		      			  
	      ( (BLK_IS_RES_BORDER(bi->flags)) && 
		(bi->level==levelMg+1) ));
  }

    // theta = distance to actual interface between fluid and neighbor air cell
  int getGhostPressure_(  int levelMg, 
		         float phi, float phi2, 
			 float faceDensity, // if non-negatiive: use instead of phi,phi2
      			 float *pPressure, // IN/OUT
			 float *pTheta // OUT
      )
  {
    int rcThis=0;
    float theta=0.5;

    if(!(faceDensity<0))
    {
      float f=faceDensity;
      UT_ASSERT2(!(f<0||f>1));
      theta = f;
    }
    else
    {
      // phi is distance field
      phi = - ( phi - _liqDensThr );
      phi2 = - ( phi2 - _liqDensThr );

      float eps=1e-5;

#if 1
      if(levelMg>=_mgLevels-1)
      {
	// On very coarses MG level, the liquid SDF is no longer reliably 
	// consistent with the corasened OPEN-cell tagging, so we revert to
	// first order and place the interface at the next cell face (theta=0.5)
	theta = 0.5;
      }
      else
#endif
      {

	if(phi<0.0f && phi2>0.0f)
	{
	  float tmp = phi-phi2;
	  if( fabsf(tmp) < eps )
	  {
	    theta = 0.5;
	  }
	  else
	  {
	    theta = phi/tmp;  
	  }
	}
	else if(phi<0.0f && phi2<0.0f && phi2>phi+eps)
	{
	  // Neighbor cell still 'inside' liquid SDF despite being tagged as OPEN.
	  // This may only happen in coarser multigrid levels (due to aggressive
	  // coarsening of dirichlet cells). In this case we place the location of the
	  // interface as far as possible away using theta=1 (but can't go farther than 1 
	  // without actually decreasing the matrix diagonal)
	  if(levelMg==0)
	  {
	    TRCERR(("phi=%g, phi2=%g, levelMg=%d\n",phi,phi2,levelMg));
	    rcThis=1;
	  }
	  theta = 1.0f;
	}
	else
	{
	  #if 0
	  if( fabsf(phi-phi2)>0.25 || phi>.1 )
	  {
	    TRCWARN(("phi=%g phi2=%g levelMg=%d\n",phi,phi2,levelMg)); 
	    rcThis=1;
	  }
	  #endif
	  theta = 0.5f;
	}
      }
      UT_ASSERT2(!(theta<0.0f || theta>1.0f));
    }

    // Clamp for numerical stability
    theta = std::max(theta,_ghost_pressure_min_theta);
    if( pTheta ) *pTheta = theta;

    // Calculate the ghost pressure
    if( pPressure )
      *pPressure *= -(1.0f-theta)/theta;
    return rcThis;
  }

  bool isBoundaryRelaxationBlock( int levelMg, const BlockInfo *bi ) const
  {
    return (bi->flags & BLK_BOUNDARY_ZONE) && isRelaxationBlock( levelMg, bi );           
  }
  
  SparseGrid<Vec3Float>* sg0(void) const { return _sg0; }

  const char *getName() const { return _name; }

  float *getZeroBlock( int level ) const
  {
    UT_ASSERT2(level>=0&&level<_mgLevels);
    return _sparseGrids[level].zeroBlock;
  }

  CellFlags *getGhostFlagsBlock( int level, int domBc ) const
  {
    UT_ASSERT2(level>=0&&level<_mgLevels);    
    return domBc == DBC_OPEN ?
	      	_sparseGrids[level].ghostFlagsBlockDirichlet :
	      	_sparseGrids[level].ghostFlagsBlockNeumann;
  }

  LongInt getNumActCells( int levelMg=0 ) const { return _nActCells[levelMg]; }
  LongInt getNumFluidCells( int levelMg=0 ) const 
  {
    UT_ASSERT2(levelMg>=0&&levelMg<ARRAY_LENGTH(_nFluidCells));
    LongInt n = _nFluidCells[levelMg];
    UT_ASSERT0( n>=0 );  // -1=undefined
    return n; 
  }

  void setTraceDetail( int detailOn )
  {
    _psTraceDetail = detailOn;
  }

  int getTraceDetail( void ) const { return _psTraceDetail; }


  void setPSolverConf( WTURB_PSolverConf *conf ) 
  {
    _psConf = *conf;
    _mgBoundZoneWidth = 
      conf->ps_mg_bz_width;
    _mgSmIter = conf->ps_mg_sm_iter;
    _mgSmBoundIter = conf->ps_mg_sm_bound_iter;           
    
    _mgSmBlockIter = conf->ps_mg_sm_block_iter;           
    //TRCP(("_mgSmBlockIter = %d -> %d\n",conf->ps_mg_sm_block_iter,_mgSmBlockIter));

    _mgCoarsestMaxIter = conf->ps_mg_coarsest_sm_iter;
  }

  int domainBc( int iDir, int iSide ) const 
  {
    UT_ASSERT2(iDir>=0&&iDir<3&&iSide>=0&&iSide<2);
    int dbc = _domainBc[iDir][iSide];
    UT_ASSERT2(dbc!=DBC_INVALID);
    return dbc; 
  }

  int domainBcOpt( int iDir, int iSide ) const 
  {
    UT_ASSERT2(iDir>=0&&iDir<3&&iSide>=0&&iSide<2);
    int opt = _domainBcOpt[iDir][iSide];
    return opt; 
  }

  double h0( void ) const { return _h0; }

  // Statistics & Visualization

  size_t showMemoryUsage( int trcLevel=0 );

  const double *getVisSlicePos(void) const { return _visSlicePos; }

  int getSlices2D( 
	int chanId,
	int ipolType,
	double *slicePos0,
	double roiSize,
	double roiZoom,
	unsigned opt,
	int levelMg,
	BmpBitmap **pOutXZ,   // result in dataFloat[0..2] channel
	BmpBitmap **pOutXY,
	BmpBitmap **pOutZY
    ); 

  void getSumSlices( 
    	        int chan,
    		BmpBitmap *BZY, BmpBitmap *BXZ, BmpBitmap *BXY );

  
  void visualizeLastErr( PnlPanel *pnl, int chanId, 
    			 unsigned options );

  void visualizeSlices( 
	  PnlPanel *pnl, 
	  int chanId,
	  int ipolType=SBG::IP_NEAREST,
	  double *slicePos=NULL,
	  unsigned options=0,
	  int levelMg=0,
	  double thresholdCol=0,
	  int	visCellFlagsMask=0,
	  SparseGrid<float> *sgVis0=NULL,
	  std::function<void(
	    float val /*original value*/, MtStat *valStat,
	    int flags, int blockFlags, int level,
	    const Vec3Int& ipos, // at level 0
	    int &colorInOut )> 
	      setCellColorFunc={}
	  ) ;

  void visualizeBlockFlags( PnlPanel *pnl,
			    unsigned flagRed,
		            unsigned flagGreen,
			    unsigned flagBlue,
			    int levelMg = 0 );
  //
  // Cell visualization
  //
  enum
  {
    VISCELL_OPT_INTERPO_VEL = (1<<0),
    VISCELL_OPT_VELO_DIFF = (1<<1)
  };

  void visualizeCells( PnlPanel *pnl,
		       int opt=0,
		       int x0=-1, int y0=-1, int z0=-1, int level0=-1,
		       int roiSz=-1
			);

  int visGetCellColor( double g, double g0,
		       uint16_t cellFlags, uint16_t blockFlags,
		       int level, int level0,
		       unsigned options );

  void saveVisualizationBitmap( BmpBitmap *B, char *title );

  void setVisConf(   WTURB_VisualizationConf *visConf,
      		     int visPsOpt,
      		     int visPsMg,
		     int visTopology )
  {
    if(visConf) _visConf = *visConf;
    _visPsOpt = visPsOpt;
    _visPsMg = visPsMg;
    _visTopology = visTopology;
  }

  void setVisRoiBox( double *pos, // box center x/sx,y/sy,z/sz
      		     double size // box size/min(sx,sy,sz) 
		     )
  {
    TRC(("%s: %g %g %g %d\n",UT_FUNCNAME,pos[0],pos[1],pos[2],size));
    for(int k=0;k<3;k++) 
    {
      _visSlicePos[k] = _visSliceRoiPos[k] = pos[k];
      _visSliceRoiSize[k] = size;
    }
  }

  void getVisRoiBox( Vec3Float *pos, float *size )
  {
    if(pos) *pos = Vec3Float( _visSliceRoiPos[0],
			      _visSliceRoiPos[1],
			      _visSliceRoiPos[2] );
    if(size) *size = _visSliceRoiSize[0];				
  }

  void setVisOutDir( int visSaveAsBitmap, char *visOutDir, int doNameSubdir=0  )
  {
    TRC(("%s: %d '%s'\n",
	  UT_FUNCNAME,visSaveAsBitmap,visOutDir?visOutDir:""));
    _visSaveAsBitmap = visSaveAsBitmap;
    strcpy(_visOutDir,visOutDir?visOutDir:"");
    if(doNameSubdir) 
    {
      char tmp[UT_MAXPATHLEN+1];
      sprintf(tmp,"/MSBG_%s/",getName());
      strcat(_visOutDir,tmp);
    }
  }

  // Testing
  int runUnitTests( void );

  int testInterpolation( void );

  int testSignedDistanceField( void );

  int testBlockIterator( void );

  int testProcessChannel( void );

  void tstCreateRefinementMap( int refinementTyp,
			  int *blockmap,
			  double x0, double y0, double z0 );

  int test( int refinementTyp, int psType, int useMultigrid ); 
    
  enum
  {
    /*MIN_LEVEL = 0,
    MAX_LEVEL = MSBG_LEVELS-1,*/
    INVALID_LEVEL = 255
  };

  typedef struct
  {
    SparseGrid<Vec3Float> *vec3_1[MSBG_MAXRESLEVELS];
    SparseGrid<Vec3Float> *velocityAir[MSBG_MAXRESLEVELS];
    SparseGrid<Vec3Float> *vec3_2[MSBG_MAXRESLEVELS];  // difference for FLIP
    SparseGrid<Vec3Float> *velocityAirDiff;  // difference for FLIP
    SparseGrid<Vec3Float> *velocityAvg;  
    SparseGrid<Vec3Float> *vec3_3[MSBG_MAXRESLEVELS]; 
    SparseGrid<Vec3Float> *vec3_4[MSBG_MAXRESLEVELS]; 
    
    SparseGrid<FaceDensity> *faceDensity[MSBG_MAXRESLEVELS]; 
    SparseGrid<float> *float1[MSBG_MAXRESLEVELS],
    		      *float6[MSBG_MAXRESLEVELS];

    SparseGrid<float> *faceArea[3][MSBG_MAXRESLEVELS],
    		      *faceCoeff[3][MSBG_MAXRESLEVELS];

    SparseGrid<float> *densityDiff,
                      *sootDiff,
		      *heatDiff;
    SparseGrid<float> *float2[MSBG_MAXRESLEVELS],
		      *float4[MSBG_MAXRESLEVELS],
		      *float8[MSBG_MAXRESLEVELS],
		      *divergenceAdj,
		      *curvature,
		      *heat;

    SparseGrid<PSFloat> *divergence[MSBG_MAXRESLEVELS],
    			*diagonal[MSBG_MAXRESLEVELS],
			*pressure[MSBG_MAXRESLEVELS],
			*cgP[MSBG_MAXRESLEVELS],
		        *cgQ[MSBG_MAXRESLEVELS],
			*floatTmpPS[MSBG_MAXRESLEVELS],
			*pressureOld[MSBG_MAXRESLEVELS];

    SparseGrid<MassDensity> *massDensity[MSBG_MAXRESLEVELS];

    SparseGrid<CellFlags> *cellFlags[MSBG_MAXRESLEVELS];   // cell types
    SparseGrid<CellFlags> *cellFlagsTmp[MSBG_MAXRESLEVELS];  

    // Auxuillaries for CG solver
    SparseGrid<float> *float7[MSBG_MAXRESLEVELS],
		      *float5[MSBG_MAXRESLEVELS],
		      *float3[MSBG_MAXRESLEVELS],
		      *floatTmp3;

    // distance to fine-coarse interface    
    SparseGrid<uint16_t> *distFineCoarse,
      			 *genUint16[MSBG_MAXRESLEVELS],
			 *genUint16_2;

    SparseGrid<uint8_t> *genUint8,
      			*genUint8_2,
			*uint8Tmp,
			*uint8Tmp2;

    double kinEnergyInBand,   // Kinetic energy (laplacian-of-velocity) in this 
	   kinEnergyInBandTm,
    	   kinEnergyInBandOld,		      // level's frequency band
	   kinEnergyInBandOldTm;
    double kinEnergyInBandDiss;

    // Auxillary info
    //
    void *auxNullPtr[MSBG_MAXRESLEVELS];
    void **channelPointers[ MSBG_MAX_CHANNELS ];

    float auxInvScale;  // 1/(1<<level)

    SBG::Block<CellFlags> *constBlockFlagsNonexVoid,
	                  *constBlockFlagsNonexSolid;
    SBG::Block<float> *constBlockFaceCoeffFullAir,
		      *constBlockFaceCoeffFullLiq;
    SBG::Block<PSFloat> *constBlockDiagFullAir,
		        *constBlockDiagFullLiq;

    size_t zeroBlockSz,
	   zeroBlockPSFloatSz;
    float *zeroBlock,
	  *constOneBlock,
	  *constFaceCoeffAirBlock,
	  *sdfInfiniteBlock;
    PSFloat *zeroBlockPSFloat;
    uint16_t *sdfInfiniteBlockUint16;
    Vec3Float *zeroBlockVec3;
    CellFlags *ghostFlagsBlockNeumann,
	      *ghostFlagsBlockDirichlet;
  }
  Level;


  float getFloatValue( int chan, Vec4f pos )
  {
    CellAccessor ca = getCellAccessor( pos );
    return getFloatChannel(chan,ca.level)->getValueInBlock0(ca.bid,ca.vid);    
  }

  template< int doClampPos=0 >
  float getFloatValue( int chan, Vec4i ipos )
  {
    if constexpr ( doClampPos )
    {
      ipos = max(ipos,_mm_setzero_si128());
      ipos = min(ipos,_v4iDomMax);
    }
    CellAccessor ca = getCellAccessor0( ipos );
    return getFloatChannel(chan,ca.level)->getValueInBlock0(ca.bid,ca.vid);    
  }

  float _obstaclesNbDist,
	_nbDistCSDF;

  inline CellFlags getCellFlagsGen( SparseGrid<CellFlags> *sgFlags, 
      		             int x, int y, int z )
  {
    int x_=x, y_=y, z_=z, isClipped, iDir, iSide;
    sgFlags->clipGridCoords( x_, y_, z_, isClipped, iDir, iSide );
    if( isClipped )
    {
      UT_ASSERT2(iDir>=0&&iDir<3 &&iSide>=0&&iSide<2);
      return CELL_OUT_OF_DOMAIN | 
	 ( domainBc(iDir,iSide) == DBC_OPEN ? CELL_VOID : CELL_SOLID );
    }
    else return sgFlags->getValueGen_(x,y,z);
  }

  CellFlags getCellFlags( const Vec4f& pos )
  {
    CellAccessor ca = getCellAccessor( pos );
    SparseGrid<CellFlags> *sg = _sparseGrids[ca.level].cellFlags[0];	
    return sg->getValueInBlock0(ca.bid,ca.vid);
  }

  CellFlags getCellFlags( float x0, float y0, float z0 )
  {
    CellAccessor ca;
    getCellAccessor( x0, y0, z0, &ca);
    SparseGrid<CellFlags> *sg = _sparseGrids[ca.level].cellFlags[0];	
    return sg->getValueInBlock(ca.bid,ca.vid);
  }

  inline int isInRange( const Vec4f& pos, const Vec4f& min, Vec4f& max )
  {
    return !( _mm_movemask_ps( pos < min || pos > max )  );
  }

  inline int isInDomainRange( const Vec4f& pos )
  {
    return isInRange( pos, _domMin, _domMax );
  }

  CellFlags getCellFlags( Vec3Float pos )
  {
    return getCellFlags(pos[0],pos[1],pos[2]);
  }

  CellFlags getCellFlags( float *pos )
  {
    return getCellFlags(pos[0],pos[1],pos[2]);
  }

  bool isUint8Channel( int ch ) const 
  { 
    return (((int64_t)1)<<((int64_t)ch)) & ( (((int64_t)1)<<CH_UINT8_1) | 
	                 	  	     (((int64_t)1)<<CH_UINT8_2) |
	                 	  	     (((int64_t)1)<<CH_UINT8_TMP ) |
	                 	             (((int64_t)1)<<CH_UINT8_TMP_2 ) ); 
  }

  bool isUint16Channel( int ch ) const 
  { 
    return (((int64_t)1)<<((int64_t)ch)) & ( (((int64_t)1)<<CH_CELL_FLAGS) | 
	                 	  	     (((int64_t)1)<<CH_CELL_FLAGS_TMP) |
	                 	  	     (((int64_t)1)<<CH_UINT16_1 ) |
	                 	             (((int64_t)1)<<CH_UINT16_2 ) ); 
  }

  bool isPSFloatChannel( int ch ) const 
  { 
    return (((int64_t)1)<<((int64_t)ch)) & ( (((int64_t)1)<<CH_PRESSURE) |
					     (((int64_t)1)<<CH_PRESSURE_OLD) |
	                 	  	     (((int64_t)1)<<CH_CG_Q) |
	                 	  	     (((int64_t)1)<<CH_CG_P ) |
	                 	  	     (((int64_t)1)<<CH_DIVERGENCE) |
	                 	  	     (((int64_t)1)<<CH_DIAGONAL ) |
	                 	             (((int64_t)1)<<CH_FLOAT_TMP_PS ) ); 
  }

  int isFlagsChannel( int ch ) const { return isUint16Channel(ch); }

  bool isVecChannel( int ch ) const 
  { 
    return (((int64_t)1)<<ch) & 
      		( (((int64_t)1)<<CH_VEC3_1) | 
		  (((int64_t)1)<<CH_VEC3_3) | 
		  (((int64_t)1)<<CH_VEC3_2) | 
		  (((int64_t)1)<<CH_VELOCITY_AVG) | 
		  (((int64_t)1)<<CH_VELOCITY_AIR) | 
		  (((int64_t)1)<<CH_VELOCITY_AIR_DIFF) | 
		  (((int64_t)1)<<CH_VEC3_4) | 
		  (((int64_t)1)<<CH_FACE_DENSITY) );
  }
  
  SparseGrid<uint16_t> *getDistFineCoarseChannel( int level=0 ) const
  {
    if(!(level>=0&&level<_mgLevels)) return NULL;
    UT_ASSERT2( ! ( (_protectedRead & (((int64_t)1)<<((int64_t)CH_DIST_FINECOARSE)) ) ||
	            (_protectedWrite & (((int64_t)1)<<((int64_t)CH_DIST_FINECOARSE)) )));
    return _sparseGrids[level].distFineCoarse;
  }

  SparseGrid<CellFlags> *getFlagsChannel0(int level) const
  {
    return _sparseGrids[level].cellFlags[0];
  }

  SparseGrid<CellFlags> *getFlagsChannel0( int level, int levelMg ) const
  {
    return _sparseGrids[level].cellFlags[levelMg<_nLevels?levelMg:0];	  
  }

  void resetEMPTYFlags( std::vector<int> *blockList );

  SparseGrid<float> *getFaceAreaChannel( int iDir, 
      					     int level, int levelMg=0 ) const
  {
    UT_ASSERT2(iDir>=0&&iDir<3);
    if(!(level>=0&&level<_mgLevels)) return NULL;
    return _sparseGrids[level].faceArea[iDir][levelMg<_nLevels?levelMg:0];	  

  }

  SparseGrid<float> *getFaceCoeffChannel( int iDir, 
      					     int level, int levelMg=0 ) const
  {
    UT_ASSERT2(iDir>=0&&iDir<3);
    if(!(level>=0&&level<_mgLevels)) return NULL;
    return _sparseGrids[level].faceCoeff[iDir][levelMg<_nLevels?levelMg:0];	  

  }

  SparseGrid<uint8_t> *getUint8Channel( int chanId, int level, 
      					  int levelMg=0 ) const
  {
    if(!(level>=0&&level<_mgLevels)) return NULL;
    UT_ASSERT2(chanId>=0&&chanId<MSBG_MAX_CHANNELS);

    UT_ASSERT2( ! ( (_protectedRead & (((int64_t)1)<<chanId) ) ||
	            (_protectedWrite & (((int64_t)1)<<chanId) )));

    int levelMgAct = levelMg<_nLevels?levelMg:0;
    return (SparseGrid<uint8_t> *)
      		_sparseGrids[level].channelPointers[chanId][levelMgAct];
  }

  void setUint8Channel( SparseGrid<uint8_t> *sg,
      int chanId, int level, 
      int levelMg=0 )
  {
    if(!(level>=0&&level<_mgLevels)) return;
    UT_ASSERT2(chanId>=0&&chanId<MSBG_MAX_CHANNELS);

    UT_ASSERT2( ! ( (_protectedRead & (((int64_t)1)<<chanId) ) ||
	            (_protectedWrite & (((int64_t)1)<<chanId) )));

    int levelMgAct = levelMg<_nLevels?levelMg:0;
    _sparseGrids[level].channelPointers[chanId][levelMgAct] = sg;
  }


  SparseGrid<CellFlags> *getFlagsChannel( int chanId, int level, 
      					  int levelMg=0 ) const
  {
    #if 1
    if(!(level>=0&&level<_mgLevels)) return NULL;
    #else
    UT_ASSERT2( (level>=0&&level<_mgLevels) );
    #endif
    UT_ASSERT2(chanId>=0&&chanId<MSBG_MAX_CHANNELS);

    UT_ASSERT2( ! ( (_protectedRead & (((int64_t)1)<<chanId) ) ||
	            (_protectedWrite & (((int64_t)1)<<chanId) )));

    int levelMgAct = levelMg<_nLevels?levelMg:0;
    return (SparseGrid<uint16_t> *)
      		_sparseGrids[level].channelPointers[chanId][levelMgAct];

  }

  SparseGrid<CellFlags> *getUint16Channel( int chanId, int level, 
      					    int levelMg=0 ) const
  {
    return getFlagsChannel(chanId,level,levelMg);
  }

  void **getChannelAddr( int chanId, int level=0, int levelMg=0 );
  
  SparseGrid<float> *getFloatChannel( int chanId, int level, int levelMg=0 )
    const
  {
    #if 1
    if(!(level>=0&&level<_mgLevels)) return NULL;
    #else
    UT_ASSERT2( (level>=0&&level<_mgLevels) );
    #endif

    UT_ASSERT2(chanId>=0&&chanId<MSBG_MAX_CHANNELS);

    #ifdef SOLVE_PRESSURE_DOUBLE
    UT_ASSERT2( UT_IMPLIES( chanId, !isPSFloatChannel( chanId )));
    #endif

    UT_ASSERT2( ! ( (_protectedRead & (((int64_t)1)<<((int64_t)chanId)) ) ||
	            (_protectedWrite & (((int64_t)1)<<((int64_t)chanId)) )));

    int levelMgAct = levelMg<_nLevels?levelMg:0;
    return (SparseGrid<float> *)
      		_sparseGrids[level].channelPointers[chanId][levelMgAct];
  }


  void *getGenChannel( int chanId, int level, int levelMg=0 )
    const
  {
    if(!(level>=0&&level<_mgLevels)) return NULL;
    UT_ASSERT2(chanId>=0&&chanId<MSBG_MAX_CHANNELS);
    UT_ASSERT2( ! ( (_protectedRead & (((int64_t)1)<<((int64_t)chanId)) ) ||
	            (_protectedWrite & (((int64_t)1)<<((int64_t)chanId) ) )));
    int levelMgAct = levelMg<_nLevels?levelMg:0;
    return _sparseGrids[level].channelPointers[chanId][levelMgAct];
  }

  SparseGrid<PSFloat> *getPSFloatChannel( int chanId, int level, int levelMg=0 )
    const
  {
    #ifdef SOLVE_PRESSURE_DOUBLE
    UT_ASSERT2( UT_IMPLIES( chanId, isPSFloatChannel(chanId) ));
    #endif
    return (SparseGrid<PSFloat> *)getGenChannel( chanId, level, levelMg );
  }

  float *getFloatBlockDataPtr( int chan, int bid, int level, int levelMg,
      			  int doAlloc, int doZeroDataOnInit ) 
  {
    UT_ASSERT2((doAlloc==0||doAlloc==1) && 
	       (doZeroDataOnInit==0||doZeroDataOnInit==1));

    SparseGrid<float> *sg=getFloatChannel(chan,level,levelMg);
    
    return sg ? sg->getBlockDataPtr(bid,doAlloc,doZeroDataOnInit) : NULL;
  }

  SparseGrid<Vec3Float> *getVecChannel( int chanId, int level, int levelMg=0 )
    const
  {
    if(!(level>=0&&level<_mgLevels)) return NULL;
    
    UT_ASSERT2( ! ( (_protectedRead & (((int64_t)1)<<((int64_t)chanId)) ) ||
	            (_protectedWrite & (((int64_t)1)<<((int64_t)chanId) ) )));

    int levelMgAct = levelMg<_nLevels?levelMg:0;
    return (SparseGrid<Vec3Float> *)
      		_sparseGrids[level].channelPointers[chanId][levelMgAct];
  }
  
  void protectChannel( int iChan, bool protRead=1, bool protWrite=1 );

  /*-----------------------------------------------------------------------*/
  /* 								           */
  /*-----------------------------------------------------------------------*/
  template <typename Data_T>
  void syncChannelsDataGenCount( int chSrc, int levelMgSrc, 
      			         int chDst, int levelMgDst )
  {
    TRC3(("%s %d %d %d %d\n",UT_FUNCNAME,chSrc,levelMgSrc, chDst,levelMgDst));
    UT_ASSERT0( ABS(levelMgSrc-levelMgDst) == 1 );
    for(int level=levelMgDst; level<_nLevels; level++)
    {
      if(level>=levelMgSrc)
      {
        getChannel<Data_T>(this,chDst,level,levelMgDst)->setDataGenCount( 
	   getChannel<Data_T>(this,chSrc,level,levelMgSrc)->getDataGenCount());
      }
    }
  }

  /*-----------------------------------------------------------------------*/
  /* 								           */
  /*-----------------------------------------------------------------------*/
  void resetChannel( int ch=-1, int ignoreNotFound=FALSE, int levelMg=-1,
		     int opt=0 /* &1=reset, 
				  &2=release physically,
				  &4=prepareDataAccess-read 
				  &8=prepareDataAccess-write */ );

  /*-----------------------------------------------------------------------*/
  /* 								           */
  /*-----------------------------------------------------------------------*/
  void prepareDataAccess( int ch, int level=-1, int levelMg=-1,
		          unsigned accType = SBG::ACC_READ|SBG::ACC_WRITE )
  {
    UT_ASSERT0(level==-1);
    UT_ASSERT0(accType==(SBG::ACC_READ|SBG::ACC_WRITE));
    resetChannel( ch, 0, levelMg, (accType & SBG::ACC_READ ? 4 : 0 ) |
			            (accType & SBG::ACC_WRITE ? 8 : 0 ) );	
  }


  void multiplyLaplacianMatrixOpt(  
			  unsigned options, 
			  int levelMg,
			  int chanX,
			  int chanB,
			  float diffusion,
			  int chanDstY,
			  long double* pDotprod,
			  std::vector<int> *blockList=NULL			  
			  );

  void ImultiplyLaplacianMatrixDenseLevel(  
			    unsigned options, 
			    int levelMg,
			    int chanX0,
			    int chanB,
			    int chanTmp,
			    int nJacobiIter,
			    int chanDstY0
			    );

  template< int opmode, // LAPL_RELAX, LAPL_MULTIPLY
	    int bsx, 
	    int SIMDW, 
	    typename VecNpfs, typename VecNf, typename VecNui >
  void IprocessBlockLaplacian(
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
		      PSFloat *pMaxScaledResid = NULL
		      );

  template< int opmode /* LAPL_RELAX, LAPL_MULTIPLY */ >
  void processBlockLaplacian(
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
		      PSFloat *pMaxScaledResid = NULL
		    );


  void relaxBlock(
      		    int tid,	// Thread number (-1=init buffers)
		    unsigned options,
    		    int bid, int bx, int by, int bz,
		    BlockInfo *bi,
		    int levelMg,
		    int chanX,
		    int chanB,
		    int chanDstY, 
		    int nIter
		    );

  void ImultiplyLaplacianMatrixOpt(  
			    unsigned options, 
			    int levelMg,
			    int chanX,
			    int chanB,
			    float kap,
			    int chanDstY,
			    long double* pDotprod,
			    std::vector<int> *blockList=NULL			  
			    );

  void multiplyLaplacianMatrixRef_(  unsigned options,
      				    int levelMg,
				    int chanX,
				    int chanB,
				    double omega,
				    double kap,
				    int chanDstY,
				    long double* pDotprod, 
			    	    std::vector<int> *blockList=NULL			  
				    );

  int verifyMatrixSymmetryRef( int levelMg, std::vector<int> *blocklist=NULL );

  int verifyMatrixSymmetry( int matTyp );

  LongInt relaxBlockList(  std::vector<int> *blocks,
      			  int levelMg,
			  int iIter,
			  int chanX,  // source channel
			  int chanB,  // right hand side
			  unsigned options,
			  int chanDstY // destination channel				
			  );

  void relax(  int boundaryZoneOnly,
      		  int reverseOrder,
      			     int levelMg,     					    
			     int chanX,
			     int chanB,
			     int chanTmp,
			     int nIter,
			     unsigned options=0
			     );

  void setChannel( double a, int iChanDstY,
      		    int levelMax=-1, int levelMg=0, CellFlags cellMask=0,
		    std::vector<int> *blockList=NULL,
		    unsigned skipBlockFlags=0
		    );

  template< typename Data_T >
  void setChannelZero(  int iChan, int levelMg=0 );

  void checkChannelNonFluidZero( int chan, int levelMg );

  template< typename Data_T >
  void setChannelConstBlock( SBG::ConstVal constVal, int iChanDstY, 
      			     int levelMg=0);
 
  void sqrtChannel( int iChanDst );

  void copyChannel( int chSrc, int chDst, 
      		    int includeDomBorder=1, int levelMax=-1, int levelMg=0, 
		    unsigned options=0,
                    std::vector<int> *blockList=NULL
		    );

  void copyChannelFromMSBG( const MSBG::MultiresSparseGrid *msgSrc,
    			    int chSrc, int chDst );

  void scalarChannelFromToVecComp( 
    			int dir,  // direction 1=vec2scalar, -1=scalar2vec
			int iVecComp,
    			int chScalar, int chVec );

  void multChannel(  int iChanX1, int iChanX2,
    		     int iChanDstY );

  void scaleChannel(  int iChanX, double a, double b,
    		      int iChanDstY,
		      std::vector<int> *blockList=NULL );

  int verifyConsistency( unsigned opt );

  void visualizeResidual( 
    			PnlPanel *pnl,
    			int chRes,  // Residual
			int chVisTmp, // dest. channel for visualization
			int levelMg
			);

  double compareChannel( int chanA, int chanB, int levelMg=0, int doVisu=0, 
      			 int doOnlyLevel0=FALSE, int chVisDiff=CH_NULL,
			 double errThr = 1e-5 );

  template< typename Float_T >
  long double dotProdChannel_( int ch1, int ch2, int levelMg=0, 
      			       std::vector<int> *blockList=NULL, int *pErrCode=NULL );
  inline long double dotProdChannel( int ch1, int ch2, int levelMg=0, 
      				     std::vector<int> *blockList=NULL, int *pErrCode=NULL )
  {
    UT_ASSERT2( !!isPSFloatChannel(ch1) == !!isPSFloatChannel(ch2) );
    return isPSFloatChannel(ch1) ?
      dotProdChannel_<PSFloat>( ch1, ch2, levelMg, blockList, pErrCode ) :
      dotProdChannel_<float>( ch1, ch2, levelMg, blockList, pErrCode );

  }

  void absChannel(  int iChan );
  void checkFloatChannel(  int iChan, int levelMg );

  Vec3Double computeSyntheticVelocity( 
    	          Vec3Double pos /* position in grid coords [0..smax] */ );
  
  void fillTestPatternChannel( int typ,
			       int chanDst, 
			       int levelMg=0,
			       double rndSeed=0,
			       unsigned options=0
			       );

  void AXPBYChannel(  double a, int iChanX, 
      		      double b, int iChanDstY,
		      int levelMax = 9999, 
		      int levelMg = 0,
		      std::vector<int> *blockList = NULL
		      );

  void AXPBYChannelCombined(  double a, double b,    						
			      int iChanQ,
			      int iChanX,
			      int iChanP,
		              std::vector<int> *blockList = NULL
			);

  void MultiplyLaplacianMatrix(  int iChanX,
    				 int iChanDstY,
				 long double *pAlpha);

  void scaleByCellVolume( int chan );

  template< typename Float_T >
  double maxAbsChannel_(  int iChan, int level0=9999, int levelMg=0,
      			 int doScaleByCellVolume=FALSE,  
			 Vec3Int *outPosMax=NULL,
			 int *outLevelMax=NULL,
		         std::vector<int> *blockList=NULL );

  inline double maxAbsChannel(  int iChan, int level0=9999, int levelMg=0,
      			 int doScaleByCellVolume=FALSE,  
			 Vec3Int *outPosMax=NULL,
			 int *outLevelMax=NULL,
		         std::vector<int> *blockList=NULL )
  {
    return isPSFloatChannel(iChan) ? 
      maxAbsChannel_<PSFloat>( iChan, level0, levelMg, doScaleByCellVolume,
	  		       outPosMax, outLevelMax, blockList ) :
      maxAbsChannel_<float>( iChan, level0, levelMg, doScaleByCellVolume,
	  		       outPosMax, outLevelMax, blockList );
  }

  void statChannel(  int iChan, 
		      int doCountFluidCellsOnly,
		      double zeroEps,
		      MtStat *stat, // OUT
		      std::vector<int> *blockList=NULL );

  void normalizeChannel(  int chan );

  int _nLevels;  // Number of resolution levels
  Level _sparseGrids[MSBG_MAX_LEVELS];
  //void *channels[MSBG_MAX_CHANNELS][MSBG_MAX_LEVELS];
  SparseGrid<Vec3Float> *_sg0;

  // Block lists (precomputed aux.)
  std::vector<int> 
  		   _blocksExisting, 
		   *_blocksActive,
		   _blocksFineCoarsePerLevel[MSBG_MAXRESLEVELS],
		   _blocksFluid[MSBG_MAXRESLEVELS],  
		   _blocksNonFrozen[MSBG_MAXRESLEVELS],  
		   _blocksRelax[MSBG_MAXRESLEVELS],  
		   _blocksBoundaryRelax[MSBG_MAXRESLEVELS],
  		   _blocksAdaptiveRelax[MSBG_MAXRESLEVELS][MSBG_AR_MAX_ITER],
		   _blocksValue[MSBG_MAXRESLEVELS];

  int _mgArActNumIter[MSBG_MAXRESLEVELS];

  void sortBlockListMorton( std::vector<int>& blockList );

  int getNumLevels( void ) const { return _nLevels; }
  int getNumMgLevels( void ) const { return _mgLevels; }
  bool isDenseLevel( int levelMg ) const { return levelMg >= _nLevels; } 

  int nBlocks( int level=0 ) const 
  { 
    return (level < _nLevels ? _nBlocks : 1); 
  }

  LongInt nActCells( int levelMg = 0) const
  {
    return _nActCells[levelMg];
  }

  double gridSpacing( int level=0 ) const
  {
    return (1<<level)*_dx0;
  }

  float scaleToNormCoords() const { return _scaleToNormCoords; }

  int dTransRes() { return _dTransRes; }

  TransResType getTransResType(
    		int level0,   // level of current cell
    		int x2, int y2, int z2, // coords of neighbor cell
		int levelMg = 0
    		)
  {
    if( levelMg>=_nLevels ) 
      return FINE_FINE; // dense coarse grids

    // Get neighbor block coords
    int bsxShift = _sg0->bsxLog2() - level0;
    UT_ASSERT2(bsxShift>=0);
    int bx2 = x2 >> bsxShift,
	by2 = y2 >> bsxShift,
	bz2 = z2 >> bsxShift;

    if( !(_sg0->blockCoordsInRange(bx2,by2,bz2))) return FINE_FINE;

    LongInt bid = _sg0->getBlockIndex(bx2, by2, bz2);
    UT_ASSERT2(bid>=0&&bid<_sg0->nBlocks());
    int level = MAX( _blockmap[0][bid].level, levelMg);

    //UT_ASSERT2(level>=0&&level<_nLevels);

    if(level==level0 || level>=_nLevels) 
    {
      return FINE_FINE;
    }
    else if(level==level0+1)
    {
      return FINE_COARSE;
    }
    else if(level==level0-1)
    {
      return COARSE_FINE;
    }
    else
    {
      TRCERR(("L=%d l=%d l2=%d pos2=%d,%d,%d\n",
	    levelMg,level0,level,x2,y2,z2));
      //return INVALID;
      UT_ASSERT0(FALSE);
    }
  }

  int getTimeLevel( int bid ) const
  {
  //UT_ASSERT2( bid>=0 && bid < _nBlocks);
    return _blockTimeLevels[bid];
  }

  int getBlockLevel( int bid, int levelMg=0 ) const
  {
  //UT_ASSERT2( bid>=0 && bid < _nBlocks);
    return levelMg < _nLevels ? 
      MAX( _blockmap[0][bid].level, levelMg ) : 
      levelMg;
  }

  int getBlockLevel0( int bid ) const
  {
    UT_ASSERT2( bid>=0 && bid < _nBlocks);
    return _blockmap[0][bid].level; 
  }

  void getBlockLevel( int bid, int levelMg,
      		     int& level0, 
		     int& level ) 
    		     const
  {
    if(levelMg<_nLevels)
    {
      UT_ASSERT2( bid>=0 && bid < _nBlocks);
      level0 = _blockmap[0][bid].level;
      level = MAX( level0, levelMg );
    }
    else
    {
      level0 = 0;
      level = levelMg;
    }  
  }

  // 
  // 'pseudo-private' members sparing the hassle with accessor functions
  //
  int	_doDumpSimulationState;

  MtStat _velMagStatLiq,
  	 _velMagStatLiqPrev;
  MtStat _velMagStatLiqFrame,
  	 _velMagStatLiqFramePrev;

  PSFloat _pressureMin,
  	  _pressureMax;

  float _avgMassDensity0;

  int _chAdvectorVelocity;

  WTURB_PSolverConf _psConf;
  WTURB_MiscConf _miscConf;

  float _ghostFloatSpecVals3x3x3[3*3*3];     

  const Vec4f _v4fStencil7[7] = 
	  { Vec4f( 0.f, 0.f, 0.f, 0.f),
	    Vec4f(-1.f, 0.f, 0.f, 0.f), Vec4f(1.f, 0.f, 0.f, 0.f),
	    Vec4f( 0.f,-1.f, 0.f, 0.f), Vec4f(0.f, 1.f, 0.f, 0.f),
	    Vec4f( 0.f, 0.f,-1.f, 0.f), Vec4f(0.f, 0.f, 1.f, 0.f) };

  const Vec4i _v4iStencil7[7] = 
	  { Vec4i( 0, 0, 0, 0),
	    Vec4i(-1, 0, 0, 0), Vec4i(1, 0, 0, 0),
	    Vec4i( 0,-1, 0, 0), Vec4i(0, 1, 0, 0),
	    Vec4i( 0, 0,-1, 0), Vec4i(0, 0, 1, 0) };

  const Vec4f _v4fOffsMAC[4] = 
  	   { Vec4f(0.0f, 0.5f, 0.5f, 0.0f),  // x-face
	     Vec4f(0.5f, 0.0f, 0.5f, 0.0f),  // y-face
	     Vec4f(0.5f, 0.5f, 0.0f, 0.0f),  // z-face
	     Vec4f(0.5f, 0.5f, 0.5f, 0.0f)   // cell centered	   
	   };

  int _chDivAdjValid;

  float _ghost_pressure_min_theta,
  	_gp_theta_out_of_dom_liq,
	_gp_theta_out_of_dom_air;

  int 	_liqDropOn,
	_liqDropInitvelPressureCorrected;
  int   _liqDropMult;
  float _liqDropMultDisp,_liqDropMultDispVel;
  unsigned _liqDropOpt;
  float _liqDropDist, _liqDropDens, _liqDropBlendVel;
  int 	_liqDropKill;
  float _liqDropKillDist, _liqDropKillDens, _liqDropKillLifetime;
  int 	_liqDropRejoin;
  float _liqDropRejoinDist, _liqDropRejoinBlendVel, _liqDropRejoinBlendWidth;
  int	_liqDropDragOn;
  float _liqDropDragStrength;

  float _liqNbWidth,
  	_liqNbWidthAir;

  int _adaptiveOpenAirBc;

  float _phiParticleLevel[MSBG_MAXRESLEVELS],
  	_phiParticleLevelMin, _phiParticleLevelMax;
  float _phiParticleLevelAir[MSBG_MAXRESLEVELS],
  	_phiParticleLevelAirMin, _phiParticleLevelAirMax;

  SparseGrid<float> *_sgCoarseDistToRefineRoi;  
  int _inflowOptions;

  int _refineAtObst;
  int getRefineAtObst(void) const { return _refineAtObst; }

  double _scaleToNormCoords;
  double _totalDensitySum;

  double _velMagMax,
  	 _velMagMaxNonOutliers;
  double getVelMagMaxNonOutliers(void) const { return _velMagMaxNonOutliers; }

  double _enstrophy,
	 _enstrophyInitial,
	 _enstrophyOld,
	 _enstrophyOldTm,
	 _enstrophyDissExp,
	 _enstrophyTimeSum;
  double _energy,
	 _energyInitial,
	 _energyOld,
	 _energyOldTm,
	 _energyDiss,
	 _energyInjected;

  long double _totalLiquidVolume,
	      _totalParticleMass,
	      _totalParticleMass0;

  // Statistics & Diagnostics
  LongInt  _nActCells[MSBG_MAX_LEVELS],
	   _nFluidCells[MSBG_MAX_LEVELS],
	   _nSolidCells[MSBG_MAX_LEVELS],
	   _nVoidCells[MSBG_MAX_LEVELS];
  LongInt  _nActCellsBlkFluid;

  int _lastErrNum,
      _lastErrNumVisualized,
      _lastErrX,
      _lastErrY,
      _lastErrZ,
      _lastErrLevel,
      _lastErrLevelMG;

  // debug

  QuadSplineWeights *getQuadBSplineWeights( int i )
  {
    return &_quadSplineWgtStaggered_[i];
  }

  private:

  // Domain border offsets:  Define enough ghost cells 
  // to accomodate a stencil of radius 1
  enum { DB = 1 };

  char	    _name[80];

  int64_t	_validForInterpolation,
	        _protectedRead,
	        _protectedWrite;
  int _options;

  int _totalSteps;
  double _totalTime,
  	 _totalTimePrev,
  	 _totalTimeFrame,
	 _totalTimeFramePrev;
  double _dt,
  	 _dtPrev;
  
  int 	 _nBlocks;
  float   _goffVel[3][3];

  float _faceCoeffFullAir,
        _faceCoeffFullLiq;
  PSFloat _diagFullAir,
	  _diagFullLiq;

  double _dx0,  // grid spacing at finest resolution
	 _h0;
 
  float _gravityAlpha;
  Vec3Float _gravityDir;

  int _coarseDistInterfaceLevel;

  void setBlockLevel( int bid, int level )
  {
    UT_ASSERT0( bid>=0 && bid<_nBlocks);
    UT_ASSERT0( level>=0 && (level<_nLevels || (level==INVALID_LEVEL)));
    _blockmap[0][bid].level = level;
  }

  // 
  // Domain boundary conditions  
  // X:left-right, Ybottom-top, Z:front-back
  //
  int _domainBc[3][2],      
      _domainBcOpt[3][2],
      _domainBcFlags[3][2];

  Vec4f _domMin,
        _domMax;

  double _domainBcInflow[3];

  float _obsVelocity[4];

  int _dTransRes;

  int   _mgOn,
	_mgUseFaceArea;
  int   _mgLevels;
  int 	_mgRestrictTyp,
  	_mgBoundZoneWidth;
  float _mgSmOmega,
  	_mgSmOmegaDefault,
  	_mgSmOmegaSched1,
	_mgSmOmegaSched2;
  int 	_mgSmType,
	_mgSmIter,
	_mgSmBlockIter,
	_mgSmRelaxType,
	_mgSmBoundIter,
	_mgVcycPerCgIter;
  int 	_mgSxyzMin,
	_mgCoarsestMaxIter;

  double _mgMaxRInitial, 
  	 _mgMaxRScaledInitial,
	 _mgCgActResid,
	 _mgCgActResidHistMean,
	 _mgCgActResidHistThr;
  int    _mgCgActIter;   
  double _mgArResid0;
  double _mgCgAvgConvRate,
  	 _mgCgAvgConvRateMovAvg;

  int   _psTraceDetail;

  //MultiresSparseGrid* _mgMsbg[MSBG_MAX_LEVELS];

  uint32_t *_blockListTmp;
  LongInt  _blockListTmpMaxLen;

  // Visualization
  
  void visualizeDenseField( 
    	PnlPanel *pnl, 
	float **velField,  // dense grid data 
	uint8_t *obstacles,
	double *slicePos,
	int doVisDistField
      );
  
  int	_visDcmpVel;

  int _visPsOpt,
      _visPsMg,
      _visTopology;
    
  double _visSlicePosAndRoi[3+3+3],
	 *_visSlicePos,
	 *_visSliceRoiPos,
	 *_visSliceRoiSize;

  WTURB_VisualizationConf _visConf;

  int    _visSaveAsBitmap;
  char   _visOutDir[UT_MAXPATHLEN+1];

  // 
  // Misc. auxillary data
  //    
  QuadSplineWeights _quadSplineWgtCellCenter,
		    _quadSplineWgtStaggered_[3];

  // Linear downsampling kernel (cell centered & staggered)
  float _krnDownSample4x4x4[4][4*4*4];

  // Auxillary data
  int _bsx0,
      _bsx0Log2;
  Vec4i _v4iDomMax;
  //Vec4i _auxBlkDims;

  int	_blkStridY,
	_blkStridZ;

  int	   _hasNullspace;

  friend class BlockIterator;

};  // MultiresSparseGrid

/*-----------------------------------------------------------------------------------*/
/* 										       */
/*-----------------------------------------------------------------------------------*/
template< typename SimdType >
void setBcSpecValsSIMD( uint16_t *dataFlags, uint16_t mskNeumann, uint16_t mskDirichlet,
			SimdType& X ) = delete;

template<> inline
void setBcSpecValsSIMD( uint16_t *dataFlags, 
		uint16_t mskNeumann, uint16_t mskDirichlet, 
		Vec4f& X )
{
  Vec8us M8us = Vec8us().load( dataFlags  );
  Vec4ui Mi = _mm_unpacklo_epi16(M8us, _mm_set1_epi16(0));  
  X = select(  ( Mi & mskNeumann ) != 0, SPEC_FLOAT_NEUMANN, X );
  X = select(  ( Mi & mskDirichlet ) != 0, SPEC_FLOAT_DIRICHLET , X );
}

template<> inline
void setBcSpecValsSIMD( uint16_t *dataFlags, 
		uint16_t mskNeumann, uint16_t mskDirichlet, 
		Vec8f& X )
{
  Vec8us M8us = Vec8us().load( dataFlags  );
  Vec4ui MiL = _mm_unpacklo_epi16(M8us, _mm_set1_epi16(0)),
	 MiH = _mm_unpackhi_epi16(M8us, _mm_set1_epi32(0));
  Vec8i Mi = Vec8i( MiL, MiH );
  X = select(  ( Mi & mskNeumann ) != 0, SPEC_FLOAT_NEUMANN, X );
  X = select(  ( Mi & mskDirichlet ) != 0, SPEC_FLOAT_DIRICHLET , X );
}


float lineFractionInside(float phi_left, float phi_right);
float areaFractionInside(float phi_bl, float phi_br, float phi_tl, float phi_tr) ;

static void permuteCoords( const int *perm, int& x, int& y, int& z )
{
  int xyz[3];
  xyz[perm[0]] = x;
  xyz[perm[1]] = y;
  xyz[perm[2]] = z;
  x = xyz[0];
  y = xyz[1];
  z = xyz[2];
}

static const int dirNeigh_[6][4] = 
	{ {1,0,0, 1}, {-1,0,0, -1},
	  {0,1,0, 1}, {0,-1,0, -1},
	  {0,0,1, 1}, {0,0,-1, -1} };

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
template<typename T>
inline T max3( const T& a, const T& b, const T& c )
{
  T largest = a;
  if (largest < b) largest = b;
  if (largest < c) largest = c;
  return largest;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
inline 
float distToBoxSq( Vec4f bmin,  // box min-corner
    		   Vec4f bmax,
		   Vec4f pos )
{
  bmin -= pos;
  pos -= bmax;

  bmin = max(max(bmin,0),pos);

  bmin *= bmin;

  return vfget_x(bmin)+vfget_y(bmin)+vfget_z(bmin);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
inline 
double distToBoxSq( double x0, double y0, double z0,  // box min-corner
    			 double x1, double y1, double z1,  // box max corner
			 double px, double py, double pz  // point
			 )
{
  double dx = max3<double>(x0-px,0,px-x1),
  	 dy = max3<double>(y0-py,0,py-y1),
	 dz = max3<double>(z0-pz,0,pz-z1);
  return dx*dx+dy*dy+dz*dz;
}

} // namespace MSBG

#endif  // __cplusplus

#endif // MSBG_H
