/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#ifndef MSBGCELL_H
#define MSBGCELL_H

// 
// Global compile time control switches
//
//
//
//
//#define DBG_VISUALIZE_CELLS
#if 0

#define DBG_CELL_X 96
#define DBG_CELL_Y 36
#define DBG_CELL_Z 96
#define DBG_CELL_LEVEL 0
#define DBG_CELL_LEVEL_MG 0

#define IS_DBG_CELL( levelMg, level, x, y, z ) \
  ((levelMg)==DBG_CELL_LEVEL_MG && (level)==DBG_CELL_LEVEL &&  \
  (x)==DBG_CELL_X && (y)==DBG_CELL_Y && (z)==DBG_CELL_Z) \

#else

#define IS_DBG_CELL( levelMg, level, x, y, z ) UT_ASSERT0(FALSE)

#endif

#if 0
#define DBG_BREAK_CELL(L,l,x,y,z) \
if( (L)==DBG_CELL_LEVEL_MG && (l)==DBG_CELL_LEVEL && \
    (x)==DBG_CELL_X && (y)==DBG_CELL_Y && (z)==DBG_CELL_Z)\
{ \
  TRCERR(("DBG\n"));\
}
#endif


#ifdef SOLVE_PRESSURE_DOUBLE
typedef double PSFloat; 
typedef Vec4d Vec4psf;
static inline Vec4psf vecf2vecpsf( const Vec4f& x ) { return _mm256_cvtps_pd( x ); }
static inline Vec4d vecpsf2vec4d( const Vec4psf& x ) { return x; }
#else
typedef float PSFloat;
typedef Vec2f Vec2psf;
typedef Vec4f Vec4psf;
#define vecf2vecpsf(x) x
static inline Vec4d vecpsf2vec4d( const Vec4psf& x ) { return _mm256_cvtps_pd( x ); }
#endif


#define G2P_FLIP_ALPHA_USE_GRID_LEVEL
#define RTA_REUSE_BLOCKMAP
#define RENDERDENS_SORT_BLOCKS_MORTON

//#define RSURF_8_BIT
#define RSURF_16_BIT

#define RQUANT_SINGLE_PREC

#if 1
#define PARTICLE_QUANT_VEL
#define PARTICLE_QUANT_POS
#endif

#define PARTICLE_IDX_64BIT
#define MG_FACE_COEFFS_GALERKIN
#define FACE_COEFF_AND_DIAG_CONST_BLOCKS
#define REGULARIZE_INTERIOR_FACE_DENSITIES
#define PARTICLES_MIN_SPLATWEIGHTSUM (1.0E-8f)  // Minimum weight sum considered as valid sample
#define PARTICLES_CONSTRAIN_SMALL_SPLATWEIGHTS 1e-6f
#define PARTICLES_CONSTRAIN_SMALL_SPLATVALS 1e-7f
#define PHYS_AIR_DENSITY	1.25f
#define PHYS_WATER_DENSITY	1000.f
#define PHYS_AIR_VISCOSITY	1.8e-5f  	/* dynamic viscosity */
#define PHYS_WATER_SURFACE_TENSION 0.073f	
#define SPLAT_AIR_VELOCITY_WITH_PARTICLE_MASS
#ifdef UT_ASSERT_LEVEL_3
#define CHECK_NONFLUID_ZERO_PRE
#define CHECK_NONFLUID_ZERO_POST
#endif
#define FACE_DENSITY_16_BIT 
#define USE_EXT_BLOCK_REFS
#define VCL_COMPR_TILES_WITH_PARALLEL_MMX
#define CONST_FLAGS_UNIFORM_BLOCKS
#define RELAX_BLOCK_GAUSS_SEIDEL
#define RELAX_BLOCKS_RED_BLACK
#define MG_CONST_PROLONGATION_ALPHA 1.0f
#define MG_CONST_PROLONGATION
#define MG_CONST_RESTRICTION
#define DOWNSAMPLE_FACEDENS_FACE_ORIENTED
#define AIR_PROC_AMPL_HF_WITH_DISTLIQ
#define TRACK_DROPLETS_WITH_SUBSTEPS
#define SWAP_PARTICLES_TO_FILE_ASYNC
#define MSBG_DISCLAIM_MEMORY 0
#define CONSTRAIN_VELO_WITH_OBST_PROJECTION
#define MSBG_USE_COARSE_SDF
#define ADV_CONSTRAIN_VELO_NEAR_OBST
#define MSBG_WARN_PRESSURE_LIMIT 1e6
#define INITIAL_STEP_DT_FROM_VMAX
#define FAST_DISCONTINUOUS_FINECOARSE_INTERPOLATION
#ifdef TEST_VORT_VIS
#define CONSTRAIN_VELO_COARSE_FINE
#endif
#define LIST_BASED_RELAXATION
#define MSBG_BOUNDZONE_INCL_RESTRANS
#define USE_LAPL_OPTIMIZED
#define REDISTANCE_INTERIOR_OBSTACLES_SDF
#define FILL_TINY_DIRICHLET_HOLES
#ifdef VOLUME_WEIGHTED_MG_TRANSFER_OPERATORS
#define USE_CELL_VOLUME_FRACTIONS
#endif
#define GHOST_PRESSURE_DENSE_MG_LEVELS
#define SPLATTING_KERNEL_EXACTLY_SPHERICAL
#define MSBG_BLOCK_RELAXATION
#define MSBG_PARTICLE_DEFER_SPLIT 0
#define MSBG_PARTICLE_DEFER_MERGE 3
#define MSBG_SIDE_OBSTACLES
#define FINE_GRAIN_THREAD_BARRIER
#define FORCE_INLINE_KERNELS
#define TIMER_MAX_LEVEL 7

//
//
//

#define MSBG_MAX_RESOLUTION 32768 // ~200 increments per cell at 24 bit quant

#define DT_DEFAULT 0.1

#define MSBG_BSXLOG2_0 4  // bsx0=16 

#define MAX_MAX_CFL	20
//#define MAX_MAX_CFL	16

#define MSBG_MAXRESLEVELS 4 

#define MSBG_MAX_TIME_LEVELS 5  // Max. number of RTA time levels

#define PARTICLES_MAX_HEAT 5000.0f
#define PARTICLES_MAX_DROPLET_RADIUS 4.0f
//#define PARTICLES_MIN_DROPLET_RADIUS 0.025f
#define PARTICLES_MIN_DROPLET_RADIUS 0.01f

#if defined( PER_PARTICLE_MASS_DENSITY ) 
#define MASS_DENSITY_GAUSS_APPROX_16
#endif

//

#define MSBG_EXTRAPOL_FOR_ADV_MINDIST_L1 7
#define MSBG_EXTRAPOL_FOR_ADV_MAXDIST_L1( cfl ) \
  (MAX( MSBG_EXTRAPOL_FOR_ADV_MINDIST_L1, \
	ceil( (3./MT_SQRT3)*(cfl) ) )+3) 

//
//
//
#define MSBG_OBST_IPOL_OPT \
 		( SBG::OPT_IPCORNER |  /* obstacle SDF is sampled at grid nodes */\
		  SBG::OPT_IPCHKBLK | \
		  SBG::OPT_IPEXTRAPOL | \
		  SBG::OPT_IPBC_NEUMAN )

#define MSBG_CSDF_IPOL_OPT \
	(OPT_IPCHKBLK | OPT_IPEXTRAPOL | OPT_IPBC_NEUMAN)

#define MSBG_MIN_OBST_GRAD_MAG 1e-4

//
//

#define PARTICLES_WITH_DENSITY

#ifdef __cplusplus
namespace MSBG
{

// Grid cells for FLIP particles

typedef struct 
{
  Vec4f quantSum_;	// [vx,vy,vz,dens]
  Vec4f v4fWeightSum_;
}
ParticleCell;

const float SPEC_FLOAT_NEUMANN = -1e20f,
      	    SPEC_FLOAT_DIRICHLET = 1e20f;

} // namespace MSBG

#endif  // __cplusplus
#endif // MSBGCELL_H
