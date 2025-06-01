/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#ifndef PNS_H
#define PNS_H

#define PNS_DOUBLE_PREC

#ifdef __cplusplus
namespace CurlNoise2
{

  Vec4f curlNoise4D( Vec4f pos );

  template<int DO_RET_POT_VAL=false, int S_CURVE_TYP=3 /*1=linear,3=quintic*/, 
    	   int DO_DEBUG = false>
    Vec4f curlNoise4D_( Vec4f pos, Vec4f *pPotValue=NULL );

  enum
  {
    OPT_VALUE = 1<<0,
    OPT_DX = 1<<1,  
    OPT_DY = 1<<2,
    OPT_DZ = 1<<3
  };

Vec4f curlNoiseSum4D( 
		const Vec4f& pos0,  	 // position (4D)
		float freq_sp,   // frequency (spatial)
		float freq_tm,	 // frequency (time)
		float alpha,
		float beta_sp,
		float beta_tm,
		float seed,
		float octaves,     // IN: ocatves 
		int oct_leadin,
		int lf_min_oct,	// optional cutoff for LF component (out_lf)
		const Vec4f& ampl,
		float ampl_hf_gain, // high frequency boost
		float ampl_hf_min_oct,
		const Vec4f *modfield,  // Modulating scalar field and derivatives f,fx,fy,fz
		unsigned options,
		Vec4f *out_lf  // Optional out: LF filtered (
    	);
}
#endif
#ifdef __cplusplus
namespace PNS
{

float genFractalScalar3D(
    	int sx, int sy, int sz,  // grid resolution
    	int x, int y, int z, // grid coords
	double seed,
    	double fscale,  // Feature scale in voxels
	double minscale,
	double octaves_leadin,
	double alpha,   // Fractal amplitude persistence
	int do_turb	// Apply fractal 'turbulence' (abs(noise))
	);

void genFractalScalarField(  
        int dim,
    	int sx, int sy, int sz,
	double seed,
    	double fscale,  // Feature scale in voxels
	double minscale,
	double octaves_leadin,
	double alpha,   // Fractal amplitude persistence
	int do_turb,	// Apply fractal 'turbulence' (abs(noise))
	float *field    // OUT: generated field
	);

namespace VVN
{
void fractal_vec3( 
    	      int dim,
	      double *pos, 
	      double octaves, 
	      int octaves_leadin,  // see http://www.planetside.co.uk/forums/index.php/topic,17566.msg170871.html#msg170871
              double lacunarity,
	      double alpha,  // roughness
	      int turb_typ,  // 0=fbm, 1=abs, 2=ridged
	      //float ridge_offs, float ridge_mult, float ridge_exp,

	      float *outvec	// 3d vector output
    	   );

  void noise3_vec3_( double pos[3],
		     float *vec3_out );
}

namespace simd
{
Vec4f noise3d_simd4_i( const Vec4d &px, const Vec4d &py, const Vec4d &pz  );
Vec4f fractal_simd4_i( 
    	      int dim, 
	      const Vec4d& px_, const Vec4d& py_, const Vec4d& pz_, 
	      double octaves, 
	      int octaves_leadin,  // see http://www.planetside.co.uk/forums/index.php/topic,17566.msg170871.html#msg170871
              double lacunarity,
	      double alpha,  // roughness
	      int turb_typ  // 0=fbm, 1=abs, 2=ridged
	      //float ridge_offs, float ridge_mult, float ridge_exp,
    	   );
}
}
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define PNS_MAX_OCTAVES	32

void PNS_noise3_vec3_test(double pos[3], float *out);

int PNS_Init2(int rseed);
int PNS_Init3(int rseed);

float PNS_simdtst_noise3_v4( double *p );
float PNS_donw_noise3d( double *p );

int pns_rand_init(int seed, void *rnds);
int pns_rand(void *rnds);

void PNS_noise3d_simd4_f( const float *px, const float *py, const float *pz, int n,
    			float *result );

void PNS_noise3d_simd8( const double *px, const double *py, const double *pz, int n,
    			float *result );

double PNS_PerlinNoise4D( double *pos, double alpha, double beta, double octaves );

typedef struct
{
  double octaves;
  int octaves_leadin;
  float noisevals[PNS_MAX_OCTAVES];
  int n_noisevals;
}
PNS_Spectrum;

typedef struct
{
  int typ;
  double ampl,		// Amplitude
	 offs,		// Value offset
    	 fscale,	// Feature scale 
         alpha,		// Amplitude decay factor
  	 beta,		// Frequency decay factor
	 coff;		// Coordinate offset
  int oct_leadin;
  double oct_leadin_pow;
}
PNS_Params;

int PNS_Init(int rseed);

float PNS_Noise4D( double* pos, int improveSmoothness );

double PNS_Noise3D( double* pos );

double PNS_PerlinNoise2D(double x,double y,double alpha,double beta,double octaves,
    			 int typ);

double PNS_PerlinNoise3D(double x,double y,double z,double alpha, double beta,
    		     double octaves,
		     int do_pyroclastic,
		     int typ);

double PNS_Fractal_2( 
    	      int dim,
	      double *pos, 
	      int do_comp_spectrum,	// opmode
	      int do_comp_valsum,	      
	      // Optional IN/OUT of <n_octaves> noise values for sepctral   
	      // texture blending as of [Neyret2003]
	      PNS_Spectrum *spectrum,   
	      double octaves, 
	      int octaves_leadin,  // see http://www.planetside.co.uk/forums/index.php/topic,17566.msg170871.html#msg170871
              double lacunarity,
	      // 
	      double alpha,
	      int turb_typ,  // 0=fbm, 1=abs, 2=ridged
	      float ridge_offs, float ridge_mult, float ridge_exp
    		);

double PNS_Fractal3D(double x,double y,double z,double alpha, double beta,
    		     double octaves, 
		     int octaves_leadin,  // see http://www.planetside.co.uk/forums/index.php/topic,17566.msg170871.html#msg170871
		     int turb_typ,  // 0=fbm, 1=abs, 2=ridged
		     float ridge_offs, float ridge_mult, float ridge_exp, 
		     int do_multifract,
		     int ntyp);

double PNS_Fractal(double *pos, int dim,  // position, dimension
    		     double alpha, double beta,
    		     double octaves, 
		     int octaves_leadin,  // see http://www.planetside.co.uk/forums/index.php/topic,17566.msg170871.html#msg170871
		     int turb_typ,  // 0=fbm, 1=abs, 2=ridged
		     float ridge_offs, float ridge_mult, float ridge_exp, 
		     int do_multifract,
		     int ntyp,
		     double seed, double seed_tm  // seed (for spatial & time dimension)
    			);

#if 1
double PNS_Fractal_warp(double *pos, int dim,  // position, dimension
    			double *W,
    		     double alpha, double beta,
    		     double octaves, 
		     int octaves_leadin,  // see http://www.planetside.co.uk/forums/index.php/topic,17566.msg170871.html#msg170871
		     int turb_typ,  // 0=fbm, 1=abs, 2=ridged
		     float ridge_offs, float ridge_mult, float ridge_exp, 
		     int do_multifract,
		     //double pa_mem, double pa_ampl, // 'structure memory' & amplitude for flow noise 'pseudo advection' see Ebert. p 384 ff
		     int ntyp,
		     double seed, double seed_tm  // seed (for spatial & time dimension)
    			);
#endif

double PNS_Fractal_SIMD_f(double *pos, int dim,  // position, dimension
    		     double alpha, double beta,
    		     double octaves, 
		     int octaves_leadin,  // see http://www.planetside.co.uk/forums/index.php/topic,17566.msg170871.html#msg170871
		     int turb_typ,  // 0=fbm, 1=abs, 2=ridged
		     float ridge_offs, float ridge_mult, float ridge_exp, 
		     int do_multifract,
		     int ntyp);


#define PNS_NTYP_NONE	0
#define PNS_NTYP_FAST	1	// ATTENTION: DO NOT CHANGE ! 
#define PNS_NTYP_IMPROVED 2	// ATTENTION: DO NOT CHANGE !

#ifdef __cplusplus
}  // extern "C"
#endif

#endif /* PNS_H */

