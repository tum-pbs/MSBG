/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#ifndef GRID_H
#define GRID_H

#ifdef __cplusplus
#include <vectorclass/vectorclass.h>
#endif

#define GRID_GEN		(1<<0)
#define GRID_IPCORNER		(1<<1)
#define GRID_FLOAT		(1<<2)
#define GRID_GAUSS		(1<<4)
#define GRID_LAPLACE		(1<<5)
#define GRID_VERBOSE		(1<<6)
#define GRID_FREE		(1<<7)
#define GRID_FASTMIPMAP		(1<<8)
#define GRID_FTRIANGLE		(1<<9)
#define GRID_NO_BND_CHK 	(1<<10)
#define GRID_CELLCENTERED	(1<<11)
#define GRID_TEST		(1<<12)

#define GRID_IPNEAREST		1
#define GRID_IPLINEAR		2
#define GRID_IPCUBIC		3
#define GRID_IPCUBIC_MONO	4
#define GRID_IPCUBIC_MONO_2	5
#define GRID_IPQUADRATIC	6
#define GRID_IPFRACTAL_MONO	7
#define GRID_IPWENO4		8
#define GRID_IPQUINTIC		9


#define GRID_BDEFAULT	0
#define GRID_BCLIP	1


#define G3D_MAX_PYRA_LEVEL 24
#define G3D_MAX_FCHAN	32

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  LongInt x0,y0,z0,
          sx,sy,sz;
}
G3D_Box;

typedef struct 
{
  struct G3D_GridAdm	*grid0;
  int			 n_level,
			 laplace_valid;
  unsigned		 opt;
  float	                *krn;
  int		         krn_r;
  double		 res_0,
			 res_factor,
			 res_act;
  struct G3D_GridAdm    *gaussian[G3D_MAX_PYRA_LEVEL],
			*laplacian[G3D_MAX_PYRA_LEVEL];
}
G3D_Pyramid;

struct G3D_GridAdm
{
  LongInt sx,sy,sz,
	  sxy, n,
	  smax,
	  smin,
	  dim[4],
	  tsx,tsy,tsz; // Tiling 
  double  ws_pos[4],    // world space position of grid origin
	  ws_dx,ws_dx_r;  // world space voxel size and reciprocal
  float   soff[3];       // Grid coordinate offset (staggering)
  float	*dataFloat[G3D_MAX_FCHAN];
  unsigned char dataFloatIsExt[G3D_MAX_FCHAN];
  uint16_t *dataUshort[G3D_MAX_FCHAN];
  G3D_Pyramid pyramid;
};

#define G3D_SetGridOffset( G, ox, oy, oz )\
{\
  (G)->soff[0] = ox;\
  (G)->soff[1] = oy;\
  (G)->soff[2] = oz;\
}

#define G3D_SetGridBorder( G, db0, db1 )\
{\
  (G)->db0 = db0;\
  (G)->db1 = db1;\
}

#define G3D_GRID_OFFS(G) ((G)->soff)

typedef struct G3D_GridAdm *G3D_Grid;

G3D_Grid G3D_CreateGrid0( int sx, int sy, int sz,
    			  unsigned opt,
			  char *mmTag);
#define G3D_CreateGrid( sx, sy, sz, opt )\
	G3D_CreateGrid0( sx, sy, sz, opt,\
	    		 MM_MK_TMP_UID())

int G3D_DeleteGrid( G3D_Grid *p_adm );
void G3D_SetBorders(G3D_Grid G, float *data, int db, float val);
void G3D_CopyOutBorders(G3D_Grid G, float *field);

void G3D_SetWorldSpaceInfo( G3D_Grid G, double *pos, double dx );

float *G3D_GetFloatChannel0( G3D_Grid grid, int i, unsigned opt,
    			     char *mmTag );
#define G3D_GetFloatChannel( grid, i, opt ) \
  G3D_GetFloatChannel0( grid, i, opt, MM_MK_TMP_UID())

uint16_t *G3D_GetUshortChannel( G3D_Grid grid, int i, unsigned opt );
void G3D_FreeUshortChannel( G3D_Grid grid, int i );
int G3D_SetFloatChannel( G3D_Grid grid, int i, float *data );
void G3D_FreeFloatChannel( G3D_Grid grid, int i );

int G3D_CreatePyramid( G3D_Grid G,  	    // target bitmap
		      int lmin,int lmax,    // min,max level (-1=auto)
		      int chan,	int ichan,  // data channel
		      float *kern_down,  // optional: 5x5x5 kernel
		      float *kern_up,    // optional
		      int btyp,
    		      unsigned opt );

int G3D_FreePyramid( G3D_Grid grid );

int G3D_UpSamplePyramid( G3D_Grid *G, int chan, int ichan,
    		        int l_from, int l_to );

int G3D_UpSamplePyraLevel( 
    	G3D_Grid bmpLo, float *dataLo, uint16_t *dataLoUs,
    	G3D_Grid bmpHi, float *dataHi, uint16_t *dataHiUs,
	float *dataLapl, 
	int oxh, int oyh, int ozh,
        float *kernel,  // optional 5x5 upsampling kernel
	int btyp // boundary handling 
	);

int G3D_GetBurtAdelsonKernel(
    		double a, // IN: smoothness/shape parameter  
			 //   0.6=neg.lobes, 0.5=triangle, 
			 //   0.4=gauss, 0.3=broader gauss
			 //   0 = box
			 //   see also [Burt1983]
    		int krn_radius,  // kernel radius: 2 for gauss, 1 for triangle,...
    		double *krn_data  // OUT: kernel 
		);

void G3D_GetGausslikeKernel4x4x4(
    		// even-sized kernel for cell centered filtering 
		// See also: "magic kernel" http://johncostella.webs.com/magic/
    		float *krn_data  // OUT: kernel (normalized wsum=1)
		);

void G3D_DownsampleBlock5x5x5( 
			 G3D_Grid grid,
			 int x1, int x2,	// target range
			 int y1, int y2,
			 int z1, int z2,
			 float *data,

			 G3D_Grid gridHi,
			 float *dataHi,

			 double *kern5x5x5
			 );
#ifdef __cplusplus

void G3D_ReDistanceFastSweeping( G3D_Grid G, float *phi,
    				float phiInfinity, float h,
				unsigned options=0 );

void G3D_InterpolateVec( 
		  G3D_Grid G,		// Input grid
		  Vec4f pos,		// Input position x,y,z
		  int dim,		// data dimension
		  const float **data,		// data pointers
		  int method,		// Interpolation method
		  unsigned opt,		// Options
		  float *p_out		// OUT: interpolated output value
		);

static inline 
Vec8f InterpolLinearVec8( 
		Vec8f *F, 
		Vec8f const & U ) 
{
  return  (U*F[1]+(1.0f-U)*F[0]);
}

static inline 
Vec8f InterpolCubicVec8( 
    		Vec8f *F, 
    		Vec8f const & U, 
		Vec8f const & U2, 
		Vec8f const & U3 )
{
  Vec8f D1 = (F[2]-F[0])*0.5f,
	D2 = (F[3]-F[1])*0.5f,
	D0 = F[2]-F[1];

  Vec8f A0 = F[1],
	A1 = D1,
	A2 = 3.f*D0 - 2.f*D1 - D2,
	A3 = D1 + D2 - 2.0f*D0;

  return A3*U3 + A2*U2 + A1*U + A0;
}

static inline 
Vec8f InterpolCubicMono2Vec8( 
    		Vec8f *F, 
    		Vec8f const & U, 
		Vec8f const & U2, 
		Vec8f const & U3 )
{
  Vec8f D0 = F[1]-F[0],
	D1 = F[2]-F[1],
	D2 = F[3]-F[2];

  Vec8f M0 = select( D1*D0 > MT_NUM_EPS_F, 
      		     6.0f/(3./(D0) + 3./(D1)),
      		     0.0f );

  Vec8f M1 = select( D2*D1 > MT_NUM_EPS_F, 
      		     6.0f/(3./(D1) + 3./(D2)),
      		     0.0f );

  Vec8f
    S = D1,
    A = 3.0f*S - 2.0f*M0 - M1,
    B = M0 + M1 - 2.0f*S,

    G = F[1] + M0*U + A*U2 + B*U3;

  return G;
}

static inline
void wenoParabola( Vec8f f1, Vec8f f2, Vec8f f3, Vec8f x, // IN
    		     Vec8f *w, Vec8f *g // OUT 
		     )
{
  Vec8f d = (f3-f1)*0.5f,	// 1st derivative at 0
	    dd = (f1-f2*2.0f+f3),	// 2nd derivative
	    S = d*(d+dd) + dd*dd*(4.f/3.f);  // Smoothness over x=[0,1]
   *w = Vec8f(2.0f-x)/((S+1e-6f)*(S+1e-6f));	// Relative weight
   *g = f2 + (d+0.5f*dd*x)*x;	// Value of parabola at x
}

static inline 
Vec8f InterpolWENO4Vec8( 
    		Vec8f *F, 
    		Vec8f const & U ) 
{
  Vec8f w1,p1;
  wenoParabola( F[0], F[1], F[2], U,
		&w1, &p1 );
  Vec8f w2,p2;
  wenoParabola( F[3], F[2], F[1], 1-U,
		&w2, &p2 );
  return (w1*p1 + w2*p2) / (w1+w2);  
}

#endif

int G3D_Interpolate2( 
		  G3D_Grid G,		// Input grid
		  const double *goff,	// Grid offset of interpolated quantity samples (.5,.5,.5=center)
		  double *p0,		// Input position x,y,z
		  const float *data,		// data field to be interpolated
		  uint8_t *obstacles,   // optional obstacles
		  int method,		// Interpolation method
		  unsigned opt,		// Options
		  double *p_out,	// OUT: interpolated output value
		  double *p_neigh_min, // MIN/MAX of interpolation neighborhood
		  double *p_neigh_max,
		  double *p_neigh_bound	// Interpolation neighborhood touches obstacle or boundary
		);

void G3D_InterpolateLinear( 
		  G3D_Grid G,		// Input grid
		  double *p0,		// Input position x,y,z
		  float *data,
		  uint16_t *data_us,
		  double *p_out		// OUT: interpolated output value
		);

int G3D_Interpolate( 
		  G3D_Grid G,		// Input grid
		  double *p0,		// Input position x,y,z
		  float *data,
		  uint16_t *data_ushort,
		  int method,		// Interpolation method
		  unsigned opt,		// Options
		  double *p_out		// OUT: interpolated output value
		);

int G3D_AdvectFieldEulerian(
    		const double dt, 
		// velocity field driving the advection
		const float* velx, const float* vely,  const float* velz,
		float velMaxMag,
		// advected field
		float* oldField, float* newField, 
		int xres, int yres, int zres,
		float dtFactor  // 'safety'factor for CFL based time step
		);

int G3D_AdvectAllFieldsSemiLagrangian(
    		const double dt, 
		// Grid resolution (co-located velocities at cell centers !)
                int xres, int yres, int zres, 
		// Advector velocity field
		const float* velx, const float* vely,  const float* velz,
		// advected fields
		int numFields,
		const float **oldFields, float **newFields, 
		// misc. parameters
		double advectMaxVel,	// max. allowed step (1=grid-cell size)
		int traceType,	// 1=euler, 2=midpoint-rule
		int ipType, 	// GRID_IP*
		int ipTypeVel
		);

int G3D_AdvectFieldSemiLagrangian2(
    		const double dt, 
		// velocity field driving the advection
		const float* velx, const float* vely,  const float* velz,
		// advected field
		float* oldField, float* newField, 
		uint8_t *obstacles,
		int xres, int yres, int zres, 
		// misc. parameters
		double advectMaxVel,	// max. allowed step (1=grid-cell size)
		int traceType,	// 1=euler, 2=midpoint-rule
		int ipType, 	// 1=linear, 2=cubic
		int ipTypeVel,
		int doMacCormack,
		int doClampZero
		);

int G3D_AdvectFieldSemiLagrangian(
    		const double dt, 
		// velocity field driving the advection
		const float* velx, const float* vely,  const float* velz,
		// advected field
		float* oldField, float* newField, 
		int xres, int yres, int zres, 
		// misc. parameters
		double advectMaxVel,	// max. allowed step (1=grid-cell size)
		int traceType,	// 1=euler, 2=midpoint-rule
		int ipType, 	// 1=linear, 2=cubic
		int ipTypeVel,
		int doMacCormack,
		int doClampZero
		);

void G3D_AdvectFieldSemiLagrangeCubic(
		 G3D_Grid G,
		 double dt, 
		 const float* velx,   // velocity field (sampled at cell centers)
		 const float* vely,  
		 const float* velz,
		 const float* oldField,   // density field (sampled at cell ceneters)
		 float* newField,  // OUT
		 unsigned options );

int G3D_TestAdvectFieldSemiLagrangeCubic(void);
int G3D_TestAdvectFieldSemiLagrange(void);

#define G3D_HAS_PYRAMID(G) ((G)->pyramid.n_level>0)

#define G3D_GET_ROI( G, roi, x0,y0,z0, sx,sy,sz ) \
{\
  if(roi)\
  {\
    x0=roi->x0; y0=roi->y0; z0=roi->z0;\
    sx=roi->sx; sy=roi->sy; sz=roi->sz;\
  }\
  else\
  {\
    x0=0;y0=0;z0=0;\
    sx=G->sx;sy=G->sy;sz=G->sz;\
  }\
  UT_ASSERT_FATAL(G3D_IN_RANGE(G,x0,y0,z0)&&\
      G3D_IN_RANGE(G,x0+sx-1,y0+sy-1,z0+sz-1));\
}

#define G3D_DIM_EQUAL( b1_, b2_ ) \
  ( ((b1_)->sx == (b2_)->sx) && ((b1_)->sy == (b2_)->sy) && ((b1_)->sz == (b2_)->sz) )

#define G3D_IN_RANGE(grid, x,y,z) \
  G3D_IN_BOX(0,0,0, (grid)->sx-1,(grid)->sy-1,(grid)->sz-1, x, y, z)

#define G3D_IN_BOX(xmin,ymin,zmin, xmax,ymax,zmax, x,y,z) \
 ((x)>=(xmin)&&(x)<=(xmax)&&(y)>=(ymin)&&(y)<=(ymax)&&z>=(zmin)&&(z)<=(zmax))

#ifdef UT_ASSERT_LEVEL_2
#define G3D_INDEX(G,x,y,z) \
  	( G3D_IN_RANGE(G,x,y,z) ? \
		((x)+(y)*((LongInt)G->sx)+(z)*(LongInt)(G->sxy)) :\
		UT_THROW_ASSERTION() )
#else
#define G3D_INDEX(G,x,y,z) ((x)+(y)*((LongInt)G->sx)+(z)*(LongInt)(G->sxy))
#endif

#define G3D_INDEX0( grid_sx, grid_sxy, x, y, z ) \
  ((x)+(y)*((LongInt)(grid_sx))+(z)*(LongInt)(grid_sxy))

#define G3D_IDX2COORDS( G, idx, /*IN*/\
    		    x, y, z /*OUT*/)\
{ \
  z = (LongInt)(idx) / (LongInt)(G->sxy); \
  LongInt rem_ = (LongInt)(idx) % (LongInt)(G->sxy); \
  y = rem_ / (LongInt)(G->sx); \
  x = rem_ % (LongInt)(G->sx); \
}

#define G3D_IDX2COORDS0( sx, sxy, idx, /*IN*/\
    		    x, y, z /*OUT*/)\
{ \
  z = (idx) / (sxy); \
  x = (idx) % (sxy); \
  y = x / (sx); \
  x = x % (sx); \
}



#define G3D_CLIP_CONST(x1,y1,z1, x2,y2,z2, /* IN: region */\
    		      x,y,z  /* OUT: coordantes to clip */\
    		     )\
{\
  if(unlikely((x)>(x2))) x = x2; \
  if(unlikely((x)<(x1))) x = x1; \
  if(unlikely((y)>(y2))) y = y2; \
  if(unlikely((y)<(y1))) y = y1; \
  if(unlikely((z)>(z2))) z = z2; \
  if(unlikely((z)<(z1))) z = z1; \
}

#ifdef __cplusplus
} // extern C
#endif

#ifdef __cplusplus
template<typename T>
double G3D_GetGaussianKernel(
    		double sigma,
    		int width,
    		T *kernel,  // OUT: kernel [width*width*width]
		int verbose=0
		);
#endif

#endif /* GRID_H */

