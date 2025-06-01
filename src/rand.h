/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#ifndef RAND_H
#define RAND_H
/*=========================================================================
 *
 * 	P R N G	(Pseudo Random Number Generator)
 *
 * =======================================================================*/
#ifdef __cplusplus
extern "C" {
#endif

#define RND_CHECK		(1<<1)
#define RND_VERBOSE		(1<<2)

#define RND_MAX_NAME_LEN 	64
#define RND_MAX_STATE_SZ	1024 /* TODO */
#define RND_EMAGIC		1346146
typedef unsigned char RND_State[RND_MAX_STATE_SZ];
typedef int (*RND_SeedFunc)( RND_State *rnds, int seed );
typedef double (*RND_GetFunc)( RND_State *rnds );

#define RND_Get( rndStream ) \
  (rndStream)->gen->rndGet( &(rndStream)->state )

#define RND_IsInit( rndStream ) \
  ((rndStream)->emagic == RND_EMAGIC)

typedef struct RND_Generator_s
{
  char name[RND_MAX_NAME_LEN+1];
  RND_SeedFunc rndSeed;
  RND_GetFunc rndGet;
}
RND_Generator;

typedef struct
{
  RND_Generator *gen;
  RND_State      state;
  int		 emagic;
}
RND_Stream;

typedef struct
{
  int	  m;
  double *L,
	 *vd,
	 *vx;  
}
RND_MultNormState;

typedef struct
{
  int	  n,	// number of varaibles (dimensions)
	  m;	// number of constraints 
  double *L,	// mxn constrains matrix
	 *LCT,
	 *Y,*Y2,*YL,*YE;  
  RND_MultNormState *rmnsE;
}
RND_MultNormLECState;

RND_Generator *RND_GetGenerator( char *name );

int 	RND_Seed_numc( RND_State *rnds, int seed );
double 	RND_Get_numc( RND_State *rnds );

int 	RND_Seed_knuth( RND_State *rndsIn, int seed );
double 	RND_Get_knuth( RND_State *rnds );


void RND_ShuffleIntArrayOld(int *v, int n, 
    		         RND_Generator *prng, RND_State *rnds );

void RND_ShuffleIntArray(int *v, int n, 
    			 RND_Stream *rnds );

void RND_ShuffleIntArray2(int *v, int n_total, 
    			  int n_shuffle,
    			 RND_Stream *rnds );

int RND_Seed(  RND_Stream *rs, int seed );
int RND_Init(  RND_Stream *rs, RND_Generator *gen, int seed );
int RND_Test(void);
// Get random integer between 0 and n-1					   
int RND_GetInt(  RND_Stream *rs, int n );
double RND_GetNorm( RND_Stream *rs, double mean, double sigma );
int RND_GetMultNorm( 
    		RND_Stream *rs, 
    		int	m,	
			// dimension
    		double *mCorrCov, 
			// correlation or covariance matrix
			// if mCorrCov is correlation matrix
			// then multiply result by std. dev. vector
			// before adding mean vector
		double *vOut,    
			// output vector 
			// -> add mean and multiply variances 
		RND_MultNormState **pState,
			// IN/OUT (optional): Cholesky Decomposition
			// for repeated sampling with same correlation
			// structure- ( *pChol==NULL -> init and return
			// *pChol != NULL -> re-use
			// Memory must be freed by caller
		unsigned opt
			// options
			);
int RND_DelMultNormState( RND_MultNormState *s );
int RND_Cholesky( 
 		double *m,		
  		int n,		
  		double *ld);
double RND_GetTruncNorm( RND_Stream *rs, double m, double sd,
    			 double left, double right,
			 int	*p_num_iter );

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int RND_GetMultNormLEC(
    	// Linear Equality Constrained Sampling (LEC)  
    	// Transform unconstrained sample X0 from N(mu,S) to
	// a sample X subject to Linear Equality Constraint
	// A * X = B  (See [Rue2002], GMRFlib documentation)
			RND_Stream *rs,
				//
			int	n,
				// number of dimensions
    			double *S,
				// IN: Covariance Matrix (nxn)
			int	m,
				// Number of constraints
			double *A,
				// Matrix (mxn) of m linear constraints 
				// m = number of constraints
				// Required: Rank(A) = m 
			double *B,
				// Right-Hand-Side Vector (mx1)
				// for constrains A*X = B
		        double *E,
				// Optional mxm Covariance Matrix for 
				// stochastic (smooth) conditioning
			double *X,
				// IN/OUT: Original unconstrained 
				// sample (nx1) from N(0,S)
				// On output: Sample from constrained 
				// distribution
			RND_MultNormLECState **pState
				// IN/OUT: (optional)
				// State memory for efficient 
				// repeated calculation
				// Must be deallocated by caller
			);

int RND_DelMultNormLECState( RND_MultNormLECState *s );

#ifdef __cplusplus
}
#endif

#endif /* RAND_H */

