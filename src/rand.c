/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdarg.h>
#include <ctype.h>
#include <limits.h>
#include "stdmath.h"
#include <float.h>
#include "globdef.h"
#include "mtool.h"
#include "rand.h"

/*=========================================================================
 *
 *
 * 	R N D	(Pseudo Random Number Generator)
 *
 *
 * =======================================================================*/

RND_Generator 
  RND_Generator_numc = { "numc", RND_Seed_numc, RND_Get_numc },
  RND_Generator_knuth = { "knuth", RND_Seed_knuth, RND_Get_knuth };

RND_Generator *RND_Generators[] = 
{
  &RND_Generator_numc,
  &RND_Generator_knuth,
/*  { "libc", RND_Seed_libc, RND_Get_libc },
  { "knuth", RND_Seed_knuth, RND_Get_knuth };*/
};

RND_Generator *defaultGen = &RND_Generator_numc;

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
RND_Generator *RND_GetGenerator( char *name )
{
  int i,found=FALSE;
  RND_Generator *rndg=NULL;

  for(i=0;i<ARRAY_LENGTH(RND_Generators);i++)
  {
    rndg = RND_Generators[i]; 
    if(!stricmp(rndg->name,name)) 
    {
      found = TRUE;
      break;
    }
  }
  return found ? rndg : NULL;
}

/*=========================================================================
 *
 * 	numc 	(standard from numerical recieps)
 *
 * =======================================================================*/
/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int RND_Seed_numc( RND_State *rnds, int seed )
{
  MtRndSSeed((MtRndStream *)rnds, seed );
  return 0;
}
/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
double RND_Get_numc( RND_State *rnds )
{
  return (double)MtRndSGet((MtRndStream *)rnds);
}
/*=========================================================================
 *
 * 	knuth	(cf. Knuth 1981)
 *
 * =======================================================================*/

typedef struct
{
  unsigned int m,a,c,x;
}
RND_State_knuth;
/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int RND_Seed_knuth( RND_State *rndsIn, int seed )
{
  RND_State_knuth *rnd=(RND_State_knuth *)rndsIn;

  rnd->m = (1<<31);	/* cf. Knuth 1981 */
  rnd->a = 7*7*7*7*7;
  rnd->c = 0;
  rnd->x = 19+seed;
  return 0;
}
/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
double RND_Get_knuth( RND_State *rnds )
{
  RND_State_knuth *rnd=(RND_State_knuth *)rnds;

  rnd->x = ( rnd->a * rnd->x + rnd->c ) % rnd->m;
  return (double)(rnd->x)/(double)(rnd->m);
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void RND_ShuffleIntArrayOld(int *v, int n, 
    			 RND_Generator *prng, RND_State *rnds )
{
  int i,j,tmp;

    if( n > 1 )
    {
      for(i=0;i<n;i++)
      {
	j = n * prng->rndGet(rnds);
	if( j <= 0 ) j = 0;
	if( j > n - 1 ) j = n - 1;
	tmp = v[i];
	v[i] = v[j];
	v[j] = tmp;  
      }  
    }

}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void RND_ShuffleIntArray(int *v, int n, 
    			 RND_Stream *rnds )
{
  int i,j,tmp;

    if( n > 1 )
    {
      for(i=0;i<n;i++)
      {
	j = n * RND_Get(rnds);
	if( j <= 0 ) j = 0;
	if( j > n - 1 ) j = n - 1;
	tmp = v[i];
	v[i] = v[j];
	v[j] = tmp;  
      }  
    }

}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void RND_ShuffleIntArray2(int *v, int n_total, 
    			  int n_shuffle,
    			 RND_Stream *rnds )
{
  int i,j,tmp;

    if( n_shuffle > 1 )
    {
      for(i=0;i<n_shuffle;i++)
      {
	j = n_total * RND_Get(rnds);
	if( j <= 0 ) j = 0;
	if( j > n_total - 1 ) j = n_total - 1;
	tmp = v[i];
	v[i] = v[j];
	v[j] = tmp;  
      }  
    }

}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int RND_Init(  RND_Stream *rs, RND_Generator *gen, int seed )
{
  rs->gen = gen ? gen : defaultGen;
  rs->gen->rndSeed( &rs->state, seed );
  rs->emagic = RND_EMAGIC;
  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int RND_Seed(  RND_Stream *rs, int seed )
{
  rs->gen->rndSeed( &rs->state, seed );
  return 0;
}

/*-------------------------------------------------------------------------*/
/* Get random integer between 0 and n-1					   */
/*-------------------------------------------------------------------------*/
int RND_GetInt(  RND_Stream *rs, int n )

{
  int i = ((n)*(rs)->gen->rndGet(&(rs)->state));
  if((i>=n&&(n!=0))||i<0) 
  {
    TRCERR(("Assertion: random number %d out of bound %d\n",
	  i,n));
    exit(1);
  }
  return i;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
double RND_GetNorm( RND_Stream *rs, double mean, double sigma )
{
  static int    iset = 0;
  static double  gset;
  double         fac, rsq, v1, v2;
  
  if( iset == 0 )
  {
    do
    {
      v1 = 2.0 * RND_Get(rs) - 1.0;
      v2 = 2.0 * RND_Get(rs) - 1.0;
      rsq = v1*v1 + v2*v2;
    }
    while(( rsq >= 1.0 )||( rsq == 0.0 ));
    fac = sqrt( -2.0*log( rsq ) / rsq );
    gset = v1 * fac;
    iset = 1;
    return sigma * ( v2 * fac ) + mean;
  }
  else
  {
    iset = 0;
    return sigma * gset + mean;
  }
}

#if 0
int dFactorCholesky (double *A, int n)
{
  int i,j,k,nskip;
  double sum,*a,*b,*aa,*bb,*cc,*recip;
  dAASSERT (n > 0 && A);
  nskip = dPAD (n);
  recip = (double*) ALLOCA (n * sizeof(double));
  aa = A;
  for (i=0; i<n; i++) {
    bb = A;
    cc = A + i*nskip;
    for (j=0; j<i; j++) {
      sum = *cc;
      a = aa;
      b = bb;
      for (k=j; k; k--) sum -= (*(a++))*(*(b++));
      *cc = sum * recip[j];
      bb += nskip;
      cc++;
    }
    sum = *cc;
    a = aa;
    for (k=i; k; k--, a++) sum -= (*a)*(*a);
    if (sum <= REAL(0.0)) return 0;
    *cc = dSqrt(sum);
    recip[i] = dRecip (*cc);
    aa += nskip;
  }
  return 1;
}
#endif

double squared_norm
( double *v,		/* The vector */
  int d,		/* Distance between elements */
  int n			/* Number of elements */
)
{
  double s;
  int i;

  s = 0;

  for (i = n; i>0; i--)
  { s += *v * *v;
    v += d;
  }

  return s;
}

double inner_product
( double *v1,		/* First vector */
  int d1,		/* Distance between elements of first vector */
  double *v2,		/* Second vector */
  int d2,		/* Distance between elements of second vector */
  int n			/* Number of elements in the vectors */
)
{
  double s;
  int i;

  s = 0;

  for (i = n; i>0; i--)
  { s += *v1 * *v2;
    v1 += d1;
    v2 += d2;
  }

  return s;
}

int cholesky
( double *m,		/* The matrix */
  int n,		/* Number of rows and columns of matrix */
  double *ld		/* Place to store log of determinant, or zero */
)
{
  double *r, *p;
  double s;
  int i, j;

  if (ld) *ld = 0;

  r = m;

  for (i = 0; i<n; i++)
  {
    s = r[i] - squared_norm(r+i-1,-1,i);

    if (s<1e-100)
    { return 0;
    }
    else
    { if (ld) *ld += log(s);
      s = sqrt(s);
      r[i] = s;
      p = r;
      for (j = i+1; j<n; j++)
      { p += n;
        p[i] = (p[i] - inner_product(r+i-1,-1,p+i-1,-1,i)) / s;
      }
    }

    r += n;
  }

  return 1;
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int RND_Cholesky( 
 double *m,		/* The matrix */
  int n,		/* Number of rows and columns of matrix */
  double *ld		/* Place to store log of determinant, or zero */
  )
{
  return cholesky(m,n,ld);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int RND_DelMultNormState( RND_MultNormState *s )
{
  if(s)
  {
    FREEMEM(s->L); FREEMEM(s->vd); FREEMEM(s->vx);
    FREEMEM(s);
  }
  return 0;
}

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
			)
{

  RND_MultNormLECState *st=NULL;
  int rcThis=0,i,j;
  double *T=NULL,*T2=NULL;

  if(pState)
  {
    st = *pState;
  }

  if(!st)
  {
    // compute state variables
    ALLOCOBJ0(st);
    ALLOCARR( st->LCT, n*m);
    ALLOCARR( st->L, m*n );
    ALLOCARR( st->Y2, m );
    ALLOCARR( st->YL, n );
    ALLOCARR( st->YE, m );

    MtMatCopy_fast(st->L,A,m,n);

    // T = A * S
    ALLOCARR( T,n*m );
    MtMatMult_d( T, A, S,  m, n, n );

    // T2 = T * A'
    ALLOCARR( T2,m*m );
    MtMatMult_d( T2, T, A,  m, n, -m );

    if(E)
    {
      // T2 = T2 + E
      for(i=0;i<m;i++) 
	for(j=0;j<m;j++)
	  MT_MATEL(T2,m,i,j) += MT_MATEL(E,m,i,j);
    }

    // T2 = T2^(-1)
    MtMatInverse( T2, m );

    // T = A' * T2
    MtMatMult_d( T, A, T2,  -n, m, m );

    //LCT = S * T
    MtMatMult_d( st->LCT, S, T,  n, n, m );
  }

  if( E )
  {
    // Sample from noise for stochastic (soft) constraint
    RND_GetMultNorm( rs, m, E, st->YE, &st->rmnsE, 0 );
    for(i=0;i<m;i++) B[i] -= st->YE[i];   	

  }

  // Apply transformation

  MtMatMult_d( st->Y2, st->L, X, m, n, 1 ); 
  for(i=0;i<m;i++) st->Y2[i] -= B[i];
  // YL = LCT * Y2
  MtMatMult_d( st->YL, st->LCT, st->Y2, n, m, 1 );
  // Y0 = Y0 - YL
  for(i=0;i<n;i++) X[i] -= st->YL[i];

//rcCatch:
  if(pState)
  {
    *pState = st;
  }
  else
  {
    RND_DelMultNormLECState(st);
    st=NULL;
  }
  FREEMEM(T); FREEMEM(T2);
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int RND_DelMultNormLECState( RND_MultNormLECState *s )
{
  if(s)
  {
    RND_DelMultNormState(s->rmnsE);
    FREEMEM(s->YE);
    FREEMEM(s->L);
    FREEMEM(s->LCT);
    FREEMEM(s->Y);
    FREEMEM(s->Y2);
    FREEMEM(s->YL);
    FREEMEM(s);
  }
  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
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
			)
{
#undef MAXDIM
#define MAXDIM 20
  int rcThis=0, dynAlloc=FALSE,k,i,j,L_valid=FALSE;
  double L_buf[MAXDIM*MAXDIM],vd_buf[MAXDIM],vx_buf[MAXDIM],
	 *L=NULL,*vd=NULL,*vx=NULL;
  RND_MultNormState *state=NULL;

  if(pState)
  {
    state = *pState;
    if(state) L_valid = TRUE;
  }

  if(!state && (pState || m>MAXDIM))
  {
    ALLOCOBJ0(state);
    ALLOCARR(state->L,m*m); 
    ALLOCARR( state->vd, m ); ALLOCARR( state->vx, m );
    state->m = m;
    dynAlloc=TRUE;
  }
  else if(!L_valid)
  {
    L = L_buf; vd = vd_buf; vx = vx_buf;
  }

  if(state)
  {
    m = state->m;
    L = state->L; vd = state->vd; vx = state->vx;
  }

  // Method: Cholesky decomposition of correlation matrix
  // cf. Uebersax JS. MVN program for random multivariate 
  // normal numbers. 2006. 

  if(!L_valid)
  {
    MtMatCopy_fast(L,mCorrCov,m,m);
#if 1
    cholesky(L,m,vd);
#else
    MtMatCholesky(L,m,vd);
#endif
    // clear upper triangle
    for(i=0;i<m;i++) 
      for(j=i+1;j<m;j++) MT_MATEL( L,m,i,j)=0;

    // Optionally verify the cholesky decomposition
    // mCorrCov =? LL'
    if(opt & RND_CHECK)
    {
      double f,sum,esum=0,err;
      MtStat estat;

      MT_STAT_INIT(&estat);

      for(i=0;i<m;i++)
	for(j=0;j<m;j++)
	{
	  f = MT_MATEL(mCorrCov,m,i,j);
	  for(sum=0,k=0;k<m;k++)
	    sum += MT_MATEL(L,m,i,k) *
	           MT_MATEL(L,m,j,k);
	  err = fabs(f-sum);
	  MT_STAT_UPDATE(&estat,err);
	  esum += err*err;
	}

      MT_STAT_RESULT(&estat);

      if(opt & RND_VERBOSE)
      {
	TRC1(("MLR: (%dx%d) Matrix Cholesky verification:\n",
	      m,m));
	TRC1(("Elem-Err: "));MtStatPrint(&estat);
	TRC1(("Total error: %g\n",
	      sqrt(esum)/(double)(m*m)));
      }

      if(estat.max>0.001)
      {
	TRCERR(("Numerical verification error (e=%g)\n",
	      estat.max));
	if(estat.max>0.01) rcThis=MI_ENUMERICAL;
      }
    }
  }


  if(pState) *pState = state;

  for(k=0;k<m;k++)
  {
    vx[k] = RND_GetNorm(rs,0,1);
  }

  MtMatMult_d( vOut, L, vx, m, m, 1 );

//rcCatch:
  if(dynAlloc && !pState)
  {
    FREEMEM(L); FREEMEM(vd); FREEMEM( vx );
  }
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
double RND_GetTruncNorm( RND_Stream *rs, double m, double sd,
    			 double left, double right,
			 int	*p_num_iter )
{
  //
  // Truncated normal sampling from x ~ N(mean,sigma) with 
  // left < x < right
  //
  // This method is based on paper "Simulation of truncated normal variables"
  // by Robert (1995)
  double low,up,u,prob,z,x;
  int i;

  // Standardize
  low = (left - m)/sd;
  up = (right - m)/sd;

  u =1; prob = 0; i=0;
  while(u>prob)
  {
    // z follows uniform distribution on [low, up]
    z = RND_Get(rs)*(up-low)+low;
    // compute the acceptance probability    
    if (low > 0)
    {
      prob = exp((low*low-z*z)/2);
    }
    else if( up < 0 )
    {
      prob = exp((up*up-z*z)/2);      
    }
    else 
    {
      prob = exp(-(z*z/2));
    }
    u = RND_Get(rs);
    i++;
  }
  
  x = m + z *sd;

  if(p_num_iter) *p_num_iter=i;
  return x;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int RND_Test(void)
{
  int m,nMc,seed,H[10000],i,j,rcThis=0;
  RND_Stream rs;  
  RND_Generator *prng=NULL; 
  double	u;
  
  USE(u);
  seed = time(NULL);
  prng = &RND_Generator_numc;
  RND_Init(&rs,prng,seed);

  for(m=2;m<=10;m++)
  {
    nMc=100000;
    for(i=0;i<m;i++) H[i] = 0;
    for(i=0;i<nMc;i++)
    {
      j = RND_GetInt(&rs,m);
      if(j<0||j>m-1)
      {
	TRCERRR(("index %d out of bound (max=%d)\n",j,m-1),MI_ERROR);
      }
      H[j]++;
    }
    for(i=0;i<m;i++) printf("%6.2f ",100.*H[i]/(double)nMc);
    printf("\n");
  }
rcCatch:
  return rcThis;
}

