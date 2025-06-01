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
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include "globdef.h"
#include "mtool.h"
#include "util.h"
#include "grid.h"

//#define DO_ALIGNED_MALLOC

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void G3D_GetGausslikeKernel4x4x4(
    		// even-sized kernel for cell centered filtering 
		// See also: "magic kernel" http://johncostella.webs.com/magic/
    		float *krn_data  // OUT: kernel (normalized wsum=1)
		)
{
  double krn1d[4] = {1.,3.,3.,1.},
	 krn3d[4*4*4];

  TRC3(("%s\n",UT_FUNCNAME));

  int dx,dy,dz,m=4,
      d0=-m/2+1,d1=m/2;
  double wsum=0;
  for(dz=d0;dz<=d1;dz++)
  {
    for(dy=d0;dy<=d1;dy++)
    {
      for(dx=d0;dx<=d1;dx++)
      {
	double w = krn1d[dx-d0]*krn1d[dy-d0]*krn1d[dz-d0];
	wsum += w;
	krn3d[ G3D_INDEX0( m,m*m, dx-d0, dy-d0, dz-d0 ) ] = w;
	TRC3(("%7.4f ",w));
      }
      TRC3(("\n"));
    }
    TRC3(("\n"));
  }

  UT_ASSERT0(wsum>MT_NUM_EPS);

  double wsum2=0;
  int k;
  for(k=0;k<m*m*m;k++) 
  {
    krn_data[k] = krn3d[k]/ wsum;
    wsum2 += krn_data[k];
  }

  TRC3(("wsum=%g\n",wsum2));
  UT_ASSERT0(fabs(wsum2-1)<1e-6);

}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int G3D_GetBurtAdelsonKernel(
    		double a, // IN: smoothness/shape parameter  
			 //   0.6=neg.lobes, 0.5=triangle, 
			 //   0.4=gauss, 0.3=broader gauss
			 //   0 = box
			 //   see also [Burt1983]
    		int krn_radius,  // kernel radius: 2 for gauss, 1 for triangle,...
    		double *krn_data  // OUT: kernel 
		)
{
  int rcThis=0,dx,dy,dz,r=krn_radius;
  double w0[3]={0},w,wsum;

  TRC3(("G3D_GetBurtAdelsonKernel\n"));

  w0[0] = a;
  w0[1] = 1./4;
  w0[2] = 1./4 - a/2.;

  wsum=0;
  for(dz=-r;dz<=r;dz++)
  {
    for(dy=-r;dy<=r;dy++) 
    {
      for(dx=-r;dx<=r;dx++)
      {
	w = w0[ABS(dx)]*w0[ABS(dy)]*w0[ABS(dz)];
	krn_data[G3D_INDEX0((2*r+1),((2*r+1)*(2*r+1)),dx+r,dy+r,dz+r)] = w;
	wsum += w;
	TRC3(("%7.4f ",w));
      }
      TRC3(("\n"));
    }
    TRC3(("\n"));
  }
  TRC3(("wsum=%g\n",wsum));

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static 
void G3D_IUpsampleGaussBlock5x5x5( 
    			     G3D_Grid grid,
    			     int x1, int x2,	// target range (hiRes)
			     int y1, int y2,
			     int z1, int z2,
			     float *data,
			     float *dataLapl,

			     G3D_Grid gridLo,
			     float *dataLo,

			     double *kern1
			     )
{
  int x,y,z;
  int xlmax=gridLo->sx-1, ylmax=gridLo->sy-1, zlmax=gridLo->sz-1;
  double *kern0=kern1+2;

  for(z=z1;z<=z2;z++)
    for(y=y1;y<=y2;y++)
      for(x=x1;x<=x2;x++)
      {
	int dx,dy,dz,xl,yl,zl;
	double sum=0;
	for(dz=-2;dz<=2;dz++)
	{
	  int z2=z+dz;
	  if(z2&1) continue;
	  zl = z2/2;  
	  MT_CLAMP(zl,0,zlmax);
	  for(dy=-2;dy<=2;dy++)
	  {
	    int y2=y+dy;
	    if(y2&1) continue;
	    yl = y2/2;
	    MT_CLAMP(yl,0,ylmax);
	    LongInt ii=G3D_INDEX(gridLo,0,yl,zl);
	    double wgt0=kern0[dy]*kern0[dz];
	    for(dx=-2;dx<=2;dx++)
	    {
	      int x2=x+dx;
	      if(x2&1) continue;
	      xl = x2/2;
	      MT_CLAMP(xl,0,xlmax);
	      double wgt = kern0[dx]*wgt0;
	      sum += wgt*dataLo[ii+xl];
	    }
	  }
	}
	//TRCP(("wsum=%g ",wsum));
        LongInt ii=G3D_INDEX(grid,x,y,z);
	double g = sum*8;
	if(dataLapl) 
	{
	  dataLapl[ii] = data[ii] - g;
	}
	else
	{
	  data[ii] = g;
	}
      }
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int G3D_UpSamplePyraLevelFast( 
    	G3D_Grid gridLo, float *dataLo, uint16_t *dataLoUs,
    	G3D_Grid gridHi, float *dataHi, uint16_t *dataHiUs,
	float *dataLapl, 
	int oxh, int oyh, int ozh,
        float *kernel,  // optional 5x5 upsampling kernel
	int btyp // boundary handling 
	)
{
  double krn1[5] = { 0.05, 0.25, 0.4, 0.25, 0.05 },
    	 kern[5*5*5]; 

  int rcThis=0,r=2,
      xlmax=gridLo->sx-1,ylmax=gridLo->sy-1,zlmax=gridLo->sz-1,
      xhmax=gridHi->sx-1,yhmax=gridHi->sy-1,zhmax=gridHi->sz-1;
  UtTimer tm;
  int	do_time=0;
  G3D_Grid T=NULL;

  UT_ASSERT(kernel==NULL); 
  UT_ASSERT(r==2);

  if(dataHiUs||dataLoUs)
  {
    UT_ASSERT(!(dataLo||dataHi));
  }
  if(dataLo||dataHi)
  {
    UT_ASSERT(!(dataLoUs||dataHiUs));
  }

  if(btyp==GRID_BDEFAULT) btyp=GRID_BCLIP;

  TRC3(("UpSamplePyraLevelFast: (%dx%dx%d) -> (%dx%dx%d)\n",
	gridLo->sx,gridLo->sy,gridLo->sz, gridHi->sx,gridHi->sy,gridHi->sz));

  UT_ASSERT(((xhmax+1) <= 2*(xlmax+1))&&((yhmax+1) <= 2*(ylmax+1))&&
      	    ((zhmax+1) <= 2*(zlmax+1))   );

  if((do_time = (LogLevel>=2)&&(gridHi->n>10000000))==TRUE)
  {
    TIMER_START(&tm);
  }

  int dx,dy,dz,k=0;
  for(dz=-2;dz<=2;dz++)
    for(dy=-2;dy<=2;dy++)
      for(dx=-2;dx<=2;dx++)
	kern[k++] = krn1[2+dx]*krn1[2+dy]*krn1[2+dz];

  //
  // Setup per-thread-task grid
  //
  // Each block should be 
  //   1. Small enough to fit into per-CPU-core cache
  //   2. Large enough to amortize thread-setup costs
  //

  {
    G3D_Grid G=gridHi;
    LongInt dx=MAX(MIN(G->sx,64),1),	
	    dy=MAX(MIN(G->sy/4,64),1),
	    dz=MAX(G->sz/16,1);		

    T=G3D_CreateGrid((G->sx + (dx-1)) / dx,
      		     (G->sy + (dy-1)) / dy,
		     (G->sz + (dz-1)) / dz, 0);

    T->tsx=dx; T->tsy=dy; T->tsz=dz;
  }

  int iTask;
  #pragma omp parallel for schedule(dynamic) 
  for(iTask=0;iTask<T->n;iTask++)
  {
    int x,y,z;
    G3D_IDX2COORDS(T,iTask,x,y,z);
    int x1 = x * T->tsx,
	x2 = MIN(x1 + T->tsx-1, gridHi->sx-1),
	y1 = y * T->tsy,
	y2 = MIN(y1 + T->tsy-1, gridHi->sy-1),
	z1 = z * T->tsz,
	z2 = MIN(z1 + T->tsz-1, gridHi->sz-1);

    G3D_IUpsampleGaussBlock5x5x5(
			      gridHi, 
			      x1,x2, y1,y2, z1,z2, 
			      dataHi, 
			      dataLapl,
			      gridLo, dataLo, krn1);

  }

  if(do_time)
  {
    TIMER_STOP(&tm);
    TRC3(("CPU = %.1f seconds (%.0f pix/sec)\n", 
	  (double)(TIMER_DIFF_MS(&tm)/1000.0), 
	  (gridHi->n)/(double)(TIMER_DIFF_MS(&tm)/1000.0)));
  }

rcCatch:
    G3D_DeleteGrid(&T);
    return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int G3D_UpSamplePyraLevel( 
    	G3D_Grid bmpLo, float *dataLo, uint16_t *dataLoUs,
    	G3D_Grid bmpHi, float *dataHi, uint16_t *dataHiUs,
	float *dataLapl, 
	int oxh, int oyh, int ozh,
        float *kernel,  // optional 5x5 upsampling kernel
	int btyp // boundary handling 
	)
{
  double krn1[5] = { 0.05, 0.25, 0.4, 0.25, 0.05 }; 
  int rcThis=0,r=2,xh0,yh0,zh0,
      xlmax=bmpLo->sx-1,ylmax=bmpLo->sy-1,zlmax=bmpLo->sz-1,
      xhmax=bmpHi->sx-1,yhmax=bmpHi->sy-1,zhmax=bmpHi->sz-1;
  UtTimer tm;
  int	do_time=0;

  UT_ASSERT(kernel==NULL); 
  UT_ASSERT(r==2);

  if(dataHiUs||dataLoUs)
  {
    UT_ASSERT(!(dataLo||dataHi));
  }
  if(dataLo||dataHi)
  {
    UT_ASSERT(!(dataLoUs||dataHiUs));
  }

  if(btyp==GRID_BDEFAULT) btyp=GRID_BCLIP;

  TRC3(("UpSamplePyraLevel: (%dx%dx%d) -> (%dx%dx%d)\n",
	bmpLo->sx,bmpLo->sy,bmpLo->sz, bmpHi->sx,bmpHi->sy,bmpHi->sz));

  UT_ASSERT(((xhmax+1) <= 2*(xlmax+1))&&((yhmax+1) <= 2*(ylmax+1))&&
      	    ((zhmax+1) <= 2*(zlmax+1))   );

  if((do_time = (LogLevel>=2)&&(bmpHi->n>10000000))==TRUE)
  {
    TIMER_START(&tm);
  }

  #pragma omp parallel for private(xh0,yh0,zh0) \
  			   schedule(dynamic)
  for(zh0=0;zh0<bmpHi->sz;zh0++)
  {
    for(yh0=0;yh0<bmpHi->sy;yh0++)
    {
      for(xh0=0;xh0<bmpHi->sx;xh0++)
      {
	int dx,dy,dz, 
	    xh,yh,zh, xh2,yh2,zh2, 
	    xl,yl,zl, xl2,yl2,zl2;
	LongInt po;
	double wsum,gsum,g=0;

	xh = xh0-oxh; xl = xh/2;
	yh = yh0-oyh; yl = yh/2; 
	zh = zh0-ozh; zl = zh/2; 
	{
	  wsum=gsum=0;
	  for(dz=-r;dz<=r;dz++)
	    for(dy=-r;dy<=r;dy++)
	      for(dx=-r;dx<=r;dx++)
	      {
		double wgt = krn1[2+dx]*krn1[2+dy]*krn1[2+dz];

		xh2 = xh-dx; 
		yh2 = yh-dy; 
		zh2 = zh-dz;

		xl2=xh2/2; 
		yl2=yh2/2; 
		zl2=zh2/2;

//		if(xh2&1 || yh2&1 || zh2&1) continue;
		if(2*xl2!=xh2 || 2*yl2!=yh2 || 2*zl2!=zh2) continue;            

#if 1
		G3D_CLIP_CONST(0,0,0, bmpLo->sx-1,bmpLo->sy-1,bmpLo->sz-1,
		    		xl2,yl2,zl2);
#else
		if(!(G3D_IN_RANGE(bmpLo,xl2,yl2,zl2))) continue;
#endif

		po = G3D_INDEX(bmpLo,xl2,yl2,zl2);
		g = dataLo ? dataLo[po] : dataLoUs[po];
		gsum += wgt*g;
		wsum += wgt;
	      }
	  g = gsum/wsum;
	  po = G3D_INDEX(bmpHi,xh0,yh0,zh0);
	  if(dataLapl)
	  {
	    float gh = dataHi?dataHi[po]:dataHiUs[po];
	    dataLapl[po] = gh-g;
	  }
	  else
	  {
	    if(dataHi)
	      dataHi[po]=g;
	    else
	      dataHiUs[po]=g;
	  }
	}
      }
    }
  }

  if(do_time)
  {
    TIMER_STOP(&tm);
    TRC3(("CPU = %.1f seconds (%.0f pix/sec)\n", 
	  (double)(TIMER_DIFF_MS(&tm)/1000.0), 
	  (bmpHi->n)/(double)(TIMER_DIFF_MS(&tm)/1000.0)));
  }

rcCatch:
    return rcThis;
}

#if 0
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int G3D_UpSamplePyramid( G3D_Grid *G, int chan, int ichan,
    		        int l_from, int l_to )
{
  int rcThis=0,l,rc;

  for(l=l_from-1;l>=l_to;l--)
  {
    G3D_Grid GL = G->pyramid.gaussian[l+1],
	      GH = G->pyramid.gaussian[l];
    UT_ASSERT(chan==GRID_FLOAT);
    float *dataHi=GH->dataFloat[ichan],
	  *dataLo=GL->dataFloat[ichan];
    
    rc = G3D_UpSamplePyraLevelFast( 
		      gridLo, dataLo, NULL,
		      gridHi, dataHi, NULL,
		      NULL,  
		      0, 0, 0, NULL,
		      GRID_BDEFAULT );
    if(rc) raiseRc(rc);
  }
  
rcCatch:
  return rcThis;
}
#endif

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static
int G3D_ICreateLaplacianPyramid( 
    		      G3D_Grid G,  	 
		      int chan, int ichan,		 
		      float *kern,
		      unsigned opt )
{
  int rcThis=0,rc,l;
  G3D_Pyramid *P=&G->pyramid;
  G3D_Grid gridHi, gridLo, gridLap;
  float *dataLap;
  float *dataLo, *dataHi;

  TRC(("G3D_ICreateLaplacianPyramid(%d/%d)\n",chan,ichan));

  UT_ASSERT(!(opt & GRID_FTRIANGLE));
  UT_ASSERT(!(opt & GRID_FASTMIPMAP));
  UT_ASSERT(!(opt & GRID_CELLCENTERED));

  for(l=P->n_level-1; l>=0; l--)
  {
    gridLo = l<P->n_level-1 ? P->gaussian[l+1] : NULL;
    gridHi = P->gaussian[l];
    UT_ASSERT(gridHi);
    gridLap = P->laplacian[l];
    if(!gridLap || (!G3D_DIM_EQUAL(gridLap,gridHi)))
    {
      G3D_DeleteGrid(&gridLap);
      P->laplacian[l] = gridLap = G3D_CreateGrid(gridHi->sx,gridHi->sy,gridHi->sz,GRID_GEN);
    }
    UT_ASSERT(gridLap);

    UT_ASSERT(chan==GRID_FLOAT);
    UT_ASSERT(kern==NULL);

    dataHi = G3D_GetFloatChannel(gridHi,ichan,0);
    dataLap = G3D_GetFloatChannel(gridLap,ichan,0);

    UT_ASSERT(dataHi);  
    UT_ASSERT(dataLap);

    if(!gridLo)
    {
      // max-level: simply copy gaussian
      LongInt ii;
      #pragma omp parallel for
      for(ii=0;ii<gridHi->n;ii++) dataLap[ii] = dataHi[ii];
    }
    else
    {
      dataLo = G3D_GetFloatChannel(gridLo,ichan,0);
      UT_ASSERT((dataLo));  
    
      rc = G3D_UpSamplePyraLevelFast( 
	  		gridLo, dataLo, NULL,
	  		gridHi, dataHi, NULL,
	  		dataLap,  
	  		0, 0, 0, kern,
	  		GRID_BDEFAULT );
      if(rc) raiseRc(rc);

#if 0	
      {
        float *dataLapRef = G3D_GetFloatChannel(gridLap,ichan+1,0);
	rc = G3D_UpSamplePyraLevel( 
			gridLo, dataLo, NULL,
			gridHi, dataHi, NULL,
			dataLapRef,  
			0, 0, 0, kern,
			GRID_BDEFAULT );
	if(rc) raiseRc(rc);

	int x,y,z;
	for(z=0;z<gridHi->sz;z++)
	  for(y=0;y<gridHi->sy;y++)
	    for(x=0;x<gridHi->sx;x++)
	{
	  LongInt i=G3D_INDEX(gridHi,x,y,z);
	  if(fabs(dataLap[i]-dataLapRef[i])>0.001)
	  {
	    printf("dataHi(l=%d:%d,%d,%d): %g %g\n",l,x,y,z,dataLap[i],dataLapRef[i]);
	  }
	}
      }
#endif

    }
  }

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static inline
void G3D_IDownsampleGaussLine( int x1, int x2, 		// target range 
    			     LongInt stride,
		             float *data,             // start of target line/column

			     int minHi, int maxHi, 	// max. source range
    			     LongInt strideHi,
    		             float *dataHi,		// start of line/column 

			     double *kern5		// 5 kernel coefficients 
			     )
{
  double *kernCent=kern5+2;
  int i;

  for(i=x1;i<=x2;i++)
  {
    int k;
    double sum=0;
    for(k=-2;k<=2;k++)
    {
      int j=2*i+k;
      if(unlikely(j<minHi)) j=minHi;
      if(unlikely(j>maxHi)) j=maxHi;
      sum += kernCent[k]*dataHi[j*strideHi];     
    }
    data[i*stride] = sum;
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static 
void G3D_IDownsampleBlock4x4x4( 
    			     G3D_Grid grid,
    			     int x1, int x2,	// target range
			     int y1, int y2,
			     int z1, int z2,
			     float *data,

			     G3D_Grid gridHi,
			     float *dataHi,

			     double *kern4x4x4
			     )
{
  int x,y,z;
  int xhmax=gridHi->sx-1, yhmax=gridHi->sy-1, zhmax=gridHi->sz-1;

  for(z=z1;z<=z2;z++)
  for(y=y1;y<=y2;y++)
  for(x=x1;x<=x2;x++)
  {
    int dx,dy,dz,xh2,yh2,zh2,
	xh=2*x,yh=2*y,zh=2*z;
    double sum=0;
    int k=0;
    for(dz=-1;dz<=2;dz++)
    {
      zh2=zh+dz;
      MT_CLAMP(zh2,0,zhmax);
      for(dy=-1;dy<=2;dy++)
      {
	yh2=yh+dy;
	MT_CLAMP(yh2,0,yhmax);
	LongInt ii=G3D_INDEX(gridHi,0,yh2,zh2);
	for(dx=-1;dx<=2;dx++)
	{
	  double wgt = kern4x4x4[k++];
	  xh2=xh+dx;  
	  MT_CLAMP(xh2,0,xhmax);
	  sum += wgt*dataHi[ii+xh2];
	  //if(x>x1+3&&y>y1+3&&z>z1+3) TRCP(("%7.4f ",wgt));
	}
	  //if(x>x1+3&&y>y1+3&&z>z1+3) TRCP(("\n"));
      }
	  //if(x>x1+3&&y>y1+3&&z>z1+3) TRCP(("\n"));
    }
    data[G3D_INDEX(grid,x,y,z)]=sum;
  }
}



/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static 
void G3D_IDownsampleGaussBlock5x5x5( 
    			     G3D_Grid grid,
    			     int x1, int x2,	// target range
			     int y1, int y2,
			     int z1, int z2,
			     float *data,

			     G3D_Grid gridHi,
			     float *dataHi,

			     double *kern5x5x5
			     )
{
  int x,y,z;
#if 0
  // X-direction
  for(z=z1;z<=z2;z++)
    for(y=y1;y<=y2;y++)
     G3D_IDownsampleGaussLine(
	 x1,x2,1, data+G3D_INDEX(grid,0,y,z),
	 0,gridHi->sx-1,1, dataHi+G3D_INDEX(gridHi,0,2*y,2*z), kern5);

  // Y-direction
  for(z=z1;z<=z2;z++)
    for(x=x1;x<=x2;x++)
     G3D_IDownsampleGaussLine(
	 y1,y2,grid->sx, data+G3D_INDEX(grid,x,0,z),
	 0,gridHi->sy-1,gridHi->sx, dataHi+G3D_INDEX(gridHi,2*x,0,2*z), kern5);

  // Z-direction
  for(y=y1;y<=y2;y++)
    for(x=x1;x<=x2;x++)
     G3D_IDownsampleGaussLine(
	 z1,z2,grid->sxy, data+G3D_INDEX(grid,x,y,0),
	 0,gridHi->sz-1,gridHi->sxy, dataHi+G3D_INDEX(gridHi,2*x,2*y,0), kern5);
#else
  int xhmax=gridHi->sx-1, yhmax=gridHi->sy-1, zhmax=gridHi->sz-1;

  for(z=z1;z<=z2;z++)
    for(y=y1;y<=y2;y++)
      for(x=x1;x<=x2;x++)
      {
	int dx,dy,dz,xh2,yh2,zh2,xh,yh,zh;
	xh = 2*x;
	yh = 2*y; 
	zh = 2*z; 
	double sum=0;
	int k=0;
#if 0
	for(dz=-2;dz<=2;dz++)
	  for(dy=-2;dy<=2;dy++)
	    for(dx=-2;dx<=2;dx++)
	    {
	      double wgt = kern5x5x5[k++];
	      xh2=xh+dx; yh2=yh+dy; zh2=zh+dz;
	      G3D_CLIP_CONST(0,0,0, xhmax,yhmax,zhmax, xh2,yh2,zh2);
	      sum += wgt*dataHi[G3D_INDEX(gridHi,xh2,yh2,zh2)];
	    }
#else
	for(dz=-2;dz<=2;dz++)
	{
	  zh2=zh+dz;
	  MT_CLAMP(zh2,0,zhmax);
	  for(dy=-2;dy<=2;dy++)
	  {
	    yh2=yh+dy;
	    MT_CLAMP(yh2,0,yhmax);
	    LongInt ii=G3D_INDEX(gridHi,0,yh2,zh2);
	    for(dx=-2;dx<=2;dx++)
	    {
	      double wgt = kern5x5x5[k++];
	      xh2=xh+dx;  
	      MT_CLAMP(xh2,0,xhmax);
	      sum += wgt*dataHi[ii+xh2];
	    }
	  }
	}
#endif
	data[G3D_INDEX(grid,x,y,z)]=sum;
      }
#endif
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void G3D_DownsampleBlock5x5x5( 
    			     G3D_Grid grid,
    			     int x1, int x2,	// target range
			     int y1, int y2,
			     int z1, int z2,
			     float *data,

			     G3D_Grid gridHi,
			     float *dataHi,

			     double *kern5x5x5
			     )

{
  G3D_IDownsampleGaussBlock5x5x5( grid, x1,x2, y1,y2, z1,z2, data,
      				  gridHi, dataHi,
				  kern5x5x5 );      				  
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static 
void G3D_IDownsampleGaussBlock( 
    			     G3D_Grid grid,
    			     int x1, int x2,	// target range
			     int y1, int y2,
			     int z1, int z2,
			     float *data,

			     G3D_Grid gridHi,
			     float *dataHi,

			     int kern_r,
			     double *kern
			     )
{
  int x,y,z,r=kern_r;
  int xhmax=gridHi->sx-1, yhmax=gridHi->sy-1, zhmax=gridHi->sz-1;

  for(z=z1;z<=z2;z++)
    for(y=y1;y<=y2;y++)
      for(x=x1;x<=x2;x++)
      {
	int dx,dy,dz,xh2,yh2,zh2,xh,yh,zh;
	xh = 2*x;
	yh = 2*y; 
	zh = 2*z; 
	double sum=0;
	int k=0;
	for(dz=-r;dz<=r;dz++)
	{
	  zh2=zh+dz;
	  MT_CLAMP(zh2,0,zhmax);
	  for(dy=-r;dy<=r;dy++)
	  {
	    yh2=yh+dy;
	    MT_CLAMP(yh2,0,yhmax);
	    LongInt ii=G3D_INDEX(gridHi,0,yh2,zh2);
	    for(dx=-r;dx<=r;dx++)
	    {
	      double wgt = kern[k++];
	      xh2=xh+dx;  
	      MT_CLAMP(xh2,0,xhmax);
	      sum += wgt*dataHi[ii+xh2];
	    }
	  }
	}
	data[G3D_INDEX(grid,x,y,z)]=sum;
      }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static 
void G3D_IDownsampleMipmapBlock( 
    			     G3D_Grid grid,
    			     int x1, int x2,	// target range
			     int y1, int y2,
			     int z1, int z2,
			     float *data,

			     G3D_Grid gridHi,
			     float *dataHi,

			     double *kern5x5x5
			     )
{
  int x,y,z;
  int xhmax=gridHi->sx-1, yhmax=gridHi->sy-1, zhmax=gridHi->sz-1;

  for(z=z1;z<=z2;z++)
    for(y=y1;y<=y2;y++)
      for(x=x1;x<=x2;x++)
      {
	int dx,dy,dz,xh2,yh2,zh2,xh,yh,zh;
	xh = 2*x;
	yh = 2*y; 
	zh = 2*z; 
	double sum=0;
	for(dz=0;dz<=1;dz++)
	{
	  zh2=zh+dz;
	  MT_CLAMP(zh2,0,zhmax);
	  for(dy=0;dy<=1;dy++)
	  {
	    yh2=yh+dy;
	    MT_CLAMP(yh2,0,yhmax);
	    LongInt ii=G3D_INDEX(gridHi,0,yh2,zh2);
	    for(dx=0;dx<=1;dx++)
	    {
	      xh2=xh+dx;  
	      MT_CLAMP_MAX(xh2,xhmax);
	      sum += dataHi[ii+xh2];
	    }
	  }
	}
	data[G3D_INDEX(grid,x,y,z)]=sum*0.125;
      }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int G3D_DownSamplePyraLevelFast( 
	  G3D_Grid gridLo, float *dataLo, uint16_t *dataLoUs, 
	  G3D_Grid gridHi, float *dataHi, uint16_t *dataHiUs, 
	  int oxh, int oyh, int ozh, // optional gridHi offset 0, or -1 (for MCM)  
          float *kernel,  // optional 5x5 downsampling kernel
	  int btyp,	 // boundary handling
	  unsigned opt
	  )
{  
  // 
  // TODO: optimize  
  //   - separable 1D loops over cache/core blocked sub-grids
  //   - SIMD
  //
  double krn1[5] = { 0.05, 0.25, 0.4, 0.25, 0.05 },
  	 kern[5*5*5]; 
  int rcThis=0,rc;
  G3D_Grid T=NULL;

  UT_ASSERT(kernel==NULL);

  UT_ASSERT(oxh==0&&oyh==0&&ozh==0);

  UtTimer tm;
  int	do_time=0;

  if(dataHiUs||dataLoUs)
  {
    UT_ASSERT(!(dataLo||dataHi));
  }
  if(dataLo||dataHi)
  {
    UT_ASSERT(!(dataLoUs||dataHiUs));
  }

  UT_ASSERT(dataHi||dataHiUs);
  UT_ASSERT(dataLo||dataLoUs);

  if(btyp==GRID_BDEFAULT) btyp=GRID_BCLIP;

  {
    int     xlmax=gridLo->sx-1,ylmax=gridLo->sy-1,zlmax=gridLo->sz-1,
	    xhmax=gridHi->sx-1,yhmax=gridHi->sy-1,zhmax=gridHi->sz-1;
    UT_ASSERT(((xhmax+1) <= 2*(xlmax+1))&&((yhmax+1) <= 2*(ylmax+1))&&
	((zhmax+1) <= 2*(zlmax+1)));
  }


  if(opt & GRID_FTRIANGLE)
  {
    // simple 3x3x3 triangle kernel 
    rc = G3D_GetBurtAdelsonKernel(0.5,1,kern);
    if(rc) raiseRc(rc);
  }
  else if( opt & GRID_CELLCENTERED )
  {
    float tmp[4*4*4];
    G3D_GetGausslikeKernel4x4x4(tmp);
    int k;
    for(k=0;k<ARRAY_LENGTH(tmp);k++) kern[k] = tmp[k];
  }
  else
  {
    int dx,dy,dz,k=0;
    for(dz=-2;dz<=2;dz++)
      for(dy=-2;dy<=2;dy++)
	for(dx=-2;dx<=2;dx++)
	  kern[k++] = krn1[2+dx]*krn1[2+dy]*krn1[2+dz];
  }


  do_time = (LogLevel>=2)&&
		(gridHi->n>1000000000);

  if(do_time)
  {
    TRC3(("G3D_DownSamplePyraLevelFast (opt=%d): (%dx%dx%d) -> (%dx%dx%d)\n",
	opt,gridHi->sx,gridHi->sy,gridHi->sz,gridLo->sx,gridLo->sy,gridLo->sz));
  }

  if(do_time)
  {
    TIMER_START(&tm);
  }

  UT_ASSERT(btyp==GRID_BCLIP);

  //
  // Setup per-thread-task grid
  //
  // Each block should be 
  //   1. Small enough to fit into per-CPU-core cache
  //   2. Large enough to amortize thread-setup costs
  //

  {
    G3D_Grid G=gridLo;
    LongInt dx=MAX(MIN(G->sx,64),1),	
	    dy=MAX(MIN(G->sy/4,64),1),
	    dz=MAX(G->sz/16,1);		

    T=G3D_CreateGrid((G->sx + (dx-1)) / dx,
      		     (G->sy + (dy-1)) / dy,
		     (G->sz + (dz-1)) / dz, 0);

    T->tsx=dx; T->tsy=dy; T->tsz=dz;
  }

#if 1
  int iTask;
  #pragma omp parallel for schedule(dynamic) 
  for(iTask=0;iTask<T->n;iTask++)
  {
    int x,y,z;
    G3D_IDX2COORDS(T,iTask,x,y,z);
    int x1 = x * T->tsx,
	x2 = MIN(x1 + T->tsx-1, gridLo->sx-1),
	y1 = y * T->tsy,
	y2 = MIN(y1 + T->tsy-1, gridLo->sy-1),
	z1 = z * T->tsz,
	z2 = MIN(z1 + T->tsz-1, gridLo->sz-1);
#if  0
    if(do_time)
    {
      TRCP(("thread %d: iTask=%d/%d x=%d-%d, y=%d-%d, z=%d-%d\n",
	    omp_get_thread_num(),iTask,T->n,x1,x2,y1,y2,z1,z2));
    }
#endif

    if(opt & GRID_FASTMIPMAP)
    {
      G3D_IDownsampleMipmapBlock(
		      gridLo, 
		      x1,x2, y1,y2, z1,z2, 
		      dataLo, 
		      gridHi, dataHi, kern);
    }
    else if(opt & GRID_FTRIANGLE)
    {
      G3D_IDownsampleGaussBlock(
		      gridLo, 
		      x1,x2, y1,y2, z1,z2, 
		      dataLo, 
		      gridHi, dataHi, 
		      1, kern );
    }
    else if( opt & GRID_CELLCENTERED )
    {
      G3D_IDownsampleBlock4x4x4(
		      gridLo, 
		      x1,x2, y1,y2, z1,z2, 
		      dataLo, 
		      gridHi, dataHi, kern);
    }
    else
    {
      G3D_IDownsampleGaussBlock5x5x5(
		      gridLo, 
		      x1,x2, y1,y2, z1,z2, 
		      dataLo, 
		      gridHi, dataHi, kern);
    }
  }
#else
    G3D_IDownsampleGaussBlock(
			      gridLo, 
			      0,gridLo->sx-1, 
			      0,gridLo->sy-1, 
			      0,gridLo->sz-1, 			     
			      dataLo, 
			      gridHi, dataHi, kern);
#endif
  if(do_time)
  {
    TIMER_STOP(&tm);
    TRC(("G3D_DownSamplePyraLevelFast: CPU = %.1f seconds (%.0f pix/sec)\n", 
	  (double)(TIMER_DIFF_MS(&tm)/1000.0), 
	  (gridHi->n)/(double)(TIMER_DIFF_MS(&tm)/1000.0)));
  }

rcCatch:
  G3D_DeleteGrid(&T);
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int G3D_DownSamplePyraLevel( 
	  G3D_Grid bmpLo, float *dataLo, uint16_t *dataLoUs, 
	  G3D_Grid bmpHi, float *dataHi, uint16_t *dataHiUs, 
	  int oxh, int oyh, int ozh, // optional bmpHi offset 0, or -1 (for MCM)  
          float *kernel,  // optional 5x5 downsampling kernel
	  int btyp	 // boundary handling
	  )
{  
  // 
  // TODO: optimize  
  //   - separable 1D loops over cache/core blocked sub-grids
  //   - SIMD
  //
  double krn1[5] = { 0.05, 0.25, 0.4, 0.25, 0.05 }; 
  int rcThis=0,r=2,xl,yl,zl,
      xlmax=bmpLo->sx-1,ylmax=bmpLo->sy-1,zlmax=bmpLo->sz-1,
      xhmax=bmpHi->sx-1,yhmax=bmpHi->sy-1,zhmax=bmpHi->sz-1;

  UT_ASSERT(kernel==NULL);
  UtTimer tm;
  int	do_time=0;

  if(dataHiUs||dataLoUs)
  {
    UT_ASSERT(!(dataLo||dataHi));
  }
  if(dataLo||dataHi)
  {
    UT_ASSERT(!(dataLoUs||dataHiUs));
  }

  UT_ASSERT(dataHi||dataHiUs);
  UT_ASSERT(dataLo||dataLoUs);

  if(btyp==GRID_BDEFAULT) btyp=GRID_BCLIP;

  r = 2;

  if((do_time = (LogLevel>=2)&&
		(bmpHi->n>10000000)));

  if(do_time)
  {
    TRC(("G3D_DownSamplePyraLevel: (%dx%dx%d) -> (%dx%dx%d)\n",
	bmpHi->sx,bmpHi->sy,bmpHi->sz,bmpLo->sx,bmpLo->sy,bmpLo->sz));
  }

  UT_ASSERT(((xhmax+1) <= 2*(xlmax+1))&&((yhmax+1) <= 2*(ylmax+1))&&
      	    ((zhmax+1) <= 2*(zlmax+1)));

  if(do_time) 
  {
    TIMER_START(&tm);
  }

  UT_ASSERT(btyp==GRID_BCLIP);
  UT_ASSERT(r==2);

  #pragma omp parallel for private(xl,yl,zl) \
  			   schedule(dynamic)
  for(zl=0;zl<bmpLo->sz;zl++)
  {
    for(yl=0;yl<bmpLo->sy;yl++)
    {
      for(xl=0;xl<bmpLo->sx;xl++)
      {
	int dx,dy,dz,xh2,yh2,zh2,xh,yh,zh;
	LongInt po;
	double wsum,gsum,g;
	xh = 2*xl+oxh;
	yh = 2*yl+oyh; 
	zh = 2*zl+ozh; 
	wsum=gsum=0;
	for(dz=-r;dz<=r;dz++)
	  for(dy=-r;dy<=r;dy++)
	    for(dx=-r;dx<=r;dx++)
	    {
	      double wgt = krn1[2+dx]*krn1[2+dy]*krn1[2+dz];
	      xh2=xh+dx; yh2=yh+dy; zh2=zh+dz;

//	      if(!(G3D_IN_RANGE(bmpHi,xh2,yh2,zh2))) continue;
//	      wgt=1./(dx+0.001);

	      G3D_CLIP_CONST(0,0,0, bmpHi->sx-1,bmpHi->sy-1,bmpHi->sz-1,
		  	     xh2,yh2,zh2);	      

	      po = G3D_INDEX(bmpHi,xh2,yh2,zh2);
	      g = dataHi ? dataHi[po] : dataHiUs[po]; 
	      gsum += wgt*g;
	      wsum += wgt;
	    }
	g = gsum/wsum;
	po = G3D_INDEX(bmpLo,xl,yl,zl);
	if(dataLo)
	  dataLo[po] = g;
	else 
	  dataLoUs[po] = g;
      }
    }
  }
  
  if(do_time)
  {
    TIMER_STOP(&tm);
    TRC(("CPU = %.1f seconds (%.0f pix/sec)\n", 
	  (double)(TIMER_DIFF_MS(&tm)/1000.0), 
	  (bmpHi->n)/(double)(TIMER_DIFF_MS(&tm)/1000.0)));
  }

rcCatch:
    return rcThis;
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static
int G3D_ICreateGaussianPyramid( G3D_Grid G,  	  	// target bitmap
    			       int lmin,
    			       int lmax,
		      	       int chan, int ichan,	// data channel
		      	       float *kern,
    		      	       unsigned opt )
{
  int rcThis=0,l,sxl,syl,szl,done,rc;  
  G3D_Pyramid *P=&G->pyramid;
  G3D_Grid gridHi,gridLo;
  float *dataLo, *dataHi;
  int sxyzMin=0;

  if(lmax<-1)
  {
    sxyzMin = -lmax;
    lmax = -1;
  }

  TRC(("G3D_ICreateGaussianPyramid(%d/%d,%p)\n",chan,ichan,kern));

  if(lmin==-1) lmin=0;

  P->gaussian[0] = G;

  for(l=lmin+1,done=FALSE;!done&&(lmax!=0);l++)
  {
    gridHi = P->gaussian[l-1];

    //
    // calculate new bitmap dimensions
    //
    sxl = (gridHi->sx+1) / 2;
    syl = (gridHi->sy+1) / 2;
    szl = (gridHi->sz+1) / 2;
    if(  ((sxl<=1 && syl<=1 && szl<=1 ) ||  (lmax>0&&l>=lmax)))
    {
      lmax = l;
      done = TRUE;
    }

    if( sxyzMin )
    {
      if( MIN(MIN(sxl,syl),szl) < sxyzMin )
      {
	lmax=l-1;
	break;
      }
    }

    // 
    // prepare bitmap & data channels
    //
    gridLo = P->gaussian[l];
    if(!gridLo || (!(gridLo->sx==sxl&&gridLo->sy==syl&&gridLo->sz==szl)))
    {
      G3D_DeleteGrid(&gridLo);
      P->gaussian[l] = gridLo = G3D_CreateGrid(sxl,syl,szl,GRID_GEN);
    }
    UT_ASSERT(gridLo);
    
    if((opt&GRID_VERBOSE)||LogLevel>=2)
    {
      if(l==lmin+1)
      {
	TRC(("G3D Pyramid\n"));
	TRC(("-----------------------------\n"));
        TRC(("level %2d: (%dx%dx%d)\n",
	      l-1,(int)gridHi->sx,(int)gridHi->sy,(int)gridHi->sz));
      }

      TRC(("  level %2d: (%dx%dx%d)\n",
	    l,(int)gridLo->sx,(int)gridLo->sy,(int)gridLo->sz));
    }

    if(chan != GRID_GEN)
    {
      UT_ASSERT(chan==GRID_FLOAT);
      UT_ASSERT(kern==NULL);
      dataHi = G3D_GetFloatChannel(gridHi,ichan,0);
      dataLo = G3D_GetFloatChannel(gridLo,ichan,0);
      UT_ASSERT((dataHi&&dataLo));  

#if 1

      // Downsampling
      #if 0

      rc = G3D_DownSamplePyraLevel(
		gridLo,dataLo, NULL,
		gridHi,dataHi, NULL, 
		0,0,0,
		kern, GRID_BDEFAULT );
      #else

      rc = G3D_DownSamplePyraLevelFast(
			gridLo,dataLo, NULL,
			gridHi,dataHi, NULL, 
			0,0,0,
			kern, GRID_BDEFAULT,
			opt );
      #endif

      if(rc) raiseRc(rc);    

#else
      {
        float *dataLoRef = G3D_GetFloatChannel(gridLo,ichan+1,0);
	rc = G3D_DownSamplePyraLevel(
		    gridLo,dataLoRef, NULL,
		    gridHi,dataHi, NULL, 
		    0,0,0,
		    kern, GRID_BDEFAULT );
	if(rc) raiseRc(rc);

	rc = G3D_DownSamplePyraLevelFast(
		    gridLo,dataLo, NULL,
		    gridHi,dataHi, NULL, 
		    0,0,0,
		    kern, GRID_BDEFAULT );
	if(rc) raiseRc(rc);

	int x,y,z;
	for(z=0;z<gridLo->sz;z++)
	  for(y=0;y<gridLo->sy;y++)
	    for(x=0;x<gridLo->sx;x++)
	{
	  LongInt i=G3D_INDEX(gridLo,x,y,z);
	  if(fabs(dataLo[i]-dataLoRef[i])>0.001)
	  {
	    printf("dataLo(%d,%d,%d): %g %g\n",x,y,z,dataLo[i],dataLoRef[i]);
	  }
	}
      }

#endif



      if( opt & GRID_FREE )
      {
	G3D_FreeFloatChannel(gridHi,ichan);
      }
    }
  }

  P->n_level = lmax+1;

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int G3D_FreePyramidChannels( G3D_Grid G,  	  // target bitmap
		      	    int chan,	int ichan,	  // data channel
    		            unsigned opt )
{
  int l,rcThis=0;
  G3D_Pyramid *P=&G->pyramid;


  // free the gaussian levels (except level 0)
  if(opt & GRID_GAUSS)
  {
    for(l=0; l<P->n_level; l++)
    {
      G3D_Grid bmp=P->gaussian[l];
      if(bmp && bmp!=G)
      {
        UT_ASSERT(chan==GRID_FLOAT);
	G3D_FreeFloatChannel(bmp,ichan);
      }
    }
  }
  else
  {
    UT_ASSERT(FALSE);//TODO
  }

rcCatch:
  return rcThis;
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int G3D_CreatePyramid( G3D_Grid G,  	    // target bitmap
		      int lmin,int lmax,    // min,max level (-1=auto)
		      int chan,	int ichan,  // data channel
		      float *kern_down,  // optional: 5x5x5 kernel
		      float *kern_up,    // optional
		      int btyp,
    		      unsigned opt )
{
  int rcThis=0,rc;
  G3D_Pyramid *P=&G->pyramid;

  UT_ASSERT(chan==GRID_GEN||chan==GRID_FLOAT);
  UT_ASSERT(kern_down==NULL&&kern_up==NULL);

  P->res_0 = 1;
  P->res_act = 1;
  P->res_factor = 2;
  P->krn = NULL;
  P->krn_r = 0;
  P->laplace_valid = 0;

  P->opt = opt;
  P->grid0 = G;

  rc = G3D_ICreateGaussianPyramid(G,lmin,lmax,chan,ichan,
      				 kern_down,opt);
  if(rc) raiseRc(rc);
  
  if(opt&GRID_LAPLACE)
  {
    UT_ASSERT(lmin==-1);
    rc = G3D_ICreateLaplacianPyramid(G,chan,ichan,kern_up,opt);
    if(rc) raiseRc(rc);
    if(!(opt&GRID_GAUSS))
    {
      G3D_FreePyramidChannels(G,chan,ichan,GRID_GAUSS);
    }
  }

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
float *G3D_GetFloatChannel0( G3D_Grid grid, int i, unsigned opt,
    			     char *mmTag )
{
  int rcThis=0;
  size_t szNeeded = grid->n*sizeof(grid->dataFloat[0][0]);

  UT_ASSERT0(opt==0);
  UT_ASSERT0(i>=0&&i<G3D_MAX_FCHAN);
  if(!grid->dataFloat[i] && 
     !grid->dataFloatIsExt[i])
  {
#ifdef DO_ALIGNED_MALLOC
    UT_ASSERT_FATAL(FALSE); // TODO: ensure that under all circumstances the correct variant of FREEMEM/FREEMEM_ALIGNED is called ...*/ 
    ALLOCMEM_ALIGNED_(grid->dataFloat[i], float, szNeeded, CPU_CACHELINE_SZ); 
#else
    ALLOCMEM_TAG( grid->dataFloat[i], szNeeded, mmTag);  
#endif
  }
rcCatch:
  return rcThis==0 ? grid->dataFloat[i] : NULL;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
uint16_t *G3D_GetUshortChannel( G3D_Grid grid, int i, unsigned opt )
{
  int rcThis=0;
  size_t szNeeded = grid->n*sizeof(grid->dataUshort[0][0]);

  UT_ASSERT0(i>=0&&i<G3D_MAX_FCHAN);
  if(!grid->dataUshort[i])
  {
    ALLOCMEM( grid->dataUshort[i], szNeeded);  
  }
rcCatch:
  return rcThis==0 ? grid->dataUshort[i] : NULL;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int G3D_SetFloatChannel( G3D_Grid grid, int i, float *data )
{
  int rcThis=0;

  UT_ASSERT0(i>=0&&i<G3D_MAX_FCHAN);
  UT_ASSERT_FATAL(grid->dataFloat[i]==NULL||
                  grid->dataFloatIsExt[i] );            

  grid->dataFloat[i] = data;
  grid->dataFloatIsExt[i] = data ? TRUE : FALSE;

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void G3D_FreeUshortChannel( G3D_Grid grid, int i )
{
  UT_ASSERT0(i>=0&&i<G3D_MAX_FCHAN);
  if(grid)
  {
    FREEMEM(grid->dataUshort[i]);
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void G3D_FreeFloatChannel( G3D_Grid grid, int i )
{
  UT_ASSERT0(i>=0&&i<G3D_MAX_FCHAN);
  if(grid)
  {
    if(!grid->dataFloatIsExt[i] )
    {
#ifdef DO_ALIGNED_MALLOC
    UT_ASSERT_FATAL(FALSE); // TODO: ensure that under all circumstances the correct variant of FREEMEM/FREEMEM_ALIGNED is called ...*/ 
    FREEMEM_ALIGNED(grid->dataFloat[i]); 
#else
      FREEMEM(grid->dataFloat[i]);
#endif
    }
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int G3D_FreePyramid( G3D_Grid grid )
{
  int i;
  if(grid)
  {
    G3D_Pyramid *pyra = &grid->pyramid;
    for(i=0;i<pyra->n_level;i++)
    {
      G3D_Grid grid2=pyra->gaussian[i];
      if(grid2!=grid) // special case for top-level
      {
	G3D_DeleteGrid(&grid2);
      }
      pyra->gaussian[i] = NULL;
      G3D_DeleteGrid(&pyra->laplacian[i]);
    }
    pyra->n_level = 0;
  }
  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void G3D_SetWorldSpaceInfo( G3D_Grid G, double *pos, double dx )
{
  int k;
  for(k=0;k<3;k++) G->ws_pos[k]=pos[k];
  G->ws_dx = dx;
  G->ws_dx_r = 1./dx;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int G3D_DeleteGrid( G3D_Grid *p_adm )
{
  int rcThis=0;
  G3D_Grid adm=p_adm?*p_adm:NULL;

  if(adm)
  {
    int i;

    G3D_FreePyramid(adm);

    for(i=0;i<G3D_MAX_FCHAN;i++) 
    {
      G3D_FreeFloatChannel(adm,i);
      //if(adm->dataFloat[i]) FREEMEM(adm->dataFloat[i]);
      if(adm->dataUshort[i]) FREEMEM(adm->dataUshort[i]);
    }
    FREEMEM(adm);
    *p_adm = NULL;
  }

//rcCatch: 
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
G3D_Grid G3D_CreateGrid0( int sx, int sy, int sz,
    			 unsigned opt,
			 char *memOwnerId )
{
  int rcThis=0,k;
  G3D_Grid adm=NULL;

  ALLOCOBJ0_TAG(adm,memOwnerId);

  adm->sx = sx;
  adm->sy = sy;
  adm->sz = sz;
  adm->sxy = (LongInt)sx*(LongInt)sy;
#if 0
  adm->sx_plus_1 = adm->sx+1;
  adm->sxy_plus_1 = adm->sxy+1;
  adm->
#endif
  adm->n = (LongInt)sx*(LongInt)sy*(LongInt)sz;
  adm->smax = MAX(MAX(sx,sy),sz);
  adm->smin = MIN(MIN(sx,sy),sz);
  adm->dim[0] = sx;
  adm->dim[1] = sy;
  adm->dim[2] = sz;

  for(k=0;k<3;k++) adm->soff[k] = 0;

  adm->pyramid.n_level=0;

//rcCatch:
  if(rcThis)
  {
    G3D_DeleteGrid(&adm);
  }
  return adm;
}




