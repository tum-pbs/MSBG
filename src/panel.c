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
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "globdef.h"
#include "util.h"
#include "mtool.h"
#include "plot.h"
#include "bitmap.h"
#include "MGWin.h"
#include "panel.h"

MSBG_NAMESPACE_BEGIN

#define MAX_VERT_PER_POLY 10000

#define SNAP_TO_EDGE (1<<0)
#define HIGHLIGHT (1<<1)
#define DISPLAY_INFO (1<<2)

#define PNL_VISU_SLICES_MAXCHAN 3


PnlPanel *globActPanel = NULL;

struct 
{
  PnlPanel	*pnlInput,
  		*pnlOutput,
		*pnlTool,
		*pnlTool2,
		*pnlCalib,
		*pnlAux[100];
  int		 nPnlAux,nPnlAuxMax;
} 
globGwx;


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static void PnlVisualize3DFieldSlices1( 
		PnlPanel *pnl,		// Output panel
    		float *field, 		// 3D field
		float *slabZY, float *slabXZ, float *slabXY,  // optional slabs
		LongInt sx, LongInt sy, LongInt sz,  // resolution
		double x0, double y0, double z0,	// slab positions 
		MtStat *gstat,		// optional: global statistics
		unsigned options,
		BmpBitmap **pOutBmp
		)
{
  BmpBitmap *B = BmpNewBitmap( sx+sz, sy+sz, BMP_GEN);
  UT_ASSERT_FATAL(B);
  //UT_ASSERT_FATAL(!(options&PNL_PROVIDE_SLABSTATS));
  int doRGB = options & PNL_RGB;
  if(doRGB)
    BmpGetRGBChannel(B,BMP_CLEAR);
  else  
    BmpGetGreyChannel(B,BMP_CLEAR);
  MtStat gstat0;
  char chbuf[1024];
  LongInt nVoxels = sx*sy*sz;

  if(!field)
  {
    UT_ASSERT_FATAL(slabXY&&slabZY&&slabXZ);
  }

  if(!gstat && field)
  {
    gstat = &gstat0;
    UtTimer tm;
    TIMER_START(&tm);

    MtGetStat_f(nVoxels,field,gstat);

    TIMER_STOP(&tm);
    TRC3(("PnlVisualize3DFieldSlices/MtGetStat:  CPU = %.4f sec (%.0f voxels/sec)\n", 
	  (double)(TIMER_DIFF_MS(&tm)/1000.0), 
	  (sx*sy*sz)/(double)(TIMER_DIFF_MS(&tm)/1000.0)));
    TRC(("PnlVisualize3DField('%s') stats: %s\n", pnl->title,
	    (MT_STAT_SPRINT(chbuf,gstat),chbuf)));
  }
  else gstat=NULL;


  int x,y,z,
      xb0,yb0;
  float zmax,zmin;
  double zsum,zsum2;
  MtStat statXY,statZY,statXZ;

  // XY
  z = sz*z0;
  xb0 = 0; yb0 = 0;

  zmax=-1e20; zmin=1e20; zsum=0; zsum2=0;
  #pragma omp parallel
  {
    float my_zmax=-1e20, my_zmin=1e20;
    double my_zsum=0, my_zsum2=0;
    #pragma omp  for private(x,y)
    for(y=0; y<sy; y++)
    for(x=0; x<sx; x++)
    {
      LongInt ii=MT_GXYZ(sx,sx*sy, x,y,z),
	      xb = xb0+x, 
	      yb = yb0+sy-1-y,
	      jj=BMP_XYOFF(B,xb,yb);
      UT_ASSERT_FATAL(BMP_IN_RANGEB(xb,yb,B));
      if(doRGB)
      {
	B->data[jj] = field ? ((uint32_t*)field)[ii] : 
	  		      ((uint32_t*)slabXY)[MT_GXY(sx,x,y)];
      }
      else
      {
	double f = field ? field[ii] : slabXY[MT_GXY(sx,x,y)];
	B->dataGrey[jj] = f;
        my_zsum += f;
        my_zsum2 += f*f;
        if(f>my_zmax) my_zmax=f;
        if(f<my_zmin) my_zmin=f;
      }
    }
    #pragma omp critical
    {
      zmax=MAX(zmax,my_zmax);
      zmin=MIN(zmin,my_zmin);
      zsum+=my_zsum;
      zsum2+=my_zsum2;
    }
    MT_STAT_INIT(&statXY);
    MT_STAT_UPDATE_N(&statXY,nVoxels,zmin,zmax,zsum,zsum2);
    MT_STAT_RESULT(&statXY);
  }
#if 1
  // ZY
  x = sx*x0;
  xb0 = sx; yb0 = 0;

  zmax=-1e20; zmin=1e20; zsum=0; zsum2=0;
  #pragma omp parallel
  {
    float my_zmax=-1e20, my_zmin=1e20;
    double my_zsum=0, my_zsum2=0;
    #pragma omp  for private(y,z)
    for(y=0; y<sy; y++)
    for(z=0; z<sz; z++)
    {
      LongInt ii=MT_GXYZ(sx,sx*sy, x,y,z),
	      xb = xb0+z, 
	      yb = yb0+sy-1-y,
	      jj=BMP_XYOFF(B,xb,yb);
      UT_ASSERT_FATAL(BMP_IN_RANGEB(xb,yb,B));
      if(doRGB)
      {
	B->data[jj] = field ? ((uint32_t*)field)[ii] :
	  		      ((uint32_t*)slabZY)[MT_GXY(sz,z,y)];
      }
      else
      {            
	double f = field ? field[ii] : slabZY[MT_GXY(sz,z,y)];
	B->dataGrey[jj] = f;
	my_zsum += f;
	my_zsum2 += f*f;
	if(f>my_zmax) my_zmax=f;
	if(f<my_zmin) my_zmin=f;
      }
    }
    #pragma omp critical
    {
      zmax=MAX(zmax,my_zmax);
      zmin=MIN(zmin,my_zmin);
      zsum+=my_zsum;
      zsum2+=my_zsum2;
    }
    MT_STAT_INIT(&statZY);
    MT_STAT_UPDATE_N(&statZY,nVoxels,zmin,zmax,zsum,zsum2);
    MT_STAT_RESULT(&statZY);
  }

  // XZ
  y = sy*y0;
  xb0 = 0; yb0 = sy;
  
  zmax=-1e20; zmin=1e20; zsum=0; zsum2=0;
  #pragma omp parallel
  {
    float my_zmax=-1e20, my_zmin=1e20;
    double my_zsum=0, my_zsum2=0;
    #pragma omp  for private(x,z)
    for(z=0; z<sz; z++)
    for(x=0; x<sx; x++)
    {
      LongInt ii=MT_GXYZ(sx,sx*sy, x,y,z),
	      xb = xb0+x, 
	      yb = yb0+z,
	      jj=BMP_XYOFF(B,xb,yb);
      UT_ASSERT_FATAL(BMP_IN_RANGEB(xb,yb,B));
      if(doRGB)      
      {
	B->data[jj] = field ? ((uint32_t*)field)[ii] :
	  		      ((uint32_t*)slabXZ)[MT_GXY(sx,x,z)];
      }
      else
      {
        double f = field ? field[ii] : slabXZ[MT_GXY(sx,x,z)];
	B->dataGrey[jj] = f;
	my_zsum += f;
	my_zsum2 += f*f;
	if(f>my_zmax) my_zmax=f;
	if(f<my_zmin) my_zmin=f;
      }
    }
    #pragma omp critical
    {
      zmax=MAX(zmax,my_zmax);
      zmin=MIN(zmin,my_zmin);
      zsum+=my_zsum;
      zsum2+=my_zsum2;
    }
    MT_STAT_INIT(&statXZ);
    MT_STAT_UPDATE_N(&statXZ,nVoxels,zmin,zmax,zsum,zsum2);
    MT_STAT_RESULT(&statXZ);
  }
#endif

  if(doRGB)
  {
    PnlShowBitmap2(pnl,B);
  }
  else
  {
    PnlShowBitmap3(pnl,B,0,BMP_GREY,0,NULL,0);

  if(!(options&PNL_NOTEXT))
  {
    float wx0,wy0;

  {
    MtStat *s=&statXY;
    sprintf(chbuf,"XY,%.1f: %g - %g",z0,s->min,s->max);
    GWsettxtl(60,1.0,1,GWkrgb(255,255,255),GWkrgb(0,0,0),(char*)""); 
#ifdef TEXTPOS_LOWER
    PnlBmpCoordsInv2(pnl, 0,sy-1, &wx0,&wy0);
#else
    PnlBmpCoordsInv2(pnl, 0, 0+10, &wx0,&wy0);
#endif
    GWputtxtl(wx0,wy0,chbuf);
  }

  {
    MtStat *s=&statZY;
    sprintf(chbuf,"ZY,%.1f: %g - %g",x0,s->min,s->max);
    GWsettxtl(60,1.0,1,GWkrgb(255,255,255),GWkrgb(0,0,0),(char*)""); 
#ifdef TEXTPOS_LOWER
    PnlBmpCoordsInv2(pnl, sx-1, sy-1, &wx0,&wy0);
#else
    PnlBmpCoordsInv2(pnl, sx-1, 0+10, &wx0,&wy0);
#endif
    GWputtxtl(wx0,wy0,chbuf);
  }

  {
    MtStat *s=&statZY;
    sprintf(chbuf,"XZ,%.1f,Z: %g - %g",y0,s->min,s->max);
    GWsettxtl(60,1.0,1,GWkrgb(255,255,255),GWkrgb(0,0,0),(char*)""); 
#ifdef TEXTPOS_LOWER
    PnlBmpCoordsInv2(pnl, 0,sy+sz-1, &wx0,&wy0);
#else
    PnlBmpCoordsInv2(pnl, 0,sy-1+10, &wx0,&wy0);
#endif
    GWputtxtl(wx0,wy0,chbuf);
  }

  if(gstat)
  {
    sprintf(chbuf,"Global: %g - %g (%g %g)",
	gstat->min,gstat->max,gstat->avg,gstat->var);
    GWsettxtl(60,1.0,1,GWkrgb(255,255,255),GWkrgb(0,0,0),(char*)""); 
    PnlBmpCoordsInv2(pnl, sx-1,sy+sz-1, &wx0,&wy0);
    GWputtxtl(wx0,wy0,chbuf);
  }

  }
  }

  if(pOutBmp)
  {
    *pOutBmp = B;
  }
  else
  {
    BmpDeleteBitmap(&B);
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void PnlVisualize3DFieldSlices( 
		PnlPanel *pnl,		// Output panel
    		float **field, 		// 3D field (nChan channels)
		float **slabZY, float **slabXZ, float **slabXY,  // optional slabs
		int nChan,		// Number of data channels 
		LongInt sx, LongInt sy, LongInt sz,  // resolution
		
		double x0, double y0, double z0,	// slab positions 

		// Pointer to stats array:  global, slabZY0, slabXZ0, slabXY0,  //channel 0 
		//			    global, slabZY1, slabXZ1, slabXY1,  //channel 1 ...						    
		MtStat *stats,	
		unsigned options,
		BmpBitmap **pOutBmp
		)
{
  if(nChan==1&&(!(options&PNL_PROVIDE_SLABSTATS)))
    return PnlVisualize3DFieldSlices1(
			pnl,field?field[0]:NULL,
			slabZY?slabZY[0]:NULL, slabXZ?slabXZ[0]:NULL, slabXY?slabXY[0]:NULL,
			sx,sy,sz, x0,y0,z0, stats, options,
			pOutBmp );
  int k,idir;
  UtTimer tm;

  MtStat statGlob_[PNL_VISU_SLICES_MAXCHAN],
    	 statXY_[PNL_VISU_SLICES_MAXCHAN],
	 statZY_[PNL_VISU_SLICES_MAXCHAN],
	 statXZ_[PNL_VISU_SLICES_MAXCHAN];
  MtStat *statGlob[PNL_VISU_SLICES_MAXCHAN]={NULL},
    	 *statXY[PNL_VISU_SLICES_MAXCHAN]={NULL},
	 *statZY[PNL_VISU_SLICES_MAXCHAN]={NULL},
	 *statXZ[PNL_VISU_SLICES_MAXCHAN]={NULL};

  MtStat *statSlab;

  UT_ASSERT0(nChan>=1&&nChan<=3);

  TIMER_START(&tm);
  int db=2;
  BmpBitmap *B = BmpNewBitmap( 2*db + 3*MAX(sx,sz), 
			       2*db + 2*sy+sz,
			       BMP_GEN);
  UT_ASSERT_FATAL(B);

  int doRGB = options & PNL_RGB;
  if(doRGB)
    BmpGetRGBChannel(B,BMP_CLEAR);
  else  
    BmpGetGreyChannel(B,BMP_CLEAR);

  char chbuf[1024];

  if(!field)
  {
    UT_ASSERT_FATAL(slabXY&&slabZY&&slabXZ);
  }

#if 0
  if(!stats && field)
  {
    TIMER_START(&tm);
    for(k=0;k<nChan;k++)
    {
      statGlob[k] = &statGlob_[k];
      MtGetStat_f(nVoxels,field[k],statGlob[k]);
      TRC(("PnlVisualize3DField('%s') stats: %s\n", pnl->title,
	    (MT_STAT_SPRINT(chbuf,statGlob[k]),chbuf)));
    }
    TIMER_STOP(&tm);
    TRC3(("PnlVisualize3DFieldSlices/MtGetStat:  CPU = %.4f sec (%.0f voxels/sec)\n", 
	  (double)(TIMER_DIFF_MS(&tm)/1000.0), 
	  (sx*sy*sz)/(double)(TIMER_DIFF_MS(&tm)/1000.0)));
  }
#endif

  for(k=0;k<nChan;k++) 
  {
    if(stats) statGlob_[k] = stats[4*k+0]; 
    statGlob[k] = &statGlob_[k];
    if(stats) statZY_[k] = stats[4*k+1]; 
    statZY[k] = &statZY_[k];
    if(stats) statXZ_[k] = stats[4*k+2]; 
    statXZ[k] = &statXZ_[k];
    if(stats) statXY_[k] = stats[4*k+3]; 
    statXY[k] = &statXY_[k];
  }

  for(idir=0;idir<3;idir++)
  {
    int xb0=0,yb0=0;
    float zmax[3],zmin[3];
    double zsum[3],zsum2[3];
    float **slab;
    LongInt nSamples;

    int i0,i1,i2,s1,s2;
    switch(idir)
    {
      case 0: // X - ZY
	statSlab=&statZY_[0]; slab=slabZY; 
	i0=sx*x0; s1=sz; s2=sy; yb0=sz+db+sy+db; break;		      
      case 1: // Y - XZ
	statSlab=&statXZ_[0]; slab=slabXZ;
	i0=sy*y0; s1=sx; s2=sz; yb0=0; break;	
      case 2: // Z - XY
	statSlab=&statXY_[0]; slab=slabXY;
	i0=sz*z0; s1=sx; s2=sy; yb0=sz+db;break;	
      default: UT_ASSERT0(FALSE); break;
    }

    for(k=0;k<ARRAY_LENGTH(zmax);k++)
    {
      zmax[k]=-1e20; zmin[k]=1e20; zsum[k]=0; zsum2[k]=0;
    }
    nSamples=0;
    #pragma omp parallel
    {
      int k;
      float my_zmax[3], my_zmin[3];
      double my_zsum[3], my_zsum2[3];
      for(k=0;k<nChan;k++)
      {
	my_zmax[k]=-1e20; my_zmin[k]=1e20; my_zsum[k]=0; my_zsum2[k]=0;
      }    
      #pragma omp  for private(i1,i2) reduction(+:nSamples)
      for(i2=0; i2<s2; i2++)
      for(i1=0; i1<s1; i1++)
      {
	LongInt ii=0;
	int yb=0;
	switch(idir)
	{
	  case 0: 
	    ii = MT_GXYZ(sx,sx*sy, i0,i2,i1); 
	    yb = yb0+s2-1-i2;
	    break;
	  case 1: 
	    ii = MT_GXYZ(sx,sx*sy, i1,i0,i2); 
	    yb = yb0+i2;
	    break;
	  case 2: 
	    ii = MT_GXYZ(sx,sx*sy, i1,i2,i0); 
	    yb = yb0+s2-1-i2;
	    break;
	  default: break;	  
	}
	for(k=0;k<nChan;k++)
	{
	  int xb = xb0+i1+k*(s1+db);
	  LongInt jj = BMP_XYOFF(B,xb,yb);
	  UT_ASSERT_FATAL(BMP_IN_RANGEB(xb,yb,B));
	  double f;
	  if(doRGB)
	  {
	    f = field ? ((uint32_t*)field[k])[ii] : 
	  		      ((uint32_t*)slab[k])[MT_GXY(s1,i1,i2)];
	    /*if(xb>=57&&xb<=59 && yb>=130&&yb<=140)
	    {
	      f = 255;
	    }*/
	    B->data[jj] = f;
	  }
	  else
	  {
	    f = field ? field[k][ii] : slab[k][MT_GXY(s1,i1,i2)];
	    B->dataGrey[jj] = f;
	  }
	  my_zsum[k] += f;
	  my_zsum2[k] += f*f;
	  if(f>my_zmax[k]) my_zmax[k]=f;
	  if(f<my_zmin[k]) my_zmin[k]=f;
	}
	nSamples++;
      } // omp for
      #pragma omp critical
      {
	for(k=0;k<nChan;k++)
	{
	  zmax[k]=MAX(zmax[k],my_zmax[k]);
	  zmin[k]=MIN(zmin[k],my_zmin[k]);
	  zsum[k]+=my_zsum[k];
	  zsum2[k]+=my_zsum2[k];
	}
      }
    }  // omp parallel
    if(!stats)
    {
      for(k=0;k<nChan;k++)
      {
	MT_STAT_INIT(&statSlab[k]);
	MT_STAT_UPDATE_N(&statSlab[k],nSamples,zmin[k],zmax[k],zsum[k],zsum2[k]);
	MT_STAT_RESULT(&statSlab[k]);
      }
    }
  }

  PnlClear(pnl);
  if(doRGB)
  {
    PnlShowBitmap2(pnl,B);
  }
  else
  {
    PnlShowBitmap3(pnl,B,0,BMP_GREY,0,NULL,0);
  }
  if(!(options & PNL_NOTEXT))
  {
    float wx0,wy0;
    
    for(idir=0;idir<3;idir++)
    {
      for(k=0;k<nChan;k++)
      {
	MtStat *s=NULL;
	int xb,yb;
	switch(idir)
	{
	  case 0: // X - ZY
	    s=statZY[k]; sprintf(chbuf,"ZY,%.1f: %g - %g",x0,s->min,s->max);
	    xb = k*(sz+db);
#ifdef TEXTPOS_LOWER
	    yb = B->sy-1;
#else
	    yb = B->sy-1-sy+10;
#endif
	    break;
	  case 1: // Y - XZ
	    s=statXZ[k]; sprintf(chbuf,"XZ,%.1f: %g - %g",y0,s->min,s->max);
	    xb = k*(sx+db);
#ifdef TEXTPOS_LOWER
	    yb = sz-1;
#else
	    yb = 0+10;
#endif
	    break;
	  case 2: // Z - XY
	    s=statXY[k]; sprintf(chbuf,"XY,%.1f: %g - %g",z0,s->min,s->max);
	    xb = k*(sx+db);;
#ifdef TEXTPOS_LOWER
	    yb = sz+sy-1;
#else
	    yb = sz-1+10;
#endif
	    break;
	  default: UT_ASSERT0(FALSE); break;
	}
	GWsettxtl(50,1.0,1,GWkrgb(255,255,255),GWkrgb(0,0,0),(char*)""); 
	PnlBmpCoordsInv2(pnl, xb,yb, &wx0,&wy0);
	GWputtxtl(wx0,wy0,chbuf);
	char *tmp=MM_strdup(pnl->title);
	UT_ASSERT0(tmp);
	MiStripTrailingCRLF(tmp,0);
	TRC(("V3DFS(t=%g): '%s' -> %s , avg=%g\n",pnl->totalTime,tmp,chbuf,s->avg)); 
	MM_free(tmp);
      }
    }
#if 0
    if(gstat)
    {
      int k;
      for(k=0;k<nChan;k++)
      {
	sprintf(chbuf,"Global: %g - %g (%g %g)",
	    gstat[k].min,gstat[k].max,gstat[k].avg,gstat[k].var);
	GWsettxtl(60,1.0,1,GWkrgb(255,255,255),GWkrgb(0,0,0),(char*)""); 
	PnlBmpCoordsInv2(pnl, nChan*(sx+db)-1,db+7*k, &wx0,&wy0);
	GWputtxtl(wx0,wy0,chbuf);
      }
    }
#endif
  }

  if(pOutBmp)
  {
    *pOutBmp=B;
  }
  else
  {
    BmpDeleteBitmap(&B);
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
PnlPanel *PnlGetActPanel(void)
{
  return globActPanel;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void PnlDrawRect( PnlPanel *pnl, int x0, int y0, int rx, int ry,
    		    int cr, int cg, int cb, int mode )
{
    float fx1,fx2,fy1,fy2;

    PnlSelect( pnl, 0 );
    PnlBmpCoordsInv(pnl,x0,y0,&fx1,&fy1);
    PnlBmpCoordsInv(pnl,x0+rx-1, y0+ry-1,&fx2,&fy2);
    GWsetpen(-1,1,0,7);
    GWrect(fx1,fy1,fx2,fy2);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void PnlDrawPixel( PnlPanel *pnl, int x0, int y0,
    		    int cr, int cg, int cb, int mode )
{
    float fx1,fx2,fy1,fy2;

    PnlSelect( pnl, 0 );
    PnlBmpCoordsInv(pnl,x0,y0,&fx1,&fy1);
    PnlBmpCoordsInv(pnl,x0+1, y0+1,&fx2,&fy2);
    GWsetpen(GWkrgb(cr,cg,cb),1,0,4);
//    GWsetpen(-1,1,1,4);
    GWrect(fx1,fy1,fx2,fy2);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void PnlClear( PnlPanel *pnl )
{
    PnlSelect( pnl, 0 );
    GWclear(-1);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void PnlDrawLine( PnlPanel *pnl, int x1, int y1, int x2, int y2,
    		    int cr, int cg, int cb, int mode )
{
    float fx1,fx2,fy1,fy2;

    PnlSelect( pnl, 0 );
    PnlBmpCoordsInv(pnl,x1,y1,&fx1,&fy1);
    PnlBmpCoordsInv(pnl,x2,y2,&fx2,&fy2);
    GWsetpen(-1,1,0,7);
    GWline(fx1,fy1,fx2,fy2);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
PnlPanel *PnlActPanel()
{
  return globActPanel;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void PnlSelect( PnlPanel *pnl, unsigned options )
{
  GWselect( pnl->win );

  if( options & PNL_MAXIMIZE )
  {
    GWshowwn( pnl->win, 5 );
  }
  if( options & PNL_RESTORE )
  {
    GWshowwn( pnl->win, 3 );
    GWarrange( 5 );
  }
  if( options & PNL_ACTIVATE )
  {
    GWshowwn( pnl->win, 4 );
  }
  globActPanel = pnl;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void PnlBmpCoords( PnlPanel *pnl, float wx, float wy,
    			  int *px, int *py )
{
  *px = wx/pnl->fx + pnl->bmp_roi.x0; 
  *py = (pnl->wsy-wy)/pnl->fy + pnl->bmp_roi.y0; 
}
void PnlBmpCoordsInv( PnlPanel *pnl, int px, int py,
    			     float *pwx, float *pwy )
{
  *pwx = (px-pnl->bmp_roi.x0)*pnl->fx;  
  *pwy = pnl->wsy - (py - pnl->bmp_roi.y0)*pnl->fy; 
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void PnlBmpCoords2( PnlPanel *pnl, float wx, float wy,
    			  int *px, int *py )
{
  *px = pnl->x02 + wx/pnl->fx2; 
  *py = pnl->y02 - wy/pnl->fy2; 
}
void PnlBmpCoordsInv2( PnlPanel *pnl, float px, float py,
    			     float *pwx, float *pwy )
{
  *pwx = (px-pnl->x02)*pnl->fx2;
  *pwy = (pnl->y02 - py)*pnl->fy2;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void PnlFree( PnlPanel *pnl )
{
  PnlFreeBitmap( pnl );
  if(pnl->bitmapTmp) 
  {
    BmpCloseBitmap( pnl->bitmapTmp );
    pnl->bitmapTmp = NULL;
  }
  if( pnl->bmpGw ) 
  {
    MGWdelbmp( pnl->bmpGw );
    pnl->bmpGw = 0;
  }
  if( pnl->bmpGw2 ) 
  {
    MGWdelbmp( pnl->bmpGw2 );
    pnl->bmpGw2 = 0;
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void PnlInit( PnlPanel *pnl, int id, int win)
{
  float fx1,fy1,fx0,fy0;
  pnl->id = id;
  pnl->win = win;
  pnl->undoCount = 0;
  pnl->bitmap = NULL;
  pnl->bitmap2 = NULL;
  pnl->bitmapSave = NULL;
  pnl->bitmapTmp = NULL;
  pnl->totalTime = -1;
  BMP_SET_RECT(&pnl->bmp_roi,0,0,0,0);
  pnl->bmpGw = 0;
  pnl->bmpGw2 = 0;
  pnl->bmpGw_rr = 1;
  pnl->bmpGw2_rr = 1;
  strcpy(pnl->bmpFileName,"default.bmp");
  GWselect( win );
  GWgetwn( &fx0, &fy0, &fx1, &fy1 );
  pnl->wx0=fx0; pnl->wy0=fy0;
  pnl->wx1=fx1; pnl->wy1=fy1;
  PlScreenInit( &pnl->plScreen, 0, 0, fx1-fx0, fy1-fy0 );
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
PnlPanel *PnlCreate2( char *title, int wsx, int wsy, unsigned opt )
{
  PnlPanel *pnl=NULL;
  int win,i,j, rcThis=0;
  float	    fx0,fy0,fx1,fy1;

  if(!wsx) wsx = 3800; 
  if(!wsy) wsy = 2300;

  win = GWopenxl(0,wsx,wsy,19,0,10,title);
  if(win==0) TRCERRR(("GWopenx() failed\n"),MI_ERROR);

  ALLOCOBJ0(pnl);

  strcpy(pnl->title,title?title:"");

  pnl->wsx=wsx; pnl->wsy=wsy;

  PnlInit(pnl, win, win );

  PnlSelect( pnl, 0 );

  GWshowwn( pnl->win, 10 );

  GWgetwn( &fx0, &fy0, &fx1, &fy1 );
  PlScreenInit( &pnl->plScreen, 0, 0, fx1-fx0, fy1-fy0 );

    for(i=0;i<2;i++) for(j=0;j<2;j++)
    {
      PlScreen *pls= &pnl->plScreen;
      PlWindowInit( &pnl->plwin[i*2+j], 
	  &pnl->plScreen, j*50, (1-i)*50, 
	  45*(float)pls->sy/(float)pls->sx,45 );
      PlPlotInit(&pnl->plplot[i*2+j],
	  	 &pnl->plwin[i*2+j],0,1,0,1,0);
    }
    PlWindowInit( &pnl->plwin0, &pnl->plScreen, 
	0,0,100,100);
    PlPlotInit(&pnl->plplot0,
	  	 &pnl->plwin0,0,1,0,1,0);

      PlWindowInit( &pnl->plwinIso, 
	  &pnl->plScreen, 0, 0, 
	  100*(float)(pnl->plScreen).sy/
	    (float)(pnl->plScreen).sx,100 );

  PlColor( GWkrgb( 70,70,70 ));
  GWsetpen(-1,1,0,4);


rcCatch:
  return pnl;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
PnlPanel *PnlCreate( char *title )
{
  return PnlCreate2(title,0,0,0);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int PnlDelete( PnlPanel *pnl )
{
  if(pnl)
  {
    GWclose(pnl->win);
    PnlFree(pnl);
    FREEMEM(pnl);
  }
  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void PnlSetTmpBitmap( PnlPanel *pnl, BmpBitmap *bmp )
{
  if( pnl->bitmapTmp )
  {
    BmpCloseBitmap( pnl->bitmapTmp );
    pnl->bitmapTmp = NULL;
  }
  pnl->bitmapTmp = bmp;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void PnlFreeBitmap( PnlPanel *pnl )
{
  int rc;

  if( pnl->bitmap )
  {
    GWselect(pnl->win);
    GWclear(-1 );
    rc = BmpCloseBitmap( pnl->bitmap );
    if( rc ) TRCERR(("BmpCloseBitmap() failed %d\n"));
    pnl->bitmap = NULL;
  } 
  if( pnl->bitmapSave )
  {
    rc = BmpCloseBitmap( pnl->bitmapSave );
    pnl->bitmapSave = NULL;
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void PnlSetBitmap( PnlPanel *pnl, BmpBitmap *bmp, BmpBitmap *bmpSave,
    		       char *fname )
{
  int	rc;
  float f,wsx,wsy,x1,y1,x2,y2;
 
  PnlFreeBitmap( pnl );

  if( bmp ) pnl->bitmap = bmp;
  pnl->bitmapSave = bmpSave;

  if( bmp )
  {
    BmpInitBitmapWindow( &pnl->bmw0, bmp,
	0,0,bmp ? bmp->sx:0,bmp ? bmp->sy:0 );

    rc = GWselect(pnl->win);

    GWgetwn( &x1, &y1, &x2, &y2 );
    wsx = x2 - x1; wsy = y2 - y1;
    pnl->wsx=wsx; pnl->wsy=wsy;

    if( bmp )
    {
#if 0    
      if( bmp->sx > bmp->sy )
      {
	f = (float)wsx/bmp->sx;
      }
      else
      {
	f = (float)wsy/bmp->sy;
      }
      /*    f *= 2;*/
#endif

      f = MIN( (float)wsx/(float)bmp->sx, (float)wsy/(float)bmp->sy );

      TRC(("win:%fx%f bmp:%dx%d f=%f\n",wsx,wsy,bmp->sx,bmp->sy));
      pnl->fx = f; pnl->x0 = 0;
      pnl->fy = f; pnl->y0 = bmp->sy;

      PnlUpdateBitmapRoi(pnl,NULL);
    }

    if( fname )
    {
      strncpy( pnl->bmpFileName, fname, MI_MAX_FNAME_LEN );
    }
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void PnlUpdateBitmapRoi( PnlPanel *pnl, BmpRectangle *roi )
{
  int	rc;
  float f,wsx,wsy,x1,y1,x2,y2;
  BmpRectangle roi0;
  if(!roi)
  {
    BMP_SET_RECT(&roi0,0,0,pnl->bitmap->sx,pnl->bitmap->sy);
    roi=&roi0;
  }
  rc = GWselect(pnl->win);
  pnl->bmp_roi = *roi;
  GWgetwn( &x1, &y1, &x2, &y2 );
  wsx = x2 - x1; wsy = y2 - y1;
  pnl->wsx=wsx; pnl->wsy=wsy; 
  f = MIN( (float)wsx/(float)roi->sx, (float)wsy/(float)roi->sy );
  pnl->fx = f;
  pnl->fy = f; 
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void PnlErase( PnlPanel *pnl )
{
  GWselect(pnl->win);
  GWerase(0,0);
}

/*-------------------------------------------------------------------------*/   
/* 									   */
/*-------------------------------------------------------------------------*/
static int PnlIDisplayBmp( PnlPanel *pnl, BmpBitmap *bmp0, int force )
{
  BmpBitmap *bmp = NULL;
  int rc,rcThis=0;
  float fx,fy;
  BmpRectangle *roi=&pnl->bmp_roi;

  if(!bmp0) 
  {
    bmp = pnl->bitmap;
    if(roi->sx*roi->sy<=0)
    {
      if(bmp) 
      {
	BMP_SET_RECT(roi,0,0,bmp->sx,bmp->sy);
      }
    }
  }
  else
  {
    bmp = bmp0;
    if(roi->sx*roi->sy<=0)
    {
      BMP_SET_RECT(roi,0,0,bmp->sx,bmp->sy);
    }
    bmp->GWnr = 0;
    rc = GWselect(pnl->win);
    bmp->GWnr = MGWmakebmp( bmp, 0, bmp->sx, bmp->sy, 24, 255.0, BMP_GEN, 0, 0, 0,
			    roi, 
			    &pnl->bmpGw_rr );
    if(!bmp->GWnr )
    TRCERRR(("MGWmakebmp(%d x %d) failed\n",bmp->sx,bmp->sy),MI_ERROR);
    GWsetbmp(bmp->GWnr,0,0,1,0,1);
  }

  if(!bmp) raiseRc(0);

  rc = GWselect(pnl->win);

  if( bmp->doRegenerate || force )
  {
    GWerase(0,0);
    rc = MGWmakebmp( bmp, bmp->GWnr, bmp->sx, bmp->sy, 24, 255.0,
		     BMP_GEN, 0, 0, 0, roi, &pnl->bmpGw_rr );
    if((!bmp0 &&!bmp->GWnr) || rc != bmp->GWnr )
      TRCERRR(("MGWmakebmp(%d x %d) failed\n",bmp->sx,bmp->sy),MI_ERROR);
    bmp->doRegenerate = FALSE;
  }

  if( bmp0 )
  {
    float f,wsx,wsy,x1,y1,x2,y2;
    GWgetwn( &x1, &y1, &x2, &y2 );
    wsx = x2 - x1; wsy = y2 - y1;
    pnl->wsx=wsx; pnl->wsy=wsy;
    f = MIN( (float)wsx/(float)roi->sx, (float)wsy/(float)roi->sy );
    fx = fy = f;
  }
  else
  {
    fx = pnl->fx; fy = pnl->fy;
  }

  rc = GWsetbmp(bmp->GWnr,roi->sx*fx,roi->sy*fy,1,0,1);
  rc = GWputbmp(bmp->GWnr, 0, 0, 1 );

  TRC3(("PnlDisplayBmp pnl=%d bmp='%s'\n",
	pnl->id,bmp->name?bmp->name:"-"));
rcCatch:  
  return rcThis;
}


/*-------------------------------------------------------------------------*/   
/* 									   */
/*-------------------------------------------------------------------------*/
int PnlDisplayBmp( PnlPanel *pnl, int force )
{
  return PnlIDisplayBmp( pnl, NULL, force );
}


/*-------------------------------------------------------------------------*/   
/* 									   */
/*-------------------------------------------------------------------------*/
int PnlDisplayBmp0( PnlPanel *pnl, BmpBitmap *bmp, int force )
{
  return PnlIDisplayBmp( pnl, bmp, force );
}

/*-------------------------------------------------------------------------*/   
/* 									   */
/*-------------------------------------------------------------------------*/
int PnlIDisplayBmp2( PnlPanel *pnl, BmpBitmap *bmp, float x0, float y0,
    		    float wsxp, float wsyp,
    		    int force )
{
  float f,wsx,wsy,x1,y1,x2,y2,bx0=x0,by0=y0,fx,fy;
  int gwBmp,rc,rcThis=0;

  GWselect( pnl->win );
  GWgetwn( &x1, &y1, &x2, &y2 );
  wsx = x2 - x1; wsy = y2 - y1;
  f = MIN( (float)wsx/(float)bmp->sx, (float)wsy/(float)bmp->sy );
  if(fabs(wsxp)>MT_REAL_EPSILON || fabs(wsyp)>MT_REAL_EPSILON) 
    f *= MAX(wsxp,wsyp)/100.;
  fx=fy=f;

  GWerase(0,0);

    gwBmp = GWmakebmp( 0, bmp->sx, bmp->sy, 24, bmp->data );
    if(!gwBmp )
      TRCERRR(("GWmakebmp(%d x %d) failed\n",bmp->sx,bmp->sy),MI_ERROR);


    printf("=====> gwBmp=%d\n",gwBmp);

    rc = GWsetbmp(gwBmp,bmp->sx*fx,bmp->sy*fy,1,0,1);
    if( rc != gwBmp ) TRCERRR(("GWsetbmp() failed %d\n",rc),MI_ERROR);
    rc = GWputbmp(gwBmp, bx0, by0, 1 );
rcCatch:
    return rcThis;
}

/*-------------------------------------------------------------------------*/   
/* 									   */
/*-------------------------------------------------------------------------*/
int PnlDisplayBmp2( PnlPanel *pnl, BmpBitmap *bmp, float x0, float y0,
    		    float wsx, float wsy, /* in percent */
    		    int force )
{
  return PnlIDisplayBmp2(pnl,bmp,x0,y0,wsx,wsy,force);    
}

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int PnlShowBitmap( PnlPanel *pnl, BmpBitmap *bmp )
{
  float f,wsx,wsy,x1,y1,x2,y2,bx0=0,by0=0,fx,fy;
  int gwBmp,rc,rcThis=0;

/*  GWselect(pnl->win);
  GWerase(0,0);*/
  if( pnl->bmpGw ) 
  {
    MGWdelbmp( pnl->bmpGw );
    pnl->bmpGw = 0;
  }
  GWselect(pnl->win);
  GWerase(0,0);

  GWgetwn( &x1, &y1, &x2, &y2 );
  wsx = x2 - x1; wsy = y2 - y1;
  f = MIN( (float)wsx/(float)bmp->sx, (float)wsy/(float)bmp->sy );
  fx=fy=f;
  pnl->fx = pnl->fy = f; pnl->x0 = 0; pnl->y0 = bmp->sy;

  gwBmp = MGWmakebmp( bmp, 0, bmp->sx, bmp->sy, 24, 255.0, BMP_GEN, 0, 0,
      		      0, NULL, &pnl->bmpGw_rr );
  if(!gwBmp )
    TRCERRR(("MGWmakebmp(%d x %d) failed\n",bmp->sx,bmp->sy),MI_ERROR);
  pnl->bmpGw = gwBmp;  
  rc = GWsetbmp(gwBmp,bmp->sx*fx,bmp->sy*fy,1,0,1);
  if( rc != gwBmp ) TRCERRR(("GWsetbmp() failed %d\n",rc),MI_ERROR);
  rc = GWputbmp(gwBmp, bx0, by0, 1 );
rcCatch:
    return rcThis;

}

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int PnlShowBitmap3( PnlPanel *pnl, BmpBitmap *bmp, 
    		    float gmax, //0: do norm
    		    int chan, int ichan,
		    BmpRectangle *roi,
		    unsigned opt )
{
  float f,wsx,wsy,x1,y1,x2,y2,bx0=0,by0=0,fx,fy;
  int gwBmp,rc,rcThis=0;
  BmpRectangle roi0;

  if(!roi)
  {
    roi=&roi0;
    BMP_SET_RECT(roi,0,0,bmp->sx,bmp->sy);
  }

  GWselect(pnl->win);
  GWerase(0,0);
  if( pnl->bmpGw2 ) 
  {
    MGWdelbmp( pnl->bmpGw2 );
    pnl->bmpGw2 = 0;
  }
  GWselect(pnl->win);
  GWerase(0,0);

  if(opt & PNL_CLEAR) GWclear(-1);

  pnl->bitmap2 = bmp;
  GWgetwn( &x1, &y1, &x2, &y2 );
  wsx = x2 - x1; wsy = y2 - y1;
  f = MIN( (float)wsx/(float)roi->sx, (float)wsy/(float)roi->sy );
  fx=fy=f;
  pnl->fx2 = pnl->fy2 = f; pnl->x02 = 0; pnl->y02 = roi->sy;

  gwBmp = MGWmakebmp( bmp, 0, bmp->sx, bmp->sy, 24, gmax, chan, ichan,
      		      gmax<MT_NUM_EPS ? TRUE:FALSE, opt & PNL_IGNAN, roi,
      		      &pnl->bmpGw2_rr );
  if(!gwBmp )
  {
    pnl->bmpGw2 = 0;
    TRCERRR(("MGWmakebmp(%d x %d) failed\n",roi->sx,roi->sy),MI_ERROR);
  }
  pnl->bmpGw2 = gwBmp;  
  rc = GWsetbmp(gwBmp,roi->sx*fx,roi->sy*fy,1,0,1);
  if( rc != gwBmp ) TRCERRR(("GWsetbmp() failed %d\n",rc),MI_ERROR);
  rc = GWputbmp(gwBmp, bx0, by0, 1 );

rcCatch:
    return rcThis;

}

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int PnlShowBitmap2( PnlPanel *pnl, BmpBitmap *bmp )
{
  return PnlShowBitmap3(pnl,bmp,255.0,BMP_GEN,0,NULL,0);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int PnlLoadBitmap( char *fnameIn0, PnlPanel *pnl, BmpBitmap *bmpIn,
    		   CoConverter *colConv,
    			    unsigned options, int resizePerc  )
{
  int		rcThis=0,rc,optGreyOnly=options&PNL_GREY_ONLY;
  unsigned	bmpOptClipboard=0;
  BmpBitmap     *bmp = NULL,*bmpSave=NULL;
  char		fname[MI_MAX_FNAME_LEN+1],*fnameIn = fnameIn0;

  if(!bmpIn)
  {
    if( options & PNL_CLIPBOARD ) 
    {
      if(!fnameIn) fnameIn = "_clipboard_";   
    }
    else
    {
      if(!fnameIn) fnameIn = pnl->bmpFileName;   

      if(strcmp(UtFileSuffix(fnameIn),"")==0)
	UtSetFilenameSuffix( fname, fnameIn, ".png" ); 
      else
	strcpy(fname,fnameIn);
    }
  }

  /*-----------------------------------------------------------------------*/
  /* free existing bitmap						   */
  /*-----------------------------------------------------------------------*/
  PnlFreeBitmap( pnl );

  /*-----------------------------------------------------------------------*/
  /* load new bitmap							   */
  /*-----------------------------------------------------------------------*/
  if(!bmpIn)
  {
    int win = pnl->win;
    bmp = BmpCreateBitmap2( fname, NULL, 0, 0, 0, fname,  colConv,
			     BMP_RGB | bmpOptClipboard |
			     (optGreyOnly ? BMP_GREY_ONLY : 0));
    if( resizePerc && bmp )
    {
	UT_ASSERT(FALSE); // deprecated
    }

    GWselect(win);
    bmp->GWnr = MGWmakebmp( bmp, 0, bmp->sx, bmp->sy, 24, 255.0, BMP_GEN,0,0,0,
			    NULL, &pnl->bmpGw_rr );
    if(!bmp->GWnr )
    TRCERRR(("MGWmakebmp(%d x %d) failed\n",bmp->sx,bmp->sy),MI_ERROR);
    rc = GWsetbmp(bmp->GWnr,0,0,1,0,1);

  }
  else 
  {
    bmp = bmpIn;
    
    if(!bmp->GWnr) 
    {
      bmp->GWnr = MGWmakebmp( bmp, 0, bmp->sx, bmp->sy, 24, 255.0,BMP_GEN,0,0,0,
	  		      NULL, &pnl->bmpGw_rr);
      if(!bmp->GWnr )
         TRCERRR(("MGWmakebmp(%d x %d) failed\n",bmp->sx,bmp->sy),MI_ERROR);
      rc = GWsetbmp(bmp->GWnr,0,0,1,0,1);
    }
  }

  if(!bmp) raiseRc(MI_ERROR);

  if( options & PNL_CREATE_BMP_SAVE )
  {
    bmpSave = BmpCreateBitmap( NULL, bmp, 0, 0, 0, NULL, 0 );
    if(!bmpSave) raiseRc(MI_ERROR);
  }

  PnlSetBitmap( pnl, bmp, bmpSave, NULL );

  PnlDisplayBmp( pnl, 0 );

  if( fnameIn0 && !bmpIn )
  {
    strcpy( pnl->bmpFileName, fname );
  }

rcCatch:
  if( rcThis ) 
  {
    if( bmp && ! bmpIn ) BmpCloseBitmap( bmp );
    if(bmpSave) BmpCloseBitmap( bmpSave );
  }
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int PnlSaveBitmap( char *fnameIn0, PnlPanel *pnl, CoConverter *cc )
{
  int 	rc,rcThis=0;
  char	*fnameIn = fnameIn0, 
  	 fname[MI_MAX_FNAME_LEN+1];

  if(! pnl->bitmap )
  {
    TRCERRR(("no bitmap\n"),MI_EINVARG);
  }

  if(!fnameIn) 
  {
    fnameIn = pnl->bmpFileName;   
  }
  
  strcpy(fname,fnameIn);
/*
  if(strcmp(UtFileSuffix(fname),"")==0)
  {
    setFilenameSuffix( fname, fname, ".bmp");
  }
*/  

  rc = BmpSaveBitmap( pnl->bitmap, pnl->win, fname, BMP_NULL, cc, 0 );
  if(rc) 
  {
    TRCERRR(("BmpSaveBitmap() failed %d",rc),rc);
  }

  if( fnameIn0 )
  {
    strcpy( pnl->bmpFileName, fname );
  }

/*  printf("image saved as '%s'\n",fname );*/

rcCatch:
  return rcThis;
}

/*====================================================================*/
/*							              */
/* 		Misc						      */
/*							              */
/*====================================================================*/

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/

#define MGW_DUMMY	9999999
#define MGW_MAXSZ	(3000*3000)
int MGWmakebmp(BmpBitmap *bmp, 
    	       int NB, int IW, int IH, int IBC, float gmax,
	       int chan, int ichan, int do_norm, int do_ignore_nan,
	       BmpRectangle *roi,
	       int *p_rr )
{
  BmpBitmap *bmp2=NULL;
  int	rc=0,x,y,x2,y2,c,g,rr=1,has_roi=roi?TRUE:FALSE,do_rgb_us=0;
  LongInt IW2,IH2,IW0=IW,IH0=IH;
  float sf= gmax > MT_NUM_EPS ? 255.0/gmax : 1;
  int   *IBITS=bmp->data;
  BmpRectangle roi0;
  //if(glob->noGwBmp) return MGW_DUMMY;
  if(chan == (BMP_RGB|BMP_USHORT)) 
  {
    chan=BMP_USHORT;
    do_rgb_us=TRUE;
  }

  if(!roi)
  {
    roi=&roi0;
    BMP_SET_RECT(roi,0,0,bmp->sx,bmp->sy);
  }
  IW = roi->sx; IH = roi->sy;
  if((LongInt)IW*(LongInt)IH>MGW_MAXSZ || (!IBITS) || (chan!=BMP_GEN)||has_roi)
  {
    MtStat st={0};
    float f,span=1;
    BmpData bmpData={0};
    float vRgb[3],vLuv[3];
    USE(vRgb);USE(vLuv);

#if 0    
    if((LongInt)IW*IH>MGW_MAXSZ) rr = 4;
    if((LongInt)IW*IH>8*MGW_MAXSZ) rr = 8;
#else
    for(rr=1; 
	( ((LongInt)(IW/rr)*(LongInt)(IH/rr)) >= MGW_MAXSZ )  && (rr<=1024); 
	rr*=2);
#endif
    bmpData.pFloat = NULL;
      
    IW2 = IW/rr; IH2 = IH/rr;
    
    TRC3(("MGWmakebmp: rr=%d,%dx%d -> %dx%d\n",rr,IW,IH,IW2,IH2));

    if(chan!=BMP_GEN)
    {
      rc = BmpGetDataChannel(bmp,chan,ichan,0,&bmpData);
      UT_ASSERT0(rc==0);
    }

    bmp2=BmpNewBitmap(IW2,IH2,BMP_RGB);
    UT_ASSERT0(bmp2);

    if(p_rr) *p_rr = rr;

    if(do_norm)
    {
      UT_ASSERT0(bmpData.pFloat);
      MT_STAT_INIT(&st);
      BMP_LOOP_XY(bmp2,x2,y2)
      {
        x=x2*rr+roi->x0;y=y2*rr+roi->y0;
        BMP_CLIP_CONST(0,0,IW0-1,IH0-1,x,y);
	f = BMP_DPIXEL(bmp,bmpData.pFloat,x,y);
	if(do_ignore_nan && isnan(f)) continue;
	MT_STAT_UPDATE(&st,f);
      }
      MT_STAT_RESULT(&st);
      span = st.max-st.min; if(span<MT_NUM_EPS) span=1;
    }

    #pragma omp for private(x2,y2,x,y,g,c,vLuv,vRgb)    
    BMP_LOOP_XY(bmp2,x2,y2)
    {
      x=x2*rr+roi->x0;y=y2*rr+roi->y0;
      BMP_CLIP_CONST(0,0,IW0-1,IH0-1,x,y);
      if(do_rgb_us)
      {
        LongInt ixy=BMP_XYOFF(bmp,x,y);
	int k,crgb[3];
	for(k=0;k<3;k++) crgb[k] = bmp->dataUshort[k][ixy]/256;
	c = BMP_MKRGB(crgb[0],crgb[1],crgb[2]);
      }
      else if(bmpData.pFloat)
      {
	float f = BMP_DPIXEL(bmp,bmpData.pFloat,x,y);
	if(do_ignore_nan && isnan(f))
	{
	  c = BMP_MKRGB(PLRGB(0x700),0,0);
	}
	else 
	{
	  if(do_norm)
	  {
	    g = 255.0*(f-st.min)/span;
	  }
	  else
	  {
	    g = BMP_DPIXEL(bmp,bmpData.pFloat,x,y)*sf;
	  }
	  if(g<0)g=0;
	  if(g>255)g=255;
	  c = BMP_MKRGB(g,g,g);
	}
      }
      else if(bmpData.pUshort)
      {
	g = (float)(BMP_DPIXEL(bmp,bmpData.pUshort,x,y))*sf;
	if(g<0)g=-g;
	if(g>255)g=255;
	c = BMP_MKRGB(g,g,g);
      }
      else if(bmpData.pUbyte)
      {
	g = (float)(BMP_DPIXEL(bmp,bmpData.pUbyte,x,y))*sf;
	if(g<0)g=-g;
	if(g>255)g=255;
	c = BMP_MKRGB(g,g,g);
      }
      else if(bmp->data)
      {
	c = BMP_PIXEL(bmp,data,x,y);
      }
      else if(bmp->dataUV&&bmp->dataGrey)
      {
	UT_ASSERT0(FALSE);
      }
      else if(bmp->dataGrey)
      {
	g = BMP_PIXEL(bmp,dataGrey,x,y)*sf;
	if(g<0)g=-g;
	if(g>255)g=255;
	c = BMP_MKRGB(g,g,g);
      }
      else 
      {
	UT_ASSERT0(FALSE);
      }
      BMP_PIXEL(bmp2,data,x2,y2) = c;      
    }
    IBITS = bmp2->data;
    IW=IW2;IH=IH2;
  }
  rc = GWmakebmp(NB,IW,IH,IBC,IBITS);
rcCatch:
  if(bmp2) BmpCloseBitmap(bmp2);
  return rc;
}
int MGWsetbmp(BmpBitmap *bmp,
    	      int NB, float W, float H, int MX, int ITR, int IOF)
{
  if(NB==MGW_DUMMY) return 1;
  return GWsetbmp( NB,W,H,MX,ITR,IOF);
}
int MGWputbmp(BmpBitmap *bmp,
    	      int NB, float X, float Y, int IBK)
{
  if(NB==MGW_DUMMY) return 1;
  return GWputbmp(NB,X,Y,IBK);
}
int MGWdelbmp(int NM)
{
  if(NM==MGW_DUMMY) return MGW_DUMMY;
  return GWdelbmp(NM);  
}

/*-------------------------------------------------------------------------*/
PnlPanel *MiGetAuxPanel2(int *pi, const char *title)
/*-------------------------------------------------------------------------*/
{

  int i,rcThis=0,nPnlAux_=globGwx.nPnlAux;
  char chbuf[1000];
  PnlPanel *pnl=NULL;
  i = pi ? *pi : -1;
  if(i==-1) 
  {
    i = globGwx.nPnlAux;
    nPnlAux_++;
  }
  if(i>ARRAY_LENGTH(globGwx.pnlAux))
  {
    TRCERRR(("panel id %d out of range\n",i),1);
  }
  if(i>globGwx.nPnlAuxMax) globGwx.nPnlAuxMax=i;
  globGwx.nPnlAux = nPnlAux_;
  pnl = globGwx.pnlAux[i];
  
  if(!pnl)
  {
    if(!title)
    {
      sprintf(chbuf,"AUX%d\n",i);
      title=chbuf;
    }
    pnl = PnlCreate((char*)title);
    if(!pnl) raiseRc(MI_ERROR);
    globGwx.pnlAux[i] = pnl;
  }
#if 1
  else
  {
    if(pi && *pi!=-1 && title && strcmp(pnl->title,title)!=0) 
    {
      int pi2=-1;
      PnlPanel *pnl2 = MiGetAuxPanel2( &pi2, title );
      *pi = pi2;
      return pnl2;
    }
  }
#endif

rcCatch:
  if(pi) *pi = i;
  return pnl;
}

/*-------------------------------------------------------------------------*/
PnlPanel *MiGetAuxPanel(int *pi)
/*-------------------------------------------------------------------------*/
{
  return MiGetAuxPanel2(pi,NULL);
}

MSBG_NAMESPACE_END

