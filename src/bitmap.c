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
#include <float.h>
#include <limits.h>
#include "stdmath.h"
#include "globdef.h"
#include "util.h"
#include "mtool.h"
#include "bitmap.h"
#include "MGWin.h"

MSBG_NAMESPACE_BEGIN

//#define BMP_WITH_GRWIN

#define BMP_FORMAT_ID_U48 	"U480"
#define BMP_U48_CCID_LEN	8

#define BMP_EYE_FREE 0xDEADBEAA
#define BMP_EYE_VALID 0xABCD1234

#define XYZRGB_GAMMA 0

static char minifont[] = 
\
"111110"
"100010"
"100010"
"100010"
"111110"
"000000"
\
"001000"
"001000"
"001000"
"001000"
"001000"
"000000"
\
"111110"
"000010"
"111110"
"100000"
"111110"
"000000"
\
"111100"
"000010"
"111100"
"000010"
"111100"
"000000"
\
"010000"
"010000"
"010100"
"011110"
"000100"
"000000"
\
"111110"
"100000"
"111110"
"000010"
"111110"
"000000"
\
"111110"
"100000"
"111110"
"100010"
"111110"
"000000"
\
"111110"
"000010"
"000100"
"001000"
"010000"
"000000"
\
"111110"
"100010"
"111110"
"100010"
"111110"
"000000"
\
"111110"
"100010"
"111110"
"000010"
"111110"
"000000"
\
"000000" /*'.'*/
"000000"
"000000"
"000000"
"001000"
"000000"
\
"000000" /*'+'*/
"001000"
"011100"
"001000"
"001000"
"000000"
\
"000000" /*'-'*/
"000000"
"011110"
"000000"
"000000"
"000000"
\
"111110" /*'E'*/
"100000"
"111110"
"100000"
"111111"
"000000"
\
"111110" /*'unknown'*/
"111110"
"111110"
"111110"
"111110"
"000000";

//
// global defines
//

#define GP_HDL_BORD_CLIP 1

//
// global data
//

#if 1
float bmp_dflt_pyra_kern_5x5[5*5] = 
{
  // standard near-gaussian convolution kernel (cf. Burt 1983)
  0.0025,    0.0125,    0.0200,    0.0125,    0.0025, 
  0.0125,    0.0625,    0.1000,    0.0625,    0.0125,
  0.0200,    0.1000,    0.1600,    0.1000,    0.0200, 
  0.0125,    0.0625,    0.1000,    0.0625,    0.0125, 
  0.0025,    0.0125,    0.0200,    0.0125,    0.0025 
};
#else

float bmp_dflt_pyra_kern_5x5[5*5] = 
{  // XXX
  // standard near-gaussian convolution kernel (cf. Burt 1983)
  0.0008,    0.004,    0.0100,    0.0125,    0.0025, 
  0.004,     0.0300,    0.1000,    0.0625,    0.0125,
  0.0100,    0.1000,    0.1600,    0.1000,    0.0200, 
  0.0125,    0.0625,    0.1000,    0.0625,    0.0125, 
  0.0025,    0.0125,    0.0200,    0.0125,    0.0025 
};

#endif

CoGlobalConst coGlobConst =
{
	1,   /* rgbClippingStrength */
	0.3,	/* grayweight_red */
	0.5,	/* grayweight_green */
	0.2,	/* grayweight_blue */
	2.22	/* gamma */
};

// 
// Forward declarations
//
void BmpGradientBackDiff( BmpBitmap *bmp, float *data,
    			  int x, int y, 
    			   double h, double	*vg );
static int BmpIFreePyramid( BmpBitmap *bmp );
static void BmpIFreeGreyChannel( BmpBitmap *bmp );

/*-------------------------------------------------------------------------*/
/* dummies for functions not yet implemeted for 64 bit		           */
/*-------------------------------------------------------------------------*/
#if 0
int MGWdelbmp( int x )
{
  TRCERR(("not yet implemented for 64 bit\n"));
  return 0;
}
#endif
#if 0
int BmpReadPNG(FILE *fp, CoConverter *cc, BmpBitmap **pBmp)
{
  TRCERR(("not yet implemented for 64 bit\n"));
  return MI_ENOTIMPL;
}
int BmpWritePNG(FILE *fp, BmpBitmap *bmp, CoConverter *cc,
    		unsigned opt )
{
  TRCERR(("not yet implemented for 64 bit\n"));
  return MI_ENOTIMPL;
}

void SmpInitPool()
{
  TRCERR(("not yet implemented for 64 bit\n"));
}
#endif

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpClipRectangle( 
    	BmpRectangle *A, // A: IN: rectangle to be clipped against
    	BmpRectangle *B  // B: IN/OUT: clipped rectangle
          )
{
  int rcThis=0,
      xmin=A->x0,ymin=A->y0,
      xmax=xmin+A->sx-1,ymax=ymin+A->sy-1,
      x1,y1,x2,y2,x3,y3,x4,y4;
  x1=B->x0; y1=B->y0; 
  BMP_CLIP_CONST(xmin,ymin,xmax,ymax,x1,y1);
  x2=B->x0+B->sx-1; y2=B->y0; 
  BMP_CLIP_CONST(xmin,ymin,xmax,ymax,x2,y2);
  x3=B->x0; y3=B->y0+B->sy-1; 
  BMP_CLIP_CONST(xmin,ymin,xmax,ymax,x3,y3);
  x4=B->x0+B->sx-1; y4=B->y0+B->sy-1; 
  BMP_CLIP_CONST(xmin,ymin,xmax,ymax,x4,y4);     
  BMP_SET_RECT(B,x1,y1,x2-x1+1,y3-y1+1);	  
//rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void BmpInitBitmapWindow( BmpBitmapWindow *bmw, 
    		     	 BmpBitmap *bmp,
    	                 int x0, int y0, int sx, int sy )
{
  bmw->bmp = bmp;
  bmw->x0 = x0;
  bmw->y0 = y0;
  bmw->sx = sx;
  bmw->sy = sy;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MiDeleteRegionInfo( BmpRegionInfo *ri)
{
  if( ri ) 
  {
    if( ri->riOrig ) MiDeleteRegionInfo(ri->riOrig);
    if( ri->mask ) BmpDeleteRegionMask( ri->mask );
    if( ri->bmpReg ) BmpCloseBitmap( ri->bmpReg );
    MM_free( ri );
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int CoClipRgb(float *vRgb)
{
  int k,ne=0;
  for(k=0;k<3;k++)
  {
    if(isnan(vRgb[k])||vRgb[k]<-1||vRgb[k]>2)
    {
      TRCERR(("invalid pixel RGB value %g\n",vRgb[k]));
      vRgb[k]=0;
      ne++;
    }
  }
  for(k=0;k<3;k++)
  {
    if(vRgb[k]<0) { vRgb[k]=0; ne++; } 
    if(vRgb[k]>1) { vRgb[k]=1; ne++; } 
  }
  return ne;
}

/*-------------------------------------------------------------------------
 * clip "excess" RGB coordinates (e.g. after brightness adjustment etc.)   
 * to the 0..255 range, thereby trying to conserve brightness and color    
 * value as much as possible.
 *
 * Method: "distribute" overflow of one axis evenly (wrt. to RGB perceived 
 * brightness constans) onto other axis until all axis are in range		
 *------------------------------------------------------------------------*/
void CoClipRGB( int *pr, int *pg, int  *pb, int max )
{

#define CLIP_RGB_2( x, fx_, y, fy_, dx, dy, max, strength ) \
  { \
    x += (strength)*dx; y += (strength)*dy; \
    if( x > max ) \
    { \
      e = fx_*(x-max); y += (strength)*(e/fy_); x = max; \
    }  \
    if( y > max ) \
    { \
      e = fy_*(y-max); x += (strength)*(e/fx_); y = max; \
      if( x > max ) x = max; \
    } \
  }

#define CLIP_RGB_3( r_, fr_, g_, fg_, b_, fb_, max_ ) \
  { \
    if( r_ > max_ ) \
    { \
      d = fr_ * (r_ -max_); q = b_ > 0 ? g_/b_ : g_; \
      dy = d/(fg_*q+fb_); dx = q*dy;	\
      r_ = max_; \
      CLIP_RGB_2( g_, fg_, b_, fb_, dx, dy, max_, \
	  coGlobConst.rgbClippingStrength );     \
    } \
  }

  float r = *pr, g=*pg, b=*pb;
  float	q,fr=coGlobConst.grayweight_red, 
  	  fg=coGlobConst.grayweight_green,
	  fb=coGlobConst.grayweight_blue,
	  dx,dy,d,e;

  /*-----------------------------------------------------------------------*/
  /* check r,g,b > max ?						   */
  /*-----------------------------------------------------------------------*/
  CLIP_RGB_3( r, fr, g, fg, b, fb, max );
  CLIP_RGB_3( g, fg, r, fr, b, fb, max );
  CLIP_RGB_3( b, fb, r, fr, g, fg, max );

  /*-----------------------------------------------------------------------*/
  /* check r,g,b < min (<=> -x > -max)					   */
  /*-----------------------------------------------------------------------*/
  max = 0; r = -r; g = -g; b = -b;
  CLIP_RGB_3( r, fr, g, fg, b, fb, max );
  CLIP_RGB_3( g, fg, r, fr, b, fb, max );
  CLIP_RGB_3( b, fb, r, fr, g, fg, max );
  r = -r; g = -g; b = -b;

  *pr=r; *pg=g; *pb=b;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int  CoConvRgb2Dummy( float *vLxy, float *vRgb )
{
  MT_VEC3_SET( vLxy, vRgb[0], 2*vRgb[1]-1, 2*vRgb[2]-1 );
  return 0;
}
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int CoConvDummy2Rgb( float *vRgb, float *vLxy )
{
  MT_VEC3_SET( vRgb, vLxy[0], (vLxy[1]+1)/2, (vLxy[2]+1)/2 );
  return 0;
}

/*=========================================================================
 *
 *   H C L  (Hue Croma Luminance)  color space  
 *
 *   as usggested by dcecar@home.com in comp.graphics.algorithms
 *
 *=========================================================================*/   

typedef struct tag_hcl
{
    int h; // hue
    int c; // cromma
    int l; // luminance
}HCL;

// this is what we arbitrarly set the number of distinct hues there are! (must
// be a multiple of 6)
#define HUENUM 510

// given number above this is 1*HPE6/6, 2*HPE6/6, 3*HPE6/6 etc.
#define HPE1 85
#define HPE2 170
#define HPE3 255
#define HPE4 340
#define HPE5 425
#define HPE6 510

// Hues range (0 - HUENUM-1)
#define HUEMAX 509

// max value for rgb components (also max for luminance)
#define RGBMAX 255


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int  CoConvRgb2Hcl( float *vLxy, float *vRgb )
{
  int r, g, b, cmin, cmax, rcThis=0;
  int ch, dr, dg, db;
  HCL hcl;

  r = vRgb[0]*255; g = vRgb[1]*255; b = vRgb[2]*255; 

  cmin = MIN(MIN(r, g), b);
  cmax = MAX(MAX(r, g), b);
  hcl.c = cmax-cmin; // chromma // global_ret_2
  hcl.l = (cmax+cmin)/2; // luminance // global_ret_3

  // acromatic case
  if (hcl.c == 0)
  {
   hcl.h = 0; // or some undefined value
   raiseRc( MI_OK );
  }

  ch = hcl.c/2;
  if (cmax == r)
  {
   dg = ((cmax-g)*HPE1 + ch)/hcl.c;
   db = ((cmax-b)*HPE1 + ch)/hcl.c;
   hcl.h = HPE1+db-dg;
  }
  else if (cmax == g)
  {
   dr = ((cmax-r)*HPE1 + ch)/hcl.c;
   db = ((cmax-b)*HPE1 + ch)/hcl.c;
   hcl.h = HPE3 + dr-db;
  }
  else
  {
   dr = ((cmax-r)*HPE1 + ch)/hcl.c;
   dg = ((cmax-g)*HPE1 + ch)/hcl.c;
   hcl.h = HPE5 + dg-dr;
  }

  if (hcl.h < 0) hcl.h += HUENUM;
  else if (hcl.h > HUEMAX) hcl.h -= HUENUM;

rcCatch:
  MT_VEC3_SET( vLxy, (float)hcl.l/255.0, 
      		     (float)hcl.h/(((float)HUENUM)/2.) - 1, 
      		     (float)hcl.c/127.0 - 1 );
  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int  CoConvHcl2Rgb( float *vRgb, float *vLxy )
{
  int 	r, g, b, cmax, cmin, oog;
  int 	dl;
  HCL	hcl;

  hcl.l = vLxy[0]*255.;
  hcl.h = (vLxy[1]+1)*((float)HUENUM/2.);
  hcl.c = (vLxy[2]+1)*127.;

  //
  // Do a little checking, make sure we clamp cromma if outside (for given
  // luminance).
  // Hue and Luminance wont be affected (which is where most
  // of the visual info is)
  //

  if (hcl.l < 128)
  {
   dl = 2*hcl.l;
   if (hcl.c > dl)
   {
    cmax = dl;
    cmin = 0;
    hcl.c = dl;
   }
   else
   {
    cmax = hcl.l+hcl.c/2;
    cmin = cmax-hcl.c;
   }
  }
  else
  {
   dl = 2*(255-hcl.l);
   if (hcl.c > dl)
   {
    cmax = 255;
    cmin = 255-dl;
    hcl.c = dl;
   }
   else
   {
    cmax = hcl.l+hcl.c/2;
    cmin = cmax-hcl.c;
   }
  }

  //
  // Now calculate the rgb components
  //

  if (hcl.h < HPE2)
  {
   r = cmax;
   if (hcl.h < HPE1)
   {
    g = cmin;
    b = cmin + (hcl.c*(HPE1-hcl.h)/HPE1);
   }
   else
   {
    b = cmin;
    g = cmin + (hcl.c*(hcl.h-HPE1)/HPE1);
   }
  }
  else if (hcl.h < HPE4)
  {
   g = cmax;
   if (hcl.h < HPE3)
   {
    b = cmin;
    r = cmin + (hcl.c*(HPE3-hcl.h)/HPE1);
   }
   else
   {
    r = cmin;
    b = cmin + (hcl.c*(hcl.h-HPE3)/HPE1);
   }
  }
  else
  {
   b = cmax;
   if (hcl.h < HPE5)
   {
    r = cmin;
    g = cmin + (hcl.c*(HPE5-hcl.h)/HPE1);
   }
   else
   {
    g = cmin;
    r = cmin + (hcl.c*(hcl.h-HPE5)/HPE1);
   }
  } 

  MT_VEC3_SET( vRgb, (float)r/255., (float)g/255., (float)b/255. );
  CO_CHECK_RGB_F( vRgb, oog );
  return oog ? CO_OUT_OF_GAMUT : 0; 

}


CoConverter 
  CoConvDummy = { "DUMMY", 0, NULL,
  		  CoConvDummy2Rgb, CoConvRgb2Dummy },
  CoConvHCL = { "HCL", 0, NULL,
	CoConvHcl2Rgb, CoConvRgb2Hcl };
		  
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BmpRegionInfo *MiCreateRegionInfo( int regionSize,
    					 int nFeatures,
				   int contextRegionSize,	 	
    				         unsigned options )
{
  int		rcThis = 0;
  BmpRegionInfo *padm = NULL;
  size_t	sz;

  sz = sizeof(BmpRegionInfo);

  padm  = MM_malloc(sz,MM_DUMMYTAG);
  if(!padm)
  {
    TRCERR(("MM_malloc(%ld) failed\n",(long)sz));
    raiseRc( MI_ENOMEM );
  }
  padm->mask = NULL;
  padm->bmpReg = NULL;

  padm->rContext = contextRegionSize;
  padm->mask = BmpCreateRegionMask( regionSize, 0 );
  if(! padm->mask )
  {
    TRCERR(("BmpCreateRegionMask failed\n",(long)sz));
    raiseRc( MI_ERROR );
  }
  
  padm->colInfo.bmpUV = NULL;
  padm->bmpMatchCNH = NULL;

  BMP_SET_RECT( &padm->window, 0, 0, 0, 0 );

  padm->riOrig = NULL;

  TRC(("MiCreateRegionInfo regSize=%d nFeatures=%d -> %p\n",
	(int)regionSize,(int)nFeatures,padm));
rcCatch:
  if( rcThis )
  {
    if( padm ) MiDeleteRegionInfo( padm );
    exit(1);
  }
  return padm;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
 BmpRegionMask *BmpCreateRegionMask( int regionSize,
    				      unsigned options )
{
  unsigned char	*pmask;
  BmpRegionMask *padm=NULL;
  int s,x,y,r,n;
  size_t  sz;

  r = regionSize;
  s = r/2;
  sz = sizeof(BmpRegionMask) + sizeof(unsigned char)*r*r;

  padm  = MM_malloc(sz,MM_DUMMYTAG);
  if(!padm)
  {
    TRCERR(("MM_malloc(%ld) failed\n",(long)sz));
    exit(1);
  }

  pmask = (unsigned char *)padm+sizeof(BmpRegionMask);
  padm->mask= pmask;
  padm->regSize = regionSize;

  n=0;
  for(x=0;x<r;x++) 
  {
    for(y=0;y<r;y++)
    {
      if(sqrt((x-s)*(x-s)+(y-s)*(y-s))<=s+0.5)
      {
	pmask[x+y*r]=1; n++;
      }
      else
	pmask[x+y*r]=0;
    }
  }

  /**** special case r=5 => make fat circle */
  if(r==5)
  {
    pmask[0]=pmask[4]=pmask[20]=pmask[24]=1;
    n+=4;
  }

  padm->nSetPixels = n;

  if(LogLevel >= 3)
  {
    printf("--------> r=%d\n",r); 
    for(x=0;x<r;x++) 
    {
      for(y=0;y<r;y++) printf( pmask[x+y*r] ? "*" : "." );
      printf("\n");
    }
  }

  TRC(("BmpCreateRegionMask regSize=%d npix=%d\n",
	regionSize,padm ? padm->nSetPixels: -1));

  return padm;
}
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void BmpDeleteRegionMask( BmpRegionMask *pmask )
{
  if( pmask ) MM_free( pmask );
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void BmpPrintRegion( BmpRegionInfo *regInfo )
{
  int   r,x,y,rx,ry,xmod;
  BmpGreylevel	 *ppSrc;
  unsigned char *pmask;
  BmpBitmap *bmp=regInfo->bmp;
  BmpRegionMask  *regMask;
  BmpRscanInfo   *rsi = &regInfo->rscanInfo;

  regMask = regInfo->mask;
  r = regMask->regSize;
  ppSrc = bmp->dataGrey + rsi->bmpOffs;
  xmod = rsi->xmod;
  rx=rsi->rx; ry=rsi->ry;
  pmask = regMask->mask;

  printf("REGION(bmp=%dx%d) x0,y0=%d,%d r=%d x1,y1=%d,%d rx,ry=%dx%d\n",
      (int)bmp->sx,(int)bmp->sy,
      regInfo->x0,regInfo->y0,r,regInfo->x1,regInfo->y1,rx,ry);

  printf("Y/X> "); for(x=0;x<rx;x++) printf("%3d ",x); printf("\n");
  for(y=0;y<ry;y++, ppSrc+=xmod ) 
  {
    printf("%3d> ",y);
    for(x=0;x<rx;x++, ppSrc++)
    {
      if( *pmask++ ) 
	printf("%3d ",(int) *ppSrc);
      else
	printf("    ");
    }
    printf("\n");
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BmpPolygon *BmpCreatePolygon( int nVert, int *vertices )
{
  int *pv;
  BmpPolygon *poly=NULL;
  size_t sz;

  sz = sizeof(BmpPolygon) + sizeof(int)*nVert*2;

  poly  = MM_malloc(sz,MM_DUMMYTAG);
  if(!poly)
  {
    TRCERR(("MM_malloc(%ld) failed\n",(long)sz));
    exit(1);
  }

  pv = (int *)((unsigned char *)poly+sizeof(BmpPolygon));
  poly->vertices = pv;
  poly->nVert = nVert;
  memcpy( poly->vertices, vertices, sizeof(int)*nVert*2 );

  TRC(("BmpCreatePolygon %d\n",nVert));

  return poly;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpDrawLine0( BmpBitmap *bmp, int channel, 
    		  int x1, int y1, int x2, int y2,
		  int cr, int cg, int cb )
{
  int rcThis=0, sx=x2-x1+1, sy=y2-y1+1, tmp,x,y;
  float g=(cr+cg+cb)/3;

  BMP_CLIP_CONST(0,0,bmp->sx-1,bmp->sy-1,x1,y1);
  BMP_CLIP_CONST(0,0,bmp->sx-1,bmp->sy-1,x2,y2);

  if( x2 < x1 ) SWAP( x1, x2, tmp );
  if( y2 < y1 ) SWAP( y1, y2, tmp );

  if( sy <= 1 )
  {
    if(channel==BMP_GREY)
      for(x=x1;x<=x2;x++) BMP_DPIXEL( bmp, bmp->dataGrey, x, y1 ) = g;
    else
      for(x=x1;x<=x2;x++) BMP_SETPIXEL( bmp, x, y1, cr, cg, cb );
  }
  if( sx <= 1 )
  {
    if(channel==BMP_GREY)
      for(y=y1;y<=y2;y++) BMP_DPIXEL( bmp, bmp->dataGrey, x1, y ) = g;
    else
      for(y=y1;y<=y2;y++) BMP_SETPIXEL( bmp, x1, y, cr, cg, cb );
  }
rcCatch:
  return rcThis;
}


#if 0
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpDrawLine2( BmpBitmap *bmp, int chan,
    		   int x0,  int y0, 
		   int x1,  int y1, 
		  int col ) 
// https://rosettacode.org/wiki/Xiaolin_Wu%27s_line_algorithm#C
{
  UT_ASSERT0(chan==BMP_RGB);
#define ipart_(X) ((int)(X))
#define round_(X) ((int)(((float)(X))+0.5f))
#define fpart_(X) (((float)(X))-(float)ipart_(X))
#define rfpart_(X) (1.0f-fpart_(X))

#define plot_(x,y,d)\
  {\
    if(BMP_IN_RANGEB(x,y,bmp))\
    {\
      float f = MIN(d,1.0f);\
      int col2 = BMP_MKRGB(\
	MAX(MIN((int)((float)BMP_RED(col) * f), 255),0),\
	MAX(MIN((int)((float)BMP_GREEN(col) * f),255),0),\
	MAX(MIN((int)((float)BMP_BLUE(col) * f),255),0));\
      BMP_SETPIXEL_RGB(bmp,((int)(x)),((int)(y)),col2);\
    }\
  }

  unsigned int tmp;

  int steep = ABS(y1 - y0) > ABS(x1 - x0);

  if(steep)
  {
    SWAP(x0, y0, tmp)
    SWAP(x1, y1, tmp)
  }
  if (x0 > x1) 
  {
    SWAP(x0, x1, tmp)
    SWAP(y0, y1, tmp)
  }

  float dx = x1 - x0,
        dy = y1 - y0;

  float gradient = fabsf(dx)<MT_NUM_EPS ? 1.0f : dy/dx;

// handle first endpoint
  float  
    xend = round_(x0),
    yend = y0 + gradient * (xend - x0),
    xgap = rfpart_(x0 + 0.5f),
    xpxl1 = xend, // this will be used in the main loop
    ypxl1 = ipart_(yend);

  if (steep) {
      plot_(ypxl1,   xpxl1, rfpart_(yend) * xgap);
      plot_(ypxl1+1, xpxl1,  fpart_(yend) * xgap);
  }
  else {
      plot_(xpxl1, ypxl1  , rfpart_(yend) * xgap);
      plot_(xpxl1, ypxl1+1,  fpart_(yend) * xgap);
  }
    

  // handle second endpoint
  xend = round_(x1);
  yend = y1 + gradient * (xend - x1);
  xgap = fpart_(x1 + 0.5f);
  float xpxl2 = xend, //this will be used in the main loop
	ypxl2 = ipart_(yend);
  if (steep) {
      plot_(ypxl2  , xpxl2, rfpart_(yend) * xgap);
      plot_(ypxl2+1.f, xpxl2,  fpart_(yend) * xgap); 
  }
  else {
      plot_(xpxl2, ypxl2,  rfpart_(yend) * xgap);
      plot_(xpxl2, ypxl2+1.f, fpart_(yend) * xgap); }

  float intery = yend + gradient; // first y-intersection for the main loop

// main loop
  if(steep) 
  {
      for (int x = (int)xpxl1 + 1; x <= (int)xpxl2 - 1; x++ )
      {
	plot_(ipart_(intery)  , x, rfpart_(intery));
	plot_(ipart_(intery)+1, x,  fpart_(intery));
	intery = intery + gradient;
      }
  } 
  else 
  { 
      for (int x = (int)xpxl1 + 1; x<= (int)xpxl2 - 1; x++ )
      {
	plot_(x, ipart_(intery),  rfpart_(intery));
	plot_(x, ipart_(intery)+1, fpart_(intery));
	intery = intery + gradient;
      }
  }

  return 0;
}
#endif


#if 1
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpDrawLine2( BmpBitmap *bmp, int chan,
    		   int x1,  int y1, 
		   int x2,  int y2, 
		  int col ) 
// https://rosettacode.org/wiki/Xiaolin_Wu%27s_line_algorithm#C
{
  UT_ASSERT0(chan==BMP_RGB);
#define ipart_(X) ((int)(X))
#define round_(X) ((int)(((float)(X))+0.5f))
#define fpart_(X) (((float)(X))-(float)ipart_(X))
#define rfpart_(X) (1.0f-fpart_(X))

#if 0
#define plot_(x,y,d)\
  {\
    if(BMP_IN_RANGEB(x,y,bmp))\
    {\
      float f = /*sqrtf*/(MIN((d),1.0f));\
      if(f>0.5f)\
      {\
	f=1.0f;\
      int col2 = BMP_MKRGB(\
	MAX(MIN((int)((float)BMP_RED(col) * f), 255),0),\
	MAX(MIN((int)((float)BMP_GREEN(col) * f),255),0),\
	MAX(MIN((int)((float)BMP_BLUE(col) * f),255),0));\
      BMP_SETPIXEL_RGB(bmp,x,y,col);\
      }\
    }\
  }
#else
#define plot_(x,y,d)\
  {\
    if(BMP_IN_RANGEB(x,y,bmp))\
    {\
      float f = /*sqrtf*/(MIN((d),1.0f));\
      if(f>0.5f)\
      {\
      BMP_SETPIXEL_RGB(bmp,x,y,col);\
      }\
    }\
  }

#endif
  unsigned int tmp;

  float dx = (float)x2 - (float)x1;
  float dy = (float)y2 - (float)y1;
  if ( fabs(dx) > fabs(dy) ) {
    if ( x2 < x1 ) {
      SWAP(x1, x2, tmp);
      SWAP(y1, y2, tmp);
    }
    float gradient = fabsf(dx)>MT_NUM_EPS ? dy / dx : 1.0f;
    float xend = round_(x1);
    float yend = y1 + gradient*(xend - x1);
    float xgap = rfpart_(x1 + 0.5f);
    int xpxl1 = xend;
    int ypxl1 = ipart_(yend);
    plot_(xpxl1, ypxl1, rfpart_(yend)*xgap);
    plot_(xpxl1, ypxl1+1, fpart_(yend)*xgap);
    float intery = yend + gradient;
 
    xend = round_(x2);
    yend = y2 + gradient*(xend - x2);
    xgap = fpart_(x2+0.5);
    int xpxl2 = xend;
    int ypxl2 = ipart_(yend);
    plot_(xpxl2, ypxl2, rfpart_(yend) * xgap);
    plot_(xpxl2, ypxl2 + 1, fpart_(yend) * xgap);
 
    int x;
    for(x=xpxl1+1; x < xpxl2; x++) {
      plot_(x, ipart_(intery), rfpart_(intery));
      plot_(x, ipart_(intery) + 1, fpart_(intery));
      intery += gradient;
    }
  } else {
    if ( y2 < y1 ) {
      SWAP(x1, x2, tmp);
      SWAP(y1, y2, tmp);
    }
    float gradient = fabsf(dy)>MT_NUM_EPS ? dx / dy : 1.0f;
    float yend = round_(y1);
    float xend = x1 + gradient*(yend - y1);
    float ygap = rfpart_(y1 + 0.5f);
    int ypxl1 = yend;
    int xpxl1 = ipart_(xend);
    plot_(xpxl1, ypxl1, rfpart_(xend)*ygap);
    plot_(xpxl1 + 1, ypxl1, fpart_(xend)*ygap);
    float interx = xend + gradient;
 
    yend = round_(y2);
    xend = x2 + gradient*(yend - y2);
    ygap = fpart_(y2+0.5);
    int ypxl2 = yend;
    int xpxl2 = ipart_(xend);
    plot_(xpxl2, ypxl2, rfpart_(xend) * ygap);
    plot_(xpxl2 + 1, ypxl2, fpart_(xend) * ygap);
 
    int y;
    for(y=ypxl1+1; y < ypxl2; y++) {
      plot_(ipart_(interx), y, rfpart_(interx));
      plot_(ipart_(interx) + 1, y, fpart_(interx));
      interx += gradient;
    }
  }
  return 0;
}
#endif

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpDrawLine( BmpBitmap *bmp, int chan,
    		 int x1, int y1, int x2, int y2, int col, 
    		 int endp ) 
{
  int rcThis=0, x, y, sx=bmp->sx,sy=bmp->sy;
  int dx, dy;
  int p;
  float colf[3]={0};

  UT_ASSERT(chan==BMP_RGB||chan==BMP_XYZ||chan==BMP_GREY);

  if(chan==BMP_XYZ)
  {
    colf[0]=BMP_RED(col);
    colf[1]=BMP_GREEN(col);
    colf[2]=BMP_BLUE(col);
  }

  dx = x2 - x1;
  dy = y2 - y1;


  /* choose variable to iterate on */
  if (abs(dx)>abs(dy)) {
    int dx2, dy2, cpy;

    /* sort by x */
    if (x1 > x2) {
      int t;
      t = x1; x1 = x2; x2 = t;
      t = y1; y1 = y2; y2 = t;
    }

    dx = abs(dx);
    dx2 = dx*2;
    dy = y2 - y1;

    if (dy<0) {
      dy = -dy;
      cpy = -1;
    } else {
      cpy = 1;
    }
    dy2 = dy*2;
    p = dy2 - dx;


    y = y1;
    for(x=x1; x<x2-1; x++) {
      if (p<0) {
        p += dy2;
      } else {
        y += cpy;
        p += dy2-dx2;
      }
      if(y>0&&y<sy&&x+1>0&&x+1<sx)
      {
	if(chan==BMP_RGB)
	{
	  BMP_SETPIXEL_RGB(bmp, x+1, y, col);
	}
	else if(chan==BMP_GREY)
	{
	  BMP_SETPIXEL_GREY(bmp, x+1, y, col);
	}
	else
	{
	  BMP_SETPIXEL_XYZ(bmp, x+1, y, colf);
	}
      }
    }
  } else {
    int dy2, dx2, cpx;

    /* sort bx y */
    if (y1 > y2) {
      int t;
      t = x1; x1 = x2; x2 = t;
      t = y1; y1 = y2; y2 = t;
    }

    dy = abs(dy);
    dx = x2 - x1;
    dy2 = dy*2;

    if (dx<0) {
      dx = -dx;
      cpx = -1;
    } else {
      cpx = 1;
    }
    dx2 = dx*2;
    p = dx2 - dy;

    x = x1;

    for(y=y1; y<y2-1; y++) {
      if (p<0) {
        p  += dx2;
      } else {
        x += cpx;
        p += dx2-dy2;
      }
      if(x>0&&(x+1<sx)&&y>0&&y<sy)
      {
	if(chan==BMP_RGB)
	{
          BMP_SETPIXEL_RGB(bmp, x+1, y, col);
	}
	else if(chan==BMP_GREY)
	{
	  BMP_SETPIXEL_GREY(bmp, x+1, y, col);
	}
	else
	{
          BMP_SETPIXEL_XYZ(bmp, x+1, y, colf);
	}
      }
    }
  }

  if(chan==BMP_RGB)
  {
  if (endp) {
    if(BMP_IN_RANGEB(x1,y1,bmp))
      BMP_SETPIXEL_RGB(bmp, x1, y1, col);
    if(BMP_IN_RANGEB(x2,y2,bmp))
      BMP_SETPIXEL_RGB(bmp, x2, y2, col);
  } else {
    if (x1 != x2 || y1 != y2)
      if(BMP_IN_RANGEB(x1,y1,bmp))
        BMP_SETPIXEL_RGB(bmp, x1, y1, col);
  }
  if (endp&2) {
    if(BMP_IN_RANGEB(x2,y2,bmp))
      BMP_SETPIXEL_RGB(bmp, x2, y2, BMP_MKRGB(255,255,255));
  }
  }

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpDrawArrow( BmpBitmap *bmp, int chan,
    		  int x1, int y1, int x2, int y2, int col, 
    		  float head_size ) 
{
  static int auxValid=0;
  static float
    arr_ss1,arr_cs1,
    arr_ss2,arr_cs2;

  if(!auxValid)
  {
     arr_ss1 = sinf(MT_DEG2RAD(30)),
     arr_cs1 = cosf(MT_DEG2RAD(30)),
     arr_ss2 = sinf(MT_DEG2RAD(-30)),
     arr_cs2 = cosf(MT_DEG2RAD(-30));
     auxValid = TRUE;
  }

  float dx=-(x2-x1),
	dy=-(y2-y1),
	dx2,dy2;
  
  BmpDrawLine(bmp,chan, x1,y1,x2,y2, col, FALSE);
#if 1
  if(dx*dx+dy*dy>0.0000001)
  {	 
    MT_VEC2_ROTATE2(dx2,dy2, dx,dy, arr_ss1, arr_cs1);
    
    BmpDrawLine(bmp,chan, 
	        x2,y2, x2+dx2*head_size,y2+dy2*head_size, 
		col, TRUE );

    MT_VEC2_ROTATE2(dx2,dy2, dx,dy, arr_ss2, arr_cs2);


    BmpDrawLine(bmp,chan,
	        x2,y2, x2+dx2*head_size,y2+dy2*head_size, 
		col, TRUE );	  
  }
#endif
  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpDrawCircle( BmpBitmap *bmp, int chan,
    		  int x0, int y0, int dmax,
		  int crgb )
{
  int x2,y2,dyf,dx,
      cr=BMP_RED(crgb), cg=BMP_GREEN(crgb), cb=BMP_BLUE(crgb);
  
  for(dx=-dmax;dx<=dmax;dx++)
  {
    dyf = sqrt(dmax*dmax-dx*dx);

    x2 = x0+dx;
    y2 = (int)(dyf+0.5)+y0;
    if(BMP_IN_RANGEB(x2,y2,bmp))
    {	
      BMP_SETPIXEL( bmp, x2, y2, cr, cg, cb );
    }
    y2 = y0-(int)(dyf+0.5);
    if(BMP_IN_RANGEB(x2,y2,bmp))
    {	
      BMP_SETPIXEL( bmp, x2, y2, cr, cg, cb );
    }
  }
  BMP_MODIFIED(bmp);
  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpDrawCross( BmpBitmap *bmp, int chan,
    		  int x0, int y0, int dmax,
		  int crgb )
{
  int cr=BMP_RED(crgb), cg=BMP_GREEN(crgb), cb=BMP_BLUE(crgb);
  BmpDrawLine0(bmp,chan,x0,y0-dmax,x0,y0+dmax,cr,cg,cb);
  BmpDrawLine0(bmp,chan,x0-dmax,y0,x0+dmax,y0,cr,cg,cb);
  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpDrawRect( BmpBitmap *bmp, int channel, 
    		 int x1, int y1, int sx, int sy,
		 int cr, int cg, int cb)
{
  int rcThis=0, x2=x1+sx-1, y2=y1+sy-1;

  if(BMP_IN_RANGEB(x1,y1,bmp)&&BMP_IN_RANGEB(x2,y2,bmp))
  {
    BmpDrawLine0( bmp, channel, x1, y1, x2, y1, cr, cg, cb );
    BmpDrawLine0( bmp, channel, x2, y1, x2, y2, cr, cg, cb );
    BmpDrawLine0( bmp, channel, x2, y2, x1, y2, cr, cg, cb );
    BmpDrawLine0( bmp, channel, x1, y2, x1, y1, cr, cg, cb );
  }

/*rcCatch:*/
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpDrawFilledRect( BmpBitmap *bmp, int channel, 
    		 int x1, int y1, int sx, int sy,
		 int cr, int cg, int cb)
{
  int rcThis=0,x,y,crgb=BMP_MKRGB(cr,cg,cb);
  float colf[3]={0};

  if(channel==BMP_XYZ)
  {
    colf[0]=cr;
    colf[1]=cg;
    colf[2]=cb;
  }

  for(y=y1;y<y1+sy;y++)
    for(x=x1;x<x1+sx;x++)
    {
      if(!BMP_IN_RANGEB(x,y,bmp)) continue;
      if(channel==BMP_RGB)
      {
	BMP_SETPIXEL_RGB(bmp,x,y,crgb);
      }
      else if(channel==BMP_GREY)
      {
	BMP_DPIXEL( bmp, bmp->dataGrey, x, y ) = cg;
      }
      else
      {
	BMP_SETPIXEL_XYZ(bmp, x, y, colf);
      }
    }
  BMP_MODIFIED(bmp);

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void BmpDeletePolygon( BmpPolygon *poly )
{
  if( poly ) MM_free( poly );
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpSaveBitmapGWBMP( BmpBitmap *bmp, int gwWin, char *filename )
{
  int rcThis=0;
#ifndef BMP_WITH_GRWIN
  TRCERRR(("not implemented for 64 bit.\n"),MI_ERROR);
#else
  int rc;

  MiTranspathCyg2Win(filename);
  if( gwWin > 0 )
    rc = GWselect( gwWin );
  TRC1(("saving bitmap as %dx%d as '%s'\n",
	bmp->sx,bmp->sy,filename));
  rc = GWsavebmp( bmp->GWnr, filename, strlen(filename) );
  if(rc==0) TRCERRR(("GWsavebmp '%s' failed\n",filename),MI_ERROR);
#endif
rcCatch:
  return(rcThis);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpSaveBitmapPNG( BmpBitmap *bmp, char *fname, CoConverter *cc,
    			unsigned opt )
{
  int rcThis=0,ok,use_16_bit = (opt&BMP_16BIT) ? TRUE : FALSE,opt2;
  FILE *fp=NULL;
  
  fp = fopen(fname,"wb");
  if(!fp) 
    TRCERRR(("can't open file '%s' errno=%d \n",fname,errno),MI_ERROR);

  opt2 = use_16_bit && (bmp->dataGrey || bmp->dataUV || bmp->dataFloat[0] ||
			   ((bmp->dataUshort0 || bmp->dataUshort[0]) && 
			    (opt & BMP_USHORT))) ? BMP_16BIT : 0;
  if(opt&BMP_GREY_ONLY) opt2 |= BMP_GREY_ONLY;


  if((opt & BMP_XYZ)&&(opt & BMP_RGB)) opt2 |= BMP_RGB|BMP_XYZ;
  if((opt & BMP_USHORT)&&(opt & BMP_RGB)) opt2 |= BMP_RGB|BMP_USHORT;
  if(opt&BMP_COMPR) opt2 |= BMP_COMPR;
  if(opt&BMP_ALPHA) opt2 |= BMP_ALPHA;

  ok = BmpWritePNG( fp, bmp, cc, opt2 );
  if(!ok) raiseRc(MI_ERROR);

rcCatch:
  if(fp) fclose(fp);
  return(rcThis);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpSuffixToFormat( char *suff )
{
  if( UtStrCaseCmp(suff,"bmp") == 0 )
  {
    return BMP_WINBMP;
  }
  else if( UtStrCaseCmp( suff, "tif" ) == 0 ||
           UtStrCaseCmp( suff, "tiff" ) == 0)
  {
    return BMP_TIF;
  }
  else if( UtStrCaseCmp( suff, "u48" ) == 0 )
  {
    return BMP_U48;
  }
  else if( UtStrCaseCmp( suff, "png" ) == 0 )
  {
    return BMP_PNG;
  }
  else if( UtStrCaseCmp( suff, "jpg" ) == 0 )
  {
    return BMP_JPG;
  }
  else
  {
    return BMP_NULL;
  }
}
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
char *BmpFormatToSuffix( int format )
{
  switch(format)
  {
    case BMP_NULL: 
    default:
      return ".";
    case BMP_TIF:
      return ".tif";
    case BMP_WINBMP:
      return ".bmp";
    case BMP_PNG:
      return ".png";
    case BMP_U48:
      return ".u48";
    case BMP_JPG:
      return ".jpg";
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BmpBitmap *BmpLoadBitmap( char *fname, unsigned opt )
{
  return BmpCreateBitmap2( fname, NULL, 0, 0, 0, fname,  NULL, opt);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpSaveBitmap( BmpBitmap *bmp, 
    		   int gwNr,
    		   char *filenameIn, 
		   int format,
		   CoConverter *cc,
		   unsigned opt )
{
  int rcThis=0;
  char filename[MI_MAX_FNAME_LEN+1];

  if(!bmp) TRCERRR(("no bitmap !\n"),MI_ERROR);

  strcpy( filename, filenameIn );
  if( format == BMP_NULL )
  {
    format = BmpSuffixToFormat( UtFileSuffix(filename) );    
  }
  if(format == BMP_NULL )  
  {
    /* default: winbmp */
    format = BMP_WINBMP;
  }

  UtSetFilenameSuffix(filename,filename,BmpFormatToSuffix(format));
  
  switch(format)
  {
    case BMP_U48:
      TRCERRR(("deprecated format 'U48'"),MI_ERROR);
      break;
    case BMP_PNG:
      return BmpSaveBitmapPNG( bmp, filename, cc, BMP_16BIT | opt);
      break;
    case BMP_JPG:
      return BmpSaveBitmapJPG( bmp, filename, 0, opt);
      break;
      break;
    default:
      TRCERRR(("unkown bitmap format %d\n",format),MI_EINVARG);
      break;
  }
rcCatch:  
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpDeleteBitmap( BmpBitmap **p_bmp )
{
  int rcThis=0;
  if(*p_bmp) 
  {
    BmpCloseBitmap(*p_bmp);
    *p_bmp=NULL;
  }
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpCloseBitmap( BmpBitmap *bmp )
{
  // deprecated function: use BmpDeleteBitmap() instead
  int i;
  if(bmp)
  {
    BmpDeleteBitmap(&bmp->stereo.bmpRight);
    BmpDeleteBitmap(&bmp->stereo.bmpAnaglyph);
    BmpIFreePyramid(bmp);
    if(bmp->GWnr) MGWdelbmp(bmp->GWnr);
    if(bmp->dataUshort0) MM_free(bmp->dataUshort0);
    if(bmp->dataPtr) MM_free(bmp->dataPtr);
    if(bmp->dataAlpha) MM_free(bmp->dataAlpha);
    BmpIFreeGreyChannel(bmp);
    FREEMEM(bmp->dataValid);
    if(bmp->dataUV) 
    {
      FREEMEM(bmp->dataUV);
    }
    if(bmp->dataZbuf) MM_free(bmp->dataZbuf);
    if(bmp->data) 
    {
      FREEMEM(bmp->data);
    }
    for(i=0;i<BMP_MAX_FCHAN;i++) 
    {
      if(bmp->dataFloat[i]) FREEMEM(bmp->dataFloat[i]);
    }
    for(i=0;i<BMP_MAX_DCHAN;i++) 
    {
      FREEMEM(bmp->dataDbl[i]);
    }
    for(i=0;i<BMP_MAX_UCHAN;i++) 
    {
      FREEMEM(bmp->dataUshort[i]);
    }
    for(i=0;i<BMP_MAX_UCHAN;i++) 
    {
      FREEMEM(bmp->dataUbyte[i]);
    }
    for(i=0;i<BMP_MAX_ICHAN;i++) 
      if(bmp->dataInt[i]) MM_free(bmp->dataInt[i]);
    if(bmp->dataFeat) MtDelTab(bmp->dataFeat);
    if(bmp->dataDFeat) MtDelMat(bmp->dataDFeat);
    for(i=0;i<3;i++) 
      if(bmp->dataN[i]) MM_free(bmp->dataN[i]);
    if(bmp->dataZ) MM_free(bmp->dataZ);
    if(bmp->dataRho) MM_free(bmp->dataRho);
    bmp->eyecatcher = BMP_EYE_FREE;
    MM_free(bmp);
  }  
  return(0);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpGetDataChannelOld(BmpBitmap *B, int chan, int ichan,
    		      float **p_dataFloat, uint16_t **p_dataUshort )
{
  int rcThis=0;
  float *data=NULL;
  uint16_t *dataUs=NULL;
  BMP_GET_CHANNEL(B,chan,ichan,data,dataUs);
  if(p_dataFloat) *p_dataFloat=data;
  if(p_dataUshort) *p_dataUshort=dataUs;
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpSwapDataPtr( BmpBitmap *B1, BmpBitmap *B2, 
    		    int chan, int ichan )
{
  int rcThis=0;
  
  UT_ASSERT(BMP_DIM_EQUAL(B1,B2));

  if(chan==BMP_GREY)
  {
    SWAP_PTR(B1->dataGrey,
	     B2->dataGrey);
  }
  else if(chan==BMP_Z)
  {
    SWAP_PTR(B1->dataZ,
	     B2->dataZ);
  }
  else if(chan==BMP_FLOAT)
  {
    SWAP_PTR(B1->dataFloat[ichan],  
	     B2->dataFloat[ichan]);
  }
  else if(chan==BMP_USHORT)
  {
    SWAP_PTR(B1->dataUshort[ichan],  
	     B2->dataUshort[ichan]);
  }
  else if(chan==BMP_UBYTE)
  {
    SWAP_PTR(B1->dataUbyte[ichan],  
	     B2->dataUbyte[ichan]);
  }
  else
  {
    UT_ASSERT(FALSE);
  }
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpGetDataPtrAddr(BmpBitmap *B, int chan, int ichan,
    		       unsigned opt,
    		       BmpData *data )
{
  int rcThis=0;
  BmpData data0;
  if(!data) data=&data0;
  memset(data,0,sizeof(*data));
  data->chan = chan;
  data->ichan = ichan;
  switch(chan)
  {
    case BMP_Z: data->ppFloat = &B->dataZ; break;
    case BMP_GREY: data->ppFloat = &B->dataGrey; break;
    case BMP_FLOAT: 
      UT_ASSERT(ichan>=0&&ichan<BMP_MAX_FCHAN);
      data->ppFloat = &B->dataFloat[ichan];
      break;
    case BMP_USHORT: 
      UT_ASSERT(ichan>=0&&ichan<BMP_MAX_UCHAN);
      data->ppUshort = &B->dataUshort[ichan];
      break;
    case BMP_UBYTE: 
      UT_ASSERT(ichan>=0&&ichan<BMP_MAX_UCHAN);
      data->ppUbyte = &B->dataUbyte[ichan];
      break;
    case BMP_VALID: 
      UT_ASSERT(ichan==0);
      data->ppUbyte = &B->dataValid;
      break;
    default:
      TRCERRR(("invalid data channel %d %d\n",chan,ichan),MI_ERROR);
      break;
  }
  data->pFloat = data->ppFloat ? *data->ppFloat : NULL;
  data->pUshort = data->ppUshort ? *data->ppUshort : NULL;
  data->pUbyte = data->ppUbyte ? *data->ppUbyte : NULL;

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpCopyDataChannel( BmpBitmap *S, int chan_s, int ichan_s, // Source
    			BmpBitmap *D, int chan_d, int ichan_d, // Destination
			unsigned opt )
{
  LongInt ii;
  int rc,rcThis=0;
  BmpData src,dst;
  UT_ASSERT(S&&D);
  UT_ASSERT(BMP_DIM_EQUAL(S,D));
  rc = BmpGetDataChannel(S,chan_s,ichan_s,0,&src);
  if(rc) raiseRc(rc);
  rc = BmpGetDataChannel(D,chan_d,ichan_d,0,&dst);
  if(rc) raiseRc(rc);
  if(src.pFloat)
  {
    UT_ASSERT(dst.pFloat);
    #pragma omp parallel for private(ii)
    for(ii=0;ii<S->sx*S->sy;ii++) dst.pFloat[ii] = src.pFloat[ii];
  }
  else if(src.pUshort)
  {
    UT_ASSERT(dst.pUshort);
    #pragma omp parallel for private(ii)
    for(ii=0;ii<S->sx*S->sy;ii++) dst.pUshort[ii] = src.pUshort[ii];
  }
  else
  {
    UT_ASSERT(FALSE);
  }
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpGetDataChannel(BmpBitmap *B, int chan, int ichan,
    		       unsigned opt,
    		       BmpData *data )
{
  int rcThis=0,rc;
  BmpData data0;  
  if(!data) data=&data0;

  memset(data,0,sizeof(*data));
  data->pFloat = NULL;
  data->pUshort = NULL;
  data->pUbyte = NULL;

  data->chan = chan;
  data->ichan = ichan;

  if( opt & BMP_CHECK )
  {
    rc = BmpGetDataPtrAddr(B,chan,ichan,opt,data);
    if(rc) raiseRc(rc);
    if(!(data->pFloat||data->pUshort||data->pUbyte)) 
      raiseRc(MI_ENOTFOUND);
    raiseRc(MI_OK);
  }

  switch(chan)
  {
    case BMP_Z: 
      if((data->pFloat =  BmpGetZChannel(B,opt))==NULL) 
	raiseRc(MI_ERROR);
      break;
    case BMP_GREY: 
      if((data->pFloat =  BmpGetGreyChannel(B,opt))==NULL) 
	raiseRc(MI_ERROR);
      break;
    case BMP_ALPHA: 
      if((data->pFloat =  BmpGetAlphaChannel(B,opt))==NULL) 
	raiseRc(MI_ERROR);
      break;
    case BMP_FLOAT: 
      if((data->pFloat =  BmpGetFloatChannel(B,ichan,opt))==NULL) 
	raiseRc(MI_ERROR);
      break;
    case BMP_USHORT: 
      if((data->pUshort =  BmpGetUshortChannel(B,ichan,opt))==NULL) 
	raiseRc(MI_ERROR);
      break;
    case BMP_UBYTE: 
      if((data->pUbyte =  BmpGetUbyteChannel(B,ichan,opt))==NULL) 
	raiseRc(MI_ERROR);
      break;
    case BMP_VALID: 
      if((data->pUbyte =  BmpGetValidChannel(B,opt))==NULL) 
	raiseRc(MI_ERROR);
      break;
    default:
      TRCERRR(("invalid data channel %d %d\n",chan,ichan),MI_ERROR);
      break;
  }

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BmpBitmap *BmpNewBitmap0( int sx_, int sy_, unsigned opt,
    			  char *mmTag )
{
  BmpBitmap *bmp = NULL;
  size_t	szNeeded;
  int		i;
  LongInt sx=sx_,sy=sy_;

  TRC3(("BmpNewBitmap: %dx%d opt=%x\n",(int)sx,(int)sy,(unsigned)opt));
  if( opt == 0 ) opt |= BMP_RGB;
//  ALLOCMEM( bmp,sizeof(BmpBitmap));
  ALLOCOBJ0_TAG(bmp, mmTag);

  bmp->optOrig = opt;
  bmp->name[0] = '\0';
  bmp->data = NULL;
  bmp->doRegenerate = TRUE;
  bmp->validChannels = 0;
  bmp->sx = sx;
  bmp->sy = sy;
  bmp->vsx=sx;
  bmp->vsy=sy;
  bmp->vx0=bmp->vy0=0;
  bmp->poffMax = sx*sy-1;
  bmp->GWnr = 0;
  bmp->data = NULL;
  bmp->dataGrey = NULL;
  bmp->dataUV = NULL;
  bmp->dataAlpha = NULL;
  bmp->dataZbuf = NULL;
  bmp->dataUshort0 = NULL;
  bmp->dataPtr = NULL;
  bmp->dataFeat = NULL;
  bmp->dataDFeat = NULL;
  bmp->dataZ = NULL;
  bmp->dataValid = NULL;
  MT_STAT_INIT(&bmp->stat);
  for(i=0;i<3;i++) bmp->dataN[i]=NULL;
  bmp->colConv = NULL;
  for(i=0;i<BMP_MAX_FCHAN;i++) bmp->dataFloat[i]=NULL;
  bmp->nFloatChan=-1;
  for(i=0;i<BMP_MAX_UCHAN;i++) bmp->dataUshort[i]=NULL;
  for(i=0;i<BMP_MAX_UCHAN;i++) bmp->dataUbyte[i]=NULL;
  for(i=0;i<BMP_MAX_ICHAN;i++) bmp->dataInt[i]=NULL;
  bmp->pyramid.n_level = 0;
  for(i=0;i<BMP_MAX_PYRA_LEVEL;i++)
  {
    bmp->pyramid.gaussian[i] = NULL;
    bmp->pyramid.laplacian[i] = NULL;
  }
  bmp->pyramid.laplace_valid = 0;
  bmp->pyramid.krn = NULL;
  bmp->eyecatcher = BMP_EYE_VALID;
  bmp->guid = 0;
  bmp->icu.header = NULL;
  bmp->icu.smpScale = 1.0;
  bmp->icu.smpOffs = 0.0;

  if( opt & BMP_RGB )
  {
    szNeeded = (sx*(LongInt)sy+1)*sizeof(int);
    ALLOCMEM_TAG( bmp->data, szNeeded, mmTag);
    /* set end-marker for fast loop access */	
    bmp->data[0] = bmp->data[(LongInt)sx*sy] = -1;  
    if( opt & BMP_CLEAR ) 
    {
      memset( bmp->data, 0, szNeeded );
    }
  }
  if( opt & BMP_ALPHA )
  {
    BmpGetAlphaChannel( bmp, opt );
  }
  if( opt & BMP_N )
  {
    BmpGetNChannel( bmp, opt );
  }
  if( opt & BMP_GREY )
  {
    BmpGetGreyChannel( bmp, opt );
  }
  if( opt & BMP_Z )
  {
    BmpGetZChannel( bmp, opt );
  }
  if( opt & BMP_UV )
  {
    BmpGetUVChannel( bmp, opt );
  }
  if( opt & BMP_UBYTE )
  {
    BmpGetUbyteChannel( bmp, 0, opt );
  }
  if( opt & BMP_XYZ )
  {
    for(i=0;i<3;i++) BmpGetFloatChannel0( bmp, i, opt, mmTag );
  }
  if( opt & BMP_USHORT )
  {
    szNeeded = (sx*sy)*sizeof(unsigned short);
    ALLOCMEM_TAG( bmp->dataUshort0, szNeeded, mmTag);
    if( opt & BMP_CLEAR ) memset( bmp->dataUshort0, 0, szNeeded );
  }

  if(bmp) 
  {
    BMP_MODIFIED(bmp);
  }

  UT_ASSERT_FATAL(bmp);
  return bmp;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BmpBitmap *BmpCreateBitmap2( char *inFileName, BmpBitmap *bmpSource,
    			int sx0, int sy0, int win, 
			char *name, CoConverter *colConv, 
			unsigned options )
{
  unsigned options2,channelOpt;
  int sx,sy,rcThis = 0;
  char	*pname=NULL;
  FILE	*inFile=NULL;
  BmpBitmap *bmp = NULL;

  channelOpt = options & ( BMP_RGB | BMP_GREY | BMP_ALPHA | 
      			   BMP_UV | BMP_USHORT );
  if( ! channelOpt ) options |= BMP_RGB;
  options2 = options;

#if 1
  options |= BMP_16BIT;
#endif

#if 0
  options2 = 0;
  if( options & BMP_ALPHA ) options2 |= BMP_ALPHA;
#endif

  if( name ) pname = name;

  if( inFileName )
  {
    char *suff = UtFileSuffix( inFileName );
    USE(suff);

    if(!pname) pname = inFileName;
    inFile = fopen(inFileName,"rb");
    if(!inFile) 
    {
      TRCERR(("ERROR: can't open input file '%s'\n", inFileName ));
      raiseRc( MI_ERROR );
    }
    printf("loading image '%s' ...\n",inFileName);
  }

  if( inFile )
  {
    char *suff = UtFileSuffix( inFileName );

    if( UtStrCaseCmp(suff,"bmp") == 0 )
    {
      UT_ASSERT0(FALSE);
    }
    else if( UtStrCaseCmp( suff, "u48" ) == 0 )
    {
      TRCERRR(("deprecated format U48\n"),MI_ERROR);
      #if 0
      rc = BmpReadBitmapU48( inFile, inFileName, &bmp );
      if(rc) TRCERRR(("BmpReadBitmapU48() failed %d\n"),rc);
      #endif
      sx = bmp->sx; sy = bmp->sy;
    }
    else
    {
      TRCERRR(("unknown file type for '%s'",inFileName),MI_EINVARG);
    }

    if(inFile) 
    {
      TRC1(("bit-depth=%d\n",bmp->bitDepth));
      fclose(inFile);
    }
    inFile = NULL;
  }  

  else if(!bmp)
  {
    if( bmpSource )
    {
      /*-------------------------------------------------------------------*/
      /* copy from surce bitmap						   */
      /*-------------------------------------------------------------------*/
      sx = bmpSource->sx; sy = bmpSource->sy;
      if(!bmpSource->data) 
      {
	FLAG_RESET(options2,BMP_RGB);
	if(bmpSource->dataGrey) 
	{
	  options2 |= BMP_GREY;
	}
	else
	{
	  UT_ASSERT(FALSE);
	}
      }
      bmp = BmpNewBitmap( sx, sy, options2 );
      if(!bmp) raiseRc(MI_ERROR);
      pname = bmpSource->name;
    }
    else
    {
      sx = sx0; sy = sy0;
      bmp = BmpNewBitmap( sx, sy, options2 );
      if(!bmp) raiseRc(MI_ERROR);
    } 
  }
  if(!pname) pname = "unnamed bitmap";

  strcpy( bmp->name, pname );
  if( win > 0 )
  {
#ifdef BMP_WITH_GRWIN
    rc = GWselect(win);

    bmp->GWnr = MGWmakebmp( bmp, 0, bmp->sx, bmp->sy, 24, 255.0, 
			    BMP_GEN, 0, 0, 0, NULL, NULL );
    if(!bmp->GWnr )
    TRCERRR(("MGWmakebmp(%d x %d) failed\n",bmp->sx,bmp->sy),MI_ERROR);
    rc = GWsetbmp(bmp->GWnr,0,0,1,0,1);
#endif
  }
  else
  {
    bmp->GWnr = 0;
  }

rcCatch:
  if( rcThis )
  {
    BmpCloseBitmap(bmp);
    bmp = NULL;
  }  
  if(inFile) fclose(inFile);
  return bmp;
}	
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpCreateGreyscaleBitmap( BmpBitmap *bmp )
{
  register int *pp,i,c,cr,cg,cb,rcThis=0;
  BmpGreylevel *ppg;

  if(!bmp->dataGrey)
  {
    bmp->dataGrey = 
      MM_malloc( (bmp->sx*bmp->sy+1)*sizeof(BmpGreylevel),
	  MM_DUMMYTAG );  
    if(!bmp->dataGrey ) raiseRc( MI_ENOMEM );
    /* set end-marker for fast loop access */	
    bmp->dataGrey[0] = bmp->dataGrey[bmp->sx*bmp->sy] = 0;  
  }
  for( i=0, pp= bmp->data, ppg=bmp->dataGrey; i<bmp->sx*bmp->sy;i++ )
  {
    c = *pp++;
    cr = BMP_RED(c); cg=BMP_GREEN(c); cb=BMP_BLUE(c);
    *ppg++ = RGB_TO_GREY(cr,cg,cb);
  }

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpRgbToGrey( BmpBitmap *bmpSrc, BmpBitmap *bmp, 
    		  CoConverter *colConv, unsigned opt )
{
  int	ic,icr,icg,icb,i,*pp=NULL,rcThis=0;
  float vRgb[3],vLuv[3],lum;
  BmpGreylevel *ppg=NULL;

  if(!bmp || !bmpSrc) TRCERRR(("unexpected NULL pointer"),MI_ERROR);

  if( ! BMP_DIM_EQ( bmpSrc, bmp ) )
    TRCERRR(("bitmap dimension differ"),MI_ERROR);

  pp = bmpSrc->data;
  
  ppg = BmpGetGreyChannel( bmp, 0 );
  if(!ppg) TRCERRR(("unexpected NULL pointer"),MI_ERROR);

  bmp->colConv = colConv;

  CLOCK_START();

  for(i=0;i<bmp->sx*bmp->sy;i++, pp++)
  {
    BMP_PP_GETPIXEL( pp, ic, icr, icg, icb );
    CO_INT2RGB( vRgb, icr, icg, icb );
    colConv->fromRgb( vLuv, vRgb  ); 
    lum = vLuv[0]*255; if(lum>255) lum=255; if(lum<0) lum=0;
    *ppg++ = lum;
  }
  
  CLOCK_STOP(); 
  TRC(("BmpRgbToLuv: CPU = %d ms (%d pixels)\n", 
    (int)CLOCK_DIFF_MS(),bmp->sx*bmp->sy));

rcCatch:
  return rcThis;

}

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpGreyToRgb( BmpBitmap *bmpSrc, BmpBitmap *bmp, 
    		  CoConverter *colConv, unsigned opt )
{
  int	icr,icg,icb,i,*pp=NULL,rcThis=0;
  float vRgb[3],vLuv[3],lum;
  BmpGreylevel *ppg=NULL;

  if(!bmp || !bmpSrc) TRCERRR(("unexpected NULL pointer"),MI_ERROR);

  if( ! BMP_DIM_EQ( bmpSrc, bmp ) )
    TRCERRR(("bitmap dimension differ"),MI_ERROR);

  ppg = bmpSrc->dataGrey;
  if(!(ppg)) TRCERRR(("no Grey channel"),MI_ERROR); 
  
  pp = bmp->data;
  if(!pp) TRCERRR(("unexpected NULL pointer"),MI_ERROR);

  CLOCK_START();

  for(i=0;i<bmp->sx*bmp->sy;i++)
  {
    lum = *ppg++;
    vLuv[0] = MIN(lum/255.0,255.0);
    vLuv[1] = vLuv[2] = 0;
    colConv->toRgb( vRgb, vLuv  ); 
    CO_RGB2INT( vRgb, icr,icg,icb );
/*    printf("%d %f %f %f\n",icuv,vLuv[0],vLuv[1],vLuv[2]);*/
    *pp++ = BMP_MKRGB(icr,icg,icb);    
  }

  BMP_MODIFIED( bmp );
  
  CLOCK_STOP(); 
  TRC(("BmpRgbToLuv: CPU = %d ms (%d pixels)\n", 
    (int)CLOCK_DIFF_MS(),bmp->sx*bmp->sy));

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpRgbToLuv_orig( BmpBitmap *bmpSrc, BmpBitmap *bmp, 
    		 CoConverter *colConv, unsigned opt )
{
  int	sx,sy,doProgr,i,j,*pp=NULL,rcThis=0;
  BmpGreylevel *ppl=NULL;
  BmpUV *ppuv=NULL;
  UtProgrInd pgi;

  if(!bmp || !bmpSrc) TRCERRR(("unexpected NULL pointer"),MI_ERROR);

  if( ! BMP_DIM_EQ( bmpSrc, bmp ) )
    TRCERRR(("bitmap dimension differ"),MI_ERROR);

  pp = bmpSrc->data;
  
  if( opt & BMP_GREY )
  {
    ppl = BmpGetGreyChannel( bmp, 0 );
    if(!ppl) TRCERRR(("unexpected NULL pointer"),MI_ERROR);
  }

  if( opt & BMP_UV )
  {
    ppuv = BmpGetUVChannel( bmp, 0 );
    if(!ppuv) TRCERRR(("unexpected NULL pointer"),MI_ERROR);
  }

  bmp->colConv = colConv;
  
  CLOCK_START();
  sx = bmp->sx; sy = bmp->sy;
  UT_PROGRIND_INIT( &pgi );
  doProgr = sy*sx > 500000;

  #pragma omp parallel for private(i,j)    
  for(i=0;i<sy;i++)
  {
    /*if(doProgr) 
    {
      UT_PROGRIND( &pgi, i, 0, bmp->sy, "completed: ", NULL );
    }*/
    for(j=0;j<sx;j++)
    {
      int ic,icr,icg,icb;
      float vRgb[3],vLuv[3],lum;
      LongInt ii=BMP_XYOFF(bmp,j,i);
      BMP_PP_GETPIXEL( pp+ii, ic, icr, icg, icb );
      CO_INT2RGB( vRgb, icr, icg, icb );
      colConv->fromRgb( vLuv, vRgb  ); 
      if( ppl ) 
      {
	lum = vLuv[0]*255; if(lum>255) lum=255; if(lum<0) lum=0;
	*(ppl+ii) = lum;
      }
      if( ppuv )
      {
	BMP_MKUV((ppuv+ii),vLuv);
      }
    }
  }

  BMP_MODIFIED( bmp );

  CLOCK_STOP(); 
  TRC(("BmpRgbToLuv: CPU = %d ms (%d pixels)\n", 
    (int)CLOCK_DIFF_MS(),bmp->sx*bmp->sy));

rcCatch:
  return rcThis;
}

#if 0
/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
typedef struct
{
  unsigned opt;
  BmpBitmap *bmp,*bmpSrc;
  CoConverter *colConv;
}
BmpRgbToLuv_Ctx;

typedef struct
{
  int y1,y2;
}
BmpRgbToLuv_Task;

int BmpRgbToLuv_SMP_2( BmpBitmap *bmpSrc, BmpBitmap *bmp, 
    		 CoConverter *colConv, unsigned opt )
{
  int	rc,sx,sy,doProgr,ic,icr,icg,icb,i,j,iu,iv,*pp=NULL,rcThis=0;
  float vRgb[3],vLuv[3],lum;
  BmpGreylevel *ppl=NULL;
  BmpUV *ppuv=NULL;
  UtProgrInd pgi;
  BmpRgbToLuv_Task *tsk;
  BmpRgbToLuv_Ctx context;
  UtTimer   tm;

  if(!bmp || !bmpSrc) TRCERRR(("unexpected NULL pointer"),MI_ERROR);

  if( ! BMP_DIM_EQ( bmpSrc, bmp ) )
    TRCERRR(("bitmap dimension differ"),MI_ERROR);

  pp = bmpSrc->data;
  
  if( opt & BMP_GREY )
  {
    ppl = BmpGetGreyChannel( bmp, 0 );
    if(!ppl) TRCERRR(("unexpected NULL pointer"),MI_ERROR);
  }

  if( opt & BMP_UV )
  {
    ppuv = BmpGetUVChannel( bmp, 0 );
    if(!ppuv) TRCERRR(("unexpected NULL pointer"),MI_ERROR);
  }

  bmp->colConv = colConv;
  
  TIMER_START(&tm);
  sx = bmp->sx; sy = bmp->sy;
  UT_PROGRIND_INIT( &pgi );
  doProgr = sy*sx > 500000;

  context.opt = opt;
  context.colConv = colConv;
  context.bmpSrc = bmpSrc;
  context.bmp = bmp;

  SmpInitPool();

  SMP_TASK_LOOP( 0, sy, 
      		 &context, 
		 BmpIRgbToLuv_Range_2, 
		 tsk )
  {
    tsk->y1 = SMP_TASK_I0(tsk); 
    tsk->y2 = tsk->y1+SMP_TASK_N(tsk)-1;

    SMP_SUBMIT( tsk );
  }

  SMP_RESULT_LOOP(rc);
  if(rc) TRCERRR(("SMP_RESULT_ALL() failed %d\n",rc),MI_ERROR);

  BMP_MODIFIED( bmp );

  TIMER_STOP(&tm);
  TRC1(("BmpRgbToLuv: CPU = %d ms (%d pixels)\n", 
    (int)TIMER_DIFF_MS(&tm),bmp->sx*bmp->sy));

rcCatch:
  return rcThis;
}
#endif

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
typedef struct 
{
  int y1,y2;
  unsigned opt;
  BmpBitmap *bmp,*bmpSrc;
  CoConverter *colConv;
}
BmpIRgbToLuv_Task;

int BmpIRgbToLuv_Range( void *t0 )
{
  BmpIRgbToLuv_Task *t=(BmpIRgbToLuv_Task*)t0;
  int	i,ic,icr,icg,icb,j,*pp=NULL;
  float vRgb[3],vLuv[3],lum;
  BmpGreylevel *ppl=NULL;
  BmpUV *ppuv=NULL;
  size_t offs=t->y1*t->bmp->sx;
  BmpBitmap *bmp=t->bmp,*bmpSrc=t->bmpSrc;

  pp = bmpSrc->data;  
  if( t->opt & BMP_GREY )
  {
    ppl = BmpGetGreyChannel( bmp, 0 );
  }
  if( t->opt & BMP_UV )
  {
    ppuv = BmpGetUVChannel( bmp, 0 );
  }

  if(ppuv) ppuv += offs;
  if(ppl) ppl +=offs;
  if(pp) pp +=offs;
  
  for(i=t->y1;i<=t->y2;i++)
  {
    for(j=0;j<bmp->sx;j++,pp++)
    {
      BMP_PP_GETPIXEL( pp, ic, icr, icg, icb );
      CO_INT2RGB( vRgb, icr, icg, icb );
      t->colConv->fromRgb( vLuv, vRgb  ); 
      if( ppl ) 
      {
	lum = vLuv[0]*255; if(lum>255) lum=255; if(lum<0) lum=0;
	*ppl++ = lum;
      }
      if( ppuv )
      {
        BMP_MKUV(ppuv,vLuv);
	ppuv++;
      }
    }
  }


  return 0;
}

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpRgbToLuv( BmpBitmap *bmpSrc, BmpBitmap *bmp, 
    		     CoConverter *colConv, unsigned opt )
{
  return BmpRgbToLuv_orig( bmpSrc, bmp, colConv, opt );
#if 0
  int	sx,sy,doProgr,i,*pp=NULL,rcThis=0,rc;
  BmpGreylevel *ppl=NULL;
  BmpUV *ppuv=NULL;
  UtProgrInd pgi;
  UtTimer   tm;
  BmpIRgbToLuv_Task  *tsk1,*tsk2;

  if(!bmp || !bmpSrc) TRCERRR(("unexpected NULL pointer"),MI_ERROR);

  if( ! BMP_DIM_EQ( bmpSrc, bmp ) )
    TRCERRR(("bitmap dimension differ"),MI_ERROR);

  pp = bmpSrc->data;
  
  if( opt & BMP_GREY )
  {
    ppl = BmpGetGreyChannel( bmp, 0 );
    if(!ppl) TRCERRR(("unexpected NULL pointer"),MI_ERROR);
  }

  if( opt & BMP_UV )
  {
    ppuv = BmpGetUVChannel( bmp, 0 );
    if(!ppuv) TRCERRR(("unexpected NULL pointer"),MI_ERROR);
  }

  bmp->colConv = colConv;
  
  TIMER_START(&tm);
  sx = bmp->sx; sy = bmp->sy;
  UT_PROGRIND_INIT( &pgi );
  doProgr = sy*sx > 500000;

  SmpInitPool();

  for(i=0;i<sy;i+=2)
  {
    if(doProgr) 
    {
      UT_PROGRIND( &pgi, i, 0, bmp->sy, "completed: ", NULL );
    }
#if 0
    {
      int dbg=1;
      while(dbg)
      {
        ;
      }
    }
#endif
    SMP_TASK(NULL,NULL,BmpIRgbToLuv_Range,tsk1);
    tsk1->y1=i; tsk1->y2=i;
    tsk1->bmp = bmp; tsk1->colConv = colConv;
    tsk1->bmpSrc = bmpSrc; tsk1->opt = opt;
    SMP_SUBMIT(tsk1);

    if(i<sy-1)
    {
      SMP_TASK(NULL,NULL,BmpIRgbToLuv_Range,tsk2);
      tsk2->y1=i+1; tsk2->y2=i+1;    
      tsk2->bmp = bmp; tsk2->colConv = colConv;
      tsk2->bmpSrc = bmpSrc; tsk2->opt = opt;
      SMP_SUBMIT(tsk2);
    } else tsk2=NULL;

    SMP_RESULT( tsk1, rc ); if(rc) raiseRc(rc);
    SMP_FREE( tsk1 );

    if(tsk2)
    {
      SMP_RESULT( tsk2, rc ); if(rc) raiseRc(rc);
      SMP_FREE( tsk2 );
    }
  }

  BMP_MODIFIED( bmp );

  TIMER_STOP(&tm);
  TRC1(("BmpRgbToLuv: CPU = %d ms (%d pixels)\n", 
    (int)TIMER_DIFF_MS(&tm),bmp->sx*bmp->sy));

rcCatch:
  return rcThis;
#endif
}

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpLuvToRgb( BmpBitmap *bmpSrc, BmpBitmap *bmp, 
    		 CoConverter *colConv, unsigned opt )
{
  int	j,icr,icg,icb,i,*pp=NULL,rcThis=0,doProgr;
  float vRgb[3],vLuv[3],lum;
  BmpGreylevel *ppl=NULL;
  BmpUV *ppuv=NULL;
  UtProgrInd pgi;

  if(!bmp || !bmpSrc) TRCERRR(("unexpected NULL pointer"),MI_ERROR);

  if( ! BMP_DIM_EQ( bmpSrc, bmp ) )
    TRCERRR(("bitmap dimension differ"),MI_ERROR);

  bmpSrc->colConv = colConv;

  ppuv = bmpSrc->dataUV;
  ppl = bmpSrc->dataGrey;

  if(!(ppuv||ppl)) TRCERRR(("no LUV channel"),MI_ERROR); 
  
  pp = bmp->data;
  if(!pp) 
  {
    BmpGetRGBChannel(bmp,0);
    pp = bmp->data;
  }
  if(!pp)  TRCERRR(("unexpected NULL pointer"),MI_ERROR);

  UT_PROGRIND_INIT( &pgi );
  doProgr = bmp->sy*bmp->sx > 500000;

  CLOCK_START();
  for(i=0;i<bmp->sy;i++)
  {
    if(doProgr) 
    {
      UT_PROGRIND( &pgi, i, 0, bmp->sy, "LUV -> RGB: ", NULL );
    }
    for(j=0;j<bmp->sx;j++)
    {
      if(ppuv)
      {
        BMP_PP_GETPIXEL_UV(ppuv,vLuv);
      }
      else
      {
	MT_VEC3_SET(vLuv,0,0,0);
      }
      lum = *ppl++;
      vLuv[0] = MIN(lum/255.0,255.0);
      colConv->toRgb( vRgb, vLuv  ); 
      CO_RGB2INT( vRgb, icr,icg,icb );
      /*    printf("%d %f %f %f\n",icuv,vLuv[0],vLuv[1],vLuv[2]);*/
      *pp++ = BMP_MKRGB(icr,icg,icb);    
    }
    if(ppuv) ppuv++;
  }

  BMP_MODIFIED( bmp );

  CLOCK_STOP(); 
  TRC(("BmpLuvToRgb: CPU = %d ms (%d pixels)\n", 
    (int)CLOCK_DIFF_MS(),bmp->sx*bmp->sy));

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpGreyToRgbFast( BmpBitmap *bmpSrc, BmpBitmap *bmp )
{
  int	icr,icg,icb,i,*pp=NULL,rcThis=0;
  float lum;
  BmpGreylevel *ppl=NULL;

  if(!bmp || !bmpSrc) TRCERRR(("unexpected NULL pointer"),MI_ERROR);

  if( ! BMP_DIM_EQ( bmpSrc, bmp ) )
    TRCERRR(("bitmap dimension differ"),MI_ERROR);

  bmpSrc->colConv = NULL;

  ppl = bmpSrc->dataGrey;
  if(!(ppl)) TRCERRR(("no Grey channel"),MI_ERROR); 
  
  pp = BmpGetRGBChannel( bmp, 0 );
  if(!pp) TRCERRR(("unexpected NULL pointer"),MI_ERROR);

  for(i=0;i<bmp->sx*bmp->sy;i++)
  {
      lum = *ppl++;
      if(lum<0) lum = 0; 
      if(lum>255) lum=255;
      icr = lum; icg = lum; icb = lum;
      *pp++ = BMP_MKRGB(icr,icg,icb);    
  }

  BMP_MODIFIED( bmp );

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpSaturationGammaCorrect( 
    	BmpBitmap *bs, BmpBitmap *bt, BmpRectangle *win,
	float gamma, float alpha, CoConverter *cc, unsigned opt )
{
  int rcThis=0,x,y,*pps,*ppt,x0,y0,sx,sy,rs,rt,ic,icr,icg,icb;
  
  float vRgb[3],vLuv[3],s,u,v,s2,u0,v0;
  UtProgrInd pgi;

  BMP_GET_WIN( bs, win, x0, y0, sx, sy );

  BMP_RSCAN_INIT( bs, bs->data, x0, y0, sx, rs, pps ); 
  BMP_RSCAN_INIT( bt, bt->data, x0, y0, sx, rt, ppt ); 
  UT_PROGRIND_INIT( &pgi );

  /*** get reference white point */
  MT_VEC3_SET(vRgb,1,1,1);
  cc->fromRgb(vLuv,vRgb);
  u0 = vLuv[1]; v0=vLuv[2];

  TRC(("white point = %f %f\n",(float)u0,(float)v0));

  BMP_RSCAN_LOOP_Y( y,sy, (pps+=rs, ppt+=rt))
  {
    UT_PROGRIND( &pgi, y, 0, sy-1, "completed: ", NULL );

    for(x=0;x<sx;x++,pps++,ppt++)
    {
      BMP_PP_GETPIXEL( pps, ic, icr, icg, icb );
      CO_INT2RGB( vRgb, icr, icg, icb );

#if 0
      CoConvHCL->fromRgb( vLuv, vRgb  ); 
      s = MAX( vLuv[2]+1, 0) ;
      vLuv[2] = pow(s,gamma)-1;
      CoConvHCL->toRgb( vRgb, vLuv  ); 
#endif

      cc->fromRgb( vLuv, vRgb  ); 

      u = vLuv[1]-u0; v = vLuv[2]-v0;
      s = sqrt(u*u+v*v);
      if( s > MT_REAL_EPSILON )
      {
        s2 = alpha * pow(s,gamma);
	vLuv[1] = u0+u * s2/s; vLuv[2] = v0+v * s2/s;
      }
      
      cc->toRgb( vRgb, vLuv  ); 

      CO_RGB2INT( vRgb, icr,icg,icb );
      *ppt = BMP_MKRGB(icr,icg,icb);    
    }
  }
  BMP_MODIFIED( bt );

/*rcCatch:*/
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpSaturationLogCorrect( 
    	BmpBitmap *bs, BmpBitmap *bt, BmpRectangle *win,
	float a, int nIter, float gamma, 
	CoConverter *cc, unsigned opt )
{
  int rcThis=0,x,y,*pps,*ppt,x0,y0,sx,sy,rs,rt,ic,icr,icg,icb,i,
  	doGamma=FALSE;
  
  float vRgb[3],vLuv[3],s,u,v,s2,u0,v0;
  UtProgrInd pgi;

  BMP_GET_WIN( bs, win, x0, y0, sx, sy );

  BMP_RSCAN_INIT( bs, bs->data, x0, y0, sx, rs, pps ); 
  BMP_RSCAN_INIT( bt, bt->data, x0, y0, sx, rt, ppt ); 
  UT_PROGRIND_INIT( &pgi );

  /*** get reference white point */
  MT_VEC3_SET(vRgb,1,1,1);
  cc->fromRgb(vLuv,vRgb);
  u0 = vLuv[1]; v0=vLuv[2];

  if(!nIter) nIter =1;
  doGamma = fabs(gamma-1.0)>MT_REAL_EPSILON;

  TRC(("white point = %f %f\n",(float)u0,(float)v0));

  BMP_RSCAN_LOOP_Y( y,sy, (pps+=rs, ppt+=rt))
  {
    UT_PROGRIND( &pgi, y, 0, sy-1, "completed: ", NULL );

    for(x=0;x<sx;x++,pps++,ppt++)
    {
      BMP_PP_GETPIXEL( pps, ic, icr, icg, icb );
      CO_INT2RGB( vRgb, icr, icg, icb );

#if 0
      CoConvHCL->fromRgb( vLuv, vRgb  ); 
      s = MAX( vLuv[2]+1, 0) ;
      vLuv[2] = pow(s,gamma)-1;
      CoConvHCL->toRgb( vRgb, vLuv  ); 
#endif

      cc->fromRgb( vLuv, vRgb  ); 

      u = vLuv[1]-u0; v = vLuv[2]-v0;

      if( doGamma )
      {
	s = sqrt(u*u+v*v);
	if( s > MT_REAL_EPSILON )
	{
	  s2 = pow(s,gamma);
	  u = u * s2/s; v = v * s2/s;
	}
      }

      for(i=0;i<nIter;i++)
      {
	s = sqrt(u*u+v*v);
	if( s > MT_REAL_EPSILON )
	{
	  s2 = log(1+s/a)*a;
	  u = u * s2/s; v = v * s2/s;
	}
      }

      vLuv[1] = u0+u; vLuv[2] = v0 + v;
      
      cc->toRgb( vRgb, vLuv  ); 

      CO_RGB2INT( vRgb, icr,icg,icb );
      *ppt = BMP_MKRGB(icr,icg,icb);    
    }
  }
  BMP_MODIFIED( bt );

/*rcCatch:*/
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpConvert16BitGrey( BmpBitmap *bs,
		         float alpha, float gamma,
			 float beta,
    			 CoConverter *cc, unsigned opt )
/* y = ( alpha * (x+beta) ) ^ gamma */
{
  int rcThis=0,x,y,x0,y0,sx,sy,rs,*pp,cr,cg,cb;
  BmpGreylevel *ppg;
  float g,gOut;
  BmpRectangle *win=NULL;

  BMP_GET_WIN( bs, win, x0, y0, sx, sy );

  if(!bs->dataGrey)
  {
    TRCERRR(("no grey channel for bitmap"),MI_ERROR);
  }

  pp = BmpGetRGBChannel( bs, 0 );

  BMP_RSCAN_INIT( bs, bs->dataGrey, x0, y0, sx, rs, ppg ); 

  BMP_RSCAN_LOOP( ppg, rs, x, y, sx, sy)
  {
    g = *ppg + beta;
    if(g<0) g = 0;
    gOut = pow(g*alpha,gamma);
    if(gOut>255) gOut = 255;
    if(gOut<0) gOut=0;
    *ppg = gOut;
    cr = cg = cb = gOut;
    *pp++ = BMP_MKRGB( cr,cg,cb );
  }
  BMP_MODIFIED( bs );

rcCatch:
  return rcThis;
}


#if 0
/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpColorClip( 
    	BmpBitmap *bs, BmpBitmap *bt, BmpRectangle *win,
	float	sMin, sMax, /* saturation */
	float   hMin, hMax, /* hue */
	float a, int nIter, float gamma, 
	CoConverter *cc, unsigned opt )
{
  int rcThis=0,x,y,*pps,*ppt,x0,y0,sx,sy,rs,rt,ic,icr,icg,icb,i,
  	doGamma=FALSE;
  
  float vRgb[3],vLuv[3],s,u,v,s2,u0,v0,s0=sMin,df=0.03,f,h,h0=hMin;
  UtProgrInd pgi;

  BMP_GET_WIN( bs, win, x0, y0, sx, sy );

  BMP_RSCAN_INIT( bs, bs->data, x0, y0, sx, rs, pps ); 
  BMP_RSCAN_INIT( bt, bt->data, x0, y0, sx, rt, ppt ); 
  UT_PROGRIND_INIT( &pgi );

  /*** get reference white point */
  MT_VEC3_SET(vRgb,1,1,1);
  cc->fromRgb(vLuv,vRgb);
  u0 = vLuv[1]; v0=vLuv[2];

  if(!nIter) nIter =1;
  doGamma = fabs(gamma-1.0)>MT_REAL_EPSILON;

  TRC(("white point = %f %f\n",(float)u0,(float)v0));

  BMP_RSCAN_LOOP_Y( y,sy, (pps+=rs, ppt+=rt))
  {
    UT_PROGRIND( &pgi, y, 0, sy-1, "completed: ", NULL );

    for(x=0;x<sx;x++,pps++,ppt++)
    {
      BMP_PP_GETPIXEL( pps, ic, icr, icg, icb );
      CO_INT2RGB( vRgb, icr, icg, icb );

      CoConvHCL->fromRgb( vLuv, vRgb  ); 
      s = MAX( vLuv[2]+1, 0) ;
      
      /* saturation dependent clipping strength */
      if( s < s0 - df )
	f = 0;
      else if( s < s0 + df )
	f = (s-s0)/(2*df) + 0.5;
      else
	f = 1;

      /* clip hue */
      h = vLuv[1];
      
      if( f > MT_REAL_EPSILON )
      {
        MtClipInterval(h,h0,h0+(1-f),1);
      }

      vLuv[1] = h;
      
      CoClipLuvToRgbChk( CoConvHCL, vLuv[0], vLuv, vRgb, &icr, &icg, &icb );
      *ppt = BMP_MKRGB(icr,icg,icb);    
    }
  }
  BMP_MODIFIED( bt );

/*rcCatch:*/
  return rcThis;
}
#endif 

#if 0
/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpUVLintrans( BmpBitmap *bs, BmpBitmap *bt, BmpRectangle *win,
 		   float *trans2D, CoConverter *cc, unsigned opt )
{
  int rcThis,x,y,*pps,*ppt,x0,y0,sx,sy,rs,rt,ir,ig,ib,icu,icv,
  		doDstUV=FALSE;
  float f,vRgb[3],vLuv[3], *t=trans2D, vTmp[3];
  UtProgrInd pgi;
  BmpUV *pptUV = NULL, *ppsUV = NULL;
  BmpGreylevel *pptG = NULL, *ppsG = NULL;

  BMP_GET_WIN( bs, win, x0, y0, sx, sy );

  if( doDstUV )
  {
    pptUV = BmpGetUVChannel( bt, 0 );
    pptG = BmpGetGreyChannel( bt, 0 );
  }

  ppsUV = bs->dataUV;
  ppsG = bs->dataGrey;
  if( !(ppsUV && ppsG ) )
  {
    TRCERRR(("no LUV source channels \n"),MI_ERROR);
  }

  BMP_RSCAN_INIT( bs, bs->data, x0, y0, sx, rs, pps ); 
  BMP_RSCAN_INIT( bt, bt->data, x0, y0, sx, rt, ppt ); 
  UT_PROGRIND_INIT( &pgi );

  BMP_RSCAN_LOOP_Y( y,sy, (pps+=rs, ppt+=rt))
  {
    UT_PROGRIND( &pgi, y, 0, sy-1, "completed: ", NULL );

    for(x=0;x<sx;x++,pps++,ppt++)
    {
	/* get source LUV */
        BMP_PP_GETPIXEL_LUV(ppsG,ppsUV,vLuv);
	ppsG++; ppsUV++;

        /* transform UV chroma values */
        


	/* transform luminance */
	f = a[0]*vLuv[0] + b[0];
	if(f<0) f = 0; if(f>=1) f = 1;
	vLuv[0] = f;

	/* rotate UV components */
	if( (cc != &CoConvHCL) && doRot )
	{
	  MT_VEC2_ROTATE(vTmp, vLuv+1, rsin, rcos );
	  vLuv[1] = vTmp[0]; vLuv[2] = vTmp[1];
	}

	/* transform chroma channel U */
	f = a[1]*vLuv[1] + b[1];
	if( cc == &CoConvHCL )
	{
	  /*** HUE => rotate modulo [-1,1] */
	  if( f<=-1 || f>=1 ) f = CO_HUE_MOD( f );
	}
	else
	{
	  if(f<=-1) f = -1; if(f>=1) f = 1;
	}
	vLuv[1] = f;

	/* transform chroma channel V */
	f = a[2]*vLuv[2] + b[2];
	if(f<=-1) f = -1; if(f>=1) f = 1;
	vLuv[2] = f;

	/* write dest. LUV */
	if( doDstUV )
	{
	  CO_UV2INT( vLuv, icu, icv );
	  *pptUV++ = BMP_MKUV( icu, icv );
	  *pptG++ = vLuv[0] * 255.;
	}

        /* LUV -> RGB */ 
        cc->toRgb( vRgb, vLuv  ); 
        CO_RGB2INT( vRgb, ir, ig, ib );
        *ppt = BMP_MKRGB( ir, ig, ib );
    }
    if( doDstUV )
    {
      pptUV += rt; pptG += rt;
      ppsUV += rs; ppsG += rs;
    }
  }
  BMP_MODIFIED( bt );

rcCatch:
  return rcThis;
}

#endif

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpLocalColorStretch( 
    	BmpBitmap *bs, BmpBitmap *bt, BmpRectangle *win,
	int r, float s,
	CoConverter *cc, unsigned opt )
{
  int rcThis=0;
#if 0
  int rcThis=0,x,y,*pps,*ppt,x0,y0,sx,sy,rs,rt,ic,icr,icg,icb;
  
  float vRgb[3],vLuv[3],s,u,v,s2,u0,v0;
  UtProgrInd pgi;

  BMP_GET_WIN( bs, win, x0, y0, sx, sy );

  BMP_RSCAN_INIT( bs, bs->data, x0, y0, sx, rs, pps ); 
  BMP_RSCAN_INIT( bt, bt->data, x0, y0, sx, rt, ppt ); 
  UT_PROGRIND_INIT( &pgi );

  /*** get reference white point */
  MT_VEC3_SET(vRgb,1,1,1);
  cc->fromRgb(vLuv,vRgb);
  u0 = vLuv[1]; v0=vLuv[2];


  TRC(("white point = %f %f\n",(float)u0,(float)v0));

  BMP_RSCAN_LOOP_Y( y,sy, (pps+=rs, ppt+=rt))
  {
    UT_PROGRIND( &pgi, y, 0, sy-1, "completed: ", NULL );

    for(x=0;x<sx;x++,pps++,ppt++)
    {
      BMP_PP_GETPIXEL( pps, ic, icr, icg, icb );
      CO_INT2RGB( vRgb, icr, icg, icb );

      cc->fromRgb( vLuv, vRgb  ); 
      if( cc == &CoConvHCL )
      {	
        s = MAX( vLuv[2]+1, 0) ;
      }
      else
      {
        u = vLuv[1]-u0; v = vLuv[2]-v0;
	s = sqrt(u*u+v*v);

	if( s > MT_REAL_EPSILON )
	{
/*	  s2 = pow(s,gamma);*/
	  s2 = MtClipInterval(s,a,b,c,0);
	  u = u * s2/s; v = v * s2/s;
	}

        vLuv[1] = u0+u; vLuv[2] = v0 + v;
      }

      cc->toRgb( vRgb, vLuv  ); 

      CO_RGB2INT( vRgb, icr,icg,icb );
      *ppt = BMP_MKRGB(icr,icg,icb);    
    }
  }
  BMP_MODIFIED( bt );
#endif
/*rcCatch:*/
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BmpAlphaVal *BmpGetAlphaChannel( BmpBitmap *bmp, unsigned opt )
{
  size_t szNeeded = (bmp->sx*bmp->sy)*sizeof(BmpAlphaVal);
  if(!bmp->dataAlpha)
  {
    ALLOCMEM( bmp->dataAlpha, szNeeded);  
  }
  if( opt & BMP_CLEAR ) 
  {
    memset( bmp->dataAlpha, 0, szNeeded );
  }

  return bmp->dataAlpha;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void BmpFreeValidChannel( BmpBitmap *bmp )
{
  if(bmp)
  {
    FREEMEM(bmp->dataValid);
  }
}

void BmpFreeAlphaChannel( BmpBitmap *bmp )
{
  if(bmp)
  {
    FREEMEM(bmp->dataAlpha);
  }
}

void BmpFreeGreyChannel( BmpBitmap *bmp )
{
  if(bmp)
  {
    FREEMEM(bmp->dataGrey);
  }
}

void BmpFreeZChannel( BmpBitmap *bmp )
{
  if(bmp)
  {
    FREEMEM(bmp->dataZ);
  }
}

void BmpFreeFloatChannel( BmpBitmap *bmp, int i )
{
  UT_ASSERT0(i>=0&&i<BMP_MAX_FCHAN);
  if(bmp)
  {
    FREEMEM(bmp->dataFloat[i]);
  }
}

void BmpFreeRgbChannel( BmpBitmap *bmp )
{
  if(bmp)
  {
    FREEMEM(bmp->data);
  }
}

void BmpFreeDblChannel( BmpBitmap *bmp, int i )
{
  if(bmp)
  {
    FREEMEM(bmp->dataDbl[i]);
  }
}

void BmpFreeUshortChannel( BmpBitmap *bmp, int i )
{
  UT_ASSERT0(i>=0&&i<BMP_MAX_UCHAN);
  if(bmp)
  {
    FREEMEM(bmp->dataUshort[i]);
  }
}

void BmpFreeUbyteChannel( BmpBitmap *bmp, int i )
{
  UT_ASSERT0(i>=0&&i<BMP_MAX_UCHAN);
  if(bmp)
  {
    FREEMEM(bmp->dataUbyte[i]);
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BmpGreylevel *BmpGetGreyChannel( BmpBitmap *bmp, unsigned opt )
{
  size_t szNeeded = (bmp->sx*bmp->sy)*sizeof(BmpGreylevel);
  if(!bmp->dataGrey)
  {
    ALLOCMEM( bmp->dataGrey, szNeeded);  
    bmp->dataGreyIsExt = FALSE;
  }
  if( opt & BMP_CLEAR ) 
  {
    memset( bmp->dataGrey, 0, szNeeded );
  }

  return bmp->dataGrey;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static void BmpIFreeGreyChannel( BmpBitmap *bmp )
{
  if(!bmp->dataGreyIsExt)
  {
    FREEMEM(bmp->dataGrey);
  }
  bmp->dataGreyIsExt=FALSE;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BmpGreylevel *BmpSetGreyChannel( BmpBitmap *bmp, BmpGreylevel *dataGrey,
    				 unsigned opt )
{
  size_t szNeeded = (bmp->sx*bmp->sy)*sizeof(BmpGreylevel);
  
  BmpIFreeGreyChannel(bmp);

  bmp->dataGrey = dataGrey;
  
  if(opt & BMP_EXTERN) bmp->dataGreyIsExt = TRUE;

  if( opt & BMP_CLEAR ) 
  {
    memset( bmp->dataGrey, 0, szNeeded );
  }

  return bmp->dataGrey;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void **BmpGetPtrChannel( BmpBitmap *bmp, unsigned opt )
{
  size_t szNeeded = (bmp->sx*bmp->sy)*sizeof(void *);
  if(!bmp->dataPtr)
  {
    ALLOCMEM( bmp->dataPtr, szNeeded);  
  }
  if( opt & BMP_CLEAR ) 
  {
    memset( bmp->dataPtr, 0, szNeeded );
  }

  return bmp->dataPtr;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int *BmpGetRGBChannel( BmpBitmap *bmp, unsigned opt )
{
  size_t szNeeded = (bmp->sx*bmp->sy)*sizeof(int);
  if(!bmp->data)
  {
    ALLOCMEM( bmp->data, szNeeded);  
  }
  if( opt & BMP_CLEAR ) 
  {
    memset( bmp->data, 0, szNeeded );
  }

  return bmp->data;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BmpUV *BmpGetUVChannel( BmpBitmap *bmp, unsigned opt )
{
  size_t szNeeded = (bmp->sx*bmp->sy)*sizeof(BmpUV);
  if(!bmp->dataUV)
  {
    ALLOCMEM( bmp->dataUV, szNeeded);  
  }
  if( opt & BMP_CLEAR ) 
  {
    memset( bmp->dataUV, 0, szNeeded );
  }

  return bmp->dataUV;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BmpZval *BmpGetZbufChannel( BmpBitmap *bmp, unsigned opt )
{
  size_t szNeeded = (bmp->sx*bmp->sy)*sizeof(BmpZval);
  if(!bmp->dataZbuf)
  {
    ALLOCMEM( bmp->dataZbuf, szNeeded);  
  }
  if( opt & BMP_CLEAR ) 
  {
    memset( bmp->dataZbuf, 0, szNeeded );
  }

  return bmp->dataZbuf;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
double *BmpGetDblChannel( BmpBitmap *bmp, int i, unsigned opt )
{
  int rcThis=0;
  size_t szNeeded = (bmp->sx*bmp->sy)*sizeof(double);

  UT_ASSERT(i<BMP_MAX_DCHAN);
  if(!bmp->dataDbl[i])
  {
    ALLOCMEM( bmp->dataDbl[i], szNeeded);  
  }
  if( opt & BMP_CLEAR ) 
  {
    LongInt ii;
    for(ii=0;ii<bmp->sx*bmp->sy;ii++) bmp->dataDbl[i][ii]=0;
  }
rcCatch:
  return rcThis==0 ? bmp->dataDbl[i] : NULL;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
float *BmpGetFloatChannel0( BmpBitmap *bmp, int i, unsigned opt,
    			    char *mmTag )
{
  int rcThis=0;
  size_t szNeeded = (bmp->sx*bmp->sy)*sizeof(float);

  UT_ASSERT(i<BMP_MAX_FCHAN);
  if(!bmp->dataFloat[i])
  {
    ALLOCMEM_TAG( bmp->dataFloat[i], szNeeded, mmTag);  
  }
  if( opt & BMP_CLEAR ) 
  {
    memset( bmp->dataFloat[i], 0, szNeeded );
  }
rcCatch:
  return rcThis==0 ? bmp->dataFloat[i] : NULL;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpGetFloatChannels( BmpBitmap *bmp, int n, unsigned opt )
{
  int rcThis=0,i;
  UT_ASSERT(n<=BMP_MAX_FCHAN);
  for(i=0;i<n;i++) 
  {
    float *fp=BmpGetFloatChannel(bmp,i,opt);
    if(!fp) raiseRc(MI_ENOMEM);
  }
  bmp->nFloatChan = n;
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpGetXYZChannel( BmpBitmap *bmp, unsigned opt )
{
  int rcThis=0,k;
  for(k=0;k<3;k++) 
  {
    float *pf=BmpGetFloatChannel(bmp,k,opt);
    if(!pf) raiseRc(MI_ENOMEM);
  }

rcCatch:
  return rcThis;
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
uint16_t *BmpGetUshortChannel( BmpBitmap *bmp, int i, unsigned opt )
{
  int rcThis=0;
  size_t szNeeded = (bmp->sx*bmp->sy)*sizeof(bmp->dataUshort[0][0]);

  UT_ASSERT(sizeof(bmp->dataUshort[0][0])==2);
  UT_ASSERT(sizeof(bmp->dataUshort[0][0])==sizeof(uint16_t));
  UT_ASSERT(i<BMP_MAX_UCHAN);
  if(!bmp->dataUshort[i])
  {
    ALLOCMEM( bmp->dataUshort[i], szNeeded);  
  }
  if( opt & BMP_CLEAR ) 
  {
    memset( bmp->dataUshort[i], 0, szNeeded );
  }
rcCatch:
  return rcThis==MI_OK ? bmp->dataUshort[i] : NULL;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
uint8_t *BmpGetUbyteChannel( BmpBitmap *bmp, int i, unsigned opt )
{
  int rcThis=0;
  size_t szNeeded = (bmp->sx*bmp->sy)*sizeof(bmp->dataUbyte[0][0]);

  UT_ASSERT(sizeof(bmp->dataUbyte[0][0])==1);
  UT_ASSERT(sizeof(bmp->dataUbyte[0][0])==sizeof(uint8_t));
  UT_ASSERT(i<BMP_MAX_UCHAN);
  if(!bmp->dataUbyte[i])
  {
    ALLOCMEM( bmp->dataUbyte[i], szNeeded);  
  }
  if( opt & BMP_CLEAR ) 
  {
    memset( bmp->dataUbyte[i], 0, szNeeded );
  }
rcCatch:
  return rcThis==MI_OK ? bmp->dataUbyte[i] : NULL;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
MtTab *BmpGetFeatureMap( BmpBitmap *bmp, int dim )
{
  MtTab *f=bmp->dataFeat;

  if( !f || (f->nCol != dim))
  {
    MtDelTab(f); f=NULL;
    f = MtCreateTab( bmp->sx*bmp->sy, dim );
    bmp->dataFeat = f;
  }
  return f;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
MtMat *BmpGetDFeatureMap( BmpBitmap *bmp, int dim )
{
  MtMat *f=bmp->dataDFeat;

  if( !f || (f->nCol != dim))
  {
    MtDelMat(f); f=NULL;
    f = MtCreateMat( bmp->sx*bmp->sy, dim );
    bmp->dataDFeat = f;
  }
  return f;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
float *BmpGetZChannel( BmpBitmap *bmp, unsigned opt )
{
  if(!bmp->dataZ)
  {
    if(opt&BMP_CLEAR)
    {
      ALLOCARR0(bmp->dataZ,bmp->sx*bmp->sy);
    }
    else
    {
      ALLOCARR(bmp->dataZ,bmp->sx*bmp->sy);
    }
  }
  return bmp->dataZ;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
unsigned char *BmpGetValidChannel( BmpBitmap *bmp, unsigned opt )
{
  if(!bmp->dataValid)
  {
    if(opt&BMP_CLEAR)
    {
      ALLOCARR0(bmp->dataValid,bmp->sx*bmp->sy);
    }
    else
    {
      ALLOCARR(bmp->dataValid,bmp->sx*bmp->sy);
    }
  }
  return bmp->dataValid;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
float *BmpGetRhoChannel( BmpBitmap *bmp, unsigned opt )
{
  if(!bmp->dataRho)
  {
    if(opt&BMP_CLEAR)
    {
      ALLOCARR0(bmp->dataRho,bmp->sx*bmp->sy);
    }
    else
    {
      ALLOCARR(bmp->dataRho,bmp->sx*bmp->sy);
    }
  }
  return bmp->dataRho;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpGetNChannel( BmpBitmap *bmp, unsigned opt )
{
  int i,rcThis=0;
  for(i=0;i<3;i++)
  {
    if(!bmp->dataN[i])
    {
      if(opt&BMP_CLEAR)
      {
        ALLOCARR0(bmp->dataN[i],bmp->sx*bmp->sy);
      }
      else
      {
        ALLOCARR(bmp->dataN[i],bmp->sx*bmp->sy);
      }
    }
  }
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int *BmpGetIntChannel( BmpBitmap *bmp, int i, unsigned opt )
{
  size_t szNeeded = (bmp->sx*bmp->sy)*sizeof(int);
  if(!bmp->dataInt[i])
  {
    ALLOCMEM( bmp->dataInt[i], szNeeded);  
  }
  if( opt & BMP_CLEAR ) 
  {
    memset( bmp->dataInt[i], 0, szNeeded );
  }

  return bmp->dataInt[i];
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void BmpSetBitmapPixels( BmpBitmap *bmp, unsigned channel, unsigned val )
{
  int i;
  if( channel & BMP_RGB )
  {
    memset( bmp->data, val, (bmp->sx*bmp->sy)*sizeof(int));
  }
  else if( channel & BMP_GREY )
  {
    memset( bmp->dataGrey, val, (bmp->sx*bmp->sy)*sizeof(BmpGreylevel));
  }
  else if( channel & BMP_UV )
  {
    memset( bmp->dataUV, val, (bmp->sx*bmp->sy)*sizeof(BmpUV));
  }
  else if( channel & BMP_ALPHA )
  {
    for(i=0;i<bmp->sx*bmp->sy;i++) bmp->dataAlpha[i]=val;
//    memset( bmp->dataAlpha, val, (bmp->sx*bmp->sy)*sizeof(BmpAlphaVal));
  }
  else if( channel & BMP_USHORT )
  {
    memset( bmp->dataUshort0, val, 
	(bmp->sx*bmp->sy)*sizeof(*bmp->dataUshort0));
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void BmpCopyAlphaToZbuf( BmpBitmap *src, BmpBitmap *dst,
    			 BmpRectangle *srcWin, BmpRectangle *dstWin, 
			 int zval )
{
  int i,j,ch;
  BmpRectangle srcWin0,dstWin0;
  int		srcOffs,srcMod,dstOffs,dstMod;
  BmpAlphaVal  *ppaSrc;
  BmpZval	*ppzDst;

  if( ! srcWin )
  {
    srcWin = &srcWin0;
    BMP_SET_RECT( srcWin, 0, 0, src->sx, src->sy );
  }
  if( ! dstWin )
  {
    dstWin = &dstWin0;
    BMP_SET_RECT( dstWin, 0, 0, dst->sx, dst->sy );
  }

  RSCAN_INFO( src->sx, srcWin->x0, srcWin->y0, srcWin->sx, srcWin->sy,
      		srcOffs, srcMod );
  RSCAN_INFO( dst->sx, dstWin->x0, dstWin->y0, dstWin->sx, dstWin->sy,
      		dstOffs, dstMod );

  ppaSrc = src->dataAlpha+srcOffs; ppzDst = dst->dataZbuf+dstOffs;
  
  for(i=0;i<srcWin->sy;i++, ppaSrc += srcMod, ppzDst += dstMod )
  {
    for(j=0;j<srcWin->sx;j++)
    {
      ch = *ppaSrc++;
      if(ch == 255) 
      {
	*ppzDst = zval;
      }
      ppzDst++;
    }
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BmpBitmap *BmpResizeBitmap( BmpBitmap *bSrc, int sx, int sy, unsigned opt )
{
  int i,j,x,y,rcThis=0,cr,cg,cb,c,doGrey=FALSE,k;
  BmpBitmap *bmp = NULL;
  BmpSamplingInfo si;
  float	dx,dy;
  BmpGreylevel gl;
  MtTab *fs=NULL,*ft=NULL;
  
  if( sx > bSrc->sx || sy > bSrc->sy )
  {
    TRCERRR(("not yet implemented"),MI_ENOTIMPL);
  }

  if( opt & BMP_GREY ) doGrey = TRUE;

  bmp= BmpNewBitmap( sx, sy, BMP_RGB | (doGrey ? BMP_GREY : 0));

  UT_ASSERT_FATAL(bSrc->dataDFeat==NULL);

  if((fs=bSrc->dataFeat)!=NULL)
  {
    ft = BmpGetFeatureMap(bmp,fs->nCol);
    if(!ft) raiseRc(MI_ERROR);
  }

  dx = (float)bSrc->sx / (float)sx;
  dy = (float)bSrc->sy / (float)sy;

  /* TODO: interpolation mode */
  
  BmpInitSampling( bSrc, NULL, dx, dy, 0, sx*sy, 0, &si );
  for(i=0;i<si.ny;i++) for(j=0;j<si.nx;j++)
  {
    BMP_GET_SAMPLE( &si, j, i, x, y );
    BMP_GETPIXEL( bSrc, x, y, c, cr, cg, cb ); 
    BMP_SETPIXEL( bmp, j, i, cr, cg, cb );
    if(doGrey)
    {
      gl = BMP_GETPIXEL_GREY( bSrc, x, y ); 
      BMP_SETPIXEL_GREY( bmp, j, i, gl );
    }

    if(ft)
      for(k=0;k<ft->nCol;k++)
	  ft->rows[BMP_GETPIXEL_OFFS(bmp,j,i)][k] =
	    fs->rows[BMP_GETPIXEL_OFFS(bSrc,x,y)][k];

  }
  BMP_MODIFIED(bmp);

rcCatch: 
  return bmp;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpPutBitmap( BmpBitmap *bmpFrom, 
    			  int x0, int y0,
			  int chan_in, int ichan_in,
			  int chan_out, int ichan_out, 
			  unsigned opt,
			  BmpBitmap *bmpTo )
{
  int rc,rcThis=0,x,y;
  LongInt sx=bmpFrom->sx,sy=bmpFrom->sy;
  BmpBitmap *B=bmpFrom,*O=bmpTo;
  BmpData dataO,dataB;

  UT_ASSERT(O!=B);
  UT_ASSERT(BMP_IN_RANGEB(x0,y0,bmpTo)&&
            BMP_IN_RANGEB(x0+sx-1,y0+sy-1,bmpTo));
 
  rc = BmpGetDataChannel(B,chan_in,ichan_in,0,&dataB);
  if(rc) raiseRc(rc);
  rc = BmpGetDataChannel(O,chan_out,ichan_out,0,&dataO);
  if(rc) raiseRc(rc);

  UT_ASSERT((dataB.pFloat&&dataO.pFloat) || 
            (dataB.pUshort&&dataO.pUshort)); //TODO: other channel types

  #pragma omp parallel for private(y,x)
  for(y=0;y<sy;y++)
  {
    int y2=y+y0;
    if(dataB.pFloat)
    {
      for(x=0;x<sx;x++)
      {
	int x2=x+x0;
	BMP_DPIXEL(O,dataO.pFloat,x2,y2) = BMP_DPIXEL(B,dataB.pFloat,x,y);
      }
    }
    else if(dataB.pUshort)
    {
      for(x=0;x<sx;x++)
      {
	int x2=x+x0;
	BMP_DPIXEL(O,dataO.pUshort,x2,y2) = BMP_DPIXEL(B,dataB.pUshort,x,y);
      }
    }
  }


rcCatch:
  
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpDrawText( BmpBitmap *bmp, int x0, int y0, char *str, int zoom )
{
  /* TODO ! */
  int i,len = strlen(str),xf,yf,xb,yb,x,y;
  LongInt ixy;
  char *pf;
  float g;

  for(i=0;i<len;i++)
  {
    int ch=str[i],
	k = ch>='0'&&ch<='9' ? (tolower((int)(str[i]))-'0') :
	    ch=='.' ? 10 :
	    ch=='+' ? 11 :
	    ch=='-' ? 12 :
	    ch=='E'||ch=='e' ? 13 :
	    14;

    pf = minifont + k*6*6;
      
    for(y=0;y<6*zoom;y++)
      for(x=0;x<6*zoom;x++)
      {
	xf = x/zoom; yf = y/zoom;
	xb = x+x0; yb = y +y0;
	int ch = *(pf + yf*6 + xf);
	ixy=BMP_XYOFF(bmp,xb,yb);
	if(ch=='1') 
	  g=255.0;
	else
	  g=0;
	if(bmp->data)
	  bmp->data[ixy]=BMP_MKRGB((int)(g),(int)(g),(int)(g));
	else if(bmp->dataGrey)
	  bmp->dataGrey[ixy]=g;
/*
	if(ch=='1') 
	{
	  BMP_SETPIXEL(bmp,xb,yb,255,255,255);
	}
	else
	{
	  BMP_SETPIXEL(bmp,xb,yb,0,0,0);
	}*/
      }
    x0 += 6*zoom;
  }
  BMP_MODIFIED(bmp);
  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BmpBitmap *BmpMakeBorderPaddedBitmap( BmpBitmap *bmpIn, 
    				      int borderSz,
    				      int borderExpTyp,
				      int channel )
{
  int rcThis=0,x1,y1,x2,y2,k,x,y,xc,yc,poSrc,poDst;
  BmpBitmap *bmp=NULL;

  if( !((channel & BMP_RGB) ||
      (channel & BMP_GREY) || 
      (channel & BMP_XYZ) ))
  {
    TRCERRR(("not implemented for color channel=%d\n",channel),
	MI_ERROR);
  }

  if(borderExpTyp==0) borderSz=0;

  bmp = BmpNewBitmap(bmpIn->sx + 2*borderSz,
      			 bmpIn->sy + 2*borderSz,channel);
  UT_ASSERT(bmpIn);

  x1=borderSz; y1=borderSz;
  x2=bmp->sx-borderSz-1; y2=bmp->sy-borderSz-1;
  UT_ASSERT(x2>x1&&y2>y1);
  for(y=0;y<bmp->sy;y++) for(x=0;x<bmp->sx;x++)
  {
    xc=x; yc=y;
    if(borderExpTyp==1)
    {
      BMP_CLIP_CONST(x1,y1,x2,y2, xc,yc);
    }
    else if(borderExpTyp==2)
    {
      xc=x; yc=y;
      do
      {
        if(xc<x1) xc = x1 + (x1-xc); 
        if(xc>x2) xc = x2 - (xc-x2); 
        if(yc<y1) yc = y1 + (y1-yc); 
        if(yc>y2) yc = y2 - (yc-y2); 
      }
      while(xc<x1||xc>x2||yc<y1||yc>y2);
    }

    poSrc = BMP_GETPIXEL_OFFS(bmpIn,xc-x1,yc-y1);
    poDst = BMP_GETPIXEL_OFFS(bmp,x,y);

    if(channel & BMP_RGB )
    {
      bmp->data[poDst] = bmpIn->data[poSrc];
    }
    if(channel & BMP_GREY )
    {
      bmp->dataGrey[poDst] = bmpIn->dataGrey[poSrc];
    }
    if(channel & BMP_XYZ )
    {
      for(k=0;k<3;k++)
	bmp->dataFloat[k][poDst] = bmpIn->dataFloat[k][poSrc];
    }
  }
  BMP_MODIFIED(bmp);
rcCatch:
  return bmp;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpCopyBitmapPixels2( BmpBitmap *src, BmpRectangle *srcWin, 
    			  BmpBitmap *dst, BmpRectangle *dstWin,
			  unsigned channel )
{
  int i,rcThis=0,dx,dy,poSrc,poDst,k,j;
  BmpRectangle srcWin0,dstWin0;
  int		srcOffs,srcMod,dstOffs,dstMod,*ppSrc,*ppDst;
  int		lineSz;
  BmpAlphaVal  *ppaSrc, *ppaDst;

  if( ! channel ) channel = BMP_RGB;

  if( ! srcWin )
  {
    srcWin = &srcWin0;
    BMP_SET_RECT( srcWin, 0, 0, src->sx, src->sy );
  }
  if( ! dstWin )
  {
    dstWin = &dstWin0;
    BMP_SET_RECT( dstWin, 0, 0, dst->sx, dst->sy );
  }

  if(!( srcWin->sx == dstWin->sx && srcWin->sy == dstWin->sy ))
  {
    TRCERRR((
      "invalid args for bitmap copy: src='%s':%dx%d dst='%s':%dx%d \n",
	src->name,srcWin->sx,srcWin->sy,dst->name,
	dstWin->sx,dstWin->sy),
	MI_EINVARG);
  }

  RSCAN_INFO( src->sx, srcWin->x0, srcWin->y0, srcWin->sx, srcWin->sy,
      		srcOffs, srcMod );
  RSCAN_INFO( dst->sx, dstWin->x0, dstWin->y0, dstWin->sx, dstWin->sy,
      		dstOffs, dstMod );

  srcMod = src->sx; dstMod = dst->sx;

  if( channel & BMP_RGB )
  {
    lineSz = srcWin->sx * sizeof( *(src->data) );
    ppSrc = src->data+srcOffs; ppDst = dst->data+dstOffs;
    for(i=0;i<srcWin->sy;i++, ppSrc += srcMod, ppDst += dstMod )
    {
      memcpy(ppDst, ppSrc, lineSz );
    }
  }
  if( channel & BMP_ALPHA )
  {
    lineSz = srcWin->sx * sizeof( *(src->dataAlpha) );
    ppaSrc = src->dataAlpha+srcOffs; ppaDst = dst->dataAlpha+dstOffs;
    for(i=0;i<srcWin->sy;i++, ppaSrc += srcMod, ppaDst += dstMod )
    {
      for(j=0;j<srcWin->sx;j++) ppaDst[j] = ppaSrc[j];
//      memcpy(ppaDst, ppaSrc, lineSz );
    }
  }
  if(channel & BMP_GREY )
  {
    BmpGreylevel  *ppgSrc, *ppgDst;

    lineSz = srcWin->sx * sizeof( *(src->dataGrey) );
    ppgSrc = src->dataGrey+srcOffs; ppgDst = dst->dataGrey+dstOffs;
    for(i=0;i<srcWin->sy;i++, ppgSrc += srcMod, ppgDst += dstMod )
    {
      memcpy(ppgDst, ppgSrc, lineSz );
    }
    dst->colConv = src->colConv;
  }
  if(channel & BMP_Z )
  {
    float  *ppgSrc, *ppgDst;

    lineSz = srcWin->sx * sizeof( *(src->dataZ) );
    ppgSrc = src->dataZ+srcOffs; ppgDst = dst->dataZ+dstOffs;
    for(i=0;i<srcWin->sy;i++, ppgSrc += srcMod, ppgDst += dstMod )
    {
      memcpy(ppgDst, ppgSrc, lineSz );
    }
  }
  if(channel & BMP_UV )
  {
    BmpUV  *ppSrc, *ppDst;

    lineSz = srcWin->sx * sizeof( *(src->dataUV) );
    ppSrc = src->dataUV+srcOffs; ppDst = dst->dataUV+dstOffs;
    for(i=0;i<srcWin->sy;i++, ppSrc += srcMod, ppDst += dstMod )
    {
      memcpy(ppDst, ppSrc, lineSz );
    }
    dst->colConv = src->colConv;
  }
  if(channel & BMP_XYZ )
  {
    for(dy=0;dy<dstWin->sy;dy++)
      for(dx=0;dx<dstWin->sx;dx++)
      {
	poDst = BMP_GETPIXEL_OFFS(dst,dstWin->x0+dx,
	    				 dstWin->y0+dy );
	poSrc = BMP_GETPIXEL_OFFS(src,srcWin->x0+dx,
	    				 srcWin->y0+dy );
	for(k=0;k<3;k++)
	  dst->dataFloat[k][poDst] = 
	  	src->dataFloat[k][poSrc];
      }
  }

  BMP_MODIFIED( dst );
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpCopyBitmapPixels( BmpBitmap *src, BmpBitmap *dest )
{
  int rcThis=0;

  if(!( src->sx == dest->sx && src->sy == dest->sy ))
  {
    TRCERRR((
      "invalid args for bitmap copy: src='%s':%p:%dx%d dest='%s':%p:%dx%d \n",
	src->name,src->data,src->sx,src->sy,dest->name,dest->data,
	dest->sx,dest->sy),
	MI_EINVARG);
  }

  if(dest->data&&src->data)
    memcpy( dest->data, src->data, 
      sizeof(int)*(src->sx*src->sy+1));

  if( src->dataGrey && dest->dataGrey )
  {
    memcpy( dest->dataGrey, src->dataGrey, 
	sizeof(BmpGreylevel)*(src->sx*src->sy+1));
  }

  BMP_MODIFIED( dest );
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpPutToDisplayRGB( BmpBitmap *src, 	// source bitmap
    		 	 unsigned chan,		// source channel
		 	 int iChan,		// optional channel number
    		 	 BmpBitmap *dst, 	// destination bitmap
		 	 int x0, int y0         // target coordinates
		 	)
{
  int rcThis=0,x,y;

  BmpGetRGBChannel(dst,0);
//  UT_ASSERT(dst->data);
//
#pragma omp parallel for private(y,x)
  for(y=0;y<src->sy;y++)
    for(x=0;x<src->sx;x++)
    {
      float g,vRgb[3];
      int k,icr,icg,icb,
	  xt = x0+x, yt = y0+y;
      if(!BMP_IN_RANGEB(xt,yt,dst)) continue;
      if(chan==BMP_GREY)
      {
	g = BMP_GETPIXEL_GREY(src,x,y);
	if(g<0) g=0;
	if(g>255) g=255;
	icr=icg=icb=g;
	BMP_SETPIXEL(dst,xt,yt,icr,icg,icb);
      }
      else if(chan==BMP_RGB)
      {
	dst->data[BMP_GETPIXEL_OFFS(dst,xt,yt)] =
	  src->data[BMP_GETPIXEL_OFFS(src,x,y)];
      }
      else if(chan==(BMP_RGB|BMP_XYZ))
      {
	for(k=0;k<3;k++) vRgb[k]=BMP_PIXEL(src,dataFloat[k],x,y);
	dst->data[BMP_GETPIXEL_OFFS(dst,xt,yt)] =
	  BMP_MKRGB((int)vRgb[0],(int)vRgb[1],(int)vRgb[2]);
      }
      else if(chan==BMP_XYZ)
      {
	for(k=0;k<3;k++) vRgb[k]=BMP_PIXEL(src,dataFloat[k],x,y);
	dst->data[BMP_GETPIXEL_OFFS(dst,xt,yt)] =
	  BMP_MKRGB((int)(255*vRgb[0]),(int)(255*vRgb[1]),(int)(255*vRgb[2]));
      }
    }

rcCatch:
  return rcThis;

}

/*-------------------------------------------------------------------------*/   		
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpSubtractBitmapUV( BmpBitmap *bs, BmpBitmap *bt, /* source/target bitmaps */
        	       int ch )			     /* source/target channels */
{
  int 		i,j,rcThis=0;
  BmpUV		*ppuvDst,*ppuvSrc;
  float		vLuvSrc[3]={}, vLuvDst[3]={};
  if( ch != BMP_UV ) TRCERRR(("not implemented"),MI_ENOTIMPL);
  UtProgrInd pgi;

  if(! (ppuvSrc = bs->dataUV) ) TRCERRR(("no UV channel\n"),MI_ERROR);
  if(! (ppuvDst = bt->dataUV) ) TRCERRR(("no UV channel\n"),MI_ERROR);

  UT_PROGRIND_INIT( &pgi );

  for(i=0;i<bs->sy;i++)
  {
    UT_PROGRIND( &pgi, i, 0, bs->sy, "completed: ", NULL );
    for(j=0;j<bs->sx;j++,ppuvSrc++,ppuvDst++)
    {
      BMP_PP_GETPIXEL_UV(ppuvDst,vLuvDst);
      BMP_PP_GETPIXEL_UV(ppuvSrc,vLuvSrc);
      vLuvDst[0] = vLuvSrc[0];
      vLuvDst[1] -= vLuvSrc[1];
      vLuvDst[2] -= vLuvSrc[2];
      BMP_MKUV(ppuvDst,vLuvDst);    
    }
  }

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/   		
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpDiffBitmaps( BmpBitmap *b1, BmpBitmap *b2, int channel,
    		    MtStat *errStat )
{
  int 		i,j,rcThis=0,c;
  float		d;
  int		*pp1, *pp2, cr1,cg1,cb1,cr2,cg2,cb2;
  if( channel != BMP_RGB ) TRCERRR(("not implemented"),MI_ENOTIMPL);
//  UtProgrInd pgi;

  if(!BMP_DIM_EQ(b1,b2)) 
    TRCERRR(("bitmap dimension mismatch"),MI_EINVARG);

  pp1 = b1->data; pp2 = b2->data;
  MT_STAT_INIT( errStat );
//  UT_PROGRIND_INIT( &pgi );

  for(i=0;i<b1->sy;i++)
  {
//    UT_PROGRIND( &pgi, i, 0, b1->sy, "completed: ", NULL );
    for(j=0;j<b1->sx;j++,pp1++,pp2++)
    {
      BMP_PP_GETPIXEL( pp1, c, cr1, cg1, cb1 );
      BMP_PP_GETPIXEL( pp2, c, cr2, cg2, cb2 );
      d = fabs( cr1-cr2 ) + fabs( cg1-cg2 ) + fabs( cb1-cb2 );
      MT_STAT_UPDATE( errStat, d );
    }
  }
  MT_STAT_RESULT( errStat );
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpGetGradient( int type,
    			BmpBitmap *bmp, float *data,	
    			int x, int y, 
    			double h, double *vg )
{
  int rcThis=0;
  switch(type)
  {
    case BMP_GRAD_SOEBEL :
      BmpGradientSoebel(bmp,data,x,y,h,vg);
      break;
    case BMP_GRAD_CDIFF :
      BmpGradientCentDiff(bmp,data,x,y,h,vg);
      break;
    case BMP_GRAD_FDIFF :
      BmpGradientForwDiff(bmp,data,x,y,h,vg);
      break;
    case BMP_GRAD_BDIFF :
      BmpGradientBackDiff(bmp,data,x,y,h,vg);
      break;
    default:
      TRCERRR(("invalid gradient type '%d'\n",type),MI_EINVARG);
      break;
  }
rcCatch:
  return rcThis;
}
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void BmpGradientSoebel( BmpBitmap *bmp, float *data,	
    			int x, int y, 
    			double h, double *vg )
{
    int px,nx,py,ny;
    double Ipp,Ipc,Ipn,Icp,Icn,Inp,Inc,Inn,IppInn,IpnInp;
    double c1 = (double)(0.25*(2-MT_SQRT2)), 
	   c2 = (double)(0.5f*(MT_SQRT2-1));

    px=x-1; if (px<0) px=0;
    nx=x+1; if (nx==bmp->sx) nx--;
    py=y-1; if (py<0) py=0;
    ny=y+1; if (ny==bmp->sy) ny--;

    Ipp=data[BMP_XYOFF(bmp,px,py)] ;
    Ipc=data[BMP_XYOFF(bmp,px,y)] ;
    Ipn=data[BMP_XYOFF(bmp,px,ny)] ;
    Icp=data[BMP_XYOFF(bmp,x,py)] ;
    Icn=data[BMP_XYOFF(bmp,x,ny)] ;
    Inp=data[BMP_XYOFF(bmp,nx,py)] ;
    Inc=data[BMP_XYOFF(bmp,nx,y)] ;
    Inn=data[BMP_XYOFF(bmp,nx,ny)] ;
    IppInn = c1*(Inn-Ipp);
    IpnInp = c1*(Ipn-Inp);
    vg[0]  = ((double)(IppInn-IpnInp-c2*Ipc+c2*Inc))/h;
    vg[1]  = ((double)(IppInn+IpnInp-c2*Icp+c2*Icn))/h;
  
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void BmpGradientCentDiff( BmpBitmap *bmp, float *data,
    			  int x, int y, 
    			   double h, double	*vg )
{
    int x0=x,y0=y,x1,y1,x2,y2,xmax=bmp->sx-1,ymax=bmp->sy-1;
    float *bmp_data = data;

    x1=MAX(x0-1,0); x2=MIN(x0+1,xmax);
    vg[0] = 0.5*(bmp_data[BMP_XYOFF(bmp,x2,y0)] -
	          bmp_data[BMP_XYOFF(bmp,x1,y0)]) / h;	

    y1=MAX(y0-1,0); y2=MIN(y0+1,ymax);
    vg[1] = 0.5*(bmp_data[BMP_XYOFF(bmp,x0,y2)] -
	          bmp_data[BMP_XYOFF(bmp,x0,y1)]) / h;	  
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void BmpGradientForwDiff( BmpBitmap *bmp, float *data,
    			  int x, int y, 
    			   double h, double	*vg )
{
    int x0=x,y0=y,x1,y1,x2,y2,xmax=bmp->sx-1,ymax=bmp->sy-1;
    float *bmp_data = data;

    x1=x0; x2=MIN(x0+1,xmax);
    vg[0] = (bmp_data[BMP_XYOFF(bmp,x2,y0)] -
	          bmp_data[BMP_XYOFF(bmp,x1,y0)]) / h;	

    y1=y0; y2=MIN(y0+1,ymax);
    vg[1] = (bmp_data[BMP_XYOFF(bmp,x0,y2)] -
	          bmp_data[BMP_XYOFF(bmp,x0,y1)]) / h;	  
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void BmpGradientBackDiff( BmpBitmap *bmp, float *data,
    			  int x, int y, 
    			   double h, double	*vg )
{
    int x0=x,y0=y,x1,y1,x2,y2;
    float *bmp_data = data;

    x1=MAX(x0-1,0); x2=x0;
    vg[0] = (bmp_data[BMP_XYOFF(bmp,x2,y0)] -
	          bmp_data[BMP_XYOFF(bmp,x1,y0)]) / h;	

    y1=MAX(y0-1,0); y2=y0;
    vg[1] = (bmp_data[BMP_XYOFF(bmp,x0,y2)] -
	          bmp_data[BMP_XYOFF(bmp,x0,y1)]) / h;	  
}



/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
double BmpGetLumaGradient( BmpBitmap *bmp, int radius,
    			   int x0, int y0,
    			   int doClip,
    			   double *pDirection )
{
  register int x,y,r=radius,rcThis=0;
  register BmpGreylevel *pp,lum;
  double lumOut=0,wsum;

  register int sumX,sumY,*ppmX,*ppmY;
  int	       boff,bmod;

  static int sobel_X[3*3] = {  -1, 0, 1, 
    			       -2, 0, 2,
			       -1, 0, 1  };

  static int sobel_Y[3*3] = {   1,  2,  1, 
    				0,  0,  0,
			       -1, -2, -1  };

  static int sobel5_X[5*5] = {  1, 2, 0, -2, -1,
    				4, 8, 0, -8, -4,
				6, 12, 0, -12, -6,
				4, 8, 0, -8, -4,
				1, 2, 0, -2, -1 };

  static int sobel5_Y[5*5] = {-1,-4,-6,-4,-1,
    				-2,-8,-12,-8,-2,
				  0,0,0,0,0,
				  2,8,12,8,2,
				  1,4,6,4,1 };

  if(doClip)
  {
    BMP_CLIP_REGION( bmp,r/2,x0,y0,r,r );
  }

  sumX = 0; sumY = 0;

  if(r==3)
  {
    ppmX = sobel_X; ppmY = sobel_Y;
    wsum=512;	// max. kernel sum = 4
  }
  else if(r==5)
  {
    ppmX = sobel5_X; ppmY = sobel5_Y;
    wsum=8192;	// max. sum=48
  }
  else 
  {
    TRCERRR(("invalid radius\n"),MI_ENOTIMPL);
  }

  RSCAN_INFO( bmp->sx, x0, y0, r, r, boff, bmod );

  RSCAN_LOOP( pp, x, y, bmp->dataGrey, boff, r, r, bmod )
  {
    lum = *pp;
    sumX += (*ppmX++) * lum; sumY += (*ppmY++) * lum;
  }

  lumOut = (ABS(sumX) + ABS(sumY))/wsum;

//  lumOut = sqrt(sumX*sumX+sumY*sumY) / 361;

  if(pDirection) 
  {
    *pDirection = atan2((double)sumY,(double)sumX);
  }
rcCatch:
  return lumOut;
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BmpBitmap *BmpGradientSobelXYM( BmpBitmap *bmp, BmpRectangle *roi,
    				float *dataChannel,
    				int radius )
/*
 *  on return: dataFloat[0]: Gradient X-component
 *	       dataFloat[1]: Gradient Y-component
 *	       dataFloat[2]: Gradient Magnitude
 *
 */ 
{
  int	r=radius,i,k,rcThis=0;
  register int x,y,r2=r/2,x1,y1,rx,ry,x2,y2,xn,yn,dx,dy,wx,wy;
  register BmpGreylevel lum;
  double wsum,divX,divY,divM;
  BmpBitmap *bmpOut=NULL;
  float *pgX,*pgY,*pgM,*data = dataChannel ? dataChannel : bmp->dataGrey;
  register int *ppmX,*ppmY;
  register float sumX,sumY,wsumX,wsumY;
  int	      po,*ppmX0,*ppmY0;
  MtStat	statX,statY,statM;
  static int sobel_X[3*3] = {  -1, 0, 1, 
    			       -2, 0, 2,
			       -1, 0, 1  };

  static int sobel_Y[3*3] = {   1,  2,  1, 
    				0,  0,  0,
			       -1, -2, -1  };

  static int sobel5_X[5*5] = {  1, 2, 0, -2, -1,
    				4, 8, 0, -8, -4,
				6, 12, 0, -12, -6,
				4, 8, 0, -8, -4,
				1, 2, 0, -2, -1 };

  static int sobel5_Y[5*5] = {-1,-4,-6,-4,-1,
    				-2,-8,-12,-8,-2,
				  0,0,0,0,0,
				  2,8,12,8,2,
				  1,4,6,4,1 };

  BMP_GET_WIN( bmp, roi, x1, y1, rx, ry );
  x2=x1+rx-1; y2=y1+ry-1;

  bmpOut = BmpNewBitmap(bmp->sx,bmp->sy,BMP_GEN);
  if(!bmpOut) raiseRc(MI_ERROR);

  pgX = BmpGetFloatChannel(bmpOut,0,0);
  pgY = BmpGetFloatChannel(bmpOut,1,0);
  pgM = BmpGetFloatChannel(bmpOut,2,0);

  if(!(pgX&&pgY&&pgM)) 
    TRCERRR(("no data channels\n"),MI_ERROR);

  if(r==3)
  {
    ppmX0 = sobel_X; ppmY0 = sobel_Y;
  }
  else if(r==5)
  {
    ppmX0 = sobel5_X; ppmY0 = sobel5_Y;
  }
  else 
  {
    TRCERRR(("invalid radius\n"),MI_ENOTIMPL);
  }

  for(wsum=0,k=0;k<r*r;k++) wsum+=fabs(ppmX0[k]);

  MT_STAT_INIT(&statX);
  MT_STAT_INIT(&statY);
  MT_STAT_INIT(&statM);

  for(y=y1;y<=y2;y++)
    for(x=x1;x<=x2;x++)
    {
      sumX = sumY = wsumX = wsumY = 0;
      ppmX = ppmX0; ppmY = ppmY0;
      for(dy=-r2;dy<=r2;dy++)
	for(dx=-r2;dx<=r2;dx++)
	{
	  wx = *ppmX++; wy = *ppmY++;
	  xn=x+dx; yn=y+dy;
	  BMP_CLIP_CONST(x1,y1,x2,y2, xn,yn);
//	  if(!BMP_IN_RANGE(xn,yn,x1,y1,x2,y2)) continue;
	  lum = data[BMP_GETPIXEL_OFFS(bmp,xn,yn)];
	  sumX += lum*wx; sumY += lum*wy;
	  wsumX += wx; wsumY += wy;
	}

      /*
      BMP_CLIP_REGION( bmp,r/2,x0,y0,r,r );
      sumX = 0; sumY = 0;
      ppmX = ppmX0; ppmY = ppmY0;
      RSCAN_INFO( bmp->sx, x0, y0, r, r, boff, bmod );
      RSCAN_LOOP( pp, ix, iy, data, boff, r, r, bmod )
      {
	lum = *pp;
	sumX += (*ppmX++) * lum; sumY += (*ppmY++) * lum;
      }
      */
      po = BMP_GETPIXEL_OFFS(bmp,x,y);
/*      
      if(fabs(wsumX)>MT_REAL_EPSILON)
        sumX /= fabs(wsumX); 
      if(fabs(wsumY)>MT_REAL_EPSILON)
        sumY /= fabs(wsumY); 
*/
      pgX[po] = sumX;
      pgY[po] = sumY;
      pgM[po] = sqrt(sumX*sumX + sumY*sumY);

      MT_STAT_UPDATE(&statX,pgX[po]);
      MT_STAT_UPDATE(&statY,pgY[po]);
      MT_STAT_UPDATE(&statM,pgM[po]);
    }
  MT_STAT_RESULT(&statX);
  MT_STAT_RESULT(&statY);
  MT_STAT_RESULT(&statM);

  // normalize
  divX = statX.var; if(divX < MT_REAL_EPSILON) divX = 1;
  divY = statY.var; if(divY < MT_REAL_EPSILON) divY = 1;
  divM = statM.max; if(divM < MT_REAL_EPSILON) divM = 1;

  for(i=0;i<bmp->sx*bmp->sy;i++) 
  {
    pgX[i] = ( pgX[i] - statX.avg ) / divX;
    pgY[i] = ( pgY[i] - statY.avg ) / divY;
    pgM[i] = pgM[i] / divM;
  }
  

rcCatch:
  if(rcThis) 
  {
    BmpCloseBitmap(bmpOut);
    bmpOut=NULL;
  }
  return bmpOut;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
typedef struct
{
  int r;
}
MwfCb_Soebel_info;

static BmpGreylevel MwfCb_SoebelGrey( int boff, int bmod, BmpGreylevel *bdata, 
    			 void *info  )
{
  register int x,y,r=3;
  register BmpGreylevel *pp,lum;

  register int sumX,sumY,*ppmX,*ppmY;
  int	       lumOut;

  static int convMask_Sobel_X[3*3] = {  -1, 0, 1, 
    				 	-2, 0, 2,
				 	-1, 0, 1  };

  static int convMask_Sobel_Y[3*3] = {   1,  2,  1, 
    				  	 0,  0,  0,
				 	-1, -2, -1  };

  sumX = 0; sumY = 0;
  ppmX = convMask_Sobel_X; ppmY = convMask_Sobel_Y;

  RSCAN_LOOP( pp, x, y, bdata, boff, r, r, bmod )
  {
    lum = *pp;
    sumX += (*ppmX++) * lum; sumY += (*ppmY++) * lum;
  }

  lumOut = ABS(sumX) + ABS(sumY);
  if(lumOut > 255 ) lumOut = 255;
  if(lumOut < 0 ) lumOut = 0;
  
  return lumOut;
}

/*-------------------------------------------------------------------------*/
/* 	MwfCb_Smooth		    				           */
/*-------------------------------------------------------------------------*/
typedef struct
{
  int r;
  float *kernel;
  float	sumw;
}
MwfCb_Smooth_Info;

static int MwfCb_Smooth( int boff, int bmod, BmpAlphaVal *bdata, 
    			      void *pinfo  )
{
  register int x,y;
  register BmpAlphaVal *pp;
  MwfCb_Smooth_Info *info = (MwfCb_Smooth_Info *)pinfo;
  register int lumOut,lum,r=info->r;
  float		sum,sumw,kf,*kp=info->kernel;


  sum = 0, sumw=0;

  RSCAN_LOOP( pp, x, y, bdata, boff, r, r, bmod )
  {
    lum = *pp;
    kf = *kp++;
    sum += (float)lum * kf; 
  }

  lumOut = sum/info->sumw;
  if(lumOut > 255 ) lumOut = 255;
  if(lumOut < 0 ) lumOut = 0;
  
  return lumOut;
}

BmpUV globTmpUV;

static BmpUV MwfCb_SmoothUV( int boff, int bmod, BmpUV *bdata, 
    			          void *pinfo  )
{
  int rcThis=0;
  register int x,y;
  register BmpUV *pp;
  MwfCb_Smooth_Info *info = (MwfCb_Smooth_Info *)pinfo;
  register int r=info->r;
  float		sum1,sum2,sumw,kf,*kp=info->kernel,vLuv[3];


  sum1 = sum2 = 0, sumw=0;

  RSCAN_LOOP( pp, x, y, bdata, boff, r, r, bmod )
  {
    BMP_PP_GETPIXEL_UV(pp,vLuv);
    kf = *kp++;
    sum1 += vLuv[1] * kf; 
    sum2 += vLuv[2] * kf; 
  }

  vLuv[1] = sum1/info->sumw;
  vLuv[2] = sum2/info->sumw;
  BMP_MKUV( &globTmpUV, vLuv );

  UT_ASSERT(FALSE);  //TODO (return value?)
rcCatch:
  return globTmpUV;
}

static int MwfCb_SmoothRGB( int boff, int bmod, int *bdata, 
    			          void *pinfo  )
{
  register int x,y,icrgb,icr,icg,icb;
  register int *pp;
  MwfCb_Smooth_Info *info = (MwfCb_Smooth_Info *)pinfo;
  register int r=info->r;
  float		sum1,sum2,sum3,sumw,kf,*kp=info->kernel;

  sum1 = sum2 = sum3 = 0, sumw=0;

  RSCAN_LOOP( pp, x, y, bdata, boff, r, r, bmod )
  {
    BMP_PP_GETPIXEL(pp,icrgb,icr,icg,icb);
    kf = *kp++;
    sum1 += icr * kf; 
    sum2 += icg * kf; 
    sum3 += icb * kf; 
  }

  icr = sum1/info->sumw;
  icg = sum2/info->sumw;
  icb = sum3/info->sumw;

  return BMP_MKRGB(icr,icg,icb);
}

static BmpGreylevel MwfCb_SmoothGrey( int boff, int bmod, 
    					   BmpGreylevel *bdata, 
    			      	  	   void *pinfo  )
{
  register int x,y;
  register BmpGreylevel *pp;
  MwfCb_Smooth_Info *info = (MwfCb_Smooth_Info *)pinfo;
  register int lumOut,lum,r=info->r;
  float		sum,sumw,kf,*kp=info->kernel;


  sum = 0, sumw=0;

  RSCAN_LOOP( pp, x, y, bdata, boff, r, r, bmod )
  {
    lum = *pp;
    kf = *kp++;
    sum += (float)lum * kf; 
  }

  lumOut = sum/info->sumw;
  
//  printf("%f ",(float)lumOut);

  return lumOut;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpMovWinFilterAlpha( BmpBitmap *src, BmpBitmap *dst, 
    		     int srcChannel, int dstChannel,
		     BmpRectangle *srcWin,
        	     BmpRectangle *dstWin,
		     int d, 
		     BmpMovWinFunctype movWinCallback,
		     void *info )
{
  int r=d,x,y,sx,sy,rcThis=0,
  	csx,csy,u,v,mwSrcOffs,mwSrcMod,dstOffs,dstMod;
  BmpRectangle srcWin0,dstWin0;
  BmpAlphaVal *bmpData=NULL, *bmpDestData=NULL, *ppDest;

  /**** check args */

  UT_ASSERT(FALSE); // TODO: changed alpha-type from char to float
  
  if( ! src || ! dst )
  {
    TRCERRR(("bitmap NULL pointer"),MI_EINVARG);
  }

  if( src == dst )
  {
    TRCERRR((
	  "can't apply MOVWIN filter to asme source and target bitmap %p",
	  src),MI_EINVARG);
  }

  if( (srcChannel & BMP_RGB) ||
      (dstChannel & BMP_RGB) ||
      (srcChannel & BMP_GREY) ||
      (dstChannel & BMP_GREY) )
  {
    TRCERRR(("channel not yet implemented\n"),MI_ENOTIMPL);
  }

  if( ! srcWin )
  {
    srcWin = &srcWin0; BMP_SET_RECT( srcWin,0,0,src->sx,src->sy );
  }
  if( ! dstWin )
  {
    dstWin = &dstWin0; BMP_SET_RECT( dstWin,0,0,dst->sx,dst->sy );
  }

  if(!( srcWin->sx == dstWin->sx && srcWin->sy == dstWin->sy ))
  {
    TRCERRR((
      "different window dimensions: src=%dx%d dst=%dx%d \n",
	srcWin->sx,srcWin->sy, dstWin->sx,dstWin->sy),
	MI_EINVARG);
  }

  sx = srcWin->sx; sy = srcWin->sy;

  if( sx < r || sy < r ) 
  {
    TRCERRR(("subwindow too small (%d)\n",r),MI_EINVARG);
  }

  /**** prepare RSCAN info */
  RSCAN_INFO( dst->sx, dstWin->x0, dstWin->y0, dstWin->sx, dstWin->sy,
      		dstOffs, dstMod );

  if( srcChannel & BMP_ALPHA )
    bmpData = src->dataAlpha;

  if( srcChannel & BMP_ALPHA )
    bmpDestData = dst->dataAlpha;

  csx = src->sx; csy = src->sy;

  CLOCK_START();

  /**** scan destination bitmap window and for each pixel call scanner 
        for moving rxr-source window */

  RSCAN_LOOP( ppDest, x, y, bmpDestData, dstOffs, sx, sy, dstMod )
  { 
    u = x-r/2; v = y-r/2;
    BMP_CLIP_RECT_MP( csx, csy, u, v, r, r );
    RSCAN_INFO( src->sx, u, v, r, r, mwSrcOffs, mwSrcMod );

    *ppDest = movWinCallback( mwSrcOffs, mwSrcMod, bmpData, info ); 
  }

  CLOCK_STOP();
  TRC(("filter processing time = %d ms (%.1f pixels/sec)\n",
	(int)CLOCK_DIFF_MS(), 
	1000.0*(float)(sx*sy)/(float)CLOCK_DIFF_MS()));

rcCatch:

  BMP_MODIFIED( dst );

  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpMovWinFilterUV( BmpBitmap *src, BmpBitmap *dst, 
		       BmpRectangle *srcWin,
        	       BmpRectangle *dstWin,
		       int d, 
		       BmpMovWinUVFunctype movWinCallback,
		       void *info )
{
  int r=d,x,y,sx,sy,rcThis=0,xs0,ys0,
  	csx,csy,u,v,mwSrcOffs,mwSrcMod,dstOffs,dstMod;
  BmpRectangle srcWin0,dstWin0;
  UtProgrInd pgi;

  BmpUV *bmpData=NULL, *bmpDestData=NULL, *ppDest;

  /**** check args */

  if( ! src || ! dst )
  {
    TRCERRR(("bitmap NULL pointer"),MI_EINVARG);
  }

  if( src == dst )
  {
    TRCERRR((
	  "can't apply MOVWIN filter to same source and target bitmap %p",
	  src),MI_EINVARG);
  }

  if( ! srcWin )
  {
    srcWin = &srcWin0; BMP_SET_RECT( srcWin,0,0,src->sx,src->sy );
  }
  if( ! dstWin )
  {
    dstWin = &dstWin0; BMP_SET_RECT( dstWin,0,0,dst->sx,dst->sy );
  }

  if(!( srcWin->sx == dstWin->sx && srcWin->sy == dstWin->sy ))
  {
    TRCERRR((
      "different window dimensions: src=%dx%d dst=%dx%d \n",
	srcWin->sx,srcWin->sy, dstWin->sx,dstWin->sy),
	MI_EINVARG);
  }

  sx = srcWin->sx; sy = srcWin->sy;

  if( sx < r || sy < r ) 
  {
    TRCERRR(("subwindow too small (%d)\n",r),MI_EINVARG);
  }

  /**** prepare RSCAN info */
  RSCAN_INFO( dst->sx, dstWin->x0, dstWin->y0, dstWin->sx, dstWin->sy,
      		dstOffs, dstMod );
  bmpData = src->dataUV;
  bmpDestData = dst->dataUV;
  csx = src->sx; csy = src->sy;
  xs0 = srcWin->x0; ys0 = srcWin->y0;

  UT_PROGRIND_INIT( &pgi );

  CLOCK_START();

  /**** scan destination bitmap window and for each pixel call scanner 
        for moving rxr-source window */

  RSCAN_LOOP( ppDest, x, y, bmpDestData, dstOffs, sx, sy, dstMod )
  { 
    UT_PROGRIND( &pgi, y, 0, sy, "filter completed: ", NULL );
    
    u = x-r/2+xs0; v = y-r/2+ys0;
    BMP_CLIP_RECT_MP( csx, csy, u, v, r, r );
    RSCAN_INFO( csx, u, v, r, r, mwSrcOffs, mwSrcMod );

    *ppDest = movWinCallback( mwSrcOffs, mwSrcMod, bmpData, info ); 
  }

  CLOCK_STOP();
  TRC(("filter processing time = %d ms (%.1f pixels/sec)\n",
	(int)CLOCK_DIFF_MS(), 
	1000.0*(float)(sx*sy)/(float)CLOCK_DIFF_MS()));

rcCatch:

  BMP_MODIFIED( dst );

  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpMovWinFilterRGB( BmpBitmap *src, BmpBitmap *dst, 
		       BmpRectangle *srcWin,
        	       BmpRectangle *dstWin,
		       int d, 
		       BmpMovWinRGBFunctype movWinCallback,
		       void *info )
{
  int r=d,x,y,sx,sy,rcThis=0,xs0,ys0,
  	csx,csy,u,v,mwSrcOffs,mwSrcMod,dstOffs,dstMod;
  BmpRectangle srcWin0,dstWin0;
  UtProgrInd pgi;

  int *bmpData=NULL, *bmpDestData=NULL, *ppDest;

  /**** check args */

  if( ! src || ! dst )
  {
    TRCERRR(("bitmap NULL pointer"),MI_EINVARG);
  }

  if( src == dst )
  {
    TRCERRR((
	  "can't apply MOVWIN filter to same source and target bitmap %p",
	  src),MI_EINVARG);
  }

  if( ! srcWin )
  {
    srcWin = &srcWin0; BMP_SET_RECT( srcWin,0,0,src->sx,src->sy );
  }
  if( ! dstWin )
  {
    dstWin = &dstWin0; BMP_SET_RECT( dstWin,0,0,dst->sx,dst->sy );
  }

  if(!( srcWin->sx == dstWin->sx && srcWin->sy == dstWin->sy ))
  {
    TRCERRR((
      "different window dimensions: src=%dx%d dst=%dx%d \n",
	srcWin->sx,srcWin->sy, dstWin->sx,dstWin->sy),
	MI_EINVARG);
  }

  sx = srcWin->sx; sy = srcWin->sy;

  if( sx < r || sy < r ) 
  {
    TRCERRR(("subwindow too small (%d)\n",r),MI_EINVARG);
  }

  /**** prepare RSCAN info */
  RSCAN_INFO( dst->sx, dstWin->x0, dstWin->y0, dstWin->sx, dstWin->sy,
      		dstOffs, dstMod );
  bmpData = src->data;
  bmpDestData = dst->data;
  csx = src->sx; csy = src->sy;
  xs0 = srcWin->x0; ys0 = srcWin->y0;

  UT_PROGRIND_INIT( &pgi );

  CLOCK_START();

  /**** scan destination bitmap window and for each pixel call scanner 
        for moving rxr-source window */

  RSCAN_LOOP( ppDest, x, y, bmpDestData, dstOffs, sx, sy, dstMod )
  { 
    UT_PROGRIND( &pgi, y, 0, sy, "filter completed: ", NULL );
    
    u = x-r/2+xs0; v = y-r/2+ys0;
    BMP_CLIP_RECT_MP( csx, csy, u, v, r, r );
    RSCAN_INFO( csx, u, v, r, r, mwSrcOffs, mwSrcMod );

    *ppDest = movWinCallback( mwSrcOffs, mwSrcMod, bmpData, info ); 
  }

  CLOCK_STOP();
  TRC(("filter processing time = %d ms (%.1f pixels/sec)\n",
	(int)CLOCK_DIFF_MS(), 
	1000.0*(float)(sx*sy)/(float)CLOCK_DIFF_MS()));

rcCatch:

  BMP_MODIFIED( dst );

  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpMovWinFilterGrey( BmpBitmap *src, BmpBitmap *dst, 
		         BmpRectangle *srcWin,
        	         BmpRectangle *dstWin,
		         int d, 
		         BmpMovWinGreyFunctype movWinCallback,
		         void *info )
{
  int r=d,x,y,sx,sy,rcThis=0,xs0,ys0,
  	csx,csy,u,v,mwSrcOffs,mwSrcMod,dstOffs,dstMod;
  BmpRectangle srcWin0,dstWin0;
  UtProgrInd pgi;

  BmpGreylevel *bmpData=NULL, *bmpDestData=NULL, *ppDest;

  /**** check args */

  if( ! src || ! dst )
  {
    TRCERRR(("bitmap NULL pointer"),MI_EINVARG);
  }

  if( src == dst )
  {
    TRCERRR((
	  "can't apply MOVWIN filter to same source and target bitmap %p",
	  src),MI_EINVARG);
  }

  if( ! srcWin )
  {
    srcWin = &srcWin0; BMP_SET_RECT( srcWin,0,0,src->sx,src->sy );
  }
  if( ! dstWin )
  {
    dstWin = &dstWin0; BMP_SET_RECT( dstWin,0,0,dst->sx,dst->sy );
  }

  if(!( srcWin->sx == dstWin->sx && srcWin->sy == dstWin->sy ))
  {
    TRCERRR((
      "different window dimensions: src=%dx%d dst=%dx%d \n",
	srcWin->sx,srcWin->sy, dstWin->sx,dstWin->sy),
	MI_EINVARG);
  }

  sx = srcWin->sx; sy = srcWin->sy;

  if( sx < r || sy < r ) 
  {
    TRCERRR(("subwindow too small (%d)\n",r),MI_EINVARG);
  }

  /**** prepare RSCAN info */
  RSCAN_INFO( dst->sx, dstWin->x0, dstWin->y0, dstWin->sx, dstWin->sy,
      		dstOffs, dstMod );
  bmpData = src->dataGrey;
  bmpDestData = dst->dataGrey;
  csx = src->sx; csy = src->sy;
  xs0 = srcWin->x0; ys0 = srcWin->y0;

  UT_PROGRIND_INIT( &pgi );

  CLOCK_START();

  /**** scan destination bitmap window and for each pixel call scanner 
        for moving rxr-source window */

  RSCAN_LOOP( ppDest, x, y, bmpDestData, dstOffs, sx, sy, dstMod )
  { 
    UT_PROGRIND( &pgi, y, 0, sy, "filter completed: ", NULL );

    u = x-r/2+xs0; v = y-r/2+ys0;
    BMP_CLIP_RECT_MP( csx, csy, u, v, r, r );
    RSCAN_INFO( csx, u, v, r, r, mwSrcOffs, mwSrcMod );

    *ppDest = movWinCallback( mwSrcOffs, mwSrcMod, bmpData, info ); 
  }

  CLOCK_STOP();
  TRC(("filter processing time = %d ms (%.1f pixels/sec)\n",
	(int)CLOCK_DIFF_MS(), 
	1000.0*(float)(sx*sy)/(float)CLOCK_DIFF_MS()));

rcCatch:

  BMP_MODIFIED( dst );

  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpMovWinFilterFloat( BmpBitmap *src, BmpBitmap *dst, 
		         BmpRectangle *srcWin,
        	         BmpRectangle *dstWin,
		         int d, 
			 int iChan,
		         BmpMovWinGreyFunctype movWinCallback,
		         void *info )
{
  int r=d,x,y,sx,sy,rcThis=0,xs0,ys0,
  	csx,csy,u,v,mwSrcOffs,mwSrcMod,dstOffs,dstMod;
  BmpRectangle srcWin0,dstWin0;
  UtProgrInd pgi;

  BmpGreylevel *bmpData=NULL, *bmpDestData=NULL, *ppDest;

  /**** check args */

  if( ! src || ! dst )
  {
    TRCERRR(("bitmap NULL pointer"),MI_EINVARG);
  }

  if( src == dst )
  {
    TRCERRR((
	  "can't apply MOVWIN filter to same source and target bitmap %p",
	  src),MI_EINVARG);
  }

  if( ! srcWin )
  {
    srcWin = &srcWin0; BMP_SET_RECT( srcWin,0,0,src->sx,src->sy );
  }
  if( ! dstWin )
  {
    dstWin = &dstWin0; BMP_SET_RECT( dstWin,0,0,dst->sx,dst->sy );
  }

  if(!( srcWin->sx == dstWin->sx && srcWin->sy == dstWin->sy ))
  {
    TRCERRR((
      "different window dimensions: src=%dx%d dst=%dx%d \n",
	srcWin->sx,srcWin->sy, dstWin->sx,dstWin->sy),
	MI_EINVARG);
  }

  sx = srcWin->sx; sy = srcWin->sy;

  if( sx < r || sy < r ) 
  {
    TRCERRR(("subwindow too small (%d)\n",r),MI_EINVARG);
  }

  /**** prepare RSCAN info */
  RSCAN_INFO( dst->sx, dstWin->x0, dstWin->y0, dstWin->sx, dstWin->sy,
      		dstOffs, dstMod );
  bmpData = src->dataFloat[iChan];
  bmpDestData = dst->dataFloat[iChan];
  csx = src->sx; csy = src->sy;
  xs0 = srcWin->x0; ys0 = srcWin->y0;

  UT_PROGRIND_INIT( &pgi );

  CLOCK_START();

  /**** scan destination bitmap window and for each pixel call scanner 
        for moving rxr-source window */

  RSCAN_LOOP( ppDest, x, y, bmpDestData, dstOffs, sx, sy, dstMod )
  { 
    UT_PROGRIND( &pgi, y, 0, sy, "filter completed: ", NULL );

    u = x-r/2+xs0; v = y-r/2+ys0;
    BMP_CLIP_RECT_MP( csx, csy, u, v, r, r );
    RSCAN_INFO( csx, u, v, r, r, mwSrcOffs, mwSrcMod );

    *ppDest = movWinCallback( mwSrcOffs, mwSrcMod, bmpData, info ); 
  }

  CLOCK_STOP();
  TRC(("filter processing time = %d ms (%.1f pixels/sec)\n",
	(int)CLOCK_DIFF_MS(), 
	1000.0*(float)(sx*sy)/(float)CLOCK_DIFF_MS()));

rcCatch:

  BMP_MODIFIED( dst );

  return rcThis;
}



/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpFilterSoebel( BmpBitmap *src, BmpBitmap *dst, 
        	     int channel,
    		     BmpRectangle *srcWin, BmpRectangle *dstWin, 
    		     int r)
{
  int rcThis=0,rc;
  MwfCb_Soebel_info info;

  info.r = r;

  if( channel & BMP_GREY )
  {
    if( ! src->dataGrey )
    {
      TRCERRR(("bitmap has no Grey channels\n"),MI_ERROR);
    }
    rc = BmpMovWinFilterGrey( src, dst, srcWin, dstWin, r,
			      MwfCb_SoebelGrey, &info );
    if(rc) raiseRc(rc);
  }
  else TRCERRR(("invalid bitmap chanell %d\n",channel ),MI_EINVARG);

rcCatch:
  return rcThis;
}
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BmpBitmap *BmpMakeKernel(int r, float s, float shrp, unsigned opt )
{
  int rcThis=0;
  BmpBitmap *bmp = NULL;
  
  bmp = BmpNewBitmap(r,r,BMP_GREY);
  if(!bmp) raiseRc(MI_ERROR);

  if(s < MT_REAL_EPSILON )
  {
    int k,m=r*r;
    for(k=0;k<bmp->sx*bmp->sy;k++) bmp->dataGrey[k]=1.0/(double)m;
  }
  else
  {
    rcThis = BmpFilterSmooth( NULL, bmp, BMP_GREY, 0, NULL, NULL, r, s, shrp, 
      		              opt | BMP_FKERN );
  }

rcCatch:
  return bmp;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpFilterSmooth( BmpBitmap *src, BmpBitmap *dst, 
        	     	  int channel, int iChannel,
    			  BmpRectangle *srcWin, BmpRectangle *dstWin, 
			  int r, float s, float shrp, unsigned opt )
{
  int			i,j,rc,rcThis=0,dstIs0;
  float			*kp,dist,g,sum=0,sumg;
  MwfCb_Smooth_Info info={0};

#define KRN_DIST2( x1,y1, x2,y2 ) \
  sqrt(((x2)-(x1))*((x2)-(x1))+((y2)-(y1))*((y2)-(y1)))

  if(opt & BMP_FKERN )
  {
    kp = dst->dataGrey;
    info.kernel = kp;
  }
  else
  {
    info.r = r;
    info.sumw = 0;
    ALLOCMEM(info.kernel,sizeof(float)*r*r);
    kp = info.kernel;
  }

  printf("Kernel:\n");

  if( opt & BMP_FGAUSS )
  {
    for(sum=0,i=0;i<r;i++) 
    {
      for(j=0;j<r;j++)
      {
	dist = KRN_DIST2( j,i, r/2,r/2 );
	*kp = MtGauss( s*4.0*dist/(float)r );
	sum += *kp;
	printf("%5.3f ",*kp);
	kp++;
      }
      printf("\n");
    }
  }
  if( opt & BMP_FFLAT )
  {
    for(sum=0,i=0;i<r;i++) 
    {
      for(j=0;j<r;j++)
      {
	*kp = 1;
	sum += *kp;
	printf("%5.3f ",*kp);
	kp++;
      }
      printf("\n");
    }
  }
  else if( opt & BMP_FCONE )
  {
    for(sum=0,i=0;i<r;i++) 
    {
      for(j=0;j<r;j++)
      {
	dist = KRN_DIST2( j,i, r/2,r/2 );
	*kp = 1-dist/pow((float)r,s);
	sum += *kp;
	printf("%5.3f ",*kp);
	kp++;
      }
      printf("\n");
    }
  }
  if( opt & BMP_FHAT )
  {
    for(sum=0,i=0;i<r;i++) for(j=0;j<r;j++)
    {
      dist = KRN_DIST2( j,i, r/2,r/2 );
      g = MtGauss( s*4.0*dist/(float)r );
      sum += g;
      *kp = -g;
      kp++;
    }
    sumg = sum;
    kp = info.kernel;
    for(sum=0,i=0;i<r;i++) 
    {
      for(j=0;j<r;j++)
      {
	dstIs0 = (j==r/2 && i==r/2);
	if(dstIs0) *kp = sumg*shrp - *kp;
	printf("%5.3f ",*kp);
	sum+=*kp;
	kp++;
      }
      printf("\n");
    }
  }
  if( opt & BMP_FLAPLACE )
  {
    TRCERRR(("not implemented\n"),MI_ENOTIMPL);
#if 0
    for(sum=0,i=0;i<r;i++) for(j=0;j<r;j++)
    {
      dist = KRN_DIST2( j,i, r/2,r/2 );
      g = MtGauss( s*4.0*dist/(float)r );
      sum += g;
      *kp = -g;
      kp++;
    }
    sumg = sum;
    kp = info.kernel;
    for(sum=0,i=0;i<r;i++) 
    {
      for(j=0;j<r;j++)
      {
	dstIs0 = (j==r/2 && i==r/2);
	if(dstIs0) *kp = sumg*shrp /*+ *kp*/;
	printf("%5.3f ",*kp);
	sum+=*kp;
	kp++;
      }
      printf("\n");
    }
#endif
  }

  if(opt & BMP_FCENTER)
  {
    kp = info.kernel;
    for(i=0;i<r;i++) for(j=0;j<r;j++)
    {
      if((j==r/2 && i==r/2)) 
      {
	sum -= *kp; *kp=0;
      }
      kp++;
    }
  }

  info.sumw = sum;
  printf("sum=%f\n",sum);

  if(opt & BMP_FKERN )
  {
    for(i=0;i<r*r;i++) info.kernel[i] /= info.sumw;
    info.kernel = NULL;
    raiseRc( MI_OK );
  }

  if( channel & BMP_UV )
  {
    if( ! src->dataUV || ! dst->dataUV )
    {
      TRCERRR(("bitmap has no UV channels\n"),MI_ERROR);
    }
    dst->colConv = src->colConv;
    rc = BmpMovWinFilterUV( src, dst, srcWin, dstWin, r,
			    MwfCb_SmoothUV, &info );
    if(rc) raiseRc(rc);

  }

  if( channel & BMP_FLOAT )
  {
    if( ! src->dataFloat[iChannel] )
    {
      TRCERRR(("no data channel\n"),MI_ERROR);
    }
    rc = BmpMovWinFilterFloat( src, dst, srcWin, dstWin, r,
				iChannel,
	  			MwfCb_SmoothGrey, &info );
    if(rc) raiseRc(rc);
  }

  if( channel & BMP_GREY )
  {
    if( ! src->dataGrey )
    {
      TRCERRR(("bitmap has no Grey channels\n"),MI_ERROR);
    }
    rc = BmpMovWinFilterGrey( src, dst, srcWin, dstWin, r,
	  			MwfCb_SmoothGrey, &info );
    if(rc) raiseRc(rc);
  }

  if( channel & BMP_RGB )
  {
    if( ! src->data )
    {
      TRCERRR(("bitmap has no Grey channels\n"),MI_ERROR);
    }
    rc = BmpMovWinFilterRGB( src, dst, srcWin, dstWin, r,
	  			MwfCb_SmoothRGB, &info );
    if(rc) raiseRc(rc);
  }
  

  if( channel & BMP_ALPHA )
  {
    rc = BmpMovWinFilterAlpha( src, dst, channel, channel, 
			  srcWin, dstWin, r,
			  
			  MwfCb_Smooth, &info );
    if(rc) raiseRc(rc);
  }

  FREEMEM(info.kernel);

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/

/*
   Direct fourier transform
*/
int DFT(int dir,int m,double *x1,double *y1)
{
   long i,k;
   double arg;
   double cosarg,sinarg;
   double *x2=NULL,*y2=NULL;

   x2 = MM_malloc(m*sizeof(double),MM_DUMMYTAG);
   y2 = MM_malloc(m*sizeof(double),MM_DUMMYTAG);
   if (x2 == NULL || y2 == NULL)
      return(FALSE);

   for (i=0;i<m;i++) {
      x2[i] = 0;
      y2[i] = 0;
      arg = - dir * 2.0 * 3.141592654 * (double)i / (double)m;
      for (k=0;k<m;k++) {
         cosarg = cos(k * arg);
         sinarg = sin(k * arg);
         x2[i] += (x1[k] * cosarg - y1[k] * sinarg);
         y2[i] += (x1[k] * sinarg + y1[k] * cosarg);
      }
   }

   /* Copy the data back */
   if (dir == 1) {
      for (i=0;i<m;i++) {
         x1[i] = x2[i] / (double)m;
         y1[i] = y2[i] / (double)m;
      }
   } else {
      for (i=0;i<m;i++) {
         x1[i] = x2[i];
         y1[i] = y2[i];
      }
   }

   MM_free(x2);
   MM_free(y2);
   return(TRUE);
}

#if 0
/*-------------------------------------------------------------------------*/   		
/* 								           */
/*-------------------------------------------------------------------------*/
void BmpFilterFFT ( float *data, int nx0, int ny0, int InverseTransform )
{
  int 	  nx,ny,i,j,kx,ky,nxy;
  double f,*real=NULL, *imag=NULL,
  	 *realBuf = NULL,
	 *imagBuf = NULL;

  /**** round dimensions to power of 2 */
  GET_NEAREST_POW2( nx0, nx, kx, 1<<20 );  
  GET_NEAREST_POW2( ny0, ny, ky, 1<<20 );  
  nxy = MAX(nx,ny);
	 
  /**** allocate temp. arrays */
  ALLOCMEM( imag, sizeof(double)*nx*ny );
  ALLOCMEM( real, sizeof(double)*nx*ny );
  ALLOCMEM( realBuf, sizeof(double)*nxy );
  ALLOCMEM( imagBuf, sizeof(double)*nxy );

  /*** transform rows */
  for (i=0;i<ny;i++) 
  {
    for (j=0;j<nx;j++) 
    {
      /*** zero pad input to power of 2 */
      if(( i < ny0 ) && ( j < nx0 ))
      {
	f = MT_MATEL(data,nx0,i,j);
      }
      else f = 0;
      MT_MATEL(real,nx,i,j) = f;
      MT_MATEL(imag,nx,i,j) = 0;
    }
    MtFFT( 1, kx, real+nx*i, imag+nx*i ); 
  }

  /*** transform columns */
  for (j=0;j<nx;j++) 
  {
    for (i=0;i<ny;i++) 
    {
      realBuf[i] = MT_MATEL( real,nx,i,j );
      imagBuf[i] = MT_MATEL( imag,nx,i,j );
    }
    MtFFT( 1, ky, realBuf, imagBuf ); 
    for(i=0;i<ny;i++)
    {
      MT_MATEL( data, nx, i, j ) = 
	sqrt( realBuf[i]*realBuf[i]+ imagBuf[i]*imagBuf[i] );
    }
  }

  FREEMEM(realBuf); FREEMEM(imagBuf);  
  FREEMEM( imag ); FREEMEM( real );
}
#endif

/*-------------------------------------------------------------------------*/   		
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpFilterFFT( BmpBitmap *bs, BmpBitmap *bt, /* source/target bitmaps */ 
        	  int chs, int cht,		/* source/target channels */
    	          BmpRectangle *win,		/* source window */ 
		  int xt, int yt )		/* target point	*/
{
  int 	  nx1,ny1,x,y,xb,yb,x0,y0,nx0,ny0,nx,ny,i,j,kx,ky,nxy,rsSrc,rsDst,
  	  db;
  BmpBitmap *bmpTmp = NULL;
  MtStat stat;
  BmpGreylevel *ppg;
  int		*pp,col;
  float		 gl,f;
  double *real=NULL, *imag=NULL,
  	 *realBuf = NULL,
	 *imagBuf = NULL;

  BMP_GET_WIN( bs, win, x0, y0, nx0, ny0 );

  /**** add border for fourier window "fade-off" */
  db = MAX(2,nx0/10); nx1=nx0+db; ny1=ny0+db;

  /**** round dimensions to power of 2 */
  MT_NEAREST_POW2( nx1, nx, kx, 1<<20 );  
  MT_NEAREST_POW2( ny1, ny, ky, 1<<20 );  
  nxy = MAX(nx,ny);

  /**** allocate temp. bitmap and copy source to it */
  bmpTmp = BmpNewBitmap( nx, ny, BMP_GREY | BMP_CLEAR );
  xb = (nx - nx0)/2;  yb = (ny - ny0)/2;

  BMP_RSCAN_INIT( bs, bs->data, x0, y0, nx0, rsSrc, pp ); 
  BMP_RSCAN_INIT( bmpTmp, bmpTmp->dataGrey, xb, yb, nx0, rsDst, ppg ); 
  BMP_RSCAN_LOOP2( pp, ppg, rsSrc, rsDst, x, y, nx0, ny0 )
  {
    BMP_RGB24_TO_GREY( *pp, gl );
    *ppg = gl;
  }

#if 1  
  /**** TEST: copy  bitmap to destination bitmap */
  BMP_RSCAN_INIT( bmpTmp, bmpTmp->dataGrey, 0, 0, nx, rsSrc, ppg ); 
  BMP_RSCAN_INIT( bt, bt->data, xt+nx, yt, nx, rsDst, pp ); 
  BMP_RSCAN_LOOP2( pp, ppg, rsDst, rsSrc, x, y, nx, ny )
  {
    col = (*ppg);
    *pp = BMP_MKRGB(col,col,col);
  } 
  BMP_MODIFIED( bt );

#endif

  /**** perform 2-dim FFT on tmp. bitmap */
  ALLOCMEM( imag, sizeof(double)*nx*ny );
  ALLOCMEM( real, sizeof(double)*nx*ny );
  ALLOCMEM( realBuf, sizeof(double)*nxy );
  ALLOCMEM( imagBuf, sizeof(double)*nxy );

  /*** transform rows */
  for (i=0;i<ny;i++) 
  {
    for (j=0;j<nx;j++) 
    {
      MT_MATEL(real,nx,i,j) = *BMP_GET_PP_GREY(bmpTmp,j,i);
      MT_MATEL(imag,nx,i,j) = 0;
    }
    MtFFT( 1, kx, real+nx*i, imag+nx*i ); 
  }

  /*** transform columns */
  for (j=0;j<nx;j++) 
  {
    for (i=0;i<ny;i++) 
    {
      realBuf[i] = MT_MATEL( real,nx,i,j );
      imagBuf[i] = MT_MATEL( imag,nx,i,j );
    }
    MtFFT( 1, ky, realBuf, imagBuf ); 
    for(i=0;i<ny;i++)
    {
      *BMP_GET_PP_GREY(bmpTmp,j,i) = 
	sqrt( realBuf[i]*realBuf[i]+ imagBuf[i]*imagBuf[i] );
    }
  }

  /**** copy transformed tmp bitmap to destination bitmap */
  MT_STAT_INIT( &stat );

  BMP_RSCAN_INIT( bmpTmp, bmpTmp->dataGrey, 0, 0, nx, rsSrc, ppg ); 
  BMP_RSCAN_INIT( bt, bt->data, xt, yt, nx, rsDst, pp ); 
  BMP_RSCAN_LOOP2( pp, ppg, rsDst, rsSrc, x, y, nx, ny )
  {
    f = *ppg;
    col = 255.0 * sqrt(f); if(col>255) col=255;
    *pp = BMP_MKRGB(col,col,col);

    MT_STAT_UPDATE( &stat, f );

  } 
  BMP_MODIFIED( bt );

  MT_STAT_RESULT( &stat );
  MtStatPrint( &stat );

  /*** sample 3x3 submatrix (LO quadrant) */
  {
    /*BmpSamplingInfo si;*/
    int u,v,kx = 3, ky=3, mx = nx/2, my=ny/2, npx;
    float sum=0, dx=mx/(float)kx, dy=my/(float)ky;
    MtStat stat;

    TRC(("mx=%d my=%d\n",mx,my));

    MT_STAT_INIT(&stat);
    for(y=0;y+dy<my;y+=dy) 
    {
      for(x=0;x+dx<mx;x+=dx)
      {
	sum = 0; npx=0;
	for(v=y;v<y+dy;v++) for(u=x;u<x+dx;u++) 
	{
	  f = 255.* (*BMP_GET_PP_GREY(bmpTmp,u,v));
	  MT_STAT_UPDATE( &stat, f );
	  sum += f;
	  npx++;
	}
	sum = sum/(float)npx;
	TRC(("%3.0f ",(float)sum));
      }
      TRC(("\n"));
    }
    MT_STAT_RESULT( &stat );
    MtStatPrint(&stat);

  }

  /*** free resources */
  FREEMEM(realBuf); FREEMEM(imagBuf);  
  FREEMEM( imag ); FREEMEM( real );
  BmpCloseBitmap( bmpTmp );

  return MI_OK;
}

/*-------------------------------------------------------------------------*/   		
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpFilterFFT2( BmpBitmap *bs, 		/* source bitmaps */ 
		   BmpBitmap **pbt, 		/* transformed bitmap */
        	   int chs,		        /* source channel */
    	           BmpRectangle *win,		/* source window */ 
		   int fftDir  )		
{
  int 	  nx1,ny1,x,y,xb,yb,x0,y0,nx0,ny0,nx,ny,i,j,kx,ky,nxy,rsSrc,rsDst,
  	  db,j2,i2,nx2,ny2;
  BmpBitmap *bmpTmp = NULL;
  BmpGreylevel *ppg,*ppg2;
  int		*pp;
  float		gl;
  double f,*real=NULL, *imag=NULL,
  	 *realBuf = NULL,
	 *imagBuf = NULL;

  BMP_GET_WIN( bs, win, x0, y0, nx0, ny0 );

  /**** add border for fourier window "fade-off" */
//  db = MAX(2,nx0/10); 
  db=0;
  
  nx1=nx0+db; ny1=ny0+db;

  /**** round dimensions to power of 2 */
  MT_NEAREST_POW2( nx1, nx, kx, 1<<20 );  
  MT_NEAREST_POW2( ny1, ny, ky, 1<<20 );  
  nxy = MAX(nx,ny);

  /**** allocate temp. bitmap and copy source to it */
  bmpTmp = BmpNewBitmap( nx, ny, BMP_GREY | BMP_CLEAR );
  xb = (nx - nx0)/2;  yb = (ny - ny0)/2;

  BMP_RSCAN_INIT( bmpTmp, bmpTmp->dataGrey, xb, yb, nx0, rsDst, ppg ); 

  if( chs == BMP_RGB )
  {
    BMP_RSCAN_INIT( bs, bs->data, x0, y0, nx0, rsSrc, pp ); 
    BMP_RSCAN_LOOP2( pp, ppg, rsSrc, rsDst, x, y, nx0, ny0 )
    {
      BMP_RGB24_TO_GREY( *pp, gl );
      *ppg = gl;
    }
  }
  else if( chs == BMP_GREY )
  {
    BMP_RSCAN_INIT( bs, bs->dataGrey, x0, y0, nx0, rsSrc, ppg2 ); 
    BMP_RSCAN_LOOP2( ppg2, ppg, rsSrc, rsDst, x, y, nx0, ny0 )
    {
      *ppg = *ppg2;
    }
  }

  /**** perform 2-dim FFT on tmp. bitmap */
  ALLOCMEM( imag, sizeof(double)*nx*ny );
  ALLOCMEM( real, sizeof(double)*nx*ny );
  ALLOCMEM( realBuf, sizeof(double)*nxy );
  ALLOCMEM( imagBuf, sizeof(double)*nxy );

  /*** transform rows */
  for (i=0;i<ny;i++) 
  {
    for (j=0;j<nx;j++) 
    {
      MT_MATEL(real,nx,i,j) = *BMP_GET_PP_GREY(bmpTmp,j,i);
      MT_MATEL(imag,nx,i,j) = 0;
    }
    MtFFT( fftDir, kx, real+nx*i, imag+nx*i ); 
  }

  /*** transform columns */
  nx2=nx/2; ny2=ny/2;
  for (j=0;j<nx;j++) 
  {
    for (i=0;i<ny;i++) 
    {
      realBuf[i] = MT_MATEL( real,nx,i,j );
      imagBuf[i] = MT_MATEL( imag,nx,i,j );
    }
    MtFFT( fftDir, ky, realBuf, imagBuf ); 
    for(i=0;i<ny;i++)
    {
       if(i<ny2)
       {
	 if(j<nx2)
	 {
	   j2 = j+nx2; i2 = i+ny2;	// Q00
	 }
	 else
	 {
	   j2 = j-nx2; i2 = i+ny2; //Q01
	 }
       }
       else
       {
	 if(j<nx2)
	 {
	   j2 = j+nx2; i2 = i-ny2;	// Q10
	 }
	 else
	 {
	   j2 = j-nx2; i2 = i-ny2;     //Q11
	 }
	}

       f = sqrt( realBuf[i]*realBuf[i]+ imagBuf[i]*imagBuf[i] );
      *BMP_GET_PP_GREY(bmpTmp,j2,i2) = f;
    }
  }

  /*** free resources */
  FREEMEM(realBuf); FREEMEM(imagBuf);  
  FREEMEM( imag ); FREEMEM( real );
  
  *pbt = bmpTmp;

  return MI_OK;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpInitSampling( BmpBitmap *bmp, BmpRectangle *win, 
    		     float dx, float dy, 
		        /* (avg.) pixel distance in sample grid */
		     int regionSize, /*TODO*/
		        /* leave <regionSize>/2 pixels clearance border */
		     int maxSamples, 
		     	/* if not zero: limit number of samples */
		     unsigned opt,
		        /* BMP_SMP_JITTER (default: uniform) */
		     BmpSamplingInfo *si )
{
  int nx,ny,n,x0,y0,sx,bsx,bsy,sy,dbx=0,dby=0;
  int minMaxSamples=10;

  if( maxSamples && maxSamples < minMaxSamples ) 
    maxSamples = minMaxSamples;

  bsx=bmp->sx; bsy=bmp->sy;
  BMP_GET_WIN( bmp, win, x0, y0, sx, sy );
  if(regionSize)
  {
    dbx=dby=regionSize/2;
    x0 += dbx; y0 += dby;
    sx -= dbx; sy -= dby; 

    if(x0>bsx-1 || y0>bsy-1 || sx > bsx || sy>bsy)
    {
      TRCERR(("invalid bitmap sampling paramaters (bmp=%dx%d)\n",
	    bsx,bsy));
      x0=bsx-1;y0=bsy-1;sx=bsx;sy=bsy;
    }
  }
  
  if( maxSamples )
  {
    dx = dy = sqrt( (float)(sx*sy) / (float)maxSamples );
    nx = (float)sx/dx; ny=(float)sy/dy;
  }
  else
  {
    nx = (float)sx/dx; ny=(float)sy/dy;
  }

  if(nx<1) nx=1; 
  if(ny<1) ny=1;
  n=nx*ny;

  TRC(("InitSampling: %dx%d (max=%d) -> grid %dx%d d=%f,%f\n",
	sx,sy,maxSamples,nx,ny,(float)dx,(float)dy));
  
  si->x0=x0; si->y0=y0;
  si->x2=x0+sx-1; si->y2=y0+sy-1;
  si->x1=x0+dx/2; si->y1=y0+dy/2; 

  si->n=n; si->nx=nx; si->ny=ny;
  si->dx=dx; si->dy=dy;

  if( si->dx<=1 && si->dy<=1 ) FLAG_RESET( opt, BMP_SMP_JITTER);

  if( opt & BMP_SMP_JITTER )
  {
    si->doJittered = TRUE;
    MtRndSSeed( &si->rand, 12345 );
  }
  else si->doJittered = FALSE;

  si->nSmpXY = 0;
  si->pSmpXY = NULL;
  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpCloseSamplingInfo(BmpSamplingInfo *si)
{
  if(si->pSmpXY)
  {
    FREEMEM(si->pSmpXY);
    si->nSmpXY = 0;
  }
  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
MtTab *BmpSampleBitmapOld( BmpBitmap *bmp, int nmax,
    		     	unsigned opt )
{
  int rcThis=0,x,y,k,i,j,dim;
  MtTab *smp=NULL;
  float *vf,vLuv[3];
  BmpSamplingInfo si;

  BmpInitSampling( bmp, NULL, 0, 0, 0, nmax, 0, &si );
 
  dim=0;
  if(opt & BMP_GREY) dim++;
  if(opt & BMP_UV) dim+=2;
  if(opt & BMP_XYZ) dim+=3;
  UT_ASSERT(dim>0&&dim<=4);

  smp=MtCreateTab(si.n,dim);
  UT_ASSERT(smp);

  smp->nRow=0;
  for(i=0;i<si.ny;i++) for(j=0;j<si.nx;j++)
  {
    BMP_GET_SAMPLE( &si, j, i, x, y );
    if(smp->nRow>=smp->nRowMax) 
    {
      TRCERR(("Assertion violation\n"));
      continue;
    }
    vf = smp->rows[smp->nRow++];

    if(opt & BMP_XYZ)
      BMP_GETPIXEL_XYZ(bmp,x,y,vLuv)
    else	
      BMP_GETPIXEL_LUV(bmp,x,y,vLuv);

    if(opt & BMP_GREY) *vf++ = vLuv[0];
    if(opt & BMP_UV) 
    {
      for(k=1;k<3;k++) 
      {
	*vf++ = vLuv[k];
        UT_ASSERT(isfinite(vLuv[k]));
      }
    }
    if(opt & BMP_XYZ) for(k=0;k<3;k++) *vf++ = vLuv[k]; 
  }
 
rcCatch:
  return smp;
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
MtTab *BmpSampleBitmap( BmpBitmap *bmp, int nmax, float sd,
    			CoConverter *cc,
    		     	unsigned opt )
{
  int rcThis=0,x,y,k,i,j,dim;
  MtTab *smp=NULL;
  float *vf,vLuv[3],vRgb[3];
  BmpSamplingInfo si;

  if(sd > MT_REAL_EPSILON)
  {
    BmpInitSampling( bmp, NULL, sd, sd, 0, 0, 0, &si );
  }
  else
  {
    BmpInitSampling( bmp, NULL, 0, 0, 0, nmax, 0, &si );
  }
 
  dim=0;
  if(opt & BMP_GREY) dim++;
  if(opt & BMP_UV) dim+=2;
  if(opt & BMP_XYZ) dim+=3;
  if(opt & BMP_RGB) dim+=3;
  if(opt & BMP_LUV) dim+=3;
  UT_ASSERT(dim>0&&dim<=4);

  smp=MtCreateTab(si.n,dim);
  UT_ASSERT(smp);

  smp->nRow=0;
  for(i=0;i<si.ny;i++) for(j=0;j<si.nx;j++)
  {
    BMP_GET_SAMPLE( &si, j, i, x, y );
    if(smp->nRow>=smp->nRowMax) 
    {
      TRCERR(("Assertion violation\n"));
      continue;
    }
    vf = smp->rows[smp->nRow++];

    if(opt & BMP_RGB)
    {
      BMP_GETPIXEL_RGB_F(bmp,x,y,vf);
    }
    else if((opt & BMP_LUV) && cc)
    {
      BMP_GETPIXEL_RGB_F(bmp,x,y,vRgb);
      cc->fromRgb(vf,vRgb);
    }
    else
    {
      if(opt & BMP_XYZ)
	BMP_GETPIXEL_XYZ(bmp,x,y,vLuv)
      else	
	BMP_GETPIXEL_LUV(bmp,x,y,vLuv);

      if(opt & BMP_GREY) *vf++ = vLuv[0];
      if(opt & BMP_UV) for(k=1;k<3;k++) *vf++ = vLuv[k];
      if(opt & BMP_XYZ) for(k=0;k<3;k++) *vf++ = vLuv[k]; 
    }
  }
 
rcCatch:
  return smp;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpIWriteWinBmp( BmpBitmap *bmp, FILE *fp )
{
  int rcThis=0;
  return rcThis;
}


/*-------------------------------------------------------------------------*/   		
/* 								           */
/*-------------------------------------------------------------------------*/
BmpBitmap *BmpCombineRGB( BmpBitmap *b1, BmpBitmap *b2, BmpBitmap *b3,
    		   	  float *cr, float *cg, float *cb, unsigned opt )
{
  int 		i,j,rcThis=0,c;
  int		*pp1, *pp2, *pp3, *pp,
  		cr1,cg1,cb1,cr2,cg2,cb2,cr3,cg3,cb3,icr,icg,icb;
  BmpBitmap	*bmp=NULL;
  UtProgrInd pgi;

  if(!BMP_DIM_EQ(b1,b2) || (b3 && !BMP_DIM_EQ(b2,b3))) 
    TRCERRR(("bitmap dimension mismatch"),MI_EINVARG);

  bmp = BmpNewBitmap(b1->sx,b1->sy,BMP_RGB);
  if(!bmp) TRCERRR(("BmpNewBitmap() failed\n"),MI_ERROR);

  pp1 = b1->data; pp2 = b2->data;  
  pp3 = b3 ? b3->data : NULL;
  pp = bmp->data;

  UT_PROGRIND_INIT( &pgi );

  for(i=0;i<b1->sy;i++)
  {
    UT_PROGRIND( &pgi, i, 0, b1->sy, "completed: ", NULL );
    for(j=0;j<b1->sx;j++,pp1++,pp2++,pp3 = pp3 ? pp3+1:NULL)
    {
      BMP_PP_GETPIXEL( pp1, c, cr1, cg1, cb1 );
      BMP_PP_GETPIXEL( pp2, c, cr2, cg2, cb2 );
      if(pp3)
      {
        BMP_PP_GETPIXEL( pp3, c, cr3, cg3, cb3 );
      }
      else
      {
	cr3=0;cg3=0;cb3=0;
      }
      icr = (int)((cr[0]*cr1 + cr[1]*cr2 + cr[2]*cr3)+.5);
      icg = (int)((cg[0]*cr1 + cg[1]*cr2 + cg[2]*cr3)+.5);
      icb = (int)((cb[0]*cr1 + cb[1]*cr2 + cb[2]*cr3)+.5);
      *pp++ = BMP_MKRGB(icr,icg,icb);    
    }
  }

  BMP_MODIFIED(bmp);

rcCatch:
  if(rcThis)
  {
    BmpCloseBitmap(bmp);
    bmp=NULL;
  }
  return bmp;
}


/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpCorrectVertLumiBands( 
    	BmpBitmap *bs, BmpBitmap *bt,
	int  y1, int y2, int d, 
	CoConverter *cc, unsigned opt )
{
  int rc,rcThis=0,x,y,icr,icg,icb,k,i,j,m;
  float vRgb[3],vLuv[3],vLuv0[3],sum,sumw,w,lum0;
  UtProgrInd pgi;
  float *refline = NULL, *reflineMAV=NULL, *mask=NULL, span, span0;
  MtStat st;

  if(!BMP_DIM_EQ(bs,bt)) 
    TRCERRR(("bitmap dimension mismatch"),MI_EINVARG);

  m = bs->sx;

  if(d%2==0) d++;

  span0 = 100/255.;

  /*** get LUV channels */
  if( ! bs->dataGrey || ! bs->dataUV )
  {
    rc = BmpRgbToLuv(bs, bs, cc, BMP_GREY|BMP_UV);
    if(rc) raiseRc(rc);
  }

  /*** get reference point & line */
/*
  BMP_GETPIXEL_LUV( bs, x0, y0, vLuv0 );
  TRC1(("luminance reference point at %d,%d = %d\n",
      x0,y0,255*vLuv0[0]));
*/

  ALLOCARR( refline, bs->sx*2 );
  ALLOCARR( reflineMAV, bs->sx*2 );
  ALLOCARR( mask, bs->sx );

  /*** build smoothing mask */
  TRC1(("smoothing delay = %d\n",d));
  for(k=0;k<d;k++)
  {
    mask[k] = MtGauss( 4*fabs(k-d/2)/(float)d );
    printf("%f ",mask[k]);
  }
  printf("\n");
 
  /**** find reference lines */
  if( y1 < 0 )
  {
    y1 = -1; y2 = -1;
    for(y=0;y<bs->sy;y++)
    {
      MT_STAT_INIT( &st );
      for(x=0;x<bs->sx;x++)
      {
	BMP_GETPIXEL_LUV( bs, x, y, vLuv );
	MT_STAT_UPDATE( &st, vLuv[0] );
	refline[x] = vLuv[0];
      }
      MT_STAT_RESULT( &st );
      span = st.max - st.min;
      TRC(("%6d> %6d %6d\n", y, (int)(span*255), (int)(st.var*255)));
      if( span < span0 && st.min > 100./255.)
      {
	if( y1 < 0 )
	{
	  y1 = y;
	}
	else
	{
	  y2 = y;
	}
      }    
      lum0 = st.avg;
      vLuv0[0] = lum0;
    }
    TRC1(("referene strip: y1 = %d y2 = %d\n",y1,y2));
  }
  if(y1 < 0 || y2 < 0 )
  {
    TRCERRR(("could not found reference lines with luminance diff < %d\n",
	  (int)span0*255),MI_ERROR);
  }

  MtVecConst( refline, m, 0 );
  for(k=0,y=y1;y<=y2;y++,k++)
  {
    for(x=0;x<bs->sx;x++)
    {
      BMP_GETPIXEL_LUV( bs, x, y, vLuv );
      refline[x] += vLuv[0];
    }   
  }
  MtVecScale_f( refline, m, 1./(float)k);

  MT_STAT_INIT( &st );
  for(x=0;x<bs->sx;x++)
  {
/*    printf("%d ",255.0*refline[x]);*/
    MT_STAT_UPDATE( &st, 255.*refline[x] );
  }
/*  printf("\n");*/
  MT_STAT_RESULT( &st );
  MtStatPrint(&st);
  lum0 = st.avg/255.;

/*  
  MtVecMAE( reflineMAV, refline, bs->sx, beta );
*/
  UT_PROGRIND_INIT( &pgi );
  for(x=0;x<bs->sx;x++)
  {
    UT_PROGRIND( &pgi, x, 0, bs->sx-1, "smoothing: ", NULL );
    for(k=0,sum=0,sumw=0,j=x-d/2;j<=x+d/2;j++)
    {
      i=j;
      if( i < 0 ) i = -i;
      if( i >= bs->sx ) i =  bs->sx - (i - bs->sx + 1);
      w = mask[k++];
      sumw += w;
      sum += refline[i] * w;
    }    
    reflineMAV[x] = sum/sumw;
  }

  /*** adjust bitmap */
  UT_PROGRIND_INIT( &pgi );
  for(y=0;y<bs->sy;y++)
  {
    UT_PROGRIND( &pgi, y, 0, bs->sy-1, "completed: ", NULL );
    for(x=0;x<bs->sx;x++)
    {
      BMP_GETPIXEL_LUV( bs, x, y, vLuv );
      vLuv[0] += (lum0 - reflineMAV[x]);
      cc->toRgb( vRgb, vLuv  ); 
      CO_RGB2INT( vRgb, icr,icg,icb );
      BMP_SETPIXEL( bt, x, y, icr, icg, icb ); 
    }
  }
  BMP_MODIFIED( bt );

rcCatch:

  FREEMEM( refline );
  FREEMEM( mask );
  FREEMEM( reflineMAV );
  return rcThis;
}


/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpEnhanceDiff( BmpBitmap *bs, BmpBitmap *bt, float alpha,
    		    CoConverter *cc )
{
  int rcThis=0,x,y,icr,icg,icb;
  UtProgrInd pgi;
  float vLuv[3],vRgb[3],lum,ls,lt;

  if(!BMP_DIM_EQ(bs,bt)) 
    TRCERRR(("bitmap dimension mismatch"),MI_EINVARG);

/*    bt = bt + alpha*(bs-bt) */

  if(!bs->dataGrey)
  {
    BmpRgbToLuv( bs, bs, cc, BMP_GREY );
  }
  if(!bt->dataGrey)
  {
    BmpRgbToLuv( bt, bt, cc, BMP_GREY );
  }

  UT_PROGRIND_INIT( &pgi );
  MT_VEC3_SET( vLuv, 0, 0, 0 );
  for(y=0;y<bs->sy;y++)
  {
    UT_PROGRIND( &pgi, y, 0, bs->sy-1, "completed: ", NULL );
    for(x=0;x<bs->sx;x++)
    {
      ls = BMP_GETPIXEL_GREY(bs,x,y);
      lt = BMP_GETPIXEL_GREY(bt,x,y);
      lum = lt + alpha*(ls-lt);
      if(lum<0) lum = 0; 
      if(lum>255) lum = 255;
      vLuv[0] = lum/255.;
      cc->toRgb( vRgb, vLuv  ); 
      CO_RGB2INT( vRgb, icr,icg,icb );
      BMP_SETPIXEL( bt, x, y, icr, icg, icb ); 
    }
  }
  BMP_MODIFIED(bt);
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
int BmpEnhanceDiff2( BmpBitmap *bs, BmpBitmap *bt, 
    		     BmpRectangle *win, /* target window */
		     int dx1, int dx2,	/* background win offset/width */
    		     float alpha,
    		     CoConverter *cc )
{
  int rcThis=0,x,y,icr,icg,icb,i,j,k;
  UtProgrInd pgi;
  float vLuv[3],vRgb[3],lum,ls,l0,sumw,w,df,df0;
  float *col0=NULL,*smask2=NULL,*smask1=NULL;

  if(!BMP_DIM_EQ(bs,bt)) 
    TRCERRR(("bitmap dimension mismatch"),MI_EINVARG);

  ALLOCARR( col0, 2*win->sy );
  ALLOCARR( smask1, 2*dx1 );
  ALLOCARR( smask2, 2*dx2 );
/*    bt = bt + alpha*(bs-bt) */

  if(!bs->dataGrey)
  {
    BmpRgbToLuv( bs, bs, cc, BMP_GREY );
  }
  if(!bt->dataGrey)
  {
    BmpRgbToLuv( bt, bt, cc, BMP_GREY );
  }

  MT_VEC3_SET(vRgb,1,1,1);
  cc->fromRgb(vLuv,vRgb);

  /*** build smoothing masks */
  for(k=0;k<dx1;k++)
  {
    smask1[k] = MtGauss( 1.3*(float)k/(float)dx1 );
    printf("%f ",smask1[k]);
  }
  printf("\n");
  for(k=0;k<dx2;k++)
  {
    smask2[k] = MtGauss( 1.3*(float)k/(float)dx2 );
    printf("%f ",smask2[k]);
  }
  printf("\n");

  /**** get background ref. column */
  for(y=win->y0;y<win->y0+win->sy;y++)
  {
    i = y-win->y0;
    col0[i] = 0;
    sumw = 0;
    for(j=0;j<dx1;j++)
    {
      x = win->x0-j-1;
      w = smask1[j];
      col0[i] += w * (BMP_GETPIXEL_GREY(bs,x,y));
      sumw += w;
    }
    for(j=0;j<dx2;j++)
    {
      x = win->x0+win->sx+j;
      w = smask2[j];
      col0[i] += w * (BMP_GETPIXEL_GREY(bs,x,y));
      sumw += w;
    }
    col0[i] /= sumw;
  }

  /*** make difference image */
  UT_PROGRIND_INIT( &pgi );
  for(y=win->y0;y<win->y0+win->sy;y++)
  {
    UT_PROGRIND( &pgi, y, win->y0, win->y0+win->sy, "completed: ", NULL );
    i = y-win->y0;
    l0 = col0[i];
    for(x=win->x0;x<win->x0+win->sx;x++)
    {
      ls = BMP_GETPIXEL_GREY(bs,x,y);
      df0 = ls-l0;
      df = MtSign(df0)*(df0/3)*(df0/3);
      if((fabs(df/df0)) > alpha) df = alpha*df0;
      lum = ls+df;
/*      lum = ls + alpha*df;*/
      if(lum<0) lum = 0; 
      if(lum>255) lum = 255;
      vLuv[0] = lum/255.;
      cc->toRgb( vRgb, vLuv  ); 
      CO_RGB2INT( vRgb, icr,icg,icb );
      BMP_SETPIXEL( bt, x, y, icr, icg, icb ); 
    }
  }
  BMP_MODIFIED(bt);

rcCatch:
  FREEMEM( col0 );
  FREEMEM( smask1 );
  FREEMEM( smask2 );
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpMakeXYZ( BmpBitmap *bs, BmpBitmap *bt, int colmod,
    		int linearGreyMap )
{
  int i,x,y,rcThis=0,cr,cg,cb,c,k,doRGB2RGB=0,doLUV2LUV=0;
  float	*px,*py,*pz,vXYZ[3],*pg=NULL;

  USE(k);
  MT_VEC3_SET(vXYZ,0,0,0);

  if( ! BMP_DIM_EQ( bs, bt ) )
    TRCERRR(("bitmap dimension differ"),MI_ERROR);

  if(colmod == BMP_RGB)
  {
    doRGB2RGB = TRUE;
    if(linearGreyMap)
    {
      pg = bt->dataGrey;
      if(!pg)
        TRCERRR(("no greay channel\n"),MI_EINVARG);
    }
  }
  else
  {
    doLUV2LUV = TRUE;
  }

  for(i=0;i<3;i++) BmpGetFloatChannel(bt,i,0);

  px = bt->dataFloat[0];
  py = bt->dataFloat[1];
  pz = bt->dataFloat[2];

  for(y=0;y<bt->sy;y++) for(x=0;x<bt->sx;x++)
  {
    if( doLUV2LUV )
    {
      BMP_GETPIXEL_LUV( bs, x, y, vXYZ ); 
      vXYZ[0] = 255.*vXYZ[0];
      vXYZ[1] = 255.*(vXYZ[1]+1)/2.;
      vXYZ[2] = 255.*(vXYZ[2]+1)/2.;
    }
    else if( doRGB2RGB )
    {
      BMP_GETPIXEL( bs, x, y, c, cr, cg, cb ); 
      MT_VEC3_SET( vXYZ, cr, cg, cb);
#if XYZRGB_GAMMA
      for(k=0;k<3;k++) vXYZ[k] /= 255.0;
      CoGammaRgbToLinear( vXYZ );
      for(k=0;k<3;k++) vXYZ[k] *= 255.0;
#endif
      if(linearGreyMap)
      {
        *pg++ = RGB_TO_GREY(vXYZ[0],vXYZ[1],vXYZ[2]);
      }
    }
    *px++ = vXYZ[0]; *py++ = vXYZ[1]; *pz++ =vXYZ[2];
  }

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpXYZToRGB( BmpBitmap *bs, BmpBitmap *bt, int colmod,
    		 float *meanXYZ,
    		 CoConverter *cc )
{
  int k,rcThis=0,cr,cg,cb,x,y;
  float *pxyz[3],vXYZ[3],vRgb[3];

  if( ! BMP_DIM_EQ( bs, bt ) )
    TRCERRR(("bitmap dimension differ"),MI_ERROR);

  for(k=0;k<3;k++) 
    if(!(pxyz[k] = bs->dataFloat[k])) 
      TRCERRR(("no float channel\n"),MI_ERROR);

  for(y=0;y<bs->sy;y++) for(x=0;x<bs->sx;x++)
  {
    for(k=0;k<3;k++) vXYZ[k] = *(pxyz[k])++;
    if(meanXYZ) for(k=0;k<3;k++) vXYZ[k] += meanXYZ[k];
    if( colmod == BMP_RGB )
    {
      for(k=0;k<3;k++) vXYZ[k] /= 255.0;
#if XYZRGB_GAMMA
      CoGammaRgbFromLinear( vXYZ );
#endif
      cr = vXYZ[0]*255; cg = vXYZ[1]*255; cb = vXYZ[2]*255; 
    }
    else if( colmod == BMP_LUV )
    {
      vXYZ[0] = vXYZ[0]/255;
      vXYZ[1] = 2.*vXYZ[1]/255.-1;
      vXYZ[2] = 2.*vXYZ[2]/255.-1.;
      cc->toRgb(vRgb,vXYZ);
      CO_RGB2INT( vRgb, cr,cg,cb );  
    }
    else TRCERR(("assertion failed %d\n",colmod));

    BMP_SETPIXEL( bt, x, y, cr, cg, cb );
  }

rcCatch:
  return rcThis;
}



/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpXYZToLUV( BmpBitmap *bs, BmpBitmap *bt, int colmod,
    		 CoConverter *cc )
{
  int k,rcThis=0,cr,cg,cb,x,y;
  float *pxyz[3],vXYZ[3],vRgb[3];

  if( ! BMP_DIM_EQ( bs, bt ) )
    TRCERRR(("bitmap dimension differ"),MI_ERROR);

  for(k=0;k<3;k++) 
    if(!(pxyz[k] = bs->dataFloat[k])) 
      TRCERRR(("no float channel\n"),MI_ERROR);

  for(y=0;y<bs->sy;y++) for(x=0;x<bs->sx;x++)
  {
    for(k=0;k<3;k++) vXYZ[k] = *(pxyz[k])++;
    if( colmod == BMP_RGB )
    {
      cr = vXYZ[0]; cg = vXYZ[1]; cb = vXYZ[2]; 
      MT_VEC3_SET( vRgb, cr/255., cg/255., cb/255.);
#if XYZRGB_GAMMA
      CoGammaRgbFromLinear( vRgb );
#endif
      cc->fromRgb(vXYZ,vRgb);
    }
    else if( colmod == BMP_LUV )
    {
      vXYZ[0] = vXYZ[0]/255;
      vXYZ[1] = 2.*vXYZ[1]/255.-1;
      vXYZ[2] = 2.*vXYZ[2]/255.-1.;
    }
    else TRCERR(("assertion failed %d\n",colmod));

    BMP_SETPIXEL_LUV(bt,x,y,vXYZ);
  }

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpRGB2XYZ( BmpBitmap *bs, BmpBitmap *bt, CoConverter *cc, unsigned opt )
{
  int rcThis = 0,i,x,y,cr,cg,cb,c;
  float	vXYZ[3],vRgb[3],*px,*py,*pz;
  if(!bs || !bt) TRCERRR(("unexpected NULL pointer"),MI_ERROR);

  if( ! BMP_DIM_EQ( bs, bt ) )
    TRCERRR(("bitmap dimension differ"),MI_ERROR);

  for(i=0;i<3;i++) BmpGetFloatChannel(bt,i,0);
  px = bt->dataFloat[0];
  py = bt->dataFloat[1];
  pz = bt->dataFloat[2];

  for(y=0;y<bt->sy;y++) for(x=0;x<bt->sx;x++)
  {
      BMP_GETPIXEL( bs, x, y, c, cr, cg, cb ); 
      MT_VEC3_SET( vRgb, cr/255., cg/255., cb/255.);
      cc->fromRgb(vXYZ,vRgb);
      *px++ = vXYZ[0]; *py++ = vXYZ[1]; *pz++ =vXYZ[2];  
  }
  BMP_MODIFIED(bt);
  
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpXYZ2RGB( BmpBitmap *bs, BmpBitmap *bt, CoConverter *cc, unsigned opt )
{
  int rcThis = 0,x,y,cr,cg,cb;
  float	vXYZ[3],vRgb[3],*px,*py,*pz;
  if(!bs || !bt) TRCERRR(("unexpected NULL pointer"),MI_ERROR);

  if( ! BMP_DIM_EQ( bs, bt ) )
    TRCERRR(("bitmap dimension differ"),MI_ERROR);

  BmpGetRGBChannel(bt,0);
  px = bs->dataFloat[0];
  py = bs->dataFloat[1];
  pz = bs->dataFloat[2];

  for(y=0;y<bt->sy;y++) for(x=0;x<bt->sx;x++)
  {
      vXYZ[0] = *px++; vXYZ[1] = *py++;vXYZ[2] = *pz++;  
      cc->toRgb(vRgb,vXYZ);
      CO_RGB2INT( vRgb, cr,cg,cb );
      BMP_SETPIXEL( bt, x, y, cr, cg, cb ); 
  }
  BMP_MODIFIED(bt);
  
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpDownSamplePyraLevelFast( 
	  BmpBitmap *bmpLo, float *dataLo, BmpUshort *dataLoUs, uint8_t *dataLoUb,
	  BmpBitmap *bmpHi, float *dataHi, BmpUshort *dataHiUs, uint8_t *dataHiUb,
	  int oxh, int oyh, // optional bmpHi offset 0, or -1 (for MCM)  
          float *kernel5x5,  // optional 5x5 downsampling kernel
	  int btyp	 // boundary handling
	  )
{
#define FAST_DOWNSAMPLE_LINE(gtype, dataLo, dataHi)\
    {\
    gtype *pDataLo,*pDataHi,*p;\
    float g0,g1,g2,g3,g4,\
	  g5,g6,g7,g8,g9,\
	  g10,g11,g12,g13,g14,\
	  g15,g16,g17,g18,g19,\
	  g20,g21,g22,g23,g24;\
\
      pDataLo = dataLo + BMP_XYOFF(bmpLo,xl,yl);\
      pDataHi = dataHi + BMP_XYOFF(bmpHi,xh-r,yh-r);\
      while(xh<=xhmax-4)\
      {\
	p = pDataHi;\
\
	g0=*p++;  g1=*p++;  g2=*p++;  g3=*p++;  g4=*p++; p+=rmod;\
	g5=*p++;  g6=*p++;  g7=*p++;  g8=*p++;  g9=*p++; p+=rmod;\
	g10=*p++; g11=*p++; g12=*p++; g13=*p++; g14=*p++; p+=rmod;\
	g15=*p++; g16=*p++; g17=*p++; g18=*p++; g19=*p++; p+=rmod;\
	g20=*p++; g21=*p++; g22=*p++; g23=*p++; g24=*p++; \
\
	g = w0*(g0+g4+g20+g24) +\
	    w1*(g1+g3+g21+g23+g5+g9+g15+g19) +\
	    w2*(g2+g22+g10+g14) +\
	    w4*(g6+g8+g16+g18) +\
	    w5*(g7+g17+g11+g13) +\
	    w8*(g12);\
\
	*pDataLo++ = g;\
	pDataHi += 2;\
	xh +=2;\
	xl++;\
      }\
    }
  
  int rcThis=0,r=2,rx,ry,xl,yl,roff,rmod,
      xlmax=bmpLo->sx-1,ylmax=bmpLo->sy-1,
      xhmax=bmpHi->sx-1,yhmax=bmpHi->sy-1;
  float *krn=kernel5x5 ? kernel5x5 : BmpGetDfltPyraConvKernel5x5(),*pkrn;
  double wsum,w0,w1,w2,w3,w4,w5,w6,w7,w8;
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

  UT_ASSERT(dataHi||dataHiUs||dataHiUb);
  UT_ASSERT(dataLo||dataLoUs||dataLoUb);

  if(btyp==BMP_BDEFAULT) btyp=BMP_BCLIP;

  r = 2;

  TRC3(("BmpDownSamplePyraLevelFast: (%dx%d) -> (%dx%d)\n",
	bmpHi->sx,bmpHi->sy,bmpLo->sx,bmpLo->sy));

  wsum=0; pkrn=krn;
  for(ry=-r;ry<=r;ry++) for(rx=-r;rx<=r;rx++) wsum += *pkrn++;
  UT_ASSERT(pkrn-krn==(2*r+1)*(2*r+1));

  // normalized quadrant-symmetrical kernel weights
  w0=krn[0]/wsum; w1=krn[1]/wsum; w2=krn[2]/wsum;
  w3=krn[5]/wsum; w4=krn[6]/wsum; w5=krn[7]/wsum;
  w6=krn[10]/wsum; w7=krn[11]/wsum; w8=krn[12]/wsum;

  TRC3(("w0=%g w1=%g w2=%g \n" 
	"w3=%g w4=%g w5=%g \n" 
	"w6=%g w7=%g w8=%g \n",w0,w1,w2,w3,w4,w5,w6,w7,w8));

  UT_ASSERT(fabs(w2-w6)<MT_NUM_EPS && 
            fabs(w7-w5)<MT_NUM_EPS &&
            fabs(w3-w1)<MT_NUM_EPS);

  roff = BMP_XYOFF(bmpHi,r,r);
  rmod = bmpHi->sx-(2*r+1);

  UT_ASSERT(((xhmax+1) <= 2*(xlmax+1))&&((yhmax+1) <= 2*(ylmax+1)));

  if((do_time = (LogLevel>=2)&&
		((LongInt)bmpHi->sx*bmpHi->sy>250000000))==TRUE)
  {
    TIMER_START(&tm);
  }

  //
  //
  //
#pragma omp parallel for private(xl,yl)
  for(yl=0;yl<bmpLo->sy;yl++)
  {
    int dx,dy,xh2,yh2,xh,yh,linedone,fetch;
    LongInt po;
    double wsum,gsum,g=0;
    float *pkrn,wgt;

    yh = yl*2+oyh; xh=oxh;
    linedone = FALSE;
    for(xl=0;xl<bmpLo->sx;
	)
    {
      xh = 2*xl+oxh;
      if((!(BMP_IN_RANGEB(xh-r,yh-r,bmpHi)&&
            BMP_IN_RANGEB(xh+r,yh+r,bmpHi)))||linedone)
      {
	// 
	// handle border case 
	//
	pkrn=krn; wsum=gsum=0;
	for(dy=-r;dy<=r;dy++)
	  for(dx=-r;dx<=r;dx++)
	  {
	    wgt = *pkrn++;
	    xh2=xh+dx; yh2=yh+dy;
	    fetch = TRUE;
            if(btyp==BMP_BCLIP)
	    {
	      if(!(BMP_IN_RANGEB(xh2,yh2,bmpHi))) continue;
	    }
            else if(btyp==BMP_BCONST)
	    {
	      BMP_CLIP_CONST(0,0,xhmax,yhmax,xh2,yh2);
	    }
            else if(btyp==BMP_BZERO)
	    {
	      if(!(BMP_IN_RANGEB(xh2,yh2,bmpHi)))
	      {
	        g = 0;
	        fetch = FALSE;
	      }
	    }
	    else 
	    {
	      TRCERR(("unexpected boundary type=%d\n",btyp));
	      fetch = FALSE;
	    }

	    if(fetch)
	    {
	      po = BMP_XYOFF(bmpHi,xh2,yh2);
	      g = dataHi ? dataHi[po] : (dataHiUs ? dataHiUs[po] : 
		  				    dataHiUb[po]);
	    }
	    gsum += wgt*g;
	    wsum += wgt;
	  }
	g = gsum/wsum;
	po = BMP_XYOFF(bmpLo,xl,yl);
	if(dataLo)
	  dataLo[po] = g;
	else if(dataLoUs)
	  dataLoUs[po] = g;
	else
	  dataLoUb[po] = g;
	xl++; xh+=2;
	continue;
      }

      // 
      // non-border case: process fast line
      //
      if(dataHi)
      {
	FAST_DOWNSAMPLE_LINE(float, dataLo, dataHi );
      }
      else if(dataHiUs )
      {
	FAST_DOWNSAMPLE_LINE(BmpUshort, dataLoUs, dataHiUs );
      }
      else
      {
	FAST_DOWNSAMPLE_LINE(uint8_t, dataLoUb, dataHiUb );
      }
      linedone = TRUE;
    }
  }
  
  if(do_time)
  {
    TIMER_STOP(&tm);
    TRC(("CPU = %.1f seconds (%.0f pix/sec)\n", 
	  (double)(TIMER_DIFF_MS(&tm)/1000.0), 
	  (bmpHi->sx*bmpHi->sy)/(double)(TIMER_DIFF_MS(&tm)/1000.0)));
  }

rcCatch:
    return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpDownSamplePyraLevel( BmpBitmap *bmp0,
    		       BmpBitmap *bmpLo, float *dataLo, double *dataLo_d,
    		       BmpBitmap *bmpHi, float *dataHi, double *dataHi_d,
		       int	doPresVariance )
{
  int rcThis=0,x,y,i,
  	x2max=bmpHi->sx-1,y2max=bmpHi->sy-1,
	r,krn_r=bmp0?bmp0->pyramid.krn_r:2,nPix;
  float *krn=bmp0?bmp0->pyramid.krn:BmpGetDfltPyraConvKernel5x5();
  double g,g0,res_factor=2;
  double eClp,gsumHi,gsumLo,gsumHi2,gsumLo2;
  MtStat stLo,stHi;

  TRC3(("BmpDownSamplePyraLevel: (%dx%d) -> (%dx%d)\n",
	bmpHi->sx,bmpHi->sy,bmpLo->sx,bmpLo->sy));
  
  stHi.var=0;

  if(!((dataLo||dataLo_d) && (dataHi||dataHi_d)))
  {
    TRCERRR(("no bitmap data\n"),MI_ERROR);
  }

  if(doPresVariance)
  {
    MT_STAT_INIT(&stHi);
    gsumHi = gsumHi2 = 0;
    nPix = bmpHi->sx*bmpHi->sy;
    for(i=0;i<nPix;i++)
    {
      g = dataHi ? dataHi[i] : dataHi_d[i];
      gsumHi += g; gsumHi2 +=g*g;
    }
    MT_STAT_UPDATE_N(&stHi,nPix,-1,-1,gsumHi,gsumHi2);
    MT_STAT_RESULT(&stHi);
  }

  MT_STAT_INIT(&stLo); 
  gsumLo = gsumLo2 = 0;
  r = krn_r;

  nPix = bmpLo->sx*bmpLo->sy;

#pragma omp parallel for private(y,x,g) \
  reduction(+:gsumLo,gsumLo2)
  for(y=0;y<bmpLo->sy;y++)
  {
    int dx,dy,po,x2,y2,x0,y0;
    double wsum,gsum;
    float *pkrn,wgt;
    for(x=0;x<bmpLo->sx;x++)
    {
      x0=x*res_factor; y0=y*res_factor;
      pkrn=krn; wsum=gsum=0;
      for(dy=-r;dy<=r;dy++)
	for(dx=-r;dx<=r;dx++)
	{
          wgt = *pkrn++;
          x2=x0+dx; y2=y0+dy;
#if GP_HDL_BORD_CLIP
	  if(!(BMP_IN_RANGE(x2,y2, 0,0,x2max,y2max))) continue;
#else
	  BMP_CLIP_CONST(0,0,x2max,y2max, x2,y2);
#endif
          po = BMP_GETPIXEL_OFFS(bmpHi,x2,y2);
	  g = dataHi ? dataHi[po] : dataHi_d[po];
	  gsum += wgt*g;
	  wsum += wgt;
	}
      //UT_ASSERT(wsum > MT_REAL_EPSILON);
      g = gsum/wsum;
      if(dataLo)
        dataLo[BMP_GETPIXEL_OFFS(bmpLo,x,y)] = g;
      else
        dataLo_d[BMP_GETPIXEL_OFFS(bmpLo,x,y)] = g;
      if(doPresVariance)
      {
	gsumLo += g; gsumLo2 += g*g;
      }
    }
  }
  BMP_MODIFIED(bmpLo);
  if( doPresVariance )
  {
    int i,nClp;
    float sHi,sLo,mHi,mLo,f;
    MT_STAT_UPDATE_N(&stLo,nPix,-1,-1,gsumLo,gsumLo2);
    MT_STAT_RESULT(&stLo); 

    sLo = stLo.var; sHi = stHi.var;
    mLo = stLo.avg; mHi = stHi.avg;
    f = sLo > MT_REAL_EPSILON ? sHi/sLo : 1;

    if(!(f>0.99&&f<1.2))
    {
      TRC1(("**** WARNING: unexpected variance factor: %g\n",f));
    }

    nClp=0;eClp=0;
    for(i=0;i<bmpLo->sx*bmpLo->sy;i++)
    {
      g0 = dataLo ? dataLo[i] : dataLo_d[i];
      g = f*(g0-mLo)+mLo;
      if(g<0)
      {
//	printf("%6.2f > %6.2f\n",dataLo[i],g); 
	g=0;nClp++;eClp+=fabs(g-g0);
      }
      if(g>255)
      {
//	printf("%6.2f > %6.2f\n",dataLo[i],g); 
	g=255;nClp++;eClp+=fabs(g-g0);
      }
      if(dataLo) 
	dataLo[i] = g;
      else
	dataLo_d[i] = g;
    }
    eClp /= (double)nPix;
    if(eClp > 2) 
    {
      TRC1(("**** WARNING: Clipped %d of %d pixels (err/pix=%g)\n",
	    nClp,nPix,eClp));
    }
  }


rcCatch:
    return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpIReducePyramid( BmpBitmap *bmp0,
    		       int level,
		       unsigned opt )
{
  int rcThis=0,rc,i,nChanDone=0;
  int sx,sy, doPresVar = 0;
  BmpBitmap *bmpParent = bmp0->pyramid.gaussian[level-1],
	    *bmp=NULL;
  float *dataHi, *dataLo;

  // Recursion end ?
  sx = bmpParent->sx;
  sy = bmpParent->sy;
/*  if(sx==1||sy==1) 
  {
    raiseRc(MI_OK);
  }*/
  if(sx==1&&sy==1) 
  {
    raiseRc(MI_OK);
  }

  sx = (sx+1) / 2;
  sy = (sy+1) / 2;

  // New Gaussian
  bmp = BmpNewBitmap(sx,sy,opt);
  if(!bmp) raiseRc(MI_ERROR);

  bmp->vsx = MAX(bmpParent->vsx/2,1);
  bmp->vsy = MAX(bmpParent->vsy/2,1);
  bmp->vx0 /= 2;
  bmp->vy0 /= 2;

#if 1
  {
    int k,sx2,sy2;
    MT_NEAREST_POW2( bmp->vsx, sx2, k, 1<<28 );  
    MT_NEAREST_POW2( bmp->vsy, sy2, k, 1<<28 );  
    UT_ASSERT(sx2==bmp->vsx);
    UT_ASSERT(sy2==bmp->vsy);
  }
#endif

  BmpCloseBitmap(bmp0->pyramid.gaussian[level]);
  bmp0->pyramid.gaussian[level] = bmp;
  bmp0->pyramid.n_level = level+1;
  bmp->pyramid.res_act = bmpParent->pyramid.res_act * 2;

  TRC(("Bmp pyramid: %2d %7.1f   %dx%d (%dx%d)\n",
	level,bmp->pyramid.res_act,sx,sy,bmp->vsx,bmp->vsy));

  if(opt & BMP_GREY)
  {
    dataHi = bmpParent->dataGrey;
    dataLo = BmpGetGreyChannel(bmp,0);
    rc = BmpDownSamplePyraLevelFast( 
		bmp, dataLo, NULL, NULL,
		bmpParent, dataHi, NULL, NULL,
		0, 0, NULL, BMP_BDEFAULT );
    if(rc) raiseRc(MI_ERROR);
    nChanDone++;
  }
  if(opt & BMP_ALPHA)
  {
    dataHi = bmpParent->dataAlpha;
    dataLo = BmpGetAlphaChannel(bmp,0);
    rc = BmpDownSamplePyraLevel( 
		bmp0, bmp, dataLo, NULL, 
		bmpParent, dataHi, NULL,
		doPresVar );
    if(rc) raiseRc(MI_ERROR);
    nChanDone++;
  }
  if(opt & BMP_Z)
  {
    dataHi = bmpParent->dataZ;
    dataLo = BmpGetZChannel(bmp,0);
    rc = BmpDownSamplePyraLevel( 
		bmp0, bmp, dataLo, NULL, 
		bmpParent, dataHi, NULL,
		doPresVar );
    if(rc) raiseRc(MI_ERROR);
    nChanDone++;
  }
  if(opt & BMP_XYZ)
  {
    for(i=0;i<3;i++)
    {
      dataHi = bmpParent->dataFloat[i];
      dataLo = BmpGetFloatChannel(bmp,i,0);
      rc = BmpDownSamplePyraLevel( 
	  	bmp0, bmp, dataLo, NULL, 
		bmpParent, dataHi, NULL, 
	  	doPresVar );
      if(rc) raiseRc(MI_ERROR);
      nChanDone++;
    }
  }
  if(!nChanDone&&(!(opt&BMP_GEN)))
  {
    TRCERRR(("no bitmap channel selected"),MI_ERROR);
  }

  rc = BmpIReducePyramid( bmp0, level+1, opt );
  if(rc) raiseRc(MI_ERROR);

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
//static inline 
float BmpLaplacianGain_FUNC1( 
    	float x, // input laplacian coefficient
	float s, // x scale factor
	float g, // overall gain factor
	float p, // power exponent
	float l, // lower limit of x*s for nonlinear enhancement
	float h  // upper limit of x*s for nonlinear enhancement
    )
{
  // see [Dippel2002]
  int sig;
  float y,xs,xs_;

  xs = x * s;


  if(xs<0)
  {
    sig = -1;
    xs_ = -xs;
  }
  else 
  {
    xs_ = xs;
    sig = 1;
  }

  if(xs_<l)
  {
    y = g*h*(xs/l)*powf(l/h,p);
  }
  else if(xs_<h)
  {
    y = g*h*sig*powf(xs_/h,p);
  }
  else
  {
    y = g*x;
  }
  return y/s;
}

typedef struct
{
  float *lut_table;
  double lut_scale,lut_offs;
  int    lut_n;	  // lookup table half size
  double  s,  // input scale factor
	  g,  // overall gain factor
	  p,  // exponent
	  l;  // lower limit for nonlinear enhancement (upper=1)
}
BmpDetailGainFunc;

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpDeleteDetailGainFunc( BmpDetailGainFunc **p_adm )
{
  int rcThis=0;
  BmpDetailGainFunc *adm=*p_adm;

  if(adm)
  {
    FREEMEM(adm->lut_table);
    FREEMEM(adm);
    *p_adm = NULL;
  }

//rcCatch: 
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static 
float BmpComputeDetailGainFunc(BmpDetailGainFunc *adm, double s, float x)
{
  if(adm->lut_table)
  {
    int i = x*adm->lut_scale+adm->lut_offs;
    if(i>=0&&i<2*adm->lut_n)
    {
      return adm->lut_table[i];
    }
  }

  return 
    BmpLaplacianGain_FUNC1(
	x*(s<MT_NUM_EPS?adm->s:s),
	1,adm->g,adm->p,adm->l,1.0) / adm->s;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BmpDetailGainFunc *BmpCreateDetailGainFunc(
		      int typ,	// function type
		      double s,	// input scale factor
		      double g,	// overall gain factor
		      double p,	// exponent
		      double l,	// lower limit for nonlinear enhancement (upper=1)
		      int lut_size // lookup table size (half)
    		      )
{
  int rcThis=0,i;
  BmpDetailGainFunc *adm=NULL;
  float *lut_table=NULL;

  UT_ASSERT(typ==1); // TODO: support other functions
  ALLOCOBJ0(adm);
  adm->lut_n = lut_size;
  adm->s = s;
  UT_ASSERT(s>0);
  adm->g = g;
  adm->p = p;
  adm->l = l;
  adm->lut_scale = adm->lut_n*adm->s;
  adm->lut_offs = adm->lut_n;

  if(adm->lut_n>0)
  {
    ALLOCARR(lut_table,adm->lut_n*2);
    for(i=0;i<2*adm->lut_n;i++)
    {
      double xs = (i/(double)adm->lut_n-1.0),
	     f =  BmpComputeDetailGainFunc(adm,1,xs);
      UT_ASSERT(isfinite(f));
      lut_table[i]=f;
    }
    adm->lut_table = lut_table; lut_table = NULL;
  }

rcCatch:
  FREEMEM(lut_table);
  if(rcThis)
  {
    BmpDeleteDetailGainFunc(&adm);
  }
  return adm;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpTestDetailGainFunc(int n_iter, BmpBitmap **pB)
{
  int rcThis=0,lut_n,i;
  double s,g,p,l,dx,x,xmax,h=1;
  MtStat es,xs;
  BmpDetailGainFunc *dgf=NULL;
  BmpBitmap *B=BmpNewBitmap(1024,768,BMP_RGB|BMP_CLEAR);

  UT_ASSERT(B);
#if 0
  s=0.123;
  g=1.8;
  p=0.77;
  l=0.002;
  lut_n = 50000;
#endif

#if 1
  s=0.1;
  g=1;
  p=0.5;
  l=0.002;
  lut_n = 100000;
#endif

  printf("lut_n = %d, (1/n=%g)\n",lut_n,1/(double)lut_n); 

  dgf = BmpCreateDetailGainFunc( 1, 
      				 s, g, p, l, lut_n ); 
  if(!dgf) raiseRc(MI_ERROR);

  xmax=5;

  dx = (2*xmax)/(double)n_iter;

  MT_STAT_INIT(&es);
  MT_STAT_INIT(&xs);
  for(i=-n_iter/2;i<=n_iter/2;i++)
  {
    x = i*dx;	
    double e,f0=BmpComputeDetailGainFunc(dgf,0,x),
	   f=BmpLaplacianGain_FUNC1(x,s,g,p,l,h);

    UT_ASSERT(isfinite(f0));
    e = f>MT_NUM_EPS ? 100*fabs((f0/f)-1) : 0;
   if(e>10) 
      TRC1(("%g -> %g %g  -> e=%g%%\n",x,f0,f,e));
    MT_STAT_UPDATE(&xs,x);
    MT_STAT_UPDATE(&es,e);

    {
      double x_min=-xmax,x_span=2*xmax,
	     f_min=-20,f_span=2*abs(f_min);
      int ix=B->sx*(x-x_min)/x_span,
	  iy=B->sy*(1-(f-f_min)/f_span),
	  iy0=B->sy*(1-(f0-f_min)/f_span),
	  iy1=B->sy*(1-(x-f_min)/f_span);
      
      if(BMP_IN_RANGEB(ix,iy,B))
        BMP_PIXEL(B,data,ix,iy) = BMP_MKRGB(255,50,255);
      if(BMP_IN_RANGEB(ix,iy0,B))
        BMP_PIXEL(B,data,ix,iy0) = BMP_MKRGB(50,255,50);
      if(BMP_IN_RANGEB(ix,iy1,B))
        BMP_PIXEL(B,data,ix,iy1) = BMP_MKRGB(50,255,50);
    }
  }
  MT_STAT_RESULT(&es);
  MT_STAT_RESULT(&xs);

  printf("x-stats:   ");MtStatPrint(&xs);
  printf("err-stats: ");MtStatPrint(&es);

rcCatch:
  BmpDeleteDetailGainFunc(&dgf);
  if(rcThis)
  {
    BmpDeleteBitmap(&B);
  }
  *pB = B;
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpUpSamplePyraLevelFast2( 
    	BmpBitmap *bmpLo, float *dataLo, 
    	BmpBitmap *bmpHi, float *dataHi, 
	float *dataLapl, 
	uint8_t *dataAlpha, int alphaMax,
	float lap_f, float lap_g, float lap_gs,	// laplacian float & gamma
	int oxh, int oyh,
        float *kernel5x5,  // optional 5x5 downsampling kernel
	int btyp, // boundary handling 
	double *pLaplSig	 // OUT: (optional) laplacian standard deviation
	)
{
#define FAST_UPSAMPLE_LINE2(gtype, dataLo, dataHi, dataLapl, dataAlpha )\
  {\
    gtype *pDataLo,*pDataHi,*p,*p0,*p1,*p2;\
    float *pDataLapl=NULL;\
    uint8_t *pDataAlpha=NULL;\
    int	i=0;\
\
    xl = xh/2;\
    pDataLo = dataLo + BMP_XYOFF(bmpLo,xl-1,yl-1);\
    pDataHi = dataHi + BMP_XYOFF(bmpHi,xh0,yh0);\
    if(dataLapl) pDataLapl = dataLapl + BMP_XYOFF(bmpHi,xh0,yh0);\
    if(dataAlpha) pDataAlpha = dataAlpha + BMP_XYOFF(bmpHi,xh0,yh0);\
\
    while(xh<=xhmax-r+oxh)\
    {\
      p = pDataLo;\
      p0 = p; p+=sxl;\
      p1 = p; p+=sxl;\
      p2 = p;	 \
\
      if(xh&1)\
      {\
	pDataLo++;\
	if(yh&1)\
	{\
	  /* odd/odd*/\
	  g = w4*(p1[1]+p1[2]+p2[1]+p2[2]);\
	}\
	else\
	{\
	  /*odd/even*/\
	  g = w1*(p0[1]+p0[2]+p2[1]+p2[2])+\
	      w7*(p1[1]+p1[2]);\
	}\
      }\
      else\
      {\
	if(yh&1)\
	{\
	  /*  even/odd*/\
	  g = w3*(p1[0]+p1[2]+p2[0]+p2[2])+\
	      w5*(p1[1]+p2[1]);\
	}\
	else\
	{\
	  /*  even/even*/\
	  g = w0*(p0[0]+p0[2]+p2[0]+p2[2])+\
	      w2*(p0[1]+p2[1]+p1[0]+p1[2])+\
	      w8*p1[1];\
	}\
      }\
      if(pDataAlpha)\
      {\
	if(do_lap_g)\
          pDataHi[i] = g + \
	    (((float)(pDataAlpha[i]))/alphaMax) * \
		/*BmpLaplacianGain_FUNC1(pDataLapl[i],lap_gs,lap_f,lap_g,0.002,1);*/\
		BmpComputeDetailGainFunc(dgf,0,pDataLapl[i]);\
	else\
          pDataHi[i] = g + \
		lap_f*(((float)(pDataAlpha[i]))/alphaMax) * (pDataLapl[i]);\
      }\
      else if(pDataLapl) \
      {\
	lap = pDataHi[i] - g;\
        pDataLapl[i] = lap;\
	lap_sum += lap; \
	lap_sum2 += lap*lap;\
      }\
      else \
      {\
        pDataHi[i] = g;\
      }\
      xh++; xh0++;\
      i++;\
    }\
  }
  
  int rcThis=0,dx,dy,r=2,xh0,yh0,do_pgi,dyy,yy1,yy2,yy0,sy,
      sxl=bmpLo->sx,do_lap_g=fabs(lap_g-1)>MT_NUM_EPS,
      xlmax=bmpLo->sx-1,ylmax=bmpLo->sy-1,
      xhmax=bmpHi->sx-1,yhmax=bmpHi->sy-1;
  float *krn=kernel5x5 ? kernel5x5 : BmpGetDfltPyraConvKernel5x5(),*pkrn;
  double wsum,w0,w1,w2,w3,w4,w5,w6,w7,w8,lap_sum2,lap_sum,lap_mean,lap_var;
  UtTimer tm;
  UtProgrInd pgi;
  BmpDetailGainFunc *dgf=NULL;
  int	do_time=0;

  UT_ASSERT(r==2);

//  UT_ASSERT(fabs(lap_g-1)<MT_NUM_EPS);

  UT_ASSERT(dataHi);

  if(dataAlpha)
  {
    UT_ASSERT(dataHi);
  }

  if(btyp==BMP_BDEFAULT) btyp=BMP_BCLIP;

  TRC(("UpSamplePyraLevelFast2: (%dx%d) -> (%dx%d)  lap_f=%g/%g(do_g=%d)\n",
	bmpLo->sx,bmpLo->sy,bmpHi->sx,bmpHi->sy,lap_f,lap_g,do_lap_g));

  UT_ASSERT(((xhmax+1) <= 2*(xlmax+1))&&((yhmax+1) <= 2*(ylmax+1)));

  if(do_lap_g)
  {
    dgf = BmpCreateDetailGainFunc( 1, 
			lap_gs, lap_f, lap_g, 0.002, 
			(double)bmpHi->sx*bmpHi->sy>500000?100000:0);
    if(!dgf) raiseRc(MI_ERROR);
  }

  wsum=0; pkrn=krn;
  for(dy=-r;dy<=r;dy++) for(dx=-r;dx<=r;dx++) wsum += *pkrn++;
  UT_ASSERT(pkrn-krn==(2*r+1)*(2*r+1));

  // normalized quadrant-symmetrical kernel weights
  w0=krn[0]/wsum; w1=krn[1]/wsum; w2=krn[2]/wsum;
  w3=krn[5]/wsum; w4=krn[6]/wsum; w5=krn[7]/wsum;
  w6=krn[10]/wsum; w7=krn[11]/wsum; w8=krn[12]/wsum;

  TRC(("w0=%g w1=%g w2=%g \n" 
	"w3=%g w4=%g w5=%g \n" 
	"w6=%g w7=%g w8=%g \n",w0,w1,w2,w3,w4,w5,w6,w7,w8));

  UT_ASSERT(fabs(w2-w6)<MT_NUM_EPS && 
            fabs(w7-w5)<MT_NUM_EPS &&
            fabs(w3-w1)<MT_NUM_EPS);

  // normalize subset of weights 
  wsum = 4*w4; 
  w4 /= wsum;	  
  wsum = 4*w1+2*w7;
  w1 /= wsum; w7 /=wsum;
  wsum = 4*w3+2*w5;
  w3 /= wsum; w5 /= wsum;
  wsum = 4*w0+4*w2+w8;
  w0 /= wsum; w2 /= wsum; w8 /= wsum;

  if((do_pgi = ((double)bmpHi->sx*bmpHi->sy>50000000)))
  {
      UT_PROGRIND_INIT2( &pgi, 0.02, 2 );
  }

  if((do_time = (LogLevel>=2)&&
		((LongInt)bmpHi->sx*bmpHi->sy>10000000))==TRUE)
  {
    TIMER_START(&tm);
  }

  //
  //
  //
  sy = bmpHi->sy;
  dyy = MAX(sy/20,1);
  lap_sum = lap_sum2 = 0;

  for(yy0=0; yy0<sy; yy0+=dyy)
  {
    if(do_pgi)
    {
      UT_PROGRIND( &pgi, yy0, 0, sy-1, "upsampling", NULL );
      if(UT_PROGRIND_LAPSED(&pgi))
      {
	UT_USER_BREAK(NULL,MI_ECANCEL);	
      }
    }
    yy1=yy0; yy2=MIN(yy0+dyy-1,sy-1);

#pragma omp parallel for private(xh0,yh0)\
    			 reduction(+:lap_sum,lap_sum2)
  for(yh0=yy1;yh0<=yy2;yh0++)
  {
    int dx,dy,xh,yh,xh2,yh2,xl,yl,xl1,yl1,xl2,yl2,linedone,fetch;
    LongInt po;
    double wsum,gsum,g=0,lap;
    float *pkrn,wgt;

    yh = yh0-oyh;
    yl = yh/2; 
    linedone = FALSE;
    for(xh0=0;xh0<bmpHi->sx;
	)
    {
      xh = xh0-oxh;
      xl1=(xh-r)/2; yl1=(yh-r)/2;
      xl2=(xh+r)/2; yl2=(yh+r)/2;

      if((!(BMP_IN_RANGE(xl1,yl1, 0,0,xlmax,ylmax)&&
	   BMP_IN_RANGE(xl2,yl2, 0,0,xlmax,ylmax))) || linedone)	
      {
	// 
	// handle border case 
	//
	pkrn=krn; wsum=gsum=0;
	for(dy=-r;dy<=r;dy++)
	  for(dx=-r;dx<=r;dx++)
	  {
	    wgt = *pkrn++;
	    xh2 = xh-dx; yh2 = yh-dy;
	    xl2=xh2/2; yl2=yh2/2;
	    if(2*xl2!=xh2 || 2*yl2!=yh2) continue;
	    fetch = TRUE;
            if(btyp==BMP_BCLIP)
	    {
	      if(!(BMP_IN_RANGEB(xl2,yl2,bmpLo))) continue;
	    }
            else if(btyp==BMP_BCONST)
	    {
	      BMP_CLIP_CONST(0,0,xlmax,ylmax, xl2,yl2);
	    }
            else if(btyp==BMP_BZERO)
	    {
	      if(!(BMP_IN_RANGEB(xl2,yl2,bmpLo))) 
	      {
	        g = 0;
	        fetch = FALSE;
	      }
	    }
	    else 
	    {
	      TRCERR(("unexpected boundary type=%d\n",btyp));
	      fetch = FALSE;
	    }

	    if(fetch)
	    {
	      po = BMP_GETPIXEL_OFFS(bmpLo,xl2,yl2);
	      g = dataLo[po];
	    }
	    gsum += wgt*g;
	    wsum += wgt;
	  }
	g = gsum/wsum;
	po = BMP_GETPIXEL_OFFS(bmpHi,xh0,yh0);
#if 0
	if(dataAlpha)
	{
	  dataHi[po] = g + (dataAlpha[po]/(float)255.0)*dataLapl[po];
	}
	else if(dataLapl)
	{
	  float gh = dataHi[po];
	  dataLapl[po] = gh-g;
	}
	else
	{
	  dataHi[po]=g;
	}
#else
      if(dataAlpha)
      {
	if(do_lap_g)
          dataHi[po] = g + 
	    (((float)(dataAlpha[po]))/alphaMax) * 
		BmpComputeDetailGainFunc(dgf,0,dataLapl[po]);
	else
          dataHi[po] = g + 
		lap_f*(((float)(dataAlpha[po]))/alphaMax) * dataLapl[po];
      }
      else if(dataLapl) 
      {
        lap = dataHi[po] - g;
	dataLapl[po] = lap;
	lap_sum += lap;
	lap_sum2 += lap*lap;
      }
      else 
      {
        dataHi[po] = g;
      }
#endif
	xh0++;
      }
      else
      {
	// 
	// non-border case process fast line
	//
	FAST_UPSAMPLE_LINE2(float, dataLo, dataHi, dataLapl, dataAlpha );
	linedone = TRUE;
      }
    }
  }
  }

  if(do_time)
  {
    TIMER_STOP(&tm);
    TRC(("CPU = %.1f seconds (%.0f pix/sec)\n", 
	  (double)(TIMER_DIFF_MS(&tm)/1000.0), 
	  (bmpHi->sx*bmpHi->sy)/(double)(TIMER_DIFF_MS(&tm)/1000.0)));
  }

  MT_SUM2MEANVAR(lap_sum,lap_sum2,((double)bmpHi->sx*bmpHi->sy),
      		 lap_mean,lap_var);
  if(pLaplSig)
  {
    *pLaplSig = sqrt(lap_var);
  }

rcCatch:
  BmpDeleteDetailGainFunc(&dgf);
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpUpSamplePyraLevelFast( 
    	BmpBitmap *bmpLo, float *dataLo, BmpUshort *dataLoUs,
    	BmpBitmap *bmpHi, float *dataHi, BmpUshort *dataHiUs,
	float *dataLapl, 
	int oxh, int oyh,
        float *kernel5x5,  // optional 5x5 upsampling kernel
	int btyp // boundary handling 
	)
{
#define FAST_UPSAMPLE_LINE(gtype, dataLo, dataHi, dataLapl )\
  {\
    gtype *pDataLo,*pDataHi,*p,*p0,*p1,*p2;\
    float *pDataLapl=NULL;\
\
    xl = xh/2;\
    pDataLo = dataLo + BMP_XYOFF(bmpLo,xl-1,yl-1);\
    pDataHi = dataHi + BMP_XYOFF(bmpHi,xh0,yh0);\
    if(dataLapl) pDataLapl = dataLapl + BMP_XYOFF(bmpHi,xh0,yh0);\
\
    while(xh<=xhmax-r+oxh)\
    {\
      p = pDataLo;\
      p0 = p; p+=sxl;\
      p1 = p; p+=sxl;\
      p2 = p;	 \
\
      if(xh&1)\
      {\
	pDataLo++;\
	if(yh&1)\
	{\
	  /* odd/odd*/\
	  g = w4*(p1[1]+p1[2]+p2[1]+p2[2]);\
	}\
	else\
	{\
	  /*odd/even*/\
	  g = w1*(p0[1]+p0[2]+p2[1]+p2[2])+\
	      w7*(p1[1]+p1[2]);\
	}\
      }\
      else\
      {\
	if(yh&1)\
	{\
	  /*  even/odd*/\
	  g = w3*(p1[0]+p1[2]+p2[0]+p2[2])+\
	      w5*(p1[1]+p2[1]);\
	}\
	else\
	{\
	  /*  even/even*/\
	  g = w0*(p0[0]+p0[2]+p2[0]+p2[2])+\
	      w2*(p0[1]+p2[1]+p1[0]+p1[2])+\
	      w8*p1[1];\
	}\
      }\
      if(pDataLapl) \
        *pDataLapl++ = *pDataHi++ - g;\
      else \
        *pDataHi++ = g;\
      xh++; xh0++;\
    }\
  }
  
  int rcThis=0,dx,dy,r=2,xh0,yh0,
      sxl=bmpLo->sx,
      xlmax=bmpLo->sx-1,ylmax=bmpLo->sy-1,
      xhmax=bmpHi->sx-1,yhmax=bmpHi->sy-1;
  float *krn=BmpGetDfltPyraConvKernel5x5(),*pkrn;
  double wsum,w0,w1,w2,w3,w4,w5,w6,w7,w8;
  UtTimer tm;
  int	do_time=0;

  UT_ASSERT(kernel5x5==NULL); // TODO
  UT_ASSERT(r==2);

  if(dataHiUs||dataLoUs)
  {
    UT_ASSERT(!(dataLo||dataHi));
  }
  if(dataLo||dataHi)
  {
    UT_ASSERT(!(dataLoUs||dataHiUs));
  }

  if(btyp==BMP_BDEFAULT) btyp=BMP_BCLIP;

  TRC(("UpSamplePyraLevelFast: (%dx%d) -> (%dx%d)\n",
	bmpLo->sx,bmpLo->sy,bmpHi->sx,bmpHi->sy));

  UT_ASSERT(((xhmax+1) <= 2*(xlmax+1))&&((yhmax+1) <= 2*(ylmax+1)));

  wsum=0; pkrn=krn;
  for(dy=-r;dy<=r;dy++) for(dx=-r;dx<=r;dx++) wsum += *pkrn++;
  UT_ASSERT(pkrn-krn==(2*r+1)*(2*r+1));

  // normalized quadrant-symmetrical kernel weights
  w0=krn[0]/wsum; w1=krn[1]/wsum; w2=krn[2]/wsum;
  w3=krn[5]/wsum; w4=krn[6]/wsum; w5=krn[7]/wsum;
  w6=krn[10]/wsum; w7=krn[11]/wsum; w8=krn[12]/wsum;

  TRC(("w0=%g w1=%g w2=%g \n" 
	"w3=%g w4=%g w5=%g \n" 
	"w6=%g w7=%g w8=%g \n",w0,w1,w2,w3,w4,w5,w6,w7,w8));

  UT_ASSERT(fabs(w2-w6)<MT_NUM_EPS && 
            fabs(w7-w5)<MT_NUM_EPS &&
            fabs(w3-w1)<MT_NUM_EPS);

  // normalize subset of weights 
  wsum = 4*w4; 
  w4 /= wsum;	  
  wsum = 4*w1+2*w7;
  w1 /= wsum; w7 /=wsum;
  wsum = 4*w3+2*w5;
  w3 /= wsum; w5 /= wsum;
  wsum = 4*w0+4*w2+w8;
  w0 /= wsum; w2 /= wsum; w8 /= wsum;

  if((do_time = (LogLevel>=2)&&
		((LongInt)bmpHi->sx*bmpHi->sy>10000000))==TRUE)
  {
    TIMER_START(&tm);
  }

  //
  //
  //
#pragma omp parallel for private(xh0,yh0)
  for(yh0=0;yh0<bmpHi->sy;yh0++)
  {
    int dx,dy,xh,yh,xh2,yh2,xl,yl,xl1,yl1,xl2,yl2,linedone,fetch;
    LongInt po;
    double wsum,gsum,g=0;
    float *pkrn,wgt;

    yh = yh0-oyh;
    yl = yh/2; 
    linedone = FALSE;
    for(xh0=0;xh0<bmpHi->sx;
	)
    {
      xh = xh0-oxh;
      xl1=(xh-r)/2; yl1=(yh-r)/2;
      xl2=(xh+r)/2; yl2=(yh+r)/2;

      if((!(BMP_IN_RANGE(xl1,yl1, 0,0,xlmax,ylmax)&&
	   BMP_IN_RANGE(xl2,yl2, 0,0,xlmax,ylmax))) || linedone)	
      {
	// 
	// handle border case 
	//
	pkrn=krn; wsum=gsum=0;
	for(dy=-r;dy<=r;dy++)
	  for(dx=-r;dx<=r;dx++)
	  {
	    wgt = *pkrn++;
	    xh2 = xh-dx; yh2 = yh-dy;
	    xl2=xh2/2; yl2=yh2/2;
	    if(2*xl2!=xh2 || 2*yl2!=yh2) continue;
	    fetch = TRUE;
            if(btyp==BMP_BCLIP)
	    {
	      if(!(BMP_IN_RANGEB(xl2,yl2,bmpLo))) continue;
	    }
            else if(btyp==BMP_BCONST)
	    {
	      BMP_CLIP_CONST(0,0,xlmax,ylmax, xl2,yl2);
	    }
            else if(btyp==BMP_BZERO)
	    {
	      if(!(BMP_IN_RANGEB(xl2,yl2,bmpLo))) 
	      {
	        g = 0;
	        fetch = FALSE;
	      }
	    }
	    else 
	    {
	      TRCERR(("unexpected boundary type=%d\n",btyp));
	      fetch = FALSE;
	    }

	    if(fetch)
	    {
	      po = BMP_GETPIXEL_OFFS(bmpLo,xl2,yl2);
	      g = dataLo ? dataLo[po] : dataLoUs[po];
	    }
	    gsum += wgt*g;
	    wsum += wgt;
	  }
	g = gsum/wsum;
	po = BMP_GETPIXEL_OFFS(bmpHi,xh0,yh0);
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
	xh0++;
      }
      else
      {
	// 
	// non-border case process fast line
	//
	if(dataHi)
	{
	  FAST_UPSAMPLE_LINE(float, dataLo, dataHi, dataLapl );
	}
	else
	{
	  FAST_UPSAMPLE_LINE(BmpUshort, dataLoUs, dataHiUs, dataLapl );
	}
	linedone = TRUE;
      }
    }
  }

  if(do_time)
  {
    TIMER_STOP(&tm);
    TRC(("CPU = %.1f seconds (%.0f pix/sec)\n", 
	  (double)(TIMER_DIFF_MS(&tm)/1000.0), 
	  (bmpHi->sx*bmpHi->sy)/(double)(TIMER_DIFF_MS(&tm)/1000.0)));
  }

rcCatch:
    return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpUpSamplePyraLevel( BmpBitmap *bmp0,
    		           BmpBitmap *bmpLo, 
			   float *dataLo, double *dataLo_d,
    		           BmpBitmap *bmpHi, 
			   float *dataHi, double *dataHi_d,
		           double res_factor_f, 
			   float *dataLapl, 
          		   float *kernel5x5  // optional 5x5 downsampling kernel

			   )
{

  int rcThis=0,x,y,
  	xmaxLo=bmpLo->sx-1,ymaxLo=bmpLo->sy-1,
	r,krn_r=2;
  float *krn = (kernel5x5!=NULL ? kernel5x5 : bmp0->pyramid.krn);

  if(!((dataLo||dataLo_d) && (dataHi||dataHi_d)))
  {
    TRCERRR(("no bitmap data\n"),MI_ERROR);
  }

  if(res_factor_f<MT_REAL_EPSILON) res_factor_f = 2;
  r = krn_r;
#pragma omp parallel for private(y,x)    
  for(y=0;y<bmpHi->sy;y++)
  {
    int x2,y2,dx,dy,po,x2h,y2h;
    float wgt,*pkrn;
    double wsum,gsum,g;
    for(x=0;x<bmpHi->sx;x++)
    {
      pkrn=krn; wsum=gsum=0;
      for(dy=-r;dy<=r;dy++)
	for(dx=-r;dx<=r;dx++)
	{
          wgt = *pkrn++;
	  x2h = x-dx; y2h = y-dy;
          x2=x2h/2; y2=y2h/2;
	  if(2*x2!=x2h || 2*y2!=y2h) continue;
/*	  printf("%+d %+d  >> %d, %d  -  %d %d\n",
	      dx,dy,2*x2,x2h,2*y2,y2h);*/
//          if(!BMP_IN_RANGE(x2,y2,0,0,xmaxLo,ymaxLo)) continue;
/*	  printf("%+d %+d  >> %d, %d  -  %d %d\n",
	      dx,dy,2*x2,x2h,2*y2,y2h);*/
#if GP_HDL_BORD_CLIP
	  if(!(BMP_IN_RANGE(x2,y2,0,0, xmaxLo,ymaxLo))) continue;
#else
	  BMP_CLIP_CONST(0,0,xmaxLo,ymaxLo, x2,y2);
#endif
          po = BMP_GETPIXEL_OFFS(bmpLo,x2,y2);
	  g = dataLo ? dataLo[po] : dataLo_d[po];
	  gsum += wgt*g;
	  wsum += wgt;
	}
      //UT_ASSERT(wsum > MT_REAL_EPSILON);
      g = gsum/wsum;
      po = BMP_GETPIXEL_OFFS(bmpHi,x,y);
      
      if(dataLapl)
      {
        g = dataHi ? BMP_LAPL_SCALE*dataHi[po]-BMP_LAPL_SCALE*g : 
	  	     BMP_LAPL_SCALE*dataHi_d[po]-BMP_LAPL_SCALE*g;
	dataLapl[po] = g;
      }
      else
      {
	UT_ASSERT0(isfinite(g)&&!isnan(g));
	if(dataHi)
	  dataHi[po] = g;
	else
	  dataHi_d[po] = g;
      }
    }
  }
  BMP_MODIFIED(bmpHi);

rcCatch:
    return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BmpBitmap *BmpIExpandPyramid( BmpBitmap *bmp0,
    			      BmpBitmap *bmpSrc,
			      BmpBitmap *bmpDst,
    		              int levelDst,
			      	// target level
			      float res_factor,
		              unsigned opt )
{
  int rcThis=0,rc,i,nChanDone=0;
  int sx,sy,level=levelDst;
  float *dataHi, *dataLo;
  BmpBitmap *bmp=NULL,*bmpThis;

  if(bmp0)
  {
    UT_ASSERT(levelDst>=0&&levelDst<bmp0->pyramid.n_level-1);
    bmpSrc = bmpSrc ? bmpSrc : bmp0->pyramid.gaussian[level+1],
    bmpThis = bmp0->pyramid.gaussian[level];
     sx = bmpThis->sx; 
     sy = bmpThis->sy;
  }
  else
  {
    sx = bmpSrc->sx * res_factor;
    sy = bmpSrc->sy * res_factor;
  }

  bmp = BmpNewBitmap(sx,sy,opt);
  if(!bmp) raiseRc(MI_ERROR);

  TRC(("Bmp pyramid EXPAND %d,%dx%d -> %d,%dx%d\n",
	level+1,bmpSrc->sx,bmpSrc->sy,level,sx,sy));

  if(opt & BMP_GREY)
  {
    dataHi = bmp->dataGrey;
    dataLo = BmpGetGreyChannel(bmpSrc,0);
    rc = BmpUpSamplePyraLevelFast( bmpSrc, dataLo, NULL,
				   bmp, dataHi, NULL,
				   NULL, 0, 0, NULL, BMP_BCONST);
    if(rc) raiseRc(MI_ERROR);
    nChanDone++;
  }
  if(opt & BMP_Z)
  {
    UT_ASSERT(FALSE);
  }
  if(opt & BMP_XYZ)
  {
    for(i=0;i<3;i++)
    {
      dataHi = bmp->dataFloat[i];
      dataLo = BmpGetFloatChannel(bmpSrc,i,0);
      rc = BmpUpSamplePyraLevel( bmp0, bmpSrc, dataLo, NULL,
	  				bmp, dataHi, NULL,
	  		      	        res_factor,0,NULL);
      if(rc) raiseRc(MI_ERROR);
      nChanDone++;
    }
  }
  if(!nChanDone)
  {
    TRCERRR(("no bitmap channel selected"),MI_ERROR);
  }

rcCatch:
  return bmp;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BmpBitmap *BmpExpandPyramidOld( BmpBitmap *bmp0,
    			BmpBitmap *bmpFrom,
    		      int levelFrom,
		      int levelTo,
		      unsigned opt
		      )
{
  int i,rcThis=0;
  BmpBitmap *bmp=bmpFrom ? bmpFrom : bmp0->pyramid.gaussian[levelFrom],
	    *bmp2;
  for(i=levelFrom;i>levelTo;i--)
  {
    bmp2 = BmpIExpandPyramid(bmp0,bmp,NULL,i-1,2,opt);
    if(!bmp2) raiseRc(MI_ERROR);
    if(i!=levelFrom) BmpCloseBitmap(bmp);
    bmp = bmp2;
  }
rcCatch:
  return bmp;
}



/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpGetLaplacianPyramid( BmpBitmap *bmp0, unsigned opt )
{
  BmpPyramid *py=&bmp0->pyramid;
  int rcThis=0,i,maxlevel=py->n_level-1,po,k;
  BmpBitmap  *bmpG=NULL,	// Gaussian current
	     *bmpGLE=NULL,	// Expanded lower gaussian
	     *bmpH=NULL;	// Laplacian

  if(py->n_level<=0)
    TRCERRR(("Missing gaussian pyramid\n"),MI_ERROR);

/*  if((opt & BMP_GREY) || !(opt & BMP_XYZ))
    TRCERRR(("not implemented\n"),MI_ENOTIMPL);*/

  for(i=maxlevel;i>=0;i--)
  {
    bmpG = py->gaussian[i];
    BmpCloseBitmap(py->laplacian[i]); py->laplacian[i]=NULL;
    py->laplacian[i] = bmpH = BmpNewBitmap(bmpG->sx,bmpG->sy,opt);
    if(!bmpH) raiseRc(MI_ERROR);
    if(i<maxlevel)
    {
      BmpCloseBitmap(bmpGLE);bmpGLE=NULL;
      bmpGLE = BmpExpandPyramidOld(bmp0,NULL,i+1,i,opt);
      if(!bmpGLE) raiseRc(MI_ERROR);
    }
    else bmpGLE = NULL;

    // Lowest level -> L = G
#pragma omp parallel for private(po,k)    
    for(po=0;po<bmpH->sx*bmpH->sy;po++)
    {
      if(opt & BMP_GREY) 
      {
	// multiply with BMP_LAPLACE to avoid numerical instabilities
	// due to very small differences
	bmpH->dataGrey[po] = 
	  bmpGLE ? BMP_LAPL_SCALE * bmpG->dataGrey[po] - 
	  	     BMP_LAPL_SCALE * bmpGLE->dataGrey[po] :
	  	   BMP_LAPL_SCALE * bmpG->dataGrey[po];
	  //UT_ASSERT(!isnan(bmpH->dataGrey[po]));
      }
      if(opt & BMP_Z) 
      {
	bmpH->dataZ[po] = 
	  bmpGLE ? BMP_LAPL_SCALE * bmpG->dataZ[po] - 
	  	     BMP_LAPL_SCALE * bmpGLE->dataZ[po] :
	  	   BMP_LAPL_SCALE * bmpG->dataZ[po];
	  //UT_ASSERT(!isnan(bmpH->dataZ[po]));
      }
      if( opt & BMP_XYZ )
      {
	for(k=0;k<3;k++) 
	{
	  bmpH->dataFloat[k][po] = 
	    bmpGLE ? 
	      BMP_LAPL_SCALE * bmpG->dataFloat[k][po] - 
	        BMP_LAPL_SCALE * bmpGLE->dataFloat[k][po] :
	      BMP_LAPL_SCALE * bmpG->dataFloat[k][po];
	  //UT_ASSERT(!isnan(bmpH->dataFloat[k][po]));
	}
      }
    }      	
  }

rcCatch:
  BmpCloseBitmap(bmpGLE);bmpGLE=NULL;
  if(!rcThis) bmp0->pyramid.laplace_valid = TRUE;
  return rcThis;
}

#if 0
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpGetLaplacianPyramid2( BmpBitmap *bmp0, unsigned opt )
{
  BmpPyramid *py=&bmp0->pyramid;
  int rc,rcThis=0,i,maxlevel=py->n_level-1,po;
  BmpBitmap  *bmpG=NULL,	// Gaussian current
	     *bmpL=NULL,
	     *bmpH=NULL;	// Laplacian

  if(py->n_level<=0)
    TRCERRR(("Missing gaussian pyramid\n"),MI_ERROR);

  UT_ASSERT(opt==BMP_GREY);

  for(i=maxlevel;i>=0;i--)
  {
    bmpG = py->gaussian[i];
    BmpCloseBitmap(py->laplacian[i]); py->laplacian[i]=NULL;
    py->laplacian[i] = bmpH = BmpNewBitmap(bmpG->sx,bmpG->sy,opt);
    if(!bmpH) raiseRc(MI_ERROR);
    if(i<maxlevel)
    {
      bmpL = py->gaussian[i+1];
      rc = BmpUpSamplePyraLevel( bmp0, 
			bmpL, bmpL->dataGrey, NULL, 
			bmpH, bmpH->dataGrey, NULL,
			0,
			bmpH->dataGrey );
      if(rc) raiseRc(MI_ERROR);
    }
    else 
    {
      // laplacian at coarsest level = gaussian
      for(po=0;po<bmpH->sx*bmpH->sy;po++)
      {
        bmpH->dataGrey[po] = BMP_LAPL_SCALE * bmpG->dataGrey[po];
      }      	    
    }
  }

rcCatch:
  if(!rcThis) bmp0->pyramid.laplace_valid = TRUE;
  return rcThis;
}
#endif
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
BmpBitmap *BmpReconPyramid( BmpBitmap *bmp0, unsigned opt )
{
  int i,ih,rcThis=0,levelFrom=bmp0->pyramid.n_level-1,po,k;
  BmpBitmap *bmp= bmp0->pyramid.gaussian[levelFrom],
	    *bmp2,*bmpL;
  if(!bmp0->pyramid.laplace_valid) 
    TRCERRR(("no laplacian pyramid\n"),MI_ERROR);
  for(i=levelFrom;i>0;i--)
  {
    ih = i-1;
    // Expand to next level
    bmp2 = BmpIExpandPyramid(bmp0,bmp,NULL,ih,2,opt);
    if(!bmp2) raiseRc(MI_ERROR);
    if(i!=levelFrom) BmpCloseBitmap(bmp);
    bmp = bmp2;
    // Add laplacian
    bmpL = bmp0->pyramid.laplacian[ih];
    if(!bmpL) raiseRc(MI_ERROR);
    for(po=0;po<bmp->sx*bmp->sy;po++)
    {
      if(opt & BMP_GREY) 
	bmp->dataGrey[po] += bmpL->dataGrey[po] / BMP_LAPL_SCALE;
      if(opt & BMP_Z) 
	bmp->dataZ[po] += bmpL->dataZ[po] / BMP_LAPL_SCALE;
      if( opt & BMP_XYZ )
	for(k=0;k<3;k++) 
	  bmp->dataFloat[k][po] 
	    += 1*bmpL->dataFloat[k][po] / BMP_LAPL_SCALE; 
      }
    }

rcCatch:
  return bmp;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpGetBurtAdelsonKernel5x5(
    		float a, // IN: smoothness/shape parameter  
			 //   0.6=neg.lobes, 0.5=triangle, 
			 //   0.4=gauss, 0.3=broader gauss
			 //   0 = box
			 //   see also [Burt1983]
    		float *kernel5x5  // OUT: 5x5 kernel 
		)
{
  int rcThis=0,dx,dy,r=2;
  double w0[3]={0},w,wsum;

  TRC3(("BmpGetBurtAdelsonKernel5x5() a=%g\n",a));

  if(a>MT_NUM_EPS)
  {
    w0[0] = a;
    w0[1] = 1./4;
    w0[2] = 1./4 - a/2.;
  }

  wsum=0;
  for(dy=-r;dy<=r;dy++) 
  {
    for(dx=-r;dx<=r;dx++)
    {
      if(a>MT_NUM_EPS)
      {
        w = w0[ABS(dx)]*w0[ABS(dy)];
      }
      else
      {
	w = dx>=0&&dy>=0&&dx<=1&&dy<=1 ? 0.25 : 0;
      }
      kernel5x5[MT_GXY(5,dx+r,dy+r)] = w;
      wsum += w;
      TRC(("%7.4f ",w));
    }
    TRC(("\n"));
  }
  TRC(("wsum=%g\n",wsum));

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
float *BmpGetDfltPyraConvKernel5x5(void)
{
  return bmp_dflt_pyra_kern_5x5;

#if 0 // XXX on-the-fly initialization not thread safe !
  float *pkrn,*krn=bmp_dflt_pyra_kern_5x5;
  
  if(krn[0]<0)
  {
    int dx,dy,r=2;
    float w0_buf[100],*w0=w0_buf+r,a_kern=0.4;
    double wsum,w,wsum2;
    // Create convolution kernel (cf. Burt 1983)
    w0[0] = a_kern;
    w0[-1] = w0[1] = 1./4;
    w0[-2] = w0[2] = 1./4 - a_kern/2.;
    
    TRC(("BmpGetDfltPyraConvKernel5x5:\n"));

    pkrn = krn; wsum=0;
    for(dy=-r;dy<=r;dy++) 
      for(dx=-r;dx<=r;dx++)
      {
	w = w0[dx]*w0[dy];
	wsum += w;
	*pkrn++ = w;
      }
    pkrn = krn;
    wsum2=0;
    for(dy=-r;dy<=r;dy++) 
    {
      for(dx=-r;dx<=r;dx++)
      {
	w = *pkrn / wsum;
	*pkrn++ = w;
	wsum2+=w;
	TRC(("  %7.4f ",w));
      }
      TRC(("\n"));
    }
//    TRC(("sum=%g\n",wsum2));
  }
  return krn;
#endif
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpCreatePyramidOld( BmpBitmap *bmp0,
    		      float	res_factor,
		      	// Resolution reduction factor (default: 2)
		      float	a_kern,
		        // Generating kernel 'a' parameter (-> Burt 1983)
		      int maxLevel,
    		      unsigned opt )
{
  int rcThis=0,rc,r,i,dx,dy;
  BmpPyramid *pyra=&bmp0->pyramid;
  float	 *w0_buf=NULL,*w0,*pkrn,w;
  double   wsum;

  USE(i);

  // Parameter defaults
  if(res_factor < MT_REAL_EPSILON) res_factor = 2;
  if(a_kern<MT_REAL_EPSILON) a_kern = 0.4;

  // Parameter plausis
  if((int)res_factor!=2) 
    TRCERRR(("not implemented\n"),MI_ENOTIMPL);

  // Free old pyramid 
  BmpIFreePyramid(bmp0);

  pyra->krn_r = r = (int)(res_factor+0.5);
  UT_ASSERT(r==2);
  ALLOCARR(w0_buf,2*r+1);
  w0 = w0_buf+r;
  pyra->chan = BMP_GEN;
  pyra->ichan = 0;
  pyra->n_level=0;
  pyra->res_0 = 1;
  pyra->res_act = 1;
  pyra->res_factor = res_factor;
  pyra->opt = opt;
  pyra->bmp0 = bmp0;
  TRC(("Creating image pyramid: (r=%g)\n",res_factor));

  // Create convolution kernel (cf. Burt 1983)
  w0[0] = a_kern;
  w0[-1] = w0[1] = 1./4;
  w0[-2] = w0[2] = 1./4 - a_kern/2.;
/*  TRC1(("Generating Kernel a=%g w0= \n"));
  for(i=-r;i<=r;i++) printf("%7.4f ",w0[i]); printf("\n");*/
  ALLOCARR(pyra->krn,(2*r+1)*(2*r+1));

  pkrn = pyra->krn; wsum=0;
  for(dy=-r;dy<=r;dy++) 
    for(dx=-r;dx<=r;dx++)
    {
      w = w0[dx]*w0[dy];
      wsum += w;
      *pkrn++ = w;
    }
  pkrn = pyra->krn;

  for(dy=-r;dy<=r;dy++) 
  {
    for(dx=-r;dx<=r;dx++)
    {
      w = *pkrn / wsum;
      *pkrn++ = w;
      TRC(("%7.4f ",w));
    }
    TRC(("\n"));
  }
  
  {
    int k;
    MT_NEAREST_POW2( bmp0->sx, bmp0->vsx, k, 1<<28 );  
    MT_NEAREST_POW2( bmp0->sy, bmp0->vsy, k, 1<<28 );  
  }

//  pyra->krn = MtMakeGaussKernel(pyra->krn_r,0.5);
  TRC(("Bmp pyramid: %2d %7.1f   %dx%d (%dx%d)\n",
	pyra->n_level,pyra->res_act,bmp0->sx,bmp0->sy,
	bmp0->vsx,bmp0->vsy));

  pyra->gaussian[pyra->n_level] = bmp0;

  // Recursively create pyramid levels
  rc = BmpIReducePyramid( bmp0, pyra->n_level+1, opt );
  if(rc) raiseRc(rc);

  // Optionally create the laplacian pyramid
  if(opt & BMP_LAPLACE)
  {
    rc = BmpGetLaplacianPyramid(bmp0,opt);
    if(rc) raiseRc(MI_ERROR);
  }
  
rcCatch:
  FREEMEM(w0_buf);
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpGetPyramid( BmpBitmap *bmp, unsigned opt )
{
  int rc,rcThis=0;
  UT_ASSERT( opt == (BMP_GREY | BMP_GAUSS) ); // others not yet supported
  if (! BmpHasPyramid( bmp, opt ))
  {
    rc = BmpCreatePyramidOld( bmp, -1, -1, -1, opt );
    UT_ASSERT(rc==0);
  }
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpFindPyramidLevel( BmpBitmap *bmp0, int md0, int sxy0,
    			 int *pLevel )
{
  int rcThis=0,level=-1,i,sxy,nlev=bmp0->pyramid.n_level;
  BmpBitmap *bmp;
  for(i=nlev-1;i>=0;i--)
  {
    bmp = bmp0->pyramid.gaussian[i]; 
    sxy = MIN(bmp->sx,bmp->sy);
    if(sxy0>0 && md0<0)
    {
      if(sxy >= sxy0) 
      {
	level = i;
	break;
      }
    }
    else if(md0>0 && sxy0<0)
    {
      if( bmp->pyramid.res_act<= md0  )
      {
	level = i;
	break;
      }
    }
    else TRCERRR(("Assertion failed\n"),MI_ERROR);    
  }
  
rcCatch:
  if(level<0)
  {
    TRCERR(("Pyramid level (max=%d) not found for md0=%d, sxy0=%d\n",
	  nlev,sxy0));
    rcThis = MI_ENOTFOUND;
  }
  *pLevel = level;
  return rcThis;
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpFreePyramid( BmpBitmap *bmp, unsigned opt )
{
  int rcThis=0;
  UT_ASSERT(opt==0); // TODO
  rcThis= BmpIFreePyramid( bmp );
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static int BmpIFreePyramid( BmpBitmap *bmp )
{
  int i,rc;
  if(bmp)
  {
    BmpPyramid *pyra = &bmp->pyramid;
    FREEMEM(pyra->krn); pyra->krn = NULL;
    for(i=0;i<pyra->n_level;i++)
    {
      BmpBitmap *bmp2=pyra->gaussian[i];
      if(bmp2!=bmp)
      {
	rc= BmpCloseBitmap(bmp2);
	if(rc) TRCERR(("BmpCloseBitmap() failed %d\n",rc));
	pyra->gaussian[i] = NULL;
      }
      BmpCloseBitmap(pyra->laplacian[i]); pyra->laplacian[i] = NULL;
    }
    pyra->n_level = 0;
  }
  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpHasPyramid( BmpBitmap *bmp, int opt )
{
  BmpPyramid *py=&bmp->pyramid;
  if(py->n_level<=0) return FALSE;
  if(opt & BMP_LAPLACE) return py->laplace_valid;
  if(opt & BMP_GAUSS) return TRUE;
  TRCERR(("unexpected options %d\n",opt));
  return FALSE;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpGetPixStats( BmpBitmap *B, 
    		    float *data,
    		    int x1, int y1, 
    		    int x2, int y2,
		    MtStat *stat )
{
  int rcThis=0,x,y;
  float zmax,zmin;
  double zsum,zsum2,nPix;

  if(x1==0&&y1==0&&x2==0&&y2==0)
  {
    x2=B->sx-1; y2=B->sy-1;
  }
  UT_ASSERT(BMP_IN_RANGEB(x1,y1,B)&&BMP_IN_RANGEB(x2,y2,B));
  
  TRC(("gathering pixel stats (%dx%d) %d,%d - %d,%d\n",
	B->sx,B->sy,x1,y1,x2,y2));

  zmax=-1e20; zmin=1e20; zsum=0; zsum2=0;
#pragma omp parallel
  {
    float my_zmax=-1e20, my_zmin=1e20;
    double my_zsum=0, my_zsum2=0;
#pragma omp for private(y,x) 
    for(y=y1;y<=y2;y++)
    {
      double z;
      for(x=x1;x<=x2;x++)
      {
	z = data[BMP_XYOFF(B,x,y)];
	my_zsum += z;
	my_zsum2 += z*z;
	if(z>my_zmax) my_zmax=z;
	if(z<my_zmin) my_zmin=z;
      }
    }
#pragma omp critical
    {
      zmax=MAX(zmax,my_zmax);
      zmin=MIN(zmin,my_zmin);
      zsum+=my_zsum;
      zsum2+=my_zsum2;
    }
  }
  nPix=(x2-x1+1)*(double)(y2-y1+1);
  MT_STAT_INIT(stat);
  MT_STAT_UPDATE_N(stat,nPix,zmin,zmax,zsum,zsum2);
  MT_STAT_RESULT(stat);

  if(!isfinite(zsum))
  {
    TRCERRR(("invalid floating point number %g\n",zsum),
	MI_ENUMERICAL);
  }
 
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpNormalizeChannel( BmpBitmap *I, int chan, int ichan )
{
  int rc,rcThis=0;
  LongInt ii;
  BmpData data;
  MtStat stat;

  BmpGetDataChannel(I,chan,ichan,0,&data);
  UT_ASSERT(data.pFloat); 

  rc = BmpGetPixStats(I,data.pFloat,0,0,0,0,&stat);
  if(rc) raiseRc(rc);

  if(stat.span<MT_NUM_EPS)
  {
    //TRCWARN(("zero variance image\n"));
    stat.span=1;
  }

  #pragma omp parallel for private(ii)
  for(ii=0;ii<I->sx*I->sy;ii++)
  {
    double g = (data.pFloat[ii]-stat.min)/stat.span;
    if(!isfinite(g))
    {
      UT_ASSERT0(FALSE);
    }
    if(g<0||g>1)
    {
      UT_ASSERT0(FALSE);
      if(g<0)g=0;
      if(g>1)g=1;
    }
    data.pFloat[ii] = g;
  }
rcCatch:
  return rcThis;
}

MSBG_NAMESPACE_END

