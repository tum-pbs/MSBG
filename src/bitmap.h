/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#ifndef BITMAP_H
#define BITMAP_H

MSBG_NAMESPACE_BEGIN

#ifdef __cplusplus
extern "C" {
#endif

#define RGB_TO_GREY( r,g,b) (coGlobConst.grayweight_red*(r)+ \
    			     coGlobConst.grayweight_green*(g)+ \
			     coGlobConst.grayweight_blue*(b))

#define CO_RGB2INT( vRgb, icr, icg, icb ) \
{ \
  icr = vRgb[0]*255.; icg = vRgb[1]*255.; icb = vRgb[2]*255.; \
  CoClipRGB( &icr, &icg, &icb, 255 ); \
}

#define CO_CLIP_RGB_F(vRgb,oog) \
{ \
  oog=0; \
  if((vRgb)[0]>1) {(vRgb)[0]=1;oog=1;} if((vRgb)[0]<0) {(vRgb)[0]=0;oog=1;}\
  if((vRgb)[1]>1) {(vRgb)[1]=1;oog=1;} if((vRgb)[1]<0) {(vRgb)[1]=0;oog=1;}\
  if((vRgb)[2]>1) {(vRgb)[2]=1;oog=1;} if((vRgb)[2]<0) {(vRgb)[2]=0;oog=1;}\
}

#define CO_CHECK_RGB_F(vRgb,oog) \
{ \
  oog = ((vRgb)[0]>1)||((vRgb)[0]<0)||\
	((vRgb)[1]>1)||((vRgb)[1]<0)||\
	((vRgb)[2]>1)||((vRgb)[2]<0);\
}

#define CO_OUT_OF_GAMUT		9

typedef struct
{
  /* !!! Attention: keep consistent with default values in color.c */
  float	rgbClippingStrength,
  	grayweight_red,
	grayweight_green,
	grayweight_blue,
	gamma;
  double transMatRGB2LXY[3*3];
}
CoGlobalConst;

extern CoGlobalConst coGlobConst;

void CoClipRGB( int *pr, int *pg, int  *pb, int max );

#define BMP_VALID	(1<<0)
#define BMP_RGB 	(1<<1)
#define BMP_ALPHA 	(1<<2)
#define BMP_GREY	(1<<3)
#define BMP_IGNAN	(1<<4)
#define BMP_CLEAR       (1<<5)
#define BMP_VERBOSE	(1<<6)
#define BMP_DOUBLE	(1<<7)
#define BMP_USHORT	(1<<8)
#define BMP_32BIT	(1<<9)
#define BMP_UV		(1<<10)
#define BMP_DEST_LUV	(1<<11)
#define BMP_DEST_NORGB  (1<<12)
#define BMP_GAUSS	(1<<13)
#define BMP_EXTERN	(1<<14)
#define BMP_UBYTE	(1<<15)
#define BMP_Z		(1<<16)
#define BMP_N		(1<<17)
#define BMP_COMPR	(1<<18)
#define BMP_INPLACE	(1<<19)
#define BMP_16BIT	(1<<20)
#define BMP_GREY_ONLY	(1<<21)
#define BMP_NONORM	(1<<22)
#define BMP_XYZ		(1<<23)
#define BMP_LUV		(1<<24)
#define BMP_LAPLACE	(1<<26)
#define BMP_FREE	(1<<27)
#define BMP_GEN		(1<<28)
#define BMP_CHECK	(1<<29)
#define BMP_FLOAT	(1<<30)
#define BMP_NORM	(1<<31)

#define BMP_IPOL_NEAREST 0
#define BMP_IPOL_LINEAR	 1
#define BMP_IPOL_CUBIC	 2

#define BMP_IPOL_CORNER	(1<<1)

#define BMP_LAPL_SCALE	1000
//#define BMP_LAPL_SCALE  1

#define PNH_FULL		(1<<1)
#define PNH_INCLCENT	(1<<2)
#define PNH_MAXCHAN		10

#define BMP_MAX_PYRA_LEVEL	64

// boundary types
// do not change values !
#define BMP_BDEFAULT	0
#define BMP_BCONST	1	// constant extrapolation
#define BMP_BREFL	2	// reflect at edge
#define BMP_BWRAP	3	// wrap-around edge
#define BMP_BZERO	4	// set zero
#define BMP_BCLIP	5	// omit out-of-bund pixels
#define BMP_BANTIREFL	6	// anti-reflective (C1-contin.)

// kernel types
#define  BMP_KERN_BOX 		1
#define  BMP_KERN_GAUSS 	2
#define  BMP_KERN_TRIANGLE 	3
#define  BMP_KERN_HAT 		4

// file types
#define BMP_NULL	0
#define BMP_WINBMP	1
#define BMP_PNG		2
#define BMP_U48		3
#define BMP_JPG		4
#define BMP_MBM		5
#define BMP_TIF		6

#define BMP_GRAD_SOEBEL 1
#define BMP_GRAD_CDIFF	2
#define BMP_GRAD_FDIFF  3
#define BMP_GRAD_BDIFF  4

#define BMP_SMP_JITTER	(1<<9)

#define BMP_FFLAT	(1<<19)
#define BMP_FKERN	(1<<22)
#define BMP_FCONE	(1<<14)
#define BMP_FGAUSS	(1<<13)
#define BMP_FHAT	(1<<15)
#define BMP_FLAPLACE    (1<<16)
#define BMP_FCENTER	(1<<27)

#define BMP_GET_CHANNEL( B, chan, ichan, /* IN */\
    			   data, dataUs /* OUT */)\
{\
  data = NULL;\
  dataUs = NULL;\
  switch(chan)\
  {\
    case BMP_Z: data =  BmpGetZChannel(B,0); break;\
    case BMP_GREY: data =  BmpGetGreyChannel(B,0); break;\
    case BMP_FLOAT: data =  BmpGetFloatChannel(B,ichan,0); break;\
    case BMP_USHORT: dataUs =  BmpGetUshortChannel(B,ichan,0); break;\
    default:\
      UT_ASSERT(FALSE); \
      break;\
  }\
}

#define BMP_FREE_CHANNEL( B, chan, ichan ) \
{\
  switch(chan)\
  {\
    case BMP_Z: BmpFreeZChannel(B); break;\
    case BMP_GREY: BmpFreeGreyChannel(B); break;\
    case BMP_FLOAT: BmpFreeFloatChannel(B,ichan); break;\
    case BMP_USHORT: BmpFreeUshortChannel(B,ichan); break;\
    default:\
      UT_ASSERT(FALSE); \
      break;\
  }\
}

#define BMP_DATA_COPY(dst,i, src,j) \
{ \
  if((dst).pFloat) \
    ((dst).pFloat)[i] = (src).pFloat[j]; \
  else \
    ((dst).pUshort)[i] = ((src).pUshort)[j]; \
}

#define BMP_DATA_SET(dst,i, a) \
{ \
  if((dst).pFloat) \
    ((dst).pFloat)[i] = a; \
  else \
    ((dst).pUshort)[i] = a; \
}

#define BMP_DATA_GET(dst,i, af, au) \
{ \
  if((dst).pFloat) \
    af = ((dst).pFloat)[i]; \
  else \
    au  = ((dst).pUshort)[i]; \
}

#define BMP_GUID(bmp) ((bmp)->guid)

#define BMP_MODIFIED( bmp ) \
{ \
  (bmp)->doRegenerate = TRUE; \
}

#define BMP_LOOP_SXY(sx,sy,x,y) \
  for(y=0;y<(sy);y++) \
    for(x=0;x<(sx);x++) \

#define BMP_LOOP_XY(bmp,x,y)\
  BMP_LOOP_SXY((bmp)->sx,(bmp)->sy,x,y)


// Get bilinear patch with neuman boundary condition
// ATTENTION: upper left point (x,y) assumed already in range ! 
#define BMP_GET_BILP( sx, sy, dataZ, \
    		      x, y, \
    		      z1, z2, z3, z4 ) /* scanline order */\
{ \
  LongInt ii;\
\
  if((x)<(sx)-1&&(y)<(sy)-1)\
  { /* fast case */\
    ii = MT_GXY(sx,x,y);  \
    z1 = (dataZ)[ii];\
    ii++; \
    z2 = (dataZ)[ii];\
    ii += sx-1; \
    z3 = (dataZ)[ii];\
    ii++; \
    z4 = (dataZ)[ii];\
  }\
  else\
  {\
    LongInt i2,i3,i4;\
    z1 = (dataZ)[MT_GXY(sx,x,y)];\
    if((x)<(sx)-1)\
    {\
      i2 = MT_GXY(sx,x+1,y);\
      i3 = MT_GXY(sx,x,y-1);\
      i4 = MT_GXY(sx,x+1,y-1);\
      z2 = (dataZ)[i2];\
      z3 = 2*z1-(dataZ)[i3];\
      z4 = 2*z2-(dataZ)[i4];\
    }\
    else if((y)<(sy)-1)\
    {\
      i2 = MT_GXY(sx,x-1,y);\
      i3 = MT_GXY(sx,x,y+1);\
      i4 = MT_GXY(sx,x-1,y+1);\
      z2 = 2*z1-(dataZ)[i2];\
      z3 = (dataZ)[i3];\
      z4 = 2*z3-(dataZ)[i4];\
    }\
    else\
    {\
      i2 = MT_GXY(sx,x-1,y);\
      i3 = MT_GXY(sx,x,y-1);\
      i4 = MT_GXY(sx,x-1,y-1);\
      z2 = 2*z1-(dataZ)[i2];\
      z3 = 2*z1-(dataZ)[i3];\
      z4 = 2*z1-(dataZ)[i4];\
    }\
  }\
}

#define BMP_EVAL_BILP( z1, z2, z3, z4, /*scanline order */\
    		       u, v ) \
  ((1-(u))*(1-(v))*(z1) + (u)*(1-(v))*(z2) + \
  (1-(u))*(v)*(z3) + (u)*(v)*(z4)) 


#define BMP_RECT_AREA( r ) \
  (((r)->sx<0||(r)->sy<0) ? 0 : (r)->sx*(r)->sy)
  
#define BMP_RECT_IN_RECT(r1,r2) \
  (BMP_POINT_IN_RECT((r1)->x0,(r1)->y0,(r2)) && \
   BMP_POINT_IN_RECT((r1)->x0+(r1)->sx-1,(r1)->y0+(r1)->sy-1,(r2)))

#define BMP_POINT_IN_RECT(x,y,r)\
  BMP_IN_RANGE(x,y,(r)->x0,(r)->y0,(r)->x0+(r)->sx-1,(r)->y0+(r)->sy-1)

#define BMP_SET_RECT( r, x_, y_, sx_, sy_ ) \
{ \
  (r)->x0 = x_; (r)->y0 = y_; (r)->sx = sx_; (r)->sy=sy_; \
  (r)->x2 = (x_) + (sx_)-1; (r)->y2 = (y_) + (sy_)-1; \
}

#define BMP_SET_RECT2( r, x_, y_, x2_, y2_ ) \
{ \
  (r)->x0 = x_; (r)->y0 = y_; (r)->x2 = x2_; (r)->y2 = y2_;\
  (r)->sx = (x2_)-(x_)+1; (r)->sy = (y2_)-(y_)+1;\
}  
  
#define BMP_GET_RECT( r, x_, y_, sx_, sy_ ) \
{ \
  x_ = (r)->x0; y_=(r)->y0; sx_=(r)->sx; sy_=(r)->sy; \
}

#define BMP_GET_RECT2( r, x_, y_, sx_, sy_, x2_, y2_ ) \
{ \
  x_ = (r)->x0; y_=(r)->y0; sx_=(r)->sx; sy_=(r)->sy; \
  x2_ = (x_)+(sx_)-1; y2_=(y_)+(sy_)-1; \
}

#define BMP_GET_WIN( bmpSrc, psrcWin, x0_, y0_, sx_, sy_ ) \
{ \
  if( ! psrcWin ) \
  { \
    x0_ = 0; y0_ = 0; sx_=(bmpSrc)->sx; sy_=(bmpSrc)->sy; \
  } \
  else \
  { \
    BMP_GET_RECT( psrcWin,  x0_, y0_, sx_, sy_ ); \
  } \
}

#define BMP_GET_WIN2( bmpSrc, psrcWin, x0_, y0_, sx_, sy_, x2, y2 ) \
{ \
  BMP_GET_WIN(bmpSrc, psrcWin, x0_, y0_, sx_, sy_); \
  x2 = x0_+sx_-1;   y2 = y0_+sy_-1; \
}

#define BMP_SET_BITMAP_WINDOW( bmw_, bmp_, x0_, y0_, sx_, sy_ ) \
{ \
  (bmw_)->bmp = bmp_; (bmw_)->x0 = x0_; (bmw_)->y0 = y0_; \
  (bmw_)->sx = sx_; (bmw_)->sy = sy_; \
}

#define BMP_DIM_EQUAL( b1_, b2_ ) \
  ( ((b1_)->sx == (b2_)->sx) && ((b1_)->sy == (b2_)->sy) )

#define BMP_XOR_PIXEL( bmp, x, y, cr, cg, cb ) \
{ \
  bmp->data[(x)+bmp->sx*(y)] ^= BMP_MKRGB(cr,cg,cb); \
} 

#define BMP_GETPIXEL( bmp, x, y, c, r, g, b ) \
{ \
  c = bmp->data[(x)+bmp->sx*(y)]; \
  r = BMP_RED( c ); \
  g = BMP_GREEN( c ); \
  b = BMP_BLUE( c ); \
}

#define BMP_PIXEL_VALID(B,x,y) (BMP_PIXEL(B,dataValid,x,y)==0)

#define BMP_PIXEL(bmp,data_,x,y) (bmp)->data_[BMP_XYOFF(bmp,x,y)]

#define BMP_DPIXEL(bmp,data_,x,y) (data_)[BMP_XYOFF(bmp,x,y)]

#define BMP_GETPIXEL0( pp, bsx, x, y ) (*(pp+(x)+((LongInt)(bsx))*(y)))

#define BMP_GETPIXEL_GREY( bmp, x, y ) \
  (bmp->dataGrey[(x)+bmp->sx*(y)])

#define BMP_RGB2GREY( r,g,b) (0.3*(r)+0.59*(g)+0.11*(b))

#define BMP_GET_PP( bmp, x_, y_ ) (bmp->data+(x_)+bmp->sx*(y_)) 
#define BMP_GET_PP_GREY( bmp, x, y ) (bmp->dataGrey+(x)+bmp->sx*(y)) 
#define BMP_GET_PP_ALPHA( bmp, x, y ) (bmp->dataAlpha+(x)+bmp->sx*(y))
#define BMP_GET_PP_ZBUF( bmp, x, y ) (bmp->dataZbuf+(x)+bmp->sx*(y))
#define BMP_GET_PP_USHORT( bmp, x, y ) (bmp->dataUshort0+(x)+bmp->sx*(y)) 
#define BMP_GET_PP_UV( bmp, x, y ) (bmp->dataUV+(x)+bmp->sx*(y)) 
#define BMP_GETPIXEL_OFFS( bmp, x, y ) ((x)+(bmp)->sx*(y))
#define BMP_XYOFF(bmp,x,y) BMP_GETPIXEL_OFFS(bmp,x,y)
#define BMP_GETPIXEL_OFFSET( bmp_sx, x, y ) ((x)+((LongInt)(bmp_sx))*(y))
#define BMP_POFF2XY( bmp, po, x, y ) \
{ \
  y = ((int)(po)) / (bmp)->sx; \
  x = (po) - y*(bmp)->sx; \
}

#define BMP_IDX2COORDS( B, idx, /*IN*/\
    		        x, y  /*OUT*/)\
{ \
  y = (idx) / (LongInt)((B)->sx); \
  x = (idx) % (LongInt)((B)->sx); \
}

#define BMP_XYZ2RGB(vRgb,vXYZ) \
{ \
      vRgb[0] = vXYZ[0]/255; \
      vRgb[1] = vXYZ[1]/255; \
      vRgb[2] = vXYZ[2]/255; \
}

#define BMP_RGB2XYZ(vXYZ,vRgb) \
{ \
      vXYZ[0] = vRgb[0]*255; \
      vXYZ[1] = vRgb[1]*255; \
      vXYZ[2] = vRgb[2]*255; \
}

#define BMP_XYZ2LUV(vLuv,vXYZ) \
{ \
      vLuv[0] = vXYZ[0]/255; \
      vLuv[1] = 2.*vXYZ[1]/255.-1; \
      vLuv[2] = 2.*vXYZ[2]/255.-1.; \
}

#define BMP_LUV2XYZ(vXYZ,vLuv) \
{ \
      vXYZ[0] = 255.*vLuv[0]; \
      vXYZ[1] = 255.*(vLuv[1]+1)/2.; \
      vXYZ[2] = 255.*(vLuv[2]+1)/2.;  \
}

#define BMP_GETPIXEL_XYZ(bmp,x,y,vXYZ) \
{ \
  int p_ = BMP_GETPIXEL_OFFS(bmp,x,y); \
  vXYZ[0] = bmp->dataFloat[0][p_]; \
  vXYZ[1] = bmp->dataFloat[1][p_]; \
  vXYZ[2] = bmp->dataFloat[2][p_]; \
}

#define BMP_GETPIXELO_XYZ(bmp,po,vXYZ) \
{ \
  vXYZ[0] = bmp->dataFloat[0][po]; \
  vXYZ[1] = bmp->dataFloat[1][po]; \
  vXYZ[2] = bmp->dataFloat[2][po]; \
}

#define BMP_GETPIXEL_XYZ_LAP(bmp,x,y,vXYZ) \
{ \
  int p_ = BMP_GETPIXEL_OFFS(bmp,x,y); \
  vXYZ[0] = (double)bmp->dataFloat[0][p_] / (double)BMP_LAPL_SCALE; \
  vXYZ[1] = (double)bmp->dataFloat[1][p_] / (double)BMP_LAPL_SCALE; \
  vXYZ[2] = (double)bmp->dataFloat[2][p_] / (double)BMP_LAPL_SCALE; \
}

#define BMP_SETPIXEL_XYZ(bmp,x,y,vXYZ) \
{ \
  int p_ = BMP_GETPIXEL_OFFS(bmp,x,y); \
  bmp->dataFloat[0][p_] = vXYZ[0]; \
  bmp->dataFloat[1][p_] = vXYZ[1]; \
  bmp->dataFloat[2][p_] = vXYZ[2]; \
}

#define	BMP_GET_NEXTPIXEL( pp, c, cr, cg, cb ) \
{ \
  c = *pp++; \
  cr = BMP_RED( c ); \
  cg = BMP_GREEN( c ); \
  cb = BMP_BLUE( c ); \
}

#define	BMP_PP_GETPIXEL( pp, c, cr, cg, cb ) \
{ \
  c = *(pp); \
  cr = BMP_RED( c ); \
  cg = BMP_GREEN( c ); \
  cb = BMP_BLUE( c ); \
}

#define	BMP_PP_GETPIXEL_UV( pp, vLuv ) \
{ \
  vLuv[1] = (pp)->cu;\
  vLuv[2] = (pp)->cv;\
}

#define BMP_PP_GETPIXEL_LUV( ppg_, ppuv_, vLuv_ ) \
{ \
  vLuv_[0] = *(ppg_)/255.; \
  BMP_PP_GETPIXEL_UV( (ppuv_), vLuv_ ); \
}

#define BMP_GETPIXEL_LUV( bmp_, x_, y_, vLuv_ ) \
{ \
  BMP_PP_GETPIXEL_LUV( \
      BMP_GET_PP_GREY( (bmp_), (x_), (y_)), \
      BMP_GET_PP_UV( (bmp_), (x_), (y_)), \
      vLuv_ ); \
}

#define BMP_GETPIXEL_UV( bmp_, x_, y_, vLuv_ ) \
{ \
  BMP_PP_GETPIXEL_UV( \
      BMP_GET_PP_UV( (bmp_), (x_), (y_)), vLuv_);\
}
#define BMP_SETPIXEL_GREY(bmp_,x_,y_,l_) \
{ \
  *(BMP_GET_PP_GREY( bmp_, x_, y_ )) = l_; \
}

#define BMP_SETPIXEL_GREY_C(bmp_,x_,y_,l_) \
{ \
  *(BMP_GET_PP_GREY( bmp_, x_, y_ )) = MIN(MAX(l_,0),255); \
}

#define BMP_SETPIXEL_ALPHA(bmp_,x_,y_,l_) \
{ \
  *(BMP_GET_PP_ALPHA( bmp_, x_, y_ )) = l_; \
}

#define BMP_SETPIXEL_LUV( bmp_, x_, y_, vLuv_ ) \
{ \
  *(BMP_GET_PP_GREY( bmp_, x_, y_ )) = (vLuv_)[0]*255.; \
  BMP_MKUV((BMP_GET_PP_UV( bmp_, x_, y_ )),vLuv_); \
}

#define BMP_SETPIXEL_UV( bmp_, x_, y_, vLuv_ ) \
{ \
  BMP_MKUV((BMP_GET_PP_UV( bmp_, x_, y_ )),vLuv_); \
}

#define BMP_SET_NEXTPIXEL( pp, cr, cg, cb ) \
{ \
  *pp++ = BMP_MKRGB(cr,cg,cb); \
}

#define BMP_SET_NEXTPIXEL_GREY( ppg, l ) \
{ \
  *ppg++ = l; \
}

#define BMP_MAX_FCHAN		256
#define BMP_MAX_DCHAN		16
#define BMP_MAX_ICHAN		16
#define BMP_MAX_UCHAN		256
#define BMP_MAX_NAME_LEN	256

#define RSCAN_LOOP( pp_, x_, y_, bmpBase_, bmpOffs_, rx_, ry_, xmod_ ) \
      for( pp_=(bmpBase_)+(bmpOffs_), y_= 0; \
	   y_ < ry_; \
	   y_++, pp_+=xmod_ ) \
	for(x_=0; x_ < rx_ ;x_++, pp_++)

#define RSCAN_LOOP2( pp_, pp2_, x_, y_, bmpBase_, bmpBase2_, bmpOffs_, \
    		     rx_, ry_, xmod_ ) \
      for( pp_=(bmpBase_)+(bmpOffs_), pp2_=(bmpBase2_)+(bmpOffs_), y_= 0; \
	   y_ < ry_; \
	   y_++, pp_+=xmod_, pp2_+=xmod_ ) \
	for(x_=0; x_ < rx_ ;x_++, pp_++, pp2_++ )

#define BMP_RED( c ) 	((c) & 0xFF )
#define BMP_GREEN( c )	(((c) & 0xFF00 )>>8)
#define BMP_BLUE( c )   (((c) & 0xFF0000 )>>16)
#define BMP_MKRGB( cr, cg, cb ) \
 (0x02000000 | (((int)(cb))<<16) | (((int)(cg))<<8) | ((int)(cr)))

#define BMP_U( c ) ((c)->cu)
#define BMP_V( c ) ((c)->cv)
#define BMP_MKUV( c, vLuv ) { (c)->cu=vLuv[1]; (c)->cv=vLuv[2]; }

#define BMP_SETPIXEL( bmp, x, y, cr, cg, cb ) \
{ \
  bmp->data[(x)+bmp->sx*(y)] = BMP_MKRGB(cr,cg,cb); \
} 

#define BMP_SETPIXEL_RGB( bmp, x, y, crgb ) \
{ \
  bmp->data[(x)+bmp->sx*(y)] = crgb; \
} 

#define BMP_SETPIXEL_RGB_F( bmp, x, y, vRgb ) \
{ \
  int icr,icg,icb;\
  CO_RGB2INT(vRgb,icr,icg,icb);\
  bmp->data[(x)+bmp->sx*(y)] = BMP_MKRGB(icr,icg,icb); \
} 

#define CO_INT2RGB( vRgb, icr, icg, icb ) \
{ \
  vRgb[0]=(float)(icr)/255.0; \
  vRgb[1]=(float)(icg)/255.0; \
  vRgb[2]=(float)(icb)/255.0; \
}

#define BMP_GETPIXEL_RGB_F( bmp, x, y, vRgb ) \
{ \
  int ic,icr,icg,icb;\
  BMP_GETPIXEL(bmp,x,y,ic,icr,icg,icb);\
  CO_INT2RGB(vRgb,icr,icg,icb);\
} 

#define BMP_MIX_ALPHA_RGB( pp_, cr_, cg_, cb_, alpha_ ) \
{ \
  if( alpha_ < 255 ) \
  { \
    int c_ = *(pp_); \
    cr_ = (((alpha_)*(cr_) + (255-(alpha_))*BMP_RED(c_))>>8); \
    cg_ = (((alpha_)*(cg_) + (255-(alpha_))*BMP_GREEN(c_))>>8); \
    cb_ = (((alpha_)*(cb_) + (255-(alpha_))*BMP_BLUE(c_))>>8); \
  } \
}

#define BMP_TRANSCOLOR_PIXEL( bmp, x, y, cr, cg, cb ) \
{ \
  int *pp_=BMP_GET_PP( bmp, x, y ), lum8_, lh_, \
  	cr_=cr, cg_=cg, cb_=cb; \
  BMP_RGB24_TO_GREY( *pp_, lum8_ );  \
  TRANS_COLOR_SIMPLE0( cr_, cg_, cb_, lum8_, lh_ ); \
  *pp_ = BMP_MKRGB(cr_,cg_,cb_); \
}

#define BMP_PIXEL_IN_RECTANGLE( x0, y0, x, y, sx, sy ) \
( (x0) >= (x) && (x0) < (x)+(sx) && (y0)>=(y) && (y0)<(y)+(sy) )

#define BMP_RECT_CANONICALIZE( x1, y1, x2, y2, sx, sy ) \
{ \
  if( x1 > x2 || y1 > y2 ) \
  { \
    SWAP( x1, x2, sx ); \
    SWAP( y1, y2, sx ); \
  } \
  sx = x2-x1+1; \
  sy = y2-y1+1; \
}

#define BMP_DIM_EQ(b1,b2) ((b1)->sx==(b2)->sx && (b1)->sy==(b2)->sy)

#define BMP_RGB24_TO_GREY( c24_, gl_ ) \
{ \
  gl_ = RGB_TO_GREY( BMP_RED(c24_),BMP_GREEN(c24_), BMP_BLUE(c24_)); \
}

#define BMP_IN_RANGE(x,y,xmin,ymin,xmax,ymax) \
 ((x)>=(xmin)&&(x)<=(xmax)&&(y)>=(ymin)&&(y)<=(ymax))

#define BMP_IN_RANGER(x,y,r) \
 ((x)>=((r)->x0)&&(x)<=((r)->x2)&&\
  (y)>=((r)->y0)&&(y)<=((r)->y2))

#define BMP_IN_RANGEB(x,y,bmp) \
  BMP_IN_RANGE(x,y,0,0,((bmp)->sx-1),((bmp)->sy-1))


#define BMP_CLIP_CONST(x1,y1,x2,y2, /* IN: region */\
    		      x,y  /* OUT: coordantes to clip */\
    		     )\
{\
  if((x)>(x2)) x = x2; \
  if((x)<(x1)) x = x1; \
  if((y)>(y2)) y = y2; \
  if((y)<(y1)) y = y1; \
}

#define BMP_CLIP_WRAP(x1,y1,x2,y2, /* IN: region */\
    		      x,y  /* OUT: coordantes to clip */\
    		     )\
{\
  if((x)>(x2)) x = (x1)+(x)-(x2)-1; \
  if((x)<(x1)) x = (x2)+(x)-(x1)+1; \
  if((y)>(y2)) y = (y1)+(y)-(y2)-1; \
  if((y)<(y1)) y = (y2)+(y)-(y1)+1; \
}

#define BMP_CLIP_MOD(x1,y1,x2,y2, /* IN: region */\
    		      x,y  /* OUT: coordantes to clip */\
    		     )\
{\
  if((x)>(x2)) x = (x1) + ((x)-(x1)) % ((x2)-(x1)+1); \
  if((x)<(x1)) x = (x2) - ((x1)-(x)-1) % ((x2)-(x1)+1); \
  if((y)>(y2)) y = (y1) + ((y)-(y1)) % ((y2)-(y1)+1); \
  if((y)<(y1)) y = (y2) - ((y1)-(y)-1) % ((y2)-(y1)+1); \
}

#define BMP_CLIP_REFL(x1,y1,x2,y2, /* IN: region */\
    		      xc,yc  /* OUT: coordantes to clip */\
    		      )\
{\
  while((xc)<(x1)||(xc)>(x2)||(yc)<(y1)||(yc)>(y2)) \
  { \
    while((xc)<(x1)) xc = (x1) + ((x1)-(xc)); \
    while((xc)>(x2)) xc = (x2) - ((xc)-(x2)); \
    while((yc)<(y1)) yc = (y1) + ((y1)-(yc)); \
    while((yc)>(y2)) yc = (y2) - ((yc)-(y2)); \
  } \
}

#define BMP_CLIP_RECT_MP( csx, csy, x0, y0, sx, sy ) \
{ \
  if( x0 < 0 ) { x0 = 0; } \
  if( y0 < 0 ) { y0 = 0; } \
  if( (x0)+(sx) > (csx) ) { x0 = (csx) - (sx); } \
  if( (y0)+(sy) > (csy) ) { y0 = (csy) - (sy); } \
}

#define BMP_CLIP_MP( \
    	x0, y0, sx, sy, /* IN: target window */\
    	x1, y1, rx, ry  /* IN/OUT: window to be move-clipped */) \
	/* Assumtion: rx<=sx, ry<=sy */ \
{ \
  if( x1 < x0 ) { x1 = x0; } \
  if( y1 < y0 ) { y1 = y0; } \
  if( (x1)+(rx) > (x0)+(sx) ) { x1 = (x0)+(sx) - (rx); } \
  if( (y1)+(ry) > (y0)+(sy) ) { y1 = (y0)+(sy) - (ry); } \
}

#define CLIP_REGION_MP( bmw, \
    		     x0_, y0_, r, \
    		     x1, y1, rx, ry ) \
{ \
  rx = r; ry = r; \
  x1 = x0_ - rx/2; y1 = y0_ - ry/2; \
  if( x1 < (bmw)->x0 ) { x1 += (bmw)->x0 - x1; } \
  if( y1 < (bmw)->y0 ) { y1 += (bmw)->y0 - y1; } \
  if( x1+rx-1 >= (bmw)->x0+(bmw)->sx ) \
    { x1 -= x1+rx - ((bmw)->x0+(bmw)->sx); } \
  if( y1+ry-1 >= (bmw)->y0+(bmw)->sy ) \
    { y1 -= y1+ry - ((bmw)->y0+(bmw)->sy); } \
}

#define BMP_CLIP_RECT2( x0,y0,sx,sy, x1, y1, rx, ry ) \
{ \
  if( x1 < x0 ) { rx -= x0 - x1; x1 = x0; } \
  if( y1 < y0 ) { ry -= y0 - y1; y1 = y0;} \
  if( x1+rx-1 >= x0+sx ) \
    { rx -= x1+rx - (x0+sx);} \
  if( y1+ry-1 >= y0+sy ) \
    { ry -= y1+ry - (y0+sy);} \
}

#define BMP_CLIP_RECT0( r1_, r2_ ) \
  BMP_CLIP_RECT2( (r1_)->x0,(r1_)->y0,(r1_)->sx,(r1_)->sy, \
      		  (r2_)->x0,(r2_)->y0,(r2_)->sx,(r2_)->sy) 


/*
 *  Shift-Clip rectangular region x1,y1,rx,ry to bitmap <bmp>
 *  considering a maximal offset (each direction)
 *  of <pomax> pixels.
 *  Clipping is performed by shifting the region location.
 *  Assumption: rx,ry < bmp->sx,bmp->sy
 *  TODO: clip also rx,ry if region size greater than bitmap 
 */
#define BMP_CLIP_REGION( \
	bmp,     /* IN: bitmap */ \
	pomax,   /* IN: max. pixel offset */ \
	x1, y1,  /* IN/OUT: upper left corner of (shifted) region */ \
	rx, ry   /* IN/OUT: pixel width/height of  region) */ \
	) \
{ \
  int sxmax, symax; \
  sxmax=bmp->sx-2*(pomax); symax=bmp->sy-2*(pomax); \
  if(rx>sxmax) rx=sxmax; \
  if(ry>symax) ry=symax; \
  BMP_CLIP_MP(pomax,pomax,sxmax,symax,x1,y1,rx,ry ); \
}

/* deprecated */
#define BMP_CLIP_REGION_OLD( \
	    bmp,     /* IN: bitmap */ \
	    roi,     /* IN: region window (NULL=bitmap) */\
            pomax,   /* IN: max. pixel offset */ \
	    x1, y1,  /* OUT: upper left corner of shift-clipped region */ \
	    rx, ry   /* OUT: size of shift-clipped region) */ \
	    ) \
{ \
  int sxmax, symax; \
  BMP_GET_WIN( bmp, roi, x1, y1, rx, ry ); \
  sxmax=bmp->sx-2*pomax; symax=bmp->sy-2*pomax; \
  if(rx>sxmax) rx=sxmax; if(ry>symax) ry=symax; \
  BMP_CLIP_MP(pomax,pomax,sxmax,symax,x1,y1,rx,ry ); \
}


/* utility info for fast bitmap rectangle scan */

#define RSCAN_INFO( bsx, rx0, ry0, rsx, rsy, \
    			bmpOffs, xmod ) \
{ \
  bmpOffs =  bsx*(ry0) + (rx0); \
  xmod =  bsx - (rsx); \
}

#define BMP_RSCAN_INFO( bmp, x1, y1, rx_, ry_, \
    			rsi ) \
{ \
  (rsi)->rx = rx_; \
  (rsi)->ry = ry_; \
  (rsi)->bmpOffs =  ((bmp)->sx)*(y1) + (x1); \
  (rsi)->xmod =  (bmp)->sx - (rx_); \
}

#define BMP_RSCAN_INIT( bmp, pBmpData, x0, y0, sx_, rsmod, pp ) \
{ \
  pp = pBmpData + (bmp)->sx*(y0)+(x0); \
  rsmod = (bmp)->sx - (sx_); \
}

#define BMP_RSCAN_LOOP( pp1_, xmod1_, x_, y_, rx_, ry_ ) \
  for( y_= 0; y_ < ry_; y_++, pp1_+=xmod1_ ) \
	for(x_=0; x_ < rx_ ;x_++, pp1_++ )

#define BMP_RSCAN_LOOP2( pp1_, pp2_, xmod1_, xmod2_, x_, y_, rx_, ry_ ) \
  for( y_= 0; y_ < ry_; y_++, pp1_+=xmod1_, pp2_+=xmod2_ ) \
	for(x_=0; x_ < rx_ ;x_++, pp1_++, pp2_++ )

#define BMP_RSCAN_LOOP_Y( y_, sy_, upd_ ) \
	for(y_=0;y_<sy_;y_++, upd_ )

#define BMP_GET_SAMPLE( si_, ix_, iy_, x_, y_ ) \
{ \
  x_ = (si_)->x1 + (ix_) * (si_)->dx; \
  y_ = (si_)->y1 + (iy_) * (si_)->dy; \
  if( (si_)->doJittered ) \
  { \
    float r_ = 0.5*(MtRndSGet(&(si_)->rand) - 0.5); \
    x_ = (float)(x_)+(si_)->dx*r_; y_ = (float)(y_)+(si_)->dy*r_;\
    if((int)x_ < (si_)->x0) x_=(si_)->x0; \
    if((int)x_ > (si_)->x2) x_=(si_)->x2; \
    if((int)y_ < (si_)->y0) y_=(si_)->y0; \
    if((int)y_ > (si_)->y2) y_=(si_)->y2; \
  }\
}

#define BMP_GET_SAMPLE_INV( si_, x_, y_, ix_, iy_ ) \
{ \
  ix_ = ((x_) - (si_)->x1)/(si_)->dx; \
  iy_ = ((y_) - (si_)->y1)/(si_)->dy; \
  if((ix_)<0) ix_=0; \
  if((ix_)>=(si_)->nx) ix_ = (si_)->nx - 1; \
  if((iy_)<0) iy_=0; \
  if((iy_)>=(si_)->ny) iy_ = (si_)->ny - 1; \
}

#define BmpSetColConv( bmp_, cc_ ) \
{ \
  (bmp_)->colConv = cc_; \
}

#define BmpGetValidChannels( bmp_ ) (bmp_)->validChannels
#define BmpSetValidChannels( bmp_, v_ ) \
{ \
  (bmp_)->validChannels = v_; \
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
typedef   struct
{
  int chan,
      ichan;

  // data pointers
  uint8_t  *pUbyte;
  uint16_t *pUshort;
  float    *pFloat;

  // pointer addresses
  uint8_t  **ppUbyte;
  uint16_t **ppUshort;
  float    **ppFloat;
}
BmpData;

typedef float BmpGreylevel;

#define	BMP_GREYLEVEL_MAX	((BmpGreylevel)255)

typedef float BmpAlphaVal;
typedef unsigned char BmpZval;
typedef uint16_t BmpUshort;  // XXX do not change ! (must be uint16_t)

typedef struct
{
  int x0,y0,
      x2,y2;
  LongInt sx,sy;
}
BmpRectangle;

typedef struct
{
  float cu,
	cv;
}
BmpUV;

typedef struct
{
  struct BmpBitmap_s	*bmp0;
  int			 n_level,
			 chan, ichan,		// data channel
			 laplace_valid;
  unsigned		 opt;
  double		 res_0,
			 res_factor,
			 res_act;
  float	                *krn;
  int		         krn_r;
  struct BmpBitmap_s    *gaussian[BMP_MAX_PYRA_LEVEL],
			*laplacian[BMP_MAX_PYRA_LEVEL];
}
BmpPyramid;

/*
 * generic conversion between 3-dim color spaces 
 *
 *    for luma-chroma (Lxy) spaces:
 *       L normalized to [0..1]
 *	 x,y normalized to [-1..1]
 */
#define CO_OK			0
#define CO_CLIPPED		10
#define CO_HTAB_SIZE		20000
#define CO_CONV_MAX_NAME_LEN	128
#define CO_MAX_MGCS		1024
#define CO_CLIP_UV( vLxy_, vmin_, vmax_ ) \
{ \
  if((vLxy_)[1]<(vmin_)) (vLxy_)[1]=(vmin_); \
  if((vLxy_)[1]>(vmax_)) (vLxy_)[1]=(vmax_);\
  if((vLxy_)[2]<(vmin_)) (vLxy_)[2]=(vmin_); \
  if((vLxy_)[2]>(vmax_)) (vLxy_)[2]=(vmax_);\
}

int CoConvRgb2Hcl( float *vLxy, float *vRgb );
int CoConvHcl2Rgb( float *vRgb, float *vLxy );

typedef int (*CoConverterFunc)( float *triIn, 
    			        float *triOut ); 					      

typedef int (*CoConverterInitFunc)( void );

typedef struct
{
  char		  	name[CO_CONV_MAX_NAME_LEN+1];
  int		  	initDone;
  CoConverterInitFunc 	init;
  CoConverterFunc 	toRgb,
  		  	fromRgb;
  float 		minGamutCone[CO_MAX_MGCS+1],
  			rMaxMGC;
  int			nMGC;

  struct
  {
    float   **table;
    float    *centroids;
    int	      nSlices, nRad;
    float     rMax;
  }
  gamut;

}
CoConverter;

extern CoConverter CoConvDummy;

typedef struct BmpBitmap_s
{
  char		 	name[BMP_MAX_NAME_LEN+1];
  unsigned int		eyecatcher,guid;
  int		 	GWnr,
  		 	doRegenerate;
  LongInt		sx,sy,poffMax,
  			tsx,tsy;	/* optional: tile dimensions*/
  int			bitDepth;
  int 			vsx,vsy,vx0,vy0;  /* virtual size and pos. */
  unsigned		validChannels,
			optOrig;
  CoConverter		*colConv;	/* color space converter */
  int     		*data;
  BmpGreylevel  	*dataGrey;
  int			 dataGreyIsExt;
  BmpAlphaVal		*dataAlpha;
  BmpZval		*dataZbuf;
  unsigned short 	*dataUshort0;
  void		       **dataPtr;
  BmpUV			*dataUV;
  float			*dataFloat[BMP_MAX_FCHAN];
  int			 nFloatChan;
  double		*dataDbl[BMP_MAX_FCHAN];
  BmpUshort		*dataUshort[BMP_MAX_UCHAN];
  uint8_t		*dataUbyte[BMP_MAX_UCHAN];
  int			*dataInt[BMP_MAX_ICHAN];
  MtTab			*dataFeat;	/* generic pixel features */
  MtMat			*dataDFeat;
  float			*dataN[3];	 /* surface normals (nx,ny,nz)*/
  float			*dataZ;	 	 /* surface height */
  float			*dataRho;	 /* Albedo */
  unsigned char		*dataValid;
  MtStat		 stat;
  BmpPyramid		 pyramid;
  struct
  {  
    struct BmpBitmap_s    *bmpLeft,	// always points to parent
      			  *bmpRight,
			  *bmpAnaglyph;
  }
  stereo;
  
  struct
  {
    void	*header;	// ISIS cube header
    float	smpOffs,
		smpScale;
  }
  icu;
}
BmpBitmap;

typedef struct
{
  float  avg,var;
  int	 isValid;
}
BmpHistogram;

typedef struct
{
  BmpBitmap     *bmp;
  int 	      x0,y0,sx,sy;
}
BmpBitmapWindow;

typedef struct
{
  int    	 nVert;
  int		*vertices;
}
BmpPolygon;

typedef struct
{
  int nSetPixels,
      regSize;
  unsigned char *mask;
}
BmpRegionMask;

typedef struct
{
  int   bmpOffs,
        xmod,
        rx,ry;
  void *data;
}
BmpRscanInfo;

typedef struct BmpRegionInfo_s
{
  BmpBitmap	*bmp,*bmpReg,*bmpMatchCNH;
  BmpRegionMask *mask;
  BmpRscanInfo   rscanInfo;
  BmpRectangle	 window;
  int		 x0,y0,  /* midpoint */
  		 x1,y1;  /* left upper point of region to scan */
  int		 rContext;		 
  struct
  {
    BmpBitmap	*bmpUV;
    int		 x0,y0;
  }
  colInfo;
  struct BmpRegionInfo_s *riOrig;
}
BmpRegionInfo;

typedef int (*BmpMovWinFunctype)( int boff, int bmod, 
    				  BmpAlphaVal *bdata, 
				  void *info );

typedef BmpUV (*BmpMovWinUVFunctype)( int boff, int bmod, 
    				      BmpUV *bdata, 
				      void *info );

typedef int (*BmpMovWinRGBFunctype)( int boff, int bmod, 
    				      int *bdata, 
				      void *info );

typedef BmpGreylevel (*BmpMovWinGreyFunctype)( int boff, int bmod, 
    				      BmpGreylevel *bdata, 
				      void *info );
typedef struct
{
  int nMax,n,nx,ny,ix,iy,x0,y0,x2,y2;
  float	x,y,x1,y1,dx,dy;
  MtRndStream rand;
  int doJittered;

  int  nSmpXY;
  int *pSmpXY;
}
BmpSamplingInfo;

#define BMP_BITMAP_RAW_MAGIC 0x47111174
typedef struct
{
  uint32_t magic,
	   width,
	   height,
	   n_channels,
	   bit_depth;
}
BmpBitmapRawHeader;

//
int BmpWritePNG(FILE *fp, BmpBitmap *bmp, CoConverter *cc,
    		unsigned opt );

int BmpGetPixelNeighborhood(
    	LongInt sx, LongInt sy, float *data, int x, int y, int r, 
    	float *pn // OUT: pixel neighborhood (2*r+1)x(2*r+1)
		  //      in scanline order
	);

int BmpResampleBitmap( 
	  BmpBitmap *A,		// Input bitmap
	  int chan, int ichan,
	  double px, double py,	// output pixel size in mult. of A pixels
	  double ox, double oy,	// offset of B in A (in A pixels)
	  int	ip_method,	// interpolation method
	  BmpBitmap *B		// Output bitmap
    	);

int BmpResampleBitmapPixel( 
	  BmpBitmap *A,		// Input bitmap
	  float *data,
	  double *data_dbl,      // alternatively
	  double px, double py,	// output pixel size in mult. of A pixels
	  double ox, double oy,	// offset of B in A (in A pixels)
	  int	ip_method,	// interpolation method
	  int   btyp,
	  int   xb, int yb,
	  double *p_f_out	// interpolated output value
    	);

int BmpClipRectangle( 
    	BmpRectangle *A, // A: IN: rectangle to be clipped against
    	BmpRectangle *B  // B: IN/OUT: clipped rectangle
          );

int BmpNormalizeChannel( BmpBitmap *bmp, int chan, int ichan );

int BmpTestDetailGainFunc(int n_iter, BmpBitmap **pB);

float BmpLaplacianGain_FUNC1( 
    	float x, // input laplacian coefficient
	float s, // x scale factor
	float g, // overall gain factor
	float p, // power exponent
	float l, // lower limit of x*s for nonlinear enhancement
	float h  // upper limit of x*s for nonlinear enhancement
    );

int BmpUpSamplePyramid( BmpBitmap *B, int chan, int ichan,
    		        int l_from, int l_to );

int BmpUpSampleBitmap( BmpBitmap *B, int chan, int ichan, 
    		       int uplevel,
		       BmpBitmap **pBH // OUT: upsampled bitmap 
		       );

int BmpBlurGaussFast( BmpBitmap *B, int chan, int ichan, int level0,
    		      int btyp );

int BmpFreePyramid( BmpBitmap *bmp, unsigned opt );

int BmpFreePyramidChannels( BmpBitmap *B,  	  // target bitmap
		      	    int chan,	int ichan,	  // data channel
    		            unsigned opt );

int BmpCreatePyramid( BmpBitmap *B,  	  // target bitmap
		      int lmin,int lmax,  // min,max level (-1=auto)
		      int chan,	int ichan,	  // data channel
		      float *kern5x5_down,  // optional: 5x5 kernel
		      float *kern5x5_up,  // optional
		      int btyp,
    		      unsigned opt );

int BmpDisplayPyramid( BmpBitmap *B, 
    		       int chan, int ichan,
    		       char *save_fname,
    		       unsigned opt,
		       BmpBitmap **pO   // output bitmap 
		       );

int BmpAllocStereo( BmpBitmap *bmp, unsigned opt );
int BmpMakeAnaglyph( 
    	BmpBitmap *L, BmpBitmap *R,   
		// IN: left&right eye views: channels: XYZ as 32bit RGB 
	unsigned opt,
    	BmpBitmap *A 	// out: anaglyph
	);

BmpBitmap *BmpCreateBitmap2( char *inFileName, BmpBitmap *bmpSource,
    			int sx0, int sy0, int win, 
			char *name, CoConverter *colConv, 
			unsigned options );
#define BmpCreateBitmap( inFileName, bmpSource, sx0, sy0, win, name,\
    			 options ) \
  BmpCreateBitmap2( inFileName, bmpSource, sx0, sy0, win, name, NULL, \
    			 options ) 

BmpBitmap *BmpNewBitmap0( int sx, int sy, unsigned opt,
    			  char *mmTag );
#define BmpNewBitmap( sx, sy, opt )\
  BmpNewBitmap0( sx, sy, opt, MM_MK_TMP_UID())

int BmpDeleteBitmap( BmpBitmap **p_bmp );
int BmpCloseBitmap( BmpBitmap *bmp );
int BmpSaveBitmap( BmpBitmap *bmp, 
    		   int gwNr,
    		   char *filenameIn, 
		   int format,
		   CoConverter *cc,
		   unsigned opt );
int BmpSaveBitmapJPG( BmpBitmap *bmp, char *fname, int quality,
    			unsigned opt );
int BmpSaveBitmapGWBMP( BmpBitmap *bmp, int gwWin, char *filename );
int BmpLoadBitmap2( char *filename, 
    		       int chan, int ichan,     // main data channel
    		       unsigned opt,
    		       BmpBitmap **pBmp );
BmpBitmap *BmpLoadBitmap( char *fname, unsigned opt );
int BmpSaveBitmapEXR( BmpBitmap *B, char *filename, 
    		      int chan, int ichan,
    		      unsigned opt );
int BmpLoadBitmapEXR( char *filename, 
    		      unsigned opt,  
    		      BmpBitmap **pBmp );
int BmpLoadBitmapPNG( char *fname, unsigned chan,
    			     BmpBitmap **pBmp);
int BmpSaveBitmapPNG( BmpBitmap *bmp, char *fname, CoConverter *cc,
    			unsigned opt );

int BmpLoadBitmapTIFF( char *filename, 
    		       int chan, int ichan,    
    		       unsigned opt,
    		       BmpBitmap **pBmp );

int BmpSaveBitmapTIFF( BmpBitmap *bmp, char *fname, unsigned chan );
int BmpSaveBitmapTIFF2( BmpBitmap *bmp, char *fname, int chan, int ichan,
    			int bpp_out,  // output bits per pixel
    			unsigned opt );

int BmpSaveBitmapISIS( BmpBitmap *B,
    		       char *filename, 
    		       int chan, int ichan,
		       int tile_sx, int tile_sy,
		       unsigned opt );  

int BmpLoadBitmapISIS( char *filename, 
    		       int chan, int ichan,     // main data channel
		       int make_mask,	// 1= dataValid mask (0=ok,non-zero=error) 
		       			// 2= alpha-mask (UINT16_MAX=ok,0=invalid)
		       int ichan_mask,	
		       double valid_min, double valid_max,
		       int do_norm,
		       int specpix,	// special pixel handling: 0=set hi/lo to max/min and treat as valid	
		       unsigned opt,
    		       BmpBitmap **pBmp );

int BmpLoadBitmapRAW( char *filename, 
    		       int chan, int ichan,     // main data channel
		       unsigned opt,
    		       BmpBitmap **pBmp );

int BmpSaveBitmapRAW( BmpBitmap *B,
    		       char *filename, 
    		       int chan, int ichan, int nchan,
		       unsigned opt );  

int BmpSaveBitmapMBM( BmpBitmap *bmp, char *fname, unsigned chan );
int BmpLoadBitmapMBM( char *fname, unsigned chan,
    			     BmpBitmap **pBmp);
int BmpSaveBitmapUCF(char *fname, BmpBitmap *bmp, CoConverter *cc,
    		     unsigned opt );
float *BmpGetZChannel( BmpBitmap *bmp, unsigned opt );
int BmpGetNChannel( BmpBitmap *bmp, unsigned opt );
BmpAlphaVal *BmpGetAlphaChannel( BmpBitmap *bmp, unsigned opt );
void BmpFreeValidChannel( BmpBitmap *bmp );
void BmpFreeAlphaChannel( BmpBitmap *bmp );
void BmpFreeGreyChannel( BmpBitmap *bmp );	
void BmpFreeZChannel( BmpBitmap *bmp );
void BmpFreeDblChannel( BmpBitmap *bmp, int i );
void BmpFreeFloatChannel( BmpBitmap *bmp, int i );
void BmpFreeRgbChannel( BmpBitmap *bmp );
void BmpFreeUshortChannel( BmpBitmap *bmp, int i );
void BmpFreeUbyteChannel( BmpBitmap *bmp, int i );
BmpZval *BmpGetZbufChannel( BmpBitmap *bmp, unsigned opt );
unsigned char *BmpGetValidChannel( BmpBitmap *bmp, unsigned opt );
BmpGreylevel *BmpGetGreyChannel( BmpBitmap *bmp, unsigned opt );
BmpGreylevel *BmpSetGreyChannel( BmpBitmap *bmp, BmpGreylevel *dataGrey,
    				 unsigned opt );
double *BmpGetDblChannel( BmpBitmap *bmp, int i, unsigned opt );

float *BmpGetFloatChannel0( BmpBitmap *bmp, int i, unsigned opt,
    			    char *mmTag );
#define BmpGetFloatChannel( bmp, i, opt ) \
  BmpGetFloatChannel0(bmp,i,opt,MM_MK_TMP_UID())

int BmpGetFloatChannels( BmpBitmap *bmp, int n, unsigned opt );
int BmpGetXYZChannel( BmpBitmap *bmp, unsigned opt );
uint8_t *BmpGetUbyteChannel( BmpBitmap *bmp, int i, unsigned opt );
uint16_t *BmpGetUshortChannel( BmpBitmap *bmp, int i, unsigned opt );
void **BmpGetPtrChannel( BmpBitmap *bmp, unsigned opt );
int *BmpGetIntChannel( BmpBitmap *bmp, int i, unsigned opt );
BmpUV *BmpGetUVChannel( BmpBitmap *bmp, unsigned opt );
int *BmpGetRGBChannel( BmpBitmap *bmp, unsigned opt );
MtTab *BmpGetFeatureMap( BmpBitmap *bmp, int dim );
MtMat *BmpGetDFeatureMap( BmpBitmap *bmp, int dim );
int BmpCreateGreyscaleBitmap( BmpBitmap *bmp );
int BmpCopyDataChannel( BmpBitmap *S, int chan_s, int ichan_s, // Source
    			BmpBitmap *D, int chan_d, int ichan_d, // Destination
			unsigned opt );
int BmpGetDataPtrAddr(BmpBitmap *B, int chan, int ichan,
    		       unsigned opt,
    		       BmpData *data );
int BmpGetDataChannel(BmpBitmap *B, int chan, int ichan,
    		       unsigned opt,
    		       BmpData *data );

int BmpSwapDataPtr( BmpBitmap *B1, BmpBitmap *B2, 
    		    int chan, int ichan );

int BmpCopyBitmapPixels( BmpBitmap *src, BmpBitmap *dest );

int BmpCopyBitmapPixels2( BmpBitmap *src, BmpRectangle *srcWin,
    			  BmpBitmap *dst, BmpRectangle *dstWin,
			  unsigned channel );

int BmpPutToDisplayRGB( BmpBitmap *src, 	// source bitmap
    		 	 unsigned chan,		// source channel
		 	 int iChan,		// optional channel number
    		 	 BmpBitmap *dst, 	// destination bitmap
		 	 int xt, int yt         // target coordinates
		 	);

void BmpSetBitmapPixels( BmpBitmap *bmp, unsigned channel, unsigned val );

void BmpCopyAlphaToZbuf( BmpBitmap *bmpAlphaSrc, BmpBitmap *bmpZbufDst,
    			 BmpRectangle *srcWin, BmpRectangle *dstWin, 
			 int zval );

BmpRegionInfo *MiCreateRegionInfo( int regionSize,
    					 int nFeatures,
				   int contextRegionSize,	 	
    				         unsigned options );

void BmpDeleteRegionMask( BmpRegionMask *pmask );

BmpRegionMask *BmpCreateRegionMask( int regionSize,
    				      unsigned options );

void BmpPrintRegion( BmpRegionInfo *regInfo );

void MiDeleteRegionInfo( BmpRegionInfo *ri);

BmpPolygon *BmpCreatePolygon( int nVert, int *vertices );

int BmpFillPolygon( BmpBitmap *bmp, int channel, 
    		    BmpPolygon *poly, int val );

void BmpDeletePolygon( BmpPolygon *poly );

int BmpColorShift( BmpBitmapWindow *bmwIn, BmpBitmapWindow *bmwOut,
		   CoConverter *cc,
    		   float cx0, float cy0, float crx, float cry,
		   float cx2, float cy2 );

void BmpInitBitmapWindow( BmpBitmapWindow *bmw, 
    		     	 BmpBitmap *bmp,
    	                 int x0, int y0, int sx, int sy );

int BmpGetGradient( int type,
    			BmpBitmap *bmp, float *data,	
    			int x, int y, 
    			double h, double *vg );

void BmpGradientSoebel( BmpBitmap *bmp, float *data,	
    			int x, int y, 
    			double h, double *vg );

void BmpGradientCentDiff( BmpBitmap *bmp, float *data,
    			  int x, int y, 
    			   double h, double	*vg );

void BmpGradientForwDiff( BmpBitmap *bmp, float *data,
    			  int x, int y, 
    			   double h, double	*vg );

double BmpGetLumaGradient( BmpBitmap *bmp, int r,
    			   int x0, int y0,
    			   int doClip,
    			   double *pDirection );

int BmpFilterSoebel( BmpBitmap *src, BmpBitmap *dst, 
        	     int channel,
    		     BmpRectangle *srcWin, BmpRectangle *dstWin, 
    		     int r);

int BmpFilterSmooth( BmpBitmap *src, BmpBitmap *dst, 
        	     	  int channel, int iChannel,
    			  BmpRectangle *srcWin, BmpRectangle *dstWin, 
			  int r, float s, float shrp, unsigned opt );

int BmpFilterFFT( BmpBitmap *bs, BmpBitmap *bt,	/* source/target bitmaps */ 
        	  int chs, int cht,		/* source/target channels */
    	          BmpRectangle *win,		/* source window */ 
		  int xt, int yt );		/* target point	*/

int BmpFilterFFT2( BmpBitmap *bs, 		/* source bitmaps */ 
		   BmpBitmap **pbt, 		/* transformed bitmap */
        	   int chs,		        /* source channel */
    	           BmpRectangle *win,		/* source window */ 
		   int dir );		

int BmpTransformHistogram( BmpBitmap *bmpIn, BmpBitmap *bmpOut, 
    			   BmpBitmap *bmpRef, int maxSmp,
			   CoConverter *colorConverter,
			   unsigned	opt );

int BmpTransformHistogramOld( BmpBitmap *bmpIn, BmpBitmap *bmpOut, 
    			      BmpBitmap *bmpRef, int maxSmp,
			      CoConverter *colorConverter );

int BmpInitSampling( BmpBitmap *bmp, BmpRectangle *win, 
    		     float dx, float dy, 
		        /* (avg.) pixel distance in sample grid */
		     int regionSize, /*TODO*/
		        /* leave <regionSize>/2 pixels clearance border */
		     int maxSamples, 
		     	/* if not zero: limit number of samples */
		     unsigned opt,
		        /* BMP_JITTER (default: uniform) */
		     BmpSamplingInfo *si );

int BmpCloseSamplingInfo(BmpSamplingInfo *si);

MtTab *BmpSampleBitmapOld( BmpBitmap *bmp, int nmax,
    		     	unsigned opt );

MtTab *BmpSampleBitmap( BmpBitmap *bmp, int nmax, float sd,
    			CoConverter *cc,
    		     	unsigned opt );
int BmpRgbToLuv( BmpBitmap *bmpSrc, BmpBitmap *bmp, 
    		 CoConverter *colConv, unsigned opt );

int BmpLuvToRgb( BmpBitmap *bmpSrc, BmpBitmap *bmp, 
    		 CoConverter *colConv, unsigned opt );

int BmpRgbToGrey( BmpBitmap *bmpSrc, BmpBitmap *bmp, 
    		  CoConverter *colConv, unsigned opt );

int BmpGreyToRgb( BmpBitmap *bmpSrc, BmpBitmap *bmp, 
    		  CoConverter *colConv, unsigned opt );

int BmpGreyToRgbFast( BmpBitmap *bmpSrc, BmpBitmap *bmp );

BmpBitmap *BmpResizeBitmap( BmpBitmap *bSrc, int sx, int sy, unsigned opt );

int BmpPutBitmap( BmpBitmap *bmpFrom, 
    			  int x0, int y0,
			  int chan_in, int ichan_in,
			  int chan_out, int ichan_out, 
			  unsigned opt,
			  BmpBitmap *bmpTo );

BmpBitmap *BmpMakeBorderPaddedBitmap( BmpBitmap *bmpIn, 
    				      int borderSz,
    				      int borderExpTyp,
				      int channel );

int BmpSubtractBitmapUV( BmpBitmap *bs, BmpBitmap *bt, int ch );			    
int BmpDiffBitmaps( BmpBitmap *b1, BmpBitmap *b2, int channel,
    		    MtStat *errStat );

int BmpGetPixStats2( BmpBitmap *B, 
    		     int chan, int ichan,
    		    int x1, int y1, 
    		    int x2, int y2,
		    unsigned opt,
		    MtStat *stat );

int BmpGetPixStats( BmpBitmap *B, 
    		    float *data,
    		    int x1, int y1, 
    		    int x2, int y2,
		    MtStat *stat );

BmpBitmap *BmpCombineRGB( BmpBitmap *b1, BmpBitmap *b2, BmpBitmap *b3,
    		   	  float *cr, float *cg, float *cb, unsigned opt );

int BmpLuvLintrans( BmpBitmap *bs, BmpBitmap *bt, BmpRectangle *win,
 		    float *vLuvFactor, float *vLuvOffset,
		    float uvRot, float uflip, float vflip, 
		    float lalpha, float lgamma,
		    void (*nonLinFunc)( float *x, int dim ),
		    float logf,
    		    CoConverter *cc, unsigned opt );

int BmpLuvLintrans2( BmpBitmap *bs, BmpBitmap *bt, BmpRectangle *win,
 		    float *vcu, float *vcv,
		    void (*nonLinFunc)( float *x, int dim ),
    		    CoConverter *cc, float logf, unsigned opt );

int BmpSaturationLogCorrect( 
    	BmpBitmap *bs, BmpBitmap *bt, BmpRectangle *win,
	float a, int nIter, float gamma, 
	CoConverter *cc, unsigned opt );

int BmpSaturationGammaCorrect( 
    	BmpBitmap *bs, BmpBitmap *bt, BmpRectangle *win,
	float gamma, float alpha, CoConverter *cc, unsigned opt );
int BmpSaturationClip( 
    	BmpBitmap *bs, BmpBitmap *bt, BmpRectangle *win,
	float a, float b, float c, 
	CoConverter *cc, unsigned opt );
int BmpSaturationInvert( 
    	BmpBitmap *bs, BmpBitmap *bt, BmpRectangle *win,
	CoConverter *cc, unsigned opt );
int BmpLocalColorStretch( 
    	BmpBitmap *bs, BmpBitmap *bt, BmpRectangle *win,
	int r, float s,
	CoConverter *cc, unsigned opt );

int BmpDrawCircle( BmpBitmap *bmp, int chan,
    		  int x0, int y0, int dmax,
		  int crgb );

int BmpDrawLine2( BmpBitmap *bmp, int chan,
    		   int x1,  int y1, 
		   int x2,  int y2, 
		  int col ); 

int BmpDrawLine( BmpBitmap *bmp, int chan,
    		 int x1, int y1, int x2, int y2, int col, 
    		 int endp ); 

int BmpDrawArrow( BmpBitmap *bmp, int chan,
    		  int x1, int y1, int x2, int y2, int col, 
    		  float head_size );

int BmpDrawText( BmpBitmap *bmp, int x0, int y0, char *str, int psz );

int BmpAddNoise( 
    	BmpBitmap *bs, BmpBitmap *bt, BmpRectangle *win,
	float sigma, CoConverter *cc, unsigned opt );

int BmpCorrectVertLumiBands( 
    	BmpBitmap *bs, BmpBitmap *bt,
	int  x0, int y0, int dly, 
	CoConverter *cc, unsigned opt );

int BmpDrawCross( BmpBitmap *bmp, int chan,
    		  int x0, int y0, int dmax,
		  int crgb );

int BmpDrawRect( BmpBitmap *bmp, int channel, 
    		 int x1, int y1, int sx, int sy,
		 int cr, int cg, int cb);

int BmpDrawFilledRect( BmpBitmap *bmp, int channel, 
    		 int x1, int y1, int sx, int sy,
		 int cr, int cg, int cb);

int BmpEnhanceDiff( BmpBitmap *bs, BmpBitmap *bt, float alpha,
    		    CoConverter *cc );
int BmpEnhanceDiff2( BmpBitmap *bs, BmpBitmap *bt, 
    		     BmpRectangle *win, /* target window */
		     int ox, int dx,	/* background win offset/width */
    		     float alpha,
    		     CoConverter *cc );
int BmpConvert16BitGrey( BmpBitmap *bs,
		         float alpha, float gamma, float beta,
    			 CoConverter *cc, unsigned opt );
BmpBitmap *BmpMakeKernel(int r, float s, float shrp, unsigned opt );
int BmpMakeXYZ( BmpBitmap *bs, BmpBitmap *bt, int colmod,
    		int linearGreyMap );
int BmpXYZToRGB( BmpBitmap *bs, BmpBitmap *bt, int colmod,
    		 float *meanXYZ,
    		 CoConverter *cc );
int BmpXYZ2RGB( BmpBitmap *bs, BmpBitmap *bt, CoConverter *cc, unsigned opt );
int BmpRGB2XYZ( BmpBitmap *bs, BmpBitmap *bt, CoConverter *cc, unsigned opt );
int BmpXYZToLUV( BmpBitmap *bs, BmpBitmap *bt, int colmod,
    		 CoConverter *cc );
BmpBitmap *BmpGradientSobelXYM( BmpBitmap *bmp, BmpRectangle *roi,
    				float *dataChannel,
    				int radius );
int BmpCreatePyramidOld( BmpBitmap *bmp0,
    		      float	res_factor,
		      	// Resolution reduction factor (default: 2)
		      float	a_kern,
		        // Generating kernel 'a' parameter (-> Burt 1983)
		      int maxLevel,
    		      unsigned opt );

float *BmpGetDfltPyraConvKernel5x5(void);

int BmpGetBurtAdelsonKernel5x5(
    		float a, // IN: smoothness/shape parameter  
			 //   0.6=neg.lobes, 0.5=triangle, 
			 //   0.4=gauss, 0.3=broader gauss
			 //   see also [Burt1983]
    		float *kernel5x5  // OUT: 5x5 kernel 
		);

int BmpGetPyramid( BmpBitmap *bmp, unsigned opt );

int BmpHasPyramid( BmpBitmap *bmp, int opt );

int BmpFindPyramidLevel( BmpBitmap *bmp0, int md0, int sxy0,
    			 int *pLevel );

BmpBitmap *BmpIExpandPyramid( BmpBitmap *bmp0,
    			      BmpBitmap *bmpSrc,
			      BmpBitmap *bmpDst,
    		              int levelDst,
			      	// target level
			      float res_factor,
		              unsigned opt );

BmpBitmap *BmpExpandPyramidOld( BmpBitmap *bmp0,
    			BmpBitmap *bmpFrom,
    		      int levelFrom,
		      int levelTo,
		      unsigned opt
		      );

int BmpDownSamplePyraLevelMaxZ( 
    		       BmpBitmap *bmpLo, float *dataLo, 
    		       BmpBitmap *bmpHi, float *dataHi );

int BmpDownSamplePyraLevelFast( 
	  BmpBitmap *bmpLo, float *dataLo, BmpUshort *dataLoUs, uint8_t *dataLoUb,
	  BmpBitmap *bmpHi, float *dataHi, BmpUshort *dataHiUs, uint8_t *dataHiUb,
	  int oxh, int oyh, // optional bmpHi offset 0, or -1 (for MCM)  
          float *kernel5x5,  // optional 5x5 downsampling kernel
	  int btyp	 // boundary handling
	  );

int BmpDownSamplePyraLevel( BmpBitmap *bmp0,
    		       BmpBitmap *bmpLo, float *dataLo, double *dataLo_d,
    		       BmpBitmap *bmpHi, float *dataHi, double *dataHi_d,
		       int	doPresVariance );

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
	);

int BmpUpSamplePyraLevelFast( 
    	BmpBitmap *bmpLo, float *dataLo, BmpUshort *dataLoUs,
    	BmpBitmap *bmpHi, float *dataHi, BmpUshort *dataHiUs,
	float *dataLapl, 
	int oxh, int oyh,
        float *kernel5x5,  // optional 5x5 upsampling kernel
	int btyp // boundary handling 
	);

int BmpUpSamplePyraLevel( BmpBitmap *bmp0,
    		           BmpBitmap *bmpLo, 
			   float *dataLo, double *dataLo_d,
    		           BmpBitmap *bmpHi, 
			   float *dataHi, double *dataHi_d,
		           double res_factor_f, 
			   float *dataLapl, 
          		   float *kernel5x5  // optional 5x5 downsampling kernel
			   );

BmpBitmap *BmpReconPyramid( BmpBitmap *bmp0, unsigned opt );

#ifdef __cplusplus
}
#endif

MSBG_NAMESPACE_END

#endif /* BITMAP_H */

