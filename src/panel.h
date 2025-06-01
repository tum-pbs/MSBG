/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#ifndef PANEL_H
#define PANEL_H

MSBG_NAMESPACE_BEGIN

#define PNL_MAXIMIZE		(1<<0)
#define PNL_RESTORE		(1<<1)
#define PNL_ACTIVATE		(1<<2)
#define PNL_CREATE_BMP_SAVE	(1<<3)
#define PNL_HIGHLIGHT		(1<<4)
#define PNL_XOR			(1<<5)
#define PNL_CLIPBOARD		(1<<6)	
#define PNL_GREY_ONLY		(1<<7)
#define PNL_IGNAN		(1<<8)
#define PNL_RGB			(1<<9)
#define PNL_PROVIDE_SLABSTATS	(1<<10)
#define PNL_NOTEXT		(1<<11)
#define PNL_CLEAR		(1<<12)

#define GWK_LMOUSE	-1
#define GWK_RMOUSE	-2
#define GWK_INS		45
#define GWK_DEL		46
#define GWK_TAB		9
#define GWK_PLUS	187
#define GWK_MINUS	189
#define GWK_F1		112
#define GWK_F2		113
#define GWK_F3		114
#define GWK_F4		115
#define GWK_F5		116
#define GWK_F6		117
#define GWK_F7		118
#define GWK_F8		119
#define GWK_F9		120
#define GWK_COMMA	188
#define GWK_DOT		190
#define GWK_ENTER	13
#define GWK_ROOF	220

#define GWK_PGUP	33
#define GWK_PGDOWN	34
#define GWK_HOME 	36
#define GWK_END		35
#define GWK_INSERT	45

#define GWK_UP		38
#define GWK_DOWN	40
#define GWK_LEFT	37
#define GWK_RIGHT	39

#define GW_WHITE	19

//
typedef struct
{
  int	    	 id,
		 win,
		 bmpGw, bmpGw_rr,
		 bmpGw2, bmpGw2_rr,
  		 undoCount;
  BmpRectangle   bmp_roi;
  BmpBitmap   	*bitmap,
		*bitmap2,
  		*bitmapSave,
		*bitmapTmp;
  BmpBitmapWindow bmw0;
  float		 wx0,wy0,wx1,wy1;
  float	    	 x0,y0,fx,fy,wsx,wsy;
  float	    	 x02,y02,fx2,fy2;
  PlScreen   	 plScreen;
  PlWindow   	 plwin[4],plwin0,plwin1,plwinIso;
  PlPlot	 plplot[4],plplot0;
  char		 bmpFileName[MI_MAX_FNAME_LEN+1];
  char		 title[1024];
  double 	 totalTime;
}
PnlPanel;


#ifdef __cplusplus
extern "C" {
#endif

void PnlVisualize3DFieldSlices( 
		PnlPanel *pnl,		// Output panel
    		float **field, 		// 3D field 
		float **slabZY, float **slabXZ, float **slabXY,  // optional slabs
		int nChan,		// Number of data channels 
		LongInt sx, LongInt sy, LongInt sz,  // resolution
		double x0, double y0, double z0,	// slab positions 
		MtStat *gstat,		// optional: global statistics
		unsigned options,
		BmpBitmap **pOutBmp
		);

PnlPanel *PnlActPanel();
void PnlClear( PnlPanel *pnl );
void PnlFree( PnlPanel *pnl );
void PnlInit( PnlPanel *pnl, int id, int win);
void PnlFreeBitmap( PnlPanel *pnl );
void PnlUpdateBitmapRoi( PnlPanel *pnl, BmpRectangle *roi );
void PnlSetBitmap( PnlPanel *pnl, BmpBitmap *bmp, BmpBitmap *bmpSave,
    		       char *fname );
int PnlDisplayBmp( PnlPanel *pnl, int force );
int PnlDisplayBmp0( PnlPanel *pnl, BmpBitmap *bmp, int force );
int PnlDisplayBmp2( PnlPanel *pnl, BmpBitmap *bmp, float x0, float y0,
    		    float wsx, float wsy, /* in percent */
    		    int force );
void PnlDrawRect( PnlPanel *pnl, int x0, int y0, int rx, int ry,
    		    int cr, int cg, int cb, int mode );
void PnlDrawPixel( PnlPanel *pnl, int x0, int y0,
    		    int cr, int cg, int cb, int mode );
void PnlDrawLine( PnlPanel *pnl, int x1, int y1, int x2, int y2,
    		    int cr, int cg, int cb, int mode );
void PnlBmpCoordsInv( PnlPanel *pnl, int px, int py,
    			     float *pwx, float *pwy );
void PnlBmpCoords( PnlPanel *pnl, float wx, float wy,
    			  int *px, int *py );
void PnlSelect( PnlPanel *pnl, unsigned options );
int PnlLoadBitmap( char *fnameIn0, PnlPanel *pnl, BmpBitmap *bmpIn,
    		   CoConverter *colConv,
    			    unsigned options, int resizePerc  );
int PnlSaveBitmap( char *fnameIn0, PnlPanel *pnl, CoConverter *cc );
int PnlCapturePolygon( PnlPanel *pnl, char *msg, BmpPolygon **ppoly );
void PnlErase( PnlPanel *pnl );
PnlPanel *PnlGetActPanel(void);
void PnlSetTmpBitmap( PnlPanel *pnl, BmpBitmap *bmp );
int PnlShowBitmap( PnlPanel *pnl, BmpBitmap *bmp );
int PnlShowBitmap2( PnlPanel *pnl, BmpBitmap *bmp );
int PnlShowBitmap3( PnlPanel *pnl, BmpBitmap *bmp, 
    		    float gmax, //0: do norm
    		    int chan, int ichan,
		    BmpRectangle *roi,
		    unsigned opt );
PnlPanel *PnlCreate( char *title );
PnlPanel *PnlCreate2( char *title, int wsx, int wsy, unsigned opt );
int PnlDelete( PnlPanel *pnl );
void PnlBmpCoords2( PnlPanel *pnl, float wx, float wy,
    			  int *px, int *py );
void PnlBmpCoordsInv2( PnlPanel *pnl, float px, float py,
    			     float *pwx, float *pwy );

// external
PnlPanel *MiGlobalPnlOutput( void );
PnlPanel *MiGlobalPnlTool( void );
PnlPanel *MiGlobalPnlInput( void );
PnlPanel *MiGlobalPnlRef( void );
PnlPanel *MiGetAuxPanel(int *pi);
PnlPanel *MiGetAuxPanel2(int *pi, const char *title);

#ifdef __cplusplus
}
#endif


#ifdef __cplusplus

template<typename func_T>
void PnlPlotFunc( PnlPanel *pnl, 
          	 double x1, double x2, 
		 double y1, double y2, 
       		 func_T func, 
		 int nStep=256,
		 int doClear=TRUE, int col=-1 )
		    
{
  int i;
  PlPlot pl; PlWindow plw;
  double dx,x;

  if(col<0) col=GWkrgb(255,255,200);

  PnlSelect(pnl,PNL_ACTIVATE);
  if(doClear) PnlClear(pnl);

  PlWindowInit(&plw,&pnl->plScreen,0,0,100,100);
  PlPlotInit(&pl,&plw,x1,x2,y1,y2,0);
  PlColor(GWkrgb(90,90,90)); 
  PlPlotAxes( &pl, PL_AXDEFAULT, 15, 5 );

  dx = (x2-x1)/nStep;
  for(i=0;i<nStep;i++)
  {
    x = x1+dx*i;
    double y = func(x);
    PlPlotDot2(&pl,x,y,col); 
  }
}

#endif


MSBG_NAMESPACE_END

#endif /* PANEL_H */

