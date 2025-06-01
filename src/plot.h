/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#ifndef PLOT_H
#define PLOT_H

MSBG_NAMESPACE_BEGIN

#ifdef __cplusplus
extern "C" {
#endif

#define PL_LINE_START  (1<<1)
 
#define PL_AUTOSPC_X  (1<<1)
#define PL_AUTOSPC_Y  (1<<2)
#define PL_TEXT_X     (1<<3)
#define PL_TEXT_Y     (1<<4)
#define PL_LOGSCALE_X (1<<5)
#define PL_LOGSCALE_Y (1<<6)
#define PL_FRAME      (1<<7)

#define PL_LINES      (1<<1)
#define PL_DOTS	      (1<<2)	
#define PL_NOAXES     (1<<3)
#define PL_STEP       (1<<4)
#define PL_FATDOTS    (1<<5)	
#define PL_FIRST      (1<<6)
#define PL_LAST	      (1<<7)

#define PL_AXDEFAULT  ( PL_AUTOSPC_X | PL_AUTOSPC_Y | PL_TEXT_X |\
                        PL_TEXT_Y | PL_FRAME ) 

#define PlSetPen( col_, style_, width_, mode_ ) \
  GWsetpen(col_,style_,width_,mode_);

#define PLRGB( col )   \
  GWkrgb(  PlColNib2Byte(((col)&0xF00)>>8),\
      			PlColNib2Byte(((col)&0x0F0)>>4),\
      			PlColNib2Byte((col)&0x00F))

typedef double (* PlTransFunc)( double x );

typedef struct
{
  int x0,y0,
      sx,sy;
}
PlScreen;

typedef struct
{
  PlScreen  *screen;
  int        x0,y0,
             sx,sy;
}
PlWindow;

typedef struct
{
  PlWindow   *window; 
  void	     *usrInfo;
          
  double      xScale, yScale,     /* scale */
              xOff, yOff,         /* offset */
              xMin, xMax,   
              yMin, yMax,
              xSpan,
              ySpan;
             
  int         flags;

  PlTransFunc transX, transY,
              transRevX, transRevY;                       
}
PlPlot;

typedef struct
{
  int pxOld, pyOld, 
      first;
}
PlLineHd;

void 	PlPlotLineInit( PlLineHd *plHd ) ;
int 	PlPlotLine( PlPlot *p, double x, double y, PlLineHd *plHd ); 
int PlPlotLine2( PlPlot *p, float x1, float y1, 
    			    float x2, float y2 ); 
int PlPlotPolyLine( PlPlot *p, double *px, double *py, int n );
int     PlPlotDot( PlPlot *p, double x, double y );
int PlPlotDot2( PlPlot *p, double x, double y, int col ); 
int PlPlotCross( PlPlot *p, double x, double y, 
    			   double sx, double sy );
int PlPlotCrossX( PlPlot *p, double x, double y, 
    			   double sx, double sy );
int PlPlotPixel( PlPlot *p, double x, double y, int col ); 
int PlPlotText( PlPlot *p, double x, double y, int col, char *str );
int 	PlPlotRect( PlPlot *p, double x, double y, 
    			   double sx, double sy );
int PlPlotRect2( PlPlot *p, float x1, float y1, 
    			    float x2, float y2 ) ;

int PlPlotCircle( PlPlot *p, double x, double y, double r );
int PlPlotArrow( PlPlot *p, float x1, float y1, 
    			    float x2, float y2,
			    float head_size
			    ) ;
int PlPlotFilledRect( PlPlot *p, double x, double y, 
    		      double sx, double sy,
		      int col );
int     PlPlotDotRect( PlPlot *p, double x, double y, int sx, int sy ); 
int     PlPlotGraph( PlPlot *p, MtVector vx, MtVector vy, int n,int style ); 
int 	PlPlotInit( PlPlot *p, PlWindow *w, double xMin, double xMax, 
        	    double yMin, double yMax, int flags );
int     PlPlotAxes( PlPlot *p, int flags, double xSpace, double ySpace );
int     PlPlotCursor( PlPlot *p, double x, double y, int flags );
int 	PlWindowInit( PlWindow *w, PlScreen *s, double xp0, double yp0, 
        	      double spx, double spy );
int     PlWindowClear( PlWindow *plw, int col );
int 	PlScreenInit( PlScreen *s, int x0, int y0, int sx, int sy );
void 	PlColor( int color ); 
void 	PlRefreshScreen( void ); 
void    PlClearScreen( void );
void    PlPlotReverse( PlPlot *p, int px, int py, double *x, double *y );
int 	PlColNib2Byte( int a );

#ifdef __cplusplus
}  // extern "C"
#endif

MSBG_NAMESPACE_END

#endif
