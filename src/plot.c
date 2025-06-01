/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include "globdef.h"
#include "mtool.h"
#include "plot.h"

static int roundupTable[10] = { 1, 2, 2, 5, 5, 5, 5, 10, 10, 10 };
  
/********************************************************************/
/*								    */
/********************************************************************/
static double PlITransFuncIdent( double x )
{
return x;
}

#if 1
/********************************************************************/
/*								    */
/********************************************************************/
static double PlITransFuncLog( double x )
{
return log( x );
}

/********************************************************************/
/*								    */
/********************************************************************/
static double PlITransFuncExp( double x )
{
return exp( x );
}
#endif
      
/********************************************************************/
/*								    */
/********************************************************************/
int PlScreenInit( PlScreen *s, int x0, int y0, int sx, int sy )
{
  s->x0 = x0;
  s->y0 = y0;
  s->sx = sx;
  s->sy = sy;
  return SUCCESS;
}

/********************************************************************/
/*								    */
/********************************************************************/
int PlWindowInit( PlWindow *w, PlScreen *s, double xp0, double yp0, 
                 double spx, double spy )
{
  /* printf("PlWindowInit %7.1f %7.1f %7.1f %7.1f\n",xp0,yp0,spx,spy); */
  w->screen = s;
  w->x0 = ( xp0 / 100. ) * s->sx;
  w->y0 = ( yp0 / 100. ) * s->sy;
  w->sx = ( spx / 100. ) * s->sx;
  w->sy = ( spy / 100. ) * s->sy;
  return SUCCESS;  
}

/********************************************************************/
/*								    */
/********************************************************************/
int PlWindowClear( PlWindow *w, int col )
{
  int px=w->x0,py=w->y0,psx=w->sx,psy=w->sy;

  GWsrect( px, py, px+psx-1, py+psy-1, col );
   
  return SUCCESS;
}


/********************************************************************/
/*								    */
/********************************************************************/
int PlPlotInit( PlPlot *p, PlWindow *w, double xMin, double xMax, 
                double yMin, double yMax, int flags )
{
  PlScreen  *s;     

  if( flags & PL_LOGSCALE_X ) 
  {
    p->transX = PlITransFuncLog;
    p->transRevX = PlITransFuncExp;
  }  
  else
  {
    p->transX = PlITransFuncIdent;
    p->transRevX = PlITransFuncIdent;
  }  

  if( flags & PL_LOGSCALE_Y ) 
  {
    p->transY = PlITransFuncLog;
    p->transRevY = PlITransFuncExp;
  }  
  else
  {
    p->transY = PlITransFuncIdent;
    p->transRevY = PlITransFuncIdent;
  }  
  
  xMin = p->transX( xMin );
  xMax = p->transX( xMax );
  yMin = p->transY( yMin );
  yMax = p->transY( yMax );
  
  s = w->screen;
  p->window = w;
  p->xMin =  xMin;
  p->xMax =  xMax;
  p->yMin =  yMin;
  p->yMax =  yMax;

  p->xSpan = p->xMax - p->xMin;
  p->ySpan = p->yMax - p->yMin;
  
  p->xScale = p->xSpan > MT_REAL_EPSILON ? 
    		((double) ( w->sx )) / p->xSpan : 0;
  p->xOff = - p->xScale * p->xMin + w->x0 + s->x0; 
    
  p->yScale = p->ySpan > MT_REAL_EPSILON ? 
    		((double) ( w->sy )) / p->ySpan : 0;
  p->yOff = - p->yScale * p->yMin + w->y0 + s->y0;
  
  p->flags = flags;
  
            
  return SUCCESS;
}

/********************************************************************/
/*								    */
/********************************************************************/
int PlPlotCursor( PlPlot *p, double x, double y, int flags )
{
  int        px, py1, py2;

  px = x * p->xScale + p->xOff;
  py1 = p->yMin * p->yScale + p->yOff;
  py2 = p->yMax * p->yScale + p->yOff;
  GWline( px, py1, px, py2 );
  return SUCCESS;
}

/********************************************************************/
/*								    */
/********************************************************************/
int PlPlotAxes( PlPlot *p, int flags, double xSpace, double ySpace )
{
  int        i, px1, py1, px2, py2, dpx, dpy, mod;
  double     ySpan, xN, yN, step, tmp, xSpc, ySpc, x, y;
  char       str[256];
  
  
  dpx = 0;
  dpy = 0;
  
  if( flags & PL_AUTOSPC_X ) 
  {
    xN = xSpace;
    step = p->xSpan / xN;
    tmp = pow( 10., floor( log10( step )));
    if(tmp>MT_NUM_EPS)
    {
      i = (floor)(( step / tmp ) + .5 ) ;
      i = MIN(MAX(i,0),ARRAY_LENGTH(roundupTable)-1);
      xSpc = roundupTable[i-1] * tmp; 
    }
    else xSpc = xSpace;
  }
  else
  {
    xSpc = xSpace;
  }

  if( flags & PL_AUTOSPC_Y ) 
  {  
    yN = ySpace;
    ySpan = p->transRevY( p->yMax ) - p->transRevY( p->yMin );
    step = ySpan / yN;
    tmp = pow( 10., floor( log10( step )));
    if(tmp>MT_NUM_EPS)
    {
      i = (floor)(( step / tmp ) + .5 ) ;
      i = MIN(MAX(i,0),ARRAY_LENGTH(roundupTable)-1);
      ySpc = roundupTable[i-1] * tmp; 
    }
    else ySpc = ySpace;
  }
  else
  {
    ySpc = ySpace;
  }

  py1 = p->yMin * p->yScale + p->yOff;
  py2 = p->yMax * p->yScale + p->yOff;

  GWsettxtl(14,1.0,1,GWkrgb(100,100,100),GWkrgb(0,0,0)," ");

  for( mod=0, x = ((floor)( p->xMin / xSpc )+0 ) * xSpc ; 
       x <= p->xMax;
       x += xSpc, mod++ )
  {
       px1 =  x * p->xScale + p->xOff;
       GWline( px1, py1-dpy, px1, py2+dpy );
       if( flags & PL_TEXT_X )
       {
         sprintf( str, "%5.5g", (double)x ); 
         if( mod % 2 == 0 )
         {
           GWputtxtl( px1-8, py1-8, str );
         }  
       }
  }

  px1 = p->xMin * p->xScale + p->xOff;
  px2 = p->xMax * p->xScale + p->xOff;

  for( y = ((floor)( p->transRevY( p->yMin ) / ySpc )+0 ) * ySpc ; 
       y <= p->transRevY( p->yMax );
       y += ySpc )
  {
       py1 = p->transY( y ) * p->yScale + p->yOff;
       GWline( px1+dpx, py1, px2-dpx, py1 );
       if( flags & PL_TEXT_Y )
       {
         sprintf( str, "%5.5g", (double)y ); 
         GWputtxtl( px1, py1, str );
       }
  }
  
  if( flags & PL_FRAME )
  {
    py1 = p->yMin * p->yScale + p->yOff;
    py2 = p->yMax * p->yScale + p->yOff;
    px1 = p->xMin * p->xScale + p->xOff;
    px2 = p->xMax * p->xScale + p->xOff;
    
    GWline( px1, py1, px2, py1 );
    GWline( px2, py1, px2, py2 );
    GWline( px2, py2, px1, py2 );
    GWline( px1, py2, px1, py1 );

  }

return SUCCESS;
}

/********************************************************************/
/*								    */
/********************************************************************/
void PlPlotLineInit( PlLineHd *plHd ) 
{
  plHd->pxOld=plHd->pyOld=0;
  plHd->first = TRUE;
}

/********************************************************************/
/*								    */
/********************************************************************/
int PlPlotLine( PlPlot *p, double x, double y, PlLineHd *plHd ) 
{
  int px,py;

  x = p->transX( x );
  y = p->transY( y );
  
  px = x * p->xScale + p->xOff;
  py = y * p->yScale + p->yOff;
   
  if( plHd->first )
  {
    plHd->first = FALSE; 
  }
  else
  {
    GWline( plHd->pxOld, plHd->pyOld, px, py );    
  }

  plHd->pxOld = px;
  plHd->pyOld = py;
  
  return SUCCESS;
}

/********************************************************************/
/*								    */
/********************************************************************/
int PlPlotPolyLine( PlPlot *p, double *ppx, double *ppy, int n ) 
{
  int px,py,i;
  float *gwPoints=NULL,x,y,*gwp;

  ALLOCARR(gwPoints,n*2);
  gwp = gwPoints;

  for(i=0;i<n;i++)
  {
    x = *ppx++; y = *ppy++;
    x = p->transX( x ); 
    y = p->transY( y );  
    px = x * p->xScale + p->xOff;
    py = y * p->yScale + p->yOff;
    *gwp++ = px; *gwp++=py;
  }

  GWpolylin( gwPoints, n );    

  FREEMEM(gwPoints);

  return SUCCESS;
}

/********************************************************************/
/*								    */
/********************************************************************/
int PlPlotLine2( PlPlot *p, float x1, float y1, 
    			    float x2, float y2 ) 
{
  int px1,py1,px2,py2;

  x1 = p->transX( x1 ); y1 = p->transY( y1 );
  x2 = p->transX( x2 ); y2 = p->transY( y2 );  

  px1 = x1 * p->xScale + p->xOff; py1 = y1 * p->yScale + p->yOff;
  px2 = x2 * p->xScale + p->xOff; py2 = y2 * p->yScale + p->yOff;   
  
  GWline( px1, py1, px2, py2 );    

  
  return SUCCESS;

}

/********************************************************************/
/*								    */
/********************************************************************/
  int PlPlotDot( PlPlot *p, double x, double y ) 
{
  int px,py;

  x = p->transX( x );
  y = p->transY( y );
  
  px = x * p->xScale + p->xOff;
  py = y * p->yScale + p->yOff;

  GWline( px, py, px, py );    

  return SUCCESS;
}

/********************************************************************/
/*								    */
/********************************************************************/
int PlPlotDot2( PlPlot *p, double x, double y, int col ) 
{
  int px,py;

  x = p->transX( x );
  y = p->transY( y );
  
  px = x * p->xScale + p->xOff;
  py = y * p->yScale + p->yOff;

  GWsetpxl( px, py, col );    

  return SUCCESS;
}


/********************************************************************/
/*								    */
/********************************************************************/
int PlPlotPixel( PlPlot *p, double x, double y, int col ) 
{
  int px,py;

  x = p->transX( x );
  y = p->transY( y );
  
  px = x * p->xScale + p->xOff;
  py = y * p->yScale + p->yOff;

  GWsetpxl( px, py, col );    

  return SUCCESS;
}

/********************************************************************/
/*								    */
/********************************************************************/
int PlPlotText( PlPlot *p, double x, double y, int col, char *str )
{
  int px,py;

  x = p->transX( x );
  y = p->transY( y );
  
  px = x * p->xScale + p->xOff;
  py = y * p->yScale + p->yOff;

  if(col>=0)
  {
    GWsettxtl(0,1.0,1,col,GWkrgb(0,0,0)," ");
  }

  GWputtxtl( px, py, str );
   
  return SUCCESS;
}
/********************************************************************/
/*								    */
/********************************************************************/
int PlPlotRect( PlPlot *p, double x, double y, 
    			   double sx, double sy )
{
  int px,py,psx,psy;

  x = p->transX( x );
  y = p->transY( y );
 
  px = x * p->xScale + p->xOff;
  py = y * p->yScale + p->yOff;
  psx=sx*p->xScale;
  psy=sy*p->yScale;

  GWrect( px, py, px+psx-1, py+psy-1 );
   
  return SUCCESS;
}

/********************************************************************/
/*								    */
/********************************************************************/
int PlPlotRect2( PlPlot *p, float x1, float y1, 
    			    float x2, float y2 ) 
{
  int px1,py1,px2,py2;

  x1 = p->transX( x1 ); y1 = p->transY( y1 );
  x2 = p->transX( x2 ); y2 = p->transY( y2 );  

  px1 = x1 * p->xScale + p->xOff; py1 = y1 * p->yScale + p->yOff;
  px2 = x2 * p->xScale + p->xOff; py2 = y2 * p->yScale + p->yOff;   
  
  GWrect( px1, py1, px2, py2 );    

  
  return SUCCESS;

}

/********************************************************************/
/*								    */
/********************************************************************/
int PlPlotCross( PlPlot *p, double x, double y, 
    			   double sx, double sy )
{
  int px,py,psx,psy;

  x = p->transX( x );
  y = p->transY( y );
 
  px = x * p->xScale + p->xOff;
  py = y * p->yScale + p->yOff;
  psx=sx*p->xScale;
  psy=sy*p->yScale;

  GWline( px-sx, py, px+sx, py );    
  GWline( px, py-sy, px, py+sy );    
 
  return SUCCESS;
}

/********************************************************************/
/*								    */
/********************************************************************/
int PlPlotCrossX( PlPlot *p, double x, double y, 
    			   double sx, double sy )
{
  int px,py,psx,psy;

  x = p->transX( x );
  y = p->transY( y );
 
  px = x * p->xScale + p->xOff;
  py = y * p->yScale + p->yOff;
  psx=sx*p->xScale;
  psy=sy*p->yScale;

  GWline( px-sx, py-sy, px+sx, py+sy );    
  GWline( px+sx, py-sy, px-sx, py+sy );    
 
  return SUCCESS;
}

/********************************************************************/
/*								    */
/********************************************************************/
int PlPlotArrow( PlPlot *p, float x1, float y1, 
    			    float x2, float y2,
			    float head_size
			    ) 
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
  PlPlotLine2(p,x1,y1,x2,y2);
  if(dx*dx+dy*dy>0.0000001)
  {	 
    MT_VEC2_ROTATE2(dx2,dy2, dx,dy, arr_ss1, arr_cs1);
    PlPlotLine2(p, x2,y2,x2+dx2*head_size,y2+dy2*head_size);
    MT_VEC2_ROTATE2(dx2,dy2, dx,dy, arr_ss2, arr_cs2);
    PlPlotLine2(p,x2,y2,x2+dx2*head_size,y2+dy2*head_size);	  
  }

  return SUCCESS;
}


/********************************************************************/
/*								    */
/********************************************************************/
int PlPlotFilledRect( PlPlot *p, double x, double y, 
    		      double sx, double sy,
		      int col )
{
  int px,py,psx,psy;

  x = p->transX( x );
  y = p->transY( y );
 
  px = x * p->xScale + p->xOff;
  py = y * p->yScale + p->yOff;
  psx=sx*p->xScale;
  psy=sy*p->yScale;

  GWsrect( px, py, px+psx-1, py+psy-1, col );
   
  return SUCCESS;
}

/********************************************************************/
/*								    */
/********************************************************************/
int PlPlotCircle( PlPlot *p, double x, double y, 
    		  double r )
{
  int px,py,psx,psy;

  x = p->transX( x );
  y = p->transY( y );
 
  px = x * p->xScale + p->xOff;
  py = y * p->yScale + p->yOff;
  psx=r*p->xScale;
  psy=r*p->yScale;

  GWellipse( px-psx, py-psy, px+psx-1, py+psy-1 );
   
  return SUCCESS;
}

/********************************************************************/
/*								    */
/********************************************************************/
int PlPlotDotRect( PlPlot *p, double x, double y, int sx, int sy ) 
{
  int px,py;

  x = p->transX( x );
  y = p->transY( y );
  
  px = x * p->xScale + p->xOff;
  py = y * p->yScale + p->yOff;
   
  GWline( px-sx/2, py-sy/2, px+sx/2, py-sy/2 );    
  GWline( px+sx/2, py-sy/2, px+sx/2, py+sy/2 );    
  GWline( px+sx/2, py+sy/2, px-sx/2, py+sy/2 );    
  GWline( px-sx/2, py+sy/2, px-sx/2, py-sy/2 );    

  return SUCCESS;
}

/********************************************************************/
/*								    */
/********************************************************************/
int PlPlotGraph( PlPlot *p, MtVector vx, MtVector vy, int n,
                 int style ) 
{
  int 	    i;
  PlLineHd  line;
  
  switch( style )
  {
  case PL_LINES:
    PlPlotLineInit( &line );
    for(i=0;i<n;i++)
      PlPlotLine( p, vx[i], vy[i], &line );  
    break;
  case PL_DOTS:  
    for(i=0;i<n;i++)
      PlPlotDot( p, vx[i], vy[i] );  
    break;
  case PL_FATDOTS:  
    for(i=0;i<n;i++)
      PlPlotDotRect( p, vx[i], vy[i], 2, 2 );  
    break;
  default:
    break;
  }      
  return SUCCESS;
}

/********************************************************************/
/*								    */
/********************************************************************/
void PlPlotReverse( PlPlot *p, int px, int py, double *x, double *y ) 
{
  *x = p->transRevX( ( px - p->xOff ) / p->xScale ); 
  *y = p->transRevY( ( py - p->yOff ) / p->yScale ); 
}

/********************************************************************/
/*								    */
/********************************************************************/
void PlColor( int color ) 
{
  GWsetpen(color,1,0,4);
  GWcolor( color, 3 );
  GWcolor( color, 7 );
  GWcolor( color, 1 );
  GWsetpen(color,1,0,4);
}

/********************************************************************/
/*								    */
/********************************************************************/
void PlRefreshScreen( void ) 
{
  //GWrefresh( );
}

/********************************************************************/
/*								    */
/********************************************************************/
void PlClearScreen( void ) 
{
  GWclear(-1 );
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int PlColNib2Byte( int a )
{
  int b = (a+1)*16;
  if(b>255) b=255;
  return b;
}

