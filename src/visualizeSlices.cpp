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
#include <float.h>
#include <limits.h>
#include <math.h>
#include <tbb/tbb.h>
#include <vectorclass/vectorclass.h>
#include "vectorclass_util.h"
#include "globdef.h"
#include "mtool.h"
#include "util.h"
#include "thread.h"
#include "plot.h"
#include "bitmap.h"
#include "panel.h"
#include "msbg.h"

namespace MSBG
{

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::saveVisualizationBitmap( BmpBitmap *B, char *title )
{
  if(!_visSaveAsBitmap) return;
  int rc;

  char path[UT_MAXPATHLEN+1];
  strcpy(path,_visOutDir);
  UT_ASSERT0(*_visOutDir);
  strcat(path,"/");
  char *chp=path+strlen(path);
  for(int i=0;i<(int)strlen(title);i++)
  {
  int ch = title[i];
    if(!isalnum(ch)) ch='_';
    *chp++=ch;
  }
  *chp++=0;

  {
    /*UtTimer tm;
    TIMER_START(&tm);*/
    UtCreatePath(path);
    /*TIMER_STOP(&tm);
    TRCP(("UtCreatePath CPU=%g sec\n",TIMER_DIFF_MS(&tm)/1000.));*/
  }

  char chbuf[100];
  sprintf(chbuf,"/frame_%06d",_totalSteps);
  UtTimer tm;
  TIMER_START(&tm);

  if(!B->data)
  {
    BmpNormalizeChannel(B,BMP_GREY,0);
    BmpGetRGBChannel(B,0);
    #pragma omp parallel for
    for(LongInt i=0;i<B->sx*B->sy;i++) 
    {
      int ig = B->dataGrey[i]*255.0;
      MT_CLAMP(ig,0,255);
      B->data[i] = BMP_MKRGB(ig,ig,ig);
    }
  }

  if(_visSaveAsBitmap & 2 )
  {
    // Lossless
    strcat(chbuf,".png"); strcat(path,chbuf);
    TRC(("Saving output image: '%s'\n",path));
    rc = BmpSaveBitmapPNG( B, path, NULL, BMP_PNG );
  }
  else
  {
    strcat(chbuf,".jpg"); strcat(path,chbuf);
    TRC(("Saving output image: '%s'\n",path));
    // No noticeable JPG artefacts with quality=100
    rc = BmpSaveBitmapJPG( B, path, 100, BMP_RGB );
  }
  
  if(rc) 
  {
    TRCERR(("BmpSaveBitmap('%s') failed %d",path,rc));
  }
  TIMER_STOP(&tm);
  TRC(("%s: CPU %.2f sec, %.0f pixels/sec)\n","BmpSaveBitmap",
      (double)TIMER_DIFF_MS(&tm)/1000.,
      B->sx*(double)B->sy/(double)(TIMER_DIFF_MS(&tm)/1000.0)));

}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::visualizeLastErr( PnlPanel *pnl, int chanId, 
    					   unsigned options )
{
  int x=_lastErrX,
      y=_lastErrY,
      z=_lastErrZ,
      level=_lastErrLevel,
      levelMg=_lastErrLevelMG;

  TRCERR(("%s: '%s' chanId=%d pos=%d,%d,%d, level=%d, levelMg=%d \n",
	UT_FUNCNAME,pnl->title,chanId,x,y,z,level,levelMg));

  Vec3Float pos0 = Vec3Float(x,y,z);
  pos0 = (pos0+0.5f)*(1<<level);
  Vec3Int ipos0(pos0.x, pos0.y, pos0.z);
  _sg0->clipGridCoords(ipos0.x,ipos0.y,ipos0.z);
  double slicePos[100] = 
    { (double)(-ipos0.x), (double)(-ipos0.y), (double)(-ipos0.z) };

  visualizeSlices( pnl, chanId, IP_NEAREST, slicePos, options, levelMg, 0,0,NULL,
    [&]( float f, MtStat *fstat,
	int flags, int blockFlags, int level, const Vec3Int& ipos,
	int &colorOut ) -> void
    {
      if(ipos.x==ipos0.x && ipos.y==ipos0.y && ipos.z==ipos0.z )
      {
        Vec3Int c(BMP_RED(colorOut), BMP_GREEN(colorOut), BMP_BLUE(colorOut));
	for(int k=0;k<3;k++)
	{
	  c[k] = MAX(MIN((int)(c[k]*3),255),0);
	}
	//colorOut = BMP_MKRGB(255,255,255); //BMP_MKRGB(c.x,c.y,c.z);
	colorOut = BMP_MKRGB(c.x,c.y,c.z);
      }           
    });		  
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void MultiresSparseGrid::visualizeSlices( 
    	PnlPanel *pnl, 
	int chanId,
	int ipolType,
	double *slicePos0,   // slice-pos(x,y,z) [ roi(x0,y0,z0,sx,sy,sz,zoom) ] 
	unsigned opt,
	int levelMg,
	double thresholdCol,
	int visCellFlagsMask,
	SparseGrid<float> *sgVis0,
	std::function<void(
	  float val, MtStat *valStat,
	  int flags, int blockFlags, int level,
	  const Vec3Int& ipos, 
	  int &colorInOut /* in/out */
	  )> setCellColorFunc
    ) 
{
//#define VERIFY_IPOS
  if(sgVis0) chanId = CH_FLOAT_2;
  if( opt & VIS_OPT_COLOR_RGB ) opt |= VIS_OPT_NONRM;
  int isVecChan = isVecChannel(chanId) || (opt & VIS_OPT_AIR_VELOCITY),
  	doVisReslevel=FALSE,
	doVisGradient = opt & VIS_OPT_GRADIENT;
  int nChan = isVecChan || doVisGradient ? 3 : 1;
  double *slicePos = slicePos0 ? slicePos0 : _visSlicePos;
  int dBorder=0;
  int zf = 1,
      doZoom = FALSE;

  if(ipolType!=IP_NEAREST && chanId!=CH_DIST_FINECOARSE)
  {
    if(!(_validForInterpolation & (((int64_t)1)<<chanId)))
    {
      prepare(chanId);
    }
  }

  if(  _visConf.tm_mod  && (_totalSteps % _visConf.tm_mod) != 0 ) 
  {
    TRC(("%s: _totalSteps=%d _visConf.tm_mod=%d -> skip.\n",
	  UT_FUNCNAME,_totalSteps,_visConf.tm_mod));
    return;
  }

#if 0
  opt |= VIS_OPT_DOMAIN_BC;

  dBorder = (opt & VIS_OPT_DOMAIN_BC) && 
    		!(opt & VIS_OPT_SUM_SLICES) ? 1 : 0;
#endif

  UtTimer tm;

  TIMER_START(&tm);

  double slicePosMax[20];
  if( opt & VIS_OPT_ABSMAX_POS )
  {
    UT_ASSERT0(chanId);
    UT_ASSERT0(!sgVis0);
    Vec3Int iposMax;
    int levelMax;
    double valMax;
    valMax = maxAbsChannel(chanId,-1,levelMg,0,
			    &iposMax,&levelMax);
    TRC(("%s chan=%d. valMax= %g posMax=%d:%d,%d,%d\n",UT_FUNCNAME,
	  valMax,levelMax,iposMax.x,iposMax.y,iposMax.z));
    Vec3Float pos0 = Vec3Float(iposMax.x,iposMax.y,iposMax.z);
    pos0 = (pos0+0.5f)*(1<<levelMax);
    Vec3Int ipos0(pos0.x, pos0.y, pos0.z);
    _sg0->clipGridCoords(ipos0.x,ipos0.y,ipos0.z);
    for(int k=0;k<3;k++) slicePosMax[k] = -ipos0[k];
    slicePos = slicePosMax;
  }

  SparseGrid<Vec3Float> *sg0 = _sparseGrids[0].vec3_1[0];

  if(!(slicePos[0]<1)) slicePos[0] *= 1.0/sg0->sx();
  if(!(slicePos[1]<1)) slicePos[1] *= 1.0/sg0->sy();
  if(!(slicePos[2]<1)) slicePos[2] *= 1.0/sg0->sz();

  Vec3Float roiPos(0,0,0),roiSize(1,1,1);
  //double roiZoom=1;
  /*if(opt&VIS_OPT_USE_ROI_)
  {
    dBorder = 0;
    for(int k=0;k<3;k++) 
    {
      roiPos[k] = slicePos0[3+k]; 
      roiSize[k] = slicePos0[3+3+k];
    }
    roiZoom = slicePos0[9];
    TRC(("using roi: %g,%g,%g %gx%gx%g zoom=%g\n",
	  roiPos.x,roiPos.y,roiPos.z,
	  roiSize.x,roiSize.y,roiSize.z,
	  roiZoom ));
  }*/

  if( opt & VIS_OPT_ZOOM )
  {
    doZoom = TRUE;
    zf = 6;
  }

  if(opt & VIS_OPT_RESLEV)
  {
    doVisReslevel = TRUE;
  }

  if(levelMg>0)
  {
    UT_ASSERT0(ipolType==IP_NEAREST);
  }

  BmpBitmap *BXZ=NULL,*BXY=NULL,*BZY=NULL;
  BXZ = BmpNewBitmap((sg0->sx()+2*dBorder)*zf, 
                     (sg0->sz()+2*dBorder)*zf, 
      		     BMP_XYZ|BMP_GREY|BMP_CLEAR);
  BXY = BmpNewBitmap((sg0->sx()+2*dBorder)*zf, 
      		     (sg0->sy()+2*dBorder)*zf, 
      		     BMP_XYZ|BMP_GREY|BMP_CLEAR);
  BZY = BmpNewBitmap((sg0->sz()+2*dBorder)*zf, 
                     (sg0->sy()+2*dBorder)*zf, 
      		     BMP_XYZ|BMP_GREY|BMP_CLEAR);

  for(int k=0;k<2;k++)
  {
    BmpGetUshortChannel(BXZ,k,BMP_CLEAR);
    BmpGetUshortChannel(BXY,k,BMP_CLEAR);
    BmpGetUshortChannel(BZY,k,BMP_CLEAR);
  }
#ifdef VERIFY_IPOS
  for(int k=2;k<5;k++)
  {
    BmpGetUshortChannel(BXZ,k,BMP_CLEAR);
    BmpGetUshortChannel(BXY,k,BMP_CLEAR);
    BmpGetUshortChannel(BZY,k,BMP_CLEAR);
  }
#endif

  //
  // Get slices
  //
  int 
    ix0 = slicePos[0] < 0 ? round(-slicePos[0]) : 
		(int)(sg0->sx()*(double)slicePos[0])*zf,
    iy0 = slicePos[1] < 0 ? round(-slicePos[1]) : 
    		(int)(sg0->sy()*(double)slicePos[1])*zf,
    iz0 = slicePos[2] < 0 ? round(-slicePos[2]) : 
		(int)(sg0->sz()*(double)slicePos[2])*zf;

  TRC(("%s iSlice=%d,%d,%d\n",UT_FUNCNAME,ix0,iy0,iz0));

  BmpRectangle roiXZ={0},roiXY={0},roiZY={0};

  if( !sg0->inRange(ix0,iy0,iz0) )
  {
    TRCERR(("visualizeSlices: slice coords %d,%d,%d, out of grid range (%dx%dx%d)\n",
	  ix0,iy0,iz0,sg0->sx(),sg0->sy(),sg0->sz()));
    goto cleanup;
  }

  {
    BmpBitmap *B = BZY;
    BMP_SET_RECT( &roiZY, 
		  B->sx*roiPos[0], B->sy*roiPos[1],
		  round(B->sx*roiSize[0]), round(B->sy*roiSize[1]));	  	
  }
  {
    BmpBitmap *B = BXZ;
    BMP_SET_RECT( &roiXZ, 
		  B->sx*roiPos[0], B->sy*roiPos[1],
		  round(B->sx*roiSize[0]), round(B->sy*roiSize[1]));	  	
  }
  {
    BmpBitmap *B = BXY;
    BMP_SET_RECT( &roiXY, 
		  B->sx*roiPos[0], B->sy*roiPos[1],
		  round(B->sx*roiSize[0]), round(B->sy*roiSize[1]));	  	
  }


#if 0
  {
    // Visualize 1D slice
    PlPlot plp; PlWindow *plw=&pnl->plwin[0];

    float tmin = ( _visConf.roi_box_pos[2] - 0.5*_visConf.roi_box_size ) * sg0->sz(),
	  tmax = ( _visConf.roi_box_pos[2] + 0.5*_visConf.roi_box_size ) * sg0->sz();

    float vmin=-.125,vmax=1.25;

    PlPlotInit( &plp, plw, tmin, tmax, vmin, vmax , 0 );
    PnlSelect( pnl, 0 );
    {
      PnlErase( pnl );
      PnlClear( pnl );
    }
    GWsetpen( 0, 1, 0, 4);
    PlColor(GWkrgb(70,70,70)); 
    PlPlotAxes(&plp,PL_AXDEFAULT|PL_FRAME,10,10);
    PlColor(GWkrgb(200,200,200)); 

    PlLineHd plhd;
    PlPlotLineInit(&plhd);
    PlColor(PLRGB(0xFF6));
    for(int i=0;i<tmax;i++)
    {
      float t=i+0.5f;
      Vec4f pos(ix0+0.5f,iy0+0.5f,t,0.0f);
      float val =  interpolateScalarFast2( chanId, pos, 									
	  				   OPT_FAST_DISCONTINUOUS_FINECOARSE );
      #if 0
      PlPlotDot2(&plp,t,val,GWkrgb(200,200,200));
      #else
      PlPlotLine(&plp,t,val,&plhd);
      #endif
    }

    goto cleanup;
  }
#endif

  if((opt&VIS_OPT_SUM_SLICES))
  {
    getSumSlices(chanId,
		 BZY,BXZ,BXY);
  }
  else
  {

  #pragma omp parallel for schedule( dynamic )
  for(int iz=-dBorder*zf; iz<(sg0->sz()+dBorder)*zf; iz++) 
  for(int iy=-dBorder*zf; iy<(sg0->sy()+dBorder)*zf; iy++)
  for(int ix=-dBorder*zf; ix<(sg0->sx()+dBorder)*zf; ix++)  
  {
    int yIsSlice = (iy==iy0),
	zIsSlice = (iz==iz0),
	xIsSlice = (ix==ix0);

    if (xIsSlice || yIsSlice || zIsSlice )
    {
      CellFlags cellFlags=0;
      uint16_t blockFlags=0;
      Vec3Float data;
      int level,tlevel=0;

      /*if(xIsSlice && iz==614 && iy==149)
      {
	TRCERR(("DBG\n"));
      }*/

      if( doZoom )
      {
	UT_ASSERT0( !(opt & VIS_OPT_AIR_VELOCITY ));
	UT_ASSERT0(ipolType != IP_NEAREST);
	UT_ASSERT0(!sgVis0 );

	int ipolOpt = OPT_IPCHKBLK|
		      OPT_IPEXTRAPOL|
		      OPT_IPBC_NEUMAN;

	Vec4f pos(ix,iy,iz,0.0f);
	pos = (pos+0.5f)*1.0f/zf;
        CellAccessor ca = getCellAccessor( pos );
	SparseGrid<CellFlags> *sgFlags = getFlagsChannel(CH_CELL_FLAGS,ca.level,levelMg);
	cellFlags = sgFlags ? sgFlags->getValueInBlock0(ca.bid,ca.vid) : 0;	
	BlockInfo *bi = getBlockInfo(ca.bid,levelMg);
	UT_ASSERT2(bi->level==ca.level);
	blockFlags = bi->flags;
	level = bi->level;
	tlevel = _blockTimeLevels[ca.bid];

	if( doVisGradient )
	{
	  Vec4f grad;
	  UT_ASSERT0(FALSE);
	  interpolateWithDerivs(IP_BSPLINE_QUAD,chanId,pos,0,0,NULL,&grad);
	  //interpolateWithDerivs(IP_BSPLINE_CUBIC,chanId,pos,0,NULL,&grad);
	  //Vec4f grad = getScalarFieldDerivs(GRAD_QUAD_SPLINE,chanId,pos,0);
	  data = Vec3Float(grad);
	}
	else
	{
	  float pos_[3] = {vfget_x(pos),vfget_y(pos),vfget_z(pos)};
	  interpolate( ipolType, chanId, pos_, ipolOpt, data ); 
	}
      }
      else
      {
	UT_ASSERT0(!doVisGradient);

      int ixc=ix,
	  iyc=iy,
	  izc=iz;
      int clipped=0, clipDir=-1, clipFace=-1;

      sg0->clipGridCoords(ixc,iyc,izc,
	  		  clipped, clipDir, clipFace );
      int bx,by,bz;
      sg0->getBlockCoords(ixc,iyc,izc, bx,by,bz);
      UT_ASSERT0(sg0->blockCoordsInRange(bx,by,bz));
      LongInt bid = sg0->getBlockIndex(bx,by,bz);
      UT_ASSERT2(bid>=0&&bid<sg0->nBlocks());

      int levelBlk,
	  isGhostBlock=FALSE;

      if(sgVis0)
      {
        levelBlk = sgVis0->getBlock(bid) ? 0 :_nLevels-1; 
      }
      else
      {
	levelBlk = getBlockLevel( bid, levelMg );
      }

      if(levelMg < _nLevels)
      {
	if( opt & VIS_OPT_GHOST_FINE )
	{
	  BlockInfo *bi = getBlockInfo( bid, levelMg );
	  if( bi->flags & BLK_COARSE_FINE )
	  {
	    UT_ASSERT2(bi->level>0);
	    UT_ASSERT2(levelBlk==bi->level);
	    isGhostBlock = TRUE;
	  }
	}
      }

      level = levelBlk;
      tlevel = _blockTimeLevels[bid];

      if(levelMg>0 && !isGhostBlock)
      {
	UT_ASSERT2(!sgVis0);
	level = MAX(level,levelMg);
      }

      if(isGhostBlock)
      {
	UT_ASSERT2(levelBlk>0&&level>0);
	levelBlk -= 1;
	level -= 1;
      }
      
      //UT_ASSERT0(!(level<0||level>_nLevels-1));	      
//      Level *lev = &_sparseGrids[level];
      int ix1=ixc>>level,
	    iy1=iyc>>level,
	    iz1=izc>>level;

      SparseGrid<CellFlags> 
	*sgFlags = getFlagsChannel(CH_CELL_FLAGS,level,levelMg);

//      CellFlags  cellFlags = lev->cellFlags[0]->getValue(ix1,iy1,iz1);
      
      cellFlags = sgFlags && sgFlags->hasData() ? sgFlags->getValueGen_(ix1,iy1,iz1) : 0;
      blockFlags = levelMg < _nLevels ? getBlockInfo(bid,levelMg)->flags : 0;

      int ipolOpt = OPT_IPCHKBLK|
		    OPT_IPEXTRAPOL|
		    OPT_IPBC_NEUMAN;

      float pos[4] = { (float)ixc+0.5f, (float)iyc+0.5f, (float)izc+0.5f, 0.0f };

      //if(xIsSlice||yIsSlice||zIsSlice)
      if(sgVis0)
      {
	float f = ipolType==IP_NEAREST ? sgVis0->getValue(ixc,iyc,izc) :
	  	  ipolType==IP_LINEAR ? sgVis0->interpolate<IP_LINEAR>(pos,ipolOpt) :
		  0;
	data[0]=data[1]=data[2] = f;
      }
      else
      {
	if(isUint8Channel(chanId))
	{
	  SparseGrid<uint8_t> *sg = getUint8Channel(chanId,level,levelMg);
	  uint8_t val = sg && sg->hasData() ? sg->getValueGen_(ix1,iy1,iz1) : 0;
	  data[0]=data[1]=data[2] = (opt & VIS_OPT_RENDERDENS) ? renderDensToFloat_(val) : val;
	}
	else if(isUint16Channel(chanId))
	{
	  SparseGrid<uint16_t> *sg = getUint16Channel(chanId,level,levelMg);
	  uint16_t val = sg && sg->hasData() ? sg->getValueGen_(ix1,iy1,iz1) : 0;
	  data[0]=data[1]=data[2] = (opt & VIS_OPT_RENDERDENS) ? renderDensToFloat_(val) : val;
	}
	else
	{
	  if(opt & VIS_OPT_AMBIENTOCCLUSION)
	  {
	#if 0
	    data[0] = data[1] = data[2] = 
	      getAmbientOcclusion( chanId, Vec4f(pos[0],pos[1],pos[2],0.0f),1,0,1,0);
	#else
	    UT_ASSERT0(FALSE);

	#endif
	  }
	  else
	  {

	  switch(ipolType)
	  {
	    case IP_NEAREST:
	      if(chanId==CH_FACE_DENSITY)
	      {
		FaceDensity fd = getChannel<FaceDensity>(this,chanId,level,levelMg)->
		    getValueGen_(ix1,iy1,iz1);
		for(int k=0;k<3;k++) data[k] = FACE_DENSITY_TO_FLOAT( fd[k] );
	      }
	      else if(chanId==CH_FACE_AREA||chanId==CH_FACE_COEFF)
	      {
		int iDir = opt & VIS_OPT_DIR_X ? 0 :
			   opt & VIS_OPT_DIR_Y ? 1 :
			   opt & VIS_OPT_DIR_Z ? 2 :
			   -1;
		UT_ASSERT2(iDir>=0);
		SparseGrid<float> *sg = chanId==CH_FACE_AREA ?
		  getFaceAreaChannel(iDir,level,levelMg) :
		  getFaceCoeffChannel(iDir,level,levelMg);
		data[0] = sg ? sg->getValueGen_(ix1,iy1,iz1) : 0.0f;
	      }
	      else
	      {
		if(levelMg>0||isGhostBlock)
		{
		  if(isVecChan)
		  {
		    data = getVecChannel(chanId,level,levelMg)->getValueGen_(ix1,iy1,iz1);
		  }
		  else
		  {
		    SparseGrid<float> *sg = getFloatChannel(chanId,level,levelMg);
		    data[0] = sg ? sg->getValueGen_(ix1,iy1,iz1) : 0.0f;
		  }
		}
		else
		{
		  interpolate<IP_NEAREST>( chanId, pos, 0, data, 
		      			    opt & VIS_OPT_LEVEL0_ONLY ); 
		}
	      }
	      break;
	    case IP_LINEAR: interpolate<IP_LINEAR>( chanId, pos, ipolOpt, data ); break;
	    case IP_CUBIC_MONO_2: interpolate<IP_CUBIC_MONO_2>( chanId, pos, ipolOpt, data ); break;
	    case IP_CUBIC: interpolate<IP_CUBIC>( chanId, pos, ipolOpt, data ); break;
	    case IP_WENO4: interpolate<IP_WENO4>( chanId, pos, ipolOpt, data ); break;
	    default:
	      UT_ASSERT0(FALSE);
	      break;
	  }


#if 0
	  if(chanId==CH_FLOAT_8)
	  {
	    float f = data[0];
	#if 1
	    MT_CLAMP(f,-10,10);	    

//	    f = 5*floorf(f/5);
//	    f = 10*floorf(f/10);
	    f = levelMg>0 ? f : 1*floorf(f/1);
	#endif

	    data[0] = data[1] = data[2] = f;
	  }
	  else 
#endif
	  }
	}
      }

      }

      if( opt & VIS_OPT_DISTFIELD )
      {
	float nb = thresholdCol;
	USE(nb);
	float f = data[0];
	#if 1
	f = std::min(std::max(f,-nb),nb);
	#else
	f = f<=-SDF_INFINITY ? -100 :
	  f<SDF_INFINITY ? f:
	  100;
	#endif
	data[0] = f;
      }

      int levelAct = opt & VIS_OPT_TIMELEVEL ? tlevel :
		     level;

      if(yIsSlice)
      {
	BmpRectangle *roi=&roiXZ;
	int ix2=ix-roi->x0+dBorder*zf,
	    iz2=iz-roi->y0+dBorder*zf;	
	if(BMP_IN_RANGE(ix2,iz2, 0,0,roi->sx-1,roi->sy-1))  
	{
	  LongInt i=BMP_XYOFF(BXZ,ix+dBorder*zf,iz+dBorder*zf);
	  for(int k=0;k<nChan;k++)  BXZ->dataFloat[k][i] = data[k];  
	  BXZ->dataGrey[i] = levelAct; 
	  BXZ->dataUshort[0][i] = cellFlags; 
	  BXZ->dataUshort[1][i] = blockFlags;
#ifdef VERIFY_IPOS
	  BXZ->dataUshort[2][i] = ix; 
	  BXZ->dataUshort[3][i] = iy;
	  BXZ->dataUshort[4][i] = iz; 
#endif
	}
      }
      if(zIsSlice)
      {
	BmpRectangle *roi=&roiXY;
	int ix2=ix-roi->x0+dBorder*zf,
	    iy2=iy-roi->y0+dBorder*zf;
	if(BMP_IN_RANGE(ix2,iy2, 0,0,roi->sx-1,roi->sy-1))
	{
	  LongInt i=BMP_XYOFF(BXY,ix+dBorder*zf,iy+dBorder*zf);
	  for(int k=0;k<nChan;k++) BXY->dataFloat[k][i] = data[k];
	  BXY->dataGrey[i] = levelAct;
	  BXY->dataUshort[0][i] = cellFlags; 
	  BXY->dataUshort[1][i] = blockFlags;
#ifdef VERIFY_IPOS
	  BXY->dataUshort[2][i] = ix; 
	  BXY->dataUshort[3][i] = iy;
	  BXY->dataUshort[4][i] = iz; 
#endif
	}
      }
      if(xIsSlice)
      {	
	BmpRectangle *roi=&roiZY;
	int iz2=iz-roi->x0+dBorder*zf,
	    iy2=iy-roi->y0+dBorder*zf;
	if(BMP_IN_RANGE(iz2,iy2, 0,0,roi->sx-1,roi->sy-1))
	{
	  LongInt i=BMP_XYOFF(BZY,iz+dBorder*zf,iy+dBorder*zf);
	  for(int k=0;k<nChan;k++) 
	    BZY->dataFloat[k][i] = data[k];
	  BZY->dataGrey[i] = levelAct;
	  BZY->dataUshort[0][i] = cellFlags;
	  BZY->dataUshort[1][i] = blockFlags;
#ifdef VERIFY_IPOS
	  BZY->dataUshort[2][i] = ix; 
	  BZY->dataUshort[3][i] = iy;
	  BZY->dataUshort[4][i] = iz; 
#endif
	}
      }    
    }
  }

  }


  // 
  // Post process slices
  //
  MtStat gstat_[4*4];
  // Color according to refinement level
  for(int k=0;k<3;k++)  // Slices
  {
    BmpBitmap *bitmaps[3] = { BZY, BXZ, BXY };   
    BmpRectangle *rois[3] = { &roiZY, &roiXZ, &roiXY },
		 *roi = rois[k];
    BmpBitmap *B=bitmaps[k];
    for(int j=0;j<nChan;j++)  // grey channels
    {
      MtStat *gstat = gstat_+j*4+1+k; 
      int igchan = j;
					
      MtGetStat2_f(B->sx*B->sy, B->dataFloat[igchan],
	  -PCELL_INF_QUANTITY, PCELL_INF_QUANTITY, gstat);
      if(gstat->span<MT_NUM_EPS) gstat->span=1;

      #pragma omp parallel for
      for(int iby=0;iby<B->sy;iby++)
      for(int ibx=0;ibx<B->sx;ibx++)
      {
	LongInt i = BMP_XYOFF(B,ibx,iby);
	int level = B->dataGrey[i];	
	unsigned cellFlags = B->dataUshort[0][i],
		 blockFlags = B->dataUshort[1][i];
	int obst = (cellFlags & CELL_SOLID);
	USE(cellFlags);
	USE(blockFlags);

	Vec3Int igpos;
	switch(k)
	{
	  case 0: // ZY
	    igpos.x = ix0;
	    igpos.y = iby + roi->y0 - dBorder*zf;
	    igpos.z = ibx + roi->x0 - dBorder*zf;
	    break;
	  case 1: //XZ
	    igpos.x = ibx + roi->x0 - dBorder*zf;
	    igpos.y = iy0;
	    igpos.z = iby + roi->y0 - dBorder*zf;
	    break;
	  case 2: //XY
	    igpos.x = ibx + roi->x0 - dBorder*zf;
	    igpos.y = iby + roi->y0 - dBorder*zf;
	    igpos.z = iz0; 
	    break;	    
	  default:
	    UT_ASSERT0(FALSE);
	}

	#ifdef VERIFY_IPOS
	{
	  int xRef=B->dataUshort[2][i],
	      yRef=B->dataUshort[3][i],
	      zRef=B->dataUshort[4][i];
	  UT_ASSERT0( igpos.x==xRef &&
	      	      igpos.y==yRef &&
		      igpos.z==zRef );
	}
	#endif

	//UT_ASSERT0(level>=0&&level<_nLevels);

	int cr=0,cg=0,cb=0;
	if(level==0)
	{
	  cr=255; cg=200; cb=50;
	  if(obst) { cr=220; cg=220; cb=220; }
	}
	else if(level==1)
	{
	  cr=80; cg=255; cb=150; 
	  if(obst) { cr=150; cg=150; cb=150; }
	}
	else if(level==2)
	{
	  cr=170; cg=170; cb=255;
//	  cr=100; cg=50; cb=255;
	  if(obst) { cr=80; cg=80; cb=80; }
	}
	double g = B->dataFloat[igchan][i];
	if(!(std::isfinite(g)))
	{
	  TRCERR(("Invalid floating point number %g\n",g));
	  TRCERR(("pnl='%s' chan=%d ipol=%d opt=%d level=%d levelMg=%d," 
		  "pos(L0)=%d,%d,%d, flags=%d, bflags=%d\n",
		pnl ? pnl->title : "-",
		chanId,ipolType,opt,level,levelMg,igpos.x,igpos.y,igpos.z,
		cellFlags,blockFlags));
	  g = 0;
	}

	double threshColAct = thresholdCol;
	float g0 = g;

	if(!( g0<PCELL_INF_QUANTITY )) g=0;

	if(!(opt & VIS_OPT_NONRM) && g0<PCELL_INF_QUANTITY)
	{
	  g = (g-gstat->min)/gstat->span;	  
	  threshColAct = (thresholdCol-gstat->min)/gstat->span;
	  UT_ASSERT_NR(!(g<-0.01||g>1.01));
	  MT_CLAMP(g,0,1);
	  //if(j==1) g = sqrtf(g);
	  //g = 0.2+0.8*g;
	}

	MT_CLAMP(g,0,1);

	int color=0;

	if(opt&VIS_OPT_THRESHCOL && g0<PCELL_INF_QUANTITY)
	{
	  Vec3Float col(g,g,g);
	  if(g>threshColAct) 
	  {
	    col.y *= 0.85;
	  }
	  else
	  {
	    col.x *= 0.85;
	    col.z *= 0.85;
	  }
	  color = BMP_MKRGB(255*col.x,255*col.y,255*col.z); 
	}
	else if(opt & VIS_OPT_TEST)
	{
#if 0
	  // 
	  // Color level
	  //
	  color = BMP_MKRGB(  (int)(g*cr), (int)(g*cg), (int)(g*cb) );
#elif 0
#elif 1
	  //
	  // Black/White
	  //
	  if(chanId==CH_CELL_FLAGS)
	  {
	   int isPartialSolid = (cellFlags&CELL_PARTIAL_SOLID);
	   color = 
	     CELL_IS_LIQUID(cellFlags) && !isPartialSolid ? BMP_MKRGB( 0,0,255 ) :
	     CELL_IS_LIQUID(cellFlags) && isPartialSolid ? BMP_MKRGB( 128,0,255 ) :

	     (cellFlags & CELL_AIR) && !isPartialSolid ? BMP_MKRGB( 200,200,0 ) :
	     (cellFlags & CELL_AIR) && isPartialSolid ? BMP_MKRGB( 100,100,0 ):

	     (cellFlags & CELL_VOID) && !isPartialSolid ? BMP_MKRGB( 0,255,0 ) :
	     (cellFlags & CELL_VOID) && isPartialSolid ? BMP_MKRGB( 100,210,50 ):
	     cellFlags & CELL_SOLID ? BMP_MKRGB( 255,0,0 ) :
	     BMP_MKRGB(0,0,0);

	   if( cellFlags & CELL_BOUNDARY_ZONE )
	   {
	     float rgb[3] = { BMP_RED(color)*0.75f, 
	       		      BMP_GREEN(color)*0.75f, 
	       		      BMP_BLUE(color)*0.75f};
	     color = BMP_MKRGB(rgb[0],rgb[1],rgb[2]);	     
	   }

	   if(cellFlags & CELL_OUT_OF_DOMAIN)
	   {
	     float rgb[3] = { BMP_RED(color)*0.25f, 
	       		      BMP_GREEN(color)*0.25f, 
	       		      BMP_BLUE(color)*0.25f};
	     color = BMP_MKRGB(rgb[0],rgb[1],rgb[2]);	     
	   }

	    if(visCellFlagsMask)
	    {
	#if 0
	      color = 
		BMP_MKRGB( 255*!!(cellFlags&CELL_SOLID),
			   255*!!(cellFlags&visCellFlagsMask),
			   opt & VIS_OPT_BLKBORDER ? (255*!!(cellFlags&CELL_BLK_BORDER )): 128
			  );	      
	#else
	      color = cellFlags & visCellFlagsMask ? BMP_MKRGB(255,255,255) : BMP_MKRGB(0,0,0);
	#endif
	    }
	  }
	  else 
	  {
	    color = visGetCellColor( g, g0, cellFlags, blockFlags, 
	      doVisReslevel ? level : -1, level,
	      opt  );
	  }


	  if( cellFlags & CELL_FIXED )
	  {
	     float rgb[3] = { (float)BMP_RED(color), 
	       		      (float)BMP_GREEN(color), 
	       		      (float)BMP_BLUE(color)};

	     for(int k=0;k<3;k++) rgb[k] = 0.5*rgb[k] + 0.5*127;
	     color = BMP_MKRGB(rgb[0],rgb[1],rgb[2]);	     
	  }
#endif
     
	}
	else if(opt & VIS_OPT_BLOCKS)
	{
	  double b = level<_nLevels ? 
	    100+(_nLevels-1-level)*(255-100)/_nLevels : 255;

	    #if 0
	    color = BMP_MKRGB( b*!!(cellFlags&CELL_SOLID),
			       b*!!(blockFlags&BLK_BOUNDARY_ZONE),
			       b*!!(cellFlags&CELL_BOUNDARY_ZONE));
	    #endif

	    #if 0
	    color = BMP_MKRGB( b*!!(cellFlags&CELL_SOLID),
			       b*!!(blockFlags&BLK_NO_OBSTACLES),
			       b*!!(blockFlags&BLK_ONLY_OBSTACLES));					   
	    #endif

	    #if 1
	    color = BMP_MKRGB( b*!!(cellFlags&CELL_SOLID),
			       b*!!(blockFlags&BLK_ONLY_FLUID),
			       b*!!(BLK_IS_RES_BORDER(blockFlags)));					   
	    #endif


/*	  uint32_t color2 = (blockFlags & (BLK_RES_BORDER)) == 0 ?
		  BMP_MKRGB(255*g,255*g,255*g) :
		  BMP_MKRGB( (int)((blockFlags&BLK_RES_BORDER)?240:0)*g,
			     (int)((level==1)?240:0)*g,
			     (int)((level==2)?240:0)*g);

	  color ^= color2;*/

	}
	else
	{
	 // MT_CLAMP(g,0,1);
	 // color = BMP_MKRGB((int)(g*255),(int)(g*255),(int)(g*255));

	  double g = B->dataFloat[igchan][i];

	  if( !(g0<PCELL_INF_QUANTITY)  /*|| 
	      ( cellFlags & cellEmptyValFlag )*/) 			     			    			     
	  {
	    color = VIS_COL_INVALID_QUANTITY;
	  }
	  else
	  {

	    if((opt & VIS_OPT_NONRM) && thresholdCol>1e-10)
	    {
	      g = g/thresholdCol;
	      MT_CLAMP(g,0,1);
	    }
	    else
	    {
	      g = (g-gstat->min)/gstat->span;
	    }

	    if(opt&VIS_OPT_DISPLAY_GAMMA)
	    {
	      g = sqrtf(g);
	    }

	    if(opt & VIS_OPT_DISTFIELD)
	    {
	      color =
		g0<0 ? BMP_MKRGB( MAX((g)*255,10), MAX((g)*255,10), 0 ) : 
		      	    BMP_MKRGB( 0, MAX((1-g)*255,10), MAX((1-g)*255,10) );
	    }
	    else
	    {
	      Vec3Int c(g*255, g*255, g*255);
	      if( opt & VIS_OPT_COLOR_PHASES )
	      {
		if( cellFlags & CELL_VOID )
		{
		  c = 0;
		}
		else if( cellFlags & CELL_SOLID )
		{
		  c.x = 112;
		  c.y = 70;
		  c.z = 42;
		}
		else if( CELL_IS_AIR(cellFlags))
		{
		  c.x = MIN((c.x+70),255);
		  c.y = MIN((c.y+70),255);
		  c.z = 0;
		}
		else if( CELL_IS_LIQUID(cellFlags))
		{
		  c.x=0;
		  c.y = MIN((c.y+70),255);
		  c.z = MIN((c.z+70),255);
		}
	      }
	      color = /*CELL_IS_AIR( cellFlags ) ? 
		BMP_MKRGB( 0, g*255, g*255 ) :*/
		BMP_MKRGB( c.x, c.y, c.z );
	    }
#if 0
	    UT_ASSERT_NR(!(g<-0.01||g>1.01));
	    MT_CLAMP(g,0,1);

	    float vHCL[3] = { (float)g,   // brightness [0,1]
			      0.5f, 	// Hue [-1,1]
			      -1 },	// Saturation [-1,1]
		  vRGB[3];
	#if 0
	    float hue = ( cellFlags & CELL_SOLID ) ? 0 : 
			( cellFlags & CELL_DOM_BORDER ) ? 0.25 :	    		
						 0.75;
	    vHCL[1] = hue;
	#endif

	    CoConvHcl2Rgb( vRGB, vHCL );
	    color = BMP_MKRGB( vRGB[0]*255, vRGB[1]*255, vRGB[2]*255 );
#endif
	  }
	}


	/*if(opt&VIS_OPT_DISTFIELD)
	{
	  int ig = BMP_RED(color);
	  color = SDF_SIGN(g0)<0 ? BMP_MKRGB(ig,0,0) : BMP_MKRGB(0,ig,0);
	}*/

	if( setCellColorFunc )
	{
	  setCellColorFunc( g0, gstat, cellFlags, blockFlags, level,igpos, 
			    color );
	}

	if(doVisReslevel)
	{
	  int c[3] = { BMP_RED(color), BMP_GREEN(color), BMP_BLUE(color) };

	  if(level>=2) c[2] = MIN(c[2]+64,255);
	  else if(level>=1) c[1] = MIN(c[1]+64,255);	  
	  color = BMP_MKRGB(c[0],c[1],c[2]);	     
	}

	((uint32_t*)B->dataFloat[igchan])[i] = color;
            
      }
    }
  }

  /*for(int i=0;i<BXZ->sx*BXZ->sy;i++)
    ((uint32_t*)BXZ->dataFloat[0])[i] = BMP_MKRGB(200,50,250);*/

  {
    BmpBitmap *B=NULL;

    pnl->totalTime = _totalTime; 

    PnlVisualize3DFieldSlices( pnl, NULL, 
			   BZY->dataFloat, BXZ->dataFloat, BXY->dataFloat,
			   nChan, BXY->sx, BXY->sy, BZY->sx, 
			   slicePos[0],slicePos[1],slicePos[2], gstat_, 
			   PNL_RGB | PNL_PROVIDE_SLABSTATS, 
			   &B 
			   );

    pnl->totalTime = -1;

    saveVisualizationBitmap(B,pnl->title);
    
    BmpDeleteBitmap(&B);
  }

  TIMER_STOP(&tm);
  TRC(("%s: CPU %.2f sec, %.0f voxels/sec)\n",UT_FUNCNAME,
      (double)TIMER_DIFF_MS(&tm)/1000.,
      _nActCells[0]/(double)(TIMER_DIFF_MS(&tm)/1000.0)));
cleanup:
  BmpDeleteBitmap(&BXZ);
  BmpDeleteBitmap(&BZY);
  BmpDeleteBitmap(&BXY);
}

/*-------------------------------------------------------------------------------------*/
/* 										       */
/*-------------------------------------------------------------------------------------*/
void MultiresSparseGrid::visualizeBlockFlags( 
    	PnlPanel *pnl,
	unsigned flagRed,
	unsigned flagGreen,
	unsigned flagBlue,
	int levelMg )
{
  visualizeSlices( pnl, CH_CELL_FLAGS,IP_NEAREST,NULL,VIS_OPT_TEST,levelMg,0,0,NULL
    ,[&]( float f, MtStat *fstat,
	int flags, int blockFlags, int level, const Vec3Int& ipos,
	int &colorOut ) -> void
    {
      Vec3Int rgb=0;
      if(blockFlags & flagRed) rgb[0] = 255;
      if(blockFlags & flagGreen) rgb[1] = 255;
      if(blockFlags & flagBlue) rgb[2] = 255;
      colorOut = BMP_MKRGB( (int)rgb.x, (int)rgb.y, (int)rgb.z );
    } );
}

} // namespace MSBG

