/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/
//#define USE_OPENEXR

#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <jpeglib.h>
#ifdef USE_OPENEXR
#include <openexr/ImathBox.h>
#include <openexr/ImfRgbaFile.h>
#include <openexr/ImfArray.h>
#include <openexr/half.h>
#include <openexr/ImfAttribute.h>
#include <openexr/ImfInputFile.h>
#include <openexr/ImfOutputFile.h>
#include <openexr/ImfChannelList.h>
#include <openexr/ImfFrameBuffer.h>
#include <openexr/half.h>
#endif
#include "globdef.h"
#include "util.h"
#include "mtool.h"
#include "bitmap.h"

MSBG_NAMESPACE_BEGIN

/*=========================================================================
 *
 *
 *   OpenEXR
 *
 *
 * =======================================================================*/

#ifdef USE_OPENEXR
using namespace OPENEXR_IMF_NAMESPACE;
using namespace std;
using namespace Imf;
using namespace Imath;
#endif

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpSaveBitmapEXR( BmpBitmap *B, char *filename, 
    		      int chan, int ichan,
    		      unsigned opt )
{
  // supported options:
  //
  // BMP_GREY 	:	save only one greylevel channel
  // BMP_XYZ 	:	save three channels dataFloat[0..2] as RGB image
  // BMP_32BIT 	:	use 32 bit precision
  // BMP_16BIT 	:	use 16 bit precsision (use 'HALF' formats for floats)
  // BMP_COMPR  :	use compression
  //

  int rcThis=0;
#ifdef USE_OPENEXR
  int k,rc;    
  int bit_depth = opt & BMP_32BIT ? 32 : 16,
      do_float = opt & BMP_FLOAT,
      do_compr = opt & BMP_COMPR,
      do_grey = opt & BMP_GREY;
  Imf::RgbaOutputFile*   fout_rgba=NULL;
  Imf::OutputFile*   fout=NULL;
  BmpData bmpData;

  //int do_grey_uint16 = (opt & BMP_16BIT) && (opt & BMP_GREY) && (!(opt & BMP_FLOAT));

  USE(bit_depth);
  UT_ASSERT( opt & BMP_FLOAT );  // unfortunately OpenEXR only supports float formats

  TRC(("%s: name='%s' opt=%d\n",UT_FUNCNAME,filename,opt));
  UT_ASSERT(do_grey || (opt&BMP_XYZ));

  rc = BmpGetDataChannel(B,chan,ichan,0,&bmpData);
  UT_ASSERT(rc==0);

  if(do_grey)
  {
    if(do_float)
    {
      UT_ASSERT(bmpData.pFloat);    
    }
    else
    {
      UT_ASSERT(bmpData.pUshort);    
    }
  }
  else
  {
    for(k=0;k<3;k++) UT_ASSERT(B->dataFloat[k]);
  }

#if 0
  {
    int nProc = MAX(UtGetNumPhysProc( ) / 1,1);
    Imf::setGlobalThreadCount (nProc);
    TRC(("%s: using %d OpenEXR worker threads\n",
	    UT_FUNCNAME,nProc));
// Imf::setGlobalThreadCount (0);
  }
#endif
#if 0
  if(do_grey_uint16)
  {
    UT_ASSERT( FALSE );  // unfortunately OpenEXR only supports float formats
    const char *chan_name = "Y";
    Header header (B->sx, B->sy); 
    header.channels().insert (chan_name, Channel (UINT));
    try 
    {
      fout = new OutputFile(filename, header);
    }
    catch (const std::exception &e)
    {
      fout = NULL;
      TRCERRR(("%s: %s\n",UT_FUNCNAME,e.what()),MI_ERROR);
    }

    uint16_t *pixels = bmpData.pUshort;

    FrameBuffer frameBuffer; // 5
    frameBuffer.insert (chan_name, // name // 6
    			Slice (UINT, // type // 7
    			(char *) pixels, // base // 8
    			sizeof (*pixels) * 1, // xStride// 9
    			sizeof (*pixels) * B->sx)); // yStride// 10
    try 
    {
      fout->setFrameBuffer (frameBuffer); // 16
      fout->writePixels (B->sy);
    }
    catch (std::exception &e)
    {
      TRCERRR(("%s: %s\n",UT_FUNCNAME,e.what()),MI_ERROR);
    }
  }
  else
#endif
  {
    // use simplified RGBA interface 
    // Open file 
    //
    try 
    {
      fout_rgba = new Imf::RgbaOutputFile( 
	  			filename, B->sx, B->sy, 
	  			Imf::WRITE_RGB, 
				1, Imath::V2f (0, 0), 1,
				Imf::INCREASING_Y,
				do_compr ? 
					#if 0
					  Imf::PIZ_COMPRESSION : 
					#else
					  Imf::ZIP_COMPRESSION : 
					#endif
					Imf::NO_COMPRESSION

				);
    }
    catch (const std::exception &e)
    {
      fout_rgba = NULL;
      TRCERRR(("%s: %s\n",UT_FUNCNAME,e.what()),MI_ERROR);
    }
    //
    // allocate the scanline buffer
    //
    const LongInt buf_max_pixels = B->sx*B->sy,   // 250*ONE_MB,  // must transfer all in one chunk: TODO: solve OpenEXR problem with chunked transfer ! 
	  	  buf_sy = MIN( buf_max_pixels / B->sx, B->sy );

    TRC(("allocating line buffer of %dx%d pixels\n",B->sx,buf_sy));
	  	
    UT_ASSERT_FATAL(buf_sy*(LongInt)B->sx<2*(LongInt)ONE_GB);  // OpenEXR limitation (?)
    Imf::Array<Imf::Rgba> pixels (buf_sy * B->sx);

    BmpRectangle dw;
    BMP_SET_RECT2(&dw,0,0,B->sx-1,B->sy-1);

    LongInt y0=0, n_chunk=0;
    while ( dw.y0 <= dw.y2 )
    {
      // fill buffer 
      int out_sy = MIN(buf_sy,dw.y2+1-dw.y0);
       TRC(("transferring %d lines %d - %d of %d\n",
	     out_sy,y0,y0+out_sy-1,B->sy));
      #pragma omp parallel for
      for(int y=0;y<out_sy;y++)
        for(int x=0;x<B->sx;x++)
	{
	  int x3=x,y3=y+y0;
	  UT_ASSERT_FATAL(BMP_IN_RANGEB(x3,y3,B));
	  LongInt ii=BMP_XYOFF(B,x3,y3),
		  jj=MT_GXY(B->sx,x,y);
	  float col[4];
	  for(int k=0;k<3;k++) col[k]=B->dataFloat[k][ii];  
	  Imf::Rgba* rgba = pixels + jj;
	  rgba->r = col[0];
	  rgba->g = col[1];
	  rgba->b = col[2];
	  rgba->a = 1;
	}
      // write buffer to file
      try 
      {
	fout_rgba->setFrameBuffer (&pixels[0], 1, B->sx);
	//TRCP(("out_sy=%d\n",out_sy));
	
	UT_ASSERT_FATAL(n_chunk<1); // TODO investigate problem with repeated call to writePixels() 
	fout_rgba->writePixels (out_sy);
      }
      catch (std::exception &e)
      {
	TRCERRR(("%s: %s\n",UT_FUNCNAME,e.what()),MI_ERROR);
      }
    
      dw.y0 += buf_sy;
      y0 += buf_sy;
      n_chunk++;
    }
  }
#else  // USE_OPENEXR
  UT_ASSERT(FALSE);
#endif
rcCatch:
#ifdef USE_OPENEXR
  if(fout) delete fout;
  if(fout_rgba) delete fout;
#endif
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpLoadBitmapEXR( char *filename, 
    		      unsigned opt,  // BMP_XYZ | BMP_32BIT 
    		      BmpBitmap **pBmp )
{
  int rcThis=0;
  BmpBitmap *B=NULL;
#ifdef USE_OPENEXR

  TRC(("%s: name='%s' opt=%d\n",UT_FUNCNAME,filename,opt));

  Imf::RgbaInputFile*   fin=NULL;
  UT_ASSERT(opt & (BMP_XYZ));

  {
    #if 0
    int nProc = MAX(UtGetNumPhysProc( ) / 1,1);
    TRC(("%s: using %d OpenEXR worker threads\n",
	  UT_FUNCNAME,nProc));
    Imf::setGlobalThreadCount (nProc);
    #endif

    //
    // Open file 
    //
    UtTimer tm;
    TIMER_START(&tm);

    try 
    {
      fin = new Imf::RgbaInputFile( filename );
    }
    catch (const std::exception &e)
    {
      fin = NULL;
      TRCERRR(("%s: %s\n",UT_FUNCNAME,e.what()),MI_ERROR);
    }

    TIMER_STOP(&tm);

    TRC(("CPU (Open EXR) %.1f sec.\n", (double)(TIMER_DIFF_MS(&tm)/1000.0))); 

    //
    // Read header
    //
    const Imf::Header& header = fin->header();

    BmpRectangle dw;
    BMP_SET_RECT2(&dw,  
		header.dataWindow().min.x,  header.dataWindow().min.y,
		header.dataWindow().max.x,  header.dataWindow().max.y
	);

    int fWidth = dw.sx,
	fHeight = dw.sy;
    const Imf::RgbaChannels channels = fin->channels();
    int fChannels = channels & Imf::WRITE_A ? 4 : 3;

    TRCP(("width=%d, height=%d, channels=%d\n",
	  fWidth,fHeight,fChannels));

    UT_ASSERT(fChannels==3);   // no alpha

    //
    // Create destination bitmap
    //
    B = BmpNewBitmap(fWidth,fHeight,BMP_XYZ);

    //
    // Read pixels (scanline mode)
    //
    const LongInt buf_max_pixels = B->sx*B->sy,   // 250*ONE_MB,  // must transfer all in one chunk: TODO: solve OpenEXR problem with chunked transfer ! 
	  	  buf_sy = MIN( buf_max_pixels / B->sx, B->sy );

    TRC(("allocating line buffer of %dx%d pixels\n",B->sx,buf_sy));
	  	
    UT_ASSERT_FATAL(buf_sy*(LongInt)B->sx<2*(LongInt)ONE_GB);  // OpenEXR limitation

    Imf::Array<Imf::Rgba> pixels (buf_sy * B->sx);
    LongInt y0=0;
    while ( dw.y0 <= dw.y2 )
    {
      int in_sy = MIN(buf_sy,dw.y2+1-dw.y0);
      TRC(("transferring %d lines %d - %d of %d\n",
	   in_sy,y0,y0+in_sy-1,B->sy));
      try 
      {
	fin->setFrameBuffer (&pixels[0], 1, B->sx);

	TIMER_START(&tm);

	fin->readPixels (dw.y0, dw.y0+in_sy-1);

	TIMER_STOP(&tm);

	TRCP((" fin->readPixels %.1f sec. (%.0f pixels/sec)\n", 
	      (double)(TIMER_DIFF_MS(&tm)/1000.0), 
	      (B->sx*B->sy)/(double)(TIMER_DIFF_MS(&tm)/1000.0)));

      }
      catch (std::exception &e)
      {
	TRCERRR(("%s: %s\n",UT_FUNCNAME,e.what()),MI_ERROR);
      }

      // transfer pixels
      #pragma omp parallel for
      for(int y=0;y<in_sy;y++)
        for(int x=0;x<B->sx;x++)
	{
	  int x3=x,y3=y+y0;
	  UT_ASSERT_FATAL(BMP_IN_RANGEB(x3,y3,B));
	  LongInt ii=BMP_XYOFF(B,x3,y3),
		  jj=MT_GXY(B->sx,x,y);

	  Imf::Rgba* src = pixels + jj;
	  float col[4];

	  col[0] = src->r;
	  col[1] = src->g;
	  col[2] = src->b;
	  col[3] = fChannels==4 ? (float)src->a : 1;

	  for(int k=0;k<3;k++) B->dataFloat[k][ii]=col[k];  
	}

      dw.y0 += buf_sy;
      y0 += buf_sy;
    }
  }
#else	// USE_OPENEXR
  UT_ASSERT(FALSE);
#endif
rcCatch:
#ifdef USE_OPENEXR
  if(fin) delete fin;
#endif
  if(rcThis)
  {
    BmpDeleteBitmap(&B);
  }
  *pBmp=B;

  return rcThis;
}

/*=========================================================================
 *
 *
 *   JPG
 *
 *   needs separate source file due to conflicts with headers 
 *   (windows.h needed by grwin.h etc.)
 *
 * =======================================================================*/

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpSaveBitmapJPG( BmpBitmap *bmp, char *fname, 
    			int quality,  // [0..100] 0=use-default
    			unsigned opt )
{
  int rcThis=0;
  LongInt i;
  FILE *outfile=NULL;
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  unsigned char *raw_image=NULL;

  TRC(("BmpSaveBitmapJPG '%s' q=%d opt=%d\n",fname,quality,opt));

  // Assumed input: RGB in dataFloat[0..2]  data range [0..1]
  UT_ASSERT(opt & (BMP_XYZ|BMP_RGB));
  //for(k=0;k<3;k++) UT_ASSERT(bmp->dataFloat[k]);
  
  {
    int bytes_per_pixel=3;
    ALLOCARR_(raw_image,unsigned char,bmp->sx*bmp->sy*bytes_per_pixel);
    int doParallel=TRUE;
    USE(doParallel);
    #pragma omp parallel for if(doParallel)
    for(i=0;i<bmp->sx*bmp->sy;i++)
    {
      int k;
      if(opt&BMP_RGB)
      {
	int c = bmp->data[i];
	raw_image[i*3+0]=BMP_RED(c);
	raw_image[i*3+1]=BMP_GREEN(c);
	raw_image[i*3+2]=BMP_BLUE(c);
      }
      else for(k=0;k<3;k++) 
      {
	int c=bmp->dataFloat[k][i];
	MT_CLAMP(c,0,255);
	raw_image[i*3+k]=c;
      }
    }

    /* this is a pointer to one row of image data */
    JSAMPROW row_pointer[1];    
    outfile = fopen(fname,"wb");
    if(!outfile) 
      TRCERRR(("can't open output file '%s'\n",fname),MI_ERROR);

    cinfo.err = jpeg_std_error( &jerr );
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo, outfile);

    /* Setting the parameters of the output file here */
    cinfo.image_width = bmp->sx;	
    cinfo.image_height = bmp->sy;
    cinfo.input_components = bytes_per_pixel;
    cinfo.in_color_space = JCS_RGB; /* or JCS_GRAYSCALE for grayscale images */

    /* default compression parameters, we shouldn't be worried about these */
    jpeg_set_defaults( &cinfo );
    /* set quality (0..100) */
    jpeg_set_quality (&cinfo, quality ? quality : 100, TRUE);
    /* Now do the compression .. */
    jpeg_start_compress( &cinfo, TRUE );
    /* like reading a file, this time write one row at a time */
    while( cinfo.next_scanline < cinfo.image_height )
    {
      row_pointer[0] = raw_image + cinfo.next_scanline * cinfo.image_width *  cinfo.input_components;
      jpeg_write_scanlines( &cinfo, row_pointer, 1 );
    }

    /* similar to read file, clean up after we're done compressing */
    jpeg_finish_compress( &cinfo );
    jpeg_destroy_compress( &cinfo );
  }

rcCatch:
  if(outfile) fclose(outfile);
  FREEMEM(raw_image);
  return(rcThis);
}

MSBG_NAMESPACE_END
