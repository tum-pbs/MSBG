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
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <stdarg.h>
#include <setjmp.h>
#include <zlib.h>
#include <png.h>
#include "globdef.h"
#include "mtool.h"
#include "util.h"
#include "bitmap.h"

typedef int            BOOL_;
typedef int32_t       LONG_;
typedef uint32_t      DWORD_;

typedef png_color PALETTE;
typedef struct 
MYIMAGE_s 
{
  int32_t    width;
  int32_t    height;
  unsigned int    pixdepth;
  unsigned int    palnum;
  BOOL_    topdown;
  BOOL_    alpha;
  DWORD_   rowbytes;
  LongInt imgbytes;
  PALETTE *palette;
  unsigned char    **rowptr;
  unsigned char    *bmpbits;
  png_color_8 sigbit;
} 
MYIMAGE;

static int  interlace = 0;
//static int  complevel = 6;
static int  complevel0 = 0;
static int  filters   = 0;

static void printfxx(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  fflush(stderr);
  va_end(ap);
}

static void my_png_error(png_structp png_ptr, png_const_charp message)
{
  printfxx("ERROR(libpng): %s - %s\n", message,
      (char *)png_get_error_ptr(png_ptr));
  longjmp(png_jmpbuf(png_ptr), 1);
}


static void my_png_warning(png_structp png_ptr, png_const_charp message)
{
  printfxx("WARNING(libpng): %s - %s\n", message,
	       (char *)png_get_error_ptr(png_ptr));
}


static void myimage_free(MYIMAGE *img)
{
  MM_free(img->palette);
  MM_free(img->rowptr);
  MM_free(img->bmpbits);
}

static void myimage_init(MYIMAGE *img)
{
  img->palette = NULL;
  img->rowptr  = NULL;
  img->bmpbits = NULL;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static BOOL_ myimage_alloc(MYIMAGE *img)
{
  unsigned char *bp, **rp;
  LONG_ n;

  if (img->palnum > 0) 
  {
    img->palette = MM_malloc((size_t)img->palnum * sizeof(PALETTE),MM_DUMMYTAG);
    if (img->palette == NULL) { myimage_init(img); return FALSE; }
  } 
  else 
  {
    img->palette = NULL;
  }
  
  img->rowbytes = ((DWORD_)img->width * img->pixdepth + 31) / 32 * 4;
  img->imgbytes = img->rowbytes * (LongInt)img->height;
  img->rowptr   = MM_malloc((size_t)img->height * sizeof(unsigned char *),MM_DUMMYTAG);
  img->bmpbits  = MM_malloc((size_t)img->imgbytes,MM_DUMMYTAG);

  if (img->rowptr == NULL || img->bmpbits == NULL) 
  {
    myimage_free(img); myimage_init(img); return FALSE;
  }

  n  = img->height;
  rp = img->rowptr;
  bp = img->bmpbits;

  if (img->topdown) 
  {
    while (--n >= 0) 
    {
      *(rp++) = bp;
      bp += img->rowbytes;
    }
  } 
  else 
  {	
    bp += img->imgbytes;
    while (--n >= 0) 
    {
      ((DWORD_ *)bp)[-1] = 0;
      bp -= img->rowbytes;
      *(rp++) = bp;
    }
  }

  return TRUE;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int BmpWritePNG(FILE *fp, BmpBitmap *bmp, CoConverter *cc,
    		unsigned opt )
{
  TRC(("BmpWritePNG()\n"));
	int img_allocated=FALSE,use_bmp_ushort=FALSE;
	png_structp png_ptr;
	png_infop info_ptr;
	int color_type,crgb,cr,cg,cb,nOog=0,k;
	int interlace_type;
	MYIMAGE	img0,*img=NULL;
	int bit_depth,channels,x,y,gli;
	float gl,vRgba[4];
	char *fn=bmp->name;
	unsigned char	*rp;
        UtProgrInd pgi;
	int	doPgi=FALSE;

//	printf("BmpWritePNG: opt = %d\n",opt);

	/*-----------------------------------------------------------------*/
	/* Setup IMG from BMP structure					   */
	/*-----------------------------------------------------------------*/
	if( (opt & BMP_16BIT ))
	{
	  if((bmp->dataUV&&(!(opt&BMP_GREY_ONLY)))||
	    	  ((opt&BMP_RGB)&&(opt&BMP_XYZ))||
	    	  ((opt&BMP_RGB)&&(opt&BMP_USHORT))		  
		  )
	  {
	    //
	    // 16 Bit RGB
	    //
	    bit_depth = 16;
	    channels = 3;
	    color_type = PNG_COLOR_TYPE_RGB;
	  }
	  else if(bmp->dataGrey)
	  {
	    //
	    // 16 Bit Grey
	    //
	    bit_depth = 16;
	    channels = 1;
	    color_type = PNG_COLOR_TYPE_GRAY;
	  }
	  else if(bmp->dataUshort0)
	  {
	    //
	    // 16 Bit Grey (unsigned short)
	    //
	    bit_depth = 16;
	    channels = 1;
	    color_type = PNG_COLOR_TYPE_GRAY;
	    use_bmp_ushort=TRUE;
	  }
	  else
	  {
	    TRCERR(("inufficient data channel(s) for PNG 16 bit\n"));
	    return 0;
	  }
	}
	else
	{
	  //
	  // 8 Bit RGB
	  //
	  bit_depth = 8;
	  channels = 3;
	  color_type = PNG_COLOR_TYPE_RGB;
	}

	if((opt & BMP_ALPHA) && bmp->dataAlpha)
	{
	  color_type |= PNG_COLOR_MASK_ALPHA;
	  channels+=1;
	  UT_ASSERT0(channels<=4);
	}

	img = &img0;
	myimage_init(img);
	img->width    = (LONG_)bmp->sx;
	img->height   = (LONG_)bmp->sy;
	img->pixdepth = (unsigned int)bit_depth * channels;
	img->palnum   = (img->pixdepth <= 8) ? (1 << img->pixdepth) : 0;
	img->topdown  = TRUE;
	img->alpha    = FALSE;

	if (!myimage_alloc(img)) {
		TRCERR(("out of memory\n"));
		return 0;
	}
	img_allocated = TRUE;

	if(bmp->sy>=10000)
	{
          UT_PROGRIND_INIT( &pgi );
          UT_PROGRIND_CONF( &pgi, 0.05, 1 );
	  //doPgi = TRUE;
	}

	/*-----------------------------------------------------------------*/
	/* Convert data from BMP to IMG					   */
	/*-----------------------------------------------------------------*/
	for(y=0;y<bmp->sy;y++)
	{
	  if(doPgi)
	  {
            UT_PROGRIND( &pgi, y, 0, bmp->sy, "saving PNG ", NULL );
	  }
	  rp = img->rowptr[y];
	  for(x=0;x<bmp->sx;x++)
	  {
	    if(bit_depth==8)
	    {
	      BMP_GETPIXEL(bmp,x,y,crgb,cr,cg,cb);
	      rp[2] = cr;
	      rp[1] = cg;
	      rp[0] = cb;
	      rp += 3;
	    }
	    else if(channels==1)
	    {
	      if(use_bmp_ushort)
	      {
		gli = bmp->dataUshort0[ BMP_GETPIXEL_OFFS(bmp,x,y)];
	      }
	      else
	      {
		gl = BMP_GETPIXEL_GREY(bmp,x,y);
		gli = (int)(gl*256.0);
		if(gli>65535) gli = 65535;
	      }
	      rp[0] = gli>>8;
	      rp[1] = gli&0xFF;
	      rp+=2;
	    }
	    else if(channels>=3)
	    {
	      if(opt&BMP_USHORT)
	      {
		for(k=0;k<3;k++)
		  vRgba[k]=BMP_PIXEL(bmp,dataUshort[k],x,y)/((double)UINT16_MAX);
	      }
	      else if(opt&BMP_XYZ)
	      {
		for(k=0;k<3;k++)
		  vRgba[k]=BMP_PIXEL(bmp,dataFloat[k],x,y)/255.0;
	      }
	      else
	      {
		UT_ASSERT0(FALSE);
	      }
	      
	      if(color_type & PNG_COLOR_MASK_ALPHA) 
	      {
		for(k=3;k>=1;k--) vRgba[k]=vRgba[k-1];
		vRgba[0]=BMP_PIXEL(bmp,dataAlpha,x,y);
		if( (!isfinite(vRgba[0]))||vRgba[0]<0||vRgba[0]>1)
		{
	          TRCERR(("Assertion failed\n"));
		  goto error_abort;
	        }
	      }

	      for(k=channels-1;k>=0;k--)
	      {
		gl = vRgba[k];
		gli = (int)(gl*65535);
		if(gli>65535) gli = 65535;
		if(gli<0) gli = 0;
		rp[0] = gli>>8;
		rp[1] = gli&0xFF;
		rp+=2;
	      }
	    }
	    else
	    {
	      TRCERR(("unsupported PNG format\n"));
	      goto error_abort;
	    }
	  }
	}

	if(nOog)
	{
	  TRC1(("WARNING: %d of %d pixels out of RGB color gamut \n",
		nOog,bmp->sx*bmp->sy));
	}

	png_ptr = png_create_write_struct( PNG_LIBPNG_VER_STRING, fn,
	                                     my_png_error, my_png_warning );
	if( png_ptr==NULL ){
	  	TRCERR(("out of memory\n"));
		return 0;
	}
	info_ptr = png_create_info_struct( png_ptr );
	if( info_ptr==NULL ){
		png_destroy_write_struct( &png_ptr, NULL );
	  	TRCERR(("out of memory\n"));
		return 0;
	}
	if( 
	    #ifdef MIMP_ON_LINUX
	    setjmp(png_jmpbuf(png_ptr)) 
	    #else
	    //  setjmp(png_ptr->jmpbuf) 
	    setjmp(png_jmpbuf(png_ptr)) 
	    #endif
	    )
	{
		png_destroy_write_struct( &png_ptr, &info_ptr );
		return 0;
	}
	png_init_io( png_ptr, fp );
	png_set_compression_level( png_ptr, 
	    			  (opt & BMP_COMPR) ? 3 : complevel0 );
	if( filters!=0 )
		png_set_filter( png_ptr, PNG_FILTER_TYPE_BASE, filters );


	png_set_compression_mem_level( png_ptr, MAX_MEM_LEVEL );
	interlace_type = (interlace)? PNG_INTERLACE_ADAM7 : PNG_INTERLACE_NONE ;

	png_set_IHDR( png_ptr, info_ptr, img->width, img->height, bit_depth,
	              color_type, interlace_type, PNG_COMPRESSION_TYPE_DEFAULT,
	              PNG_FILTER_TYPE_DEFAULT );

	png_write_info( png_ptr, info_ptr );


	if( color_type == PNG_COLOR_TYPE_RGB ||
	    color_type == (PNG_COLOR_TYPE_RGB | PNG_COLOR_MASK_ALPHA)) 
	  png_set_bgr( png_ptr );

	png_write_image( png_ptr, img->rowptr );

	png_write_end( png_ptr, info_ptr );
	png_destroy_write_struct( &png_ptr, &info_ptr );


	if(img_allocated) myimage_free(img);
	return TRUE;

error_abort:
	if(img_allocated) myimage_free(img);
	return FALSE;

}

