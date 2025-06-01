/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#include "common_headers.h"
#include "render.h"

namespace RDR
{
  // 
  // Very simple demo raymarching renderer for MSBG/SBG grids
  //
  BmpBitmap *RaymarchRenderer::render( void ) const
  {
    BmpBitmap *bmpOut=BmpNewBitmap(_width,_height,BMP_RGB);

    constexpr RealType ISO_LEVEL   = 0.5f;
    RealType step_size   = _h;
    constexpr RealType MAX_DIST    = 2.0f;
    int   max_steps   = int(MAX_DIST / step_size);

    const RealType invW = 1.0f / RealType(_width);
    const RealType invH = 1.0f / RealType(_height);
    const RealType aspect = RealType(_width) / RealType(_height);

    ThrRunParallel<int>( _height, nullptr,
      [&](int &tls, int tid, size_t y) 
      {
	for(int x=0; x<_width; ++x)
	{
	    // --- Generate primary ray ---
	    RealType ndcX = ( (x + 0.5f) * invW - 0.5f ) * 2.0f * aspect;
	    RealType ndcY = ( 0.5f - (y + 0.5f) * invH ) * 2.0f;

	    Vec3 rayDir = normalize(_forward * _focalLen + _right * ndcX + _up * ndcY);
	    Vec3 rayPos = _camPos;

	    // --- Ray–march ---
	    RealType prevDensity = sampleDensity(rayPos);
	    Vec3  prevPos     = rayPos;

	    RealType traveled = 0.0f;
	    bool  hit      = false;
	    Vec3  hitPos{};

	    for(int step=0; step<max_steps; ++step)
	    {
		traveled += step_size;
		rayPos += rayDir * step_size;
		RealType density = sampleDensity(rayPos);

		if(prevDensity < ISO_LEVEL && density >= ISO_LEVEL)
		{
		    // Linear interpolation to approximate surface
		    RealType t = (ISO_LEVEL - prevDensity) / (density - prevDensity);
		    hitPos  = prevPos + (rayPos - prevPos) * t;
		    hit = true;
		    break;
		}

		prevDensity = density;
		prevPos     = rayPos;
	    }

	    // --- Shading / write pixel ---

	    int colOut = BMP_MKRGB(0,0,0);
	    if(hit)
	    {
		Vec3 normal = normalize(sampleDensityGradient(hitPos));
		Vec3 lightDir = normalize(_sunPos - hitPos);
		RealType diff = std::max(dot(normal, lightDir), (RealType)0.0f);

		Vec3 col = _surfaceColor * diff; // simple Lambert
		col.x = std::clamp(col.x, (RealType)0.0f, (RealType)1.0f);
		col.y = std::clamp(col.y, (RealType)0.0f, (RealType)1.0f);
		col.z = std::clamp(col.z, (RealType)0.0f, (RealType)1.0f);

		 colOut = BMP_MKRGB( uint8_t(col.x * 255.0f + 0.5f),
		    	                 uint8_t(col.y * 255.0f + 0.5f),
		                         uint8_t(col.z * 255.0f + 0.5f) );
	    }

	    int idx = BMP_XYOFF(bmpOut,x,y);
	    UT_ASSERT0(idx>=0&&idx<bmpOut->sx*bmpOut->sy);
	    bmpOut->data[idx] = colOut;
	}
    }
  );

    return bmpOut;
  }

} // namespace

