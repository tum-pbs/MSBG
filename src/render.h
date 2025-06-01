/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#pragma once
#include <vector>
#include <array>
#include <cmath>
#include <cstdint>
#include <algorithm>

namespace RDR
{

  typedef double RealType;


struct Vec3 {
    RealType x, y, z;

    // Constructors
    constexpr Vec3() : x(0), y(0), z(0) {}
    constexpr Vec3(RealType xx, RealType yy, RealType zz) : x(xx), y(yy), z(zz) {}

    // Arithmetic
    Vec3 operator + (const Vec3& rhs) const { return {x+rhs.x, y+rhs.y, z+rhs.z}; }
    Vec3 operator - (const Vec3& rhs) const { return {x-rhs.x, y-rhs.y, z-rhs.z}; }
    Vec3 operator * (RealType s) const         { return {x*s,   y*s,   z*s}; }
    Vec3 operator / (RealType s) const         { return {x/s,   y/s,   z/s}; }

    Vec3& operator += (const Vec3& rhs) { x+=rhs.x; y+=rhs.y; z+=rhs.z; return *this; }
    Vec3& operator *= (RealType s)        { x*=s; y*=s; z*=s; return *this; }
};

inline Vec3 operator * (RealType s, const Vec3& v) { return v * s; }
inline RealType dot(const Vec3& a, const Vec3& b)  { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline Vec3  cross(const Vec3& a, const Vec3& b){
    return { a.y*b.z - a.z*b.y,
             a.z*b.x - a.x*b.z,
             a.x*b.y - a.y*b.x };
}
inline RealType length(const Vec3& v){ return std::sqrt(dot(v,v)); }
inline Vec3  normalize(const Vec3& v){
    RealType len = length(v);
    return len>1e-20f ? v / len : Vec3(0,1,0);
}

class RaymarchRenderer {
public:
    RaymarchRenderer(int w, int h, SBG::SparseGrid<uint16_t> *sgDensity)
        : _width(w), _height(h), _sgDens(sgDensity)
    {
        _res = _sgDens->sxyzMax();
        _h = 1.0/(double)_res;
        setCamera({0,0,-5}, {0,0,0}, 1.0f);
        setSunLight({10,10,-10});
        setSurfaceColor({1,1,1});
    }

    void setCamera(const Vec3& position,
                   const Vec3& lookAt,
                   float       focalLength /* >0 == zoom factor */)
    {
        _camPos   = position;
        _lookAt   = lookAt;
        _focalLen = focalLength;
        // pre–compute camera basis (simple look-at with world-up = Y+)
        const Vec3 worldUp{0,1,0};
        _forward = normalize(lookAt - position);
        _right   = normalize(cross(_forward, worldUp));
        _up      = cross(_right, _forward);
    }

    void setSunLight(const Vec3& sunPos){ _sunPos = sunPos; }
    void setSurfaceColor(const Vec3& color){ _surfaceColor = color; }


    float sampleDensity(const Vec3& p) const
    {
      Vec4f pos = _res * Vec4f(p.x,p.y,p.z,0);
      float f = 0;
      //if(_sgDens->inRange( pos ))
      {
        Vec4i ipos = truncate_to_int(floor(pos));	
    	if( _sgDens->inRange(ipos)) 
	{
	  int bid = _sgDens->getBlockIndex(_sgDens->getBlockCoords( ipos ));
	  // Skip fast over empty blocks
	  if(_sgDens->isEmptyBlock(bid)) return 0.f;

	  #if 0
	  constexpr float eps=1e-1;
	  if( std::max( std::max( std::min( fabsf(p.x-0.f), fabsf(p.x-1.f) ),
				  std::min( fabsf(p.y-0.f), fabsf(p.y-1.f) ) ),
			std::min( fabsf(p.z-0.f), fabsf(p.z-1.f) ) ) 
	      < eps )
	  {	  
	    return .9f;  
	  }
	  #endif

	  f = _sgDens->interpolateUint16ToFloatFast( pos, 0 );
	  RENDER_DENS_TO_FLOAT_UI16(f);
	  UT_ASSERT2(std::isfinite(f) && !( f<0.0f || f>1.0f));
	}
      }
      return f;
    }

    Vec3 sampleDensityGradient(const Vec3& p) const
    {
      float outVal=0;
      Vec4f outGrad(0);
      Vec4f pos = _res * Vec4f(p.x,p.y,p.z,0);
      if(_sgDens->inRange( pos ))
      {
	_sgDens->interpolateWithDerivs<0>( SBG::IP_BSPLINE_CUBIC, pos, 0, 
					      &outVal, &outGrad); 
	//RENDER_DENS_TO_FLOAT_UI16(outVal);
	RENDER_DENS_TO_FLOAT_UI16(outGrad); 
      }
      return -1*Vec3(vfget_x(outGrad),vfget_y(outGrad),vfget_z(outGrad));
    }

    /**
     * Render the scene into an RGB byte buffer of size width*height*3.
     * Pixels are stored row-major, top-to-bottom, left-to-right.
     */
    BmpBitmap *render(void) const;

    int width()  const { return _width;  }
    int height() const { return _height; }

private:
    int  _width  = 0;
    int  _height = 0;

    RealType _h;
    int	_res;

    //
    SBG::SparseGrid<uint16_t> *_sgDens;

    // camera
    Vec3 _camPos{};
    Vec3 _lookAt{};
    RealType _focalLen = 1.0f;
    Vec3 _forward{}, _right{}, _up{};

    // lighting & material
    Vec3 _sunPos{};
    Vec3 _surfaceColor{};
};

}

