/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#ifndef MTOOL2_H
#define MTOOL2_H

#include <algorithm>
#include <limits>
#include "util.h"

MSBG_NAMESPACE_BEGIN

namespace MT
{
  inline float sqf( const float& x ) { return x*x; }

  inline double relativeError( double a, 
      		   double b, 
		   double epsilon=std::numeric_limits<double>::epsilon() )
  {      

    //TRCP(("%s : %g %g %g\n",UT_FUNCNAME,a,b,epsilon));

    // https://floating-point-gui.de/errors/comparison/
    double absA = std::abs(a),
	      absB = std::abs(b),
              diff = std::abs(a - b);

    double minNormal = std::numeric_limits<double>::min(),
	      maxValue = std::numeric_limits<double>::max();

    //TRCP(("%s : %g %g %g %g %g\n",UT_FUNCNAME,absA,absB,diff,minNormal,maxValue));


    UT_ASSERT0(std::numeric_limits<double>::epsilon() <= epsilon);
    UT_ASSERT0(epsilon < 1);

    if (a == b) 
    { // shortcut, handles infinities
        //TRCP(("%s : a==b\n", UT_FUNCNAME));
        return 0;
    } 
    else if (a == 0 || b == 0 || (absA + absB < minNormal))
    {
        // a or b is zero or both are extremely close to it

        // relative error is less meaningful here
        //TRCP(("%s : absA,absB < min: %g %g\n", UT_FUNCNAME, absA + absB,epsilon*minNormal));
        return diff < (epsilon * minNormal) ? 0 : 1;
    } 
    else 
    { // use relative error
        double e = diff / std::min( (absA + absB), maxValue );
        //TRCP(("%s : %g\n", UT_FUNCNAME, e));
	return e;
    }
  }
      
  template<typename FloatType=double, typename AccuType=double > class Stats
  {
    public:

    Stats( void ) 
    {
      reset();
    };

    void update( FloatType val )
    {
      UT_ASSERT2_NF( std::isfinite( val ));
      _n++;
      _min = std::min( _min, val );
      _max = std::max( _max, val );
      _sum += val;
      _sumSq += val*val;
    }

    void update( LongInt n, FloatType min, FloatType max, 
		 AccuType sum, AccuType sumSq )
    {
      UT_ASSERT2_NF( UT_IMPLIES(n>0,!(max<min)) ); 

      _n += n; 
      
      _min = std::min( _min, min );
      _max = std::max( _max, max );

      _sum += sum; 
      _sumSq += sumSq; 
    }

    void update( const Stats& stat2 )
    {
      update( stat2._n, stat2._min, stat2._max, stat2._sum, stat2._sumSq );
    }

    void result( void )
    {
      UT_ASSERT2_NF( std::isfinite( _sum ));
      UT_ASSERT2_NF( std::isfinite( _sumSq ) && !(_sumSq<0));
      AccuType nInv = _n>0 ? 1.0/_n : 0;
      _avg = _n > 0  ? (_sum * nInv ) : 0; 
      _sig =  _n > 0 ? ( (_sumSq * nInv ) - _sum*nInv * _sum*nInv ) : 0; 
      if(_sig < 0) _sig = 0; 
      _sig = sqrt(_sig); 
      _span = _max-_min;       
    }

    void reset()
    {
      _n =
      _sum = 
      _sumSq = 
      _span = 
      _avg = 
      _sig = 0;
      _min = std::numeric_limits<FloatType>::max();
      _max = std::numeric_limits<FloatType>::lowest();
    }

    LongInt _n;

    AccuType _sum,
	     _sumSq;

    FloatType _min, 
	      _max,
	      _span,
	      _avg,
	      _sig;
  };

} // namespace MT

MSBG_NAMESPACE_END

#endif

