#ifndef STDMATH_H
#define STDMATH_H

#ifdef MI_WITH_ICC

#include <mathimf.h>

#else

#include <math.h>

#endif

// isnan() implementation for MSVC
#ifdef MI_WITH_MSVC
#if 0
#include <iostream>
#include <cmath>

int isnan(double x) { 
  return std::isnan(x);
}

int isinf(double x) { 
  return std::isinf(x);
}
#endif
#define isnan(x) _isnan(x)
#endif

#endif  /* STDMATH_H */

