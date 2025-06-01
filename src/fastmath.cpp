/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/
#include <memory>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <vectorclass/vectorclass.h>
#include "globdef.h"
#include "fastmath.h"

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int FMA_FastRoundLog2(double x) 
{
  return FMA_FAST_ROUND_LOG2(x);
}

int FMA_FastIntLog2(double x) 
{
  return FMA_FAST_INT_LOG2(x);
}

float FMA_FastApproxLog2( float x )
{
  return FMA_FAST_APPROX_LOG2(x);
}

