/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/
/* Coherent noise function over 1, 2 or 3 dimensions */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "globdef.h"
#include "fastmath.h"
#include "util.h"
#include "pnoise.h"
#include "mtool.h"
#include "rand.h"


#ifdef PNS_DOUBLE_PREC
typedef double pns_grad_t;
typedef LongInt pns_int_t;
#else
typedef float pns_grad_t;
typedef int pns_int_t;
#endif

/*=========================================================================
 *
 *
 *   Improved Noise http://mrl.nyu.edu/~perlin/noise/
 *
 *
 * =======================================================================*/
 
   static pns_int_t pp[512];
   static pns_int_t permutation[] = { 151,160,137,91,90,15,
   131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
   190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
   88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
   77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
   102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
   135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
   5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
   223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
   129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
   251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
   49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
   138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180
   };
   
int improvedNoiseInit()
{ 
  int i;
  for (i=0; i < 256 ; i++) pp[256+i] = pp[i] = permutation[i]; 
  return 0;
}

static  double fade(double t) { return t * t * t * (t * (t * 6 - 15) + 10); }
static double lerp(double t, double a, double b) { return a + t * (b - a); }
static double grad(pns_int_t hash, double x, double y, double z) 
{
      pns_int_t h = hash & 15;                      // CONVERT LO 4 BITS OF HASH CODE
      double u = h<8 ? x : y,                 // INTO 12 GRADIENT DIRECTIONS.
             v = h<4 ? y : h==12||h==14 ? x : z;
      return ((h&1) == 0 ? u : -u) + ((h&2) == 0 ? v : -v);
}

double improvedNoise3( double x, double y, double z)
{
      pns_int_t X = (pns_int_t)(floor(x)) & 255,                  // FIND UNIT CUBE THAT
          Y = (pns_int_t)(floor(y)) & 255,                  // CONTAINS POINT.
          Z = (pns_int_t)(floor(z)) & 255;
      x -= floor(x);                                // FIND RELATIVE X,Y,Z
      y -= floor(y);                                // OF POINT IN CUBE.
      z -= floor(z);
      double u = fade(x),                                // COMPUTE FADE CURVES
             v = fade(y),                                // FOR EACH OF X,Y,Z.
             w = fade(z);

      pns_int_t A = pp[X  ]+Y, AA = pp[A]+Z, AB = pp[A+1]+Z,      // HASH COORDINATES OF
          B = pp[X+1]+Y, BA = pp[B]+Z, BB = pp[B+1]+Z;      // THE 8 CUBE CORNERS,

      return lerp(w, lerp(v, lerp(u, grad(pp[AA  ], x  , y  , z   ),  // AND ADD
                                     grad(pp[BA  ], x-1, y  , z   )), // BLENDED
                             lerp(u, grad(pp[AB  ], x  , y-1, z   ),  // RESULTS
                                     grad(pp[BB  ], x-1, y-1, z   ))),// FROM  8
                     lerp(v, lerp(u, grad(pp[AA+1], x  , y  , z-1 ),  // CORNERS
                                     grad(pp[BA+1], x-1, y  , z-1 )), // OF CUBE
                             lerp(u, grad(pp[AB+1], x  , y-1, z-1 ),
                                     grad(pp[BB+1], x-1, y-1, z-1 ))));
}

/*=========================================================================
 *
 *
 *   Original Perlin Implementation
 *
 *
 * =======================================================================*/


#if 1
#define _B_ 0x1000
#define BM 0xfff
#else
#define _B_ 0x100
#define BM 0xff
#endif

#ifdef PNS_DOUBLE_PREC
  #ifdef MI_WITH_64BIT
  #define N 0x80000000
  #else
  #define N 0x10000000
  #endif
#else
  #define N 0x1000
#endif

#ifdef MI_WITH_64BIT
#define CORR_OFFSET_X 12345678.98765	
#define CORR_OFFSET_Y 98765432.12345
#define CORR_OFFSET_Z 31417549.14133
#else
#define CORR_OFFSET_X 123.98765	
#define CORR_OFFSET_Y 987.12345
#define CORR_OFFSET_Z 314.14133
#endif

#define s_curve(t) ( t * t * (3. - 2. * t) )
#define s_curve2(t) (t * t * t * (t * (t * 6 - 15) + 10))

#define LERP(t, a, b) ( a + t * (b - a) )
#define setup(i,b0,b1,r0,r1)\
        t = vec[i] + N;\
        b0 = ((pns_int_t)t) & BM;\
        b1 = (b0+1) & BM;\
        r0 = t - (pns_int_t)t;\
        r1 = r0 - 1.;
#define at2(rx,ry) ( rx * q[0] + ry * q[1] )
#define at3(rx,ry,rz) ( rx * q[0] + ry * q[1] + rz * q[2] )
#define at4(rx,ry,rz,rw) ( rx * q[0] + ry * q[1] + rz * q[2] + rw * q[3] )

void normalize3(pns_grad_t v[3]);
void normalize2(pns_grad_t v[2]);

static pns_int_t p[_B_ + _B_ + 2];
static pns_grad_t g4[_B_ + _B_ + 2][4];
static pns_grad_t g3[_B_ + _B_ + 2][3];
static pns_grad_t g2[_B_ + _B_ + 2][2];
static pns_grad_t g1[_B_ + _B_ + 2];

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void PNS_noise3_vec3_test(double pos[3], float *out)
{
  int k;

  for(k=0;k<3;k++)
  {
//    double vec[3] = {pos[0]+k*173,pos[1]+k*333,pos[2]+k*912};
    double vec[3] = {pos[0]+k*173.7,pos[1]+k*333.12,pos[2]+k*912.71};
    pns_int_t bx0, bx1, by0, by1, bz0, bz1, b00, b10, b01, b11;
    pns_grad_t rx0, rx1, ry0, ry1, rz0, rz1, sx, sy, sz, a, b, c, d, u, v;
    double t;
    pns_int_t i, j;
    pns_grad_t *q;

    //IACA_START

    setup(0, bx0,bx1, rx0,rx1);
    setup(1, by0,by1, ry0,ry1);
    setup(2, bz0,bz1, rz0,rz1);

    i = p[ bx0 ];
    j = p[ bx1 ];

    b00 = p[ i + by0 ];
    b10 = p[ j + by0 ];
    b01 = p[ i + by1 ];
    b11 = p[ j + by1 ];

    sx  = s_curve(rx0);
    sy = s_curve(ry0);
    sz = s_curve(rz0);

    q = g3[ b00 + bz0 ] ; u = at3(rx0,ry0,rz0);   
    q = g3[ b10 + bz0 ] ; v = at3(rx1,ry0,rz0);
    a = LERP(sx, u, v);
    q = g3[ b01 + bz0 ] ; u = at3(rx0,ry1,rz0);
    q = g3[ b11 + bz0 ] ; v = at3(rx1,ry1,rz0);
    b = LERP(sx, u, v);

    c = LERP(sy, a, b);

    q = g3[ b00 + bz1 ] ; u = at3(rx0,ry0,rz1);
    q = g3[ b10 + bz1 ] ; v = at3(rx1,ry0,rz1);
    a = LERP(sx, u, v);

    q = g3[ b01 + bz1 ] ; u = at3(rx0,ry1,rz1);
    q = g3[ b11 + bz1 ] ; v = at3(rx1,ry1,rz1);
    b = LERP(sx, u, v);

    d = LERP(sy, a, b);

    //IACA_END

    out[k] = LERP(sz, c, d);
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/

double noise1(double arg, int typ)
{
  if(typ==PNS_NTYP_IMPROVED)
  {
    return improvedNoise3(arg,0,0);
  }
  else if(typ==PNS_NTYP_FAST)
  {
   pns_int_t bx0, bx1;
   double rx0, rx1, sx, t, u, v, vec[1];

   vec[0] = arg;

   setup(0,bx0,bx1,rx0,rx1);

   sx = s_curve(rx0);
   u = rx0 * g1[ p[ bx0 ] ];
   v = rx1 * g1[ p[ bx1 ] ];

   return(LERP(sx, u, v));
  }
  else
  {
    fprintf(stderr,"invalid noise type %d\n",typ);
    return 0;
  }
}

double noise2(double vec[2], int typ)
{
  if(typ==PNS_NTYP_IMPROVED)
  {
    return improvedNoise3(vec[0],vec[1],0);
  }
  else if(typ==PNS_NTYP_FAST)
  {
   pns_int_t bx0, bx1, by0, by1, b00, b10, b01, b11;
   double rx0, rx1, ry0, ry1, sx, sy, a, b, t, u, v;
   pns_int_t i, j;

   pns_grad_t *q;

   setup(0, bx0,bx1, rx0,rx1);
   setup(1, by0,by1, ry0,ry1);

   i = p[ bx0 ];
   j = p[ bx1 ];

   b00 = p[ i + by0 ];
   b10 = p[ j + by0 ];
   b01 = p[ i + by1 ];
   b11 = p[ j + by1 ];

   sx = s_curve(rx0);
   sy = s_curve(ry0);

   q = g2[ b00 ] ; u = at2(rx0,ry0);
   q = g2[ b10 ] ; v = at2(rx1,ry0);
   a = LERP(sx, u, v);

   q = g2[ b01 ] ; u = at2(rx0,ry1);
   q = g2[ b11 ] ; v = at2(rx1,ry1);
   b = LERP(sx, u, v);

   return LERP(sy, a, b);
  }
  else
  {
    fprintf(stderr,"invalid noise type %d\n",typ);
    return 0;
  }
}


#if 0

#define gw( q, w, p_ ) \
{\
  pns_grad_t *p=&(p_)[0]; \
  q[0] = w[0]*p[0] + w[1]*p[1] + w[2]*p[2];\
  q[1] = w[3]*p[0] + w[4]*p[1] + w[5]*p[2];\
  q[2] = w[6]*p[0] + w[7]*p[1] + w[8]*p[2];\
}

double noise3_warp(double vec[3], double *W, int typ)
{
   pns_int_t bx0, bx1, by0, by1, bz0, bz1, b00, b10, b01, b11;
   pns_grad_t rx0, rx1, ry0, ry1, rz0, rz1, sx, sy, sz, a, b, c, d, u, v;
   double t;
   pns_int_t i, j;
   pns_grad_t q[3];

   setup(0, bx0,bx1, rx0,rx1);
   setup(1, by0,by1, ry0,ry1);
   setup(2, bz0,bz1, rz0,rz1);

   i = p[ bx0 ];
   j = p[ bx1 ];

   b00 = p[ i + by0 ];
   b10 = p[ j + by0 ];
   b01 = p[ i + by1 ];
   b11 = p[ j + by1 ];

   sx  = s_curve(rx0);
   sy = s_curve(ry0);
   sz = s_curve(rz0);
   
   gw(q,W,g3[ b00 + bz0 ]) ; u = at3(rx0,ry0,rz0);   
   gw(q,W,g3[ b10 + bz0 ]) ; v = at3(rx1,ry0,rz0);
   a = LERP(sx, u, v);
   gw(q,W,g3[ b01 + bz0 ]) ; u = at3(rx0,ry1,rz0);
   gw(q,W,g3[ b11 + bz0 ]) ; v = at3(rx1,ry1,rz0);
   b = LERP(sx, u, v);

   c = LERP(sy, a, b);

   gw(q,W,g3[ b00 + bz1 ]) ; u = at3(rx0,ry0,rz1);
   gw(q,W,g3[ b10 + bz1 ]) ; v = at3(rx1,ry0,rz1);
   a = LERP(sx, u, v);

   gw(q,W,g3[ b01 + bz1 ]) ; u = at3(rx0,ry1,rz1);
   gw(q,W,g3[ b11 + bz1 ]) ; v = at3(rx1,ry1,rz1);
   b = LERP(sx, u, v);

   d = LERP(sy, a, b);
   return LERP(sz, c, d);
}

#endif

double noise3(double vec[3], int typ)
{
/*  if(typ==PNS_NTYP_NONE)  // testing
  {
    return PNS_Noise3DOrigSIMD( vec );
  }
  else*/ 
  if(typ==PNS_NTYP_IMPROVED)
  {
#if 1
    return improvedNoise3(vec[0],vec[1],vec[2]);
#else
#if 1
    float px[8],py[8],pz[8];

    px[0]=vec[0];
    py[0]=vec[1];
    pz[0]=vec[2];

    //TRCP(("%g -> %g\n",vec[0],px[0]));

    float res[8];
#if 1
    PNS_noise3d_simd4_f(px,py,pz,1,
		      res);
#else
    PNS_noise3d_simd8(px,py,pz,1,
		      res);
#endif

//    TRCP(("%g %g %g -> %g\n",vec[0],vec[1],vec[2],res[0]));

    return res[0];
#else
    return PNS_simdtst_noise3_v4(vec);
#endif
#endif
  }
  else if(typ==PNS_NTYP_FAST)
  {
   pns_int_t bx0, bx1, by0, by1, bz0, bz1, b00, b10, b01, b11;
   pns_grad_t rx0, rx1, ry0, ry1, rz0, rz1, sx, sy, sz, a, b, c, d, u, v;
   double t;
   pns_int_t i, j;
   pns_grad_t *q;

   //IACA_START

   setup(0, bx0,bx1, rx0,rx1);
   setup(1, by0,by1, ry0,ry1);
   setup(2, bz0,bz1, rz0,rz1);

   i = p[ bx0 ];
   j = p[ bx1 ];

   b00 = p[ i + by0 ];
   b10 = p[ j + by0 ];
   b01 = p[ i + by1 ];
   b11 = p[ j + by1 ];

   sx  = s_curve(rx0);
   sy = s_curve(ry0);
   sz = s_curve(rz0);
   
   q = g3[ b00 + bz0 ] ; u = at3(rx0,ry0,rz0);   
   q = g3[ b10 + bz0 ] ; v = at3(rx1,ry0,rz0);
   a = LERP(sx, u, v);
   q = g3[ b01 + bz0 ] ; u = at3(rx0,ry1,rz0);
   q = g3[ b11 + bz0 ] ; v = at3(rx1,ry1,rz0);
   b = LERP(sx, u, v);

   c = LERP(sy, a, b);

   q = g3[ b00 + bz1 ] ; u = at3(rx0,ry0,rz1);
   q = g3[ b10 + bz1 ] ; v = at3(rx1,ry0,rz1);
   a = LERP(sx, u, v);

   q = g3[ b01 + bz1 ] ; u = at3(rx0,ry1,rz1);
   q = g3[ b11 + bz1 ] ; v = at3(rx1,ry1,rz1);
   b = LERP(sx, u, v);

   d = LERP(sy, a, b);
   
   //IACA_END
   
   return LERP(sz, c, d);
  }
  else
  {
    fprintf(stderr,"invalid noise type %d\n",typ);
    return 0;
  }

}

double PNS_Noise3D( double* pos )
{
  return noise3(pos,PNS_NTYP_FAST);
}

float noise4(double vec[4], 
    	     int improveSmoothness )
{
   pns_int_t bx0, bx1, by0, by1, bz0, bz1, bw0, bw1;
   pns_grad_t rx0, rx1, ry0, ry1, rz0, rz1, rw0, rw1, 
	  sx, sy, sz, sw, a, b, c, d, e, t, u, v, f;
   pns_int_t i, j;
   pns_grad_t *q;

//   IACA_START
   setup(0, bx0,bx1, rx0,rx1);
   setup(1, by0,by1, ry0,ry1);
   setup(2, bz0,bz1, rz0,rz1);
   setup(3, bw0,bw1, rw0,rw1);

   i = p[ bx0 ];
   j = p[ bx1 ];

   int 
   b00 = p[ i + by0 ],
   b10 = p[ j + by0 ],
   b01 = p[ i + by1 ],
   b11 = p[ j + by1 ],
   b000 = p[ b00 + bz0 ],
   b100 = p[ b10 + bz0 ],
   b010 = p[ b01 + bz0 ],
   b110 = p[ b11 + bz0 ],
   b001 = p[ b00 + bz1 ],
   b101 = p[ b10 + bz1 ],
   b011 = p[ b01 + bz1 ],
   b111 = p[ b11 + bz1 ];

   if(improveSmoothness)
   { // 5th order polynomial as in Perlin "Improving Noise"
     sx  = s_curve2(rx0);
     sy = s_curve2(ry0);
     sz = s_curve2(rz0);
     sw = s_curve2(rw0);
   }
   else
   {
     sx  = s_curve(rx0);
     sy = s_curve(ry0);
     sz = s_curve(rz0);
     sw = s_curve(rw0);
   }
   
   q = g4[ b000 + bw0 ] ; u = at4(rx0,ry0,rz0,rw0);   
   q = g4[ b100 + bw0 ] ; v = at4(rx1,ry0,rz0,rw0);
   a = LERP(sx, u, v);
   q = g4[ b010 + bw0 ] ; u = at4(rx0,ry1,rz0,rw0);
   q = g4[ b110 + bw0 ] ; v = at4(rx1,ry1,rz0,rw0);
   b = LERP(sx, u, v);
   c = LERP(sy, a, b);
   q = g4[ b001 + bw0 ] ; u = at4(rx0,ry0,rz1,rw0);
   q = g4[ b101 + bw0 ] ; v = at4(rx1,ry0,rz1,rw0);
   a = LERP(sx, u, v);
   q = g4[ b011 + bw0 ] ; u = at4(rx0,ry1,rz1,rw0);
   q = g4[ b111 + bw0 ] ; v = at4(rx1,ry1,rz1,rw0);
   b = LERP(sx, u, v);
   d = LERP(sy, a, b);
   e = LERP(sz, c, d);

   q = g4[ b000 + bw1 ] ; u = at4(rx0,ry0,rz0,rw1);   
   q = g4[ b100 + bw1 ] ; v = at4(rx1,ry0,rz0,rw1);
   a = LERP(sx, u, v);
   q = g4[ b010 + bw1 ] ; u = at4(rx0,ry1,rz0,rw1);
   q = g4[ b110 + bw1 ] ; v = at4(rx1,ry1,rz0,rw1);
   b = LERP(sx, u, v);
   c = LERP(sy, a, b);
   q = g4[ b001 + bw1 ] ; u = at4(rx0,ry0,rz1,rw1);
   q = g4[ b101 + bw1 ] ; v = at4(rx1,ry0,rz1,rw1);
   a = LERP(sx, u, v);
   q = g4[ b011 + bw1 ] ; u = at4(rx0,ry1,rz1,rw1);
   q = g4[ b111 + bw1 ] ; v = at4(rx1,ry1,rz1,rw1);
   b = LERP(sx, u, v);
   d = LERP(sy, a, b);
   f = LERP(sz, c, d);

   //IACA_END
   return LERP(sw, e, f);
}

float PNS_Noise4D( double* pos, int improveSmoothness )
{
  return noise4(pos, improveSmoothness);
}

void normalize2(pns_grad_t v[2])
{
   double s;

   s = sqrt(v[0] * v[0] + v[1] * v[1]);
   v[0] = v[0] / s;
   v[1] = v[1] / s;
}

void normalize3(pns_grad_t v[3])
{
   double s;

   s = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
   v[0] = v[0] / s;
   v[1] = v[1] / s;
   v[2] = v[2] / s;
}

void normalize4(pns_grad_t v[4])
{
   double s;

   s = sqrt((double)v[0] * (double)v[0] + 
            (double)v[1] * (double)v[1] + 
	    (double)v[2] * (double)v[2] + 
	    (double)v[3] * (double)v[3]);
   v[0] = v[0] / s;
   v[1] = v[1] / s;
   v[2] = v[2] / s;
   v[3] = v[3] / s;
}

#ifndef MI_WITH_64BIT
#define PNS_RAND_OLD
#endif

int pns_rand_init(int seed, void *rnds)
{
  int rcThis=0;
#ifdef PNS_RAND_OLD
  srand(seed);
  return 0;
#else
  {
    srand(seed);
    int rc = RND_Init((RND_Stream*)rnds,RND_GetGenerator("knuth"),seed);
    if(rc) raiseRc(rc);
  }
#endif
rcCatch:
  return rcThis;
}

int pns_rand(void *rnds)
{
#ifdef PNS_RAND_OLD
  return rand();
#else
  int i = RND_GetInt(rnds,32768),
      j = (uint32_t)(rand()) & (32768-1);
  int k=j*i;
  return k & (32768-1);
#endif
}

int PNS_Init(int rseed)
{
  UtTimer tm;

  TIMER_START(&tm);

   pns_int_t i, j, k;
   RND_Stream rnds;
 
//   TRCP(("RAND_MAX = %ld\n",(long)RAND_MAX));

   pns_rand_init(rseed,&rnds); 
   improvedNoiseInit();

   for (i = 0 ; i < _B_ ; i++) {
      p[i] = i;		// TODO*: use better PRNG than 'rand()'
      g1[i] = (double)((pns_rand(&rnds) % (_B_ + _B_)) - _B_) / _B_;

      for (j = 0 ; j < 2 ; j++)
         g2[i][j] = (double)((pns_rand(&rnds) % (_B_ + _B_)) - _B_) / _B_;
      normalize2(g2[i]);

      for (j = 0 ; j < 3 ; j++)
         g3[i][j] = (double)((pns_rand(&rnds) % (_B_ + _B_)) - _B_) / _B_;
      normalize3(g3[i]);
      //TRCP(("g3[%d] = %g %g %g\n",i,g3[i][0],g3[i][1],g3[i][2]));
   }

   while (--i) {
      k = p[i];
      j = pns_rand(&rnds);
      UT_ASSERT0(j>=0&&j<32768);
      j = j % _B_;
      UT_ASSERT0(j>=0&&j<32768);
      p[i] = p[j];
      p[j] = k;
   }

#if 0
   for (i = 0 ; i < _B_ ; i++) 
   {
     TRCP(("p_[%d]=%d\n",i,p[i]));
   }
#endif

   // do the 4D gradients separately for compatibility of g3 with old seed 
   for (i = 0 ; i < _B_ ; i++) 
   {  
     for (j = 0 ; j < 4 ; j++)
	 g4[i][j] = (double)((pns_rand(&rnds) % (_B_ + _B_)) - _B_) / _B_;
      normalize4(g4[i]);
    }

   // fill wrap around
   for (i = 0 ; i < _B_ + 2 ; i++) {
      p[_B_ + i] = p[i];
      g1[_B_ + i] = g1[i];
      for (j = 0 ; j < 2 ; j++)
         g2[_B_ + i][j] = g2[i][j];
      for (j = 0 ; j < 3 ; j++)
         g3[_B_ + i][j] = g3[i][j];
      for (j = 0 ; j < 4 ; j++)
         g4[_B_ + i][j] = g4[i][j];
   }

   PNS_Init2(rseed);

   PNS_Init3(rseed);

   TIMER_STOP(&tm);

  TRC(("PNS_Init-> CPU = %.0f seconds\n", 
	  (double)(TIMER_DIFF_MS(&tm)/1000.0)));

   return 0;
}

/* --- My harmonic summing functions - PDB --------------------------*/

/*
   In what follows "alpha" is the weight when the sum is formed.
   Typically it is 2, As this approaches 1 the function is noisier.
   "beta" is the harmonic scaling/spacing, typically 2.
*/

double PerlinNoise1D(double x,double alpha,double beta,double octaves,
    		     int typ)
{
   int i,n=(int)octaves;
   double val,rem,sum = 0, alpha_inv=1./alpha;
   double p,scale = 1;

   p = x;
   for (i=0;i<n;i++) {
      val = noise1(p,typ);
      sum += val * scale;
      scale *= alpha_inv;
      p *= beta;
   }
   if((rem = octaves - n) > 1.0E-12)
   {
      val = noise1(p,typ);
      sum += rem * val * scale;
   }
   return(sum);
}

double PNS_PerlinNoise2D(
    	double x,double y,double alpha,double beta,
    	double octaves,  // allow non-integer octaves for seamless LOD blending
	int typ
	)
{
   int i,n=(int)octaves;
   double val,sum = 0, rem, alpha_inv=1./alpha;
   double p[2],scale = 1;

   p[0] = x+CORR_OFFSET_X;
   p[1] = y+CORR_OFFSET_Y;
   for (i=0;i<n;i++) {
      val = noise2(p,typ);
      sum += val * scale;
      scale *= alpha_inv;
      p[0] *= beta;
      p[1] *= beta;
   }

   if((rem = octaves - n) > 1.0E-12)
   {
      val = noise2(p,typ);
      sum += rem * val * scale;
   }

   return(sum);
}

double PNS_PerlinNoise3D(double x,double y,double z,double alpha, double beta,
    		     double octaves,
		     int do_pyroclastic,
		     int typ)
{
   int i,n=(int)octaves;
   float val,sum = 0,rem, alpha_inv=1./alpha,scale=1;
   double p[3];

   p[0] = x;
   p[1] = y;
   p[2] = z;
   for (i=0;i<n;i++) {
      val = noise3(p,typ);
      if(do_pyroclastic) val = fabs(val);
      sum += val * scale;
      scale *= alpha_inv;
      p[0] *= beta;
      p[1] *= beta;
      p[2] *= beta;
   }
   if((rem = octaves - n) > 1.0E-12)
   {
      val = noise3(p,typ);
      if(do_pyroclastic) val = fabs(val);
      sum += rem * val * scale;
   }
   return(sum);
}

double PNS_Fractal3D(double x,double y,double z,double alpha, double beta,
    		     double octaves, 
		     int octaves_leadin,  // see http://www.planetside.co.uk/forums/index.php/topic,17566.msg170871.html#msg170871
		     int turb_typ,  // 0=fbm, 1=abs, 2=ridged
		     float ridge_offs, float ridge_mult, float ridge_exp, 
		     int do_multifract,
		     int ntyp)
{
   octaves += octaves_leadin;
   int i,n=(int)octaves;
   float val,sum = 0,rem, alpha_inv=1./alpha,scale=1,prev=1;
   double p[3];

   if(turb_typ==-2)
   {
     UT_ASSERT0(FALSE);
   }
   
   p[0] = x;
   p[1] = y;
   p[2] = z;
   for (i=0;i<n;i++) {
      val = noise3(p,ntyp); 
      if(turb_typ==1)
      {
	val = fabs(val);
      }
      else if(turb_typ==2)
      {
	UT_ASSERT0(FALSE);
      }
      sum += val * scale * prev;

#if 0
      TRCP(("oct %d of %d+%d -> =%g\n",i,(int)octaves,octaves_leadin)); 
#endif

      if(do_multifract) prev=val;
      if(!(i<octaves_leadin)) scale *= alpha_inv;
      p[0] *= beta;
      p[1] *= beta;
      p[2] *= beta;
   }
   if((rem = octaves - n) > 1.0E-12)
   {
      val = noise3(p,ntyp); 
      if(turb_typ==1)
      {
	val = fabs(val);
      }
      else if(turb_typ==2)
      {
	UT_ASSERT0(FALSE);
      }
      sum += rem * val * scale;
#if 0
      TRCP(("oct %d of %d+%d -> =%g\n",i,(int)octaves,octaves_leadin)); 
#endif

   }
   return(sum);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
double PNS_Fractal(double *pos, int dim,  // position, dimension
    		     double alpha, double beta,
    		     double octaves, 
		     int octaves_leadin,  // see http://www.planetside.co.uk/forums/index.php/topic,17566.msg170871.html#msg170871
		     int turb_typ,  // 0=fbm, 1=abs, 2=ridged
		     float ridge_offs, float ridge_mult, float ridge_exp, 
		     int do_multifract,
		     //double pa_mem, double pa_ampl, // 'structure memory' & amplitude for flow noise 'pseudo advection' see Ebert. p 384 ff
		     int ntyp,
		     double seed, double seed_tm  // seed (for spatial & time dimension)
    			)
{
   octaves += octaves_leadin;
   int i,n=(int)octaves,k;
   float val,sum = 0,rem, alpha_inv=1./alpha,scale=1,prev=1;
   double p[4];

   if(!(dim>=2&&dim<=4))
   {
     UT_ASSERT_NR(FALSE);
     return 0;
   }

   p[0] = pos[0] + seed;
   p[1] = pos[1] + seed*0.7;
   if(dim>=3) p[2] = pos[2];
   if(dim>=4) p[3] = pos[3] + seed_tm*3.141;

   if(turb_typ==-2)
   {
     UT_ASSERT0(FALSE);
   }

   for (i=0;i<n;i++) {
      
     val = dim==2 ? noise2(p,ntyp) :
           dim==3 ? noise3(p,ntyp) :
           dim==4 ? PNS_Noise4D(p,0) : 0;

      if(turb_typ==1)
      {
	val = fabs(val);
      }
      else if(turb_typ==2)
      {
	UT_ASSERT0(FALSE);
      }
      sum += val * scale * prev;

      if(do_multifract) prev=val;
      if(!(i<octaves_leadin)) scale *= alpha_inv;
      for(k=0;k<dim;k++) p[k] *= beta;
   }
   if((rem = octaves - n) > 1.0E-12)
   {
     val = dim==2 ? noise2(p,ntyp) :
           dim==3 ? noise3(p,ntyp) :
           dim==4 ? PNS_Noise4D(p,0) : 0;
      if(turb_typ==1)
      {
	val = fabs(val);
      }
      else if(turb_typ==2)
      {
	UT_ASSERT0(FALSE);
      }
      sum += rem * val * scale;
   }
   return(sum);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
double PNS_Fractal_2( 
    	      int dim,
	      double *pos, 
	      int do_comp_spectrum,	// opmode
	      int do_comp_valsum,	      
	      // Optional IN/OUT of <n_octaves> noise values for sepctral   
	      // texture blending as of [Neyret2003]
	      PNS_Spectrum *spectrum,   
	      double octaves, 
	      int octaves_leadin,  // see http://www.planetside.co.uk/forums/index.php/topic,17566.msg170871.html#msg170871
              double lacunarity,
	      // 
	      double alpha,
	      int turb_typ,  // 0=fbm, 1=abs, 2=ridged
	      float ridge_offs, float ridge_mult, float ridge_exp
    		)
{
   int i,k;
   float val,sum = 0,rem, alpha_inv=1./alpha,scale=1;
   double p[4];

   if(!(dim>=2&&dim<=4))
   {
     UT_ASSERT_NR(FALSE);
     return 0;
   }

   if(do_comp_valsum && turb_typ==-2)
   {
     UT_ASSERT0(FALSE);
   }

   if(pos) 
   {
     for(k=0;k<dim;k++) p[k]=pos[k];
   }
   else
   {
     UT_ASSERT_FATAL(spectrum);
   }

   if(do_comp_spectrum)
   {
     if(spectrum)
     {
       spectrum->octaves = octaves;
       spectrum->octaves_leadin = octaves_leadin;
     }
   }
   else
   {
     octaves = spectrum->octaves;
     octaves_leadin = spectrum->octaves_leadin;
   }

   octaves += octaves_leadin;

   int n=(int)octaves;

   if((rem = octaves - n) > 1.0E-12) n++;

   if(spectrum)
   {
     if(do_comp_spectrum)
     {
       spectrum->n_noisevals = n;
     }
     else
     {
       UT_ASSERT0(spectrum->n_noisevals==n);
     }
   }

   for (i=0;i<n;i++) 
   { 
     if(do_comp_spectrum)
     {
       val = dim==2 ? noise2(p,PNS_NTYP_FAST) :
             dim==3 ? noise3(p,PNS_NTYP_FAST) :
             dim==4 ? PNS_Noise4D(p,0) : 0;
       if(spectrum) spectrum->noisevals[i] = val;
       for(k=0;k<dim;k++) p[k] *= lacunarity;
     }
     else
     {
       val = spectrum->noisevals[i];
     }

     if(!do_comp_valsum) continue;

    if(turb_typ==1)
    {
      val = fabs(val);
    }
    else if(turb_typ==2)
    {
      UT_ASSERT0(FALSE);
    }

    if(i>=(int)octaves) val*=rem;
    sum += val * scale;
    if(!(i<octaves_leadin)) scale *= alpha_inv;
   }

   return(sum);
}

