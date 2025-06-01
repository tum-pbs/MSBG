#!/bin/bash

make -f ../makefile \
  OBJE=o\
  CPP_FLAGS="-std=gnu++17"\
  CFLAGS_OPT="-O3 -DNDEBUG -mtune=native -march=native" \
  CFLAGS_PROF= \
  CFLAGS_OMP=-fopenmp \
  CFLAGS_TBB="" \
  LIB_OMP="-lgomp"\
  CC="gcc" \
  LD="g++  -o" \
  LDFLAGS_BW=\
  AR="ar" \
  ISPC="ispc"\
  CFLAGS_BW="-m64 -DMI_WITH_64BIT"\
  -j \
  $1

