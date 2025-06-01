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
#include <stdlib.h>
#include <tbb/tbb.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <errno.h>
#include <unistd.h>

#include "msbg.h"
#include "msbg_demo.h"

//#include <openvdb/openvdb.h>

int nMaxThreads=-1;

int main(int argc, char **argv) 
{
  int multiProc=0,
      blockSize0=16,
      testCase=0;

  //
  // Open log file 
  //
  LogLevel = 2; 
  {
    char *logFile = (char*)"msbg.log";
    int rc = UtOpenLogFile(logFile,FALSE,FALSE);
    if(rc)
    {
      fprintf(stderr,"can't open log file '%s'\n",logFile);
      exit(1);
    }

    {
      time_t rawtime;
      struct tm * timeinfo;
      time ( &rawtime );
      timeinfo = localtime ( &rawtime );
      TRC1(("local time: %s",asctime (timeinfo)));
    }
  }

  //
  // Parse command line parameters
  //
  char filename[UT_MAXPATHLEN];
  strcpy(filename,"c:/tmp/msbgtest.vdb");
  int resolution = 768;

  extern float camPos[3];
  extern float camLookAt[3];
  extern float camLight[3];
  extern float camZoom;
  extern int camRes[2];

  {
    char c;
    while((c = getopt(argc,argv,"l:i:a:b:u:r:f:c:jv:")) != EOF)
    {
      switch(c)
      {
        case 'a':
	   {
	     UT::StrTok t(optarg,","); 
	     if(t.n()!=7)
	     {
	       TRCERR(("7 sub-parameters expected for option '-a'\n"));
	       exit(1);
	     }
	     for(int i=0;i<3;i++) camPos[i] = t.getFloat(i);
	     for(int i=0;i<3;i++) camLookAt[i] = t.getFloat(3+i);
	     camZoom = t.getFloat(6);
	   }
           break;
        case 'l':
	   {
	     UT::StrTok t(optarg,","); 
	     if(t.n()!=3)
	     {
	       TRCERR(("3 sub-parameters expected for option '-3'\n"));
	       exit(1);
	     }
	     for(int i=0;i<3;i++) camLight[i] = round(t.getFloat(i));
	   }
           break;
        case 'i':
	   {
	     UT::StrTok t(optarg,","); 
	     if(t.n()!=2)
	     {
	       TRCERR(("2 sub-parameters expected for option '-i'\n"));
	       exit(1);
	     }
	     for(int i=0;i<2;i++) camRes[i] = round(t.getFloat(i));
	   }
           break;

        case 'u':
           multiProc = atoi( optarg );
	   nMaxThreads = multiProc;
           break;
	case 'f':
	  strcpy(filename,optarg);
	  break;
	case 'b':
	  blockSize0 = atoi( optarg );
          UT_ASSERT0( blockSize0==16 || blockSize0==32 );
	  break;
	case 'c':
	  testCase = atoi( optarg );
	  break;
	case 'r':
	  resolution = atoi( optarg );
	  break;
	case 'v':
	  LogLevel = atoi( optarg );
	  break;
      }   
    }
  }

  //
  // Initialize multi threading
  //
  {
    int nMaxOmpThreadsAct,
	nMaxTBBThreadsAct,
	nMaxThreadsAct;

    ThrInit();

    ThrGetMaxNumberOfThreads( &nMaxOmpThreadsAct, 
			      &nMaxTBBThreadsAct );
    nMaxThreadsAct = MAX( nMaxOmpThreadsAct, nMaxTBBThreadsAct );

    if(nMaxThreads>0)
    {
      TRCP(("Setting max. TBB threads %d\n",nMaxThreads));
      ThrSetMaxNumberOfTBBThreads( nMaxThreads );
    }

    if(nMaxThreads<=0)
    {
      nMaxThreads = nMaxThreadsAct;
    }
    else
    {
      nMaxThreads = MIN( nMaxThreadsAct, nMaxThreads );
    }

    UT_ASSERT0(nMaxThreads<=THR_MAX_THREADS);
    
    if(nMaxThreads<=0) nMaxThreads=nMaxThreadsAct;

    TRCP(("Using %d of maximum number of threads OMP=%d TBB=%d \n",nMaxThreads,
	  nMaxOmpThreadsAct, nMaxTBBThreadsAct ));

    //TRCP(("omp_proc_bind() = %d\n",omp_get_proc_bind()));

    UT_ASSERT0(nMaxThreads<=MI_MAX_THREADS); // TODO tmpBlkDataFBuf etc.
  }


  {
    int rc =PNS_Init(1234);
    if(rc) 
    {
      TRCERR(("PNS_Init() failed %d\n",rc));
      exit(1);
    }
  }

  //
  //
  //

  resolution = ALIGN(resolution,blockSize0);
  int sx = resolution, sy = sx, sz = sx;

  TRCP(("opcode=%d fname=%s %dx%dx%d\n",testCase,filename,sx,sy,sz)); 

  const char *basePointsFile = testCase==2 ? 
    				  "../data/bun_zipper.ply" :
				  "../data/bun_zipper_res2.ply";
  msbg_test_sparse(testCase, basePointsFile, blockSize0, sx, sy, sz);
  
  return 0;
}

