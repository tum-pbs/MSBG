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
#include <omp.h>
#include <tbb/tbb.h>
//#include <tbb/global_control.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <vectorclass/vectorclass.h>
#include "globdef.h"
#include "mtool.h"
#include "fastmath.h"
#include "util.h"
#include "thread.h"



int thrMaxAllowedParallelism,
    thrIsRunningParallel;

static LONG volatile thrIdxTaskNext=0,
     	             thrNumTasks,
		     thrNumTasks0;



#define THR_CHUNK_SIZE 10
#if 0
static int thrNumThreads;
static ThrTaskFunc thrTaskFunc;
static ThrThreadAdm thrAdmTable[THR_MAX_THREADS];

/*-------------------------------------------------------------------------*/
/* 								           */
/*-------------------------------------------------------------------------*/
static long unsigned int ThrThreadLoopFunc( void *arg )
{
  int tid = *(int*)arg;
  UT_ASSERT0(tid>=0&&tid<THR_MAX_THREADS);
  ThrThreadAdm &thr = thrAdmTable[tid];

  for(;;)
  {
    thr.evtStart.waitReset();
  
    if(thrTaskFunc)
    {
      int iTask;
      while((iTask = 
	 InterlockedIncrement((LONG volatile *)&thrIdxTaskNext) - 1) < thrNumTasks)      
      {
	#ifndef THR_CHUNK_SIZE
	thrTaskFunc( tid, iTask );
	#else
	for(int j=0;j<THR_CHUNK_SIZE;j++)
	{
	  int iTaskAct = iTask*THR_CHUNK_SIZE+j;
	  UT_ASSERT2( iTaskAct>=0 && iTask<thrNumTasks0);	   
	  thrTaskFunc( tid, iTaskAct );
	}
	if(iTask==thrNumTasks-1)
	{
	  for(int j=thrNumTasks*THR_CHUNK_SIZE;j<thrNumTasks0;j++)
	  {
	    int iTaskAct = j;
	    UT_ASSERT2( iTaskAct>=0 && iTask<thrNumTasks0);	   
	    thrTaskFunc( tid, iTaskAct );
	  }
	}
	#endif
      }
    }
    
    thr.evtFinish.set();        
  }
};

/*---------------------------------------------------------------------*/
/* 								       */
/*---------------------------------------------------------------------*/
static void ThrWaitThreads()
{
  for(int i=1;i<thrNumThreads;i++)
  {
    ThrThreadAdm &thr = thrAdmTable[i];
    thr.evtFinish.waitReset();
  }
}

/*---------------------------------------------------------------------*/
/* 								       */
/*---------------------------------------------------------------------*/
static void ThrStartThreads()
{
  for(int i=1;i<thrNumThreads;i++)
  {
    ThrThreadAdm &thr = thrAdmTable[i];
    thr.evtStart.set();
  }
}

/*---------------------------------------------------------------------*/
/* 								       */
/*---------------------------------------------------------------------*/
void ThrInitThreadPool( int nThreads ) 
{
  if(thrNumThreads) return;
  TRCP(("*** ThrInitThreadPool %d\n",nThreads));
  thrNumThreads = nThreads;
  UT_ASSERT0( sizeof( thrAdmTable[0]) >= CPU_CACHELINE_SZ );

  for(int tid=0; tid<thrNumThreads;tid++)
  {
    ThrThreadAdm &thr = thrAdmTable[tid];
    thr.tid = tid;
    
    if(tid==0) continue; // reserve 0 for main thread

    thr.evtStart.init(false);
    thr.evtFinish.init(false);
    
    UT_ASSERT0(false);  // TODO replace CreateThread by _beginThreadEx
    thr.hdl = 
      CreateThread(nullptr, 0, &ThrThreadLoopFunc, (void*)&thr.tid, 0, nullptr);
  }
}

/*---------------------------------------------------------------------*/
/* 								       */
/*---------------------------------------------------------------------*/
void ThrRunParallelNew( int nTasks,
    	std::function<void( int tid, int iTask )> func, void *dummy ) 
{
  thrNumTasks0 = nTasks;
#ifndef THR_CHUNK_SIZE
  thrNumTasks = nTasks;
#else
  thrNumTasks = nTasks/THR_CHUNK_SIZE;
#endif

  thrIdxTaskNext = 0;

  thrTaskFunc = func;

  UT_ASSERT0( sizeof( thrAdmTable[0]) >= CPU_CACHELINE_SZ );

  ThrStartThreads();

  ThrWaitThreads();
}
#endif

#if 0
static unsigned long __stdcall ThrThreadFunc( void *vdata )
{
  return 0;
}
#endif

namespace THR
{

void getMaxNumberOfThreads( int *pNumThreadsOMP,
    			      int *pNumThreadsTBB )
{
  if(pNumThreadsOMP) *pNumThreadsOMP = omp_get_max_threads();
//  *pNumThreadsTBB = tbb::task_scheduler_init::default_num_threads();
  if(pNumThreadsTBB) *pNumThreadsTBB = tbb::this_task_arena::max_concurrency();
}

} // namespace THR


void ThrInit(void)
{
  thrMaxAllowedParallelism = 0;
  thrIsRunningParallel = 0;
}


void ThrSetMaxNumberOfTBBThreads( int nMaxThreads )
{
 thrMaxAllowedParallelism=nMaxThreads;
 static tbb::global_control c(tbb::global_control::max_allowed_parallelism, nMaxThreads);
 //static tbb::task_scheduler_init init(1);
}

void ThrGetMaxNumberOfThreads( int *pNumThreadsOMP,
    			       int *pNumThreadsTBB )
{
  return THR::getMaxNumberOfThreads( pNumThreadsOMP, pNumThreadsTBB );
}

