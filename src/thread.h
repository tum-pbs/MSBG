/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#ifndef THREAD_H
#define THREAD_H

using namespace MSBG_NAMESPACE;

//#define THR_STATISTICS

extern int nMaxThreads;
extern int thrMaxAllowedParallelism;
extern int thrIsRunningParallel;

#define THR_MAX_THREADS	128

#ifdef __cplusplus
#include <memory>
#include <atomic>
#include <thread>
#include <tbb/tbb.h>
#include <functional>
//#include <thread>
//
//========================================================================
//
//
//  	ThrThreadLocals
//
//
//========================================================================

template <typename Type >
class ThrThreadLocals
{
public:

  //! Default constructor
  ThrThreadLocals( 
      std::function<void( Type *tls, int idx )> initFunc={},
      std::function<void( Type *tls, int idx )> reduceFunc={}
      ) 
    :  _reduceFunc(reduceFunc)
  {
    for(int i=0;i<THR_MAX_THREADS;i++) 
    {
      Type *tls = &_locals[i].slot;
      *tls = {};  // Initialize members to zero 
      //memset(tls,0,sizeof(Type));
      if(initFunc) initFunc( tls, i);
    }
  }

  //! Destructor
  ~ThrThreadLocals()
  {
  };

  void reduce()
  {
    UT_ASSERT0(_reduceFunc);
    for(int i=0;i<THR_MAX_THREADS;i++) 
    {
      _reduceFunc( &_locals[i].slot, i);
    }
  };

  int size() const { return THR_MAX_THREADS; }

  Type *operator[](int idx) 
  {
    UT_ASSERT2(idx>=0 && idx<THR_MAX_THREADS);
    return &_locals[idx].slot;
  }

private:

  std::function<void( Type *tls, int idx )> _reduceFunc;

  struct // align to avoid false cache sharing
    __attribute__((aligned(CPU_CACHELINE_SZ)))
  {
    Type slot;
  }
  _locals[THR_MAX_THREADS];
};

#if 0
//========================================================================
//
//
// 	ThrEvent
//
//
//========================================================================
class ThrEvent
{
  public:

  /*---------------------------------------------------------------------*/ 
  /* 								         */
  /*---------------------------------------------------------------------*/
  ThrEvent( ) : _valid(false) {}

  /*-----------------------------------------------------------------------*/
  /* 								           */
  /*-----------------------------------------------------------------------*/
  void init(bool state)
  { 
    _hdl = CreateEvent(NULL,	/* default security attribs */
			  FALSE,	/* no manual-reset */
			  state ? TRUE : FALSE,	/* initial state*/
			  NULL);	/* unnamed */
    if(_hdl == NULL)
    {
      UT_ASSERT0(false);
    }
    else
    {
      _valid=TRUE;
    }
  }

  /*---------------------------------------------------------------------*/ 
  /* 								         */
  /*---------------------------------------------------------------------*/
  ~ThrEvent()
  {
    if(_valid)
    {
      //TRCP(("~ThrEvent '%p'\n",_hdl));
      BOOL ok = CloseHandle( _hdl );
      if(!ok)
      {
	TRCERR(("CloseHandle() failed (errno=%d)\n",
	      (int)GetLastError()));
      }
    }
  }

  /*---------------------------------------------------------------------*/ 
  /* 								         */
  /*---------------------------------------------------------------------*/
  void set( void )
  {  
    BOOL ok = SetEvent( _hdl );
    if(!ok)
    {
      TRCERR(("SetEvent() failed (errno=%d)\n",
	  (int)GetLastError()));
    }
  }

  /*---------------------------------------------------------------------*/ 
  /* 								         */
  /*---------------------------------------------------------------------*/
  int waitReset( int timeoutMs=-1 )
  {
    if(timeoutMs==1) timeoutMs = INFINITE;
      
    DWORD rc = WaitForSingleObject( _hdl, timeoutMs);
    if(rc==WAIT_OBJECT_0)
    {
      return (MI_OK);
    }
    else if(rc==WAIT_TIMEOUT)
    {
      return (MI_ETIMEOUT);
    }
    else
    {
      TRCERR(("WaitForSingleObject() failed (errno=%d)\n",
	  (int)GetLastError()));
      UT_ASSERT0(FALSE);
    }
  }

  private:

  int _valid;
  HANDLE  _hdl;
};

//========================================================================
//
//
// 	ThrThreadPool
//
//
//========================================================================
typedef std::function<void(int iThread, int iTask)> ThrTaskFunc;

typedef struct // align to avoid false cache sharing
  __attribute__((aligned(CPU_CACHELINE_SZ)))
{
  int tid;
  HANDLE hdl;
  ThrEvent evtStart,
	   evtFinish;
}
ThrThreadAdm;

void ThrInitThreadPool( int nThreads = THR_MAX_THREADS );
#endif

void ThrRunParallelNew( int nTasks,
    	std::function<void( int tid, int iTask )> func, void *dummy=nullptr ); 

//========================================================================
//
//
// 	ThrRunParallel
//
//
//========================================================================
template < typename TLS, typename bodyFuncType, int doUseBlocklist >
void ThrIRunParallel( size_t _nTasks, 
    		std::vector<int> *blockList,
    		std::function<void( TLS &tls, int tid )> _initFunc,
    		const bodyFuncType &_bodyFunc,  // avoid heavyweight std::function
		std::function<void( TLS &tls, int tid )> _reduceFunc={},
		int doInterlocked = FALSE,
		int doSerial = FALSE )
{
  struct // align to avoid false cache sharing
    __attribute__((aligned(CPU_CACHELINE_SZ)))
  {
    TLS slot;
  }
  _locals[THR_MAX_THREADS] = {};

  //UT_ASSERT0(!doSerial);

  //if(_reduceFunc||_initFunc)
  {
    for(int i=0;i<THR_MAX_THREADS;i++) 
    {
      TLS &tls = _locals[i].slot;
      //tls = {};  // Initialize members to zero 
      if(_initFunc) _initFunc( tls, i);
    }
  }

  /*if(doSerial && !(doSerial==2))
  {
    TRCWARN(("%s: doSerial = TRUE in parallel loop\n",UT_FUNCNAME));
  }*/

  if(thrMaxAllowedParallelism==1) doSerial = TRUE;

  size_t nTasksAct = doUseBlocklist ? blockList->size() :
    		    		      _nTasks;

  if(doSerial)
  {
    for(size_t iTask=0;iTask<nTasksAct;iTask++)
    {
      int tid = 0;
      TLS &tls = _locals[tid].slot;
      size_t iTaskAct;
      if constexpr ( doUseBlocklist )
      {
        iTaskAct = (*blockList)[iTask];
      }
      else
      {
        iTaskAct = iTask;
      }
      _bodyFunc( tls, tid, iTaskAct );
    }
  }
  else
  {
    #if 0
    UtGetAverageCPULoad(NULL);  // call once at begin of interval 
    #endif

    thrIsRunningParallel = TRUE;

    tbb::parallel_for( tbb::blocked_range<size_t>(0,nTasksAct), 
    [&](const tbb::blocked_range<size_t>& range) 
    {
      int tid = tbb::this_task_arena::current_thread_index();
      //TRCP(("tid=%d -> %d-%d/%d\n",tid,range.begin(),range.end(),_nTasks));
      TLS &tls = _locals[tid].slot;
      for(size_t i=range.begin(); i !=range.end(); ++i)
      {
        size_t j;
	if constexpr ( doUseBlocklist )
	{
	  j = (*blockList)[i];
	}
	else
	{
	  j = i;
	}

	_bodyFunc( tls, tid, j );
      }
    } );

    thrIsRunningParallel = FALSE;

    #if 0 
    double deltaTimeMs;
    float avgCpuLoad = UtGetAverageCPULoad(&deltaTimeMs);
    if(!(avgCpuLoad<0))
    {
      TRC(("%s: Avg. CPU load was %.2fs x %g%%\n", UT_FUNCNAME, 
	    deltaTimeMs/1000.0, avgCpuLoad));
    }
    #endif
  }


  if(_reduceFunc)
  {
    for(int i=0;i<THR_MAX_THREADS;i++) 
    {
      _reduceFunc( _locals[i].slot, i);
    }
  }
}

template < typename TLS, typename bodyFuncType >
void ThrRunParallel( 
                std::vector<int> *blocklist, 
    		std::function<void( TLS &tls, int tid )> _initFunc,
    		const bodyFuncType &_bodyFunc,  // avoid heavyweight std::function
		std::function<void( TLS &tls, int tid )> _reduceFunc={},
		int doInterlocked = FALSE,
		int doSerial = FALSE )
{
  return ThrIRunParallel<TLS, bodyFuncType, 1>( 0, blocklist, _initFunc, _bodyFunc, _reduceFunc,
      			  doInterlocked, doSerial );
}

template < typename TLS, typename bodyFuncType >
void ThrRunParallel( size_t _nTasks, 
    		std::function<void( TLS &tls, int tid )> _initFunc,
    		const bodyFuncType &_bodyFunc,  // avoid heavyweight std::function
		std::function<void( TLS &tls, int tid )> _reduceFunc={},
		int doInterlocked = FALSE,
		int doSerial = FALSE )
{
  return ThrIRunParallel<TLS, bodyFuncType, 0>( _nTasks, NULL, _initFunc, _bodyFunc, _reduceFunc,
      			  doInterlocked, doSerial );
}


//========================================================================
//
//
// 	T e s t i n g
//
//
//========================================================================

// https://programmer.help/blogs/c-binds-lambda-expressions-to-callback-functions.html

typedef struct
{
  void *arg;
  int tid;
}
start_routine_arg_t;

static volatile LONG nThreadsActive;

template <typename Fx> long unsigned int start_routine_t(void* arg)
{
    start_routine_arg_t *sra = (start_routine_arg_t*)arg;
    Fx* f = (Fx*)sra->arg;
    int tid = sra->tid;
    (void)(*f)(tid);

    InterlockedDecrement( (LONG volatile *)&nThreadsActive );
    return 0;
}

#if 0
template<class Function>
inline
void ThrRunParallel_(Function f, int n = THR_MAX_THREADS)
{
  HANDLE threads[THR_MAX_THREADS];
  start_routine_arg_t args[THR_MAX_THREADS];

  nThreadsActive = n;
    
  for(int tid=0;tid<n;tid++)
  {
    start_routine_arg_t *sra = &args[tid];
    sra->arg = &f;
    sra->tid = tid;
//    threads[tid] = CreateThread(nullptr, 0, start_routine_t<decltype(f)>, sra, 0, nullptr);
//    pthread_create(&threads[tid], nullptr, start_routine_t<decltype(f)>, sra);
    //TRCP(("tid=%d -> h=%p\n",tid,(void*)threads[tid]));
  }

  while( nThreadsActive ) 
  {
    _mm_pause();
  }

/*  for(int tid=0;tid<n;tid++)
  {
    pthread_join(threads[tid], nullptr);
  }*/
}
#endif

#if 0
//========================================================================
//
//
// 	ThrRunParallel
//
//
//========================================================================

// https://programmer.help/blogs/c-binds-lambda-expressions-to-callback-functions.html

typedef struct
{
  void *arg;
  int tid;
}
start_routine_arg_t;

template <typename Fx> void *start_routine_t(void* arg)
{
    start_routine_arg_t *sra = (start_routine_arg_t*)arg;
    Fx* f = (Fx*)sra->arg;
    int tid = sra->tid;
    (void)(*f)(tid);
    return nullptr;
}

template<class Function>
inline
void ThrRunParallel(Function f, int n = THR_MAX_THREADS)
{
  pthread_t threads[THR_MAX_THREADS];
  start_routine_arg_t args[THR_MAX_THREADS];

  for(int tid=0;tid<n;tid++)
  {
    start_routine_arg_t *sra = &args[tid];
    sra->arg = &f;
    sra->tid = tid;
    pthread_create(&threads[tid], nullptr, start_routine_t<decltype(f)>, sra);
  }

  for(int tid=0;tid<n;tid++)
  {
    pthread_join(threads[tid], nullptr);
  }

}

#if 0
template<class Function>
inline
void parallel_fun(int n, Function f)
{
    std::thread threads[THR_MAX_THREADS];
    for (int i=0; i<n; ++i)
    {
      threads[i] = std::thread(f);
      i++;
    }

    for (int i=0; i<n;i++)
    {
      if(threads[i].joinable())
      threads[i].join();
    }
}
#endif

template< typename Func >
void ThrRunParallel( int nThreads, const Func& func )
{
  if(nThreads==-1)
  {
    nThreads = THR_MAX_THREADS;
  }

  HANDLE th[THR_MAX_THREADS];

  for( int i=0; i<nThreads; i++ )
  {
    int arg=i;
    func( i );
    //th[i] = CreateThread(nullptr, 0, func, &arg, 0, nullptr);
  }
}
#endif

//========================================================================
//
//
//
//
//========================================================================

inline float ThrAtomicMax(float* pf, float f)
{
    union
    {
	volatile LONG iOld;
	float fOld;
    };
    union
    {
	LONG iNew;
	float fNew;
    };

    for(;;)
    {
	iOld = *(int32_t*)pf;   // current old value
	fNew = std::max(fOld,f);
	if(InterlockedCompareExchange( (volatile LONG*)pf, iNew, iOld) == iOld)
	    return fNew; 
    }
}

inline float ThrAtomicMin(float* pf, float f)
{
    union
    {
	volatile LONG iOld;
	float fOld;
    };
    union
    {
	LONG iNew;
	float fNew;
    };

    for(;;)
    {
	iOld = *(int32_t*)pf;   // current old value
	fNew = std::min(fOld,f);
	if(InterlockedCompareExchange( (volatile LONG*)pf, iNew, iOld) == iOld)
	    return fNew; 
    }
}

inline uint16_t ThrAtomicMinU16(uint16_t* pf, uint16_t f)
{
  union
  {
    volatile int16_t iOld;
    uint16_t fOld;
  };
  union
  {
    int16_t iNew;
    uint16_t fNew;
  };


  for(;;)
  {
    iOld = *(int16_t*)pf;   // current old value
    fNew = std::min(fOld,f);
    if(InterlockedCompareExchange16( (volatile int16_t*)pf, iNew, iOld) == iOld)
	return fNew; 
  }
}


inline float ThrAtomicAdd(float* pf, float f)
{
    union
    {
	volatile LONG iOld;
	float fOld;
    };
    union
    {
	LONG iNew;
	float fNew;
    };

    for(;;)
    {
	iOld = *(int32_t*)pf;   // current old value
	fNew = fOld + f;
	if(InterlockedCompareExchange( (volatile LONG*)pf, iNew, iOld) == iOld)
	    return fNew; 
    }
}


inline uint16_t ThrAtomicAddU16(uint16_t* pf, uint16_t f)
{
    union
    {
	volatile int16_t iOld;
	uint16_t fOld;
    };
    union
    {
	int16_t iNew;
	uint16_t fNew;
    };

    for(;;)
    {
	iOld = *(int16_t*)pf;   // current old value
	fNew = MIN((uint32_t)fOld + (uint32_t)f, UINT16_MAX);
	if(InterlockedCompareExchange16( (volatile int16_t*)pf, iNew, iOld) == iOld)
	    return fNew; 
    }
}



inline float ThrAtomicMulAdd(float* pf, float a, float c)  // *pf = (*pf) * a + c
{
    union
    {
	volatile LONG iOld;
	float fOld;
    };
    union
    {
	LONG iNew;
	float fNew;
    };

    for(;;)
    {
	iOld = *(int32_t*)pf;   // current old value
	fNew = a*fOld + c;
	if(InterlockedCompareExchange( (volatile LONG*)pf, iNew, iOld) == iOld)
	    return fNew; 
    }
}


struct qtBarrier
	{
	 volatile LONG c1;
	 volatile LONG c2;
	 qtBarrier() { c1 = 0; c2 = 0; }

	 int enter(LONG iThread, LONG nThreads)
	 {
	  if(iThread)
	  {
	   // indicate this thread reached barrier
	   InterlockedIncrement((volatile LONG*)&c1);
	   while(c1)  _mm_pause();
	   InterlockedIncrement((volatile LONG*)&c2);
	   while(c2) _mm_pause();
	   return 1;
	  }
	  else
	  {
	   while((c1+1)!=nThreads) _mm_pause();
	   return 0;
	  }
	 }

	 void release(int nThreads)
	 {
	   c1 = 0;
	   while((c2+1)!=nThreads) _mm_pause();
	   c2 = 0;
	 }
	};



#if 0
class ThrBarrier
   {
   public:
       inline ThrBarrier(uint32_t numThreads) : m_threadCount(numThreads)
       {
           m_allThreadsReachedSync = 0;

	   InitializeCriticalSection( &m_numThreadsReachedSyncLock );
	   SetCriticalSectionSpinCount( &m_numThreadsReachedSyncLock, 100000 );

           m_numThreadsReachedSync = 0;
           m_syncNumber = 0;
	   m_taCount = 0;
       }
   
       ~ThrBarrier() { 
	    DeleteCriticalSection( &m_numThreadsReachedSyncLock );       
       }

       int  enter( int tid=-1 )
       {
	   int taCount;
           EnterCriticalSection( &m_numThreadsReachedSyncLock);
	   taCount = m_taCount++;

           uint32_t syncNr = m_syncNumber;

           m_numThreadsReachedSync++;

	   TRC1(("BARTST#%d %d/%d sn=%d n=%d a=%d\n", m_taCount, tid, m_threadCount, m_syncNumber, 
		 m_numThreadsReachedSync, m_allThreadsReachedSync ));

           if (m_numThreadsReachedSync == m_threadCount)
           {
	     return 0;
           }
           else
           {
               if ( (m_syncNumber == syncNr) && (m_numThreadsReachedSync == 1) )
               {
		 m_allThreadsReachedSync = 0;
                 //  ResetEvent(m_allThreadsReachedSync);
               }

               while (m_syncNumber == syncNr)
               {
		   int syncNumber = m_syncNumber;
                   LeaveCriticalSection( &m_numThreadsReachedSyncLock);

		   {
		     for(;;)
		     {
       			

			int allThreadsReachedSync = // Atomic read with membarrier
			  InterlockedCompareExchange( &m_allThreadsReachedSync, 0,0);

			TRC1(("BARTST#%d %d/%d W sn=%d/%d a=%d\n", taCount, tid, m_threadCount, syncNumber, syncNr, m_allThreadsReachedSync ));

			if(allThreadsReachedSync) break;
			_mm_pause();			
			Sleep(0);
			//YieldProcessor();
		     }
		   }

                   EnterCriticalSection( &m_numThreadsReachedSyncLock);
       	           taCount = m_taCount++;
        }
               LeaveCriticalSection( &m_numThreadsReachedSyncLock);
           }
	   return 1;
       }

       inline void release(int tid=-1)
       {
	 int taCount=m_taCount;
	 InterlockedIncrement(&m_allThreadsReachedSync);
	 int syncNumber = ++m_syncNumber;
	 m_numThreadsReachedSync = 0;
	 LeaveCriticalSection( &m_numThreadsReachedSyncLock);
	 TRC1(("BARTST#%d %d/%d R sn=%d\n",taCount,tid,m_threadCount,syncNumber));
       }

   private:
       uint32_t m_threadCount,
		m_taCount;
       CRITICAL_SECTION m_numThreadsReachedSyncLock;
       volatile LONG m_allThreadsReachedSync;
       uint32_t m_numThreadsReachedSync;
       uint32_t m_syncNumber;
   };
#endif


#if 0
class ThrBarrier
   {
   public:
       inline ThrBarrier(uint32_t numThreads) : m_threadCount(numThreads)
       {
           m_allThreadsReachedSync = CreateEvent( NULL, TRUE, FALSE, NULL );
           InitializeCriticalSection( &m_numThreadsReachedSyncLock);
           m_numThreadsReachedSync = 0;
           m_syncNumber = 0;
       }
   
       ~ThrBarrier() { 
	    DeleteCriticalSection( &m_numThreadsReachedSyncLock );
	    BOOL ok = CloseHandle( m_allThreadsReachedSync );
	    UT_ASSERT0(ok);

       
       }
#if 0
       inline void threadSync()
       {
           EnterCriticalSection( &m_numThreadsReachedSyncLock);
           uint32_t syncNr = m_syncNumber;
           m_numThreadsReachedSync++;
           if (m_numThreadsReachedSync == m_threadCount)
           {
               SetEvent(m_allThreadsReachedSync);
               m_syncNumber++;
               m_numThreadsReachedSync = 0;
               LeaveCriticalSection( &m_numThreadsReachedSyncLock);
           }
           else
           {
               if ( (m_syncNumber == syncNr) && (m_numThreadsReachedSync == 1) )
               {
                   ResetEvent(m_allThreadsReachedSync);
               }
               while (m_syncNumber == syncNr)
               {
                   LeaveCriticalSection( &m_numThreadsReachedSyncLock);
                   WaitForSingleObject(m_allThreadsReachedSync, 100);
                   EnterCriticalSection( &m_numThreadsReachedSyncLock);
               }
               LeaveCriticalSection( &m_numThreadsReachedSyncLock);
           }
       }
#else
       int  enter()
       {
           EnterCriticalSection( &m_numThreadsReachedSyncLock);
           uint32_t syncNr = m_syncNumber;
           m_numThreadsReachedSync++;
           if (m_numThreadsReachedSync == m_threadCount)
           {
	     return 0;
           }
           else
           {
               if ( (m_syncNumber == syncNr) && (m_numThreadsReachedSync == 1) )
               {
                   ResetEvent(m_allThreadsReachedSync);
               }
               while (m_syncNumber == syncNr)
               {
                   LeaveCriticalSection( &m_numThreadsReachedSyncLock);
                   WaitForSingleObject(m_allThreadsReachedSync, 100);
                   EnterCriticalSection( &m_numThreadsReachedSyncLock);
               }
               LeaveCriticalSection( &m_numThreadsReachedSyncLock);
           }
	   return 1;
       }

       inline void release()
       {
	 SetEvent(m_allThreadsReachedSync);
	 m_syncNumber++;
	 m_numThreadsReachedSync = 0;
	 LeaveCriticalSection( &m_numThreadsReachedSyncLock);
       }

#endif


   private:
       uint32_t m_threadCount;
       CRITICAL_SECTION m_numThreadsReachedSyncLock;
       HANDLE m_allThreadsReachedSync;
       uint32_t m_numThreadsReachedSync;
       uint32_t m_syncNumber;
   };

#endif

#if 0
//---------------------------------------------------------
  180 // pthread implementation without using pthread's barrier
  181 //---------------------------------------------------------
  182 
  183 class Barrier
  184 {
  185 public:
  186     inline Barrier(__uint32 numThreads) : m_threadCount(numThreads)
  187     {
  188         pthread_cond_init( &m_allThreadsReachedSync, NULL);
  189         pthread_mutex_init( &m_numThreadsReachedSyncLock, NULL);
  190         m_numThreadsReachedSync = 0;
  191         m_syncNumber = 0;
  192     }
  193 
  194     ~Barrier()
  195     {
  196         pthread_cond_destroy( &m_allThreadsReachedSync);
  197         pthread_mutex_destroy( &m_numThreadsReachedSyncLock);
  198     }
  199 
  200     inline void threadSync()
  201     {
  202         pthread_mutex_lock( &m_numThreadsReachedSyncLock);
  203         __uint32 syncNr = m_syncNumber;
  204         m_numThreadsReachedSync++;
  205         if (m_numThreadsReachedSync == m_threadCount)
  206         {
  207             m_syncNumber++;
  208             pthread_cond_signal( &m_allThreadsReachedSync);
  209             m_numThreadsReachedSync = 0;
  210         }
  211         else
  212         {
  213             while (syncNr == m_syncNumber)
  214                 pthread_cond_wait( &m_allThreadsReachedSync, &m_numThreadsReachedSyncLock);
  215         }
  216         pthread_mutex_unlock( &m_numThreadsReachedSyncLock);
  217     }
  218 private:
  219     __uint32 m_threadCount;
  220     pthread_mutex_t m_numThreadsReachedSyncLock;
  221     pthread_cond_t m_allThreadsReachedSync;
  222     __uint32 m_numThreadsReachedSync;
  223     __uint32 m_syncNumber;
  224 };
#endif

typedef volatile int32_t ThrSpinLockType
	 __attribute__((aligned(CPU_CACHELINE_SZ)));

inline void ThrSpinLock( ThrSpinLockType *lock, int64_t *nLockCollisions=NULL  )
{
  for(;;)
  {
    if( InterlockedCompareExchange((volatile LONG*)lock, 1, 0 ) == 0 )
    {
	break;
    }
    if(nLockCollisions) *nLockCollisions += 1;
    Sleep(0);
  }
}

inline void ThrSpinUnlock( ThrSpinLockType *lock )
{    
  InterlockedExchange( (LONG volatile *)lock, 0 );
}

/*-------------------------------------------------------------------------*/
/* 							                   */
/*-------------------------------------------------------------------------*/
// Full memory barrier
inline void ThrMemFence( void )
{
  __sync_synchronize();  
}

#if 1
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
class ThrSpinBarrier
{
  // 
  // Efficient threads barrier optimized for multi/many core architectures
  // Overhead: about 2000 CPU cycles on 60 CPU cores
  // https://software.intel.com/en-us/forums/intel-many-integrated-core/topic/392587
public:
  
  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
  ThrSpinBarrier( int nMaxThreads=THR_MAX_THREADS ) : 
    _nMaxThreads(nMaxThreads)
  {
    reset();
  }

  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
  void reset(void)
  {
    _counter=0;
    _passed=0;
    #ifdef THR_STATISTICS
    for(int i=0;i<ARRAY_LENGTH(_nWait);i++)
    {
      _nWait[i] = 0;
    }
    #endif
  }

#if 0 // test
  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
  int enter( int me=0, int nThreadsIn=0  )
  {
    int nThreads = _nMaxThreads;
    UT_ASSERT2(me>=0&&me<nThreads);
    if(nThreads==1) return 0;

    int passedOld = _passed;

    LONG counter = InterlockedIncrement( (LONG volatile *)&_counter);

    //TRCP(("enter %d/%d count=%d passed=%d\n",me,nThreads,counter,passedOld));  // UUUU

    //TRC(("THR:%d enter counter=%d\n",me,counter));

    UT_ASSERT0(counter>=1&&counter<=nThreads);
    if( counter == nThreads  )
    {
      return 0;
    }

    // wait
    for(int nWait=0;;nWait++)
    {
      int passed = // Atomic read with membarrier
	InterlockedCompareExchange((LONG volatile*)&_passed,0,0);
      if(passed != passedOld) break;

      static constexpr std::int32_t LOOPS_BEFORE_YIELD = 16;
      int count=1;
      if (count <= LOOPS_BEFORE_YIELD) 
      {
	  //machine_pause(count);
	  int delay=count;
	  while (delay-- > 0) { _mm_pause(); }

	  // Pause twice as long the next time.
	  count *= 2;
      } 
      else 
      {
	  // Pause is so long that we might as well yield CPU to scheduler.
	  //yield();

	#if 0
	SwitchToThread();
	#else
	std::this_thread::yield();
	#endif

      }
    }

    __sync_synchronize();

    return 1;
  }
#endif


#if 1 // original
  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
  int enter( int me=0, int nThreadsIn=0  )
  {
    int nThreads = _nMaxThreads;
    UT_ASSERT2(me>=0&&me<nThreads);
    if(nThreads==1) return 0;

    int passedOld = _passed;


    LONG counter = InterlockedIncrement( (LONG volatile *)&_counter);

    //TRCP(("enter %d/%d count=%d passed=%d\n",me,nThreads,counter,passedOld));  // UUUU

    //TRC(("THR:%d enter counter=%d\n",me,counter));

    UT_ASSERT0(counter>=1&&counter<=nThreads);
    if( counter == nThreads  )
    {
      return 0;
    }

    uint32_t h = passedOld ^ me;
    UT_FAST_HASH_32(h);
    int rnds = h+1;
    for(int nWait=0;;nWait++)
    {
      int passed = // Atomic read with membarrier
	InterlockedCompareExchange((LONG volatile*)&_passed,0,0);
      if(passed != passedOld) break;

#if 0

#if 0
      double u = UtFastRand(&rnds);
      //TRCP(("#%d u=%g\n",me,u));
      for(volatile int i=0;i<u*20;i++) _mm_pause();
#endif
#ifdef THR_STATISTICS
      InterlockedIncrement64((volatile LONGLONG *)&_nWait[me*8]);
#endif
    
      if(nWait%10==0)
      {
	BOOL ok = SwitchToThread();
	if(!ok)
	{
	  Sleep(0);
	}
      }
      
      /*if(nWait%100000==0)
      {
	Sleep(1);
      }*/

      /*if(nWait%10000==0)
      {
	Sleep(1);
      }*/

      YieldProcessor();
      _mm_pause();

#else
      Sleep(0); 
      continue;

      double u = UtFastRand(&rnds);
      ////TRCP(("#%d u=%g\n",me,u));
      //
      
      for(volatile int i=0;i<u*200;i++);

      #if 1
      if(nWait>200)
      {
	Sleep(0);
      }
      #endif

      /*if(nWait>20)
      {
        SwitchToThread();
	//Sleep(1);
	//TRCP(("WAIT %d -> %d\n",me,nWait));
      }*/
      _mm_pause();
#endif


    }

    __sync_synchronize();

    return 1;
  }
#endif

  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
  // Release method to be called by the master, i.e. the one to whom
  // calling enter() returned zero
  void release( int tid=-1 )
  {
    //TRCP(("release %d\n",tid));  // UUUU

    _counter = 0;
    __sync_synchronize();	
    InterlockedIncrement( (LONG volatile *)&_passed);
  }


  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
//private:
public:
  int _nMaxThreads;
  volatile LONG _counter,
	        _passed;
#ifdef THR_STATISTICS
  volatile LONGLONG _nWait[8*THR_MAX_THREADS];
#endif
};
#endif


#if 0
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
class ThrSpinBarrier
{
  // 
  // Efficient threads barrier optimized for multi/many core architectures
  // Overhead: about 2000 CPU cycles on 60 CPU cores
  // https://software.intel.com/en-us/forums/intel-many-integrated-core/topic/392587
public:
  
  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
  ThrSpinBarrier( int nMaxThreads=THR_MAX_THREADS ) : 
    _nMaxThreads(nMaxThreads),
    _counter(NULL)
  {
    //TRCP(("%s\n",UT_FUNCNAME));
    UT_ASSERT0(sizeof(*_counter)*8>=CPU_CACHELINE_SZ);
    UT_ASSERT0(sizeof(LONGLONG)==sizeof(int64_t));
    size_t szNeeded = sizeof(*_counter)*8*nMaxThreads;
    ALLOCMEM_ALIGNED_(_counter,int64_t,szNeeded,CPU_CACHELINE_SZ);
    memset((void*)_counter,0,szNeeded);
  }

  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
  ~ThrSpinBarrier( void )
  {
    //TRCP(("%s\n",UT_FUNCNAME));
    FREEMEM_ALIGNED(_counter);
  }

  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
  // Enter barrier and wait until all other threads are inside too.
  // It is guaranteed that the master thread (id 0) is the last one
  // to enter the barrier. In this case, enter() returns zero 
  // and at this point all other threads are waiting for the master to 
  // call the release() method in order for the whole team to proceed.
  int enter( int me, int nThreads )
  {
    UT_ASSERT2(nThreads<=_nMaxThreads);
    UT_ASSERT0(UT_ISPOW2(nThreads));  // is this limitation necessary ?
    UT_ASSERT2(me>=0&&me<nThreads);
    volatile int64_t *counter = _counter;
    int64_t myCount = counter[ me*8 ],
	    next = 1;
    
    if(nThreads>1)
    {
      while( next < nThreads - me )
      {
	if ((me & next) != 0) break;
	while( counter[(me+next)*8] <= myCount ) _mm_pause();
	next *= 2;
      }
    }

    if( me==0 ) 
    {
      // I'm the master thread: all others are now waiting for me
      // In the follwoing code, the master thread may want to 
      // modify shared data much like within a mutex protected critical section,
      // so we issue a memory barrier here (really necessary on x86 ?)
      __sync_synchronize();  
      return 0;
    }

    // Signal I'm done
    counter[me*8] += 1;

    // Waiting for master thread
    while (counter[0] <= myCount) _mm_pause();
    return 1;
  }

  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
  // Release method to be called by the master, i.e. the one to whom
  // calling enter() returned zero
  void release( void )
  {
    // Master signals done
    UT_ASSERT2(omp_get_thread_num()==0); 
    volatile int64_t *counter = _counter;

#if 0
    counter[0] += 1;
#else
    // Use interlocked increment to generate full memory barrier here
    // (not necessary for the non-master increments)
    InterlockedIncrement64( (LONGLONG volatile *)&counter[0]);
#endif
  }


  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
private:
  int64_t _nMaxThreads;
  volatile int64_t *_counter;
};

#endif

#if 0
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
class ThrSpinBarrier
{
  // 
  // Efficient threads barrier optimized for multi/many core architectures
  // Overhead: about 2000 CPU cycles on 60 CPU cores
  // https://software.intel.com/en-us/forums/intel-many-integrated-core/topic/392587
public:
  
  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
  ThrSpinBarrier( int nMaxThreads=THR_MAX_THREADS ) : 
    _nMaxThreads(nMaxThreads),
    _counter(NULL),
    _nWait(0),
    _nWaitMaster(0)
  {
    //TRCP(("%s\n",UT_FUNCNAME));
    UT_ASSERT0(sizeof(*_counter)*8>=CPU_CACHELINE_SZ);
    UT_ASSERT0(sizeof(LONGLONG)==sizeof(int64_t));
    size_t szNeeded = sizeof(*_counter)*8*nMaxThreads;
    ALLOCMEM_ALIGNED_(_counter,int64_t,szNeeded,CPU_CACHELINE_SZ);
    memset((void*)_counter,0,szNeeded);
  }

  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
  ~ThrSpinBarrier( void )
  {
    //TRCP(("%s\n",UT_FUNCNAME));
    //TRCP(("_nWait=%d, _nWaitMaster=%d\n",(int)_nWait,_nWaitMaster));
    FREEMEM_ALIGNED(_counter);
  }


  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
  // Enter barrier and wait until all other threads are inside too.
  // It is guaranteed that the master thread (id 0) is the last one
  // to enter the barrier. In this case, enter() returns zero 
  // and at this point all other threads are waiting for the master to 
  // call the release() method in order for the whole team to proceed.
  int enter( int me, int nThreads )
  {
    UT_ASSERT2(nThreads<=_nMaxThreads);
    UT_ASSERT0(UT_ISPOW2(nThreads));  // is this limitation necessary ?
    UT_ASSERT2(me>=0&&me<nThreads);
    volatile int64_t *counter = _counter;
    int64_t myCount = counter[ me*8 ],
	    next = 1;
    
    if(nThreads>1)
    {
      while( next < nThreads - me )
      {
	if ((me & next) != 0) break;
	while( counter[(me+next)*8] <= myCount ) 
	{
          InterlockedIncrement64( (LONGLONG volatile *)&_nWait);
	  _mm_pause();
	}
	next *= 2;
      }
    }

    if( me==0 ) 
    {
      // I'm the master thread: all others are now waiting for me
      // In the follwoing code, the master thread may want to 
      // modify shared data much like within a mutex protected critical section,
      // so we issue a memory barrier here (really necessary on x86 ?)
      __sync_synchronize();  
      return 0;
    }

    // Signal I'm done
    counter[me*8] += 1;

    // Waiting for master thread
    while (counter[0] <= myCount) 
    {
      InterlockedIncrement64( (LONGLONG volatile *)&_nWaitMaster);
      _mm_pause();
    }
    return 1;
  }

  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
  // Release method to be called by the master, i.e. the one to whom
  // calling enter() returned zero
  void release( void )
  {
    // Master signals done
    UT_ASSERT2(omp_get_thread_num()==0); 
    volatile int64_t *counter = _counter;

#if 0
    counter[0] += 1;
#else
    // Use interlocked increment to generate full memory barrier here
    // (not necessary for the non-master increments)
    InterlockedIncrement64( (LONGLONG volatile *)&counter[0]);
#endif
  }


  /*-----------------------------------------------------------------------*/
  /* 									   */
  /*-----------------------------------------------------------------------*/
private:
  int64_t _nMaxThreads;
  volatile int64_t *_counter;
public:
  volatile LONGLONG _nWait,
	            _nWaitMaster;
};

#endif

#endif // __cplusplus

#ifdef __cplusplus
extern "C" {
#endif

  void ThrInit(void);

void ThrGetMaxNumberOfThreads( int *pNumThreadsOMP,
    			      int *pNumThreadsTBB );

void ThrSetMaxNumberOfTBBThreads( int nMaxThreads );

#ifdef __cplusplus
}
#endif

#endif // THR_H

