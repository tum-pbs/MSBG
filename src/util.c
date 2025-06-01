/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#ifdef MIMP_ON_LINUX
#define _LARGEFILE64_SOURCE
#endif

#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>
#ifdef MIMP_ON_LINUX
#include <sys/mman.h> 
#endif
#include <math.h>
#include <float.h>
#include <string.h>
#include <pthread.h>
#include <errno.h>
#include <stdarg.h>
#include <stddef.h>
#include <ctype.h>
#include <limits.h>
#include "globdef.h"
#include "util.h"

/**** trace */

int LogLevel = 1;
int UtCheckLevel = 1;
char LogSaveLoc[256];
static UtMutex trcLock;
static int trcLockIsInitialized=0;

clock_t globClockStart,
	globClockStop;
float   globClockDiffMs;

#ifdef MIMP_ON_WINDOWS
static double hrTimerFreq;
static int hrTimerFreqValid=0;
#endif

static FILE *LogFile=NULL,
	    *LogFile2=NULL;  
static char LogFileName[UT_MAXPATHLEN+1]={0};

int UtGlobDbgCheckCount=0,
    UtGlobDbgErrCode=0;
UtGlobDbgErrInfoS UtGlobDbgErrInfo;

UtTimer UtGlobTimer;
int UtGlobTimerActive=0;

static char tempDir[UT_MAXPATHLEN+1] = "c:/tmp";

#if 0
int GWsetmsg( char *txt)
{
  printf("gwmsg: %s\n",txt);
  return 1;
}
#endif
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void UtSetTempDir( char *tmpDir )
{
  TRC(("%s: '%s'\n",UT_FUNCNAME,tmpDir?tmpDir:""));
  strcpy(tempDir,tmpDir?tmpDir:"");
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
char *UtGetTempDir( void ) 
{ 
  return tempDir; 
}

#ifdef MIMP_ON_LINUX

/*-------------------------------------------------------------------------------------*/
/* 									               */
/*-------------------------------------------------------------------------------------*/
float UtGetAverageCPULoad( double *deltaTimeMilliSec, UtTimer *tm )
{
  return -1;
}

#else

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static __int64 FileTimeToInt64 ( FILETIME ft )
{
    ULARGE_INTEGER foo;

    foo.LowPart = ft.dwLowDateTime;
    foo.HighPart = ft.dwHighDateTime;

    return ( foo.QuadPart );
}

/*-------------------------------------------------------------------------------------*/
/* 									               */
/*-------------------------------------------------------------------------------------*/
float UtGetAverageCPULoad( double *deltaTimeMilliSec, UtTimer *tm )
{
  FILETIME IdleTime, KernelTime, UserTime;
  
  static uint64_t PrevTotal_ = 0, *pPrevTotal;
  static uint64_t PrevIdle_ = 0, *pPrevIdle;
  static uint64_t PrevUser_ = 0, *pPrevUser;

  if(tm)
  {
    pPrevTotal = &tm->PrevTotal;
    pPrevIdle = &tm->PrevIdle;
    pPrevUser = &tm->PrevUser;
  }
  else
  {
    pPrevTotal = &PrevTotal_;
    pPrevIdle = &PrevIdle_;
    pPrevUser = &PrevUser_;

  }


  uint64_t ThisTotal;
  uint64_t ThisIdle, ThisKernel, ThisUser;
  uint64_t TotalSinceLast, IdleSinceLast, UserSinceLast;
    // GET THE KERNEL / USER / IDLE times.  
    // And oh, BTW, kernel time includes idle time
  
  BOOL ok = GetSystemTimes( & IdleTime, & KernelTime, & UserTime);

  if(!ok)
  {
    TRCERR(("GetSystemTimes() failed. errno=%d\n",
	  (int)GetLastError()));
    return MI_ERROR;
  }

  ThisIdle = FileTimeToInt64(IdleTime);
  ThisKernel = FileTimeToInt64 (KernelTime);
  ThisUser = FileTimeToInt64 (UserTime);

  ThisTotal = ThisKernel + ThisUser;
  TotalSinceLast = ThisTotal - *pPrevTotal;
  if(deltaTimeMilliSec) *deltaTimeMilliSec = TotalSinceLast/10000.;
  if(TotalSinceLast<=0) return -1;
  IdleSinceLast = ThisIdle - *pPrevIdle;
  UserSinceLast = ThisUser - *pPrevUser;
  double Headroom;
  Headroom =  (double)IdleSinceLast / (double)TotalSinceLast ;
  double Load;
  Load = 1.0 - Headroom;
  Headroom *= 100.0;  // to make it percent
  Load *= 100.0;  // percent

  *pPrevTotal = ThisTotal;
  *pPrevIdle = ThisIdle;
  *pPrevUser = ThisUser;

  return Load;
}

#endif


/*-------------------------------------------------------------------------------------*/
/* 									               */
/*-------------------------------------------------------------------------------------*/
uint64_t UtGetCurrentThreadId(void)
{
  #ifdef MIMP_ON_WINDOWS
  DWORD id = GetCurrentThreadId();
  #else
  unsigned long id = (unsigned long)pthread_self();
  #endif

  return id;
}

/*-------------------------------------------------------------------------------------*/
/* 									               */
/*-------------------------------------------------------------------------------------*/
static inline uintptr_t UtAlignUp(uintptr_t sz, size_t alignment) 
{
  UT_ASSERT0(alignment != 0);
  uintptr_t mask = alignment - 1;
  if ((alignment & mask) == 0) 
  {  // power of two?
    return ((sz + mask) & ~mask);
  }
  else 
  {
    return (((sz + mask)/alignment)*alignment);
  }
}
static void* UtAlignUpPtr(void* p, size_t alignment) 
{
  return (void*)UtAlignUp((uintptr_t)p, alignment);
}

static inline uintptr_t UtAlignDown(uintptr_t sz, size_t alignment) 
{
  UT_ASSERT0(alignment != 0);
  uintptr_t mask = alignment - 1;
  if ((alignment & mask) == 0) 
  { // power of two?
    return (sz & ~mask);
  }
  else 
  {
    return ((sz / alignment) * alignment);
  }
}
static void* UtAlignDownPtr(void* p, size_t alignment) 
{
  return (void*)UtAlignDown((uintptr_t)p, alignment);
}


/*-------------------------------------------------------------------------------------*/
/* 									               */
/*-------------------------------------------------------------------------------------*/
static void* UtPageAlign(
    		int isConservative, void* addr, size_t size, size_t* newsize) 
{
  UT_ASSERT0(addr != NULL && size > 0);
  if (newsize != NULL) *newsize = 0;
  if (size == 0 || addr == NULL) return NULL;

  size_t sysPageSize = UtGetSystemPageSize();

  void* start = (isConservative ? UtAlignUpPtr(addr,sysPageSize)
    : UtAlignDownPtr(addr,sysPageSize));
  void* end = (isConservative ? UtAlignDownPtr((uint8_t*)addr + size,sysPageSize)
    : UtAlignUpPtr((uint8_t*)addr + size,sysPageSize));
  ptrdiff_t diff = (uint8_t*)end - (uint8_t*)start;
  if (diff <= 0) return NULL;

  UT_ASSERT0((isConservative && (size_t)diff <= size) || 
      		(!isConservative && (size_t)diff >= size));
  if (newsize != NULL) *newsize = (size_t)diff;
  return start;
}



/*-------------------------------------------------------------------------------------*/
/* 									               */
/*-------------------------------------------------------------------------------------*/
int UtMemDisclaim( void *addr, size_t size,
    		   int doReleasePhysical,
		   int isPageAligned
		   )
{
  #ifdef MIMP_ON_WINDOWS

  int rcThis=0;
  size_t csize;
  void* start;

  if(!isPageAligned)
  {
    start = UtPageAlign(TRUE, addr, size, &csize);
  }
  else
  {
    #ifdef UT_ASSERT_LEVEL_2
    {
      size_t psz = UtGetSystemPageSize();
      uintptr_t ptr = (uintptr_t)addr; 
      UT_ASSERT0( size == ALIGN(size,psz) );
      UT_ASSERT0( ptr == ALIGN(ptr,psz) );
    }
    #endif

    start = addr;
    csize = size;
  }
  if (csize == 0) return MI_OK; 

  void* p = VirtualAlloc(start, csize, MEM_RESET, PAGE_READWRITE);

  if(p==NULL || p!=start)
  {
    TRCERRR(("VirtualAlloc(MEM_RESET) failed. errno=%d\n",
      (int)GetLastError()),MI_ERROR );
  }

  if(doReleasePhysical)
  {
    if(p == start && start != NULL) 
    {
      BOOL ok = VirtualUnlock(start,csize); 
      USE(ok);
      #if 0
      if(!ok)
      {
	TRCERR(("VirtualUnlock(%p %p) failed. errno=%d\n",
	      start,(void*)(uintptr_t)csize,
	      (int)GetLastError()));
      }
      #endif
    }
  }

rcCatch:
  return rcThis;

#else
  
  TRCERR(("%s not supported on non-Windows platform\n",UT_FUNCNAME));
  return MI_ERROR;

#endif
}


/*-------------------------------------------------------------------------------------*/
/* 									               */
/*-------------------------------------------------------------------------------------*/
int UtMemProtect( void *addr, size_t size,
    		  int doProtectWrite )
{
  size_t sysPageSz = UtGetSystemPageSize();
  UT_ASSERT0( ALIGN( ((uint64_t)addr), ((uint64_t)sysPageSz)) == (uint64_t)addr);
  UT_ASSERT0( ALIGN( ((uint64_t)size), ((uint64_t)sysPageSz)) == (uint64_t)size);

  int rcThis = 0;

#ifdef MIMP_ON_WINDOWS

  DWORD 
    protFlagsPrev =0,
    protFlags = doProtectWrite ? PAGE_READONLY : PAGE_READWRITE;

  BOOL ok = VirtualProtect(addr, size, protFlags, &protFlagsPrev );

  TRC(("UtMemProtect %d -> %d, rc=%d\n",protFlagsPrev,protFlags,ok));

  if(!ok)
  {
    rcThis = MI_ERROR;
    TRCERR(("VirtualProtect(%p,%d) failed. errno=%d\n",addr,size,
	  (int)GetLastError()));
  }

  return ok ? 0 : MI_ERROR;

#else

  int protFlags = doProtectWrite ? PROT_READ : (PROT_READ | PROT_WRITE);
  if (mprotect(addr, size, protFlags) != 0)
  {
    rcThis = MI_ERROR;
    TRCERR(("mprotect(%p,%zu) failed. errno=%d\n", addr, size, errno));
  }

#endif

  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void *UtAllocMemNuma( size_t sizeBytes )
{
  UT_ASSERT0(FALSE); // TODO: pin affinity of threads to processors
  return NULL;
#if 0
LPVOID WINAPI VirtualAllocExNuma(
  _In_     HANDLE hProcess,
  _In_opt_ LPVOID lpAddress,
  _In_     SIZE_T dwSize,
  _In_     DWORD  flAllocationType,
  _In_     DWORD  flProtect,
  _In_     DWORD  nndPreferred
);


#if 0
  ULONG HighestNodeNumber;
  if (!GetNumaHighestNodeNumber (&HighestNodeNumber))
  {
    TRCERR(("GetNumaHighestNodeNumber() failed. errno=%d\n",
	  (int)GetLastError()));
    return NULL;
  }

  if(HighestNodeNumber == 0)
  {
    TRCERR(("Not aNUMA system\n"));
    return NULL;
  }
#endif

  DWORD processorNumber = GetCurrentProcessorNumber();
  UCHAR NodeNumber=0xFF;

  UT_ASSERT0(THR_MAX_THREADS<=64); //  GetCurrentProcessorNumberEx() for more than 64 threads

  if (!GetNumaProcessorNode (processorNumber, &NodeNumber))
  {
    TRCERR(("GetNumaProcessorNode(%d) failed. errno=%d\n",
	  (int)processorNumber, (int)GetLastError()));
    return NULL;
  }

  UT_ASSERT0(NodeNumber!=0xFF);

  void *buffer = VirtualAllocExNuma(
      GetCurrentProcess(),
      NULL,
      sizeBytes,
      MEM_RESERVE | MEM_COMMIT,
      PAGE_READWRITE,
      NodeNumber
  );

  if(!buffer)
  {
    TRCERR(("VirtualAllocExNuma(%.0f,%d) failed. errno=%d\n",
      (double)sizeBytes,(int)NodeNumber,
      (int)GetLastError()));
    return NULL;
  }

  TRCP(("%s : proc=%d NUMA-Node=%d size=%.0f\n",UT_FUNCNAME,
	//omp_get_thread_num(),
	(int)processorNumber,(int)NodeNumber,(double)sizeBytes));

  memset(buffer,0,sizeBytes);

  return buffer;
#endif
}

#if 0
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void UtFreeMemNuma( void **ptr )
{
  if(*ptr)
  {
    BOOL ok = VirtualFree(*ptr, 0, MEM_RELEASE);
    if(!ok)
    {
      TRCERR(("VirtualFree(%p) failed. errno=%d\n",ptr,
	    (int)GetLastError()));
    }
    else *ptr=NULL;
  }
}
#endif

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
char *UtThousandSep( int64_t x, char *outbuf )
{
  LongInt n;
  int i=0;
  char tmpbuf[100];
  static char outbuf0[100];

  if(!outbuf) 
  {
    UT_ASSERT0(!omp_in_parallel());
    outbuf = outbuf0;
  }

  char *chp=tmpbuf;
  n=x;
  for(;;)
  {
    *chp++ = '0'+(n%10);
    i++;
    n/=10;	
    if(n<=0) break;
    if((i>2)&&(i%3)==0) *chp++='.';
  }

  int l=chp-tmpbuf;
  for(i=0;i<l;i++)
  {
    outbuf[i]=tmpbuf[l-1-i];
  }
  outbuf[l]=0;
  return outbuf;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtFileLastModified( char *path, 		// IN
    			time_t *pTimeModified   // OUT 
			)
{
  int rcThis=0;
  struct stat file_stat;
  int err = stat(path, &file_stat);
  if(err!=0)
  {
    TRCERRR(("stat('%s') failed rc=%d. errno=%d\n",
	  path, err, errno),MI_ERROR );    
  }
  *pTimeModified = file_stat.st_mtime;
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtTestFileRead(char *fname)
{
  int rcThis=0,rc;
  size_t maxsz=3*((size_t)ONE_GB);
  UtTimer tm;

  USE(rc);

  TIMER_START(&tm);
  void *buffer=NULL;

  ALLOCMEM(buffer,maxsz);

#if 0
  //
  // Win ReadFile
  //
  DWORD actsz; 
  HANDLE hFile = CreateFile(fname,GENERIC_READ|GENERIC_WRITE, 
  		FILE_SHARE_READ,NULL,OPEN_ALWAYS,FILE_ATTRIBUTE_NORMAL,NULL); 

  if (hFile == INVALID_HANDLE_VALUE) 
  {
    TRCERRR(("can't open file '%s'\n",fname),MI_EIO);
  }

  ReadFile(hFile,buffer,maxsz,&actsz,NULL); 

  CloseHandle(hFile); 
  hFile=INVALID_HANDLE_VALUE;

  TIMER_STOP(&tm);
  TRCP(("CPU (ReadFile) = %.4f sec\n", 
	(double)(TIMER_DIFF_MS(&tm)/1000.0)));

  TRCP(("%lu bytes read\n",(unsigned long)actsz));
#endif

#if 0
  //
  // open/ead
  //
  UT_ASSERT_FATAL(FALSE);
  size_t fsize=2334523392;
  TRCP(("size=%.0f\n",(double)fsize));

  TIMER_START(&tm);
  rc = UtReadWriteBinData(0,buffer,fsize,fname,0);
  UT_ASSERT(rc==0);
  TIMER_STOP(&tm);
  TRCP(("CPU (UtReadWriteBinData) = %.4f sec\n", 
	(double)(TIMER_DIFF_MS(&tm)/1000.0)));
#endif

rcCatch:
  FREEMEM(buffer);
  return rcThis;
}


#if 0
/*-------------------------------------------------------------------------*/
/* spinlocks							           */
/*-------------------------------------------------------------------------*/

typedef LONG UtSpinlock;
#define SMP_SPIN_COUNT	10
#define SMP_SPIN_LOCKED 1
#define SMP_SPIN_UNLOCKED 0

#define SmpSpinlockInit( sl_, val_ ) \
{ \
  *((LONG volatile *)(sl_)) = val_; \
} 

#define SmpSpinLock( sl_ ) \
{ \
  if (InterlockedCompareExchange((LONG volatile *)(sl_), 1, 0)) \
  { \
    int i_; \
    for (i_=0;;i_++) { \
      if ((*((LONG volatile *)(sl_))==0) && \
	  !InterlockedCompareExchange((LONG volatile *)(sl_), 1, 0)) break; \
      if(i_ > SMP_SPIN_COUNT) \
      { \
	i_ = 0; \
	Sleep(0); \
      } \
      YieldProcessor();; \
    } \
  } \
}

#define SmpSpinUnlock( sl_ ) \
  InterlockedExchange ((LONG volatile *)(sl_), 0);

void *UtCreateSpinlock(  )
{
  UtSpinlock *s=NULL;
  ALLOCOBJ0(s);
  SmpSpinlockInit(s,SMP_SPIN_UNLOCKED);
  return s;
}

void UtAcquireSpinlock( void *lock )
{
  static const unsigned int YIELD_ITERATION = 100; 
  static const unsigned int MAX_SLEEP_ITERATION = 40; 

  int m_iterations = 0;
  while(TRUE)
  {
    if(InterlockedCompareExchange((LONG volatile *)(lock), 1, 0) == 0)
    {
      break;			
    }

    // spin wait to acquire 
    while(*(LONG volatile *)(lock) != 0)
    {
      if((m_iterations >= YIELD_ITERATION))
      {
	if(m_iterations + YIELD_ITERATION >= MAX_SLEEP_ITERATION)
	  Sleep(0);

	if(m_iterations >= YIELD_ITERATION && m_iterations < MAX_SLEEP_ITERATION)
	{
	  m_iterations = 0;
	  SwitchToThread();
	}
      }
      // Yield processor on multi-processor but if on single processor then give other thread the CPU
      m_iterations++;
      YieldProcessor(/*no op*/); 
    }				
  }
}

void UtReleaseSpinlock( void *lock )
{
  InterlockedExchange((LONG volatile *)(lock), 0);
}
#endif

#ifdef MIMP_ON_WINDOWS

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
static int WinPrivilege(TCHAR* pszPrivilege, BOOL bEnable)
{
  int rcThis=0;
  HANDLE      hToken;
  TOKEN_PRIVILEGES tp;
  BOOL       status;
  DWORD      error;

  OpenProcessToken(GetCurrentProcess(), TOKEN_ADJUST_PRIVILEGES | TOKEN_QUERY, &hToken);
  LookupPrivilegeValue(NULL, pszPrivilege, &tp.Privileges[0].Luid);
  tp.PrivilegeCount = 1;
  tp.Privileges[0].Attributes = SE_PRIVILEGE_ENABLED;
  status = AdjustTokenPrivileges(hToken, FALSE, &tp, 0, (PTOKEN_PRIVILEGES)NULL, 0);
  error = GetLastError();
  if(error)
  {
    TRCERRR(("AdjustTokenPrivileges() failed. errno=%d\n",
      (int)GetLastError()),MI_ERROR );
  }
  CloseHandle(hToken);
rcCatch:
  return rcThis;
}

#endif

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
uint32_t UtAtomicCompareExchange( uint32_t volatile *dest, 
    			          uint32_t exchange, uint32_t comparand,
				  unsigned opt )
{
  return InterlockedCompareExchange((LONG volatile *)(dest), 
      				     exchange, comparand);
}
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
uint32_t UtAtomicExchange( uint32_t volatile *dest, 
    			          uint32_t exchange,
				  unsigned opt )
{
  return InterlockedExchange((LONG volatile *)(dest), exchange);
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void *UtAllocLargePageMem( size_t size )
{
#ifdef MIMP_ON_WINDOWS

  int rcThis=0;
  void *ptr=NULL;
  size_t pagesz = GetLargePageMinimum();

  UT_ASSERT_FATAL(pagesz > 4*ONE_KB &&
      		  !(((pagesz) & ((pagesz) - 1))));  // pow-of-2
  
  TRCP(("GetLargePageMinimum = %d\n",(int)pagesz));

  size_t actsz = (size + pagesz - 1) & (~(pagesz-1));  
  
  TRCP(("UtAllocLargePageMem: reqsize=%.0f, pagesz=%.0f, actsize=%.0f\n",
	(double)size,(double)pagesz,(double)actsz)); 

  WinPrivilege(TEXT("SeLockMemoryPrivilege"), TRUE);

  ptr = VirtualAlloc( NULL, actsz,
//      			MEM_COMMIT|MEM_RESERVE |

      			MEM_COMMIT |
      			  MEM_LARGE_PAGES, 
			PAGE_READWRITE );
  if(!ptr)
  {
    TRCERRR(("VirtualAlloc() failed. errno=%d\n",
      (int)GetLastError()),MI_ERROR );
  }

rcCatch:
  return ptr;

#else

  UT_ASSERT_NR(FALSE);
  return NULL;
 
#endif

}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void UtFreeLargePageMem( void *address )
{
#ifdef MIMP_ON_WINDOWS

  BOOL ok = VirtualFree(address, 0, MEM_RELEASE);
  if(!ok)
  {
    TRCERR(("VirtualFree(%p) failed. errno=%d\n",address,
	  (int)GetLastError()));
  }

#else

  UT_ASSERT_NR(FALSE);
 
#endif
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
/* Parse S into tokens separated by characters in DELIM.
   If S is NULL, the saved pointer in SAVE_PTR is used as
   the next starting point. For example:
   char s[] = "-abc-=-def";
   char *sp;
   x = strtok_r(s, "-", &sp); // x = "abc", sp = "=-def"
   x = strtok_r(NULL, "-=", &sp); // x = "def", sp = NULL
   x = strtok_r(NULL, "=", &sp); // x = NULL
// s = "abc\0-def\0"
*/
char *UtStrtok(char *s,
    	       const char *delim,
    	       char **save_ptr)
{
  char *token;

  if (s == NULL)
    s = *save_ptr;

  /* Scan leading delimiters. */
  s += strspn (s, delim);
  if (*s == '\0')
    return NULL;

  /* Find the end of the token. */
  token = s;
  s = strpbrk (token, delim);
  if (s == NULL)
    /* This token finishes the string. */
    *save_ptr = strchr (token, '\0');
  else
  {
    /* Terminate the token and make *SAVE_PTR point past it. */
    *s = '\0';
    *save_ptr = s + 1;
  }
  return token;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtSleep( unsigned millisec )
{
  Sleep(millisec);

  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtParseFloat( char *tok, double *pdbl )
{
  int rcThis=0;
  char *endptr;

  *pdbl = strtod(tok,&endptr);
  // TODO: also check errno (over,underflow etc. see 'man strtod')
  if(endptr==tok)
  {
    TRCERRR(("invalid number format: '%s'\n",tok),MI_ESYNTAX);
  }

rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtParseInt( char *tok, unsigned opt, LongInt *plong )
{
  int rcThis=0;
  char *endptr;
  long long ll;

  ll = strtoll(tok,&endptr,10);
  if((endptr==tok)||(endptr&&*endptr!='\0'))
  {
    if(opt&UT_NOERRTRC) raiseRc(MI_ESYNTAX);
    TRCERRR(("invalid number format: '%s'\n",tok),MI_ESYNTAX);
  }
  if(ll==LLONG_MAX||ll==LLONG_MIN)
  {
    if(opt&UT_NOERRTRC) raiseRc(MI_ESYNTAX);
    TRCERRR(("overflow in integr conversion: %s\n",tok),
	MI_ESYNTAX);
  }
  *plong = (LongInt)ll;
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtParseFloatArray( double *f, int *pLen, char *tok_in, char *delim )
{
  int i=0,rcThis=0,len,maxlen=0,rc;
  double fval;
  char *chp=NULL, *tok;
  chp = MM_strdup(tok_in);
  if(!chp) TRCERRR(("MM_strdup() failed\n"),MI_ERROR);
  if(pLen) maxlen=*pLen;
  for(tok=strtok(chp,delim),i=0;
      tok;
      tok=strtok(NULL,delim),i++)
  {
    if(pLen && i>=maxlen) break;
    rc = UtParseFloat(tok,&fval);
    if(rc) TRCERRR(("UtParseFloat() failed %d (string='%s')\n",tok),
	MI_ESYNTAX);
    f[i] = fval;
//    UT_ASSERT(isfinite(f[i]));
    UT_ASSERT(!isnan(f[i]));
  }
rcCatch:
  len = i;
  if(pLen) *pLen = len;
  if(chp) MM_free(chp);
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtGetFileSize( char *fpath, size_t *p_size )
{
  int rcThis=0,rc,fd=-1;
  *p_size =0;
  
  fd = open( fpath, O_RDONLY| O_BINARY ,
      		    S_IRUSR | S_IWUSR );
  if (fd < 0)
  {
    TRCERRR(("can't open file '%s' for reading. errno=%d\n",
	  fpath,(int)errno),MI_EIO);
  }

  off64_t fsize = lseek64(fd, 0, SEEK_END);
  if(fsize<0)
  {
    TRCERRR(("lseek64 failed for file '%s' errno=%d\n",
	  fpath,(int)errno),MI_EIO);
  }

  *p_size = fsize;

rcCatch:

  if(fd!=-1) 
  {
    rc = close(fd);
    fd=-1;
    if(rc!=0)
    {
      TRCERR(("close() failed for file '%s' errno=%d\n",fpath,(int)errno));
      rcThis=MI_EIO;
    }
  }

  return rcThis;
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtReadWriteBinData2( int do_write, 
    			     const void *data_ptr, size_t data_sz, 
    			     char *fname,
			     int doUnbuffered )
{
  #ifndef MIMP_ON_WINDOWS

  return UtReadWriteBinData( do_write, data_ptr, data_sz, fname, -1, 0 );

  #else

  int rcThis=0;
  LongInt szio,file_offs_total=0;
  LongInt trans_sz,total_sz,chunk_sz;
  BOOL ok;

  TRC(("%s: write=%d data=%p sz=%.0f name='%s' %d\n",UT_FUNCNAME,
	do_write,data_ptr,(double)data_sz,fname,doUnbuffered));

  if(doUnbuffered)
  {
    // https://docs.microsoft.com/en-us/windows/win32/fileio/file-buffering
    uintptr_t ptr = (uintptr_t)data_ptr;
    size_t psz=4096;
    UT_ASSERT0( data_sz == ALIGN(data_sz,psz) );
    UT_ASSERT0( ptr == ALIGN(ptr,psz) );
  }

  HANDLE f = CreateFileA(
      fname,
      do_write ? GENERIC_WRITE : GENERIC_READ,
      0,
      NULL,
      do_write ? CREATE_ALWAYS : OPEN_ALWAYS,

      FILE_ATTRIBUTE_NORMAL 
        | (doUnbuffered ? (FILE_FLAG_NO_BUFFERING|FILE_FLAG_WRITE_THROUGH ) : 0)
	,
      
      NULL
      );

  if( f == INVALID_HANDLE_VALUE )
  {
    TRCERRR(("%s: Error %s file '%s' ernno=%d\n",UT_FUNCNAME,
	  do_write?"writing":"reading",
	  fname,(int)GetLastError()),1);
  }
    
  //
  // Transfer data (in chunks to avoid possible >2GB problems)
  //
  trans_sz = 0;
  total_sz = data_sz;
  chunk_sz = 2*(LongInt)(ONE_GB);

  while(trans_sz<total_sz)
  {
    szio = MIN(chunk_sz,total_sz-trans_sz);
    UT_ASSERT(((LongInt)szio)<=2*ONE_GB);

    TRC(("transfer chunk=%.0f total=%.0f trans=%.1f%%\n", 
	  (double)szio,(double)total_sz,
	  100.0*(double)trans_sz/(double)total_sz)); 

    DWORD szio_act;

    if(do_write)
    {
      ok = WriteFile( f, ((unsigned char*)data_ptr)+trans_sz, szio, &szio_act, NULL ); 
    }
    else
    {
      ok = ReadFile( f, ((unsigned char*)data_ptr)+trans_sz, szio, &szio_act, NULL ); 
    }

    if(!ok || szio != (LongInt)szio_act )
    {
      TRCERRR(("%s: Error %s file '%s' %ld %ld ernno=%d\n",UT_FUNCNAME,
	    do_write?"writing":"reading",
	    fname,(long)szio,(long)szio_act,(int)GetLastError()),1);
    }

    trans_sz += szio_act;
    file_offs_total += szio_act;
  }        

  TRC(("transferred total of %.0f bytes\n", (double)file_offs_total)); 

  //
  // Close/cleanup
  //

rcCatch:
  {
    ok = CloseHandle(f);
    if(!ok)
    {
      TRCERR(("%s: close() failed for file '%s' errno=%d\n",UT_FUNCNAME,
	    fname,(int)GetLastError()));
      rcThis=1;
    }
  }

  return rcThis;

  #endif
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtReadWriteBinData( int do_write, 
    			     const void *data_ptr, size_t data_sz, 
    			     char *fname, int fd,
			     unsigned opt )
{
  int rc,rcThis=0,fd_supplied=TRUE;
  LongInt szio,szio_act,file_offs_total=0;
  LongInt trans_sz,total_sz,chunk_sz;

  TRC(("UtReadWriteBinData: write=%d data=%p sz=%.0f name='%s'\n",
	do_write,data_ptr,(double)data_sz,fname));

  if(!data_ptr) raiseRc(0);	// ignore NULL-pointers

  //
  // Open file
  //
  if(fd==-1)
  {
    fd_supplied = FALSE;

    fd = open( fname, 
	       O_BINARY | ( do_write ? O_WRONLY| O_CREAT : 
				       O_RDONLY), 
	       S_IRUSR | S_IWUSR );

    if (fd < 0)
    {
      TRCERRR(("can't open file '%s' for %s. errno=%d\n",
	    fname,do_write?"writing":"reading",
	    (int)errno),1);
    }
  }

  //
  // Transfer data (in chunks to avoid possible >2GB problems)
  //
  trans_sz = 0;
  total_sz = data_sz;
  chunk_sz = 128*1024*1024;

  while(trans_sz<total_sz)
  {
    szio = MIN(chunk_sz,total_sz-trans_sz);
    UT_ASSERT(((LongInt)szio)<2000000000);
    
    TRC(("transfer chunk=%d total=%.0f trans=%.1f%%\n", 
	  (int)szio,(double)total_sz,
	  100.0*(double)trans_sz/(double)total_sz)); 

    szio_act = do_write ? 
      		write( fd, ((unsigned char*)data_ptr)+trans_sz, szio ) :
		read( fd, ((unsigned char*)data_ptr)+trans_sz, szio );
    if( szio_act != szio ) 
    {
      TRCERRR(("Error %s file '%s' %ld %ld ernno=%d\n",
	    do_write?"writing":"reading",
	    fname,(long)szio,(long)szio_act,(int)errno),1);
    }
    trans_sz += szio_act;
    file_offs_total += szio_act;
  }        

  TRC(("transferred total of %.0f bytes\n", (double)file_offs_total)); 
  
  //
  // Close/cleanup
  //

rcCatch:
  if(!fd_supplied && fd!=-1) 
  {
    rc = close(fd);
    fd=-1;
    if(rc!=0)
    {
      TRCERR(("close() failed for file '%s' errno=%d\n",fname,(int)errno));
      rcThis=1;
    }
  }
  
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
size_t UtGetSystemPageSize( void )
{
  #ifdef MIMP_ON_WINDOWS

  static size_t systemPageSize=0;
  if(!systemPageSize)
  {
    SYSTEM_INFO sSysInfo;
    GetSystemInfo(&sSysInfo);
    systemPageSize = sSysInfo.dwPageSize;
  }
  return systemPageSize;

  #else

  return getpagesize();

  #endif
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void UtGetSystemMemSize( size_t *p_phys_total, 
    			 size_t *p_phys_avail,
			 size_t *p_commit_avail,  // max. per process 
			 size_t *p_virtual_avail )

{
  #ifdef MIMP_ON_WINDOWS

  MEMORYSTATUSEX status={0};
  BOOL ok;
  status.dwLength = sizeof(status);
  ok = GlobalMemoryStatusEx( &status );
  if( ! ok )
  {
    TRCERR(("GlobalMemoryStatusEx(%p) failed. errno=%d\n",
	  (int)GetLastError()));    
    return;
  }    

  TRC(("%s: phys=%.1f/%.1f swap=%.1f virt=%.1f  (in MB)\n",
	UT_FUNCNAME,
	(double)((size_t)status.ullAvailPhys)/ONE_MB,
	(double)((size_t)status.ullTotalPhys)/ONE_MB,
	(double)((size_t)status.ullAvailPageFile)/ONE_MB,
	(double)((size_t)status.ullAvailVirtual)/ONE_MB	
	));

  if(p_phys_total) *p_phys_total = (size_t)status.ullTotalPhys;
  if(p_phys_avail) *p_phys_avail = (size_t)status.ullAvailPhys;
  if(p_commit_avail ) *p_commit_avail = (size_t)status.ullAvailPageFile;
  if(p_virtual_avail ) *p_virtual_avail = (size_t)status.ullAvailVirtual;

  #else

  size_t pages = sysconf(_SC_PHYS_PAGES);
  size_t page_size = sysconf(_SC_PAGE_SIZE);
  size_t phys_total = pages * page_size;

  if(p_phys_total) *p_phys_total = phys_total; 
  if(p_phys_avail) *p_phys_avail = 0;
  if(p_commit_avail ) *p_commit_avail = 0;
  if(p_virtual_avail ) *p_virtual_avail = 0;

  //TRCERR(("%s not availbale on non-windows platform\n"));

  #endif
}

/*--------------------------------------------------------------------*/
/* 								      */
/*--------------------------------------------------------------------*/
static int LogIErr( int isWarn, const char *format, va_list arg_ptr )
{
  int     rc,isInfo=0;
  char    buf[512],buf2[512];

  if(isWarn==2)
  {
    isInfo = 1;
    isWarn = 0;
  }


#if 0  // Should be already locked by LogLock() via TRC... macros 
  static UtMutex lock;
  static int lock_is_initialized=FALSE;

  if(!lock_is_initialized)
  {
    UtMutexInit(&lock);
    lock_is_initialized = TRUE;
  }
#endif

  sprintf(buf2,"%s",isInfo ? "INFO:" : isWarn ? "WARNING:":"ERROR:");
  rc = vsnprintf( buf, sizeof(buf)-1, format, arg_ptr );
  MiStripTrailingCRLF( buf, 0 );


#if 0
  UtMutexLock( &lock );
#endif

  #ifdef MIMP_ON_WINDOWS

     // Windows color attributes 
   // http://git.vger.kernel.narkive.com/FqeWCXtw/need-your-help-with-mingw-issue-17-color-options-don-t-work-produce-garbage
  WORD outcol = isInfo ? FOREGROUND_INTENSITY | FOREGROUND_GREEN :
		isWarn ? FOREGROUND_INTENSITY | FOREGROUND_GREEN | FOREGROUND_RED :
		         FOREGROUND_INTENSITY | FOREGROUND_RED;

  CONSOLE_SCREEN_BUFFER_INFO ConsoleInfo;
  GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &ConsoleInfo);
  WORD attrOld = ConsoleInfo.wAttributes;

  SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), outcol );

  #else

  fprintf( stdout, 
      isInfo ? "\x1b[32m" : /*ANSI_COLOR_GREEN */
      isWarn ? "\x1b[33m" : /*ANSI_COLOR_YELLOW */ 
      	       "\x1b[31m" );  /*ANSI_COLOR_RED */
  #endif

  if(isInfo)
    fprintf(stdout,"%s %s\n", buf2, buf);   
  else
    fprintf(stdout,"*** %s %s %s\n", buf2, buf,LogSaveLoc);   

  #ifdef MIMP_ON_WINDOWS

  SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), attrOld );

  #else

  fprintf( stdout, "\x1b[0m" /*ANSI_COLOR_RESET */ );

  #endif

#if 0
  UtMutexUnlock( &lock );
#endif

  if(LogFile)
  {
    if(isInfo)
      fprintf(LogFile,"%s %s\n", buf2, buf);   
    else
      fprintf(LogFile,"*** %s %s %s\n",buf2, buf,LogSaveLoc);   
    fflush(LogFile);
  }
  if(LogFile2)
  {
    if(isInfo)
      fprintf(LogFile2,"%s %s\n", buf2, buf);   
    else
      fprintf(LogFile2,"*** %s %s %s\n",buf2, buf,LogSaveLoc);   
    fflush(LogFile2);
  }

  return rc;
}

/*--------------------------------------------------------------------*/
/* 								      */
/*--------------------------------------------------------------------*/
void LogLock(void)
{
  if(!trcLockIsInitialized)
  {
    UtMutexInit(&trcLock);
    trcLockIsInitialized = TRUE;
  }
  UtMutexLock( &trcLock );
}

void LogUnlock(void)
{
  if(!trcLockIsInitialized)
  {
    fprintf(stderr,"ERROR: LogUnlock(): internal error\n");
  }
  UtMutexUnlock( &trcLock );
}

/*--------------------------------------------------------------------*/
/* 								      */
/*--------------------------------------------------------------------*/
int LogThrowError( const char*file, int line )
{
  LogLock(); 
  sprintf(LogSaveLoc,"  [%-36.36s %d]",file,line); 
  LogErr("Assertion failed (fatal).");  
  LogUnlock();   
  exit(1);
  return 0;
}

/*--------------------------------------------------------------------*/
/* 								      */
/*--------------------------------------------------------------------*/
int LogErr( const char *format, ...)
{
  int rc;
  va_list arg_ptr;
  va_start (arg_ptr, format);
  rc = LogIErr(0,format,arg_ptr);
  va_end (arg_ptr);
  return rc;
}

/*--------------------------------------------------------------------*/
/* 								      */
/*--------------------------------------------------------------------*/
int LogWarn( const char *format, ...)
{
  int rc;
  va_list arg_ptr;
  va_start (arg_ptr, format);
  rc = LogIErr(1,format,arg_ptr);
  va_end (arg_ptr);
  return rc;
}

/*--------------------------------------------------------------------*/
/* 								      */
/*--------------------------------------------------------------------*/
int LogInfo( const char *format, ...)
{
  int rc;
  va_list arg_ptr;
  va_start (arg_ptr, format);
  rc = LogIErr(2,format,arg_ptr);
  va_end (arg_ptr);
  return rc;
}

#if 0
/*--------------------------------------------------------------------*/
/* 								      */
/*--------------------------------------------------------------------*/
int LogPrefix( const char *prefix, const char *format, ...)
{
  int rc;
  va_list arg_ptr;
  va_start (arg_ptr, format);
  rc = LogIPrefix(prefix,format,arg_ptr);
  va_end (arg_ptr);
  return rc;
}
#endif


#if 1
/*--------------------------------------------------------------------*/
/* 								      */
/*--------------------------------------------------------------------*/
int LogPrnLog( const char *format, ...)
{
  va_list arg_ptr;
  int     rc;

  LogLock();
 
  va_start (arg_ptr, format);
  rc = vfprintf( stdout, format, arg_ptr );
  va_end (arg_ptr);


  if(LogFile)
  {
    va_start (arg_ptr, format);
    rc = vfprintf( LogFile, format, arg_ptr );
    va_end (arg_ptr);
    fflush(LogFile);
  }
  if(LogFile2)
  {
    va_start (arg_ptr, format);
    rc = vfprintf( LogFile2, format, arg_ptr );
    va_end (arg_ptr);
    fflush(LogFile2);
  }

  LogUnlock();
  
  return rc;
}
#else
/*--------------------------------------------------------------------*/
/* 								      */
/*--------------------------------------------------------------------*/
int LogPrnLog( const char *format, ...)
{
  va_list arg_ptr;
  int     rc;

  LogLock();
 
  va_start (arg_ptr, format);
  rc = vfprintf( stdout, format, arg_ptr );
  if(LogFile)
  {
    rc = vfprintf( LogFile, format, arg_ptr );
    fflush(LogFile);
  }
  if(LogFile2)
  {
    rc = vfprintf( LogFile2, format, arg_ptr );
    fflush(LogFile2);
  }
  va_end (arg_ptr);

  LogUnlock();
  
  return rc;
}
#endif



#if 0
int LogPrnLogTST( const char *format )
{
  printf("LogPrnLogTST: %s\n",format);
  return 0;
}
#endif
/*--------------------------------------------------------------------*/
/* 								      */
/*--------------------------------------------------------------------*/
int LogLog( const char *format, ...)
{
  va_list arg_ptr;
  int     rc=0;
 
  LogLock();
  
  if(LogFile)
  {
    va_start (arg_ptr, format);
    rc = vfprintf( LogFile, format, arg_ptr );
    va_end (arg_ptr);
    fflush(LogFile);
  }
  if(LogFile2)
  {
    va_start (arg_ptr, format);
    rc = vfprintf( LogFile2, format, arg_ptr );
    va_end (arg_ptr);
    fflush(LogFile2);
  }

  LogUnlock();

  return rc;
}
                             
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtOpenLogFile(char *fname, int do_copy_old, int do_append)
{
  int rcThis=0,rc;
  char fname_old[1024];

  strcpy(LogFileName,fname);

  if(do_copy_old&&(!do_append))
  {
    UT_ASSERT(strlen(fname)+20<sizeof(fname_old));
    sprintf(fname_old,"%s.old",fname);  
    rc = UtCopyFile(fname,fname_old,0);
    USE(rc);
  }

  LogFile = fopen(fname,do_append?"at":"wt");
  if(!LogFile)
  {
    TRCERRR(("can't open trace file '%s'\n",fname),MI_ERROR);
  }
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtOpenLogFile2(char *fname)
{
  int rcThis=0,rc;
  TRC(("%s: %s\n",UT_FUNCNAME,fname));

  UT_ASSERT0(LogFile);
  UT_ASSERT0(*LogFileName);

  if(!LogFile2)
  {
    // Copy existing trace file
    rc = UtCopyFile(LogFileName,fname,0);
    UT_ASSERT(rc==0);

    // Open for append
    LogFile2 = fopen(fname,"at");
    if(!LogFile2)
    {
      TRCERRR(("can't open trace file '%s'\n",fname),MI_ERROR);
    }
  }
rcCatch:
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtCloseLogFile()
{
  if(LogFile) 
  {
    fclose(LogFile);
    LogFile = NULL;
    *LogFileName = 0;
  }
  if(LogFile2) 
  {
    fclose(LogFile2);
    LogFile2 = NULL;
  }
  return MI_OK;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtCloseLogFile2()
{
  if(LogFile2) 
  {
    fclose(LogFile2);
    LogFile2=NULL;
  }
  return MI_OK;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtProgrIndDisplay( UtProgrInd *pgi, float percCompleted, char *chbuf0, 
    		     int (*callback)( char *msg, int l) )
{
#define PGI_MIN_SPEED 0.00001
  char chbuf[1000];
  double dMs,remMs,dPerc,actSpeed=0;

  strcpy(chbuf,chbuf0);
  int len = strlen(chbuf); 
  if(pgi->tmStepSec>0)
  {
    double eps = 1./1000000000.0;
    int dSec;
    time_t tmAct = time(NULL);
    dSec = (double)difftime(tmAct,pgi->tmStampSec);

    // Calculate new speed every second
    if(dSec>=1)
    {
      TIMER_STOP_INT(&pgi->tm);
      dMs = TIMER_DIFF_MS(&pgi->tm);
      actSpeed = (double)percCompleted / ((double)dMs);
      pgi->actSpeed = actSpeed;
    }
    else actSpeed = pgi->actSpeed;
    remMs = actSpeed > eps ? 
      (1.0-percCompleted) / actSpeed : 0;
    sprintf(chbuf+len," %.0f%% remaining: %.0f sec",
	(double)(100.*percCompleted), 
	actSpeed < eps ? - 1: (remMs/1000.)); 

    // Update percentage display
    callback( chbuf, 0 );

    // Trigger event each d seconds
    if(dSec>=pgi->tmStepSec)
    {
      pgi->tmStampSec = tmAct;
      pgi->trig = TRUE;
    }
    return 0;
  }
  else if(pgi->doTime)
  {
    TIMER_STOP(&pgi->tm);
    dMs = TIMER_DIFF_MS(&pgi->tm);
    TIMER_START(&pgi->tm);
    dPerc = percCompleted - pgi->percCompletedOld;
    actSpeed = dMs > 1 ? (float)dPerc / (float)dMs : -1;
    if(actSpeed < 0) 
      actSpeed = (pgi)->actSpeed;
    if(actSpeed > 0)
    {
      remMs = (1.0-percCompleted) / actSpeed;
      sprintf(chbuf+len,"%.0f%% remaining: %d sec",
	  (float)(100.*percCompleted), (int)(remMs/1000.)); 
      pgi->actSpeed = actSpeed;
    }
  }
  else
  {
    sprintf(chbuf+len,"%.0f%%",(float)(100.*percCompleted)); 
  }
  callback( chbuf, 0 );
  pgi->trig = TRUE;
  return 0;
}


#if 0
/*-------------------------------------------------------------------------------------*/
/* 										       */
/*-------------------------------------------------------------------------------------*/
int UtGetKeyPressed( int *keyOut, int *codeOut )
{
  HANDLE input_handle = GetStdHandle(STD_INPUT_HANDLE);
  DWORD events = 0;			// how many events took place
  INPUT_RECORD input_record;	// a record of input events
  DWORD input_size = 1;		// how many characters to read

  BOOL peeked = PeekConsoleInput(input_handle, &input_record, input_size, &events);

  if(peeked && input_record.EventType == KEY_EVENT)
  { // PeekConsoleInput succeeded and a key was pressed, so set and return keypress.
    *keyOut = input_record.Event.KeyEvent.wVirtualKeyCode;
    *codeOut = input_record.Event.KeyEvent.dwControlKeyState;
    return MI_OK;
  }
  return MI_EAGAIN;	
}
#endif

/*-------------------------------------------------------------------*/
/* 								     */			
/*-------------------------------------------------------------------*/
char *MiHexdump2Str( unsigned char *memaddr, size_t memlen, 
		      char *buf, int buflen )
{
  int             i,j,byte;
  unsigned char        *ptr;
  char          hexdigits[16] = { '0','1','2','3',
                                    '4','5','6','7',
                                    '8','9','A','B',
                                    'C','D','E','F' };  
  ptr = (unsigned char *)memaddr;
  for( i=0, j=0;
       ( i < (int)memlen ) && ( j < buflen - 4 ); 
       i++ )
  {
    byte = *ptr++;
    buf[j++] = hexdigits[byte >> 4];
    buf[j++] = hexdigits[byte & 15];
    if(( (i+1) & 3 ) == 0 ) buf[j++] = ' ';  
  }
  
  buf[j] = '\0';
  return buf;
}

/*-------------------------------------------------------------------*/
/* 								     */			
/*-------------------------------------------------------------------*/
size_t MiStripTrailingCRLF( char *str, size_t slen ) 
{
  if( ! slen ) slen = strlen(str);
  int i =  slen - 1; 
  
  if( (i>=0) && (str[i] == '\n' ))
  {
    str[i] = '\0'; slen = i;
    i--;
    if((i>=0)&&( str[i] == '\r' )) 
    { 
      str[i] = '\0'; slen = i;
    }  
  }  
  return slen;
}

/*-------------------------------------------------------------------*/
/* 								     */
/*-------------------------------------------------------------------*/
int UtWriteCRLF( FILE *fp )
{
#define CHCR   '\x0d'
#define CHLF   '\x0a'
  if( putc (CHCR,fp) != CHCR ) return MI_ERROR;
  if( putc (CHLF,fp) != CHLF ) return MI_ERROR;
  return 0;
}

/*-------------------------------------------------------------------*/
/* 								     */
/*-------------------------------------------------------------------*/
int UtReadCRLF( FILE *fp )
{
  if( fgetc (fp) != CHCR ) return MI_ERROR;
  if( fgetc (fp) != CHLF ) return MI_ERROR;
  return 0;
}

/*-------------------------------------------------------------------*/
/* 								     */			
/*-------------------------------------------------------------------*/
size_t UtStripTrailingWS( char *str, size_t slen ) 
{
  if( ! slen ) slen = strlen(str);
  int i; 
  
  for(i=slen-1;i>=0;i--)
  {
    if(!((str[i]==' ' || str[i]=='\t'))) break;
  }
  str[i+1] = '\0';
  return i+1;
}


/*--------------------------------------------------------------------*/
/* 								      */			
/*--------------------------------------------------------------------*/
int UtSimpleRegexMatch(char	*string,
			char	*pattern,
			int      accept_sub_string,
			int      case_insensitive)
{
    for (;; pattern++, string++ )
    {
        if (*pattern == '\0')
        {
	    if (accept_sub_string)
		return (TRUE);
            else if (*string == '\0')
                return (TRUE);
            else
                return (FALSE);
        }
        if ((*string == '\0') && (*pattern != '*'))
            return (FALSE);

        if (*pattern == '*')
        {
            pattern++;
            if (*pattern == '\0')
            {
                return (TRUE);
            }
            for (;;)
            {
                if (UtSimpleRegexMatch (string, pattern,
				      accept_sub_string, case_insensitive))
                    return (TRUE);

                if (*string == '\0')
                    return (FALSE);

                string++;
            }
        }

        if (*pattern == '?')
            continue;

        if (*pattern == '\\')
        {
            pattern++;
            if (*pattern == '\0') return (FALSE);
        }

	if (case_insensitive)
	{
	    if (toupper ((int)*pattern) != toupper ((int)*string))
		return (FALSE);
	}
	else
	{
	    if (*pattern != *string)
		return (FALSE);
	}
    }

    return (FALSE);
}


/*--------------------------------------------------------------------*/
/*                                                        	      */
/*--------------------------------------------------------------------*/			           	     	
char *UtGetNextLineFromString( char **textBegin, char *textEnd )
{
  register char  ch,
                  *pEnd = *textBegin, *pBegin = *textBegin;
  register int	   state = 0, len = textEnd - *textBegin;                  
  
  if( ! len ) return NULL;
  while( len-- )
  {
    ch = *pEnd++;
    switch( state )
    {
      case 0:
        if( ch == '\r' ) 
          state = 1;
        if( ch == '\n' )
        {
          *(pEnd-1) = '\0';
          *textBegin = pEnd;
           return pBegin;
        }   
        break;
      case 1:
        if( ch == '\n' ) 
        {
          *(pEnd-1) = *(pEnd-2) = '\0';
          *textBegin = pEnd;
           return pBegin;
        }   
        else
        {
          state = 0;  
        }  
        break;
      default:
        break;              
    }    
  }

  if( pEnd - pBegin > 0 )
  {
    *textBegin = pEnd;
    return pBegin;
  }

  return NULL;      
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtGetTextLine( FILE *fi, char *line, size_t maxLineLen, unsigned opt,
    		   char **pline, int *plen )
{
  int rcThis=0,len=0;
  char *lp=NULL;

  for(;;)
  {
    lp = fgets(line,maxLineLen,fi);
    if(feof(fi)) raiseRc( MI_EOF );   
    if( ! lp ) raiseRc( MI_ERROR );
    len = strlen(line);
    MiStripTrailingCRLF(line,len);
    lp =line;
    if( opt & UT_ALL ) break;
    if(line[0]=='#') continue;
    while(*lp == ' '||*lp=='\t') lp++;
    if(lp[0]=='\0') continue;
    break;
  }

rcCatch:
  *pline = lp;
  if(plen) *plen = len;
  return rcThis;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
char *UtFileSuffix( char *fname )
{
  int i,len=strlen(fname);

  for(i=len-1;i>=0;i--)
  {
    if(fname[i]=='.') return fname+i+1; 
  }
  return fname+len;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
char *UtSetFilenameSuffix( char *fnameCompl, 
    			    char *fname, 
			    char *suffix )
{
  int i,len=strlen(fname), slen = strlen(suffix);
  
  if( fnameCompl != fname )
    strcpy(fnameCompl, fname);

  for(i=len-1;i>=0;i--)
  {
    if(fnameCompl[i]=='.') 
    { 
      fnameCompl[i] = 0; 
      len = i+1;
      break; 
    }
  }
  
  if( len <= 4 || strcmp(fname+len-slen,suffix)!=0 )
  {
    strcat( fnameCompl, suffix );
  }

  return fnameCompl;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
char *MiTranspathCyg2Win( char *fname )
{
  int i;
  for(i=0;i<strlen(fname);i++) 
    if(fname[i]=='/') fname[i]='\\';
  return fname;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
char *UtPathCyg2Win( char *pin )
{
  int i,len;
  char dletter, *chp, *pin0=pin;

  if(((len=strlen(pin)) >=11 ) && 
      ( sscanf(pin,"/cygdrive/%c",&dletter) == 1 ))
  {
    pin+=11; 
    chp = pin0;
    *chp++ = dletter; *chp++=':';
    for(i=0;i<len-11;i++) chp[i] = pin[i];
    len -= 9;
  }

  for(i=0;i<len;i++) 
    if(pin0[i]=='/') pin0[i]='\\';
  pin0[i]= '\0';

  return pin0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
char *UtPathWin2Cyg( char *pin, int maxlen )
{
  int i,len,len2,rcThis=0;
  char dletter, *pin0=pin;

  USE(rcThis);

  if(((len=strlen(pin)) >=3 ) && 
      ( sscanf(pin,"%c:\\",&dletter) == 1 ))
  {
    len2 = len+9;
    if(len>maxlen) 
      TRCERRR(("string buffer too short for file path '%s'",pin),
	  MI_EOVERFLOW);
    for(i=0;i<len-2;i++) pin[len2-i-1] = pin[len-i-1];
    memcpy(pin,"/cygdrive/",10);
    pin[10]=dletter;
    pin[len2] = '\0';
  }
  else len2=len;

  for(i=0;i<len2;i++) 
    if(pin0[i]=='\\') pin0[i]='/';
rcCatch:
  return pin0;
}



/*-------------------------------------------------------------------*/
/* 								     */			
/*-------------------------------------------------------------------*/
int UtShellExec( char *command, char *params, unsigned options,
    		 int *pRetCode )
{

  #ifdef MIMP_ON_WINDOWS

  int rcThis = 0;  

  TRC(("MiShellExecute: cmd='%s' params='%s'\n",command,params));

  if(pRetCode) *pRetCode = -1;
 
  if( options & UT_WAIT )
  {
    SHELLEXECUTEINFO 	si;
    BOOL sysok;
    DWORD dwCode;
    
    memset( &si, 0, sizeof(si));
    si.lpVerb = "open";
    si.cbSize = sizeof(si);
    si.lpFile = command;
    si.lpParameters = params;
    si.nShow = SW_SHOWNOACTIVATE;
    si.fMask = SEE_MASK_NOCLOSEPROCESS;

    sysok = (int)ShellExecuteEx( &si);
    if( ! sysok )
    {
      TRCERRR(("ShellExecute('%s','%s') failed. hInstApp=%p, errno=%d\n",
	    command, params, (void *)si.hInstApp,
	    (int)GetLastError()),
	    MI_ERROR );    
    }    

    // Wait for process to finish
    TRC1(("waiting for process to finish ...\n"));
    WaitForSingleObject(si.hProcess, INFINITE);
    GetExitCodeProcess(si.hProcess, &dwCode);

    if(pRetCode) *pRetCode = dwCode;

    if(dwCode==STILL_ACTIVE)
    {
      TRCERRR(("Process still active\n"),MI_EINTERNAL);
    }
    
    TRC1(("-> Finished with code %d\n",(int)dwCode)); 
  }
  else
  {
    void *sysrc = (void*)
      			ShellExecute( NULL,	// hwnd 
			      "open", 	// lpOperation
			      command,  // lpFile
			      params,   // lpParameters
			      NULL,	// lpDirectory
			      SW_SHOWNOACTIVATE // nShowCommand
			      );
    if( (int)(ptrdiff_t)sysrc <= 32 )
    {
      TRCERRR(("ShellExecute('%s','%s') failed. sysrc=%lu\n",
	    command, params, sysrc),MI_ERROR );    
    }

    TRC(("MiShellExecute: return code =%lu\n",sysrc));
  }
   
rcCatch:        
  return rcThis;      

#else

  UT_ASSERT_NR(FALSE);
  return MI_ERROR;

#endif
}

#if 0  // deprecated because too dangerous wrt accidentally deleting whole dirs
/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtDeleteDirectory(char *dirnameIn)
{
  int rcThis=0,
      noRecycleBin=TRUE,
      len = strlen(lpszDir);
  char *pszFrom = NULL;
  
  ALLOCARR(pszFrom,len+2);
  strcpy(pszFrom,lpszDir);
  pszFrom[len] = 0;
  pszFrom[len+1] = 0;

  UtPathCyg2Win(pszFrom);
  
  TRC(("UtDeleteDirectory '%s'\n",pszFrom));

  SHFILEOPSTRUCT fileop;
  fileop.hwnd   = NULL;    // no status display
  fileop.wFunc  = FO_DELETE;  // delete operation
  fileop.pFrom  = pszFrom;  // source file name as double null terminated string
  fileop.pTo    = NULL;    // no destination needed
  fileop.fFlags = FOF_NOCONFIRMATION|FOF_SILENT;  // do not prompt the user
  
  if(!noRecycleBin)
    fileop.fFlags |= FOF_ALLOWUNDO;

  fileop.fAnyOperationsAborted = FALSE;
  fileop.lpszProgressTitle     = NULL;
  fileop.hNameMappings         = NULL;

  {
	    TRCP(("Delete directory '%s' (y/n) ?\n",pszFrom));
	    fflush(stdout);
	    if(getch()!='y') 
	    {
	      TRCWARN(("operation aborted by user\n"));
	      raiseRc(MI_ECANCEL);
	    }
  }

  int ret = SHFileOperation(&fileop);
  if(ret) 
  {
    TRCERRR(("UtDeleteDirectory: SHFileOperation failed for '%s', errno=%d\n",
	pszFrom,(int)GetLastError()),MI_ERROR);
  }

rcCatch:
  FREEMEM(pszFrom);
  return rcThis;
}
#endif

/*-------------------------------------------------------------------*/
/* 								     */			
/*-------------------------------------------------------------------*/
int UtGetCurDir( char *pbuf, size_t *plen )
{
  int 	        rcThis = MI_OK;  

  #ifdef MIMP_ON_WINDOWS

  DWORD		len = 0;
  LPTSTR        lpBuffer = pbuf;
  len = GetCurrentDirectory( (DWORD)((*plen)-1), lpBuffer );
  if(!len)
  {
    TRCERRR(("GetCurrentDirectory() failed. errno=%d\n",
      (int)GetLastError()),MI_ERROR );
  }
  pbuf[len]='\0';

  #else

  size_t len=0;
  size_t lenMax = (*plen)-1;
  char *pbuf_ = getcwd(pbuf, lenMax);
  if(!pbuf_)
  {
    TRCERRR(("getcwd() failed. errno=%d\n",errno),MI_ERROR );
  }
  len = strlen(pbuf);

  #endif


rcCatch:  

  *plen = len;
      
  return rcThis;      
}

/*-------------------------------------------------------------------*/
/* 								     */			
/*-------------------------------------------------------------------*/
int UtSetCurDir( char *pbuf )
{
  int 	        rcThis = MI_OK;  

  #ifdef MIMP_ON_WINDOWS

  BOOL		ok;
  LPTSTR        lpBuffer = pbuf;

  ok = SetCurrentDirectory( lpBuffer );
  if(!ok)
  {
    TRCERRR(("SetCurrentDirectory('%s') failed. errno=%d\n",
      pbuf,(int)GetLastError()),MI_ERROR );
  }

  #else

  int result = chdir(pbuf);
  if (result != 0) 
  {
    TRCERRR(("chdir('%s') failed. errno=%d\n",
      pbuf,errno),MI_ERROR );
  }

  #endif

rcCatch:  
 
  return rcThis;      
}

/*-------------------------------------------------------------------*/
/* 								     */			
/*-------------------------------------------------------------------*/
int UtDiffFiles( char *fn1, char *fn2 )
{
  int 	        sysrc, rcThis = 0;  
  char		syscmd[2*UT_MAXPATHLEN+100];

  sprintf(syscmd,"diff %s %s\n",fn1,fn2);
  sysrc = system(syscmd);
  TRC(("system command: '%s' -> %d\n",syscmd,sysrc));
  if(sysrc!=0) raiseRc(MI_EDIFFER); // TODO: check other errors

rcCatch:        
  TRC(("UtDiffFiles '%s' '%s' -> %d\n",fn1,fn2,rcThis));
  return rcThis;      
}

/*-------------------------------------------------------------------*/
/* 								     */			
/*-------------------------------------------------------------------*/
int UtDeleteFile( char *fname )
{
  int rcThis=0,sysrc=0;

  #ifdef MIMP_ON_WINDOWS

  sysrc = (int)DeleteFile( fname );
  if( sysrc == 0 )
  {
    int syerr=(int)GetLastError();

    if(syerr==ERROR_FILE_NOT_FOUND)
    {
      TRC(("UtDeleteFile: '%s' does not exist.\n",fname));
      raiseRc(MI_ENOTFOUND);
    }
    else if(syerr==ERROR_ACCESS_DENIED)
    {
      TRCERRR(("UtDeleteFile: '%s' access denied.\n",fname),
	  MI_ERROR);
    }
    else
    {
      TRCERRR(("DeleteFile('%s') failed. errno=%d\n",
	    fname, GetLastError()),MI_ERROR );    
    }
  }  			
 
  #else

  sysrc = unlink( fname );
  if(sysrc)
  {
    TRCERRR(("unlink('%s') failed. errno=%d\n",
	  fname,errno),MI_ERROR );    
  }

  #endif

rcCatch:        
  TRC(("UtDeleteFile '%s' -> %d\n",fname,rcThis));
  return rcThis;      
}

/*-------------------------------------------------------------------*/
/* 								     */			
/*-------------------------------------------------------------------*/
int UtCopyFile( char *src, char *dst, int doFailIfExists )
{
  int rcThis = 0;

  #ifdef MIMP_ON_WINDOWS
 
  int sysrc = (int)CopyFile( src, dst, doFailIfExists );
  if( sysrc == 0 )
  {
    int syerr=(int)GetLastError();

    if(syerr==ERROR_FILE_NOT_FOUND)
    {
      TRC(("UtCopyFile() source file does not exist\n"));
      raiseRc( MI_ENOTFOUND );
    }
    else
    {
      TRCERRR(("CopyFile('%s','%s') failed. errno=%d\n",
	    src, dst, syerr),MI_ERROR );    
    }
  }  			
   
  #else

  {
    char cmd[4096];
    #if 0
    sprintf( cmd, "[ -e src ]; then cp -p %s \'%s\' \'%s\'; fi", doFailIfExists?"":"-f", src, dst);
    #else
    sprintf( cmd, "cp -p %s \'%s\' \'%s\'", doFailIfExists?"":"-f", src, dst);
    #endif

    int sysrc = system( cmd );
    if(sysrc)
    {
      TRCERRR(("system(%s) failed sysrc=%d. errno=%d\n",cmd,sysrc,errno),MI_ERROR );    
    }
  }

  #endif

rcCatch:        
  TRC(("UtCopyFile '%s' -> '%s'\n",src,dst));
  return rcThis;      
}

/*-------------------------------------------------------------------*/
/* 								     */
/*-------------------------------------------------------------------*/
char *UtReadFileToString( char *path )
{
  int rcThis=0,ch,i;
  FILE *fp=NULL;
  char *buf=NULL;
  size_t fsize;

  fp = fopen(path,"r");
  if(!fp)
  {
    TRCERRR(("can't open log-file '%s'\n",path),MI_ERROR);
  }

  fseek( fp, 0, SEEK_END );  /* get file size */
  fsize = ftell( fp );
  rewind(fp);

  buf = MM_malloc(fsize+1,MM_DUMMYTAG);
  if(!buf)
  {
    TRCERRR(("MM_malloc(%ld) failed\n",(long)(fsize+1)),MI_ENOMEM);
  }

  for(i=0;i<fsize;)
  {
    ch=fgetc(fp);
    if(ch==EOF) break;
    buf[i++]=ch;
  }
  if(i!=fsize)
  {
    TRCERRR(("error reading file '%s'\n",path),MI_EIO);
  }
  buf[i]=0;

rcCatch:
  if(fp) fclose(fp);
  if(rcThis)
  {
    if(buf) MM_free(buf);
    buf=NULL;
  }
  return buf;
}

/*-------------------------------------------------------------------*/
/* 								     */			
/*-------------------------------------------------------------------*/
int UtCheckFileExists( char *path )
{
  int exists = FALSE;
  FILE *fp = fopen(path,"r");
  if( fp ) 
  {
    exists=TRUE;
    fclose(fp);
  }
  return exists ? MI_OK : MI_ENOTFOUND;
}

#if 0
int DeleteDirectoryAndAllSubfolders(char *dir)
{
    SHFILEOPSTRUCT fo = {0};

     fo.wFunc = FO_DELETE;
     fo.pFrom = chbuf;
     fo.fFlags = FOF_SILENT | FOF_NOERRORUI | FOF_NOCONFIRMATION;



#if 0
    WCHAR szDir[MAX_PATH+1];  // +1 for the double null terminate

    StringCchCopy(szDir, MAX_PATH, wzDirectory);
    int len = lstrlenW(szDir);
    szDir[len+1] = 0; // double null terminate for SHFileOperation

    // delete the folder and everything inside
    fos.wFunc = FO_DELETE;
    fos.pFrom = szDir;
    fos.fFlags = FOF_NO_UI;
#endif
    return SHFileOperation( &fos );
}
#endif

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtCreatePath( char *pathIn )
{
  char *path2=NULL,*path;
  int 	rcThis = 0,i,posOld,len,alreadyExists;  

  char sep;

  #ifdef MIMP_ON_WINDOWS
  path2 = MM_strdup(UtPathCyg2Win(pathIn));
  sep = '\\';
  #else
  path2 = MM_strdup(pathIn);
  sep = '/';
  #endif

/*  UtParseFilePath(path2,&path,NULL,NULL);*/
  path = path2;
  len = strlen(path);

  TRC(("UtCreatePath: '%s'\n",path));

  alreadyExists = TRUE;
  for(i=0,posOld=0;i<=len;i++)
  {
    if(((path[i]==sep || path[i]=='\0')) && (i - posOld > 1))
    {
      // skip root part 
      if(i>0 && path[i-1]==':') continue;

      path[i] = '\0';

      TRC3(("CreateDirectory('%s')\n",path));

      #ifdef MIMP_ON_WINDOWS

      int sysrc = (int)CreateDirectory( path, NULL );
      if( sysrc == 0 )
      {
	int eno = GetLastError();
	if(eno!=ERROR_ALREADY_EXISTS)
	{
	  TRCERRR(("CreateDirectory('%s') failed. errno=%d\n",
		path, eno),MI_ERROR );    
	}
      }
      else
      {
	alreadyExists = FALSE;
      }

      #else
      
      int status = 0;
      {
	struct stat     st;
	if (stat(path, &st) != 0)
	{
	  /* directory does not exist. eexist for race condition */
	  alreadyExists = FALSE;
	  mode_t mode = 0777; 
	  if (mkdir(path, mode) != 0 && errno != EEXIST)
	      status = -1;
	}
	else if (!S_ISDIR(st.st_mode))
	{
	  errno = ENOTDIR;
	  status = -1;
	}
      }

      if(status!=0)
      {
	TRCERRR(("CreateDirectory('%s') failed. errno=%d\n",
	      path, errno),MI_ERROR );    
      }

      #endif

      path[i] = sep;
    }
  }
  if(path2) MM_free(path2);
  if(alreadyExists) rcThis=MI_EEXISTS;
rcCatch:
  return rcThis;
}


/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtIsAbsFilePath( char *fullpath )
{
  char path[UT_MAXPATHLEN+1];
  strcpy(path,fullpath);
  if(strlen(path)&&(path[0]=='/'||path[0]=='\\')) return TRUE;
  char *dir=NULL,*file=NULL,*suff=NULL;
  (void)UtParseFilePath(path,&dir,&file,&suff);
  if(strlen(dir)==0) return FALSE;
  if(strlen(dir)>=2&&dir[1]==':') return TRUE;
  return FALSE;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtParseFilePath( char *fullPath, 
    		     char **pDir, char **pFile, char **pSuffix )
{
  int j,len=strlen(fullPath);
  char *chp = fullPath;

  for(j=strlen(chp)-1;j>=0 && chp[j]!='/'&&chp[j]!='\\';j--);
  if(j>=0) 
  {
    chp[j]=0; 
    if(pDir) *pDir = chp;
  }
  else 
  {
    if(pDir) *pDir = chp+len;
  }
  j++;
  if(pFile) *pFile = chp+j;

  if(pSuffix && pFile) *pSuffix = chp = UtFileSuffix(*pFile);
  if(pFile && chp && strlen(*pFile)>0 && *(chp-1)=='.') *(chp-1) = '\0';

  return 0;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtStripLastDirSep( char *path )
{
  int l = UtStripTrailingWS(path,strlen(path));
  if(l>0 && (path[l-1]=='/'||path[l-1]=='\\'))
  {
    l--;
    path[l] = '\0';
  }
  return l;
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
#ifdef MI_ON_NATWIN

#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

struct timezone 
{
  int  tz_minuteswest; /* minutes W of Greenwich */
  int  tz_dsttime;     /* type of dst correction */
};

int gettimeofday(struct timeval *tv, struct timezone *tz)
{
  FILETIME ft;
  unsigned __int64 tmpres = 0;
  static int tzflag;

  if (NULL != tv)
  {
    GetSystemTimeAsFileTime(&ft);

    tmpres |= ft.dwHighDateTime;
    tmpres <<= 32;
    tmpres |= ft.dwLowDateTime;

    /*converting file time to unix epoch*/
    tmpres -= DELTA_EPOCH_IN_MICROSECS; 
    tmpres /= 10;  /*convert into microseconds*/
    tv->tv_sec = (long)(tmpres / 1000000UL);
    tv->tv_usec = (long)(tmpres % 1000000UL);
  }

  if (NULL != tz)
  {
    if (!tzflag)
    {
      _tzset();
      tzflag++;
    }
    tz->tz_minuteswest = _timezone / 60;
    tz->tz_dsttime = _daylight;
  }

  return 0;
}
#endif

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtHTimerStart(UtHTimer *tm)
{
  #ifdef MIMP_ON_WINDOWS

  int rcThis=0;
  LARGE_INTEGER lpCount;
  BOOL		ok;

  if(!hrTimerFreqValid)
  {
    LARGE_INTEGER lpFrequency; 
    ok =  QueryPerformanceFrequency( &lpFrequency );
    if(!ok)
    {
      TRCERRR(("QueryPerformanceFrequency() failed. errno=%d\n",
       (int)GetLastError()),MI_ERROR );
    }
    hrTimerFreq = (double)lpFrequency.QuadPart;
    TRC1(("high resolution timer frequency = %g\n",
	  hrTimerFreq));
    hrTimerFreqValid = TRUE;
  }

  ok = QueryPerformanceCounter(&lpCount);
  if(!ok)
  {
    TRCERRR(("QueryPerformanceCounter() failed. errno=%d\n",
	  (int)GetLastError()),MI_ERROR );
  }

  tm->hrStart = (double)lpCount.QuadPart;

rcCatch:
  return rcThis;

  #else

  UT_ASSERT_NR(FALSE);
  return MI_ERROR;

  #endif
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
int UtHTimerStop(UtHTimer *tm)
{

  #ifdef MIMP_ON_WINDOWS
  
  int rcThis=0;
  LARGE_INTEGER lpCount; 
  BOOL		ok;

  ok = QueryPerformanceCounter(&lpCount);
  if(!ok)
  {
    TRCERRR(("QueryPerformanceCounter() failed. errno=%d\n",
	  (int)GetLastError()),MI_ERROR );
  }

  tm->hrStop = (double)lpCount.QuadPart;

  tm->hrDiffSec = ( tm->hrStop - tm->hrStart ) / hrTimerFreq;

rcCatch:
  return rcThis;

  #else

  UT_ASSERT_NR(FALSE);
  return MI_ERROR;

  #endif
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
void UtDeleteFSM( UtFSM *fsm )
{
  if(fsm)
  {
    FREEMEM(fsm->data);
    FREEMEM(fsm->stack);
    FREEMEM(fsm);
  }
}

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/
UtFSM *UtCreateFSM( size_t maxBlocks, size_t blockSize )
{
  int rcThis=0;
  size_t i,memalign=8; /* TODO */
  UtFSM *fsm=NULL;
  blockSize = ALIGN( blockSize, memalign );
  ALLOCOBJ0( fsm );
  ALLOCMEM( fsm->data, blockSize*maxBlocks );
  ALLOCARR( fsm->stack, maxBlocks );
  fsm->sp = 0;
  //#pragma omp parallel for
  for(i=0;i<maxBlocks;i++)
  {
    fsm->stack[i] = (char*)fsm->data+blockSize*(maxBlocks-i-1);
//    UT_FFREE( fsm, (char*)fsm->data+blockSize*i );
  }
  fsm->sp = maxBlocks;

/*rcCatch:*/
  if(rcThis)
  {
    UtDeleteFSM( fsm );
    fsm = NULL;
  }
  return fsm;
}


