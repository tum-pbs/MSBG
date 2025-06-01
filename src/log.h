/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#ifndef TRACE_H
#define TRACE_H

#ifdef __cplusplus
extern "C" {
#endif

int LogErr( const char *format, ...);
int LogWarn( const char *format, ...);
int LogInfo( const char *format, ...);
int LogLog( const char *format, ...);
int LogPrnLog( const char *format, ...);
int LogThrowError( const char*file, int line );
void LogLock(void);
void LogUnlock(void);

extern int LogLevel;
extern char LogSaveLoc[];

#define TRC_LEVEL LogLevel

#define TRC(x) { if(unlikely( LogLevel >= 2 )) { LogLog x; }}
#define TRC3(x) { if(unlikely( LogLevel >= 3 )) { LogLog x; }}

//#define TRC1(x) { if( LogLevel >= 1 ) { LogPrnLog x; }}  /* TODO */
#define TRC1(x) { if( LogLevel >= 1 ) { LogLog x; }}  /* TODO */

#define TRCP(x) { if( LogLevel >= 1 ) { LogPrnLog x; }}
#define TRCLOG(x) { if( LogLevel >= 1 ) { LogLog x; }}
#define TRCL(tl_,x) { if( LogLevel >= tl_ ) { LogLog x; }}

#define TRCWARN(x) \
{ if( LogLevel >= 1 ) { \
  LogLock(); \
  sprintf(LogSaveLoc,"  <%-36.36s %d>",__FILE__,__LINE__); LogWarn x;  \
  LogUnlock(); }}

#define TRCINFO(x) \
{ if( LogLevel >= 1 ) { \
  LogLock(); \
  sprintf(LogSaveLoc,"  <%-36.36s %d>",__FILE__,__LINE__); LogInfo x;  \
  LogUnlock(); }}

#define TRCERR(x) \
{ if( LogLevel >= 1 ) { \
  LogLock(); \
  sprintf(LogSaveLoc,"  <%-36.36s %d>",__FILE__,__LINE__); LogErr x;  \
  LogUnlock(); }}

#define TRCERRR(x,rc) \
{ if( LogLevel >= 1 ) { \
  LogLock(); \
  sprintf(LogSaveLoc,"  <%-36.36s %d>",__FILE__,__LINE__); LogErr x; \
  LogUnlock(); \
  } \
  raiseRc( rc ) }

#ifdef __cplusplus
#define TRCERRT(x,rc) \
{ if( LogLevel >= 1 ) { \
  LogLock(); \
  sprintf(LogSaveLoc,"  <%-36.36s %d>",__FILE__,__LINE__); LogErr x; \
  LogUnlock(); } \
  throw( rc ); }
#endif

#ifdef __cplusplus
}
#endif

#endif /* TRACE_H */

