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
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#define MIMP_NOREDEFINE_GWX
#include "globdef.h"
#include "util.h"
#include "gwx.h"

#define GW_OK 1

int globGwxImpl=1;

//#ifdef MIMP_ON_LINUX

int GWsavebmp(int NB, char *FN, int l) { return 0; }
int GWmove2(float X, float Y) { return 0; }
int GWline2(float X, float Y) { return 0; }
int GWcapvec(float X1, float Y1, float *X2, float *Y2, char *TEXT, int l) { return 0; }
int GWcaplin(float *X1, float *Y1, float *X2, float *Y2, char *TEXT, int l) { return 0; }
int GWclose(int NW) { return 0; }
int GWclear(int K) { return 0; }
int GWrefresh(void) { return 0; }
int GWpolylin(float *POINTS, int N) { return 0; }
int GWpolygon(float *POINTS, int N, int MF) { return 0; }
int GWquitx(int MQ) { return 0; }
int GWsetogn(int ISRC, int IDST) { return 0; }
int GWellipse(float X1, float Y1, float X2, float Y2) { return 0; }
int GWsrect(float X1, float Y1, float X2, float Y2, int K) { return 0; }
int GWsetpxl(float X, float Y, int K) { return 0; }
int GWline(float X1, float Y1, float X2, float Y2) { return 0; }
int GWrect(float X1, float Y1, float X2, float Y2) { return 0; }
int GWmakebmp(int NB, int IW, int IH, int IBC, int *IBITS) { return 0; }
int GWquit(void) { return 0; }
int GWerase(int N, int LRF) { return 0; }
int GWclipimg(float X1, float Y1, float X2, float Y2) { return 0; }
int GWshowwn(int NW, int IS) { return 0; }
int GWinput(char *TITLE, int l1, char *TXT, int l2) { return 0; }
int GWarrange(int M) { return 0; }
int GWgetpen(int *IPC, int *IPS, int *IPW, int *MX) { return 0; }
int GWcolor(int K, int IDST) { return 0; }
int GWopenx(int NW, int IW, int IH, int IFC, int IBC, int M, char *FN, int l) { return 0; }
int GWcappnt(float *X, float *Y, char *TEXT, int l) { return 0; }
int GWcaprect(float *X1, float *Y1, float *X2, float *Y2, char *TEXT, int l) { return 0; }
int GWdelbmp(int NM) { return 0; }
int GWevent(float *X, float *Y) { *X = 0; *Y = 0; return 0; }
int GWfiledlg(char *FN, int l) { return 0; }
int GWgetwn(float *X1, float *Y1, float *X2, float *Y2) { return 0; }
int GWkybrd(int *ICH, int *NCNT, int *IFLG, int M) { return 0; }
int GWputbmp(int NB, float X, float Y, int IBK) { return 0; }
int GWputtxt(float X, float Y, char *TXT, int l) { return 0; }
int GWmsgbox(char *TXT, int l) { return 0; }
int GWselect(int NW) { return 0; }
int GWsetbmp(int NB, float W, float H, int MX, int ITR, int IOF) { return 0; }
int GWsetpen(int IPC, int IPS, int IPW, int MX) { return 0; }
int GWsettxt(float H, float A, int IO, int K, int KB, char *FACE, int l) { return 0; }
int GWkrgb(int IR, int IG, int IB) { return 0; }
int GWinitx(int IRB, int IX, int IY, int IW, int IH, int MA, int MM, 
		int MZ, int ND) { return 0; }
int GWsetmsg(char *TXT, int l) { return 0; }

//#endif

int GWXsavebmp(int NB, char *FN, int l)
{
  return globGwxImpl ? GW_OK :  GWsavebmp( NB, FN,  l);
}

int GWXline2(float X, float Y)
{
  return globGwxImpl ? GW_OK :  GWline2( X,  Y);
}

int GWXmove2(float X, float Y)
{
  return globGwxImpl ? GW_OK :  GWmove2( X,  Y);
}

int GWXcapvec(float X1, float Y1, float *X2, float *Y2, char *TEXT, int l)
{
  if(globGwxImpl)
  {
    return 0;
  }
  else
  {
    return GWcapvec( X1,  Y1, X2, Y2, TEXT,  l);
  }  
}

int GWXcaplin(float *X1, float *Y1, float *X2, float *Y2, char *TEXT, int l)
{
  if(globGwxImpl)
  {
    return 0;
  }
  else
  {
    return GWcaplin(X1, Y1, X2, Y2, TEXT,  l);
  }  
}

int GWXclose(int NW)
{
  return globGwxImpl ? GW_OK : GWclose( NW);
}

int GWXclear(int K)
{
  return globGwxImpl ? GW_OK : GWclear(K);
}

int GWXrefresh(void)
{
  return globGwxImpl ? GW_OK : GWrefresh();
}

int GWXpolylin(float *POINTS, int N)
{
  return globGwxImpl ? GW_OK : GWpolylin( POINTS, N );
}

int GWXpolygon(float *POINTS, int N, int MF)
{
  return globGwxImpl ? GW_OK : GWpolygon(POINTS,  N,  MF);
}

int GWXquitx(int MQ)
{
  return globGwxImpl ? GW_OK : GWquitx(MQ);
}

int GWXsetogn(int ISRC, int IDST)
{
  return globGwxImpl ? GW_OK :  GWsetogn( ISRC,  IDST);
}

int GWXellipse(float X1, float Y1, float X2, float Y2)
{
  return globGwxImpl ? GW_OK :  GWellipse( X1,  Y1,  X2,  Y2);
}

int GWXsrect(float X1, float Y1, float X2, float Y2, int K)
{
  return globGwxImpl ? GW_OK : GWsrect( X1,  Y1,  X2,  Y2,  K);
}

int GWXsetpxl(float X, float Y, int K)
{
  return globGwxImpl ? GW_OK :  GWsetpxl( X,  Y,  K);
}

int GWXline(float X1, float Y1, float X2, float Y2)
{
  return globGwxImpl ? GW_OK :  GWline( X1,  Y1,  X2,  Y2);
}

int GWXrect(float X1, float Y1, float X2, float Y2)
{
  return globGwxImpl ? GW_OK :  GWrect( X1,  Y1,  X2,  Y2);
}

int GWXmakebmp(int NB, int IW, int IH, int IBC, int *IBITS)
{
  return globGwxImpl ? 1 :  GWmakebmp( NB,  IW,  IH,  IBC, IBITS); 
}

int GWXquit(void)
{
  return globGwxImpl ? GW_OK : GWquit();
}

int GWXerase(int N, int LRF)
{
  return globGwxImpl ? GW_OK : GWerase( N,  LRF);
}

int GWXclipimg(float X1, float Y1, float X2, float Y2)
{
  return globGwxImpl ? GW_OK : GWclipimg( X1,  Y1,  X2,  Y2);
}

int GWXshowwn(int NW, int IS)
{
  return globGwxImpl ? GW_OK : GWshowwn( NW,  IS);
}

int GWXinput(char *TITLE, int l1, char *TXT, int l2)
{
  if(globGwxImpl)
  {
    *TXT=0;
    return 0;
  }
  else
  {
    return GWinput(TITLE,  l1, TXT,  l2);
  }
}

int GWXarrange(int M)
{
  return globGwxImpl ? GW_OK : GWarrange( M);
}

int GWXgetpen(int *IPC, int *IPS, int *IPW, int *MX)
{
  if(globGwxImpl)
  {
    *IPC = 0;
    *IPS = 0;
    *IPW = 1;
    *MX = 0;
    return GW_OK;
  }
  else
  {
    return GWgetpen(IPC, IPS, IPW, MX);
  }
}

int GWXcolor(int K, int IDST)
{
  return globGwxImpl ? GW_OK : GWcolor( K,  IDST);
}

int GWXopenx(int NW, int IW, int IH, int IFC, int IBC, int M, char *FN, int l)
{
  return globGwxImpl ? GW_OK : 
	 GWopenx( NW,  IW,  IH,  IFC,  IBC,  M, FN,  l);
}


int GWXcappnt(float *X, float *Y, char *TEXT, int l)
{
  if(globGwxImpl)
  {
    return 0;
  }
  else
  {
    return GWcappnt(X, Y, TEXT, l);
  }
}


int GWXcaprect(float *X1, float *Y1, float *X2, float *Y2, char *TEXT, int l)
{
  if(globGwxImpl)
  {
    return 0;
  }
  else
  {
    return GWcaprect(X1, Y1, X2, Y2, TEXT, l);
  }
}

int GWXdelbmp(int NM)
{
  return globGwxImpl ? GW_OK : GWdelbmp( NM);
}

int GWXevent(float *X, float *Y)
{
  if(globGwxImpl)
  {
    *X = 0;
    *Y = 0;
    return GW_OK;
  }
  else
  {
    return GWevent(X, Y);
  }
}

int GWXfiledlg(char *FN, int l)
{
  if(globGwxImpl)
  {
    return GW_OK;
  }
  else
  {
    return GWfiledlg(FN,  l);
  }
}

int GWXgetwn(float *X1, float *Y1, float *X2, float *Y2)
{
  if(globGwxImpl)
  {
    *X1=100;
    *Y1=100;
    *X2=200;
    *Y2=200;
    return GW_OK;
  }
  else return GWgetwn(X1, Y1, X2, Y2);
}

int GWXkybrd(int *ICH, int *NCNT, int *IFLG, int M)
{
  if(globGwxImpl)
  {
    *ICH =0;
    *NCNT = 0;
    *IFLG = 0;
    return 0;
  }
  else
  {
    return GWkybrd(ICH, NCNT, IFLG,  M);
  }
}

int GWXputbmp(int NB, float X, float Y, int IBK)
{
  return globGwxImpl ? GW_OK : GWputbmp( NB,  X,  Y,  IBK);
}

int GWXputtxt(float X, float Y, char *TXT, int l)
{
  return globGwxImpl ? GW_OK : GWputtxt( X,  Y, TXT,  l);
}

int GWXmsgbox(char *TXT, int l)
{
  if(globGwxImpl)
  {
    char chbuf[1000];
    if(l==0) l=strlen(TXT);
    int len = MIN(ARRAY_LENGTH(chbuf)-1,l);
    memcpy(chbuf,TXT,len);
    chbuf[len]=0;
    printf("%s (yes/no)",chbuf);fflush(stdout);
    int ch=getch();
    if(ch=='y') return 1;
    else if(ch=='n') return -1;
    else return 0;
  }
  else
  {
    return GWmsgbox(TXT,l);
  }
}

int GWXselect(int NW)
{
  return globGwxImpl ? GW_OK : GWselect(NW);
}

int GWXsetbmp(int NB, float W, float H, int MX, int ITR, int IOF)
{
  return globGwxImpl ? GW_OK : GWsetbmp( NB,  W,  H,  MX,  ITR,  IOF);
}

int GWXsetpen(int IPC, int IPS, int IPW, int MX)
{
  return globGwxImpl ? GW_OK :  GWsetpen( IPC,  IPS,  IPW,  MX);
}

int GWXsettxt(float H, float A, int IO, int K, int KB, char *FACE, int l)
{
  return globGwxImpl ? GW_OK : GWsettxt( H,  A,  IO,  K,  KB,  FACE,  l);

}

int GWXkrgb(int IR, int IG, int IB)
{
  return globGwxImpl ? 0 : GWkrgb( IR,  IG,  IB);
}  

int GWXinitx(int IRB, int IX, int IY, int IW, int IH, int MA, int MM, 
		int MZ, int ND)
{
  return globGwxImpl ? GW_OK : 
     GWinitx( IRB,  IX,  IY,  IW,  IH,  MA,  MM, 
		 MZ,  ND);

}

int GWXsetmsg(char *TXT, int l)
{
  if(globGwxImpl==0)
  {
    return GWsetmsg(TXT,  l);
  }
  else
  {
    char chbuf[1000];
    if(l==0) l=strlen(TXT);
    int len = MIN(ARRAY_LENGTH(chbuf)-1,l);
    memcpy(chbuf,TXT,len);
    chbuf[len]=0;
#if 1
    {
      int i,len=strlen(TXT);
      for(i=0;i<len;i++) putchar('\b');
      for(i=0;i<len;i++) putchar(TXT[i]); 
      fflush(stdout);
    }
#else
    TRCP(("GWXsetmsg ==> '%s'\n",chbuf));
#endif
    return GW_OK;
  }
}

