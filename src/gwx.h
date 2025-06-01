/******************************************************************************
 *
 * Copyright 2025 Bernhard Braun 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 ******************************************************************************/

#ifndef GWX_H
#define GWX_H

#ifdef __cplusplus
extern "C" {
#endif

extern int globGwxImpl;

int GWXsavebmp(int NB, char *FN, int l);
int GWXline2(float X, float Y);  
int GWXmove2(float X, float Y);
int GWXcapvec(float X1, float Y1, float *X2, float *Y2, char *TEXT, int l);
int GWXcaplin(float *X1, float *Y1, float *X2, float *Y2, char *TEXT, int l);
int GWXclose(int NW);
int GWXclear(int K);
int GWXrefresh(void);
int GWXpolylin(float *POINTS, int N);
int GWXpolygon(float *POINTS, int N, int MF);  
int GWXquitx(int MQ);
int GWXsetogn(int ISRC, int IDST);
int GWXellipse(float X1, float Y1, float X2, float Y2);  
int GWXsrect(float X1, float Y1, float X2, float Y2, int K);
int GWXsetpxl(float X, float Y, int K);
int GWXline(float X1, float Y1, float X2, float Y2);
int GWXrect(float X1, float Y1, float X2, float Y2);
int GWXmakebmp(int NB, int IW, int IH, int IBC, int *IBITS);
int GWXquit(void);
int GWXerase(int N, int LRF);
int GWXclipimg(float X1, float Y1, float X2, float Y2);
int GWXshowwn(int NW, int IS);
int GWXinput(char *TITLE, int l1, char *TXT, int l2);
int GWXarrange(int M);
int GWXgetpen(int *IPC, int *IPS, int *IPW, int *MX);
int GWXcolor(int K, int IDST);
int GWXopenx(int NW, int IW, int IH, int IFC, int IBC, int M, char *FN, int l);
int GWXcappnt(float *X, float *Y, char *TEXT, int l);
int GWXcaprect(float *X1, float *Y1, float *X2, float *Y2, char *TEXT, int l);
int GWXdelbmp(int NM);
int GWXevent(float *X, float *Y);
int GWXfiledlg(char *FN, int l);
int GWXgetwn(float *X1, float *Y1, float *X2, float *Y2);
int GWXkybrd(int *ICH, int *NCNT, int *IFLG, int M);
int GWXputbmp(int NB, float X, float Y, int IBK);
int GWXputtxt(float X, float Y, char *TXT, int l);
int GWXmsgbox(char *TXT, int l);
int GWXselect(int NW);
int GWXsetbmp(int NB, float W, float H, int MX, int ITR, int IOF);
int GWXsetpen(int IPC, int IPS, int IPW, int MX);
int GWXsettxt(float H, float A, int IO, int K, int KB, char *FACE, int l);
int GWXkrgb(int IR, int IG, int IB);
int GWXinitx(int IRB, int IX, int IY, int IW, int IH, int MA, int MM, 
		int MZ, int ND);
int GWXsetmsg(char *TXT, int l);

#ifdef __cplusplus
}
#endif

#endif /* GWX_H */

