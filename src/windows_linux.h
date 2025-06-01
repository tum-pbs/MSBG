#ifdef MIMP_ON_WINDOWS
#include <windows.h>
#else

// 
// Stub for use with MIMP_ON_LINUX 
//
//
#ifndef WINDOWS_LINUX_STUB_H
#define WINDOWS_LINUX_STUB_H

#include <unistd.h>
#include <stdint.h>

#define getch getchar
#define stricmp(s1,s2) strcasecmp(s1,s2)

#define GetLastError() (errno)

#define O_BINARY 0 // https://stackoverflow.com/questions/2266992/no-o-binary-and-o-text-flags-in-linux

typedef int32_t INT;
typedef int32_t LONG;
typedef uint32_t DWORD;

typedef struct tagPOINT {
  LONG x;
  LONG y;
} POINT, *PPOINT, *NPPOINT, *LPPOINT;

typedef DWORD COLORREF;
typedef DWORD* LPCOLORREF;


#define MEMORY_ALLOCATION_ALIGNMENT 16

typedef struct __attribute__((aligned(16)))
{
  unsigned char data[16];
}
SLIST_HEADER;

typedef struct __attribute__((aligned(16)))
{
  unsigned char data[16];
}
SLIST_ENTRY;

#define InterlockedIncrement(ptr) __sync_add_and_fetch( ptr, 1 )
#define InterlockedDecrement(ptr) __sync_sub_and_fetch( ptr, 1 )
#define InterlockedExchange(ptr,val) __sync_lock_test_and_set(ptr, val )
#define InterlockedCompareExchange(ptr,newval,oldval) __sync_val_compare_and_swap( ptr, oldval, newval )
#define InterlockedCompareExchange16(ptr,newval,oldval) __sync_val_compare_and_swap( ptr, oldval, newval )
#define Sleep(n) ((n)>0 ? sleep(n) : sched_yield())

inline int itoa(int value, char *sp, int radix)
{
    char tmp[16];// be careful with the length of the buffer
    char *tp = tmp;
    int i;
    unsigned v;

    int sign = (radix == 10 && value < 0);    
    if (sign)
        v = -value;
    else
        v = (unsigned)value;

    while (v || tp == tmp)
    {
        i = v % radix;
        v /= radix;
        if (i < 10)
          *tp++ = i+'0';
        else
          *tp++ = i + 'a' - 10;
    }

    int len = tp - tmp;

    if (sign) 
    {
        *sp++ = '-';
        len++;
    }

    while (tp > tmp)
        *sp++ = *--tp;

    return len;
}
#endif
#endif

