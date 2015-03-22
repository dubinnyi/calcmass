#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "isotopes.h"

#define DEBUG 0
#define IFDEBUG if(DEBUG)

void Loop_Copy(formula_t *dst, formula_t *src)
{
  int e;
  for(e=0;e<NI+1;e++)
    (*dst)[e]=(*src)[e];
};


void Loop_Init(formula_t *ft, limits_t *lim)
{
  int e;
  (*ft)[0]=0;
  IFDEBUG printf("Loop_Init: ");
  for(e=1;e<NI;e++){
    (*ft)[e]=(*lim)[e].min;
    IFDEBUG printf("%2d ",(*ft)[e]);
  }
  IFDEBUG printf("\n");
  (*ft)[NI]=0;
};

void Loop_Next(formula_t *ft, limits_t *lim)
{
  int e;
  (*ft)[1]++;
  for(e=1; e<NI&&(*ft)[e]>(*lim)[e].max; e++){
    (*ft)[e]=(*lim)[e].min;
    (*ft)[e+1]++;
  }
  IFDEBUG {
    printf("Loop_Next: ");
    for(e=1;e<NI;e++)printf("%2d ",(*ft)[e]);
    printf("\n");
  }
};


int Loop_End(formula_t *ft, limits_t *lim)
{
  int e,flag=(*ft)[NI]==0;
//  for(e=1; e<NI&&flag; e++)
//    flag=flag && (*ft)[e]>=(*lim)[e].max;
  IFDEBUG { 
    printf("Loop_End : ");
    for(e=1;e<NI;e++)printf("%2d ",(*ft)[e]);
    printf(" flag=%d\n",flag);
  }
  return flag;
};
