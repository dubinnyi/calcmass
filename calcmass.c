#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "isotopes.h"

#define DEBUG 0
#define IFDEBUG if(DEBUG)

void get_formula(char * fc, formula_t fi);
double Print_Formula_Mass(formula_t fi, mass_t mt);
void Print_Formula(formula_t fi, char* str);
void Print_Limits(limits_t lim, char* str);

#define N_Best 30
#define MAX_DIFF 1.0
#define N_Chem 512
fmass_t CHEM_LIST[N_Chem];
fmass_t *chem = CHEM_LIST;

int main(int argc, char** argv)
{
  int i;
  formula_t fi,  *pfi=&fi;
  limits_t  lim, *plim=&lim;
  double mass, Target_Mass, *TM=&Target_Mass;
  char fs[1024];
  mass_t MT;

  printf("\ncalcmass by Maxim A. Dubinnyi\n", NI);
  printf("\nDatabase of %d known isotopes loaded:\n", NI);
  for(i=0;i<NI;i++)
    printf("%s(%d)%s",Isotopes[i].name, Isotopes[i].imass, i==NI-1?"\n":", ");

  if(argc<4){
    printf(
      "Too few arguments found (%d).\nUsage:\n"
      "calcmass -mass [type] [Chemical formula] \n"
      "calcmass -chem [type] [Chemical limits]  [mass]\n\n"
      "[type] is one of M, M+, M-, MH+, MNa+, MCa+, MH-\n"
      "[Chemical formula] is in the form C10,N6,O2,H13,Cl\n"
      "[Chemical limits]  is in the form C8-10,N6,O1-2,H1-13,Cl\n",argc)
      ;
    exit(0);
  }

  if(strcmp(argv[2],"M")==0)        MT=MASS_M;
  else if(strcmp(argv[2],"MH+")==0) MT=MASS_POS_H;
  else if(strcmp(argv[2],"M+")==0)  MT=MASS_POS;
  else if(strcmp(argv[2],"M-")==0)  MT=MASS_NEG;
  else if(strcmp(argv[2],"MNa+")==0)MT=MASS_POS_Na;
  else if(strcmp(argv[2],"MK+")==0) MT=MASS_POS_K;
  else if(strcmp(argv[2],"MCa+")==0)MT=MASS_POS_Ca;
  else if(strcmp(argv[2],"MH-")==0) MT=MASS_NEG_H;
  else{
    printf("Unknown mass type: \'%s\'. Exiting.\n",argv[2]);
    return(-1);
  }
  
  if(strcmp(argv[1],"-mass")==0){
    Parse_Formula(argv[3],fi);
    Print_Formula(fi,fs);
    printf("\nParsed formula: \'%s\'\n",fs);
    mass=Print_Formula_Mass(fi,MT);
    return(0);
  }else if(strcmp(argv[1],"-chem")==0){
    Parse_Limits(argv[3],lim);
    Print_Limits(lim,fs);
    printf("\nParsed limits: \'%s\'\n",fs);
    if(argc<5){
      printf("Please enter the target mass\n");
      return(-1);
    }else{
      float fm, *pfm=&fm;
      /*printf("argv[4] is \'%s\'\n",argv[4]);*/
      i=sscanf(argv[4],"%f", pfm);
      Target_Mass=fm;
      printf("Target mass is %10.6f\n",Target_Mass);
    }
    int i=0,count=0;
    chem = CHEM_LIST;
    for(
	Loop_Init(pfi,plim); 
	Loop_End (pfi,plim);
	Loop_Next(pfi,plim)){
      Loop_Copy(&(chem->f),pfi);
      chem->mass=Calc_Mass(chem->f,MT);
      chem->mass_diff=chem->mass-Target_Mass;
      i++;chem++;count++;
      if(i>=N_Chem){
	qsort(CHEM_LIST,N_Chem,sizeof(fmass_t),Mass_Diff);
	i=N_Best;
	chem=CHEM_LIST+i;
      }
      /*
      Print_Formula(*pfi,fs);
      printf("%15s: %10.6f  %10.6f\n",fs, mass, mass-Target_Mass);
      */
    }
    qsort(CHEM_LIST,i,sizeof(fmass_t),Mass_Diff);
    int n_chem=i;
    printf("Chemical formulas considered: %d\n",count);

    int n_act=n_chem<N_Best?n_chem:N_Best;
    
    int n_print;
    for(n_print=0;n_print<n_act;n_print++){
      if(fabs(CHEM_LIST[n_print].mass_diff)>=MAX_DIFF)break;
    }
    
    printf("Best %d formulas with mass differense less then %8.6f:\n",n_print, MAX_DIFF);
    for(i=0;i<n_print;i++){
      Print_Formula(CHEM_LIST[i].f,fs);
      printf("%25s: %10.6f  %10.6f\n",fs, CHEM_LIST[i].mass, CHEM_LIST[i].mass_diff);
    }
    return(0);
  }else{
    printf("Unknown argument: %s", argv[1]);
    return(-1);
  }

};

double Print_Formula_Mass(formula_t fi, mass_t mt)
{
  int m,e;
  printf("\nFormula:\n");
  double mass=0.0, me;
  for(e=0;e<NI;e++){
    me= Isotopes[e].dmass*fi[e];
    mass+=me;
    if(fi[e])
      printf("%2s %2d %9.6f %11.6f\n", 
	Isotopes[e].name, fi[e], Isotopes[e].dmass, me);
  }
  printf("------------------------------\n");
/*  for(m=0;m<NM;m++){*/
  printf("    %10s: %11.6f\n",Mass[mt].name,Calc_Mass(fi,mt));
/*  }*/
  return(mass);
};

void Print_Formula(formula_t fi, char* str)
{
   int e,n,flag=0;
   char *s=str;
   for(e=1;e<NI;e++){
     if(fi[e]>=1 && flag){
       n=sprintf(s,"%s",",");
       if(n>0){s+=n;} else {*s=0;break;};
     }
     if(fi[e]>=1){
       n=sprintf(s,"%s",Isotopes[e].name);
       if(n>0){s+=n;} else {*s=0;break;};
       flag=1;
     }
     if(fi[e]>=2){
       n=sprintf(s,"%d",fi[e]);
       if(n>0){s+=n;} else {*s=0;break;};
     }
   }
   if(fi[0]>0){
     n=sprintf(s,"%s","-");
     if(n>0){s+=n;} else {*s=0;return;};
   }
   if(fi[0]<0){
     n=sprintf(s,"%s","+");
     if(n>0){s+=n;} else {*s=0;return;};
   }
   if(fi[0]<=-2 || fi[0]>=2){
     n=sprintf(s,"%d",fi[0]*(fi[0]>0?1:-1));
     if(n>0){s+=n;} else {*s=0;return;};
   }
};

void Print_Limits(limits_t lim, char* str)
{
   int e,n,flag=0;
   char *s=str;
   for(e=1;e<NI;e++){
     if(lim[e].max && flag){
       n=sprintf(s,"%s",",");
       if(n>0){s+=n;} else {*s=0;break;};
     }
     if(lim[e].max){
       n=sprintf(s,"%s",Isotopes[e].name);
       if(n>0){s+=n;} else {*s=0;break;};
       flag=1;
       n=sprintf(s,"%d-",lim[e].min);
       if(n>0){s+=n;} else {*s=0;break;};
       n=sprintf(s,"%d",lim[e].max);
       if(n>0){s+=n;} else {*s=0;break;};
     }
   }
};
