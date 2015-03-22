#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "isotopes.h"

#define DEBUG 0
#define IFDEBUG if(DEBUG)

int issign(char c)
{
  if(c=='+' || c=='-') return(1);
  return(0);
};

void Parse_Formula(char * fc, formula_t fi)
{
  char elem[3],digit[30],c,d,*fd;
  int i=0,e=0,ne=0,ni,id,it,E=-1;
  int parse_sign=0;

  for(e=0;e<NI;e++)fi[e]=0;

  IFDEBUG printf("\nParsing formula string: \'%s\'\n\n", fc);
  while(*fc)
  {
    /* Skip some characters */
    it=0;
    while(*fc && !isalpha(*fc) && !isdigit(*fc) && !issign(*fc))
       {fc++;it++;};

    IFDEBUG printf("%d trailing characters skipped\n",it);

    /* get element name */
    i=0;
   
    parse_sign=0;
    while(*fc)
    {
      c=*fc;
      if(!isalpha(c) && !issign(c) )break;
      IFDEBUG printf("Parsing character \'%c\'\n",c);
      fc++;
      if(issign(c)){
	strcpy(elem,"e");
	i=1;
	if(c=='+')parse_sign=-1; /* Electrons are SUBTRACTED */
	if(c=='-')parse_sign=+1; /* Electrons are ADDED      */
	IFDEBUG printf("Molecule charge is parsed\n",c);
	IFDEBUG printf("Electrons are added with sign %d\n",parse_sign);
	break;
      }
      if(i==0){
	elem[i]=toupper(c);
	i++;
      }else if(i==1){
	elem[i]=tolower(c);
	i++;
      }else if(i==2){
	elem[i]=0;
	i++;
	break;
      }
    };
    elem[i]=0;
    
    IFDEBUG printf("Element name is \'%s\'\n",elem);
    
    E=-1;
    for(e=0;e<NI;e++)
      if(strcmp(elem,Isotopes[e].name)==0)
	E=e;

    if(E==-1)
      printf("Unknown elemet: \'%s\'\n",elem);
    else
      IFDEBUG printf("Element ID is:  %d\n",E);
    
    /* get integer - number of elements */
    fd=digit;
    while(*fc){
      d=*fc;
      if(!isdigit(d))break;
      fc++;
      IFDEBUG printf("Parsing character \'%c\'\n",d);
      *fd++=d;
    }
    *fd=0;
    if(strlen(digit)==0)id=1;
    else id=atoi(digit);
    if(parse_sign)id*=parse_sign;

    IFDEBUG printf("Element counting is %d\n", id);

    if(E!=-1)fi[E]=id;
    IFDEBUG printf("\n");
  }
  
};



void Parse_Limits(char * fc, limits_t lim)
{

  char elem[3],digit[30],c,d,*fd;
  int i=0,e=0,ne=0,ni,id,it,E=-1,tmpi;

  for(e=0;e<NI;e++){
    lim[e].min=0;
    lim[e].max=0;
  }

  IFDEBUG printf("\nParsing limits string: \'%s\'\n\n", fc);
  while(*fc)
  {
    /* Skip some characters */
    it=0;
    while(*fc && !isalpha(*fc) && !isdigit(*fc))
       {fc++;it++;};

    IFDEBUG printf("%d trailing characters skipped\n",it);

    /* get element name */
    i=0;
   
    while(*fc)
    {
      c=*fc;
      if( !isalpha(c) )break;
      IFDEBUG printf("Parsing character \'%c\'\n",c);
      fc++;
      if(i==0){
	elem[i]=toupper(c);
	i++;
      }else if(i==1){
	elem[i]=tolower(c);
	i++;
      }else if(i==2){
	elem[i]=0;
	i++;
	break;
      }
    };
    elem[i]=0;
    
    IFDEBUG printf("Element name is \'%s\'\n",elem);
    
    E=-1;
    for(e=0;e<NI;e++)
      if(strcmp(elem,Isotopes[e].name)==0)
	E=e;

    if(E==-1)
      printf("Unknown elemet: \'%s\'\n",elem);
    else
      IFDEBUG printf("Element ID is:  %d\n",E);
    
    /* get integer - lower of element count */
    fd=digit;
    while(*fc){
      d=*fc;
      if(!isdigit(d))break;
      fc++;
      IFDEBUG printf("Parsing character \'%c\'\n",d);
      *fd++=d;
    }
    *fd=0;
    if(strlen(digit)==0)id=1;
    else id=atoi(digit);

    IFDEBUG printf("Element lower limit is %d\n", id);

    if(E!=-1)lim[E].min=id;
    IFDEBUG printf("\n");
    
    if(*fc=='-'){
      fc++;
      fd=digit;
      while(*fc){
	d=*fc;
	if(!isdigit(d))break;
	fc++;
	IFDEBUG printf("Parsing character \'%c\'\n",d);
	*fd++=d;
      }
      *fd=0;
      if(strlen(digit)==0)id=1;
      else id=atoi(digit);
    }
    IFDEBUG printf("Element upper limit is %d\n", id);
    if(E!=-1)lim[E].max=id;
    if(lim[E].max<lim[E].min){
      tmpi=lim[E].max;
      lim[E].max=lim[E].min;
      lim[E].min=tmpi;
      IFDEBUG printf("Upper and lower limits are changed\n", id);
    }
  }
  
};


