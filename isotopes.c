#include <stdio.h>
#include <math.h>
#include "isotopes.h"

static int Mass_Init_Flag;

isotope_t Isotopes[NI] =
{
  { "e",  0,  0.000549 },  /*  0  ELECTRON    */
  { "H",  1,  1.007825 },  /*  1  PROTON      */
  { "C", 12, 12.000000 },  /*  2  CARBON      */
  { "N", 14, 14.003074 },  /*  3  NITROGEN    */
  { "O", 16, 15.994915 },  /*  4  OXIGEN      */
  { "P", 31, 30.973763 },  /*  5  PHOSPHORUS  */
  { "S", 32, 31.972072 },  /*  6  SERUM       */ 
  {"Na", 23, 22.989770 },  /*  7  SODIUM      */ 
  {"Ca", 40, 39.962591 },  /*  8  CALCIUM     */ 
  { "K", 39, 38.963708 },  /*  9  POTASSIUM   */ 
  
  { "B", 11, 11.009305 },  /* 10  BORON       */ 
  { "F", 19, 18.998403 },  /* 11  FLUORINE    */ 
  {"Cl", 35, 34.968853 },  /* 12  CHLORINE    */ 
  {"Br", 79, 78.918336 },  /* 13  BROMINE     */ 
};

mass_type_t Mass[NM]={
/* name, elem1,count1,elem2,count2,elem3,count3,madd,charge */
 {"[M] ",    0,     0,    0,     0,    0,     0, 0.0, +0},  /* MASS_M:      No changes      */
 {"[M+H]+",  1,     1,    0,    -1,    0,     0, 0.0, +1},  /* MASS_POS_H:    1H(+1), e(-1) */
 {"[M+Na]+", 7,     1,    0,    -1,    0,     0, 0.0, +1},  /* MASS_POS_Na: 23Na(+1), e(-1) */
 {"[M+Ca-H]+",8,    1,    0,    -1,    1,    -1, 0.0, +1},  /* MASS_POS_Ca: 40Ca(+1), e(-1) H(-1)*/
 {"[M-H]-",  1,    -1,    0,    +1,    0,     0, 0.0, -1},  /* MASS_NEG_H:    1H(-1), e(+1) */
 {"[M]+",    0,    -1,    0,     0,    0,     0, 0.0, +1},  /* MASS_POS:      e(-1) */
 {"[M]-",    0,    +1,    0,     0,    0,     0, 0.0, -1},  /* MASS_NEG:      e(+1) */
 {"[M+K]+",  9,     1,    0,    -1,    0,     0, 0.0, +1},  /* MASS_POS_K:  39K(+1),  e(-1) */
 {"[M+2H]++",1,     2,    0,    -2,    0,     0, 0.0, +2},  /* MASS_POS_2H:  1H(+2), e(-2) */
 {"[M+H+Na]++",1,   1,    0,    -2,    7,     1, 0.0, +2},  /* MASS_POS_2H:  1H(+2), e(-2) */
 {"[M+2Na]++",7,    2,    0,    -2,    0,     0, 0.0, +2},  /* MASS_POS_2H:  1H(+2), e(-2) */
};

void mass_type_init()
{
  int m;
  for(m=0;m<NM;m++){
    Mass[m].madd=
      Isotopes[Mass[m].elem1].dmass*Mass[m].count1+
      Isotopes[Mass[m].elem2].dmass*Mass[m].count2+
      Isotopes[Mass[m].elem3].dmass*Mass[m].count3;
    Mass[m].charge=(Mass[m].elem1==0?(-Mass[m].count1):0)+
                   (Mass[m].elem2==0?(-Mass[m].count2):0)+
                   (Mass[m].elem3==0?(-Mass[m].count3):0);
  };
  Mass_Init_Flag=1;
};

double Calc_Mass_M(formula_t fi)
{
  if(!Mass_Init_Flag)mass_type_init();
  return(Calc_Mass(fi, MASS_M));
};

double Calc_Mass(formula_t fi, mass_t M)
{
  int e;
  double mass=0.0, me;

  if(!Mass_Init_Flag)mass_type_init();
  for(e=0;e<NI;e++)mass+=Isotopes[e].dmass*fi[e];
  mass+=Mass[M].madd;
  return(mass);
};

int Calc_Charge(formula_t fi, mass_t M)
{
  int e;
  double me;
  double charge=0.0;

  if(!Mass_Init_Flag)mass_type_init();
  charge=-fi[0];
  charge+=Mass[M].charge;
  return(charge);
};

double Calc_Mass_Charge_Ratio(formula_t fi, mass_t M)
{
  double mass;
  int    charge;

  mass =   Calc_Mass(fi,M);
  charge = Calc_Charge(fi, M);
  return(mass/charge);
}

int Mass_Diff(const void *fm1, const void *fm2)
{
  return( fabs(((fmass_t*)fm1)->mass_diff)>fabs(((fmass_t*)fm2)->mass_diff) );
};
