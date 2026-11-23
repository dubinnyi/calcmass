#ifndef ISOTOPES_H_INCLUDED
#define ISOTOPES_H_INCLUDED

#define NI 14 /* Number of isotopes */

typedef struct
{
  char   name[3];
  int    imass;
  double dmass;
  int    valence;
  int    v_unambig;
} isotope_t;

extern isotope_t Isotopes[NI];

typedef int formula_t[NI+1];
typedef struct {int min, max;} int_limits_t;
typedef int_limits_t limits_t[NI+1];

typedef struct {
  formula_t f;
  double mass;
  double mass_diff;
  double mass_diff_ppm;
  double RDB;  /* Ring and Double Bond counter */
               /* Neutral formula should have integer and non-negative RDB value */
               /* half-integer and negative RDB means broken valence rules */
} fmass_t;

int Mass_Diff(const void *fm1, const void *fm2);

#define NM 11 /* Number of mass types */

typedef enum {
  MASS_M,
  MASS_POS_H,
  MASS_POS_Na,
  MASS_POS_Ca,
  MASS_NEG_H,
  MASS_POS,
  MASS_NEG,
  MASS_POS_K,
  MASS_POS_2H,
  MASS_POS_HNa,
  MASS_POS_2Na,
} mass_t;

typedef struct {
  char name[20];
  int elem1, count1;
  int elem2, count2;
  int elem3, count3;
  double madd;
  int    charge;
} mass_type_t;

extern mass_type_t Mass[NM];

extern int Mass_Init_Flag;

double Calc_Mass(formula_t fi, mass_t M);
int    Calc_RDB(formula_t fi, double *rdb);
double Calc_Charge(formula_t fi, mass_t M);
double Calc_Mass_Charge_Ratio(formula_t fi, mass_t M);
double Calc_Mass_M(formula_t fi);
void   mass_type_init(void);

void Parse_Formula(const char *fc, formula_t fi);
void Parse_Limits(const char *fc, limits_t lim);

void Loop_Copy(formula_t *dst, formula_t *src);
void Loop_Init(formula_t *ft, limits_t *lim);
void Loop_Next(formula_t *ft, limits_t *lim);
int  Loop_End(formula_t *ft, limits_t *lim);

#endif /* ISOTOPES_H_INCLUDED */
