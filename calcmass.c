#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "isotopes.h"

double Print_Formula_Mass(formula_t fi, mass_t mt);
void Print_Formula(formula_t fi, char *str, size_t str_size);
void Print_Limits(limits_t lim, char *str, size_t str_size);

#define N_Best 30
#define MAX_DIFF 0.1
#define MAX_DIFF_PPM 20.0
#define N_Chem 512

static fmass_t CHEM_LIST[N_Chem];
static fmass_t *chem = CHEM_LIST;

static void print_usage(int argc)
{
  printf(
    "Too few arguments found (%d).\nUsage:\n"
    "calcmass -mass [type] [Chemical formula]\n"
    "calcmass -chem [type] [Chemical limits] [mass]\n"
    "calcmass --chem-rdb [type] [Chemical limits] [mass]\n\n"
    "[type] is one of M, M+, M-, MH+, MNa+, MK+, MCa+, MH-,\n"
    "                 M2H++ (MH++, MH2+, MH+2), MHNa++, M2Na++\n"
    "[Chemical formula] is in the form C10,N6,O2,H13,Cl or C10H13N6O2Cl\n"
    "[Chemical limits]  is in the form C8-10,N6,O1-2,H1-13,Cl\n", argc);
}

int main(int argc, char **argv)
{
  formula_t fi, *pfi = &fi;
  limits_t lim, *plim = &lim;
  double Target_Mass;
  char fs[1024];
  mass_t MT;
  int rdb_flag = 0;

  mass_type_init();

  printf("\ncalcmass by Maxim A. Dubinnyi\n");
  printf("Version 2025-03-01\n");

  if (argc < 4) {
    print_usage(argc);
    return 1;
  }

  if (strcmp(argv[2], "M") == 0)           MT = MASS_M;
  else if (strcmp(argv[2], "MH+") == 0)    MT = MASS_POS_H;
  else if (strcmp(argv[2], "M+") == 0)     MT = MASS_POS;
  else if (strcmp(argv[2], "M-") == 0)     MT = MASS_NEG;
  else if (strcmp(argv[2], "MNa+") == 0)   MT = MASS_POS_Na;
  else if (strcmp(argv[2], "MK+") == 0)    MT = MASS_POS_K;
  else if (strcmp(argv[2], "MCa+") == 0)   MT = MASS_POS_Ca;
  else if (strcmp(argv[2], "MH-") == 0)    MT = MASS_NEG_H;
  else if (strcmp(argv[2], "M2H++") == 0)  MT = MASS_POS_2H;
  else if (strcmp(argv[2], "MH++") == 0)   MT = MASS_POS_2H;
  else if (strcmp(argv[2], "MH2+") == 0)   MT = MASS_POS_2H;
  else if (strcmp(argv[2], "MH+2") == 0)   MT = MASS_POS_2H;
  else if (strcmp(argv[2], "MHNa++") == 0) MT = MASS_POS_HNa;
  else if (strcmp(argv[2], "M2Na++") == 0) MT = MASS_POS_2Na;
  else {
    fprintf(stderr, "Unknown mass type: '%s'. Exiting.\n", argv[2]);
    return 1;
  }

  printf("Mass type is '%s' (%d), charge= %d, madd=%9.6f\n",
         Mass[MT].name, MT, Mass[MT].charge, Mass[MT].madd);

  if (strcmp(argv[1], "-mass") == 0 || strcmp(argv[1], "--mass") == 0) {
    Parse_Formula(argv[3], fi);
    Print_Formula(fi, fs, sizeof(fs));
    double rdb;
    Calc_RDB(fi, &rdb);
    printf("\nParsed formula: '%s', RDB = %4.1f\n", fs, rdb);
    (void)Print_Formula_Mass(fi, MT);
    return 0;
  } else if (strcmp(argv[1], "-chem") == 0 || strcmp(argv[1], "--chem") == 0 ||
             strcmp(argv[1], "-chem-rdb") == 0 || strcmp(argv[1], "--chem-rdb") == 0) {
    if (strcmp(argv[1], "-chem-rdb") == 0 || strcmp(argv[1], "--chem-rdb") == 0) {
      /* skip negative and half-integer RDB values */
      rdb_flag = 1;
    } else {
      rdb_flag = 0;
    }

    Parse_Limits(argv[3], lim);
    Print_Limits(lim, fs, sizeof(fs));
    printf("\nParsed limits: '%s'\n", fs);

    if (argc < 5) {
      fprintf(stderr, "Please enter the target mass\n");
      return 1;
    }

    char *endptr = NULL;
    errno = 0;
    Target_Mass = strtod(argv[4], &endptr);
    if (errno != 0 || endptr == argv[4] || *endptr != '\0' ||
        !isfinite(Target_Mass) || Target_Mass <= 0.0) {
      fprintf(stderr, "Invalid target mass: '%s'\n", argv[4]);
      return 1;
    }
    printf("Target mass is %10.6f\n", Target_Mass);
    printf("RDB check flag is %d\n", rdb_flag);

    int i = 0;
    long count = 0;
    chem = CHEM_LIST;
    for (Loop_Init(pfi, plim); Loop_End(pfi, plim); Loop_Next(pfi, plim)) {
      Loop_Copy(&(chem->f), pfi);
      chem->mass          = Calc_Mass_Charge_Ratio(chem->f, MT);
      chem->mass_diff     = chem->mass - Target_Mass;
      chem->mass_diff_ppm = (chem->mass_diff) / Target_Mass * 1000000.0;
      (void)Calc_RDB(chem->f, &(chem->RDB));
      double rdb = chem->RDB;
      if ((!rdb_flag) ||
          (rdb_flag && (rdb >= -0.01) && fabs(round(rdb) - rdb) <= 0.01)) {
        i++;
        chem++;
        count++;
      }
      if (i >= N_Chem) {
        qsort(CHEM_LIST, N_Chem, sizeof(fmass_t), Mass_Diff);
        i = N_Best;
        chem = CHEM_LIST + i;
      }
    }

    if (i > 1) {
      qsort(CHEM_LIST, (size_t)i, sizeof(fmass_t), Mass_Diff);
    }
    int n_chem = i;
    printf("Chemical formulas considered: %ld\n", count);

    int n_act = n_chem < N_Best ? n_chem : N_Best;
    if (n_act == 0) {
      printf("No matching chemical formulas found.\n");
      return 0;
    }

    int n_print;
    for (n_print = 0;
         n_print < n_act &&
         fabs(CHEM_LIST[n_print].mass_diff) <= MAX_DIFF &&
         fabs(CHEM_LIST[n_print].mass_diff_ppm) <= MAX_DIFF_PPM;
         n_print++) {
    }

    if (n_print == 0) n_print = 1;

    printf("Best %d formulas with mass difference less than %8.6f Da and %7.3f ppm",
           n_print, MAX_DIFF, MAX_DIFF_PPM);
    if (rdb_flag) {
      printf(" and integer non-negative RDB values:\n");
    } else {
      printf(":\n");
    }
    printf("%25s  %10s  %12s %11s  %4s\n", "Formula", "Theory", "Mass diff", "Mass diff", "RDB");
    for (int idx = 0; idx < n_print; idx++) {
      Print_Formula(CHEM_LIST[idx].f, fs, sizeof(fs));
      printf("%25s: %7.6f  %8.3f mDa %7.3f ppm  %4.1f\n", fs,
             CHEM_LIST[idx].mass,
             CHEM_LIST[idx].mass_diff * 1000.0,
             CHEM_LIST[idx].mass_diff_ppm,
             CHEM_LIST[idx].RDB);
    }
    if (!rdb_flag) {
      printf("use --chem-rdb instead of --chem\nto ignore formulas with NEGATIVE or HALF-INTEGER RDB values\n");
    }
    return 0;
  } else {
    fprintf(stderr, "Unknown argument: %s\n", argv[1]);
    return 1;
  }
}

double Print_Formula_Mass(formula_t fi, mass_t mt)
{
  int e;
  double mass = 0.0;
  double charge = Calc_Charge(fi, mt);

  printf("\nFormula:\n");
  for (e = 0; e < NI; e++) {
    double me = Isotopes[e].dmass * fi[e];
    mass += me;
    if (fi[e]) {
      printf("%2s %2d %9.6f %11.6f\n",
             Isotopes[e].name, fi[e], Isotopes[e].dmass, me);
    }
  }
  printf("------------------------------\n");
  printf("    %10s: (M)   %11.6f\n", Mass[mt].name, Calc_Mass(fi, mt));
  if (charge == 0.0) {
    printf("    %10s: (M/Z) N/A (neutral)\n", Mass[mt].name);
  } else {
    printf("    %10s: (M/Z) %11.6f\n", Mass[mt].name, Calc_Mass_Charge_Ratio(fi, mt));
  }
  return mass;
}

void Print_Formula(formula_t fi, char *str, size_t str_size)
{
  int e, n, flag = 0;
  char *s = str;
  size_t rem = str_size;

  if (rem == 0) return;
  *s = '\0';

  for (e = 1; e < NI; e++) {
    if (fi[e] >= 1 && flag) {
      n = snprintf(s, rem, ",");
      if (n < 0 || (size_t)n >= rem) return;
      s += n;
      rem -= (size_t)n;
    }
    if (fi[e] >= 1) {
      n = snprintf(s, rem, "%s", Isotopes[e].name);
      if (n < 0 || (size_t)n >= rem) return;
      s += n;
      rem -= (size_t)n;
      flag = 1;
    }
    if (fi[e] >= 2) {
      n = snprintf(s, rem, "%d", fi[e]);
      if (n < 0 || (size_t)n >= rem) return;
      s += n;
      rem -= (size_t)n;
    }
  }
  if (fi[0] > 0) {
    n = snprintf(s, rem, "-");
    if (n < 0 || (size_t)n >= rem) return;
    s += n;
    rem -= (size_t)n;
  } else if (fi[0] < 0) {
    n = snprintf(s, rem, "+");
    if (n < 0 || (size_t)n >= rem) return;
    s += n;
    rem -= (size_t)n;
  }
  if (fi[0] <= -2 || fi[0] >= 2) {
    n = snprintf(s, rem, "%d", fi[0] > 0 ? fi[0] : -fi[0]);
    if (n < 0 || (size_t)n >= rem) return;
  }
}

void Print_Limits(limits_t lim, char *str, size_t str_size)
{
  int e, n, flag = 0;
  char *s = str;
  size_t rem = str_size;

  if (rem == 0) return;
  *s = '\0';

  for (e = 1; e < NI; e++) {
    if (lim[e].max && flag) {
      n = snprintf(s, rem, ",");
      if (n < 0 || (size_t)n >= rem) return;
      s += n;
      rem -= (size_t)n;
    }
    if (lim[e].max) {
      n = snprintf(s, rem, "%s%d-%d", Isotopes[e].name, lim[e].min, lim[e].max);
      if (n < 0 || (size_t)n >= rem) return;
      s += n;
      rem -= (size_t)n;
      flag = 1;
    }
  }
}
