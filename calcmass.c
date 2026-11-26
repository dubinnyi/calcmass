/**
 * @file calcmass.c
 * @brief CLI entry point and output formatting for the calcmass mass-spectrometry tool.
 *
 * Provides two primary command-line modes:
 * 1. `-mass` / `--mass [-s|--s] <adduct> <formula>`
 *    Calculates monoisotopic mass, m/z ratio, elemental breakdown, and RDB for a formula.
 * 2. `-chem` / `--chem` / `-chem-rdb` / `--chem-rdb [-s|--s] <adduct> <limits> <target_mass>`
 *    Enumerates candidate formulas within the given element bounds, retains the closest
 *    matches within `MAX_DIFF` (Da) and `MAX_DIFF_PPM` (ppm), and prints up to `N_Best`
 *    results sorted by absolute mass error.
 *
 * Passing `-s` or `--s` anywhere on the command line enables short output mode
 * (only the M/Z line for `-mass`, or the results table starting from `Formula` for `-chem`).
 */

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

/** @brief Maximum number of top-ranked candidate formulas retained and displayed. */
#define N_Best 30
/** @brief Maximum absolute mass difference (Da) for printing search hits. */
#define MAX_DIFF 0.1
/** @brief Maximum relative mass difference (ppm) for printing search hits. */
#define MAX_DIFF_PPM 20.0
/** @brief Capacity of the working candidate buffer flushed via qsort down to N_Best. */
#define N_Chem 512

static fmass_t CHEM_LIST[N_Chem];
static fmass_t *chem = CHEM_LIST;

/**
 * @brief Print command-line usage instructions to stdout.
 *
 * @param[in] argc Actual argument count received by `main()`, reported in the diagnostic line.
 */
static void print_usage(int argc)
{
  printf(
    "Too few arguments found (%d).\nUsage:\n"
    "calcmass [-s|--s] -mass [type] [Chemical formula]\n"
    "calcmass [-s|--s] -chem [type] [Chemical limits] [mass]\n"
    "calcmass [-s|--s] --chem-rdb [type] [Chemical limits] [mass]\n\n"
    "Options:\n"
    "  -s, --s          Short output (M/Z line only for -mass; table only for -chem)\n\n"
    "[type] is one of M, M+, M-, MH+, MNa+, MK+, MCa+, MH-,\n"
    "                 M2H++ (MH++, MH2+, MH+2), MHNa++, M2Na++\n"
    "[Chemical formula] is in the form C10,N6,O2,H13,Cl, C10H13N6O2Cl, or Ca(OH)2\n"
    "[Chemical limits]  is in the form C8-10,N6,O1-2,H1-13,Cl\n", argc);
}

/**
 * @brief Program entry point: parses CLI arguments and dispatches `-mass` or `-chem` modes.
 *
 * Recognizes optional `-s` / `--s` (short output) at any position in `argv`, leaving
 * positional arguments:
 * - `pos[0]`: Mode flag: `-mass`/`--mass`, `-chem`/`--chem`, or `-chem-rdb`/`--chem-rdb`.
 * - `pos[1]`: Ion/adduct name (e.g., `M`, `MH+`, `MNa+`, `M2H++`).
 * - `pos[2]`: Chemical formula (in `-mass` mode) or element limits (in `-chem` modes).
 * - `pos[3]`: Required in `-chem` modes; strictly positive finite target mass/m/z in Da.
 *
 * @param[in] argc Number of command-line arguments (at least 3 positional args for `-mass`, 4 for `-chem`).
 * @param[in] argv Array of NUL-terminated command-line argument strings.
 * @return `0` on success; `1` on invalid arguments, unknown mode/adduct, parse error, or invalid target mass.
 */
int main(int argc, char **argv)
{
  formula_t fi, *pfi = &fi;
  limits_t lim, *plim = &lim;
  double Target_Mass;
  char fs[1024];
  mass_t MT;
  int rdb_flag = 0;
  int short_flag = 0;
  int npos = 0;
  int rc = 0;

  mass_type_init();

  char **pos = NULL;
  if (argc > 1) {
    pos = (char **)malloc((size_t)(argc - 1) * sizeof(*pos));
    if (pos == NULL) {
      fprintf(stderr, "Memory allocation failed\n");
      return 1;
    }
  }

  for (int a = 1; a < argc; a++) {
    if (strcmp(argv[a], "-s") == 0 || strcmp(argv[a], "--s") == 0) {
      short_flag = 1;
    } else {
      pos[npos++] = argv[a];
    }
  }

  if (!short_flag) {
    printf("\ncalcmass by Maxim A. Dubinnyi\n");
    printf("Version 2026-11-26\n");
  }

  if (npos < 3) {
    print_usage(argc);
    free(pos);
    return 1;
  }

  if (strcmp(pos[1], "M") == 0)           MT = MASS_M;
  else if (strcmp(pos[1], "MH+") == 0)    MT = MASS_POS_H;
  else if (strcmp(pos[1], "M+") == 0)     MT = MASS_POS;
  else if (strcmp(pos[1], "M-") == 0)     MT = MASS_NEG;
  else if (strcmp(pos[1], "MNa+") == 0)   MT = MASS_POS_Na;
  else if (strcmp(pos[1], "MK+") == 0)    MT = MASS_POS_K;
  else if (strcmp(pos[1], "MCa+") == 0)   MT = MASS_POS_Ca;
  else if (strcmp(pos[1], "MH-") == 0)    MT = MASS_NEG_H;
  else if (strcmp(pos[1], "M2H++") == 0)  MT = MASS_POS_2H;
  else if (strcmp(pos[1], "MH++") == 0)   MT = MASS_POS_2H;
  else if (strcmp(pos[1], "MH2+") == 0)   MT = MASS_POS_2H;
  else if (strcmp(pos[1], "MH+2") == 0)   MT = MASS_POS_2H;
  else if (strcmp(pos[1], "MHNa++") == 0) MT = MASS_POS_HNa;
  else if (strcmp(pos[1], "M2Na++") == 0) MT = MASS_POS_2Na;
  else {
    fprintf(stderr, "Unknown mass type: '%s'. Exiting.\n", pos[1]);
    free(pos);
    return 1;
  }

  if (!short_flag) {
    printf("Mass type is '%s' (%d), charge= %d, madd=%9.6f\n",
           Mass[MT].name, MT, Mass[MT].charge, Mass[MT].madd);
  }

  if (strcmp(pos[0], "-mass") == 0 || strcmp(pos[0], "--mass") == 0) {
    if (npos > 3) {
      fprintf(stderr, "Too many arguments for -mass mode (unexpected: '%s')\n", pos[3]);
      free(pos);
      return 1;
    }
    if (Parse_Formula(pos[2], fi) != 0) {
      free(pos);
      return 1;
    }
    if (short_flag) {
      double charge = Calc_Charge(fi, MT);
      if (charge == 0.0) {
        printf("%s: (M)   %11.6f\n", Mass[MT].name, Calc_Mass(fi, MT));
      } else {
        printf("%s: (M/Z) %11.6f\n", Mass[MT].name, Calc_Mass_Charge_Ratio(fi, MT));
      }
      free(pos);
      return 0;
    }
    Print_Formula(fi, fs, sizeof(fs));
    double rdb;
    Calc_RDB(fi, &rdb);
    printf("\nParsed formula: '%s', RDB = %4.1f\n", fs, rdb);
    (void)Print_Formula_Mass(fi, MT);
    free(pos);
    return 0;
  } else if (strcmp(pos[0], "-chem") == 0 || strcmp(pos[0], "--chem") == 0 ||
             strcmp(pos[0], "-chem-rdb") == 0 || strcmp(pos[0], "--chem-rdb") == 0) {
    if (npos > 4) {
      fprintf(stderr, "Too many arguments for %s mode (unexpected: '%s')\n", pos[0], pos[4]);
      free(pos);
      return 1;
    }
    if (strcmp(pos[0], "-chem-rdb") == 0 || strcmp(pos[0], "--chem-rdb") == 0) {
      /* skip negative and half-integer RDB values */
      rdb_flag = 1;
    } else {
      rdb_flag = 0;
    }

    if (Parse_Limits(pos[2], lim) != 0) {
      free(pos);
      return 1;
    }
    if (!short_flag) {
      Print_Limits(lim, fs, sizeof(fs));
      printf("\nParsed limits: '%s'\n", fs);
    }

    if (npos < 4) {
      fprintf(stderr, "Please enter the target mass\n");
      free(pos);
      return 1;
    }

    char *endptr = NULL;
    errno = 0;
    Target_Mass = strtod(pos[3], &endptr);
    if (errno != 0 || endptr == pos[3] || *endptr != '\0' ||
        !isfinite(Target_Mass) || Target_Mass <= 0.0) {
      fprintf(stderr, "Invalid target mass: '%s'\n", pos[3]);
      free(pos);
      return 1;
    }
    if (!short_flag) {
      printf("Target mass is %10.6f\n", Target_Mass);
      printf("RDB check flag is %d\n", rdb_flag);
    }

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
    if (!short_flag) {
      printf("Chemical formulas considered: %ld\n", count);
    }

    int n_act = n_chem < N_Best ? n_chem : N_Best;
    if (n_act == 0) {
      printf("No matching chemical formulas found.\n");
      free(pos);
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

    if (!short_flag) {
      printf("Best %d formulas with mass difference less than %8.6f Da and %7.3f ppm",
             n_print, MAX_DIFF, MAX_DIFF_PPM);
      if (rdb_flag) {
        printf(" and integer non-negative RDB values:\n");
      } else {
        printf(":\n");
      }
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
    if (!short_flag && !rdb_flag) {
      printf("use --chem-rdb instead of --chem\nto ignore formulas with NEGATIVE or HALF-INTEGER RDB values\n");
    }
    rc = 0;
  } else {
    fprintf(stderr, "Unknown argument: %s\n", pos[0]);
    rc = 1;
  }

  free(pos);
  return rc;
}

/**
 * @brief Print an elemental mass breakdown, adduct mass `(M)`, and `(M/Z)` to stdout.
 *
 * Prints each non-zero isotope contribution in `fi`, followed by the adduct-adjusted
 * total mass `(M)` and mass-to-charge ratio `(M/Z)` (or `"N/A (neutral)"` when net
 * charge is 0).
 *
 * @param[in] fi Composition vector with initialized entries in `[0, NI)`.
 * @param[in] mt Adduct type selector; must be a valid `mass_t` value in `[0, NM)`.
 * @return Sum of elemental masses in `fi` (in Da) before adding `Mass[mt].madd`.
 */
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

/**
 * @brief Format a `formula_t` composition vector as a comma-separated string.
 *
 * Emits elements in `Isotopes[]` order (`1 .. NI-1`), omitting count `1` (e.g.,
 * `"C10,H13,N6,O2,Cl"`), and appends any explicit charge encoded in `fi[0]`
 * (`-` for `fi[0] > 0`, `+` for `fi[0] < 0`, with magnitude suffix when `|fi[0]| >= 2`).
 * Output is safely truncated to fit within `str_size` and always NUL-terminated
 * when `str_size > 0`.
 *
 * @param[in]  fi       Composition vector with initialized entries in `[0, NI)`.
 * @param[out] str      Destination character buffer; must be non-NULL when `str_size > 0`.
 * @param[in]  str_size Total capacity of `str` in bytes (including NUL terminator).
 */
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

/**
 * @brief Format a `limits_t` bounds array as a comma-separated range string.
 *
 * Emits `"<Symbol><min>-<max>"` for each element `e` in `[1, NI)` whose `max != 0`
 * (e.g., `"C8-10,H1-13,N6-6,O1-2,Cl1-1"`). Output is safely truncated to `str_size`
 * and always NUL-terminated when `str_size > 0`.
 *
 * @param[in]  lim      Element bounds array with initialized entries in `[1, NI)`.
 * @param[out] str      Destination character buffer; must be non-NULL when `str_size > 0`.
 * @param[in]  str_size Total capacity of `str` in bytes (including NUL terminator).
 */
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
