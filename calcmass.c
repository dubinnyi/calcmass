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

/** @brief RDB filter configuration applied during candidate enumeration. */
typedef struct {
  int    require_int; /**< Non-zero for `--chem-rdb`: require integer neutral-molecule RDB. */
  int    has_min;     /**< Non-zero if `min_val` lower bound is enforced. */
  double min_val;     /**< Inclusive lower RDB bound (default `0.0` for `--chem-rdb`). */
  int    has_max;     /**< Non-zero if `max_val` upper bound is enforced. */
  double max_val;     /**< Inclusive upper RDB bound. */
  int    strict_min;  /**< Non-zero (`--rdb-strict`): evaluate bounds against `RDB_min` (S=2, P=3) only. */
} rdb_filter_t;

/**
 * @brief Parse a finite floating-point number from a CLI argument string.
 *
 * @param[in]  s   NUL-terminated string to parse.
 * @param[out] out Destination for the parsed value.
 * @return `0` on success; `-1` if malformed, out of range, or non-finite.
 */
static int parse_double_arg(const char *s, double *out)
{
  char *endptr = NULL;
  errno = 0;
  double val = strtod(s, &endptr);
  if (errno != 0 || endptr == s || *endptr != '\0' || !isfinite(val)) {
    return -1;
  }
  *out = val;
  return 0;
}

/**
 * @brief Format an RDB range `[rdb_min, rdb_max]` into a human-readable string.
 *
 * Formats unambiguous values (`rdb_min == rdb_max`) as `"%4.1f"` (e.g., `"10.0"`,
 * `" 8.0"`) and variable-valence ranges as `"%.1f..%.1f"` (e.g., `"1.0..3.0"`).
 *
 * @param[in]  rdb_min  Lower RDB bound (minimum valences, e.g., S=2, P=3).
 * @param[in]  rdb_max  Upper RDB bound (maximum valences, e.g., S=6, P=5).
 * @param[out] buf      Destination character buffer.
 * @param[in]  buf_size Capacity of `buf` in bytes.
 */
static void format_rdb(double rdb_min, double rdb_max, char *buf, size_t buf_size)
{
  if (fabs(rdb_max - rdb_min) < 1e-6) {
    snprintf(buf, buf_size, "%4.1f", rdb_min);
  } else {
    snprintf(buf, buf_size, "%.1f..%.1f", rdb_min, rdb_max);
  }
}

/**
 * @brief Test whether a formula's neutral-molecule RDB interval `[rdb_min, rdb_max]` satisfies `flt`.
 *
 * When `flt->require_int` is set (`--chem-rdb`), rejects half-integer RDB values
 * (parity is identical across `[rdb_min, rdb_max]` because sulfur/phosphorus valence
 * steps are always even). In normal mode, range bounds pass if `[rdb_min, rdb_max]`
 * overlaps `[min_val, max_val]`. With `strict_min` (`--rdb-strict`), both bounds
 * are tested strictly against `rdb_min` (organic baseline S=2, P=3).
 *
 * @param[in] rdb_min Lower RDB bound of the neutral formula.
 * @param[in] rdb_max Upper RDB bound of the neutral formula.
 * @param[in] flt     Active RDB filter configuration.
 * @return `1` if the candidate passes the filter; `0` otherwise.
 */
static int rdb_passes(double rdb_min, double rdb_max, const rdb_filter_t *flt)
{
  if (flt->require_int) {
    if (fabs(round(rdb_min) - rdb_min) > 0.01) return 0;
  }

  if (flt->strict_min) {
    if (flt->has_min && rdb_min < flt->min_val - 1e-6) return 0;
    if (flt->has_max && rdb_min > flt->max_val + 1e-6) return 0;
  } else {
    if (flt->has_min && rdb_max < flt->min_val - 1e-6) return 0;
    if (flt->has_max && rdb_min > flt->max_val + 1e-6) return 0;
  }

  return 1;
}

/**
 * @brief Print command-line usage instructions to stdout.
 *
 * @param[in] argc Actual argument count received by `main()`, reported in the diagnostic line.
 */
static void print_usage(int argc)
{
  printf(
    "Too few arguments found (%d).\nUsage:\n"
    "calcmass [options] -mass [type] [Chemical formula]\n"
    "calcmass [options] -chem [type] [Chemical limits] [mass]\n"
    "calcmass [options] --chem-rdb [type] [Chemical limits] [mass]\n\n"
    "Options:\n"
    "  -s, --s          Short output (M/Z line only for -mass; table only for -chem)\n"
    "  --rdb-min <val>  Minimum allowed neutral-molecule RDB (default 0.0 for --chem-rdb)\n"
    "  --rdb-max <val>  Maximum allowed neutral-molecule RDB\n"
    "  --rdb-strict     Apply RDB bounds strictly to organic baseline (S=2, P=3)\n\n"
    "[type] is one of M, M+, M-, MH+, MNa+, MK+, MCa+, MH-,\n"
    "                 M2H++ (MH++, MH2+, MH+2), MHNa++, M2Na++\n"
    "[Chemical formula] is in the form C10,N6,O2,H13,Cl, C10H13N6O2Cl, or 'Ca(OH)2'\n"
    "                   (quote formulas containing parentheses in the shell)\n"
    "[Chemical limits]  is in the form C8-10,N6,O1-2,H1-13,Cl\n"
    "[mass]             is a positive number or +/- expression, e.g. '362.1003 - 242.0576'\n", argc);
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
  char rdb_buf[64];
  mass_t MT;
  int short_flag = 0;
  int npos = 0;
  int rc = 0;

  rdb_filter_t rdb_flt = { 0, 0, 0.0, 0, 0.0, 0 };

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
    } else if (strcmp(argv[a], "--rdb-strict") == 0) {
      rdb_flt.strict_min = 1;
    } else if (strcmp(argv[a], "--rdb-min") == 0) {
      if (a + 1 >= argc) {
        fprintf(stderr, "Missing value for --rdb-min\n");
        free(pos);
        return 1;
      }
      if (parse_double_arg(argv[++a], &rdb_flt.min_val) != 0) {
        fprintf(stderr, "Invalid --rdb-min value: '%s'\n", argv[a]);
        free(pos);
        return 1;
      }
      rdb_flt.has_min = 1;
    } else if (strcmp(argv[a], "--rdb-max") == 0) {
      if (a + 1 >= argc) {
        fprintf(stderr, "Missing value for --rdb-max\n");
        free(pos);
        return 1;
      }
      if (parse_double_arg(argv[++a], &rdb_flt.max_val) != 0) {
        fprintf(stderr, "Invalid --rdb-max value: '%s'\n", argv[a]);
        free(pos);
        return 1;
      }
      rdb_flt.has_max = 1;
    } else {
      pos[npos++] = argv[a];
    }
  }

  if (rdb_flt.has_min && rdb_flt.has_max && rdb_flt.min_val > rdb_flt.max_val) {
    fprintf(stderr, "--rdb-min (%.3f) cannot exceed --rdb-max (%.3f)\n",
            rdb_flt.min_val, rdb_flt.max_val);
    free(pos);
    return 1;
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
    double rdb_min, rdb_max;
    (void)Calc_RDB_Range(fi, &rdb_min, &rdb_max);
    format_rdb(rdb_min, rdb_max, rdb_buf, sizeof(rdb_buf));
    printf("\nParsed formula: '%s', RDB = %s\n", fs, rdb_buf);
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
      rdb_flt.require_int = 1;
      if (!rdb_flt.has_min) {
        rdb_flt.has_min = 1;
        rdb_flt.min_val = 0.0;
      }
      if (rdb_flt.has_max && rdb_flt.min_val > rdb_flt.max_val) {
        fprintf(stderr, "--rdb-min (%.3f) cannot exceed --rdb-max (%.3f)\n",
                rdb_flt.min_val, rdb_flt.max_val);
        free(pos);
        return 1;
      }
    }

    int rdb_active = rdb_flt.require_int || rdb_flt.has_min || rdb_flt.has_max;

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

    if (Parse_Mass_Expr(pos[3], &Target_Mass) != 0 || Target_Mass <= 0.0) {
      fprintf(stderr, "Invalid target mass: '%s'\n", pos[3]);
      free(pos);
      return 1;
    }
    if (!short_flag) {
      printf("Target mass is %10.6f\n", Target_Mass);
      printf("RDB check flag is %d\n", rdb_active);
    }

    int i = 0;
    long count = 0;
    chem = CHEM_LIST;
    for (Loop_Init(pfi, plim); Loop_End(pfi, plim); Loop_Next(pfi, plim)) {
      Loop_Copy(&(chem->f), pfi);
      chem->mass          = Calc_Mass_Charge_Ratio(chem->f, MT);
      chem->mass_diff     = chem->mass - Target_Mass;
      chem->mass_diff_ppm = (chem->mass_diff) / Target_Mass * 1000000.0;
      (void)Calc_RDB_Range(chem->f, &(chem->RDB_min), &(chem->RDB_max));
      if (rdb_passes(chem->RDB_min, chem->RDB_max, &rdb_flt)) {
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
      if (rdb_active) {
        int sep = 0;
        printf(" (neutral RDB filter: ");
        if (rdb_flt.require_int) { printf("integer"); sep = 1; }
        if (rdb_flt.has_min) { printf("%smin=%.1f", sep ? ", " : "", rdb_flt.min_val); sep = 1; }
        if (rdb_flt.has_max) { printf("%smax=%.1f", sep ? ", " : "", rdb_flt.max_val); sep = 1; }
        if (rdb_flt.strict_min) { printf("%sstrict S=2/P=3", sep ? ", " : ""); }
        printf("):\n");
      } else {
        printf(":\n");
      }
    }
    printf("%25s  %10s  %12s %11s  %4s\n", "Formula", "Theory", "Mass diff", "Mass diff", "RDB");
    for (int idx = 0; idx < n_print; idx++) {
      Print_Formula(CHEM_LIST[idx].f, fs, sizeof(fs));
      format_rdb(CHEM_LIST[idx].RDB_min, CHEM_LIST[idx].RDB_max, rdb_buf, sizeof(rdb_buf));
      printf("%25s: %7.6f  %8.3f mDa %7.3f ppm  %s\n", fs,
             CHEM_LIST[idx].mass,
             CHEM_LIST[idx].mass_diff * 1000.0,
             CHEM_LIST[idx].mass_diff_ppm,
             rdb_buf);
    }
    if (!short_flag && !rdb_flt.require_int) {
      printf("use --chem-rdb instead of --chem\nto ignore formulas with NEGATIVE or HALF-INTEGER neutral RDB values\n");
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
