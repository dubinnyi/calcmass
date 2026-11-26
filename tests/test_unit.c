/**
 * @file tests/test_unit.c
 * @brief Unit tests for calcmass core modules (isotopes.c, parse.c, multiloop.c).
 */

#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <unistd.h>
#include <fcntl.h>

#include "../isotopes.h"

static int g_pass = 0;
static int g_fail = 0;

/* Temporarily suppress expected parser diagnostic messages on stderr during negative tests. */
static int silence_stderr_begin(void)
{
  fflush(stderr);
  int saved = dup(STDERR_FILENO);
  if (saved < 0) return -1;
  int devnull = open("/dev/null", O_WRONLY);
  if (devnull >= 0) {
    (void)dup2(devnull, STDERR_FILENO);
    close(devnull);
  }
  return saved;
}

static void silence_stderr_end(int saved)
{
  if (saved < 0) return;
  fflush(stderr);
  (void)dup2(saved, STDERR_FILENO);
  close(saved);
}

#define ASSERT_TRUE(expr) do { \
  if (expr) { \
    g_pass++; \
  } else { \
    fprintf(stderr, "FAIL [%s:%d]: %s\n", __FILE__, __LINE__, #expr); \
    g_fail++; \
  } \
} while (0)

#define ASSERT_EQ_INT(actual, expected) do { \
  int _a = (actual), _e = (expected); \
  if (_a == _e) { \
    g_pass++; \
  } else { \
    fprintf(stderr, "FAIL [%s:%d]: %s == %d, expected %d\n", \
            __FILE__, __LINE__, #actual, _a, _e); \
    g_fail++; \
  } \
} while (0)

#define ASSERT_NEAR(actual, expected, tol) do { \
  double _a = (actual), _e = (expected); \
  if (fabs(_a - _e) <= (tol)) { \
    g_pass++; \
  } else { \
    fprintf(stderr, "FAIL [%s:%d]: %s == %.12f, expected %.12f (diff %.3e > %.3e)\n", \
            __FILE__, __LINE__, #actual, _a, _e, fabs(_a - _e), (double)(tol)); \
    g_fail++; \
  } \
} while (0)

static void test_adducts_and_masses(void)
{
  const double me  = 0.000548579909;
  const double mC  = 12.0;
  const double mH  = 1.007825031898;
  const double mN  = 14.00307400425;
  const double mO  = 15.99491461926;
  const double mNa = 22.98976928195;
  const double mCa = 39.962590863;
  const double mK  = 38.963706486;

  mass_type_init();
  ASSERT_EQ_INT(Mass_Init_Flag, 1);

  /* Idempotent re-init */
  mass_type_init();
  ASSERT_NEAR(Mass[MASS_M].madd,       0.0,                 1e-12);
  ASSERT_NEAR(Mass[MASS_POS_H].madd,   mH - me,             1e-12);
  ASSERT_NEAR(Mass[MASS_NEG_H].madd,  -mH + me,             1e-12);
  ASSERT_NEAR(Mass[MASS_POS_Na].madd,  mNa - me,            1e-12);
  ASSERT_NEAR(Mass[MASS_POS_K].madd,   mK - me,             1e-12);
  ASSERT_NEAR(Mass[MASS_POS_Ca].madd,  mCa - mH - me,       1e-12);
  ASSERT_NEAR(Mass[MASS_POS].madd,    -me,                  1e-12);
  ASSERT_NEAR(Mass[MASS_NEG].madd,     me,                  1e-12);
  ASSERT_NEAR(Mass[MASS_POS_2H].madd,  2.0 * (mH - me),     1e-12);
  ASSERT_NEAR(Mass[MASS_POS_HNa].madd, mH + mNa - 2.0 * me, 1e-12);
  ASSERT_NEAR(Mass[MASS_POS_2Na].madd, 2.0 * (mNa - me),    1e-12);

  formula_t fi;
  ASSERT_EQ_INT(Parse_Formula("C16H17N3O7", fi), 0);
  double m_neutral = 16.0 * mC + 17.0 * mH + 3.0 * mN + 7.0 * mO;
  ASSERT_NEAR(Calc_Mass_M(fi),                        m_neutral,                       1e-10);
  ASSERT_NEAR(Calc_Charge(fi, MASS_M),                0.0,                             1e-12);
  ASSERT_NEAR(Calc_Mass_Charge_Ratio(fi, MASS_M),     m_neutral,                       1e-10);
  ASSERT_NEAR(Calc_Mass_Charge_Ratio(fi, MASS_NEG_H), m_neutral - mH + me,             1e-10);
  ASSERT_NEAR(Calc_Mass_Charge_Ratio(fi, MASS_POS_2H),(m_neutral + 2.0*(mH - me))/2.0, 1e-10);

  /* Explicit charge in formula */
  ASSERT_EQ_INT(Parse_Formula("C10H14+2", fi), 0);
  ASSERT_EQ_INT(fi[0], -2);
  ASSERT_NEAR(Calc_Charge(fi, MASS_M),            2.0,                         1e-12);
  ASSERT_NEAR(Calc_Mass_Charge_Ratio(fi, MASS_M), (10.0*mC + 14.0*mH - 2*me)/2.0, 1e-10);
}

static void test_rdb_ranges(void)
{
  formula_t fi;
  double rmin = 0.0, rmax = 0.0, rsingle = 0.0;

  /* Unambiguous organic: C16H17N3O7 -> 10.0 */
  ASSERT_EQ_INT(Parse_Formula("C16H17N3O7", fi), 0);
  ASSERT_EQ_INT(Calc_RDB_Range(fi, &rmin, &rmax), 1);
  ASSERT_NEAR(rmin, 10.0, 1e-9);
  ASSERT_NEAR(rmax, 10.0, 1e-9);
  ASSERT_EQ_INT(Calc_RDB(fi, &rsingle), 1);
  ASSERT_NEAR(rsingle, 10.0, 1e-9);

  /* Explicit ion charge does not alter RDB */
  ASSERT_EQ_INT(Parse_Formula("C6H6+2", fi), 0);
  ASSERT_EQ_INT(Calc_RDB_Range(fi, &rmin, &rmax), 1);
  ASSERT_NEAR(rmin, 4.0, 1e-9);
  ASSERT_NEAR(rmax, 4.0, 1e-9);

  /* Cysteine (C3H7NO2S): S=2 -> 1.0, S=6 -> 3.0 */
  ASSERT_EQ_INT(Parse_Formula("C3H7NO2S", fi), 0);
  ASSERT_EQ_INT(Calc_RDB_Range(fi, &rmin, &rmax), 0);
  ASSERT_NEAR(rmin, 1.0, 1e-9);
  ASSERT_NEAR(rmax, 3.0, 1e-9);

  /* DMSO ((CH3)2SO): S=2 -> 0.0, S=6 -> 2.0 */
  ASSERT_EQ_INT(Parse_Formula("(CH3)2SO", fi), 0);
  ASSERT_EQ_INT(Calc_RDB_Range(fi, &rmin, &rmax), 0);
  ASSERT_NEAR(rmin, 0.0, 1e-9);
  ASSERT_NEAR(rmax, 2.0, 1e-9);

  /* SF6 (hypervalent S=6): S=2 -> -2.0, S=6 -> 0.0 */
  ASSERT_EQ_INT(Parse_Formula("SF6", fi), 0);
  ASSERT_EQ_INT(Calc_RDB_Range(fi, &rmin, &rmax), 0);
  ASSERT_NEAR(rmin, -2.0, 1e-9);
  ASSERT_NEAR(rmax,  0.0, 1e-9);

  /* Impossible H8S: S=2 -> -3.0, S=6 -> -1.0 */
  ASSERT_EQ_INT(Parse_Formula("H8S", fi), 0);
  ASSERT_EQ_INT(Calc_RDB_Range(fi, &rmin, &rmax), 0);
  ASSERT_NEAR(rmin, -3.0, 1e-9);
  ASSERT_NEAR(rmax, -1.0, 1e-9);

  /* Three sulfur atoms (C3H10S3): base=-1.0, S=2..6 shifts by 0..+6 -> [-1.0, 5.0] */
  ASSERT_EQ_INT(Parse_Formula("C3H10S3", fi), 0);
  ASSERT_EQ_INT(Calc_RDB_Range(fi, &rmin, &rmax), 0);
  ASSERT_NEAR(rmin, -1.0, 1e-9);
  ASSERT_NEAR(rmax,  5.0, 1e-9);

  /* Phosphorus: H3PO4 -> [0.0, 1.0], PF5 -> [-1.0, 0.0] */
  ASSERT_EQ_INT(Parse_Formula("H3PO4", fi), 0);
  ASSERT_EQ_INT(Calc_RDB_Range(fi, &rmin, &rmax), 0);
  ASSERT_NEAR(rmin, 0.0, 1e-9);
  ASSERT_NEAR(rmax, 1.0, 1e-9);

  ASSERT_EQ_INT(Parse_Formula("PF5", fi), 0);
  ASSERT_EQ_INT(Calc_RDB_Range(fi, &rmin, &rmax), 0);
  ASSERT_NEAR(rmin, -1.0, 1e-9);
  ASSERT_NEAR(rmax,  0.0, 1e-9);

  /* Combined P + S (C2H7O2PS): base=-0.5, P adds 0.5..1.5, S adds 0..2 -> [0.0, 3.0] */
  ASSERT_EQ_INT(Parse_Formula("C2H7O2PS", fi), 0);
  ASSERT_EQ_INT(Calc_RDB_Range(fi, &rmin, &rmax), 0);
  ASSERT_NEAR(rmin, 0.0, 1e-9);
  ASSERT_NEAR(rmax, 3.0, 1e-9);

  /* Even-electron fragment ion composition: CH3 -> 0.5, NH4 -> -0.5 */
  ASSERT_EQ_INT(Parse_Formula("CH3", fi), 0);
  ASSERT_EQ_INT(Calc_RDB_Range(fi, &rmin, &rmax), 1);
  ASSERT_NEAR(rmin, 0.5, 1e-9);
  ASSERT_NEAR(rmax, 0.5, 1e-9);

  ASSERT_EQ_INT(Parse_Formula("NH4", fi), 0);
  ASSERT_EQ_INT(Calc_RDB_Range(fi, &rmin, &rmax), 1);
  ASSERT_NEAR(rmin, -0.5, 1e-9);
  ASSERT_NEAR(rmax, -0.5, 1e-9);
}

static void test_parsers(void)
{
  formula_t fi;
  limits_t lim;

  /* Formula positive cases */
  ASSERT_EQ_INT(Parse_Formula("C2H5OH", fi), 0);
  ASSERT_EQ_INT(fi[1], 2); /* C */
  ASSERT_EQ_INT(fi[2], 6); /* H */
  ASSERT_EQ_INT(fi[4], 1); /* O */

  ASSERT_EQ_INT(Parse_Formula("Ca(OH)2", fi), 0);
  ASSERT_EQ_INT(fi[8], 1); /* Ca */
  ASSERT_EQ_INT(fi[4], 2); /* O */
  ASSERT_EQ_INT(fi[2], 2); /* H */

  ASSERT_EQ_INT(Parse_Formula("((CH2)2O)3H2O", fi), 0);
  ASSERT_EQ_INT(fi[1], 6);  /* C */
  ASSERT_EQ_INT(fi[2], 14); /* H */
  ASSERT_EQ_INT(fi[4], 4);  /* O */

  ASSERT_EQ_INT(Parse_Formula("c10, h13, cl", fi), 0);
  ASSERT_EQ_INT(fi[1], 10);
  ASSERT_EQ_INT(fi[2], 13);
  ASSERT_EQ_INT(fi[12], 1);

  /* Formula negative cases (suppress expected stderr diagnostics) */
  int saved_err = silence_stderr_begin();
  int r_null = Parse_Formula(NULL, fi);
  int r_empty = Parse_Formula("", fi);
  int r_ws = Parse_Formula("  ,  ", fi);
  int r_si = Parse_Formula("C10Si", fi);
  int fi_c_after_err = fi[1];
  int r_foo = Parse_Formula("Foo", fi);
  int r_unclosed = Parse_Formula("Ca(OH2", fi);
  int r_unmatched = Parse_Formula("CaOH)2", fi);
  int r_emptyparen = Parse_Formula("Ca()2", fi);
  int r_leadnum = Parse_Formula("10CH4", fi);
  int r_badchar = Parse_Formula("C10@H14", fi);
  int r_ovf1 = Parse_Formula("C99999999999999", fi);
  int r_ovf2 = Parse_Formula("C2000000000C2000000000", fi);
  int r_ovf3 = Parse_Formula("(C1000000)1000000", fi);
  int r_plusonly = Parse_Formula("+", fi);
  int r_minusonly = Parse_Formula("-2", fi);
  int r_zeroatoms = Parse_Formula("C0", fi);
  int r_zeromult = Parse_Formula("(OH)0", fi);
  silence_stderr_end(saved_err);

  ASSERT_EQ_INT(r_null, -1);
  ASSERT_EQ_INT(r_empty, -1);
  ASSERT_EQ_INT(r_ws, -1);
  ASSERT_EQ_INT(r_si, -1);
  ASSERT_EQ_INT(fi_c_after_err, 0); /* zeroed on parse failure */
  ASSERT_EQ_INT(r_foo, -1);
  ASSERT_EQ_INT(r_unclosed, -1);
  ASSERT_EQ_INT(r_unmatched, -1);
  ASSERT_EQ_INT(r_emptyparen, -1);
  ASSERT_EQ_INT(r_leadnum, -1);
  ASSERT_EQ_INT(r_badchar, -1);
  ASSERT_EQ_INT(r_ovf1, -1);
  ASSERT_EQ_INT(r_ovf2, -1);
  ASSERT_EQ_INT(r_ovf3, -1);
  ASSERT_EQ_INT(r_plusonly, -1);
  ASSERT_EQ_INT(r_minusonly, -1);
  ASSERT_EQ_INT(r_zeroatoms, -1);
  ASSERT_EQ_INT(r_zeromult, -1);

  /* Limits positive cases */
  ASSERT_EQ_INT(Parse_Limits("C8-10,N6,O1-2,H1-13,Cl", lim), 0);
  ASSERT_EQ_INT(lim[1].min, 8);  ASSERT_EQ_INT(lim[1].max, 10);
  ASSERT_EQ_INT(lim[3].min, 6);  ASSERT_EQ_INT(lim[3].max, 6);
  ASSERT_EQ_INT(lim[4].min, 1);  ASSERT_EQ_INT(lim[4].max, 2);
  ASSERT_EQ_INT(lim[2].min, 1);  ASSERT_EQ_INT(lim[2].max, 13);
  ASSERT_EQ_INT(lim[12].min, 1); ASSERT_EQ_INT(lim[12].max, 1);

  /* Reversed bounds swap */
  ASSERT_EQ_INT(Parse_Limits("C10-8", lim), 0);
  ASSERT_EQ_INT(lim[1].min, 8);
  ASSERT_EQ_INT(lim[1].max, 10);

  /* Limits negative cases */
  saved_err = silence_stderr_begin();
  int rl_null = Parse_Limits(NULL, lim);
  int rl_empty = Parse_Limits("", lim);
  int rl_si = Parse_Limits("C1-5,Si0-2", lim);
  int lim_c_after_err = lim[1].max;
  int rl_nolow = Parse_Limits("C-10", lim);
  int rl_nohigh = Parse_Limits("C8-", lim);
  int rl_twohyph = Parse_Limits("C8-10-12", lim);
  int rl_num = Parse_Limits("10-20", lim);
  int rl_ovf = Parse_Limits("C0-999999999999999", lim);
  silence_stderr_end(saved_err);

  ASSERT_EQ_INT(rl_null, -1);
  ASSERT_EQ_INT(rl_empty, -1);
  ASSERT_EQ_INT(rl_si, -1);
  ASSERT_EQ_INT(lim_c_after_err, 0); /* zeroed on parse failure */
  ASSERT_EQ_INT(rl_nolow, -1);
  ASSERT_EQ_INT(rl_nohigh, -1);
  ASSERT_EQ_INT(rl_twohyph, -1);
  ASSERT_EQ_INT(rl_num, -1);
  ASSERT_EQ_INT(rl_ovf, -1);

  /* Target mass arithmetic expression parser */
  double mval = 0.0;
  ASSERT_EQ_INT(Parse_Mass_Expr("362.1003", &mval), 0);
  ASSERT_NEAR(mval, 362.1003, 1e-12);

  ASSERT_EQ_INT(Parse_Mass_Expr("362.1003 - 242.0576", &mval), 0);
  ASSERT_NEAR(mval, 120.0427, 1e-12);

  ASSERT_EQ_INT(Parse_Mass_Expr("362.1003-242.0576", &mval), 0);
  ASSERT_NEAR(mval, 120.0427, 1e-12);

  /* Unicode minus U+2212 (\xE2\x88\x92) and en-dash U+2013 (\xE2\x80\x93) */
  ASSERT_EQ_INT(Parse_Mass_Expr("362.1003 \xE2\x88\x92 242.0576", &mval), 0);
  ASSERT_NEAR(mval, 120.0427, 1e-12);

  ASSERT_EQ_INT(Parse_Mass_Expr("362.1003 \xE2\x80\x93 242.0576", &mval), 0);
  ASSERT_NEAR(mval, 120.0427, 1e-12);

  ASSERT_EQ_INT(Parse_Mass_Expr("242.0576 + 120.0427", &mval), 0);
  ASSERT_NEAR(mval, 362.1003, 1e-12);

  ASSERT_EQ_INT(Parse_Mass_Expr("362.1003 - 242.0576 + 18.010565", &mval), 0);
  ASSERT_NEAR(mval, 138.053265, 1e-12);

  ASSERT_EQ_INT(Parse_Mass_Expr("-242.0576 + 362.1003", &mval), 0);
  ASSERT_NEAR(mval, 120.0427, 1e-12);

  ASSERT_EQ_INT(Parse_Mass_Expr("362.1003 - 1.0e-3", &mval), 0);
  ASSERT_NEAR(mval, 362.0993, 1e-12);

  ASSERT_EQ_INT(Parse_Mass_Expr(NULL, &mval), -1);
  ASSERT_EQ_INT(Parse_Mass_Expr("", &mval), -1);
  ASSERT_EQ_INT(Parse_Mass_Expr("   ", &mval), -1);
  ASSERT_EQ_INT(Parse_Mass_Expr("362.1003 -", &mval), -1);
  ASSERT_EQ_INT(Parse_Mass_Expr("362.1003 + + 1", &mval), -1);
  ASSERT_EQ_INT(Parse_Mass_Expr("362.1003 242.0576", &mval), -1);
  ASSERT_EQ_INT(Parse_Mass_Expr("abc", &mval), -1);
  ASSERT_EQ_INT(Parse_Mass_Expr("nan", &mval), -1);
  ASSERT_EQ_INT(Parse_Mass_Expr("inf", &mval), -1);
  ASSERT_EQ_INT(Parse_Mass_Expr("0x10", &mval), -1);
}

static void test_multiloop_and_comparator(void)
{
  limits_t lim;
  formula_t fi, copy;
  ASSERT_EQ_INT(Parse_Limits("C1-3,H0-2,O1-2", lim), 0);

  int visited = 0;
  for (Loop_Init(&fi, &lim); Loop_End(&fi, &lim); Loop_Next(&fi, &lim)) {
    ASSERT_TRUE(fi[1] >= 1 && fi[1] <= 3);
    ASSERT_TRUE(fi[2] >= 0 && fi[2] <= 2);
    ASSERT_TRUE(fi[4] >= 1 && fi[4] <= 2);
    Loop_Copy(&copy, &fi);
    ASSERT_EQ_INT(copy[1], fi[1]);
    ASSERT_EQ_INT(copy[2], fi[2]);
    ASSERT_EQ_INT(copy[4], fi[4]);
    visited++;
  }
  /* 3 * 3 * 2 = 18 combinations */
  ASSERT_EQ_INT(visited, 18);

  /* INT_MAX boundary in multiloop must terminate cleanly without signed overflow */
  memset(&lim, 0, sizeof(lim));
  lim[1].min = INT_MAX - 1;
  lim[1].max = INT_MAX;
  visited = 0;
  for (Loop_Init(&fi, &lim); Loop_End(&fi, &lim); Loop_Next(&fi, &lim)) {
    ASSERT_TRUE(fi[1] == INT_MAX - 1 || fi[1] == INT_MAX);
    visited++;
  }
  ASSERT_EQ_INT(visited, 2);

  /* Strict weak ordering and deterministic tie-breaking of Mass_Diff */
  fmass_t a, b, c, d;
  memset(&a, 0, sizeof(a));
  memset(&b, 0, sizeof(b));
  memset(&c, 0, sizeof(c));
  memset(&d, 0, sizeof(d));
  a.mass_diff = -0.001;
  b.mass_diff =  0.002;
  c.mass_diff =  0.001;
  d.mass_diff =  0.001;
  d.f[1] = 5;

  ASSERT_EQ_INT(Mass_Diff(&a, &b), -1);
  ASSERT_EQ_INT(Mass_Diff(&b, &a),  1);
  ASSERT_EQ_INT(Mass_Diff(&a, &c),  0);
  ASSERT_EQ_INT(Mass_Diff(&a, &a),  0);
  ASSERT_EQ_INT(Mass_Diff(&a, &d), -1);
  ASSERT_EQ_INT(Mass_Diff(&d, &a),  1);
}

int main(void)
{
  test_adducts_and_masses();
  test_rdb_ranges();
  test_parsers();
  test_multiloop_and_comparator();

  printf("Unit tests: %d passed, %d failed\n", g_pass, g_fail);
  return (g_fail == 0) ? 0 : 1;
}
