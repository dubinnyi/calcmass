/**
 * @file isotopes.c
 * @brief Built-in isotope/adduct tables and mass, charge, m/z, and RDB calculations.
 *
 * Defines the static monoisotopic element database (`Isotopes[]`), adduct definitions
 * (`Mass[]`), and core numerical functions used both for single-formula mass calculation
 * and for scoring candidate formulas against a target mass.
 *
 * Note: Adduct mass deltas (`Mass[m].madd`) are lazily initialized via `mass_type_init()`
 * on first calculation if not already called explicitly at startup.
 */

#include <stdio.h>
#include <math.h>
#include "isotopes.h"

#define False 0
#define True  1

/* Monoisotopic neutral atomic masses from AME2020 (Wang et al., Chin. Phys. C 45, 030003, 2021)
 * and electron relative atomic mass from CODATA 2022 (in unified atomic mass units / Da). */
isotope_t Isotopes[NI] =
{
  { "e",  0,  0.000548579909, 0, False },  /*  0  ELECTRON    */
  { "C", 12, 12.000000000000, 4, True  },  /*  1  CARBON      */
  { "H",  1,  1.007825031898, 1, True  },  /*  2  PROTON      */
  { "N", 14, 14.00307400425,  3, True  },  /*  3  NITROGEN    */
  { "O", 16, 15.99491461926,  2, True  },  /*  4  OXYGEN      */
  { "P", 31, 30.97376199768,  5, False },  /*  5  PHOSPHORUS  */
  { "S", 32, 31.9720711744,   6, True  },  /*  6  SULFUR      */
  {"Na", 23, 22.98976928195,  1, True  },  /*  7  SODIUM      */
  {"Ca", 40, 39.962590863,    2, True  },  /*  8  CALCIUM     */
  { "K", 39, 38.963706486,    1, True  },  /*  9  POTASSIUM   */

  { "B", 11, 11.009305166,    3, True  },  /* 10  BORON       */
  { "F", 19, 18.99840316207,  1, True  },  /* 11  FLUORINE    */
  {"Cl", 35, 34.968852694,    1, True  },  /* 12  CHLORINE    */
  {"Br", 79, 78.918337574,    1, True  },  /* 13  BROMINE     */
};

int Mass_Init_Flag = 0;

mass_type_t Mass[NM] = {
/* name, elem1,count1,elem2,count2,elem3,count3, madd,charge */
 {"[M]",       0,   0,    0,     0,    0,     0,  0.0,  0},  /* MASS_M:       Neutral molecule    */
 {"[M+H]+",    2,   1,    0,    -1,    0,     0,  0.0, +1},  /* MASS_POS_H:   1H(+1), e(-1)       */
 {"[M+Na]+",   7,   1,    0,    -1,    0,     0,  0.0, +1},  /* MASS_POS_Na:  23Na(+1), e(-1)     */
 {"[M+Ca-H]+", 8,   1,    0,    -1,    2,    -1,  0.0, +1},  /* MASS_POS_Ca:  40Ca(+1), e(-1), H(-1) */
 {"[M-H]-",    2,  -1,    0,    +1,    0,     0,  0.0, -1},  /* MASS_NEG_H:   1H(-1), e(+1)       */
 {"[M]+",      0,  -1,    0,     0,    0,     0,  0.0, +1},  /* MASS_POS:     e(-1)               */
 {"[M]-",      0,  +1,    0,     0,    0,     0,  0.0, -1},  /* MASS_NEG:     e(+1)               */
 {"[M+K]+",    9,   1,    0,    -1,    0,     0,  0.0, +1},  /* MASS_POS_K:   39K(+1), e(-1)      */
 {"[M+2H]++",  2,   2,    0,    -2,    0,     0,  0.0, +2},  /* MASS_POS_2H:  1H(+2), e(-2)       */
 {"[M+H+Na]++",2,   1,    0,    -2,    7,     1,  0.0, +2},  /* MASS_POS_HNa: 1H(+1), 23Na(+1), e(-2) */
 {"[M+2Na]++", 7,   2,    0,    -2,    0,     0,  0.0, +2},  /* MASS_POS_2Na: 23Na(+2), e(-2)     */
};

/**
 * @brief Precompute exact mass offsets (`madd`) for every entry in `Mass[]`.
 *
 * Computes `Mass[m].madd` from the three (isotope index, signed count) pairs
 * defined in each `mass_type_t` entry and sets `Mass_Init_Flag = 1`.
 * Safe to call multiple times (idempotent).
 */
void mass_type_init(void)
{
  int m;
  for (m = 0; m < NM; m++) {
    Mass[m].madd =
      Isotopes[Mass[m].elem1].dmass * Mass[m].count1 +
      Isotopes[Mass[m].elem2].dmass * Mass[m].count2 +
      Isotopes[Mass[m].elem3].dmass * Mass[m].count3;
  }
  Mass_Init_Flag = 1;
}

/**
 * @brief Calculate the neutral monoisotopic mass of a chemical formula.
 *
 * Convenience wrapper equivalent to `Calc_Mass(fi, MASS_M)`.
 * Lazily invokes `mass_type_init()` if not yet initialized.
 *
 * @param[in] fi Composition vector of length at least `NI`; elements `[0, NI)`
 *               must be initialized (non-negative atom counts; index 0 is electron count).
 * @return Exact monoisotopic mass in Daltons (Da).
 */
double Calc_Mass_M(formula_t fi)
{
  if (!Mass_Init_Flag) mass_type_init();
  return Calc_Mass(fi, MASS_M);
}

/**
 * @brief Calculate the exact monoisotopic mass of a formula with an adduct applied.
 *
 * Sums `Isotopes[e].dmass * fi[e]` for all `e` in `[0, NI)` (including any explicit
 * electron adjustment at `fi[0]`) and adds the adduct offset `Mass[M].madd`.
 * Lazily invokes `mass_type_init()` if not yet initialized.
 *
 * @param[in] fi Composition vector with initialized entries in `[0, NI)`.
 * @param[in] M  Adduct type selector; must be a valid `mass_t` value in `[0, NM)`.
 * @return Total monoisotopic mass of the species in Daltons (Da), not divided by charge.
 */
double Calc_Mass(formula_t fi, mass_t M)
{
  int e;
  double mass = 0.0;

  if (!Mass_Init_Flag) mass_type_init();
  for (e = 0; e < NI; e++) mass += Isotopes[e].dmass * fi[e];
  mass += Mass[M].madd;
  return mass;
}

/**
 * @brief Calculate Ring and Double Bond (RDB) equivalents for a chemical formula.
 *
 * Uses the standard formula: `RDB = 1 + sum_{e=1}^{NI-1} ((valence_e - 2) / 2) * count_e`.
 * Index 0 (electrons) is intentionally skipped so explicit ion charge does not skew RDB.
 * Lazily invokes `mass_type_init()` if not yet initialized.
 *
 * @param[in]  fi  Composition vector with initialized entries in `[0, NI)`.
 * @param[out] rdb Pointer to receive the computed RDB value; must be non-NULL.
 * @return 1 (`True`) if every present element has an unambiguous valence (`v_unambig != 0`);
 *         0 (`False`) if any present element has variable valence (e.g., phosphorus).
 */
int Calc_RDB(formula_t fi, double *rdb)
{
  *rdb = 1.0;
  int unambig_flag = True;
  int e;

  if (!Mass_Init_Flag) mass_type_init();
  /* Start at e = 1 to skip electrons (index 0) */
  for (e = 1; e < NI; e++) {
    (*rdb) += (Isotopes[e].valence - 2.0) / 2.0 * fi[e];
    if (fi[e])
      unambig_flag &= Isotopes[e].v_unambig;
  }
  return unambig_flag;
}

/**
 * @brief Calculate the total net electric charge of a formula and adduct.
 *
 * Combines the intrinsic charge implied by explicit electrons in the formula
 * (`-fi[0]`) with the adduct charge (`Mass[M].charge`).
 * Lazily invokes `mass_type_init()` if not yet initialized.
 *
 * @param[in] fi Composition vector with initialized `fi[0]` (excess electrons).
 * @param[in] M  Adduct type selector; must be a valid `mass_t` value in `[0, NM)`.
 * @return Signed net charge as a floating-point integer (e.g., `0.0`, `+1.0`, `-1.0`, `+2.0`).
 */
double Calc_Charge(formula_t fi, mass_t M)
{
  double charge = 0.0;

  if (!Mass_Init_Flag) mass_type_init();
  charge = -fi[0];
  charge += Mass[M].charge;
  return charge;
}

/**
 * @brief Calculate the mass-to-charge ratio (m/z) of a formula with an adduct.
 *
 * Divides `Calc_Mass(fi, M)` by `fabs(Calc_Charge(fi, M))`. If the net charge is
 * zero (e.g., neutral molecule with `MASS_M`), returns the undivided monoisotopic
 * mass directly to avoid division by zero.
 *
 * @param[in] fi Composition vector with initialized entries in `[0, NI)`.
 * @param[in] M  Adduct type selector; must be a valid `mass_t` value in `[0, NM)`.
 * @return Theoretical m/z value in Da/e, or neutral mass in Da when net charge is 0.
 */
double Calc_Mass_Charge_Ratio(formula_t fi, mass_t M)
{
  double mass = Calc_Mass(fi, M);
  double charge = Calc_Charge(fi, M);

  if (charge == 0.0)
    return mass;
  return mass / fabs(charge);
}

/**
 * @brief Standard `qsort` comparator for `fmass_t` entries by absolute mass error.
 *
 * Establishes a strict weak ordering on `fabs(mass_diff)` so the closest
 * theoretical formulas to the target mass sort first.
 *
 * @param[in] fm1 Pointer to the first `fmass_t` element; must be non-NULL.
 * @param[in] fm2 Pointer to the second `fmass_t` element; must be non-NULL.
 * @return `-1` if `|fm1->mass_diff| < |fm2->mass_diff|`,
 *         `+1` if `|fm1->mass_diff| > |fm2->mass_diff|`,
 *          `0` if equal.
 */
int Mass_Diff(const void *fm1, const void *fm2)
{
  double d1 = fabs(((const fmass_t *)fm1)->mass_diff);
  double d2 = fabs(((const fmass_t *)fm2)->mass_diff);

  if (d1 < d2) return -1;
  if (d1 > d2) return 1;
  return 0;
}
