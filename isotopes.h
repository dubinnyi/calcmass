/**
 * @file isotopes.h
 * @brief Core data structures, constants, and public API for calcmass.
 *
 * Declares the built-in monoisotopic element table, supported mass-spectrometry
 * ion/adduct types, chemical formula and search-limit representations, mass/m/z/RDB
 * calculators, string parsers, and the multi-dimensional formula iterator.
 *
 * Array layout convention for formula_t and limits_t:
 * - Index 0        : Electron count (`e`). Positive = extra electrons (negative charge);
 *                    negative = removed electrons (positive charge).
 * - Indices 1..NI-1: Atom counts for chemical elements in `Isotopes[]` order.
 * - Index NI       : Overflow/carry sentinel slot used by the multiloop iterator.
 */

#ifndef ISOTOPES_H_INCLUDED
#define ISOTOPES_H_INCLUDED

/** @brief Number of entries in the built-in Isotopes[] table (electron + 13 elements). */
#define NI 14

/**
 * @brief Monoisotopic mass and valence metadata for a single element or particle.
 */
typedef struct
{
  char   name[3];    /**< NUL-terminated element symbol (e.g., "C", "Cl", "e"). */
  int    imass;      /**< Nominal integer mass number (0 for electron). */
  double dmass;      /**< Exact monoisotopic mass in Daltons (Da). */
  int    valence;    /**< Standard chemical valence used in RDB calculations. */
  int    v_unambig;  /**< Non-zero if valence is unambiguous; 0 for variable-valence elements (e.g., P). */
} isotope_t;

/** @brief Built-in table of supported isotopes/elements indexed by element ID [0, NI). */
extern isotope_t Isotopes[NI];

/**
 * @brief Chemical composition vector indexed by element ID [0, NI), plus sentinel at [NI].
 * @see isotopes.h file documentation for index 0 (electrons) and index NI (sentinel).
 */
typedef int formula_t[NI+1];

/** @brief Inclusive lower and upper atom count bounds for a single element. */
typedef struct {
  int min; /**< Minimum atom count (invariant after parsing: 0 <= min <= max). */
  int max; /**< Maximum atom count (invariant after parsing: max >= min). */
} int_limits_t;

/** @brief Per-element search bounds indexed by element ID [0, NI), plus unused slot [NI]. */
typedef int_limits_t limits_t[NI+1];

/**
 * @brief Candidate chemical formula evaluated against a target mass/m/z.
 */
typedef struct {
  formula_t f;            /**< Evaluated chemical formula. */
  double mass;            /**< Theoretical m/z (or neutral mass for MASS_M) in Da. */
  double mass_diff;       /**< Signed mass error (theory - target) in Da. */
  double mass_diff_ppm;   /**< Signed relative mass error in parts per million (ppm). */
  double RDB;             /**< Ring and Double Bond equivalent. Neutral closed-shell molecules
                               normally have an integer, non-negative RDB; half-integer or
                               negative values indicate unsatisfied valence rules. */
} fmass_t;

/** @brief Number of supported ion/adduct types in Mass[]. */
#define NM 11

/**
 * @brief Identifiers for supported neutral and ionized adduct species.
 */
typedef enum {
  MASS_M,        /**< Neutral molecule [M] (charge 0). */
  MASS_POS_H,    /**< Protonated ion [M+H]+ (charge +1). */
  MASS_POS_Na,   /**< Sodiated ion [M+Na]+ (charge +1). */
  MASS_POS_Ca,   /**< Calciated ion [M+Ca-H]+ (charge +1). */
  MASS_NEG_H,    /**< Deprotonated ion [M-H]- (charge -1). */
  MASS_POS,      /**< Radical cation [M]+ (charge +1). */
  MASS_NEG,      /**< Radical anion [M]- (charge -1). */
  MASS_POS_K,    /**< Potassiated ion [M+K]+ (charge +1). */
  MASS_POS_2H,   /**< Doubly protonated ion [M+2H]++ (charge +2). */
  MASS_POS_HNa,  /**< Doubly charged proton/sodium adduct [M+H+Na]++ (charge +2). */
  MASS_POS_2Na,  /**< Doubly sodiated ion [M+2Na]++ (charge +2). */
} mass_t;

/**
 * @brief Definition of an ion/adduct modification relative to neutral molecule M.
 *
 * Up to three isotope adjustments (elem1..3 with signed count1..3) are combined
 * into the cached mass offset `madd` by mass_type_init().
 */
typedef struct {
  char name[20];          /**< Human-readable adduct label (e.g., "[M+H]+"). */
  int elem1, count1;      /**< First adduct isotope index [0, NI) and signed count. */
  int elem2, count2;      /**< Second adduct isotope index [0, NI) and signed count. */
  int elem3, count3;      /**< Third adduct isotope index [0, NI) and signed count. */
  double madd;            /**< Precomputed monoisotopic mass delta in Da (initialized by mass_type_init()). */
  int    charge;          /**< Intrinsic net integer charge contributed by the adduct. */
} mass_type_t;

/** @brief Table of adduct definitions indexed by mass_t [0, NM). */
extern mass_type_t Mass[NM];

/** @brief Non-zero once mass_type_init() has populated Mass[m].madd for all adducts. */
extern int Mass_Init_Flag;

/* --- Mass, charge, RDB, and sorting (isotopes.c) --- */

/** @brief Precompute adduct mass offsets Mass[m].madd for all ion types. */
void   mass_type_init(void);

/** @brief Calculate the monoisotopic mass (Da) of formula `fi` with adduct `M`. */
double Calc_Mass(formula_t fi, mass_t M);

/** @brief Calculate the neutral monoisotopic mass (Da) of formula `fi` (MASS_M). */
double Calc_Mass_M(formula_t fi);

/** @brief Calculate Ring and Double Bond (RDB) equivalents for formula `fi`. */
int    Calc_RDB(formula_t fi, double *rdb);

/** @brief Calculate net electric charge of formula `fi` combined with adduct `M`. */
double Calc_Charge(formula_t fi, mass_t M);

/** @brief Calculate m/z ratio (or neutral mass if net charge is 0) of `fi` with adduct `M`. */
double Calc_Mass_Charge_Ratio(formula_t fi, mass_t M);

/** @brief qsort comparator ordering fmass_t entries by ascending |mass_diff|. */
int    Mass_Diff(const void *fm1, const void *fm2);

/* --- Formula and search-limit parsing (parse.c) --- */

/** @brief Parse a chemical formula string into composition vector `fi`. */
void Parse_Formula(const char *fc, formula_t fi);

/** @brief Parse an element-range limits string into bounds array `lim`. */
void Parse_Limits(const char *fc, limits_t lim);

/* --- Multi-dimensional formula enumeration (multiloop.c) --- */

/** @brief Copy an entire formula_t vector (including sentinel slot). */
void Loop_Copy(formula_t *dst, formula_t *src);

/** @brief Initialize iterator state `ft` to the lower bounds in `lim`. */
void Loop_Init(formula_t *ft, limits_t *lim);

/** @brief Advance iterator state `ft` to the next combination within `lim`. */
void Loop_Next(formula_t *ft, limits_t *lim);

/** @brief Check whether the iterator `ft` still holds a valid combination. */
int  Loop_End(formula_t *ft, limits_t *lim);

#endif /* ISOTOPES_H_INCLUDED */
