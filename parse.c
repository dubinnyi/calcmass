/**
 * @file parse.c
 * @brief Parsers for chemical formula strings and element-range search limits.
 *
 * Supports both comma-separated notation (e.g., "C10,N6,O2,H13,Cl", "C8-10,H1-13,Cl")
 * and compact standard chemical notation (e.g., "C10H13N6O2Cl", "C2H5OH").
 * Repeated elements in a formula are accumulated, and unknown element symbols are
 * reported to stderr and ignored without aborting.
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "isotopes.h"

#define DEBUG 0
#define IFDEBUG if(DEBUG)

/**
 * @brief Test whether a character represents an ionic charge sign.
 *
 * @param[in] c Character to inspect.
 * @return Non-zero (1) if `c` is `'+'` or `'-'`; 0 otherwise.
 */
static int issign(char c)
{
  return (c == '+' || c == '-');
}

/**
 * @brief Consume a non-negative decimal integer from the front of a string pointer.
 *
 * Advances `*pfc` past all consecutive ASCII digit characters. Excess digits beyond
 * the internal 31-character buffer are still consumed from `*pfc` to keep parsing
 * synchronized, and parsed values exceeding `INT_MAX` are clamped to `INT_MAX`.
 *
 * @param[in,out] pfc         Address of the current cursor in the input string;
 *                            must be non-NULL and `*pfc` must point to a NUL-terminated string.
 * @param[in]     default_val Value returned when `*pfc` does not start with a digit.
 * @return Parsed non-negative integer clamped to `INT_MAX`, or `default_val` if no digits were present.
 */
static int parse_uint(const char **pfc, int default_val)
{
  char digit[32];
  size_t dlen = 0;
  const char *fc = *pfc;

  while (*fc && isdigit((unsigned char)*fc)) {
    if (dlen < sizeof(digit) - 1) {
      digit[dlen++] = *fc;
    }
    fc++;
  }
  digit[dlen] = '\0';
  *pfc = fc;

  if (dlen == 0) {
    return default_val;
  }

  long val = strtol(digit, NULL, 10);
  if (val > INT_MAX) val = INT_MAX;
  return (int)val;
}

/**
 * @brief Parse a chemical formula string into a `formula_t` composition vector.
 *
 * Zero-initializes all `NI + 1` entries of `fi` before parsing. Recognizes:
 * - 1- or 2-letter element symbols in `Isotopes[]` (e.g., `C`, `H`, `Cl`, `Na`),
 *   separated by delimiters (such as commas/spaces) or uppercase boundaries.
 * - Optional integer counts (omitted count defaults to 1).
 * - Repeated occurrences of the same element, whose counts are summed (e.g., `C2H5OH` -> `C2,H6,O`).
 * - Charge suffixes `'+'` and `'-'` with optional multiplicity (e.g., `+`, `+2`, `-`),
 *   mapped to electron count `fi[0]` (`+` subtracts electrons, `-` adds electrons).
 *
 * Unknown element symbols are logged to `stderr` and skipped.
 *
 * @param[in]  fc NUL-terminated formula string to parse; must be non-NULL.
 * @param[out] fi Destination composition vector; overwritten completely on entry.
 */
void Parse_Formula(const char *fc, formula_t fi)
{
  char elem[3];
  int i, e, id, E;
  int parse_sign;

  for (e = 0; e <= NI; e++) fi[e] = 0;

  IFDEBUG printf("\nParsing formula string: '%s'\n\n", fc);
  while (*fc) {
    /* Skip delimiter characters */
    while (*fc && !isalpha((unsigned char)*fc) &&
           !isdigit((unsigned char)*fc) && !issign(*fc)) {
      fc++;
    }
    if (!*fc) break;

    /* Get element name (or charge sign) */
    i = 0;
    parse_sign = 0;
    while (*fc) {
      unsigned char uc = (unsigned char)*fc;
      if (!isalpha(uc) && !issign(*fc)) break;
      if (issign(*fc)) {
        if (i > 0) break; /* Finish preceding element before parsing charge */
        fc++;
        elem[0] = 'e';
        i = 1;
        if (uc == '+') parse_sign = -1; /* Electrons are SUBTRACTED */
        if (uc == '-') parse_sign = +1; /* Electrons are ADDED      */
        IFDEBUG printf("Molecule charge '%c' parsed (sign %d)\n", (char)uc, parse_sign);
        break;
      }
      /* Start a new element on an uppercase letter after at least one letter */
      if (i > 0 && isupper(uc)) break;
      if (i >= 2) break;

      fc++;
      if (i == 0) {
        elem[i++] = (char)toupper(uc);
      } else {
        elem[i++] = (char)tolower(uc);
      }
    }
    elem[i] = '\0';

    IFDEBUG printf("Element name is '%s'\n", elem);

    E = -1;
    for (e = 0; e < NI; e++) {
      if (strcmp(elem, Isotopes[e].name) == 0) {
        E = e;
        break;
      }
    }

    if (E == -1) {
      fprintf(stderr, "Unknown element: '%s'\n", elem);
    } else {
      IFDEBUG printf("Element ID is: %d\n", E);
    }

    /* Get element count */
    id = parse_uint(&fc, 1);
    if (parse_sign) id *= parse_sign;

    IFDEBUG printf("Element count is %d\n", id);

    if (E != -1) fi[E] += id;
  }
}

/**
 * @brief Parse an element-range specification string into a `limits_t` bounds array.
 *
 * Zero-initializes all `NI + 1` entries of `lim` before parsing. Each element entry
 * may specify:
 * - No count (e.g., `Cl`)      : fixed range `[1, 1]`.
 * - Single count (e.g., `N6`)  : fixed range `[6, 6]`.
 * - Range (e.g., `C8-10`)      : inclusive range `[8, 10]`. If bounds are reversed
 *                                (e.g., `C10-8`), they are automatically swapped so
 *                                `lim[E].min <= lim[E].max` always holds.
 *
 * Unknown element symbols are logged to `stderr` and skipped. Repeated entries for
 * the same element overwrite the earlier range for that element.
 *
 * @param[in]  fc  NUL-terminated limits string (e.g., `"C8-10,N6,O1-2,H1-13,Cl"`);
 *                 must be non-NULL.
 * @param[out] lim Destination bounds array; overwritten completely on entry.
 */
void Parse_Limits(const char *fc, limits_t lim)
{
  char elem[3];
  int i, e, id, E, tmpi;

  for (e = 0; e <= NI; e++) {
    lim[e].min = 0;
    lim[e].max = 0;
  }

  IFDEBUG printf("\nParsing limits string: '%s'\n\n", fc);
  while (*fc) {
    /* Skip delimiter characters */
    while (*fc && !isalpha((unsigned char)*fc) && !isdigit((unsigned char)*fc)) {
      fc++;
    }
    if (!*fc) break;

    /* Get element name */
    i = 0;
    while (*fc) {
      unsigned char uc = (unsigned char)*fc;
      if (!isalpha(uc)) break;
      if (i > 0 && isupper(uc)) break;
      if (i >= 2) break;

      fc++;
      if (i == 0) {
        elem[i++] = (char)toupper(uc);
      } else {
        elem[i++] = (char)tolower(uc);
      }
    }
    elem[i] = '\0';

    IFDEBUG printf("Element name is '%s'\n", elem);

    E = -1;
    for (e = 0; e < NI; e++) {
      if (strcmp(elem, Isotopes[e].name) == 0) {
        E = e;
        break;
      }
    }

    if (E == -1) {
      fprintf(stderr, "Unknown element: '%s'\n", elem);
    } else {
      IFDEBUG printf("Element ID is: %d\n", E);
    }

    /* Get lower limit of element count */
    id = parse_uint(&fc, 1);
    IFDEBUG printf("Element lower limit is %d\n", id);

    if (E != -1) lim[E].min = id;

    if (*fc == '-') {
      fc++;
      id = parse_uint(&fc, 1);
    }
    IFDEBUG printf("Element upper limit is %d\n", id);

    if (E != -1) {
      lim[E].max = id;
      if (lim[E].max < lim[E].min) {
        tmpi = lim[E].max;
        lim[E].max = lim[E].min;
        lim[E].min = tmpi;
        IFDEBUG printf("Upper and lower limits swapped\n");
      }
    }
  }
}


