/**
 * @file parse.c
 * @brief Parsers for chemical formula strings and element-range search limits.
 *
 * Supports comma/whitespace-separated notation (e.g., "C10,N6,O2,H13,Cl", "C8-10,H1-13,Cl"),
 * compact standard chemical notation (e.g., "C10H13N6O2Cl", "C2H5OH"), and nested
 * parenthesized groups in formulas (e.g., "Ca(OH)2", "(CH3)2SO").
 *
 * All syntax errors, unknown element symbols, malformed ranges, and integer overflows
 * are reported to `stderr` and cause the parser to return `-1` rather than silently
 * producing incorrect compositions or bounds.
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <math.h>

#include "isotopes.h"

#define DEBUG 0
#define IFDEBUG if(DEBUG)

/** @brief Maximum allowed nesting depth for parenthesized formula groups. */
#define MAX_PAREN_DEPTH 32

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
 * @brief Test whether a character is an allowed token separator.
 *
 * @param[in] c Character to inspect.
 * @return Non-zero if `c` is a comma or ASCII space/tab; 0 otherwise.
 */
static int isdelim(char c)
{
  return (c == ',' || c == ' ' || c == '\t');
}

/**
 * @brief Parse a non-negative decimal integer with overflow checking.
 *
 * Advances `*pfc` past all consecutive ASCII digit characters. If `*pfc` does not
 * start with a digit, sets `*out = default_val` and returns `0`. If the parsed
 * value exceeds `INT_MAX`, returns `-1`.
 *
 * @param[in,out] pfc         Address of the current cursor in the input string.
 * @param[in]     default_val Value stored in `*out` when no digits are present.
 * @param[out]    out         Destination for the parsed integer.
 * @return `0` on success; `-1` on integer overflow.
 */
static int parse_uint_checked(const char **pfc, int default_val, int *out)
{
  const char *fc = *pfc;
  long long acc = 0;
  int has_digits = 0;

  while (*fc && isdigit((unsigned char)*fc)) {
    has_digits = 1;
    acc = acc * 10 + (*fc - '0');
    if (acc > INT_MAX) {
      /* Advance past remaining digits to keep cursor consistent */
      while (*fc && isdigit((unsigned char)*fc)) fc++;
      *pfc = fc;
      return -1;
    }
    fc++;
  }

  *pfc = fc;
  *out = has_digits ? (int)acc : default_val;
  return 0;
}

/**
 * @brief Parse an element symbol starting at `*pfc` and resolve its index in `Isotopes[]`.
 *
 * Consumes an initial alphabetic character (normalized to uppercase) followed by
 * any immediately following lowercase alphabetic characters (so multi-letter invalid
 * tokens like `"Foo"` are reported intact rather than split).
 *
 * @param[in,out] pfc       Address of the current cursor (must point at an alpha char).
 * @param[out]    elem      Buffer receiving the normalized NUL-terminated symbol.
 * @param[in]     elem_size Capacity of `elem` in bytes (must be >= 4).
 * @return Element index in `[1, NI)` if found in `Isotopes[]`, or `-1` if unknown.
 */
static int parse_element(const char **pfc, char *elem, size_t elem_size)
{
  const char *fc = *pfc;
  size_t len = 0;

  if (*fc && isalpha((unsigned char)*fc)) {
    if (len + 1 < elem_size) {
      elem[len++] = (char)toupper((unsigned char)*fc);
    }
    fc++;
    while (*fc && islower((unsigned char)*fc)) {
      if (len + 1 < elem_size) {
        elem[len++] = *fc;
      }
      fc++;
    }
  }
  elem[len] = '\0';
  *pfc = fc;

  IFDEBUG printf("Element name is '%s'\n", elem);

  for (int e = 1; e < NI; e++) {
    if (strcmp(elem, Isotopes[e].name) == 0) {
      IFDEBUG printf("Element ID is: %d\n", e);
      return e;
    }
  }
  return -1;
}

/**
 * @brief Safely add a signed delta to an `int` accumulator with overflow detection.
 *
 * @param[in,out] target Accumulator to update.
 * @param[in]     delta  Signed value to add.
 * @return `0` on success; `-1` if `*target + delta` would overflow `int`.
 */
static int safe_add_int(int *target, long long delta)
{
  long long sum = (long long)(*target) + delta;
  if (sum > INT_MAX || sum < INT_MIN) {
    return -1;
  }
  *target = (int)sum;
  return 0;
}

/**
 * @brief Recursive helper for `Parse_Formula()` supporting nested parenthesized groups.
 *
 * @param[in,out] pfc        Address of the current cursor in the formula string.
 * @param[in,out] acc        Composition vector accumulating counts for the current group.
 * @param[in]     depth      Current parenthesis nesting depth (`0` at top level).
 * @param[out]    parsed_any Set to `1` if at least one valid token was parsed in this group.
 * @return `0` on success; `-1` on any syntax error, unknown element, or overflow.
 */
static int parse_formula_group(const char **pfc, formula_t acc, int depth, int *parsed_any)
{
  *parsed_any = 0;

  if (depth > MAX_PAREN_DEPTH) {
    fprintf(stderr, "Parentheses nested too deeply in formula (max %d)\n", MAX_PAREN_DEPTH);
    return -1;
  }

  while (**pfc) {
    while (**pfc && isdelim(**pfc)) {
      (*pfc)++;
    }
    if (!**pfc) {
      if (depth > 0) {
        fprintf(stderr, "Unclosed '(' in formula\n");
        return -1;
      }
      break;
    }

    char c = **pfc;

    if (c == ')') {
      if (depth == 0) {
        fprintf(stderr, "Unmatched ')' in formula\n");
        return -1;
      }
      (*pfc)++;
      return 0;
    }

    if (c == '(') {
      (*pfc)++;
      formula_t sub;
      for (int e = 0; e <= NI; e++) sub[e] = 0;
      int sub_any = 0;

      if (parse_formula_group(pfc, sub, depth + 1, &sub_any) != 0) {
        return -1;
      }
      if (!sub_any) {
        fprintf(stderr, "Empty parentheses '()' in formula\n");
        return -1;
      }

      int mult = 1;
      if (parse_uint_checked(pfc, 1, &mult) != 0) {
        fprintf(stderr, "Integer overflow in group multiplier in formula\n");
        return -1;
      }

      for (int e = 0; e < NI; e++) {
        if (sub[e] != 0) {
          long long scaled = (long long)sub[e] * (long long)mult;
          if (scaled > INT_MAX || scaled < INT_MIN || safe_add_int(&acc[e], scaled) != 0) {
            fprintf(stderr, "Integer overflow for element '%s' in formula\n", Isotopes[e].name);
            return -1;
          }
        }
      }
      *parsed_any = 1;
      continue;
    }

    if (issign(c)) {
      int sign = (c == '+') ? -1 : +1; /* '+' removes electrons, '-' adds electrons */
      (*pfc)++;

      int cnt = 1;
      if (parse_uint_checked(pfc, 1, &cnt) != 0) {
        fprintf(stderr, "Integer overflow in charge specification in formula\n");
        return -1;
      }

      IFDEBUG printf("Molecule charge '%c' parsed (delta e = %lld)\n",
                     c, (long long)sign * (long long)cnt);
      if (safe_add_int(&acc[0], (long long)sign * (long long)cnt) != 0) {
        fprintf(stderr, "Integer overflow in net charge in formula\n");
        return -1;
      }
      *parsed_any = 1;
      continue;
    }

    if (isalpha((unsigned char)c)) {
      char elem[16];
      int E = parse_element(pfc, elem, sizeof(elem));
      if (E == -1) {
        fprintf(stderr, "Unknown element: '%s'\n", elem);
        return -1;
      }

      int cnt = 1;
      if (parse_uint_checked(pfc, 1, &cnt) != 0) {
        fprintf(stderr, "Integer overflow for element '%s' in formula\n", elem);
        return -1;
      }

      IFDEBUG printf("Element count is %d\n", cnt);
      if (safe_add_int(&acc[E], (long long)cnt) != 0) {
        fprintf(stderr, "Integer overflow for element '%s' in formula\n", elem);
        return -1;
      }
      *parsed_any = 1;
      continue;
    }

    fprintf(stderr, "Unexpected character '%c' in formula\n", c);
    return -1;
  }

  if (depth > 0) {
    fprintf(stderr, "Unclosed '(' in formula\n");
    return -1;
  }
  return 0;
}

/**
 * @brief Parse a chemical formula string into a `formula_t` composition vector.
 *
 * Zero-initializes all `NI + 1` entries of `fi` before parsing. Recognizes:
 * - 1- or 2-letter element symbols in `Isotopes[]` (e.g., `C`, `H`, `Cl`, `Na`),
 *   separated by commas/spaces or uppercase boundaries.
 * - Optional non-negative integer counts (omitted count defaults to 1).
 * - Nested parenthesized groups with optional multipliers (e.g., `Ca(OH)2`, `(CH3)2SO`).
 * - Repeated occurrences of the same element, whose counts are summed (`C2H5OH` -> `C2,H6,O`).
 * - Charge suffixes `'+'` and `'-'` with optional multiplicity (e.g., `+`, `+2`, `++`, `-`),
 *   mapped to electron count `fi[0]` (`+` subtracts electrons, `-` adds electrons).
 *
 * @param[in]  fc NUL-terminated formula string to parse; must be non-NULL.
 * @param[out] fi Destination composition vector; overwritten completely on entry.
 * @return `0` on success; `-1` on empty input, unknown element, syntax error, or overflow.
 */
int Parse_Formula(const char *fc, formula_t fi)
{
  for (int e = 0; e <= NI; e++) fi[e] = 0;

  if (fc == NULL) {
    fprintf(stderr, "Empty chemical formula\n");
    return -1;
  }

  IFDEBUG printf("\nParsing formula string: '%s'\n\n", fc);

  const char *cursor = fc;
  int parsed_any = 0;
  if (parse_formula_group(&cursor, fi, 0, &parsed_any) != 0) {
    for (int e = 0; e <= NI; e++) fi[e] = 0;
    return -1;
  }
  if (!parsed_any) {
    fprintf(stderr, "Empty chemical formula\n");
    for (int e = 0; e <= NI; e++) fi[e] = 0;
    return -1;
  }

  int has_atoms = 0;
  for (int e = 1; e < NI; e++) {
    if (fi[e] > 0) {
      has_atoms = 1;
      break;
    }
  }
  if (!has_atoms) {
    fprintf(stderr, "Chemical formula must contain at least one atom\n");
    for (int e = 0; e <= NI; e++) fi[e] = 0;
    return -1;
  }
  return 0;
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
 * Missing bounds (e.g., `C-10`, `C8-`), extra hyphens (`C8-10-12`), unexpected
 * characters, unknown element symbols, and integer overflows are reported to
 * `stderr` and cause an error return. Repeated entries for the same element
 * overwrite the earlier range for that element.
 *
 * @param[in]  fc  NUL-terminated limits string (e.g., `"C8-10,N6,O1-2,H1-13,Cl"`).
 * @param[out] lim Destination bounds array; overwritten completely on entry.
 * @return `0` on success; `-1` on empty input, unknown element, malformed range, or overflow.
 */
int Parse_Limits(const char *fc, limits_t lim)
{
  for (int e = 0; e <= NI; e++) {
    lim[e].min = 0;
    lim[e].max = 0;
  }

  if (fc == NULL) {
    fprintf(stderr, "Empty element limits\n");
    return -1;
  }

  IFDEBUG printf("\nParsing limits string: '%s'\n\n", fc);

  int parsed_any = 0;
  while (*fc) {
    while (*fc && isdelim(*fc)) {
      fc++;
    }
    if (!*fc) break;

    if (!isalpha((unsigned char)*fc)) {
      fprintf(stderr, "Unexpected character '%c' in limits\n", *fc);
      goto fail;
    }

    char elem[16];
    int E = parse_element(&fc, elem, sizeof(elem));
    if (E == -1) {
      fprintf(stderr, "Unknown element: '%s'\n", elem);
      goto fail;
    }

    if (*fc == '-') {
      fprintf(stderr, "Missing lower limit before '-' for element '%s'\n", elem);
      goto fail;
    }

    int low = 1;
    if (parse_uint_checked(&fc, 1, &low) != 0) {
      fprintf(stderr, "Integer overflow in lower limit for element '%s'\n", elem);
      goto fail;
    }
    IFDEBUG printf("Element lower limit is %d\n", low);

    int high = low;
    if (*fc == '-') {
      fc++;
      if (!isdigit((unsigned char)*fc)) {
        fprintf(stderr, "Missing upper limit after '-' for element '%s'\n", elem);
        goto fail;
      }
      if (parse_uint_checked(&fc, 0, &high) != 0) {
        fprintf(stderr, "Integer overflow in upper limit for element '%s'\n", elem);
        goto fail;
      }
      if (*fc == '-') {
        fprintf(stderr, "Unexpected '-' in limits for element '%s'\n", elem);
        goto fail;
      }
    }
    IFDEBUG printf("Element upper limit is %d\n", high);

    if (high < low) {
      int tmp = high;
      high = low;
      low = tmp;
      IFDEBUG printf("Upper and lower limits swapped\n");
    }

    lim[E].min = low;
    lim[E].max = high;
    parsed_any = 1;
  }

  if (!parsed_any) {
    fprintf(stderr, "Empty element limits\n");
    goto fail;
  }
  return 0;

fail:
  for (int e = 0; e <= NI; e++) {
    lim[e].min = 0;
    lim[e].max = 0;
  }
  return -1;
}

/**
 * @brief Advance `*pp` past ASCII spaces/tabs and UTF-8 non-breaking spaces (U+00A0).
 */
static void skip_expr_ws(const char **pp)
{
  const char *p = *pp;
  while (*p) {
    if (*p == ' ' || *p == '\t') {
      p++;
    } else if ((unsigned char)p[0] == 0xC2 && (unsigned char)p[1] == 0xA0) {
      p += 2;
    } else {
      break;
    }
  }
  *pp = p;
}

/**
 * @brief Consume a `+` or `-` operator (including Unicode minus U+2212 and en-dash U+2013).
 *
 * @param[in,out] pp Address of the current cursor in the expression string.
 * @return `+1` for `+`, `-1` for `-`/`−`/`–`, or `0` if no operator is present.
 */
static int consume_mass_op(const char **pp)
{
  const char *p = *pp;
  if (*p == '+') {
    *pp = p + 1;
    return +1;
  }
  if (*p == '-') {
    *pp = p + 1;
    return -1;
  }
  /* UTF-8 U+2212 MINUS SIGN (E2 88 92) and U+2013 EN DASH (E2 80 93) */
  if ((unsigned char)p[0] == 0xE2 &&
      (((unsigned char)p[1] == 0x88 && (unsigned char)p[2] == 0x92) ||
       ((unsigned char)p[1] == 0x80 && (unsigned char)p[2] == 0x93))) {
    *pp = p + 3;
    return -1;
  }
  return 0;
}

/**
 * @brief Parse a single unsigned finite decimal floating-point term at `*pp`.
 *
 * Rejects leading signs, hex floats (`0x...`), `inf`, and `nan`.
 *
 * @param[in,out] pp  Address of the current cursor; advanced past the parsed number on success.
 * @param[out]    out Destination for the parsed value.
 * @return `0` on success; `-1` if no valid unsigned decimal number starts at `*pp`.
 */
static int parse_unsigned_term(const char **pp, double *out)
{
  skip_expr_ws(pp);
  const char *p = *pp;

  if (!isdigit((unsigned char)*p) &&
      !(*p == '.' && isdigit((unsigned char)p[1]))) {
    return -1;
  }
  if (p[0] == '0' && (p[1] == 'x' || p[1] == 'X')) {
    return -1;
  }

  char *endptr = NULL;
  errno = 0;
  double val = strtod(p, &endptr);
  if (errno != 0 || endptr == p || !isfinite(val)) {
    return -1;
  }

  *pp = endptr;
  *out = val;
  return 0;
}

/**
 * @brief Parse a numeric target mass or additive `+`/`-` arithmetic expression.
 *
 * Evaluates a left-to-right sum/difference of finite decimal numbers, supporting
 * ASCII `+`/`-` as well as UTF-8 Unicode minus (`−`, U+2212) and en-dash (`–`, U+2013)
 * with optional whitespace (e.g., `"362.1003 - 242.0576"`, `"362.1003 − 242.0576"`,
 * `"242.0576 + 120.0427"`). Useful for direct neutral-loss mass searches in MS/MS.
 *
 * @param[in]  expr NUL-terminated expression string; must be non-NULL.
 * @param[out] out  Destination for the evaluated result; must be non-NULL.
 * @return `0` on success; `-1` on empty input, syntax error, or non-finite result.
 */
int Parse_Mass_Expr(const char *expr, double *out)
{
  if (expr == NULL || out == NULL) {
    return -1;
  }

  const char *p = expr;
  skip_expr_ws(&p);
  if (!*p) {
    return -1;
  }

  int sign = consume_mass_op(&p);
  if (sign == 0) sign = +1;

  double term = 0.0;
  if (parse_unsigned_term(&p, &term) != 0) {
    return -1;
  }

  double total = sign * term;
  for (;;) {
    skip_expr_ws(&p);
    if (!*p) break;

    int op = consume_mass_op(&p);
    if (op == 0) {
      return -1;
    }
    if (parse_unsigned_term(&p, &term) != 0) {
      return -1;
    }
    total += op * term;
    if (!isfinite(total)) {
      return -1;
    }
  }

  *out = total;
  return 0;
}
