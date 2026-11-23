#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "isotopes.h"

#define DEBUG 0
#define IFDEBUG if(DEBUG)

static int issign(char c)
{
  return (c == '+' || c == '-');
}

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


