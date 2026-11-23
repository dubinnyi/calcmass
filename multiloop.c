/**
 * @file multiloop.c
 * @brief Mixed-radix multi-dimensional iterator over chemical element bounds.
 *
 * Implements an odometer-style counter across all elements `1 .. NI-1` bounded
 * by a `limits_t` structure. Intended to drive a standard `for` loop:
 * @code
 *   for (Loop_Init(&fi, &lim); Loop_End(&fi, &lim); Loop_Next(&fi, &lim)) { ... }
 * @endcode
 *
 * Slot `(*ft)[NI]` serves as an overflow carry sentinel: it starts at 0 and
 * increments to 1 once all combinations in the Cartesian product have been visited.
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "isotopes.h"

#define DEBUG 0
#define IFDEBUG if(DEBUG)

/**
 * @brief Copy all `NI + 1` slots of a `formula_t` vector from `src` to `dst`.
 *
 * @param[out] dst Destination formula pointer; must be non-NULL.
 * @param[in]  src Source formula pointer; must be non-NULL.
 */
void Loop_Copy(formula_t *dst, formula_t *src)
{
  int e;
  for (e = 0; e < NI + 1; e++)
    (*dst)[e] = (*src)[e];
}

/**
 * @brief Initialize iterator state `ft` to the first combination in `lim`.
 *
 * Sets electron count `(*ft)[0] = 0`, each element `e` in `[1, NI)` to
 * `(*lim)[e].min`, and the completion sentinel `(*ft)[NI] = 0`.
 *
 * @param[out] ft  Iterator state to initialize; must be non-NULL.
 * @param[in]  lim Element bounds; must be non-NULL with `(*lim)[e].min <= (*lim)[e].max`
 *                 for every `e` in `[1, NI)`.
 */
void Loop_Init(formula_t *ft, limits_t *lim)
{
  int e;
  (*ft)[0] = 0;
  IFDEBUG printf("Loop_Init: ");
  for (e = 1; e < NI; e++) {
    (*ft)[e] = (*lim)[e].min;
    IFDEBUG printf("%2d ", (*ft)[e]);
  }
  IFDEBUG printf("\n");
  (*ft)[NI] = 0;
}

/**
 * @brief Advance iterator state `ft` to the next combination in lexicographical order.
 *
 * Increments element 1 (`C`) and propagates carries through `1 .. NI-1` whenever
 * `(*ft)[e] > (*lim)[e].max`, resetting rolled-over positions to `(*lim)[e].min`.
 * When the last element (`NI - 1`) overflows, the sentinel `(*ft)[NI]` becomes non-zero.
 *
 * @param[in,out] ft  Iterator state previously initialized with `Loop_Init()`; must be non-NULL.
 * @param[in]     lim Element bounds used for roll-over checks; must be non-NULL and unchanged
 *                    since `Loop_Init()`, with `(*lim)[e].min <= (*lim)[e].max`.
 */
void Loop_Next(formula_t *ft, limits_t *lim)
{
  int e;
  (*ft)[1]++;
  for (e = 1; e < NI && (*ft)[e] > (*lim)[e].max; e++) {
    (*ft)[e] = (*lim)[e].min;
    (*ft)[e + 1]++;
  }
  IFDEBUG {
    printf("Loop_Next: ");
    for (e = 1; e < NI; e++) printf("%2d ", (*ft)[e]);
    printf("\n");
  }
}

/**
 * @brief Check whether the iterator `ft` still points to a valid (unexhausted) combination.
 *
 * Despite its name, this function acts as a loop-continuation predicate (returns
 * truthy while the loop should keep running, and 0 once exhausted).
 *
 * @param[in] ft  Iterator state initialized with `Loop_Init()`; must be non-NULL.
 * @param[in] lim Unused parameter kept for uniform signature with `Loop_Init`/`Loop_Next`.
 * @return `1` while combinations remain (`(*ft)[NI] == 0`); `0` once enumeration is complete.
 */
int Loop_End(formula_t *ft, limits_t *lim)
{
  int flag = ((*ft)[NI] == 0);
  (void)lim;
  IFDEBUG {
    int e;
    printf("Loop_End : ");
    for (e = 1; e < NI; e++) printf("%2d ", (*ft)[e]);
    printf(" flag=%d\n", flag);
  }
  return flag;
}
