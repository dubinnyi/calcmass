#!/bin/sh
# End-to-end CLI integration tests for calcmass.
# Usage: sh tests/test_cli.sh [path/to/calcmass]

set -u

BIN="${1:-./calcmass}"
PASS=0
FAIL=0

assert_ok_contains() {
  desc="$1"
  expected="$2"
  shift 2
  out="$("$BIN" "$@" 2>&1)"
  rc=$?
  if [ "$rc" -ne 0 ]; then
    printf "FAIL: %s (expected exit 0, got %d)\n%s\n" "$desc" "$rc" "$out" >&2
    FAIL=$((FAIL + 1))
    return
  fi
  case "$out" in
    *"$expected"*)
      PASS=$((PASS + 1))
      ;;
    *)
      printf "FAIL: %s (expected substring '%s')\nActual output:\n%s\n" "$desc" "$expected" "$out" >&2
      FAIL=$((FAIL + 1))
      ;;
  esac
}

assert_ok_not_contains() {
  desc="$1"
  forbidden="$2"
  shift 2
  out="$("$BIN" "$@" 2>&1)"
  rc=$?
  if [ "$rc" -ne 0 ]; then
    printf "FAIL: %s (expected exit 0, got %d)\n%s\n" "$desc" "$rc" "$out" >&2
    FAIL=$((FAIL + 1))
    return
  fi
  case "$out" in
    *"$forbidden"*)
      printf "FAIL: %s (did not expect substring '%s')\nActual output:\n%s\n" "$desc" "$forbidden" "$out" >&2
      FAIL=$((FAIL + 1))
      ;;
    *)
      PASS=$((PASS + 1))
      ;;
  esac
}

assert_err_contains() {
  desc="$1"
  expected="$2"
  shift 2
  out="$("$BIN" "$@" 2>&1)"
  rc=$?
  if [ "$rc" -eq 0 ]; then
    printf "FAIL: %s (expected non-zero exit, got 0)\n%s\n" "$desc" "$out" >&2
    FAIL=$((FAIL + 1))
    return
  fi
  case "$out" in
    *"$expected"*)
      PASS=$((PASS + 1))
      ;;
    *)
      printf "FAIL: %s (expected error '%s')\nActual output:\n%s\n" "$desc" "$expected" "$out" >&2
      FAIL=$((FAIL + 1))
      ;;
  esac
}

# --- 1. Exact masses and adducts (-mass / --mass) ---
assert_ok_contains "neutral H2O short"         "18.010565"  -s --mass M H2O
assert_ok_contains "MH+ H2O short"             "19.017841"  -s --mass MH+ H2O
assert_ok_contains "MH- C16H17N3O7 short"      "362.099373" -s --mass MH- C16H17N3O7
assert_ok_contains "M2H++ C16H17N3O7"          "182.560601" -s --mass M2H++ C16H17N3O7
assert_ok_contains "MH++ alias"                "182.560601" -s --mass MH++ C16H17N3O7
assert_ok_contains "MH2+ alias"                "182.560601" -s --mass MH2+ C16H17N3O7
assert_ok_contains "MH+2 alias"                "182.560601" -s --mass MH+2 C16H17N3O7
assert_ok_contains "MNa+ H2O"                  "40.999785"  -s --mass MNa+ H2O
assert_ok_contains "MK+ H2O"                   "56.973723"  -s --mass MK+ H2O
assert_ok_contains "MCa+ H2O"                  "56.964782"  -s --mass MCa+ H2O
assert_ok_contains "parentheses Ca(OH)2"       "73.968070"  -s --mass M "Ca(OH)2"
assert_ok_contains "parentheses (CH3)2SO"      "78.013936"  -s --mass M "(CH3)2SO"
assert_ok_contains "nested ((CH2)2O)3H2O"      "150.089209" -s --mass M "((CH2)2O)3H2O"

# --- 2. RDB range display for S and P ---
assert_ok_contains "cysteine RDB range -mass"  "RDB = 1.0..3.0"  --mass M C3H7NO2S
assert_ok_contains "SF6 RDB range -mass"       "RDB = -2.0..0.0" --mass M SF6
assert_ok_contains "H3PO4 RDB range -mass"     "RDB = 0.0..1.0"  --mass M H3PO4
assert_ok_contains "unambiguous RDB -mass"     "RDB = 10.0"      --mass MH- C16H17N3O7

# --- 3. RDB search modes & filters ---
assert_ok_contains "chem-rdb keeps SF6 (S=6)" \
  "-2.0..0.0" -s --chem-rdb M S0-1,F0-6 145.962490
assert_ok_contains "chem-rdb --rdb-strict rejects SF6" \
  "No matching chemical formulas found." -s --chem-rdb --rdb-strict M S1,F6 145.962490
assert_ok_contains "chem-rdb --rdb-strict keeps cysteine" \
  "C3,H7,N,O2,S" -s --chem-rdb --rdb-strict M C3,H7,N1,O2,S1 121.019750
assert_ok_contains "chem-rdb min/max keeps C15H17NO7" \
  "C15,H17,N,O7" -s --chem-rdb --rdb-min 8 --rdb-max 13 MH+ C0-30,H0-60,N0-5,O0-20 324.1067
assert_ok_not_contains "chem-rdb min/max excludes RDB 4" \
  "C10,H17,N3,O9" -s --chem-rdb --rdb-min 8 --rdb-max 13 MH+ C0-30,H0-60,N0-5,O0-20 324.1067
assert_ok_not_contains "chem-rdb min/max excludes RDB 17" \
  "C22,H13,N,O2" -s --chem-rdb --rdb-min 8 --rdb-max 13 MH+ C0-30,H0-60,N0-5,O0-20 324.1067
assert_ok_contains "chem with --rdb-min keeps half-integer hit" \
  "C13,H15,N4,O6" -s --chem --rdb-min 8 --rdb-max 9 MH+ C0-30,H0-60,N0-5,O0-20 324.1067
assert_ok_not_contains "chem with --rdb-min excludes negative RDB" \
  "H23,N2,O17" -s --chem --rdb-min 0 MH+ C0-30,H0-60,N0-5,O0-20 324.1067
assert_ok_contains "target mass ASCII subtraction" \
  "C4,H8,O4" -s --chem-rdb M C0-10,H0-20,N0-6,O0-8 "362.1003 - 242.0576"
assert_ok_contains "target mass Unicode minus subtraction" \
  "C4,H8,O4" -s --chem-rdb M C0-10,H0-20,N0-6,O0-8 "362.1003 − 242.0576"
assert_ok_contains "target mass addition" \
  "C16,H17,N3,O7" -s --chem-rdb MH- C10-30,H0-30,N0-8,O0-10 "242.0576 + 120.0427"

# --- 4. Negative / error handling tests ---
assert_err_contains "too few args"             "Too few arguments" --mass M
assert_err_contains "unknown adduct"           "Unknown mass type" --mass MLi+ H2O
assert_err_contains "unknown mode"             "Unknown argument"  --foo M H2O
assert_err_contains "unknown element formula"  "Unknown element"   --mass M C10Si
assert_err_contains "unclosed paren"           "Unclosed '('"      --mass M "Ca(OH2"
assert_err_contains "malformed limits C-10"    "Missing lower limit" --chem MH+ C-10 100.0
assert_err_contains "malformed limits C8-"     "Missing upper limit" --chem MH+ C8- 100.0
assert_err_contains "bad target mass"          "Invalid target mass" --chem MH+ C1-5 abc
assert_err_contains "negative target mass"     "Invalid target mass" --chem MH+ C1-5 -10
assert_err_contains "negative expr result"     "Invalid target mass" --chem M C0-10 "242.0576 - 362.1003"
assert_err_contains "trailing op in mass expr" "Invalid target mass" --chem M C0-10 "362.1003 -"
assert_err_contains "bad --rdb-min"            "Invalid --rdb-min"   --chem --rdb-min abc M C1-2 12.0
assert_err_contains "rdb-min > rdb-max"        "cannot exceed"       --chem --rdb-min 10 --rdb-max 5 M C1-2 12.0
assert_err_contains "extra arg in -mass"       "Too many arguments"  --mass M H2O extra
assert_err_contains "extra arg in -chem"       "Too many arguments"  --chem M C1-2 12.0 extra

printf "CLI tests: %d passed, %d failed\n" "$PASS" "$FAIL"
[ "$FAIL" -eq 0 ]
