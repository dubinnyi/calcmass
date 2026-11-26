# calcmass

A lightweight C11 command-line utility for high-resolution mass spectrometry (HRMS): exact monoisotopic mass and *m/z* calculation, Ring and Double Bond (RDB) unsaturation analysis, and molecular formula determination from accurate mass measurements.

> **Note:** The codebase was audited, bug-fixed, and documented with **Gemini 4 Neon (High)**.

## Features

- **Exact mass & *m/z* calculation** (`-mass` / `--mass`) for neutral molecules and 11 common ESI/MS ion adducts (including singly and doubly charged species).
- **Brute-force formula determination** (`-chem` / `--chem`) over user-specified elemental composition ranges, ranked by mass error in `mDa` and `ppm`.
- **RDB / unsaturation filtering** (`-chem-rdb` / `--chem-rdb`) to automatically reject chemically invalid formulas with negative or half-integer RDB values.
- **Short script-friendly output** (`-s` / `--s`) printing only the *m/z* line for `--mass` or the results table starting from `Formula` for `--chem` / `--chem-rdb`.
- **Flexible, strict formula input**: accepts comma-separated (`C10,H13,N6,O2,Cl`), compact standard notation (`C10H13N6O2Cl`, `C2H5OH`), and nested parenthesized groups (`Ca(OH)2`, `(CH3)2SO`), failing fast with non-zero exit codes on unknown elements, malformed ranges, or integer overflows.
- **Zero external dependencies**: builds with standard GCC and `libm`.

---

## Building

### Requirements

- A C11-compatible compiler (`gcc`)
- `make`

### Compile

```bash
git clone https://github.com/dubinnyi/calcmass
cd calcmass
make
```

This builds the `./calcmass` executable in the project root. To remove build artifacts:

```bash
make clean
```

---

## Usage

```bash
./calcmass [-s|--s] --mass     <type> <formula>
./calcmass [-s|--s] --chem     <type> <limits> <target_mass>
./calcmass [-s|--s] --chem-rdb <type> <limits> <target_mass>
```

- Short single-dash mode flags (`-mass`, `-chem`, `-chem-rdb`) are also accepted.
- `-s` / `--s` can be placed anywhere in the command line to suppress banners and metadata, leaving only the essential result.

### Supported Ion / Adduct Types (`<type>`)

| CLI argument | Species | Charge | Notes |
| :--- | :--- | :---: | :--- |
| `M` | `[M]` | `0` | Neutral molecule (useful for neutral losses) |
| `M+` | `[M]+` | `+1` | Radical cation (`-e`) |
| `M-` | `[M]-` | `-1` | Radical anion (`+e`) |
| `MH+` | `[M+H]+` | `+1` | Protonated molecule (ESI+) |
| `MH-` | `[M-H]-` | `-1` | Deprotonated molecule (ESI−) |
| `MNa+` | `[M+Na]+` | `+1` | Sodium adduct |
| `MK+` | `[M+K]+` | `+1` | Potassium adduct |
| `MCa+` | `[M+Ca-H]+` | `+1` | Calcium adduct |
| `M2H++` | `[M+2H]++` | `+2` | Doubly protonated (aliases: `MH++`, `MH2+`, `MH+2`) |
| `MHNa++` | `[M+H+Na]++` | `+2` | Doubly charged H/Na adduct |
| `M2Na++` | `[M+2Na]++` | `+2` | Doubly sodiated adduct |

### Built-in Elements & Monoisotopic Masses (AME2020 / CODATA 2022)

`C` (`12.000000000000`), `H` (`1.007825031898`), `N` (`14.00307400425`), `O` (`15.99491461926`), `P` (`30.97376199768`), `S` (`31.9720711744`), `Na` (`22.98976928195`), `Ca` (`39.962590863`), `K` (`38.963706486`), `B` (`11.009305166`), `F` (`18.99840316207`), `Cl` (`34.968852694`), `Br` (`78.918337574`), and electron `e` (`0.000548579909`).

### Formula and Limits Syntax

- **Formula**: `C16H17N3O7`, `C16,H17,N3,O7`, `C2H5OH`, `Ca(OH)2`, or `(CH3)2SO` (repeated elements are summed; nested parentheses are supported).
- **Limits**: `<Elem><min>-<max>` or fixed count `<Elem><n>` / `<Elem>` (implies `1`), e.g. `C0-30,H0-60,N0-5,O0-20` or `C8-10,N6,O1-2,H1-13,Cl`.
- **Search window**: prints up to the top 30 hits within `0.1 Da` and `20 ppm` of `<target_mass>` (or the single closest match if none fall within the window).

---

## Examples

### 1. Exact Mass Calculation for a Known Formula (`--mass`)

Calculate the theoretical mass, *m/z*, elemental breakdown, and RDB of `C16H17N3O7` in negative mode (`MH-`):

```bash
./calcmass --mass MH- C16H17N3O7
```

Full output:

```text
calcmass by Maxim A. Dubinnyi
Version 2026-11-26
Mass type is '[M-H]-' (4), charge= -1, madd=-1.007276

Parsed formula: 'C16,H17,N3,O7', RDB = 10.0

Formula:
 C 16 12.000000  192.000000
 H 17  1.007825   17.133026
 N  3 14.003074   42.009222
 O  7 15.994915  111.964402
------------------------------
        [M-H]-: (M)    362.099373
        [M-H]-: (M/Z)  362.099373
```

With `-s` / `--s` (short output):

```bash
./calcmass -s --mass MH- C16H17N3O7
```

```text
[M-H]-: (M/Z)  362.099373
```

---

### 2. Formula Search Without RDB Filtering (`--chem -s`)

Suppose an ESI+ spectrum gives a molecular ion `[M+H]+` at *m/z* `324.1067`, and we search across `C0-30,H0-60,N0-5,O0-20`:

```bash
./calcmass -s --chem MH+ C0-30,H0-60,N0-5,O0-20 324.1067
```

```text
                  Formula      Theory     Mass diff   Mass diff   RDB
               H23,N2,O17: 324.106949     0.249 mDa   0.767 ppm  -9.5
            C13,H15,N4,O6: 324.106436    -0.264 mDa  -0.816 ppm   8.5
             C15,H17,N,O7: 324.107778     1.078 mDa   3.327 ppm   8.0
              C12,H19,O10: 324.105098    -1.602 mDa  -4.942 ppm   3.5
            C16,H13,N5,O3: 324.109116     2.416 mDa   7.454 ppm  13.0
            C3,H21,N3,O14: 324.109629     2.929 mDa   9.037 ppm  -5.0
            C10,H17,N3,O9: 324.103756    -2.944 mDa  -9.085 ppm   4.0
            C18,H15,N2,O4: 324.110458     3.758 mDa  11.596 ppm  12.5
               C5,H23,O15: 324.110971     4.271 mDa  13.179 ppm  -5.5
             C22,H13,N,O2: 324.101905    -4.795 mDa -14.794 ppm  17.0
            C6,H19,N4,O11: 324.112309     5.609 mDa  17.306 ppm  -0.5
            C7,H19,N2,O12: 324.101075    -5.625 mDa -17.354 ppm  -0.5
             C20,H11,N4,O: 324.100562    -6.138 mDa -18.937 ppm  17.5
             C21,H13,N3,O: 324.113138     6.438 mDa  19.865 ppm  17.0
```

Notice that the closest raw mass matches include chemically impossible compositions such as `H23,N2,O17` (RDB `-9.5`) and radical/odd-valence formulas like `C13,H15,N4,O6` (RDB `8.5`).

---

### 3. Filtering by Unsaturation / RDB (`--chem-rdb -s`)

Running the same query with `--chem-rdb` discards negative and half-integer RDB entries, keeping only valid neutral closed-shell molecular formulas:

```bash
./calcmass -s --chem-rdb MH+ C0-30,H0-60,N0-5,O0-20 324.1067
```

```text
                  Formula      Theory     Mass diff   Mass diff   RDB
             C15,H17,N,O7: 324.107778     1.078 mDa   3.327 ppm   8.0
            C16,H13,N5,O3: 324.109116     2.416 mDa   7.454 ppm  13.0
            C10,H17,N3,O9: 324.103756    -2.944 mDa  -9.085 ppm   4.0
             C22,H13,N,O2: 324.101905    -4.795 mDa -14.794 ppm  17.0
             C21,H13,N3,O: 324.113138     6.438 mDa  19.865 ppm  17.0
```

All nonsensical hits are eliminated, leaving only 5 chemically plausible formulas within 20 ppm (`C15,H17,N,O7` at +3.33 ppm, `C16,H13,N5,O3` at +7.45 ppm, `C10,H17,N3,O9` at −9.09 ppm, etc.).

---

### 4. Multi-Stage MS(*n*) Fragmentation Interpretation (ESI−)

When a single precursor *m/z* still admits several valid formulas, combining MS2, neutral loss, and MS3 measurements resolves the ambiguity.

#### Step A — MS2 precursor ion (`[M-H]-` = `362.1003`)

```bash
./calcmass -s --chem-rdb MH- C10-30,H0-30,N0-8,O0-10 362.1003
```

```text
                  Formula      Theory     Mass diff   Mass diff   RDB
            C17,H13,N7,O3: 362.100711     0.411 mDa   1.135 ppm  15.0
            C16,H17,N3,O7: 362.099373    -0.927 mDa  -2.559 ppm  10.0
                C28,H13,N: 362.097523    -2.777 mDa  -7.669 ppm  23.0
             C21,H17,N,O5: 362.103396     3.096 mDa   8.551 ppm  14.0
             C22,H13,N5,O: 362.104734     4.434 mDa  12.244 ppm  19.0
            C11,H17,N5,O9: 362.095351    -4.949 mDa -13.668 ppm   6.0
            C10,H17,N7,O8: 362.106584     6.284 mDa  17.355 ppm   6.0
            C23,H13,N3,O2: 362.093500    -6.800 mDa -18.779 ppm  19.0
```

#### Step B — Neutral loss (`M` = `120.0427`, i.e., `362.1003 − 242.0576`)

Using adduct type `M` (neutral molecule) identifies candidate neutral losses:

```bash
./calcmass -s --chem-rdb M C0-10,H0-20,N0-6,O0-8 120.0427
```

```text
                  Formula      Theory     Mass diff   Mass diff   RDB
                 C4,H8,O4: 120.042259    -0.441 mDa  -3.676 ppm   1.0
                 C5,H4,N4: 120.043596     0.896 mDa   7.465 ppm   6.0
```

#### Step C — MS3 fragment ion (`[M-H]-` = `242.0576`)

```bash
./calcmass -s --chem-rdb MH- C0-20,H0-20,N0-6,O0-6 242.0576
```

```text
                  Formula      Theory     Mass diff   Mass diff   RDB
             C12,H9,N3,O3: 242.057115    -0.485 mDa  -2.005 ppm  10.0
               C17,H9,N,O: 242.061137     3.537 mDa  14.614 ppm  14.0
              C7,H9,N5,O5: 242.053092    -4.508 mDa -18.624 ppm   6.0
```

Combining the MS3 fragment (`C12H9N3O3`, −2.01 ppm) and neutral loss (`C4H8O4`, −3.68 ppm) gives **`C16H17N3O7`** (`C12H9N3O3 + C4H8O4 = C16H17N3O7`), uniquely confirming the precursor formula (consistent with the additional observed `H2O` loss to `344.0895` and `C5H10O4` loss to `228.0420`).

---

## Project Structure

- `calcmass.c` — CLI entry point, candidate ranking, and formatted output.
- `isotopes.h` / `isotopes.c` — Element table, adduct definitions, mass/*m/z*/RDB calculators, and sort comparator.
- `parse.c` — Chemical formula and elemental range parsers.
- `multiloop.c` — Multi-dimensional iterator over elemental search limits.
- `Makefile` — Build configuration (`gcc`, C11, strict warnings).
- `IsotopeMassTable.txt`, `Isotope_Table_NIST.txt` — Reference isotope mass tables.

## Author

Dr. Maxim A. Dubinnyi

ORCID https://orcid.org/0000-0002-5463-0169

Shemyakin-Ovchinnikov Institute of bioorganic chemistry, 
Moscow, Russia (http://www.ibch.ru/en)
