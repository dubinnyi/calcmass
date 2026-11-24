# calcmass

A lightweight C11 command-line utility for high-resolution mass spectrometry (HRMS): exact monoisotopic mass and *m/z* calculation, Ring and Double Bond (RDB) unsaturation analysis, and molecular formula determination from accurate mass measurements.

> **Note:** The codebase was audited, bug-fixed, and documented with **Gemini 4 Neon (High)**.

## Features

- **Exact mass & *m/z* calculation** (`-mass` / `--mass`) for neutral molecules and 11 common ESI/MS ion adducts (including singly and doubly charged species).
- **Brute-force formula determination** (`-chem` / `--chem`) over user-specified elemental composition ranges, ranked by mass error in `mDa` and `ppm`.
- **RDB / unsaturation filtering** (`-chem-rdb` / `--chem-rdb`) to automatically reject chemically invalid formulas with negative or half-integer RDB values.
- **Short script-friendly output** (`-s` / `--s`) printing only the *m/z* line for `--mass` or the results table starting from `Formula` for `--chem` / `--chem-rdb`.
- **Flexible formula input**: accepts both comma-separated (`C10,H13,N6,O2,Cl`) and standard compact chemical notation (`C10H13N6O2Cl`, `C2H5OH`).
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

### Built-in Elements & Monoisotopic Masses

`C` (12.000000), `H` (1.007825), `N` (14.003074), `O` (15.994915), `P` (30.973763), `S` (31.972072), `Na` (22.989770), `Ca` (39.962591), `K` (38.963708), `B` (11.009305), `F` (18.998403), `Cl` (34.968853), `Br` (78.918336), and electron `e` (0.000549).

### Formula and Limits Syntax

- **Formula**: `C16H17N3O7`, `C16,H17,N3,O7`, or `C2H5OH` (repeated elements are summed).
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
Version 2026-11-24
Mass type is '[M-H]-' (4), charge= -1, madd=-1.007276

Parsed formula: 'C16,H17,N3,O7', RDB = 10.0

Formula:
 C 16 12.000000  192.000000
 H 17  1.007825   17.133025
 N  3 14.003074   42.009222
 O  7 15.994915  111.964405
------------------------------
        [M-H]-: (M)    362.099376
        [M-H]-: (M/Z)  362.099376
```

With `-s` / `--s` (short output):

```bash
./calcmass -s --mass MH- C16H17N3O7
```

```text
[M-H]-: (M/Z)  362.099376
```

---

### 2. Formula Search Without RDB Filtering (`--chem -s`)

Suppose an ESI+ spectrum gives a molecular ion `[M+H]+` at *m/z* `324.1067`, and we search across `C0-30,H0-60,N0-5,O0-20`:

```bash
./calcmass -s --chem MH+ C0-30,H0-60,N0-5,O0-20 324.1067
```

```text
                  Formula      Theory     Mass diff   Mass diff   RDB
               H23,N2,O17: 324.106954     0.254 mDa   0.784 ppm  -9.5
            C13,H15,N4,O6: 324.106437    -0.263 mDa  -0.811 ppm   8.5
             C15,H17,N,O7: 324.107780     1.080 mDa   3.332 ppm   8.0
              C12,H19,O10: 324.105101    -1.599 mDa  -4.934 ppm   3.5
            C16,H13,N5,O3: 324.109116     2.416 mDa   7.454 ppm  13.0
            C3,H21,N3,O14: 324.109633     2.933 mDa   9.049 ppm  -5.0
            C10,H17,N3,O9: 324.103758    -2.942 mDa  -9.077 ppm   4.0
            C18,H15,N2,O4: 324.110459     3.759 mDa  11.598 ppm  12.5
               C5,H23,O15: 324.110976     4.276 mDa  13.193 ppm  -5.5
             C22,H13,N,O2: 324.101905    -4.795 mDa -14.795 ppm  17.0
            C6,H19,N4,O11: 324.112312     5.612 mDa  17.315 ppm  -0.5
            C7,H19,N2,O12: 324.101079    -5.621 mDa -17.343 ppm  -0.5
             C20,H11,N4,O: 324.100562    -6.138 mDa -18.938 ppm  17.5
             C21,H13,N3,O: 324.113138     6.438 mDa  19.864 ppm  17.0
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
             C15,H17,N,O7: 324.107780     1.080 mDa   3.332 ppm   8.0
            C16,H13,N5,O3: 324.109116     2.416 mDa   7.454 ppm  13.0
            C10,H17,N3,O9: 324.103758    -2.942 mDa  -9.077 ppm   4.0
             C22,H13,N,O2: 324.101905    -4.795 mDa -14.795 ppm  17.0
             C21,H13,N3,O: 324.113138     6.438 mDa  19.864 ppm  17.0
```

All nonsensical hits are eliminated, leaving only 5 chemically plausible formulas within 20 ppm (`C15,H17,N,O7` at +3.33 ppm, `C16,H13,N5,O3` at +7.45 ppm, `C10,H17,N3,O9` at −9.08 ppm, etc.).

---

### 4. Multi-Stage MS(*n*) Fragmentation Interpretation (ESI−)

When a single precursor *m/z* still admits several valid formulas, combining MS2, neutral loss, and MS3 measurements resolves the ambiguity.

#### Step A — MS2 precursor ion (`[M-H]-` = `362.1003`)

```bash
./calcmass -s --chem-rdb MH- C10-30,H0-30,N0-8,O0-10 362.1003
```

```text
                  Formula      Theory     Mass diff   Mass diff   RDB
            C17,H13,N7,O3: 362.100712     0.412 mDa   1.138 ppm  15.0
            C16,H17,N3,O7: 362.099376    -0.924 mDa  -2.552 ppm  10.0
                C28,H13,N: 362.097523    -2.777 mDa  -7.669 ppm  23.0
             C21,H17,N,O5: 362.103398     3.098 mDa   8.556 ppm  14.0
             C22,H13,N5,O: 362.104734     4.434 mDa  12.245 ppm  19.0
            C11,H17,N5,O9: 362.095354    -4.946 mDa -13.659 ppm   6.0
            C10,H17,N7,O8: 362.106587     6.287 mDa  17.363 ppm   6.0
            C23,H13,N3,O2: 362.093501    -6.799 mDa -18.777 ppm  19.0
```

#### Step B — Neutral loss (`M` = `120.0427`, i.e., `362.1003 − 242.0576`)

Using adduct type `M` (neutral molecule) identifies candidate neutral losses:

```bash
./calcmass -s --chem-rdb M C0-10,H0-20,N0-6,O0-8 120.0427
```

```text
                  Formula      Theory     Mass diff   Mass diff   RDB
                 C4,H8,O4: 120.042260    -0.440 mDa  -3.665 ppm   1.0
                 C5,H4,N4: 120.043596     0.896 mDa   7.464 ppm   6.0
```

#### Step C — MS3 fragment ion (`[M-H]-` = `242.0576`)

```bash
./calcmass -s --chem-rdb MH- C0-20,H0-20,N0-6,O0-6 242.0576
```

```text
                  Formula      Theory     Mass diff   Mass diff   RDB
             C12,H9,N3,O3: 242.057116    -0.484 mDa  -2.000 ppm  10.0
               C17,H9,N,O: 242.061138     3.538 mDa  14.616 ppm  14.0
              C7,H9,N5,O5: 242.053094    -4.506 mDa -18.615 ppm   6.0
```

Combining the MS3 fragment (`C12H9N3O3`, −2.00 ppm) and neutral loss (`C4H8O4`, −3.67 ppm) gives **`C16H17N3O7`** (`C12H9N3O3 + C4H8O4 = C16H17N3O7`), uniquely confirming the precursor formula (consistent with the additional observed `H2O` loss to `344.0895` and `C5H10O4` loss to `228.0420`).

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
