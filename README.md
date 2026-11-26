# calcmass

A lightweight C11 command-line utility for high-resolution mass spectrometry (HRMS): exact monoisotopic mass and *m/z* calculation, Ring and Double Bond (RDB) unsaturation analysis, and molecular formula determination from accurate mass measurements.

> **Note:** The codebase was audited, bug-fixed, and documented with **Gemini 4 Neon (High)**.

## Features

- **Exact mass & *m/z* calculation** (`-mass` / `--mass`) for neutral molecules and 11 common ESI/MS ion adducts (including singly and doubly charged species).
- **Brute-force formula determination** (`-chem` / `--chem`) over user-specified elemental composition ranges, ranked by mass error in `mDa` and `ppm`.
- **Neutral-molecule RDB / unsaturation filtering** (`-chem-rdb` / `--chem-rdb`, `--rdb-min`, `--rdb-max`, `--rdb-strict`) computing exact `[RDB_min, RDB_max]` intervals for variable-valence elements (`S`: 2..6, `P`: 3..5) and rejecting negative or half-integer neutral-molecule RDB values.
- **Short script-friendly output** (`-s` / `--s`) printing only the *m/z* line for `--mass` or the results table starting from `Formula` for `--chem` / `--chem-rdb`.
- **Flexible, strict formula input**: accepts comma-separated (`C10,H13,N6,O2,Cl`), compact standard notation (`C10H13N6O2Cl`, `C2H5OH`), and nested parenthesized groups (`'Ca(OH)2'`, `'(CH3)2SO'`), failing fast with non-zero exit codes on unknown elements, malformed ranges, or integer overflows.
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

This builds the `./calcmass` executable in the project root. Useful `make` targets:

```bash
make test       # Build and run C unit tests + CLI integration tests
make debug      # Build with -O0 -g and ASan/UBSan sanitizers (requires libasan/libubsan)
make install    # Install to /usr/local/bin (override with PREFIX=...)
make clean      # Remove executables and object files
```

---

## Usage

```bash
./calcmass [options] --mass     <type> <formula>
./calcmass [options] --chem     <type> <limits> <target_mass>
./calcmass [options] --chem-rdb <type> <limits> <target_mass>
```

- Short single-dash mode flags (`-mass`, `-chem`, `-chem-rdb`) are also accepted.
- Options can be placed anywhere on the command line:
  - `-h`, `--help` — show usage help and exit with status `0`.
  - `-s`, `--s` — suppress banners and metadata, printing only the *m/z* line (`--mass`) or the results table (`--chem`).
  - `--ppm <val>` — maximum relative mass error in ppm for `--chem` / `--chem-rdb` (default: `20.0`).
  - `--da <val>` — maximum absolute mass error in Da for `--chem` / `--chem-rdb` (default: `0.1`).
  - `--top <n>` — maximum number of top-ranked formulas to display (default: `30`).
  - `--rdb-min <val>`, `--rdb-max <val>` — inclusive lower/upper bounds on **neutral-molecule RDB** (default lower bound for `--chem-rdb` is `0.0`).
  - `--rdb-strict` — apply RDB bounds strictly against the organic baseline (`S=2`, `P=3`) instead of accepting any valid valence in `[S=2..6, P=3..5]`.

### RDB Convention (Important Note for Xcalibur / MassLynx Users)

`calcmass` always separates the ionizing adduct (`+H+`, `-H+`, `+Na+`, `+2H2+`, etc.) from the **neutral molecule `M`** and computes RDB strictly for the **neutral formula `M`**:

$$\mathrm{RDB}(M) = 1 + \frac{1}{2}\sum_{e} (v_e - 2)\, n_e$$

This differs from vendor tools such as **Thermo Xcalibur** or **Waters MassLynx**, which report the **ionic formula** (including the ionizing proton/cation) and therefore display **half-integer RDB** for even-electron ions (`[M+H]+`, `[M-H]-`):

| Species for `M = C16H17N3O7` | `calcmass` (neutral `M` convention) | Xcalibur / MassLynx (ionic formula convention) |
| :--- | :--- | :--- |
| Neutral `M` | `C16H17N3O7`, **RDB = `10.0`** | `C16H17N3O7`, RDB = `10.0` |
| Protonated `[M+H]+` | `C16H17N3O7`, **RDB = `10.0`** | `[C16H18N3O7]+`, RDB = `9.5` |
| Deprotonated `[M-H]-` | `C16H17N3O7`, **RDB = `10.0`** | `[C16H16N3O7]-`, RDB = `10.5` |
| Sodiated `[M+Na]+` | `C16H17N3O7`, **RDB = `10.0`** | `[C16H17N3O7Na]+`, RDB = `9.5` |

**Why `calcmass` uses neutral-molecule RDB:**
- RDB stays an intrinsic structural property of the molecule (rings + $\pi$-bonds), identical across `MH+`, `MH-`, `MNa+`, and `M2H++` spectra — matching the neutral-formula convention commonly used in metabolomics, natural-product, and organic-chemistry publications (as well as SIRIUS and Agilent MassHunter neutral tables).
- In MS(*n*) fragmentation analysis, precursor, product ion, and neutral loss formulas balance directly without proton bookkeeping (`C16H17N3O7 = C12H9N3O3 + C4H8O4`).
- Consequently, `--chem-rdb` always selects **integer, non-negative RDB** (closed-shell neutral molecules) even when searching `MH+` or `MH-`. If you specifically want to inspect odd-electron/radical neutral compositions (half-integer neutral RDB), use `--chem` without `-rdb`.
- For variable-valence elements — **sulfur** (`S`, valence `2, 4, 6`) and **phosphorus** (`P`, valence `3, 5`) — RDB is reported as an interval `RDB_min..RDB_max` (e.g., `1.0..3.0` for cysteine `C3H7NO2S`, `-2.0..0.0` for `SF6`, `0.0..1.0` for `H3PO4`). Parity never changes because valence steps are even (`±2`). By default `--chem-rdb` keeps a formula if any valid valence assignment has `RDB >= 0` (`RDB_max >= 0`); pass `--rdb-strict` to enforce `RDB_min >= 0` (divalent `S=2`, trivalent `P=3`).

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

- **Formula**: `C16H17N3O7`, `C16,H17,N3,O7`, `C2H5OH`, `'Ca(OH)2'`, or `'(CH3)2SO'` (repeated elements are summed; nested parentheses are supported).
  > **Shell note:** Always enclose formulas containing parentheses in single quotes (e.g., `./calcmass --mass M '(CH3)2SO'`), since `(` and `)` are shell metacharacters in `bash`/`sh`/`zsh` and will trigger a shell syntax error if left unquoted.
- **Limits**: `<Elem><min>-<max>` or fixed count `<Elem><n>` / `<Elem>` (implies `1`), e.g. `C0-30,H0-60,N0-5,O0-20` or `C8-10,N6,O1-2,H1-13,Cl`.
- **Target mass (`<target_mass>`)**: accepts either a single positive number (`120.0427`) or an additive `+`/`-` expression (including Unicode minus `−` and en-dash `–`), e.g., `'362.1003 - 242.0576'` or `'242.0576 + 120.0427'` (quote expressions containing spaces in the shell).
- **Search window**: prints up to `--top` hits (default `30`) within `--da` (default `0.1 Da`) and `--ppm` (default `20 ppm`) of `<target_mass>` (or the single closest match if none fall within the window).

---

## Examples

### 1. Exact Mass Calculation for a Known Formula (`--mass`)

Calculate the theoretical mass, *m/z*, elemental breakdown, and RDB of the deazaflavin bioluminescence cofactor **ActH** (archaeal cofactor F0, 7,8-didemethyl-8-hydroxy-5-deazariboflavin, `C16H17N3O7`) [\[2\]](#references) in negative mode (`MH-`):

```bash
./calcmass --mass MH- C16H17N3O7
```

Full output:

```text
calcmass by Maxim A. Dubinnyi
Version 2026-11-27
Mass type is '[M-H]-' (4), charge= -1, madd=-1.007276

Parsed formula: 'C16,H17,N3,O7', RDB = 10.0

Formula:
 C  16 12.000000  192.000000
 H  17  1.007825   17.133026
 N   3 14.003074   42.009222
 O   7 15.994915  111.964402
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

### 2. Formula Search Without RDB Filtering (`--chem -s`) — Compound AsLn7 [\[1\]](#references)

In the structure elucidation of **Compound AsLn7** — a natural luciferin analogue from the bioluminescent Siberian earthworm *Fridericia heliota* consisting of the tyrosine-derived chromophore **CompX** conjugated with γ-aminobutyric acid (GABA) [\[1\]](#references) — ESI+ HRMS gives a protonated molecular ion `[M+H]+` at *m/z* `324.1067`. Searching across `C0-30,H0-60,N0-5,O0-20` without RDB filtering:

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

Notice that the two closest raw mass matches are chemically impossible or odd-electron neutral compositions: `H23,N2,O17` (RDB `-9.5`) and `C13,H15,N4,O6` (RDB `8.5`).

---

### 3. Filtering by Unsaturation / RDB (`--chem-rdb -s`) — Compound AsLn7 [\[1\]](#references)

Running the same query for **AsLn7** with `--chem-rdb` discards negative and half-integer neutral-molecule RDB values, leaving only valid closed-shell formulas:

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

All nonsensical hits are eliminated, and the true formula of **AsLn7** — **`C15H17NO7`** (RDB `8.0`, +1.078 mDa / +3.33 ppm) — ranks #1 as the sole candidate within 5 ppm (verified with `--ppm 5`) [\[1\]](#references).

---

### 4. Multi-Stage MS(*n*) Fragmentation Interpretation — Deazaflavin Cofactor (ActH / ActS) [\[2\]](#references)

In the identification of the *Henlea* sp. bioluminescence activators **ActH** (deazaflavin cofactor F0, `C16H17N3O7`) and its novel 3′-O-sulfated derivative **ActS** (`C16H17N3O10S`) [\[2\]](#references), full-scan and MS2 masses alone admit multiple valid closed-shell formulas. Combining multi-stage MS(*n*) daughter ions and neutral losses resolves the ambiguity (data from [\[2\]](#references), Fig. 2 and ESI Tables S1–S2, Figs. S1–S3).

#### Step A — ActS MS2 ion (`[M-H]-` = `362.10031`, ESI Table S2 of [\[2\]](#references))

In ESI− mode, the ActS molecular ion `[M-H]-` at *m/z* `442.06003` yields a single MS2 daughter ion at *m/z* `362.10031` via neutral loss of `79.95972 Da` (`SO3`, exact `79.956815 Da`, −2.91 mDa / −36.3 ppm; preferred over `HPO3` at `79.966331 Da`, +6.61 mDa / +82.7 ppm, and confirmed by absence of ¹H–³¹P NMR splitting [\[2\]](#references)). Searching the `362.10031` ion (`ActH` core):

```bash
./calcmass -s --chem-rdb --ppm 15 MH- C5-30,H0-30,N0-10,O0-10 362.10031
```

```text
                  Formula      Theory     Mass diff   Mass diff   RDB
            C17,H13,N7,O3: 362.100711     0.401 mDa   1.107 ppm  15.0
            C16,H17,N3,O7: 362.099373    -0.937 mDa  -2.586 ppm  10.0
            C5,H17,N9,O10: 362.102561     2.251 mDa   6.217 ppm   2.0
                C28,H13,N: 362.097523    -2.787 mDa  -7.697 ppm  23.0
             C21,H17,N,O5: 362.103396     3.086 mDa   8.523 ppm  14.0
            C12,H13,N9,O5: 362.096688    -3.622 mDa -10.002 ppm  11.0
             C22,H13,N5,O: 362.104734     4.424 mDa  12.217 ppm  19.0
            C11,H17,N5,O9: 362.095351    -4.959 mDa -13.696 ppm   6.0
```

At MS2 alone, `C17H13N7O3` (+1.11 ppm) actually lies closer to `362.10031` than `C16H17N3O7` (−2.59 ppm).

#### Step B — Neutral loss (`M` = `362.10031 − 242.05761` = `120.04270`) and MS3 fragment (`242.05761`)

Using adduct type `M` with inline mass subtraction identifies the neutral loss, while `MH-` identifies the MS3 fragment ion:

```bash
./calcmass -s --chem-rdb M C0-10,H0-20,N0-6,O0-8 '362.10031 - 242.05761'
```

```text
                  Formula      Theory     Mass diff   Mass diff   RDB
                 C4,H8,O4: 120.042259    -0.441 mDa  -3.676 ppm   1.0
                 C5,H4,N4: 120.043596     0.896 mDa   7.465 ppm   6.0
```

```bash
./calcmass -s --chem-rdb MH- C0-20,H0-20,N0-10,O0-6 242.05761
```

```text
                  Formula      Theory     Mass diff   Mass diff   RDB
             C12,H9,N3,O3: 242.057115    -0.495 mDa  -2.046 ppm  10.0
               C,H9,N9,O6: 242.060303     2.693 mDa  11.124 ppm   2.0
               C8,H5,N9,O: 242.054429    -3.181 mDa -13.140 ppm  11.0
               C17,H9,N,O: 242.061137     3.527 mDa  14.573 ppm  14.0
              C7,H9,N5,O5: 242.053092    -4.518 mDa -18.665 ppm   6.0
```

#### Step C — Ribityl neutral loss (`M` = `362.10031 − 228.04200` = `134.05831`) and aromatic core (`228.04200`)

```bash
./calcmass -s --chem-rdb M C0-10,H0-20,N0-6,O0-8 '362.10031 - 228.04200'
```

```text
                  Formula      Theory     Mass diff   Mass diff   RDB
                C5,H10,O4: 134.057909    -0.401 mDa  -2.993 ppm   1.0
                 C6,H6,N4: 134.059246     0.936 mDa   6.984 ppm   6.0
```

```bash
./calcmass -s --chem-rdb MH- C0-20,H0-20,N0-10,O0-6 228.04200
```

```text
                  Formula      Theory     Mass diff   Mass diff   RDB
             C11,H7,N3,O3: 228.041465    -0.535 mDa  -2.348 ppm  10.0
                 H7,N9,O6: 228.044653     2.653 mDa  11.632 ppm   2.0
               C7,H3,N9,O: 228.038779    -3.221 mDa -14.123 ppm  11.0
               C16,H7,N,O: 228.045487     3.487 mDa  15.293 ppm  14.0
              C6,H7,N5,O5: 228.037442    -4.558 mDa -19.988 ppm   6.0
```

Although `C17H13N7O3` had a smaller MS2 mass error (+1.11 ppm), matching it to the MS3 daughters (`C12H9N3O3`, `C11H7N3O3`) requires `C5H4N4` (+7.47 ppm) and `C6H6N4` (+6.98 ppm) neutral losses — exceeding the 5 ppm threshold. **`C16H17N3O7`** (`C12H9N3O3 + C4H8O4 = C11H7N3O3 + C5H10O4 = C16H17N3O7`, also consistent with the `H2O` loss `362.10031 → 344.08948`) is the **only** formula for which the parent ion, all neutral losses, and all MS3 daughter ions simultaneously agree within 5 ppm (for ActH itself in ESI−, Full MS `362.09832` gives the same unique self-consistent solution; see Table S1 of [\[2\]](#references)).

#### Step D — Positive-mode MS4 unambiguous anchor (`[M+H]+` = `132.04555`, Fig. S1 of [\[2\]](#references))

In ESI+ mode, following the deepest fragmentation path `364.12003 → 230.05746 → 187.05197 → 132.04555` reaches a low-mass MS4 fragment with a **single unambiguous formula** across `C0-20,H0-40,N0-10,O0-20`:

```bash
./calcmass -s --chem-rdb MH+ C0-20,H0-40,N0-10,O0-20 132.04555
```

```text
                  Formula      Theory     Mass diff   Mass diff   RDB
                C8,H5,N,O: 132.044390    -1.160 mDa  -8.783 ppm   7.0
```

Tracing back to the root (`C8H5NO → C10H6N2O2 → C11H7N3O3 → C16H17N3O7`) independently confirms **ActH** = `C16H17N3O7` and **ActS** = `C16H17N3O10S` (`ActH + SO3`) [\[2\]](#references).

---

## Project Structure

- `calcmass.c` — CLI entry point, threshold/RDB filtering, candidate ranking, and formatted output.
- `isotopes.h` / `isotopes.c` — Element table (AME2020/CODATA 2022 masses and valence ranges), adduct definitions, mass/*m/z*/RDB calculators, and sort comparator.
- `parse.c` — Strict chemical formula (including nested parentheses), elemental range, and target-mass expression parsers.
- `multiloop.c` — Multi-dimensional iterator over elemental search limits.
- `tests/test_unit.c`, `tests/test_cli.sh` — C unit tests and CLI integration test suite.
- `Makefile` — Build, test, sanitizer, and install configuration (`gcc`/`clang`, C11, strict warnings).
- `IsotopeMassTable.txt`, `Isotope_Table_NIST.txt` — Reference isotope mass tables.

## References

1. Dubinnyi MA, Tsarkova AS, Petushkov VN, *et al.* (2015) Novel Peptide Chemistry in Terrestrial Animals: Natural Luciferin Analogues from the Bioluminescent Earthworm *Fridericia heliota*. *Chem. Eur. J.* **21**:3942–3947. https://doi.org/10.1002/chem.201406498
2. Petushkov VN, Vavilov MV, Ivanov IA, *et al.* (2023) Deazaflavin cofactor boosts earthworms *Henlea* bioluminescence. *Org. Biomol. Chem.* **21**:415–427. https://doi.org/10.1039/D2OB01946A

## Author

Dr. Maxim A. Dubinnyi

ORCID https://orcid.org/0000-0002-5463-0169

Shemyakin-Ovchinnikov Institute of bioorganic chemistry, 
Moscow, Russia (http://www.ibch.ru/en)
