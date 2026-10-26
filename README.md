# A polynomial counterexample of degree 29,984

This folder contains an exact polynomial, its certificate, and the code needed
to check the stored evidence or recompute the interval bounds. It accompanies
**Tom Adamczewski and Matthew J. Colbrook, _From boundary geometry to
counterexamples in Smale's mean value problem_.**

The supplied polynomial $p\in\mathbb Q(i)[z]$ has degree $29\,984$,
satisfies $p(0)=0$ and $p'(0)=1$, and obeys

~~~text
1.00000049595801 < min_{p'(c)=0} |p(c)/c| < 1.00000049595804.
~~~

All 29,983 critical points are simple, and exactly one attains the minimum.
Every critical-point quotient therefore has modulus greater than one.
Read [proof.pdf](proof.pdf) for the mathematical argument; its LaTeX sources
are included.

## Start here: three kinds of check

Run commands from this folder. Python 3.10 or later is required; no third-party
Python package or `pip` installation is needed.

1. **File integrity:** check the files against the supplied hash manifest.
   This detects missing or changed files; it does not prove the mathematics.

   ~~~sh
   python3 verify_integrity.py
   ~~~

2. **Stored certificate:** check the supplied inputs and certificate records.

   ~~~sh
   python3 check_stored_certificate.py
   ~~~

   On Windows, use `py -3 check_stored_certificate.py`. This route needs no
   native compiler. It writes three audit JSON files. Its `PASS` result means
   that the stored-record checks passed: **it does not recompute the interval
   calculations**.

3. **Fresh native verification:** build the verifier and recompute the interval
   enclosures, using the Linux configuration described below.

   ~~~sh
   python3 certify.py --jobs 4
   ~~~

   This includes the local arithmetic checks, all 16 exhaustive partitions,
   refinement of 316 disks, and the final audits. It replaces the stored
   certificate logs and summaries.

The stored checker checks exact input serialization, all 29,983 disk records,
rational containing-disk disjointness, strict printed comparisons,
accepted-disk containment, normalization in all 16 partitions, refined bounds,
and uniqueness of the minimum. It also checks the stored model-anchor
comparison.

## The exact polynomial

[polynomial.json](polynomial.json) contains every parameter and coefficient
array in the finite definition. No external coefficient table is needed.
Every decimal or scientific-notation numerical string denotes an exact rational
number. Read these strings as exact data, not binary floating-point values.

[parameters.txt](parameters.txt) is the corresponding exact rational input for
the native verifier; `make_parameters.py` checks that the two agree.

| Mathematical array | JSON field | Complex coefficients |
| --- | --- | ---: |
| $d_{\ell k}$ in $R_0$ | `correction_coefficients` | 1,024 |
| $\eta_j$ in the cubic remainder | `free_coefficients` | 4 |
| $x_{\ell k}$ in $X/(z-a)^2$ | `extra_correction_coefficients` | 2,048 |
| $t_{\ell k}$ in the high-degree tail | `tail_coefficients` | 128 |

The two-dimensional arrays include explicit `ell` and `k` indices.
The formulas and finite recurrences in [proof.tex](proof.tex) determine the
remaining coefficients.

The construction uses the degree-16 coordinate
$q(z)=(z^{16}+iA)/(A+1/A)$, with $A=1061603/250000$, a finite
approximation to a Blaschke factor, polynomial corrections, and a high-degree
tail with nonzero leading term proportional to $q^{1874}$. The normalization is

~~~text
p(z) = [P(a_hat + z) - P(a_hat)] / P'(a_hat),
a_hat = 1.04015231339331770677 + 0.302102329270430736141 i.
~~~

The model anchor `a` is a different point. The verifier encloses the full finite
values $P(\widehat a)$ and $P'(\widehat a)$ and proves that the latter is
nonzero. At the model anchor, the same polynomial has a critical quotient
between 0.99982426177940 and 0.99982426177944, so the independent normalization
point matters.

The disks in [critical_disks.txt](critical_disks.txt) enclose zeros of **$P'$**.
Translate their centres by $-\widehat a$, leaving their radii unchanged, to
obtain disks for the critical points of $p$.

Some descriptive JSON metadata retains the historical word `candidate`.
Those fields are preserved from the construction records; they do not determine
the verification status. The certificate and [result.json](result.json) supply
the verified result.

## Recompute the interval bounds

The validated configuration is **x86-64 Linux with GNU C++17** and an IEEE 754
binary `long double` type with **64 significand bits**. The supplied calculation
used GNU C++ 14.2.0 and MPFR 4.2.2. GMP, including its C++ library and headers,
and MPFR headers and libraries are required.

On Debian or Ubuntu:

~~~sh
sudo apt-get install build-essential python3 libgmp-dev libmpfr-dev
~~~

Run `python3 certify.py --jobs 4` for the complete calculation.
`--jobs` accepts 1 through 16; `--compiler /path/to/g++` selects the compiler.
The scripts use these exact flags:

~~~text
-O2 -std=c++17 -fno-fast-math -ffp-contract=off -frounding-math
-DN_PARAM=1874 -DBLOCKS_PARAM=64 -DTAIL_BLOCKS_PARAM=8 -DEXTRA_BLOCKS_PARAM=128
~~~

All four `-D` definitions select this example. Compiling `verify.cpp` with its
built-in defaults selects a different truncation. The scripts link MPFR and
the GMP C++ and C libraries. Basic arithmetic is enclosed outward;
transcendental endpoints use MPFR directed rounding at 192 bits.
The arithmetic backend checks its floating-point assumptions, and every
certificate stage must pass. The printed-endpoint checks are calibrated to the
supplied 64-bit-significand output. An ordinary MSVC build with a 53-bit
`long double` does not meet these requirements; Windows users can run the native
verification in a suitable Linux environment.

For an optional shorter native check:

~~~sh
python3 check_local.py
~~~

This tests 64,800 exact-rational disk-arithmetic images, checks the minimizing
critical disk (index 28891), and checks the model-anchor comparison (index
27158). These checks are already included in `certify.py`; running them
separately is optional. They do not replace the exhaustive critical-point count.

A recorded exhaustive run took approximately 19 minutes, followed by about
32 seconds of refinement. These are historical timings, not a fresh verification
of this upload. Runtime depends on the machine and number of workers.

## Evidence and source files

- `verify.cpp`, `rigorous_interval.hpp`: interval verifier and arithmetic backend.
- `refine_disks.cpp`, `refinement_indices.txt`: refinement of 316 disks.
- `certificate.log`, `logs/part_00.log` through `logs/part_15.log`:
  exhaustive interval calculation.
- `certificate_refined.log` and the certificate summaries: refined minimum
  combined with the complete calculation.
- `audit_certificate.py`, `audit_exact_disks.py`, `audit_normalization.py`:
  coverage, exact-rational, and normalization checks.
- `check_base_comparison.cpp`, `base_comparison.log`: model-anchor comparison.
- `test_ball_arithmetic.cpp`, `test_containment_audit.py`: arithmetic and
  exact-containment tests.
- `result.json`: degree, normalization, strict bounds, and unique minimum.

The work-partition number printed as `sheet` in the logs is the critical-disk
index modulo 16. Local Rouché comparisons give one simple zero of $P'$ in each
accepted disk. Disjointness and the count $29\,983=\deg P'$ exhaust the
critical set.

To compile the proof, run `pdflatex proof.tex` twice.
The inputs `seed_parameters.tex` and `verified_bounds.tex` are included.
See [CITATION.cff](CITATION.cff) for citation metadata.

The core polynomial inputs and historical certificate records are preserved
byte for byte from the supplied package. `manifest.json` records the size and
SHA-256 hash of each delivered file other than the manifest itself.
Run `verify_integrity.py` before regenerating the certificate: a new calculation
legitimately replaces logs and summaries and will no longer match the delivered
snapshot. The generated audit reports are not part of that snapshot.
`.gitattributes` prevents Git from changing line endings in the hashed files.

## Checks on this package

On 26 October 2026, the stored-certificate checks passed from a separate copy
using Python 3.12 on Windows. The native local checks also passed in x86-64
Ubuntu under WSL2 with GNU C++ 13.3.0: all 64,800 arithmetic cases, the
minimizing disk, and the model-anchor comparison. The supplied proof PDF was
compiled from the included sources and visually checked.

The full exhaustive interval computation was not rerun during this packaging
check; its supplied records are included, and `certify.py` recomputes it.

## Licence

This repository is licensed under the [Apache License 2.0](LICENSE).
