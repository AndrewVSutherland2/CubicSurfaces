# CubicSurfaces

A [Magma](http://magma.maths.usyd.edu.au/) implementation of **Elsenhans–Jahnel
Algorithm 5.1** ([arXiv:1209.5591](https://arxiv.org/abs/1209.5591), *On the
inverse Galois problem for cubic surfaces*): given a finite group `G` realised as
a subgroup of the Weyl group `W(E6)` (equivalently, a prescribed Galois action on
the 27 lines), construct an **explicit smooth cubic surface over `Q`** whose
Galois action on its 27 lines is `G`.

The implementation works through the two Galois descents of the paper using
Coble's irrational invariants ("gamma coordinates") on the moduli space of
marked cubic surfaces: it builds the universal 40 invariants and 30 cubic
relations once, computes the `W(E6)` action on the 10-dimensional linear span,
descends the twisted moduli space over the splitting field of `f`, **produces a
rational point of the twist with an explicit dominant map from affine 6-space**
(the A6 cuspidal-cubic producer — no brute-force search), and converts the
recovered Clebsch invariants into a pentahedral quintic and finally a cubic
surface in `P^3`.

## Contents

| File | Description |
|------|-------------|
| `cubic_surface_resolvent_twist.m` | The core implementation (all entry points below), including the **A6 cuspidal-cubic point producer**. |
| `database_pipeline.m`             | Generator for **many** surfaces per number field — `GenerateDatabase(f, …)` (see below). |
| `examples.m`                      | Builds and certifies the four cyclic examples `C4, C5, C9, C12`. |
| `algorithm_5_1_explicit_A6_map.tex` | The note describing the explicit dominant `A^6 →` moduli map used to avoid the point search. |
| `cubic_surface_discriminant_clebsch.tex` | The Clebsch–Salmon discriminant `Δ_Cl` (EJ Lemma 2.6), used to rank surfaces by intrinsic discriminant. |
| `experiments/point_search_helpers.m` | Exploratory routines from the original heuristic point search (now superseded by the A6 map). |

## Requirements

Magma (any reasonably recent version). No other dependencies.

## Quick start

From the repository root:

```
magma
> load "cubic_surface_resolvent_twist.m";
> Pq<t> := PolynomialRing(Rationals());
> R := CubicSurfaceFromPolynomial(t^4 + t^3 + t^2 + t + 1);   // Gal(f) = C4
> R`surface;                                                  // the cubic in P^3
> LinesGaloisGroupOfCubicSurface(R`surface);                  // independent check: C4
```

`BuildUniversalCobleData` (the one-time setup, ~1–2 s) is called automatically;
pass it explicitly via the `Universal` parameter to reuse it across many calls.

## Interface

Everything is an ordinary Magma `function` (no intrinsics), so just `load` the
file. The main entry points:

#### `BuildUniversalCobleData( : InterpolationPoints := 300, Print := true)`
Builds the group- and surface-independent data: the Weyl group `W(E6)`, the
sign/basis matrix `bm`, the 30 universal cubic relations `rel3`, the degree-3
monomials `mon3`, and the `W(E6)` action. Returns a record `U`; pass it as
`Universal := U` to the functions below to avoid rebuilding it.

#### `CubicSurfaceFromPolynomial(f : Universal := false, MaxClasses := 6, RandomSlices := 300, RandomBound := 2, Print := true)`
Convenience driver. Realises `Gal(f)` on the 27 lines **without** naming a
`W(E6)`-subgroup: it finds the `W(E6)`-subgroup classes isomorphic to `Gal(f)`
and tries them in turn (a given abstract group embeds in `W(E6)` in several
non-conjugate ways, giving genuinely different twists) until the point search
succeeds. Returns `R, ci` — the result record and the index of the subgroup
class that worked.

#### `CubicSurfaceFromSubgroup(G, f : Universal := false, Point := false, RandomSlices := 300, RandomBound := 2, Print := true, ...)`
The core routine for one **named** subgroup `G < W(E6)` and a polynomial `f`
whose Galois group is isomorphic to `G`. Performs the descent over the splitting
field of `f` and the rational-point search. If you already know a rational point
of the twist you can supply it as `Point := [...]` to skip the search. Returns a
result record `R`.

#### `CubicSurfaceFromResolventPolynomial(classNo, f : Universal := false, ...)`
As above, but selects the subgroup by Magma's `SubgroupClasses(W(E6))` index
`classNo` instead of by an explicit group.

#### `RunC9Example( : HeightBound := 2, Print := true)`
Self-contained reproduction of the paper's flagship hard case — subgroup no. 73,
cyclic of order 9 (`c9_example.m` in EJ's distribution). Returns a result record.

#### `LinesGaloisGroupOfCubicSurface(gl : Tries := 12, Print := true)`
**Independent certification.** For a smooth cubic `gl = 0`, computes the Galois
group of the field of definition of its 27 lines directly (via the Fano scheme of
the surface, using a generic primitive element so the full group — not a quotient
— is returned), as a permutation group. Use it to check that the construction
really realises the prescribed group:
`IsIsomorphic(LinesGaloisGroupOfCubicSurface(R`surface), R`subgroup)`.

#### `IsSmoothCubicSurface(gl)`
Boolean test that `gl = 0` is a smooth cubic surface.

### The result record

`CubicSurfaceFromSubgroup` / `…FromPolynomial` / `…FromResolventPolynomial`
return a record with (among others) the fields:

| Field | Meaning |
|-------|---------|
| `surface`            | the cubic form in `P^3` (reduced/minimised) |
| `surface_unreduced`  | the cubic before lattice minimisation |
| `invariants`         | the Clebsch invariants `[i8, i16, i24, i32, i40]` |
| `pentahedral_quintic`| the pentahedral quintic `g(T)` |
| `point`              | the rational point found on the twist |
| `subgroup`           | the `W(E6)` subgroup used |

## Reproducing the cyclic examples

```
magma
> load "examples.m";
```

`examples.m` builds a smooth cubic surface over `Q` for each cyclic group below
and certifies it by recomputing the Galois action on the 27 lines of the *output*
surface. Each surface is smooth over `Q` and its 27-line Galois group is
isomorphic to the prescribed group. Approximate run times are on a single core.

| Group | Defining polynomial `f` | Splitting field used | runtime |
|-------|--------------------------|------------------------|---------|
| `C4`  | `t^4 + t^3 + t^2 + t + 1`                | `Q(ζ₅)` (5th cyclotomic field) | ~1–3 min\* |
| `C5`  | `t^5 + t^4 − 4t^3 − 3t^2 + 3t + 1`       | `Q(ζ₁₁)⁺`, the real cyclic quintic of conductor 11 | ~2 s |
| `C9`  | `t^9 − 19t^8 + 152t^7 − 665t^6 + 1729t^5 − 2717t^4 + 2508t^3 − 1254t^2 + 285t − 19` | the cyclic degree-9 field of conductor 19 (subfield of `Q(ζ₁₉)`) | ~2 s |
| `C12` | `t^12 + t^11 + ⋯ + t + 1`                | `Q(ζ₁₃)` (13th cyclotomic field) | ~5 s |

\* `C4` is slower only because the first two `W(E6)`-subgroup classes it tries
yield twists with no easily-found point; the third class succeeds.

You can also build any one of them directly, e.g.

```
> Pq<t> := PolynomialRing(Rationals());
> R := CubicSurfaceFromPolynomial(t^5 + t^4 - 4*t^3 - 3*t^2 + 3*t + 1);
> LinesGaloisGroupOfCubicSurface(R`surface);   // C5
```

## Number fields used

All four cyclic examples descend over an **abelian (cyclic) number field of small
prime conductor** — the splitting field `L` of the defining polynomial. Because
each `L` is Galois over `Q`, the splitting field equals the field generated by a
root.

| Group | Field `L` | degree | conductor | `disc(O_L)` |
|-------|-----------|--------|-----------|-------------|
| `C4`  | `Q(ζ₅)`                          | 4 | 5  | `5³ = 125` |
| `C5`  | `Q(ζ₁₁)⁺` (maximal real subfield) | 5 | 11 | `11⁴ = 14641` |
| `C9`  | degree-9 subfield of `Q(ζ₁₉)`    | 9 | 19 | `19⁸` |
| `C12` | `Q(ζ₁₃)`                         | 12 | 13 | `13¹¹` |

The descent itself is field-agnostic — `CubicSurfaceFromPolynomial(f)` accepts
any `f` whose Galois group embeds in `W(E6)`; these particular cyclotomic /
cyclotomic-subfield choices simply give clean, small-conductor instances of each
cyclic group.

## How it works (one paragraph)

`BuildUniversalCobleData` constructs, once, the 40 Coble gamma invariants, their
10-dimensional linear span and the 30 cubic relations, the `W(E6)` action on the
span, and **Pinkham's 6-dimensional reflection representation** of `W(E6)` on the
parameters of six points of a cuspidal cubic. For a target subgroup `G` and a
polynomial `f` with `Gal(f) ≅ G`, the code forms the semilinear descent condition
`R_{ρ(σ)}·σ(y) = y` over the splitting field `L` of `f` and solves it for the
fixed `Q`-rational 10-dimensional space (the descended twist). A rational point of
the twist is then produced **without any search**, by an explicit **dominant
rational map `Φ_ρ : A^6 ⇢ M̃_ρ`** (the A6 cuspidal-cubic producer of
`algorithm_5_1_explicit_A6_map.tex`): six points `p(t)=(t:t^3:1)` on the fixed
cuspidal cubic blow up to a marked cubic surface, Pinkham's representation is
twisted by the same cocycle `ρ`, and Coble's gammas restricted to the cusp give
the point. Because the map is dominant, a generic small integer vector `u ∈ Z^6`
lands on a smooth surface — the degenerate locus is a proper closed subset to
avoid, not a wall. The recovered Clebsch invariants `[i8 : … : i40]` give a
pentahedral quintic and finally a cubic surface in `P^3`, certified independently
from the 27 lines by `LinesGaloisGroupOfCubicSurface`.

This replaces the heuristic rational-point search of EJ's Algorithm 5.1, which
failed for the non-cyclic / small twists (every easily-accessible point of the
descended fourfold was degenerate). The cuspidal map needs no descended cubics at
all; `UseA6 := true` is the default path of `CubicSurfaceFromSubgroup`.

## Example cubic surfaces

One small cubic surface `= 0` per group, each realising the named group on its 27
lines (certified independently via `LinesGaloisGroupOfCubicSurface`):

```
C2   (t^2 - 2):                  -x^2*z + x*y*z - x*y*w - 5*x*z^2 - 8*x*z*w - 3*x*w^2 - y^2*w
                                  + 2*y*z^2 + 3*y*z*w - 12*z^3 + 12*z^2*w + 27*z*w^2 + 9*w^3
C3   (t^3+t^2-2t-1):              2*x^2*z + 4*x*y*z + 4*x*z^2 + 3*x*z*w + 2*x*w^2 + y^3 + 2*y^2*w
                                  + y*z^2 - 3*y*z*w - 3*y*w^2 + z^3 - 10*z^2*w - 14*z*w^2 - 3*w^3
C4   (t^4+t^3+t^2+t+1):           -x^2*w + x*y^2 + x*y*w + 5*x*z*w - 4*x*w^2 + y^3 + 4*y^2*z - 2*y^2*w
                                  + 2*y*z^2 + 3*y*z*w - 4*y*w^2 + 4*z^3 + z^2*w + 5*z*w^2 - 10*w^3
C5   (t^5+t^4-4t^3-3t^2+3t+1):    x^2*y + x^2*w - 2*x*y*z - 2*x*z*w + x*w^2 + 2*y^2*w - 2*y*z^2
                                  - 5*y*z*w - 2*y*w^2 - z^3 - 3*z^2*w - z*w^2
C9   (t^9-19t^8+...+285t-19):     x^2*w - x*y^2 + 2*x*y*z + x*y*w - x*z^2 + 5*x*z*w - 2*x*w^2 + 2*y^2*z
                                  - 3*y^2*w - 2*y*z^2 - 15*y*z*w + 21*y*w^2 - z^3 + 7*z^2*w - 21*z*w^2 + 16*w^3
C12  (t^12+t^11+...+t+1):         -2*x^2*z + x^2*w - x*y^2 + 4*x*y*z - 8*x*y*w - 2*x*z^2 - 9*x*z*w - 4*x*w^2
                                  - 3*y^2*z - y^2*w - 8*y*z*w - 11*y*w^2 - z^3 + 3*z^2*w - 2*z*w^2 - w^3
S3   (t^3-t-1):                  -15*x^2*w - 5*x*y*w + 5*x*z^2 - 17*x*z*w - 5*x*w^2 + 18*y^3 + 25*y^2*z
                                  + 20*y^2*w + 8*y*z^2 + 18*y*z*w - 21*y*w^2 + 4*z^3 - 16*z^2*w - 58*z*w^2 + 268*w^3
```

The cyclic groups produce very small models in seconds (`C9`, a degree-9 field,
reduces to `max|coef| = 4`); `S3` is a "large-discriminant" twist where only
larger models are currently reachable (see Status).

## Database generation (`database_pipeline.m`)

`database_pipeline.m` (load after the core file) turns the A6 producer into a
generator of *many* surfaces per number field:

```
load "cubic_surface_resolvent_twist.m";
load "database_pipeline.m";
U    := BuildUniversalCobleData();
rows := GenerateDatabase(f, U, "/tmp/cs" : PerClass := 5);
```

For `f` with `Gal(f) = G`, and **each** `W(E6)`-conjugacy class of subgroups
isomorphic to `G` and **each distinct embedding `ρ`** (up to normalizer-conjugacy
— genuinely different twists), it: generates many distinct moduli points (fast —
the A6 map, no reduction); ranks them by the **Clebsch–Salmon discriminant**
`Δ_Cl = A⁴ − 128A²B + 4096B² − 2048AC − 16384D` (the intrinsic discriminant, with
`(A,…,E) = (i8,…,i40)`; EJ Lemma 2.6, `cubic_surface_discriminant_clebsch.tex`);
reconstructs the lowest-discriminant candidates and **minimizes them in parallel
under a per-job timeout**; and keeps the smallest. Each row is `⟨class, surface,
max|coef|, z, Clebsch invariants, Δ_Cl, ρ⟩`. For example,
`GenerateDatabase(t^3+t^2-2*t-1, U, dir)` returns 15 certified `C3` cubics (5 per
class × 3 classes) in ~7 s.

### Canonical `W(E6)` subgroup labels

`WE6subgroups.txt` assigns a canonical label to each of the 350 conjugacy classes
of subgroups of `W(E6) = 27T1161` (for eventual inclusion in the LMFDB), given by
generators in `S27` (the action on the 27 lines) together with the degree `d` and
T-number `t` of the minimal-degree transitive group isomorphic to the subgroup.
`RealizeCubicSurfaces(f, n, U, subs, dir)` returns a list of `⟨label, cubic form⟩`
tuples — up to `n` per label — the labels ranging over the `W(E6)` subgroup
classes isomorphic to `Gal(f)`:

```
subs   := LoadWE6Subgroups(U, "WE6subgroups.txt");
tuples := RealizeCubicSurfaces(t^3 + t^2 - 2*t - 1, 5, U, subs, "/tmp/cs");
```

`LoadWE6Subgroups` transports each file subgroup into the 40-dimensional gamma
representation via `IsIsomorphic(27T1161, we6)` — verified to agree with the
explicit 27-line action of the generators on **all 350 classes** (so the `W(E6)`
outer automorphism does not perturb the labels). For `t^3+t^2-2t-1` this yields
the three `C3` classes `51840.b.17280.{a1,b1,c1}.a1`, each realised by small
certified cubics (max|coef| ≤ 8).

## Status and limitations

- **The point search is resolved.** The A6 cuspidal map is dominant, so a smooth
  marked surface is found for **every** twist by enumerating small `u ∈ Z^6` —
  including the non-cyclic / small groups (`C2, C3, S3, …`) where the old
  heuristic search failed entirely. Certified end to end: `C2, C3, C4` (all four
  `W(E6)` classes), `C5, C9, C12, S3`.
- **Field degree is not the blocker.** `C9` (degree-9 splitting field) reduces to
  `max|coef| = 4` in ~3 s.
- **Open: small explicit models for "large-Clebsch" twists.** The difficulty of
  producing a *small* explicit equation is governed by the twist's
  Clebsch-invariant height (set by the descent basis), not the field degree.
  Cyclic / small twists give tiny `Δ_Cl` (≤ ~30 digits) and minimize fast;
  `S3`-type twists give `Δ_Cl ≳ 60` digits, whose models are badly non-minimal
  (large spurious primes in the model discriminant), so full minimization via
  `MinimizeReduceCubicSurface` is expensive. The pipeline's parallel-timeout
  harvest still extracts whatever minimizes, but full coverage of these twists
  awaits a lower-height reconstruction or a smaller descent basis. (The first ten
  `S3` cubic fields from the LMFDB are all in this regime.)

## References

- A.-S. Elsenhans and J. Jahnel, *On the inverse Galois problem for cubic
  surfaces*, [arXiv:1209.5591](https://arxiv.org/abs/1209.5591).
- A. B. Coble, *Point sets and allied Cremona groups*, Trans. AMS **16** (1915).
