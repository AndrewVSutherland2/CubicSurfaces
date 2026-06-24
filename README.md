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
descends the twisted moduli space over the splitting field of `f`, searches the
resulting rational fourfold in `P^9` for a rational point, and converts the
recovered Clebsch invariants into a pentahedral quintic and finally a cubic
surface in `P^3`.

## Contents

| File | Description |
|------|-------------|
| `cubic_surface_resolvent_twist.m` | The implementation (all entry points below). |
| `examples.m`                      | Builds and certifies the four cyclic examples `C4, C5, C9, C12`. |
| `experiments/point_search_helpers.m` | Exploratory point-search routines (Minkowski reduction, bounded integer search, fibre search) used while investigating the harder non-cyclic twists. |

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

`BuildUniversalCobleData` constructs, once, the 40 Coble gamma invariants, the 80
signed coordinates, their 10-dimensional linear span and the 30 cubic relations,
together with the `W(E6)` action on the span. For a target subgroup `G` and a
polynomial `f` with `Gal(f) ≅ G`, the code forms the semilinear descent condition
`R_{ρ(σ)}·σ(y) = y` over the splitting field `L` of `f`, solves it for the
fixed `Q`-rational 10-dimensional space, and substitutes a reduced integral basis
of that space into the 30 universal cubics to obtain the descended twist — a
rational fourfold in `P^9`. A rational point on the twist is found by slicing it
to dimension zero with coordinate and random hyperplanes and solving exactly; its
coordinates give the Clebsch invariants `[i8 : … : i40]`, from which a pentahedral
quintic and then a cubic surface in `P^3` are reconstructed. Finally
`LinesGaloisGroupOfCubicSurface` certifies the result independently from the 27
lines of the output surface.

Two coordinate reductions keep the construction fast: the fixed space is taken in
a **saturated** LLL basis, and the 30 descended cubics are replaced by an
LLL-reduced basis of their integer coefficient lattice. Without these the cubics
carry enormous coefficients (growing with `[L:Q]`); with them the cyclic searches
run in seconds.

## Status and limitations

- **Cyclic groups work and are certified** end to end: `C4, C5, C9, C12` (above).
- The **descent is correct for every group** (verified separately), but the
  **rational-point search is heuristic** (this is the heuristic step of EJ's
  Algorithm 5.1). For several non-cyclic / very small groups (`C2, C3, S3, D4,
  C3²`, …) the easily-accessible rational points of the twist are all degenerate
  (singular surfaces), and the slicing / bounded / line searches here did not
  reach a smooth point. Elsenhans–Jahnel reached such points with a
  Minkowski-reduced `O_L¹⁰ ∩ V` integral-lattice search (and a choice of
  favourable `f`); matching them on the remaining six "hard" classes (gap
  nos. 155, 169, 177, 179, 266, 286; no. 73 = `C9` is reproduced here) is the
  natural next step.

## References

- A.-S. Elsenhans and J. Jahnel, *On the inverse Galois problem for cubic
  surfaces*, [arXiv:1209.5591](https://arxiv.org/abs/1209.5591).
- A. B. Coble, *Point sets and allied Cremona groups*, Trans. AMS **16** (1915).
