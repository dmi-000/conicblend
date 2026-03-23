# Design History — conicblend

Records of significant design decisions, algorithm comparisons, and proved
properties.  The intent is that any future contributor can understand *why*
the library looks as it does, not just what it does.

---

## 2026-03-21  Initial release

**Parametrization: cross-ratio / stereographic**

The first working implementation used the cross-ratio (Möbius) parametrization
for same-branch windows and a separate φ-unwrap for cross-branch hyperbolas.
The cross-ratio is projectively invariant but requires different code paths for
each conic type.

**Architecture: 5-point windows + smoothstep blend**

Every interior segment j is blended from two overlapping 5-point windows
(win[j-2] and win[j-1]).  Because the same function evaluates from *both sides*
of each knot, C^N continuity is exact by construction — not approximate.

---

## 2026-03-21  Lagrange fallback (commit 18cab05)

**Decision: per-window fallback, not global**

When a ConicWindow fails (non-monotone φ, degenerate conic, …), the original
design fell back to the *entire curve* using circle arcs.  This violates
locality: a single bad window caused a global change in curve type.

Replaced with `LagrangeWindow<Dim>` (degree-4 Lagrange through 5 control
points) as a *per-window* fallback.  A local invalidity now has only local
effect.  This is now the library's stated locality principle.

---

## 2026-03-21  Vertex-P₀ + φ = 2·arctan(s)  (commit 1e77fbd)

**Decision: algebraic vertex as canonical P₀**

The cross-ratio parametrization depended on the choice of reference points and
produced different P₀ for the two windows covering the same segment.  If two
overlapping windows share the *same* conic but different P₀ they evaluate the
same geometric arc via different maps, and their blend is not exactly conic
(only approximately so).

The algebraic vertex (where ∂Q/∂x = 0 in the principal frame) is a property
of the conic itself, independent of which 5 points define it.  Using it as P₀
means two windows that happen to lie on the same conic get identical maps →
their blend is exact.

**Proved property:** for locally-conic data (control points exactly on a conic),
overlapping windows share the same P₀, so the smoothstep blend evaluates the
same point at every parameter value.  C^N continuity tightens to exact conic
reproduction.

**Parabola fix:** near-parabola guard `|λ₀| < 1e-6·|λ₁|` prevents the
cross-branch detector from misfiring when both eigenvalues are near zero with
opposite floating-point signs.

---

## 2026-03-22  Cylinder solver — 2D Newton approach  (diag_cylinder2.cpp)

**Context:** To handle 3D curves with torsion better than the current
plane-projection approach, we need to find all real right circular cylinders
through 5 arbitrary 3D points (Lichtblau ICMS 2006).

**Failed approach: Sylvester resultant**

The standard algebraic approach (Lichtblau):
1. Parametrize axis as {y=ax+b, z=cx+d}; 5 distance equations linear in (b,d).
2. Eliminate (b,d) via 3×3 minors G = det([F₀;F₁;F₂]) and H = det([F₁;F₂;F₃]),
   both degree-5 polynomials in c.
3. Eliminate c via the 10×10 Sylvester resultant R(a) = Res_c(G,H).
4. Find roots of R(a) (univariate degree-20 poly), back-substitute for c, b, d.

In exact arithmetic this is correct.  In double precision it fails completely:
at a genuine root a*, the 10×10 determinant R(a*) exhibits catastrophic
cancellation.  Measured ratio R(a*)/max_a|R(a)| ≈ 1e-98 for the helix — far
beyond the noise floor of double precision.  Sign-change bisection converges to
noise-generated artifacts, not true roots.

**Chosen approach: 2D Newton on (G=0, H=0)**

Instead of eliminating c to get a univariate R(a), solve the bivariate system
G(a,c)=0, H(a,c)=0 directly using 2D Newton's method from a dense 40×40
sinh-spaced starting grid.  Catastrophic cancellation never occurs because G and
H are each cheap 3×3 scalar determinants evaluated pointwise.

Results:
- Helix (exact cylinder): knot error ~1e-16 (machine precision) for most windows
- Twisted cubic: 2 real cylinders per window found consistently
- Runtime: ~30ms per 5-point window (0.47s for 16 windows)
- 3 permutations of (x,y,z) cover all axis orientations

**Design lesson:** when the diagnosis of a numerical failure shows the algorithm
is fundamentally incompatible with the precision of the medium (not merely
buggy), switch mediums or switch algorithms rather than fighting the medium.
Here: switch from "compute a degree-20 poly determinant" to "Newton on a 2-
equation system."  An alternative that was never tried but would also have
worked: pipe to a persistent `WolframKernel` process (startup cost paid once;
NSolve handles the exact arithmetic internally).

---

## Pending / open questions

- **CylinderWindow integration** *(in progress)*: `conicblend_cylinder.hpp`
  will add a `CylinderWindow<3>` class and `cylinder_tag{}` blend_curve overload
  as a separate opt-in alongside `ConicWindow`.

- **Replace ConicWindow with CylinderWindow for 3D?** Open.  If it can be
  proved (or measured across a broad curve family) that cylinder unrolling always
  gives equal or lower interpolation error than plane projection for 3D input,
  replace and note the proof/evidence here.  Until then, both coexist.

- **Non-planar data and the plane fallback chain**: when `solve_cylinders`
  returns no real cylinders for a window, the right fallback is: plane ConicWindow
  → LagrangeWindow.  This is the 3D analogue of the existing fallback chain.
