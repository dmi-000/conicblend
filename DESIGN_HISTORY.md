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

---

## 2026-03-22  CylinderWindow integration  (conicblend_cylinder.hpp)

**Motivation: why cylinders instead of planes**

`ConicWindow<3>` projects 3D control points onto their best-fit plane and fits
a 2D conic.  This discards torsion: the plane projection of a curve with
non-zero torsion τ introduces an O(κτh³) geometric error per segment that
cannot be eliminated by adding more control points, because the error is
structural — the curve leaves the plane between knots.

A right circular cylinder is the natural next step:

1. A helix (canonical constant-κ, constant-τ curve) is a geodesic on its
   cylinder.  It unrolls to a straight line, which any interpolator
   reproduces exactly.
2. More generally, the osculating cylinder captures both curvature and torsion
   locally, whereas the osculating plane captures only curvature.
3. As the window shrinks, the best-fit cylinder through 5 nearby points
   converges to the osculating cylinder, giving an O(h⁴) or better geometric
   approximation vs O(h³) torsion leakage from the plane.

In short: planes are the right primitive for 2D curves; cylinders are the
right primitive for 3D curves with torsion.

**Architecture: separate opt-in `cylinder_tag{}`**

Integrated the 2D-Newton cylinder solver as `CylinderWindow<3>` in
`conicblend_cylinder.hpp`.  Available via `blend_curve(ctrl, times,
cylinder_tag{})`.  Kept as a separate opt-in alongside `ConicWindow`:
both coexist until evidence proves one is universally better.

**Fallback chain per window:**
1. `CylinderWindow<3>`: find best-fit cylinder → unroll → `ConicWindow<2>` or
   `LagrangeWindow<2>` → re-roll.
2. `ConicWindow<3>` (no real cylinder found): project to best-fit plane.
3. `LagrangeWindow<3>` (both above fail): at `blend_curve` level.

**Geodesic degeneracy (helix on its own cylinder)**

A helix is a geodesic on its cylinder, so it unrolls to a straight line in
(r·φ, z) space.  `ConicWindow<2>` correctly detects collinearity (via its
PCA eigenvalue guard) and returns `valid_=false`.  We fall back to
`LagrangeWindow<2>`, which reproduces a straight line exactly (degree-4
Lagrange is exact for degree-1 data).  This gives near-machine-precision
helix interpolation (~2e-15 error).

**Cylinder selection: monotone-φ + geodesic-first sort**

For 5 points on a helix, up to 6 distinct exact-fit cylinders exist (all with
knot_err ≈ 1e-16).  A naïve sort by φ-span picks a non-monotone oblique
cylinder (span=2.5) rather than the true helix cylinder (span=7.2), leading to
errors of ~1.8 (10× *worse* than conic).  Sorting by φ-span among monotone
cylinders alone is also insufficient: an oblique monotone cylinder (span=5.4)
is preferred over the true helix cylinder.

**Fix:** Two-level sort in `cyl_solve`:
1. Monotone-φ cylinders first.  A non-monotone φ sequence means the control
   points' angular projections onto the cylinder fold back as t increases
   monotonically — a non-injective parametrization.  The cylinder is simply a
   poor frame for this arc; it does not mean the 3D curve is degenerate.
2. Among monotone cylinders: geodesic cylinders (`lam_ratio < 1e-4`, where
   `lam_ratio = lam_min/lam_max` of the 2D unrolled scatter) come first,
   because they give *exact* interpolation via `LagrangeWindow<2>`.
3. Among non-geodesic cylinders: sort by φ-span ascending (smaller = more
   compact = more stable).

**Measured results (demo_cylinder.cpp):**

| Curve | n | conic max_dev | cylinder max_dev | improvement |
|-------|---|---------------|-----------------|-------------|
| Helix | 8 | 9.6e-02 | 1.9e-15 | ~5×10¹³× |
| Helix | 12 | 1.1e-02 | 1.9e-15 | ~6×10¹²× |
| Helix | 20 | 1.1e-02 | 2.4e-15 | ~5×10¹²× |
| Twisted cubic | 8 | 1.3e-03 | 2.5e-05 | 53× |
| Twisted cubic | 12 | 4.7e-04 | 2.7e-06 | 175× |
| Twisted cubic | 20 | 1.0e-04 | 1.7e-07 | 601× |

**Proved property:** for a curve that lies exactly on a cylinder (e.g., helix),
`CylinderWindow<3>` achieves near-machine-precision interpolation, independent
of n (number of control points).  The geodesic-first sort is essential for this.

## Pending / open questions

- **Replace ConicWindow with CylinderWindow for 3D?** Open.  The measured data
  shows cylinder is always better for the two test curves.  Needs broader testing
  across curve families before replacing.  Document proof/evidence here when done.

- **Non-planar data and the plane fallback chain**: when `cyl_solve` returns no
  real cylinders for a window, the current fallback is ConicWindow<3>.  This is
  the 3D analogue of the existing fallback chain.
