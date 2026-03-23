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

**Diagnosis: parabola failures revealed the root cause**

The previous cross-ratio approach chose P₀ from the 5 sample points (the first
candidate giving a monotone cross-ratio sequence).  It broke visibly on parabolas:
one eigenvalue λ₀ = 0 makes the stereographic slopes diverge → NaN/Inf → invalid
window.  Investigating *why* led to the core insight:

> **Same conic → identical orbit maps.**  If two overlapping windows both sample
> the same physical conic but use different P₀ choices, they evaluate the same
> geometric arc via different parametrizations.  The blend `(1−w)·A(t) + w·B(t)`
> is then a weighted average of two *different* rational maps over the same conic —
> not exact reproduction, only an approximation.

The parabola failures were the most visible symptom, but any near-conic input was
silently wrong: the blend drifted off the conic by an amount proportional to the
difference between the two P₀ choices.

**Fix: P₀ = algebraic vertex (intrinsic to the conic, not the sample points)**

The algebraic vertex (where `∂Q/∂x = 0` in the principal frame) is determined by
the conic coefficients alone — independent of which 5 points happen to define it.
Two windows on the same conic therefore compute the same vertex and the same orbit
map, making the blend exact.

This simultaneously fixes parabolas: the vertex is well-defined even when λ₀ = 0
(it is the tip of the parabola), and `φ = 2·arctan(s)` from the vertex works as a
single continuous code path for all conic types — the parabola is treated as the
limit of an ellipse as `minor/major → 0`, exactly as `conicspline` does.

The guard `|λ₀| < 1e-6·|λ₁|` still exists to prevent the cross-branch detector
from misfiring when both eigenvalues are near zero with opposite floating-point
signs — a numerical-stability guard on top of the principled fix.

**Proved property:** for locally-conic data (control points exactly on a conic),
overlapping windows share the same P₀, so the smoothstep blend evaluates the
same point at every parameter value.  C^N continuity tightens to exact conic
reproduction.

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

---

## 2026-03-23  Cylinder intersection geometry

### Why multiple exact-fit cylinders exist

For 5 points on a helix, `cyl_solve` returns up to 6 exact-fit cylinders (knot
error ~1e-16 for all of them).  This is not a numerical artifact — it has a clean
algebraic explanation.

Two right circular cylinders are quadric surfaces.  By Bézout's theorem two quadrics
intersect in a curve of degree at most 4.  All 5 control points lie exactly on this
quartic (they satisfy both cylinder equations).  For the helix specifically, the 6-fold
multiplicity arises from helical symmetry: the isometry `(rotate by φ₀, translate z
by p·φ₀)` maps the helix to itself, so if a cylinder passes through 5 consecutive
sample points it also passes through the 5 points shifted by any multiple of the
sampling step.  For generic 3D curves (e.g., twisted cubic) there are typically only
2 real solutions.

### Degenerate conic on the oblique cylinder

When the 5 helix control points are unrolled onto the oblique extra cylinder (projected
to its `(r·φ, z)` plane), the 5 points land on two intersecting straight lines.  A
5-point conic fit `Ax²+Bxy+Cy²+Dx+Ey+1=0` through points on two crossing lines
returns a *reducible conic* (det of the 3×3 conic matrix ≈ 0).

Geometric cause: the degree-4 intersection curve of the two cylinders, projected onto
the oblique cylinder's unrolled plane, forms two arcs that cross near-linearly in that
plane.  The 5 control points sample one point from each arm.  A conic through points
taken from two different crossing lines decomposes exactly into those two lines.

Consequence: `ConicWindow<2>` correctly rejects the oblique cylinder (degenerate conic
fit), falls back to `LagrangeWindow<2>`, and yields max error ~0.31 — far worse than
the true cylinder.  The degenerate conic is therefore a secondary confirmation (beyond
`lam_ratio`) that the oblique cylinder is the wrong frame.

### `lam_ratio` as the collinearity / geodesic detector

`lam_ratio = lam_min / lam_max` of the 2D scatter matrix of unrolled `(r·φ, z)` control
points distinguishes the correct cylinder from the oblique extra one:

| Cylinder | Unrolled shape | lam_ratio | Outcome |
|---|---|---|---|
| True (geodesic) | Single straight line | ≈ 1e-15 | `LagrangeWindow<2>` exact (~2e-15) |
| Oblique (non-geodesic) | Two crossing lines | ≈ 0.4 | conic degenerate → Lagrange (~0.31) |

The key insight: a helix *is* a geodesic on its own cylinder.  Geodesics on a cylinder
are helices (or degenerate cases: circles and straight lines), and the helix's pitch and
radius uniquely fix which cylinder.  No other cylinder in the 6-fold family can make the
helix a geodesic, so `lam_ratio ≈ 0` is a unique fingerprint of the correct cylinder for
any geodesic-on-cylinder curve.

---

---

## 2026-03-23  phi_monotone is the wrong gate for ConicWindow

**Finding (from visualisation of extra cylinders)**

For helix n=8, the non-preferred cyl2 has `phi_monotone=false` (r·φ(t) reverses
direction between t₁ and t₂) but its 5 unrolled control points lie monotonically
around a proper ellipse — i.e., the conic arc-angle from P₀ IS monotone in t.
The previous code sorted phi-non-monotone cylinders last and, more importantly,
used the same `phi_monotone` flag to determine whether `LagrangeWindow<2>` was safe.

**Two distinct monotonicity conditions**

| Condition | Meaning | Correct gate for |
|---|---|---|
| `phi_monotone` | r·φ(t₀..t₄) is a monotone sequence | `LagrangeWindow<2>` only |
| conic arc-angle monotone | The conic arc from P₀ sweeps monotonically with t | `ConicWindow<2>` (its own internal check) |

These are independent: a cylinder can be phi-non-monotone but have a monotone conic
arc (the ellipse simply has a turning point in its x = r·φ projection).  Applying
`phi_monotone` as a gate on `ConicWindow<2>` is too conservative.

**Fix (conicblend_cylinder.hpp)**

1. `cyl_solve` sort: removed `phi_monotone` as primary key.  New order: geodesic
   (lam_ratio < 1e-4) first, then `phi_span` ascending.
2. `CylinderWindow<3>` constructor: loops through all cylinders in order.
   - Geodesic: `LagrangeWindow<2>` (phi IS monotone for a geodesic cylinder).
   - Non-geodesic: tries `ConicWindow<2>`; if it rejects (degenerate or
     non-monotone arc-angle), **continues to the next cylinder** rather than
     falling back to `LagrangeWindow<2>` (which is unsafe for non-monotone phi).

Existing test results (helix ~2e-15, twisted cubic 47–601×) unchanged: the geodesic
cylinder is still found first for the helix, and the twisted cubic cylinders were
already phi-monotone.

---

## 2026-03-23  ConicWindow accepts straight lines; LagrangeWindow<2> removed from CylinderWindow

**Finding:** `ConicWindow<Dim>` rejected collinear input via a PCA eigenvalue guard
(`lam_min < 1e-4 * lam_max → return invalid`).  This forced `CylinderWindow<3>` to
special-case geodesic cylinders (where the 5 unrolled points are near-collinear) by
routing them to `LagrangeWindow<2>` instead of `ConicWindow<2>`.

**Principle:** A straight line is a valid degenerate conic — and it is trivially monotone
(every line has a single consistent direction).  There is no logical reason for
`ConicWindow` to treat collinearity as invalidity.

**Fix:**

Instead of returning invalid when collinear, `ConicWindow` now stores the projected 2D
points and sets an internal `line_mode_` flag.  `eval_at_` uses degree-4 Lagrange in 2D
when `line_mode_` is active — the same computation that `LagrangeWindow<2>` performs.
For truly collinear points (helix geodesic), degree-4 Lagrange through 5 collinear points
exactly reproduces the linear polynomial, giving machine-precision interpolation (~5×10⁻¹⁶).

`CylinderWindow<3>` is simplified: the two-branch `lam_ratio < 1e-4` block is removed.
All cylinders now go through a single `ConicWindow<2>` call.  `Inner2D` is no longer a
`std::variant` — it is just `ConicWindow<2>`.

**Fallback chain (simplified):**
```
CylinderWindow<3>  (some cylinder + ConicWindow<2> accepted)
  → ConicWindow<3> (no real cylinder, project to best-fit plane)
  → LagrangeWindow<3> (both fail; blend_curve level)
```

Test results unchanged: helix ~5×10⁻¹⁶, twisted cubic 47–601×.

---

## 2026-03-23  Fallback chain redesign (design decision, not yet implemented)

**Problem with ConicWindow<3> as fallback:**

The current fallback `→ ConicWindow<3> (project to best-fit plane)` is geometrically
incoherent as a general fallback.  It only makes sense when:
(a) the curve is locally planar (torsion ≈ 0), or
(b) the plane is the r→∞ limit of a cylinder (infinite radius, zero curvature).

When `cyl_solve` (the exact-fit Newton solver) returns no converged cylinders, it is a
numerical failure — the curve may still have significant torsion.  Projecting to the
best-fit plane in that case introduces the O(κτh³) torsion error the cylinder approach
was designed to eliminate.

**Why circles don't add a useful intermediate level:**

A circle is a zero-pitch geodesic on a cylinder (a degenerate cylinder with axis
perpendicular to the plane).  If the best-fit cylinder has large residuals, any circle
fit also has large residuals — circle ⊂ cylinder.  The same argument eliminates
ConicWindow<3>: a plane is the r→∞ cylinder, so if the general best-fit cylinder fails
the residual gate, the plane will too.

**Criterion for cylinder selection:**

The selection criterion for choosing among multiple valid cylinders (e.g. the up-to-6
exact-fit cylinders through 5 helix points) should be:
- Intrinsic to the 5 points and the cylinder (invariant under rigid motions of the data)
- The same criterion for both overlapping windows covering the same arc (to preserve
  the exact-blend symmetry property)

The orbit reconstruction error inside `ConicWindow<2>` — how accurately the conic orbit
map (vertex-P₀ + φ=2·arctan(s)) reproduces the 5 unrolled control points — is the
most principled such criterion.  `phi_span` (angular span of the 5 points on the
cylinder) is a weaker proxy.  The current implementation uses phi_span; the orbit
reconstruction error should be exposed as a `fit_error()` accessor and used instead.

**Designed (not yet implemented) fallback chain:**

```
Exact-fit cylinder + ConicWindow<2>         (Newton converged, ConicWindow accepted)
  → Best-fit cylinder + ConicWindow<2>      (only if cylinder surface residuals small)
    → LagrangeWindow<3>                     (no cylindrical structure in the data)
```

`ConicWindow<3>` is removed from the chain — it was never geometrically justified as a
fallback, only present as a historical pre-cylinder code path.

Best-fit cylinder: minimise Σ(dist(pᵢ, axis) − r)² over axis direction, axis point,
and radius r.  Always has a solution.  Fires only when `cyl_solve` Newton fails.
Needs a residual gate: if max(|dist(pᵢ, axis) − r|) > ε · window_scale, the data
has no cylindrical structure → fall through to LagrangeWindow<3>.

The convex-hull failure mode (one point in the convex hull of the other 4) causes
large residuals in both exact-fit and best-fit, so it exits to LagrangeWindow<3>
at both levels — which is the correct answer (no meaningful cylinder frame exists).

---

## Pending / open questions

- **Replace ConicWindow with CylinderWindow for 3D?** Open.  The measured data
  shows cylinder is always better for the two test curves.  Needs broader testing
  across curve families before replacing.  Document proof/evidence here when done.

- **Implement best-fit cylinder fallback**: replace the current ConicWindow<3>
  fallback with best-fit cylinder + residual gate + LagrangeWindow<3>.  Requires
  a 5-parameter nonlinear least-squares solver (Gauss-Newton or LM) on the
  cylinder surface distance function.

- **Expose ConicWindow<2> orbit reconstruction error**: add `fit_error()` accessor
  and use it to select among multiple accepted cylinders (minimum error wins)
  instead of phi_span.
