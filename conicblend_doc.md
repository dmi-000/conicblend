# conicblend — C^N 3D Curve Interpolation from Circle Arcs

## What it does

Given a sequence of **n control points** in ℝ³ and their **parameter values** (times),
`conicblend` constructs a smooth 3D parametric curve that passes exactly through all
of them.  Smoothness is C^N — position, velocity, acceleration, and higher derivatives
are continuous everywhere, not just at the control points.

Each local segment is built from **3-point circle arcs**: any three non-collinear
points in ℝ³ define a unique circumscribed circle lying in the plane of those points.
Adjacent windows overlap and are blended with a smoothstep weight, giving exact
continuity without solving a global linear system.

Torsion (out-of-plane twist) is handled **implicitly**: adjacent windows live in
slightly different planes, so the blended curve naturally follows the Frenet binormal
rotation without any explicit torsion estimation.

---

## Core theory: overlapping 3-point windows + smoothstep blend

Every **interior segment** `ctrl[j] → ctrl[j+1]` is covered by two overlapping
**3-point window functions**:

```
Window A (index j−1):  ctrl[j−1], ctrl[j],   ctrl[j+1]
Window B (index j)  :  ctrl[j],   ctrl[j+1],  ctrl[j+2]
```

Each window is a circle arc through its three control points.  The blended segment is:

```
curve(t) = (1 − w(s)) · A(t)  +  w(s) · B(t)
```

where `s = (t − t_j) / (t_{j+1} − t_j) ∈ [0,1]` and `w = smoothstep(s, N)`.

### Why this gives exact C^N continuity

At knot `ctrl[j]` (`s = 0`), `w = 0` and `dᵏw/dsᵏ = 0` for `k = 1 … N`, so
the blended curve reduces to `A(t)` on both sides of the junction — the **same
function instance** approaches from both the ending segment `j−1` and the starting
segment `j`.  No matching condition is required: continuity is an algebraic identity.

Concretely: segment `j−1` ends at `s=1` where `w=1` → blend = win[j-1](t_{j}`) = `ctrl[j]`;
segment `j` starts at `s=0` where `w=0` → blend = `win[j-1](t_j)` = `ctrl[j]`.
Both sides evaluate the **same window** at the **same parameter** → identical
position, velocity, and N higher derivatives.

### Smoothstep

`smoothstep(s, N)` is the unique polynomial satisfying:
- `w(0) = 0`, `w(1) = 1`
- `dᵏw/dsᵏ = 0` at `s = 0` and `s = 1` for `k = 1 … N`

| N | Degree | Continuity | Formula |
|---|--------|------------|---------|
| 1 | 3      | C¹         | 3s² − 2s³ |
| 2 | 5      | C²         | 10s³ − 15s⁴ + 6s⁵ |
| 3 | 7      | C³         | 35s⁴ − 84s⁵ + 70s⁶ − 20s⁷ |

Default is N=2 (C², quintic weight).

---

## Circle window construction

### Circumscribed circle (`fc::circumcircle`)

Three non-collinear 3D points `p0, p1, p2` define a unique circle:

```
a = p1 − p0,   b = p2 − p0
A = a·a,  B = a·b,  C = b·b,  D = AC − B²   (= |a×b|²  by Lagrange identity)
s = C(A−B)/(2D),  t = A(C−B)/(2D)
center = p0 + s·a + t·b
```

The in-plane basis is `e1 = a/|a|`, normal `= (a×b)/|a×b|`, `e2 = normal × e1`.
A point on the circle: `center + R·cos(φ)·e1 + R·sin(φ)·e2`.

**Collinearity guard**: `D = |a×b|² < ε·A·C` signals near-collinear points.
The threshold `ε = 1e-14` gives ~1e-7 relative tolerance on the radius computation.
A `std::invalid_argument` is thrown; the blend pipeline interprets this as a signal
to request more control points.

### Angle unwrapping (`fc::CircleWindow`)

Each control point maps to an angle: `φᵢ = atan2((pᵢ − center)·e2, (pᵢ − center)·e1)`.

The raw angles `φ0, φ1, φ2` are unwrapped so the sequence is monotone:
- `d1 = wrap(φ1 − φ0)`, `d2 = wrap(φ2 − φ0)`  where `wrap` maps to (−π, π]
- If `d1` and `d2` have opposite signs, the arc crosses the `±π` branch cut:
  adjust `d2 += ±2π` to bring it to the same side as `d1`.
- Store `phi1 = phi0 + d1`, `phi2 = phi0 + d2`.

**Angular interpolation**: `φ(t)` is the unique quadratic Lagrange polynomial through
`(t0, phi0), (t1, phi1), (t2, phi2)`.  This is smooth, exactly interpolating, and
uses only elementary arithmetic — no transcendental evaluation at each t.
Degree 2 (three knots), vs. degree 4 for the conic-window cross-ratio parametrisation.

**Window evaluation**: `CircleWindow(t) = center + R·cos(φ(t))·e1 + R·sin(φ(t))·e2`.

**Invariance**: `φ = atan2(...)` is a Euclidean metric quantity.
- *Similarity-invariant* (rotation, uniform scale, translation): similarities map circles
  to circles; the new angles satisfy `φ_i' = φ_i + rotation_offset` (same constant for
  all i), so `φ'(t) = φ(t) + const` and the arc evaluates to the transformed point.
- *Not affinely invariant*: a non-uniform stretch maps the circle to an ellipse, but
  the window fits a new circumscribed circle to the stretched points — a geometrically
  different object.  The angle parametrisation is defined from the new center, with no
  algebraic relationship to the original angles.  Achieving affine invariance for
  3-point windows would require fitting the optimal local conic (ellipse) instead of a
  circle, which is exactly the 5-point conic window.

**Near-collinear fallback**: when `D = AC − B² < 1e-10·AC` (radius > ~10⁵ × chord),
`CircleWindow` silently switches to **quadratic Lagrange position interpolation** in
place of the circle arc.  No exception is thrown; C^N continuity is unaffected because
the proof depends only on the same window evaluating from both sides of each knot,
regardless of window type.  This is a local fallback (per window); by contrast,
`ConicWindow` failure causes the entire `blend_curve` call to fall back to circle windows.

---

## Conic window construction (`conicblend.hpp`)

For a 5-point window `ctrl[i..i+4]`, the conic is constructed as follows.

### Best-fit 2D plane

For nD ≥ 3, the 5 points are projected to their best-fit 2D plane via PCA
(covariance matrix eigendecomposition).  For Dim=2, this step is a no-op.

### SVD conic fit

The 5×6 design matrix `M[k] = [x²,x·y,y²,x,y,1]` is formed from the 2D projected
coordinates (centered for numerical stability).  The conic coefficients `[A,B,C,D,E,F]`
are the smallest eigenvector of `M^T·M` — the null-space of M.

### Cross-ratio parametrisation (affinely invariant)

The conic `Ax²+Bxy+Cy²+Dx+Ey+F=0` is parametrised via the rational stereographic map
from the base point P₀:

```
slope s  →  u(s) = -(L_e + M_e·s) / (A_e + C_e·s²)
             x = x₀ + u,   y = y₀ + s·u
```

where `A_e, C_e` are eigenvalues of the quadratic form matrix (principal frame),
`L_e = 2A_e·x₀ + D_e`, `M_e = 2C_e·y₀ + E_e`.

**Cross-ratio invariant** — Instead of the direct stereographic angle `φ = 2·arctan(s)`
(which is NOT invariant under non-orthogonal transforms), the 5 slopes
`s₀, s₁, s₂, s₃, s₄` from P₀ to each control point are encoded as cross-ratio values:

```
q_i = (s_i − s₀)(s₂ − s₄) / [(s_i − s₄)(s₂ − s₀)]
```

where `s₀, s₂, s₄` are the slopes at the first, middle, and last control points,
used as anchor values.  The `q_i` are then normalised to `r_i = q_i / (q_i + 1) ∈ [0,1]`
so that all 5 values are finite (`r₀ = 0`, `r₄ = 1` always).

A degree-4 Lagrange polynomial (C^∞) interpolates `r(t)` at the 5 knots.
At evaluation: `r → q = r/(1−r) → s` via the Möbius inversion:

```
s = [q·s₄·(s₂−s₀) − s₀·(s₂−s₄)] / [q·(s₂−s₀) − (s₂−s₄)]
```

then the rational map evaluates `x, y` from `s`.

**Why this achieves exact affine invariance**: under any affine transform
`T = (x,y) → (ax+by, cx+dy)`, all slopes from P₀ undergo the _same_ Möbius map
`s → (c+ds)/(a+bs)`.  The cross-ratio `q_i` is invariant under any Möbius transform
of the slopes — so `q_i^new = q_i^old`, `r_i^new = r_i^old`, and the Lagrange
polynomial `r(t)` is identical before and after any affine transform.
The inversion `r → q → s` then reconstructs the correct (Möbius-transformed)
slopes for the new conic.  Error: floating-point noise only (~2e-13).

### Why P₀ selection gives machine-epsilon blends

P₀ is chosen from among the 5 control points (trying order `{0, 4, 1, 3, 2}`) as
the first candidate for which the cross-ratio values `r₀, …, r₄` are monotone.
All 5-point windows fitting the same physical conic compute the same cross-ratio
sequence → the same Lagrange polynomial → identical orbit functions → the blend
stays exactly on the conic for locally-conic input data.

### Fallback policy

If any step fails (degenerate anchor slopes, non-monotone cross-ratio sequence,
orbit reconstruction error > 1e-3), `ConicWindow::valid()` returns false and the
entire `blend_curve` call falls back to circle windows.

---

## API

Two headers; both live in namespace `fc`, both are single-include C++17.

### Circle windows (`conicblend_circle.hpp`)

#### Types

```cpp
fc::Vec3            // 3D vector (= VecN<3>); access via [0],[1],[2]
fc::Circle3         // circumscribed circle: {center, radius, e1, e2}
fc::BlendResult     // {std::vector<Vec3> pts, std::vector<double> times}

// nD variants (template parameter Dim ≥ 2):
fc::VecN<Dim>
fc::CircleND<Dim>
fc::BlendResultND<Dim>
```

#### Functions

```cpp
// Fit circumscribed circle to three non-collinear points (any Dim ≥ 2).
// Throws std::invalid_argument if points are (near-)collinear.
fc::CircleND<Dim> fc::circumcircle<Dim>(VecN<Dim> p0, VecN<Dim> p1, VecN<Dim> p2);

// Smoothstep blend weight. N in {1,2,3}. Clamps s to [0,1].
double fc::smoothstep(double s, int N = 2);

// Build and evaluate the blended curve.
// ctrl.size() == times.size() >= 4; times strictly increasing.
// Returns pts_per_seg points per interior segment; final point included once.
// Blended segments j = 1…n-3; curve passes through ctrl[1]…ctrl[n-2].
fc::BlendResult fc::blend_curve(                        // 3D, no template syntax
    std::vector<Vec3>   const& ctrl,
    std::vector<double> const& times,
    int pts_per_seg = 60, int smooth_N = 2);

fc::BlendResultND<Dim> fc::blend_curve<Dim>(            // nD template
    std::vector<VecN<Dim>> const& ctrl,
    std::vector<double>    const& times,
    int pts_per_seg = 60, int smooth_N = 2);
```

### Conic windows (`conicblend.hpp`)

Requires `conicblend_circle.hpp` (included automatically).

#### Types

```cpp
fc::ConicWindow<Dim>   // 5-point conic window; call valid() before use
```

#### Function

```cpp
// Build and evaluate the blended curve using 5-point conic windows.
// ctrl.size() == times.size() >= 6; times strictly increasing.
// Blended segments j = 2…n-4; curve passes through ctrl[2]…ctrl[n-3].
// If any window fails (non-monotone phi, bad fit), falls back to circle windows.
fc::BlendResultND<Dim> fc::blend_curve<Dim>(
    std::vector<VecN<Dim>> const& ctrl,
    std::vector<double>    const& times,
    fc::conic_tag,
    int pts_per_seg = 60, int smooth_N = 2);

fc::BlendResult fc::blend_curve(                        // 3D convenience overload
    std::vector<Vec3>   const& ctrl,
    std::vector<double> const& times,
    fc::conic_tag,
    int pts_per_seg = 60, int smooth_N = 2);
```

#### Window-type tag dispatch

A tag argument selects the window type at compile time with zero runtime overhead:

```cpp
// Circle windows (3-point, minimum n=4):
auto r = fc::blend_curve(ctrl, times);                    // untagged 3D
auto r = fc::blend_curve(ctrl, times, fc::circle_tag{}); // tagged

// Conic windows (5-point, minimum n=6) — requires conicblend.hpp:
#include "conicblend.hpp"
auto r = fc::blend_curve(ctrl, times, fc::conic_tag{});

// nD template forms:
auto r = fc::blend_curve<4>(ctrl, times, fc::circle_tag{});
auto r = fc::blend_curve<4>(ctrl, times, fc::conic_tag{});
```

Calling `blend_curve` with `fc::conic_tag{}` without including `conicblend.hpp`
gives a clear compile error — no silent fallback.

#### Embedded / no-exceptions builds

Define `FC_NO_EXCEPTIONS` before including either header to replace all
`std::invalid_argument` throws with `std::terminate()`:

```cpp
#define FC_NO_EXCEPTIONS
#include "conicblend_circle.hpp"   // or conicblend.hpp
```

#### Build

Single-header C++17, no dependencies:

```bash
# Direct compilation (Homebrew clang++ on macOS, or g++ on Linux)
clang++ -std=c++17 -O2 -I. -o demo demo.cpp
clang++ -std=c++17 -O2 -I. -o demo_conic demo_conic.cpp

# CMake
cmake -B build && cmake --build build
# Targets: demo, demo_nd, demo_conic
```

### N-dimensional generalisation

Both headers are fully nD via `template<int Dim>` (Dim ≥ 2):

**Circle windows** use Gram-Schmidt instead of the cross product to complete the in-plane basis:
```
e2 = (b − (b·e1)·e1) / |b − (b·e1)·e1|
```
The circumcenter formula (Cramer's rule on scalar quantities A, B, C, D) is purely in terms of dot products — dimension-independent.

**Conic windows** use a Dim×Dim covariance eigendecomposition to find the best-fit 2D plane, then project to 2D for the conic fit.

```cpp
#include "conicblend.hpp"

using V4 = fc::VecN<4>;

// 4D ellipse embedded in ℝ⁴
std::vector<V4>     ctrl(n);
std::vector<double> times(n);
// ... fill ctrl and times ...

auto r_circle = fc::blend_curve<4>(ctrl, times, fc::circle_tag{});
auto r_conic  = fc::blend_curve<4>(ctrl, times, fc::conic_tag{});
```

`Dim=2` reproduces 2D arcs; `Dim=3` matches the 3D convenience overloads to ~1e-15.

### Python prototype (`frenet_blend_proto.py`, `demo_curves.py`)

The Python prototype replicates the C++ architecture exactly and adds:
- `blend_curve_frenet(ctrl, times, ...)`: Frenet-space local blend (frozen frame per segment, approximate)
- `estimate_kappa_tau(ctrl, j)`: 5-point Menger curvature κ and torsion τ estimator
- `smoothblend_tau(ctrl, times, ...)`: τ estimated at ctrl pts, smoothblended between knots
- `demo_curves.py`: six 3D test curves with window-arc visualisation

---

## Typical usage (C++)

```cpp
// Circle windows (3-point, minimum n=4):
#include "conicblend_circle.hpp"
#include <vector>
#include <cmath>

int main()
{
    // Sample a helix: P(t) = (cos t, sin t, 0.3 t)
    int n = 20;
    std::vector<fc::Vec3>   ctrl(n);
    std::vector<double>     times(n);
    for (int i = 0; i < n; ++i) {
        double t  = 4.0 * fc::detail::PI * i / (n - 1);
        ctrl[i]   = fc::Vec3(std::cos(t), std::sin(t), 0.3 * t);
        times[i]  = t;
    }

    auto result = fc::blend_curve(ctrl, times, /*pts_per_seg=*/80, /*N=*/2);

    for (std::size_t i = 0; i < result.pts.size(); ++i) {
        auto& p = result.pts[i];
        printf("t=%.4f  x=%.6f  y=%.6f  z=%.6f\n",
               result.times[i], p[0], p[1], p[2]);
    }
}
```

```cpp
// Conic windows (5-point, minimum n=6, exact on conics):
#include "conicblend.hpp"

    auto result = fc::blend_curve(ctrl, times, fc::conic_tag{}, 80, 2);
    // Falls back to circle windows if any window fails
```

The result contains `(n−3) × pts_per_seg + 1` points covering the interior segments
`j = 1, …, n−3`.  The curve passes exactly through `ctrl[1], …, ctrl[n−2]`.

---

## Typical usage (Python)

```python
import numpy as np
from demo_curves import CircleWindow, blend_curve, smoothblend_tau

# Sample a torus knot
n    = 30
ts   = np.linspace(0, 2*np.pi, n)
ctrl = np.c_[(2 + np.cos(3*ts))*np.cos(2*ts),
             (2 + np.cos(3*ts))*np.sin(2*ts),
             np.sin(3*ts)]

pts, t_out = blend_curve(ctrl, ts, pts_per_seg=80)
tau_b, tau_c = smoothblend_tau(ctrl, ts, pts_per_seg=80)

# pts: (n_interior_segs * 80 + 1, 3) array of blended positions
# tau_b: smoothblended torsion estimate at each output point
```

---

## Mathematical properties

### Exact interpolation

The blended curve passes through every control point exactly.  At `t = times[j]`
(`s = 0`), `w = 0` and the blend equals `win[j-1](times[j]) = ctrl[j]`.
At `t = times[j+1]` (`s = 1`), `w = 1` and the blend equals `win[j](times[j+1]) = ctrl[j+1]`.

### C^N continuity (exact, not approximate)

Position and the first N derivatives are continuous everywhere.  At junction `times[j]`:
- Segment `j−1` at `s=1`: blend = `win[j-1](times[j])` = `ctrl[j]`
- Segment `j`   at `s=0`: blend = `win[j-1](times[j])` = `ctrl[j]`

Identical calls — same window, same parameter, same result.  No tolerance, no matching.

### Exact reproduction of circles and planar arcs

If all control points lie on a circle, every window computes the **same circle**.
The blend `(1−w)·A(t) + w·B(t)` lies on the chord between two points on the circle,
which is not generally on the circle — *except* when `A(t) = B(t)`.  For a uniform
parameter spacing on a circle, both windows evaluate to identical points at all `t`,
so the blend reproduces the circle to floating-point precision (error < 1e-15).

Circle windows do **not** exactly reproduce ellipses, parabolas, or hyperbolas — they
approximate them with the best locally-fitting circle arc (O(h⁴) accuracy per segment).
Conic windows (`conic_tag`) exactly reproduce all conics, including circles.

### Implicit torsion

Each `CircleWindow` lives in the plane of its three control points.  Adjacent
windows (A and B for segment j) live in different planes.  The blend
`(1−w)·A(t) + w·B(t)` traces a path that exits the plane of A and enters the plane
of B — this is the torsion of the blended curve.  No explicit binormal integration
is required: the 3D geometry of the overlapping windows encodes it.

This is equivalent to saying that the blended curve in Frenet decomposition
`P = u·T + v·N + b·B` has its B-component (torsion accumulation) carried implicitly
by the global positions of the window circle arcs, not by a separate scalar field.
See `frenet_blend_proto.py` for proof that "global Frenet-space blending" is
algebraically identical to position-space blending.

### Local support

Modifying `ctrl[k]` affects at most **three windows** (those covering it:
`wins[k-2], wins[k-1], wins[k]`) and at most **four consecutive blend segments**.
Influence radius is O(1/n) of the total parameter range — far smaller than global
splines which shift the entire curve.

### No global linear system

Each window is fitted from three points by direct closed-form algebra (Cramer's rule
for the circumcenter).  There is no n×n system.  A degenerate region (near-collinear
triplet) falls back to a linear interpolant locally without affecting the rest of the
curve.

### Invariance properties

"T-invariant" means: if `ctrl_new = T(ctrl)`, then `result_new.pts[i] = T(result_orig.pts[i])`
to the listed precision.

#### Circle windows (`circle_tag`)

| Transform class | Error bound | Why |
|---|---|---|
| Rigid / similarity (rotation, uniform scale, translation) | < 2e-14 | Similarities map circles to circles; angles shift by a constant, so φ'(t) = φ(t) + const |
| General affine (non-uniform scale, shear) | O(stretch factor) | Affine maps a circle to an ellipse; the window re-fits a circle to the stretched points — a different object |

The angle `φ = atan2(...)` is defined by the Euclidean metric (dot products with `e1, e2`),
so only metric-preserving transforms (similarities) can be invariant.

#### Conic windows (`conic_tag`)

| Transform class | Error bound | Mechanism |
|---|---|---|
| Rigid / similarity (rotation, uniform scale, translation) | < 2e-14 | Orthogonal frame — principal-frame rotation is equivariant |
| General affine (non-uniform scale, shear, reflection) | < 3e-13 | Cross-ratio of slopes is Möbius-invariant; affine transforms Möbius-transform all slopes by the same map |
| Projective (perspective with ε depth gradient) | O(ε · scale) ≈ 8e-9 at ε=10% | Projective ≈ affine locally; residual from inter-point denominator variation |

**Affine invariance — proof sketch**: Under affine `T = (x,y) → (ax+by, cx+dy)`, the
stereographic slope from P₀ to any other conic point transforms as the Möbius map
`f: s → (c+ds)/(a+bs)`.  All 5 slopes undergo the _same_ `f`, so the cross-ratio
`q_i = CR(s₀, s₂; s_i, s₄)` is preserved exactly (cross-ratio is Möbius-invariant).
The normalized values `r_i = q_i/(q_i+1)` are preserved → the Lagrange polynomial
`r(t)` is identical → the inversion `r → q → s` restores the Möbius-transformed
slopes → the rational map gives the affinely-transformed orbit point.

**What true projective invariance would require**: the slopes s_i are projective
parameters on the conic (via the rational parametrization from P₀).  Under a
projective transform, the new slope from T(P₀) to T(Pᵢ) is NOT a fixed Möbius map
of the old slope (the denominator `gx+hy+k` differs per point).  Full projective
invariance would require computing the projective parameter in homogeneous coordinates
throughout — a significant redesign with no practical benefit, given the O(ε)
sub-nanometre error already achieved for typical perspective strengths.

### Convergence

For a smooth curve in ℝ³, the approximation error per segment is O(h⁴) where h is
the control-point spacing, matching the order of natural cubic splines.  The helix
at n=20 achieves max error ~5e-3; doubling n roughly halves the error at low n and
gives O(h²) convergence once the angular steps are small.

---

## Design principles

- **No global smoothing**: the curve passes through every control point exactly.
- **No matching conditions**: C^N continuity follows algebraically from the
  overlapping-window construction, not from solving a system of equations.
- **Exact 3D, no projection**: three non-collinear points in ℝ³ always define a
  unique circle in their shared plane.  No SVD, no dimensionality reduction needed.
- **Torsion implicit**: the binormal rotation between adjacent window planes naturally
  encodes torsion.  No explicit Frenet integration required.
- **Collinearity = signal**: a near-collinear triplet means the local curve is nearly
  straight and the circumcircle radius diverges.  This is geometrically correct and
  the library signals it cleanly (exception in C++; `LinearWindow` fallback in Python).
- **Mirror of conicspline in 3D**: `conicblend` replaces 5-point conic windows with
  3-point circle windows, reducing the window size by 2 (less context, simpler fit)
  while gaining native 3D support.  The smoothstep blend and C^N proof are identical.

---

## Comparison with conicspline

| Property | conicspline (2D) | conicblend circle | conicblend conic |
|---|---|---|---|
| Window size | 5 points | 3 points | 5 points |
| Window type | Conic arc | Circle arc | Conic arc |
| Fitting | SVD null-space (5×6) | Cramer's rule (2×2) | SVD null-space (5×6) + best-fit plane |
| Orbit parameter | phi = 2·arctan(s) | angle via atan2 | cross-ratio r_i = q_i/(q_i+1) |
| Invariance | similarity | similarity | **affine** (exact); projective (~8e-9 at 10%) |
| Fallback | Natural cubic spline | exception (C++) | circle windows |
| Torsion | N/A (2D) | Implicit via tilting planes | Implicit via tilting planes |
| Exact reproduction | Conics | Circles and circular arcs | Conics (ellipse, hyperbola, parabola) |
| Ambient dimension | 2D only | nD (any Dim ≥ 2) | nD (any Dim ≥ 2) |
| Minimum n | 6 | 4 | 6 |
| Language | Python | C++17 header-only | C++17 header-only |

---

## File layout

| File | Role |
|---|---|
| `conicblend_circle.hpp` | 3-pt circle windows: namespace `fc`, `template<int Dim>`, `VecN`/`CircleND`/`CircleWindow`/`blend_curve`; 3D aliases; `circle_tag`/`conic_tag` stub; `FC_NO_EXCEPTIONS` |
| `conicblend.hpp` | 5-pt conic windows: `ConicWindow<Dim>`, cross-ratio orbit, `SymEig<N>` Jacobi SVD, `Pchip5`, `best_fit_plane<Dim>`; `blend_curve(..., conic_tag{})` |
| `demo.cpp` | 13 tests (3D circle): helix, torus knot, C^N, collinearity, exact circle, input validation, n=4, large arc, clockwise, tilted circle, non-uniform times, N=1/N=3, tag dispatch |
| `demo_nd.cpp` | 6 tests (nD circle): Dim=2 circle, Dim=3 helix, Dim=4 torus, collinearity, near-collinear fallback, input validation |
| `demo_conic.cpp` | 9 tests (5-pt conic): tilted ellipse (~2e-14), parabola (~3.5e-15), 3D helix (~1.3e-15), 4D ellipse (~2.3e-14), input validation, used_conic_flag, C^N continuity, similarity invariance (~2e-14), **affine invariance (~2e-13)** |
| `CMakeLists.txt` | CMake build: targets `demo`, `demo_nd`, `demo_conic` |
| `frenet_blend_proto.py` | Python prototype: position-space and Frenet-space blend comparison |
| `demo_curves.py` | Six 3D test curves with window-arc visualisation and τ coloring |
| `plot_blend.py` | Visualisation of helix/torus knot τ: smoothblended, convergence, analytic |
