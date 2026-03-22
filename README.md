# conicblend

C^N-continuous parametric curve interpolation using local conic-section windows.
Header-only C++17, no dependencies, works in any dimension ≥ 2.

---

## What it does

Given a sequence of control points and their parameter values, `conicblend` builds a
smooth curve that passes **exactly** through every control point.  Smoothness is C^N —
position and the first N derivatives are continuous everywhere, not just at knots.

Each local segment is built from a best-fit conic arc (ellipse, parabola, or hyperbola)
through 5 consecutive control points.  Adjacent windows overlap and are blended with a
smoothstep weight, giving exact C^N continuity without solving any global system.

**Key properties:**
- **Exact interpolation** — curve passes through all interior control points
- **C^N continuity** — exact, algebraic, no tolerance (default N=2, quintic blend)
- **Affine invariant** — stretching or shearing the control points gives the same curve,
  stretched (~2e-13 error); similarity-invariant for circle windows
- **Exact conic reproduction** — if input lies on a conic, output is machine-epsilon exact
- **nD support** — works in ℝ² through ℝⁿ via `template<int Dim>`
- **Header-only C++17** — single `#include`, no build step, no dependencies
- **Local support** — changing one control point affects at most 4 consecutive segments

---

## Quick start

```cpp
#include "conicblend.hpp"   // or conicblend_circle.hpp for 3-point circle windows
#include <vector>
#include <cmath>

int main()
{
    // Sample a helix: P(t) = (cos t, sin t, 0.3 t)
    int n = 20;
    std::vector<fc::Vec3>   ctrl(n);
    std::vector<double>     times(n);
    for (int i = 0; i < n; ++i) {
        double t  = 4.0 * M_PI * i / (n - 1);
        ctrl[i]   = fc::Vec3(std::cos(t), std::sin(t), 0.3 * t);
        times[i]  = t;
    }

    // 3-point circle windows (conicblend_circle.hpp, minimum n=4)
    auto r_circle = fc::blend_curve(ctrl, times, 80, 2);

    // 5-point conic windows (conicblend.hpp, minimum n=6, affine-invariant)
    auto r_conic  = fc::blend_curve(ctrl, times, fc::conic_tag{}, 80, 2);

    // r.pts  — std::vector<Vec3>, dense output points
    // r.times — corresponding parameter values
}
```

Compile:
```bash
clang++ -std=c++17 -O2 -I. -o my_program my_program.cpp
```

---

## Two headers

| Header | Window | Points | Min n | Exact for | Invariance |
|---|---|---|---|---|---|
| `conicblend_circle.hpp` | Circle arc | 3 | 4 | Circles | Similarity |
| `conicblend.hpp` | Conic arc | 5 | 6 | All conics | **Affine** |

`conicblend.hpp` automatically includes `conicblend_circle.hpp` and falls back to
circle windows if any conic window fails.

### Tag dispatch

```cpp
auto r = fc::blend_curve(ctrl, times);                    // circle (3D, no tag)
auto r = fc::blend_curve(ctrl, times, fc::circle_tag{}); // circle (explicit)
auto r = fc::blend_curve(ctrl, times, fc::conic_tag{});  // conic (requires conicblend.hpp)
```

### nD template form

```cpp
auto r = fc::blend_curve<4>(ctrl4d, times, fc::conic_tag{});  // 4D
auto r = fc::blend_curve<2>(ctrl2d, times);                   // 2D
```

---

## Build and test

```bash
# Run tests directly
clang++ -std=c++17 -O2 -I. -o demo      demo.cpp      && ./demo
clang++ -std=c++17 -O2 -I. -o demo_nd   demo_nd.cpp   && ./demo_nd
clang++ -std=c++17 -O2 -I. -o demo_conic demo_conic.cpp && ./demo_conic

# Or via CMake
cmake -B build && cmake --build build
ctest --test-dir build
```

Tests cover: exact interpolation, C^N continuity, similarity and affine invariance,
collinearity guards, near-collinear fallback, nD (Dim=2/3/4), input validation.

---

## Documentation

Full documentation — architecture, mathematical properties, invariance proofs, API
reference, comparison with NURBS — is in [`conicblend_doc.md`](conicblend_doc.md).

---

## License

MIT — see [LICENSE](LICENSE).
