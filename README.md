# cfsem-cpp

A **C-ABI shim library** over the [`cfsem`](https://github.com/cfs-energy/cfsem-py)
Rust crate, so C and C++ code can link directly against it.

Built with [Corrosion](https://github.com/corrosion-rs/corrosion), which compiles
the Rust crates as part of the normal CMake build.

## About cfsem

[`cfsem`](https://docs.rs/cfsem) is a library of quasi-steady electromagnetics
formulas: Biot-Savart fields, vector potentials, and inductances for linear and
circular current filaments, Grad-Shafranov operators, volumetric and
boundary-element sources, and related math (e.g. elliptic integrals). Its stated
goals are library-level formulas with literature references and units,
quantitative unit tests, and performance; the Rust core is parallelized with
`rayon` and also ships as a Python package.

For linear filaments it offers two evaluation paths, both exposed here:

- **Direct summation** — every observation point against every filament,
  O(N·M). Exact up to floating point, and the reference the other path is
  checked against.
- **Hierarchical (Barnes-Hut)** — a source tree with an acceptance angle
  `theta`, trading a controllable approximation error for scaling to large
  problems. There is no guaranteed error bound, so validate `theta` against the
  direct functions for your geometry.

Near a filament, cfsem treats each source as a finite-thickness conductor of
radius `wire_radius` rather than an ideal line current, so fields stay bounded
on and inside the conductor. Codes that regularize this singularity differently
will agree with cfsem in the far field but not for self-field or on-filament
evaluation; compare on well-separated source/observation configurations.

This repo exposes only the linear-filament subset. The rest of cfsem is
reachable from Rust or Python directly.

## Data layout

Coordinate arrays use the **interleaved 3×N** layout `[x0,y0,z0, x1,y1,z1, …]`,
matching an Armadillo `arma::Mat<double>(3, N)` in memory (column-major). A
consumer can pass `mat.memptr()` directly — no transpose at the boundary.
cfsem works in structure-of-arrays internally, so each call de-interleaves into
temporaries (an O(N) copy).

## Layout

| Path          | What it is                                                              |
| ------------- | ----------------------------------------------------------------------- |
| `cfsem-py/`   | Rust crate `cfsem` — git submodule, upstream, **untouched**             |
| `cfsem-capi/` | Rust crate that exposes `#[no_mangle] extern "C"` wrappers → `libcfsem_capi.a` |
| `include/`    | `cfsem.h` (C ABI) and `cfsem.hpp` (`namespace cfsem` C++ wrapper)        |
| `main.cpp`    | smoke test (not the deliverable)                                        |

`cfsem-capi` depends on `cfsem` with **default features only** (no `python`), so
no pyo3 / Python runtime is linked.

## Requirements

- **Rust ≥ 1.85** — `cfsem` uses edition 2024. (`rustup update` if older.)
- **CMake ≥ 3.22** and a C++17 compiler.

## Install and build

```bash
# Clone with the cfsem-py submodule
git clone --recursive https://github.com/kolbbond/cfsem-cpp.git
# (existing clone) git submodule update --init

mkdir build && cd build
cmake ..
make
./cfsem-cpp        # smoke test; asserts its own tolerances
```

The first build compiles `cfsem` and its dependencies (`faer`, `nalgebra`,
`rayon`), so it takes a couple of minutes; later builds are incremental.

## Linking from another project

This repo exposes a CMake `INTERFACE` target, `cfsem_cpp`, that bundles the Rust
static library and the headers. Add this repo (submodule or `add_subdirectory`)
and link it:

```cmake
add_subdirectory(cfsem-cpp)
target_link_libraries(<consumer> PRIVATE cfsem_cpp)   # gets libcfsem_capi.a + include/
```

Then in C++ (interleaved 3×N arrays; see `cfsem.h` for the full contract):

```cpp
#include "cfsem.hpp"

double k = cfsem::ellipk(0.1);

double b_out[3];  // interleaved Bx,By,Bz at the single observation point
int rc = cfsem::flux_density_linear_filament(
    rs_obs, /*n_obs=*/1, rs_fil, drs_fil, ifil, wire_radius, /*n_fil=*/1, b_out);
// rc == 0 on success; cfsem::vector_potential_linear_filament has the same shape.
```

…or from C via `#include "cfsem.h"` and the `cfsem_*` functions.

## Current C API

| Function | Returns |
| -------- | ------- |
| `cfsem_ellipk(m)`, `cfsem_ellipe(m)` | complete elliptic integrals K(m), E(m) (scalar) |
| `cfsem_flux_density_linear_filament(...)` | B (T), direct sum |
| `cfsem_vector_potential_linear_filament(...)` | A (V·s/m), direct sum |
| `cfsem_{flux_density,vector_potential}_linear_filament_matrix(...)` | per-source B/A response matrices, row-major `(n_obs, n_fil)` |
| `cfsem_inductance_linear_filaments(...)` | mutual inductance (H) from all sources to each target segment |
| `cfsem_inductance_linear_filaments_matrix(...)` | dense mutual-inductance matrix (H), row-major `(n_src, n_tgt)` |
| `cfsem_body_force_density_linear_filament(...)` | J×B body force density (N/m³) |
| `cfsem_flux_density_linear_filament_hierarchical(..., theta, par, skip, ...)` | B (T), Barnes-Hut |
| `cfsem_vector_potential_linear_filament_hierarchical(..., theta, par, skip, ...)` | A (V·s/m), Barnes-Hut |
| `cfsem_near_field_interaction_map_linear_filament(...)` | CSC pattern of the near-field (direct) interactions |
| `cfsem_inductance_linear_filaments_sparse_csc(...)` | mutual inductance (H) on a caller-supplied CSC pattern |
| `cfsem_lf_tree_*` | persistent Barnes-Hut tree: build once, update currents, evaluate A/B |

The `*_hierarchical` variants add `theta` (acceptance angle; smaller = more
accurate, slower), `par` (nonzero = parallel), and `skip` (`CFSEM_SKIP_*`), which
omits near- or far-field interactions; a full solve equals near-only plus
far-only. Together with the interaction map and sparse inductance, this supports
a near/far split: an explicit sparse near-field M plus a far-field-only tree
solve. The `cfsem_lf_tree_*` handle avoids rebuilding the tree when only
currents change.

Status codes: `0` ok · `1` cfsem error (e.g. length mismatch) · `2` null pointer
· `3` panic caught at the FFI boundary · `4` buffer too small (required length
written back) · `5` argument out of range.

## Adding a function to the C API

Exposing another `cfsem` function takes coordinated edits:

1. **`cfsem-capi/src/lib.rs`** — add a `#[no_mangle] pub extern "C"` wrapper,
   wrapped in `catch_unwind`. Array functions follow the interleaved-3×N +
   element-count convention used by `cfsem_flux_density_linear_filament` (the
   shared `linear_filament_field` helper de-interleaves into the structure-of-
   arrays form cfsem wants); output buffers are caller-allocated, return is an
   `int` status code.
2. **`include/cfsem.h`** — declare the matching C signature (signatures are
   unchecked across the FFI boundary; keep them exactly in sync).
3. **`include/cfsem.hpp`** — add the `namespace cfsem` inline forwarder.
4. **`main.cpp`** — exercise it; an export with no smoke test doesn't exist.
