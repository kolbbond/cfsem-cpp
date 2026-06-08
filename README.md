# cfsem-cpp

A **C-ABI shim library** over the [`cfsem`](https://github.com/kolbbond/cfsem-py)
Rust crate (quasi-steady electromagnetics: filamentized Biot-Savart,
Grad-Shafranov, etc.). It wraps `cfsem`'s Rust functions in `extern "C"` so C and
C++ code — eventually the `goose` libraries — can link directly against it.

Built with [Corrosion](https://github.com/corrosion-rs/corrosion), which compiles
the Rust crates as part of the normal CMake build.

## Role and scope

cfsem is intended as an **alternative / reference B- and A-field calculator** for
goose — an *exact direct-summation* Biot-Savart backend to cross-check goose's
fast-multipole (rat MLFMM) results, and for small problems.

Two things follow from that:

- It is **O(N·M)** (every observation point against every filament), not an FMM,
  so it is a validation/reference path, not a drop-in performance replacement.
- cfsem and rat regularize the near-filament singularity **differently** (cfsem
  models a finite-thickness conductor via `wire_radius`; rat van Lanen uses an
  `eps` softening). They agree in the **far field** (observation well separated
  from sources) but **not** for self-field / on-filament evaluation. Validate on
  separated source/observation configurations.

## Data layout

Coordinate arrays use the **interleaved 3×N** layout `[x0,y0,z0, x1,y1,z1, …]`,
matching an Armadillo `arma::Mat<double>(3, N)` in memory (column-major). goose
can pass `mat.memptr()` directly — no transpose at the boundary.

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
./cfsem-cpp        # smoke test: prints ellipk + a Biot-Savart B-field
```

The first build compiles `cfsem` and its dependencies (`faer`, `nalgebra`,
`rayon`), so it takes a couple of minutes; later builds are incremental.

## Linking from another project (e.g. goose)

This repo exposes a CMake `INTERFACE` target, `cfsem_cpp`, that bundles the Rust
static library and the headers. Add this repo (submodule or `add_subdirectory`)
and link it:

```cmake
add_subdirectory(cfsem-cpp)
target_link_libraries(goose PRIVATE cfsem_cpp)   # gets libcfsem_capi.a + include/
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
| `cfsem_ellipk(m)` | complete elliptic integral K(m) (scalar) |
| `cfsem_flux_density_linear_filament(...)` | B-field (T) from linear filaments |
| `cfsem_vector_potential_linear_filament(...)` | A-field (V·s/m) from linear filaments |

Status codes: `0` ok · `1` cfsem error (e.g. length mismatch) · `2` null pointer
· `3` panic caught at the FFI boundary.

## Adding a function to the C API

Exposing another `cfsem` function takes three coordinated edits:

1. **`cfsem-capi/src/lib.rs`** — add a `#[no_mangle] pub extern "C"` wrapper,
   wrapped in `catch_unwind`. Array functions follow the interleaved-3×N +
   element-count convention used by `cfsem_flux_density_linear_filament` (the
   shared `linear_filament_field` helper de-interleaves into the structure-of-
   arrays form cfsem wants); output buffers are caller-allocated, return is an
   `int` status code.
2. **`include/cfsem.h`** — declare the matching C signature (signatures are
   unchecked across the FFI boundary; keep them exactly in sync).
3. **`include/cfsem.hpp`** — add the `namespace cfsem` inline forwarder.
