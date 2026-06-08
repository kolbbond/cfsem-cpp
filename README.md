# cfsem-cpp

A **C-ABI shim library** over the [`cfsem`](https://github.com/kolbbond/cfsem-py)
Rust crate (quasi-steady electromagnetics: filamentized Biot-Savart,
Grad-Shafranov, etc.). It wraps `cfsem`'s Rust functions in `extern "C"` so C and
C++ code — eventually the `goose` libraries — can link directly against it.

Built with [Corrosion](https://github.com/corrosion-rs/corrosion), which compiles
the Rust crates as part of the normal CMake build.

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

Then in C++:

```cpp
#include "cfsem.hpp"
double k = cfsem::ellipk(0.1);
```

…or from C via `#include "cfsem.h"` and the `cfsem_*` functions.

## Adding a function to the C API

Exposing another `cfsem` function takes three coordinated edits:

1. **`cfsem-capi/src/lib.rs`** — add a `#[no_mangle] pub extern "C"` wrapper that
   rebuilds the Rust arguments from raw pointers and calls the `cfsem` function.
   Array arguments follow the pointer-plus-length convention used by
   `cfsem_flux_density_linear_filament` (component arrays + element counts, output
   buffers caller-allocated, `int` status return).
2. **`include/cfsem.h`** — declare the matching C signature (signatures are
   unchecked across the FFI boundary; keep them exactly in sync).
3. **`include/cfsem.hpp`** — add the `namespace cfsem` inline forwarder.
