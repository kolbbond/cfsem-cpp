//! C-ABI shim over the `cfsem` crate.
//!
//! Every function here is `#[no_mangle] pub extern "C"` so a C or C++ linker can
//! find it by its unmangled symbol name. Each one is a thin wrapper: rebuild the
//! Rust argument types from raw C pointers, call the safe `cfsem` function, and
//! map the result to a C-friendly return value.
//!
//! Naming convention: every exported symbol is prefixed `cfsem_`. The companion
//! C++ header (`include/cfsem.hpp`) re-wraps these in `namespace cfsem`, dropping
//! the prefix.
//!
//! Array convention (see `cfsem_flux_density_linear_filament`): coordinate
//! components are passed as separate `*const f64` pointers plus element counts;
//! results are written into caller-allocated `*mut f64` buffers; the return value
//! is an `int` status code (0 = success).

use std::slice;

/// (dimensionless) Complete elliptic integral of the first kind, parameter `m`.
///
/// Scalar smoke-test wrapper over [`cfsem::math::ellipk`].
#[no_mangle]
pub extern "C" fn cfsem_ellipk(m: f64) -> f64 {
    cfsem::math::ellipk(m)
}

/// Biot-Savart magnetic flux density from many straight current filaments at many
/// observation points. Thin wrapper over
/// [`cfsem::physics::linear_filament::flux_density_linear_filament`].
///
/// # Arguments
///
/// Observation points (each array length `n_obs`):
/// * `xp`, `yp`, `zp`  — (m) coordinates.
///
/// Filament segments (each array length `n_fil`):
/// * `xfil`, `yfil`, `zfil` — (m) segment start coordinates.
/// * `dlx`, `dly`, `dlz`    — (m) segment length deltas (end = start + delta).
/// * `ifil`                 — (A) segment current.
/// * `wire_radius`          — (m) conductor (half-)thickness.
///
/// Outputs (caller-allocated, each length `n_obs`):
/// * `bx`, `by`, `bz` — (T) flux density components, overwritten in place.
///
/// # Returns
///
/// * `0` — success.
/// * `1` — the underlying `cfsem` call reported an error (e.g. length mismatch).
/// * `2` — a required pointer was null.
///
/// # Safety
///
/// Every pointer must be non-null and point to at least the stated number of
/// `f64` elements. Output buffers must not alias the inputs or each other.
#[no_mangle]
pub unsafe extern "C" fn cfsem_flux_density_linear_filament(
    xp: *const f64,
    yp: *const f64,
    zp: *const f64,
    n_obs: usize,
    xfil: *const f64,
    yfil: *const f64,
    zfil: *const f64,
    dlx: *const f64,
    dly: *const f64,
    dlz: *const f64,
    ifil: *const f64,
    wire_radius: *const f64,
    n_fil: usize,
    bx: *mut f64,
    by: *mut f64,
    bz: *mut f64,
) -> i32 {
    // Defensive null check — cheap, and turns the most common caller mistake into
    // a clean error code instead of undefined behavior.
    let in_ptrs = [
        xp,
        yp,
        zp,
        xfil,
        yfil,
        zfil,
        dlx,
        dly,
        dlz,
        ifil,
        wire_radius,
    ];
    if in_ptrs.iter().any(|p| p.is_null()) || bx.is_null() || by.is_null() || bz.is_null() {
        return 2;
    }

    // Rebuild the crate's tuple-of-slices arguments from the raw pointers.
    let xyzp = (
        slice::from_raw_parts(xp, n_obs),
        slice::from_raw_parts(yp, n_obs),
        slice::from_raw_parts(zp, n_obs),
    );
    let xyzfil = (
        slice::from_raw_parts(xfil, n_fil),
        slice::from_raw_parts(yfil, n_fil),
        slice::from_raw_parts(zfil, n_fil),
    );
    let dlxyzfil = (
        slice::from_raw_parts(dlx, n_fil),
        slice::from_raw_parts(dly, n_fil),
        slice::from_raw_parts(dlz, n_fil),
    );
    let ifil = slice::from_raw_parts(ifil, n_fil);
    let wire_radius = slice::from_raw_parts(wire_radius, n_fil);
    let out = (
        slice::from_raw_parts_mut(bx, n_obs),
        slice::from_raw_parts_mut(by, n_obs),
        slice::from_raw_parts_mut(bz, n_obs),
    );

    match cfsem::physics::linear_filament::flux_density_linear_filament(
        xyzp,
        xyzfil,
        dlxyzfil,
        ifil,
        wire_radius,
        out,
    ) {
        Ok(()) => 0,
        Err(_) => 1,
    }
}
