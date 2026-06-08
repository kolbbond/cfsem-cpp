//! C-ABI shim over the `cfsem` crate.
//!
//! Every exported function is `#[no_mangle] pub extern "C"`, prefixed `cfsem_`,
//! and declared in `include/cfsem.h`. The companion C++ header `cfsem.hpp`
//! re-wraps these in `namespace cfsem`.
//!
//! ## Data layout — interleaved 3xN
//!
//! Coordinate arrays use the **interleaved** layout `[x0,y0,z0, x1,y1,z1, ...]`,
//! which matches an Armadillo `arma::Mat<double>(3, N)` in memory (column-major).
//! Downstream C++ (e.g. goose) can therefore pass `mat.memptr()` directly with no
//! transpose. Internally `cfsem` wants structure-of-arrays, so each wrapper
//! de-interleaves into temporary component vectors (an O(N) copy per call —
//! acceptable for a validation/reference backend) and re-interleaves the result.
//!
//! ## Panic safety
//!
//! A Rust `panic!` unwinding across the FFI boundary is undefined behavior, so
//! every wrapper runs its body inside `catch_unwind` and maps a caught panic to a
//! status code rather than letting it cross into C/C++.
//!
//! ## Status codes
//!
//! * `0` — success.
//! * `1` — the underlying `cfsem` call reported an error (e.g. length mismatch).
//! * `2` — a required pointer was null.
//! * `3` — a panic was caught at the FFI boundary.

// The C ABI dictates flat argument lists; bundling them into a struct would
// change the contract C/C++ callers depend on. This lint doesn't apply here.
#![allow(clippy::too_many_arguments)]

use std::panic::{catch_unwind, AssertUnwindSafe};
use std::slice;

/// Signature shared by the linear-filament field functions in `cfsem`
/// (`flux_density_linear_filament`, `vector_potential_linear_filament`): identical
/// argument shape, different physics.
type LinFilamentFn = fn(
    (&[f64], &[f64], &[f64]),             // observation points (xs, ys, zs)
    (&[f64], &[f64], &[f64]),             // filament start coords
    (&[f64], &[f64], &[f64]),             // filament segment deltas
    &[f64],                               // currents
    &[f64],                               // wire radii (regularization scale)
    (&mut [f64], &mut [f64], &mut [f64]), // output (xs, ys, zs)
) -> Result<(), &'static str>;

/// (dimensionless) Complete elliptic integral of the first kind, parameter `m`.
/// Scalar smoke-test wrapper over [`cfsem::math::ellipk`].
#[no_mangle]
pub extern "C" fn cfsem_ellipk(m: f64) -> f64 {
    cfsem::math::ellipk(m)
}

/// Biot-Savart magnetic flux density (T) from many straight current filaments at
/// many observation points. Wraps
/// [`cfsem::physics::linear_filament::flux_density_linear_filament`].
///
/// See [`linear_filament_field`] for the argument and return-code contract.
///
/// # Safety
///
/// All pointers must be non-null and sized as documented in `cfsem.h`.
#[no_mangle]
pub unsafe extern "C" fn cfsem_flux_density_linear_filament(
    rs_obs: *const f64,
    n_obs: usize,
    rs_fil: *const f64,
    drs_fil: *const f64,
    ifil: *const f64,
    wire_radius: *const f64,
    n_fil: usize,
    b_out: *mut f64,
) -> i32 {
    catch_unwind(AssertUnwindSafe(|| {
        linear_filament_field(
            rs_obs,
            n_obs,
            rs_fil,
            drs_fil,
            ifil,
            wire_radius,
            n_fil,
            b_out,
            cfsem::physics::linear_filament::flux_density_linear_filament,
        )
    }))
    .unwrap_or(3)
}

/// Magnetic vector potential (V·s/m) from many straight current filaments at many
/// observation points. Wraps
/// [`cfsem::physics::linear_filament::vector_potential_linear_filament`].
///
/// See [`linear_filament_field`] for the argument and return-code contract.
///
/// # Safety
///
/// All pointers must be non-null and sized as documented in `cfsem.h`.
#[no_mangle]
pub unsafe extern "C" fn cfsem_vector_potential_linear_filament(
    rs_obs: *const f64,
    n_obs: usize,
    rs_fil: *const f64,
    drs_fil: *const f64,
    ifil: *const f64,
    wire_radius: *const f64,
    n_fil: usize,
    a_out: *mut f64,
) -> i32 {
    catch_unwind(AssertUnwindSafe(|| {
        linear_filament_field(
            rs_obs,
            n_obs,
            rs_fil,
            drs_fil,
            ifil,
            wire_radius,
            n_fil,
            a_out,
            cfsem::physics::linear_filament::vector_potential_linear_filament,
        )
    }))
    .unwrap_or(3)
}

/// Shared implementation for the linear-filament field wrappers.
///
/// Validates pointers, de-interleaves the `3xN` inputs into component vectors,
/// calls the supplied `cfsem` function `f`, and re-interleaves the result into
/// `out`. Returns a status code (see module docs).
///
/// # Arguments
///
/// * `rs_obs`:  interleaved observation coords, length `3 * n_obs`.
/// * `rs_fil`:  interleaved filament segment start coords, length `3 * n_fil`.
/// * `drs_fil`: interleaved filament segment deltas (end = start + delta),
///   length `3 * n_fil`.
/// * `ifil`:    filament currents (A), length `n_fil`.
/// * `wire_radius`: conductor (half-)thickness / regularization scale (m),
///   length `n_fil`.
/// * `out`:     caller-allocated interleaved output, length `3 * n_obs`.
///
/// # Safety
///
/// Every pointer must be non-null and point to at least the stated number of
/// `f64` elements. `out` must not alias the inputs.
unsafe fn linear_filament_field(
    rs_obs: *const f64,
    n_obs: usize,
    rs_fil: *const f64,
    drs_fil: *const f64,
    ifil: *const f64,
    wire_radius: *const f64,
    n_fil: usize,
    out: *mut f64,
    f: LinFilamentFn,
) -> i32 {
    if rs_obs.is_null()
        || rs_fil.is_null()
        || drs_fil.is_null()
        || ifil.is_null()
        || wire_radius.is_null()
        || out.is_null()
    {
        return 2;
    }

    // Borrow the raw inputs as slices (interleaved for the 3xN coordinate arrays).
    let rs_obs = slice::from_raw_parts(rs_obs, 3 * n_obs);
    let rs_fil = slice::from_raw_parts(rs_fil, 3 * n_fil);
    let drs_fil = slice::from_raw_parts(drs_fil, 3 * n_fil);
    let ifil = slice::from_raw_parts(ifil, n_fil);
    let wire_radius = slice::from_raw_parts(wire_radius, n_fil);

    // De-interleave into the structure-of-arrays form cfsem expects.
    let (xp, yp, zp) = deinterleave3(rs_obs, n_obs);
    let (xfil, yfil, zfil) = deinterleave3(rs_fil, n_fil);
    let (dlx, dly, dlz) = deinterleave3(drs_fil, n_fil);

    // Component output buffers, re-interleaved into `out` after the call.
    let mut ox = vec![0.0_f64; n_obs];
    let mut oy = vec![0.0_f64; n_obs];
    let mut oz = vec![0.0_f64; n_obs];

    let result = f(
        (&xp, &yp, &zp),
        (&xfil, &yfil, &zfil),
        (&dlx, &dly, &dlz),
        ifil,
        wire_radius,
        (&mut ox, &mut oy, &mut oz),
    );
    if result.is_err() {
        return 1;
    }

    let out = slice::from_raw_parts_mut(out, 3 * n_obs);
    for i in 0..n_obs {
        out[3 * i] = ox[i];
        out[3 * i + 1] = oy[i];
        out[3 * i + 2] = oz[i];
    }
    0
}

/// Split an interleaved `[x0,y0,z0, x1,y1,z1, ...]` slice of length `3 * n` into
/// three component vectors `(xs, ys, zs)`, each length `n`.
fn deinterleave3(rs: &[f64], n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut xs = Vec::with_capacity(n);
    let mut ys = Vec::with_capacity(n);
    let mut zs = Vec::with_capacity(n);
    for i in 0..n {
        xs.push(rs[3 * i]);
        ys.push(rs[3 * i + 1]);
        zs.push(rs[3 * i + 2]);
    }
    (xs, ys, zs)
}
