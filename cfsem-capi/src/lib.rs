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
//! Downstream C++ can therefore pass `mat.memptr()` directly with no
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
//! * `4` — a caller-supplied buffer was too small; the required length was written
//!   to the corresponding `*_out` size argument.
//! * `5` — an argument was outside its valid range (e.g. an unrecognized `skip`).
//!
//! ## Variable-length output
//!
//! The near-field interaction map has a length that is only known after the tree
//! walk, so it uses a two-call protocol: call once with a null data pointer to
//! learn the length, allocate, then call again to fill. Nothing is allocated on
//! the Rust side and handed to C — every buffer stays caller-owned.

// The C ABI dictates flat argument lists; bundling them into a struct would
// change the contract C/C++ callers depend on. This lint doesn't apply here.
#![allow(clippy::too_many_arguments)]

use std::panic::{catch_unwind, AssertUnwindSafe};
use std::slice;
use std::time::Instant;

// The tree/evaluator pieces are only `pub` via their submodule paths (the
// mod-level re-exports are crate-private).
use cfsem::physics::hierarchical::evaluator::{
    eval, eval_par, scratch_len, scratch_len_par, update_summaries, EvaluationScratch,
    SourceNodeSummaries,
};
use cfsem::physics::hierarchical::kernel::{HierarchicalError, HierarchicalKernel};
use cfsem::physics::hierarchical::kernels::{
    DipoleTarget, DipoleTargets, LinearFilamentFluxDensityKernel, LinearFilamentFluxDensitySummary,
    LinearFilamentSource, LinearFilamentSources, LinearFilamentVectorPotentialKernel,
    LinearFilamentVectorPotentialSummary,
};
use cfsem::physics::hierarchical::tree::{BuildMethod, ClusterTree};
use cfsem::physics::hierarchical::Skip;

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

/// Hierarchical (Barnes-Hut) flux density (T). Approximate — `theta` trades
/// accuracy for speed, no guaranteed bound; not for safety field limits. Wraps
/// [`cfsem::physics::hierarchical::flux_density_linear_filament_hierarchical`];
/// see [`linear_filament_field_hierarchical`] for the contract.
///
/// # Safety
///
/// All pointers must be non-null and sized as documented in `cfsem.h`.
#[no_mangle]
pub unsafe extern "C" fn cfsem_flux_density_linear_filament_hierarchical(
    rs_obs: *const f64,
    n_obs: usize,
    rs_fil: *const f64,
    drs_fil: *const f64,
    ifil: *const f64,
    wire_radius: *const f64,
    n_fil: usize,
    theta: f64,
    par: i32,
    skip: i32,
    b_out: *mut f64,
) -> i32 {
    catch_unwind(AssertUnwindSafe(|| {
        let skip = match decode_skip(skip) {
            Ok(skip) => skip,
            Err(code) => return code,
        };
        linear_filament_field_hierarchical(
            rs_obs,
            n_obs,
            rs_fil,
            drs_fil,
            ifil,
            wire_radius,
            n_fil,
            theta,
            par != 0,
            b_out,
            |xyzp, xyzfil, dlxyzfil, ifil, wire_radius, method, theta, par, extra, out| match skip {
                None => cfsem::physics::hierarchical::flux_density_linear_filament_hierarchical(
                    xyzp, xyzfil, dlxyzfil, ifil, wire_radius, method, theta, par, extra, out,
                ),
                Some(skip) => {
                    cfsem::physics::hierarchical::flux_density_linear_filament_hierarchical_with_skip(
                        xyzp, xyzfil, dlxyzfil, ifil, wire_radius, method, theta, par, skip, extra,
                        out,
                    )
                }
            },
        )
    }))
    .unwrap_or(3)
}

/// Hierarchical (Barnes-Hut) vector potential (V·s/m). Approximate; same caveats
/// as [`cfsem_flux_density_linear_filament_hierarchical`]. Wraps
/// [`cfsem::physics::hierarchical::vector_potential_linear_filament_hierarchical`].
///
/// # Safety
///
/// All pointers must be non-null and sized as documented in `cfsem.h`.
#[no_mangle]
pub unsafe extern "C" fn cfsem_vector_potential_linear_filament_hierarchical(
    rs_obs: *const f64,
    n_obs: usize,
    rs_fil: *const f64,
    drs_fil: *const f64,
    ifil: *const f64,
    wire_radius: *const f64,
    n_fil: usize,
    theta: f64,
    par: i32,
    skip: i32,
    a_out: *mut f64,
) -> i32 {
    catch_unwind(AssertUnwindSafe(|| {
        let skip = match decode_skip(skip) {
            Ok(skip) => skip,
            Err(code) => return code,
        };
        linear_filament_field_hierarchical(
            rs_obs,
            n_obs,
            rs_fil,
            drs_fil,
            ifil,
            wire_radius,
            n_fil,
            theta,
            par != 0,
            a_out,
            |xyzp, xyzfil, dlxyzfil, ifil, wire_radius, method, theta, par, extra, out| match skip {
                None => cfsem::physics::hierarchical::vector_potential_linear_filament_hierarchical(
                    xyzp, xyzfil, dlxyzfil, ifil, wire_radius, method, theta, par, extra, out,
                ),
                Some(skip) => {
                    cfsem::physics::hierarchical::vector_potential_linear_filament_hierarchical_with_skip(
                        xyzp, xyzfil, dlxyzfil, ifil, wire_radius, method, theta, par, skip, extra,
                        out,
                    )
                }
            },
        )
    }))
    .unwrap_or(3)
}

/// Map a C `skip` code to the `cfsem` interaction filter.
///
/// `0` none, `1` near, `2` far, `3` both. `Ok(None)` means "evaluate everything".
fn decode_skip(skip: i32) -> Result<Option<Skip>, i32> {
    match skip {
        0 => Ok(None),
        1 => Ok(Some(Skip::Near)),
        2 => Ok(Some(Skip::Far)),
        3 => Ok(Some(Skip::Both)),
        _ => Err(5),
    }
}

/// Near-field (direct) source-target interaction pattern selected by the
/// Barnes-Hut acceptance test, in canonical CSC: rows are filament (source)
/// indices, columns are observation (target) indices, shape `(n_fil, n_obs)`.
///
/// This is the sparsity pattern of the near-field mutual-inductance matrix: pass
/// filament midpoints as `rs_obs` to classify segment-to-segment interactions,
/// then feed the result to [`cfsem_inductance_linear_filaments_sparse_csc`].
///
/// Runs the tree walk with both kernel classes skipped, so no field is evaluated —
/// this is a classification pass only. The A- and B-field filament kernels share
/// one acceptance function, so the same map is valid for both.
///
/// Two-call protocol: call with `row_indices` null to learn `*nnz_out`, allocate,
/// then call again. Returns `4` if `row_capacity < nnz` (with `*nnz_out` set).
///
/// The pattern is **not symmetric** — acceptance tests a target against a source
/// node, so `i` near `j` does not imply `j` near `i`. Symmetrizing (union of
/// `P` and `P^T`) is a modeling choice left to the caller.
///
/// The map must be rebuilt whenever geometry, `theta`, or currents change.
///
/// # Safety
///
/// All pointers must be non-null (except `row_indices` and `column_pointers`, which
/// may be null on a sizing call) and sized as documented in `cfsem.h`.
#[no_mangle]
pub unsafe extern "C" fn cfsem_near_field_interaction_map_linear_filament(
    rs_obs: *const f64,
    n_obs: usize,
    rs_fil: *const f64,
    drs_fil: *const f64,
    ifil: *const f64,
    wire_radius: *const f64,
    n_fil: usize,
    theta: f64,
    par: i32,
    row_indices: *mut usize,
    row_capacity: usize,
    column_pointers: *mut usize,
    nnz_out: *mut usize,
) -> i32 {
    catch_unwind(AssertUnwindSafe(|| {
        if rs_obs.is_null()
            || rs_fil.is_null()
            || drs_fil.is_null()
            || ifil.is_null()
            || wire_radius.is_null()
            || nnz_out.is_null()
        {
            return 2;
        }

        let rs_obs = slice::from_raw_parts(rs_obs, 3 * n_obs);
        let rs_fil = slice::from_raw_parts(rs_fil, 3 * n_fil);
        let drs_fil = slice::from_raw_parts(drs_fil, 3 * n_fil);
        let ifil = slice::from_raw_parts(ifil, n_fil);
        let wire_radius = slice::from_raw_parts(wire_radius, n_fil);

        let (xp, yp, zp) = deinterleave3(rs_obs, n_obs);
        let (xfil, yfil, zfil) = deinterleave3(rs_fil, n_fil);
        let (dlx, dly, dlz) = deinterleave3(drs_fil, n_fil);

        // Skip::Both zeroes these without evaluating a kernel; only the tree walk runs.
        let mut ox = vec![0.0_f64; n_obs];
        let mut oy = vec![0.0_f64; n_obs];
        let mut oz = vec![0.0_f64; n_obs];

        let diagnostics =
            cfsem::physics::hierarchical::vector_potential_linear_filament_hierarchical_with_skip(
                (&xp, &yp, &zp),
                (&xfil, &yfil, &zfil),
                (&dlx, &dly, &dlz),
                ifil,
                wire_radius,
                BuildMethod::LongestAxis,
                theta,
                par != 0,
                Skip::Both,
                true, // extra_diagnostics: the near map is the whole point of this call
                (&mut ox, &mut oy, &mut oz),
            );
        let Ok(diagnostics) = diagnostics else {
            return 1;
        };
        let Some(map) = diagnostics.near_field_interaction_map() else {
            return 1;
        };

        let nnz = map.row_indices.len();
        *nnz_out = nnz;

        if !column_pointers.is_null() {
            let dst = slice::from_raw_parts_mut(column_pointers, n_obs + 1);
            dst.copy_from_slice(&map.column_pointers);
        }

        if row_indices.is_null() {
            return 0; // sizing call
        }
        if row_capacity < nnz {
            return 4;
        }
        let dst = slice::from_raw_parts_mut(row_indices, nnz);
        // Widen u32 -> size_t so the C caller can hand these straight back to the
        // sparse inductance fill without a conversion pass.
        for (slot, &row) in dst.iter_mut().zip(map.row_indices.iter()) {
            *slot = row as usize;
        }
        0
    }))
    .unwrap_or(3)
}

/// Mutual inductance (H) evaluated only on a caller-supplied CSC sparsity pattern —
/// the sparse near-field `M` matrix. Wraps
/// [`cfsem::physics::linear_filament::inductance_linear_filaments_sparse_csc`].
///
/// Rows are source filaments, columns are target filaments, shape `(n_src, n_tgt)`,
/// matching the pattern from
/// [`cfsem_near_field_interaction_map_linear_filament`]. Each stored entry is
/// integrated with 3-point Gauss-Legendre over the whole target segment, so the
/// entries are segment-to-segment inductances, not pointwise vector potentials.
/// `values_out` holds one value per stored coordinate, in CSC order, and retains
/// explicit numerical zeros.
///
/// # Safety
///
/// All pointers must be non-null and sized as documented in `cfsem.h`.
#[no_mangle]
pub unsafe extern "C" fn cfsem_inductance_linear_filaments_sparse_csc(
    rs_fil_tgt: *const f64,
    drs_fil_tgt: *const f64,
    n_tgt: usize,
    rs_fil_src: *const f64,
    drs_fil_src: *const f64,
    wire_radius_src: *const f64,
    n_src: usize,
    row_indices: *const usize,
    column_pointers: *const usize,
    nnz: usize,
    par: i32,
    values_out: *mut f64,
) -> i32 {
    catch_unwind(AssertUnwindSafe(|| {
        if rs_fil_tgt.is_null()
            || drs_fil_tgt.is_null()
            || rs_fil_src.is_null()
            || drs_fil_src.is_null()
            || wire_radius_src.is_null()
            || row_indices.is_null()
            || column_pointers.is_null()
            || values_out.is_null()
        {
            return 2;
        }

        let rs_fil_tgt = slice::from_raw_parts(rs_fil_tgt, 3 * n_tgt);
        let drs_fil_tgt = slice::from_raw_parts(drs_fil_tgt, 3 * n_tgt);
        let rs_fil_src = slice::from_raw_parts(rs_fil_src, 3 * n_src);
        let drs_fil_src = slice::from_raw_parts(drs_fil_src, 3 * n_src);
        let wire_radius_src = slice::from_raw_parts(wire_radius_src, n_src);
        let row_indices = slice::from_raw_parts(row_indices, nnz);
        let column_pointers = slice::from_raw_parts(column_pointers, n_tgt + 1);
        let values_out = slice::from_raw_parts_mut(values_out, nnz);

        let (xt, yt, zt) = deinterleave3(rs_fil_tgt, n_tgt);
        let (dxt, dyt, dzt) = deinterleave3(drs_fil_tgt, n_tgt);
        let (xs, ys, zs) = deinterleave3(rs_fil_src, n_src);
        let (dxs, dys, dzs) = deinterleave3(drs_fil_src, n_src);

        let f = if par != 0 {
            cfsem::physics::linear_filament::inductance_linear_filaments_sparse_csc_par
        } else {
            cfsem::physics::linear_filament::inductance_linear_filaments_sparse_csc
        };
        let result = f(
            (&xt, &yt, &zt),
            (&dxt, &dyt, &dzt),
            (&xs, &ys, &zs),
            (&dxs, &dys, &dzs),
            wire_radius_src,
            row_indices,
            column_pointers,
            values_out,
        );
        if result.is_err() {
            return 1;
        }
        0
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

/// Hierarchical counterpart of [`linear_filament_field`]: same de-/re-interleave
/// contract, plus `theta` and `par`; tree builder fixed to
/// [`BuildMethod::LongestAxis`].
///
/// `R`/`E` are generic so one helper serves both kernels' `Diagnostics` — only
/// `Result::is_err` is read (any error → status `1`).
///
/// # Safety
///
/// Every pointer must be non-null and point to at least the documented number of
/// `f64` elements. `out` must not alias the inputs.
unsafe fn linear_filament_field_hierarchical<R, E, F>(
    rs_obs: *const f64,
    n_obs: usize,
    rs_fil: *const f64,
    drs_fil: *const f64,
    ifil: *const f64,
    wire_radius: *const f64,
    n_fil: usize,
    theta: f64,
    par: bool,
    out: *mut f64,
    f: F,
) -> i32
where
    F: Fn(
        (&[f64], &[f64], &[f64]),
        (&[f64], &[f64], &[f64]),
        (&[f64], &[f64], &[f64]),
        &[f64],
        &[f64],
        BuildMethod,
        f64,
        bool,
        bool,
        (&mut [f64], &mut [f64], &mut [f64]),
    ) -> Result<R, E>,
{
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
        BuildMethod::LongestAxis,
        theta,
        par,
        false, // extra_diagnostics: tree/near-map diagnostics not exported over the C ABI
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

/// Persistent Barnes-Hut tree over fixed linear-filament geometry.
///
/// The one-shot hierarchical wrappers rebuild the tree on every call; a solver
/// whose geometry is fixed and whose currents change per iteration pays that
/// build each time. This handle builds once and refreshes only the node
/// summaries when currents change (lazily, per kernel, on the next evaluation).
pub struct LinearFilamentTree {
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
    dlx: Vec<f64>,
    dly: Vec<f64>,
    dlz: Vec<f64>,
    wire_radius: Vec<f64>,
    ifil: Vec<f64>,
    tree: ClusterTree<f64>,
    a_summaries: Vec<LinearFilamentVectorPotentialSummary<f64>>,
    b_summaries: Vec<LinearFilamentFluxDensitySummary<f64>>,
    a_stale: bool,
    b_stale: bool,
    build_seconds: f64,
    summary_seconds: f64,
    eval_seconds: f64,
}

impl LinearFilamentTree {
    fn sources(&self) -> LinearFilamentSources<'_, f64> {
        LinearFilamentSources::new(
            (&self.x, &self.y, &self.z),
            (&self.dlx, &self.dly, &self.dlz),
            &self.wire_radius,
        )
    }
}

/// Build a persistent tree over `n_fil` filaments (same layout as
/// [`cfsem_vector_potential_linear_filament_hierarchical`], minus currents).
/// Currents start at zero; set them with [`cfsem_lf_tree_set_currents`].
/// Returns null on a null pointer, tree-construction failure or panic. Free
/// with [`cfsem_lf_tree_free`].
///
/// # Safety
///
/// All pointers must be non-null and sized as documented in `cfsem.h`.
#[no_mangle]
pub unsafe extern "C" fn cfsem_lf_tree_new(
    rs_fil: *const f64,
    drs_fil: *const f64,
    wire_radius: *const f64,
    n_fil: usize,
) -> *mut LinearFilamentTree {
    if rs_fil.is_null() || drs_fil.is_null() || wire_radius.is_null() || n_fil == 0 {
        return std::ptr::null_mut();
    }
    catch_unwind(AssertUnwindSafe(|| {
        let start = Instant::now();
        let (x, y, z) = deinterleave3(slice::from_raw_parts(rs_fil, 3 * n_fil), n_fil);
        let (dlx, dly, dlz) = deinterleave3(slice::from_raw_parts(drs_fil, 3 * n_fil), n_fil);
        let wire_radius = slice::from_raw_parts(wire_radius, n_fil).to_vec();
        let sources = LinearFilamentSources::new((&x, &y, &z), (&dlx, &dly, &dlz), &wire_radius);
        let tree = match ClusterTree::build_with_method(sources, BuildMethod::LongestAxis) {
            Ok(t) => t,
            Err(_) => return std::ptr::null_mut(),
        };
        let a_summaries =
            SourceNodeSummaries::<LinearFilamentVectorPotentialKernel<f64>>::new(tree.as_view())
                .node_summaries;
        let b_summaries =
            SourceNodeSummaries::<LinearFilamentFluxDensityKernel<f64>>::new(tree.as_view())
                .node_summaries;
        let build_seconds = start.elapsed().as_secs_f64();
        Box::into_raw(Box::new(LinearFilamentTree {
            x,
            y,
            z,
            dlx,
            dly,
            dlz,
            wire_radius,
            ifil: vec![0.0; n_fil],
            tree,
            a_summaries,
            b_summaries,
            a_stale: true,
            b_stale: true,
            build_seconds,
            summary_seconds: 0.0,
            eval_seconds: 0.0,
        }))
    }))
    .unwrap_or(std::ptr::null_mut())
}

/// Release a tree from [`cfsem_lf_tree_new`]. Null is a no-op.
///
/// # Safety
///
/// `tree` must be null or a live pointer from [`cfsem_lf_tree_new`], not used afterwards.
#[no_mangle]
pub unsafe extern "C" fn cfsem_lf_tree_free(tree: *mut LinearFilamentTree) {
    if !tree.is_null() {
        drop(Box::from_raw(tree));
    }
}

/// Replace the filament currents (A), length `n_fil`. Summaries are refreshed
/// lazily on the next evaluation of each field.
///
/// # Safety
///
/// `tree` must be live; `ifil` must hold `n_fil` values.
#[no_mangle]
pub unsafe extern "C" fn cfsem_lf_tree_set_currents(
    tree: *mut LinearFilamentTree,
    ifil: *const f64,
) -> i32 {
    if tree.is_null() || ifil.is_null() {
        return 2;
    }
    let t = &mut *tree;
    let n = t.ifil.len();
    t.ifil.copy_from_slice(slice::from_raw_parts(ifil, n));
    t.a_stale = true;
    t.b_stale = true;
    0
}

/// Hierarchical vector potential (V·s/m) at `n_obs` interleaved targets from the
/// cached tree. `theta`/`par`/`a_out` as in
/// [`cfsem_vector_potential_linear_filament_hierarchical`].
///
/// # Safety
///
/// `tree` must be live; `rs_obs`/`a_out` sized `3 * n_obs`, not aliased.
#[no_mangle]
pub unsafe extern "C" fn cfsem_lf_tree_vector_potential(
    tree: *mut LinearFilamentTree,
    rs_obs: *const f64,
    n_obs: usize,
    theta: f64,
    par: i32,
    a_out: *mut f64,
) -> i32 {
    if tree.is_null() || rs_obs.is_null() || a_out.is_null() {
        return 2;
    }
    catch_unwind(AssertUnwindSafe(|| {
        let t = &mut *tree;
        let kernel = LinearFilamentVectorPotentialKernel::<f64>::new();
        let mut summaries = std::mem::take(&mut t.a_summaries);
        let stale = t.a_stale;
        let rc = tree_eval(
            t,
            &kernel,
            &mut summaries,
            stale,
            rs_obs,
            n_obs,
            theta,
            par != 0,
            a_out,
        );
        t.a_summaries = summaries;
        if rc == 0 {
            t.a_stale = false;
        }
        rc
    }))
    .unwrap_or(3)
}

/// Hierarchical flux density (T) from the cached tree. Same contract as
/// [`cfsem_lf_tree_vector_potential`].
///
/// # Safety
///
/// `tree` must be live; `rs_obs`/`b_out` sized `3 * n_obs`, not aliased.
#[no_mangle]
pub unsafe extern "C" fn cfsem_lf_tree_flux_density(
    tree: *mut LinearFilamentTree,
    rs_obs: *const f64,
    n_obs: usize,
    theta: f64,
    par: i32,
    b_out: *mut f64,
) -> i32 {
    if tree.is_null() || rs_obs.is_null() || b_out.is_null() {
        return 2;
    }
    catch_unwind(AssertUnwindSafe(|| {
        let t = &mut *tree;
        let kernel = LinearFilamentFluxDensityKernel::<f64>::new();
        let mut summaries = std::mem::take(&mut t.b_summaries);
        let stale = t.b_stale;
        let rc = tree_eval(
            t,
            &kernel,
            &mut summaries,
            stale,
            rs_obs,
            n_obs,
            theta,
            par != 0,
            b_out,
        );
        t.b_summaries = summaries;
        if rc == 0 {
            t.b_stale = false;
        }
        rc
    }))
    .unwrap_or(3)
}

/// Wall-clock seconds of the tree build, the last summary refresh and the last
/// evaluation. Any out pointer may be null.
///
/// # Safety
///
/// `tree` must be live.
#[no_mangle]
pub unsafe extern "C" fn cfsem_lf_tree_timings(
    tree: *const LinearFilamentTree,
    build_s: *mut f64,
    summary_s: *mut f64,
    eval_s: *mut f64,
) -> i32 {
    if tree.is_null() {
        return 2;
    }
    let t = &*tree;
    for (p, v) in [
        (build_s, t.build_seconds),
        (summary_s, t.summary_seconds),
        (eval_s, t.eval_seconds),
    ] {
        if !p.is_null() {
            *p = v;
        }
    }
    0
}

/// Refresh `summaries` if `stale`, then evaluate kernel `K` at the targets.
unsafe fn tree_eval<K>(
    t: &mut LinearFilamentTree,
    kernel: &K,
    summaries: &mut [K::SourceSummary],
    stale: bool,
    rs_obs: *const f64,
    n_obs: usize,
    theta: f64,
    par: bool,
    out: *mut f64,
) -> i32
where
    K: HierarchicalKernel<
            Scalar = f64,
            Output = [f64; 3],
            SourceGeometry = LinearFilamentSource<f64>,
            SourceMoment = f64,
            TargetGeometry = DipoleTarget<f64>,
        > + Sync,
{
    if stale {
        let start = Instant::now();
        if update_summaries(
            kernel,
            t.tree.as_view(),
            t.sources(),
            t.ifil.as_slice(),
            summaries,
        ) != HierarchicalError::Ok
        {
            return 1;
        }
        t.summary_seconds = start.elapsed().as_secs_f64();
    }

    let start = Instant::now();
    let (xp, yp, zp) = deinterleave3(slice::from_raw_parts(rs_obs, 3 * n_obs), n_obs);
    let targets = DipoleTargets::new(&xp, &yp, &zp);
    let mut ox = vec![0.0_f64; n_obs];
    let mut oy = vec![0.0_f64; n_obs];
    let mut oz = vec![0.0_f64; n_obs];
    let mut scratch_values = vec![
        [0.0_f64; 3];
        if par {
            scratch_len_par(n_obs)
        } else {
            scratch_len()
        }
    ];
    let mut scratch = EvaluationScratch {
        contribution: &mut scratch_values,
    };
    let out_components = [ox.as_mut_slice(), oy.as_mut_slice(), oz.as_mut_slice()];
    let sources = t.sources();
    let moments = t.ifil.as_slice();
    let err = if par {
        eval_par(
            kernel,
            t.tree.as_view(),
            summaries,
            sources,
            targets,
            moments,
            theta,
            None,
            out_components,
            &mut scratch,
        )
    } else {
        eval(
            kernel,
            t.tree.as_view(),
            summaries,
            sources,
            targets,
            moments,
            theta,
            None,
            out_components,
            &mut scratch,
        )
    };
    if err != HierarchicalError::Ok {
        return 1;
    }

    let out = slice::from_raw_parts_mut(out, 3 * n_obs);
    for i in 0..n_obs {
        out[3 * i] = ox[i];
        out[3 * i + 1] = oy[i];
        out[3 * i + 2] = oz[i];
    }
    t.eval_seconds = start.elapsed().as_secs_f64();
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
