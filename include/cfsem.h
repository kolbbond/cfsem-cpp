/*
 * cfsem.h — C ABI for the cfsem electromagnetics library.
 *
 * These declarations mirror the `#[no_mangle] extern "C"` functions exported by
 * the cfsem-capi Rust crate (libcfsem_capi.a). Link against that static library
 * and include this header to call cfsem from C or C++.
 *
 * For an ergonomic C++ surface, prefer cfsem.hpp, which wraps these in
 * `namespace cfsem`.
 *
 * Data layout: coordinate arrays are INTERLEAVED 3xN, i.e.
 *   [x0,y0,z0, x1,y1,z1, ...], matching an Armadillo arma::Mat<double>(3, N)
 *   in memory (column-major). Pass mat.memptr() directly.
 *
 * Status codes (functions returning int):
 *   0  success
 *   1  internal cfsem error (e.g. length mismatch)
 *   2  a required pointer was NULL
 *   3  a panic was caught at the FFI boundary
 *   4  a caller-supplied buffer was too small (required length written to *_out)
 *   5  an argument was out of range (e.g. an unrecognized skip code)
 *
 * Ownership: every buffer is caller-allocated and caller-freed. Nothing returned
 * by this library needs to be freed.
 */
#ifndef CFSEM_H
#define CFSEM_H

#include <stddef.h> /* size_t */

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Complete elliptic integral of the first kind, parameter m. Dimensionless.
 */
double cfsem_ellipk(double m);

/*
 * Biot-Savart magnetic flux density (T) from many straight current filaments at
 * many observation points.
 *
 *   rs_obs      [in]  interleaved observation coords,        length 3*n_obs   [m]
 *   n_obs       [in]  number of observation points
 *   rs_fil      [in]  interleaved filament segment START,    length 3*n_fil   [m]
 *   drs_fil     [in]  interleaved filament segment DELTA,    length 3*n_fil   [m]
 *                     (segment end = start + delta)
 *   ifil        [in]  filament currents,                     length n_fil     [A]
 *   wire_radius [in]  conductor (half-)thickness / softening, length n_fil    [m]
 *   n_fil       [in]  number of filament segments
 *   b_out       [out] interleaved Bx,By,Bz, caller-allocated, length 3*n_obs  [T]
 *
 * b_out must not alias the inputs. Returns a status code (see file header).
 */
int cfsem_flux_density_linear_filament(
    const double* rs_obs, size_t n_obs,
    const double* rs_fil, const double* drs_fil,
    const double* ifil, const double* wire_radius, size_t n_fil,
    double* b_out);

/*
 * Magnetic vector potential (V*s/m) from many straight current filaments at many
 * observation points. Same argument contract as
 * cfsem_flux_density_linear_filament; a_out receives interleaved Ax,Ay,Az,
 * length 3*n_obs.
 */
int cfsem_vector_potential_linear_filament(
    const double* rs_obs, size_t n_obs,
    const double* rs_fil, const double* drs_fil,
    const double* ifil, const double* wire_radius, size_t n_fil,
    double* a_out);

/*
 * Which class of interaction to OMIT in a hierarchical solve. A full solve equals
 * the sum of the near-only and far-only solves, up to floating-point roundoff.
 */
#define CFSEM_SKIP_NONE 0 /* evaluate near and far                */
#define CFSEM_SKIP_NEAR 1 /* far-field only                       */
#define CFSEM_SKIP_FAR  2 /* near-field (direct) only             */
#define CFSEM_SKIP_BOTH 3 /* zero output, no kernel evaluation    */

/*
 * Hierarchical (Barnes-Hut) flux density (T). Same contract as the direct version
 * plus theta (acceptance angle; smaller = more accurate, slower), par (nonzero
 * = parallel), and skip (one of CFSEM_SKIP_*). APPROXIMATE — no guaranteed error
 * bound; not for safety field limits. Tree builder fixed to "longest axis".
 *
 * Skipped interactions are not evaluated, so skip=CFSEM_SKIP_NEAR is cheaper than
 * a full solve, not just a filtered one.
 */
int cfsem_flux_density_linear_filament_hierarchical(
    const double* rs_obs, size_t n_obs,
    const double* rs_fil, const double* drs_fil,
    const double* ifil, const double* wire_radius, size_t n_fil,
    double theta, int par, int skip,
    double* b_out);

/* Hierarchical (Barnes-Hut) vector potential (V*s/m). Caveats as above. */
int cfsem_vector_potential_linear_filament_hierarchical(
    const double* rs_obs, size_t n_obs,
    const double* rs_fil, const double* drs_fil,
    const double* ifil, const double* wire_radius, size_t n_fil,
    double theta, int par, int skip,
    double* a_out);

/*
 * Near-field (direct) interaction pattern chosen by the Barnes-Hut acceptance
 * test, as a canonical CSC sparsity pattern: rows are filament (source) indices,
 * columns are observation (target) indices, shape (n_fil, n_obs). Row indices are
 * sorted and unique within each column.
 *
 * This is the sparsity of the near-field mutual-inductance matrix. Pass filament
 * MIDPOINTS as rs_obs to classify segment-to-segment interactions, then hand the
 * pattern to cfsem_inductance_linear_filaments_sparse_csc.
 *
 * No field is evaluated — this is a classification walk only. The A- and B-field
 * filament kernels share one acceptance function, so one map serves both.
 *
 *   theta, par      [in]  as in the hierarchical solvers
 *   row_indices     [out] CSC row indices, length >= nnz; NULL for a sizing call
 *   row_capacity    [in]  allocated length of row_indices
 *   column_pointers [out] CSC column offsets, length n_obs+1; may be NULL
 *   nnz_out         [out] number of stored interactions (always written)
 *
 * Two-call protocol:
 *   1. row_indices = NULL  -> *nnz_out = required length
 *   2. allocate, call again -> pattern filled
 * Returns 4 if row_capacity < nnz, with *nnz_out set.
 *
 * The pattern is NOT symmetric: acceptance tests a target against a source node,
 * so i near j does not imply j near i. Symmetrizing (union of P and P^T) is a
 * modeling choice left to the caller. Rebuild the map whenever geometry, theta,
 * or currents change.
 */
int cfsem_near_field_interaction_map_linear_filament(
    const double* rs_obs, size_t n_obs,
    const double* rs_fil, const double* drs_fil,
    const double* ifil, const double* wire_radius, size_t n_fil,
    double theta, int par,
    size_t* row_indices, size_t row_capacity,
    size_t* column_pointers,
    size_t* nnz_out);

/*
 * Mutual inductance (H) evaluated only on a caller-supplied CSC pattern — the
 * sparse near-field M matrix. Rows are source filaments, columns are target
 * filaments, shape (n_src, n_tgt), matching the map above.
 *
 * Each entry is integrated with 3-point Gauss-Legendre over the WHOLE target
 * segment, so entries are segment-to-segment inductances, not pointwise vector
 * potentials sampled at a midpoint.
 *
 *   row_indices     [in]  CSC row indices, length nnz
 *   column_pointers [in]  CSC column offsets, length n_tgt+1
 *   nnz             [in]  number of stored entries
 *   values_out      [out] caller-allocated, length nnz, CSC order              [H]
 *
 * Explicit numerical zeros are retained. Returns 1 if the pattern is
 * non-canonical (unsorted or duplicate rows, decreasing column pointers, or a row
 * index >= n_src).
 */
int cfsem_inductance_linear_filaments_sparse_csc(
    const double* rs_fil_tgt, const double* drs_fil_tgt, size_t n_tgt,
    const double* rs_fil_src, const double* drs_fil_src,
    const double* wire_radius_src, size_t n_src,
    const size_t* row_indices, const size_t* column_pointers, size_t nnz,
    int par,
    double* values_out);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CFSEM_H */
