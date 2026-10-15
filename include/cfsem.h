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
 *   in memory (column-major). Pass mat.memptr() directly. Axisymmetric
 *   (circular-filament) coordinates are interleaved 2xN, [r0,z0, r1,z1, ...].
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

/*
 * Complete elliptic integral of the second kind, parameter m. Dimensionless.
 */
double cfsem_ellipe(double m);

/*
 * Mutual inductance (H) from ALL source segments to EACH target segment, i.e.
 * the column sums of cfsem_inductance_linear_filaments_matrix without storing
 * the matrix. Same 3-point Gauss-Legendre segment integral as the sparse fill.
 *
 *   rs_fil_tgt, drs_fil_tgt [in]  interleaved target start/delta, 3*n_tgt [m]
 *   rs_fil_src, drs_fil_src [in]  interleaved source start/delta, 3*n_src [m]
 *   wire_radius_src         [in]  source radii, length n_src             [m]
 *   out                     [out] caller-allocated, length n_tgt         [H]
 *
 * Returns a status code (see file header).
 */
int cfsem_inductance_linear_filaments(
    const double* rs_fil_tgt, const double* drs_fil_tgt, size_t n_tgt,
    const double* rs_fil_src, const double* drs_fil_src,
    const double* wire_radius_src, size_t n_src,
    double* out);

/*
 * Dense mutual-inductance matrix (H). Arguments as above, plus par (nonzero =
 * parallel). m_out is caller-allocated, length n_src*n_tgt, ROW-MAJOR
 * (n_src, n_tgt): entry (i, j) = m_out[i*n_tgt + j]. In Armadillo terms that is
 * arma::Mat<double>(n_tgt, n_src), i.e. the transpose. Returns 5 if
 * n_src*n_tgt overflows size_t.
 */
int cfsem_inductance_linear_filaments_matrix(
    const double* rs_fil_tgt, const double* drs_fil_tgt, size_t n_tgt,
    const double* rs_fil_src, const double* drs_fil_src,
    const double* wire_radius_src, size_t n_src,
    int par,
    double* m_out);

/*
 * Per-source flux density (T) response matrices: entry (i, j) is the field at
 * observation i due to filament j, currents applied, so summing each row
 * recovers cfsem_flux_density_linear_filament. Same inputs as that function,
 * plus par (nonzero = parallel).
 *
 *   bx_out, by_out, bz_out [out] caller-allocated, each length n_obs*n_fil,
 *                                ROW-MAJOR (n_obs, n_fil)                  [T]
 *
 * Components are separate buffers, not interleaved, so the matrices are written
 * in place without a copy. In Armadillo terms each is
 * arma::Mat<double>(n_fil, n_obs). Returns 5 if n_obs*n_fil overflows size_t.
 */
int cfsem_flux_density_linear_filament_matrix(
    const double* rs_obs, size_t n_obs,
    const double* rs_fil, const double* drs_fil,
    const double* ifil, const double* wire_radius, size_t n_fil,
    int par,
    double* bx_out, double* by_out, double* bz_out);

/* Per-source vector potential (V*s/m) response matrices; contract as above. */
int cfsem_vector_potential_linear_filament_matrix(
    const double* rs_obs, size_t n_obs,
    const double* rs_fil, const double* drs_fil,
    const double* ifil, const double* wire_radius, size_t n_fil,
    int par,
    double* ax_out, double* ay_out, double* az_out);

/*
 * J x B body force density (N/m^3) on a current density at each observation
 * point, due to the filaments' field.
 *
 *   rs_fil, drs_fil, ifil, wire_radius, n_fil  [in] filaments, as above
 *   rs_obs  [in]  interleaved observation coords,  length 3*n_obs   [m]
 *   j_obs   [in]  interleaved current density,     length 3*n_obs   [A/m^2]
 *   par     [in]  nonzero = parallel
 *   f_out   [out] interleaved fx,fy,fz, caller-allocated, 3*n_obs   [N/m^3]
 *
 * Returns a status code (see file header).
 */
int cfsem_body_force_density_linear_filament(
    const double* rs_fil, const double* drs_fil,
    const double* ifil, const double* wire_radius, size_t n_fil,
    const double* rs_obs, const double* j_obs, size_t n_obs,
    int par,
    double* f_out);

/*
 * Circular filaments: closed forms via elliptic integrals. The (r, z) functions
 * take coaxial loops about the z axis; the Cartesian ones take loops at any
 * center and orientation. Shared arguments:
 *
 *   rz_fil      [in]  interleaved loop (r, z),          length 2*n_fil   [m]
 *   r_fil       [in]  loop major radii,                 length n_fil     [m]
 *   ifil        [in]  loop currents,                    length n_fil     [A-turns]
 *   loc         [in]  interleaved 3xN loop centers,     length 3*n_fil   [m]
 *   normal      [in]  interleaved 3xN loop normals,     length 3*n_fil
 *                     (normalized internally; positive current is
 *                     right-handed about it)
 *   wire_radius [in]  conductor cross-section radius,   length n_fil     [m]
 *                     (0 = ideal filament)
 *   rz_obs      [in]  interleaved observation (r, z),   length 2*n_obs   [m]
 *   rs_obs      [in]  interleaved 3xN observation pts,  length 3*n_obs   [m]
 *   par         [in]  nonzero = parallel over observation points
 *
 * Finite on axis and for negative r. An ideal filament is singular on the loop
 * itself; with wire_radius > 0 the field blends to a finite-section model within
 * 1.5-3 wire radii of the centerline. Invalid geometry (e.g. a zero normal)
 * yields NaN, not a status code. Each returns a status code.
 */

/* Poloidal flux (Wb): total flux through the loop of radius r, = 2*pi*r*A_phi.
 * psi_out has length n_obs. */
int cfsem_flux_circular_filament(
    const double* rz_fil, const double* ifil, const double* wire_radius,
    size_t n_fil,
    const double* rz_obs, size_t n_obs,
    int par,
    double* psi_out);

/* A_phi (V*s/m); A_r = A_z = 0 by symmetry. a_phi_out has length n_obs. */
int cfsem_vector_potential_circular_filament(
    const double* rz_fil, const double* ifil, const double* wire_radius,
    size_t n_fil,
    const double* rz_obs, size_t n_obs,
    int par,
    double* a_phi_out);

/* (B_r, B_z) (T), b_out interleaved, length 2*n_obs. */
int cfsem_flux_density_circular_filament(
    const double* rz_fil, const double* ifil, const double* wire_radius,
    size_t n_fil,
    const double* rz_obs, size_t n_obs,
    int par,
    double* b_out);

/* Cartesian B (T) at rs_obs; b_out interleaved, length 3*n_obs. */
int cfsem_flux_density_circular_filament_cartesian(
    const double* r_fil, const double* ifil, const double* loc,
    const double* normal, const double* wire_radius, size_t n_fil,
    const double* rs_obs, size_t n_obs,
    int par,
    double* b_out);

/* Cartesian A (V*s/m) at rs_obs; a_out interleaved, length 3*n_obs. */
int cfsem_vector_potential_circular_filament_cartesian(
    const double* r_fil, const double* ifil, const double* loc,
    const double* normal, const double* wire_radius, size_t n_fil,
    const double* rs_obs, size_t n_obs,
    int par,
    double* a_out);

/*
 * J x B body force density (N/m^3) on current density j_obs [A/m^2], interleaved
 * 3xN, at rs_obs. f_out interleaved, length 3*n_obs.
 */
int cfsem_body_force_density_circular_filament(
    const double* r_fil, const double* ifil, const double* loc,
    const double* normal, const double* wire_radius, size_t n_fil,
    const double* rs_obs, const double* j_obs, size_t n_obs,
    int par,
    double* f_out);

/*
 * Persistent Barnes-Hut tree over fixed filament geometry. The one-shot
 * *_hierarchical functions rebuild the tree every call; for fixed geometry with
 * changing currents, build once here and only refresh currents.
 *
 *   cfsem_lf_tree_new          build over n_fil filaments (rs/drs/wire_radius as
 *                              above); currents start at zero. NULL on failure.
 *   cfsem_lf_tree_set_currents replace currents (length n_fil) [A]
 *   cfsem_lf_tree_*            evaluate at n_obs targets; theta/par/out as in the
 *                              one-shot hierarchical functions
 *   cfsem_lf_tree_timings      seconds of the build, last summary refresh, last
 *                              evaluation; any out pointer may be NULL
 *
 * A handle is not thread-safe: evaluation refreshes its cached summaries.
 */
typedef struct cfsem_lf_tree cfsem_lf_tree;

cfsem_lf_tree* cfsem_lf_tree_new(
    const double* rs_fil, const double* drs_fil,
    const double* wire_radius, size_t n_fil);

void cfsem_lf_tree_free(cfsem_lf_tree* tree);

int cfsem_lf_tree_set_currents(cfsem_lf_tree* tree, const double* ifil);

int cfsem_lf_tree_vector_potential(
    cfsem_lf_tree* tree, const double* rs_obs, size_t n_obs,
    double theta, int par, double* a_out);

int cfsem_lf_tree_flux_density(
    cfsem_lf_tree* tree, const double* rs_obs, size_t n_obs,
    double theta, int par, double* b_out);

int cfsem_lf_tree_timings(
    const cfsem_lf_tree* tree, double* build_s, double* summary_s, double* eval_s);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CFSEM_H */
