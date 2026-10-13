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
 * Hierarchical (Barnes-Hut) flux density (T). Same contract as the direct version
 * plus theta (acceptance angle; smaller = more accurate, slower) and par (nonzero
 * = parallel). APPROXIMATE — no guaranteed error bound; not for safety field
 * limits. Tree builder fixed to "longest axis".
 */
int cfsem_flux_density_linear_filament_hierarchical(
    const double* rs_obs, size_t n_obs,
    const double* rs_fil, const double* drs_fil,
    const double* ifil, const double* wire_radius, size_t n_fil,
    double theta, int par,
    double* b_out);

/* Hierarchical (Barnes-Hut) vector potential (V*s/m). Caveats as above. */
int cfsem_vector_potential_linear_filament_hierarchical(
    const double* rs_obs, size_t n_obs,
    const double* rs_fil, const double* drs_fil,
    const double* ifil, const double* wire_radius, size_t n_fil,
    double theta, int par,
    double* a_out);

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
