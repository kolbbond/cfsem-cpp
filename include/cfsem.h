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

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CFSEM_H */
