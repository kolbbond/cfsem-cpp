/*
 * cfsem.h — C ABI for the cfsem electromagnetics library.
 *
 * These declarations mirror the `#[no_mangle] extern "C"` functions exported by
 * the cfsem-capi Rust crate (libcfsem_capi.a). Link against that static library
 * and include this header to call cfsem from C or C++.
 *
 * For an ergonomic C++ surface, prefer cfsem.hpp, which wraps these in
 * `namespace cfsem`.
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
 * Observation points  (each array length n_obs): xp, yp, zp        [m]
 * Filament segments    (each array length n_fil):
 *     xfil, yfil, zfil  — segment start coordinates                [m]
 *     dlx,  dly,  dlz   — segment length deltas (end = start+delta) [m]
 *     ifil              — segment current                           [A]
 *     wire_radius       — conductor (half-)thickness                [m]
 * Outputs (caller-allocated, each length n_obs): bx, by, bz        [T]
 *
 * Returns: 0 on success, 1 on an internal cfsem error (e.g. length mismatch),
 *          2 if any required pointer is NULL.
 *
 * Output buffers must not alias the inputs or each other.
 */
int cfsem_flux_density_linear_filament(
    const double* xp, const double* yp, const double* zp, size_t n_obs,
    const double* xfil, const double* yfil, const double* zfil,
    const double* dlx, const double* dly, const double* dlz,
    const double* ifil, const double* wire_radius, size_t n_fil,
    double* bx, double* by, double* bz);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CFSEM_H */
