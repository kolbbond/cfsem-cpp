// cfsem.hpp — C++ namespace wrapper over the cfsem C ABI (cfsem.h).
//
// These are thin `inline` forwarders inside `namespace cfsem`; they add no
// behavior of their own, only ergonomics (drop the `cfsem_` prefix, get a real
// namespace). Header-only — no extra linking beyond libcfsem_capi.a.
#ifndef CFSEM_HPP
#define CFSEM_HPP

#include <cstddef> // std::size_t

#include "cfsem.h"

namespace cfsem {

// Complete elliptic integral of the first kind, parameter m. Dimensionless.
inline double ellipk(double m) { return ::cfsem_ellipk(m); }

// Biot-Savart flux density from many straight filaments. See cfsem.h for the
// full argument and return-code contract.
inline int flux_density_linear_filament(
    const double* xp, const double* yp, const double* zp, std::size_t n_obs,
    const double* xfil, const double* yfil, const double* zfil,
    const double* dlx, const double* dly, const double* dlz,
    const double* ifil, const double* wire_radius, std::size_t n_fil,
    double* bx, double* by, double* bz) {
  return ::cfsem_flux_density_linear_filament(
      xp, yp, zp, n_obs, xfil, yfil, zfil, dlx, dly, dlz, ifil, wire_radius,
      n_fil, bx, by, bz);
}

} // namespace cfsem

#endif // CFSEM_HPP
