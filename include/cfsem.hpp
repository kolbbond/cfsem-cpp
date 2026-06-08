// cfsem.hpp — C++ namespace wrapper over the cfsem C ABI (cfsem.h).
//
// These are thin `inline` forwarders inside `namespace cfsem`; they add no
// behavior of their own, only ergonomics (drop the `cfsem_` prefix, get a real
// namespace). Header-only — no extra linking beyond libcfsem_capi.a.
//
// Coordinate arrays are interleaved 3xN (see cfsem.h), matching
// arma::Mat<double>(3, N) memory, so e.g. `cfsem::flux_density_linear_filament(
// Rt.memptr(), Rt.n_cols, ...)` works directly.
#ifndef CFSEM_HPP
#define CFSEM_HPP

#include <cstddef> // std::size_t

#include "cfsem.h"

namespace cfsem {

// Complete elliptic integral of the first kind, parameter m. Dimensionless.
inline double ellipk(double m) { return ::cfsem_ellipk(m); }

// Biot-Savart flux density (T) from many straight filaments. See cfsem.h for the
// full argument and status-code contract.
inline int flux_density_linear_filament(
    const double* rs_obs, std::size_t n_obs,
    const double* rs_fil, const double* drs_fil,
    const double* ifil, const double* wire_radius, std::size_t n_fil,
    double* b_out) {
  return ::cfsem_flux_density_linear_filament(
      rs_obs, n_obs, rs_fil, drs_fil, ifil, wire_radius, n_fil, b_out);
}

// Vector potential (V*s/m) from many straight filaments. Same contract as above.
inline int vector_potential_linear_filament(
    const double* rs_obs, std::size_t n_obs,
    const double* rs_fil, const double* drs_fil,
    const double* ifil, const double* wire_radius, std::size_t n_fil,
    double* a_out) {
  return ::cfsem_vector_potential_linear_filament(
      rs_obs, n_obs, rs_fil, drs_fil, ifil, wire_radius, n_fil, a_out);
}

} // namespace cfsem

#endif // CFSEM_HPP
