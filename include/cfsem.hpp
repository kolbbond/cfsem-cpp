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

// Interaction class to omit in a hierarchical solve; see cfsem.h.
enum Skip {
  kSkipNone = CFSEM_SKIP_NONE,
  kSkipNear = CFSEM_SKIP_NEAR, // far-field only
  kSkipFar = CFSEM_SKIP_FAR,   // near-field (direct) only
  kSkipBoth = CFSEM_SKIP_BOTH,
};

// Hierarchical (Barnes-Hut) flux density (T). Approximate; theta/par/skip as in
// cfsem.h. Same status-code contract.
inline int flux_density_linear_filament_hierarchical(
    const double* rs_obs, std::size_t n_obs,
    const double* rs_fil, const double* drs_fil,
    const double* ifil, const double* wire_radius, std::size_t n_fil,
    double theta, int par, int skip,
    double* b_out) {
  return ::cfsem_flux_density_linear_filament_hierarchical(
      rs_obs, n_obs, rs_fil, drs_fil, ifil, wire_radius, n_fil, theta, par,
      skip, b_out);
}

// Hierarchical (Barnes-Hut) vector potential (V*s/m). Caveats as above.
inline int vector_potential_linear_filament_hierarchical(
    const double* rs_obs, std::size_t n_obs,
    const double* rs_fil, const double* drs_fil,
    const double* ifil, const double* wire_radius, std::size_t n_fil,
    double theta, int par, int skip,
    double* a_out) {
  return ::cfsem_vector_potential_linear_filament_hierarchical(
      rs_obs, n_obs, rs_fil, drs_fil, ifil, wire_radius, n_fil, theta, par,
      skip, a_out);
}

// Near-field (direct) interaction pattern as canonical CSC, shape
// (n_fil, n_obs). Two-call protocol: pass row_indices = nullptr to size, then
// again to fill. See cfsem.h for the full contract and the asymmetry caveat.
inline int near_field_interaction_map_linear_filament(
    const double* rs_obs, std::size_t n_obs,
    const double* rs_fil, const double* drs_fil,
    const double* ifil, const double* wire_radius, std::size_t n_fil,
    double theta, int par,
    std::size_t* row_indices, std::size_t row_capacity,
    std::size_t* column_pointers,
    std::size_t* nnz_out) {
  return ::cfsem_near_field_interaction_map_linear_filament(
      rs_obs, n_obs, rs_fil, drs_fil, ifil, wire_radius, n_fil, theta, par,
      row_indices, row_capacity, column_pointers, nnz_out);
}

// Mutual inductance (H) on a caller-supplied CSC pattern — the sparse near-field
// M matrix, shape (n_src, n_tgt). See cfsem.h.
inline int inductance_linear_filaments_sparse_csc(
    const double* rs_fil_tgt, const double* drs_fil_tgt, std::size_t n_tgt,
    const double* rs_fil_src, const double* drs_fil_src,
    const double* wire_radius_src, std::size_t n_src,
    const std::size_t* row_indices, const std::size_t* column_pointers,
    std::size_t nnz, int par,
    double* values_out) {
  return ::cfsem_inductance_linear_filaments_sparse_csc(
      rs_fil_tgt, drs_fil_tgt, n_tgt, rs_fil_src, drs_fil_src, wire_radius_src,
      n_src, row_indices, column_pointers, nnz, par, values_out);
}

} // namespace cfsem

#endif // CFSEM_HPP
