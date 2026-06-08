// Smoke test for the cfsem C-ABI shim library.
//
// Not the deliverable — it proves the linkable `cfsem_cpp` library builds, links,
// and runs, and sanity-checks the numbers against a textbook result. Coordinate
// arrays use the interleaved 3xN layout the C API expects (see cfsem.h).
// Returns nonzero if any cfsem call fails its status code.

#include <cmath>
#include <cstddef>
#include <cstdio>

#include "cfsem.hpp"

namespace {
constexpr double kMu0 = 1.25663706212e-6; // (H/m) vacuum permeability
constexpr double kPi = 3.14159265358979323846;
} // namespace

int main() {
  std::printf("Initialize cfsem C-ABI shim\n");

  // 1) Scalar path: complete elliptic integral of the first kind.
  std::printf("cfsem::ellipk(0.1) = %.6f\n", cfsem::ellipk(0.1));

  // 2) Field paths: a single long, straight filament approximating an infinite
  //    wire. Segment runs along +x from (-L/2,0,0); observe at perpendicular
  //    distance r on the +y axis. For L >> r the flux density approaches the
  //    textbook infinite-wire result B = mu0 * I / (2*pi*r), pointing in +z.
  const double length = 2.0e3; // (m) very long vs the standoff distance
  const double current = 1.0e3; // (A)
  const double r = 0.1;         // (m) perpendicular standoff

  const double rs_obs[3] = {0.0, r, 0.0};               // 1 observation point
  const double rs_fil[3] = {-0.5 * length, 0.0, 0.0};   // segment start
  const double drs_fil[3] = {length, 0.0, 0.0};         // segment delta
  const double ifil[1] = {current};
  const double wire_radius[1] = {1.0e-4};               // thin wire
  const std::size_t n_obs = 1;
  const std::size_t n_fil = 1;

  double b_out[3] = {0.0, 0.0, 0.0};
  double a_out[3] = {0.0, 0.0, 0.0};

  const int rc_b = cfsem::flux_density_linear_filament(
      rs_obs, n_obs, rs_fil, drs_fil, ifil, wire_radius, n_fil, b_out);
  const int rc_a = cfsem::vector_potential_linear_filament(
      rs_obs, n_obs, rs_fil, drs_fil, ifil, wire_radius, n_fil, a_out);

  std::printf("B (rc=%d) = (%.6e, %.6e, %.6e) T\n", rc_b, b_out[0], b_out[1],
              b_out[2]);
  std::printf("A (rc=%d) = (%.6e, %.6e, %.6e) V*s/m\n", rc_a, a_out[0], a_out[1],
              a_out[2]);

  if (rc_b != 0 || rc_a != 0) {
    std::printf("FAIL: cfsem call failed (rc_b=%d rc_a=%d)\n", rc_b, rc_a);
    return rc_b != 0 ? rc_b : rc_a;
  }

  // Analytic sanity check on the B-field magnitude (informational).
  const double b_analytic = kMu0 * current / (2.0 * kPi * r);
  const double b_mag =
      std::sqrt(b_out[0] * b_out[0] + b_out[1] * b_out[1] + b_out[2] * b_out[2]);
  const double rel_err = std::fabs(b_mag - b_analytic) / b_analytic;
  std::printf("infinite-wire check: |B|=%.6e T, analytic=%.6e T, rel_err=%.2e\n",
              b_mag, b_analytic, rel_err);
  std::printf(rel_err < 1e-2 ? "OK\n" : "WARN: B differs from analytic > 1%%\n");
  return 0;
}
