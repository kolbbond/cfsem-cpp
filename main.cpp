// Smoke test for the cfsem C-ABI shim: proves it links and the numbers are sane.
// Interleaved 3xN arrays throughout (see cfsem.h). Nonzero on any status failure.

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <vector>

#include "cfsem.hpp"

namespace {
constexpr double kMu0 = 1.25663706212e-6; // (H/m) vacuum permeability
constexpr double kPi = 3.14159265358979323846;

double mag3(const double v[3]) {
  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}
} // namespace

int main() {
  std::printf("Initialize cfsem C-ABI shim\n");

  // scalar path
  std::printf("cfsem::ellipk(0.1) = %.6f\n", cfsem::ellipk(0.1));

  // Long straight filament ~ infinite wire: for L >> r, B -> mu0*I/(2*pi*r) in
  // +z, observed at perpendicular standoff r on +y.
  const double length = 2.0e3;  // (m), >> r
  const double current = 1.0e3; // (A)
  const double r = 0.1;         // (m) perpendicular standoff

  const double rs_obs[3] = {0.0, r, 0.0};
  const double rs_fil[3] = {-0.5 * length, 0.0, 0.0}; // start
  const double drs_fil[3] = {length, 0.0, 0.0};       // delta (end = start+delta)
  const double ifil[1] = {current};
  const double wire_radius[1] = {1.0e-4};
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

  // analytic check
  const double b_analytic = kMu0 * current / (2.0 * kPi * r);
  const double b_mag = mag3(b_out);
  const double rel_err = std::fabs(b_mag - b_analytic) / b_analytic;
  std::printf("infinite-wire check: |B|=%.6e T, analytic=%.6e T, rel_err=%.2e\n",
              b_mag, b_analytic, rel_err);
  std::printf(rel_err < 1e-2 ? "OK\n" : "WARN: B differs from analytic > 1%%\n");

  // Hierarchical path: one source is trivially accepted, so discretize the wire
  // into many segments to actually build/walk the tree. Small theta -> tracks the
  // exact direct sum (the validation regime).
  constexpr std::size_t kSeg = 400;
  std::vector<double> seg_fil(3 * kSeg), seg_dfil(3 * kSeg), seg_i(kSeg),
      seg_wr(kSeg);
  const double seg = length / static_cast<double>(kSeg);
  for (std::size_t s = 0; s < kSeg; ++s) {
    seg_fil[3 * s] = -0.5 * length + static_cast<double>(s) * seg;
    seg_fil[3 * s + 1] = 0.0;
    seg_fil[3 * s + 2] = 0.0;
    seg_dfil[3 * s] = seg; // along +x
    seg_dfil[3 * s + 1] = 0.0;
    seg_dfil[3 * s + 2] = 0.0;
    seg_i[s] = current;
    seg_wr[s] = wire_radius[0];
  }

  double b_dir[3] = {0.0, 0.0, 0.0};
  double b_hier[3] = {0.0, 0.0, 0.0};
  double a_dir[3] = {0.0, 0.0, 0.0};
  double a_hier[3] = {0.0, 0.0, 0.0};
  const double theta = 0.1; // small -> near-direct accuracy

  const int rc_bd = cfsem::flux_density_linear_filament(
      rs_obs, n_obs, seg_fil.data(), seg_dfil.data(), seg_i.data(), seg_wr.data(),
      kSeg, b_dir);
  const int rc_bh = cfsem::flux_density_linear_filament_hierarchical(
      rs_obs, n_obs, seg_fil.data(), seg_dfil.data(), seg_i.data(), seg_wr.data(),
      kSeg, theta, /*par=*/0, b_hier);
  const int rc_ad = cfsem::vector_potential_linear_filament(
      rs_obs, n_obs, seg_fil.data(), seg_dfil.data(), seg_i.data(), seg_wr.data(),
      kSeg, a_dir);
  const int rc_ah = cfsem::vector_potential_linear_filament_hierarchical(
      rs_obs, n_obs, seg_fil.data(), seg_dfil.data(), seg_i.data(), seg_wr.data(),
      kSeg, theta, /*par=*/0, a_hier);

  if (rc_bd != 0 || rc_bh != 0 || rc_ad != 0 || rc_ah != 0) {
    std::printf("FAIL: hierarchical call failed (rc_bd=%d rc_bh=%d rc_ad=%d "
                "rc_ah=%d)\n",
                rc_bd, rc_bh, rc_ad, rc_ah);
    return rc_bh != 0 ? rc_bh : (rc_ah != 0 ? rc_ah : (rc_bd != 0 ? rc_bd : rc_ad));
  }

  const double b_rel = std::fabs(mag3(b_hier) - mag3(b_dir)) / mag3(b_dir);
  const double a_rel = std::fabs(mag3(a_hier) - mag3(a_dir)) / mag3(a_dir);
  std::printf("hierarchical vs direct (%zu segments, theta=%.2f): "
              "B rel_diff=%.2e, A rel_diff=%.2e\n",
              kSeg, theta, b_rel, a_rel);
  std::printf((b_rel < 1e-2 && a_rel < 1e-2)
                  ? "OK\n"
                  : "WARN: hierarchical differs from direct > 1%%\n");
  return 0;
}
