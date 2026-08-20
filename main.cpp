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
      kSeg, theta, /*par=*/0, CFSEM_SKIP_NONE, b_hier);
  const int rc_ad = cfsem::vector_potential_linear_filament(
      rs_obs, n_obs, seg_fil.data(), seg_dfil.data(), seg_i.data(), seg_wr.data(),
      kSeg, a_dir);
  const int rc_ah = cfsem::vector_potential_linear_filament_hierarchical(
      rs_obs, n_obs, seg_fil.data(), seg_dfil.data(), seg_i.data(), seg_wr.data(),
      kSeg, theta, /*par=*/0, CFSEM_SKIP_NONE, a_hier);

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

  // --- skip: near + far must partition the full solve exactly ---------------
  // Not an approximation — it is one sum split in two, so it holds to roundoff.
  double a_near[3] = {0.0, 0.0, 0.0};
  double a_far[3] = {0.0, 0.0, 0.0};
  double a_none[3] = {1.0, 1.0, 1.0}; // nonzero, to prove SKIP_BOTH overwrites
  const int rc_near = cfsem::vector_potential_linear_filament_hierarchical(
      rs_obs, n_obs, seg_fil.data(), seg_dfil.data(), seg_i.data(), seg_wr.data(),
      kSeg, theta, /*par=*/0, CFSEM_SKIP_FAR, a_near);
  const int rc_far = cfsem::vector_potential_linear_filament_hierarchical(
      rs_obs, n_obs, seg_fil.data(), seg_dfil.data(), seg_i.data(), seg_wr.data(),
      kSeg, theta, /*par=*/0, CFSEM_SKIP_NEAR, a_far);
  const int rc_both = cfsem::vector_potential_linear_filament_hierarchical(
      rs_obs, n_obs, seg_fil.data(), seg_dfil.data(), seg_i.data(), seg_wr.data(),
      kSeg, theta, /*par=*/0, CFSEM_SKIP_BOTH, a_none);
  const int rc_bad = cfsem::vector_potential_linear_filament_hierarchical(
      rs_obs, n_obs, seg_fil.data(), seg_dfil.data(), seg_i.data(), seg_wr.data(),
      kSeg, theta, /*par=*/0, /*skip=*/99, a_none);

  if (rc_near != 0 || rc_far != 0 || rc_both != 0) {
    std::printf("FAIL: skip call failed (near=%d far=%d both=%d)\n", rc_near,
                rc_far, rc_both);
    return 1;
  }
  if (rc_bad != 5) {
    std::printf("FAIL: bad skip code returned %d, expected 5\n", rc_bad);
    return 1;
  }

  double split_err = 0.0;
  for (int k = 0; k < 3; ++k) {
    split_err = std::fmax(split_err, std::fabs(a_hier[k] - (a_near[k] + a_far[k])));
  }
  const double split_rel = split_err / mag3(a_hier);
  const bool both_zeroed =
      a_none[0] == 0.0 && a_none[1] == 0.0 && a_none[2] == 0.0;
  std::printf("skip partition: |A_full - (A_near + A_far)| / |A_full| = %.2e, "
              "SKIP_BOTH zeroed = %s\n",
              split_rel, both_zeroed ? "yes" : "no");
  if (split_rel > 1e-12 || !both_zeroed) {
    std::printf("FAIL: near/far split is not a partition of the full solve\n");
    return 1;
  }
  std::printf("OK\n");

  // --- near-field interaction map + sparse M --------------------------------
  // Targets are segment midpoints, so the map classifies segment-to-segment
  // interactions: exactly the sparsity of the near-field inductance matrix.
  std::vector<double> mid(3 * kSeg);
  for (std::size_t s = 0; s < kSeg; ++s) {
    for (int k = 0; k < 3; ++k) {
      mid[3 * s + k] = seg_fil[3 * s + k] + 0.5 * seg_dfil[3 * s + k];
    }
  }

  // theta=0 would accept nothing; use a coarse angle so the map is genuinely
  // sparse rather than the full dense pattern.
  const double map_theta = 0.5;
  std::size_t nnz = 0;
  const int rc_size = cfsem::near_field_interaction_map_linear_filament(
      mid.data(), kSeg, seg_fil.data(), seg_dfil.data(), seg_i.data(),
      seg_wr.data(), kSeg, map_theta, /*par=*/0, /*row_indices=*/nullptr,
      /*row_capacity=*/0, /*column_pointers=*/nullptr, &nnz);
  if (rc_size != 0 || nnz == 0 || nnz >= kSeg * kSeg) {
    std::printf("FAIL: map sizing call rc=%d nnz=%zu\n", rc_size, nnz);
    return 1;
  }

  std::vector<std::size_t> rows(nnz), colptr(kSeg + 1);
  std::size_t nnz_check = 0;
  const int rc_small = cfsem::near_field_interaction_map_linear_filament(
      mid.data(), kSeg, seg_fil.data(), seg_dfil.data(), seg_i.data(),
      seg_wr.data(), kSeg, map_theta, /*par=*/0, rows.data(),
      /*row_capacity=*/nnz - 1, colptr.data(), &nnz_check);
  if (rc_small != 4 || nnz_check != nnz) {
    std::printf("FAIL: undersized buffer returned %d (expected 4), nnz=%zu\n",
                rc_small, nnz_check);
    return 1;
  }

  const int rc_map = cfsem::near_field_interaction_map_linear_filament(
      mid.data(), kSeg, seg_fil.data(), seg_dfil.data(), seg_i.data(),
      seg_wr.data(), kSeg, map_theta, /*par=*/0, rows.data(), nnz,
      colptr.data(), &nnz_check);
  if (rc_map != 0 || nnz_check != nnz || colptr[kSeg] != nnz) {
    std::printf("FAIL: map fill rc=%d nnz=%zu colptr[n]=%zu\n", rc_map,
                nnz_check, colptr[kSeg]);
    return 1;
  }

  // Canonical CSC: nondecreasing column pointers, sorted unique rows per column.
  for (std::size_t c = 0; c < kSeg; ++c) {
    if (colptr[c] > colptr[c + 1]) {
      std::printf("FAIL: column pointers decrease at %zu\n", c);
      return 1;
    }
    for (std::size_t e = colptr[c] + 1; e < colptr[c + 1]; ++e) {
      if (rows[e] <= rows[e - 1] || rows[e] >= kSeg) {
        std::printf("FAIL: non-canonical rows in column %zu\n", c);
        return 1;
      }
    }
  }
  std::printf("near map (theta=%.2f): nnz=%zu of %zu (%.1f%% dense), canonical "
              "CSC\n",
              map_theta, nnz, kSeg * kSeg,
              100.0 * static_cast<double>(nnz) /
                  static_cast<double>(kSeg * kSeg));
  std::printf("OK\n");

  std::vector<double> m_near(nnz, 0.0);
  const int rc_m = cfsem::inductance_linear_filaments_sparse_csc(
      seg_fil.data(), seg_dfil.data(), kSeg, seg_fil.data(), seg_dfil.data(),
      seg_wr.data(), kSeg, rows.data(), colptr.data(), nnz, /*par=*/0,
      m_near.data());
  if (rc_m != 0) {
    std::printf("FAIL: sparse inductance rc=%d\n", rc_m);
    return 1;
  }

  // Mutual inductance is symmetric, so wherever the (asymmetric) pattern happens
  // to store both (i,j) and (j,i), the two values must agree.
  double sym_rel = 0.0;
  double m_max = 0.0;
  std::size_t pairs = 0;
  for (std::size_t c = 0; c < kSeg; ++c) {
    for (std::size_t e = colptr[c]; e < colptr[c + 1]; ++e) {
      const std::size_t rrow = rows[e];
      m_max = std::fmax(m_max, std::fabs(m_near[e]));
      if (rrow == c) {
        continue;
      }
      for (std::size_t f = colptr[rrow]; f < colptr[rrow + 1]; ++f) {
        if (rows[f] == c) {
          sym_rel = std::fmax(sym_rel, std::fabs(m_near[e] - m_near[f]));
          ++pairs;
          break;
        }
      }
    }
  }
  sym_rel = m_max > 0.0 ? sym_rel / m_max : 0.0;
  std::printf("sparse M: %zu entries, max |M|=%.3e H, %zu mutual pairs, "
              "asymmetry=%.2e\n",
              nnz, m_max, pairs, sym_rel);
  if (!(m_max > 0.0) || pairs == 0 || sym_rel > 1e-12) {
    std::printf("FAIL: sparse M failed the symmetry check\n");
    return 1;
  }
  std::printf("OK\n");
  return 0;
}
