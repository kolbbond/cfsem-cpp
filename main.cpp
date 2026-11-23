// Smoke test for the cfsem C-ABI shim: proves it links and the numbers are sane.
// Interleaved 3xN arrays throughout (see cfsem.h). Nonzero on any status failure.

#include <cmath>
#include <cstddef>
#include <cstdint>
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

  // Persistent tree: same tree, summaries and walk as the one-shot path, so it
  // must match it exactly after every current update (both fields, both modes).
  constexpr std::size_t kObs = 64;
  std::vector<double> obs(3 * kObs);
  for (std::size_t j = 0; j < kObs; ++j) {
    obs[3 * j] = -0.4 * length + 0.8 * length * static_cast<double>(j) / kObs;
    obs[3 * j + 1] = r * static_cast<double>(1 + j % 5);
    obs[3 * j + 2] = 0.01 * static_cast<double>(j % 3);
  }
  cfsem::lf_tree* tree =
      cfsem::lf_tree_new(seg_fil.data(), seg_dfil.data(), seg_wr.data(), kSeg);
  if (tree == nullptr) {
    std::printf("FAIL: lf_tree_new returned NULL\n");
    return 1;
  }
  std::vector<double> ref(3 * kObs), got(3 * kObs);
  int mismatches = 0;
  for (int update = 0; update < 3; ++update) {
    for (std::size_t s = 0; s < kSeg; ++s) {
      seg_i[s] = current * std::cos(0.01 * static_cast<double>(s * (update + 1)));
    }
    if (cfsem::lf_tree_set_currents(tree, seg_i.data()) != 0) {
      std::printf("FAIL: lf_tree_set_currents\n");
      return 1;
    }
    for (int par = 0; par < 2; ++par) {
      int rc = cfsem::vector_potential_linear_filament_hierarchical(
          obs.data(), kObs, seg_fil.data(), seg_dfil.data(), seg_i.data(),
          seg_wr.data(), kSeg, theta, par, CFSEM_SKIP_NONE, ref.data());
      rc |= cfsem::lf_tree_vector_potential(tree, obs.data(), kObs, theta, par,
                                            got.data());
      mismatches += (rc != 0 || ref != got);
      rc = cfsem::flux_density_linear_filament_hierarchical(
          obs.data(), kObs, seg_fil.data(), seg_dfil.data(), seg_i.data(),
          seg_wr.data(), kSeg, theta, par, CFSEM_SKIP_NONE, ref.data());
      rc |= cfsem::lf_tree_flux_density(tree, obs.data(), kObs, theta, par,
                                        got.data());
      mismatches += (rc != 0 || ref != got);
    }
  }
  double build_s = 0.0, summary_s = 0.0, eval_s = 0.0;
  cfsem::lf_tree_timings(tree, &build_s, &summary_s, &eval_s);
  cfsem::lf_tree_free(tree);
  std::printf("persistent tree vs one-shot (3 current updates x A/B x "
              "serial/par): %d mismatches; build %.2e s, summary %.2e s, "
              "eval %.2e s\n",
              mismatches, build_s, summary_s, eval_s);
  if (mismatches != 0) {
    std::printf("FAIL: persistent tree differs from one-shot\n");
    return 1;
  }
  std::printf("OK\n");

  // E(0) = pi/2, E(1) = 1; the polynomial fit is good to ~2e-8.
  const double e0 = cfsem::ellipe(0.0);
  const double e1 = cfsem::ellipe(1.0 - 1e-15);
  std::printf("ellipe: E(0)=%.9f (pi/2=%.9f), E(1)=%.9f\n", e0, 0.5 * kPi, e1);
  if (std::fabs(e0 - 0.5 * kPi) > 1e-7 || std::fabs(e1 - 1.0) > 1e-6) {
    std::printf("FAIL: ellipe off its endpoints\n");
    return 1;
  }
  std::printf("OK\n");

  // Response matrices: row sums must reproduce the contracted field, and the
  // parallel path must match the serial one.
  {
    const std::size_t n_mat = kObs * kSeg;
    std::vector<double> mx(n_mat), my(n_mat), mz(n_mat);
    std::vector<double> px(n_mat), py(n_mat), pz(n_mat);
    std::vector<double> contracted(3 * kObs);
    double worst = 0.0;
    bool par_equal = true;
    for (int field = 0; field < 2; ++field) {
      const auto matrix = field == 0
                              ? cfsem::flux_density_linear_filament_matrix
                              : cfsem::vector_potential_linear_filament_matrix;
      const auto direct = field == 0 ? cfsem::flux_density_linear_filament
                                     : cfsem::vector_potential_linear_filament;
      int rc = matrix(obs.data(), kObs, seg_fil.data(), seg_dfil.data(),
                      seg_i.data(), seg_wr.data(), kSeg, /*par=*/0, mx.data(),
                      my.data(), mz.data());
      rc |= matrix(obs.data(), kObs, seg_fil.data(), seg_dfil.data(),
                   seg_i.data(), seg_wr.data(), kSeg, /*par=*/1, px.data(),
                   py.data(), pz.data());
      rc |= direct(obs.data(), kObs, seg_fil.data(), seg_dfil.data(),
                   seg_i.data(), seg_wr.data(), kSeg, contracted.data());
      if (rc != 0) {
        std::printf("FAIL: response-matrix call failed (field %d)\n", field);
        return 1;
      }
      par_equal = par_equal && mx == px && my == py && mz == pz;
      for (std::size_t i = 0; i < kObs; ++i) {
        double sum[3] = {0.0, 0.0, 0.0};
        for (std::size_t j = 0; j < kSeg; ++j) {
          sum[0] += mx[i * kSeg + j];
          sum[1] += my[i * kSeg + j];
          sum[2] += mz[i * kSeg + j];
        }
        double diff[3];
        for (int k = 0; k < 3; ++k) {
          diff[k] = sum[k] - contracted[3 * i + k];
        }
        worst = std::fmax(worst, mag3(diff) / mag3(&contracted[3 * i]));
      }
    }
    std::printf("B/A response matrices: row-sum vs contracted rel_diff=%.2e, "
                "par == serial: %s\n",
                worst, par_equal ? "yes" : "no");
    if (worst > 1e-10 || !par_equal) {
      std::printf("FAIL: response matrices disagree with contracted field\n");
      return 1;
    }
  }
  std::printf("OK\n");

  // Dense M: symmetric, positive diagonal, column sums equal the contracted
  // coupling, parallel matches serial.
  {
    std::vector<double> m(kSeg * kSeg), m_par(kSeg * kSeg), m_col(kSeg);
    int rc = cfsem::inductance_linear_filaments_matrix(
        seg_fil.data(), seg_dfil.data(), kSeg, seg_fil.data(), seg_dfil.data(),
        seg_wr.data(), kSeg, /*par=*/0, m.data());
    rc |= cfsem::inductance_linear_filaments_matrix(
        seg_fil.data(), seg_dfil.data(), kSeg, seg_fil.data(), seg_dfil.data(),
        seg_wr.data(), kSeg, /*par=*/1, m_par.data());
    rc |= cfsem::inductance_linear_filaments(
        seg_fil.data(), seg_dfil.data(), kSeg, seg_fil.data(), seg_dfil.data(),
        seg_wr.data(), kSeg, m_col.data());
    const int rc_overflow = cfsem::inductance_linear_filaments_matrix(
        seg_fil.data(), seg_dfil.data(), static_cast<std::size_t>(-1) / 2,
        seg_fil.data(), seg_dfil.data(), seg_wr.data(), 3, /*par=*/0, m.data());
    if (rc != 0 || rc_overflow != 5) {
      std::printf("FAIL: dense inductance rc=%d overflow rc=%d (expected 5)\n",
                  rc, rc_overflow);
      return 1;
    }
    double m_max = 0.0, asym = 0.0, col_err = 0.0;
    bool diag_pos = true;
    for (std::size_t j = 0; j < kSeg; ++j) {
      double col = 0.0;
      for (std::size_t i = 0; i < kSeg; ++i) {
        const double mij = m[i * kSeg + j];
        col += mij;
        m_max = std::fmax(m_max, std::fabs(mij));
        asym = std::fmax(asym, std::fabs(mij - m[j * kSeg + i]));
      }
      diag_pos = diag_pos && m[j * kSeg + j] > 0.0;
      col_err = std::fmax(col_err, std::fabs(col - m_col[j]) / std::fabs(m_col[j]));
    }
    std::printf("dense M (%zux%zu): asymmetry=%.2e, column-sum vs contracted "
                "rel_diff=%.2e, diag>0: %s, par == serial: %s\n",
                kSeg, kSeg, asym / m_max, col_err, diag_pos ? "yes" : "no",
                m == m_par ? "yes" : "no");
    if (asym / m_max > 1e-12 || col_err > 1e-10 || !diag_pos || m != m_par) {
      std::printf("FAIL: dense inductance checks\n");
      return 1;
    }
  }
  std::printf("OK\n");

  // Body force density must equal J x B with B from the direct evaluator.
  {
    std::vector<double> j_obs(3 * kObs), f(3 * kObs), f_par(3 * kObs),
        b(3 * kObs);
    for (std::size_t i = 0; i < kObs; ++i) {
      j_obs[3 * i] = 1.0e6 * std::sin(0.3 * static_cast<double>(i));
      j_obs[3 * i + 1] = 2.0e6;
      j_obs[3 * i + 2] = -1.0e6 * std::cos(0.2 * static_cast<double>(i));
    }
    int rc = cfsem::body_force_density_linear_filament(
        seg_fil.data(), seg_dfil.data(), seg_i.data(), seg_wr.data(), kSeg,
        obs.data(), j_obs.data(), kObs, /*par=*/0, f.data());
    rc |= cfsem::body_force_density_linear_filament(
        seg_fil.data(), seg_dfil.data(), seg_i.data(), seg_wr.data(), kSeg,
        obs.data(), j_obs.data(), kObs, /*par=*/1, f_par.data());
    rc |= cfsem::flux_density_linear_filament(
        obs.data(), kObs, seg_fil.data(), seg_dfil.data(), seg_i.data(),
        seg_wr.data(), kSeg, b.data());
    if (rc != 0) {
      std::printf("FAIL: body force rc=%d\n", rc);
      return 1;
    }
    double worst = 0.0;
    for (std::size_t i = 0; i < kObs; ++i) {
      const double* jv = &j_obs[3 * i];
      const double* bv = &b[3 * i];
      const double jxb[3] = {jv[1] * bv[2] - jv[2] * bv[1],
                             jv[2] * bv[0] - jv[0] * bv[2],
                             jv[0] * bv[1] - jv[1] * bv[0]};
      double diff[3];
      for (int k = 0; k < 3; ++k) {
        diff[k] = f[3 * i + k] - jxb[k];
      }
      worst = std::fmax(worst, mag3(diff) / mag3(jxb));
    }
    std::printf("body force vs J x B: rel_diff=%.2e, par == serial: %s\n",
                worst, f == f_par ? "yes" : "no");
    if (worst > 1e-10 || f != f_par) {
      std::printf("FAIL: body force disagrees with J x B\n");
      return 1;
    }
  }
  std::printf("OK\n");

  // --- circular filaments ---------------------------------------------------
  const double loop_r = 0.5, loop_z = 0.1, loop_i = 2.0e3; // (m, m, A-turns)
  const double rz_loop[2] = {loop_r, loop_z};
  const double r_loop[1] = {loop_r};
  const double i_loop[1] = {loop_i};
  const double thin[1] = {0.0};
  const double loc_loop[3] = {0.0, 0.0, loop_z};
  const double z_hat[3] = {0.0, 0.0, 1.0};

  // On axis: B_z = mu0 I R^2 / (2 (R^2 + dz^2)^1.5), B_r = 0. Slightly off axis,
  // B_r against the small-r series B_r = -(r/2) dB_z/dz, good to O(r^3).
  {
    const double dz = 0.3, r_near = 1.0e-3;
    const double s2 = loop_r * loop_r + dz * dz;
    const double bz_analytic =
        kMu0 * loop_i * loop_r * loop_r / (2.0 * std::pow(s2, 1.5));
    const double br_series = -0.5 * r_near * (-3.0 * dz * bz_analytic / s2);
    const double rz_axis[2] = {0.0, loop_z + dz};
    const double rz_near[2] = {r_near, loop_z + dz};
    double b_axis[2], b_near[2];
    int rc = cfsem::flux_density_circular_filament(rz_loop, i_loop, thin, 1,
                                                   rz_axis, 1, /*par=*/0, b_axis);
    rc |= cfsem::flux_density_circular_filament(rz_loop, i_loop, thin, 1,
                                                rz_near, 1, /*par=*/0, b_near);
    if (rc != 0) {
      std::printf("FAIL: circular B call\n");
      return 1;
    }
    const double rel_bz = std::fabs(b_axis[1] - bz_analytic) / bz_analytic;
    const double rel_br = std::fabs(b_near[0] - br_series) / br_series;
    std::printf("circular on-axis: Bz rel_err=%.2e, Br=%.1e; Br at r=%.0e vs "
                "series rel_err=%.2e\n",
                rel_bz, b_axis[0], r_near, rel_br);
    if (rel_bz > 1e-6 || b_axis[0] != 0.0 || rel_br > 1e-3) {
      std::printf("FAIL: circular near-axis field\n");
      return 1;
    }
  }
  std::printf("OK\n");

  // Off-axis grid: psi = 2 pi r A_phi; cartesian B = (Br cos phi, Br sin phi,
  // Bz) and A = A_phi (-sin phi, cos phi, 0); J x B matches; par == serial.
  {
    constexpr std::size_t kRz = 48;
    std::vector<double> rz(2 * kRz), xyz(3 * kRz), j(3 * kRz);
    for (std::size_t i = 0; i < kRz; ++i) {
      const double r_i = 0.05 + 1.2 * static_cast<double>(i) / kRz;
      const double z_i = -0.4 + 0.017 * static_cast<double>(i);
      const double phi = 0.13 * static_cast<double>(i);
      rz[2 * i] = r_i;
      rz[2 * i + 1] = z_i;
      xyz[3 * i] = r_i * std::cos(phi);
      xyz[3 * i + 1] = r_i * std::sin(phi);
      xyz[3 * i + 2] = z_i;
      j[3 * i] = 1.0e6;
      j[3 * i + 1] = -2.0e6 * std::cos(phi);
      j[3 * i + 2] = 0.5e6;
    }
    std::vector<double> psi(kRz), psi_p(kRz), aphi(kRz), aphi_p(kRz),
        b_rz(2 * kRz), b_rz_p(2 * kRz), b_xyz(3 * kRz), b_xyz_p(3 * kRz),
        a_xyz(3 * kRz), a_xyz_p(3 * kRz), f(3 * kRz), f_p(3 * kRz);
    int rc = 0;
    for (int par = 0; par < 2; ++par) {
      rc |= cfsem::flux_circular_filament(rz_loop, i_loop, thin, 1, rz.data(),
                                          kRz, par,
                                          par ? psi_p.data() : psi.data());
      rc |= cfsem::vector_potential_circular_filament(
          rz_loop, i_loop, thin, 1, rz.data(), kRz, par,
          par ? aphi_p.data() : aphi.data());
      rc |= cfsem::flux_density_circular_filament(
          rz_loop, i_loop, thin, 1, rz.data(), kRz, par,
          par ? b_rz_p.data() : b_rz.data());
      rc |= cfsem::flux_density_circular_filament_cartesian(
          r_loop, i_loop, loc_loop, z_hat, thin, 1, xyz.data(), kRz, par,
          par ? b_xyz_p.data() : b_xyz.data());
      rc |= cfsem::vector_potential_circular_filament_cartesian(
          r_loop, i_loop, loc_loop, z_hat, thin, 1, xyz.data(), kRz, par,
          par ? a_xyz_p.data() : a_xyz.data());
      rc |= cfsem::body_force_density_circular_filament(
          r_loop, i_loop, loc_loop, z_hat, thin, 1, xyz.data(), j.data(), kRz,
          par, par ? f_p.data() : f.data());
    }
    if (rc != 0) {
      std::printf("FAIL: circular grid calls rc=%d\n", rc);
      return 1;
    }
    double psi_err = 0.0, cart_err = 0.0, a_err = 0.0, f_err = 0.0;
    for (std::size_t i = 0; i < kRz; ++i) {
      const double r_i = rz[2 * i];
      psi_err = std::fmax(psi_err, std::fabs(psi[i] - 2.0 * kPi * r_i * aphi[i]) /
                                       std::fabs(psi[i]));
      const double c = xyz[3 * i] / r_i, s = xyz[3 * i + 1] / r_i;
      const double b_ref[3] = {b_rz[2 * i] * c, b_rz[2 * i] * s, b_rz[2 * i + 1]};
      const double a_ref[3] = {-aphi[i] * s, aphi[i] * c, 0.0};
      const double* bv = &b_xyz[3 * i];
      double diff[3];
      for (int k = 0; k < 3; ++k) {
        diff[k] = bv[k] - b_ref[k];
      }
      cart_err = std::fmax(cart_err, mag3(diff) / mag3(b_ref));
      for (int k = 0; k < 3; ++k) {
        diff[k] = a_xyz[3 * i + k] - a_ref[k];
      }
      a_err = std::fmax(a_err, mag3(diff) / mag3(a_ref));
      const double* jv = &j[3 * i];
      const double jxb[3] = {jv[1] * bv[2] - jv[2] * bv[1],
                             jv[2] * bv[0] - jv[0] * bv[2],
                             jv[0] * bv[1] - jv[1] * bv[0]};
      for (int k = 0; k < 3; ++k) {
        diff[k] = f[3 * i + k] - jxb[k];
      }
      f_err = std::fmax(f_err, mag3(diff) / mag3(jxb));
    }
    const bool par_equal = psi == psi_p && aphi == aphi_p && b_rz == b_rz_p &&
                           b_xyz == b_xyz_p && a_xyz == a_xyz_p && f == f_p;
    std::printf("circular grid: psi vs 2 pi r A_phi=%.2e, cartesian vs rotated "
                "B=%.2e A=%.2e, body force vs J x B=%.2e, par == serial: %s\n",
                psi_err, cart_err, a_err, f_err, par_equal ? "yes" : "no");
    if (psi_err > 1e-12 || cart_err > 1e-12 || a_err > 1e-12 || f_err > 1e-10 ||
        !par_equal) {
      std::printf("FAIL: circular grid consistency\n");
      return 1;
    }
  }
  std::printf("OK\n");

  // Pose: a loop centered at c with normal +x. Cyclic map local (x,y,z) -> world
  // (y,z,x) is a proper rotation taking the local axis to the normal, so world B
  // at p is (B_loc,z, B_loc,x, B_loc,y) of the coaxial field at local q, where
  // q = (p_y - c_y, p_z - c_z, p_x - c_x).
  {
    constexpr std::size_t kPose = 16;
    const double c[3] = {0.3, -0.2, 0.7};
    const double x_hat[3] = {1.0, 0.0, 0.0};
    const double rz_origin[2] = {loop_r, 0.0};
    std::vector<double> p(3 * kPose), b(3 * kPose), rz_q(2 * kPose),
        b_q(2 * kPose);
    std::vector<double> phi_q(kPose);
    for (std::size_t i = 0; i < kPose; ++i) {
      const double qr = 0.1 + 0.08 * static_cast<double>(i);
      const double qz = -0.3 + 0.05 * static_cast<double>(i);
      phi_q[i] = 0.4 * static_cast<double>(i);
      const double q[3] = {qr * std::cos(phi_q[i]), qr * std::sin(phi_q[i]), qz};
      p[3 * i] = c[0] + q[2];
      p[3 * i + 1] = c[1] + q[0];
      p[3 * i + 2] = c[2] + q[1];
      rz_q[2 * i] = qr;
      rz_q[2 * i + 1] = qz;
    }
    int rc = cfsem::flux_density_circular_filament_cartesian(
        r_loop, i_loop, c, x_hat, thin, 1, p.data(), kPose, /*par=*/0, b.data());
    rc |= cfsem::flux_density_circular_filament(rz_origin, i_loop, thin, 1,
                                                rz_q.data(), kPose, /*par=*/0,
                                                b_q.data());
    if (rc != 0) {
      std::printf("FAIL: circular pose calls rc=%d\n", rc);
      return 1;
    }
    double pose_err = 0.0;
    for (std::size_t i = 0; i < kPose; ++i) {
      const double br = b_q[2 * i], bz = b_q[2 * i + 1];
      const double b_loc[3] = {br * std::cos(phi_q[i]), br * std::sin(phi_q[i]), bz};
      const double b_ref[3] = {b_loc[2], b_loc[0], b_loc[1]};
      double diff[3];
      for (int k = 0; k < 3; ++k) {
        diff[k] = b[3 * i + k] - b_ref[k];
      }
      pose_err = std::fmax(pose_err, mag3(diff) / mag3(b_ref));
    }
    std::printf("circular pose (center %.1f,%.1f,%.1f, normal +x) vs coaxial "
                "rotated: %.2e\n",
                c[0], c[1], c[2], pose_err);
    if (pose_err > 1e-12) {
      std::printf("FAIL: circular pose\n");
      return 1;
    }
  }
  std::printf("OK\n");

  // Finite wire radius: identical to the ideal filament beyond the 3-radius
  // blend band, finite on the centerline where the ideal field diverges.
  {
    const double wr[1] = {0.01};
    const double rz_far[4] = {loop_r + 0.05, loop_z, loop_r, loop_z + 0.2};
    const double rz_on[2] = {loop_r, loop_z};
    double b_thin[4], b_wire[4], b_on[2];
    int rc = cfsem::flux_density_circular_filament(rz_loop, i_loop, thin, 1,
                                                   rz_far, 2, /*par=*/0, b_thin);
    rc |= cfsem::flux_density_circular_filament(rz_loop, i_loop, wr, 1, rz_far,
                                                2, /*par=*/0, b_wire);
    rc |= cfsem::flux_density_circular_filament(rz_loop, i_loop, wr, 1, rz_on,
                                                1, /*par=*/0, b_on);
    if (rc != 0) {
      std::printf("FAIL: circular wire radius calls rc=%d\n", rc);
      return 1;
    }
    double far_err = 0.0;
    for (int i = 0; i < 2; ++i) {
      const double d[2] = {b_wire[2 * i] - b_thin[2 * i],
                           b_wire[2 * i + 1] - b_thin[2 * i + 1]};
      far_err = std::fmax(far_err, std::hypot(d[0], d[1]) /
                                       std::hypot(b_thin[2 * i], b_thin[2 * i + 1]));
    }
    const bool on_finite = std::isfinite(b_on[0]) && std::isfinite(b_on[1]);
    std::printf("circular wire radius %.0e m: vs ideal beyond blend=%.2e, "
                "centerline |B|=%.3e T (finite: %s)\n",
                wr[0], far_err, std::hypot(b_on[0], b_on[1]),
                on_finite ? "yes" : "no");
    if (far_err > 1e-12 || !on_finite) {
      std::printf("FAIL: circular wire radius\n");
      return 1;
    }
  }
  std::printf("OK\n");

  // Mutual inductance loop -> coaxial polygon as sum of A(midpoint) . dl per
  // amp-turn. (Upstream mutual_inductance_circular_to_linear is ~0.8% low on
  // this case, so it is not exported.) Two independent references: the loop's
  // flux through the target radius (psi per amp-turn), and a dense
  // linear-filament M with the source loop discretized too. All converge as
  // O(1/N^2) in polygon resolution.
  {
    constexpr std::size_t kPoly = 720;
    const double tgt_r = 0.8, tgt_z = 0.45;
    auto polygon = [](double r, double z, std::vector<double>& start,
                      std::vector<double>& delta) {
      start.resize(3 * kPoly);
      delta.resize(3 * kPoly);
      for (std::size_t s = 0; s < kPoly; ++s) {
        const double a0 = 2.0 * kPi * static_cast<double>(s) / kPoly;
        const double a1 = 2.0 * kPi * static_cast<double>(s + 1) / kPoly;
        start[3 * s] = r * std::cos(a0);
        start[3 * s + 1] = r * std::sin(a0);
        start[3 * s + 2] = z;
        delta[3 * s] = r * (std::cos(a1) - std::cos(a0));
        delta[3 * s + 1] = r * (std::sin(a1) - std::sin(a0));
        delta[3 * s + 2] = 0.0;
      }
    };
    std::vector<double> tgt_start, tgt_delta, src_start, src_delta;
    polygon(tgt_r, tgt_z, tgt_start, tgt_delta);
    polygon(loop_r, loop_z, src_start, src_delta);

    const double one_turn[1] = {1.0};
    std::vector<double> mid(3 * kPoly), a_mid(3 * kPoly);
    for (std::size_t k = 0; k < 3 * kPoly; ++k) {
      mid[k] = tgt_start[k] + 0.5 * tgt_delta[k];
    }
    double psi_tgt = 0.0;
    const double rz_tgt[2] = {tgt_r, tgt_z};
    int rc = cfsem::vector_potential_circular_filament_cartesian(
        r_loop, one_turn, loc_loop, z_hat, thin, 1, mid.data(), kPoly,
        /*par=*/0, a_mid.data());
    rc |= cfsem::flux_circular_filament(rz_loop, one_turn, thin, 1, rz_tgt, 1,
                                        /*par=*/0, &psi_tgt);

    const std::vector<double> src_wr(kPoly, 1.0e-4);
    std::vector<double> m_col(kPoly);
    rc |= cfsem::inductance_linear_filaments(
        tgt_start.data(), tgt_delta.data(), kPoly, src_start.data(),
        src_delta.data(), src_wr.data(), kPoly, m_col.data());
    if (rc != 0) {
      std::printf("FAIL: circular-to-linear calls rc=%d\n", rc);
      return 1;
    }
    double m_lin = 0.0, m_circ = 0.0;
    for (double v : m_col) {
      m_lin += v;
    }
    for (std::size_t k = 0; k < 3 * kPoly; ++k) {
      m_circ += a_mid[k] * tgt_delta[k];
    }
    const double rel_psi = std::fabs(m_circ - psi_tgt) / psi_tgt;
    const double rel_lin = std::fabs(m_circ - m_lin) / m_lin;
    std::printf("circular->linear M=%.6e H; vs psi/I=%.6e (rel %.2e), vs "
                "polygon-polygon M=%.6e (rel %.2e)\n",
                m_circ, psi_tgt, rel_psi, m_lin, rel_lin);
    if (rel_psi > 1e-4 || rel_lin > 1e-4) {
      std::printf("FAIL: circular-to-linear inductance\n");
      return 1;
    }
  }
  std::printf("OK\n");

  // Triangle mesh: faceted cylinder carrying an azimuthal current sheet, s linear
  // in z, so the total current is s_top - s_bottom.
  {
    const double rad = 0.1, len = 0.2, i_tot = 1.0e3; // (m), (m), (A)
    struct Mesh {
      std::vector<double> nodes, s;
      std::vector<std::int64_t> tris;
      std::size_t n_node() const { return s.size(); }
      std::size_t n_tri() const { return tris.size() / 3; }
    };
    auto cylinder = [&](std::size_t nphi, std::size_t nz) {
      Mesh m;
      for (std::size_t j = 0; j <= nz; ++j) {
        const double z = -0.5 * len + len * double(j) / double(nz);
        for (std::size_t k = 0; k < nphi; ++k) {
          const double phi = 2.0 * kPi * double(k) / double(nphi);
          m.nodes.insert(m.nodes.end(),
                         {rad * std::cos(phi), rad * std::sin(phi), z});
          m.s.push_back(i_tot * double(j) / double(nz));
        }
      }
      for (std::size_t j = 0; j < nz; ++j) {
        for (std::size_t k = 0; k < nphi; ++k) {
          const auto n = [&](std::size_t jj, std::size_t kk) {
            return std::int64_t(jj * nphi + kk % nphi);
          };
          m.tris.insert(m.tris.end(), {n(j, k), n(j, k + 1), n(j + 1, k + 1)});
          m.tris.insert(m.tris.end(), {n(j, k), n(j + 1, k + 1), n(j + 1, k)});
        }
      }
      return m;
    };

    // Center field vs the finite solenoid; the inscribed facets shrink the
    // effective radius, an O((pi/nphi)^2) bias.
    const Mesh fine = cylinder(128, 32);
    const double center[3] = {0.0, 0.0, 0.0};
    double b_c[3];
    int rc = cfsem::flux_density_triangle_mesh(
        fine.nodes.data(), fine.n_node(), fine.tris.data(), fine.n_tri(),
        fine.s.data(), center, 1, /*par=*/0, b_c);
    const double bz_exact = kMu0 * (i_tot / len) * (0.5 * len) /
                            std::sqrt(0.25 * len * len + rad * rad);
    const double rel_bz = std::fabs(std::fabs(b_c[2]) - bz_exact) / bz_exact;

    // curl A vs B by central differences at an interior point.
    const double p[3] = {0.03, 0.02, 0.01}, h = 1.0e-5;
    double probe[18], a_p[18], b_p[3];
    for (int d = 0; d < 3; ++d) {
      for (int sgn = 0; sgn < 2; ++sgn) {
        double* q = probe + 3 * (2 * d + sgn);
        q[0] = p[0], q[1] = p[1], q[2] = p[2];
        q[d] += sgn ? h : -h;
      }
    }
    rc |= cfsem::vector_potential_triangle_mesh(
        fine.nodes.data(), fine.n_node(), fine.tris.data(), fine.n_tri(),
        fine.s.data(), probe, 6, /*par=*/1, a_p);
    rc |= cfsem::flux_density_triangle_mesh(
        fine.nodes.data(), fine.n_node(), fine.tris.data(), fine.n_tri(),
        fine.s.data(), p, 1, /*par=*/1, b_p);
    // dA_i/dx_d from probes (2d, 2d+1)
    const auto da = [&](int i, int d) {
      return (a_p[3 * (2 * d + 1) + i] - a_p[3 * (2 * d) + i]) / (2.0 * h);
    };
    const double curl[3] = {da(2, 1) - da(1, 2), da(0, 2) - da(2, 0),
                            da(1, 0) - da(0, 1)};
    const double rel_curl =
        std::sqrt(std::pow(curl[0] - b_p[0], 2) + std::pow(curl[1] - b_p[1], 2) +
                  std::pow(curl[2] - b_p[2], 2)) /
        mag3(b_p);

    // Hierarchical vs direct over a line through the wall, plus the near/far
    // partition and serial/parallel equality.
    constexpr std::size_t kLine = 64;
    std::vector<double> line(3 * kLine);
    for (std::size_t k = 0; k < kLine; ++k) {
      line[3 * k] = 0.25 * rad + 1.5 * rad * double(k) / double(kLine);
      line[3 * k + 1] = 0.1 * rad;
      line[3 * k + 2] = 0.3 * len * double(k) / double(kLine) - 0.1 * len;
    }
    std::vector<double> b_dir(3 * kLine), b_ser(3 * kLine), b_hier(3 * kLine),
        b_near(3 * kLine), b_far(3 * kLine), a_dir(3 * kLine),
        a_hier(3 * kLine);
    const double theta = 0.5;
    const auto* nd = fine.nodes.data();
    const auto* tr = fine.tris.data();
    const auto* sv = fine.s.data();
    const std::size_t nn = fine.n_node(), nt = fine.n_tri();
    rc |= cfsem::flux_density_triangle_mesh(nd, nn, tr, nt, sv, line.data(),
                                            kLine, 1, b_dir.data());
    rc |= cfsem::flux_density_triangle_mesh(nd, nn, tr, nt, sv, line.data(),
                                            kLine, 0, b_ser.data());
    rc |= cfsem::vector_potential_triangle_mesh(nd, nn, tr, nt, sv, line.data(),
                                                kLine, 1, a_dir.data());
    rc |= cfsem::flux_density_triangle_mesh_hierarchical(
        nd, nn, tr, nt, sv, line.data(), kLine, theta, 1, CFSEM_SKIP_NONE,
        b_hier.data());
    rc |= cfsem::flux_density_triangle_mesh_hierarchical(
        nd, nn, tr, nt, sv, line.data(), kLine, theta, 1, CFSEM_SKIP_FAR,
        b_near.data());
    rc |= cfsem::flux_density_triangle_mesh_hierarchical(
        nd, nn, tr, nt, sv, line.data(), kLine, theta, 1, CFSEM_SKIP_NEAR,
        b_far.data());
    rc |= cfsem::vector_potential_triangle_mesh_hierarchical(
        nd, nn, tr, nt, sv, line.data(), kLine, theta, 1, CFSEM_SKIP_NONE,
        a_hier.data());
    // max-norm error relative to the max field on the line
    const auto rel_max = [](const std::vector<double>& x,
                            const std::vector<double>& ref) {
      double num = 0.0, den = 0.0;
      for (std::size_t k = 0; k < x.size(); ++k) {
        num = std::fmax(num, std::fabs(x[k] - ref[k]));
        den = std::fmax(den, std::fabs(ref[k]));
      }
      return num / den;
    };
    std::vector<double> b_sum(3 * kLine);
    for (std::size_t k = 0; k < b_sum.size(); ++k) {
      b_sum[k] = b_near[k] + b_far[k];
    }
    const double rel_hier_b = rel_max(b_hier, b_dir);
    const double rel_hier_a = rel_max(a_hier, a_dir);
    const double rel_part = rel_max(b_sum, b_hier);
    const bool par_eq = b_ser == b_dir;

    if (rc != 0) {
      std::printf("FAIL: triangle-mesh field calls rc=%d\n", rc);
      return 1;
    }
    std::printf("triangle mesh (%zu tris): center |Bz| rel_err=%.2e, curl A vs "
                "B=%.2e, par == serial: %s\n",
                nt, rel_bz, rel_curl, par_eq ? "yes" : "no");
    std::printf("triangle mesh hierarchical (theta=%.2f): B rel=%.2e, A "
                "rel=%.2e, near+far vs full=%.2e\n",
                theta, rel_hier_b, rel_hier_a, rel_part);
    if (rel_bz > 1e-3 || rel_curl > 1e-6 || !par_eq || rel_hier_b > 1e-2 ||
        rel_hier_a > 1e-2 || rel_part > 1e-12) {
      std::printf("FAIL: triangle-mesh fields\n");
      return 1;
    }
    std::printf("OK\n");

    // Inductance: s^T L s / 2 vs 1/2 int K.A dS, A from the mesh's own field at
    // the degree-2 rule (2/3, 1/6, 1/6); K = n x grad s = -sum_i s_i e_i / 2area,
    // e_i the edge opposite node i taken counterclockwise about n.
    const Mesh coarse = cylinder(32, 8);
    const std::size_t cn = coarse.n_node(), ct = coarse.n_tri();
    std::vector<double> lmat(cn * cn), lmat_ser(cn * cn);
    rc = cfsem::triangle_mesh_inductance_matrix(coarse.nodes.data(), cn,
                                                coarse.tris.data(), ct, 5, 1,
                                                lmat.data());
    rc |= cfsem::triangle_mesh_inductance_matrix(coarse.nodes.data(), cn,
                                                 coarse.tris.data(), ct, 5, 0,
                                                 lmat_ser.data());
    double asym = 0.0, l_max = 0.0, e_l = 0.0;
    for (std::size_t a = 0; a < cn; ++a) {
      for (std::size_t b = 0; b < cn; ++b) {
        asym = std::fmax(asym, std::fabs(lmat[a * cn + b] - lmat[b * cn + a]));
        l_max = std::fmax(l_max, std::fabs(lmat[a * cn + b]));
        e_l += 0.5 * coarse.s[a] * lmat[a * cn + b] * coarse.s[b];
      }
    }
    double lmat_par_diff = 0.0;
    for (std::size_t k = 0; k < lmat.size(); ++k) {
      lmat_par_diff = std::fmax(lmat_par_diff, std::fabs(lmat[k] - lmat_ser[k]));
    }

    std::vector<double> qp(9 * ct), a_q(9 * ct), k_t(3 * ct), area(ct);
    const double bary[3][3] = {{2.0 / 3, 1.0 / 6, 1.0 / 6},
                               {1.0 / 6, 2.0 / 3, 1.0 / 6},
                               {1.0 / 6, 1.0 / 6, 2.0 / 3}};
    for (std::size_t t = 0; t < ct; ++t) {
      const double* v[3];
      double sv3[3];
      for (int i = 0; i < 3; ++i) {
        const auto idx = std::size_t(coarse.tris[3 * t + i]);
        v[i] = &coarse.nodes[3 * idx];
        sv3[i] = coarse.s[idx];
      }
      const double e1[3] = {v[1][0] - v[0][0], v[1][1] - v[0][1],
                            v[1][2] - v[0][2]};
      const double e2[3] = {v[2][0] - v[0][0], v[2][1] - v[0][1],
                            v[2][2] - v[0][2]};
      const double cr[3] = {e1[1] * e2[2] - e1[2] * e2[1],
                            e1[2] * e2[0] - e1[0] * e2[2],
                            e1[0] * e2[1] - e1[1] * e2[0]};
      const double twice_area = mag3(cr);
      area[t] = 0.5 * twice_area;
      for (int c = 0; c < 3; ++c) {
        double kc = 0.0;
        for (int i = 0; i < 3; ++i) {
          kc -= sv3[i] * (v[(i + 2) % 3][c] - v[(i + 1) % 3][c]);
        }
        k_t[3 * t + c] = kc / twice_area;
        for (int q = 0; q < 3; ++q) {
          qp[9 * t + 3 * q + c] = bary[q][0] * v[0][c] + bary[q][1] * v[1][c] +
                                  bary[q][2] * v[2][c];
        }
      }
    }
    rc |= cfsem::vector_potential_triangle_mesh(
        coarse.nodes.data(), cn, coarse.tris.data(), ct, coarse.s.data(),
        qp.data(), 3 * ct, /*par=*/1, a_q.data());
    double e_a = 0.0;
    for (std::size_t t = 0; t < ct; ++t) {
      for (int q = 0; q < 3; ++q) {
        for (int c = 0; c < 3; ++c) {
          e_a += 0.5 * (area[t] / 3.0) * k_t[3 * t + c] * a_q[9 * t + 3 * q + c];
        }
      }
    }

    // Status contract: bad quadrature rule (5), out-of-range node index (1).
    std::vector<std::int64_t> bad = coarse.tris;
    bad[0] = std::int64_t(cn);
    double b_bad[3];
    const int rc_quad = cfsem::triangle_mesh_inductance_matrix(
        coarse.nodes.data(), cn, coarse.tris.data(), ct, 0, 0, lmat.data());
    const int rc_idx = cfsem::flux_density_triangle_mesh(
        coarse.nodes.data(), cn, bad.data(), ct, coarse.s.data(), center, 1, 0,
        b_bad);

    if (rc != 0) {
      std::printf("FAIL: triangle-mesh inductance calls rc=%d\n", rc);
      return 1;
    }
    const double rel_e = std::fabs(e_l - e_a) / std::fabs(e_l);
    std::printf("triangle mesh L (%zux%zu): asymmetry=%.2e, par vs serial=%.2e, "
                "energy s^T L s/2=%.6e J vs 1/2 int K.A=%.6e J (rel %.2e); "
                "rc bad quad=%d, bad index=%d\n",
                cn, cn, asym / l_max, lmat_par_diff / l_max, e_l, e_a, rel_e,
                rc_quad, rc_idx);
    if (asym / l_max > 1e-14 || lmat_par_diff / l_max > 1e-12 || rel_e > 1e-2 ||
        rc_quad != 5 || rc_idx != 1) {
      std::printf("FAIL: triangle-mesh inductance\n");
      return 1;
    }
  }
  std::printf("OK\n");
  return 0;
}
