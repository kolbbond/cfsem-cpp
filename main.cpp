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
  return 0;
}
