// Smoke test for the cfsem C-ABI shim library.
//
// This is NOT the deliverable — it just proves the linkable `cfsem_cpp` library
// builds, links, and runs. The deliverable is the library + headers that goose
// (and other C/C++ code) links against. Returns nonzero if any cfsem call fails.

#include <cstddef>
#include <cstdio>

#include "cfsem.hpp"

int main() {
  std::printf("Initialize cfsem C-ABI shim\n");

  // 1) Scalar path: complete elliptic integral of the first kind.
  std::printf("cfsem::ellipk(0.1) = %.6f\n", cfsem::ellipk(0.1));

  // 2) Array path: Biot-Savart B-field from a single straight filament.
  //    Segment from (0,0,0) to (1,0,0) m carrying 1000 A; observe 0.1 m off the
  //    wire at its midpoint. Expect a B-field of order 1e-3 T pointing in -z.
  const double xp[] = {0.5};
  const double yp[] = {0.1};
  const double zp[] = {0.0};
  const std::size_t n_obs = 1;

  const double xfil[] = {0.0};
  const double yfil[] = {0.0};
  const double zfil[] = {0.0};
  const double dlx[] = {1.0};
  const double dly[] = {0.0};
  const double dlz[] = {0.0};
  const double ifil[] = {1000.0};
  const double wire_radius[] = {0.001};
  const std::size_t n_fil = 1;

  double bx[1] = {0.0};
  double by[1] = {0.0};
  double bz[1] = {0.0};

  const int rc = cfsem::flux_density_linear_filament(
      xp, yp, zp, n_obs, xfil, yfil, zfil, dlx, dly, dlz, ifil, wire_radius,
      n_fil, bx, by, bz);

  std::printf("cfsem::flux_density_linear_filament rc=%d  B=(%.6e, %.6e, %.6e) T\n",
              rc, bx[0], by[0], bz[0]);

  if (rc != 0) {
    std::printf("FAIL: cfsem call returned %d\n", rc);
    return rc;
  }
  std::printf("OK\n");
  return 0;
}
