/*
 * (C) Copyright 2026- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <netcdf.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "eckit/testing/Test.h"

#include "oops/mpi/mpi.h"
#include "oops/runs/Run.h"
#include "oops/runs/Test.h"

#include "soca/Geometry/VerticalCoordinate.h"

namespace test {

// -----------------------------------------------------------------------------

void expectSame(const std::vector<double> & a, const std::vector<double> & b) {
  EXPECT_EQUAL(a.size(), b.size());
  for (size_t k = 0; k < a.size(); ++k) EXPECT(std::fabs(a[k] - b[k]) < 1.0e-12);
}

const std::vector<double> kDepths{1.0, 3.0, 5.0, 7.0, 9.0, 11.0};

// -----------------------------------------------------------------------------

void testLevels() {
  expectSame(soca::verticalCoordinate("levels", 4, {}, 10.0), {1.0, 2.0, 3.0, 4.0});
}

void testMetres() {
  expectSame(soca::verticalCoordinate("m", kDepths.size(), kDepths, 10.0), kDepths);
}

void testLogDepth() {
  const auto c = soca::verticalCoordinate("log depth", kDepths.size(), kDepths, 10.0);
  EXPECT_EQUAL(c.size(), kDepths.size());
  for (size_t k = 0; k < c.size(); ++k) {
    EXPECT(std::fabs(c[k] - std::log(kDepths[k] + 10.0)) < 1.0e-12);
  }
}

void testCoordinateErrors() {
  EXPECT_THROWS(soca::verticalCoordinate("pressure", 6, kDepths, 10.0));
  EXPECT_THROWS(soca::verticalCoordinate("m", 6, {}, 10.0));
  EXPECT_THROWS(soca::verticalCoordinate("log depth", 6, {}, 10.0));
  EXPECT_THROWS(soca::verticalCoordinate("m", 7, kDepths, 10.0));
  EXPECT_THROWS(soca::verticalCoordinate("m", 3, {1.0, 3.0, 3.0}, 10.0));
  EXPECT_THROWS(soca::verticalCoordinate("m", 3, {-1.0, 3.0, 5.0}, 10.0));
  EXPECT_THROWS(soca::verticalCoordinate("log depth", 6, kDepths, 0.0));
  EXPECT_NO_THROW(soca::verticalCoordinate("m", 6, kDepths, 0.0));  // offset unused
}

// -----------------------------------------------------------------------------

// z* layers 2, 2, 4 m thick (centres 1, 3, 6 m), 2 x 3 columns: 0-2 full depth, 3 partial
// bottom cell, 4 shallow (layer 3 vanished, 0), 5 land (`dry`). Partial cells are left out
// above the bottom layer; bottom layer: median(6, 6, 6, 4.5) = 6.
std::vector<std::vector<double>> thickness(const double dry) {
  return {{2.0, 2.0, 2.0, 2.0, 2.0, dry},
          {2.0, 2.0, 2.0, 2.0, 1.0, dry},
          {4.0, 4.0, 4.0, 1.0, 0.0, dry}};
}

// h([time,] layer, y, x) with optional _FillValue / missing_value, written on rank 0.
void writeFile(const std::string & name, const std::vector<std::vector<double>> & h,
               const nc_type type, const bool withTime, const double * fill) {
  const eckit::mpi::Comm & comm = oops::mpi::world();
  if (comm.rank() == 0) {
    int ncid, dims[4], varid;
    nc_create(name.c_str(), NC_CLOBBER, &ncid);
    int nd = 0;
    if (withTime) nc_def_dim(ncid, "time", NC_UNLIMITED, &dims[nd++]);
    nc_def_dim(ncid, "z_l", h.size(), &dims[nd++]);
    nc_def_dim(ncid, "yh", 2, &dims[nd++]);
    nc_def_dim(ncid, "xh", 3, &dims[nd++]);
    nc_def_var(ncid, "h", type, nd, dims, &varid);
    if (fill) {
      nc_put_att_double(ncid, varid, "_FillValue", type, 1, fill);
      nc_put_att_double(ncid, varid, "missing_value", type, 1, fill);
    }
    nc_enddef(ncid);
    for (size_t k = 0; k < h.size(); ++k) {
      const size_t start[4] = {0, k, 0, 0};
      const size_t count[4] = {1, 1, 2, 3};
      nc_put_vara_double(ncid, varid, withTime ? start : start + 1,
                         withTime ? count : count + 1, h[k].data());
    }
    nc_close(ncid);
  }
  comm.barrier();
}

std::vector<double> depths(const std::string & file, const std::string & var,
                           const size_t nlevs) {
  return soca::layerDepthsFromFile(file, var, nlevs, oops::mpi::world());
}

void testFile() {
  const std::string f = "test_verticalcoord_h.nc";
  const double zero = 0.0, big = 1.0e20;

  // MOM6 history: float, with time, fill 0
  writeFile(f, thickness(0.0), NC_FLOAT, true, &zero);
  expectSame(depths(f, "h", 3), {1.0, 3.0, 6.0});
  EXPECT_THROWS(depths(f, "h", 4));          // wrong layer count
  EXPECT_THROWS(depths(f, "thickness", 3));  // no such variable

  // fill 1e20, no time
  writeFile(f, thickness(1.0e20), NC_DOUBLE, false, &big);
  expectSame(depths(f, "h", 3), {1.0, 3.0, 6.0});

  // no fill attributes: the netCDF default fill is missing
  writeFile(f, thickness(NC_FILL_FLOAT), NC_FLOAT, false, nullptr);
  expectSame(depths(f, "h", 3), {1.0, 3.0, 6.0});

  // a layer with no wet column
  auto h = thickness(0.0);
  h[1].assign(6, 0.0);
  writeFile(f, h, NC_DOUBLE, false, &zero);
  EXPECT_THROWS(depths(f, "h", 3));

  if (oops::mpi::world().rank() == 0) std::remove(f.c_str());
  EXPECT_THROWS(depths("no_such_file.nc", "h", 3));
}

// -----------------------------------------------------------------------------

class VerticalCoordinateTest : public oops::Test {
 private:
  std::string testid() const override {return "test::VerticalCoordinateTest";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();
    ts.emplace_back(CASE("soca/VerticalCoordinate/levels") {testLevels();});
    ts.emplace_back(CASE("soca/VerticalCoordinate/m") {testMetres();});
    ts.emplace_back(CASE("soca/VerticalCoordinate/logDepth") {testLogDepth();});
    ts.emplace_back(CASE("soca/VerticalCoordinate/errors") {testCoordinateErrors();});
    ts.emplace_back(CASE("soca/VerticalCoordinate/layerDepthsFromFile") {testFile();});
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::VerticalCoordinateTest tests;
  return run.execute(tests);
}
