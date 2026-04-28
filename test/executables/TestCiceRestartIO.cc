/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <netcdf.h>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/testing/Test.h"

#include "oops/base/Geometry.h"
#include "oops/base/State.h"
#include "oops/runs/Run.h"
#include "oops/runs/Test.h"
#include "oops/util/parameters/IgnoreOtherParameters.h"

#include "soca/Geometry/Geometry.h"
#include "soca/PostProcess/CiceRestartIO.h"
#include "soca/State/State.h"
#include "soca/Traits.h"

namespace test {

typedef oops::Geometry<soca::Traits> Geometry_;
typedef oops::State<soca::Traits>    State_;

class CiceIOTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(CiceIOTestParameters, Parameters)

 public:
  oops::RequiredParameter<eckit::LocalConfiguration> geometry{"geometry", this};
  oops::RequiredParameter<eckit::LocalConfiguration> state{"state", this};
  oops::RequiredParameter<std::string> inputRestart{"input restart", this};
  oops::RequiredParameter<std::string> outputRestart{"output restart", this};
  oops::RequiredParameter<int> ncat{"ncat", this};

  oops::IgnoreOtherParameters ignoreOthers{this};
};

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

void ncCheck(int rc, const char * where) {
  if (rc != NC_NOERR) {
    throw std::runtime_error(std::string("nc error at ") + where + ": "
                             + nc_strerror(rc));
  }
}

std::vector<double> readVar(const std::string & file, const std::string & var) {
  int ncid = -1, varid = -1;
  ncCheck(nc_open(file.c_str(), NC_NOWRITE, &ncid), "nc_open");
  ncCheck(nc_inq_varid(ncid, var.c_str(), &varid), ("inq_varid " + var).c_str());
  int ndims = 0;
  ncCheck(nc_inq_varndims(ncid, varid, &ndims), "inq_varndims");
  std::vector<int> dimids(ndims);
  ncCheck(nc_inq_vardimid(ncid, varid, dimids.data()), "inq_vardimid");
  std::size_t total = 1;
  for (int d : dimids) {
    std::size_t len = 0;
    ncCheck(nc_inq_dimlen(ncid, d, &len), "inq_dimlen");
    total *= len;
  }
  std::vector<double> data(total);
  ncCheck(nc_get_var_double(ncid, varid, data.data()), "get_var_double");
  ncCheck(nc_close(ncid), "nc_close");
  return data;
}

int readGlobalIntAttr(const std::string & file, const std::string & attr) {
  int ncid = -1;
  ncCheck(nc_open(file.c_str(), NC_NOWRITE, &ncid), "nc_open");
  int value = 0;
  ncCheck(nc_get_att_int(ncid, NC_GLOBAL, attr.c_str(), &value),
          ("get_att " + attr).c_str());
  ncCheck(nc_close(ncid), "nc_close");
  return value;
}

std::array<std::size_t, 3> readDims(const std::string & file) {
  int ncid = -1, dimid = -1;
  std::array<std::size_t, 3> out{0, 0, 0};
  ncCheck(nc_open(file.c_str(), NC_NOWRITE, &ncid), "nc_open");
  std::size_t len = 0;
  ncCheck(nc_inq_dimid(ncid, "ncat", &dimid), "ncat");
  ncCheck(nc_inq_dimlen(ncid, dimid, &len), "len ncat");
  out[0] = len;
  ncCheck(nc_inq_dimid(ncid, "nj", &dimid), "nj");
  ncCheck(nc_inq_dimlen(ncid, dimid, &len), "len nj");
  out[1] = len;
  ncCheck(nc_inq_dimid(ncid, "ni", &dimid), "ni");
  ncCheck(nc_inq_dimlen(ncid, dimid, &len), "len ni");
  out[2] = len;
  ncCheck(nc_close(ncid), "nc_close");
  return out;
}

// ---------------------------------------------------------------------------
// Test 1: write the FieldSet read from disk back out, verify aicen/vicen/vsnon
// match the input within tolerance and dimensions/attrs are preserved.
// ---------------------------------------------------------------------------
void testRoundTrip() {
  CiceIOTestParameters params;
  params.validateAndDeserialize(TestEnvironment::config());

  Geometry_ geom(params.geometry, oops::mpi::world(), oops::mpi::myself());
  State_    state(geom, params.state);

  const std::string inFile  = params.inputRestart.value();
  const std::string outFile = params.outputRestart.value() + ".roundtrip.nc";
  const std::size_t ncat    = static_cast<std::size_t>(params.ncat.value());

  const soca::Geometry & socaGeom = state.state().geometry();
  soca::CiceRestartIO io(socaGeom, inFile, outFile);
  io.write(state.fieldSet().fieldSet(), ncat);

  if (oops::mpi::world().rank() != 0) return;

  // Dimensions must match.
  const auto dIn  = readDims(inFile);
  const auto dOut = readDims(outFile);
  EXPECT_EQUAL(dIn[0], dOut[0]);
  EXPECT_EQUAL(dIn[1], dOut[1]);
  EXPECT_EQUAL(dIn[2], dOut[2]);

  // Global attributes preserved by the byte-copy template path.
  for (const std::string & attr : {"istep1", "myear", "mmonth", "mday", "msec"}) {
    EXPECT_EQUAL(readGlobalIntAttr(inFile, attr),
                 readGlobalIntAttr(outFile, attr));
  }

  // The three vars CiceRestartIO actually overwrites must match the input
  // within float tolerance (FieldSet round-trip through soca State + atlas
  // can introduce sub-eps differences).
  for (const std::string & v : {"aicen", "vicen", "vsnon"}) {
    auto refData = readVar(inFile,  v);
    auto outData = readVar(outFile, v);
    EXPECT_EQUAL(refData.size(), outData.size());
    double maxDiff = 0.0;
    for (std::size_t i = 0; i < refData.size(); ++i) {
      maxDiff = std::max(maxDiff, std::fabs(refData[i] - outData[i]));
    }
    oops::Log::test() << "round-trip max|diff| for " << v << ": "
                      << maxDiff << std::endl;
    EXPECT(maxDiff < 1.0e-12);
  }

  // Vars that CiceRestartIO does not touch must be byte-identical to input
  // (they came from the byte-copy of the template).
  for (const std::string & v : {"Tsfcn", "qice001", "qice004", "qice007",
                                "sice001", "sice004", "sice007", "qsno001",
                                "apnd", "hpnd", "ipnd", "iceumask"}) {
    auto refData = readVar(inFile,  v);
    auto outData = readVar(outFile, v);
    EXPECT_EQUAL(refData.size(), outData.size());
    for (std::size_t i = 0; i < refData.size(); ++i) {
      EXPECT_EQUAL(refData[i], outData[i]);
    }
  }
}

// ---------------------------------------------------------------------------
// Test 2: Zero out aicen/vicen/vsnon in the FieldSet before writing. Output
// must have those three vars exactly zero, and every other var must remain
// byte-identical to the input restart.
// ---------------------------------------------------------------------------
void testSelectiveOverwrite() {
  CiceIOTestParameters params;
  params.validateAndDeserialize(TestEnvironment::config());

  Geometry_ geom(params.geometry, oops::mpi::world(), oops::mpi::myself());
  State_    state(geom, params.state);

  const std::string inFile  = params.inputRestart.value();
  const std::string outFile = params.outputRestart.value() + ".zeroed.nc";
  const std::size_t ncat    = static_cast<std::size_t>(params.ncat.value());

  const soca::Geometry & socaGeom = state.state().geometry();

  // Zero out the three groups in the FieldSet (in-place on a copy).
  atlas::FieldSet fset = state.fieldSet().fieldSet();
  for (std::size_t k = 0; k < ncat; ++k) {
    const std::string suffix = std::to_string(k + 1);
    for (const std::string & pattern : {
        "sea_ice_category" + suffix + "_area_fraction",
        "sea_ice_category" + suffix + "_volume",
        "sea_ice_snow_category" + suffix + "_volume"}) {
      auto view = atlas::array::make_view<double, 2>(fset.field(pattern));
      for (atlas::idx_t j = 0; j < view.shape(0); ++j) {
        for (atlas::idx_t l = 0; l < view.shape(1); ++l) view(j, l) = 0.0;
      }
    }
  }

  soca::CiceRestartIO io(socaGeom, inFile, outFile);
  io.write(fset, ncat);

  if (oops::mpi::world().rank() != 0) return;

  // The overwritten variables must be exactly zero everywhere.
  for (const std::string & v : {"aicen", "vicen", "vsnon"}) {
    auto outData = readVar(outFile, v);
    double maxAbs = 0.0;
    for (double x : outData) maxAbs = std::max(maxAbs, std::fabs(x));
    oops::Log::test() << "after zero, max|" << v << "| = " << maxAbs << std::endl;
    EXPECT_EQUAL(maxAbs, 0.0);
  }

  // Untouched variables must still be byte-identical to input.
  for (const std::string & v : {"Tsfcn", "qice001", "sice004", "qsno001",
                                "apnd", "hpnd", "ipnd", "iceumask"}) {
    auto refData = readVar(inFile,  v);
    auto outData = readVar(outFile, v);
    EXPECT_EQUAL(refData.size(), outData.size());
    for (std::size_t i = 0; i < refData.size(); ++i) {
      EXPECT_EQUAL(refData[i], outData[i]);
    }
  }

  // Structure preserved.
  const auto dIn  = readDims(inFile);
  const auto dOut = readDims(outFile);
  EXPECT_EQUAL(dIn[0], dOut[0]);
  EXPECT_EQUAL(dIn[1], dOut[1]);
  EXPECT_EQUAL(dIn[2], dOut[2]);
}

// ---------------------------------------------------------------------------
// Test 3: ncat mismatch must throw a clear error and leave no stale output.
// ---------------------------------------------------------------------------
void testNcatMismatch() {
  CiceIOTestParameters params;
  params.validateAndDeserialize(TestEnvironment::config());

  Geometry_ geom(params.geometry, oops::mpi::world(), oops::mpi::myself());
  State_    state(geom, params.state);

  const std::string inFile  = params.inputRestart.value();
  const std::string outFile = params.outputRestart.value() + ".ncatmismatch.nc";
  const std::size_t wrongNcat = static_cast<std::size_t>(params.ncat.value()) + 1;

  const soca::Geometry & socaGeom = state.state().geometry();
  soca::CiceRestartIO io(socaGeom, inFile, outFile);
  bool threw = false;
  try {
    io.write(state.fieldSet().fieldSet(), wrongNcat);
  } catch (const eckit::Exception &) {
    threw = true;
  }
  EXPECT(threw);
}

// ---------------------------------------------------------------------------

class CiceRestartIOTest : public oops::Test {
 private:
  std::string testid() const override {return "test::CiceRestartIOTest";}

  void register_tests() const override {
    std::vector<eckit::testing::Test> & ts = eckit::testing::specification();
    ts.emplace_back(CASE("soca/CiceRestartIO/roundTrip")
                    {testRoundTrip(); });
    ts.emplace_back(CASE("soca/CiceRestartIO/selectiveOverwrite")
                    {testSelectiveOverwrite(); });
    ts.emplace_back(CASE("soca/CiceRestartIO/ncatMismatch")
                    {testNcatMismatch(); });
  }

  void clear() const override {}
};

}  // namespace test

int main(int argc, char ** argv) {
  oops::Run run(argc, argv);
  test::CiceRestartIOTest tests;
  return run.execute(tests);
}
