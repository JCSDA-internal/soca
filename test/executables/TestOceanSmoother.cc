/*
* (C) Copyleft 2024-2024 UCAR
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/testing/Test.h"

#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Geometry.h"
#include "oops/runs/Run.h"
#include "oops/runs/Test.h"
#include "oops/util/parameters/IgnoreOtherParameters.h"

#include "soca/Traits.h"
#include "soca/Utils/OceanSmoother.h"

namespace test {

typedef oops::Geometry<soca::Traits> Geometry_;
typedef oops::Increment<soca::Traits> Increment_;
typedef oops::State<soca::Traits> State_;

// -----------------------------------------------------------------------------

class OceanSmootherTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OceanSmootherTestParameters, Parameters)

 public :
  oops::RequiredParameter<eckit::LocalConfiguration> geometry{"geometry", this};
  oops::RequiredParameter<soca::OceanSmoother::Parameters> smoother{"oceanSmoother", this};
  oops::RequiredParameter<eckit::LocalConfiguration> output{"output", this};
  oops::RequiredParameter<eckit::LocalConfiguration> state{"state", this};

  oops::IgnoreOtherParameters ignoreOthers{this};
};

// -----------------------------------------------------------------------------

void testMultiply() {
  OceanSmootherTestParameters parameters;
  parameters.validateAndDeserialize(TestEnvironment::config());

  // read in geometry
  Geometry_ geom (parameters.geometry, oops::mpi::world(), oops::mpi::myself());

  // read in state
  State_ state(geom, parameters.state);

  // create a random increment
  const std::string varname = "sea_water_potential_temperature";
  oops::Variables vars({oops::Variable(varname)});
  util::DateTime time(20240101,0);
  Increment_ dx(geom, vars, time);
  dx.random();
  oops::Log::test() << "random: " << dx << std::endl;
  atlas::FieldSet fset = dx.fieldSet().fieldSet();

  // initialize the smoother
  const int levels = fset.field(varname).shape(1);
  soca::OceanSmoother smoother(geom.generic(), parameters.smoother, levels, state.fieldSet().fieldSet());

  // smooth it and write it out
  smoother.multiply(fset);
  dx.fromFieldSet(fset);
  oops::Log::test() << "smoothed random: " << dx << std::endl;
  dx.write(parameters.output);

  // test that smoothing a field of ones results in ones
  dx.ones();
  fset = dx.fieldSet().fieldSet();
  smoother.multiply(fset);
  dx.fromFieldSet(fset);
  oops::Log::test() << "smoothed ones: " << dx << std::endl;

  // test that smoothing a field of zeros results in zeros
  dx.zero();
  fset = dx.fieldSet().fieldSet();
  smoother.multiply(fset);
  dx.fromFieldSet(fset);
  oops::Log::test() << "smoothed zeros: " << dx << std::endl;
}

// -----------------------------------------------------------------------------

class OceanSmootherTest : public oops::Test {
 private:
  std::string testid() const override {return "test::OceanSmootherTest";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();
    ts.emplace_back(CASE("soca/OceanSmoother/testMultiply")
                    {testMultiply(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test

// -----------------------------------------------------------------------------

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::OceanSmootherTest tests;
  return run.execute(tests);
}
