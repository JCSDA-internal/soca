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

std::shared_ptr<Geometry_> geom;
std::shared_ptr<soca::OceanSmoother> smoother;

// -----------------------------------------------------------------------------

class OceanSmootherTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OceanSmootherTestParameters, Parameters)

 public :
  oops::RequiredParameter<Geometry_::Parameters_> geometry{"geometry", this};
  oops::RequiredParameter<soca::OceanSmoother::Parameters> smoother{"oceanSmoother", this};
  oops::RequiredParameter<eckit::LocalConfiguration> output{"output", this};
  oops::RequiredParameter<eckit::LocalConfiguration> state{"state", this};

  oops::IgnoreOtherParameters ignoreOthers{this};
};

// -----------------------------------------------------------------------------

// void testConstructor() {
//   OceanSmootherTestParameters parameters;
//   parameters.validateAndDeserialize(TestEnvironment::config());

// }

// -----------------------------------------------------------------------------

void testMultiply() {
  OceanSmootherTestParameters parameters;
  parameters.validateAndDeserialize(TestEnvironment::config());

  // read in geometry
  geom.reset(new Geometry_(parameters.geometry, oops::mpi::world(), oops::mpi::myself()));

  // read in state
  State_ state(*geom, parameters.state);

  // create smoother
  EXPECT(geom.get() != nullptr);
  smoother.reset(new soca::OceanSmoother(geom->generic(), parameters.smoother, 25, state.fieldSet().fieldSet()));

  // create a random increment
  oops::Variables vars({oops::Variable("tocn")});
  util::DateTime time(20240101,0);
  Increment_ dx(*geom, vars, time);
  dx.random();
  oops::Log::test() << "random: " << dx << std::endl;

  // smooth it
  atlas::FieldSet fset = dx.fieldSet().fieldSet();
  smoother->multiply(fset);
  dx.fromFieldSet(fset);
  oops::Log::test() << "smoothed random: " << dx << std::endl;
  dx.write(parameters.output);

  // smoothing a field of ones should result in ones
  dx.ones();
  fset = dx.fieldSet().fieldSet();
  smoother->multiply(fset);
  dx.fromFieldSet(fset);
  oops::Log::test() << "smoothed ones: " << dx << std::endl;

  // smoothing a field of zeros should result in zeros
  dx.zero();
  fset = dx.fieldSet().fieldSet();
  smoother->multiply(fset);
  dx.fromFieldSet(fset);
  oops::Log::test() << "smoothed zeros: " << dx << std::endl;

}

// -----------------------------------------------------------------------------

class OceanSmootherTest : public oops::Test {
 private:
  std::string testid() const override {return "test::OceanSmootherTest";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();
    // ts.emplace_back(CASE("soca/OceanSmoother/testConstructor")
    //                 {testConstructor(); });
    ts.emplace_back(CASE("soca/OceanSmoother/testMultiply")
                    {testMultiply(); });
  }

  void clear() const override {}
};
}  // namespace test

// -----------------------------------------------------------------------------

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::OceanSmootherTest tests;
  return run.execute(tests);
}