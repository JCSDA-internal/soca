
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MOM5CICE5_TEST_OP_OBS

#include <iostream>
#include <iomanip>
#include <cmath>

#include <boost/test/unit_test.hpp>
#include "oops/runs/Run.h"
#include "model/Traits.h"
#include "model/instantiateCovarFactory.h"
#include "model/instantiateObsFactory.h"

#include "test/base/TestSuiteOpObsFixture.h"

namespace mom5cice5 {

/*!
 *  TestEnv initializes the test suite for the MOM5CICE5 model as an OOPS application
 */
class TestEnv: public oops::Run {
  public:
    TestEnv() :
        oops::Run(
            (const int) boost::unit_test::framework::master_test_suite().argc,
            (const char**) boost::unit_test::framework::master_test_suite().argv) {
      setup(config());
    }

    virtual ~TestEnv() {
    }

  private:
    // Initialize datas for the whole test suite
    void setup(const eckit::Configuration & fullConfig) {
      instantiateObsFactory();
      instantiateCovarFactory();
      test::TestSuiteOpObsFixture<mom5cice5::Traits>::getInstance().setup(
          fullConfig);
    }
};
}

// -----------------------------------------------------------------------------
typedef mom5cice5::TestEnv TestEnv_;
BOOST_GLOBAL_FIXTURE(TestEnv_);

/*!
 *  Run the generic tests for the MOM5CICE5 model
 */

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(tl_ad, test::TestSuiteOpObsFixture<mom5cice5::Traits>)
#include "test/base/TestSuiteOpObsTLAD.h"
BOOST_AUTO_TEST_SUITE_END()
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(tl, test::TestSuiteOpObsFixture<mom5cice5::Traits>)
#include "test/base/TestSuiteOpObsTL.h"
BOOST_AUTO_TEST_SUITE_END()
// -----------------------------------------------------------------------------
