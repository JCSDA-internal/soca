
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SOCA_TEST_OP_OBS

#include <iostream>
#include <iomanip>
#include <cmath>

#include <boost/test/unit_test.hpp>
#include "src/Run/Run.h"
#include "src/Traits.h"
#include "src/instantiateCovarFactory.h"
#include "src/instantiateObsFactory.h"

#include "test/base/TestSuiteOpObsFixture.h"

namespace soca {

/*!
 *  TestEnv initializes the test suite for the SOCA model as an OOPS application
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
      test::TestSuiteOpObsFixture<soca::Traits>::getInstance().setup(
          fullConfig);
    }
};
}

// -----------------------------------------------------------------------------
typedef soca::TestEnv TestEnv_;
BOOST_GLOBAL_FIXTURE(TestEnv_);

/*!
 *  Run the generic tests for the SOCA model
 */

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(tl_ad, test::TestSuiteOpObsFixture<soca::Traits>)
#include "test/base/TestSuiteOpObsTLAD.h"
BOOST_AUTO_TEST_SUITE_END()
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(tl, test::TestSuiteOpObsFixture<soca::Traits>)
#include "test/base/TestSuiteOpObsTL.h"
BOOST_AUTO_TEST_SUITE_END()
// -----------------------------------------------------------------------------
