
#include "model/Traits.h"
#include "model/Run/Run.h"
#include "test/interface/Increment.h"

int main(int argc,  char ** argv) {
  soca::Run run(argc, argv);
  test::Increment<soca::Traits> tests;
  run.execute(tests);
  return 0;
};

