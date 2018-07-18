
#include "src/Traits.h"
#include "src/Run/Run.h"
#include "test/interface/Variables.h"

int main(int argc,  char ** argv) {
  soca::Run run(argc, argv);
  test::Variables<soca::Traits> tests;
  run.execute(tests);
  return 0;
};

