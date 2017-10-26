
#include "model/Traits.h"
#include "oops/runs/Run.h"
#include "test/interface/Variables.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::Variables<mom5cice5::Traits> tests;
  run.execute(tests);
  return 0;
};

