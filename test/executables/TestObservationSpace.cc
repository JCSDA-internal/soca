
#include "model/Traits.h"
#include "oops/runs/Run.h"
#include "test/interface/ObservationSpace.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::ObservationSpace<mom5cice5::Traits> tests;
  run.execute(tests);
  return 0;
};

