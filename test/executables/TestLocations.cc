
#include "model/Traits.h"
#include "model/Run/Run.h"
#include "test/interface/Locations.h"

int main(int argc,  char ** argv) {
  model::Run run(argc, argv);
  test::Locations<soca::Traits> tests;
  run.execute(tests);
  return 0;
};

