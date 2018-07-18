
#include "src/Traits.h"
#include "src/Run/Run.h"
#include "test/interface/Model.h"

int main(int argc,  char ** argv) {
  soca::Run run(argc, argv);
  test::Model<soca::Traits> tests;
  run.execute(tests);
  return 0;
};
