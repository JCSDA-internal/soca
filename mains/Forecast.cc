
#include "src/Traits.h"
#include "oops/runs/Forecast.h"
#include "src/Run/Run.h"

int main(int argc,  char ** argv) {
  soca::Run run(argc, argv);
  oops::Forecast<soca::Traits> fc;
  run.execute(fc);
  return 0;
};
