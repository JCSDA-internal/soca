
#include "model/Traits.h"
#include "oops/runs/Forecast.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::Forecast<soca::Traits> fc;
  run.execute(fc);
  return 0;
};
