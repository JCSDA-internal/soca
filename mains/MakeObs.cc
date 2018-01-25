
//#include "model/instantiateObsFactory.h"
#include "model/Traits.h"
#include "oops/runs/MakeObs.h"
#include "model/Run/Run.h"

int main(int argc,  char ** argv) {
  soca::Run run(argc, argv);
  oops::MakeObs<soca::Traits> mkobs;
  run.execute(mkobs);
  return 0;
};
