#include "oops/runs/Run.h"
#include "soca/Traits.h"
#include "soca/Model/OceanIceEmulator/test/test_OceanIceFFNN.h"


int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  soca::test_OceanIceFFNN testoceaniceemul;
  return run.execute(testoceaniceemul);
}