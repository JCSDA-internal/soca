
#ifndef SOCA_RUN_H
#define SOCA_RUN_H

#include "oops/runs/Run.h"

namespace soca {

// -----------------------------------------------------------------------------

class Run : public oops::Run {
public:
  Run(int, char **);
  ~Run();
};

// -----------------------------------------------------------------------------

}  // namespace soca

#endif // SOCA_RUN_H
