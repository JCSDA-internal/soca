
#include <iostream>

#include "model/Variables.h"
#include "model/Fortran.h"

// -----------------------------------------------------------------------------
namespace mom5cice5 {
// -----------------------------------------------------------------------------

void Variables::print(std::ostream & os) const {
  int nv;
  //int nl;
  mom5cice5_var_info_f90(keyVar_, nv);
  os << nv;
  //if (nl == 1) os << " with LBC";
  //ASSERT(nl == 0 || nl == 1);
}

// -----------------------------------------------------------------------------

}  // namespace mom5cice5
