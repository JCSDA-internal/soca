
#include <iostream>

#include "model/Variables/Variables.h"
#include "model/Fortran.h"

// -----------------------------------------------------------------------------
namespace soca {
  // -----------------------------------------------------------------------------

  void Variables::print(std::ostream & os) const {
    int nv;
    //int nl;
    soca_var_info_f90(keyVar_, nv);
    os << nv;
    //if (nl == 1) os << " with LBC";
    //ASSERT(nl == 0 || nl == 1);
  }

  // -----------------------------------------------------------------------------

}  // namespace soca
