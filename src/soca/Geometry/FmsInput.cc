/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Geometry/FmsInput.h"

// -----------------------------------------------------------------------------
namespace soca {
  // -----------------------------------------------------------------------------
  FmsInput::FmsInput(const eckit::mpi::Comm & comm,
                     const eckit::Configuration & conf)
    : comm_(comm), conf_(conf) {

    // Get the file name of the fms namelist
    inputnml_orig_ = conf_.getString("mom6_input_nml");

    // "input.nml" is protected
    ASSERT(inputnml_orig_ != "input.nml");

    if ( comm_.rank() == 0 ) {
      // Get serial # of original file
      inputnml_sn_ = getFileSN(inputnml_orig_);
    }
    comm_.broadcast(inputnml_sn_, 0);

    // Check that inputnml_sn_ is a valid serial #
    ASSERT(inputnml_sn_ > 0);
    comm_.barrier();
  }
  // -----------------------------------------------------------------------------
  FmsInput::FmsInput(const FmsInput & other)
    : comm_(other.comm_),
      conf_(other.conf_),
      inputnml_sn_(other.inputnml_sn_) {
  }
  // -----------------------------------------------------------------------------
  FmsInput::~FmsInput() {
    inputnml_sn_ = -1;
  }

  // -----------------------------------------------------------------------------
  void FmsInput::updateNameList() {
    int testlink = -1;
    if ( comm_.rank() == 0 ) {
      // Get serial # of current input.nml
      int inputnml_current_sn = getFileSN("input.nml");

      // If the file's serial # don't match, update the input.nml link
      if ( inputnml_sn_ != inputnml_current_sn ) {
        // Remove input.nml (link or file) if it exist
        struct stat buffer;
        int status = stat("input.nml", &buffer);
        if ( status == 0 ) { remove("input.nml"); }

        // Create new symbolic link
        char* target = const_cast<char*>(inputnml_orig_.c_str());
        char linkname[] = "input.nml";
        testlink = symlink(target, linkname);
        if ( testlink != 0 ) { comm_.abort(); }
      }
    }
    comm_.barrier();
  }
  // -----------------------------------------------------------------------------
  int FmsInput::getFileSN(const std::string & filename) {
    // Get serial # of file
    struct stat buffer;
    char* filename_char = const_cast<char*>(filename.c_str());
    int status = stat(filename_char, &buffer);
    if ( status == 0 ) {
      return buffer.st_ino;
    } else {
      return -1;
    }
  }
}  // namespace soca
