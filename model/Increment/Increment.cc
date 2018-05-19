/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/Increment/Increment.h"

#include <algorithm>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Variables.h"
#include "oops/generic/UnstructuredGrid.h"
#include "util/DateTime.h"
#include "util/Logger.h"
#include "model/ModelBiasIncrement.h"
#include "model/Covariance/ErrorCovariance.h"
#include "model/Fields/Fields.h"
#include "model/Geometry/Geometry.h"
#include "model/State/State.h"
#include "util/Duration.h"
#include "ufo/GeoVaLs.h"
#include "ioda/Locations.h"

using oops::Log;

namespace soca {

  // -----------------------------------------------------------------------------
  /// Constructor, destructor
  // -----------------------------------------------------------------------------
  Increment::Increment(const Geometry & resol, const oops::Variables & vars,
		       const util::DateTime & vt)
    : fields_(new Fields(resol, vars, vt)), stash_()
  {
    fields_->zero();
    Log::trace() << "Increment constructed." << std::endl;
  }
  // -----------------------------------------------------------------------------
  Increment::Increment(const Geometry & resol, const Increment & other)
    : fields_(new Fields(*other.fields_, resol)), stash_()
  {
    Log::trace() << "Increment constructed from other." << std::endl;
  }
  // -----------------------------------------------------------------------------
  Increment::Increment(const Increment & other, const bool copy)
    : fields_(new Fields(*other.fields_, copy)), stash_()
  {
    Log::trace() << "Increment copy-created." << std::endl;
  }
  // -----------------------------------------------------------------------------
  Increment::Increment(const Increment & other)
    : fields_(new Fields(*other.fields_)), stash_()
  {
    Log::trace() << "Increment copy-created." << std::endl;
  }
  // -----------------------------------------------------------------------------
  Increment::~Increment() {
    Log::trace() << "Increment destructed" << std::endl;
  }
  // -----------------------------------------------------------------------------
  void Increment::activateModel() {
    // Should get variables from model. YT
     const std::vector<std::string> vv{
          "cicen",
	  "hicen",
	  "hsnon",
	  "tsfcn",
	  "qsnon",
	  "sicnk",
	  "qicnk",
	  "socn",
	  "tocn",
	  "ssh",
	  "hocn"	    
	};
    oops::Variables vars(vv);    
    //eckit::LocalConfiguration modelvars;
    //modelvars.set("variables", "tl");
    //oops::Variables vars(modelvars);
    // Should get variables from model. YT
    stash_.reset(new Fields(*fields_, vars));
    swap(fields_, stash_);
    ASSERT(fields_);
    ASSERT(stash_);
    Log::trace() << "TLM ... NO IMPLEMENTED" << std::endl;
  }
  // -----------------------------------------------------------------------------
  void Increment::deactivateModel() {
    swap(fields_, stash_);
    *fields_ = *stash_;
    stash_.reset();
    ASSERT(fields_);
    ASSERT(!stash_);
    Log::trace() << "Increment deactivated for TLM" << std::endl;
  }
  // -----------------------------------------------------------------------------
  /// Basic operators
  // -----------------------------------------------------------------------------
  void Increment::diff(const State & x1, const State & x2) {

    ASSERT(this->validTime() == x1.validTime());
    ASSERT(this->validTime() == x2.validTime());
    Log::debug() << "Increment:diff incr " << *fields_ << std::endl;
    Log::debug() << "Increment:diff x1 " << x1.fields() << std::endl;
    Log::debug() << "Increment:diff x2 " << x2.fields() << std::endl;
    fields_->diff(x1.fields(), x2.fields());
  }
  // -----------------------------------------------------------------------------
  Increment & Increment::operator=(const Increment & rhs) {
    *fields_ = *rhs.fields_;
    return *this;
  }
  // -----------------------------------------------------------------------------
  Increment & Increment::operator+=(const Increment & dx) {
    ASSERT(this->validTime() == dx.validTime());
    *fields_ += *dx.fields_;
    return *this;
  }
  // -----------------------------------------------------------------------------
  Increment & Increment::operator-=(const Increment & dx) {
    ASSERT(this->validTime() == dx.validTime());
    *fields_ -= *dx.fields_;
    return *this;
  }
  // -----------------------------------------------------------------------------
  Increment & Increment::operator*=(const double & zz) {
    *fields_ *= zz;
    return *this;
  }
  // -----------------------------------------------------------------------------
  void Increment::zero() {
    fields_->zero();
  }
  // -----------------------------------------------------------------------------
  void Increment::dirac(const eckit::Configuration & config) {
    fields_->dirac(config);
    Log::trace() << "Increment dirac initialized" << std::endl;    
  }  
  // -----------------------------------------------------------------------------
  void Increment::zero(const util::DateTime & vt) {
    fields_->zero(vt);
  }
  // -----------------------------------------------------------------------------
  void Increment::axpy(const double & zz, const Increment & dx,
                       const bool check) {
    ASSERT(!check || this->validTime() == dx.validTime());
    fields_->axpy(zz, *dx.fields_);
  }
  // -----------------------------------------------------------------------------
  void Increment::accumul(const double & zz, const State & xx) {
    fields_->axpy(zz, xx.fields());
  }
  // -----------------------------------------------------------------------------
  void Increment::schur_product_with(const Increment & dx) {
    fields_->schur_product_with(*dx.fields_);
  }
  // -----------------------------------------------------------------------------
  double Increment::dot_product_with(const Increment & other) const {
    return dot_product(*fields_, *other.fields_);
  }
  // -----------------------------------------------------------------------------
  void Increment::random() {
    fields_->random();
  }
  /// Interpolate to observation location
  // -----------------------------------------------------------------------------
  void Increment::interpolateTL(const ioda::Locations & locs, const oops::Variables & vars, ufo::GeoVaLs & cols) const {
    Log::debug() << "Increment::interpolateTL fields in" << *fields_ << std::endl;
    fields_->interpolateTL(locs, vars, cols);
    Log::debug() << "Increment::interpolateTL " << cols << std::endl;    
  }
  // -----------------------------------------------------------------------------
  void Increment::interpolateAD(const ioda::Locations & locs, const oops::Variables & vars, const ufo::GeoVaLs & cols) {
    Log::debug() << "Increment::interpolateAD gom " << cols << std::endl;
    Log::debug() << "Increment::interpolateAD fields in" << *fields_ << std::endl;    
    fields_->interpolateAD(locs, vars, cols);
  }
  // -----------------------------------------------------------------------------
  /// Define and convert to/from unstructured grid
  // -----------------------------------------------------------------------------
  void Increment::define(oops::UnstructuredGrid & ug) const {
    //fields_->define(ug);
  }  
  void Increment::convert_to(oops::UnstructuredGrid & ug) const {
    Log::debug() << "Increment::convert_to" << std::endl;    
    fields_->convert_to(ug);
  }
  // -----------------------------------------------------------------------------
  void Increment::convert_from(const oops::UnstructuredGrid & ug) {
    fields_->convert_from(ug);
  }
  // -----------------------------------------------------------------------------
  /// I/O and diagnostics
  // -----------------------------------------------------------------------------
  void Increment::read(const eckit::Configuration & files) {
    fields_->read(files);
  }
  // -----------------------------------------------------------------------------
  void Increment::write(const eckit::Configuration & files) const {
    fields_->write(files);
  }
  // -----------------------------------------------------------------------------
  void Increment::print(std::ostream & os) const {
    os << std::endl << "  Valid time: " << validTime();
    os << *fields_;
  }
  // -----------------------------------------------------------------------------

}  // namespace soca
