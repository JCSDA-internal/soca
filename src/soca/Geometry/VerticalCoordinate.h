/*
 * (C) Copyright 2026- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

namespace eckit {
namespace mpi {
class Comm;
}
}

namespace soca {

/// Vertical coordinate for vertical localization, one value per layer, top down:
///   "levels"     layer index 1..nlevs
///   "m"          layer-centre depth (m)
///   "log depth"  ln(depth + offset)
std::vector<double> verticalCoordinate(const std::string & units, const size_t nlevs,
                                       const std::vector<double> & layerDepths,
                                       const double logDepthOffset);

/// One depth profile for the grid from the layer thicknesses `var` ([time,] layer, y, x) in
/// netCDF `file`. Depth of layer k: median layer-centre depth over columns where layers k and
/// k+1 are wet (all wet columns for the bottom layer). Wet: not missing and > 0.
/// Read on rank 0, broadcast to all ranks.
std::vector<double> layerDepthsFromFile(const std::string & file, const std::string & var,
                                        const size_t nlevs, const eckit::mpi::Comm & comm);

}  // namespace soca
