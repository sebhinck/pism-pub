// Copyright (C) 2011, 2012, 2013, 2014, 2015 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <gsl/gsl_math.h>

#include "PAIndexForcing.hh"
#include "base/util/PISMVars.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/pism_const.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/pism_options.hh"
#include "base/util/error_handling.hh"

namespace pism {
namespace atmosphere {

IndexForcing::IndexForcing(IceGrid::ConstPtr g)
  : AtmosphereModel(g),
    m_air_temp_snapshot(m_sys, "air_temp_snapshot") {

  // allocate IceModelVecs for storing temperature and precipitation fields:

  // create mean annual ice equivalent precipitation rate (before separating
  // rain, and before melt, etc. in SurfaceModel)
  m_precipitation.create(m_grid, "precipitation", WITHOUT_GHOSTS);
  m_precipitation.set_attrs("climate_state",
                          "mean annual ice-equivalent precipitation rate",
                          "m s-1",
                          ""); // no CF standard_name ??
  m_precipitation.metadata().set_string("glaciological_units", "m year-1");
  m_precipitation.write_in_glaciological_units = true;
  m_precipitation.set_time_independent(true);

  m_air_temp.create(m_grid, "air_temp", WITHOUT_GHOSTS);
  m_air_temp.set_attrs("climate_state",
                     "mean annual near-surface (2 m) air temperature",
                     "K",
                     "");
  m_air_temp.set_time_independent(true);

  // initialize metadata for "air_temp_snapshot"
  m_air_temp_snapshot.set_string("pism_intent", "diagnostic");
  m_air_temp_snapshot.set_string("long_name",
                               "snapshot of the near-surface air temperature");
  m_air_temp_snapshot.set_string("units", "K");
  
  process_options();
  
}

void IndexForcing::process_options()
{
  
    std::string climate_option_prefix = "IndexForcing_climate";
  
    options::String climate_file(climate_option_prefix + "_file",
                         "Specifies a file with boundary conditions");
    if (climate_file.is_set()) {
      m_climate_file = climate_file;
      m_log->message(2,
                 "  - Reading boundary conditions from '%s'...\n",
                 m_climate_file.c_str());
    } else {
      // find PISM input file to read data from:
      bool do_regrid; int start;   // will be ignored
      find_pism_input(m_climate_file, do_regrid, start);

      m_log->message(2,
                 "  - Option %s_file is not set. Trying the input file '%s'...\n",
                 climate_option_prefix.c_str(), m_climate_file.c_str());
    }

    
    std::string index_option_prefix = "IndexForcing_index";
  
    options::String index_file(index_option_prefix + "_file",
                         "Specifies a file with boundary conditions");
    if (index_file.is_set()) {
      m_index_file = index_file;
      m_log->message(2,
                 "  - Reading boundary conditions from '%s'...\n",
                 m_index_file.c_str());
    } else {
      // find PISM input file to read data from:
      bool do_regrid; int start;   // will be ignored
      find_pism_input(m_index_file, do_regrid, start);

      m_log->message(2,
                 "  - Option %s_file is not set. Trying the input file '%s'...\n",
                 index_option_prefix.c_str(), m_index_file.c_str());
    }
    
    options::Integer period(index_option_prefix + "_period",
                            "Specifies the length of the climate index data period (in years)", 0);
    if (period.value() < 0.0) {
      throw RuntimeError::formatted("invalid %s_period %d (period length cannot be negative)",
                                    option_prefix.c_str(), period.value());
    }
    m_period = (unsigned int)period;

    options::Integer ref_year(index_option_prefix + "_reference_year",
                              "Boundary condition reference year", 0);
    if (ref_year.is_set()) {
      m_reference_time = units::convert(Model::m_sys, ref_year, "years", "seconds");
    } else {
      m_reference_time = 0;
    }
  
}



IndexForcing::~IndexForcing(pism::IceGrid::ConstPtr g) {

  // Do nothing....
  
}

void IndexForcing::init()
{
  init_data();
}

void IndexForcing::init_data()
{
  
  
  
}



}
}

