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
#include "base/util/io/PIO.hh"

namespace pism {
namespace atmosphere {

IndexForcing::IndexForcing(IceGrid::ConstPtr g)
  : AtmosphereModel(g),
    m_air_temp_snapshot(m_sys, "air_temp_snapshot") {
      
  m_option_prefix = "-atmosphere_index";

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

  m_air_temp.create(m_grid, "air_temp", WITHOUT_GHOSTS);
  m_air_temp.set_attrs("climate_state",
                     "mean annual near-surface (2 m) air temperature",
                     "K",
                     "");

  // initialize metadata for "air_temp_snapshot"
  m_air_temp_snapshot.set_string("pism_intent", "diagnostic");
  m_air_temp_snapshot.set_string("long_name",
                               "snapshot of the near-surface air temperature");
  m_air_temp_snapshot.set_string("units", "K");
  
  process_options();
  
  m_T_0.set_n_records(12);
  m_T_0.create(g, "airtemp_0");
  m_T_0.set_attrs("climate_forcing", "air temperature at t0", "K", "");

  m_T_1.set_n_records(12);
  m_T_1.create(g, "airtemp_1");
  m_T_1.set_attrs("climate_forcing", "air temperature at t1", "K", "");
  
  m_P_0.set_n_records(12);
  m_P_0.create(g, "precip_0");
  m_P_0.set_attrs("climate_forcing", "precipitation at t0", "m s-1", "");
  
  m_P_1.set_n_records(12);
  m_P_1.create(g, "precip_1");
  m_P_1.set_attrs("climate_forcing", "precipitation at t1", "m s-1", "");
  
  m_h_0.create(m_grid, "usurf_0", WITHOUT_GHOSTS);
  m_h_0.set_attrs("climate_state",
                          "surface elevation at t0",
                          "m",
                          ""); // no CF standard_name ??
  m_h_0.set_time_independent(true);
  
  m_h_1.create(m_grid, "usurf_1", WITHOUT_GHOSTS);
  m_h_1.set_attrs("climate_state",
                          "surface elevation at t1",
                          "m",
                          ""); // no CF standard_name ??
  m_h_1.set_time_independent(true);
  
  m_temp_lapse_rate = 5.; // temp lapse rate [K/km]
  m_precip_lapse_rate = 1.; // precip_lapse rate [m/year/km]
  m_precip_thresh_height = 5. ; // threshold height for precip [km]
}

void IndexForcing::process_options()
{
  
  std::string climate_option_prefix = m_option_prefix + "_climate";

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
  
  std::string index_option_prefix = m_option_prefix + "_index";

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
				  index_option_prefix.c_str(), period.value());
  }
  m_period = (unsigned int)period;

  options::Integer ref_year(index_option_prefix + "_reference_year",
			    "Boundary condition reference year", 0);
  if (ref_year.is_set()) {
    m_reference_time = units::convert(m_sys, ref_year, "years", "seconds");
  } else {
    m_reference_time = 0;
  }
  
  m_index = new Timeseries(*m_grid, "glac_index", m_config->get_string("time_dimension_name"));
  m_index->metadata().set_string("units", "");
  m_index->metadata().set_string("long_name", "glacial index");
  m_index->dimension_metadata().set_string("units", m_grid->ctx()->time()->units_string());
 
}



IndexForcing::~IndexForcing() {

  // Do nothing....
  
}

void IndexForcing::init()
{
  m_temp_lapse_rate = options::Real("-temp_lapse_rate", "Elevation lapse rate for the temperature," "in K per km", m_temp_lapse_rate);
  m_temp_lapse_rate = units::convert(m_sys, m_temp_lapse_rate, "K/km", "K/m");

  m_precip_lapse_rate = options::Real("-precip_lapse_rate", "elevation lapse rate for the surface mass balance," "in m/year per km", m_precip_lapse_rate);
  m_precip_lapse_rate = units::convert(m_sys, m_precip_lapse_rate, "m/year / km", "m/s / m");

  m_precip_thresh_height = options::Real("-precip_thresh_height", "Threshold height for precipitation," "km", m_precip_thresh_height);
  m_precip_thresh_height = units::convert(m_sys, m_precip_thresh_height, "km", "m");

  init_data();
}

void IndexForcing::init_data()
{
  
  m_T_0.init(m_climate_file, 0, 1);
  m_T_1.init(m_climate_file, 0, 1);
  m_P_0.init(m_climate_file, 0, 1);
  m_P_1.init(m_climate_file, 0, 1);
  
  
  m_log->message(2,
	    "  reading %s data from forcing file %s...\n",
	    m_index->short_name.c_str(), m_index_file.c_str());

  PIO nc(m_grid->com, "netcdf3");
  nc.open(m_index_file, PISM_READONLY);
  {
    m_index->read(nc, *m_grid->ctx()->time(), *m_grid->ctx()->log());
  }
  nc.close();
  
  //Maybe already calculate the forcing fields at sea level here Instead of repeating the same calculation over and over again!
  
  m_log->message(2,
	    "  reading surface elevation at t0 & t1 data from forcing file %s...\n",
	     m_index_file.c_str());
  
  bool do_regrid = false;
  int start = 0;  
  
  if (do_regrid) {
    m_h_0.regrid(m_climate_file, CRITICAL);
    m_h_1.regrid(m_climate_file, CRITICAL);
  } else {
    m_h_0.read(m_climate_file, start); // fails if not found!
    m_h_1.read(m_climate_file, start);
  }
  
  
}

void IndexForcing::add_vars_to_output_impl(const std::string& keyword, std::set< std::string >& result)
{
  result.insert("precipitation");
  result.insert("air_temp");

  if (keyword == "big") {
    result.insert("air_temp_snapshot");
  }
}

void IndexForcing::define_variables_impl(const std::set< std::string >& vars, const PIO& nc, IO_Type nctype)
{
  if (set_contains(vars, "air_temp_snapshot")) {
    std::string order = m_grid->ctx()->config()->get_string("output_variable_order");
    io::define_spatial_variable(m_air_temp_snapshot, *m_grid, nc, nctype, order, false);
  }

  if (set_contains(vars, "precipitation")) {
    m_precipitation.define(nc, nctype);
  }

  if (set_contains(vars, "air_temp")) {
    m_air_temp.define(nc, nctype);
  }
}

void IndexForcing::write_variables_impl(const std::set< std::string >& vars, const PIO& nc)
{
  if (set_contains(vars, "air_temp_snapshot")) {
    IceModelVec2S tmp;
    tmp.create(m_grid, "air_temp_snapshot", WITHOUT_GHOSTS);
    tmp.metadata() = m_air_temp_snapshot;

    temp_snapshot(tmp);

    tmp.write(nc);
  }

  if (set_contains(vars, "precipitation")) {
    m_precipitation.write(nc);
  }

  if (set_contains(vars, "air_temp")) {
    m_air_temp.write(nc);
  }
}

void IndexForcing::begin_pointwise_access()
{
  m_T_0.begin_access();
  m_T_1.begin_access();
  m_P_0.begin_access();
  m_P_1.begin_access();
  m_h_0.begin_access();
  m_h_1.begin_access();
  m_surface->begin_access();
  m_precipitation.begin_access();
  m_air_temp.begin_access();
}

void IndexForcing::end_pointwise_access()
{
  m_T_0.end_access();
  m_T_1.end_access();
  m_P_0.end_access();
  m_P_1.end_access();
  m_h_0.end_access();
  m_h_1.end_access();
  m_surface->end_access();
  m_precipitation.end_access();
  m_air_temp.end_access();
}

MaxTimestep IndexForcing::max_timestep_impl(double t)
{
  (void) t;
  return MaxTimestep();
}

void IndexForcing::mean_annual_temp(IceModelVec2S& result)
{
  result.copy_from(m_air_temp);
}

void IndexForcing::temp_snapshot(IceModelVec2S& result)
{
  result.copy_from(m_air_temp);
}

void IndexForcing::mean_precipitation(IceModelVec2S& result)
{
  result.copy_from(m_precipitation);
}

void IndexForcing::update_impl(double my_t, double my_dt)
{

  m_t  = m_grid->ctx()->time()->mod(my_t, 1);
  m_dt = my_dt;

 // m_log->message(5,
	//	" Update at t: '%s', dt:'%s'...\n",
	//	m_grid->ctx()->time()->date(my_t).c_str(), m_grid->ctx()->time()->date(my_dt).c_str());

  m_surface = m_grid->variables().get_2d_scalar("surface_altitude");
  
  m_t_index = m_grid->ctx()->time()->mod(my_t - m_reference_time, m_period);
  
  double index = (*m_index)(m_t_index);
  
  m_T_0.update(m_t, m_dt);
  m_T_0.average(m_t, m_dt);
  m_T_1.update(m_t, m_dt);
  m_T_1.average(m_t, m_dt);
  m_P_0.update(m_t, m_dt);
  m_P_0.average(m_t, m_dt);
  m_P_1.update(m_t, m_dt);
  m_P_1.average(m_t, m_dt);
  
  IceModelVec::AccessList list;
  list.add(m_air_temp);
  list.add(m_precipitation);
  list.add(*m_surface);
  list.add(m_T_0);
  list.add(m_T_1);
  list.add(m_P_0);
  list.add(m_P_1);
  list.add(m_h_0);
  list.add(m_h_1);
  
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    m_air_temp(i, j) = compute_T_ij(m_T_0(i,j), m_T_1(i,j), m_h_0(i,j), m_h_1(i,j), (*m_surface)(i,j), index);
    m_precipitation(i, j) = compute_P_ij(m_P_0(i,j), m_P_1(i,j), m_h_0(i,j), m_h_1(i,j), (*m_surface)(i,j), index);
  }
  
}

void IndexForcing::init_timeseries(const std::vector< double >& ts)
{
  m_surface = m_grid->variables().get_2d_scalar("surface_altitude");
  
  int N = ts.size();
  m_ts_index.resize(N);
  m_ts_times.resize(N);
  
  for(unsigned int k=0; k<N; k++){
    m_ts_index[k] = (*m_index)(m_grid->ctx()->time()->mod(ts[k] - m_reference_time, m_period));
    m_ts_times[k] = ts[k];
  }
  
  m_T_0.init_interpolation(ts);
  m_T_1.init_interpolation(ts);
  m_P_0.init_interpolation(ts);
  m_P_1.init_interpolation(ts);

}

void IndexForcing::temp_time_series(int i, int j, std::vector< double >& result)
{
  std::vector<double> T0(m_ts_times.size()),
		      T1(m_ts_times.size());
		      
  m_T_0.interp(i, j, T0);
  m_T_1.interp(i, j, T1);
		      
  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = compute_T_ij(T0[k], T1[k], m_h_0(i,j), m_h_1(i,j), (*m_surface)(i,j), m_ts_index[k]);
  }
}

void IndexForcing::precip_time_series(int i, int j, std::vector< double >& result)
{
  std::vector<double> P0(m_ts_times.size()),
		      P1(m_ts_times.size());
  
  m_P_0.interp(i, j, P0);
  m_P_1.interp(i, j, P1);
		      
  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = compute_P_ij(P0[k], P1[k], m_h_0(i,j), m_h_1(i,j), (*m_surface)(i,j), m_ts_index[k]);
  }
}


double IndexForcing::compute_T_ij(double T0, double T1, double h0, double h1, double h, double index)
{
  double T0_sl = applyLapseRateT(T0, h0, 0.0), // do at init?!
    T1_sl = applyLapseRateT(T1, h1, 0.0),
    T_sl = (T1_sl - T0_sl) * index + T0_sl,
    T = applyLapseRateT(T_sl, 0.0, h);
   
  return(T);
}

double IndexForcing::compute_P_ij(double P0, double P1, double h0, double h1, double h, double index)
{
  double result = applyLapseRateP((P0 + (P1 - P0) * index), 0.0, h); //h_ref is ignored anyway!
  return(std::max(0.0, result));
}

double IndexForcing::applyLapseRateT(double T, double h_ref, double h)
{
  double result = T - m_temp_lapse_rate * (h - h_ref);
  return(result);
}

double IndexForcing::applyLapseRateP(double P, double h_ref, double h)
{
  (void) h_ref;
  double result = P * exp(m_precip_lapse_rate * std::max(0.0, (h - m_precip_thresh_height)));
  return(result);
}





}
}

