%{
#include "coupler/atmosphere/YearlyCycle.hh"
#include "coupler/atmosphere/YearlyCycle.hh"
#include "coupler/atmosphere/Anomaly.hh"
#include "coupler/atmosphere/ConstantPIK.hh"
#include "coupler/atmosphere/CosineYearlyCycle.hh"
#include "coupler/util/ScalarForcing.hh"
#include "coupler/atmosphere/Delta_P.hh"
#include "coupler/atmosphere/Delta_T.hh"
#include "coupler/atmosphere/Frac_P.hh"
#include "coupler/atmosphere/GivenClimate.hh"
#include "coupler/atmosphere/LapseRates.hh"
#include "coupler/atmosphere/OrographicPrecipitation.hh"
#include "coupler/atmosphere/Paleo_precip.hh"
#include "coupler/atmosphere/SeariseGreenland.hh"
#include "coupler/atmosphere/WeatherStation.hh"
#include "coupler/atmosphere/Uniform.hh"
%}

%shared_ptr(pism::atmosphere::AtmosphereModel)
%include "coupler/AtmosphereModel.hh"

%shared_ptr(pism::Diag<pism::atmosphere::AtmosphereModel>)
%template(_AtmosphereDiag) pism::Diag<pism::atmosphere::AtmosphereModel>;

%shared_ptr(pism::atmosphere::YearlyCycle)
%shared_ptr(pism::Diag<pism::atmosphere::YearlyCycle>)
%shared_ptr(pism::atmosphere::PA_mean_july_temp)
%template(_YearlyCycleDiag) pism::Diag<pism::atmosphere::YearlyCycle>;
%rename(AtmosphereYearlyCycle) pism::atmosphere::YearlyCycle;
%include "coupler/atmosphere/YearlyCycle.hh"
  
%shared_ptr(pism::atmosphere::Anomaly)
%rename(AtmosphereAnomaly) pism::atmosphere::Anomaly;
%include "coupler/atmosphere/Anomaly.hh"

%shared_ptr(pism::atmosphere::CosineYearlyCycle)
%rename(AtmosphereCache) pism::atmosphere::CosineYearlyCycle;
%include "coupler/atmosphere/CosineYearlyCycle.hh"

%shared_ptr(pism::atmosphere::Delta_P)
%rename(AtmosphereDeltaP) pism::atmosphere::Delta_P;
%include "coupler/atmosphere/Delta_P.hh"

%shared_ptr(pism::atmosphere::Delta_T)
%rename(AtmosphereDeltaT) pism::atmosphere::Delta_T;
%include "coupler/atmosphere/Delta_T.hh"

%shared_ptr(pism::atmosphere::Frac_P)
%rename(AtmosphereFraceP) pism::atmosphere::Frac_P;
%include "coupler/atmosphere/Frac_P.hh"

%shared_ptr(pism::atmosphere::Given)
%rename(AtmosphereGiven) pism::atmosphere::Given;
%include "coupler/atmosphere/GivenClimate.hh"

%shared_ptr(pism::atmosphere::LapseRates)
%rename(AtmosphereLapseRates) pism::atmosphere::LapseRates;
%include "coupler/atmosphere/LapseRates.hh"

%shared_ptr(pism::atmosphere::OrographicPrecipitation)
%rename(AtmosphereOrographicPrecipitation) pism::atmosphere::OrographicPrecipitation;
%include "coupler/atmosphere/OrographicPrecipitation.hh"

%shared_ptr(pism::atmosphere::PIK)
%rename(AtmospherePIK) pism::atmosphere::PIK;
%include "coupler/atmosphere/ConstantPIK.hh"

%shared_ptr(pism::atmosphere::PaleoPrecip)
%rename(AtmospherePaleoPrecip) pism::atmosphere::PaleoPrecip;
%include "coupler/atmosphere/Paleo_precip.hh"

%shared_ptr(pism::atmosphere::SeaRISEGreenland)
%rename(AtmosphereSeaRISEGreenland) pism::atmosphere::SeaRISEGreenland;
%include "coupler/atmosphere/SeariseGreenland.hh"

%shared_ptr(pism::atmosphere::WeatherStation)
%rename(AtmosphereWeatherStation) pism::atmosphere::WeatherStation;
%include "coupler/atmosphere/WeatherStation.hh"

%shared_ptr(pism::atmosphere::Uniform)
%rename(AtmosphereUniform) pism::atmosphere::Uniform;
%include "coupler/atmosphere/Uniform.hh"
