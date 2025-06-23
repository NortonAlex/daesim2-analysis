from daesim2_analysis.utils import calculate_soilTheta_z
from typing_extensions import Self
from dataclasses import dataclass
from pandas import to_datetime
from pandas import DataFrame
from dataclasses import field
from daesim.climate import *
from dataclasses import field
import numpy as np

@dataclass(frozen=True)
class ForcingData:
    SiteX                   : ClimateModule
    sowing_dates            : list[date]
    harvest_dates           : list[date]
    df                      : DataFrame
    
    start_doy_f             : int = field(init=False)
    start_year_f            : int = field(init=False)
    nrundays_f              : int = field(init=False)

    time_nday_f             : list[int] = field(init=False)
    time_doy_f              : list[float] = field(init=False)
    time_year_f             : list[int] = field(init=False)
    time_index_f            : list[int] = field(init=False)

    _soilTheta_z            : np.ndarray = field(init=False)

    Climate_doy_f           : np.ndarray = field(init=False)
    Climate_year_f          : np.ndarray = field(init=False)
    Climate_airTempCMin_f   : np.ndarray = field(init=False)
    Climate_airTempCMax_f   : np.ndarray = field(init=False)
    Climate_airTempC_f      : np.ndarray = field(init=False)
    Climate_solRadswskyb_f  : np.ndarray = field(init=False)
    Climate_solRadswskyd_f  : np.ndarray = field(init=False) 
    Climate_airPressure_f   : np.ndarray = field(init=False)
    Climate_airRH_f         : np.ndarray = field(init=False)
    Climate_airU_f          : np.ndarray = field(init=False)
    Climate_airCO2_f        : np.ndarray = field(init=False)
    Climate_airO2_f         : np.ndarray = field(init=False)
    Climate_soilTheta_z_f   : np.ndarray = field(init=False) 
    Climate_nday_f          : np.ndarray = field(init=False)

    sowing_days             : np.ndarray = field(init=False)
    sowing_years            : np.ndarray = field(init=False)
    # harvest_days            : np.ndarray = field(init=False)
    harvest_years           : np.ndarray = field(init=False)
    reset_days              : list[np.uint64] = field(init=False)
    zero_crossing_indices   : list[int] = field(default_factory=lambda: [4, 5, 6])

    def set_starts(s: Self):
        object.__setattr__(s, 'start_doy_f', s.df['DOY'].values[0])
        object.__setattr__(s, 'start_year_f', s.df['Year'].values[0])
        object.__setattr__(s, 'nrundays_f', s.df.index.size)

    def time_descretisation(s: Self):
        time_nday_f, time_doy_f, time_year_f = s.SiteX.time_discretisation(s.start_doy_f, s.start_year_f, nrundays=s.nrundays_f)
        time_doy_f = [time_doy_f[i]+0.5 for i in range(len(time_doy_f))]
        object.__setattr__(s, 'time_nday_f', time_nday_f)
        object.__setattr__(s, 'time_doy_f', time_doy_f)
        object.__setattr__(s, 'time_year_f', time_year_f)
        object.__setattr__(s, 'time_index_f', to_datetime(s.df['Date'].values))
        object.__setattr__(s, '_soilTheta_z', calculate_soilTheta_z(s.df))

    def validate_time(s: Self):
        try: s.SiteX.validate_event_dates(s.sowing_dates, s.time_index_f, event_name='Sowing')
        except Exception as e: raise ValueError('Sowing date before the earliest date in forcing dataframe')
        try: s.SiteX.validate_event_dates(s.harvest_dates, s.time_index_f, event_name='Harvest')
        except Exception as e: raise ValueError('Sowing date before the earliest date in forcing dataframe')

    def generate_forcing_inputs(s: Self):
        Climate_doy_f = interp_forcing(s.time_nday_f, s.time_doy_f, kind='pconst') #, fill_value=(time_doy[0],time_doy[-1]))
        Climate_year_f = interp_forcing(s.time_nday_f, s.time_year_f, kind='pconst') #, fill_value=(time_year[0],time_year[-1]))
        Climate_airTempCMin_f = interp1d(s.time_nday_f, s.df['Minimum temperature'].values)
        Climate_airTempCMax_f = interp1d(s.time_nday_f, s.df['Maximum temperature'].values)
        Climate_airTempC_f = interp1d(s.time_nday_f, (s.df['Minimum temperature'].values+s.df['Maximum temperature'].values)/2)
        Climate_solRadswskyb_f = interp1d(s.time_nday_f, 10*(s.df['Global Radiation'].values-s.df['Diffuse Radiation'].values))
        Climate_solRadswskyd_f = interp1d(s.time_nday_f, 10*s.df['Diffuse Radiation'].values)
        Climate_airPressure_f = interp1d(s.time_nday_f, 100*s.df['Pressure'].values)
        Climate_airRH_f = interp1d(s.time_nday_f, s.df['Relative Humidity'].values)
        Climate_airU_f = interp1d(s.time_nday_f, s.df['Uavg'].values) 
        Climate_airCO2_f = interp1d(s.time_nday_f, s.df['Atmospheric CO2 Concentration (bar)'].values)
        Climate_airO2_f = interp1d(s.time_nday_f, s.df['Atmospheric O2 Concentration (bar)'].values)
        Climate_soilTheta_z_f = interp1d(s.time_nday_f, s._soilTheta_z, axis=0)  # Interpolates across timesteps, handles all soil layers at once
        Climate_nday_f = interp1d(s.time_nday_f, s.time_nday_f)  # nday represents the ordinal day-of-year plus each simulation day (e.g. a model run starting on Jan 30 and going for 2 years will have nday=30+np.arange(2*365))

        object.__setattr__(s, 'Climate_doy_f', Climate_doy_f)
        object.__setattr__(s, 'Climate_year_f', Climate_year_f)
        object.__setattr__(s, 'Climate_airTempCMin_f', Climate_airTempCMin_f)
        object.__setattr__(s, 'Climate_airTempCMax_f', Climate_airTempCMax_f)
        object.__setattr__(s, 'Climate_airTempC_f', Climate_airTempC_f)
        object.__setattr__(s, 'Climate_solRadswskyb_f', Climate_solRadswskyb_f)
        object.__setattr__(s, 'Climate_solRadswskyd_f', Climate_solRadswskyd_f)
        object.__setattr__(s, 'Climate_airPressure_f', Climate_airPressure_f)
        object.__setattr__(s, 'Climate_airRH_f', Climate_airRH_f)
        object.__setattr__(s, 'Climate_airU_f', Climate_airU_f)
        object.__setattr__(s, 'Climate_airCO2_f', Climate_airCO2_f)
        object.__setattr__(s, 'Climate_airO2_f', Climate_airO2_f)
        object.__setattr__(s, 'Climate_soilTheta_z_f', Climate_soilTheta_z_f)
        object.__setattr__(s, 'Climate_nday_f', Climate_nday_f)

    inputs = property(
        lambda s: [
            s.Climate_solRadswskyb_f,
            s.Climate_solRadswskyd_f,
            s.Climate_airTempCMin_f,
            s.Climate_airTempCMax_f,
            s.Climate_airPressure_f,
            s.Climate_airRH_f,
            s.Climate_airCO2_f,
            s.Climate_airO2_f,
            s.Climate_airU_f,
            s.Climate_soilTheta_z_f,
            s.Climate_doy_f,
            s.Climate_year_f
        ]
    )

    def find_sowing_and_harvest_days(s: Self):
        sowing_steps_f = s.SiteX.find_event_steps(s.sowing_dates, s.time_index_f)
        harvest_steps_f = s.SiteX.find_event_steps(s.harvest_dates, s.time_index_f)

        time_axis = s.time_nday_f[sowing_steps_f[0]:harvest_steps_f[-1]+2]
        time_index = s.time_index_f[sowing_steps_f[0]:harvest_steps_f[-1]+2]

        sowing_steps_itax = s.SiteX.find_event_steps(s.sowing_dates, time_index)
        harvest_steps_itax = s.SiteX.find_event_steps(s.harvest_dates, time_index)
        reset_days_itax = s.SiteX.find_event_steps(s.sowing_dates + s.harvest_dates, time_index)

        object.__setattr__(s, 'time_axis', time_axis)
        object.__setattr__(s, 'time_index', time_index)
        object.__setattr__(s, 'sowing_days', np.floor(s.Climate_doy_f(time_axis[sowing_steps_itax])))
        object.__setattr__(s, 'sowing_years', np.floor(s.Climate_year_f(time_axis[sowing_steps_itax])))
        object.__setattr__(s, 'harvest_dates', np.floor(s.Climate_doy_f(time_axis[harvest_steps_itax])))
        object.__setattr__(s, 'harvest_years', np.floor(s.Climate_year_f(time_axis[harvest_steps_itax])))
        object.__setattr__(s, 'reset_days', list(time_axis[sorted(reset_days_itax)]))


    def __post_init__(s: Self):
        s.set_starts()
        s.time_descretisation()
        s.validate_time()
        s.generate_forcing_inputs()
        s.find_sowing_and_harvest_days()