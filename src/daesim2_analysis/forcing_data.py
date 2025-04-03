from daesim2_analysis.utils import calculate_soilTheta_z
from typing_extensions import Self
from dataclasses import dataclass
from pandas import to_datetime
from pandas import DataFrame
from dataclasses import field
from daesim.climate import *
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
        s.Climate_doy_f = interp_forcing(s.time_nday_f, s.time_doy_f, kind="pconst") #, fill_value=(time_doy[0],time_doy[-1]))
        s.Climate_year_f = interp_forcing(s.time_nday_f, s.time_year_f, kind="pconst") #, fill_value=(time_year[0],time_year[-1]))
        s.Climate_airTempCMin_f = interp1d(s.time_nday_f, s.df["Minimum temperature"].values)
        s.Climate_airTempCMax_f = interp1d(s.time_nday_f, s.df["Maximum temperature"].values)
        s.Climate_airTempC_f = interp1d(s.time_nday_f, (s.df["Minimum temperature"].values+s.df["Maximum temperature"].values)/2)
        s.Climate_solRadswskyb_f = interp1d(s.time_nday_f, 10*(s.df["Global Radiation"].values-s.df["Diffuse Radiation"].values))
        s.Climate_solRadswskyd_f = interp1d(s.time_nday_f, 10*s.df["Diffuse Radiation"].values)
        s.Climate_airPressure_f = interp1d(s.time_nday_f, 100*s.df["Pressure"].values)
        s.Climate_airRH_f = interp1d(s.time_nday_f, s.df["Relative Humidity"].values)
        s.Climate_airU_f = interp1d(s.time_nday_f, s.df["Uavg"].values)
        s.Climate_airCO2_f = interp1d(s.time_nday_f, s.df["Atmospheric CO2 Concentration (bar)"].values)
        s.Climate_airO2_f = interp1d(s.time_nday_f, s.df["Atmospheric O2 Concentration (bar)"].values)
        s.Climate_soilTheta_z_f = interp1d(s.time_nday_f, s._soilTheta_z, axis=0)  # Interpolates across timesteps, handles all soil layers at once
        s.Climate_nday_f = interp1d(s.time_nday_f, s.time_nday_f)

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
   


    def __post_init__(s: Self):
        s.set_starts()
        s.time_descretisation()
        s.validate_time()
        s.generate_forcing_inputs()









