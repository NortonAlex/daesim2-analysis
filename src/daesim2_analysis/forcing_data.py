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
    SiteX           : ClimateModule
    sowing_dates    : list[date]
    harvest_dates   : list[date]
    df              : DataFrame
    
    start_doy_f     : int = field(init=False)
    start_year_f    : int = field(init=False)
    nrundays_f      : int = field(init=False)

    time_nday_f     : list[int] = field(init=False)
    time_doy_f      : list[float] = field(init=False)
    time_year_f     : list[int] = field(init=False)
    time_index_f    : list[int] = field(init=False)

    _soilTheta_z    : np.ndarray = field(init=False)

    inputs          : list[np.ndarray] = field(init=False)

    def __post_init__(s: Self):
        object.__setattr__(s, 'start_doy_f', s.df['DOY'].values[0])
        object.__setattr__(s, 'start_year_f', s.df['Year'].values[0])
        object.__setattr__(s, 'nrundays_f', s.df.index.size)

        time_nday_f, time_doy_f, time_year_f = s.SiteX.time_discretisation(s.start_doy_f, s.start_year_f, nrundays=s.nrundays_f)
        time_doy_f = [time_doy_f[i]+0.5 for i in range(len(time_doy_f))]
        object.__setattr__(s, 'time_nday_f', time_nday_f)
        object.__setattr__(s, 'time_doy_f', time_doy_f)
        object.__setattr__(s, 'time_year_f', time_year_f)
        object.__setattr__(s, 'time_index_f', to_datetime(s.df['Date'].values))
        object.__setattr__(s, '_soilTheta_z', calculate_soilTheta_z(s.df))

        Climate_doy_f = interp_forcing(time_nday_f, time_doy_f, kind="pconst") #, fill_value=(time_doy[0],time_doy[-1]))
        Climate_year_f = interp_forcing(time_nday_f, time_year_f, kind="pconst") #, fill_value=(time_year[0],time_year[-1]))
        Climate_airTempCMin_f = interp1d(time_nday_f, s.df["Minimum temperature"].values)
        Climate_airTempCMax_f = interp1d(time_nday_f, s.df["Maximum temperature"].values)
        Climate_airTempC_f = interp1d(time_nday_f, (s.df["Minimum temperature"].values+s.df["Maximum temperature"].values)/2)
        Climate_solRadswskyb_f = interp1d(time_nday_f, 10*(s.df["Global Radiation"].values-s.df["Diffuse Radiation"].values))
        Climate_solRadswskyd_f = interp1d(time_nday_f, 10*s.df["Diffuse Radiation"].values)
        Climate_airPressure_f = interp1d(time_nday_f, 100*s.df["Pressure"].values)
        Climate_airRH_f = interp1d(time_nday_f, s.df["Relative Humidity"].values)
        Climate_airU_f = interp1d(time_nday_f, s.df["Uavg"].values)
        Climate_airCO2_f = interp1d(time_nday_f, s.df["Atmospheric CO2 Concentration (bar)"].values)
        Climate_airO2_f = interp1d(time_nday_f, s.df["Atmospheric O2 Concentration (bar)"].values)
        Climate_soilTheta_z_f = interp1d(time_nday_f, s._soilTheta_z, axis=0)  # Interpolates across timesteps, handles all soil layers at once
        Climate_nday_f = interp1d(time_nday_f, time_nday_f)

        inputs = [ 
            Climate_solRadswskyb_f,
            Climate_solRadswskyd_f,
            Climate_airTempCMin_f,
            Climate_airTempCMax_f,
            Climate_airPressure_f,
            Climate_airRH_f,
            Climate_airCO2_f,
            Climate_airO2_f,
            Climate_airU_f,
            Climate_soilTheta_z_f,
            Climate_doy_f,
            Climate_year_f
        ]

        object.__setattr__(s, 'inputs', inputs)







