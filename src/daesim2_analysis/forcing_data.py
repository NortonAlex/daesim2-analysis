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



