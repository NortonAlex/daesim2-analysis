from daesim2_analysis.utils import calculate_soilTheta_z
from typing_extensions import Self
from dataclasses import dataclass
from pandas import to_datetime
from pandas import DataFrame
from dataclasses import field
from daesim.climate import *
from dataclasses import field
import numpy as np
import re

@dataclass(frozen=True)
class ForcingData:
    SiteX                           : ClimateModule
    sowing_dates                    : list[date]
    harvest_dates                   : list[date]
    df                              : DataFrame
    df_type_1_cols                  : list[str] = field(default_factory=lambda: [
                                            'Global Radiation',
                                            'Diffuse Radiation',
                                            'Atmospheric CO2 Concentration (bar)',
                                            'Relative Humidity',
                                            "Soil moisture 5 cm",
                                            "Soil moisture 8 cm"
                                            "Soil moisture 14 cm"
                                            "Soil moisture 22 cm"
                                            "Soil moisture 34 cm"
                                            "Soil moisture 52 cm"
                                        ]
                                    )
    nlevmlsoil                      : int = 2
    df_type_2_cols                  : list[str] = field(default_factory=lambda: ['SRAD', 'VPeff'])
    df_common_cols                  : list[str] = field(default_factory=lambda: [
                                        'Maximum temperature',
                                        'Minimum temperature',
                                        'Uavg'
                                    ]
                                    )

    df_type                         : str = '3'  #field(init=False)
    diffuse_fraction                : float = 0.2
    uniform_moisture_across_layers  : bool = False
    start_doy_f                     : int = field(init=False)
    start_year_f                    : int = field(init=False)
    end_doy_f                       : int = field(init=False)
    end_year_f                      : int = field(init=False)
    nrundays_f                      : int = field(init=False)
    time_nday_f                     : list[int] = field(init=False)
    time_doy_f                      : list[float] = field(init=False)
    time_year_f                     : list[int] = field(init=False)
    time_index_f                    : list[int] = field(init=False)
    _soilTheta_z                    : np.ndarray = field(init=False)

    Climate_doy_f                   : np.ndarray = field(init=False)
    Climate_year_f                  : np.ndarray = field(init=False)
    Climate_airTempCMin_f           : np.ndarray = field(init=False)
    Climate_airTempCMax_f           : np.ndarray = field(init=False)
    Climate_airTempC_f              : np.ndarray = field(init=False)
    Climate_solRadswskyb_f          : np.ndarray = field(init=False)
    Climate_solRadswskyd_f          : np.ndarray = field(init=False) 
    Climate_airPressure_f           : np.ndarray = field(init=False)
    Climate_airRH_f                 : np.ndarray = field(init=False)
    Climate_airU_f                  : np.ndarray = field(init=False)
    Climate_airCO2_f                : np.ndarray = field(init=False)
    Climate_airO2_f                 : np.ndarray = field(init=False)
    Climate_soilTheta_z_f           : np.ndarray = field(init=False) 
    Climate_nday_f                  : np.ndarray = field(init=False)

    sowing_days                     : np.ndarray = field(init=False)
    sowing_years                    : np.ndarray = field(init=False)
    harvest_days                    : np.ndarray = field(init=False)
    harvest_years                   : np.ndarray = field(init=False)
    reset_days                      : list[np.uint64] = field(init=False)
    zero_crossing_indices           : list[int] = field(default_factory=lambda: [4, 5, 6])

    def set_df_type(s: Self):
        if s.df_type == '0':
            # Use column names to determine df_type
            if all(col in s.df.columns for col in s.df_type_1_cols + s.df_common_cols):
                object.__setattr__(s, 'df_type', '1')
            elif all(col in s.df.columns for col in s.df_type_2_cols + s.df_common_cols):
                object.__setattr__(s, 'df_type', '2')
        elif s.df_type == '3':
            pass
        else:
            raise ValueError('The dataframe is missing columns required to build the Climate Variables \n\
                                Either Build Code to support the new config of columns provided or \n\
                                    provide the the type 1 or type 2 columns'
                            )

    def update_nlevmlsoil(s: Self):
        moisture_cols = [c for c in s.df.columns if c.lower().startswith('soil moisture')]
        object.__setattr__(s, 'nlevmlsoil', max(s.nlevmlsoil, len(moisture_cols)))


    def set_starts(s: Self):
        object.__setattr__(s, 'start_doy_f', s.df['DOY'].values[0])
        object.__setattr__(s, 'start_year_f', s.df['Year'].values[0])
        object.__setattr__(s, 'end_doy_f', s.df['DOY'].values[-1])
        object.__setattr__(s, 'end_year_f', s.df['Year'].values[-1])
        object.__setattr__(s, 'nrundays_f', s.df.index.size)

    def time_descretisation(s: Self):
        time_nday_f, time_doy_f, time_year_f = s.SiteX.time_discretisation(s.start_doy_f, s.start_year_f, end_doy=s.end_doy_f, end_year=s.end_year_f)
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
        if s.df_type == '1': # RUTHERGLEN
            moistures = s.df[[c for c in s.df.columns if c.lower().startswith('soil moisture')]]
            _soilTheta_z = np.column_stack(moistures.values)
        
        ### CALCULATE SYNTHETIC SUBSTITUTES WHEN DATA COLUMNS ABSENT

        if s.df_type == '2': # MILGADARA
            num_moisture_cols = len([c for c in s.df.columns if c.lower().startswith('soil moisture')])
            if num_moisture_cols == 0:
                _soilTheta =  0.35*np.ones(s.nrundays_f)
            else:
                soil_moisture_interp = s.df['Soil moisture'].interpolate('quadratic')
                f_soilTheta_min = 0.25
                f_soilTheta_max = 0.40

                f_soilTheta_min_mm = soil_moisture_interp.min()
                f_soilTheta_max_mm = soil_moisture_interp.max()

                f_soilTheta_norm_mm = (soil_moisture_interp.values - f_soilTheta_min_mm)/(f_soilTheta_max_mm - f_soilTheta_min_mm)
                f_soilTheta_norm = f_soilTheta_min + f_soilTheta_norm_mm * (f_soilTheta_max - f_soilTheta_min)
                _soilTheta = f_soilTheta_norm
            
            if s.uniform_moisture_across_layers:
                _soilTheta_z = np.repeat(_soilTheta[:, np.newaxis], s.nlevmlsoil, axis=1)
            else:
                _soilTheta_z0 = _soilTheta-0.06
                _soilTheta_z1 = _soilTheta+0.02
                _soilTheta_z = np.column_stack((_soilTheta_z0, _soilTheta_z1))

            ## Create synthetic data for other forcing variables
            _Rsb_Wm2 = (1-s.diffuse_fraction) * s.df["SRAD"].values * 1e6 / (60*60*24)
            _Rsd_Wm2 = s.diffuse_fraction * s.df["SRAD"].values * 1e6 / (60*60*24)

            _p = 101325*np.ones(s.nrundays_f)
            _es = s.SiteX.compute_sat_vapor_pressure_daily(s.df["Minimum temperature"].values, s.df["Maximum temperature"].values)
            _RH = s.SiteX.compute_relative_humidity(s.df["VPeff"].values/10,_es/1000)
            _RH[_RH > 100] = 100
            _CO2 = 400*(_p/1e5)*1e-6     ## carbon dioxide partial pressure (bar)
            _O2 = 209000*(_p/1e5)*1e-6   ## oxygen partial pressure (bar)

            Climate_doy_f = interp_forcing(s.time_nday_f, s.time_doy_f, kind="pconst", fill_value=(s.time_doy_f[0],s.time_doy_f[-1]))
            Climate_year_f = interp_forcing(s.time_nday_f, s.time_year_f, kind="pconst", fill_value=(s.time_year_f[0],s.time_year_f[-1]))
            Climate_airTempCMin_f = interp1d(s.time_nday_f, s.df["Minimum temperature"].values)
            Climate_airTempCMax_f = interp1d(s.time_nday_f, s.df["Maximum temperature"].values)
            Climate_airTempC_f = interp1d(s.time_nday_f, (s.df["Minimum temperature"].values+s.df["Maximum temperature"].values)/2)
            Climate_solRadswskyb_f = interp1d(s.time_nday_f, _Rsb_Wm2)
            Climate_solRadswskyd_f = interp1d(s.time_nday_f, _Rsd_Wm2)
            Climate_airPressure_f = interp1d(s.time_nday_f, _p)
            Climate_airRH_f = interp1d(s.time_nday_f, _RH)
            Climate_airU_f = interp1d(s.time_nday_f, s.df["Uavg"].values)
            Climate_airCO2_f = interp1d(s.time_nday_f, _CO2)
            Climate_airO2_f = interp1d(s.time_nday_f, _O2)
            Climate_soilTheta_f = interp1d(s.time_nday_f, _soilTheta)
            Climate_soilTheta_z_f = interp1d(s.time_nday_f, _soilTheta_z, axis=0)  # Interpolates across timesteps, handles all soil layers at once
            Climate_nday_f = interp1d(s.time_nday_f, s.time_nday_f)

        if s.df_type == '3': # ASSUMES CONVENTIONAL COLUMN NAMES ARE AVAILABLE
            ## Read soil moisture columns and sort into ascending order (uppermost layer first)
            soil_moisture_cols = [col for col in s.df.columns if col.startswith("Soil moisture interp")]
            def extract_depth(colname):
                match = re.search(r"Soil moisture (\d+)-(\d+)", colname)
                if match:
                    return int(match.group(1))  # use the start depth for sorting
                else:
                    return float('inf')  # fallback if pattern doesn't match

            soil_moisture_cols_sorted = sorted(soil_moisture_cols, key=extract_depth)
            _soilTheta_z = s.df[soil_moisture_cols_sorted].to_numpy()
            _soilTheta = np.nanmean(_soilTheta_z, axis=1)

            Climate_doy_f = interp_forcing(s.time_nday_f, s.time_doy_f, kind="pconst", fill_value=(s.time_doy_f[0],s.time_doy_f[-1]))
            Climate_year_f = interp_forcing(s.time_nday_f, s.time_year_f, kind="pconst", fill_value=(s.time_year_f[0],s.time_year_f[-1]))
            Climate_airTempCMin_f = interp1d(s.time_nday_f, s.df["Minimum temperature"].values)
            Climate_airTempCMax_f = interp1d(s.time_nday_f, s.df["Maximum temperature"].values)
            Climate_airTempC_f = interp1d(s.time_nday_f, (s.df["Minimum temperature"].values+s.df["Maximum temperature"].values)/2)
            Climate_solRadswskyb_f = interp1d(s.time_nday_f, s.df["Downwelling shortwave beam radiation"].values)
            Climate_solRadswskyd_f = interp1d(s.time_nday_f, s.df["Downwelling shortwave diffuse radiation"].values)
            Climate_airPressure_f = interp1d(s.time_nday_f, s.df["Atmospheric pressure"].values)
            Climate_airRH_f = interp1d(s.time_nday_f,  s.df["Relative humidity"].values)
            Climate_airU_f = interp1d(s.time_nday_f, s.df["Wind speed"].values)
            Climate_airCO2_f = interp1d(s.time_nday_f, s.df["Atmospheric carbon dioxide concentration"].values)
            Climate_airO2_f = interp1d(s.time_nday_f, s.df["Atmospheric oxygen concentration"].values)
            Climate_soilTheta_f = interp1d(s.time_nday_f, _soilTheta)
            Climate_soilTheta_z_f = interp1d(s.time_nday_f, _soilTheta_z, axis=0)  # Interpolates across timesteps, handles all soil layers at once
            Climate_nday_f = interp1d(s.time_nday_f, s.time_nday_f)

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
        object.__setattr__(s, 'harvest_days', np.floor(s.Climate_doy_f(time_axis[harvest_steps_itax])))
        object.__setattr__(s, 'harvest_years', np.floor(s.Climate_year_f(time_axis[harvest_steps_itax])))
        object.__setattr__(s, 'reset_days', list(time_axis[sorted(reset_days_itax)]))


    def __post_init__(s: Self):
        s.set_starts()
        s.time_descretisation()
        s.validate_time()
        s.set_df_type()
        s.update_nlevmlsoil()
        s.generate_forcing_inputs()
        s.find_sowing_and_harvest_days()
        # s.time_descretisation()
        # s.validate_time()
        # s.generate_forcing_inputs()
   
