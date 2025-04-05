from daesim2_analysis import fast_sensitivity as fastsa
from daesim.plant_1000 import PlantModuleCalculator
from daesim.utils import daesim_io_write_diag_to_nc
from SALib.sample import fast_sampler
from SALib.analyze import fast
from numpy import column_stack
from pandas import to_datetime
from pandas import DataFrame
from pandas import read_csv
from pandas import concat
import numpy as np

def load_df_forcing(paths_df_forcing: list[str])->DataFrame:
    dfs: list[DataFrame] = []
    for path_df_forcing in paths_df_forcing:
        df = read_csv(path_df_forcing)
        df['Date'] = to_datetime(df['Date'])
        dfs += [df]
    df_forcing = concat(dfs)
    df_forcing = df_forcing.drop_duplicates(subset='Date')
    df_forcing = df_forcing.sort_values(by='Date')
    df_forcing = df_forcing.reset_index(drop=True)
    df_forcing['DOY'] = df_forcing['Date'].dt.dayofyear
    df_forcing['Year'] = df_forcing['Date'].dt.year
    return df_forcing

def calculate_soilTheta_z(df: DataFrame):
    moistures = df[[c for c in df.columns if c.lower().startswith('soil moisture')]]
    # moistures = moistures.dropna(axis=1, how='all')
    return column_stack(moistures.values).T

def get_iparamsets_per_run(param_values: np.ndarray, n_processes: int):
    return [
        list(range(i, i+n_processes))
        for i
        in range(0, param_values.shape[0], n_processes)
    ]

def evaluate_paramset(
    iparamset: int,
    param_values: np.ndarray,
    PlantX: PlantModuleCalculator,
    input_data: list[np.ndarray],
    parameters_df: DataFrame,
    problem: dict,
    xsite: str,
    dir_xsite_parameters: str,
    time_index: np.ndarray,
    title: str,
    description: str
):
    
    nparamset = iparamset + 1
    paramset = param_values[iparamset]
    model_output = fastsa.update_and_run_model(paramset, PlantX, input_data, parameters_df, problem)
    Mpxi, diagnostics = model_output[0], model_output[1]

    nsigfigures = len(str(np.shape(param_values)[0]))
    filename_write = f"FAST_results_{xsite}_paramset{nparamset:0{nsigfigures}}.nc"
    daesim_io_write_diag_to_nc(
        PlantX,
        diagnostics,
        dir_xsite_parameters,
        filename_write,
        time_index,
        problem=problem,
        param_values=paramset,
        nc_attributes={'title': title, 'description': description}
    )
    return np.insert(Mpxi, 0, nparamset)