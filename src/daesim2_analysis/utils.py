from numpy import column_stack
from pandas import to_datetime
from pandas import DataFrame
from pandas import read_csv
from pandas import concat

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
