from dataclasses import dataclass
from pandas import DataFrame
from typing_extensions import Union
from typing_extensions import Self
import pandas as pd


def load_df_forcing(paths_df_forcing: list[str])->pd.DataFrame:
    dfs: list[pd.DataFrame] = []
    for path_df_forcing in paths_df_forcing:
        df = pd.read_csv(path_df_forcing)
        df['Date'] = pd.to_datetime(df['Date'])
        dfs += [df]
    df_forcing = pd.concat(dfs)
    df_forcing = df_forcing.drop_duplicates(subset='Date')
    df_forcing = df_forcing.sort_values(by='Date')
    df_forcing = df_forcing.reset_index(drop=True)
    df_forcing['DOY'] = df_forcing['Date'].dt.dayofyear
    df_forcing['Year'] = df_forcing['Date'].dt.year
    return df_forcing

class ForcingData:
    df: DataFrame

    # def __post_init__(s: Self):
    #     object.__setattr__(s, 'df', )
  
    @classmethod
    def __from_csv__(
        cls: 'ForcingData',
        location: Union[str, list[str]],
        diffuse_fraction: float,
        nevlmsoil: Union[int, None],
    )->'ForcingData':
        
        df = load_df_forcing([location]) if isinstance(location, str) else load_df_forcing([location])
        return ForcingData(df)