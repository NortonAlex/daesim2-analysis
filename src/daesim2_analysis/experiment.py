import sys
import attr
from argparse import ArgumentParser
from typing_extensions import Self
from pandas import Timestamp
from datetime import date
from os import makedirs
from os.path import join

from daesim.management import ManagementModule
from daesim.plantgrowthphases import PlantGrowthPhases
from daesim.boundarylayer import BoundaryLayerModule
from daesim.leafgasexchange2 import LeafGasExchangeModule2
from daesim.canopygasexchange import CanopyGasExchange
from daesim.plantcarbonwater import PlantModel as PlantCH2O
from daesim.plantallocoptimal import PlantOptimalAllocation
from daesim.canopylayers import CanopyLayers
from daesim.canopyradiation import CanopyRadiation
from daesim.soillayers import SoilLayers
from daesim.plant_1000 import PlantModuleCalculator
from daesim.utils import ODEModelSolver
from daesim.climate import *
from functools import partial
from daesim2_analysis.forcing_data import ForcingData
from daesim2_analysis.daesim_module_partial_load import DAESIMModulePartialLoad

def is_interactive() -> bool: return hasattr(sys, 'ps1') or sys.flags.interactive

PlantDevXPartialLoad = DAESIMModulePartialLoad(
    PlantGrowthPhases,
    phases=["germination","vegetative","spike","anthesis","grainfill","maturity"],
    gdd_requirements=[50,800,280,150,300,300],
    vd_requirements=[0,30,0,0,0,0],
    allocation_coeffs=[
        [0.2,0.1,0.7,0.0,0.0],
        [0.5,0.1,0.4,0.0,0.0],
        [0.3,0.4,0.3,0.0,0.0],
        [0.3,0.4,0.3,0.0,0.0],
        [0.1,0.1,0.1,0.7,0.0],
        [0.1,0.1,0.1,0.7,0.0]
    ],
    turnover_rates=[
        [0.001,0.001,0.001,0.0,0.0],
        [0.01,0.002,0.008,0.0,0.0],
        [0.01,0.002,0.008,0.0,0.0],
        [0.01,0.002,0.008,0.0,0.0],
        [0.033,0.016,0.033,0.0002,0.0],
        [0.10,0.033,0.10,0.0002,0.0]
    ]
)

@attr.define(frozen=True)
class Experiment:
    CLatDeg                 : float = -36.05
    CLonDeg                 : float = 146.5
    tz                      : int = 10
    crop_type               : str = "Wheat"
    sowing_dates            : list[date] = attr.Factory(lambda: [date(1971,5,11)])
    harvest_dates           : list[date] = attr.Factory(lambda: [date(1971,12,23)])
    n_processes             : int = 1
    n_samples               : int = 100
    dir_results             : str = "DAESIM_data/FAST_results"
    paths_df_forcing        : list[str] = attr.Factory(lambda: ["DAESIM_data/DAESim_forcing_data/Rutherglen_1971.csv"])
    path_parameters_file    : str = "parameters/Fast1.json"

    ClimateModule           : DAESIMModulePartialLoad = DAESIMModulePartialLoad(ClimateModule)
    ManagementModule        : DAESIMModulePartialLoad = DAESIMModulePartialLoad(ManagementModule)
    ForcingData             : DAESIMModulePartialLoad = DAESIMModulePartialLoad(ForcingData)
    PlantGrowthPhases       : DAESIMModulePartialLoad = PlantDevXPartialLoad
    CanopyLayers            : DAESIMModulePartialLoad = DAESIMModulePartialLoad(CanopyLayers, nlevmlcan=3)
    BoundayLayerModule      : DAESIMModulePartialLoad = DAESIMModulePartialLoad(BoundaryLayerModule)
    LeafExchangeModule2     : DAESIMModulePartialLoad = DAESIMModulePartialLoad(LeafGasExchangeModule2)
    CanopyRadiation         : DAESIMModulePartialLoad = DAESIMModulePartialLoad(CanopyRadiation)
    CanopyGasExchange       : DAESIMModulePartialLoad = DAESIMModulePartialLoad(CanopyGasExchange)
    PlantCH2O               : DAESIMModulePartialLoad = DAESIMModulePartialLoad(PlantCH2O)
    PlantOptimalAllocation  : DAESIMModulePartialLoad = DAESIMModulePartialLoad(PlantOptimalAllocation)
    SoilLayers              : DAESIMModulePartialLoad = DAESIMModulePartialLoad(SoilLayers)
    PlantModuleCalculator   : DAESIMModulePartialLoad = DAESIMModulePartialLoad(PlantModuleCalculator)

    xsite                   : str = attr.field(init=False)
    title                   : str = attr.field(init=False)
    description             : str = attr.field(init=False)
    dir_xsite_FAST_results  : str = attr.field(init=False)
    dir_xsite_parameters    : str = attr.field(init=False)
    path_Mpx                : str = attr.field(init=False)

    def __attrs_post_init__(self):
        xsite = '-'.join(p.split('/')[-1].split('.')[0] for p in self.paths_df_forcing)
        object.__setattr__(self, 'xsite', xsite)
        title = f'DAESIM2-Plant FAST Sensitivity Analysis {xsite}'
        object.__setattr__(self, 'title', title)
        object.__setattr__(self, 'description', title)
        dir_fast = join(self.dir_results, xsite)
        params_dir = join(dir_fast, 'parameters')
        object.__setattr__(self, 'dir_xsite_FAST_results', dir_fast)
        object.__setattr__(self, 'dir_xsite_parameters', params_dir)
        object.__setattr__(self, 'path_Mpx', join(dir_fast, 'Mpx.npy'))
        makedirs(dir_fast, exist_ok=True)
        makedirs(params_dir, exist_ok=True)
        object.__setattr__(self, 'sowing_dates', [Timestamp(d) for d in self.sowing_dates])
        object.__setattr__(self, 'harvest_dates', [Timestamp(d) for d in self.harvest_dates])

    @staticmethod
    def from_cli() -> 'Experiment':
        parser = ArgumentParser()
        g1 = parser.add_argument_group('Optimisation Arguments')
        g1.add_argument('--n_processes', type=int, required=True)
        g1.add_argument('--n_samples', type=int, required=True)
        g2 = parser.add_argument_group('File Arguments')
        g2.add_argument('--crop', type=str, required=True)
        g2.add_argument('--dir_results', type=str, required=True)
        g2.add_argument('--paths_df_forcing', type=str, required=True)
        g2.add_argument('--path_parameters_file', type=str, required=True)
        ns = parser.parse_args()
        paths = ns.paths_df_forcing.split(',')

        return Experiment(
            n_processes=ns.n_processes,
            n_samples=ns.n_samples,
            dir_results=ns.dir_results,
            paths_df_forcing=paths,
            path_parameters_file=ns.path_parameters_file
        )
