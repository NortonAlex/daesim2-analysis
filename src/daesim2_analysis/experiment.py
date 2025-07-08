import sys
import attr
from argparse import ArgumentParser
from typing_extensions import Self
from pandas import Timestamp
from datetime import date
from os import makedirs
from os.path import join
from typing_extensions import Callable

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
from daesim2_analysis.utils import *
from daesim2_analysis.parameters import Parameters

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
    sowing_dates            : list[date] = attr.Factory(lambda: [date(2018,1,1)])
    harvest_dates           : list[date] = attr.Factory(lambda: [date(2018, 12, 31)])
    n_processes             : int = 1
    n_samples               : int = 100
    dir_results             : str = "DAESIM_data/FAST_results"
    paths_df_forcing        : list[str] = attr.Factory(lambda: ["DAESIM_data/DAESim_forcing_data/Rutherglen_1971.csv"])
    path_parameters_file    : str = "parameters/Fast1.json"

    ClimateModule           : DAESIMModulePartialLoad = DAESIMModulePartialLoad(ClimateModule)
    ManagementModule        : DAESIMModulePartialLoad = DAESIMModulePartialLoad(ManagementModule)
    ForcingData             : DAESIMModulePartialLoad = DAESIMModulePartialLoad(ForcingData)
    PlantGrowthPhases       : DAESIMModulePartialLoad = PlantDevXPartialLoad
    CanopyLayers            : DAESIMModulePartialLoad = DAESIMModulePartialLoad(CanopyLayers, nlevmlcan=3, arg=1)
    BoundayLayerModule      : DAESIMModulePartialLoad = DAESIMModulePartialLoad(BoundaryLayerModule)
    LeafExchangeModule2     : DAESIMModulePartialLoad = DAESIMModulePartialLoad(LeafGasExchangeModule2)
    CanopyRadiation         : DAESIMModulePartialLoad = DAESIMModulePartialLoad(CanopyRadiation)
    CanopyGasExchange       : DAESIMModulePartialLoad = DAESIMModulePartialLoad(CanopyGasExchange)
    PlantCH2O               : DAESIMModulePartialLoad = DAESIMModulePartialLoad(PlantCH2O)
    PlantOptimalAllocation  : DAESIMModulePartialLoad = DAESIMModulePartialLoad(PlantOptimalAllocation)
    SoilLayers              : DAESIMModulePartialLoad = DAESIMModulePartialLoad(SoilLayers)
    PlantModuleCalculator   : DAESIMModulePartialLoad = DAESIMModulePartialLoad(PlantModuleCalculator)

    SiteX                   : ClimateModule = attr.field(init=False)
    ForcingDataX            : ForcingData = attr.field(init=False)
    ManagementX             : ManagementModule = attr.field(init=False)
    PlantDevX               : PlantGrowthPhases = attr.field(init=False)
    BoundaryLayerX          : BoundaryLayerModule = attr.field(init=False)
    LeafX                   : LeafGasExchangeModule2 = attr.field(init=False)
    CanopyX                 : CanopyLayers = attr.field(init=False)
    CanopyRadX              : CanopyRadiation = attr.field(init=False)
    CanopyGasExchangeX      : CanopyGasExchange = attr.field(init=False)
    SoilLayersX             : SoilLayers = attr.field(init=False)
    PlantCH2OX              : PlantCH2O = attr.field(init=False)
    PlantAllocX             : PlantOptimalAllocation = attr.field(init=False)
    PlantX                  : PlantModuleCalculator = attr.field(init=False)
    PlantXCalc              : Callable = attr.field(init=False)
    Model                   : ODEModelSolver = attr.field(init=False)
    input_data              : list = attr.field(init=False)
    parameters              : Parameters = attr.field(init=False)

    xsite                   : str = attr.field(init=False)
    title                   : str = attr.field(init=False)
    description             : str = attr.field(init=False)
    dir_xsite_FAST_results  : str = attr.field(init=False)
    dir_xsite_parameters    : str = attr.field(init=False)
    path_Mpx                : str = attr.field(init=False)



    def _setup_output_structure(s: Self):
        xsite = '-'.join(p.split('/')[-1].split('.')[0] for p in s.paths_df_forcing)
        object.__setattr__(s, 'xsite', xsite)
        title = f'DAESIM2-Plant FAST Sensitivity Analysis {xsite}'
        object.__setattr__(s, 'title', title)
        object.__setattr__(s, 'description', title)
        dir_fast = join(s.dir_results, xsite)
        params_dir = join(dir_fast, 'parameters')
        object.__setattr__(s, 'dir_xsite_FAST_results', dir_fast)
        object.__setattr__(s, 'dir_xsite_parameters', params_dir)
        object.__setattr__(s, 'path_Mpx', join(dir_fast, 'Mpx.npy'))
        makedirs(dir_fast, exist_ok=True)
        makedirs(params_dir, exist_ok=True)

    def _dates_to_timestamp(s: Self):
        object.__setattr__(s, 'sowing_dates', [Timestamp(d) for d in s.sowing_dates])
        object.__setattr__(s, 'harvest_dates', [Timestamp(d) for d in s.harvest_dates])

    def _initialise_daesim_modules(s: Self):
        SiteX = s.ClimateModule(CLatDeg=s.CLatDeg,CLonDeg=s.CLonDeg,timezone=s.tz)
        ForcingDataX = s.ForcingData(
            SiteX=SiteX,
            sowing_dates=s.sowing_dates,
            harvest_dates=s.harvest_dates,
            df=load_df_forcing(s.paths_df_forcing)
        )
        ManagementX = s.ManagementModule(
            cropType=s.crop_type,
            sowingDays=ForcingDataX.sowing_days,
            harvestDays=ForcingDataX.harvest_dates,
            sowingYears=ForcingDataX.sowing_years,
            harvestYears=ForcingDataX.harvest_years,
        )
        PlantDevX = s.PlantGrowthPhases()
        BoundaryLayerX = s.BoundayLayerModule(Site=SiteX)
        LeafX = s.LeafExchangeModule2(Site=SiteX)
        CanopyX = s.CanopyLayers()
        CanopyRadX = s.CanopyRadiation(Canopy=CanopyX)
        CanopyGasExchangeX = s.CanopyGasExchange(Leaf=LeafX, Canopy=CanopyX, CanopyRad=CanopyRadX)
        SoilLayersX = s.SoilLayers(nlevmlsoil=ForcingDataX.nlevmlsoil)
        PlantCH2OX = s.PlantCH2O(
            Site=SiteX,
            SoilLayers=SoilLayersX,
            CanopyGasExchange=CanopyGasExchangeX,
            BoundaryLayer=BoundaryLayerX
        )
        PlantAllocX = s.PlantOptimalAllocation(Plant=PlantCH2OX)
        PlantX = s.PlantModuleCalculator(
            Site=SiteX,
            Management=ManagementX,
            PlantDev=PlantDevX,
            PlantCH2O = PlantCH2OX,
            PlantAlloc=PlantAllocX
        )
        PlantXCalc = PlantX.calculate
        Model = ODEModelSolver(
            calculator=PlantXCalc,
            states_init=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            time_start=ForcingDataX.time_axis[0],
            log_diagnostics=True
        )

        input_data = [
            ODEModelSolver,
            ForcingDataX.time_axis,
            ForcingDataX.time_index,
            ForcingDataX.inputs,
            ForcingDataX.reset_days,
            ForcingDataX.zero_crossing_indices
        ]

        object.__setattr__(s, 'SiteX', SiteX)
        object.__setattr__(s, 'ForcingDataX', ForcingDataX)
        object.__setattr__(s, 'ManagementX', ManagementX)
        object.__setattr__(s, 'PlantDevX', PlantDevX)
        object.__setattr__(s, 'BoundaryLayerX', BoundaryLayerX)
        object.__setattr__(s, 'LeafX', LeafX)
        object.__setattr__(s, 'CanopyX', CanopyX)
        object.__setattr__(s, 'CanopyRadX', CanopyRadX)
        object.__setattr__(s, 'CanopyGasExchangeX', CanopyGasExchangeX)
        object.__setattr__(s, 'SoilLayersX', SoilLayersX)
        object.__setattr__(s, 'PlantCH2OX', PlantCH2OX)
        object.__setattr__(s, 'PlantAllocX', PlantAllocX)
        object.__setattr__(s, 'PlantX', PlantX)
        object.__setattr__(s, 'PlantXCalc', PlantXCalc)
        object.__setattr__(s, 'Model', Model)
        object.__setattr__(s, 'input_data', input_data)
        

    def __attrs_post_init__(s: Self):
        s._setup_output_structure()
        s._dates_to_timestamp()
        s._initialise_daesim_modules()
        object.__setattr__(s, 'parameters', Parameters.__from_file__(s.path_parameters_file))

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
        