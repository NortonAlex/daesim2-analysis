from argparse import ArgumentParser
from typing_extensions import Self
from dataclasses import dataclass
from dataclasses import field
from pandas import Timestamp
from datetime import date
from os.path import join
from os import makedirs
import sys

def is_interactive() -> bool: return hasattr(sys, 'ps1') or sys.flags.interactive

@dataclass(frozen=True)
class Args:
    CLatDeg                 : str = -36.05
    CLonDeg                 : str = 146.5
    tz                      : int = 10
    crop_type               : str = 'Wheat'
    sowing_dates            : list[date] = field(default_factory=lambda: [date(year=1971, month=5, day=11)])
    harvest_dates           : list[date] = field(default_factory=lambda: [date(year=1971, month=12, day=23)])
    n_processes             : int = 1
    n_samples               : int = 100
    dir_results             : str = 'DAESIM_data/FAST_results'
    paths_df_forcing        : list[str] = field(default_factory=lambda: ['DAESIM_data/DAESim_forcing_data/Rutherglen_1971.csv'])
    path_parameters_file    : str = 'parameters/Fast1.json'

    ## ManagementX
    phases                  : list[str] = field(default_factory=lambda: ["germination", "vegetative", "spike", "anthesis", "grainfill", "maturity"])
    gdd_requirements        : list[int] = field(default_factory=lambda: [50,800,280,150,300,300])
    vd_requirements         : list[int] = field(default_factory=lambda: [0, 30, 0, 0, 0, 0])
    allocation_coeffs       : list[list[float]] = field(
                                default_factory=lambda: [
                                    [0.2, 0.1, 0.7, 0.0, 0.0],
                                    [0.5, 0.1, 0.4, 0.0, 0.0],
                                    [0.30, 0.4, 0.30, 0.0, 0.0],
                                    [0.30, 0.4, 0.30, 0.0, 0.0],
                                    [0.1, 0.1, 0.1, 0.7, 0.0],
                                    [0.1, 0.1, 0.1, 0.7, 0.0]
                                ]
                            )
    turnover_rates          : list[list[float]] = field(
                                default_factory=lambda: [
                                    [0.001,  0.001, 0.001, 0.0, 0.0],
                                    [0.01, 0.002, 0.008, 0.0, 0.0],
                                    [0.01, 0.002, 0.008, 0.0, 0.0],
                                    [0.01, 0.002, 0.008, 0.0, 0.0],
                                    [0.033, 0.016, 0.033, 0.0002, 0.0],
                                    [0.10, 0.033, 0.10, 0.0002, 0.0]
                                ]
                            )
    
    ## LeafGasExchangeModule2
    Jmax_opt_rVcmax         : float = 0.89
    Jmax_opt_rVcmax_method  : str = 'log'

    ## SoilLayers Args
    nlevmlsoil              : int = 6
    z_max                   : float = 0.66
    z_top                   : float = 0.10
    discretise_method       : str = 'horizon'
    z_horizon               : list[float] = field(default_factory=lambda: [0.06, 0.06, 0.06, 0.10, 0.10, 0.28])
    Psi_e                   : list[float] = field(default_factory=lambda: [-1.38E-03, -1.38E-03, -1.38E-03, -1.32E-03, -2.58E-03, -0.960E-03])
    b_soil                  : list[float] = field(default_factory=lambda: [4.74, 4.74, 4.74, 6.77, 8.17, 10.73])
    K_sat                   : list[float] = field(default_factory=lambda: [29.7, 29.7, 29.7, 25.2, 13.9, 40.9])
    soilThetaMax            : list[float] = field(default_factory=lambda: [0.12, 0.12, 0.12, 0.20, 0.3, 0.4])

    ## PlantCH20X
    maxLAI                  : float = 6.5
    ksr_coeff               : int = 1500
    SLA                     : float = 0.02
    sf                      : float = 1.0
    Psi_f                   : float = -5.0

    ## PlantXModuleCalculator
    GDD_method              : str = 'nonlinear'
    GDD_Tbase               : float = 0.0
    GDD_Topt                : float = 22.5
    GDD_Tupp                : float = 35.0

    d_r_max                 : float = 2.0
    Vmaxremob               : float = 3.0
    Kmremob                 : float = 0.5
    remob_phase             : list[str] = field(default_factory=lambda:['grainfill','maturity'])
    specified_phase         : str = 'spike'
    grainfill_phase         : str = field(default_factory=lambda:['grainfill','maturity'])

    xsite                   : str = field(init=False)
    title                   : str = field(init=False)
    description             : str = field(init=False)
    dir_xsite_FAST_results  : str = field(init=False)
    dir_xsite_parameters    : str = field(init=False)
    path_Mpx                : str = field(init=False)

    def __post_init__(s: Self):
        xsite = '-'.join([path.split('/')[-1].split('.')[0] for path in s.paths_df_forcing])
        title = f'DAESIM2-Plant FAST Sensitivity Analyssis {xsite}'
        description = title
        dir_xsite_FAST_results = join(s.dir_results, xsite)
        dir_xsite_parameters = join(dir_xsite_FAST_results, 'parameters')
        path_Mpx = join(dir_xsite_FAST_results, 'Mpx.npy')

        object.__setattr__(s, 'xsite', xsite)
        object.__setattr__(s, 'title', title)
        object.__setattr__(s, 'description', description)
        object.__setattr__(s, 'dir_xsite_FAST_results', dir_xsite_FAST_results)
        object.__setattr__(s, 'dir_xsite_parameters', dir_xsite_parameters)
        object.__setattr__(s, 'path_Mpx', path_Mpx)

        makedirs(dir_xsite_FAST_results, exist_ok=True)
        makedirs(dir_xsite_parameters, exist_ok=True)

        object.__setattr__(s, 'sowing_dates', [Timestamp(date) for date in s.sowing_dates])
        object.__setattr__(s, 'harvest_dates', [Timestamp(date) for date in s.harvest_dates])


    @staticmethod
    def from_cli() -> 'Args':
        parser = ArgumentParser()
        group1 = parser.add_argument_group('Optimisation Arguments')
        group1.add_argument(
            '--n_processes',
            type=int,
            required=True,
            help='Number of processes for FAST SA'
        )
        group1.add_argument(
            '--n_samples',
            type=int,
            required=True,
            help='Number of samples to generate for FAST SA'
        )

        group2 = parser.add_argument_group('File Arguments')
        group2.add_argument(
            '--crop',
            type=str,
            required=True,
            help='Name of crop'
        )
        group2.add_argument(
            '--dir_results',
            type=str,
            required=True,
            help='Directory for storing FAST SA results'
        )
        group2.add_argument(
            '--paths_df_forcing',
            type=str,
            required=True,
            help='Comma-separated list of forcing data CSVs'
        )
        group2.add_argument(
            '--path_parameters_file',
            type=str,
            required=True,
            help='Path to parameter file'
        )

        args = parser.parse_args()

        return Args(
            n_processes=args.n_processes,
            n_samples=args.n_samples,
            dir_results=args.dir_results,
            paths_df_forcing=args.paths_df_forcing.split(','),
            path_parameters_file=args.path_parameters_file
        )



