from SALib.sample import fast_sampler
from typing_extensions import Union
from typing_extensions import Self
from dataclasses import dataclass
from dataclasses import fields
from dataclasses import field
from pandas import DataFrame
from hashlib import sha256
import json

@dataclass(frozen=True)
class Parameters:
    paths             : list[str]   = field(metadata={'label': 'Module Path'})
    modules           : list[str]   = field(metadata={'label': 'Module'})
    names             : list[str]   = field(metadata={'label': 'Name'})
    units             : list[str]   = field(metadata={'label': 'Unit'})
    init              : list[float] = field(metadata={'label': 'Initial Value'})
    min               : list[float] = field(metadata={'label': 'Min'})
    max               : list[float] = field(metadata={'label': 'Max'})
    phase_specific    : list[bool]  = field(metadata={'label': 'Phase Specific'})
    phase             : list[str]   = field(metadata={'label': 'Phase'})

    df                : DataFrame   = field(init=False)
    problem           : dict        = field(init=False)
    
    __iter__          = lambda s: iter(((f.name, f.metadata.get('label', f.name), getattr(s, f.name)) for f in fields(s) if f.init == True))
    make_df           = lambda s: DataFrame({label: value for _, label, value in s})
    make_problem      = lambda s: {
                        'num_vars': len(s.df),
                        'names': s.df['Name'].values,
                        'bounds': [[row['Min'], row['Max']] for _, row in s.df.iterrows()]
                        }
    __str__           = lambda s: str(s.df)

    lengths           = property(lambda s: {label:len(value) for _, label, value in s})
    avg_length        = property(lambda s: sum([s.lengths[name] for name in s.lengths]) / len(s.lengths))
    consistent_length = property(lambda s: s.lengths['Module Path'] - s.avg_length == 0)
    __sha256__        = property(lambda s: sha256(s.__str__().encode()).hexdigest())
    unique_id         = property(lambda s: s.__sha256__)
 
    def __post_init__(s: Self):
        if not s.consistent_length:
            raise ValueError(
                f'All lists must be of the same length, but got lengths: {s.lengths}'
            )

        object.__setattr__(s, 'df', s.make_df())
        object.__setattr__(s, 'problem', s.make_problem())

    def sample(s: Self, n: int, seed: int = 0): return fast_sampler.sample(s.problem, n, seed=seed)

    @classmethod
    def __from_json__(cls: 'Parameters', data: Union[str, dict]) -> 'Parameters':
        data = data if isinstance(data, dict) else json.loads(data)
        # Map: label -> field_name
        label_to_field = {f.metadata.get("label", f.name): f.name for f in fields(cls)}

        # Remap input keys (labels) to field names
        remapped_data = {
            label_to_field[key]: value
            for key, value in data.items()
            if key in label_to_field
        }

        return cls(**remapped_data)
    
    @classmethod
    def __from_file__(cls: 'Parameters', filepath: str) -> 'Parameters':
        with open(filepath, 'r', encoding='utf-8') as f:
            data = json.load(f)
        return cls.__from_json__(data)

def test_init()->bool:
    p = Parameters(
        paths           = ["PlantCH2O.CanopyGasExchange.Leaf", "PlantCH2O.CanopyGasExchange.Leaf", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantDev", "PlantDev", "", "", "", ""],
        modules         = ["Leaf", "Leaf", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantDev", "PlantDev", "", "", "", ""],
        names           = ["Vcmax_opt", "g1", "SLA", "maxLAI", "ksr_coeff", "Psi_f", "sf", "gdd_requirements", "gdd_requirements", "GY_FE", "GY_SDW_50", "CI", "d_r_max"],
        units           = ["mol CO2 m-2 s-1", "kPa^0.5", "m2 g d.wt-1", "m2 m-2", "g d.wt-1 m-1", "MPa", "MPa-1", "deg C d", "deg C d", "thsnd grains g d.wt spike-1", "g d.wt m-2", "-", "m"],
        init            = [60e-6, 3, 0.03, 6, 1000, -3.5, 3.5, 900, 650, 0.1, 100, 0.75, 0.5],
        min             = [30e-6, 1, 0.015, 5, 300, -8.0, 1.5, 600, 350, 0.08, 80, 0.5, 0.15],
        max             = [120e-6, 6, 0.035, 7, 5000, -1.0, 7.0, 1800, 700, 0.21, 150, 1.0, 0.66],
        phase_specific  = [False, False, False, False, False, False, False, True, True, False, False, False, False],
        phase           = [None, None, None, None, None, None, None, "vegetative", "grainfill", None, None, None, None]
    )
    return True

def test_from_json():
    example_data = {
        "Module Path": ["PlantCH2O.CanopyGasExchange.Leaf", "PlantCH2O.CanopyGasExchange.Leaf", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantDev", "PlantDev", "", "", "", ""],
        "Module": ["Leaf", "Leaf", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantDev", "PlantDev", "", "", "", ""],
        "Name": ["Vcmax_opt", "g1", "SLA", "maxLAI", "ksr_coeff", "Psi_f", "sf", "gdd_requirements", "gdd_requirements", "GY_FE", "GY_SDW_50", "CI", "d_r_max"],
        "Unit": ["mol CO2 m-2 s-1", "kPa^0.5", "m2 g d.wt-1", "m2 m-2", "g d.wt-1 m-1", "MPa", "MPa-1", "deg C d", "deg C d", "thsnd grains g d.wt spike-1", "g d.wt m-2", "-", "m"],
        "Initial Value": [60e-6, 3, 0.03, 6, 1000, -3.5, 3.5, 900, 650, 0.1, 100, 0.75, 0.5],
        "Min": [30e-6, 1, 0.015, 5, 300, -8.0, 1.5, 600, 350, 0.08, 80, 0.5, 0.15],
        "Max": [120e-6, 6, 0.035, 7, 5000, -1.0, 7.0, 1800, 700, 0.21, 150, 1.0, 0.66],
        "Phase Specific": [False, False, False, False, False, False, False, True, True, False, False, False, False],
        "Phase": [None, None, None, None, None, None, None, "vegetative", "grainfill", None, None, None, None]
    }
    p = Parameters.__from_json__(example_data)
    return True

def test_from_file():
    p = Parameters.__from_file__('parameters/Fast1.json')
    return True

def test_samples_generation():
    p: Parameters = Parameters.__from_file__('parameters/Fast1.json')
    samples = p.sample(100, 0)
    return samples.shape == (p.avg_length * 100, p.avg_length)

def t():
    return all(
        [
            test_init(),
            test_from_json(),
            test_from_file(),
            test_samples_generation()
        ]
    )

if __name__ == '__main__':
    if t(): print('all tests passed')


