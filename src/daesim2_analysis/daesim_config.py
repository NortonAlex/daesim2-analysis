from dataclasses import dataclass, field, fields
from pandas import DataFrame
from typing import Any, List, Type
from hashlib import sha256
from json import load

@dataclass(frozen=True)
class DAESIMConfig:
    """
    Container for collecting DAESIM module argument specifications.
    Each field should be a list of the same length, describing one argument across modules.
    """
    module_path    : List[str]   = field(metadata={'label': 'Module Path'})
    module_name    : List[str]   = field(metadata={'label': 'Module'})
    arg_name       : List[str]   = field(metadata={'label': 'Argument'})
    arg_type       : List[Type]  = field(metadata={'label': 'Type'})
    default        : List[Any]   = field(metadata={'label': 'Default'})
    required       : List[bool]  = field(metadata={'label': 'Required'})
    description    : List[str]   = field(metadata={'label': 'Description'})

    # Derived
    df             : DataFrame   = field(init=False)

    def __post_init__(self):
        # Build DataFrame of args
        data = {}
        for f in fields(self):
            if f.init and not f.name.startswith('_'):
                label = f.metadata.get('label', f.name)
                data[label] = getattr(self, f.name)
        object.__setattr__(self, 'df', DataFrame(data))

        # Validate equal lengths
        lengths = {col: len(vals) for col, vals in self.df.items()}
        if len(set(lengths.values())) != 1:
            raise ValueError(f"Inconsistent lengths in DAESIMConfig fields: {lengths}")

    def __iter__(self):
        """Yield tuples of (field_name, label, values) for init fields."""
        for f in fields(self):
            if f.init and not f.name.startswith('_'):
                yield (f.name, f.metadata.get('label', f.name), getattr(self, f.name))

    def make_df(self) -> DataFrame:
        """Return a copy of the arguments DataFrame."""
        return self.df.copy()

    @property
    def lengths(self) -> dict:
        """Length of each argument list by label."""
        return {label: len(col) for label, col in self.df.items()}

    @property
    def avg_length(self) -> float:
        """Average length across all columns."""
        vals = list(self.lengths.values())
        return sum(vals) / len(vals)

    @property
    def consistent_length(self) -> bool:
        """True if all lists share the same length."""
        return len(set(self.lengths.values())) == 1

    @property
    def unique_id(self) -> str:
        """SHA-256 hash of the DataFrame's string representation."""
        s = self.df.to_csv(index=False).encode()
        return sha256(s).hexdigest()

    def __str__(self) -> str:
        return str(self.df)

    def get_module_args(self, module_path: str) -> DataFrame:
        """
        Return a DataFrame of parameters for the given module path.
        """
        # Filter rows where the Module Path matches exactly
        return self.df[self.df['Module Path'] == module_path].copy()

    @classmethod
    def from_dict(cls, config: dict, required: bool = True) -> "DAESIMConfig":
        """
        Construct DAESIMConfig from a nested config dict of the form:
          { "module.path.ClassName": { arg1: val1, arg2: val2, ... }, ... }
        """
        paths = []
        modules = []
        arg_names = []
        arg_types = []
        defaults = []
        required_list = []
        descriptions = []
        for module_path, args_map in config.items():
            module_name = module_path.split('.')[-1]
            for arg, val in args_map.items():
                paths.append(module_path)
                modules.append(module_name)
                arg_names.append(arg)
                arg_types.append(type(val))
                defaults.append(val)
                required_list.append(required)
                descriptions.append(f"{arg} for {module_name}")
        return cls(
            module_path=paths,
            module_name=modules,
            arg_name=arg_names,
            arg_type=arg_types,
            default=defaults,
            required=required_list,
            description=descriptions,
        )

    @classmethod
    def from_json_dict(cls, path: str):
        return cls.from_dict(load(open(path)))


# Example instantiation
if __name__ == "__main__":
    # Define multiple args for ClimateModule and ManagementModule
    args = DAESIMConfig(
        module_path=[
            "daesim.climate.ClimateModule",
            "daesim.climate.ClimateModule",
            "daesim.management.ManagementModule",
        ],
        module_name=["ClimateModule", "ClimateModule", "ManagementModule"],
        arg_name=["CLatDeg", "CLonDeg", "cropType"],
        arg_type=[float, float, str],
        default=[-36.05, 146.5, "Wheat"],
        required=[True, True, True],
        description=[
            "Latitude of the site (deg)",
            "Longitude of the site (deg)",
            "Crop type identifier",
        ],
    )
    print(args.df)
    config = {
        "daesim.climate.ClimateModule": {
            "CLatDeg": -36.05,
            "CLonDeg": 146.5,
        },
        "daesim.management.ManagementModule": {
            "cropType": "wheat"
        }
    }
    args2 = DAESIMConfig.from_dict(config)

    # Inspect the DataFrame
    print(args2.df)

    args3 = DAESIMConfig.from_json_dict('daesim_configs/daesim_config1.json')
    print(args3)