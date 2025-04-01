from typing_extensions import Union
from typing_extensions import Self
from dataclasses import dataclass
from dataclasses import fields
from dataclasses import field
from hashlib import sha256

@dataclass(frozen=True)
class Parameters:
  paths           : list[str]
  modules         : list[str]
  names           : list[str]
  units           : list[str]
  init            : list[float]
  min             : list[float]
  max             : list[float]
  phase_specific  : list[bool]
  phase           : list[str]

  __iter__ = lambda s: ((f.name, getattr(s, f.name)) for f in fields(s))
  lengths  = property(lambda s: {name:len(value) for name, value in s})
  avg_length = property(lambda s: sum([s.lengths[name] for name in s.lengths]) / len(s.lengths))
  consistent_length = property(lambda s: s.lengths['paths'] - s.avg_length == 0)
  __str__ = lambda s: ' '.join(['-'.join([name, str(value)]) for name, value in s])
  __sha256__ = property(lambda s: sha256(s.__str__().encode()).hexdigest())
  unique_id = property(lambda s: s.__sha256__)

  def __post_init__(s: Self):
    if not s.consistent_length:
      raise ValueError(
        f'All lists must be of the same length, but got lengths: {s.lengths}'
      )

def t():
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
  print(p.unique_id)


if __name__ == '__main__':
  t()