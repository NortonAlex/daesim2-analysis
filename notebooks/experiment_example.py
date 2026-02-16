# %% [markdown]
# ## 0. Load relevant modules

# %%
import pandas as pd
from daesim2_analysis.experiment import Experiment
from daesim2_analysis.parameters import Parameters
from daesim2_analysis.daesim_config import DAESIMConfig
from daesim2_analysis.utils import load_df_forcing

# %% [markdown]
# ## 1. Load forcing data into a Pandas DataFrame

# %%
df_forcing = load_df_forcing("../DAESIM_data/DAESim_forcing_data/DAESim_forcing_Milgadara_2018.csv")

# %% [markdown]
# ## 2. DAESIM2 configuration file (e.g. from JSON or dict)

# %% [markdown]
# ### E.g. load from JSON file (best practice)

# %%
daesim_config = DAESIMConfig.from_json_dict("../daesim_configs/DAESIM1.json")  # or: DAESIMConfig.from_dict(your_config_dict)

# %% [markdown]
# ### E.g. define in-line

# %%
d_config = {
    "plantgrowthphases.PlantGrowthPhases": {
        "phases": ["germination","vegetative","spike","anthesis","grainfill","maturity"],
        "gdd_requirements": [50,800,280,150,300,300],
        "vd_requirements": [0,30,0,0,0,0],
        "allocation_coeffs": [
            [0.2,0.1,0.7,0.0,0.0],
            [0.5,0.1,0.4,0.0,0.0],
            [0.3,0.4,0.3,0.0,0.0],
            [0.3,0.4,0.3,0.0,0.0],
            [0.1,0.1,0.1,0.7,0.0],
            [0.1,0.1,0.1,0.7,0.0]
        ],
        "turnover_rates": [
            [0.001,0.001,0.001,0.0,0.0],
            [0.01,0.002,0.008,0.0,0.0],
            [0.01,0.002,0.008,0.0,0.0],
            [0.01,0.002,0.008,0.0,0.0],
            [0.033,0.016,0.033,0.0002,0.0],
            [0.10,0.033,0.10,0.0002,0.0]
        ]
    },
    "canopylayers.CanopyLayers": {
        "nlevmlcan": 3
    }
}

daesim_config2 = DAESIMConfig.from_dict(d_config)

# %% [markdown]
# ## 3. Create a `Parameters` object

# %% [markdown]
# ### E.g. load from JSON file (best practice)

# %%
parameters = Parameters.__from_file__("../parameters/Fast1.json")

# %% [markdown]
# ### E.g. define in-line

# %%
parameters = Parameters(
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

# %% [markdown]
# ### E.g. ignore `Parameters` altogether (rely on default parameters and the configuration)

# %%
parameters = Parameters.__from_file__("../parameters/paramsNone.json")

# %% [markdown]
# ## 4. Instantiate an `Experiment`

# %%
exp_programmatic = Experiment(
	xsite="TestSite",
	crop_type='Wheat',
	CLatDeg=-36.05,
	CLonDeg=148.01,
	df_forcing=df_forcing,
	parameters=parameters,
	daesim_config=daesim_config,
)

# %%
# # Verify
print(exp_programmatic.ForcingDataX.df.head())
print(exp_programmatic.parameters.df.head())
print(exp_programmatic.daesim_config.df.head())

# %%
