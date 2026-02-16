# DAESIM2 Analysis

This Python package supports the development and testing of the Dynamic Agro-Ecosystem Simulator 2 (DAESIM2), a biophysical model for agro-ecosystems. The package provides tools to support simulations, sensitivity analysis, calibration, and uncertainty analysis, helping to test, evaluate and improve model performance.

References: Taghikhah et al. (2022) https://doi.org/10.1016/j.ecolmodel.2022.109930

## 1. Overview

### 1.1  What is DAESIM2?

DAESIM2 is a modular, process-based biophysical model designed to simulate the growth, development, and resource dynamics of agro-ecosystems. Currently, the model represents plant physiological processes (e.g. photosynthesis, carbon allocation, phenological development, water uptake) in a mechanistic and extensible framework. 

### 1.2 Design Philosophy

The DAESIM2 analysis package (`daesim2-analysis`) is implemented in Python using an object-oriented architecture that separates:
 * Core biophysical processes
 * Parameter configuration
 * Experiment design
 * Forcing data handling
 * Analysis workflows

This separation ensures:
 * Ease of use — the model can be run immediately using internally defined defaults.
 * Configurability — parameters can be customised for specific sites, genotypes, plant functional types, or experimental designs.
 * Scientific workflow integration — the architecture separates physics from experiment design and supports sensitivity analysis, calibration, optimisation, and ensemble experimentation. 

## 2. Model Analysis Architecture

### 2.1 Conceptual Architecture

The `daesim2-analysis` architecture is structured into four conceptual layers:
 * Biophysical modules (core process representations)
 * Configuration (site / genotype / PFT definition)
 * Parameter overrides (for calibration and experimentation)
 * Experiment (model orchestration + forcing + I/O)

This structure allows a new user to run the model quickly with minimal setup, while also allowing an advanced user to precisely control parameter sets and experiment design. This also enables clean integration with scientific tools such as sensitivity analyses packages (e.g. SALib) or optimization tools (e.g. MCMC frameworks, SciPy's `optimize` package). It also provides a relatively robust framework with which to scale the model e.g. run across multiple sites in parallel. 

Another benefit of this architecture is that it separates core process physics from experiment design. This allows different versions structures of the `DAESIM2` model to be used within the same analysis framework. 

### 2.2 Biophysical Modules

The DAESIM2 biophysical model source code is maintained in a separate repository: [DAESIM2](https://github.com/NortonAlex/DAESIM). 

The `daesim2-analysis` repository provides the experimental and analytical framework that interfaces with the underlying model, but does not define the process-level physics itself. It serves as the primary entry point for configuring, running, analysing, and calibrating the model

Each biophysical process in `DAESIM2` is implemented as a dedicated Python class (module). Within each module:

 * The mathematical formulation of a biological or physical process.
 * Its associated parameters, defined as class attributes with default values.

These class attributes define the model’s default parameter values. This is an important concept as we consider a "model" to be the combination of mathematical formulations and the coefficients used within them. By embedding default parameters within the biophysical modules, `DAESIM2` makes this relationship transparent. 

### 2.3 Parameter Hierarchy

Parameterization follows a hierarchical override system:
```sh
Default class attributes
        ↓
Configuration file
        ↓
Parameters object
```
Later layers override earlier ones.

This structure allows:
 * Quick default runs (no configuration required)
 * Site/genotype-specific setups
 * Clean integration with sensitivity and optimisation frameworks

#### 2.3.1 Default Parameters (Defined in Model Code)

Each `DAESIM2` module is implemented as a Python class. Default parameters are defined directly in the class definition as class attributes. For example:

```sh
from dataclasses import dataclass
from typing import List

@dataclass
class ManagementModule:
    cropType: str = "Wheat"
    sowingDays: List[int] = (120,)
    harvestDays: List[int] = (330,)
```

These class attributes define the model's default parameter values. These defaults ensure the model can run without external configuration and that the code remains self-contained and transparent. Additionally, this means parameters are documented alongside the physics or biology they govern, making it much easier for testing (e.g. unit tests). 

From the `daesim2-analysis` package, users can inspect these defaults by importing the relevant `DAESIM2` module and examining its class attributes.

This is the base parameter layer. 

#### 2.3.2 Configuration File (Site / Genotype / PFT)

A configuration file allows a user to customise model parameters for e.g: 
 * A specific site (e.g. lat/lon, soil properties)
 * A genotype or functional type (e.g. physiological, developmental or structural traits)
 * A scenario or study (e.g. any combination of the above)

The configuration file is typically a JSON file with a nested dictionary structure.

```sh
{
  "management.ManagementModule": {
  	"sowingDays": (130,),
  	"harvestDays": (340,),
  }
  "plantgrowthphases.PlantGrowthPhases": {
    "GDD_emerg": 75.0,
    "GDD_flowering": 820.0,
  },
  "canopy.CanopyLayers": {
    "nlevmlcan": 3,
  },
}
```

**How it works**

 * The top-level key identifies the module using its full import path.
 * The nested dictionary contains:
   * attribute: value pairs
   * These map directly to class attributes.

Internally, the configuration loader instantiates the model modules, updates matching class attributes, and leaves unspecified attributes at their defaults. 

#### 2.3.3 The `Parameters` Object (Sampling / Calibration / Sensitivity)

The `Parameters` object provides a lightweight, flat interface for:

 * Sensitivity analysis (e.g. SALib)
 * Parameter optimisation (e.g. scipy.optimize)
 * MCMC / uncertainty analysis
 * Ensemble experiments

This layer overrides both:

 * Defaults
 * Configuration file values

In addition to parameter values, the `Parameters` object stores additional data and metadata such as parameter initial values, minimum and maximum bounds. This supports direct integration with with optimisation and sampling frameworks that require parameter vectors and constraints. Future extensions may include data such as probability distribution types (e.g. normal, log-normal, uniform) to support probabilistic sampling and Bayesian calibration workflows. 

**Why this layer exists**

Optimisation and sensitivity frameworks typically require:

 * Flat parameter vectors
 * Simple name → value mappings
 * Easy conversion between dict and array

The `Parameters` object provides exactly that. 

A common use case for this could be: A user wants to configure their model to a field trial site and assign some known plant type attributes (e.g. flowing time, crop height). Then, they want calibrate other model parameters against field trial data. In this case, they can create a configuration file specific to that site and then use the `Parameters` object to specify the parameters they want to calibrate.

Example structure:

```sh
parameters = Parameters(
        paths           = ["PlantCH2O.CanopyGasExchange.Leaf", "PlantCH2O", "PlantCH2O", "PlantDev", "PlantDev"],
        modules         = ["Leaf", "PlantCH2O", "PlantCH2O", "PlantDev", "PlantDev"],
        names           = ["Vcmax_opt", "SLA", "maxLAI", "gdd_requirements", "gdd_requirements"],
        units           = ["mol CO2 m-2 s-1", "m2 g d.wt-1", "m2 m-2", "deg C d", "deg C d"],
        init            = [60e-6, 0.03, 6, 900, 650],
        min             = [30e-6, 0.015, 5, 600, 350],
        max             = [120e-6, 0.035, 7, 1800, 700],
        phase_specific  = [False, False, False, True, True],
        phase           = [None, None, None, "vegetative", "grainfill"]
    )
```

## 3. Running DAESIM2

### 3.1 Quick Start (Minimal Example)

For a quick start a user must:

 * Instantiate Experiment
 * Provide forcing data
 * Run model

An example is provided in `./notebooks/experiment_example.py`. 

### 3.2 Configuration-Based Run

You can copy and modify the example provided in `./notebooks/experiment_example.py`. You may like to provide your own forcing data and a customised configuration (either directly in-line or with a JSON file). 

### 3.3 Calibration/Sensitivity Analysis

Discuss with developers. 

## 4. The Experiment Class

The `Experiment` class is the high-level orchestrator of `daesim2-analysis`. It brings together:

 * The biophysical model modules
 * Configuration and parameterization
 * Forcing data
 * Model integration and I/O

Conceptually:

```sh
Experiment
 ├── Loads configuration
 ├── Applies parameter overrides
 ├── Loads forcing data
 ├── Instantiates model
 ├── Runs simulation
 └── Handles outputs
```

### Creating an Experiment

An Experiment typically requires:

 * Forcing data (e.g. meteorological drivers)
 * A configuration file (optional but recommended)
 * Optional `Parameters` object (for calibration or experiments)

There are some examples at the bottom of this page. 

### Inputs
| Argument        | Type                                         | Default                                                       | Description                                                             |
|-----------------|----------------------------------------------|---------------------------------------------------------------|-------------------------------------------------------------------------|
| `xsite`         | `str`                                        | `'Milgadara_2018'`                                            | Unique site/run identifier (used in output folders).                   |
| `CLatDeg`       | `float`                                      | `-36.05`                                                      | Site latitude in decimal degrees.                                       |
| `CLonDeg`       | `float`                                      | `146.5`                                                       | Site longitude in decimal degrees.                                      |
| `tz`            | `int`                                        | `10`                                                          | Time zone offset from UTC (e.g., 10 for Australia/Sydney).              |
| `crop_type`     | `str`                                        | `'Wheat'`                                                     | Crop species or variety identifier.                                     |
| `sowing_dates`  | `List[date]`                                 | `[date(2018,1,1)]`                                            | List of sowing dates.                                                   |
| `harvest_dates` | `List[date]`                                 | `[date(2018,12,31)]`                                          | List of harvest dates.                                                  |
| `n_processes`   | `int`                                        | `1`                                                           | Number of parallel worker processes.                                    |
| `n_samples`     | `int`                                        | `100`                                                         | Number of parameter samples for FAST.                                    |
| `dir_results`   | `str`                                        | `'DAESIM_data/FAST_results'`                                  | Base directory for storing results.                                     |
| `df_forcing`    | `List[str]` or `DataFrame`                   | `['.../DAESim_forcing_Milgadara_2018.csv']`                  | Forcing CSV paths or pre-loaded DataFrame.                              |
| `parameters`    | `str` or `Parameters`                        | `'parameters/Fast1.json'`                                     | Path or object for SALib parameter definitions.                         |
| `daesim_config` | `str` or `DAESIMConfig`                      | `'daesim_configs/daesim_config1.json'`                        | Path or object for DAESIM module argument defaults.                     |

### Initialization Sequence
1. **Load forcing data**: If `df_forcing` is a path (or list of paths), it is read into a pandas `DataFrame`.
2. **Setup output structure**: Creates `dir_results/{xsite}/parameters` and defines `path_Mpx`.
3. **Convert dates**: Converts `sowing_dates` and `harvest_dates` into `pandas.Timestamp`.
4. **Load configs**:
   - `Parameters` via `Parameters.__from_file__` if a filepath is given.
   - `DAESIMConfig` via `DAESIMConfig.from_json_dict` if a filepath is given.
5. **Instantiate DAESIM modules** in module dependency order. For example:
   ```sh
   SiteX → ForcingDataX → ManagementX → PlantDevX → BoundaryLayerX →
   LeafX → CanopyX → CanopyRadX → CanopyGasExchangeX → SoilLayersX →
   PlantCH2OX → PlantAllocX → PlantX (PlantModuleCalculator)
   ```

6. Setup solver: Wraps PlantX.calculate in ODEModelSolver, setting Model and assembling input_data.



## 5. Installation

To install and run the `daesim2-analysis` package on your local machine, you must download it, install it and then set up the environment. To download, click on the green "Code" button and follow your download method of choice. Alternatively, you can clone the repository from the command line with the following:

```sh
git clone https://github.com/NortonAlex/daesim2_analysis.git
```

Then, change into the newly created directory with `cd daesim2-analysis`. 

### 5.1 Using Conda (Recommended)

Before proceeding further, ensure you have the necessary dependencies. For Mac and Linux users:

* Mac: Install Homebrew (if not already installed)
* Install Conda (via Miniconda or Anaconda)
* Initialize Conda for shell use

Once they are installed, you can proceed with the installation steps below. 

Make your conda environment called `daesim2-analysis` using the environment.yml file and the command:

```sh
conda env create --name daesim2-analysis --file environment.yml
```

Then, you can activate your new environment with the command:

```sh
conda activate daesim2-analysis
```

Now you have created the conda virtual environment to work in, you must install the package. To install the package use the command:

```sh
pip install -e .
```

Hopefully everything installed successfully. 

You can now run the code or open a jupyter notebook to start testing. To open a jupyter notebook run `jupyter-notebook` to open a Jupyter Notebook server or run `jupyter-lab` to open a JupyterLab server (recommended). From there, go into the notebooks directory and work through the examples.

Once you're finished, you can deactivate the conda environment with the command:

```sh
conda deactivate
```

### 5.2 Installation Without Conda

In some environments, Anaconda (conda) may not be available for managing Python environments. This is often the case on high-performance computing (HPC) systems such as NCI's Gadi, where Python environments must be manually configured using system modules and pip.

The following instructions outline how to install and run the `daesim2-analysis` package using a Python virtual environment (venv) instead of Anaconda. First, you must load the required python module. At the time of writing, the required version is Python 3.9.2. On NCI's Gadi, load the correct module with: 

```sh
module load python3/3.9.2
```

If you are not using NCI, check your system’s Python version:

```sh
python3 --version
```

If your system does not have Python 3.9+, install it using your package manager:

* Ubuntu/Debian: sudo apt install python3 python3-venv python3-pip
* CentOS/RHEL: sudo yum install python3 python3-venv python3-pip
* MacOS: brew install python (if using Homebrew)

Then, to download `daesim2-analysis`, clone the repository to your desired directory from the command line with:

```sh
git clone https://github.com/NortonAlex/daesim2_analysis.git
```

Then, navigate to the newly created directory with `cd daesim2-analysis`. Next, to isolate dependencies, you must create a Python virtual environment (venv). On NCI Gadi, store the venv in /g/data or /scratch to avoid exceeding home directory quotas: 

```sh
python3 -m venv /g/data/$projectid/$userid/venv-daesim2-analysis
```

For general systems, you can create the venv in your project directory:

```sh
python3 -m venv venv-daesim2-analysis
```

Activate the new venv in the current session (note the leading "dot-space" `. ` to source the file into the current shell environment, used as a POSIX-standard command instead of the typical `source` which is specific to selected shells e.g. bash, zsh):

```sh
. /g/data/$projectid/$userid/venv-daesim2-analysis/bin/activate
```

Once activated, it is also recommended to update the pip and setuptools packages:

```sh
pip install -U pip setuptools
```

With the virtual environment activated, install the required dependencies listed in requirements.txt: 

```sh
pip install -r requirements.txt
```

Finally, to install `daesim2-analysis` as an editable package (allowing for easy development and updates), run: 

```sh
pip install -e .
```

Hopefully everything installed successfully. You can verify that everything is working by running:

```sh
python -c "import daesim2_analysis; print('Installation successful')"
```

You can now run the code or open a jupyter notebook to start testing. 

Once you're finished, you can deactivate the venv with the command:

```sh
deactivate
```

### 5.3 Linking to the DAESIM2 GitHub repository

This repository is linked to the DAESIM2 biophysical model source code [DAESIM2](https://github.com/NortonAlex/DAESIM). The analysis tools are kept separate from the biophysical model code for better organization and maintainability. Keeping them independent prevents unnecessary complexity in the model codebase and ensures that the analysis tools can evolve separately without impacting core model functionality.

This separation benefits both model users and developers:

* Users can run and interact with the `DAESIM2` model without needing to understand or install additional analysis tools (including their specialized Python libraries).

* Developers can work on analysis tasks (e.g., calibration, sensitivity analysis, uncertainty analysis) without interfering with the core model implementation, also allowing for a stable reference to `DAESIM2` model versions.


#### 5.3.1 Updating the Environment and Dependencies

If the `DAESIM2` biophysical model (linked via GitHub) has been updated or if you want to reference a difference version of the model, the following steps must be followed to ensure those changes are properly integrated into `daesim2-analysis`:

1. Get [DAESIM2](https://github.com/NortonAlex/DAESIM) commit hash

	Note the GitHub commit hash of the `DAESIM2` version that you want to integrate

2. Update the dependencies in `daesim2-analysis`

	2.1 Update `environment.yml`

	Locate the `- pip:` section in `environment.yml` and update the git URL reference with the new commit hash.

  	2.2 Update `setup.cfg`

	Locate the `install_requires` section within the `setup.cfg` file and update the git URL reference with the new commit hash. 

3. Update the Conda Environment

	Option 1: Update the Existing Environment (Recommended)

	Run the following command:

	```sh
	conda env update --file environment.yml --prune
	```

	This will update the environment while keeping existing dependencies intact.

	Option 2: Recreate the Environment (If Issues Occur)

	Sometimes, conda env update does not reflect the new changes and the environment must be removed and created again:

	```sh
	conda env remove -n daesim2-analysis
	conda env create -f environment.yml
	conda activate daesim2-analysis
	```

4. Reinstall the Package

	After updating dependencies, reinstall `daesim2-analysis` to ensure the latest changes are recognized:

	```sh
	pip install -e .
	```

	To force the latest GitHub dependency installation:

	```sh
	pip install --force-reinstall git+https://github.com/NortonAlex/DAESIM@<new-commit-hash>#egg=daesim
	```

5. Verify the Update

	To ensure everything is correctly installed:

	Open a Python shell or Jupyter Notebook.

	Try importing the packages:

	```sh
	import daesim
	import daesim2_analysis
	```

	If there are issues, verify the installed package versions:

	```sh
	pip list | grep daesim
	```

	If any errors occur, consider clearing the pip cache and reinstalling:

	```sh
	pip cache purge
	pip install -e .
	```

## General guidance on working in this repository

### Repository guiding principles

Text todo: Outline how to use this repository. e.g. /data/ should only contain a test set of forcing data, do not commit large datasets to it. Similarly, do not commit anything to /results/. Keep the repo clean and tidy. 

### Branch etiquette

In general, don't push to other people's branches (things can get weird if two people work on the same branch at once).
Instead, check out the branch and then create your own branch in case you want to make any commits (see "Checking out a branch related to a merge request" below).
Then, to add your changes, make a merge request from your own branch back into the original branch.

## Support

If you need support, the first place to go is the [issue tracker] (todo: link the repository / issues).
From there, you can tag other model users to ask for help.
As a second step, reach out directly to your collaborators.

## Other helpful background

### Why do we use Anaconda?

There are benefits to using Anaconda rather than just using pip. As discussed in this article (https://pythonspeed.com/articles/conda-vs-pip/) the main benefits are portability and reproducibility.
- Portability across operating systems: Instead of installing Python in three different ways on Linux, macOS, and Windows, you can use the same environment.yml on all three.
- Reproducibility: It’s possible to pin almost the whole stack, from the Python interpreter upwards.
- Consistent configuration: You don’t need to install system packages and Python packages in two different ways; (almost) everything can go in one file, the environment.yml.

Using Conda also addresses another problem: How to deal with Python libraries that require compiled code? At some point we may like to convert some of the DAESIM2 source code into another, more computationally efficient language (e.g. C++ or Fortran) so that we can run larger simulations. If we decide to do that we will need to compile the source code in a language other than Python, in which case pip would be insufficient and Conda will be required. 

### Why do I only see `.py` files in the `notebooks` directory?

Tracking notebooks is a pain because of all the extra metadata that is in them.
As a result, we use [jupytext](https://jupytext.readthedocs.io/).
This allows us to track clean `.py` files while keeping all the power of jupyter notebooks.
You'll notice that we ignore `*.ipynb` in our `.gitignore`.

Thanks to jupytext, you can just click on the `.py` files and they will open like standard notebooks.

### Checking out a branch related to a merge request

```sh
# Fetch any changes on the remote server
git fetch origin
# If you get no output, you're up to date
# If you get output, it will show you which branches have changed
# on the remote server

# Checkout the others' branch and create a local branch to track it
# git checkout -b branch-name origin/branch-name
# for example:
git checkout -b local-source origin/local-source

# Activate your environment
conda activate daesim2-analysis

# If any DAESIM2 dependencies have changed, update your conda environment
conda env update --file environment.yml --prune

# Checkout (create) your own branch in case you want to make any commits
# (typically we use the branch name plus our initials)
git checkout -b test-notebook-an

# You're now ready to work
# E.g. by starting a notebook server and looking at the notebooks
jupyter-lab
```

## Reading and viewing

Not all of these will be neccesary for everyone, but a helpful list

### What is self

https://www.educative.io/answers/what-is-self-in-python

### General approach to coding

- [Clean Code](https://thixalongmy.haugiang.gov.vn/media/1175/clean_code.pdf)
- [Refactoring](http://silab.fon.bg.ac.rs/wp-content/uploads/2016/10/Refactoring-Improving-the-Design-of-Existing-Code-Addison-Wesley-Professional-1999.pdf)
- [Refactoring guru](refactoring.guru) resource for understanding coding patterns and how to make them better

### Scientific software

- [Research software engineering with Python](https://merely-useful.tech/py-rse/) A good resource for general practice and examples of developing software with Python

### Numerical coding

- [Numerical recipes](http://numerical.recipes/book/book.html)

### Miscellaneous

- [Basic introduction to Jupyter notebooks](https://realpython.com/jupyter-notebook-introduction/)

## Recommended tools

- An IDE (e.g. pycharm) or more light-weight editor (e.g. Sublime)
- If on Mac, homebrew
- Git
- Anaconda (or Miniconda)

## Experiment Examples

### 1. Using JSON configs + CLI
```bash
python run_FAST.py \
  --n_processes 4 \
  --n_samples 200 \
  --crop_type Wheat \
  --CLatDeg -36.05 \
  --CLonDeg 146.5 \
  --dir_results ./results \
  --paths_df_forcing data/forcings1.csv,data/forcings2.csv \
  --parameters parameters/Fast1.json \
  --daesim_config daesim_configs/DAESIM1.json

----------------------------------------------------------------------------------
#DAESIMConfig.json
{
    "plantgrowthphases.PlantGrowthPhases": {
        "phases": [
					"germination",
					"vegetative",
					"spike",
					"anthesis",
					"grainfill",
					"maturity"
				],
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

#Parameters.json
{
  "Module Path": ["PlantCH2O.CanopyGasExchange.Leaf", "PlantCH2O.CanopyGasExchange.Leaf", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantDev", "PlantDev", "", "", "", ""],
  "Module": ["Leaf", "Leaf", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantCH2O", "PlantDev", "PlantDev", "", "", "", ""],
  "Name": ["Vcmax_opt", "g1", "SLA", "maxLAI", "ksr_coeff", "Psi_f", "sf", "gdd_requirements", "gdd_requirements", "GY_FE", "GY_SDW_50", "CI", "d_r_max"],
  "Unit": ["mol CO2 m-2 s-1", "kPa^0.5", "m2 g d.wt-1", "m2 m-2", "g d.wt-1 m-1", "MPa", "MPa-1", "deg C d", "deg C d", "thsnd grains g d.wt spike-1", "g d.wt m-2", "-", "m"],
  "Initial Value": [60e-6, 3, 0.03, 6, 1000, -3.5, 3.5, 900, 650, 0.1, 100, 0.75, 0.5],
  "Min": [30e-6, 1, 0.015, 5, 300, -8.0, 1.5, 600, 350, 0.08, 80, 0.5, 0.15],
  "Max": [120e-6, 6, 0.035, 7, 5000, -1.0, 7.0, 1800, 700, 0.21, 150, 1.0, 0.66],
  "Phase Specific": [false, false, false, false, false, false, false, true, true, false, false, false, false],
  "Phase": [null, null, null, null, null, null, null, "vegetative", "grainfill", null, null, null, null]
}
```
### 2. Using Interactive Instantiation Python

```py
from datetime import date
from daesim2_analysis.experiment import Experiment

# Default interactive experiment (uses defaults defined in Experiment)
exp = Experiment()
print(exp)

# Fully customized instantiation
exp_custom = Experiment(
    xsite="Milgadara_2018",
    CLatDeg=-36.05,
    CLonDeg=146.5,
    crop_type="Wheat",
    sowing_dates=[date(2018, 1, 1)],
    harvest_dates=[date(2018, 12, 31)],
    n_processes=4,
    n_samples=200,
    dir_results="./results",
    df_forcing=[
        "data/forcings1.csv",
        "data/forcings2.csv"
    ],
    parameters="parameters/Fast1.json",
    daesim_config="daesim_configs/daesim_config1.json"
)

# Access instantiated modules and solver
print("Site lat/lon:", exp_custom.SiteX.CLatDeg, exp_custom.SiteX.CLonDeg)
```

### 3. Programmatic instantiation with pre-loaded objects

```py
import pandas as pd
from daesim2_analysis.experiment import Experiment
from daesim2_analysis.parameters import Parameters
from daesim2_analysis.daesim_config import DAESIMConfig
from daesim2_analysis.utils import load_df_forcing

# 1. Load forcing data into a DataFrame
# df_force = pd.read_csv("data/forcings1.csv", parse_dates=["Date"])
df_forcing = load_df_forcing("DAESIM_data/DAESim_forcing_data/DAESim_forcing_Milgadara_2018.csv")

# 2. Create a Parameters object (e.g. from JSON or dict)

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

daesim_config = DAESIMConfig.from_json_dict("daesim_configs/DAESIM1.json")  # or: DAESIMConfig.from_dict(your_config_dict)

# # 4. Instantiate Experiment with DataFrame, Parameters, and DAESIMConfig
exp_programmatic = Experiment(
	xsite="TestSite",
	crop_type='Wheat',
	CLatDeg=-36.05,
	CLonDeg=148.01,
	df_forcing=df_forcing,
	parameters=parameters,
	daesim_config=daesim_config,
)

# # Verify
print(exp_programmatic.ForcingDataX.df.head())
print(exp_programmatic.parameters.df.head())
print(exp_programmatic.daesim_config.df.head())
```
