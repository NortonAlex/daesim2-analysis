"""
Analysis helper functions to support DAESIM2 analysis, sensitivity, calibration, etc.
"""

from typing import Any
from netCDF4 import date2num, Dataset
from datetime import datetime, timedelta
import time
import subprocess
import numpy as np
import pandas as pd

def run_model_and_get_outputs(Plant, ODEModelSolver, time_axis, forcing_inputs, reset_days, zero_crossing_indices):
    ## Define the callable calculator that defines the right-hand-side ODE function
    PlantCalc = Plant.calculate
    
    Model = ODEModelSolver(calculator=PlantCalc, states_init=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], time_start=time_axis[0], log_diagnostics=True)

    Model.reset_diagnostics()
    
    ## Run the model solver
    res = Model.run(
        time_axis=time_axis,
        forcing_inputs=forcing_inputs,
        solver="euler",
        zero_crossing_indices=zero_crossing_indices,
        reset_days=reset_days,
    )

    # Convert the defaultdict to a regular dictionary
    _diagnostics = dict(Model.diagnostics)
    # Convert each list in the dictionary to a NumPy array
    diagnostics = {key: np.array(value) for key, value in _diagnostics.items()}

    # Convert the array to a numeric type, handling mixed int and float types
    diagnostics['idevphase_numeric'] = np.array(diagnostics['idevphase'],dtype=np.float64)
    
    # In the model idevphase can equal None but that is not useable in post-processing, so we set None values to np.nan
    diagnostics["idevphase_numeric"][diagnostics["idevphase"] == None] = np.nan

    ## Conversion notes: When _E units are mol m-2 s-1, multiply by molar mass H2O to get g m-2 s-1, divide by 1000 to get kg m-2 s-1, multiply by 60*60*24 to get kg m-2 d-1, and 1 kg m-2 d-1 = 1 mm d-1. 
    ## Noting that 1 kg of water is equivalent to 1 liter (L) of water (because the density of water is 1000 kg/m³), and 1 liter of water spread over 1 square meter results in a depth of 1 mm
    diagnostics["E_mmd"] = diagnostics["E"]*18.015/1000*(60*60*24)
    
    # Turnover rates per pool
    _tr_Leaf = np.zeros(diagnostics['t'].size)
    _tr_Root = np.zeros(diagnostics['t'].size)
    _tr_Stem = np.zeros(diagnostics['t'].size)
    _tr_Seed = np.zeros(diagnostics['t'].size)
    for it, t in enumerate(diagnostics['t']):
        if np.isnan(diagnostics['idevphase_numeric'][it]):
            tr_ = Plant.PlantDev.turnover_rates[-1]
            _tr_Leaf[it] = tr_[Plant.PlantDev.ileaf]
            _tr_Root[it] = tr_[Plant.PlantDev.iroot]
            _tr_Stem[it] = tr_[Plant.PlantDev.istem]
            _tr_Seed[it] = tr_[Plant.PlantDev.iseed]
        else:
            tr_ = Plant.PlantDev.turnover_rates[diagnostics['idevphase'][it]]
            _tr_Leaf[it] = tr_[Plant.PlantDev.ileaf]
            _tr_Root[it] = tr_[Plant.PlantDev.iroot]
            _tr_Stem[it] = tr_[Plant.PlantDev.istem]
            _tr_Seed[it] = tr_[Plant.PlantDev.iseed]
    
    diagnostics['tr_Leaf'] = _tr_Leaf
    diagnostics['tr_Root'] = _tr_Root
    diagnostics['tr_Stem'] = _tr_Stem
    diagnostics['tr_Seed'] = _tr_Seed

    # Add np.nan to the end of each array in the dictionary to represent the last time point in the time_axis (corresponds to the last time point of the state vector)
    for key in diagnostics:
        if key == "t":
            diagnostics[key] = np.append(diagnostics[key], res["t"][-1])
        else:
            diagnostics[key] = np.append(diagnostics[key], np.nan)

    
    # Add state variables to the diagnostics dictionary
    diagnostics["Cleaf"] = res["y"][0,:]
    diagnostics["Cstem"] = res["y"][1,:]
    diagnostics["Croot"] = res["y"][2,:]
    diagnostics["Cseed"] = res["y"][3,:]
    diagnostics["Bio_time"] = res["y"][4,:]
    diagnostics["VRN_time"] = res["y"][5,:]
    diagnostics["Cstate"] = res["y"][7,:]
    diagnostics["Cseedbed"] = res["y"][8,:]

    # Add forcing inputs to diagnostics dictionary
    for i,f in enumerate(forcing_inputs):
        ni = i+1
        if f(time_axis[0]).size == 1:
            fstr = f"forcing {ni:02}"
            diagnostics[fstr] = f(time_axis)
        elif f(time_axis[0]).size > 1:
            # this forcing input has levels/layers (e.g. multilayer soil moisture)
            nz = f(time_axis[0]).size
            for iz in range(nz):
                fstr = f"forcing {ni:02} z{iz}"
                diagnostics[fstr] = f(time_axis)[:,iz]

    return diagnostics

def update_attribute(obj: Any, path: str, new_value: Any) -> None:
    """
    Update the attribute specified by the path with a new value.

    Parameters:
    obj (Any): The main class instance to update.
    path (str): Dot-separated path to the attribute (e.g., "Site.temperature").
    new_value (Any): The new value to set for the attribute.
    """
    attributes = path.split('.')
    # Traverse the path except for the last attribute
    if attributes[0] == "":
        # attribute is in the parent class (i.e. not in a sub-class of the obj)
        setattr(obj, attributes[-1], new_value)
    else:
        for attr in attributes[:-1]:
            obj = getattr(obj, attr)
            # Set the new value to the last attribute in the path
        setattr(obj, attributes[-1], new_value)

def update_attribute_in_phase(obj: Any, path: str, new_value: Any, phase: str) -> None:
    """
    Update the attribute specified by the path and a specific developmental phase with a new value.

    Parameters:
    obj (Any): The main class instance to update.
    path (str): Dot-separated path to the attribute (e.g., "Site.temperature").
    new_value (Any): The new value to set for the attribute.
    phase (str): The plant developmental phase of the attribute
    """
    attributes = path.split('.')
    # Traverse the path except for the last attribute
    for attr in attributes[:-1]:
        obj = getattr(obj, attr)
    # Make a copy of the phase-specific values
    new_phase_values = getattr(obj, attributes[-1])
    # Get the phase index for the value, assuming the obj is the correct module that includes the "phases" attribute
    phase_i = obj.phases.index(phase)
    # Update only the phase-specific value
    new_phase_values[phase_i] = new_value
    # Set the new phase value list to the last attribute in the path
    setattr(obj, attributes[-1], new_phase_values)

def update_and_run_model(param_values, model_instance, input_data, param_info, salib_problem):
    """
    Updates model parameters, runs the biophysical model, and returns the outputs.

    Parameters
    ----------
    param_values : array_like
        Vector of parameter values to be updated in the model.
    model_instance : object
        Instance of the model class that includes methods for running simulations.
    input_data : tuple
        A tuple containing input data required for the model run, structured as:
        (ODEModelSolver, time_axis, forcing_inputs, reset_days, zero_crossing_indices, time_nday_f, time_doy_f, time_year_f).
    param_info : pd.DataFrame
        DataFrame with parameter metadata, including:
        - "Name": Parameter name as used in the model.
        - "Module Path": Path to the model attribute to update (e.g., "Plant.Parameters").
        - "Phase Specific": Boolean indicating whether the parameter is specific to a growth phase.
        - "Phase": (Optional) Index of the growth phase if "Phase Specific" is True.
    salib_problem : dict
        Dictionary required for SALib sensitivity analysis. Includes:
        - "num_vars": Number of parameters.
        - "names": List of parameter names.
        - "bounds": List of parameter bounds.

    Returns
    -------
    array_like
        Model outputs after running the simulation.
    """

    # Unpack input data to pass to model run function
    ODEModelSolver, time_axis, forcing_inputs, reset_days, zero_crossing_indices, time_nday_f, time_doy_f, time_year_f = input_data

    # Update the model instance with the provided parameters
    for idx, value in enumerate(param_values):
        param_name = salib_problem["names"][idx]
        param_path = param_info.loc[param_info["Name"] == param_name, "Module Path"].values[0]
        full_path = f"{param_path}.{param_name}"

        if param_info.loc[param_info["Name"] == param_name, "Phase Specific"].values[0]:
            # Handle phase-specific parameters
            phase = param_info.loc[param_info["Name"] == param_name, "Phase"].values[0]
            update_attribute_in_phase(model_instance, full_path, value, phase)
        else:
            if (param_name == "sowingDays") or (param_name == "harvestDays"):
                # Update parameters that must be defined as a list type
                if isinstance(value, list):
                    update_attribute(model_instance, full_path, value)
                else:
                    update_attribute(model_instance, full_path, [value])
            else:
                # Update regular parameters
                update_attribute(model_instance, full_path, value)

        # Make sure the solver knows about the sowing and harvest dates as well (to reset the state variables like GDD and VD)
        if (param_name == "sowingDays") or (param_name == "harvestDays"):
            # Find value of time_nday_f where time_doy_f == sowingDay and time_year_f == sowingYear.
            sowingDay, sowingYear = model_instance.Management.sowingDays, model_instance.Management.sowingYears
            sowing_nday = time_nday_f[(np.floor(time_doy_f) == sowingDay) & (np.array(time_year_f) == sowingYear)]
            
            # Find value of time_nday_f where time_doy_f == sowingDay and time_year_f == sowingYear.
            harvestDay, harvestYear = model_instance.Management.harvestDays, model_instance.Management.harvestYears
            harvest_nday = time_nday_f[(np.floor(time_doy_f) == harvestDay) & (np.array(time_year_f) == harvestYear)]
            
            # Set reset_days to be the updated sowing and harvest nday
            reset_days = [sowing_nday[0], harvest_nday[0]]

    # Run the model and get the outputs
    model_outputs = run_model_and_get_outputs(
        model_instance, ODEModelSolver, time_axis, forcing_inputs, reset_days, zero_crossing_indices
    )

    return model_outputs