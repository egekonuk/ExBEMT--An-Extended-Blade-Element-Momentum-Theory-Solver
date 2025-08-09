import neuralfoil as nf
import numpy as np
import json

def get_neuralfoil_aero(dat_file, alpha_deg, reynolds_number, model_size):
    """
    Computes airfoil aerodynamics using NeuralFoil for a single or a list of AoAs.
    Returns a JSON string of a single dictionary or a list of dictionaries.
    """
    try:
        # NeuralFoil can accept a scalar or a list/array for alpha
        aero = nf.get_aero_from_dat_file(
            dat_file,
            alpha=alpha_deg,
            Re=reynolds_number,
            model_size=model_size,
        )

        # Check if the output is for a batch (i.e., CL is a numpy array)
        if isinstance(aero['CL'], np.ndarray):
            # Batch mode: Convert the dictionary of arrays into a list of dictionaries
            # This is a much cleaner format for MATLAB's jsondecode function
            
            # First, convert all numpy arrays in the dictionary to python lists
            for key, value in aero.items():
                aero[key] = value.tolist()
            
            # Then, restructure the data
            num_points = len(aero['CL'])
            list_of_aero_dicts = [dict(zip(aero, t)) for t in zip(*aero.values())]
            
            return json.dumps(list_of_aero_dicts)
            
        else:
            # Single point mode: Convert numpy types to native Python types
            for key, value in aero.items():
                if isinstance(value, (np.float32, np.float64, np.int32, np.int64)):
                    aero[key] = float(value)
            
            return json.dumps(aero)

    except Exception as e:
        return json.dumps({"error": str(e)})