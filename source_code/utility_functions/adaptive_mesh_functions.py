import numpy as np
import os
import warnings
from typing import Union

from fluid_component import FluidComponent
from jacket_component import JacketComponent

from conductor import Conductor
from stack_component import StackComponent
from strand_component import StrandComponent
from strand_mixed_component import StrandMixedComponent
from strand_stabilizer_component import StrandStabilizerComponent
from cylindrical_helix import CylindricalHelix

def detect_quench_front(conductor:Conductor)->dict:
    """Function that detects the quench fronts comparing for each StackComponent and StrandMixedComponent the current sharing temperature and its own temperature. Where these temperature crosses each other quench fronts are identified and a local mesh adaptation (refinement/coarsening) may be necessary.
    The index of the quench fronts are stored in a dictionary (front_idx).

    Args:
        conductor (Conductor): object with all information to detect the quench fronts.

    Returns:
        dict: collection of quench fronts indices for each StackComponent and StrandMixedComponent for which the temperature crosses the current sharing temperature. Keys of the dictionary are object identifiers; values of the dictionaries are ndarrays with the quench fronts indices.
    """

    # Inialize data structure collecting index of quench fronts (if any)
    front_idx = dict()

    # Loop on StackComponent to identify the index corresponding to quench 
    # front coordinates. Can be converted into a function/method.
    for obj in conductor.inventory["StackComponent"].collection:
        
        # Evaluate the sign of the difference between the current sharign 
        # temperature and the temperature of the object.
        sign_Tcs_minus_T = np.sign(
            obj.dict_node_pt["T_cur_sharing"] - obj.dict_node_pt["temperature"]
        )
        
        # Compute the elementwise difference of array sign_Tcs_minus_T. Where 
        # this difference is > 0 a left quench front can be identified; where b 
        # < 0 a right quench front can be identified.
        sign_diff = np.sign(sign_Tcs_minus_T[1:] - sign_Tcs_minus_T[:-1])
        
        # Where the value in array sign_diff is != 0 a quench front can be 
        # identified. Get the index of values != 0 in array sign_diff
        front_idx[obj.identifier] = np.nonzero(sign_diff)[0]

    # Loop on StrandMixedComponent to identify the index corresponding to 
    # quench front coordinates. Can be converted into a function/method.
    for obj in conductor.inventory["StrandMixedComponent"].collection:
        
        # Evaluate the sign of the difference between the current sharign 
        # temperature and the temperature of the object.
        sign_Tcs_minus_T = np.sign(
            obj.dict_node_pt["T_cur_sharing"] - obj.dict_node_pt["temperature"]
        )
        
        # Compute the elementwise difference of array sign_Tcs_minus_T. Where 
        # this difference is > 0 a left quench front can be identified; where b 
        # < 0 a right quench front can be identified.
        sign_diff = np.sign(sign_Tcs_minus_T[1:] - sign_Tcs_minus_T[:-1])
        
        # Where the value in array sign_diff is != 0 a quench front can be 
        # identified. Get the index of values != 0 in array sign_diff
        front_idx[obj.identifier] = np.nonzero(sign_diff)[0]

    # Remove keys in front_idx that are not associated with qench fronts 
    # exploiting dictionary comprehension.
    front_idx = {key:value for key,value in front_idx.items() if value.size > 0}

    return front_idx

def eval_gaussian_mesh_density(conductor:Conductor, front_idx:dict)->np.ndarray:
    """Function that evaluates the mesh density according to a gaussian distribution centered around each quench front coordinate. A combination of all these mesh density is returned.

    Args:
        conductor (Conductor): object with all information to evaluate the gaussian mesh density.
        front_idx (dict): dictonary with the quench fronts coordinates for each StackComponent and StrandMixedComponent object as returned from function detect_quench_front

    Returns:
        np.ndarray: updated nparray with the gaussian mesh density centered around the quench front coordinates, to be used as criterion to locally adapt the mesh.
    """

    # Alias
    zcoord = conductor.grid_features["zcoord"]
    zcoord_gauss = conductor.grid_features["zcoord_gauss"]
    sigma = conductor.grid_features["sigma"]
    exp_lim = conductor.grid_features["exp_lim"]
    rho_mesh = conductor.grid_features["rho_mesh"]
    dz_min = conductor.grid_input["SIZMIN"]

    # Loop on each item of dictionary front_idx to evaluate the mesh density 
    # associated to each quench front coordinate.
    for value in front_idx.values():
        for idx in value:

            # Identify the coordinate of the quench front.
            z_front = (zcoord(idx) + zcoord(idx+1)) / 2

            # Compute exponent of the gaussian distribution centered in 
            # z_front; used to evaluate the mesh density
            e_gauss = (zcoord_gauss - z_front) ** 2 / (2 * sigma ** 2)
            
            # Filter the exponent of the gaussian distribution wrt the minimum 
            # accepted value.
            e_gauss = np.minimum(e_gauss,exp_lim)
            
            # Update mesh density; used to understand if the mesh should be 
            # refined, coarsened or left as it is.
            rho_mesh = np.maximum(np.exp(-e_gauss)/dz_min,rho_mesh)

    return rho_mesh
