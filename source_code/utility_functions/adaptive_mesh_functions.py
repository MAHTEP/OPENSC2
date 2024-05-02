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

# Alias for flag values
OK_MESH = 0 # no need of mesh refinement/coarsening
REFINE_MESH = 1 # need to refine the mesh
COARSE_MESH = -1 # need to coarse the mesh
HARD_NODE = True # an hard node is a node that belong to the initial mesh
# A soft node is a node added with mesh refinement. False is referred to 
# the name of the variable: hard_node_flag which is True for hard nodes and 
# False for soft nodes.
SOFT_NODE = False

def adaptive_mesh(conductor:Conductor,environment:object)->Conductor:
    """Function that manages the mesh adaptivity according to an algorithm based on quench front detection. When the temperature of an instance of class StackComponent or StrandMixedComponent crosses the current sharing temperature, a quench front is found. The regions in the nearby a quench front are typically characterized by strong gradients and needs mesh refinement. Thus, in those regions the mesh is locally refined adding (soft) nodes. Since the quench front propagates, the added soft nodes may also be removed when they are no longer needed. This process is called mesh coarsening. For more info on the algorithm contact laura.savoldi@polito.it.

    Args:
        conductor (Conductor): object with all information to carry out mesh refinement.
        environment (object): object with all the info that characterize the environment.

    Returns:
        Conductor: instance of class conductor with updated dictionaryes as below.
            * grid_input:
                * NELEMS
            * grid_features:
                * zcoord
                * N_nod
                * hard_node_flag
                * N_nod_lst
                * N_removed_node
            * dict_N_equation:
                * Total
            * dict_Step
                * SYSVAR
        Each conductor component updates the array of the solution (pressure, temperature, velocity for FluidComponent and temperature for SolidComponent) in dictionary dict_node_pt as well as attribute tau on the new mesh exploiting linear interpolation. 
    """

    # Initialize array with the "ideal" shape of the mesh density.
    conductor.grid_features["rho_mesh"] = 1. / conductor.grid_input["SIZMAX"] * np.ones(conductor.grid_input["NELEMS"])

    # Call function detect_quench_front to indentify, for each StackComponent 
    # and StrandMixedComponent, the index of the mesh that correspond to a 
    # quench front.
    front_index = detect_quench_front(conductor)

    # Check if dictionary front index is not empty:
    if front_index:
        # Front index is not empty: perform mesh adaptivity.

        # Call function eval_gaussian_mesh_density to evaluate the mesh density 
        # according to a gaussian distribution centered in each quench front 
        # identified with function detect_quench_front.
        conductor.grid_features["rho_mesh"] = eval_gaussian_mesh_density(
            conductor, front_index
        )

        # Identify regions that need mesh coarsening/refinement and build a new 
        # mesh accordingly.
        conductor.grid_features = update_mesh(conductor)

        # Interpolate conductor solution on the new mesh.
        conductor.interp_solution_on_new_mesh()
        # Update the load therm vector SYSLOD on the new mesh to correctly 
        # solve the next thermal-hydraulic time step.
        conductor.dict_Step["SYSLOD"] = conductor.update_syslod_on_new_mesh()

        # Loop on conductor component to update the angular discretization,
        # used to update the coordinates of the barycenter of each conductor 
        # component.
        for comp in conductor.inventory["all_component"].collection:
            # Interpolate the angular discretization on the new mesh.
            comp.tau = np.interp(
                conductor.grid_features["zcoord_new"],
                conductor.grid_features["zcoord"],
                comp.tau,
            )

        # Update conductor features attributes that are related to the mesh.
        (
            conductor.grid_input,
            conductor.grid_features,
            conductor.dict_N_equation,
            conductor.dict_Step,
        ) = conductor.update_cond_mesh_related_features()

        # Loop on conductor component to update the coordinates of the 
        # barycenter used to build the matrices of inductances and conductances 
        # for the electric module.
        for comp in conductor.inventory["all_component"].collection:
            (
                comp.coordinate["x"],
                comp.coordinate["y"],
                comp.coordinate["z"],
            ) = comp.update_coordinates_of_barycenter(
                conductor.grid_features["N_nod"],
                conductor.grid_features["zcoord"],
                conductor.inventory["StrandComponent"].number,
            )
        
        # Update contact perimeters calling method update_contact_perimeters
        conductor.dict_interf_peri = conductor.update_contact_perimeters(
            environment
        )

        return conductor
    else:
        if (
            conductor.inventory["StackComponent"] == 0
            and conductor.inventory["StrandMixedComponent"] == 0
        ):
            warnings.warn("There are no instances of class StackComponent and of class StrandMixedComponent, therefore the adaptive mesh cannot be exploited. The simulation is carried out with the initial mesh.")
        else:
            warnings.warn("The adaptive mesh is not activated at this thermal hydraulic time step because any of the instances of class StackComponent and StrandMixedComponent do not manifest quench front.")

        return conductor

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
            z_front = (zcoord[idx] + zcoord[idx+1]) / 2

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

def update_mesh(conductor:Conductor)->dict:
    """Function that identifies regions of the mesh that need refinement and regions of the mesh that need coarsening, performing the corresponding adjustment of the mesh.
    Regions requiring coarsening or refinement are identified by comparing the actual mesh density wrt the "ideal" mesh density computed with the eval_gaussian_mesh_density function.
    For each element of the mesh:
        * coarsening is performed by calling the coarse_mesh function;
        * refinement is performed by calling the refine_mesh function;
        * regions that do not require coarsening/refinement are treated with the set_node function.

    Args:
        conductor (Conductor): object with all information to update the mesh according to the result of the comparision of the actual mesh density wrt the "ideal" mesh density from function eval_gaussian_mesh_density.

    Raises:
        ValueError: if the updated number of nodes is larger than the maximum allowed number of nodes defined by the user (MAXNOD).
        ValueError: if flags in array mesh_quality_flag is different from -1 (COARSE_MESH), 0 (OK_MESH) and 1 (REFINE_MESH).

    Returns:
        dict: dictionary grid_features with all the info associated to the new mesh. Updated dictionary key-value pairs:
            * N_nod_new -> the total number of nodes of the new mesh.
            * zcoord_new -> the spatial discretization of the new mesh.
            * hard_node_flag_new -> the list of flags that specify whether a node of the new mesh is hard (True) or soft (False).
            * N_added_node -> total number of added (soft) node in the new mesh due to refinement needs.
            * N_removed_node -> total number of removed (soft) node in the new mesh due to coarsening needs.
    """
    
    # Alias
    cond_id = conductor.identifier
    nelems = conductor.grid_input["NELEMS"]
    nelems_refinement = conductor.grid_input["NELEMS_REFINEMENT"]
    nnode_max = conductor.grid_input["MAXNOD"]
    grid_features = conductor.grid_features
    zcoord = conductor.grid_features["zcoord"]
    nnode = conductor.grid_features["N_nod"]
    rho_mesh = conductor.grid_features["rho_mesh"]

    # Build path to input file conductor_grid.xlsx
    f_path = os.path.join(
            conductor.BASE_PATH,
            conductor.file_input["GRID_DEFINITION"]
        )

    mesh_quality_flag = np.zeros(nelems)
    # Evaluate the actual mesh density
    actual_rho_mesh = 1. / (zcoord[1:] - zcoord[:-1])
    # Evaluate the desidered number of elements in which each element of the 
    # mesh should actually be discretized. If > 1 refinement is needed.
    n_ref = np.round(rho_mesh / actual_rho_mesh)
    # Evaluate the inverse of the desidered number of elements in which each 
    # element of the mesh should actually be distretized. If > 1 coarsening is 
    # needed.
    n_coa = np.round(actual_rho_mesh / rho_mesh)

    # Set mesh_quality_flag = 1 where refinement is needed (n_ref > 1)
    mesh_quality_flag[n_ref > 1.] = REFINE_MESH
    # Compute the number of elements that needs refinement
    n_marked = np.sum(n_ref > 1.)
    
    # Total number of added nodes
    grid_features["N_added_node"].append(n_marked * (nelems_refinement - 1))

    # Compute the total number of nodes that characterize the new mesh.
    tot_node = nnode + grid_features["N_added_node"][-1]

    if tot_node >= nnode_max:
        
        raise ValueError(f"Unable to refine mesh. Total number of nodes is larger than the maximum number of nodes {nnode_max}. Plesae check sheet GRID in input file {f_path} for conductor {cond_id}.\n")

    # Set mesh_quality_flag = -1 where coarsening is needed (n_coa > 1)
    mesh_quality_flag[n_coa > 1.] = COARSE_MESH

    # Check on mesh coarsening: avoid to remove more than 2 consecutive nodes 
    # at a time.
    # Sum triplets of consecutive values of flag in mesh_quality_flag
    mesh_quality_sum = (
        mesh_quality_flag[:-2]
        + mesh_quality_flag[1:-1]
        + mesh_quality_flag[2:]
    )
    # Find index in mesh_quality_sum = -3: it means that there are at least a 
    # triplet of nodes that is going to be removed.
    indx = np.nonzero(mesh_quality_sum == -3)[0]
    # Check if there are index corresponding to mesh_quality_sum = -3
    if indx.size > 0:
        # Array indx is not empty: it means that there is at least a 
        # triplet of nodes that is going to be removed, and this should be 
        # avoided. Set central value of the triplet to 0 in mesh_quality_flag.
        
        # To get the central value of the triplet (second addend) indx should 
        # be increased by 1, since values in indx refers to the first addend in 
        # mesh_quality_flag.
        mesh_quality_flag[indx+1] = OK_MESH

    # Set the first item of zcoord_new to 0 (the first axial coordinates is 
    # always z = 0 m).
    grid_features["zcoord_new"].append(0.0)
    # Set the first item of hard_node_flag_new to 1 (the first axial coordinate 
    # is aways an hard node).
    grid_features["hard_node_flag_new"].append(HARD_NODE)
    # Set to 0 the current value of the number of nodes of the new mesh. The 
    # total number of nodes of the new spatial discretization is computed 
    # iteratively while updating the mesh.
    grid_features["N_nod_new"] = 0

    # Loop on the elements of the current mesh (old mesh).
    for jj, mesh_flag in enumerate(mesh_quality_flag):

        # Check if each element of the current mesh needs refinement, 
        # coarsenign or if it of the proper size.
        if mesh_flag == OK_MESH:
            # Elment j-th of the current mesh does not need refinement or 
            # coarsening: call function set_node.
            grid_features = set_node(
                grid_features,
                conductor.grid_input,
                jj,
            )

        elif mesh_flag == REFINE_MESH:
            # Elment j-th of the current mesh needs refinement: call function 
            # refine_mesh.
            grid_features = refine_mesh(
                grid_features,
                conductor.grid_input,
                jj,
            )
            print("Refined mesh.\n")

        elif mesh_flag == COARSE_MESH:
            # Elment j-th of the current mesh needs coarsening: call function 
            # coarse_mesh.
            grid_features = coarse_mesh(
                grid_features,
                conductor.grid_input,
                jj,
            )

        else:
            raise ValueError(f"Not valid value for mesh quality flag:\n{mesh_flag = }\n")

    return grid_features

def set_node(grid_feat:dict,grid_input:dict,jj:int)->dict:
    """Function that adds a new node in the new mesh according to the info available in the old mesh. The added node corresponds to the upper boundary of the j-th element of the mesh (at index j+1). The node is characterized as soft or hard according to the information stored in key "hard_node_value" of dictionary grid_feat.

    Args:
        grid_feat (dict): dictionary that stores all the features of the mesh.
        grid_input (dict): dictonary that stores all the input values used to build the initial mesh.
        jj (int): index that identify the present element of the mesh that is queried for coarsening or refinement.

        N.B. Input argument grid_input is not used but it is required in order to have the function set_node, refine_mesh and coarse_mesh with the same signature for future refactoring.

    Returns:
        dict: dictionary grid_feat with all the info associated to the new mesh. Updated dictionary key-value pairs:
            * N_nod_new -> the total number of nodes of the new mesh.
            * zcoord_new -> the spatial discretization of the new mesh.
            * hard_node_flag_new -> the list of flags that specify whether a node of the new mesh is hard (True) or soft (False).
    """

    # Update node counter of the new spatial discretization.
    grid_feat["N_nod_new"] += 1
    
    # Insert a node in the new mesh
    grid_feat["zcoord_new"].append(grid_feat["zcoord"][jj+1])
    # Classify the new node according to the classification used in the old 
    # mesh. Possible category: 1 = hard node, 0 = soft node.
    grid_feat["hard_node_flag_new"].append(grid_feat["hard_node_flag"][jj+1])
    
    return grid_feat

def refine_mesh(grid_feat:dict,grid_input:dict,jj:int)->dict:
    """Function that locally refines the mesh splitting the j-th element of the old mesh into N evenly spaced intervals, where N is the number of element for local mesh refinement as defined by the user for variable NELEMS_REFINEMENT of input file conductor_grid.xlsx.
    The set_node function is also called to correctly deal with the end node of the element.

    Args:
        grid_feat (dict): dictionary that stores all the features of the mesh.
        grid_input (dict): dictonary that stores all the input values used to build the initial mesh.
        jj (int): index that identify the present element of the mesh that is queried for coarsening or refinement.

        N.B. Input argument grid_input is not used but it is required in order to have the function set_node, refine_mesh and coarse_mesh with the same signature for future refactoring.

    Returns:
        dict: dictionary grid_feat with all the info associated to the new mesh. Updated dictionary key-value pairs:
            * N_nod_new -> the total number of nodes of the new mesh.
            * zcoord_new -> the spatial discretization of the new mesh.
            * hard_node_flag_new -> the list of flags that specify whether a node of the new mesh is hard (True) or soft (False).
    """

    # Alias
    nelems_refinement = grid_input["NELEMS_REFINEMENT"]
    zcoord = grid_feat["zcoord"]

    # Update the total number of nodes after local refinement.
    grid_feat["N_nod_new"] += nelems_refinement - 1
    
    # Evaluate the nodes to be added in the localli refined element.
    zz_ref_loc = np.linspace(zcoord[jj],zcoord[jj+1],nelems_refinement+1)

    # Loop to refine the mesh.
    for zz in zz_ref_loc[1:-1]:
        # Add a (soft) node to the new mesh.
        grid_feat["zcoord_new"].append(zz)
        # Specify that the new added node is soft (set flag to False).
        grid_feat["hard_node_flag_new"].append(SOFT_NODE)

    # Call function set_node to set the end node of the current interval 
    # as an hard node.
    grid_feat = set_node(grid_feat,grid_input,jj)

    return grid_feat

def coarse_mesh(grid_feat:dict,grid_input:dict,jj:int)->dict:
    """Function that locally coarsen the mesh removing nodes only if they are soft. If nodes are hard, they are not removed as they correspond to the nodes of the initial mesh.
    The set_node function is also called to correctly deal with the end node of the element.

    Args:
        grid_feat (dict): dictionary that stores all the features of the mesh.
        grid_input (dict): dictonary that stores all the input values used to build the initial mesh.
        jj (int): index that identify the present element of the mesh that is queried for coarsening or refinement.

    Returns:
        dict: dictionary grid_feat with all the info associated to the new mesh. Updated dictionary key-value pairs:
            * N_nod_new -> the total number of nodes of the new mesh.
            * zcoord_new -> the spatial discretization of the new mesh.
            * hard_node_flag_new -> the list of flags that specify whether a node of the new mesh is hard (True) or soft (False).
            * N_removed_node -> total number of removed (soft) node in the new mesh due to coarsening needs.
    """

    # Alias
    zcoord = grid_feat["zcoord"]
    hard_node_flag = grid_feat["hard_node_flag"]
    node_flag_new = grid_feat["hard_node_flag_new"][grid_feat["N_nod_new"]]
    node_flag_old = hard_node_flag[jj+1]

    if node_flag_new and node_flag_old:
        # Both nodes are hard, it is not possible to coarsen the mesh. Set the 
        # node corresponding to the upper bound of the current element of the 
        # mesh as an hard node.
        grid_feat = set_node(grid_feat,grid_input,jj)

    elif node_flag_new == False and node_flag_old:
        # Start node is soft and should be removed.
        grid_feat["zcoord_new"].append(zcoord[jj+1])
        # Mark the new node in the new mesh according to the caracterization 
        # used in the old mesh.
        grid_feat["hard_node_flag_new"].append(node_flag_old)
        # Update the counter of the removed nodes.
        grid_feat["N_removed_node"][-1] += 1
        print("Coarsened mesh.\n")

    elif node_flag_new and node_flag_old == False:
        # No action needed: the soft node at the end of the element will be 
        # removed with the next call of function set_node.
        
        # Update the counter of the removed nodes to keep track of the fact 
        # that at the next call of function set_node the soft node will be 
        # removed.
        grid_feat["N_removed_node"][-1] += 1
        print("Coarsened mesh.\n")
    
    elif node_flag_new == False and node_flag_old == False:
        # Start node is soft and should be removed.
        grid_feat["zcoord_new"].append(zcoord[jj+1])
        # Mark the new node in the new mesh according to the caracterization 
        # used in the old mesh.
        grid_feat["hard_node_flag_new"].append(node_flag_old)
        # Update the counter of the removed nodes.
        grid_feat["N_removed_node"][-1] += 1
        print("Coarsened mesh.\n")

    return grid_feat