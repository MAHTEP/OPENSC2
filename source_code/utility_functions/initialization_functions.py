import logging
import numpy as np
import os
import pandas as pd
import warnings
from typing import Union

from fluid_component import FluidComponent
from jacket_component import JacketComponent

# from conductor import Conductor
from stack_component import StackComponent
from strand_component import StrandComponent
from strand_mixed_component import StrandMixedComponent
from strand_stabilizer_component import StrandStabilizerComponent
from cylindrical_helix import CylindricalHelix

logger_discretization = logging.getLogger("opensc2Logger.discretization")


def conductor_spatial_discretization(conductor: object):
    """Function that evaluate the z component of the spatial discretization of conductor.

    Args:
        conductor (object): object whit the information to evaluate the spatial discretization and of which the spatial discretization should be evaluated.
    """

    # COMPUTE THE COORDINATES IN THE FIRST TURN
    if (
        conductor.grid_input["ITYMSH"] == 0
        or conductor.grid_input["ITYMSH"] == 2
        or abs(conductor.grid_input["XEREFI"] - conductor.grid_input["XBREFI"]) <= 1e-3
    ):
        # Consider the case of adaptive, not-initially-refined grid
        # UNIFORM SPACING
        zcoord = uniform_spatial_discretization(conductor)
    # !*LOCALLY REFINED MESH. COMPUTED ON A SINGLE TURN BASIS
    elif conductor.grid_input["ITYMSH"] == 1 or conductor.grid_input["ITYMSH"] == 3:
        zcoord = fixed_refined_spatial_discretization(
            conductor, np.zeros(conductor.grid_features["N_nod"])
        )

    check_grid_features(
        conductor.inputs["ZLENGTH"], zcoord[0], zcoord[-1], conductor.identifier
    )

    conductor.grid_features["zcoord"] = zcoord
    # end function


def uniform_spatial_discretization(conductor: object) -> np.ndarray:
    """Evaluates straight uniform spatial discretization in z direction.

    Args:
        conductor (Conductor): conductor object, has all the information to evaluate the unifrom mesh.

    Returns:
        np.ndarray: array with uniform spatial discretization along z direction of length conductor.grid_features["N_nod"] for straight conductor components.
    """
    return np.linspace(
        0.0, conductor.inputs["ZLENGTH"], conductor.grid_features["N_nod"]
    )


def uniform_angular_discretization(
    conductor: object,
    comp: Union[
        StrandMixedComponent, StrandStabilizerComponent, StackComponent
    ],
) -> np.ndarray:
    """Function that evaluates uniform angular discretization, used for helicoidal geometry.

    Args:
        conductor (Conductor): conductor object, has all the information to evaluate the unifrom mesh.
        comp (Union[StrandMixedComponent, StrandStabilizerComponent, StackComponent]): generic object of for wich the uniform angular discretization should be evaluated.

    Returns:
        np.ndarray: array with uniform angular discretization of length conductor.grid_features["N_nod"] for helicoidal components.
    """
    return np.linspace(
        0.0,
        comp.cyl_helix.winding_number * 2 * np.pi,
        conductor.grid_features["N_nod"],
    )


def fixed_refined_spatial_discretization(
    conductor: object, zcoord: np.ndarray
) -> np.ndarray:
    """Function that evaluate fixed refined spatial discretization along z direction.

    Args:
        conductor (Conductor): conductor object, has all the information to evaluate the fixed refined grid.
        zcoord (np.ndarray): array of length conductor.grid_features["N_nod"] initialized to zeros.

    Raises:
        ValueError: if dz in refined region is lower than minimum dz.
        ValueError: if dz in refined region is larger than maximum dz.

    Returns:
        np.ndarray: array of length conductor.grid_features["N_nod"] with fixed refined spatial discretization for straight conductor components.
    """
    n_elem = dict()

    # total number of elements to be used for coarse region of the mesh
    n_elem["coarse"] = conductor.grid_input["NELEMS"] - conductor.grid_input["NELREF"]
    # number of elements to be used in coarse region left to refined mesh zone
    n_elem["left"] = round(
        (conductor.grid_input["XBREFI"] - 0.0)
        / (
            conductor.inputs["ZLENGTH"]
            - (conductor.grid_input["XEREFI"] - conductor.grid_input["XBREFI"])
        )
        * n_elem["coarse"]
    )
    n_elem["right"] = n_elem["coarse"] - n_elem["left"]

    NOD_ref = conductor.grid_input["NELREF"] + 1

    # Refined zone discretization
    dx_ref = (
        conductor.grid_input["XEREFI"] - conductor.grid_input["XBREFI"]
    ) / conductor.grid_input[
        "NELREF"
    ]  # m refined zone discretization pitch
    if (dx_ref >= conductor.grid_input["SIZMIN"]) and (
        dx_ref <= conductor.grid_input["SIZMAX"]
    ):
        # refined mesh
        zcoord[n_elem["left"] : n_elem["left"] + NOD_ref] = np.linspace(
            conductor.grid_input["XBREFI"], conductor.grid_input["XEREFI"], NOD_ref
        )
    elif dx_ref < conductor.grid_input["SIZMIN"]:
        raise ValueError(
            f"ERROR: {dx_ref=} m in refined zone < {conductor.grid_input['SIZMIN']=} m!!!\n"
        )
    elif dx_ref > conductor.grid_input["SIZMAX"]:
        raise ValueError(
            f"ERROR: {dx_ref=} m in refined zone > {conductor.grid_input['SIZMAX']} m!!!\n"
        )

    if n_elem["left"] > 0:
        # Discretization of coarse region left to refined zone
        dx_try = (conductor.grid_input["XBREFI"] - 0.0) / n_elem["left"]
        dx1 = dx_ref  # dummy to not overwrite dx_ref
        ii = 0
        while (dx_try / dx1 > conductor.grid_input["DXINCRE"]) and (
            ii <= n_elem["left"]
        ):
            ii = ii + 1
            dx = dx1 * conductor.grid_input["DXINCRE"]
            zcoord[n_elem["left"] - ii] = zcoord[n_elem["left"] + 1 - ii] - dx
            dx1 = dx
            dx_try = (zcoord[n_elem["left"] - ii] - 0.0) / (n_elem["left"] - ii)

        zcoord[0 : n_elem["left"] - ii + 1] = np.linspace(
            0.0, zcoord[n_elem["left"] - ii], n_elem["left"] - ii + 1
        )

    if n_elem["right"] > 0:
        # Discretization of coarse region right to refined zone
        dx_try = (
            conductor.inputs["ZLENGTH"] - conductor.grid_input["XEREFI"]
        ) / n_elem["right"]
        dx1 = dx_ref  # dummy to not overwrite dx_ref
        ii = 0
        while (dx_try / dx1 > conductor.grid_input["DXINCRE"]) and (
            ii <= n_elem["right"]
        ):
            ii = ii + 1
            dx = dx1 * conductor.grid_input["DXINCRE"]
            zcoord[n_elem["left"] + conductor.grid_input["NELREF"] + ii] = (
                zcoord[n_elem["left"] + conductor.grid_input["NELREF"] + ii - 1] + dx
            )
            dx1 = dx
            dx_try = (
                conductor.inputs["ZLENGTH"]
                - zcoord[n_elem["left"] + conductor.grid_input["NELREF"] + ii]
            ) / (n_elem["right"] - ii)

        zcoord[
            n_elem["left"]
            + conductor.grid_input["NELREF"]
            + ii : conductor.grid_input["NELEMS"]
            + 1
        ] = np.linspace(
            zcoord[n_elem["left"] + conductor.grid_input["NELREF"] + ii],
            conductor.inputs["ZLENGTH"],
            n_elem["right"] - ii + 1,
        )
    return zcoord


def fixed_refined_angular_discretization(
    conductor: object,
    comp: Union[
        StrandMixedComponent, StrandStabilizerComponent, StackComponent
    ],
    tau: np.ndarray,
) -> np.ndarray:
    """Function that evaluate fixed refined angular discretization, used for hlicoidal geometry.

    Args:
        conductor (Conductor): conductor object, has all the information to evaluate the fixed refined grid.
        comp (Union[ StrandMixedComponent, StrandStabilizerComponent, StackComponent ]): generic object for wich the fixed refined angular discretization should be evaluated.
        tau (np.ndarray): array of length conductor.grid_features["N_nod"] initialized to zeros.

    Raises:
        ValueError: if dtau in refined region is lower than minimum dtau.
        ValueError: if dtau in refined region is larger than maximum dtau.

    Returns:
        np.ndarray: array of length conductor.grid_features["N_nod"] with fixed refined angular discretization for helicoidal conductor components.
    """
    n_elem = dict()
    n_elem["coarse"] = conductor.grid_input["NELEMS"] - conductor.grid_input["NELREF"]
    n_elem["left"] = round(
        (conductor.grid_input["XBREFI"] - 0.0)
        / (
            conductor.inputs["ZLENGTH"]
            - (conductor.grid_input["XEREFI"] - conductor.grid_input["XBREFI"])
        )
        * n_elem["coarse"]
    )
    n_elem["right"] = n_elem["coarse"] - n_elem["left"]

    n_winding = dict()
    # number of windings left to the refined region
    n_winding["left"] = conductor.grid_input["XBREFI"] / (
        2 * np.pi * comp.cyl_helix.reduced_pitch
    )
    # number of windings right to the refined region
    n_winding["right"] = (
        conductor.inputs["ZLENGTH"] - conductor.grid_input["XEREFI"]
    ) / (2 * np.pi * comp.cyl_helix.reduced_pitch)
    # number of windings in the refined region
    n_winding["ref"] = (
        conductor.grid_input["XEREFI"] - conductor.grid_input["XBREFI"]
    ) / (2 * np.pi * comp.cyl_helix.reduced_pitch)

    assert (
        comp.cyl_helix.winding_number
        - (n_winding["left"] + n_winding["ref"] + n_winding["right"])
        <= 1e-15
    )

    dtau_ref = 2 * np.pi * n_winding["ref"] / conductor.grid_input["NELREF"]

    if (dtau_ref >= conductor.grid_input["SIZMIN"]) and (
        dtau_ref <= conductor.grid_input["SIZMAX"]
    ):

        tau_beg = 2 * np.pi * n_winding["left"] + dtau_ref
        tau_end = 2 * np.pi * (n_winding["left"] + n_winding["ref"])
        tau[
            n_elem["left"] + 1 : (n_elem["left"] + 1) + (conductor.grid_input["NELREF"])
        ] = np.linspace(tau_beg, tau_end, conductor.grid_input["NELREF"] + 1)
    elif dtau_ref < conductor.grid_input["SIZMIN"]:
        raise ValueError(
            f"ERROR: {dtau_ref=} m in refined zone < {conductor.grid_input['SIZMIN']=} m!!!\n"
        )
    elif dtau_ref > conductor.grid_input["SIZMAX"]:
        raise ValueError(
            f"ERROR: {dtau_ref=} m in refined zone > {conductor.grid_input['SIZMAX']} m!!!\n"
        )

    if n_elem["left"] > 0:

        d_tau_try = 2 * np.pi * n_winding["left"] / n_elem["left"]
        d_tau1 = dtau_ref
        ii = 0
        while (d_tau_try / d_tau1 > conductor.grid_input["DXINCRE"]) and (
            ii <= n_elem["left"]
        ):
            d_tau = d_tau1 * conductor.grid_input["DXINCRE"]
            tau[n_elem["left"] - ii] = tau[n_elem["left"] + 1 - ii] - d_tau
            d_tau1 = d_tau
            d_tau_try = tau[n_elem["left"] - ii] / (n_elem["left"] - ii - 1)
            ii = ii + 1

        tau_beg = 0.0
        tau_end = tau(n_elem["left"] + 1 - ii)
        tau[1 : n_elem["left"] + 1 - ii] = np.linspace(
            tau_beg, tau_end, n_elem["left"] + 1 - ii
        )

    if n_elem["right"] > 0:

        d_tau_try = 2 * np.pi * n_winding["right"] / n_elem["right"]
        d_tau1 = dtau_ref
        ii = 0
        while (d_tau_try / d_tau1 > conductor.grid_input["DXINCRE"]) and (
            ii <= n_elem["right"]
        ):
            d_tau = d_tau1 * conductor.grid_input["DXINCRE"]
            tau[(n_elem["left"] + 1) + (conductor.grid_input["NELREF"]) + ii] = (
                tau[(n_elem["left"] + 1) + (conductor.grid_input["NELREF"]) + ii - 1]
                + d_tau
            )
            d_tau1 = d_tau
            d_tau_try = (
                2 * np.pi * comp.cyl_helix.winding_number
                - tau[(n_elem["left"] + 1) + (conductor.grid_input["NELREF"]) + ii]
            ) / (n_elem["left"] - ii - 1)
            ii = ii + 1

        tau_beg = tau[(n_elem["left"] + 1) + (conductor.grid_input["NELREF"]) + ii - 1]
        tau_end = 2 * np.pi * comp.cyl_helix.winding_number
        tau[
            (n_elem["left"] + 1) + (conductor.grid_input["NELREF"]) + ii - 1 : -1
        ] = np.linspace(tau_beg, tau_end, n_elem["right"] + 1 - ii + 1)

    return tau


def user_defined_grid(conductor: object):
    """Fuction that loads user defined spatial discretization of the generic conductor component.

    Args:
        conductor (Conductor): conductor object, has all the information to evaluate load and assign user defined grid.

    Notes: does not allow interpolation in space and or in time, i.e. only fixed spatial discretization is available.
    """
    # Build file path.
    file_path = os.path.join(conductor.BASE_PATH, conductor.file_input["EXTERNAL_GRID"])
    # Load all sheets of user defined grid auxiliary input file.
    coord_dfs = pd.read_excel(file_path, sheet_name=None)
    # Check user defined grid features.
    (
        conductor.grid_features["N_nod"],
        conductor.grid_features["zcoord"],
    ) = check_user_defined_grid(coord_dfs, conductor, file_path)
    # Evaluate the number of elements from the checked number of nodes.
    conductor.grid_input["NELEMS"] = conductor.grid_features["N_nod"] - 1
    # Assign spatial discretization coordinates to each conductor component.
    for comp in conductor.inventory["all_component"].collection:
        assign_user_defined_spatial_discretization(comp, coord_dfs[comp.identifier])
        check_grid_features(
            conductor.inputs["ZLENGTH"],
            comp.coordinate["z"][0],
            comp.coordinate["z"][-1],
            comp.identifier,
        )


def check_grid_features(zlength: float, z0: float, z1: float, identifier: str):
    """Function that cheks initial anf final coordinates of the discretization to be consistent with the input values.

    Args:
        zlength (float): straight length of the cable
        z0 (float): first element of the array with the z component of the spatial discretization.
        z1 (float): last element of the array with the z component of the spatial discretization.
        identifier (str): identifier of the object for which the check is make.

    Raises:
        ValueError: if first coordinate is smaller 0.0 (below a tolerance)
        ValueError: if first coordinate is lager than 0.0 (above a tolerance)
        ValueError: if last coordinate is smaller or larger than zlength (within a tollerance).
    """

    tol = 1e-6
    if abs(z0 - 0.0) > tol:
        message = f"{identifier = }: z0 must be 0.0; {z0 = } (m) > {0.0} (m)."
        logger_discretization.error(message)
        raise ValueError(message)

    if abs(z1 - zlength) > tol:
        message = (
            f"{identifier = }: z1 does not equal zlength!\n{z1 = } m; {zlength =} m"
        )
        logger_discretization.warning(message)
        raise ValueError(message)


def check_user_defined_grid(
    dfs: dict, conductor: object, file_path: str
) -> "tuple[int, np.ndarray]":
    """Function that checks the consistency of the user defined spatial discretization for all user defined components.

    Args:
        dfs (dict): dictionary of dataframes, each dataframes has the coordinate of the corresponding coductor component.
        conductor (Conductor): conductor object.
        file_path (str): path of the input file with user defined spatial discretization.

    Raises:
        ValueError: if the number of sheets in used defined auxiliary input file differs from the total number of conductor components defined.
        ValueError: if less than three spatial coordinates are provided (x, y, z).
        ValueError: if the number of nodes is not the same in at least one sheet of the file. The reference number of nodes is inferred from the first sheet of the file.
        ValueError: if the z component of the spatial discretization is not exactly the same for FluidComponent and JacketComponent objects.

    Returns:
        tuple[int, np.ndarray]: total number of nodes of the user defined spatial discretization; z component of the conuctor spatial discretization: it is the z component of the spatial discretization of the first deifined FluidComponent object.
    """
    if len(dfs) != conductor.inventory["all_component"].number:
        raise ValueError(
            f"The number of sheets in file {file_path} must be equal to the number of defined conductor components.\n{len(dfs)} != {conductor.inventory['all_component'].number}.\n"
        )

    comp_ref = conductor.inventory["FluidComponent"].collection[0]
    n_node_ref = check_max_node_number(
        dfs[comp_ref.identifier].shape[0], conductor, file_path, comp_ref.identifier
    )
    z_ref = dfs[comp_ref.identifier]["z [m]"].to_numpy()

    for comp in conductor.inventory["all_component"].collection[1:]:
        if dfs[comp.identifier].shape[1] < 3:
            raise ValueError(
                f"User must provide three coordinates in sheeT {comp.identifier} of input file {file_path}.\n"
            )
        if dfs[comp.identifier].shape[0] != n_node_ref:
            raise ValueError(
                f"Inconsistent number of user defined nodes. The number of nodes in sheet {comp.identifier} of file {file_path} must be equal to the one defined in sheet {comp_ref.identifier} of the same file."
            )

        if not np.allclose(
            dfs[comp.identifier]["z [m]"].to_numpy(), z_ref, equal_nan=True
        ):
            raise ValueError(
                f"User must provide the same z component of the coordinate for FluidComponent, JacketComponent and StrandComponent objects. Please check column z [m] in sheet {comp.identifier} of file {file_path}."
            )

    return n_node_ref, z_ref


def check_max_node_number(
    n_nod: int, conductor: object, file_path: str, *args: str
) -> int:
    """Function that checks if the number of nodes of the spatial discretization is lower than the maximum allowed value.

    Args:
        n_nod (int): user defined number of nodes. Different oringin according to flag ITYMSH; see file_path for further detail.
        conductor (Conductor): conductor object with useful informations to make the check.
        file_path (str): path of the file from which the total number of nodes is evaluated.\n
        * If flag ITYMSH > 0 n_nod is evaluated from the number of elements assigned to conductor in input file conductor_grid;
        * If flag ITYMSH < 0 n_nod is set equal to the number of rows of the first sheet of user defined auxiliary input file used for the spatial discretization.

    Raises:
        ValueError: raises error message if the number of nodes excheedes the maximum allowed number of nodes.

    Returns:
        (int): number of nodes used for the spatial discretization.
    """
    if n_nod > conductor.grid_input["MAXNOD"]:
        if conductor.grid_input["ITYMSH"] >= 0:
            message = f"The number of nodes should not exceed the maximum value. {n_nod = } > {conductor.grid_input['MAXNOD'] = }.\nPlease check {conductor.identifier} in file {file_path}.\n"
        else:
            message = f"The number of nodes should not exceed the maximum value. {n_nod = } > {conductor.grid_input['MAXNOD'] = }.\nPlease check sheet {args[0]} in file {file_path}.\n"
        raise ValueError(message)
    return n_nod


def assign_user_defined_spatial_discretization(
    comp: Union[
        FluidComponent,
        JacketComponent,
        StrandMixedComponent,
        StrandStabilizerComponent,
        StackComponent,
    ],
    df: pd.DataFrame,
):
    """Function that assigns the user defined coordinates in sheet comp.identifier of auxiliary input file with the spatial discretization to the corresponding conductor component.

    Args:
        comp (Union[ FluidComponent, JacketComponent, StrandMixedComponent, StrandStabilizerComponent, StackComponent, ]): generic object for which user defines spatial discretization.
        dfs (pd.DataFrame): dataframe with the coordinate of the corresponding coductor component.
    """
    comp.coordinate["x"] = df["x [m]"].to_numpy()
    comp.coordinate["y"] = df["y [m]"].to_numpy()
    comp.coordinate["z"] = df["z [m]"].to_numpy()


def straight_coordinates(
    cond: object,
    comp: Union[
        FluidComponent,
        JacketComponent,
        StrandMixedComponent,
        StrandStabilizerComponent,
        StackComponent,
    ],
    xb: float,
    yb: float,
):
    """Function that builds the straight coordinates of the barycenter of conductor components.

    Args:
        cond (object): conductor object, has all the information to evaluate load and assign user defined grid.
        comp (Union[ FluidComponent, JacketComponent, StrandMixedComponent, StrandStabilizerComponent, StackComponent, ]): generic object for which the coordinates of the barycenter should be evaluated.
        xb (float): x coordinate of the barycenter of the component (passed in this way due to the workaround in function build_coordinates_of_barycenter).
        yb (float): y coordinate of the barycenter of the component (passed in this way due to the workaround in function build_coordinates_of_barycenter).
    """
    comp.coordinate["x"] = xb * np.ones(cond.grid_features["N_nod"])
    comp.coordinate["y"] = yb * np.ones(cond.grid_features["N_nod"])
    if (
        cond.grid_input["ITYMSH"] == 0
        or cond.grid_input["ITYMSH"] == 2
        or abs(cond.grid_input["XEREFI"] - cond.grid_input["XBREFI"]) <= 1e-3
    ):
        comp.coordinate["z"] = uniform_spatial_discretization(cond)
    elif cond.grid_input["ITYMSH"] == 1 or cond.grid_input["ITYMSH"] == 3:
        comp.coordinate["z"] = fixed_refined_spatial_discretization(
            cond, zcoord=np.zeros(cond.grid_features["N_nod"])
        )


def helicoidal_coordinates(
    cond: object,
    comp: Union[
        StrandMixedComponent,
        StrandStabilizerComponent,
        StackComponent,
    ],
):
    """Function that builds the helicoidal coordinates of the barycenter of conductor components of kind StrandMixedComponent, StrandStabilizerComponent and StackComponent.

    Args:
        cond (object): conductor object, has all the information to evaluate load and assign user defined grid.
        comp (Union[ StrandMixedComponent, StrandStabilizerComponent, StackComponent, ]): generic object for which the coordinates of the barycenter should be evaluated.
    """
    comp.cyl_helix = CylindricalHelix(
        comp.inputs["X_barycenter"],
        comp.inputs["Y_barycenter"],
        cond.inputs["ZLENGTH"],
        comp.inputs["COSTETA"],
    )
    # Evaluate the angular discretization according to the value of flag ITYMSH.
    if (
        cond.grid_input["ITYMSH"] == 0
        or cond.grid_input["ITYMSH"] == 2
        or abs(cond.grid_input["XEREFI"] - cond.grid_input["XBREFI"]) <= 1e-3
    ):
        tau = uniform_angular_discretization(cond, comp)
    elif cond.grid_input["ITYMSH"] == 1 or cond.grid_input["ITYMSH"] == 3:
        tau = fixed_refined_angular_discretization(
            cond, comp, tau=np.zeros(cond.grid_features["N_nod"])
        )
    # Evalute coordinates exploiting the helix parametrization.
    (
        comp.coordinate["x"],
        comp.coordinate["y"],
        comp.coordinate["z"],
    ) = comp.cyl_helix.helix_parametrization(tau)


def build_coordinates_of_barycenter(
    cond: object,
    comp: Union[
        FluidComponent,
        JacketComponent,
        StrandMixedComponent,
        StrandStabilizerComponent,
        StackComponent,
    ],
):
    """Function that builds the coordinates of the barycenter of conductor components.

    Args:
        cond (object): conductor object, has all the information to evaluate load and assign user defined grid.
        comp (Union[ FluidComponent, JacketComponent, StrandMixedComponent, StrandStabilizerComponent, StackComponent, ]): generic object for which the coordinates of the barycenter should be evaluated.

    Raises:
        ValueError: if costetha is not equal to 1 for FluidComponent and JacketComponent.

    Note:
        If user defines only one component of type StrandMixedComponent or StrandStabilizerComponent or StackComponent, the straight spatial discretization is used, i.e function straight_coordinates is called instead of helicoidal_coordinates.
    """
    # This is a workaround not a clean solution.
    if isinstance(comp, FluidComponent):
        costheta = comp.coolant.inputs["COSTETA"]
        xb = comp.coolant.inputs["X_barycenter"]
        yb = comp.coolant.inputs["Y_barycenter"]
    else:
        costheta = comp.inputs["COSTETA"]
        xb = comp.inputs["X_barycenter"]
        yb = comp.inputs["Y_barycenter"]

    if costheta == 1:
        # Evaluate straight spatial coordinates. The coordinates are assinged
        # directly to the component in the function.
        straight_coordinates(cond, comp, xb, yb)
    else:
        if isinstance(comp, StrandComponent):
            if cond.inventory["StrandComponent"].number > 1:
                # Evaluate helicoidal spatial coordinates. The coordinates are assinged directly to the component in the function.
                helicoidal_coordinates(cond, comp)
            elif cond.inventory["StrandComponent"].number == 1:
                # Use straight spatial discretization in this case.
                # Evaluate straight spatial coordinates. The coordinates are
                # assinged directly to the component in the function.
                warnings.warn(
                    "User defined only one strand component, it is considered straight since there is no need to evaluate the inductances."
                )
                straight_coordinates(cond, comp, xb, yb)
        else:
            raise ValueError(
                r"$Cos(\theta)$"
                + f"for {comp.__class__.__name__} must be 1.0; current value {costheta = }\n"
            )

    check_grid_features(
        cond.inputs["ZLENGTH"],
        comp.coordinate["z"][0],
        comp.coordinate["z"][-1],
        comp.identifier,
    )
