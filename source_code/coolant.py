import numpy as np
import warnings
from CoolProp.CoolProp import PropsSI

from fluid_component import FluidComponentInput


class Coolant(FluidComponentInput):
    """docstring for Coolant."""

    KIND = "Coolant"

    def __init__(self, sheet, sheetOpar, dict_file_path, identifier):
        super().__init__(sheet, sheetOpar, dict_file_path, identifier)
        # Kind of the coolant; make it lowercase to be userd as alias in PropsSI.
        self.type = self.inputs["FLUID_TYPE"].lower()
        # Identifier of the coolant.
        self.identifier = f"{self.type}_{identifier.split('_')[1]}"
        # Dictionary dict_node_pt declaration.
        self.dict_node_pt = dict()
        # Dictionary dict_Gauss_pt declaration.
        self.dict_Gauss_pt = dict()
        # Empty dictionary of list to save variable time evolutions at selected spatial coordinates.
        self.time_evol = dict(
            velocity=dict(),
            pressure=dict(),
            temperature=dict(),
            total_density=dict(),
        )
        headers_inl_out = [
            "time (s)",
            "velocity_inl (m/s)",
            "pressure_inl (Pa)",
            "temperature_inl (K)",
            "total_density_inl (kg/m^3)",
            "mass_flow_rate_inl (kg/s)",
            "velocity_out (m/s)",
            "pressure_out (Pa)",
            "temperature_out (K)",
            "total_density_out (kg/m^3)",
            "mass_flow_rate_out (kg/s)",
        ]
        # Empty dictionary of list to save variable time evolutions at inlet and outlet spatial coordinates.
        self.time_evol_io = {key: list() for key in headers_inl_out}
        # Remove key FLUID_TYPE from self.inputs (it becomes attribute of object coolant); removes also for object channel.
        del self.inputs["FLUID_TYPE"]

    # End method __init__.

    def __repr__(self):
        return f"{self.__class__.__name__}(Type: {self.type}, ID: {self.ID})"

    # End method __repr__

    def eval_coolant_density_din_viscosity_gen_flow(self, pressure, temperature):
        """[summary]

        Args:
            pressure ([type]): [description]
            temperature ([type]): [description]
            fluid_mode_info ([type]): [description]
        """
        # Compute density and dynamic viscosity according to the mode (needed to perform the flow initialization)
        density = PropsSI("Dmass", "T", temperature, "P", pressure, self.type)
        viscosity = PropsSI("viscosity", "T", temperature, "P", pressure, self.type)
        return density, viscosity

    # End method eval_coolant_density_din_viscosity_gen_flow

    def eval_reynolds_from_mass_flow_rate(self, mass_flow_rate, din_viscosity):
        """[summary]

        Args:
            mass_flow_rate ([type]): [description]
            din_viscosity ([type]): [description]
            inputs ([type]): [description]
        """
        # Evaluate Reynold number for module Gen_Flow.py purpose for each fluid component
        return (
            mass_flow_rate
            * self.inputs["HYDIAMETER"]
            / (din_viscosity * self.inputs["CROSSECTION"])
        )

    # end method eval_reynolds_from_mass_flow_rate

    def _eval_nodal_pressure_temperature_velocity_initialization(self, conductor):
        # Pressure, temperaure and velocity form initial condition (numpy array since fx is a numpy array); if time is larger than zero, velocity, pressure and temperature came from problem solution after function STEP is invocked.
        # Compute pressure from inlet and outlet values by linear interpolation.
        if self.operations["MDTIN"] >= 0.0:
            # Flow direction from x = 0 to x = L.
            self.dict_node_pt["pressure"] = np.interp(
                conductor.grid_features["zcoord"],
                [0.0, conductor.inputs["ZLENGTH"]],
                [self.operations["PREINL"], self.operations["PREOUT"]],
            )
        else:
            # Flow direction from x = L to x = 0.
            self.dict_node_pt["pressure"] = np.interp(
                conductor.grid_features["zcoord"],
                [0.0, conductor.inputs["ZLENGTH"]],
                [self.operations["PREOUT"], self.operations["PREINL"]],
            )
        # End if self.operations["MDTIN"] >= 0.
        # Compute temperature from inlet and outlet valuesby linear interpolation.
        self.dict_node_pt["temperature"] = np.interp(
            conductor.grid_features["zcoord"],
            [0.0, conductor.inputs["ZLENGTH"]],
            [self.operations["TEMINL"], self.operations["TEMOUT"]],
        )
        # Compute density according to the mode (needed to compute the velocity from mass flow rate)
        self.dict_node_pt["total_density"] = PropsSI(
            "Dmass",
            "T",
            self.dict_node_pt["temperature"],
            "P",
            self.dict_node_pt["pressure"],
            self.type,
        )
        # Compute velocity form mass flow rate, the sing is determined from mass flow rate.
        self.dict_node_pt["velocity"] = self.operations["MDTIN"] / (
            np.maximum(self.inputs["CROSSECTION"], 1e-7)
            * self.dict_node_pt["total_density"]
        )

    # End method _eval_nodal_pressure_temperature_velocity_initialization

    def _eval_gauss_pressure_temperature_velocity(self, conductor):
        """[summary]

        Args:
            conductor ([type]): [description]

        Raises:
            ValueError: [description]
        """
        if bool(self.dict_node_pt):
            # Pressure, temperature and velocity directly evaluated from the value in nodal point, averaging on two consecutives nodes.
            list_average_prop = ["temperature", "pressure", "velocity"]
            # Compurte pressure, temperature and velocity in Gauss points exploiting dictionary comprehension. Remember that old keys in self.dict_Gauss_pt will be deleted and the new self.dict_Gauss_pt have only the keys in list_average_prop (i.e self.dict_Gauss_pt is cleaned and constructed from scratch).
            # N.B. valutare se posso usare np.interp per calcolare pressione e temperatura nel Gauss. Per la velocità capire se posso passare per la densità nel Gauss invertendo la formula della portata come fatto per l'inizializzazione della velocità nei nodi. In questo caso userei come portata il valore medio tra i primi due nodi. Occhio che questa funzione viene utilizzata at ogni timestep non solo all'inizializzazione.
            self.dict_Gauss_pt = {
                key: (
                    self.dict_node_pt[key][: conductor.grid_features["N_nod"] - 1]
                    + self.dict_node_pt[key][1:]
                )
                / 2.0
                for key in list_average_prop
            }
        else:
            # Empty dictionary self.dict_node_pt: remember that nodal properties must be evaluated before the Gauss one.
            raise ValueError(
                f"ERROR! dictionary dict_node_pt is empty. User must first evaluate properties in nodes and than in Gauss points.\n"
            )
        # End if bool(self.dict_node_pt)

    # End method _eval_gauss_pressure_temperature_velocity

    def eval_dimensionless_numbers(self, dict_dummy):
        # Compute Reynolds dimensionless number.
        dict_dummy["Reynolds"] = np.abs(
            self.inputs["HYDIAMETER"]
            * dict_dummy["velocity"]
            * dict_dummy["total_density"]
            / dict_dummy["total_dynamic_viscosity"]
        )
        # Compute Gruneisen dimensionless number.
        dict_dummy["Gruneisen"] = (
            dict_dummy["isobaric_expansion_coefficient"]
            / dict_dummy["isothermal_compressibility"]
            / dict_dummy["total_isochoric_specific_heat"]
            / dict_dummy["total_density"]
        )
        return dict_dummy

    # End method eval_dimensionless_numbers

    def _eval_properties_nodal_gauss(self, conductor, aliases, nodal=True):
        """
        Method that evaluate density, dynamic viscosity, enthalpy, ... of Coolant
        class objects in both nodal points and Gauss points according to nodal
        input parameter exploiting the CoolProp library.
        """
        # Properties evaluation in each nodal point
        if nodal:
            self.dict_node_pt = self._eval_properties(self.dict_node_pt, aliases)
        # Properties evaluation in each Gauss point
        else:
            # Evaluate pressure temperature and velocity in Gauss point (still to be decided where to really put this line of code)
            self._eval_gauss_pressure_temperature_velocity(conductor)
            self.dict_Gauss_pt = self._eval_properties(self.dict_Gauss_pt, aliases)
        # End if nodal

    # End method _eval_properties_nodal_gauss

    def _eval_properties(self, dict_dummy, aliases):
        """
        Method that actually evaluate density, specific_heat and thermal
        conductivity of FluidComponent class objects regardless of the location
        (nodal or Gauss points) (cdp, 09/2020)
        """
        # Evaluate density, dynamic viscosity, Gruneisen, enthalpy, isobaric specific heat, isochoric specific heat, speed of sound and thermal conductivity with MatLab functions and data base
        for prop_name, alias in aliases.items():
            dict_dummy[prop_name] = PropsSI(
                alias,
                "T",
                dict_dummy["temperature"],
                "P",
                dict_dummy["pressure"],
                self.type,
            )
        # End for prop_name, function
        # Compute Reynolds and Prandtl dimensionless number invoking method self.eval_dimensionless_numbers
        dict_dummy = self.eval_dimensionless_numbers(dict_dummy)

        return dict_dummy

    # End method _eval_properties.

    def _compute_density_and_mass_flow_rates_nodal_gauss(self, conductor, nodal=True):
        """
        Method that evaluates channel density and mass flow rate in nodal points or in Gauss point according to options value, after that solution is evaluated with function STEP, to make plots of spatial distribution and time evolution. This function is also invoked in method Conductor.Initialization. (cdp, 10/2020)
        """
        if nodal:
            self.dict_node_pt = self._compute_density_and_mass_flow_rates(
                self.dict_node_pt, conductor
            )
        else:
            # Velocity directly evaluated from the value in nodal point, averaging on two consecutives nodes.
            # Velocity in Gauss points
            self.dict_Gauss_pt["velocity"] = (
                self.dict_node_pt["velocity"][
                    : conductor.grid_features["N_nod"] - 1
                ]
                + self.dict_node_pt["velocity"][1:]
            ) / 2.0
            self.dict_Gauss_pt = self._compute_density_and_mass_flow_rates(
                self.dict_Gauss_pt, conductor
            )

    # End method _compute_density_and_mass_flow_rates_nodal_gauss

    def _compute_density_and_mass_flow_rates(self, dict_dummy, conductor):
        if conductor.cond_num_step > 0:
            # Compute density spatial distribution in nodal points or in Gauss points, only in if conductor.num_time_step > 0 since at the Initialization phase it is evaluated calling method _eval_properties_nodal_gauss
            dict_dummy["total_density"] = PropsSI(
                "Dmass",
                "T",
                dict_dummy["temperature"],
                "P",
                dict_dummy["pressure"],
                self.type,
            )
        # end if conductor.num_time_step > 0
        # Compute mass flow rate spatial distribution
        dict_dummy["mass_flow_rate"] = (
            self.inputs["CROSSECTION"]
            * dict_dummy["velocity"]
            * dict_dummy["total_density"]
        )
        return dict_dummy

    # end method _compute_density_and_mass_flow_rates

    def compute_velocity_gen_flow(
        self,
        ZLENGTH,
        channel,
        MAXITR,
        PDROP,
        RHOREF,
        VISREF,
        TOLRNC,
        velocity_guess=0.0,
        friction_guess=0.01,
        nodal=None,
    ):
        """[summary]

        Args:
            ZLENGTH ([type]): [description]
            channel ([type]): [description]
            MAXITR ([type]): [description]
            PDROP ([type]): [description]
            RHOREF ([type]): [description]
            VISREF ([type]): [description]
            TOLRNC ([type]): [description]
            velocity_guess (float, optional): [description]. Defaults to 0.0.
            friction_guess (float, optional): [description]. Defaults to 0.01.
            nodal ([type], optional): [description]. Defaults to None.

        Returns:
            [type]: [description]
        """
        leff = ZLENGTH / self.inputs["COSTETA"]
        err = 1.0
        # Initialize total friction factor with the guess value (gen_flow sub dictionary)
        channel.dict_friction_factor[nodal]["total"] = friction_guess
        # Initialize velocity with the guess value.
        velocity = velocity_guess
        ii = 1
        # While loop to evaluate the velocity
        while ii <= MAXITR and err >= TOLRNC:
            ii = ii + 1
            old_velocity = velocity
            velocity = np.sqrt(
                (self.inputs["HYDIAMETER"] * PDROP)
                / (2.0 * channel.dict_friction_factor[nodal]["total"] * RHOREF * leff)
            )
            reynolds = (RHOREF * velocity * self.inputs["HYDIAMETER"]) / VISREF
            channel.eval_friction_factor(reynolds, nodal=nodal)
            err = np.fabs(velocity - old_velocity) / np.fabs(velocity)
        # End while ii.
        if ii > MAXITR and err >= TOLRNC:
            warnings.warn(
                f"WARNING. TOLERANCE NOT ACHIEVED FOR {self.ID} IN GENFLW \
      Fuidcomponent coolant compute_velocity:\nerr = {err}"
            )
        # End if ii.
        return velocity

    # End method compute_velocity_gen_flow.

    def compute_mass_flow_with_direction(
        self, FLDIR, RHOREF, VELCT
    ):  # tested: (cdp, 06/2020)
        """[summary]

        Args:
            FLDIR ([type]): [description]
            RHOREF ([type]): [description]
            VELCT ([type]): [description]

        Returns:
            [type]: [description]
        """
        return FLDIR * self.inputs["CROSSECTION"] * RHOREF * VELCT

    # End method compute_mass_flow_with_direction
