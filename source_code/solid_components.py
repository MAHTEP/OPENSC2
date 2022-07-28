import numpy as np
import os

from UtilityFunctions.auxiliary_functions import (
    get_from_xlsx,
    load_auxiliary_files,
    build_interpolator,
    do_interpolation,
)

# Cu properties
from Properties_of_materials.copper import (
    thermal_conductivity_cu_nist,
    isobaric_specific_heat_cu_nist,
    density_cu,
    electrical_resistivity_cu_nist,
)

# NbTi properties
from Properties_of_materials.niobium_titanium import (
    thermal_conductivity_nbti,
    isobaric_specific_heat_nbti,
    density_nbti,
)

# Nb3Sn properties
from Properties_of_materials.niobium3_tin import (
    thermal_conductivity_nb3sn,
    isobaric_specific_heat_nb3sn,
    density_nb3sn,
)

# RE123 properties
from Properties_of_materials.rare_earth_123 import (
    thermal_conductivity_re123,
    isobaric_specific_heat_re123,
    density_re123,
)

# Stainless steel properties
from Properties_of_materials.stainless_steel import (
    thermal_conductivity_ss,
    isobaric_specific_heat_ss,
    density_ss,
    electrical_resistivity_ss,
)

# Glass-epoxy properties
from Properties_of_materials.glass_epoxy import (
    thermal_conductivity_ge,
    isobaric_specific_heat_ge,
    density_ge,
)

# Silver properties
from Properties_of_materials.silver import (
    thermal_conductivity_ag,
    isobaric_specific_heat_ag,
    density_ag,
    electrical_resistivity_ag,
)

# HASTELLOY - C276 properties
from Properties_of_materials.hastelloy_c276 import (
    thermal_conductivity_hc276,
    isobaric_specific_heat_hc276,
    density_hc276,
    electrical_resistivity_hc276,
)

# Solder Sn60Pb40 properties
from Properties_of_materials.solder_sn60_pb40 import (
    thermal_conductivity_sn60pb40,
    isobaric_specific_heat_sn60pb40,
    density_sn60pb40,
    electrical_resistivity_sn60pb40,
)

# Aluminium properties
from Properties_of_materials.aluminium import (
    thermal_conductivity_al,
    isobaric_specific_heat_al,
    density_al,
    electrical_resistivity_al,
)


class SolidComponents:
    def __init__(self, simulation, s_comp):

        """
        Constructor method of class SolidComponents (cdp, 11/2020)
        """

        # Questa è una bozza, quando e se si dovranno considerare altri flag come \
        # IBFUN o IALPHAB, valutare se è il caso di sfruttare un metodo per \
        # evitare di scrivere lo stesso codice più volte (cdp, 11/2020)
        if s_comp.dict_operation["IQFUN"] > 0:
            # External heating parameters given in file operation.xlsx \
            # (cdp, 11/2020)
            # The heating will be on at some times (cdp, 01/2021)
            s_comp.flag_heating = "On"
            if simulation.transient_input["IADAPTIME"] == 0:
                # Time adaptivity off and (cdp, 11/2020)
                s_comp.dict_num_step["IQFUN"] = dict(
                    ON=int(
                        s_comp.dict_operation["TQBEG"]
                        / simulation.transient_input["STPMIN"]
                    ),
                    OFF=int(
                        s_comp.dict_operation["TQEND"]
                        / simulation.transient_input["STPMIN"]
                    ),
                )
            else:
                # Time adaptivity on
                print("Still to be decided what to do there (cdp, 11/2020)\n")
        # End s_comp.dict_operation["IQFUN"].

    # end method __init__ (cdp, 11/2020)

    def eval_sol_comp_properties(self, dict_obj_inventory, nodal=True):

        """
        Method that evaluate total_density, specific_heat and thermal conductivity of SolidComponents class objects in both nodal points and Gauss points according to **options input parameter (cdp, 07/2020)
        """

        # Properties evaluation in each nodal point (cdp, 07/2020)
        if nodal:
            dict_dummy = self.dict_node_pt
            self.dict_node_pt = self.eval_properties(dict_dummy, dict_obj_inventory)
        # Properties evaluation in each Gauss point (cdp, 07/2020)
        elif nodal == False:
            dict_dummy = self.dict_Gauss_pt
            self.dict_Gauss_pt = self.eval_properties(dict_dummy, dict_obj_inventory)

    def eval_properties(self, dict_dummy, dict_obj_inventory):

        """
        Method that actually evaluate total_density, specific_heat and thermal conductivity of SolidComponents class objects regardless of the location (nodal or Gauss points) (cdp, 07/2020)
        """
        # keys = list(self.dict_input.keys())
        if self.NAME == dict_obj_inventory["MixSCStabilizer"]["Name"]:
            # STR_MIX: stabilizer and superconductor strand (cdp, 07/2020)
            # initialization (cdp, 07/2020)
            rho_num = 0.0
            cp_num = np.zeros(dict_dummy["temperature"].shape)
            kk_num = np.zeros(dict_dummy["temperature"].shape)
            rhoe_num = np.zeros(dict_dummy["temperature"].shape)
            R_stab_non_stab = self.dict_input["STAB_NON_STAB"]

            for ntype in range(self.dict_input["NUM_MATERIAL_TYPES"]):
                if ntype == 0:
                    # stabilizer properties evaluation (cdp, 07/2020)
                    if self.dict_input["ISUPERCONDUCTOR"] != "HTS":
                        # to keep into account that if 3 we know all the crossections \
                        # (cdp, 07/2020)
                        if self.dict_input["ISTABILIZER"] == "Cu":  # Cu (cdp, 07/2020)
                            rho_num = rho_num + density_cu() * R_stab_non_stab
                            cp_num = (
                                cp_num
                                + density_cu()
                                * R_stab_non_stab
                                * isobaric_specific_heat_cu_nist(
                                    dict_dummy["temperature"]
                                )
                            )
                            kk_num = (
                                kk_num
                                + R_stab_non_stab
                                * thermal_conductivity_cu_nist(
                                    dict_dummy["temperature"],
                                    dict_dummy["B_field"],
                                    self.dict_input["RRR"],
                                )
                            )
                            rhoe_num = (
                                rhoe_num
                                + R_stab_non_stab
                                * electrical_resistivity_cu_nist(
                                    dict_dummy["temperature"],
                                    dict_dummy["B_field"],
                                    self.dict_input["RRR"],
                                )
                            )
                        elif (
                            self.dict_input["ISTABILIZER"] == "Al"
                        ):  # Al (cdp, 07/2020)
                            rho_num = rho_num + density_al() * R_stab_non_stab
                            cp_num = (
                                cp_num
                                + density_al()
                                * R_stab_non_stab
                                * isobaric_specific_heat_al(dict_dummy["temperature"])
                            )
                            kk_num = kk_num + R_stab_non_stab * thermal_conductivity_al(
                                dict_dummy["temperature"]
                            )
                            rhoe_num = (
                                rhoe_num
                                + R_stab_non_stab
                                * electrical_resistivity_al(dict_dummy["temperature"])
                            )
                        else:
                            raise ValueError(
                                f"""ERROR: material corresponding to
              {list(self.dict_input.keys())[2]} = 
              {self.dict_input["ISTABILIZER"]} is not defined yet.\n"""
                            )
                    else:  # self.dict_input["ISUPERCONDUCTOR"] == 3 (cdp, 08/2020)
                        Ar_sc = 1.1088e-06
                        Ar_hs = 5.5440e-05
                        Ar_so = 7.7904e-05
                        Ar_ag = 1.1088e-06
                        Ar_tot = Ar_sc + Ar_hs + Ar_so + Ar_ag
                        A_stab = self.dict_input["CROSSECTION"] - Ar_tot
                        if self.dict_input["ISTABILIZER"] == "Cu":  # Cu (cdp, 07/2020)
                            rho_num = rho_num + density_cu() * A_stab
                            cp_num = (
                                cp_num
                                + density_cu()
                                * A_stab
                                * isobaric_specific_heat_cu_nist(
                                    dict_dummy["temperature"]
                                )
                            )
                            kk_num = kk_num + A_stab * thermal_conductivity_cu_nist(
                                dict_dummy["temperature"],
                                dict_dummy["B_field"],
                                self.dict_input["RRR"],
                            )
                            rhoe_num = (
                                rhoe_num
                                + A_stab
                                * electrical_resistivity_cu_nist(
                                    dict_dummy["temperature"],
                                    dict_dummy["B_field"],
                                    self.dict_input["RRR"],
                                )
                            )
                        elif (
                            self.dict_input["ISTABILIZER"] == "Al"
                        ):  # Al (cdp, 07/2020)
                            rho_num = rho_num + density_al() * A_stab
                            cp_num = (
                                cp_num
                                + density_al()
                                * A_stab
                                * isobaric_specific_heat_al(dict_dummy["temperature"])
                            )
                            kk_num = kk_num + A_stab * thermal_conductivity_al(
                                dict_dummy["temperature"]
                            )
                            rhoe_num = rhoe_num + A_stab * electrical_resistivity_al(
                                dict_dummy["temperature"]
                            )
                        else:
                            raise ValueError(
                                f"""ERROR: material corresponding to
              {list(self.dict_input.keys())[2]} = 
              {self.dict_input["ISTABILIZER"]} is not defined yet.\n"""
                            )
                else:  # ntype > 0 (cdp, 08/2020)
                    # superconductor properties evaluation (cdp, 07/2020)
                    if self.dict_input["ISUPERCONDUCTOR"] == "NbTi":
                        # LTS: NbTi (cdp, 07/2020)
                        rho_num = rho_num + density_nbti()
                        cp_num = cp_num + density_nbti() * isobaric_specific_heat_nbti(
                            dict_dummy["temperature"],
                            dict_dummy["B_field"],
                            dict_dummy["T_cur_sharing_min"],
                            dict_dummy["T_critical"],
                        )
                        kk_num = kk_num + thermal_conductivity_nbti(
                            dict_dummy["temperature"]
                        )
                    elif self.dict_input["ISUPERCONDUCTOR"] == "Nb3Sn":
                        # LTS: Nb3Sn (cdp, 07/2020)
                        rho_num = rho_num + density_nb3sn()
                        cp_num = (
                            cp_num
                            + density_nb3sn()
                            * isobaric_specific_heat_nb3sn(
                                dict_dummy["temperature"],
                                dict_dummy["T_cur_sharing_min"],
                                dict_dummy["T_critical"],
                                self.dict_input["Tc0m"],
                            )
                        )
                        kk_num = kk_num + thermal_conductivity_nb3sn(
                            dict_dummy["temperature"]
                        )
                    elif self.dict_input["ISUPERCONDUCTOR"] == "HTS":
                        # HTS: RE123 EU DEMO CS (SPC design, 2016) (cdp, 07/2020)
                        rho_num = rho_num + (
                            density_re123() * Ar_sc
                            + density_hc276() * Ar_hs
                            + density_ag() * Ar_ag
                            + density_sn60pb40() * Ar_so
                        )
                        cp_num = cp_num + (
                            density_re123()
                            * isobaric_specific_heat_re123(dict_dummy["temperature"])
                            * Ar_sc
                            + density_hc276()
                            * isobaric_specific_heat_hc276(dict_dummy["temperature"])
                            * Ar_hs
                            + density_sn60pb40()
                            * isobaric_specific_heat_sn60pb40(dict_dummy["temperature"])
                            * Ar_so
                            + density_ag()
                            * isobaric_specific_heat_ag(dict_dummy["temperature"])
                            * Ar_ag
                        )
                        kk_num = kk_num + (
                            thermal_conductivity_re123(dict_dummy["temperature"])
                            * Ar_sc
                            + thermal_conductivity_hc276(dict_dummy["temperature"])
                            * Ar_hs
                            + thermal_conductivity_sn60pb40(dict_dummy["temperature"])
                            * Ar_so
                            + thermal_conductivity_ag(dict_dummy["temperature"]) * Ar_ag
                        )
                        rhoe_num = rhoe_num + (
                            electrical_resistivity_hc276(dict_dummy["temperature"])
                            * Ar_hs
                            + electrical_resistivity_sn60pb40(dict_dummy["temperature"])
                            * Ar_so
                            + electrical_resistivity_ag(dict_dummy["temperature"])
                            * Ar_ag
                        )
                    elif self.dict_input["ISUPERCONDUCTOR"] == "scaling.dat":
                        # user defined scaling
                        # not implemented
                        # userscaling_margin()
                        print("Not implemented yet\n")
                    else:
                        raise ValueError(
                            f"""ERROR: material corresponding to
            {list(self.dict_input.keys())[6]} = 
            {self.dict_input["ISUPERCONDUCTOR"]} is not defined yet.\n"""
                        )
                # end if ntype (cdp, 07/2020)
            # end for ntype (cdp, 07/2020)
            # STR_MIX properties evaluation (cdp, 07/2020)
            if self.dict_input["ISUPERCONDUCTOR"] != "HTS":
                dict_dummy.update(total_density=rho_num / (1 + R_stab_non_stab))
                dict_dummy.update(
                    total_thermal_conductivity=kk_num / (1 + R_stab_non_stab)
                )
                dict_dummy.update(
                    total_electrical_resistivity=rhoe_num / (1 + R_stab_non_stab)
                )
            else:
                dict_dummy.update(
                    total_density=rho_num / self.dict_input["CROSSECTION"]
                )
                dict_dummy.update(
                    total_thermal_conductivity=kk_num / self.dict_input["CROSSECTION"]
                )
                dict_dummy.update(
                    total_electrical_resistivity=rhoe_num
                    / self.dict_input["CROSSECTION"]
                )
            # This expression is always the same, what change is the way in which \
            # cp_num and rho_num are evaluated (cdp, 07/2020)
            dict_dummy.update(total_isobaric_specific_heat=cp_num / rho_num)
        elif self.NAME == dict_obj_inventory["SuperConductor"]["Name"]:
            # STR_SC: superconductor strand (cdp, 07/2020)
            if self.dict_input["ISUPERCONDUCTOR"] == "NbTi":
                # LTS: NbTi (cdp, 07/2020)
                dict_dummy.update(total_density=density_nbti())
                dict_dummy.update(
                    total_isobaric_specific_heat=isobaric_specific_heat_nbti(
                        dict_dummy["temperature"],
                        dict_dummy["B_field"],
                        dict_dummy["T_cur_sharing_min"],
                        dict_dummy["T_critical"],
                    )
                )
                dict_dummy.update(
                    total_thermal_conductivity=thermal_conductivity_nbti(
                        dict_dummy["temperature"]
                    )
                )
            elif self.dict_input["ISUPERCONDUCTOR"] == "Nb3Sn":
                # LTS: Nb3Sn (cdp, 07/2020)
                dict_dummy.update(total_density=density_nb3sn())
                dict_dummy.update(
                    total_isobaric_specific_heat=isobaric_specific_heat_nb3sn(
                        dict_dummy["temperature"],
                        dict_dummy["T_cur_sharing_min"],
                        dict_dummy["T_critical"],
                        self.dict_input["Tc0m"],
                    )
                )
                dict_dummy.update(
                    total_thermal_conductivity=thermal_conductivity_nb3sn(
                        dict_dummy["temperature"]
                    )
                )
            elif self.dict_input["ISUPERCONDUCTOR"] == "HTS":
                # HTS: RE123 EU DEMO CS (SPC design, 2016) (cdp, 07/2020)
                Ar_sc = 1.1088e-06
                Ar_hs = 5.5440e-05
                Ar_so = 7.7904e-05
                Ar_ag = 1.1088e-06
                Ar_tot = Ar_sc + Ar_hs + Ar_so + Ar_ag
                rho_num = (
                    density_re123() * Ar_sc
                    + density_hc276() * Ar_hs
                    + density_ag() * Ar_ag
                    + density_sn60pb40() * Ar_so
                )
                # Keep in mind that specific heat is averaged on mass (cdp, 08/2020)
                cp_num = (
                    density_re123()
                    * Ar_sc
                    * isobaric_specific_heat_re123(dict_dummy["temperature"])
                    + density_hc276()
                    * Ar_hs
                    * isobaric_specific_heat_hc276(dict_dummy["temperature"])
                    + density_sn60pb40()
                    * Ar_so
                    * isobaric_specific_heat_sn60pb40(dict_dummy["temperature"])
                    + density_ag()
                    * Ar_ag
                    * isobaric_specific_heat_ag(dict_dummy["temperature"])
                )
                dict_dummy.update(total_density=rho_num / Ar_tot)
                dict_dummy.update(total_isobaric_specific_heat=cp_num / rho_num)
                dict_dummy.update(
                    total_thermal_conductivity=(
                        thermal_conductivity_re123(dict_dummy["temperature"]) * Ar_sc
                        + thermal_conductivity_hc276(dict_dummy["temperature"]) * Ar_hs
                        + thermal_conductivity_sn60pb40(dict_dummy["temperature"])
                        * Ar_so
                        + thermal_conductivity_ag(dict_dummy["temperature"]) * Ar_ag
                    )
                    / Ar_tot
                )

                dict_dummy.update(
                    total_electrical_resistivity=(
                        electrical_resistivity_hc276(dict_dummy["temperature"]) * Ar_hs
                        + electrical_resistivity_sn60pb40(dict_dummy["temperature"])
                        * Ar_so
                        + electrical_resistivity_ag(dict_dummy["temperature"]) * Ar_ag
                    )
                    / Ar_tot
                )
            elif self.dict_input["ISUPERCONDUCTOR"] == "scaling.dat":
                # user defined scaling
                # not implemented
                # userscaling_margin()
                print("Not implemented yet\n")
            else:
                raise ValueError(
                    f"""ERROR: material corresponding to
        {list(self.dict_input.keys())[3]} = 
        {self.dict_input["ISUPERCONDUCTOR"]} is not defined yet.\n"""
                )
        elif self.NAME == dict_obj_inventory["Stabilizer"]["Name"]:
            # STR_STAB: stabilizer strand (cdp, 07/2020)
            if self.dict_input["ISTABILIZER"] == "Cu":
                # Cu strand (cdp, 07/2020)
                dict_dummy.update(total_density=density_cu())
                dict_dummy.update(
                    total_isobaric_specific_heat=isobaric_specific_heat_cu_nist(
                        dict_dummy["temperature"]
                    )
                )
                dict_dummy.update(
                    total_thermal_conductivity=thermal_conductivity_cu_nist(
                        dict_dummy["temperature"],
                        dict_dummy["B_field"],
                        self.dict_input["RRR"],
                    )
                )
                dict_dummy.update(
                    total_electrical_resistivity=electrical_resistivity_cu_nist(
                        dict_dummy["temperature"],
                        dict_dummy["B_field"],
                        self.dict_input["RRR"],
                    )
                )
            elif self.dict_input["ISTABILIZER"] == "Al":
                # Al strand (cdp, 07/2020)
                dict_dummy.update(total_density=density_al())
                dict_dummy.update(
                    total_isobaric_specific_heat=isobaric_specific_heat_al(
                        dict_dummy["temperature"]
                    )
                )
                dict_dummy.update(
                    total_thermal_conductivity=thermal_conductivity_al(
                        dict_dummy["temperature"]
                    )
                )
                dict_dummy.update(
                    total_electrical_resistivity=electrical_resistivity_al(
                        dict_dummy["temperature"]
                    )
                )
            else:
                raise ValueError(
                    f"""ERROR: material corresponding to
          {list(self.dict_input.keys())[2]} = 
          {self.dict_input["ISTABILIZER"]} is not defined yet.\n"""
                )
        elif self.NAME == dict_obj_inventory["Jacket"]["Name"]:
            # Z_JKT: jacket (cdp, 07/2020)
            # initialization (cdp, 07/2020)
            self.dict_input.update(CROSSECTION=0.0)
            # self.dict_input["CROSSECTION"] = 0.0
            rho_num = 0.0
            cp_num = np.zeros(dict_dummy["temperature"].shape)
            kk_num = np.zeros(dict_dummy["temperature"].shape)

            rhoe_num = np.zeros(dict_dummy["temperature"].shape)

            if self.dict_input["NUM_MATERIAL_TYPES"] == 1:
                if self.dict_input["IMATERIAL_JK"] != "None":
                    if self.dict_input["IMATERIAL_JK"] == "steinless_steel":
                        # stainless steel (cdp, 07/2020)
                        self.dict_input.update(
                            CROSSECTION=self.dict_input["CROSSECTION_JK"]
                        )
                        rho_num = (
                            rho_num + density_ss() * self.dict_input["CROSSECTION_JK"]
                        )
                        cp_num = (
                            cp_num
                            + density_ss()
                            * isobaric_specific_heat_ss(dict_dummy["temperature"])
                            * self.dict_input["CROSSECTION_JK"]
                        )
                        kk_num = (
                            kk_num
                            + thermal_conductivity_ss(dict_dummy["temperature"])
                            * self.dict_input["CROSSECTION_JK"]
                        )
                        rhoe_num = (
                            rhoe_num
                            + electrical_resistivity_ss(dict_dummy["temperature"])
                            * self.dict_input["CROSSECTION_JK"]
                        )
                    else:
                        raise ValueError(
                            f"ERROR: material corresponding to {list(self.dict_input.keys())[3]} = {self.dict_input['IMATERIAL_JK']} is not defined yet.\n"
                        )
                if self.dict_input["IMATERIAL_IN"] != "None":
                    if self.dict_input["IMATERIAL_IN"] == "glass_epoxy":
                        # Glass-epoxy (cdp, 07/2020)
                        self.dict_input.update(
                            CROSSECTION=self.dict_input["CROSSECTION_IN"]
                        )
                        rho_num = (
                            rho_num + density_ge() * self.dict_input["CROSSECTION_IN"]
                        )
                        cp_num = (
                            cp_num
                            + density_ge()
                            * isobaric_specific_heat_ge(dict_dummy["temperature"])
                            * self.dict_input["CROSSECTION_IN"]
                        )
                        kk_num = (
                            kk_num
                            + thermal_conductivity_ge(dict_dummy["temperature"])
                            * self.dict_input["CROSSECTION_IN"]
                        )
                        # dummy value for electrical resistivity, high value since it is \
                        # an insulator (cdp, 08/2020)
                        rhoe_num = rhoe_num + 1e10 * self.dict_input["CROSSECTION_IN"]
                    else:
                        raise ValueError(
                            f"ERROR: material corresponding to {list(self.dict_input.keys())[4]} = {self.dict_input['IMATERIAL_IN']} is not defined yet.\n"
                        )
            elif self.dict_input["NUM_MATERIAL_TYPES"] == 2:
                for ntype in range(self.dict_input["NUM_MATERIAL_TYPES"]):
                    if ntype == 0:
                        if self.dict_input["IMATERIAL_JK"] == "steinless_steel":
                            # stainless steel (cdp, 07/2020)
                            self.dict_input.update(
                                CROSSECTION=self.dict_input["CROSSECTION"]
                                + self.dict_input["CROSSECTION_JK"]
                            )
                            rho_num = (
                                rho_num
                                + density_ss() * self.dict_input["CROSSECTION_JK"]
                            )
                            cp_num = (
                                cp_num
                                + density_ss()
                                * isobaric_specific_heat_ss(dict_dummy["temperature"])
                                * self.dict_input["CROSSECTION_JK"]
                            )
                            kk_num = (
                                kk_num
                                + thermal_conductivity_ss(dict_dummy["temperature"])
                                * self.dict_input["CROSSECTION_JK"]
                            )
                            rhoe_num = (
                                rhoe_num
                                + electrical_resistivity_ss(dict_dummy["temperature"])
                                * self.dict_input["CROSSECTION_JK"]
                            )
                        else:
                            raise ValueError(
                                f"""ERROR: material corresponding to
              {list(self.dict_input.keys())[3]} = 
              {self.dict_input["IMATERIAL_JK"]} is not defined yet.\n"""
                            )
                    elif ntype == self.dict_input["NUM_MATERIAL_TYPES"] - 1:
                        if self.dict_input["IMATERIAL_IN"] == "glass_epoxy":
                            # Glass-epoxy (cdp, 07/2020)
                            self.dict_input.update(
                                CROSSECTION=self.dict_input["CROSSECTION"]
                                + self.dict_input["CROSSECTION_IN"]
                            )
                            rho_num = (
                                rho_num
                                + density_ge() * self.dict_input["CROSSECTION_IN"]
                            )
                            cp_num = (
                                cp_num
                                + density_ge()
                                * isobaric_specific_heat_ge(dict_dummy["temperature"])
                                * self.dict_input["CROSSECTION_IN"]
                            )
                            kk_num = (
                                kk_num
                                + thermal_conductivity_ge(dict_dummy["temperature"])
                                * self.dict_input["CROSSECTION_IN"]
                            )
                            # dummy value for electrical resistivity, high value since it is \
                            # an insulator (cdp, 08/2020)
                            rhoe_num = (
                                rhoe_num + 1e5 * self.dict_input["CROSSECTION_IN"]
                            )
                        else:
                            raise ValueError(
                                f"""ERROR: material corresponding to
              {list(self.dict_input.keys())[4]} = 
              {self.dict_input["IMATERIAL_IN"]} is not defined yet.\n"""
                            )
                    # end if ntype (cdp, 07/2020)
                # end for ntype (cdp, 07/2020)

            # Jacket properties evaluation
            dict_dummy.update(total_density=rho_num / self.dict_input["CROSSECTION"])
            dict_dummy.update(total_isobaric_specific_heat=cp_num / rho_num)
            dict_dummy.update(
                total_thermal_conductivity=kk_num / self.dict_input["CROSSECTION"]
            )
            dict_dummy.update(
                total_electrical_resistivity=rhoe_num / self.dict_input["CROSSECTION"]
            )
        else:
            raise NameError(
                f"""ERROR: there are no objects with this name:
            {self.NAME}.\n"""
            )

        return dict_dummy

    # end Eval_properties

    def get_current(self, conductor):

        """
        ############################################################################
        #              Get_I(self, conductor)
        ############################################################################
        #
        # Method that initialize electrical current in python objects of class
        # SolidComponents.
        #
        ############################################################################
        # VARIABLE              I/O    TYPE              DESCRIPTION            UNIT
        # --------------------------------------------------------------------------
        # self                  I     object             python object of
        #                                                class SolidComponents     -
        # IOPFUN*               I     scalar integer     flag to decide how to
        #                                                evaluate current:
        #                                                == 0 -> constant;
        #                                                == 1 -> exponential decay;
        #                                                == -1 -> read from
        #                                                I_file_dummy.xlsx         -
        # IOP_TOT*              I     scalar float       total operation
        #                                                current         A
        # TAUDET*               I     scalar float                                s
        # TAUDUM*               I     scalar float                                s
        # IOP0_FRACTION§        I     scalar float      fraction of IOP0
        #                                               transported @ time = 0
        #                                               by self                    -
        # BASE_PATH*              I    string            path of folder Description
        #                                               of Components              -
        # External_current_path* I    string            name of the external
        #                                               file from which
        #                                               interpolate current value  -
        # IOP§                   O    scalar float      component current          A
        #############################################################################
        # Invoched functions/methods: Get_from_xlsx
        #
        ############################################################################
        # * IOPFUN, IOP_TOT, TAUDET, TAUDUM, BASE_PATH and External_current_path are
        # given by conductor.dict_input["IOPFUN"], conductor.IOP_TOT, conductor.dict_input["TAUDET"], conductor.dict_input["TAUDUM"],
        # conductor.BASE_PATH and conductor.file_input["EXTERNAL_CURRENT"].
        # § IOP0_FRACTION and IOP are component attributes: self.dict_node_pt["IOP"]0_FRACTION,
        #  .
        # N.B. IOP is a SolidComponents attribute so its value can be assigned
        # directly, it is a scalar floar since there is no current redistribution.
        ############################################################################
        #
        # Translated and optimized by D. Placido Polito 06/2020
        #
        ############################################################################
        """

        if conductor.dict_input["IOPFUN"] == -1:

            if conductor.cond_time[-1] == 0:
                # Build file path.
                file_path = os.path.join(
                    conductor.BASE_PATH, conductor.file_input["EXTERNAL_CURRENT"]
                )
                # Load auxiliary input file.
                current_df, flagSpecfield = load_auxiliary_files(file_path, sheetname=self.ID)
                # Build interpolator and get the interpolaion flag (space_only,time_only or space_and_time).
                self.current_interpolator, self.current_interp_flag = build_interpolator(
                    current_df, self.dict_operation["IOP_INTERPOLATION"]
                )

            # call load_user_defined_quantity on the component.
            self.dict_node_pt["IOP"] = do_interpolation(
                self.current_interpolator,
                conductor.dict_discretization["xcoord"],
                conductor.cond_time[-1],
                self.current_interp_flag,
            )

            # evaluate IOP_TOT as the sum of first value of current vector \
            # self.dict_node_pt["IOP"] of each SolidComponent. This is because there \
            # is no current redistribution along the conductor \
            # (self.dict_node_pt["IOP"] have se same value in each index).\
            # (cdp, 06/2020)
            conductor.IOP_TOT = conductor.IOP_TOT + self.dict_node_pt["IOP"][0]
            if flagSpecfield == 2:
                print("still to be decided what to do here\n")
        elif conductor.dict_input["IOPFUN"] == 0:
            self.dict_node_pt["IOP"] = (
                conductor.IOP_TOT * self.dict_operation["IOP0_FRACTION"]
            )
        elif conductor.dict_input["IOPFUN"] == 1:
            if conductor.cond_time[-1] <= conductor.dict_input["TAUDET"]:
                self.dict_node_pt["IOP"] = (
                    conductor.IOP_TOT * self.dict_operation["IOP0_FRACTION"]
                )
            elif conductor.cond_time[-1] > conductor.dict_input["TAUDET"]:
                conductor.IOP_TOT = conductor.dict_input["IOP0_TOT"] * np.exp(
                    -(conductor.cond_time[-1] - conductor.dict_input["TAUDET"])
                    / conductor.dict_input["TAUDUM"]
                )
                self.dict_node_pt["IOP"] = (
                    conductor.IOP_TOT * self.dict_operation["IOP0_FRACTION"]
                )
        # Conversion of float to float array if necessary, this avoid following \
        # error: TypeError: 'float' object is not subscriptable (cdp, 08/2020)
        self.dict_node_pt["IOP"] = np.array([self.dict_node_pt["IOP"]])

    # end Get_I

    def get_magnetic_field(self, conductor, nodal=True):
        if nodal:
            # compute B_field in each node (cdp, 07/2020)
            if self.dict_operation["IBIFUN"] < 0:  
                # cza to enable other negative \
                # (read from file) flags -> ibifun.eq.-3, see below (August 29, 2018)

                if conductor.cond_time[-1] == 0:
                    # Build file path.
                    file_path = os.path.join(
                        conductor.BASE_PATH, conductor.file_input["EXTERNAL_BFIELD"]
                    )
                    # Load auxiliary input file.
                    bfield_df,_ = load_auxiliary_files(file_path, sheetname=self.ID)
                    # Build interpolator and get the interpolaion flag (space_only,time_only or space_and_time).
                    self.bfield_interpolator, self.bfield_interp_flag = build_interpolator(
                        bfield_df, self.dict_operation["B_INTERPOLATION"]
                    )

                # call load_user_defined_quantity on the component.
                self.dict_node_pt["B_field"] = do_interpolation(
                    self.bfield_interpolator,
                    conductor.dict_discretization["xcoord"],
                    conductor.cond_time[-1],
                    self.bfield_interp_flag,
                )
                if self.dict_operation["B_field_units"] == "T/A":
                    # BFIELD is per unit of current
                    self.dict_node_pt["B_field"] = (
                        self.dict_node_pt["B_field"] * conductor.IOP_TOT
                    )
                if conductor.dict_input["IOPFUN"] != 0:
                    #### bfield e' un self e' un vettore
                    self.dict_node_pt["B_field"] = (
                        self.dict_node_pt["B_field"]
                        * conductor.IOP_TOT
                        / conductor.dict_input["IOP0_TOT"]
                    )
            elif self.dict_operation["IBIFUN"] == 0:
                self.dict_node_pt["B_field"] = np.linspace(
                    self.dict_operation["BISS"],
                    self.dict_operation["BOSS"],
                    conductor.dict_discretization["N_nod"],
                )
            elif self.dict_operation["IBIFUN"] == 1:
                self.dict_node_pt["B_field"] = np.linspace(
                    self.dict_operation["BISS"],
                    self.dict_operation["BOSS"],
                    conductor.dict_discretization["N_nod"],
                ) + conductor.IOP_TOT / conductor.dict_input["IOP0_TOT"] * np.linspace(
                    self.dict_operation["BITR"],
                    self.dict_operation["BOTR"],
                    conductor.dict_discretization["N_nod"],
                )
        elif nodal == False:
            # compute B_field in each Gauss point (cdp, 07/2020)
            self.dict_Gauss_pt["B_field"] = (
                np.abs(
                    self.dict_node_pt["B_field"][:-1] + self.dict_node_pt["B_field"][1:]
                )
                / 2.0
            )

    # end Get_B_field

    # HERE STARTS THE DEFINITION OF MODULES USEFUL TO INITIALIZE THE DRIVERS FOR \
    # THE EXTERNAL HEATING. D. Placido (06/2020)

    def get_heat(self, conductor):

        """
        Method that evaluates the external heating according to the value of flag IQUFN, thaing unto account the chosen solution method (cdp, 11/2020)
        """

        # START INITIALIZATION (cdp, 10/2020)
        if conductor.cond_num_step == 0:
            if self.dict_operation["IQFUN"] == 0:
                # Initialization is done always in the same way ragardless of the \
                # solution method: a column vector to exploit the array smart notation \
                # in Conductors class method Eval_Gauss_point. It is the only times at \
                # which this method is invoked (cdp, 11/2020)
                self.dict_node_pt["EXTFLX"] = np.zeros(
                    (conductor.dict_discretization["N_nod"], 1)
                )
            else:
                # Initialization is done always in the same way ragardless of the \
                # value of the flag IQFUN and coherently with the chosen solution \
                # algorithm (cdp, 10/2020)
                if (
                    conductor.dict_input["METHOD"] == "BE"
                    or conductor.dict_input["METHOD"] == "CN"
                ):
                    # Backward Euler or Crank-Nicolson (cdp, 10/2020)
                    self.dict_node_pt["EXTFLX"] = np.zeros(
                        (conductor.dict_discretization["N_nod"], 2)
                    )
                elif conductor.dict_input["METHOD"] == "AM4":
                    # Adams-Moulton 4 (cdp, 10/2020)
                    self.dict_node_pt["EXTFLX"] = np.zeros(
                        (conductor.dict_discretization["N_nod"], 4)
                    )
            # end if self.dict_operation["IQFUN"] (cdp, 11/2020)
        # end if conductor.cond_num_step (cdp, 10/2020)
        # END INITIALIZATION (cdp, 10/2020)
        # criteria to decide how to evaluate the external heating (cdp, 10/2020)
        if self.dict_operation["IQFUN"] > 0:
            # invoke method Q0_where to evaluate the external heating according to \
            # the function corresponding to the value of flag IQFUN. It is not \
            # necessary to distinguish according to "Method" options at this point \
            # (cdp, 10/2020)
            if (
                conductor.cond_time[-1] > self.dict_operation["TQBEG"]
                and conductor.cond_time[-1] <= self.dict_operation["TQEND"]
            ):
                self.heat0_where(conductor)
            elif (
                conductor.cond_time[-1] > self.dict_operation["TQEND"]
                and self.flag_heating == "On"
            ):
                self.dict_node_pt["EXTFLX"][:, 0] = 0.0
                self.flag_heating = "Off"
            # end if (cdp, 10/2020)
        elif self.dict_operation["IQFUN"] == -1:
            if conductor.cond_time[-1] == 0:
                # Build file path.
                file_path = os.path.join(
                    conductor.BASE_PATH, conductor.file_input["EXTERNAL_HEAT"]
                )
                # Load auxiliary input file.
                heat_df,_ = load_auxiliary_files(file_path, sheetname=self.ID)
                # Build interpolator and get the interpolaion flag (space_only,time_only or space_and_time).
                self.heat_interpolator, self.heat_interp_flag = build_interpolator(
                    heat_df, self.dict_operation["Q_INTERPOLATION"]
                )

                # compute external heating at conductor initialization calling function do_interpolation.

                self.dict_node_pt["EXTFLX"][:, 0] = do_interpolation(
                    self.heat_interpolator,
                    conductor.dict_discretization["xcoord"],
                    conductor.cond_time[-1],
                    self.heat_interp_flag,
                )
            elif conductor.cond_num_step > 0:
                if conductor.cond_num_step == 1:
                    # Store the old values only immediately after the initializzation, since after that the whole SYSLOD array is saved and there is no need to compute twice the same values.
                    self.dict_node_pt["EXTFLX"][:, 1] = self.dict_node_pt["EXTFLX"][
                        :, 0
                    ].copy()
                # end if conductor.cond_num_step (cdp, 10/2020)
                # call method load_user_defined_quantity to compute heat and overwrite the previous values.
                self.dict_node_pt["EXTFLX"][:, 0] = do_interpolation(
                    self.heat_interpolator,
                    conductor.dict_discretization["xcoord"],
                    conductor.cond_time[-1],
                    self.heat_interp_flag,
                )
            # end if conductor.cond_num_step (cdp, 10/2020)
        elif self.dict_operation["IQFUN"] == -2:
            # AM2 part to be implemented
            if conductor.cond_time[-1] == 0:
                self.dict_node_pt["EXTFLX"][:, 0] = self.user_heat_function()
            elif conductor.cond_num_step > 0:
                if conductor.cond_num_step == 1:
                    # Store the old values only immediately after the initializzation, since after that the whole SYSLOD array is saved and there is no need to compute twice the same values.
                    self.dict_node_pt["EXTFLX"][:, 1] = self.dict_node_pt["EXTFLX"][
                        :, 0
                    ].copy()
                # end if conductor.cond_num_step (cdp, 10/2020)
                # call method load_user_defined_quantity to compute heat and overwrite the previous values.
                self.dict_node_pt["EXTFLX"][:, 0] = self.user_heat_function()
            # end if conductor.cond_num_step (cdp, 10/2020)

        # end self.dict_operation["IQFUN"] (cdp, 10/2020)

    # end Get_Q

    def heat0_where(self, conductor):

        # Compute at each time step since the mesh can change
        lower_bound = np.min(
            np.nonzero(
                conductor.dict_discretization["xcoord"] >= self.dict_operation["XQBEG"]
            )
        )
        upper_bound = np.max(
            np.nonzero(
                conductor.dict_discretization["xcoord"] <= self.dict_operation["XQEND"]
            )
        )
        if self.dict_operation["IQFUN"] == 1:
            # Square wave in time and space (cdp, 11/2020)
            if (
                conductor.dict_input["METHOD"] == "BE"
                or conductor.dict_input["METHOD"] == "CN"
            ):
                # Backward Euler or Crank-Nicolson (cdp, 10/2020)
                if conductor.cond_num_step == 0:
                    # Initialization to Q0 value: this occurs when TQBEG = 0.0 s, i.e. \
                    # heating starts at the beginning of the simulation (cdp, 10/2020)
                    self.dict_node_pt["EXTFLX"][
                        lower_bound : upper_bound + 1, 0
                    ] = self.dict_operation["Q0"]
                elif conductor.cond_num_step > 0:
                    if conductor.cond_num_step == 1:
                        # Store the old values only immediately after the initializzation, \
                        # since after that the whole SYSLOD array is saved and there is no \
                        # need to compute twice the same values (cdp, 10/2020)
                        self.dict_node_pt["EXTFLX"][:, 1] = self.dict_node_pt["EXTFLX"][
                            :, 0
                        ].copy()
                    self.dict_node_pt["EXTFLX"][
                        lower_bound : upper_bound + 1, 0
                    ] = self.dict_operation["Q0"]
                # end if (cdp, 10/2020)
            elif conductor.dict_input["METHOD"] == "AM4":
                # Adams-Moulton 4 (cdp, 10/2020)
                if conductor.cond_num_step == 0:
                    self.dict_node_pt["EXTFLX"][
                        lower_bound : upper_bound + 1, 0
                    ] = self.dict_operation["Q0"]
                    for cc in range(1, 4):
                        # initialize the other columns to the same array: dummy steady \
                        # state (cdp, 10/2020)
                        self.dict_node_pt["EXTFLX"][:, cc] = self.dict_node_pt[
                            "EXTFLX"
                        ][:, 0].copy()
                    # end for cc (cdp, 10/2020)
                elif conductor.cond_num_step > 0:
                    self.dict_node_pt["EXTFLX"][:, 1:4] = self.dict_node_pt["EXTFLX"][
                        :, 0:3
                    ].copy()
                    self.dict_node_pt["EXTFLX"][
                        lower_bound : upper_bound + 1, 0
                    ] = self.dict_operation["Q0"]
                # end if (cdp, 10/2020)
            # end if conductor.dict_input["METHOD"] (cdp, 10/2020)
        # end if self.dict_operation["IQFUN"] (cdp, 11/2020)

    # end Q0_where

    def user_heat_function(self):
        """User defined function to evaluate heat load."""

        # Bad code to evaluate Joule power in joints of feeder CS3U2.
        # Electrical resistivity.
        el_res = 5e-9 # Ohm
        # Length of the joint.
        joint_length = 0.5 # m
        # Compute linear joule power
        return el_res * self.dict_node_pt["IOP"]**2/joint_length

    def jhtflx_new_0(self, conductor):  # tesded: ok (cdp, 06/2020)

        """
        ############################################################################
        #              JHTFLX_new_0(self, conductor)
        ############################################################################
        #
        # Method that initialize to zero the Joule heating flux in python objects
        # of class SolidComponents.
        #
        ############################################################################
        # VARIABLE    I/O    TYPE              DESCRIPTION                      UNIT
        # --------------------------------------------------------------------------
        # self        I      object            python object of
        #                                      class SolidComponents           -
        # xcoord*     I      np array float    conductor spatial
        #                                      discretization                  m
        # JHTFLX      O      np array float    Joule heating flux
        #                                      vector                          W/m^2
        #############################################################################
        # Invoched functions/methods: none
        #
        ############################################################################
        # * xcoord is given by conductor.xcoord
        # N.B. JHTFLX is a SolidComponents attribute so its value can be assigned
        # directly, it has the same of shape of xcoord and it is a np array.
        ############################################################################
        #
        # Author D. Placido Polito 06/2020
        #
        ############################################################################
        """

        # Method JHTFLX_new_0 starts here. (cdp, 06/2020)
        if (
            conductor.dict_input["METHOD"] == "BE"
            or conductor.dict_input["METHOD"] == "CN"
        ):
            # Backward Euler or Crank-Nicolson (cdp, 10/2020)
            if conductor.cond_time[-1] == 0:
                # Initialization (cdp, 10/2020)
                self.dict_node_pt["JHTFLX"] = np.zeros(
                    (conductor.dict_discretization["N_nod"], 2)
                )
            elif conductor.cond_time[-1] > 0:
                if conductor.cond_num_step == 1:
                    # Store the old values only immediately after the initializzation, \
                    # since after that the whole SYSLOD array is saved and there is no \
                    # need to compute twice the same values (cdp, 10/2020)
                    self.dict_node_pt["JHTFLX"][:, 1] = self.dict_node_pt["JHTFLX"][
                        :, 0
                    ].copy()
                # Update value at the current time step (cdp, 10/2020)
                self.dict_node_pt["JHTFLX"][:, 0] = 0.0
            # end if conductor.cond_time[-1] (cdp, 10/2020)
        elif conductor.dict_input["METHOD"] == "AM4":
            # Adams-Moulton 4 (cdp, 10/2020)
            if conductor.cond_time[-1] == 0:
                # Initialization (cdp, 10/2020)
                self.dict_node_pt["JHTFLX"] = np.zeros(
                    (conductor.dict_discretization["N_nod"], 4)
                )
            elif conductor.cond_time[-1] > 0:
                self.dict_node_pt["JHTFLX"][:, 1:4] = self.dict_node_pt["JHTFLX"][
                    :, 0:3
                ].copy()
                # Update value at the current time step (cdp, 10/2020)
                self.dict_node_pt["JHTFLX"][:, 0] = 0.0
                # end if conductor.cond_time[-1] (cdp, 10/2020)
        # end if conductor.dict_input["METHOD"] (cdp, 10/2020)

    # end JHTFLX_new_0

    def set_energy_counters(self, conductor):
        # tesded: ok (cdp, 06/2020)

        """
        ############################################################################
        #              Set_energy_counters(self, conductor)
        ############################################################################
        #
        # Method that initialize to zero the external energy and Joule heating in
        # python objects of class SolidComponents.
        #
        ############################################################################
        # VARIABLE    I/O    TYPE              DESCRIPTION                      UNIT
        # --------------------------------------------------------------------------
        # self        I      object            python object of
        #                                      class SolidComponents           -
        # xcoord*     I      np array float    conductor spatial
        #                                      discretization                  m
        # EEXT        O      np array float    external heating vector         MJ
        # EJHT        O      np array float    external Joule heating
        #                                      vector                          MJ
        #############################################################################
        # Invoched functions/methods: none
        #
        ############################################################################
        # * xcoord is given by conductor.xcoord
        # N.B. EEXT and EJHT are SolidComponents attributes so therir value can be
        # assigned directly they have the same of shape of xcoord and they are np
        # arrays.
        ############################################################################
        #
        # Author D. Placido Polito 06/2020
        #
        ############################################################################
        """

        # Method Set_energy_counters starts here. (cdp, 06/2020)

        if (
            conductor.dict_input["METHOD"] == "BE"
            or conductor.dict_input["METHOD"] == "CN"
        ):
            # Backward Euler or Crank-Nicolson (cdp, 10/2020)
            if conductor.cond_time[-1] == 0:
                # Initialization (cdp, 10/2020)
                self.dict_node_pt["EEXT"] = np.zeros(
                    (conductor.dict_discretization["N_nod"], 2)
                )
                self.dict_node_pt["EJHT"] = np.zeros(
                    (conductor.dict_discretization["N_nod"], 2)
                )
            elif conductor.cond_time[-1] > 0:
                if conductor.cond_num_step == 1:
                    # Store the old values only immediately after the initializzation, \
                    # since after that the whole SYSLOD array is saved and there is no \
                    # need to compute twice the same values (cdp, 10/2020)
                    self.dict_node_pt["EEXT"][:, 1] = self.dict_node_pt["EEXT"][
                        :, 0
                    ].copy()
                self.dict_node_pt["EJHT"][:, 1] = self.dict_node_pt["EJHT"][:, 0].copy()
                # Update value at the current time step (cdp, 10/2020)
                self.dict_node_pt["EEXT"][:, 0] = 0.0
                self.dict_node_pt["EJHT"][:, 0] = 0.0
            # end if conductor.cond_time[-1] (cdp, 10/2020)
        elif conductor.dict_input["METHOD"] == "AM4":
            # Adams-Moulton 4 (cdp, 10/2020)
            if conductor.cond_time[-1] == 0:
                # Initialization (cdp, 10/2020)
                self.dict_node_pt["EEXT"] = np.zeros(
                    (conductor.dict_discretization["N_nod"], 4)
                )
                self.dict_node_pt["EJHT"] = np.zeros(
                    (conductor.dict_discretization["N_nod"], 4)
                )
            elif conductor.cond_time[-1] > 0:
                self.dict_node_pt["EEXT"][:, 1:4] = self.dict_node_pt["EEXT"][
                    :, 0:3
                ].copy()
                self.dict_node_pt["EJHT"][:, 1:4] = self.dict_node_pt["EJHT"][
                    :, 0:3
                ].copy()
                # Update value at the current time step (cdp, 10/2020)
                self.dict_node_pt["EEXT"][:, 0] = 0.0
                self.dict_node_pt["EJHT"][:, 0] = 0.0
            # end if conductor.cond_time[-1] (cdp, 10/2020)
        # end if conductor.dict_input["METHOD"] (cdp, 10/2020)

    # end Set_energy_counters
