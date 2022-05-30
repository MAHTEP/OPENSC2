import numpy as np

# import sys
import warnings
from fluid_component import FluidComponentInput


class Channel(FluidComponentInput):
    """docstring for Channel."""

    KIND = "Channel"

    def __init__(self, sheet, sheetOpar, dict_file_path, identifier):
        """[summary]

        Args:
            sheet ([type]): [description]
            sheetOpar ([type]): [description]
            cond ([type]): [description]
            identifier ([type]): [description]
        """
        super().__init__(sheet, sheetOpar, dict_file_path, identifier)
        # Identifier of the channel.
        self.identifier = f"{self.KIND}_{identifier.split('_')[1]}"
        # Assign the type of channel (hole or bundle) to the attribute self.kind
        self.type = self.inputs["CHANNEL_TYPE"]
        # Delete key/value pair "CHANNEL_TYPE" from self.inputs
        del self.inputs["CHANNEL_TYPE"]
        sign = dict(forward=1.0, backward=-1.0)
        flow_dir = self.operations["FLOWDIR"].lower()
        # Assign the direction of the flow to the attribute self.flow_dir
        self.flow_dir = (flow_dir, sign[flow_dir])
        if self.flow_dir[0] != "forward" and self.flow_dir[0] != "backward":
            raise ValueError(
                f"{self.operations['FLOWDIR']} is not a valid alias for the flag FLOWDIR.\nPlease, check {self.identifier} in sheet CHANNEL of input file conductor_operation.xlsx.\n"
            )
        # Delete key/value pair "FLOWDIR" from self.inputs
        del self.operations["FLOWDIR"]

        # FRICTION FACTOR COEFFICIENT

        # Dictionary with all the methods available for the evaluation of the turbulent friction factor according to the geometry of the channel; build exploiting update() method.
        turbulent_friction_factor_methods = {
            key: self.iter_conductor_hole for key in range(100, 102)
        }
        turbulent_friction_factor_methods.update(
            {key: self.turbulent_friction_with_newton_hole for key in range(102, 110)}
        )
        turbulent_friction_factor_methods.update(
            {
                key: self.incropera_rectangular_duct_demo_turbulent_flow_hole
                for key in range(116, 118)
            }
        )
        turbulent_friction_factor_methods.update(
            {
                key: self.duct_demo_common_memo_turbulent_flow_hole
                for key in range(119, 121)
            }
        )
        turbulent_friction_factor_methods.update(
            {key: self.katheder_correlation_bundle for key in range(204, 206)}
        )
        turbulent_friction_factor_methods.update(
            {
                -99: self.user_defined_turbulent_friction_factor,
                110: self.blasius_formula,
                111: self.bhatti_shah_correlation_turbulent_flow_hole,
                112: self.lhc_magnets_recipe_hole,
                113: self.memorandum_bessette_iter_cs_spiral_hole_y2015,
                114: self.csi_hole_ls_best_fit_icec_hole_y2016,
                115: self.tronza_iter_tf_hole_correlation,
                118: self.flat_spiral_from_cfd_hole,
                121: self.haaland_equation,
                122: self.colebrook_formula,
                123: self.petukhov_correlation,
                124: self.rectangular_duct_enea_hts_cicc_trubulent_flow_hole,
                200: self.tronza_iter_tf_correlation_bundle,
                201: self.nicollet_general_correlation_bundle,
                202: self.wanner_correlation_bundle,
                203: self.kstar_correlation_bundle,
                206: self.darcy_forcheimer_porus_medium_correlation_bundle,
                207: self.lhc_poncet_correlation_bundle,
                208: self.iter_cs_correlation_bundle,
                209: self.dtt_correlation_bundle,
                210: self.demo_tf_hts_correlation_bundle,
                211: self.hts_cl_correlation_bundle,
            }
        )
        # Dictionary with all the methods available for the evaluation of the laminar friction factor according to the geometry of the channel; build exploiting update() method.
        laminar_friction_factor_methods = {
            key: self.general_laminar_correlation for key in range(100, 116)
        }
        laminar_friction_factor_methods.update(
            {key: self.general_laminar_correlation for key in range(200, 211)}
        )
        laminar_friction_factor_methods.update(
            {
                key: self.incropera_rectangular_duct_demo_laminar_flow_hole
                for key in range(116, 118)
            }
        )
        laminar_friction_factor_methods.update(
            {
                -99: self.user_defined_laminar_friction_factor,
                118: self.general_laminar_correlation,
                119: self.rectangular_duct_demo_common_memo_laminar_flow_hole,
                120: self.triangular_duct_demo_common_memo_laminar_flow_hole,
                121: self.haaland_equation,  # the same correlation of the turbulent one.
                122: self.colebrook_formula,  # the same correlation of the turbulent one.
                123: self.petukhov_correlation,  # the same correlation of the turbulent one.
                124: self.rectangular_duct_enea_hts_cicc_trubulent_flow_hole,  # the same correlation of the turbulent one.
                211: self.hts_cl_correlation_bundle,  # the same correlation of the turbulent one.
            }
        )
        # Dictionary with all the methods available fot the evaluation of the total friction factor according to the geometry of the channel; build exploiting the update() method. Remember that the keys of the last dictionary override the keys of the other if already available.
        total_friction_factor_methods = {
            key: self._eval_total_friction_factor_max for key in range(100, 116)
        }
        total_friction_factor_methods.update(
            {key: self._eval_total_friction_factor_L_lim_T for key in range(116, 121)}
        )
        total_friction_factor_methods.update(
            {key: self._eval_total_friction_factor_max for key in range(121, 125)}
        )
        total_friction_factor_methods.update(
            {key: self._eval_total_friction_factor_max for key in range(200, 212)}
        )
        total_friction_factor_methods.update(
            {
                -99: self.user_defined_friction_factor,
                106: self._eval_total_friction_factor_L_lim_T,
                118: self._eval_total_friction_factor_max,
            }
        )
        # Assign to attribute self.turbulent_friction_factor_correlation the correlation to be used to evaluyate the turbulent friction factor of the channel.
        self.turbulent_friction_factor_correlation = turbulent_friction_factor_methods[
            self.inputs["IFRICTION"]
        ]
        # Assign to attribute self.laminar_friction_factor_correlation the correlation to be used to evaluyate the laminar friction factor of the channel.
        self.laminar_friction_factor_correlation = laminar_friction_factor_methods[
            self.inputs["IFRICTION"]
        ]
        self.total_friction_factor_correlation = total_friction_factor_methods[
            self.inputs["IFRICTION"]
        ]
        # Define the list of keys and dummy values to build dictionary self.dict_friction_factor[nodal]. total is set to 0. an not to None to avoid error in multiplication between float and None in method _quantities_initialization.
        list_keys_vals = [("laminar", None), ("turbulent", None), ("total", 0.0)]
        # Declare nested dictionary attribute with laminar, turbulent and total friction factors for evaluation at gen_flow (key None) as well as in nodal (key True) and Gauss (key False) points.
        self.dict_friction_factor = {
            None: {kv[0]: kv[1] for kv in list_keys_vals},
            True: {kv[0]: kv[1] for kv in list_keys_vals},
            False: {kv[0]: kv[1] for kv in list_keys_vals},
        }

        # STEADY STATE HEAT TRANSFER COEFFICIENT

        # Dictionary with all the methods available fot the evaluation of the steady state heat transfer coefficient according to the geometry of the channel; build exploiting the update() method. Remember that the keys of the last dictionary override the keys of the other if already available.
        htc_steady_methods = {
            key: self.dittus_boelter_nusselt_lower_limit
            for key in [1, *range(119, 122)]  # * unpacks the range values.
        }
        htc_steady_methods.update(
            {211: self.correlation_211, 212: self.rectangular_duct_enea_hts_cicc_htc}
        )
        # Assign to attribute self.nusselt_correlation the correlation to be used to evaluyate the steady state heat transfer coefficient of the channel.
        self.nusselt_correlation = htc_steady_methods[
            self.inputs["Flag_htc_steady_corr"]
        ]
        # Define the dictionary for the evaluation of steady heat transfer coefficient of the channel in nodal (key True) and Gauss (key False) points.
        self.dict_htc_steady = {True: None, False: None}
        # Define the dictionary for the evaluation of nusselt number of the channel in nodal (key True) and Gauss (key False) points.
        self.dict_nusselt = {True: None, False: None}

        # For time evolution plots.
        self.time_evol = dict(friction_factor=dict())

    # End method __init__.

    def __repr__(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        return f"{self.__class__.__name__}(kind: {self.KIND}, type: {self.type}, identifier: {self.identifier})"

    # End method __repr__.

    def _eval_ft_newton(
        self, hp, dict_coeff, hd2, goverh, tol=1e-2, frict_guess=1e-2, nodal=True
    ):
        """[summary]

        Args:
            hp ([type]): [description]
            dict_coeff ([type]): [description]
            hd2 ([type]): [description]
            goverh ([type]): [description]
            tol ([type], optional): [description]. Defaults to 1e-2.
            frict_guess ([type], optional): [description]. Defaults to 1e-2.
            nodal (bool, optional): [description]. Defaults to True.
        """
        self.dict_friction_factor[nodal]["turbulent"] = frict_guess * np.ones(hp.shape)
        err = 10 * np.ones(hp.shape)
        # While loop to evaluate turbulent friction factor.
        while np.max(err) >= tol:
            f_turb_old = self.dict_friction_factor[nodal]["turbulent"]
            hpiu = hp * np.sqrt(0.5 * self.dict_friction_factor[nodal]["turbulent"])
            rhpiu = (
                dict_coeff["aa"] * hpiu ** dict_coeff["bb"] * goverh ** dict_coeff["cc"]
            )
            self.dict_friction_factor[nodal]["turbulent"] = (
                2.0 / (rhpiu - 2.5 * np.log10(hd2) - 3.75) ** 2
            )
            err = np.fabs(
                (self.dict_friction_factor[nodal]["turbulent"] - f_turb_old)
                / self.dict_friction_factor[nodal]["turbulent"]
            )
        # End while np.max(err)

    # End method eval_ft_newton.

    def user_defined_turbulent_friction_factor(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # User defined turbulent friction factor. Dummy method at the moment but may be useful in future.
        pass

    # End method user_defined_turbulent_friction_factor.

    def user_defined_laminar_friction_factor(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # User defined laminar friction factor. Dummy method at the moment but may be useful in future.
        pass

    # End method user_defined_laminar_friction_factor.

    def user_defined_friction_factor(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # User should write here friction factor value or correlation. Keep in mind that correlation should be array smart if function of reynolds; it must be an array).
        # If constant multiplication by Friction_multiplayer is done in method eval_friction_factor.
        self.dict_friction_factor[nodal]["total"] = np.ones(reynolds.shape)

    # End method user_defined_friction_factor.

    def _quantities_initialization(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.

        Returns:
            [type]: [description]
        """
        reynolds = np.fabs(reynolds)
        # Laminar friction factor initialization
        self.dict_friction_factor[nodal]["laminar"] = np.zeros(reynolds.shape)
        # Turbulent friction factor initialization
        self.dict_friction_factor[nodal]["turbulent"] = np.zeros(reynolds.shape)
        # Friction factor initialization, use ones instead one zeros to avoid overvriting the guess value if nodal = None (gen_flow)
        self.dict_friction_factor[nodal]["total"] = self.dict_friction_factor[nodal][
            "total"
        ] * np.ones(reynolds.shape)
        return reynolds

    # End method _quantities_initialization.

    def general_laminar_correlation(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # This is a numpy boolean array.
        ind = reynolds >= 1.0e-5
        # Define the dictionary with sutable correlations for laminar friction factor evaluation.
        dict_correlation = {
            True: 16.0 / reynolds[ind],  # Laminar correlation for smooth tube
            False: 16.0 / np.minimum(reynolds[ind], 2000),
        }
        # Evaluate laminar friction factor.
        self.dict_friction_factor[nodal]["laminar"][ind] = dict_correlation[
            self.inputs["IFRICTION"] != 9
        ]

    # End method general_laminar_correlation.

    def _eval_correction_diameter_hole(self, tthick=1e-3):
        """[summary]

        Args:
            tthick ([type], optional): [description]. Defaults to 1e-3.

        Returns:
            [type]: [description]
        """
        doutct = self.inputs["HYDIAMETER"] + 2 * tthick
        # cdb correction for scaling from inner diameter to outer diameter (hole only)
        corr_diam = doutct / self.inputs["HYDIAMETER"]
        return doutct, corr_diam

    # End method _eval_correction_diameter_hole.

    def iter_conductor_hole(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Dictionary with the correlation coefficients.
        # 100: 7/9 ITER spiral (average 6/8 - 8/10), from DDD 11, Magnet: Section 1. Engineering Description, December 2004, pp. 32.
        # 101: 8/10 ITER spiral, from N. Peng,L.Q. Liu, L. Serio, L.Y. Xiong, L. Zhang "Thermo-hydraulic analysis of the gradual cool-down to 80 K of the ITER toroidal field coil", Cryogenics (49), 2009, pp.402-406.
        dict_coeff = {
            100: dict(aa=0.45, bb=-0.034, cc=4.0, dd=5.0),
            101: dict(aa=0.36, bb=-0.038, cc=4.0, dd=5.0),
        }
        # Evaluate the correction diameter for scaling invoking method self._eval_correction_diameter_hole
        corr_diam = self._eval_correction_diameter_hole()[1]
        # Compute turbulent friction factor.
        self.dict_friction_factor[nodal]["turbulent"] = (
            (
                dict_coeff[self.inputs["IFRICTION"]]["aa"]
                * (reynolds / corr_diam) ** (dict_coeff[self.inputs["IFRICTION"]]["bb"])
            )
            / dict_coeff[self.inputs["IFRICTION"]]["cc"]
            / (corr_diam ** dict_coeff[self.inputs["IFRICTION"]]["dd"])
        )

    # End method iter_conductor_hole.

    def turbulent_friction_with_newton_hole(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Dictionary with the correlation coefficients.
        # 102: #w_perc1*sp_pitch(cond.ICOND) # 3.05e-3 value is correct for ITER CS
        # aa = 6.4 #!crb Coefficient for 7/9 mm spiral (bad fit)
        # 106: aa = 6.4 #!crb Coefficient for 7/9 mm spiral (bad fit)
        # 107: aa = 6.4 #!crb Coefficient for 7/9 mm spiral (bad fit)
        dict_coeff = {
            102: dict(aa=6.4, bb=0.1717, cc=-0.3428, gg=3.05e-3, tthick1=1e-3),
            103: dict(aa=11.88, bb=0.039, cc=-0.299, gg=3.0e-3, tthick1=1e-3),
            104: dict(aa=11.88, bb=0.039, cc=-0.299, gg=1.5e-3, tthick1=1.5e-3),
            105: dict(aa=11.88, bb=0.039, cc=-0.299, gg=3.05e-3, tthick1=1e-3),
            106: dict(aa=6.4, bb=0.1717, cc=-0.3428, gg=2.0e-3, tthick1=1e-3),
            107: dict(aa=6.4, bb=0.1717, cc=-0.3428, gg=3.05e-3, tthick1=1e-3),
            108: dict(aa=11.88, bb=0.039, cc=-0.299, gg=2.40e-3, tthick1=1e-3),
            109: dict(aa=11.88, bb=0.039, cc=-0.299, gg=5.30e-3, tthick1=1e-3),
        }
        goverh = (
            dict_coeff[self.inputs["IFRICTION"]]["gg"]
            / dict_coeff[self.inputs["IFRICTION"]]["tthick1"]
        )
        hp = (
            dict_coeff[self.inputs["IFRICTION"]]["tthick1"]
            / self.inputs["HYDIAMETER"]
            * reynolds
        )
        hd2 = (
            2.0
            * dict_coeff[self.inputs["IFRICTION"]]["tthick1"]
            / self.inputs["HYDIAMETER"]
        )
        # Compute turbulent friction factor calling method self._eval_ft_newton
        self._eval_ft_newton(
            hp, dict_coeff[self.inputs["IFRICTION"]], hd2, goverh, nodal=nodal
        )

    # End method turbulent_friction_with_newton_hole.

    def blasius_formula(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        self.dict_friction_factor[nodal]["turbulent"] = np.array(
            (0.316 * reynolds ** -0.25) / 4.0
        )

    # End method blasius_formula.

    def bhatti_shah_correlation_turbulent_flow_hole(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        self.dict_friction_factor[nodal]["turbulent"] = (
            0.00128 + 0.1143 * np.maximum(reynolds, 4000) ** -0.311
        )

    # End method bhatti_shah_correlation_turbulent_flow_hole.

    def lhc_magnets_recipe_hole(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Evaluate turbulent friction factor: conversion to Fanning friction factor included
        # This is a numpy boolean array
        ind = reynolds <= 2.0e3
        self.dict_friction_factor[nodal]["turbulent"][ind] = 64 / reynolds[ind] / 4.0
        # This is a numpy boolean array
        ind = (reynolds > 2.0e3) & (reynolds < 4.0e3)
        self.dict_friction_factor[nodal]["turbulent"][ind] = (
            2.6471e-3 * reynolds[ind] ** 0.3279 / 4.0
        )
        # This is a numpy boolean array
        ind = (reynolds > 4.0e3) & (reynolds < 3.0e4)
        self.dict_friction_factor[nodal]["turbulent"][ind] = (
            0.3194 * reynolds[ind] ** -0.25 / 4.0
        )
        # This is a numpy boolean array
        ind = reynolds >= 3.0e4
        self.dict_friction_factor[nodal]["turbulent"][ind] = (
            1.0 / (1.8 * np.log10(reynolds[ind]) / np.log10(10.0) - 1.64) ** 2 / 4.0
        )

    # End method lhc_magnets_recipe_hole.

    def memorandum_bessette_iter_cs_spiral_hole_y2015(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Evaluate the correction diameter for scaling invoking method self._eval_correction_diameter_hole
        corr_diam = self._eval_correction_diameter_hole()[1]
        # This is a numpy boolean array
        ind = reynolds / corr_diam < 1.5e5
        self.dict_friction_factor[nodal]["turbulent"][ind] = (
            3.160 * (reynolds[ind] / corr_diam) ** -0.249
        )
        # This is a numpy boolean array
        ind = reynolds / corr_diam >= 1.5e5
        self.dict_friction_factor[nodal]["turbulent"][ind] = (
            0.216 * (reynolds[ind] / corr_diam) ** -0.024
        )
        self.dict_friction_factor[nodal]["turbulent"] / 4.0 / corr_diam ** 5

    # End method memorandum_bessette_iter_cs_spiral_hole_y2015.

    def csi_hole_ls_best_fit_icec_hole_y2016(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Evaluate Fanning turbulent friction factor.
        self.dict_friction_factor[nodal]["turbulent"] = 0.0958 * reynolds ** -0.181

    # End method csi_hole_ls_best_fit_icec_hole_y2016.

    def tronza_iter_tf_hole_correlation(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Evaluate correction diameter: different formulation wrt the one in method self._eval_correction_diameter_hole
        corr_diam = 3.18e-3 / self.inputs["HYDIAMETER"]
        # Evaluate turbulent friction factor.
        self.dict_friction_factor[nodal]["turbulent"] = (
            0.25 * 0.3164 * (reynolds / corr_diam) ** -0.25 / corr_diam
        )

    # End method tronza_iter_tf_hole_correlation.

    def _incropera_rectangular_duct_demo_tf_coeff_hole(
        self, reynolds, flag_eval_ccc=True
    ):
        """[summary]

        Args:
            reynolds ([type]): [description]

        Returns:
            [type]: [description]
        """
        # crb For rect thickness of 0.215 mm 2.17*10^ -4 = 2.17d - 04
        Coef = 2.0 / 3.0 + 11.0 / 24.0 * 2.17e-4 / (3.28e-3) * (2 - 2.17e-4 / 3.28e-3)
        # Define the dictionary with the needed coefficients evaluated according to the value of flag_eval_ccc. CCC is required only for laminar friction factor.
        dict_coeff = {True: (reynolds * Coef, 96.0), False: reynolds * Coef}
        return dict_coeff(flag_eval_ccc)

    # End method _incropera_rectangular_duct_demo_tf_coeff_holes.

    def _incropera_rectangular_duct_demo_cs_coeff_hole(
        self, reynolds, flag_eval_ccc=True, tthick1=5e-4
    ):
        """[summary]

        Args:
            reynolds ([type]): [description]
            flag_eval_ccc (bool, optional): [description]. Defaults to True.
            tthick1 ([type], optional): [description]. Defaults to 5e-4.

        Returns:
            [type]: [description]
        """
        # Define dictionary for the geometry evaluation. If flag self.inputs["ISRECTANGULAR"] is True TTHICK1 = SIDE1 and DOUTCT = SIDE2, else use default values.
        dict_geometry = {
            True: dict(TTHICK1=self.inputs["SIDE1"], DOUTCT=self.inputs["SIDE2"]),
            False: dict(
                TTHICK1=tthick1,
                DOUTCT=self._eval_correction_diameter_hole(tthick=tthick1)[0],
            ),
        }
        # Evaluate ggg as TTHICK1/DOUTCT, values of TTHICK1 and DOUTCT according to the boolean value of flag self.inputs["ISRECTANGULAR"].
        ggg = (
            dict_geometry[self.inputs["ISRECTANGULAR"]]["TTHICK1"]
            / dict_geometry[self.inputs["ISRECTANGULAR"]]["DOUTCT"]
        )
        # Define dictionary that correct the value of ggg: key True if ggg > 1.
        dict_ggg = {True: ggg ** -1, False: ggg}
        Coef = 2.0 / 3.0 + 11.0 / 24.0 * dict_ggg[ggg > 1.0] * (
            2.0 - dict_ggg[ggg > 1.0]
        )
        # Define the dictionary with the needed coefficients evaluated according to the value of flag_eval_ccc. CCC is required only for laminar friction factor.
        dict_coeff = {
            True: (reynolds * Coef, self._eval_coefficient_ccc(ggg ** -1)),
            False: reynolds * Coef,
        }
        return dict_coeff(flag_eval_ccc)

    # End method _incropera_rectangular_duct_demo_cs_coeff_holeself.

    def _eval_coefficient_ccc(self, ggg):
        """[summary]

        Args:
            ggg ([type]): [description]

        Raises:
            ValueError: [description]

        Returns:
            [type]: [description]
        """
        # Evaluate CCC only if needed (laminar friction factor)
        if ggg < 1.0:
            # sys.error('error: outside correlation range in frioction')
            raise ValueError("error: outside correlation range in frioction")
        elif (ggg >= 1.0) and (ggg < 1.43):
            dict_coeff = dict(aa=1.0, bb=1.43, jj=57.0, kk=59.0)
        elif (ggg >= 1.43) and (ggg < 2.0):
            dict_coeff = dict(aa=1.43, bb=2.0, jj=59.0, kk=62.0)
        elif (ggg >= 2.0) and (ggg < 3.0):
            dict_coeff = dict(aa=2.0, bb=3.0, jj=62.0, kk=69.0)
        elif (ggg >= 3.0) and (ggg < 4.0):
            dict_coeff = dict(aa=3.0, bb=4.0, jj=69.0, kk=73.0)
        elif (ggg >= 4.0) and (ggg <= 8.0):
            dict_coeff = dict(aa=4.0, bb=8.0, jj=73.0, kk=82.0)
        else:
            dict_coeff = dict(aa=8.0, bb=9.0, jj=96.0, kk=96.0)
            ggg = dict_coeff["bb"]
        # End if ggg.
        # Evaluate CCC by interpolation: CCC = jj + (gg - aa)/(bb - aa)*(kk - jj)
        return np.interp(
            ggg,
            [dict_coeff["aa"], dict_coeff["bb"]],
            [dict_coeff["jj"], dict_coeff["kk"]],
        )

    # End method _eval_coefficient_ccc.

    def incropera_rectangular_duct_demo_laminar_flow_hole(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Define dictionary of methods to evaluate coefficients for the correlation according to IFRICTION VALUE
        dict_methods = {
            116: self._incropera_rectangular_duct_demo_tf_coeff_hole,
            117: self._incropera_rectangular_duct_demo_cs_coeff_hole,
        }
        # Evaluate coefficnets for the correlations.
        re_rect, CCC = dict_methods[self.inputs["IFRICTION"]](reynolds)
        # Evaluate laminar friction factor.
        self.dict_friction_factor[nodal]["laminar"] = CCC / re_rect / 4.0

    # End method incropera_rectangular_duct_demo_laminar_flow_hole.

    def incropera_rectangular_duct_demo_turbulent_flow_hole(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Define dictionary of methods to evaluate coefficients for the correlation according to IFRICTION VALUE
        dict_methods = {
            116: self._incropera_rectangular_duct_demo_tf_coeff_hole,
            117: self._incropera_rectangular_duct_demo_cs_coeff_hole,
        }
        # Evaluate coefficients for the correlations.
        re_rect = dict_methods[self.inputs["IFRICTION"]](reynolds, flag_eval_ccc=False)
        # Evaluate turbulent friction factor.
        self.dict_friction_factor[nodal]["turbulent"] = (
            0.79 * np.log10(re_rect) - 1.64
        ) ** -2

    # End method incropera_rectangular_duct_demo_turbulent_flow_hole.

    def flat_spiral_from_cfd_hole(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Evaluate turbulent friction factor.
        self.dict_friction_factor[nodal]["turbulent"] = (
            0.1687 / 4.0 * reynolds ** -0.1129
        )

    # End method flat_spiral_from_cfd_hole.

    def duct_demo_common_memo_turbulent_flow_hole(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Evaluate turbulent friction factor.
        self.dict_friction_factor[nodal]["turbulent"] = (
            0.00128 + 0.1143 * np.maximum(reynolds, 4000) ** -0.311
        )

    # End method duct_demo_common_memo_turbulent_flow_hole.

    def rectangular_duct_demo_common_memo_laminar_flow_hole(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        alpha = self.inputs["SIDE2"] / self.inputs["SIDE1"]
        # Evaluate laminar friction factor.
        self.dict_friction_factor[nodal]["laminar"] = (
            24.0
            * (
                1.0
                - 1.3553 * alpha
                + 1.9467 * alpha ** 2
                - 1.7012 * alpha ** 3
                + 0.9564 * alpha ** 4
                - 0.2537 * alpha ** 5
            )
            / np.minimum(reynolds, 2000)
        )

    # End method rectangular_duct_demo_common_memo_laminar_flow_hole.

    def triangular_duct_demo_common_memo_laminar_flow_hole(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Evaluate laminar friction factor.
        self.dict_friction_factor[nodal]["laminar"] = 13.33 / np.minimum(reynolds, 2000)

    # End method triangular_duct_demo_common_memo_laminar_flow_hole.

    def haaland_equation(self, reynolds, nodal=True):
        """Method that evaluates the turbulent friction factor for the hole with the Haaland formula.

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Evaluate turbulent friction factor: conversion to Fanning friction factor included.
        self.dict_friction_factor[nodal]["turbulent"] = (
            (
                -1.8
                * np.log10(
                    (self.inputs["Roughness"] / self.inputs["HYDIAMETER"] / 3.7) ** 1.11
                    + (6.9 / reynolds)
                )
            )
            ** -2
        ) / 4.0

    # End method haaland_equation.

    def colebrook_formula(self, reynolds, nodal=True):
        """Method that evaluates the turbulent friction factor for the hole with the Colebrook equation; initial guess from the Haaland equation.

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        err = 10.0 * np.ones(reynolds.shape)
        tol = 1e-5
        max_iter = 100
        ii = 0
        # Guess friction factor from Haaland equation
        self.haaland_equation(reynolds, nodal)
        # Assign the guess value to a dummy variable.
        fric_old = self.dict_friction_factor[nodal]["turbulent"].copy()
        # Loop to invert the Colebrook equation.
        while err.max() > tol and ii < max_iter:
            # Evaluate turbulent friction factor: conversion to Fanning friction factor included.
            self.dict_friction_factor[nodal]["turbulent"] = (
                (
                    -2.0
                    * np.log10(
                        (self.inputs["Roughness"] / self.inputs["HYDIAMETER"] / 3.7)
                        + (2.51 / reynolds / np.sqrt(fric_old))
                    )
                )
                ** -2
            ) / 4.0

            err = abs(fric_old - self.dict_friction_factor[nodal]["turbulent"])

            fric_old = self.dict_friction_factor[nodal]["turbulent"].copy()

            ii += 1
        # End while err.max().

    # End method colebrook_formula.

    def petukhov_correlation(self, reynolds, nodal=True):
        """Method that evaluates the turbulent friction factor with the Petukhov correlation.

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Evaluate turbulent friction factor: conversion to Fanning friction factor included.
        self.dict_friction_factor[nodal]["turbulent"] = (
            0.790 * np.log10(reynolds) - 1.64
        ) ** -2.0 / 4.0
        # This is a boolean numpy array.
        ind = reynolds >= 3e3 and reynolds <= 5e6
        message = {
            True: self.do_noting(),
            False: warnings.warn(
                "Reinolds outside the validity range of the Petukhov correlation: friction factor may be badly evaluated."
            ),
        }
        message[ind.all() == True]

    # End method petukhov_correlation

    def rectangular_duct_enea_hts_cicc_trubulent_flow_hole(self, reynolds, nodal=True):
        """Method that evaluates the turbulent friction factor with the Blasius formula if Reynold <= 1e4 and with another correlation if Reynolds > 1e4
        (given by prof L. Savoldy, no reference)

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Evaluate turbulent friction factor: conversion to Fanning friction factor included.
        self.blasius_formula(reynolds, nodal)
        # This is a boolean numpy array.
        self.dict_friction_factor[nodal]["turbulent"][reynolds > 1e4] = (
            2.21 * reynolds[reynolds > 1e4] ** (-0.4) / 4.0
        )

    # End method rectangular_duct_enea_hts_cicc_trubulent_flow_hole

    def tronza_iter_tf_correlation_bundle(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Evaluate the doutct invoking method self._eval_correction_diameter_hole.
        doutct = self._eval_correction_diameter_hole()[0]
        # Evaluate turbulent friction factor.
        self.dict_friction_factor[nodal]["turbulent"] = (
            2.46
            * (reynolds / self.inputs["HYDIAMETER"] * 0.25e-3) ** -0.52
            * self.inputs["HYDIAMETER"]
            / 0.25e-3
            * (
                self.inputs["CROSSECTION"]
                / (0.3 * (39.8e-3 ** 2 - doutct ** 2) * 0.25 * np.pi)
            )
            ** 2
            / 4.0
        )

    # End method tronza_iter_tf_correlation_bundle.

    def nicollet_general_correlation_bundle(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Set the corr_per value: LAURA WILL SET VALUE!
        corr_peri = 6.0 / 5.0
        # Evaluate turbulent friction factor.
        self.dict_friction_factor[nodal]["turbulent"] = (
            (0.0231 + 19.5 / (reynolds / corr_peri) ** 0.7953)
            / (self.inputs["VOID_FRACTION"] ** 0.742)
            / 4.0
        )

    # End method nicollet_general_correlation_bundle

    def wanner_correlation_bundle(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # correlation from paper "Pressure Drop of the W7-X Cable-in-Conduit Conductor" M. Wanner
        # Evaluate turbulent friction factor.
        self.dict_friction_factor[nodal]["turbulent"] = 6.1 * reynolds ** -0.51 / 4.0

    # End method wanner_correlation_bundle.

    def kstar_correlation_bundle(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Deduced from experimental data (Private Communications with KSTAR Korean people); validity range: 2500 < Re < 6500
        # Evaluate turbulent friction factor.
        self.dict_friction_factor[nodal]["turbulent"] = 0.4335 * reynolds ** -0.263

    # End method kstar_correlation_bundle.

    def katheder_correlation_bundle(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # 204: J. Li, Q. Wang, J. Li, Y. Wu, J. Qian, the EAST team, "Thermal-Hydraulic Analysis of PF Coils During Plasma Discharges on EAST", J Supercond Nov Magn, vol. 25, pp. 2033-2039, 2012, DOI 10.1007/s10948-012-1557-6
        # Fitting formula for measurements with pressurized nitrogen based on Katheder correlation.
        # 205: KATHEDER, TO BE USED WITH THE WETTED PERIMETER ACCOUNTING FOR 5/6...
        # Define the dictionary of coefficients
        dict_coeff = {204: dict(aa=0.843, bb=0.0265), 205: dict(aa=0.88, bb=0.051)}
        # Divided by 4 by rb to obtain the Fanning friction factor
        # Evaluate turbulent friction factor.
        self.dict_friction_factor[nodal]["turbulent"] = (
            self.inputs["VOID_FRACTION"] ** 0.72
            * (
                19.5 / reynolds ** dict_coeff[self.inputs["IFRICTION"]]["aa"]
                + dict_coeff[self.inputs["IFRICTION"]]["bb"]
            )
            / 4.0
        )

    # End method katheder_correlation_bundle

    def _eval_jj_and_kk(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        # Define the dictionary of coefficients.
        dict_coeff = {206: (19.6e-9, 2.42, 5.80), 209: (20.9e-9, 19.1, 4.23)}
        # Evaluate jj.
        jj = (
            dict_coeff[self.inputs["IFRICTION"]][1]
            / self.inputs["VOID_FRACTION"] ** dict_coeff[self.inputs["IFRICTION"]][2]
        )
        # Evaluate kk.
        kk = (
            dict_coeff[self.inputs["IFRICTION"]][0]
            * (self.inputs["VOID_FRACTION"] ** 3)
            / ((1 - self.inputs["VOID_FRACTION"]) ** 2)
        )
        return jj, kk

    # End method _eval_jj_and_kk

    def darcy_forcheimer_porus_medium_correlation_bundle(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # TURBULENT CORRELATION Darcy-Forcheimer equation for the flow in porous \
        # media
        # Evaluate jj and kk invocking method _eval_jj_and_kk.
        jj, kk = self._eval_jj_and_kk()
        # Evaluate turbulent friction factor.
        self.dict_friction_factor[nodal]["turbulent"] = (
            self.inputs["HYDIAMETER"] ** 2
            * self.inputs["VOID_FRACTION"]
            / 2.0
            / kk
            / reynolds
            + (self.inputs["HYDIAMETER"] * self.inputs["VOID_FRACTION"] ** 2) / 2.0 * jj
        )

    # End method darcy_forcheimer_porus_medium_correlation_bundle.

    def lhc_poncet_correlation_bundle(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Evaluate turbulent friction factor: conversion to Fanning friction factor included.
        # This is a boolean numpy array.
        ind = (reynolds >= 2e3) & (reynolds < 4e3)
        self.dict_friction_factor[nodal]["turbulent"][ind] = (
            2.6471e-3 * reynolds[ind] ** 0.3279 / 4.0
        )
        # This is a boolean numpy array.
        ind = (reynolds >= 4e3) & (reynolds < 3e4)
        self.dict_friction_factor[nodal]["turbulent"][ind] = (
            0.3194 * reynolds[ind] ** -0.25 / 4.0
        )
        # This is a boolean numpy array.
        ind = reynolds >= 3e4
        self.dict_friction_factor[nodal]["turbulent"][ind] = (
            1.0 / (1.8 * np.log10(reynolds[ind]) / np.log10(10.0) - 1.64) ** 2 / 4.0
        )

    # End method lhc_poncet_correlation_bundle.

    def iter_cs_correlation_bundle(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.

        Raises:
            ValueError: [description]
        """
        raise ValueError("corr_peri value not defined!\n")
        ### LAURA mi dira' quanto. Da ottenere per mezzo debug \
        # dal Fortran quando si testa una simulazione con IFRICTION = 208
        corr_peri = ""  # ??
        # Evaluate turbulent friction factor.
        self.dict_friction_factor[nodal]["turbulent"] = (
            (0.047 + 19.5 / (reynolds / corr_peri) ** 0.857)
            / (self.inputs["VOID_FRACTION"] ** 0.742)
            / 4.0
            * corr_peri
        )

    # End method iter_cs_correlation_bundle.

    def dtt_correlation_bundle(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Evaluate jj and kk invocking method _eval_jj_and_kk.
        jj, kk = self._eval_jj_and_kk()
        # Evaluate turbulent friction factor.
        self.dict_friction_factor[nodal]["turbulent"] = (
            self.inputs["HYDIAMETER"] ** 2
            * self.inputs["VOID_FRACTION"]
            / 2.0
            / kk
            / reynolds
            + (self.inputs["HYDIAMETER"] * self.inputs["VOID_FRACTION"] ** 2)
            / 2.0
            * jj
            * (self.inputs["HYDIAMETER"] / self.inputs["VOID_FRACTION"] / np.sqrt(kk))
            ** 0.14
            / reynolds ** 0.14
        )

    # End method dtt_correlation_bundle.

    def demo_tf_hts_correlation_bundle(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Evaluate turbulent friction factor.
        # This is a boolean numpy array.
        ind = reynolds < 1.5e3
        self.dict_friction_factor[nodal]["turbulent"][ind] = (
            47.65 * reynolds[ind] ** -0.885 / 4.0
        )
        # This is a boolean numpy array.
        ind = (reynolds >= 1.5e3) & (reynolds <= 2e5)
        self.dict_friction_factor[nodal]["turbulent"][ind] = (
            1.093 * reynolds[ind] ** -0.338 / 4.0
        )
        # This is a boolean numpy array.
        ind = reynolds > 2e5
        self.dict_friction_factor[nodal]["turbulent"][ind] = 0.0377 / 4.0

    # End method demo_tf_hts_correlation_bundle.

    def hts_cl_correlation_bundle(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Evaluate turbulent friction factor: conversion to Fanning friction factor included.
        # This is a boolean numpy array.
        ind = reynolds < 1e3
        self.dict_friction_factor[nodal]["turbulent"][ind] = (
            194.002 * reynolds[ind] ** -0.52 / 4.0
        )  # (1-0.07) errorbar up.
        # This is a boolean numpy array.
        ind = (reynolds >= 1e3) & (reynolds < 2e3)
        self.dict_friction_factor[nodal]["turbulent"][ind] = (
            4.8563 + 4.87e-4 * reynolds[ind] / 4.0
        )  # (1-0.07) errorbar up.
        # This is a boolean numpy array.
        ind = reynolds >= 2e3
        self.dict_friction_factor[nodal]["turbulent"][ind] = (
            14.5148 * reynolds[ind] ** -0.12 / 4.0
        )  # (1-0.09) errorbar up.
        # N.B FL = FT in such a way to use max function.
        # forced to be 1 bto avoid problem from input.
        self.inputs["FRICTION_MULTIPLIER"] = 1.0

    # End method hts_cl_correlation_bundle.

    def _eval_total_friction_factor_L_lim_T(self, reynolds, L_lim_T=4e3, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            L_lim_T ([type], optional): [description]. Defaults to 4e3.
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Evaluate total friction factor by interpolation.
        # This is a boolean numpy array.
        ind = reynolds <= 2.0e3
        self.dict_friction_factor[nodal]["total"][ind] = self.dict_friction_factor[
            nodal
        ]["laminar"][ind]
        # This is a boolean numpy array.
        ind = (reynolds > 2.0e3) & (reynolds < L_lim_T)
        self.dict_friction_factor[nodal]["total"][ind] = (reynolds[ind] - 2000.0) / (
            L_lim_T - 2000.0
        ) * (
            self.dict_friction_factor[nodal]["turbulent"][ind]
            - self.dict_friction_factor[nodal]["laminar"][ind]
        ) + self.dict_friction_factor[
            nodal
        ][
            "laminar"
        ][
            ind
        ]
        # This is a boolean numpy array.
        ind = reynolds >= L_lim_T
        self.dict_friction_factor[nodal]["total"][ind] = self.dict_friction_factor[
            nodal
        ]["turbulent"][ind]

    # End method _eval_total_friction_factor_L_lim_T.

    def _eval_total_friction_factor_max(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Evaluate total friction factor.
        self.dict_friction_factor[nodal]["total"] = np.maximum(
            self.dict_friction_factor[nodal]["laminar"],
            self.dict_friction_factor[nodal]["turbulent"],
        )

    # End method _eval_total_friction_factor_max.

    def eval_friction_factor(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Initialize useful quantities.
        reynolds = self._quantities_initialization(reynolds, nodal)
        # Evaluate laminar friction factor according to the correlation defined in input.
        self.laminar_friction_factor_correlation(reynolds, nodal)
        # Evaluate turbulent friction factor according to the correlation defined in input.
        self.turbulent_friction_factor_correlation(reynolds, nodal)
        # Evaluate total friction factor.
        self.total_friction_factor_correlation(reynolds, nodal)
        # Multiply the total friction factor by the friction multiplier
        self.dict_friction_factor[nodal]["total"] = (
            self.dict_friction_factor[nodal]["total"]
            * self.inputs["FRICTION_MULTIPLIER"]
        )

    # End method eval_friction_factor.

    def _initialize_nusselt_and_htc(self, reynolds, nodal=True):
        """[summary]

        Args:
            reynolds ([type]): [description]
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Initialzie nusselt.
        self.dict_nusselt[nodal] = np.zeros(reynolds.shape)
        # Initialzie steady state heat transfer coefficient.
        self.dict_htc_steady[nodal] = np.zeros(reynolds.shape)

    # End method _initialize_htc.

    def dittus_boelter(self, dict_prop, nn=0.3, nodal=True):
        """[summary]

        Args:
            dict_prop ([type]): [description]
            nn (float, optional): [description]. Defaults to 0.3.
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Evaluate Nusselt number with the Dittus Boelter correlation.
        self.dict_nusselt[nodal] = (
            0.023 * dict_prop["Reynolds"] ** 0.8 * dict_prop["Prandtl"] ** nn
        )

    # End method dittus_boelter.

    def dittus_boelter_nusselt_lower_limit(self, dict_prop, nn=0.3, nodal=True):
        """[summary]

        Args:
            dict_prop ([type]): [description]
            nn (float, optional): [description]. Defaults to 0.3.
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Define dictionary with the lower limit value for the Nusselt number
        dict_nusselt_lower_limit = { 
            1: 8.235,
            119: self._eval_nusselt_lower_limit_119(),
            120: 2.181,
            121: 4.01,
        }
        # Evaluate Nusselt number with Dittus Boelter correlation.
        self.dittus_boelter(dict_prop, nn, nodal)
        # Correct the Nusselt number with the lower limit value.
        self.dict_nusselt[nodal][
            self.dict_nusselt[nodal]
            < dict_nusselt_lower_limit[self.inputs["Flag_htc_steady_corr"]]
        ] = dict_nusselt_lower_limit[self.inputs["Flag_htc_steady_corr"]]

    # End method dittus_boelter_nusselt_lower_limit.

    def _eval_nusselt_lower_limit_119(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        # Called if self.inputs["Flag_htc_steady_corr"] = 119
        alpha = self.inputs["HYDIAMETER"] / (self.inputs["HYDIAMETER"] + 2 * 1e-3)
        NuTemp = 7.541 * (
            1.0
            - 2.610 * alpha
            + 4.970 * alpha ** 2
            - 5.119 * alpha ** 3
            + 2.702 * alpha ** 4
            - 0.548 * alpha ** 5
        )
        Nuflux = 8.235 * (
            1.0
            - 10.6044 * alpha
            + 61.1755 * alpha ** 2
            - 155.1803 * alpha ** 3
            + 176.9203 * alpha ** 4
            - 72.9236 * alpha ** 5
        )
        # Evaluate lower limit Nusselt.
        return (NuTemp + Nuflux) / 2.0

    # End method _eval_nusselt_lower_limit_119

    def correlation_211(self, dict_prop, nn=3, nodal=True):
        """[summary]

        Args:
            dict_prop ([type]): [description]
            nn (int, optional): [description]. Defaults to 3.
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Boolean numpy array.
        ind = dict_prop["Reynolds"] < 1000.0
        self.dict_nusselt[nodal][ind] = 5.0969 * dict_prop["Reynolds"][ind] ** 0.10
        # Boolean numpy array.
        ind = (dict_prop["Reynolds"] >= 1000.0) & (dict_prop["Reynolds"] < 2000.0)
        self.dict_nusselt[nodal][ind] = 10.2029 - 3.3242e-5 * dict_prop["Reynolds"][ind]
        # Boolean numpy array.
        ind = dict_prop["Reynolds"] >= 2000.0
        self.dict_nusselt[nodal][ind] = 0.0395 * dict_prop["Reynolds"][ind] ** 0.73

    # End method correlation_211.

    def rectangular_duct_enea_hts_cicc_htc(self, dict_prop, nn=0.71, nodal=True):
        """Method that evaluates the heat transfer coefficient of the rectangular duct of the ENEA HTS CICC (given by prof L. Savoldy, no reference).

        Args:
            dict_prop ([type]): [description]
            nn (float, optional): [description]. Defaults to 0.71.
            nodal (bool, optional): [description]. Defaults to True.
        """
        self.dict_nusselt[nodal] = 0.42 * dict_prop["Reynolds"] ** nn
        # End method rectangular_duct_enea_hts_cicc_htc.

    def eval_steady_state_htc(self, dict_prop, nn=0.3, nodal=True):
        """[summary]

        Args:
            dict_prop ([type]): [description]
            nn (float, optional): [description]. Defaults to 0.3.
            nodal (bool, optional): [description]. Defaults to True.
        """
        # Initiaize useful quantities.
        self._initialize_nusselt_and_htc(dict_prop["Reynolds"], nodal)
        # Evaluate Nusselt with the sutable correlation.
        self.nusselt_correlation(dict_prop, nn, nodal)
        # Evaluate steady state heat transfer coefficient: Nu*ther_cond/L
        self.dict_htc_steady[nodal] = (
            self.dict_nusselt[nodal]
            * dict_prop["total_thermal_conductivity"]
            / self.inputs["HYDIAMETER"]
        )

    # End method eval_steady_state_htc.

    def do_nothing(self):
        """Method that does nothing"""
        pass

    # End method do_nothing
