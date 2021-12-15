# Import python libraries and other functions
import numpy as np

# Function CONDJK2LB starts here
def theramal_conductivity_jk2lb(TT):
    """
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #% **********************************
    #% ***   Roberto Bonifetto    ***
    #% ***   Politecnico di Torino  ***
    #% ***  Dipartimento Energia  ***
    #% ***    May 07, 2015    ***
    #% **********************************
    #% Array-smart function for JK2LB specific heat
    #% TTT = temperature [K]
    #% COND = thermal conductivity [W/m/K]
    #% -------------------------------------------------------------------------
    #% Properties from:
    #% J. Lu, R.P. Walsh, K. Han
    #% Low temperature physical properties of a high Mn austenitic steel JK2LB
    #% Cryogenics 49, pp. 133-137, 2009.
    #%
    #% Data extracted from the graphs in the paper
    #%
    #% Source: Cryosoft
    #% Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    # Correction after test: D.Placido PoliTo 04/04/2020
    # Tested against corresponding MatLab function in temperature range [4,400] K
    #% -------------------------------------------------------------------------
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    """

    TT = np.array(TT)
    COND = np.zeros(TT.shape)  # variable initialization

    T0 = 20.0
    T1 = 300.0
    Tmax = 1000.0
    #% Polinomial of degree 3 for T <= T0
    #% Polinomial of degree 6 for T > T0
    #% Coefficients( with 95 % confidence bounds):
    p1a = -3.495e-05  # (-4.704e-05, -2.286e-05)
    p2a = 0.001794  # (0.001323, 0.002265)
    p3a = 0.06738  # (0.06205, 0.07272)
    p4a = -0.03533  # (-0.05162, -0.01904)

    p1b = 2.174e-14  # % (-1.541e-15, 4.503e-14)
    p2b = -1.451e-11  # % (-3.663e-11, 7.621e-12)
    p3b = -2.417e-10  #% (-8.404e-09, 7.92e-09)
    p4b = 2.426e-06  # % (9.55e-07, 3.898e-06)
    p5b = -0.0007596  # % (-0.0008933, -0.0006259)
    p6b = 0.1175  # % (0.1118, 0.1231)
    p7b = -0.2904  # % (-0.3743, -0.2065)

    cond250 = 9.7147
    cond300 = 10.7446
    dy = -0.0243

    TT[TT > Tmax] = Tmax

    intervals = [(TT <= T0), (TT > T0) & (TT <= T1), (TT > T1)]

    behavior = [
        lambda TT: p1a * TT ** 3 + p2a * TT ** 2 + p3a * TT + p4a,
        lambda TT: p1b * TT ** 6
        + p2b * TT ** 5
        + p3b * TT ** 4
        + p4b * TT ** 3
        + p5b * TT ** 2
        + p6b * TT ** 1
        + p7b
        + dy,
        lambda TT: cond250 + dy + (cond300 - cond250) / (300 - 250) * (TT - 250),
    ]

    COND = np.piecewise(TT, intervals, behavior)

    return COND  # end of the function


# Function CPJK2LB starts here
def isobaric_specific_heat_jk2lb(TTT):
    """
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % **********************************
    # % ***     Roberto Bonifetto      ***
    # % ***   Politecnico di Torino    ***
    # % ***    Dipartimento Energia    ***
    # % ***        May 07, 2015        ***
    # % **********************************
    # % Array-smart function for JK2LB specific heat
    # % TTT = temperature [K]
    # % CP = specific heat [J/kg/K]
    # % -------------------------------------------------------------------------
    # % Properties from:
    # % J. Lu, R.P. Walsh, K. Han
    # % Low temperature physical properties of a high Mn austenitic steel JK2LB
    # % Cryogenics 49, pp. 133-137, 2009.
    # %
    # % Data extracted from the graphs in the paper
    # %
    # % Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    # Correction after test: D.Placido PoliTo 02/04/2020
    # Tested against corresponding MatLab function in temperature range [4,400] K
    # % -------------------------------------------------------------------------
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    """

    TTT = np.array(TTT)
    T0 = 50.0
    T1 = 300.0
    Tmax = 1000.0

    CP = np.zeros(TTT.shape)  # variable initialization

    # Polinomial of degree 6 for T <= T0
    # Polinomial of degree 3 for T > T0
    # Coefficients (with 95% confidence bounds):

    p1a = -4.18e-09  #% (-1.016e-08, 1.796e-09)
    p2a = 3.04e-07  #% (-6.271e-07, 1.235e-06)
    p3a = -1.039e-05  #% (-6.585e-05, 4.508e-05)
    p4a = 0.001128  #% (-0.000453, 0.002709)
    p5a = -0.01538  #% (-0.03742, 0.006657)
    p6a = 0.3118  #% (0.1757, 0.4478)
    p7a = -0.1974  #% (-0.4711, 0.07628)
    p1b = 2.092e-05  #% (1.806e-05, 2.378e-05)
    p2b = -0.01945  #% (-0.02092, -0.01797)
    p3b = 6.149  #% (5.925, 6.373)
    p4b = -179.3  #% (-189, -169.5)

    cp299 = 479.7580
    cp300 = 479.8868
    dy = 0.5230

    TTT[TTT > Tmax] = Tmax

    intervals = [(TTT <= T0), (TTT > T0) & (TTT <= T1), (TTT > T1)]

    behavior = [
        lambda TTT: p1a * TTT ** 6
        + p2a * TTT ** 5
        + p3a * TTT ** 4
        + p4a * TTT ** 3
        + p5a * TTT ** 2
        + p6a * TTT
        + p7a,
        lambda TTT: p1b * TTT ** 3 + p2b * TTT ** 2 + p3b * TTT + p4b + dy,
        lambda TTT: cp299 + dy + (cp300 - cp299) / (300 - 299) * (TTT - 299),
    ]

    CP = np.piecewise(TTT, intervals, behavior)

    return CP  # end of the function


# Function RHOEJK2LB starts here
def electrical_resistivity_jk2lb(TTT):
    """
    ###########################################################################
    # **********************************
    # ***     Roberto Bonifetto      ***
    # ***   Politecnico di Torino    ***
    # ***    Dipartimento Energia    ***
    # ***        May 07, 2015        ***
    # **********************************
    # Array-smart function for JK2LB electrical resistivity
    # TTT = temperature [K]
    # RHOE = electrical resistivity [Ohm*m]
    # -------------------------------------------------------------------------
    # Properties from:
    # J. Lu, R.P. Walsh, K. Han
    # Low temperature physical properties of a high Mn austenitic steel JK2LB
    # Cryogenics 49, pp. 133-137, 2009.
    #
    # Data extracted from the graphs in the paper
    #
    # -------------------------------------------------------------------------
    #
    # Translation from MatLab to Python: D.Placido PoliTo 03/2020
    # Correction after test: D.Placido PoliTo 02/04/2020
    # Tested against corresponding MatLab function in temperature range [4,400] K
    ###########################################################################
    """

    TTT = np.array(TTT)
    RHOE = np.array(TTT.shape)

    T1 = 300.0
    Tmax = 1000.0
    # Polinomial of degree 6
    # Coefficients (with 95% confidence bounds):
    p1a = -2.506e-21  # (-3.319e-21, -1.693e-21)
    p2a = 2.214e-18  # (1.489e-18, 2.939e-18)
    p3a = -6.662e-16  # (-9.124e-16, -4.201e-16)
    p4a = 7.45e-14  # (3.48e-14, 1.142e-13)
    p5a = -5.208e-13  # (-3.603e-12, 2.562e-12)
    p6a = -1.466e-10  # (-2.487e-10, -4.46e-11)
    p7a = 9.109e-07  # (9.098e-07, 9.119e-07)

    cond280 = 9.7261e-07
    cond300 = 9.8883e-07

    TTT[TTT > Tmax] = Tmax

    intervals = [(TTT <= T1), (TTT > T1)]

    behavior = [
        lambda TTT: p1a * TTT ** 6
        + p2a * TTT ** 5
        + p3a * TTT ** 4
        + p4a * TTT ** 3
        + p5a * TTT ** 2
        + p6a * TTT ** 1
        + p7a,
        lambda TTT: cond280 + (cond300 - cond280) / (300.0 - 280.0) * (TTT - 280.0),
    ]

    RHOE = np.piecewise(TTT, intervals, behavior)

    return RHOE  # end of the function
