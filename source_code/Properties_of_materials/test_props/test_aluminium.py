from ..temperature_dependent_materials import TemperatureDependentMaterial
# import properties_parent_class as ppc
import numpy as np

def testAluminiumThermalConductivity():
    # Tuples of temperature in K and corresponding thermal conductiity in W/m/K
    values = [
        (2.0, 6.6510),
        (4.0, 15.12),
        (6.0, 23.52),
        (8.0, 32.25),
        (10.0, 41.96),
        (20.0, 85.39),
        (40.0, 163.0),
        (60.0, 203.75),
        (80.0, 214.47),
        (100.0, 216.76),
        (200.0, 212.24),
        (300.0, 209.33),
        (1000.0, 209.33),
    ]
    al = TemperatureDependentMaterial().aluminium()
    tol = 1e-12
    for _, val in enumerate(values):
        assert abs(val[1] - al.thermal_conductivity(val[0])) <= tol


def testAluminiumTermalConductivityBelowLB():

    temperature = 1.0
    al = TemperatureDependentMaterial().aluminium()
    assert 6.6510 == al.thermal_conductivity(temperature)


def testAluminiumTermalConductivityAboveUB():

    temperature = 1500.0
    al = TemperatureDependentMaterial().aluminium()
    assert 209.3300 == al.thermal_conductivity(temperature)

def testAluminiumIsobaricSpecificHeat():
    # Tuples of temperature in K and corresponding thermal conductiity in W/m/K
    values = [
        (2.0, 0.10384646),
        (4.0, 0.26538452),
        (6.0, 0.51508074),
        (8.0, 0.88340168),
        (10.0, 1.4008139),
        (20.0, 8.886111950000004),
        (40.0, 77.047596550000020),
        (60.0, 2.088806270176337e+02),
        (80.0, 3.584085604675565e+02),
        (100.0, 4.913183769264615e+02),
        (200.0, 7.881989008530381e+02),
        (300.0, 8.822207857778576e+02),
        (1000.0, 9.704934656408732e+02),
    ]
    al = TemperatureDependentMaterial().aluminium()
    tol = 1e-12
    for ii, val in enumerate(values):
        assert abs(val[1] - al.isobaric_specific_heat(val[0])) <= tol, f"{ii}; {val}\n"

def testAluminiumIsobaricSpecificHeatBelowLB():

    temperature = 1.0
    al = TemperatureDependentMaterial().aluminium()
    assert 0.10384646 == al.isobaric_specific_heat(temperature)

def testAluminiumIsobaricSpecificHeatAboveUB():

    temperature = 1500.0
    al = TemperatureDependentMaterial().aluminium()
    assert 9.704934656408732e+02 == al.isobaric_specific_heat(temperature)


def testAluminiumElectricalResistivity():
    # Tuples of temperature in K and corresponding thermal conductiity in W/m/K
    values = [
        (2.0, 5.76e-09),
        (4.0, 5.759997441216000e-09),
        (6.0, 5.760002558776000e-09),
        (8.0, 5.759992323680000e-09),
        (10.0, 5.760028146420000e-09),
        (20.0, 5.759215479800000e-09),
        (40.0, 5.771738084000000e-09),
        (60.0, 6.695839936000000e-09),
        (80.0, 8.868304228000000e-09),
        (100.0, 1.121755103120000e-08),
        (200.0, 2.236174823180000e-08),
        (300.0, 3.511866612080000e-08),
        (1000.0, 1.161198020995000e-07),
    ]
    al = TemperatureDependentMaterial().aluminium()
    tol = 1e-12
    for ii, val in enumerate(values):
        assert abs(val[1] - al.electrical_resistivity(val[0])) <= tol, f"{ii}; {val}\n"

def testAluminiumElectricalResistivityBelowLB():

    temperature = 1.0
    al = TemperatureDependentMaterial().aluminium()
    assert 5.76e-09 == al.electrical_resistivity(temperature)

def testAluminiumElectricalResistivityAboveUB():

    temperature = 1500.0
    al = TemperatureDependentMaterial().aluminium()
    assert 1.161198020995000e-07 == al.electrical_resistivity(temperature)

def testAluminiumDensity():
    al = TemperatureDependentMaterial().aluminium()
    
    assert 2700 == al.density()