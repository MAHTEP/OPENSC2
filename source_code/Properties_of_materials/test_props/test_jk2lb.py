from ..temperature_dependent_materials import TemperatureDependentMaterial
import os
import pandas as pd


def testDensity():
    jk2lb = TemperatureDependentMaterial().high_mn_austenitic_steinless_steel_jk2lb()

    assert 7.860e3 == jk2lb.density()
    assert "Density" == jk2lb.property_name


def testThermalConductivity():

    file_path = os.path.relpath(
        "C:/Users\Daniele.Placido\OneDrive - Politecnico di Torino\MAHTEP\OPENSC2_work\Material_properties_Matlab\JK2LB\properties.xlsx"
    )
    df = pd.read_excel(file_path, usecols=["temperature", "thermal_conductivity"])

    jk2lb = TemperatureDependentMaterial().high_mn_austenitic_steinless_steel_jk2lb()

    tol = 1e-12
    df["thermal_conductivity_py"] = jk2lb.thermal_conductivity(df["temperature"])
    assert "Thermal conductivity" == jk2lb.property_name

    df["error"] = abs(df["thermal_conductivity"] - df["thermal_conductivity_py"])
    assert df["error"].min() <= tol
    assert df["error"].max() <= tol


def testIsobaricSpecifiHeat():

    file_path = os.path.relpath(
        "C:/Users\Daniele.Placido\OneDrive - Politecnico di Torino\MAHTEP\OPENSC2_work\Material_properties_Matlab\JK2LB\properties.xlsx"
    )
    df = pd.read_excel(file_path, usecols=["temperature", "specific_heat"])

    jk2lb = TemperatureDependentMaterial().high_mn_austenitic_steinless_steel_jk2lb()

    tol = 1e-12
    df["specific_heat_py"] = jk2lb.isobaric_specific_heat(df["temperature"])
    assert "Isobaric specific heat" == jk2lb.property_name

    df["error"] = abs(df["specific_heat"] - df["specific_heat_py"])
    assert df["error"].min() <= tol
    assert df["error"].max() <= tol


def testElectricalResistivity():

    file_path = os.path.relpath(
        "C:/Users\Daniele.Placido\OneDrive - Politecnico di Torino\MAHTEP\OPENSC2_work\Material_properties_Matlab\JK2LB\properties.xlsx"
    )
    df = pd.read_excel(file_path, usecols=["temperature", "electrical_resistivity"])

    jk2lb = TemperatureDependentMaterial().high_mn_austenitic_steinless_steel_jk2lb()

    tol = 1e-12
    df["electrical_resistivity_py"] = jk2lb.electrical_resistivity(df["temperature"])
    assert "Electrical resistivity" == jk2lb.property_name

    df["error"] = abs(df["electrical_resistivity"] - df["electrical_resistivity_py"])
    assert df["error"].min() <= tol
    assert df["error"].max() <= tol