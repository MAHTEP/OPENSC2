from ..temperature_dependent_materials import TemperatureDependentMaterial
import os
import pandas as pd


def testDensity():
    inconel = TemperatureDependentMaterial().inconel()

    assert 8.170e3 == inconel.density()
    assert "Density" == inconel.property_name

def testThermalConductivity():

    file_path = os.path.relpath(
        "C:/Users\Daniele.Placido\OneDrive - Politecnico di Torino\MAHTEP\OPENSC2_work\Material_properties_Matlab/NE\properties.xlsx"
    )
    df = pd.read_excel(file_path, sheet_name="kk")

    inconel = TemperatureDependentMaterial().inconel()

    tol = 1e-12
    df["thermal_conductivity_py"] = inconel.thermal_conductivity(df["temperature"])
    assert "Thermal conductivity" == inconel.property_name

    df["error"] = abs(df["thermal_conductivity"] - df["thermal_conductivity_py"])
    assert df["error"].min() <= tol
    assert df["error"].max() <= tol


def testIsobaricSpecificHeat():

    file_path = os.path.relpath(
        "C:/Users\Daniele.Placido\OneDrive - Politecnico di Torino\MAHTEP\OPENSC2_work\Material_properties_Matlab/NE\properties.xlsx"
    )
    df = pd.read_excel(file_path, sheet_name="cp")

    inconel = TemperatureDependentMaterial().inconel()

    tol = 1e-12
    df["isobaric_specific_heat_py"] = inconel.isobaric_specific_heat(df["temperature"])
    assert "Isobaric specific heat" == inconel.property_name

    df["error"] = abs(df["isobaric_specific_heat"] - df["isobaric_specific_heat_py"])
    assert df["error"].min() <= tol
    assert df["error"].max() <= tol

def testElectricalResistivity():

    file_path = os.path.relpath(
        "C:/Users\Daniele.Placido\OneDrive - Politecnico di Torino\MAHTEP\OPENSC2_work\Material_properties_Matlab/NE\properties.xlsx"
    )
    df = pd.read_excel(file_path, sheet_name="rhoe")

    inconel = TemperatureDependentMaterial().inconel()

    tol = 1e-12
    df["electrical_resistivity_py"] = inconel.electrical_resistivity(df["temperature"])
    assert "Electrical resistivity" == inconel.property_name

    df["error"] = abs(df["electrical_resistivity"] - df["electrical_resistivity_py"])
    assert df["error"].min() <= tol
    assert df["error"].max() <= tol