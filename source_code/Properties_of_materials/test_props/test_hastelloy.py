from ..temperature_dependent_materials import TemperatureDependentMaterial
import os
import pandas as pd


def testThermalConductivity():

    file_path = os.path.relpath(
        "C:/Users\Daniele.Placido\OneDrive - Politecnico di Torino\MAHTEP\OPENSC2_work\Material_properties_Matlab\HC276\properties.xlsx"
    )
    df = pd.read_excel(
        file_path, sheet_name="k_cp", usecols=["temperature", "thermal_conductivity"]
    )

    ha = TemperatureDependentMaterial().hastelloy()

    tol = 1e-12
    df["thermal_conductivity_py"] = ha.thermal_conductivity(df["temperature"])
    assert "Thermal conductivity" == ha.property_name

    df["error"] = abs(df["thermal_conductivity"] - df["thermal_conductivity_py"])
    assert df["error"].min() <= tol
    assert df["error"].max() <= tol


def testIsobaricSpecificHeat():
    # Tuples of temperature in K and corresponding thermal conductiity in W/m/K

    file_path = os.path.relpath(
        "C:/Users\Daniele.Placido\OneDrive - Politecnico di Torino\MAHTEP\OPENSC2_work\Material_properties_Matlab\HC276\properties.xlsx"
    )
    df = pd.read_excel(
        file_path, sheet_name="k_cp", usecols=["temperature", "specific_heat"]
    )

    ha = TemperatureDependentMaterial().hastelloy()

    tol = 1e-12
    df["specific_heat_py"] = ha.isobaric_specific_heat(df["temperature"])
    assert "Isobaric specific heat" == ha.property_name

    df["error"] = abs(df["specific_heat"] - df["specific_heat_py"])
    assert df["error"].min() <= tol
    assert df["error"].max() <= tol


def testElectricalResistivity():
    file_path = os.path.relpath(
        "C:/Users\Daniele.Placido\OneDrive - Politecnico di Torino\MAHTEP\OPENSC2_work\Material_properties_Matlab\HC276\properties.xlsx"
    )
    df = pd.read_excel(
        file_path, sheet_name="rhoe", usecols=["temperature", "electrical_resistivity"]
    )

    ha = TemperatureDependentMaterial().hastelloy()

    tol = 1e-12
    df["electrical_resistivity_py"] = ha.electrical_resistivity(df["temperature"])
    assert "Electrical resistivity" == ha.property_name

    df["error"] = abs(df["electrical_resistivity"] - df["electrical_resistivity_py"])
    assert df["error"].min() <= tol
    assert df["error"].max() <= tol


def testDensity():
    ha = TemperatureDependentMaterial().hastelloy()

    assert 8.89e3 == ha.density()
    assert "Density" == ha.property_name