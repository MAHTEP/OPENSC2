from ..temperature_dependent_materials import TemperatureDependentMaterial
import os
import pandas as pd


def testDensity():
    for source in range(2):
        incoloy = TemperatureDependentMaterial().incoloy(source=source)
        assert 8.17e3 == incoloy.density()
        assert "Density" == incoloy.property_name


def testThermalConductivity():
    for source in range(2):
        incoloy = TemperatureDependentMaterial().incoloy(source=source)

        file_path = os.path.relpath(
            f"C:/Users\Daniele.Placido\OneDrive - Politecnico di Torino\MAHTEP\OPENSC2_work\Material_properties_Matlab\IN\properties_{incoloy.source_name}.xlsx"
        )
        df = pd.read_excel(
            file_path, sheet_name="kk", usecols=["temperature", "thermal_conductivity"]
        )

        tol = 1e-12
        df["thermal_conductivity_py"] = incoloy.thermal_conductivity(df["temperature"])
        assert "Thermal conductivity" == incoloy.property_name

        df["error"] = abs(df["thermal_conductivity"] - df["thermal_conductivity_py"])
        assert df["error"].min() <= tol, f"{incoloy.source_name}, {df['error'].min() = }"
        assert df["error"].max() <= tol, f"{incoloy.source_name}, {df['error'].max() = }"


def testIsobaricSpecificHeat():
    for source in range(2):
        incoloy = TemperatureDependentMaterial().incoloy(source=source)

        file_path = os.path.relpath(
            f"C:/Users\Daniele.Placido\OneDrive - Politecnico di Torino\MAHTEP\OPENSC2_work\Material_properties_Matlab\IN\properties_{incoloy.source_name}.xlsx"
        )
        df = pd.read_excel(
            file_path, sheet_name="cp", usecols=["temperature", "isobaric_specific_heat"]
        )

        tol = 1e-12
        df["isobaric_specific_heat_py"] = incoloy.isobaric_specific_heat(df["temperature"])
        assert "Isobaric specific heat" == incoloy.property_name

        df["error"] = abs(df["isobaric_specific_heat"] - df["isobaric_specific_heat_py"])
        assert df["error"].min() <= tol, f"{incoloy.source_name}, {df['error'].min() = }"
        assert df["error"].max() <= tol, f"{incoloy.source_name}, {df['error'].max() = }"

def testElectricalResistivity():
    for source in range(2):
        incoloy = TemperatureDependentMaterial().incoloy(source=source)

        file_path = os.path.relpath(
            f"C:/Users\Daniele.Placido\OneDrive - Politecnico di Torino\MAHTEP\OPENSC2_work\Material_properties_Matlab\IN\properties_{incoloy.source_name}.xlsx"
        )
        df = pd.read_excel(
            file_path, sheet_name="rhoe", usecols=["temperature", "electrical_resistivity"]
        )

        tol = 1e-12
        df["electrical_resistivity_py"] = incoloy.electrical_resistivity(df["temperature"])
        assert "Electrical resistivity" == incoloy.property_name

        df["error"] = abs(df["electrical_resistivity"] - df["electrical_resistivity_py"])
        assert df["error"].min() <= tol, f"{incoloy.source_name}, {df['error'].min() = }"
        assert df["error"].max() <= tol, f"{incoloy.source_name}, {df['error'].max() = }"


# def testAllButDensity():

#     properties = dict(
#         electrical_resistivity = ("rhoe", "Electrical resistivity"),
#         isobaric_specific_heat = ("cp", "Isobaric specific heat"),
#         thermal_conductivity = ("kk", "Thermal conductivity"),
#     )

#     for source in range(2):
#         incoloy = TemperatureDependentMaterial().incoloy(source=source)

#         for key, value in properties.items():

#             file_path = os.path.relpath(
#                 f"C:/Users\Daniele.Placido\OneDrive - Politecnico di Torino\MAHTEP\OPENSC2_work\Material_properties_Matlab\IN\properties_{incoloy.source_name}.xlsx"
#             )
#             df = pd.read_excel(
#                 file_path, sheet_name=value[0], usecols=["temperature", key]
#             )

#             tol = 1e-12
#             df[f"{key}_py"] = getattr(incoloy, key)(df["temperature"])
#             assert value[1] == incoloy.property_name

#             df["error"] = abs(df[key] - df[f"{key}_py"])
#             assert df["error"].min() <= tol, f"{incoloy.source_name}, {df['error'].min() = }"
#             assert df["error"].max() <= tol, f"{incoloy.source_name}, {df['error'].max() = }"