from ..temperature_dependent_materials import TemperatureDependentMaterial
import os
import pandas as pd


def testDensity():
    rho = {0: 1.380e3, 1: 1.420e3}
    for source in range(2):
        ka = TemperatureDependentMaterial().kapton(source=source)
        assert rho[ka.source] == ka.density()
        assert "Density" == ka.property_name

def testThermalConductivity():

    file_path = os.path.relpath(
        "C:/Users\Daniele.Placido\OneDrive - Politecnico di Torino\MAHTEP\OPENSC2_work\Material_properties_Matlab\Pl\properties.xlsx"
    )
    df = pd.read_excel(file_path, sheet_name="kk")

    ka = TemperatureDependentMaterial().kapton()

    tol = 1e-12
    df["thermal_conductivity_py"] = ka.thermal_conductivity(df["temperature"])
    assert "Thermal conductivity" == ka.property_name

    df["error"] = abs(df["thermal_conductivity"] - df["thermal_conductivity_py"])
    assert df["error"].min() <= tol
    assert df["error"].max() <= tol