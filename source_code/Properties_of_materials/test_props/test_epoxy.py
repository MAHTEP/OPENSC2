from ..temperature_dependent_materials import TemperatureDependentMaterial

def testEpoxyThermalConductivity():
    # Tuples of temperature in K and corresponding thermal conductiity in W/m/K
    values = [
        (1.0, 0.016214163102889),
        (5.0, 0.051878541100108),
        (10.0, 0.058521473865065),
        (20.0, 0.067404069784445),
        (40.0, 0.089881996122016),
        (60.0, 0.111451089311360),
        (80.0, 0.130064206488107),
        (100.0, 0.145845892207336),
        (200.0, 0.195596593609067),
        (300.0, 0.218280561900595),
        (450.0, 0.230727494849419),
        (600.0, 0.233582983276373),
        (833.0, 0.237676907406309),
    ]
    ep = TemperatureDependentMaterial().epoxy()

    tol = 1e-12
    for _, val in enumerate(values):
        assert abs(val[1] - ep.thermal_conductivity(val[0])) <= tol

    assert "Thermal conductivity" == ep.property_name


def testEpoxyIsobaricSpecofocHeat():
    ep = TemperatureDependentMaterial().epoxy()
    values = [
        (1.0, 0.018025760356040),
        (5.0, 1.679015730147560),
        (10.0, 15.449315580576910),
        (20.0, 78.367233426112350),
        (40.0, 2.518766393355338e02),
        (60.0, 4.212494928273756e02),
        (80.0, 5.745040738947638e02),
        (100.0, 7.146337257319477e02),
        (200.0, 1.289027542151594e03),
        (300.0, 1.713966357707212e03),
        (450.0, 2.163543007837275e03),
        (600.0, 2.471335211731798e03),
        (833.0, 2.790755712728986e03),
    ]
    ep = TemperatureDependentMaterial().epoxy()

    tol = 1e-12
    for _, val in enumerate(values):
        assert abs(val[1] - ep.isobaric_specific_heat(val[0])) <= tol

    assert "Isobaric specific heat" == ep.property_name


def testEpoxyDensity():
    ep = TemperatureDependentMaterial().epoxy()

    assert 1.15e3 == ep.density()