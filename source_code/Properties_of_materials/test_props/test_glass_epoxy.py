from ..temperature_dependent_materials import TemperatureDependentMaterial

def testSwicthFitEquaion():
    ge = TemperatureDependentMaterial().glass_epoxy(source=1)
    if ge.source == 0:
        assert ge._fit_cryosoft == ge._fit_equation
    elif ge.source == 1:
        assert ge._fit_nist == ge._fit_equation


def testGlassEpoxyThermalConductivity():
    tol = 1e-12
    ge = TemperatureDependentMaterial().glass_epoxy(source=1)
    assert "GlassEpoxy" == ge.material
    if ge.source == 0:
        # Tuples of temperature in K and corresponding thermal conductiity in W/m/K
        values = [
            (1.0, 0.002851314481037),
            (5.0, 0.069136714005337),
            (10.0, 0.106378073518238),
            (20.0, 0.151760633692667),
            (40.0, 0.232802052454053),
            (60.0, 0.305905559247382),
            (80.0, 0.367758133998533),
            (100.0, 0.418883095028124),
            (200.0, 0.569617247617708),
            (300.0, 0.637061363339246),
            (500.0, 0.697221925886015),
            (750.0, 0.734439058410778),
            (1000.0, 0.762085220457567),
        ]
    elif ge.source == 1:
        # Tuples of temperature in K and corresponding thermal conductiity in W/m/K
        values = [
            (4.0, 0.073219988906989),
            (5.0, 0.085428700327756),
            (10.0, 0.136084919271237),
            (20.0, 0.198172573902598),
            (40.0, 0.274015451829511),
            (60.0, 0.336644433101416),
            (80.0, 0.394183268207437),
            (100.0, 0.447723613274968),
            (150.0, 0.567734478188575),
            (200.0, 0.674129296922664),
            (250.0, 0.772083452600512),
            (300.0, 0.863637841735930),
        ]

    for ii, val in enumerate(values):
        assert abs(val[1] - ge.thermal_conductivity(val[0])) <= tol, f"{ii = }, {val[0] = }\n"

    assert "Thermal conductivity" == ge.property_name


def testGlassEpoxyIsobaricSpecofocHeat():
    ge = TemperatureDependentMaterial().glass_epoxy(source=1)
    assert "GlassEpoxy" == ge.material
    tol = 1e-12
    if ge.source == 0:
        values = [
            (1.0, 0.005678744796811),
            (5.0, 0.751328021202572),
            (10.0, 7.160971275895215),
            (20.0, 40.220237920295830),
            (40.0, 1.373565231794619e02),
            (60.0, 2.335385198627686e02),
            (80.0, 3.197573864810314e02),
            (100.0, 3.975104291616320e02),
            (200.0, 7.196152121942605e02),
            (300.0, 9.906750853160461e02),
            (500.0, 1.456864530348298e03),
            (750.0, 1.947609763158504e03),
            (1000.0, 2.369980896025295e03),
        ]
    elif ge.source == 1:
        tol = 1e-10
        # Tuples of temperature in K and corresponding thermal conductiity in W/m/K
        values = [
            (4.0, 2.015888566386214),
            (5.0, 3.571339512312817),
            (10.0, 15.355925625514057),
            (20.0, 47.171387832312426),
            (40.0, 1.153594130992956e+02),
            (60.0, 1.825304286854036e+02),
            (80.0, 2.494067294193294e+02),
            (100.0, 3.168604247190897e+02),
            (150.0, 4.892850679525539e+02),
            (200.0, 6.644291629061815e+02),
            (250.0, 8.363324610546729e+02),
            (300.0, 9.987427171612554e+02),
        ]

    for ii, val in enumerate(values):
        assert abs(val[1] - ge.isobaric_specific_heat(val[0])) <= tol, f"{ii = }, {val[0] = }\n"

    assert "Isobaric specific heat" == ge.property_name


def testGlassEpoxyDensity():
    ge = TemperatureDependentMaterial().glass_epoxy(source=0)
    assert "GlassEpoxy" == ge.material
    if ge.source == 0:
        assert 2e3 == ge.density()
    elif ge.source == 1:
        assert 1.8e3 == ge.density()
    
    assert "Density" == ge.property_name