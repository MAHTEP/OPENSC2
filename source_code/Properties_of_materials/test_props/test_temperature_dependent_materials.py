from ..temperature_dependent_materials import TemperatureDependentMaterial
import numpy as np


def testTempertaureDependentMaterial():
    prop = TemperatureDependentMaterial()
    assert "Prop" == prop.property_name


def testMaterialInstance():
    ranges = dict(
        aluminium=dict(
            electrical_resistivity=(2.0, 1000.0),  # K
            isobaric_specific_heat=(2.0, 1000.0),  # K
            thermal_conductivity=(2.0, 1000.0),  # K
        ),
        epoxy=dict(
            isobaric_specific_heat=(1.0, 833.0),  # K
            thermal_conductivity=(1.0, 833.0),  # K
        ),
        glass_epoxy=dict(
            isobaric_specific_heat=dict(
                Cryosoft=(1.0, 1000.0),  # K
                NIST=(4.0, 300.0),  # K
            ),
            thermal_conductivity=dict(
                Cryosoft=(1.0, 1000.0),  # K
                NIST=(4.0, 300.0),  # K),  # K
            ),
        ),
        hastelloy=dict(
            electrical_resistivity=(4.0, 300.0),  # K
            isobaric_specific_heat=(1.0, 350.0),  # K
            thermal_conductivity=(1.0, 350.0),  # K
        ),
        high_mn_austenitic_steinless_steel_jk2lb=dict(
            electrical_resistivity=(1.0, 1000.0),  # K
            isobaric_specific_heat=(1.0, 1000.0),  # K
            thermal_conductivity=(1.0, 1000.0),  # K
        ),
        incoloy=dict(
            thermal_conductivity=dict(Cryosoft=(1.0, 1000.0), KSTAR=(4.0, 300.0)),  # K
            isobaric_specific_heat=dict(
                Cryosoft=(4.27, 1423.0), KSTAR=(4.2, 300.0)
            ),  # K
            electrical_resistivity=dict(Cryosoft=(1.0, 300.0), KSTAR=(4.0, 300.0)),  # K
        ),
        inconel=dict(
            electrical_resistivity=(1.0, 300.0),  # K
            isobaric_specific_heat=(1, 1000.0),  # K
            thermal_conductivity=(1.0, 500.0),  # K
        ),
        kapton=dict(
            isobaric_specific_heat=dict(Cryosoft=(1, 300.0), NIST=(4.0, 300.0)),  # K
            thermal_conductivity=dict(Cryosoft=(1.0, 833.0), NIST=(4.0, 300.0)),  # K
        ),
    )
    names = dict(
        aluminium="Aluminium",
        epoxy="Epoxy",
        glass_epoxy="GlassEpoxy",
        hastelloy="Hastelloy",
        high_mn_austenitic_steinless_steel_jk2lb="HighMnAusteniticSteinlessSteelJK2LB",
        incoloy="Incoloy",
        inconel="Inconel",
        kapton="Kapton",
    )

    # method_list = [
    #     attribute
    #     for attribute in dir(TemperatureDependentMaterial)
    #     if callable(getattr(TemperatureDependentMaterial, attribute))
    #     and attribute.startswith("_") == False
    # ]

    for method, value in ranges.items():
        if (
            method != "glass_epoxy"
            and method != "incoloy"
            and method != "kapton"
        ):
            material = getattr(TemperatureDependentMaterial(), method)()
            assert names[method] == material.material
            for prop, temps in value.items():
                for ii, tt in enumerate(temps):
                    assert tt == material.temperature_ranges[prop][ii]
                # End for
            # End for
        elif method == "glass_epoxy" or method == "kapton":
            source_dict = {0: "Cryosoft", 1: "NIST"}
            for source in range(2):
                # source = 0: Cryosoft
                # sourcre = 1: NIST
                material = getattr(TemperatureDependentMaterial(), method)(
                    source=source
                )
                assert names[method] == material.material
                assert source_dict[source] == material.source_name

                for prop, temps in value.items():
                    for ii, tt in enumerate(temps[material.source_name]):
                        assert tt == material.temperature_ranges[prop][ii]
                    # End for
                # End for
            # End for

        elif method == "incoloy":

            source_dict = {0: "Cryosoft", 1: "KSTAR"}
            for source in range(2):
                # source = 0: Cryosoft
                # sourcre = 1: KSTAR
                material = getattr(TemperatureDependentMaterial(), method)(
                    source=source
                )
                assert names[method] == material.material
                assert source_dict[source] == material.source_name

                for prop, temps in value.items():
                    for ii, tt in enumerate(temps[material.source_name]):
                        assert tt == material.temperature_ranges[prop][ii]
                    # End for
                # End for
            # End for


def testConversionToNumpyArray():
    temperature = 10.0
    method_list = [
        attribute
        for attribute in dir(TemperatureDependentMaterial)
        if callable(getattr(TemperatureDependentMaterial, attribute))
        and attribute.startswith("_") == False
    ]
    for method in method_list:
        material = getattr(TemperatureDependentMaterial(), method)()
        assert True == isinstance(material._convert_to_nparray(temperature), np.ndarray)


def testTemperatureLowerBoundary():
    method_list = [
        attribute
        for attribute in dir(TemperatureDependentMaterial)
        if callable(getattr(TemperatureDependentMaterial, attribute))
        and attribute.startswith("_") == False
    ]
    for method in method_list:
        if (
            method != "glass_epoxy"
            and method != "incoloy"
            and method != "kapton"
        ):
            material = getattr(TemperatureDependentMaterial(), method)()
            for key, val in material.temperature_ranges.items():
                temperature = 0.5  # K
                temperature = material._convert_to_nparray(temperature)
                assert val[0] == material._fix_temperature_lb(
                    temperature, key
                ), f"{material.__class__.__name__}"

        elif method == "glass_epoxy" or method == "incoloy" or method == "kapton":
            for source in range(2):
                material = getattr(TemperatureDependentMaterial(), method)(
                    source=source
                )
                for key, val in material.temperature_ranges.items():
                    temperature = 0.5  # K
                    temperature = material._convert_to_nparray(temperature)
                    assert val[0] == material._fix_temperature_lb(
                        temperature, key
                    ), f"{material.__class__.__name__}"


def testTemperatureUpperBoundary():
    method_list = [
        attribute
        for attribute in dir(TemperatureDependentMaterial)
        if callable(getattr(TemperatureDependentMaterial, attribute))
        and attribute.startswith("_") == False
    ]
    for method in method_list:
        if (
            method != "glass_epoxy"
            and method != "incoloy"
            and method != "kapton"
        ):
            material = getattr(TemperatureDependentMaterial(), method)()
            for key, val in material.temperature_ranges.items():
                temperature = 2000.0  # K
                temperature = material._convert_to_nparray(temperature)
                assert val[1] == material._fix_temperature_ub(
                    temperature, key
                ), f"{material.__class__.__name__}"

        elif method == "glass_epoxy" or method == "incoloy" or method == "kapton":

            for source in range(2):
                material = getattr(TemperatureDependentMaterial(), method)(
                    source=source
                )
                for key, val in material.temperature_ranges.items():
                    temperature = 2000.0  # K
                    temperature = material._convert_to_nparray(temperature)
                    assert val[1] == material._fix_temperature_ub(
                        temperature, key
                    ), f"{material.__class__.__name__}"


def testPrelimiaryOperations():
    temperature = [ii for ii in range(10)]
    temperature = [*temperature, 2e3]

    method_list = [
        attribute
        for attribute in dir(TemperatureDependentMaterial)
        if callable(getattr(TemperatureDependentMaterial, attribute))
        and attribute.startswith("_") == False
    ]
    for method in method_list:
        if (
            method != "glass_epoxy"
            and method != "incoloy"
            and method != "kapton"
        ):
            material = getattr(TemperatureDependentMaterial(), method)()
            for key, val in material.temperature_ranges.items():
                material._preprocessing(temperature, key)
                expected = np.maximum(temperature, val[0])
                expected = np.minimum(expected, val[1])
                for ii, vv in enumerate(expected):
                    assert (
                        vv == material.temperature[ii]
                    ), f"{method = }, {vv = }, {material.temperature[ii] = }\n"
        elif method == "glass_epoxy" or method == "incoloy" or method == "kapton":
            for source in range(2):
                material = getattr(TemperatureDependentMaterial(), method)(
                    source=source
                )

                temperature = [ii for ii in range(10)]
                temperature = [*temperature, 2e3]
                material._preprocessing(temperature, key)
                expected = np.maximum(temperature, material.temperature_ranges[key][0])
                expected = np.minimum(expected, material.temperature_ranges[key][1])
                for ii, vv in enumerate(expected):
                    assert (
                        vv == material.temperature[ii]
                    ), f"{method = }, {vv = }, {material.temperature[ii] = }\n"
