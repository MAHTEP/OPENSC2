import numpy as np
from properties_of_materials.aluminium import (
    density_al6063_cryosoft, isobaric_specific_heat_al6063_cryosoft,
    thermal_conductivity_al6063_cryosoft, electrical_resistivity_al6063_cryosoft,
    thermal_conductivity_al_cryosoft, electrical_resistivity_al_cryosoft,
)
from strand_stabilizer_component import DENSITY_FUNC as SD, ISOBARIC_SPECIFIC_HEAT_FUNC as SC, THERMAL_CONDUCTIVITY_FUNC as SK, ELECTRICAL_RESISTIVITY_FUNC as SR
from jacket_component import DENSITY_FUNC as JD, ISOBARIC_SPECIFIC_HEAT_FUNC as JC, THERMAL_CONDUCTIVITY_FUNC as JK, ELECTRICAL_RESISTIVITY_FUNC as JR

T = np.array([4.5,10.,20.,40.,60.,80.,100.,120.,150.])

def test_al6063_against_compiled_fortran():
    np.testing.assert_allclose(density_al6063_cryosoft(T), 2700.0, rtol=0, atol=0)
    np.testing.assert_allclose(isobaric_specific_heat_al6063_cryosoft(T), [0.31015554,1.4300442,8.8538723,76.642509,213.02734,358.52734,483.50391,581.125,687.13672], rtol=1e-4)
    np.testing.assert_allclose(thermal_conductivity_al6063_cryosoft(T), [16.910074,41.921383,87.758591,158.80775,198.00609,216.03241,222.81277,224.14310,222.02826], rtol=3e-6)
    np.testing.assert_allclose(electrical_resistivity_al6063_cryosoft(T), [5.7601364e-9,5.7623382e-9,5.7875775e-9,6.0852514e-9,7.1350748e-9,9.1250065e-9,1.1214644e-8,1.3377372e-8,1.6862764e-8], rtol=3e-6)

def test_al6063_is_distinct_from_pure_al_r300_b17():
    b = np.full(T.shape, 17.0)
    kp = thermal_conductivity_al_cryosoft(T,b,300.0)
    rp = electrical_resistivity_al_cryosoft(T,b,300.0)
    assert np.max(kp/thermal_conductivity_al6063_cryosoft(T)) > 10.0
    assert np.max(electrical_resistivity_al6063_cryosoft(T)/rp) > 10.0

def test_dispatch_maps():
    key='al6063_cryosoft'
    assert SD[key] is density_al6063_cryosoft and JD[key] is density_al6063_cryosoft
    assert SC[key] is isobaric_specific_heat_al6063_cryosoft and JC[key] is isobaric_specific_heat_al6063_cryosoft
    assert SK[key] is thermal_conductivity_al6063_cryosoft and JK[key] is thermal_conductivity_al6063_cryosoft
    assert SR[key] is electrical_resistivity_al6063_cryosoft and JR[key] is electrical_resistivity_al6063_cryosoft
