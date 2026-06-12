# Copyright (c) 2021-2026 Cubillos & Blecic
# Pyrat Bay is open-source software under the GPL-2.0 license (see LICENSE)

import os
import re
import pytest

import numpy as np
import pyratbay.spectrum as ps
import pyratbay.constants as pc

os.chdir(pc.ROOT + 'tests')


expected_eps_spot_phoenix = np.array([
    1.00520325, 1.00450367, 1.0036623 , 1.00385527, 1.00365159,
    1.0035096 , 1.00388583, 1.00447088, 1.00488525, 1.00526333,
    1.00508391, 1.00332297, 1.00343777, 1.00435936, 1.00448354,
    1.00403631, 1.00445369, 1.00447114, 1.00476787, 1.00375426,
    1.00471423, 1.00392382, 1.0040069 , 1.00392366, 1.00423105,
    1.0039243 , 1.00467319, 1.00407956, 1.00357499, 1.00453642,
    1.00428823, 1.00367059, 1.00340318, 1.00361182, 1.00342656,
    1.00373367, 1.00369197, 1.00337278, 1.00342969, 1.00328798,
    1.00331867, 1.00385249, 1.00350894, 1.00320202, 1.00332404,
    1.00338611, 1.00336859, 1.00351523, 1.00331752, 1.00332205,
    1.00354256, 1.00327367, 1.00329523, 1.0036399 , 1.00318012,
    1.00253936, 1.00350947, 1.00326464, 1.00335974, 1.00329767,
    1.00350311, 1.00409591, 1.00333869, 1.00385898, 1.003493  ,
    1.00349184, 1.00311293, 1.00385737, 1.00340929, 1.00337596,
    1.00360386, 1.00373383, 1.00369698, 1.00373902, 1.0037085 ,
    1.00380252, 1.00351968, 1.00394672, 1.00326631, 1.00358472,
    1.00355464, 1.0023713 , 1.00344563, 1.00328757, 1.00333676,
    1.00361614, 1.00352429, 1.00327353, 1.00325995, 1.00371426,
    1.0035305 , 1.00310589, 1.00316837, 1.0035944 , 1.00346402,
    1.00323172, 1.00324688, 1.00346998, 1.00321092, 1.00340294,
])

expected_eps_spot_bin = np.array([
    1.00330436, 1.00314108, 1.00306117, 1.00298392, 1.00311167,
    1.00316829, 1.00323794, 1.00311764, 1.00303873, 1.00288806,
    1.00291419, 1.00318762, 1.00303618, 1.00284281, 1.00268267,
    1.0027962 , 1.00271896, 1.00257293, 1.00264637, 1.00260749,
    1.00258784, 1.00283469, 1.00271167, 1.00267288, 1.0027025 ,
    1.00268291, 1.00251774, 1.00249096, 1.00243169, 1.00247853,
    1.00244112, 1.00226395, 1.00228434, 1.00233399, 1.00235605,
    1.00239541, 1.0022577 , 1.00230785, 1.00223687, 1.00225782,
    1.00227445, 1.00224631, 1.00217975, 1.00217524, 1.00220036,
    1.00217279, 1.00213575, 1.00213339, 1.00217072, 1.00213685,
    1.00211011, 1.00204395, 1.0020247 , 1.00203603, 1.0019961 ,
    1.00198941, 1.00195935, 1.00195522, 1.0018906 , 1.0018349 ,
    1.00180014, 1.00180023, 1.00170936, 1.00176077, 1.0017647 ,
    1.00162694, 1.00168688, 1.00152079, 1.00160271, 1.0014927 ,
    1.00150938, 1.00146087, 1.00138523, 1.00144724, 1.00132384,
    1.00134265, 1.00133126, 1.00137874, 1.00130811, 1.00129764,
    1.00128938, 1.0013356 , 1.00129068, 1.00131883, 1.00141441,
    1.00120305, 1.00121989, 1.00134247, 1.00142237, 1.00147455,
    1.00128085, 1.00157483, 1.0012911 , 1.00123853, 1.0012291 ,
    1.00122952, 1.00126698, 1.00127082, 1.00117079, 1.00123367,
    1.00117921, 1.00118055, 1.00120497, 1.00117643, 1.00122339,
    1.00115186, 1.00116892, 1.00118469, 1.00118602, 1.00119201,
    1.0011965 , 1.00118651, 1.00125476, 1.00127934, 1.00122138,
    1.00124935, 1.00120266, 1.00126044, 1.00133831, 1.00113473,
    1.00123417, 1.00120946, 1.0011472 , 1.0013145 , 1.0010907 ,
    1.00130677, 1.00120242, 1.00114848, 1.00117381, 1.00116387,
    1.00124264, 1.00128661, 1.00124439,
])

expected_eps_spot_interpolate = np.array([
    1.00337065, 1.0031303 , 1.00309687, 1.0031635 , 1.00328219,
    1.00310897, 1.00306129, 1.00326358, 1.00308941, 1.00287629,
    1.00269815, 1.00301214, 1.00291473, 1.00263166, 1.00227876,
    1.00239583, 1.00261806, 1.00262259, 1.00252037, 1.00272443,
    1.00268466, 1.00274423, 1.00258078, 1.00270003, 1.00258397,
    1.00254398, 1.00247615, 1.0025018 , 1.00207829, 1.00247585,
    1.0023911 , 1.00235617, 1.00235615, 1.00207148, 1.00224655,
    1.0023478 , 1.00209397, 1.00230591, 1.00228111, 1.00224356,
    1.00225617, 1.0022957 , 1.00220657, 1.00219108, 1.00216296,
    1.00215691, 1.00213532, 1.00224023, 1.00232356, 1.00210137,
    1.00207962, 1.00206472, 1.00201049, 1.00202516, 1.00213961,
    1.00198756, 1.00195582, 1.00192978, 1.00193141, 1.00184147,
    1.00171137, 1.0017591 , 1.0014259 , 1.00175749, 1.0017467 ,
    1.00170553, 1.00142768, 1.00159969, 1.00153524, 1.00136582,
    1.00143727, 1.00219191, 1.00137131, 1.00146734, 1.0014491 ,
    1.00136001, 1.00146244, 1.00118037, 1.0014112 , 1.00133849,
    1.00132816, 1.0008815 , 1.00131252, 1.00129285, 1.00132418,
    1.0013758 , 1.0011903 , 1.00131594, 1.00128962, 1.00145421,
    1.00125466, 1.00194241, 1.00128671, 1.00126839, 1.00115838,
    1.00121655, 1.00124696, 1.00137633, 1.00118659, 1.00124679,
    1.00116438, 1.00115084, 1.00110659, 1.00109487, 1.00147671,
    1.00112892, 1.00116094, 1.00115295, 1.00121244, 1.00113913,
    1.00155844, 1.00114302, 1.00157691, 1.0011182 , 1.00121721,
    1.00123784, 1.00108424, 1.00120044, 1.00166204, 1.00111497,
    1.0010612 , 1.00116699, 1.00111032, 1.00220404, 1.00113549,
    1.00158751, 1.00121853, 1.00109147, 1.0010509 , 1.00106227,
    1.0012345 , 1.00111076, 1.00111989,
])

expected_eps_spot_fac = np.array([
    1.00126218, 1.00113355, 1.00108133, 1.0010702 , 1.00112892,
    1.00120265, 1.0013039 , 1.00123887, 1.00111982, 1.00103691,
    1.00106197, 1.00129889, 1.00115168, 1.00107891, 1.00097168,
    1.001011  , 1.00102268, 1.00095342, 1.00098427, 1.00091617,
    1.00098645, 1.00106262, 1.00097882, 1.0010145 , 1.0010112 ,
    1.00101553, 1.00096723, 1.00093125, 1.00090819, 1.00090913,
    1.00089085, 1.00087031, 1.00089247, 1.00090728, 1.00090396,
    1.00095526, 1.00088481, 1.00090824, 1.00089804, 1.00089436,
    1.00092145, 1.00092084, 1.00090511, 1.00090716, 1.00093673,
    1.00093389, 1.00092741, 1.00093514, 1.00095401, 1.00093971,
    1.00094097, 1.00091809, 1.00094205, 1.00094269, 1.00095209,
    1.00096305, 1.00095807, 1.00097099, 1.0009271 , 1.0009253 ,
    1.00090688, 1.00095036, 1.00086375, 1.0009029 , 1.00096969,
    1.0008865 , 1.00095493, 1.00085837, 1.00081887, 1.00081792,
    1.00084166, 1.00085033, 1.00080509, 1.00084468, 1.00079527,
    1.00077792, 1.00078855, 1.00083292, 1.00078183, 1.00078291,
    1.00077157, 1.00080614, 1.00075336, 1.00075197, 1.00081926,
    1.0007075 , 1.0007031 , 1.00078797, 1.00082923, 1.00081664,
    1.00074774, 1.00086349, 1.00076188, 1.00073286, 1.00071007,
    1.00072721, 1.00075523, 1.00072574, 1.00068664, 1.00071256,
    1.0006865 , 1.00068457, 1.00065884, 1.00063823, 1.00068279,
    1.00064943, 1.00063484, 1.00065132, 1.00065096, 1.00057408,
    1.00060716, 1.00063797, 1.00066521, 1.000691  , 1.00068261,
    1.00065137, 1.00067616, 1.00069354, 1.00076797, 1.00061448,
    1.00067185, 1.00065371, 1.00063808, 1.00078474, 1.00058792,
    1.00076536, 1.0006978 , 1.00065067, 1.00067638, 1.00067139,
    1.00073343, 1.00076594, 1.00073181,
])

def test_tls_init_phoenix():
    sed_folder = f'{pc.ROOT}tests/inputs/'
    teff = 4600.0
    tls = ps.TransitLightSource(sed_folder, teff)

    np.testing.assert_allclose(tls.teff, teff)
    np.testing.assert_allclose(tls.temps, [4000, 4500, 5000])
    np.testing.assert_allclose(tls.wl[0], 0.7)
    np.testing.assert_allclose(tls.wl[-1], 4.9899999)
    assert len(tls.wl) == 2500

    eps = tls.epsilon(4200, 0.01)
    np.testing.assert_allclose(eps[0:100], expected_eps_spot_phoenix)


def test_tls_init_bin():
    sed_folder = f'{pc.ROOT}tests/inputs/'
    teff = 4600.0
    wl = ps.constant_resolution_spectrum(0.8, 3.0, resolution=100.0)
    tls = ps.TransitLightSource(sed_folder, teff, wl, sampling='bin')

    np.testing.assert_allclose(tls.teff, teff)
    np.testing.assert_allclose(tls.temps, [4000, 4500, 5000])
    np.testing.assert_allclose(tls.wl[0], 0.8)
    np.testing.assert_allclose(tls.wl[-1], 2.994770044)
    assert len(tls.wl) == 133

    eps = tls.epsilon(4200, 0.01)
    np.testing.assert_allclose(eps, expected_eps_spot_bin)


@pytest.mark.parametrize('method', (None, 'interpolate'))
def test_tls_init_interpolate(method):
    sed_folder = f'{pc.ROOT}tests/inputs/'
    teff = 4600.0
    wl = ps.constant_resolution_spectrum(0.8, 3.0, resolution=100.0)
    if method is None:
        tls = ps.TransitLightSource(sed_folder, teff, wl)
    else:
        tls = ps.TransitLightSource(sed_folder, teff, wl, sampling='interpolate')

    np.testing.assert_allclose(tls.teff, teff)
    np.testing.assert_allclose(tls.temps, [4000, 4500, 5000])
    np.testing.assert_allclose(tls.wl[0], 0.8)
    np.testing.assert_allclose(tls.wl[-1], 2.994770044)
    assert len(tls.wl) == 133

    eps = tls.epsilon(4200, 0.01)
    np.testing.assert_allclose(eps, expected_eps_spot_interpolate)


def test_tls_init_spot_fac():
    sed_folder = f'{pc.ROOT}tests/inputs/'
    teff = 4600.0
    wl = ps.constant_resolution_spectrum(0.8, 3.0, resolution=100.0)
    tls = ps.TransitLightSource(sed_folder, teff, wl, sampling='bin')

    eps = tls.epsilon(4200, 0.01, 4800, 0.01)
    np.testing.assert_allclose(eps, expected_eps_spot_fac)


def test_tls_out_of_bounds_temp():
    sed_folder = f'{pc.ROOT}tests/inputs/'
    teff = 4600.0
    tls = ps.TransitLightSource(sed_folder, teff)

    eps = tls.epsilon(3200, 0.01)
    assert np.all(np.isnan(eps))


def test_resample_out_of_bounds_fractions():
    sed_folder = f'{pc.ROOT}tests/inputs/'
    teff = 4600.0
    tls = ps.TransitLightSource(sed_folder, teff)

    error_msg = re.escape("Unphysical spot+faculae coverage > 1")
    with pytest.raises(ValueError, match=error_msg):
        tls.epsilon(4200, 0.6, 4800, 0.6)


@pytest.mark.parametrize('f_spot', (-0.5, 1.5))
def test_tls_out_of_bounds_f_spot(f_spot):
    sed_folder = f'{pc.ROOT}tests/inputs/'
    teff = 4600.0
    tls = ps.TransitLightSource(sed_folder, teff)

    error_msg = "Out of bounds f_spot fraction, value must be between 0 and 1"
    with pytest.raises(ValueError, match=error_msg):
        tls.epsilon(4200, f_spot)


@pytest.mark.parametrize('f_fac', (-0.5, 1.5))
def test_tls_out_of_bounds_f_fac(f_fac):
    sed_folder = f'{pc.ROOT}tests/inputs/'
    teff = 4600.0
    tls = ps.TransitLightSource(sed_folder, teff)

    error_msg = "Out of bounds f_fac fraction, value must be between 0 and 1"
    with pytest.raises(ValueError, match=error_msg):
        tls.epsilon(4200, 0.01, 4500, f_fac)


def test_tls_missing_folder():
    sed_folder = f'{pc.ROOT}tests/inputs_not_here/'
    teff = 4600.0

    error_msg = "SED folder for TLS model not found"
    with pytest.raises(ValueError, match=error_msg):
        ps.TransitLightSource(sed_folder, teff)


def test_tls_missing_models():
    sed_folder = f'{pc.ROOT}tests/'
    teff = 4600.0

    error_msg = "No PHOENIX models found in TLS folder"
    with pytest.raises(ValueError, match=error_msg):
        ps.TransitLightSource(sed_folder, teff)


def test_tls_bad_teff():
    sed_folder = f'{pc.ROOT}tests/inputs'
    error_msg = re.escape(
        'Effective temperature (3200) is not in range of SED '
        'temperatures (4000.0, 5000.0)'
    )
    with pytest.raises(ValueError, match=error_msg):
        ps.TransitLightSource(sed_folder, 3200)


def test_tls_bad_sampling():
    sed_folder = f'{pc.ROOT}tests/inputs'
    teff = 4600.0

    error_msg = "Invalid wavelength sampling method"
    with pytest.raises(ValueError, match=error_msg):
        ps.TransitLightSource(sed_folder, teff, sampling='jump')

