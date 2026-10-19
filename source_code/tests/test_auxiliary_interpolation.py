import unittest

import numpy as np
import pandas as pd

from utility_functions.auxiliary_functions import (
    build_interpolator,
    do_interpolation,
)


class AuxiliaryInterpolationTests(unittest.TestCase):
    def test_linear_space_only_profile_returns_callable_interpolator(self):
        profile = pd.DataFrame(
            [
                [1.0, 0.0],
                [0.0, 10.0],
                [5.0, 20.0],
                [10.0, 30.0],
            ]
        )

        interpolator, kind = build_interpolator(profile, "linear")
        values = do_interpolation(
            interpolator,
            np.array([-1.0, 0.0, 2.5, 7.5, 10.0, 11.0]),
            123.0,
            kind,
        )

        self.assertEqual(kind, "space_only")
        np.testing.assert_allclose(
            values,
            np.array([10.0, 10.0, 15.0, 25.0, 30.0, 30.0]),
        )

    def test_linear_time_only_profile_returns_callable_interpolator(self):
        profile = pd.DataFrame(
            [
                [1.0, 0.0, 5.0, 10.0],
                [0.0, 10.0, 20.0, 30.0],
            ]
        )

        interpolator, kind = build_interpolator(profile, "linear")

        self.assertEqual(kind, "time_only")
        self.assertEqual(
            float(do_interpolation(interpolator, np.array([99.0]), 2.5, kind)),
            15.0,
        )
        self.assertEqual(
            float(do_interpolation(interpolator, np.array([99.0]), -1.0, kind)),
            10.0,
        )
        self.assertEqual(
            float(do_interpolation(interpolator, np.array([99.0]), 11.0, kind)),
            30.0,
        )


if __name__ == "__main__":
    unittest.main()
