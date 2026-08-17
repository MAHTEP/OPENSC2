import copy
import unittest

from utility_functions import checkpoint as checkpoint_module
from utility_functions.checkpoint import ContinuationProfile


class CheckpointContinuationProfileCompatibilityTests(unittest.TestCase):
    def _profile(self):
        return ContinuationProfile(
            immutable={
                "IADAPTIME": 0,
                "conductors": {
                    "COND_1": {
                        "inputs": {
                            "ZLENGTH": 10.0,
                            "METHOD": "BE",
                        },
                    },
                },
            },
            time_policy={
                "TIME_STEP": 0.025,
                "STPMIN": 0.0025,
                "TEND": 0.4,
                "conductors": {
                    "COND_1": {
                        "ELECTRIC_TIME_STEP": 0.001,
                    },
                },
            },
            drivers={
                "COND_1": {
                    "current": {
                        "source": {"kind": "canonical_input"},
                        "parameters": {
                            "I0_OP_MODE": 0,
                            "I0_OP_TOT": 12000.0,
                        },
                    },
                    "external_heat": {
                        "components": {
                            "STACK_1": {
                                "source": {"kind": "canonical_input"},
                                "parameters": {
                                    "IQFUN": 1,
                                    "Q0": 250.0,
                                },
                            },
                        },
                    },
                },
            },
        )

    def _compare(self, checkpoint_profile, runtime_profile):
        return checkpoint_module.compare_continuation_profiles(
            checkpoint_profile,
            runtime_profile,
        )

    def test_identical_profiles_are_compatible(self):
        checkpoint_profile = self._profile()
        runtime_profile = copy.deepcopy(checkpoint_profile)

        comparison = self._compare(
            checkpoint_profile,
            runtime_profile,
        )

        self.assertTrue(comparison.is_compatible)
        self.assertEqual(comparison.immutable_differences, ())
        self.assertEqual(comparison.time_policy_differences, ())
        self.assertEqual(comparison.driver_differences, ())

    def test_time_policy_and_driver_changes_are_allowed_and_reported(self):
        checkpoint_profile = self._profile()
        runtime_profile = copy.deepcopy(checkpoint_profile)
        runtime_profile.time_policy["TIME_STEP"] = 0.005
        runtime_profile.time_policy["TEND"] = 1.0
        runtime_profile.time_policy["conductors"]["COND_1"][
            "ELECTRIC_TIME_STEP"
        ] = 0.0005
        runtime_profile.drivers["COND_1"]["current"]["source"] = {
            "kind": "auxiliary_file",
            "path": "drivers/current.tsv",
            "sha256": "a" * 64,
            "size": 128,
        }
        runtime_profile.drivers["COND_1"]["current"]["parameters"][
            "I0_OP_MODE"
        ] = -1
        runtime_profile.drivers["COND_1"]["external_heat"][
            "components"
        ]["STACK_1"]["parameters"]["Q0"] = 500.0

        comparison = self._compare(
            checkpoint_profile,
            runtime_profile,
        )

        self.assertTrue(comparison.is_compatible)
        self.assertEqual(comparison.immutable_differences, ())
        self.assertEqual(
            comparison.time_policy_differences,
            (
                "time_policy.TEND",
                "time_policy.TIME_STEP",
                (
                    "time_policy.conductors.COND_1."
                    "ELECTRIC_TIME_STEP"
                ),
            ),
        )
        self.assertEqual(
            comparison.driver_differences,
            (
                (
                    "drivers.COND_1.current.parameters."
                    "I0_OP_MODE"
                ),
                "drivers.COND_1.current.source.kind",
                "drivers.COND_1.current.source.path",
                "drivers.COND_1.current.source.sha256",
                "drivers.COND_1.current.source.size",
                (
                    "drivers.COND_1.external_heat.components."
                    "STACK_1.parameters.Q0"
                ),
            ),
        )

    def test_immutable_added_missing_and_modified_values_are_blocking(self):
        checkpoint_profile = self._profile()
        runtime_profile = copy.deepcopy(checkpoint_profile)
        runtime_inputs = runtime_profile.immutable["conductors"]["COND_1"][
            "inputs"
        ]
        runtime_inputs["ZLENGTH"] = 12.0
        del runtime_inputs["METHOD"]
        runtime_inputs["NELEMS"] = 250

        comparison = self._compare(
            checkpoint_profile,
            runtime_profile,
        )

        self.assertFalse(comparison.is_compatible)
        self.assertEqual(
            comparison.immutable_differences,
            (
                "immutable.conductors.COND_1.inputs.METHOD",
                "immutable.conductors.COND_1.inputs.NELEMS",
                "immutable.conductors.COND_1.inputs.ZLENGTH",
            ),
        )

    def test_comparison_is_deterministic_and_non_mutating(self):
        checkpoint_profile = self._profile()
        runtime_profile = copy.deepcopy(checkpoint_profile)
        runtime_profile.time_policy.update(
            Z_LAST=3,
            A_FIRST=1,
        )
        runtime_profile.drivers["COND_1"]["current"]["parameters"].update(
            Z_LAST=3,
            A_FIRST=1,
        )
        checkpoint_before = copy.deepcopy(checkpoint_profile)
        runtime_before = copy.deepcopy(runtime_profile)

        first = self._compare(checkpoint_profile, runtime_profile)
        second = self._compare(checkpoint_profile, runtime_profile)

        self.assertEqual(first, second)
        self.assertEqual(
            first.time_policy_differences,
            (
                "time_policy.A_FIRST",
                "time_policy.Z_LAST",
            ),
        )
        self.assertEqual(
            first.driver_differences,
            (
                (
                    "drivers.COND_1.current.parameters."
                    "A_FIRST"
                ),
                (
                    "drivers.COND_1.current.parameters."
                    "Z_LAST"
                ),
            ),
        )
        self.assertEqual(checkpoint_profile, checkpoint_before)
        self.assertEqual(runtime_profile, runtime_before)


if __name__ == "__main__":
    unittest.main()
