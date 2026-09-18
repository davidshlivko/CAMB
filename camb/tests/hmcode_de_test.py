import os
import sys
import tempfile
import unittest

import numpy as np

try:
    import camb
except ImportError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
    import camb

from camb.nonlinear import Halofit, halofit_version_names  # type: ignore[attr-defined]


class HMcodeDarkEnergyTest(unittest.TestCase):
    @staticmethod
    def _get_pk(halofit_version, *, tabulated=False):
        pars = camb.CAMBparams(WantCls=False)
        pars.set_cosmology(H0=67.5, ombh2=0.02237, omch2=0.12, mnu=0.06)
        pars.InitPower.set_params(As=2.1e-9, ns=0.965)
        if tabulated:
            a = np.linspace(0.05, 1.0, 80)
            w = -1.0 + 0.25 * np.sin(np.pi * (1.0 - a) / 0.95) ** 2
            pars.set_dark_energy_w_a(a, w, dark_energy_model="ppf")
        else:
            pars.set_dark_energy(w=-0.9, wa=0.2, dark_energy_model="ppf")
        pars.set_matter_power(redshifts=[0.5], kmax=2.0, k_per_logint=10)
        pars.NonLinear = camb.model.NonLinear_both
        pars.NonLinearModel.set_params(halofit_version=halofit_version, HMCode_logT_AGN=7.8)
        results = camb.get_results(pars)
        return results.get_nonlinear_matter_power_spectrum()[2]

    def test_names_and_ini_round_trip(self):
        self.assertEqual(halofit_version_names["mead2020_de"], 11)
        self.assertEqual(halofit_version_names["mead2020_feedback_de"], 12)

        for version in ("mead2020_de", "mead2020_feedback_de"):
            with self.subTest(version=version), tempfile.TemporaryDirectory() as temp_dir:
                pars = camb.CAMBparams(WantCls=False)
                pars.NonLinearModel = Halofit()
                pars.NonLinearModel.set_params(halofit_version=version)
                ini_file = os.path.join(temp_dir, "params.ini")
                pars.write_ini(ini_file)
                round_tripped = camb.read_ini(ini_file)
                self.assertEqual(round_tripped.NonLinearModel.get_halofit_version(), version)

    def test_cpl_matches_legacy_modes(self):
        for legacy, exact_de in (
            ("mead2020", "mead2020_de"),
            ("mead2020_feedback", "mead2020_feedback_de"),
        ):
            with self.subTest(legacy=legacy, exact_de=exact_de):
                legacy_pk = self._get_pk(legacy)
                exact_de_pk = self._get_pk(exact_de)
                np.testing.assert_allclose(exact_de_pk, legacy_pk, rtol=5e-12, atol=0.0)

    def test_tabulated_background_changes_new_modes(self):
        for legacy, exact_de in (
            ("mead2020", "mead2020_de"),
            ("mead2020_feedback", "mead2020_feedback_de"),
        ):
            with self.subTest(legacy=legacy, exact_de=exact_de):
                legacy_pk = self._get_pk(legacy, tabulated=True)
                exact_de_pk = self._get_pk(exact_de, tabulated=True)
                self.assertTrue(np.all(np.isfinite(exact_de_pk)))
                self.assertGreater(np.max(np.abs(exact_de_pk / legacy_pk - 1.0)), 1e-3)


if __name__ == "__main__":
    unittest.main()
