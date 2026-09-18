import os
import sys
import unittest

import numpy as np

try:
    import camb
except ImportError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
    import camb

from camb.dark_energy import DarkEnergyPPF  # type: ignore[attr-defined]


class PadeDarkEnergyTest(unittest.TestCase):
    @staticmethod
    def _expected_w(a, pade_eps0, pade_lneta0):
        return -1.0 + 2.0 * pade_eps0 / (3.0 + np.exp(pade_lneta0) * (a**-3 - 1.0))

    def _assert_background_matches_pade(self, pade_eps0, pade_lneta0):
        pars = camb.CAMBparams(WantCls=False)
        pars.set_cosmology(H0=67.5, ombh2=0.02237, omch2=0.12, mnu=0.06)
        pars.set_dark_energy(pade_eps0=pade_eps0, pade_lneta0=pade_lneta0, dark_energy_model="ppf")
        background = camb.get_background(pars, no_thermo=True)
        a = np.array([1e-3, 1e-2, 0.1, 0.5, 0.9, 1.0])
        _rho, w = background.get_dark_energy_rho_w(a)
        np.testing.assert_allclose(w, self._expected_w(a, pade_eps0, pade_lneta0), rtol=2e-7, atol=2e-10)

    def test_cpl_behavior_unchanged(self):
        model = DarkEnergyPPF()
        model.set_params(w=-0.9, wa=0.2, cs2=0.8)
        self.assertFalse(model.use_tabulated_w)
        self.assertEqual(model.w, -0.9)
        self.assertEqual(model.wa, 0.2)
        self.assertEqual(model.cs2, 0.8)

    def test_pade_public_entry_points(self):
        model = DarkEnergyPPF()
        model.set_params(pade_eps0=0.3, pade_lneta0=2.0)
        self.assertTrue(model.use_tabulated_w)

        pars = camb.CAMBparams(WantCls=False)
        pars.set_dark_energy(pade_eps0=0.3, pade_lneta0=2.0, dark_energy_model="ppf")
        self.assertIsInstance(pars.DarkEnergy, DarkEnergyPPF)
        self.assertTrue(pars.DarkEnergy.use_tabulated_w)

        pars = camb.set_params(pade_eps0=0.3, pade_lneta0=2.0, dark_energy_model="ppf")
        self.assertIsInstance(pars.DarkEnergy, DarkEnergyPPF)
        self.assertTrue(pars.DarkEnergy.use_tabulated_w)

    def test_pade_background(self):
        for pade_eps0, pade_lneta0 in ((0.3, 2.0), (-0.2, 2.0), (1e-5, -8.0), (1.5, 18.0)):
            with self.subTest(pade_eps0=pade_eps0, pade_lneta0=pade_lneta0):
                self._assert_background_matches_pade(pade_eps0, pade_lneta0)

    def test_zero_epsilon_is_lcdm(self):
        self._assert_background_matches_pade(0.0, 18.0)

    def test_pade_params_are_cobaya_visible(self):
        valid = camb.get_valid_numerical_params(dark_energy_model="ppf")
        self.assertIn("pade_eps0", valid)
        self.assertIn("pade_lneta0", valid)

    def test_invalid_pade_inputs(self):
        invalid = (
            {"pade_eps0": 0.3},
            {"pade_lneta0": 2.0},
            {"pade_eps0": 0.3, "pade_lneta0": 2.0, "w": -0.9},
            {"pade_eps0": 0.3, "pade_lneta0": 2.0, "wa": 0.1},
            {"pade_eps0": 0.3, "pade_lneta0": 2.0, "use_tabulated_w": True},
            {"pade_eps0": 0.3, "pade_lneta0": 2.0, "wde_a_array": np.array([0.5, 1.0])},
            {"pade_eps0": np.nan, "pade_lneta0": 2.0},
            {"pade_eps0": 0.3, "pade_lneta0": np.inf},
            {"pade_eps0": 0.3, "pade_lneta0": 1000.0},
            {"pade_eps0": 1e308, "pade_lneta0": 2.0},
        )
        for kwargs in invalid:
            with self.subTest(kwargs=kwargs), self.assertRaises(camb.CAMBValueError):
                DarkEnergyPPF().set_params(**kwargs)


if __name__ == "__main__":
    unittest.main()
