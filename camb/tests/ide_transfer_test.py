"""Linear matter-transfer checks for exponential coupled CDM."""

import importlib.util
import unittest
from pathlib import Path

import numpy as np

import camb
from camb import model

OMB = 0.02237
OMC_EARLY = 0.12
MNU = 0.06
AS = 2.10e-9
NS = 0.9649


def ide_params(*, beta=0.0, lam=0.0, accuracy=1.0, kmax=5.0, **ide_kwargs):
    pars = camb.CAMBparams()
    pars.set_ide_cosmology(
        V0_over_rho100=0.31,
        omc_early=OMC_EARLY,
        beta=beta,
        lam=lam,
        ombh2=OMB,
        mnu=MNU,
        tau=0.0544,
        **ide_kwargs,
    )
    pars.InitPower.set_params(As=AS, ns=NS)
    pars.set_accuracy(
        AccuracyBoost=accuracy,
        lSampleBoost=accuracy,
        lAccuracyBoost=accuracy,
        DoLateRadTruncation=False,
    )
    pars.set_for_lmax(300, lens_potential_accuracy=0, nonlinear=False)
    pars.set_matter_power(
        redshifts=[0.0, 0.5, 3.0],
        kmax=kmax,
        accurate_massive_neutrino_transfers=True,
        nonlinear=False,
        silent=True,
    )
    return pars


def ordinary_params(h0, *, accuracy=1.0, kmax=5.0):
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=h0, ombh2=OMB, omch2=OMC_EARLY, mnu=MNU, tau=0.0544)
    pars.InitPower.set_params(As=AS, ns=NS)
    pars.set_accuracy(
        AccuracyBoost=accuracy,
        lSampleBoost=accuracy,
        lAccuracyBoost=accuracy,
        DoLateRadTruncation=False,
    )
    pars.set_for_lmax(300, lens_potential_accuracy=0, nonlinear=False)
    pars.set_matter_power(
        redshifts=[0.0, 0.5, 3.0],
        kmax=kmax,
        accurate_massive_neutrino_transfers=True,
        nonlinear=False,
        silent=True,
    )
    return pars


class IDETransferTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.uncoupled = camb.get_results(ide_params())
        cls.ordinary = camb.get_results(ordinary_params(cls.uncoupled.Params.H0))

    def test_lcdm_linear_matter_power_and_sigma8_parity(self):
        self.assertEqual(self.uncoupled.Params.NonLinear, model.NonLinear_none)
        for name in ("get_sigma8", "get_sigma8_0"):
            with self.subTest(output=name):
                np.testing.assert_allclose(
                    getattr(self.uncoupled, name)(),
                    getattr(self.ordinary, name)(),
                    rtol=3e-3,
                )
        np.testing.assert_allclose(self.uncoupled.get_sigmaR(8), self.ordinary.get_sigmaR(8), rtol=3e-3)

        kh_ide, z_ide, pk_ide = self.uncoupled.get_linear_matter_power_spectrum()
        kh_std, z_std, pk_std = self.ordinary.get_linear_matter_power_spectrum()
        np.testing.assert_allclose(z_ide, z_std, rtol=0, atol=0)
        # IDE modes begin at a fixed early epoch, so CAMB chooses a slightly
        # different internal k grid. Compare the same physical wavenumbers.
        mask = (kh_ide > 1e-3) & (kh_ide < 1.0)
        for ide, std in zip(pk_ide, pk_std):
            std_at_ide_k = np.exp(np.interp(np.log(kh_ide[mask]), np.log(kh_std), np.log(std)))
            np.testing.assert_allclose(ide[mask], std_at_ide_k, rtol=5e-3)

        ide_trans = self.uncoupled.get_matter_transfer_data()
        std_trans = self.ordinary.get_matter_transfer_data()
        for name in ("delta_cdm", "delta_baryon", "delta_nu", "delta_tot", "delta_nonu"):
            with self.subTest(transfer=name):
                ide = ide_trans.transfer_z(name, z_index=2)
                std = std_trans.transfer_z(name, z_index=2)
                q = ide_trans.q
                mask = (q > 1e-3) & (q < 1.0)
                std_at_ide_q = np.interp(np.log(q[mask]), np.log(std_trans.q), std)
                np.testing.assert_allclose(ide[mask], std_at_ide_q, rtol=5e-3)

    def test_coupled_transfers_match_mode_evolution(self):
        for beta, lam in ((0.07, 0.7), (0.07, -0.7), (-0.07, -0.7)):
            with self.subTest(beta=beta, lam=lam):
                result = camb.get_results(ide_params(beta=beta, lam=lam))
                transfer = result.get_matter_transfer_data()
                q = transfer.q
                indices = np.unique(np.searchsorted(q, [0.003, 0.01, 0.05, 0.2, 0.7]))
                indices = indices[indices < len(q)]
                q = q[indices]
                evolution = result.get_redshift_evolution(q, [0.5], ["delta_cdm", "delta_baryon"])
                for j, name in enumerate(("delta_cdm", "delta_baryon")):
                    transfer_delta = transfer.transfer_z(name, z_index=1)[indices] * q**2
                    np.testing.assert_allclose(transfer_delta, evolution[:, 0, j], rtol=1e-2, atol=1e-7)
                sigma8 = result.get_sigma8()
                self.assertTrue(np.all(np.isfinite(sigma8)))
                self.assertTrue(np.all(sigma8 > 0))
                _, _, pk = result.get_linear_matter_power_spectrum()
                self.assertTrue(np.all(np.isfinite(pk)))
                self.assertTrue(np.all(pk > 0))
                np.testing.assert_allclose(result.get_sigma8_0(), result.get_sigmaR(8)[-1], rtol=5e-3)

    def test_coupled_sigma8_numerical_convergence(self):
        base = camb.get_results(ide_params(beta=0.07, lam=0.7, kmax=5.0))
        refined = camb.get_results(
            ide_params(
                beta=0.07,
                lam=0.7,
                kmax=10.0,
                accuracy=2.0,
                a_init=1e-12,
                history_points=4097,
                integrate_tol=1e-11,
            )
        )
        np.testing.assert_allclose(base.get_sigma8(), refined.get_sigma8(), rtol=1e-2)

    def test_repeated_set_params_preserves_matter_transfer_request(self):
        pars = camb.set_params(
            dark_energy_model="IDEExponential",
            V0_over_rho100=0.31,
            omc_early=OMC_EARLY,
            beta=0.07,
            lam=0.7,
            ombh2=OMB,
            mnu=MNU,
            nnu=3.046,
            tau=0.0544,
        )
        pars.WantTransfer = True
        updated = camb.set_params(
            pars.copy(),
            V0_over_rho100=0.31,
            omc_early=OMC_EARLY,
            beta=0.07,
            lam=0.7,
            ombh2=OMB,
            mnu=MNU,
            nnu=3.046,
            tau=0.0544,
        )
        self.assertTrue(updated.WantTransfer)

    @unittest.skipUnless(importlib.util.find_spec("cobaya"), "Cobaya not installed")
    def test_cobaya_can_derive_sigma8_from_ide_cmb_run(self):
        from cobaya.model import get_model

        def dummy_cmb_likelihood(_self):
            # Request a scalar CMB spectrum without requiring external data.
            return 0.0 * _self.provider.get_Cl(ell_factor=True)["tt"][30]

        info = {
            "theory": {
                "camb": {
                    "path": str(Path(__file__).resolve().parents[2]),
                    "extra_args": {
                        "dark_energy_model": "IDEExponential",
                        "lens_potential_accuracy": 0,
                        "max_eta_k": 18000,
                        "DoLateRadTruncation": False,
                    },
                }
            },
            "likelihood": {"dummy_cmb": {"external": dummy_cmb_likelihood, "requires": {"Cl": {"tt": 100}}}},
            "params": {
                "V0_over_rho100": 0.31,
                "omc_early": OMC_EARLY,
                "beta": 0.07,
                "lam": 0.7,
                "ombh2": OMB,
                "mnu": MNU,
                "nnu": 3.046,
                "tau": 0.0544,
                "As": AS,
                "ns": NS,
                "H0": {"derived": True},
                "omch2": {"derived": True},
                "sigma8": {"derived": True},
            },
        }
        with get_model(info) as cobaya_model:
            posterior = cobaya_model.logposterior({}, cached=False)
            derived = dict(zip(cobaya_model.parameterization.derived_params(), posterior.derived))
        self.assertTrue(np.isfinite(posterior.logpost))
        self.assertGreater(derived["sigma8"], 0)
        self.assertGreater(derived["H0"], 0)
        self.assertGreater(derived["omch2"], 0)


if __name__ == "__main__":
    unittest.main()
