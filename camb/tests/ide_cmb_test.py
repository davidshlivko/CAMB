"""Scalar-CMB and mode-evolution regression checks for exponential coupled CDM."""

import unittest

import numpy as np

import camb
from camb import constants, model
from camb.baseconfig import CAMBError

V0 = 0.31
OMC_EARLY = 0.12
OMB = 0.02237
AS = 2.10e-9
NS = 0.9649
TAU = 0.0544
MNU = 0.06
LMAX = 1000


def ide_params(*, beta=0.0, lam=0.0, lmax=LMAX, accuracy=1.0, **kwargs):
    pars = camb.CAMBparams()
    pars.set_ide_cosmology(
        V0_over_rho100=V0,
        omc_early=OMC_EARLY,
        beta=beta,
        lam=lam,
        ombh2=OMB,
        mnu=MNU,
        tau=TAU,
        **kwargs,
    )
    pars.InitPower.set_params(As=AS, ns=NS)
    pars.set_accuracy(
        AccuracyBoost=accuracy,
        lSampleBoost=accuracy,
        lAccuracyBoost=accuracy,
        DoLateRadTruncation=False,
    )
    pars.Accuracy.TimeStepBoost = accuracy
    pars.set_for_lmax(lmax, lens_potential_accuracy=0, nonlinear=False)
    return pars


def ordinary_params(h0, *, lmax=LMAX, accuracy=1.0):
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=h0, ombh2=OMB, omch2=OMC_EARLY, mnu=MNU, tau=TAU)
    pars.InitPower.set_params(As=AS, ns=NS)
    pars.set_accuracy(
        AccuracyBoost=accuracy,
        lSampleBoost=accuracy,
        lAccuracyBoost=accuracy,
        DoLateRadTruncation=False,
    )
    pars.Accuracy.TimeStepBoost = accuracy
    pars.set_for_lmax(lmax, lens_potential_accuracy=0, nonlinear=False)
    return pars


def max_scaled_difference(actual, reference, *, floor_fraction=1e-4):
    """Compare a signed spectrum without inflating noise near its zero crossings."""
    scale = np.maximum(np.abs(reference), floor_fraction * np.max(np.abs(reference)))
    return np.max(np.abs(actual - reference) / scale)


def relative_rms_residual(residual, *terms):
    """Use the sizes of individual equation terms to scale a cancellation residual."""
    interior = slice(8, -8)

    def rms(values):
        return np.sqrt(np.mean(np.square(values[interior])))

    return rms(residual) / max(np.sqrt(sum(rms(term) ** 2 for term in terms)), 1e-30)


class IDECMBTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ide = camb.get_results(ide_params())
        cls.ordinary = camb.get_results(ordinary_params(cls.ide.Params.H0))
        cls.coupled = camb.get_results(ide_params(beta=0.08, lam=0.7, lmax=500))

    def test_lcdm_recombination_and_distances(self):
        z = np.array([0.1, 0.5, 2.0, 10.0, 100.0, 1100.0, 2500.0])
        np.testing.assert_allclose(
            self.ide.angular_diameter_distance(z),
            self.ordinary.angular_diameter_distance(z),
            rtol=3e-6,
        )
        np.testing.assert_allclose(self.ide.hubble_parameter(z), self.ordinary.hubble_parameter(z), rtol=3e-6)
        z_thermal = np.geomspace(1.0, 5000.0, 500) - 1.0
        ide = self.ide.get_background_redshift_evolution(z_thermal, vars=["x_e", "visibility"], format="array")
        ordinary = self.ordinary.get_background_redshift_evolution(
            z_thermal, vars=["x_e", "visibility"], format="array"
        )
        np.testing.assert_allclose(ide[:, 0], ordinary[:, 0], rtol=1e-3, atol=1e-6)
        visibility_scale = np.max(ordinary[:, 1])
        self.assertLess(np.max(np.abs(ide[:, 1] - ordinary[:, 1])) / visibility_scale, 1e-3)

    def test_lcdm_unlensed_lensed_and_lensing_potential(self):
        for getter in ("get_unlensed_scalar_cls", "get_lensed_scalar_cls"):
            with self.subTest(spectrum=getter):
                ide = getattr(self.ide, getter)(LMAX)[2:]
                ordinary = getattr(self.ordinary, getter)(LMAX)[2:]
                for column in (0, 1):  # TT, EE
                    self.assertLess(max_scaled_difference(ide[:, column], ordinary[:, column]), 2e-3)
                te_scale = np.sqrt(np.maximum(ordinary[:, 0] * ordinary[:, 1], 0.0))
                te_scale = np.maximum(te_scale, 1e-5 * np.max(te_scale))
                self.assertLess(np.max(np.abs(ide[:, 3] - ordinary[:, 3]) / te_scale), 2e-3)
        ide_pp = self.ide.get_lens_potential_cls(LMAX)[2:, 0]
        ordinary_pp = self.ordinary.get_lens_potential_cls(LMAX)[2:, 0]
        self.assertLess(max_scaled_difference(ide_pp, ordinary_pp), 2e-3)

    def test_coupled_spectra_are_finite_and_active(self):
        for getter in ("get_unlensed_scalar_cls", "get_lensed_scalar_cls"):
            spectra = getattr(self.coupled, getter)(500)[2:]
            self.assertTrue(np.all(np.isfinite(spectra)))
            self.assertTrue(np.all(spectra[:, 0] > 0))  # TT
            self.assertTrue(np.all(spectra[:, 1] > 0))  # EE
        pp = self.coupled.get_lens_potential_cls(500)[2:, 0]
        self.assertTrue(np.all(np.isfinite(pp)))
        self.assertTrue(np.all(pp > 0))
        ordinary_tt = self.ordinary.get_unlensed_scalar_cls(500)[2:, 0]
        coupled_tt = self.coupled.get_unlensed_scalar_cls(500)[2:, 0]
        self.assertGreater(max_scaled_difference(coupled_tt, ordinary_tt), 1e-3)

    def test_coupled_superhorizon_adiabatic_relations(self):
        # The first sampled node has a larger integration transient. Include
        # an earlier anchor and compare the later superhorizon nodes.
        z = np.array([1e8, 1e7, 5e6, 2.5e6])
        for beta in (-0.08, 0.08):
            with self.subTest(beta=beta):
                result = self.coupled if beta > 0 else camb.get_results(ide_params(beta=beta, lam=0.7, lmax=100))
                eta = result.conformal_time(z)
                k = 0.01
                self.assertLess(k * eta[-1], 0.01)
                evolution = result.get_time_evolution(
                    k, eta, vars=["delta_cdm", "delta_baryon", "delta_photon"], lAccuracyBoost=1
                )
                delta_c, delta_b, delta_g = evolution.T
                self.assertTrue(np.all(np.isfinite(evolution)))
                np.testing.assert_allclose(delta_b[1:], 0.75 * delta_g[1:], rtol=0.002, atol=1e-12)
                np.testing.assert_allclose(delta_c[1:], delta_b[1:], rtol=0.005, atol=1e-12)
                scaled = delta_g[1:] / eta[1:] ** 2
                self.assertLess(np.max(np.abs(scaled / scaled[0] - 1)), 0.005)

    def test_coupled_field_cdm_and_einstein_residuals(self):
        # These are independently differentiated mode outputs; the scalar and
        # CDM coupling terms must have opposite energy-transfer signs.
        k = 0.02
        z = np.r_[1e8, np.geomspace(5e4, 3e3, 300)]
        grho100 = 3 * (100 * 1000 / constants.c) ** 2
        names = [
            "a",
            "H",
            "delta_cdm",
            "etak",
            "delta_phi_ide",
            "delta_phi_prime_ide",
            "v_cdm_synch_ide",
            "z_synch_ide",
            "dgq_ide",
        ]
        for beta in (-0.08, 0.0, 0.08):
            with self.subTest(beta=beta):
                result = self.coupled if beta > 0 else camb.get_results(ide_params(beta=beta, lam=0.7, lmax=100))
                eta = result.conformal_time(z)
                values = result.get_time_evolution(k, eta, vars=names, lAccuracyBoost=1)[1:]
                times = eta[1:]
                a, hubble, delta_c, etak, field, field_prime, vc, z_synch, dgq = values.T
                phi = result.get_ide_history(a)["phi"]
                phi_prime = np.gradient(phi, times, edge_order=2)
                potential = grho100 * V0 * np.exp(-0.7 * phi)
                grhoc_t = grho100 * OMC_EARLY * np.exp(-beta * phi) / a

                continuity_terms = (
                    np.gradient(delta_c, times, edge_order=2),
                    k * (z_synch + vc),
                    beta * field_prime,
                )
                continuity = sum(continuity_terms)
                self.assertLess(relative_rms_residual(continuity, *continuity_terms), 0.001)

                euler_terms = (
                    np.gradient(vc, times, edge_order=2),
                    (hubble - beta * phi_prime) * vc,
                    beta * k * field,
                )
                euler = sum(euler_terms)
                self.assertLess(relative_rms_residual(euler, *euler_terms), 0.001)

                kg_terms = (
                    np.gradient(field_prime, times, edge_order=2),
                    2 * hubble * field_prime,
                    k * z_synch * phi_prime,
                    (k**2 + a**2 * 0.7**2 * potential) * field,
                    -beta * grhoc_t * delta_c,
                )
                kg = sum(kg_terms)
                self.assertLess(relative_rms_residual(kg, *kg_terms), 0.001)

                momentum_terms = (np.gradient(etak, times, edge_order=2), -0.5 * dgq)
                self.assertLess(relative_rms_residual(sum(momentum_terms), *momentum_terms), 0.001)

    def test_coupled_accuracy_convergence(self):
        # Fast low-l regression guard. Full Planck-range precision is assessed
        # separately by analysis/ide_accuracy_sweep.py and the likelihood pilot.
        tighter = camb.get_results(ide_params(beta=0.08, lam=0.7, lmax=500, accuracy=2.0))
        base = self.coupled.get_lensed_scalar_cls(500)[30:, :]
        refined = tighter.get_lensed_scalar_cls(500)[30:, :]
        for column in (0, 1):
            self.assertLess(max_scaled_difference(refined[:, column], base[:, column]), 0.02)
        te_scale = np.sqrt(np.maximum(base[:, 0] * base[:, 1], 0.0))
        self.assertLess(np.max(np.abs(refined[:, 3] - base[:, 3]) / te_scale), 0.02)

    def test_coupled_initial_epoch_and_history_convergence(self):
        earlier = camb.get_results(
            ide_params(beta=0.08, lam=0.7, lmax=500, a_init=1e-12, history_points=4097, integrate_tol=1e-11)
        )
        base_cls = self.coupled.get_lensed_scalar_cls(500)[30:]
        earlier_cls = earlier.get_lensed_scalar_cls(500)[30:]
        for column in (0, 1):
            self.assertLess(max_scaled_difference(earlier_cls[:, column], base_cls[:, column]), 0.005)
        te_scale = np.sqrt(np.maximum(base_cls[:, 0] * base_cls[:, 1], 0.0))
        self.assertLess(np.max(np.abs(earlier_cls[:, 3] - base_cls[:, 3]) / te_scale), 0.005)
        base_pp = self.coupled.get_lens_potential_cls(500)[30:, 0]
        earlier_pp = earlier.get_lens_potential_cls(500)[30:, 0]
        self.assertLess(max_scaled_difference(earlier_pp, base_pp), 0.005)

    def test_unsupported_modes_are_rejected(self):
        for name, mutate, error in (
            ("tensor", lambda p: setattr(p, "WantTensors", True), "tensor"),
            ("isocurvature", lambda p: setattr(p, "scalar_initial_condition", "initial_iso_CDM"), "adiabatic"),
            ("nonlinear", lambda p: setattr(p, "NonLinear", model.NonLinear_lens), "nonlinear"),
            ("21cm", lambda p: setattr(p, "Do21cm", True), "21 cm"),
            ("ionization perturbations", lambda p: setattr(p, "Evolve_delta_xe", True), "ionization"),
            ("late truncation", lambda p: setattr(p, "DoLateRadTruncation", True), "DoLateRadTruncation"),
        ):
            with self.subTest(mode=name):
                pars = ide_params(lmax=100)
                mutate(pars)
                with self.assertRaisesRegex(CAMBError, error):
                    camb.get_results(pars)
        with self.assertRaisesRegex(CAMBError, "nonlinear matter power"):
            self.coupled.get_nonlinear_matter_power_spectrum()
        with self.assertRaisesRegex(CAMBError, "WantTransfer"):
            self.coupled.get_matter_transfer_data()

    def test_late_initial_epoch_is_rejected_before_mode_evolution(self):
        pars = ide_params(beta=0.08, lam=0.7, lmax=2500, a_init=1e-7)
        with self.assertRaisesRegex(CAMBError, "initial scale factor is too late"):
            camb.get_results(pars)

    def test_oversized_scalar_hierarchy_returns_error(self):
        eta = self.coupled.conformal_time(1000.0)
        with self.assertRaisesRegex(CAMBError, "scalar hierarchy exceeds max_l_evolve"):
            self.coupled.get_time_evolution(0.02, [eta], lAccuracyBoost=50.0)


if __name__ == "__main__":
    unittest.main()
