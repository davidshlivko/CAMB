"""Physical and API checks for the background-only exponential IDE model."""

import unittest

import numpy as np
from scipy.integrate import solve_ivp

import camb
from camb import constants, model
from camb.baseconfig import CAMBError, CAMBValueError

V0 = 0.31
OMC = 0.12
OMB = 0.0224


def ide_background(*, beta=0.0, lam=0.0, mnu=0.0, **kwargs):
    pars = camb.CAMBparams()
    pars.set_ide_cosmology(
        V0_over_rho100=V0,
        omc_early=OMC,
        beta=beta,
        lam=lam,
        ombh2=OMB,
        mnu=mnu,
        **kwargs,
    )
    return pars, camb.get_background(pars, no_thermo=True)


class IDEBackgroundTest(unittest.TestCase):
    def test_lcdm_limit_and_massive_neutrinos(self):
        for mnu in (0.0, 0.06):
            with self.subTest(mnu=mnu):
                pars, bg = ide_background(mnu=mnu)
                self.assertTrue(np.isnan(pars.H0))
                self.assertTrue(np.isnan(pars.omch2))
                self.assertAlmostEqual(bg.Params.omch2, OMC, delta=2e-11)
                self.assertGreater(bg.Params.H0, 0)
                ordinary = camb.CAMBparams()
                ordinary.set_cosmology(H0=bg.Params.H0, ombh2=OMB, omch2=OMC, mnu=mnu)
                lcdm = camb.get_background(ordinary, no_thermo=True)
                z = np.array([0.0, 0.3, 1.0, 10.0, 1100.0])
                np.testing.assert_allclose(bg.hubble_parameter(z), lcdm.hubble_parameter(z), rtol=2e-7)
                np.testing.assert_allclose(
                    bg.angular_diameter_distance(z), lcdm.angular_diameter_distance(z), rtol=2e-7, atol=1e-9
                )
                self.assertAlmostEqual(bg.sound_horizon(1100.0), lcdm.sound_horizon(1100.0), delta=2e-5)
                self.assertAlmostEqual(bg.physical_time_a1_a2(0.0, 1.0), lcdm.physical_time_a1_a2(0.0, 1.0), delta=2e-5)
                a = np.geomspace(1e-9, 1.0, 20)
                history = bg.get_ide_history(a)
                np.testing.assert_allclose(history["phi"], 0.0, atol=2e-10)
                np.testing.assert_allclose(history["rho_c_over_rho100"] * a**3, OMC, rtol=2e-10)

    def test_uncoupled_scalar_against_independent_integration(self):
        lam = 0.9
        pars, bg = ide_background(lam=lam, mnu=0.0)
        ai = pars.DarkEnergy.a_init
        grho100 = 3 * (100 * 1000 / constants.c) ** 2
        omega_r = (bg.grhog + bg.grhornomass) / grho100

        def rhs(x, y):
            phi, u = y
            a = np.exp(x)
            matter = (OMB + OMC) / a**3
            radiation = omega_r / a**4
            potential = V0 * np.exp(-lam * phi)
            e2 = (matter + radiation + potential) / (1 - u * u / 6)
            dlnh = -0.5 * ((3 * matter + 4 * radiation) / e2 + u * u)
            return [u, -(3 + dlnh) * u + 3 * lam * potential / e2]

        a = np.geomspace(1e-6, 1.0, 18)
        ref = solve_ivp(
            rhs,
            (np.log(ai), 0.0),
            (0.0, 0.0),
            t_eval=np.log(a),
            method="DOP853",
            rtol=2e-11,
            atol=2e-13,
        )
        self.assertTrue(ref.success, ref.message)
        phi, u = ref.y
        e2 = ((OMB + OMC) / a**3 + omega_r / a**4 + V0 * np.exp(-lam * phi)) / (1 - u * u / 6)
        np.testing.assert_allclose(bg.get_ide_history(a)["phi"], phi, rtol=3e-5, atol=2e-6)
        np.testing.assert_allclose(bg.hubble_parameter(1 / a - 1) / 100, np.sqrt(e2), rtol=3e-5)

    def test_coupled_continuity_densities_and_closure(self):
        beta, lam = 0.4, 0.8
        _, bg = ide_background(beta=beta, lam=lam, mnu=0.0)
        a = np.geomspace(1e-8, 1.0, 400)
        history = bg.get_ide_history(a)
        rho_c = history["rho_c_over_rho100"]
        phi = history["phi"]
        np.testing.assert_allclose(rho_c * a**3 * np.exp(beta * phi), OMC, rtol=2e-10)
        grho100 = 3 * (100 * 1000 / constants.c) ** 2
        densities = bg.get_background_densities(a)
        np.testing.assert_allclose(densities["cdm"], grho100 * rho_c * a**4, rtol=2e-9)
        np.testing.assert_allclose(densities["baryon"], bg.grhob * a, rtol=2e-12)
        np.testing.assert_allclose(
            densities["tot"],
            sum(densities[name] for name in model.density_names if name != "tot"),
            rtol=2e-12,
        )
        self.assertAlmostEqual(bg.Params.omch2, rho_c[-1], delta=2e-10)
        self.assertAlmostEqual(bg.hubble_parameter(0), bg.Params.H0, delta=2e-6)
        self.assertAlmostEqual(densities["tot"][-1] / grho100, (bg.Params.H0 / 100) ** 2, delta=2e-7)

        # Total energy conservation uses the scalar kinetic term, with the CDM
        # transfer cancelling between the two species.
        ln_a = np.log(a)
        scalar_rho = densities["de"] / (grho100 * a**4)
        u = np.gradient(phi, ln_a, edge_order=2)
        dscalar = np.gradient(scalar_rho, ln_a, edge_order=2)
        e2 = (bg.hubble_parameter(1 / a - 1) / 100) ** 2
        residual = dscalar + e2 * u**2 - beta * u * rho_c
        region = slice(20, -20)
        scale = np.maximum(np.abs(dscalar[region]) + np.abs(beta * u[region] * rho_c[region]), 1e-7)
        self.assertLess(np.max(np.abs(residual[region]) / scale), 0.03)

        rho_total = densities["tot"] / (grho100 * a**4)
        rho_r = (bg.grhog + bg.grhornomass) / (grho100 * a**4)
        rho_b = bg.grhob / (grho100 * a**3)
        energy_residual = np.gradient(rho_total, ln_a, edge_order=2) + 4 * rho_r + 3 * rho_b + 3 * rho_c + e2 * u**2
        energy_scale = 4 * rho_r + 3 * rho_b + 3 * rho_c + e2 * u**2
        self.assertLess(np.max(np.abs(energy_residual[region]) / energy_scale[region]), 0.02)

    def test_initial_time_and_accuracy_controls(self):
        a = np.array([1e-6, 1e-3, 0.1, 1.0])
        _, baseline = ide_background(beta=0.4, lam=0.8, mnu=0.0)
        _, earlier = ide_background(beta=0.4, lam=0.8, mnu=0.0, a_init=1e-12)
        _, tighter = ide_background(beta=0.4, lam=0.8, mnu=0.0, integrate_tol=1e-11, history_points=4097)
        for other in (earlier, tighter):
            np.testing.assert_allclose(
                baseline.get_ide_history(a)["phi"], other.get_ide_history(a)["phi"], rtol=2e-4, atol=2e-6
            )
            np.testing.assert_allclose(
                baseline.hubble_parameter(1 / a - 1), other.hubble_parameter(1 / a - 1), rtol=2e-5
            )
        below = baseline.get_ide_history(np.array([1e-12, 1e-11, 1e-10]))
        np.testing.assert_allclose(below["phi"], 0, atol=0)
        np.testing.assert_allclose(below["rho_c_over_rho100"] * below["a"] ** 3, OMC, rtol=1e-12)

    def test_unsupported_paths_and_bad_inputs(self):
        pars, bg = ide_background(beta=0.2, lam=0.5)
        thermo = camb.get_background(pars)
        self.assertTrue(np.isfinite(thermo.get_derived_params()["zstar"]))
        self.assertEqual(thermo.get_derived_params()["H0"], thermo.Params.H0)
        self.assertEqual(thermo.get_derived_params()["omch2"], thermo.Params.omch2)
        with self.assertRaisesRegex(CAMBError, "theta"):
            bg.cosmomc_theta()
        with self.assertRaises(CAMBValueError):
            bg.get_ide_history([0.0, 1.0])
        with self.assertRaises(CAMBValueError):
            camb.CAMBparams().set_ide_cosmology(V0_over_rho100=-1, omc_early=OMC, beta=0, lam=0)
        pars.NonLinear = model.NonLinear_pk
        with self.assertRaisesRegex(CAMBError, "nonlinear"):
            camb.get_background(pars, no_thermo=True)

    def test_set_params_discovery_and_partial_ide_update(self):
        names = camb.get_valid_numerical_params(dark_energy_model="IDEExponential")
        self.assertTrue({"V0_over_rho100", "omc_early", "beta", "lam", "tau"} <= names)
        self.assertFalse({"H0", "omch2", "cosmomc_theta", "thetastar"} & names)
        pars = camb.set_params(
            dark_energy_model="IDEExponential",
            V0_over_rho100=V0,
            omc_early=OMC,
            beta=0.01,
            lam=0.5,
            ombh2=OMB,
            mnu=0.1,
            nnu=3.046,
            history_points=3001,
            tau=0.0544,
        )
        updated = camb.set_params(pars.copy(), beta=0.02)
        self.assertEqual(updated.DarkEnergy.beta, 0.02)
        self.assertEqual(updated.DarkEnergy.history_points, 3001)
        self.assertEqual(updated.ombh2, pars.ombh2)
        self.assertEqual(updated.omnuh2, pars.omnuh2)
        self.assertEqual(updated.Reion.optical_depth, pars.Reion.optical_depth)
        with self.assertRaisesRegex(CAMBError, "physical inputs requires"):
            camb.set_params(pars.copy(), ombh2=0.023)
        with self.assertRaisesRegex(CAMBError, "derives"):
            camb.set_params(
                dark_energy_model="IDEExponential",
                V0_over_rho100=V0,
                omc_early=OMC,
                beta=0.0,
                lam=0.0,
                H0=67.0,
            )


if __name__ == "__main__":
    unittest.main()
