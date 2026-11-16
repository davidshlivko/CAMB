Interacting dark energy
=======================

``IDEExponential`` evolves a dimensionless scalar field
:math:`\varphi=\phi/M_{\rm Pl}` with
:math:`V(\varphi)=V_0e^{-\lambda\varphi}`. The mass of every CDM particle
scales as :math:`e^{-\beta\varphi}`; baryons remain uncoupled. The physical
CDM density is

.. math::

   \frac{\rho_c(a)}{\rho_{\mathrm{crit},100}}
   =\omega_{c,\mathrm{early}}a^{-3}e^{-\beta\varphi(a)}.

The field starts at ``a_init`` (default :math:`10^{-10}`) with
:math:`\varphi=d\varphi/d\ln a=0`. Below ``a_init``, the field is frozen and
the CDM density follows its continuous :math:`a^{-3}` continuation. The model
assumes flat geometry and couples all CDM. ``V0_over_rho100`` and
``omc_early`` replace ``H0`` and ``omch2`` as inputs; CAMB derives the latter
two in its Fortran background initialization.

Use ``get_background(..., no_thermo=True)``. The input ``pars.H0`` and
``pars.omch2`` remain unresolved; read ``bg.Params.H0`` and
``bg.Params.omch2`` after calculation. The following example plots the field
and physical CDM density:

.. code-block:: python

   import camb
   import matplotlib.pyplot as plt
   import numpy as np

   pars = camb.CAMBparams()
   pars.set_ide_cosmology(
       V0_over_rho100=0.31,
       omc_early=0.12,
       beta=0.4,
       lam=0.8,
       ombh2=0.0224,
       mnu=0.06,
   )
   bg = camb.get_background(pars, no_thermo=True)
   print(f"H0 = {bg.Params.H0:.4f}, omch2 = {bg.Params.omch2:.6f}")

   a = np.geomspace(1e-8, 1, 400)
   history = bg.get_ide_history(a)
   fig, axes = plt.subplots(2, 1, sharex=True)
   axes[0].semilogx(history["a"], history["phi"])
   axes[0].set_ylabel(r"$\varphi$")
   axes[1].loglog(history["a"], history["rho_c_over_rho100"])
   axes[1].set_ylabel(r"$\rho_c/\rho_{\mathrm{crit},100}$")
   axes[1].set_xlabel("Scale factor a")
   plt.show()

``get_background_densities(a)["cdm"]`` also reflects the coupled density.
Its units are CAMB's :math:`8\pi G\rho_c a^4/c^2` in :math:`\mathrm{Mpc}^{-2}`;
``get_ide_history`` returns :math:`\rho_c/\rho_{\mathrm{crit},100}` instead.
Expansion, time, and distance methods are available. ``get_background(pars)``
also computes the RECFAST thermal and reionization history. For adiabatic
scalar modes, ``get_results(pars)`` computes unlensed and linearly lensed CMB
spectra, and ``get_transfer_functions(pars)`` provides intermediate CMB and
linear matter transfer functions used by Cobaya. For linear matter power and
:math:`\sigma_8`, call ``pars.set_matter_power(redshifts=[0], kmax=5, nonlinear=False)``
before ``get_results``. The default ``Transfer_tot`` includes
CDM, baryons, and massive neutrinos with their IDE background weights. Because
CDM moves in the IDE synchronous gauge, these transfer densities use that gauge
convention; their interpretation on very large scales needs care. Nonlinear
corrections, :math:`f\sigma_8`, source windows, tensors, vectors, isocurvature,
21 cm, and the approximate CosmoMC theta are not validated. The result's
``get_derived_params()`` includes the resolved ``H0`` and ``omch2``.

For example, a linear lensed CMB calculation can be configured through
``camb.set_params``:

.. code-block:: python

   import camb

   pars = camb.set_params(
       dark_energy_model="IDEExponential",
       V0_over_rho100=0.31,
       omc_early=0.12,
       beta=0.02,
       lam=0.3,
       ombh2=0.02237,
       mnu=0.06,
       tau=0.0544,
       As=2.10e-9,
       ns=0.9649,
       lmax=2500,
       lens_potential_accuracy=0,
       max_eta_k=18000,
   )
   results = camb.get_results(pars)
   cls = results.get_cmb_power_spectra(lmax=2500, CMB_unit="muK")
   print(results.get_derived_params()["H0"])

For a Cobaya run, add ``sigma8: {derived: true}`` to ``params``. This requests
linear matter transfers at :math:`z=0`; keep ``lens_potential_accuracy=0`` and
do not request nonlinear corrections.

The numerical settings must be checked for convergence at each coupled point
before interpreting a likelihood or chain. The pilot validation record is in
``analysis/ide_planck_pilot_results.md``.

The synchronous-gauge IDE mode diagnostics ``delta_phi_ide``,
``delta_phi_prime_ide``, ``v_cdm_synch_ide``, ``z_synch_ide``, ``dgq_ide``,
``dgrho_ide``, and ``sigma_synch_ide`` are available through
``CAMBdata.get_time_evolution``. The ``frame`` option for custom symbolic
sources has not been validated for moving CDM and is rejected for IDE.

Scalar perturbation convention
------------------------------

CDM particle number is conserved. Since
:math:`\delta\rho_c/\bar\rho_c=\delta n_c/\bar n_c-\beta\delta\varphi`, its
mass variation contributes to the density contrast even in synchronous
gauge. With CAMB's velocity convention :math:`\theta_c=k v_c` and its metric
variable :math:`z`, number conservation and the particle force law give

.. math::

   \delta_c'=-k(z+v_c)-\beta\delta\varphi',\qquad
   v_c'=-(\mathcal H-\beta\bar\varphi')v_c-\beta k\delta\varphi.

The field equation and its linear perturbation are

.. math::

   \bar\varphi''+2\mathcal H\bar\varphi'+a^2(8\pi G)V_{,\varphi}
   &=\beta(8\pi G)a^2\bar\rho_c,\\
   \delta\varphi''+2\mathcal H\delta\varphi'
   +(k^2+a^2(8\pi G)V_{,\varphi\varphi})\delta\varphi
   +kz\bar\varphi'
   &=\beta(8\pi G)a^2\bar\rho_c\delta_c.

CAMB adds the scalar density and momentum perturbations to the Einstein
sources, and includes :math:`\bar\rho_c v_c` in total momentum. The code
evolves :math:`v_c`, :math:`\delta\varphi`, and
:math:`\delta\varphi'` through every scalar hierarchy switch. Its initial
conditions set these three variables to zero at ``a_init`` while the standard
adiabatic growing mode supplies the other variables.
