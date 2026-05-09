# Physics Contract

This project is a scalar low-NA validation harness for common-path microscope/OCT model development. It is not a final physical truth simulator.

## Units

- Spatial coordinates are in micrometers.
- Wavelengths in config files are in nanometers.
- OCT k values are in radians per micrometer.
- Pupil coordinates use normalized radius rho in the objective pupil.

## Defocus

- `objective.defocus_um` is a Zernike defocus OPD coefficient in micrometers.
- Physical sample or focus displacement is not a Zernike OPD coefficient.
- Physical z displacement enters propagation through the low-NA physical defocus phase:

```text
phase(rho) = pi * z_um * NA^2 * rho^2 / wavelength_um
```

- `errors.objective_focus_shift_um` and `objective.focus_shift_um_850_vs_visible` are physical focus shifts.

## Direct OCT Models

- `hybrid_rci` is a fast diagnostic approximation. It combines a spectral axial gate with a center-wavelength through-focus RCI stack.
- `full_spectral_rci` is the preferred scalar low-NA review path. It inserts `h_RCI(x,y,z;k)` before k-linearization, windowing, dispersion phase, and FFT reconstruction.
- The full-spectral review artifact stores sampled complex `h_RCI[depth,k,y,x]`, true/measured/linearized k grids, wavelength grid, source spectrum, window vector, dispersion phase, fixed scatterer-depth phase origin, and normalized max/on-axis/integrated axial profiles.
- Neither path is vector Debye, Zemax POP, FDTD, or an internal Mie-scattering truth model.

## Common-Path Dispersion

- Common-path OCT does not retain a physically meaningful absolute quadratic spectral phase when the sample and reference arms share the same dispersive path.
- The modeled reconstruction phase is therefore the mismatch term:

```text
phi_diff(k) = (errors.dispersion_quadratic_rad - errors.reference_dispersion_quadratic_rad) * xi(k)^2
xi(k) = 2 * (k - mean(k)) / (max(k) - min(k))
```

- `errors.dispersion_quadratic_rad` is the sample-arm quadratic coefficient.
- `errors.reference_dispersion_quadratic_rad` is the reference-arm quadratic coefficient.
- If the two coefficients are equal, the common-path quadratic dispersion contribution cancels to numerical zero.
- `absolute_dispersion_phase(...)` exists only for diagnostics; direct OCT field assembly uses `differential_dispersion_phase(...)`.

## Sensitivity Roll-Off

- The scalar OCT forward model applies a finite-spectrometer sensitivity roll-off before reconstruction.
- The baseline theoretical slope follows the practical Leitgeb/Yun-style convention used in the Round 5 review plan: approximately 6.7 dB loss at half the k-sampling Nyquist imaging range.
- The modeled amplitude attenuation is:

```text
A(z) = exp(-|z_um| * alpha_total)
alpha_total = db_per_mm_to_amplitude_rolloff_per_um(theoretical_sensitivity_rolloff_db_per_mm(config))
              + max(errors.rolloff_per_um, 0)
```

- `errors.rolloff_per_um` is therefore an additional empirical error term, not the sole spectrometer roll-off source.
- `checks["sensitivity_rolloff_v_gate"]` fits the roll-off slope from forward-model `oct_raw` roll-off amplitudes and compares it with `theoretical_sensitivity_rolloff_db_per_mm(config)`.

## Phase Stability

- `checks["phase_stability_v_gate"]` is measured through the OCT forward and reconstruction path, not by injecting phase noise directly.
- The validation path simulates a point reflector with `simulate_oct_raw_direct`, adds Gaussian read noise using `errors.camera_read_noise_e / errors.photon_gain_e_per_adu`, reconstructs with `reconstruct_sd_oct`, and reports the phase standard deviation at the reflector depth over repeated draws.

## Gate Semantics

- `checks.all_pass` means every implemented machine gate passed.
- `verdict.pilot_pass` means `checks.all_pass=True` and `blocker_count=0`.
- `verdict.review_required=True` means the package must not be interpreted as a final physics validation result.
- Unsupported physics claims are counted as blockers until implemented or explicitly downgraded to warnings.

## FWHM Metrics

- Global max-over-XY axial FWHM is retained for backwards comparison.
- Central-lobe FWHM is used where mirror peaks or sidelobes could merge multiple lobes.
- On-axis and integrated axial FWHM are reported to expose focus drift and off-axis artifacts.

## Calibration Paths

Relative calibration paths in YAML configs are resolved at load time. The loader first tries paths relative to the config file directory, then the project root. Resolved configs write absolute paths so review packages can be reproduced from a clean working directory.
