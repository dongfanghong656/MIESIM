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
