# Design Summary

The source documents define a common-path SD-OCT system with a visible microscope branch and an 830-850 nm OCT branch sharing the objective/sample-side optics.

The implemented baseline keeps the following hierarchy:

- shared pupil `Pc`
- microscope channel pupil `Pm = Qm * Dm * Pc`
- OCT illumination and detection pupils `Pi`, `Pd`
- OCT two-way spatial kernel `h_RCI = ui * conj(ud)`
- effective OCT PSF from spectral weighting and FFT over `k`

Path E is the first converter under test:

`microscope spatial term M0` plus `OCT axial term gamma0` gives a separable baseline `h_pred = M0 * gamma0`.

This is intentionally only a baseline. Path B, C, A, and G can later replace `M0` with fitted pupil, reflectance/confocal, complex pupil, or joint inversion estimates while leaving direct OCT truth isolated.

## Current Error Model

The direct OCT chain now includes the first processing-error hooks from the design package:

- residual quadratic dispersion phase in k-space
- measured k-grid perturbation expressed as RMS pixel residual
- depth roll-off attenuation for off-zero scatterers
- objective focus shift and dichroic tilt as pupil-level perturbations

These are deliberately low-dimensional first-stage perturbations. They are meant to expose failure fingerprints and test boundaries before upgrading to full Jones AOI tables, vector Debye propagation, or imported POP/FDTD fields.
