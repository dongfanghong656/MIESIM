# Project Rules

This simulator follows the design documents' central separation:

1. Direct OCT truth is generated only by the OCT forward chain.
2. Microscope observations can initialize or constrain a converter, but cannot define OCT truth.
3. Any microscope-to-OCT formula must live outside `oct_forward.py` and be tested against direct output.

The current implementation targets a small but reproducible scalar baseline. It is intended to produce testable data products, not to claim a production-grade vectorial OCT solver.

## Implemented Validation Gates

- Airy sanity: clear circular pupil lateral FWHM is checked against the expected `0.514 lambda / NA` scale.
- Gaussian underfill sanity: smaller fill ratio broadens the lateral focus.
- Parseval sanity: unnormalized propagation preserves energy under orthonormal FFT.
- Dichroic sanity: non-absorbing scalar coefficients conserve power.
- OCT axial sanity: broader source bandwidth produces a narrower axial coherence gate.
- No-leak sanity: `oct_forward.py` must not import microscope or conversion modules.
- Spectral-error sanity: residual quadratic dispersion broadens the axial PSF.
- k-linearization sanity: configured pixel residuals perturb the measured k-grid while preserving monotonic sampling.
- Roll-off sanity: deeper scatterers have lower raw interferometric modulation when roll-off is enabled.

## Data Products

The validation runner writes a timestamped output folder under `outputs/validation_*` with:

- `microscope_forward/microscope_zstack.h5`
- `oct_direct/oct_raw_interferogram.h5`
- `oct_direct/oct_psf_direct.h5`
- `theory_conversion_under_test/path_e_converted_psf.h5`
- `metrics/validation_summary.json`
- `metrics/sweep_summary.csv`

These files are the first-stage test data package.

The sweep table currently covers:

- pupil fill ratio
- finite bead diameter
- visible/OCT focus shift
- residual quadratic dispersion
- k-linearization residual in pixel units
- depth roll-off

Rows marked `diagnostic` are OCT-processing perturbations rather than microscope-to-OCT conversion pass/fail judgments. Rows marked `review` indicate a baseline path E sensitivity that should be inspected before using the converted PSF as evidence.
