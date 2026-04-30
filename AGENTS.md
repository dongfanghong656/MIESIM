# Common-Path OCT Simulation Rules

This project is the executable companion to the Pro design package in:

- `C:\Users\1\Downloads\Common_Path_OCT_Simulation_Design_Deliverables`
- `C:\Users\1\Desktop\mic_to_oct_paths_rigorous_v2`

## Non-Negotiable Physics Boundary

- `cop_oct_sim/oct_forward.py` is the direct OCT truth chain.
- Direct OCT truth must be generated from shared pupil, OCT illumination/detection pupils, source spectrum, reference field, sample scatterer, spectrometer sampling, and reconstruction.
- Direct OCT code must not import `microscope_forward`, `theory_conversion`, or any microscope-to-OCT converter.
- Microscope paths A-G are measurement or constraint generators. They are not truth generators for OCT.

## Minimum Evidence Rule

Every validation run must save:

- resolved configuration and random seed
- direct OCT raw interferogram and reconstructed/direct PSF data
- microscope ideal and measured z-stack data
- converted or baseline PSF data when a converter is tested
- metrics as CSV/JSON, including pass/fail flags

## Numerical Guardrails

- Preallocate arrays in simulation loops.
- Keep large image stacks in `float32` or `complex64` unless a calculation needs `float64`.
- Do not use display-only `fftshift` inside core numerical comparisons unless the axis convention is explicitly documented.
- Keep random seeds fixed for validation data.
- Prefer single-factor sweeps before Monte Carlo combinations.

## First-Stage Scope

The first implemented stage is a rigorous low-NA scalar baseline:

- shared pupil with circular aperture, Gaussian underfill, and defocus
- scalar dichroic/Jones-compatible operator
- Fraunhofer/defocus-stack propagation
- widefield/confocal microscope z-stack generator
- independent direct OCT raw and effective PSF generator
- path E separable baseline converter as the first model under test

Vector Debye, Zemax POP import, full Jones pupil-dependent AOI, and FDTD checks are extension points, not silently implied completed features.
