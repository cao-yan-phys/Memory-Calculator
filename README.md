# Vacuum Null Memory Calculator

Calculators for vacuum nonlinear-null gravitational-wave memory modes.

This version intentionally contains only:

- angular coupling coefficients;
- displacement, spin, and CM memory evaluators;
- leading-order PN helpers used for initial-offset checks;
- a bundled `lmax=10` gamma table for the modes used by the example;
- one example using `SEOBNRv5EHM` modes from `pyseobnr`;
- one `SEOBNRv5EHM`-vs-`NRHybSur3dq8_CCE` comparison example.

The core package has no `pyseobnr` or `gwsurrogate` dependency.  The examples do.
No local research artifacts, private notebooks, or generated JSON files are
included.  The only generated artifacts kept in the public tree are the small
reference CSV/PNG files produced by the examples.

The bundled gamma table lives at
`src/vacuum_memory_modes/data/gamma_coeffs_lmax10.npz`.  The public API loads
this file automatically for the example targets, so normal use does not rerun
the SymPy/Wigner-3j coefficient generation.  All valid azimuthal modes,
including `|m|=1`, are always included in the table.  If a different `lmax` or
target is requested, the code falls back to on-the-fly generation.

## Install

```bash
pip install -e .
```

For the `SEOBNRv5EHM` example, install `pyseobnr` in the same environment.
For the `SEOBNRv5EHM`-vs-`NRHybSur3dq8_CCE` example, install both `pyseobnr`
and `gwsurrogate`.

The waveform models used by the examples are:

- [`SEOBNRv5EHM` through `pyseobnr`](https://github.com/AEI-ACR/pyseobnr)
- [`NRHybSur3dq8_CCE` through `gwsurrogate`](https://github.com/sxs-collaboration/gwsurrogate)

## Example

```bash
python examples/eob_circular_memory_demo.py
```

The example generates a nonprecessing circular `SEOBNRv5EHM` event, computes
`h_20`, `h_30`, and the leading CM modes, infers an effective 0PN `x` from the
initial `dot h_20`, and compares the initial numerical memory modes with the
corresponding leading PN formulas.

The default start is `x ~= 0.015`, implemented as
`omega_start = 0.015**1.5`.  The saved PNG is a waveform comparison of
differences from the initial sample: `h_20` and `h_30` use the displayed
real/imaginary component of `h(t)-h(t0)`, while CM panels show
`|h_CM(t)-h_CM(t0)|`.  The plotted values are normalized by `nu M/R`,
i.e. the dimensionless mode amplitudes are divided by `nu`.  Numerical results
are black solid curves and the effective-0PN waveforms are red dashed curves.

For CM modes the example prints both:

- the full Nichols leading-PN value; and
- the leading-PN value truncated to the modes actually returned by pyEOB.

This matters because `SEOBNRv5EHM` does not provide every leading PN radiative
mode, e.g. it does not provide `(3,1)`.

A reference run is included:

- `examples/output/eob_circular_memory_q2_omega0.00183712.csv`
- `examples/output/eob_circular_memory_q2_omega0.00183712.png`

## SEOBNRv5EHM-Vs-NRHybSur3dq8_CCE Example

```bash
python examples/eob_cce_h20_h30_comparison.py
```

This example loads `NRHybSur3dq8_CCE`, finds the `NRHybSur3dq8_CCE` time where
`Omega_orb ~= 0.015**1.5`, fits the initial `NRHybSur3dq8_CCE` `(2,2)` phase
to get the `SEOBNRv5EHM` `omega_start`, and then compares `h(t)-h(t0)` for
the `NRHybSur3dq8_CCE` `(2,0)` and `(3,0)` modes against `SEOBNRv5EHM`
perturbative `h_20` and spin-memory `h_30`.  The plot covers the full
`NRHybSur3dq8_CCE` time range through the final sample.  The legend uses the
concrete model names `NRHybSur3dq8_CCE` and `SEOBNRv5EHM`, the y axis uses a
positive logarithmic scale for the displayed components, and the plotted
components are normalized by `nu M/R`.

To use a different oscillatory input model, replace
`_generate_pyseobnr_positive_modes` or call `_compute_memory_from_positive_modes`
with another nonprecessing positive-`m` mode dictionary.  The memory calculation
does not depend on `SEOBNRv5EHM` specifically once those modes are supplied.

The reference output is:

- `examples/output/eob_cce_h20_h30_q2_x0.015.csv`
- `examples/output/eob_cce_h20_h30_q2_x0.015.png`
