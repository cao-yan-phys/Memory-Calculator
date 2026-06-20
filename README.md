# Vacuum Null Memory Calculator

Perturbative calculators for vacuum nonlinear-null gravitational-wave memory modes.

Includes:

- angular coupling coefficients;
- displacement, spin, and CM memory evaluators;
- leading-order PN helpers for nonprecessing quasicircular compact binaries;
- a bundled `lmax=10` gamma table for common memory-mode targets.

The bundled gamma table lives at `src/vacuum_memory_modes/data/gamma_coeffs_lmax10.npz`. The package loads this file automatically for the example targets, so normal use does not rerun the SymPy/Wigner-3j coefficient generation. If a different `lmax` or target is requested, the code falls back to on-the-fly generation.

## Install

```bash
pip install -e .
```

For the `SEOBNRv5EHM` example, install `pyseobnr` in the same environment. For the `SEOBNRv5EHM`-vs-`NRHybSur3dq8_CCE` example, install both `pyseobnr` and `gwsurrogate`. For the `FastEMRIWaveforms` example, install `fastemriwaveforms`.

The waveform models used by the examples are:

- [`SEOBNRv5EHM` through `pyseobnr`](https://github.com/AEI-ACR/pyseobnr)
- [`NRHybSur3dq8_CCE` through `gwsurrogate`](https://github.com/sxs-collaboration/gwsurrogate)
- [`FastEMRIWaveforms`](https://github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms)

## Example

```bash
python examples/seobnrv5ehm_circular_memory_demo.py
```

The example generates a nonprecessing circular `SEOBNRv5EHM` event, computes $h_{20}$, $h_{30}$, and the leading CM-memory modes, infers an effective 0PN $x$ parameter from the initial $\dot h_{20}$, and compares the initial numerical memory modes with the corresponding leading PN formulas.

The default initial PN parameter is $x_0=0.015$, implemented as `omega_start = 0.015**1.5`. $h_{20}$ and $h_{30}$ panels show the real/imaginary component of $h(t)-h(t_0)$, while CM-memory panels show $|h_{\rm CM}(t)-h_{\rm CM}(t_0)|$.

For CM-memory modes the example prints both:

- the full Nichols leading-PN value; and
- the leading-PN value truncated to the modes actually returned by `pyseobnr`.

This matters because `SEOBNRv5EHM` does not provide every leading PN radiative mode, e.g. it does not provide the $(3,1)$ mode.

Example outputs:

- `examples/output/seobnrv5ehm_circular_memory_q2_omega0.00183712.csv`
- `examples/output/seobnrv5ehm_circular_memory_q2_omega0.00183712.png`

## SEOBNRv5EHM-Vs-NRHybSur3dq8_CCE Example

```bash
python examples/seobnrv5ehm_nrhybsur3dq8_cce_h20_h30_comparison.py
```

This example loads `NRHybSur3dq8_CCE`, finds the `NRHybSur3dq8_CCE` time where $\Omega_{\rm orb}=0.015^{3/2}$, fits the initial `NRHybSur3dq8_CCE` $(2,2)$ phase to get the `SEOBNRv5EHM` `omega_start`, and then compares $h(t)-h(t_0)$ for the `NRHybSur3dq8_CCE` $(2,0)$ and $(3,0)$ modes against `SEOBNRv5EHM` perturbative $h_{20}$ and $h_{30}$.

To use a different oscillatory input model, replace `_generate_pyseobnr_positive_modes` or call `_compute_memory_from_positive_modes` with another nonprecessing positive-$m$ mode dictionary.

Example outputs:

- `examples/output/seobnrv5ehm_nrhybsur3dq8_cce_h20_h30_q2_x0.015.csv`
- `examples/output/seobnrv5ehm_nrhybsur3dq8_cce_h20_h30_q2_x0.015.png`

## FastEMRIWaveforms Example

```bash
python examples/fastemriwaveforms_emri_h20_h30_demo.py
```

This example uses `FastEMRIWaveforms` to generate a circular equatorial Kerr trajectory with mass ratio $q=10^5$ and spin $\chi=0.8$, computes perturbative $h_{20}$ and $h_{30}$ from the oscillatory modes, and compares $h(t)-h(t_0)$ with the same effective-0PN construction used in the `SEOBNRv5EHM` circular example.

Example outputs:

- `examples/output/fastemriwaveforms_emri_h20_h30_q100000_p0_100_chi0p8.csv`
- `examples/output/fastemriwaveforms_emri_h20_h30_q100000_p0_100_chi0p8.png`

## References

- Marc Favata, "Post-Newtonian corrections to the gravitational-wave memory for quasicircular, inspiralling compact binaries", [arXiv:0812.0069](https://arxiv.org/abs/0812.0069).
- David A. Nichols, "Spin memory effect for compact binaries in the post-Newtonian approximation", [arXiv:1702.03300](https://arxiv.org/abs/1702.03300).
- David A. Nichols, "Center-of-mass angular momentum and memory effect in asymptotically flat spacetimes", [arXiv:1807.08767](https://arxiv.org/abs/1807.08767).
