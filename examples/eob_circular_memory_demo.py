"""Compute h20, h30, and CM memory from a circular nonprecessing EOB event."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

MTSUN_SI = 4.925490947641266978197229498498379006e-6
DEFAULT_X_START = 0.015
DEFAULT_OMEGA_START = DEFAULT_X_START**1.5

from vacuum_memory_modes import (  # noqa: E402
    cm_strain_lo_modes,
    complete_nonprecessing_modes,
    compute_memory_modes,
    cumulative_integral,
    differentiate_modes,
    h20_lo,
    h30_spin_lo,
    infer_x_eff_from_dh20,
    k20_lo,
    phase_from_h22_lo,
    symmetric_mass_ratio,
)


def _relative_error(value: complex, reference: complex) -> float:
    reference_abs = abs(reference)
    return float(abs(value - reference) / reference_abs) if reference_abs else np.nan


def _format_complex(value: complex) -> str:
    value = complex(value)
    return f"{value.real:+.6e}{value.imag:+.6e}j"


def _x_0pn_series(t: np.ndarray, x0: float, q: float) -> np.ndarray:
    """Evolve ``x`` with ``dx/dt = (64 nu / 5) x^5`` from the first sample."""

    nu = symmetric_mass_ratio(q)
    denominator = float(x0) ** -4 - (256.0 / 5.0) * nu * (t - t[0])
    if np.any(denominator <= 0.0):
        raise ValueError("0PN x evolution reached its formal coalescence before the plot end")
    return denominator ** -0.25


def _plot_indices(t: np.ndarray, duration: float | None, max_points: int = 2500) -> np.ndarray:
    if duration is None or duration <= 0.0:
        end = len(t)
    else:
        end = int(np.searchsorted(t, t[0] + duration, side="right"))
        end = max(2, min(end, len(t)))
    step = max(1, int(np.ceil(end / max_points)))
    return np.arange(0, end, step, dtype=int)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--q", type=float, default=2.0)
    parser.add_argument("--omega-start", type=float, default=DEFAULT_OMEGA_START)
    parser.add_argument("--lmax", type=int, default=10)
    parser.add_argument("--delta-t", type=float, default=20.0, help="output spacing in units of M")
    parser.add_argument("--total-mass-solar", type=float, default=50.0)
    parser.add_argument("--output-dir", default=str(ROOT / "examples" / "output"))
    parser.add_argument("--plot-duration", type=float, default=1_250_000.0)
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()

    try:
        from pyseobnr.generate_waveform import generate_modes_opt
    except ImportError as exc:
        raise SystemExit("This example requires pyseobnr. Install the optional EOB dependency.") from exc

    t, raw_modes = generate_modes_opt(
        args.q,
        0.0,
        0.0,
        args.omega_start,
        eccentricity=0.0,
        approximant="SEOBNRv5EHM",
        settings={
            "EccIC": 0,
            "M": args.total_mass_solar,
            "dt": args.delta_t * args.total_mass_solar * MTSUN_SI,
            "lmax_nyquist": 1,
        },
    )
    eob_positive_modes = {tuple(map(int, key.split(","))): value for key, value in raw_modes.items()}
    oscillatory_modes = complete_nonprecessing_modes(eob_positive_modes)
    oscillatory_hdot = differentiate_modes(t, oscillatory_modes)

    targets = [(2, 0), (3, 0), (3, 1), (3, 3), (5, 1), (5, 3), (5, 5)]
    primary = compute_memory_modes(t, oscillatory_modes, targets, lmax=args.lmax, hdot=oscillatory_hdot)

    h20 = primary[(2, 0)]["h_displacement"]
    dh20_dt = primary[(2, 0)]["dh_displacement_dt"]
    x_eff = infer_x_eff_from_dh20(dh20_dt[0], args.q)
    h20_offset = h20_lo(args.q, x_eff)
    h20_absolute = h20 - h20[0] + h20_offset

    modes_with_h20 = dict(oscillatory_modes)
    modes_with_h20[(2, 0)] = h20_absolute
    hdot_with_h20 = dict(oscillatory_hdot)
    hdot_with_h20[(2, 0)] = dh20_dt
    with_h20 = compute_memory_modes(t, modes_with_h20, targets, lmax=args.lmax, hdot=hdot_with_h20)

    phase0 = phase_from_h22_lo(oscillatory_modes[(2, 2)][0])
    available_pn_modes = set(modes_with_h20)
    cm_lo_full = cm_strain_lo_modes(args.q, x_eff, phase0)
    cm_lo_available = cm_strain_lo_modes(args.q, x_eff, phase0, available_modes=available_pn_modes)

    h30_num = primary[(3, 0)]["h_spin_mode"][0]
    h30_lo = h30_spin_lo(args.q, x_eff)

    print("EOB circular nonprecessing memory demo")
    print(f"q = {args.q:g}, omega_start = {args.omega_start:g}")
    print(f"samples = {len(t)}, time range = [{t[0]:.3f}, {t[-1]:.3f}] M")
    print(f"pyEOB positive-m modes = {sorted(eob_positive_modes)}")
    print()
    print("Effective 0PN initial x from dot h20")
    print(f"Re(dot h20)_0 = {np.real(dh20_dt[0]):.12e}")
    print(f"x_eff = {x_eff:.12e}")
    print(f"h20_LO(x_eff) = {h20_offset:.12e}")
    print(f"h20_abs_numeric(t0) = {h20_absolute[0].real:.12e}")
    print()
    print("Initial spin-memory h30")
    print("mode        numeric                 LO PN                 rel.err")
    print(
        f"h30   {_format_complex(h30_num):>24}  {_format_complex(h30_lo):>24}"
        f"  {_relative_error(h30_num, h30_lo):.3e}"
    )
    print()
    print("Initial CM strain modes")
    print("mode        numeric                 LO PN full            rel.err    LO PN available       rel.err")
    rows = []
    for target in [(3, 1), (3, 3), (5, 1), (5, 3), (5, 5)]:
        numeric = with_h20[target]["h_cm_mode"][0]
        full = cm_lo_full[target]
        available = cm_lo_available[target]
        rows.append(
            {
                "mode": str(target),
                "numeric": numeric,
                "lo_full": full,
                "lo_available": available,
                "relerr_full": _relative_error(numeric, full),
                "relerr_available": _relative_error(numeric, available),
            }
        )
        print(
            f"{target!s:<6} {_format_complex(numeric):>24}  {_format_complex(full):>24}"
            f"  {_relative_error(numeric, full):.3e}"
            f"  {_format_complex(available):>24}  {_relative_error(numeric, available):.3e}"
        )
    print()
    print("Note: 'LO PN available' omits PN radiative moments that pyEOB did not return.")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f"eob_circular_memory_q{args.q:g}_omega{args.omega_start:g}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "quantity",
                "numeric_real",
                "numeric_imag",
                "lo_full_real",
                "lo_full_imag",
                "relerr_full",
                "lo_available_real",
                "lo_available_imag",
                "relerr_available",
            ]
        )
        writer.writerow(
            [
                "h30",
                h30_num.real,
                h30_num.imag,
                h30_lo.real,
                h30_lo.imag,
                _relative_error(h30_num, h30_lo),
                h30_lo.real,
                h30_lo.imag,
                _relative_error(h30_num, h30_lo),
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row["mode"],
                    row["numeric"].real,
                    row["numeric"].imag,
                    row["lo_full"].real,
                    row["lo_full"].imag,
                    row["relerr_full"],
                    row["lo_available"].real,
                    row["lo_available"].imag,
                    row["relerr_available"],
                ]
            )
    print(f"Saved comparison table: {csv_path}")

    if not args.no_plot:
        import matplotlib.pyplot as plt

        x_0pn = _x_0pn_series(t, x_eff, args.q)
        phase_0pn = phase0 + cumulative_integral(t, x_0pn**1.5)
        indices = _plot_indices(t, args.plot_duration)
        t_plot = t[indices] - t[0]
        x_plot = x_0pn[indices]
        phase_plot = phase_0pn[indices]

        h20_0pn = k20_lo(args.q) * x_plot
        h30_0pn = np.array([h30_spin_lo(args.q, x_value) for x_value in x_plot])
        cm_0pn = {
            target: np.array(
                [
                    cm_strain_lo_modes(args.q, x_value, phase_value)[target]
                    for x_value, phase_value in zip(x_plot, phase_plot)
                ]
            )
            for target in [(3, 1), (3, 3), (5, 1), (5, 3), (5, 5)]
        }
        nu = symmetric_mass_ratio(args.q)

        plot_specs = [
            (
                r"Re $\Delta h_{20}/(\nu M/R)$",
                np.real(h20_absolute[indices] - h20_absolute[0]) / nu,
                np.real(h20_0pn - h20_0pn[0]) / nu,
            ),
            (
                r"Im $\Delta h_{30}/(\nu M/R)$",
                np.imag(primary[(3, 0)]["h_spin_mode"][indices] - primary[(3, 0)]["h_spin_mode"][0]) / nu,
                np.imag(h30_0pn - h30_0pn[0]) / nu,
            ),
            *[
                (
                    rf"$|\Delta h_{{{target[0]}{target[1]}}}^{{CM}}|/(\nu M/R)$",
                    np.abs(with_h20[target]["h_cm_mode"][indices] - with_h20[target]["h_cm_mode"][0]) / nu,
                    np.abs(cm_0pn[target] - cm_0pn[target][0]) / nu,
                )
                for target in [(3, 1), (3, 3), (5, 1), (5, 3), (5, 5)]
            ],
        ]

        fig, axes = plt.subplots(4, 2, figsize=(11, 10), sharex=True, constrained_layout=True)
        flat_axes = axes.ravel()
        for ax, (label, numeric, effective_0pn) in zip(flat_axes, plot_specs):
            ax.plot(t_plot, numeric, color="black", linewidth=1.4, label="numeric")
            ax.plot(t_plot, effective_0pn, color="red", linestyle="--", linewidth=1.3, label="effective 0PN")
            ax.set_ylabel(label)
            ax.grid(True, alpha=0.25)
        for ax in flat_axes[len(plot_specs) :]:
            ax.set_visible(False)
        for ax in flat_axes[-2:]:
            if ax.get_visible():
                ax.set_xlabel(r"$t-t_0$ [$M$]")
        flat_axes[0].legend(loc="best", frameon=False)
        fig.suptitle(f"SEOBNRv5EHM memory-mode waveform check, q={args.q:g}")
        png_path = output_dir / f"eob_circular_memory_q{args.q:g}_omega{args.omega_start:g}.png"
        fig.savefig(png_path, dpi=180)
        plt.close(fig)
        print(f"Saved comparison plot: {png_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
