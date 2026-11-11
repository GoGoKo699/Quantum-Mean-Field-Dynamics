"""Bounded illustrative IID/RQMC comparison for a nonlinear full-chain MV model.

Not a baseline-optimality or quantum-crossover experiment. It compares particle
estimators at the same finite Euler mesh. Sobol gets Brownian bridge and Rouse
normal-mode ordering, both also available to the IID implementation. Both use
the identical sampled-input distribution, force code, payoff, and particle cost.
"""
from pathlib import Path
import json
import time
import numpy as np
from scipy.fft import idct
from scipy.special import ndtri
from scipy.stats import qmc

OUT = Path(__file__).resolve().parent
L, NSTEPS, T, KAPPA, DIFF, CENTER_STD = 32, 32, 1.0, 4.0, 0.2, 0.25
QVEC = np.array([[1.0, 0.7, 0.3], [0.2, 1.1, 0.6], [0.5, -0.4, 0.9]])
VCOEF = np.array([1.3, -0.7, 0.9])
COUNTS = [64, 128, 256, 512, 1024]
REPLICATES = 20


def gaussian_drivers(n, kind, seed):
    dim = (NSTEPS + 1) * L * 3
    if kind == "rqmc":
        u = qmc.Sobol(d=dim, scramble=True, seed=seed).random_base2(int(np.log2(n)))
        z = ndtri(np.clip(u, np.finfo(float).eps, 1.0 - np.finfo(float).eps))
    else:
        z = np.random.default_rng(seed).standard_normal((n, dim))
    z = z.reshape(n, NSTEPS + 1, L, 3)

    # Exact Gaussian equilibrium of the isolated harmonic chain, with an
    # explicitly localized Gaussian center, hence nonzero Fourier mean modes.
    lam = 2.0 * (1.0 - np.cos(np.pi * np.arange(L) / L))
    std = np.empty(L)
    std[0] = np.sqrt(L) * CENTER_STD
    std[1:] = np.sqrt(DIFF / (KAPPA * lam[1:]))
    x = idct(z[:, 0] * std[None, :, None], type=2, norm="ortho", axis=1)

    # Brownian bridge in independent Rouse normal modes. Lowest temporal and
    # spatial frequencies are given early Sobol coordinates.
    b = np.zeros((n, NSTEPS + 1, L, 3))
    b[:, -1] = np.sqrt(T) * z[:, 1]
    intervals, next_z = [(0, NSTEPS)], 2
    while intervals:
        new_intervals = []
        for left, right in intervals:
            if right - left <= 1:
                continue
            mid = (left + right) // 2
            scale = np.sqrt(T * (right - left) / (4.0 * NSTEPS))
            b[:, mid] = 0.5 * (b[:, left] + b[:, right]) + scale * z[:, next_z]
            next_z += 1
            new_intervals.extend([(left, mid), (mid, right)])
        intervals = new_intervals
    increments = idct(np.diff(b, axis=1), type=2, norm="ortho", axis=2)
    return x, increments


def run_one(n, kind, seed, field_scale=1.0):
    x, increments = gaussian_drivers(n, kind, seed)
    dt = T / NSTEPS
    initial_mc = np.cos(x @ QVEC.T).mean(axis=(0, 1))
    for step in range(NSTEPS):
        phase = x @ QVEC.T
        s, c = np.sin(phase), np.cos(phase)
        mc = c.mean(axis=(0, 1))
        # The supplied initial Gaussian is invariant under R -> -R and the
        # exact drift is odd when sine means vanish. Uniqueness preserves that
        # symmetry, so both classical methods get these known zero moments
        # for free rather than estimating avoidable finite-ensemble noise.
        ms = np.zeros(len(QVEC))
        lap = np.empty_like(x)
        lap[:, 0] = x[:, 0] - x[:, 1]
        lap[:, -1] = x[:, -1] - x[:, -2]
        lap[:, 1:-1] = 2.0 * x[:, 1:-1] - x[:, :-2] - x[:, 2:]
        nonlinear = ((s * mc - c * ms) * (field_scale * VCOEF)) @ QVEC
        x += dt * (-KAPPA * lap + nonlinear) + np.sqrt(2.0 * DIFF) * increments[:, step]
    phase = x @ QVEC.T
    c, s = np.cos(phase), np.sin(phase)
    per_chain_c, per_chain_s = c.mean(axis=1), s.mean(axis=1)
    return {
        "density_cos_q1": float(per_chain_c[:, 0].mean()),
        "intrachain_scattering_q1": float(np.mean(per_chain_c[:, 0]**2 + per_chain_s[:, 0]**2)),
        "initial_density_cos_q1": float(initial_mc[0]),
        "terminal_mean_modes": c.mean(axis=(0, 1)).tolist(),
    }


def main():
    started = time.perf_counter()
    rows = []
    for n in COUNTS:
        for kind in ["iid", "rqmc"]:
            t0 = time.perf_counter()
            values = [run_one(n, kind, 400000 + 1000*n + r) for r in range(REPLICATES)]
            row = {"n": n, "method": kind, "replicates": REPLICATES,
                   "elapsed_seconds": time.perf_counter() - t0}
            for name in ["density_cos_q1", "intrachain_scattering_q1"]:
                a = np.array([v[name] for v in values])
                row[name] = {"mean": float(a.mean()), "variance": float(a.var(ddof=1)),
                             "standard_error_of_mean": float(a.std(ddof=1)/np.sqrt(REPLICATES)),
                             "replicate_values": a.tolist()}
            rows.append(row)
            print(json.dumps({k: v for k, v in row.items() if k not in ["density_cos_q1", "intrachain_scattering_q1"]}
                  | {"density_mean": row["density_cos_q1"]["mean"],
                     "density_variance": row["density_cos_q1"]["variance"]}), flush=True)

    # Matched seeded reference verifies that the self-consistent field remains
    # nonzero and materially changes the output. It is not a certified exact mean.
    reference = {"field_on": run_one(2048, "rqmc", 700001, 1.0),
                 "field_off": run_one(2048, "rqmc", 700001, 0.0)}
    fits = {}
    for kind in ["iid", "rqmc"]:
        by_kind = [r for r in rows if r["method"] == kind]
        fits[kind] = {}
        for name in ["density_cos_q1", "intrachain_scattering_q1"]:
            fit = np.polyfit(np.log([r["n"] for r in by_kind]),
                             np.log([r[name]["variance"] for r in by_kind]), 1)
            fits[kind][name] = {"log_variance_log_n_slope": float(fit[0])}
    result = {"config": {"L": L, "d": 3*L, "K": 2*len(QVEC), "nsteps": NSTEPS,
                          "T": T, "kappa": KAPPA, "diffusion": DIFF,
                          "center_std": CENTER_STD, "q_vectors": QVEC.tolist(),
                          "exact_reflection_symmetry_used": True,
                          "nonzero_mean_coordinates": len(QVEC),
                          "v_coefficients": VCOEF.tolist(), "replicates": REPLICATES},
              "rows": rows, "reference": reference, "variance_slope_fits": fits,
              "total_elapsed_seconds": time.perf_counter() - started,
              "limitations": ["Finite time-mesh bias is not controlled by this experiment.",
                              "Finite particle bias is not certified by this experiment.",
                              "Twenty randomizations give noisy empirical variance slopes.",
                              "No quantum resource or crossover claim is tested.",
                              "No claim that these classical methods or parameter choices are optimal."]}
    (OUT / "classical_polymer_stress_results.json").write_text(json.dumps(result, indent=2))
    print(json.dumps({"variance_slope_fits": fits, "reference": reference,
                      "total_elapsed_seconds": result["total_elapsed_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
