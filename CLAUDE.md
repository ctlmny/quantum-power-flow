# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

A reproduction of the five-bus "Quantum Power Flow" paper (Feng, Zhou, Zhang, arXiv:2104.04888). It solves AC power flow with the fast-decoupled method and replaces each linear solve with a custom Harrow–Hassidim–Lloyd (HHL) circuit simulated in Qiskit 2.5 / Qiskit Aer. The goal is scientific transparency: every reproduced number, every inferred parameter, and every disagreement with the paper is documented. Treat it as a research-reproduction codebase, not a general power-flow library.

## Commands

```bash
python -m pip install -e '.[dev]'          # dev install (add ,notebook for Jupyter extras)
ruff check src tests scripts               # lint (CI runs this)
pytest -q                                  # full test suite (CI runs this on Python 3.10–3.13)
pytest tests/test_hhl.py::test_prepared_hhl_matrix_is_reused   # single test
qpf-reproduce all --output results         # regenerate every figure/table/diagnostic under results/
python scripts/infer_missing_line_data.py  # rerun the line-impedance inference alone
```

On Windows the venv lives at `.venv\Scripts\`. `mypy` is configured in `pyproject.toml` but is not run by CI.

If a change affects numerical outputs, rerun `qpf-reproduce all`. The files under `results/` are committed, and the README's results table and `results/diagnostics/reproduction_metrics.json` must agree with them.

## Architecture

The core design is a **pluggable linear-solver backend** (`src/qpf_repro/linear.py`):
- `LinearSystemFactory.prepare(matrix)` returns a `PreparedLinearSystem` with `solve(rhs)` / `solve_many(rhs_batch)`.
- `NumpyLinearSystemFactory` is the classical reference. `quantum/hhl.py::QiskitHHLFactory` is the HHL backend.
- `powerflow.solve_fast_decoupled(case, linear_system_factory=...)` calls `prepare` **once** for each of `B′` and `B″`, then reuses them every iteration. This implements the paper's idea that the expensive circuit parts (QPE and the reciprocal rotation) are built once and cached. Any change must keep these prepared objects reused across iterations; `tests/test_hhl.py` checks this.

Data flow: `cases.py` builds a `PowerFlowCase` (`models.py`). `network.py` builds the Y-bus, mismatches, the reduced `B′`/`B″` (by default `-imag(Ybus[pq,pq])`), and branch flows. `powerflow.py` iterates (it also contains a Newton–Raphson solver for comparison). `experiments/reproduce.py::reproduce_all` runs every experiment and writes outputs through `reporting.py` and `plotting.py`. `reference_data.py` holds the paper's Table I values.

HHL details (`quantum/hhl.py`, explained in `docs/algorithm_mapping.md`):
- Seven logical qubits: 4 phase, 2 system, 1 ancilla.
- `phase_scale` comes from a spectral fit so that the eigenvalues fall exactly on phase bins.
- The global phase is fixed by rotating the postselected state so that its overlap with `b` is positive real. This gives signed real updates without consulting the direct solution.

There are two execution modes, set with `HHLConfiguration.execution`:
- `"spectral_model"` (the default, used for all committed `results/`): `solve()` applies the idealized map `V diag(1/λ_bin) Vᵀ` and does **not** run the circuit.
- `"statevector"`: each solve runs the right-hand side through the circuit using a cached unitary, then postselects.

Only the second mode tests the circuit itself, and it shows effects of inexact bins (QPE leakage) that the spectral model misses. Qiskit's `phase_estimation` returns the phase register bit-reversed, so the reciprocal rotation controls on `_reverse_bits(bin)`; `tests/test_hhl.py` guards this. The precision and loading study, with its write-up, is in `analysis/` (`hhl_resource_study.py`, `REPORT.md`).

`stochastic.py` (5,000 Monte Carlo samples) does **not** simulate circuits per sample. It applies the circuit-equivalent linear map `V_eig @ diag(1/λ_phase_bin) @ V_eig.T` in a vectorized batch and verifies a few samples against full Aer circuits. Keep the distinction between this classical acceleration and an actual quantum execution explicit.

## Gotchas

- **The `configs/*.yaml` files are documentation and are not loaded by the code.** The cases are hard-coded in `src/qpf_repro/cases.py`, and solver/HHL/stochastic parameters are passed in `experiments/reproduce.py`. If you change a parameter, update both the code and the YAML.
- The YAML and the paper use 1-based bus labels (slack = bus 5). The code uses 0-based indices (`slack_bus=4`).
- The six line impedances (`0.02 + j0.20` p.u.) and the stochastic variance are **reconstructions**; the paper does not print them. See `docs/reconstruction_notes.md`. Label any new assumed value as a reconstruction, not as a paper parameter.
- Known discrepancies are deliberate and tested: the Table I iteration-2 `θ4 = -0.0340` is treated as a typo, and the stressed case converges in 32 updates where the paper reports 34. Do not "fix" code to match them.
- Do not present simulator results as evidence of a quantum speedup (`CONTRIBUTING.md`, `docs/limitations_and_interpretation.md`).
- `CONTEXT.md` defines the student-facing vocabulary used in `notebooks/reproduce_paper.ipynb`. For example, say "HHL simulation", not "quantum-computer result", and "adjustment round", not "training step". Use those terms when editing the notebook.
