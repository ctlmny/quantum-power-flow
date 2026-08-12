# Changelog

All notable changes to this reproduction are documented here.

## Unreleased

## 0.2.0 — 2026-08-13

### Added

- README results section with committed reproduction figures, a Mermaid algorithm
  diagram, and reader-oriented figure labels; `results/figures/` is now versioned and
  regenerated deterministically by the reproduction command.
- Guided student notebook with an iterative algorithm diagram, linear equations, beginner pseudocode, convergence logs, demand and quantum-precision controls, stochastic voltage distributions, result comparisons, and an explained Qiskit circuit.
- Optional notebook environment with JupyterLab and circuit-drawing support.
- Continuous integration now tests Python 3.10 and 3.11 in addition to 3.12 and
  3.13, so the declared `requires-python` range is verified end to end.

## 0.1.0 — 2026-07-19

### Added

- Modern-Qiskit implementation of HHL for the fixed fast-decoupled power-flow matrices.
- Five-bus AC network reconstruction and parameter-inference script.
- Classical fast-decoupled and Newton–Raphson baselines.
- Normal, stressed, phase-register, and 5,000-sample stochastic experiments.
- Reproduced tables, figures, circuit diagnostics, and machine-readable metrics.
- Automated regression tests and GitHub Actions CI.
- Documentation of reconstruction assumptions, discrepancies, and quantum-complexity limitations.
