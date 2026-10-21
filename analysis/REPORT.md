# Milestone 1: Reproducing the five-bus HHL power flow, then varying precision and load

All numbers below come from `analysis/hhl_resource_study.py`, written to `analysis/output/study/`. Everything ran as noiseless statevector simulation; no hardware was used.

## Environment

| Package | Version |
|---|---|
| Python | 3.13.15 (Windows 11, x64) |
| Qiskit | 2.5.0 |
| Qiskit Aer | 0.17.2 |
| NumPy | 2.5.3 |
| SciPy | 1.18.1 |
| pandas | 3.0.6 |
| matplotlib | 3.11.2 |
| PyYAML | 6.0.3 |

These are the versions the repository pins (`qiskit==2.5.0`, `qiskit-aer==0.17.2`). The recorded copy is in `analysis/output/study/environment.json`.

## 1. Reproducing the normal operating case (paper Test I)

The setup follows the repository's defaults: 4 phase qubits, exact-bin ("spectral fit") phase scale, and 6 adjustment rounds as in the paper's Table I. "Reference" means a Newton–Raphson solve converged to 1e-12.

| Method | V3 (p.u.) | V4 (p.u.) | θ3 (rad) | θ4 (rad) | Max diff vs paper after 6 rounds | Rounds to 1e-6 mismatch |
|---|---:|---:|---:|---:|---:|---:|
| Paper, Table I | 0.9948 | 1.0182 | −0.1144 | −0.0393 | — | — |
| Classical FDLF (NumPy) | 0.994860 | 1.018247 | −0.114375 | −0.039248 | 6.0e-5 | 6 |
| HHL, repo default (spectral model) | 0.994860 | 1.018247 | −0.114375 | −0.039248 | 6.0e-5 | 6 |
| **HHL, executed circuit (fixed)** | 0.994860 | 1.018247 | −0.114375 | −0.039248 | 6.0e-5 | 6 |
| HHL, executed circuit (original) | 0.996081 | 1.019350 | −0.112987 | −0.038067 | 1.4e-3 | 23 |

- Plot of the iteration traces: `figures/test_i_trace_comparison.png`
- Convergence plot: `figures/test_i_convergence.png`
- Per-round table: `tables/test_i_iteration_comparison.csv`

**What matched.** With the fix below, the executed HHL circuit agrees with classical FDLF to within about 2.5e-14 on every linear solve, and the two trajectories coincide in every round. Final values agree with the paper to its four printed decimals; the 6e-5 difference is rounding. The paper's iteration-2 value θ4 = −0.0340 is the one exception, and the repository already documents it as a suspected typo (`docs/reconstruction_notes.md`). Running `pytest` (10 tests) and `qpf-reproduce all` reproduces the committed results to within 2.4e-15.

**What differed, and why it matters.** In the repository as published, the "quantum" result is not obtained by running a circuit. `QiskitHHLPreparedSystem.solve()` multiplies by a precomputed matrix, `V · diag(1/λ_bin) · Vᵀ`. The Qiskit circuit is built only to report gate counts. So the README's reported quantum–classical difference (exactly 0.0 when rerun) compares two classical calculations. When I ran the circuit itself, the 4-phase-qubit case had **48% error on each linear solve**, and the power flow needed 23 rounds instead of 6.

## 2. Focused change: fixing the phase-register bit order and running the circuit

**Cause.** Qiskit 2.5's `phase_estimation` leaves the phase register in bit-reversed order. An eigenvalue in bin 5 (`0101`) comes out as `1010` = 10. The reciprocal rotation was controlled on `ctrl_state=logical_bin`, so the small eigenvalue (4.95) received the rotation meant for 9.9. That halves its contribution to the solution. Bin 15 (`1111`) reads the same in both orders, which is why part of the answer still looked correct.

**Change** (in `src/qpf_repro/quantum/hhl.py`):
1. The reciprocal rotation now controls on the bit-reversed bin.
2. New opt-in option `HHLConfiguration(execution="statevector")`. Each solve runs the right-hand side through the circuit, postselects phase register `|0…0⟩` with ancilla `|1⟩`, fixes the global phase, and rescales the result. The circuit's unitary is computed once and reused for every solve.
3. Two new tests. One checks that the executed circuit matches a direct solve for 2 and 4 phase qubits. The other checks that the solve agrees with a full Aer statevector run to 1e-12. The first test **fails on the original code and passes after the fix**.

The default (`execution="spectral_model"`) is unchanged, so none of the committed `results/` change. All 12 tests pass and `ruff` is clean.

**Fair comparison.** Both versions use the same case, start point, tolerance (1e-6), round limit (300), calibration and number of phase qubits. The only difference is the control state. Errors are for one representative solve; rounds are counts to converge (`tables/original_vs_fixed_circuit.csv`).

| Calibration | Phase qubits | Original circuit: solve error / rounds | Fixed circuit: solve error / rounds |
|---|---:|---|---|
| Exact bins | 2 | 48% / 23 | 1.4e-15 / **6** |
| Exact bins | 3 | 28% / **diverged** | 1.2e-15 / **6** |
| Exact bins | 4 (paper) | 48% / 23 | 2.5e-14 / **6** |
| Exact bins | 5 | 28% / **diverged** | 4.2e-14 / **6** |
| Power-of-two | 4 | 48% / 23 | 1.2% / **6** |
| Power-of-two | 6 | 83% / **diverged** | 2.3% / **7** |

At the paper's setting, the fix also raises the postselection success probability from 0.062 to 0.164. That cuts the expected number of circuit runs per power flow from 749 to 114.

**What I learned.**
1. The outer power-flow loop hides solver errors. The mismatch ΔP, ΔQ is always computed classically and exactly, so an inaccurate linear solve only slows convergence; the loop still lands on the correct operating point as long as it converges. That is why the bug went unnoticed: the final answer was correct, just reached 4× more slowly.
2. The spectral model does not predict how the real circuit behaves when eigenvalues fall between phase bins. With a scale of 16, the model predicts 0.99% error for 4, 5 and 6 phase qubits. The executed circuit gives 1.2%, 2.5% and 2.3%. More phase qubits gave *worse* results, because phase-estimation probability leaks into neighbouring bins, which the model ignores.

## 3. Varying precision and loading

**Settings swept.**
- Phase qubits: 1–6.
- Two phase-scale calibrations:
  - **Exact bins.** The scale is chosen so every eigenvalue lands exactly on a bin. This requires knowing the eigenvalues of B′, i.e. a classical eigendecomposition, so it acts as an oracle.
  - **Power-of-two.** Scale 16 here, which needs only an upper bound on the largest eigenvalue. It is identical to the notebook's `fixed_scale=16`.
- Demand multipliers 0.5–3.0 on the consuming buses, defined as in the notebook, plus the paper's stressed Test II. Both classical FDLF and Newton fail beyond about 3.03×, the voltage-collapse point.

**Metrics.**
- Accuracy against a Newton–Raphson reference.
- Rounds to reach a 1e-6 mismatch.
- Gate counts after transpiling to the basis {cx, rz, sx, x} at optimization level 1 with a fixed seed.
- Postselection success probability.
- Expected CX executions for a whole power flow: Σ over solves of (1 / success probability) × CX per circuit.

Figures: `precision_tradeoff.png`, `loading_convergence.png`, `total_quantum_cost.png`. Data: `tables/precision_loading_sweep.csv`, `tables/circuit_resources.csv`.

**Circuit size and normal-demand behaviour (power-of-two scale):**

| Phase qubits | Logical qubits | CX | Depth | Error per solve | Rounds at 1× demand | Rounds at 3× demand |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 4 | 23 | 111 | 61% | 141 | 191 |
| 2 | 5 | 81 | 311 | 18% | 24 | diverged |
| 3 | 6 | 225 | 817 | 5.8% | 8 | 156 |
| **4** | 7 | 489 | 1095 | 1.2% | **6** | **70** |
| 5 | 8 | 1405 | 2475 | 2.5% | 7 | 70 |
| 6 | 9 | 3745 | 5443 | 2.3% | 7 | 69 |
| Classical FDLF | — | — | — | — | 6 | 68 |

**Findings.**
- **Accuracy.** Every run that converged was within 2.3e-6 of the Newton reference. That meets a 1e-4 target, which matches the paper's printed precision. Precision therefore affects how many rounds are needed and whether the loop converges, not the final accuracy.
- **Convergence vs load.** Low precision costs more rounds as demand rises. At 3 phase qubits the penalty is 8 vs 6 rounds at 1×, but 156 vs 68 at 3×, and 76 vs 39 in Test II. Runs that failed to converge:
  - 2 phase qubits, power-of-two scale, at 3×;
  - 1 phase qubit, exact bins, at 2.8×, 3× and Test II.
- **Circuit cost grows exponentially with phase qubits.** The reciprocal rotation uses 2ⁿ − 1 multi-controlled RY gates, so CX count rises from 23 to 3745 across 1–6 phase qubits.

**Cheapest settings that converge within 0.0001 of the reference** (`tables/cheapest_settings_meeting_target.csv`):

| Calibration | Cheapest setting | Expected CX per power flow, 1× demand | Same, Test II |
|---|---|---:|---:|
| Power-of-two (no eigenvalue knowledge) | **4 phase qubits** at every load | 5.1e4 | 2.8e5 |
| Exact bins (needs eigenvalues) | **2 phase qubits** at every load | 8.6e3 | 5.1e4 |

- **Without eigenvalue knowledge, 4 phase qubits is the cost-optimal choice.** Three qubits uses 54% fewer CX per circuit (225 vs 489), but it needs more rounds and succeeds at postselection less often. Its expected total is 9.2e4 CX, 1.8× more. Five and six qubits cost more and do no better.
- **With exact-bin calibration, 2 phase qubits gives the exact answer.** That uses 5 instead of 7 qubits, 75 instead of 485 CX, and about 6× less total quantum work, with the same number of rounds as classical FDLF at every load. **Caveat:** this depends on two things that limit it as a general result:
  - Choosing the scale requires B′'s eigenvalues, so the classical problem is effectively already solved.
  - The eigenvalues are exactly 4.95 and 14.85, a ratio of 3. That follows from the reconstructed network, in which all six branches have the same impedance. A real network would not have such a clean spectrum.
- **A 6-round budget is not achievable at high load, quantum or not.** At 2× demand and above, even classical FDLF needs more than 6 rounds. A fixed round budget only makes sense near nominal load.

## Limitations

- **The cost figures leave out readout.** They count one successful postselection per solve. Reading out the amplitudes to 1e-4 accuracy would take on the order of 10⁸ shots per solve, far more than the circuit cost. Readout cost is the most important open question for any claim about quantum resources.
- **Not a scalable circuit construction.** Each controlled-U gate is built by classically computing `expm(2πiA/scale)` and synthesizing the dense 2-qubit unitary. The norm ‖b‖/C is also taken as known rather than estimated.
- **Noiseless only.** No noise, device connectivity or shot noise is modelled. Nothing here is evidence of a quantum speedup.

## Other reproducibility notes

- **Upcoming Qiskit 3.0 break.** `Gate.control()` raises a deprecation warning about the `annotated` argument; the default changes in Qiskit 3.0. Pinning the current behaviour or adopting annotated controls should be done before upgrading.
- **Notebook label.** The notebook calls 4 phase qubits with `fixed_scale=16` "high precision", but on B′ it has 1.2% error per solve, because 4.95 does not land exactly on a bin.
- **Untracked result file.** `results/tables/stochastic_explicit_circuit_validation.csv` is committed but no code in the repository produces it. The claim in `docs/algorithm_mapping.md` that the circuit was rerun for five stochastic samples cannot be reproduced from the current code.
- **Config files not used.** `configs/*.yaml` are documentation only; the cases are hard-coded in `cases.py`.

## How to rerun

```bash
python -m pip install -e '.[dev]'
pytest -q
python analysis/hhl_resource_study.py              # ~10 min; 6-phase-qubit circuits take longest
python analysis/hhl_resource_study.py --plots-only # redraw figures from the saved tables
```
