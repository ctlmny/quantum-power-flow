# Answers to the reflection questions

These are answers to the seven questions at the end of [`reproduce_paper.ipynb`](reproduce_paper.ipynb) (Section 10). Unless marked otherwise, every number comes from the notebook's saved outputs, and the cell numbers refer to cells in that notebook. I reran the same calculations with the current code (Qiskit 2.5.0, Qiskit Aer 0.17.2, Python 3.13.15). The voltages, errors, round counts and mismatches were unchanged. The compiled circuit depths changed by a few layers; see Question 4.

---

## 1. What happened to the bus voltages when you increased demand?

**Answer:** Voltages fell at every bus except the reference bus. The largest drops were at the buses whose demand was increased (buses 1 and 2).

**Evidence:** In Section 7, cell 22 was run with the demand slider at **"Higher demand" (×1.3)**:

| Bus | Normal demand (p.u.) | ×1.3 demand, classical (p.u.) | ×1.3 demand, HHL (p.u.) | Change |
|---|---:|---:|---:|---:|
| 1 | 0.94302 | 0.85474 | 0.85474 | −0.088 |
| 2 | 0.94446 | 0.85800 | 0.85800 | −0.086 |
| 3 | 0.99486 | 0.95345 | 0.95345 | −0.041 |
| 4 | 1.01825 | 0.99190 | 0.99190 | −0.026 |
| 5 (reference) | 1.00200 | 1.00200 | 1.00200 | 0 |

**Explanation:**
- **Buses 1 and 2** drop the most. They are the heaviest consumers relative to their connections, and they are fed only through lines from bus 5 and from each other.
- **Buses 3 and 4** drop less. Bus 4 is a small generator (injection +0.20 p.u.), which helps hold up voltage on that side of the network.
- **Bus 5** is the reference ("slack") bus, so its voltage is fixed at 1.002 p.u. by definition and cannot change.
- **The HHL simulation matched the classical solver** at the higher demand to all five displayed decimals.

At ×1.3, buses 1 and 2 fall to about 0.85–0.86 p.u. That is well below the roughly 0.95 p.u. lower limit commonly used in grid operation. So a 30% increase in demand would already be an operating problem in this network.

---

## 2. How did the mismatch change across adjustment rounds?

**Answer:** The mismatch fell steadily, by roughly a factor of 5–20 each round. It dropped below the 1e-6 stopping tolerance after **6 rounds**.

**Evidence:** Section 5, cell 15 (convergence log for the 4-phase-qubit HHL simulation, normal demand):

| Round | Maximum mismatch | Bus 3 voltage (p.u.) | Bus 4 voltage (p.u.) |
|---:|---:|---:|---:|
| 1 | 0.14348 | 1.01398 | 1.02798 |
| 2 | 0.00806 | 0.99504 | 1.01841 |
| 3 | 0.00131 | 0.99501 | 1.01826 |
| 4 | 0.00006 | 0.99486 | 1.01825 |
| 5 | 0.00001 | 0.99486 | 1.01825 |
| 6 | 0.00000 | 0.99486 | 1.01825 |

The table rounds to five decimals. At full precision the last three values are 6.5e-5, 6.7e-6 and 6.0e-7, so round 6 is the first below the 1e-6 tolerance.

**Explanation:**
- The left panel of cell 15's figure plots the mismatch on a log scale. It falls almost as a straight line, which means each round removes a roughly constant *fraction* of the remaining error. That pattern is called linear convergence, and it is what the fast-decoupled method is expected to do.
- The right panel shows the voltages settling. After round 1 overshoots (bus 3 at 1.014 p.u.), they are stable to five decimals by round 4.
- The classical solver also needed 6 rounds (cell 6).

---

## 3. Which precision level produced the smallest adjustment error?

**Answer:** **High precision (4 phase qubits)**, with an adjustment error of **0.00%**. This was also the only level that needed no more rounds than the classical solver.

**Evidence:** Section 8, cell 27:

| Precision | Phase qubits | Adjustment error (%) | Adjustment rounds |
|---|---:|---:|---:|
| Low | 2 | 25.00 | 14 |
| Medium | 3 | 24.11 | 13 |
| **High** | **4** | **0.00** | **6** |

Cell 25 (precision toggle set to 4) also prints "Representative adjustment error: 0.00%" and "Power-flow adjustment rounds: 6".

**Explanation:**
- The matrix used in this experiment (cell 11) has eigenvalues of exactly 5 and 15. With a scale of 16, four phase qubits divide the range into 16 steps of exactly 1.0. Both eigenvalues land exactly on a step (5/16 and 15/16), so nothing is rounded and the inverse is exact.
- With 2 or 3 phase qubits the steps are 4.0 or 2.0 wide. The eigenvalue 15 cannot be represented exactly, so it gets rounded and the error is large (about 25%).
- **Effect on rounds:** a less accurate correction did not give a wrong final answer; the loop just needed more rounds (14 and 13 instead of 6). This is because the outer loop recomputes the true power imbalance with exact classical arithmetic every round, so errors from one correction get fixed in the next.

**Additional note (not from the notebook):** The notebook's "adjustment error" comes from an idealized model of the HHL circuit (see the note at the end). When I simulated the circuit itself, the ranking stayed the same, but the medium-precision error is much smaller than the model says:

| Phase qubits | Notebook's model | Simulated circuit |
|---:|---:|---:|
| 2 | 25.0% | 20.4% |
| 3 | 24.1% | **4.1%** |
| 4 | 0.00% | 0.00% |

The exact circuit spreads a value that falls between two steps across neighbouring steps, rather than rounding it to one, and the idealized model ignores that.

---

## 4. What resource cost increased with precision?

**Answer:** Both the **number of qubits** and the **circuit depth** (sequential layers of operations) increased. Depth increased much faster than the qubit count.

**Evidence:** Section 8, cell 27:

| Precision | Phase qubits | Logical qubits | Compiled circuit depth |
|---|---:|---:|---:|
| Low | 2 | 5 | 102 |
| Medium | 3 | 6 | 220 |
| High | 4 | 7 | 817 |

Cell 18 shows the 2-phase-qubit circuit uses 5 qubits with a raw depth of 13 before compilation.

**Explanation:**
- **Qubits:** each extra phase qubit adds exactly one qubit. The other 3 qubits are fixed: 2 system qubits hold the 4-number imbalance, plus 1 helper (ancilla) qubit.
- **Depth:** going from 2 to 4 phase qubits multiplies depth by about **8**. The controlled-rotation step (Section 6, item 3) needs one rotation for every possible value the phase register can record, which is 2ⁿ − 1 rotations: 3, 7 and 15 for n = 2, 3, 4. Each rotation also gains a control qubit as n grows, which makes it more expensive to compile. Phase estimation itself also grows.
- **Trade-off:** higher precision makes each circuit more expensive but reduces the number of rounds (14 → 6), so it also reduces the number of times the circuit must be run. The notebook's figure title makes the same point: "Higher precision improves the adjustment but requires a deeper circuit."

**Caveats:**
- With the current code the depths are **101 / 217 / 809**, because I fixed an ordering bug in the controlled rotations (see the note at the end). The saved outputs predate that fix; the trend is unchanged.
- These depths come from compiling for the Aer simulator, which keeps some multi-qubit unitary blocks intact. Compiled to basic hardware gates, the counts are larger.

---

## 5. What do the stochastic scatter plots say about relationships between physical quantities?

**Answer:** The bus 3 and bus 4 voltages move together almost perfectly: when one rises, the other rises. The input powers are also correlated, but much less strongly. The network turns moderately correlated inputs into nearly perfectly correlated voltages.

**Evidence:**
- **Cell 29** prints "Converged samples: 5000 / 5000" and "Observed voltage correlation: **0.997**". The input powers were generated with a target correlation of **0.75** (`correlation=0.75`).
- **Cell 30, panel (b)**, "Voltage correlation": the V3-vs-V4 scatter is a very narrow diagonal band.
- **Cell 30, panel (c)**, "Input-power correlation": the bus 4 injection vs bus 3 load scatter is a wider diagonal cloud, matching the weaker 0.75 correlation.
- **Panels (a) and (d)**: the voltage histograms are narrow. Cell 31 explains this as small voltage changes for the 1% input variations.

**Explanation:**
- **Why the voltages are so tightly linked:** buses 3 and 4 are both connected to the reference bus and to each other, forming a triangle. A change at either bus shifts the power flowing through that shared triangle, so both voltages respond together.
- **Why panel (c) is diagonal:** the input-power plot shows the correlation that was *put into* the experiment. Panel (b) shows the correlation that *comes out* of the physics.
- **Reading the plots correctly:** as cell 31 notes, each dot is one complete power-flow solution, not one quantum measurement. The scatter shows variation in the physical system, not quantum randomness.
- **Strength vs size:** a strong correlation does not mean large changes. The histograms show the voltages vary only slightly for 1% load variations.

---

## 6. Why did we compare the HHL result with a classical reference?

**Answer:** The classical fast-decoupled solver is a trusted, well-understood method. Comparing against it is the only way to tell whether the HHL simulation gives the right answer, and to measure *how* wrong it is when it doesn't. HHL output cannot be trusted without that check.

**Evidence:**
- **When it works:** cell 13 shows the 4-phase-qubit HHL simulation agreeing with the classical solver at every bus, with a largest voltage difference of **7.38e-08 p.u.** Cell 16's bar chart shows the two methods give the same voltages.
- **When it doesn't:** cell 27 shows that at lower precision each HHL correction is off by **25%** and the loop needs **14 rounds instead of 6**. Without a reference, nothing would show that these corrections are inaccurate.
- **Same answer at higher demand:** cell 22 shows the two methods also agree at ×1.3 demand.
- **Notebook's own warning:** cell 32 states that "seeing a correct simulation result does not prove that a noisy quantum computer would produce the same result." The comparison is what establishes correctness.

**Explanation:**
- **Errors are hidden by the loop.** As Question 3 showed, the outer loop corrects inaccurate HHL steps. So a broken HHL solver can still produce the correct *final* voltages, just more slowly. Only comparing each step and the round count against a classical reference exposes the problem.
- **This happened in practice.** While reproducing this project I found that the HHL circuit, when actually simulated, had an ordering bug that made each correction about 48% wrong. The final voltages still came out correct, so the bug had gone unnoticed. Comparing against the classical solver step by step exposed it. (Details: [`analysis/REPORT.md`](../analysis/REPORT.md).)
- **What the reference provides.** The comparison defines "correct", gives a yardstick for choosing settings such as precision, and catches errors in the simulation code itself.

---

## 7. What additional evidence would be needed before claiming an advantage over classical computing?

**Answer:** Much more than this notebook shows. A credible claim would need a full, honest cost comparison against the best classical method on problems large enough to matter, run on real or realistically noisy hardware.

**Evidence that the notebook does *not* yet show an advantage:**
- Cell 0 and cell 32 both say the results are from a simulation, not a quantum computer, and that a correct simulation does not prove a noisy quantum computer would give the same result.
- The network has only **4 unknowns** (cell 11), which a classical computer solves essentially instantly. The quantum circuit for those 4 numbers already needs **7 qubits and depth 817** (cell 27).
- Even at the best precision, HHL needed **the same 6 rounds** as the classical solver (cells 6 and 13), so it did not reduce the work of the outer loop.
- The HHL solve is only one step inside a classical loop (cells 8 and 9). The imbalance calculation, the voltage updates and the stopping check all remain classical.

**Evidence that would be needed:**
1. **Runs on real hardware or a realistic noise model**, showing the corrections stay accurate enough for the loop to converge.
2. **The full cost of reading out the answer.** HHL leaves the solution encoded in quantum amplitudes, and getting usable numbers requires running the circuit many times. Accuracy to 4 decimals would take on the order of 10⁸ runs per correction. A real comparison must include this, plus the cost of loading the imbalance into the qubits, and of the many runs discarded because the helper qubit did not report success.
3. **Fair accounting of the circuit construction.** Here the matrix exponential used inside the circuit is computed classically and turned into gates, which does not scale to large networks. An honest comparison needs an efficient quantum method for that step.
4. **Scaling evidence on large networks.** HHL's possible advantage only appears for very large, well-conditioned systems. Results would be needed for networks with thousands of buses, compared against modern sparse classical solvers, including how performance depends on the matrix condition number.
5. **End-to-end timing.** Total wall-clock time and cost for the whole power-flow calculation, including the classical loop and the quantum–classical communication, compared with the best classical alternative.

The repository's `CONTRIBUTING.md` sets the same bar: simulator results should not be presented as speedup evidence "without a complete accounting of state preparation, Hamiltonian simulation, measurement, and classical post-processing."

---

## Note on what "HHL simulation" means in this notebook

The notebook's HHL results come from the repository's default setting (`execution="spectral_model"`). That setting computes the answer an *idealized* HHL circuit would give, using ordinary linear algebra. It builds the circuit in Cell 18 for display and gate counting, but does not simulate it. When the circuit itself is simulated (`execution="statevector"`):

- with the bit-order fix, it reproduces the notebook's high-precision result to about 1e-12;
- its medium-precision error is smaller than the model predicts (Question 3).

The investigation, the bug fix and the precision/loading study are in [`analysis/REPORT.md`](../analysis/REPORT.md).
