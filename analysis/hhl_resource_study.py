#!/usr/bin/env python3
"""Reproduce the normal five-bus case on the executed HHL circuit and study precision/loading.

Run from the repository root:

    python analysis/hhl_resource_study.py

Outputs go to ``analysis/output/study/`` (tables, figures, environment.json).

Every quantum solve in this script evolves the right-hand side through the constructed Qiskit
circuit (``HHLConfiguration(execution="statevector")``), unless it is labelled "spectral model"
(the repository's default idealized inverse map) or "original circuit" (the circuit as built
before the phase-register bit-order fix, reconstructed below for a like-for-like comparison).
"""

from __future__ import annotations

import json
import platform
import sys
import warnings
from dataclasses import dataclass, replace
from importlib.metadata import version
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from qiskit import QuantumCircuit, QuantumRegister, transpile  # noqa: E402
from qiskit.circuit.library import RYGate, UnitaryGate, phase_estimation  # noqa: E402
from scipy.linalg import expm  # noqa: E402

from qpf_repro.cases import paper_case5, stressed_paper_case5  # noqa: E402
from qpf_repro.models import PowerFlowCase, PowerFlowResult  # noqa: E402
from qpf_repro.network import fast_decoupled_matrices  # noqa: E402
from qpf_repro.powerflow import solve_fast_decoupled, solve_newton_raphson  # noqa: E402
from qpf_repro.quantum import HHLConfiguration, QiskitHHLFactory  # noqa: E402
from qpf_repro.quantum.hhl import QiskitHHLPreparedSystem  # noqa: E402
from qpf_repro.reference_data import paper_table_i_dataframe  # noqa: E402
from qpf_repro.reporting import history_dataframe  # noqa: E402

warnings.filterwarnings("ignore", category=DeprecationWarning)

OUT = Path("analysis/output/study")
TABLES, FIGURES = OUT / "tables", OUT / "figures"
TOLERANCE = 1e-6
MAX_ITERATIONS = 300
ACCURACY_TARGET = 1e-4  # p.u. / rad: the paper reports four decimals
PAPER_ROUNDS = 6  # Table I budget
PHASE_QUBITS = (1, 2, 3, 4, 5, 6)
CALIBRATIONS = {
    # Exact bins: needs the eigenvalues of B' (a classical eigendecomposition) - an oracle.
    "spectral_fit": dict(scale_strategy="spectral_fit"),
    # Needs only an upper bound on the largest eigenvalue; here scale = 16, identical to the
    # notebook's "fixed_scale=16" setting.
    "power_of_two": dict(scale_strategy="power_of_two"),
}
DEMAND_MULTIPLIERS = (0.5, 1.0, 1.5, 2.0, 2.5, 2.8, 3.0)
BASIS = ["cx", "rz", "sx", "x"]
REPRESENTATIVE_RHS = np.asarray([-0.55, -0.55, -0.95, 0.20])

# Validated categorical slots (light mode), fixed order; text stays in ink colors.
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
N_RAMP = ["#c6dbf5", "#97bdee", "#6a9fe3", "#3f81d6", "#2a63b0", "#1d4580"]  # n = 1..6

plt.rcParams.update(
    {
        "figure.dpi": 150,
        "savefig.dpi": 200,
        "font.size": 9,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "legend.frameon": False,
        "lines.linewidth": 2,
    }
)


# --------------------------------------------------------------------------------------------
# Circuit variants
# --------------------------------------------------------------------------------------------


def original_core(prepared: QiskitHHLPreparedSystem) -> QuantumCircuit:
    """The HHL core exactly as built before the fix: ``ctrl_state=logical_bin``."""

    config = prepared.configuration
    n = config.phase_qubits
    phase = QuantumRegister(n, "phase")
    system = QuantumRegister(prepared.system_qubits, "system")
    ancilla = QuantumRegister(1, "ancilla")
    core = QuantumCircuit(phase, system, ancilla, name="HHL-core-original")
    unitary = UnitaryGate(expm(2j * np.pi * prepared.matrix / prepared.calibration.phase_scale))
    qpe = phase_estimation(n, unitary)
    core.compose(qpe, list(phase) + list(system), inplace=True)
    for logical_bin in range(1, 1 << n):
        lambda_est = logical_bin * prepared.calibration.phase_scale / (1 << n)
        angle = 2.0 * np.arcsin(min(1.0, prepared.reciprocal_constant / lambda_est))
        gate = RYGate(float(angle)).control(n, ctrl_state=logical_bin)
        core.append(gate, list(phase) + [ancilla[0]])
    core.compose(qpe.inverse(), list(phase) + list(system), inplace=True)
    return core


@dataclass
class OriginalCircuitFactory:
    """Statevector-executed HHL using the pre-fix circuit."""

    configuration: HHLConfiguration
    name: str = "qiskit_hhl_original_circuit"

    def __post_init__(self) -> None:
        self._inner = QiskitHHLFactory(self.configuration)

    def prepare(self, matrix: np.ndarray, *, label: str = "A") -> QiskitHHLPreparedSystem:
        prepared = self._inner.prepare(matrix, label=label)
        if prepared.core_circuit.name != "HHL-core-original":
            prepared.core_circuit = original_core(prepared)
            prepared._core_unitary = None
        return prepared


_FACTORIES: dict[tuple[str, int, str], object] = {}


def factory(calibration: str, n: int, variant: str = "circuit"):
    """Cache factories so each circuit unitary is built once per (calibration, n, variant)."""

    key = (calibration, n, variant)
    if key not in _FACTORIES:
        execution = "spectral_model" if variant == "spectral_model" else "statevector"
        config = HHLConfiguration(phase_qubits=n, execution=execution, **CALIBRATIONS[calibration])
        if variant == "original_circuit":
            _FACTORIES[key] = OriginalCircuitFactory(config)
        else:
            _FACTORIES[key] = QiskitHHLFactory(config)
    return _FACTORIES[key]


def circuit_resources(calibration: str, n: int) -> dict[str, object]:
    matrix, _ = fast_decoupled_matrices(paper_case5())
    prepared = factory(calibration, n).prepare(matrix, label="B_prime")
    circuit = prepared.circuit_for_rhs(REPRESENTATIVE_RHS)
    compiled = transpile(circuit, basis_gates=BASIS, optimization_level=1, seed_transpiler=7)
    ops = compiled.count_ops()
    return {
        "calibration": calibration,
        "phase_qubits": n,
        "logical_qubits": prepared.logical_qubits,
        "phase_scale": prepared.calibration.phase_scale,
        "eigenvalue_bins": "/".join(map(str, sorted(set(prepared.calibration.logical_phase_bins)))),
        "max_rel_eigenvalue_error": prepared.calibration.maximum_relative_eigenvalue_error,
        "controlled_rotations": (1 << n) - 1,
        "depth": compiled.depth(),
        "cx_count": int(ops.get("cx", 0)),
        "total_gates": int(sum(ops.values())),
    }


# --------------------------------------------------------------------------------------------
# Cases and metrics
# --------------------------------------------------------------------------------------------


def demand_case(multiplier: float) -> PowerFlowCase:
    """Scale every consuming bus, as in the notebook's demand control."""

    base = paper_case5()
    injections = base.specified_power_injection.copy()
    injections[injections.real < 0] *= multiplier
    return replace(base, name=f"demand_x{multiplier:g}", specified_power_injection=injections)


def state_error(result: PowerFlowResult, reference: PowerFlowResult) -> float:
    return float(
        max(
            np.max(np.abs(result.voltage_magnitude - reference.voltage_magnitude)),
            np.max(np.abs(result.voltage_angle_rad - reference.voltage_angle_rad)),
        )
    )


def run(case, linear_factory, **kwargs) -> PowerFlowResult | None:
    try:
        return solve_fast_decoupled(case, linear_system_factory=linear_factory, **kwargs)
    except RuntimeError:
        return None  # the iteration produced a non-physical / non-finite state


def solve_statistics(result: PowerFlowResult) -> dict[str, float]:
    errors, probabilities = [], []
    for record in result.history:
        for key in ("active_linear_solve", "reactive_linear_solve"):
            meta = record.metadata[key]
            errors.append(meta.get("relative_solution_error_vs_direct", 0.0))
            if "postselection_probability_phase_zero_ancilla_one" in meta:
                probabilities.append(meta["postselection_probability_phase_zero_ancilla_one"])
    probabilities_array = np.asarray(probabilities)
    return {
        "median_linear_solve_error": float(np.median(errors)) if errors else np.nan,
        "mean_postselection_probability": (
            float(np.mean(probabilities_array)) if probabilities else np.nan
        ),
        # Expected circuit executions to obtain one postselected success per linear solve.
        # A lower bound: amplitude readout would multiply this by the shots per solve.
        "expected_circuit_runs": (
            float(np.sum(1.0 / probabilities_array)) if probabilities else np.nan
        ),
    }


# --------------------------------------------------------------------------------------------
# Part 1 - normal case reproduction (paper Test I)
# --------------------------------------------------------------------------------------------


def part1_reproduction() -> dict[str, object]:
    case = paper_case5()
    reference = solve_newton_raphson(case, tolerance=1e-12, max_iterations=50)
    paper = paper_table_i_dataframe()
    methods = {
        "Classical FDLF (numpy)": None,
        "HHL spectral model (repo default)": factory("spectral_fit", 4, "spectral_model"),
        "HHL circuit, original": factory("spectral_fit", 4, "original_circuit"),
        "HHL circuit, bit-order fixed": factory("spectral_fit", 4, "circuit"),
    }
    rows, traces, converged_runs = [], {}, {}
    paper_final = paper.iloc[-1]
    for label, linear_factory in methods.items():
        fixed = run(case, linear_factory, fixed_iterations=PAPER_ROUNDS)
        to_tol = run(case, linear_factory, tolerance=TOLERANCE, max_iterations=MAX_ITERATIONS)
        traces[label] = history_dataframe(fixed)
        converged_runs[label] = to_tol
        final = traces[label].iloc[-1]
        rows.append(
            {
                "method": label,
                "V3_pu": final["V3"],
                "V4_pu": final["V4"],
                "theta3_rad": final["theta3_rad"],
                "theta4_rad": final["theta4_rad"],
                "max_abs_diff_vs_paper_after_6": float(
                    np.max(
                        np.abs(
                            final[["V3", "V4", "theta3_rad", "theta4_rad"]].to_numpy(float)
                            - paper_final[["V3", "V4", "theta3_rad", "theta4_rad"]].to_numpy(float)
                        )
                    )
                ),
                "max_state_error_vs_reference_after_6": state_error(fixed, reference),
                "rounds_to_tolerance_1e-6": to_tol.iterations if to_tol else None,
                "state_error_at_tolerance": state_error(to_tol, reference) if to_tol else None,
                **solve_statistics(to_tol),
            }
        )
    table = pd.DataFrame(rows)
    table.to_csv(TABLES / "test_i_method_comparison.csv", index=False)

    # Iteration-by-iteration comparison of HHL (fixed circuit) against classical.
    classical_trace = traces["Classical FDLF (numpy)"]
    circuit_trace = traces["HHL circuit, bit-order fixed"]
    side = paper.rename(columns=lambda c: c if c == "iteration" else f"paper_{c}")
    for column in ("V3", "V4", "theta3_rad", "theta4_rad"):
        side[f"classical_{column}"] = classical_trace[column]
        side[f"hhl_circuit_{column}"] = circuit_trace[column]
        side[f"abs_diff_{column}"] = (circuit_trace[column] - classical_trace[column]).abs()
    side.to_csv(TABLES / "test_i_iteration_comparison.csv", index=False)

    plot_test_i(paper, traces)
    plot_test_i_convergence(converged_runs)
    return {"table": table, "reference": reference}


def plot_test_i(paper: pd.DataFrame, traces: dict[str, pd.DataFrame]) -> None:
    styles = {
        "Classical FDLF (numpy)": dict(color=BLUE, linestyle="-"),
        "HHL circuit, bit-order fixed": dict(color=ORANGE, linestyle="--"),
        "HHL circuit, original": dict(color=AQUA, linestyle=":"),
    }
    panels = [
        ("V3", "Bus 3 voltage (p.u.)"),
        ("V4", "Bus 4 voltage (p.u.)"),
        ("theta3_rad", "Bus 3 angle (rad)"),
        ("theta4_rad", "Bus 4 angle (rad)"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(9, 6), sharex=True)
    for ax, (column, title) in zip(axes.flat, panels, strict=True):
        for label, style in styles.items():
            ax.plot(traces[label]["iteration"], traces[label][column], label=label, **style)
        ax.plot(
            paper["iteration"],
            paper[column],
            "o",
            color=INK,
            markersize=6,
            markerfacecolor="white",
            markeredgewidth=1.5,
            label="Paper Table I",
            zorder=5,
        )
        ax.set_title(title, loc="left")
    for ax in axes[1]:
        ax.set_xlabel("Adjustment round")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle(
        "Normal case (Test I): fixed HHL circuit overlays classical FDLF and the paper",
        x=0.01,
        ha="left",
        fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    fig.savefig(FIGURES / "test_i_trace_comparison.png", bbox_inches="tight")
    plt.close(fig)


def plot_test_i_convergence(runs: dict[str, PowerFlowResult | None]) -> None:
    styles = {
        "Classical FDLF (numpy)": dict(color=BLUE, linestyle="-", marker="o"),
        "HHL spectral model (repo default)": dict(color=VIOLET, linestyle="-.", marker="s"),
        "HHL circuit, bit-order fixed": dict(color=ORANGE, linestyle="--", marker="^"),
        "HHL circuit, original": dict(color=AQUA, linestyle=":", marker="D"),
    }
    fig, ax = plt.subplots(figsize=(7.5, 4))
    for label, style in styles.items():
        result = runs[label]
        if result is None:
            continue
        mismatch = [record.max_mismatch for record in result.history]
        ax.semilogy(range(1, len(mismatch) + 1), mismatch, label=label, markersize=4, **style)
    ax.axhline(TOLERANCE, color=MUTED, linewidth=1, linestyle=(0, (4, 3)))
    ax.text(1, TOLERANCE * 1.4, "tolerance 1e-6", color=MUTED, fontsize=8)
    ax.set_xlabel("Adjustment round")
    ax.set_ylabel("Largest power mismatch (p.u.)")
    ax.set_title("Normal case: convergence of each linear-solver backend", loc="left")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(FIGURES / "test_i_convergence.png", bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------------------------
# Part 2 - focused change: original vs fixed circuit, same configuration grid
# --------------------------------------------------------------------------------------------


def part2_original_vs_fixed() -> pd.DataFrame:
    case = paper_case5()
    reference = solve_newton_raphson(case, tolerance=1e-12, max_iterations=50)
    matrix, _ = fast_decoupled_matrices(case)
    direct = np.linalg.solve(matrix, REPRESENTATIVE_RHS)
    rows = []
    for calibration in CALIBRATIONS:
        for n in PHASE_QUBITS:
            row: dict[str, object] = {"calibration": calibration, "phase_qubits": n}
            for variant in ("spectral_model", "original_circuit", "circuit"):
                linear_factory = factory(calibration, n, variant)
                solved = linear_factory.prepare(matrix).solve(REPRESENTATIVE_RHS).solution
                row[f"{variant}_linear_error"] = float(
                    np.linalg.norm(solved - direct) / np.linalg.norm(direct)
                )
                result = run(
                    case, linear_factory, tolerance=TOLERANCE, max_iterations=MAX_ITERATIONS
                )
                row[f"{variant}_rounds"] = (
                    result.iterations if result is not None and result.converged else None
                )
                row[f"{variant}_state_error"] = (
                    state_error(result, reference) if result is not None else None
                )
            rows.append(row)
            print(f"  part2 {calibration} n={n} done", flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(TABLES / "original_vs_fixed_circuit.csv", index=False)
    return frame


# --------------------------------------------------------------------------------------------
# Part 3 - precision x loading sweep on the fixed, executed circuit
# --------------------------------------------------------------------------------------------


def part3_sweep(resources: pd.DataFrame) -> pd.DataFrame:
    cases = {f"x{m:g}": demand_case(m) for m in DEMAND_MULTIPLIERS}
    cases["Test II (stressed)"] = stressed_paper_case5()
    rows = []
    for case_label, case in cases.items():
        multiplier = float(case_label[1:]) if case_label.startswith("x") else np.nan
        reference = solve_newton_raphson(case, tolerance=1e-12, max_iterations=100)
        classical = run(case, None, tolerance=TOLERANCE, max_iterations=MAX_ITERATIONS)
        rows.append(
            {
                "case": case_label,
                "demand_multiplier": multiplier,
                "min_voltage_pu": float(np.min(reference.voltage_magnitude)),
                "backend": "classical",
                "calibration": None,
                "phase_qubits": 0,
                "converged": bool(classical and classical.converged),
                "rounds": classical.iterations if classical else None,
                "state_error": state_error(classical, reference) if classical else None,
                "state_error_after_6": state_error(
                    run(case, None, fixed_iterations=PAPER_ROUNDS), reference
                ),
            }
        )
        for calibration in CALIBRATIONS:
            for n in PHASE_QUBITS:
                linear_factory = factory(calibration, n)
                result = run(
                    case, linear_factory, tolerance=TOLERANCE, max_iterations=MAX_ITERATIONS
                )
                six = run(case, linear_factory, fixed_iterations=PAPER_ROUNDS)
                converged = bool(result is not None and result.converged)
                stats = solve_statistics(result) if result is not None else {}
                rows.append(
                    {
                        "case": case_label,
                        "demand_multiplier": multiplier,
                        "min_voltage_pu": float(np.min(reference.voltage_magnitude)),
                        "backend": "hhl_circuit",
                        "calibration": calibration,
                        "phase_qubits": n,
                        "converged": converged,
                        "rounds": result.iterations if converged else None,
                        "state_error": state_error(result, reference) if result else None,
                        "state_error_after_6": state_error(six, reference) if six else None,
                        **stats,
                    }
                )
        print(f"  part3 {case_label} done", flush=True)
    sweep = pd.DataFrame(rows)
    merged = sweep.merge(
        resources[["calibration", "phase_qubits", "logical_qubits", "depth", "cx_count"]],
        on=["calibration", "phase_qubits"],
        how="left",
    )
    # Expected two-qubit-gate executions for the whole power flow (one success per solve).
    merged["expected_total_cx"] = merged["expected_circuit_runs"] * merged["cx_count"]
    merged["meets_target"] = merged["converged"] & (merged["state_error"] <= ACCURACY_TARGET)
    merged["meets_target_in_6_rounds"] = merged["state_error_after_6"] <= ACCURACY_TARGET
    merged.to_csv(TABLES / "precision_loading_sweep.csv", index=False)
    return merged


def cheapest_settings(sweep: pd.DataFrame) -> pd.DataFrame:
    hhl = sweep[sweep["backend"] == "hhl_circuit"]
    rows = []
    for case_label, group in hhl.groupby("case", sort=False):
        for objective, column, mask in (
            ("converge to 1e-6, error <= 1e-4", "expected_total_cx", group["meets_target"]),
            ("error <= 1e-4 within 6 rounds", "cx_count", group["meets_target_in_6_rounds"]),
        ):
            for calibration in CALIBRATIONS:
                candidates = group[mask & (group["calibration"] == calibration)]
                if candidates.empty:
                    rows.append(
                        {
                            "case": case_label,
                            "objective": objective,
                            "calibration": calibration,
                            "phase_qubits": None,
                        }
                    )
                    continue
                best = candidates.sort_values([column, "phase_qubits"]).iloc[0]
                rows.append(
                    {
                        "case": case_label,
                        "objective": objective,
                        "calibration": calibration,
                        "phase_qubits": int(best["phase_qubits"]),
                        "logical_qubits": int(best["logical_qubits"]),
                        "cx_per_circuit": int(best["cx_count"]),
                        "depth": int(best["depth"]),
                        "rounds": best["rounds"],
                        "state_error": best["state_error"],
                        "expected_total_cx": best["expected_total_cx"],
                    }
                )
    frame = pd.DataFrame(rows)
    frame.to_csv(TABLES / "cheapest_settings_meeting_target.csv", index=False)
    return frame


def plot_precision(resources: pd.DataFrame, comparison: pd.DataFrame, sweep: pd.DataFrame) -> None:
    nominal = sweep[(sweep["case"] == "x1") & (sweep["backend"] == "hhl_circuit")]
    classical_rounds = int(
        sweep[(sweep["case"] == "x1") & (sweep["backend"] == "classical")]["rounds"].iloc[0]
    )
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9))
    calibration_styles = {
        "spectral_fit": dict(
            color=ORANGE, marker="o", label="Exact-bin calibration (needs eigenvalues)"
        ),
        "power_of_two": dict(
            color=BLUE, marker="s", label="Power-of-two scale (needs only λmax bound)"
        ),
    }
    ax = axes[0]
    for calibration, style in calibration_styles.items():
        rows = comparison[comparison["calibration"] == calibration]
        ax.semilogy(rows["phase_qubits"], rows["circuit_linear_error"].clip(lower=1e-16), **style)
        ax.semilogy(
            rows["phase_qubits"],
            rows["spectral_model_linear_error"].clip(lower=1e-16),
            color=style["color"],
            linestyle=":",
            linewidth=1.2,
            marker=None,
        )
    ax.set_title("One linear solve: error vs direct", loc="left")
    ax.set_xlabel("Phase qubits")
    ax.set_ylabel("Relative error (dotted: spectral model)")

    ax = axes[1]
    for calibration, style in calibration_styles.items():
        rows = nominal[nominal["calibration"] == calibration]
        ax.plot(rows["phase_qubits"], rows["rounds"], **style)
    ax.axhline(classical_rounds, color=MUTED, linewidth=1, linestyle=(0, (4, 3)))
    ax.text(
        1, classical_rounds + 0.4, f"classical FDLF: {classical_rounds}", color=MUTED, fontsize=8
    )
    ax.set_title("Normal demand: rounds to 1e-6 mismatch", loc="left")
    ax.set_xlabel("Phase qubits")
    ax.set_ylabel("Adjustment rounds")

    ax = axes[2]
    rows = resources[resources["calibration"] == "power_of_two"]
    ax.semilogy(rows["phase_qubits"], rows["cx_count"], color=INK, marker="o", label="CX gates")
    ax.semilogy(
        rows["phase_qubits"], rows["depth"], color=MUTED, marker="s", linestyle="--", label="Depth"
    )
    for _, row in rows.iterrows():
        ax.annotate(
            f"{int(row['logical_qubits'])}q",
            (row["phase_qubits"], row["cx_count"]),
            textcoords="offset points",
            xytext=(-14, 4),
            fontsize=7,
            color=MUTED,
        )
    ax.set_title("Circuit cost (basis cx/rz/sx/x)", loc="left")
    ax.set_xlabel("Phase qubits")
    ax.set_ylabel("Count (log scale)")
    ax.legend(loc="lower right")

    handles = [axes[0].lines[0], axes[0].lines[2]]
    fig.legend(
        handles,
        [h.get_label() for h in handles],
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.36, -0.03),
    )
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(FIGURES / "precision_tradeoff.png", bbox_inches="tight")
    plt.close(fig)


def plot_loading(sweep: pd.DataFrame) -> None:
    demand = sweep[sweep["demand_multiplier"].notna()]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=False)
    for ax, calibration, title in (
        (axes[0], "power_of_two", "Power-of-two scale: rounds to converge"),
        (axes[1], "spectral_fit", "Exact-bin calibration: rounds to converge"),
    ):
        classical = demand[demand["backend"] == "classical"]
        ax.plot(
            classical["demand_multiplier"],
            classical["rounds"],
            color=INK,
            linewidth=2.5,
            label="Classical FDLF",
        )
        for n in PHASE_QUBITS:
            rows = demand[(demand["calibration"] == calibration) & (demand["phase_qubits"] == n)]
            ax.plot(
                rows["demand_multiplier"],
                rows["rounds"],
                color=N_RAMP[n - 1],
                marker="o",
                markersize=4,
                label=f"{n} phase qubits",
            )
            failed = rows[~rows["converged"]]
            ax.plot(
                failed["demand_multiplier"],
                np.full(len(failed), MAX_ITERATIONS),
                "x",
                color=N_RAMP[n - 1],
                markersize=7,
            )
        ax.set_yscale("log")
        ax.set_title(title, loc="left")
        ax.set_xlabel("Demand multiplier (consuming buses)")
        ax.set_ylabel("Adjustment rounds (log); × = did not converge")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=7, fontsize=8, bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout()
    fig.savefig(FIGURES / "loading_convergence.png", bbox_inches="tight")
    plt.close(fig)


def plot_cost(sweep: pd.DataFrame) -> None:
    hhl = sweep[(sweep["backend"] == "hhl_circuit") & sweep["converged"]]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    for ax, calibration, title in (
        (axes[0], "power_of_two", "Power-of-two scale"),
        (axes[1], "spectral_fit", "Exact-bin calibration"),
    ):
        for n in PHASE_QUBITS:
            rows = hhl[(hhl["calibration"] == calibration) & (hhl["phase_qubits"] == n)]
            rows = rows[rows["demand_multiplier"].notna()]
            ax.plot(
                rows["demand_multiplier"],
                rows["expected_total_cx"],
                color=N_RAMP[n - 1],
                marker="o",
                markersize=4,
                label=f"{n} phase qubits",
            )
        ax.set_yscale("log")
        ax.set_title(f"{title}: expected CX executions per power flow", loc="left")
        ax.set_xlabel("Demand multiplier (consuming buses)")
    axes[0].set_ylabel("CX gates × circuit runs (log)")
    axes[1].legend(loc="upper left", ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "total_quantum_cost.png", bbox_inches="tight")
    plt.close(fig)


def environment() -> dict[str, str]:
    packages = ("qiskit", "qiskit-aer", "numpy", "scipy", "pandas", "matplotlib", "PyYAML")
    env = {"python": platform.python_version(), "platform": platform.platform()}
    env.update({name: version(name) for name in packages})
    return env


def replot() -> None:
    """Redraw the sweep figures from the saved tables (``--plots-only``)."""

    resources = pd.read_csv(TABLES / "circuit_resources.csv")
    comparison = pd.read_csv(TABLES / "original_vs_fixed_circuit.csv")
    sweep = pd.read_csv(TABLES / "precision_loading_sweep.csv")
    plot_precision(resources, comparison, sweep)
    plot_loading(sweep)
    plot_cost(sweep)


def main() -> None:
    for path in (TABLES, FIGURES):
        path.mkdir(parents=True, exist_ok=True)
    if "--plots-only" in sys.argv:
        replot()
        return
    env = environment()
    (OUT / "environment.json").write_text(json.dumps(env, indent=2))
    print("environment:", env, flush=True)

    print("Part 1: normal case", flush=True)
    part1 = part1_reproduction()
    print(part1["table"].to_string(), flush=True)

    print("Circuit resources", flush=True)
    resources = pd.DataFrame(
        [circuit_resources(calibration, n) for calibration in CALIBRATIONS for n in PHASE_QUBITS]
    )
    resources.to_csv(TABLES / "circuit_resources.csv", index=False)
    print(resources.to_string(), flush=True)

    print("Part 2: original vs fixed circuit", flush=True)
    comparison = part2_original_vs_fixed()
    print(comparison.to_string(), flush=True)

    print("Part 3: precision x loading sweep", flush=True)
    sweep = part3_sweep(resources)
    best = cheapest_settings(sweep)
    print(best.to_string(), flush=True)

    plot_precision(resources, comparison, sweep)
    plot_loading(sweep)
    plot_cost(sweep)
    print(f"Wrote outputs to {OUT}", flush=True)


if __name__ == "__main__":
    main()
