"""Finite-register HHL model for the fixed fast-decoupled QPF matrices.

The implementation uses Qiskit to construct the QPE / reciprocal-rotation circuit. By default
(``execution="spectral_model"``) linear solves use the idealized inverse map
``V diag(1/lambda_bin) V^T``, which avoids repeatedly simulating the same circuit in the
5,000-sample experiment. With ``execution="statevector"`` every solve instead evolves the
right-hand side through the constructed circuit, postselects, and reads out the solution, so the
result includes finite-register QPE leakage and any circuit-construction error.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import ceil, log2
from typing import Any, Literal

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, transpile
from qiskit.circuit.library import RYGate, StatePreparation, UnitaryGate, phase_estimation
from qiskit.quantum_info import Operator
from qiskit_aer import AerSimulator
from scipy.linalg import expm

from ..linear import LinearSolveResult
from ..models import FloatArray

ScaleStrategy = Literal["spectral_fit", "power_of_two", "fixed"]
Execution = Literal["spectral_model", "statevector"]


@dataclass(frozen=True, slots=True)
class HHLConfiguration:
    phase_qubits: int = 4
    scale_strategy: ScaleStrategy = "spectral_fit"
    fixed_scale: float | None = None
    reciprocal_constant_ratio: float = 0.5
    transpile_optimization_level: int = 0
    execution: Execution = "spectral_model"

    def __post_init__(self) -> None:
        if self.phase_qubits < 1:
            raise ValueError("phase_qubits must be positive")
        if self.scale_strategy == "fixed" and (self.fixed_scale is None or self.fixed_scale <= 0):
            raise ValueError("fixed_scale must be positive for fixed scaling")
        if self.execution not in ("spectral_model", "statevector"):
            raise ValueError("execution must be 'spectral_model' or 'statevector'")


def _reverse_bits(value: int, width: int) -> int:
    return int(format(value, f"0{width}b")[::-1], 2)


@dataclass(frozen=True, slots=True)
class SpectralCalibration:
    eigenvalues: FloatArray
    phase_scale: float
    logical_phase_bins: np.ndarray
    reconstructed_eigenvalues: FloatArray
    maximum_relative_eigenvalue_error: float
    condition_number: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "eigenvalues": self.eigenvalues.tolist(),
            "phase_scale": self.phase_scale,
            "logical_phase_bins": self.logical_phase_bins.tolist(),
            "reconstructed_eigenvalues": self.reconstructed_eigenvalues.tolist(),
            "maximum_relative_eigenvalue_error": self.maximum_relative_eigenvalue_error,
            "condition_number": self.condition_number,
        }


def _spectral_fit_scale(eigenvalues: FloatArray, phase_qubits: int) -> float:
    resolution = 1 << phase_qubits
    best: tuple[tuple[float, int], float] | None = None
    for eigenvalue in np.unique(np.round(eigenvalues, 14)):
        for phase_bin in range(1, resolution):
            scale = float(eigenvalue * resolution / phase_bin)
            phases = eigenvalues / scale
            if np.max(phases) >= 1.0 - 1e-13:
                continue
            bins = np.clip(np.rint(phases * resolution).astype(int), 1, resolution - 1)
            reconstructed = bins * scale / resolution
            error = float(np.max(np.abs(reconstructed - eigenvalues) / eigenvalues))
            key = (0.0 if error < 1e-12 else error, -int(np.max(bins)))
            if best is None or key < best[0]:
                best = (key, scale)
    if best is None:
        raise RuntimeError("Could not calibrate the phase register")
    return best[1]


def _calibrate(eigenvalues: FloatArray, config: HHLConfiguration) -> SpectralCalibration:
    resolution = 1 << config.phase_qubits
    if config.scale_strategy == "spectral_fit":
        scale = _spectral_fit_scale(eigenvalues, config.phase_qubits)
    elif config.scale_strategy == "power_of_two":
        scale = float(2 ** ceil(log2(float(np.max(eigenvalues)))))
        if scale <= np.max(eigenvalues):
            scale *= 2.0
    else:
        assert config.fixed_scale is not None
        scale = float(config.fixed_scale)
    bins = np.clip(np.rint(eigenvalues / scale * resolution).astype(int), 1, resolution - 1)
    reconstructed = bins * scale / resolution
    rel_error = float(np.max(np.abs(reconstructed - eigenvalues) / eigenvalues))
    return SpectralCalibration(
        eigenvalues=eigenvalues,
        phase_scale=scale,
        logical_phase_bins=bins,
        reconstructed_eigenvalues=reconstructed,
        maximum_relative_eigenvalue_error=rel_error,
        condition_number=float(np.max(eigenvalues) / np.min(eigenvalues)),
    )


@dataclass(slots=True)
class QiskitHHLPreparedSystem:
    matrix: FloatArray
    label: str
    configuration: HHLConfiguration
    calibration: SpectralCalibration
    eigenvectors: FloatArray
    reciprocal_constant: float
    core_circuit: QuantumCircuit
    _inverse_map: FloatArray
    _stats: dict[str, Any] | None = None
    _core_unitary: np.ndarray | None = None

    @property
    def dimension(self) -> int:
        return self.matrix.shape[0]

    @property
    def system_qubits(self) -> int:
        return int(log2(self.dimension))

    @property
    def logical_qubits(self) -> int:
        return self.configuration.phase_qubits + self.system_qubits + 1

    @classmethod
    def build(
        cls, matrix: FloatArray, *, label: str, configuration: HHLConfiguration
    ) -> QiskitHHLPreparedSystem:
        matrix = np.asarray(matrix, dtype=float)
        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
            raise ValueError("matrix must be square")
        if matrix.shape[0] & (matrix.shape[0] - 1):
            raise ValueError("matrix dimension must be a power of two")
        if not np.allclose(matrix, matrix.T, atol=1e-12):
            raise ValueError("matrix must be real symmetric")
        eigenvalues, eigenvectors = np.linalg.eigh(matrix)
        if np.min(eigenvalues) <= 0:
            raise ValueError("matrix must be positive definite")
        calibration = _calibrate(eigenvalues, configuration)
        reciprocal_constant = configuration.reciprocal_constant_ratio * float(
            np.min(calibration.reconstructed_eigenvalues)
        )
        inverse_map = (
            eigenvectors @ np.diag(1.0 / calibration.reconstructed_eigenvalues) @ eigenvectors.T
        )

        phase = QuantumRegister(configuration.phase_qubits, "phase")
        system = QuantumRegister(int(log2(matrix.shape[0])), "system")
        ancilla = QuantumRegister(1, "ancilla")
        core = QuantumCircuit(phase, system, ancilla, name="HHL-core")
        unitary = UnitaryGate(expm(2j * np.pi * matrix / calibration.phase_scale), label="exp(iAt)")
        qpe = phase_estimation(configuration.phase_qubits, unitary)
        core.compose(qpe, list(phase) + list(system), inplace=True)
        resolution = 1 << configuration.phase_qubits
        for logical_bin in range(1, resolution):
            lambda_est = logical_bin * calibration.phase_scale / resolution
            ratio = min(1.0, reciprocal_constant / lambda_est)
            angle = 2.0 * np.arcsin(ratio)
            # Qiskit's ``phase_estimation`` leaves the phase register bit-reversed relative to the
            # integer bin (bin 5 = 0b0101 is read as 0b1010), so control on the reversed pattern.
            gate = RYGate(float(angle)).control(
                configuration.phase_qubits,
                ctrl_state=_reverse_bits(logical_bin, configuration.phase_qubits),
            )
            core.append(gate, list(phase) + [ancilla[0]])
        core.compose(qpe.inverse(), list(phase) + list(system), inplace=True)
        return cls(
            matrix,
            label,
            configuration,
            calibration,
            eigenvectors,
            reciprocal_constant,
            core,
            inverse_map,
        )

    def circuit_for_rhs(self, rhs: FloatArray) -> QuantumCircuit:
        rhs = np.asarray(rhs, dtype=float)
        norm = np.linalg.norm(rhs)
        if norm == 0:
            raise ValueError("rhs must be nonzero")
        circuit = QuantumCircuit(self.logical_qubits)
        system_start = self.configuration.phase_qubits
        circuit.append(
            StatePreparation(rhs / norm), range(system_start, system_start + self.system_qubits)
        )
        circuit.compose(self.core_circuit, range(self.logical_qubits), inplace=True)
        return circuit

    def execute_statevector(self, rhs: FloatArray) -> tuple[FloatArray, dict[str, float]]:
        """Run ``rhs`` through the HHL circuit and return the postselected, rescaled solution.

        The core circuit's unitary is computed once and cached. The phase register ``|0...0>``
        and ancilla ``|1>`` are postselected, the global phase is fixed so the overlap with ``b``
        is positive real, and the amplitudes are rescaled by ``|b| / C``.
        """

        rhs = np.asarray(rhs, dtype=float)
        norm = float(np.linalg.norm(rhs))
        if norm == 0:
            raise ValueError("rhs must be nonzero")
        if self._core_unitary is None:
            self._core_unitary = Operator(self.core_circuit).data
        phase_qubits = self.configuration.phase_qubits
        unit_rhs = rhs / norm
        initial = np.zeros(1 << self.logical_qubits, dtype=complex)
        system_indices = np.arange(self.dimension) << phase_qubits
        initial[system_indices] = unit_rhs
        final = self._core_unitary @ initial
        selected = final[(1 << (phase_qubits + self.system_qubits)) | system_indices]
        success_probability = float(np.sum(np.abs(selected) ** 2))
        overlap = np.vdot(unit_rhs, selected)
        if abs(overlap) > 0:
            selected = selected * np.conj(overlap) / abs(overlap)
        solution = selected.real * norm / self.reciprocal_constant
        return np.asarray(solution, dtype=float), {
            "postselection_probability_phase_zero_ancilla_one": success_probability,
            "imaginary_solution_norm_after_phase_fix": float(np.linalg.norm(selected.imag)),
        }

    def solve(self, rhs: FloatArray) -> LinearSolveResult:
        rhs = np.asarray(rhs, dtype=float)
        if rhs.shape != (self.dimension,):
            raise ValueError(f"rhs must have shape ({self.dimension},)")
        execution_metadata: dict[str, float] = {}
        if self.configuration.execution == "statevector":
            solution, execution_metadata = self.execute_statevector(rhs)
        else:
            solution = self._inverse_map @ rhs
        direct = np.linalg.solve(self.matrix, rhs)
        residual = float(np.linalg.norm(self.matrix @ solution - rhs))
        norm_solution = np.linalg.norm(solution)
        norm_direct = np.linalg.norm(direct)
        fidelity = (
            1.0
            if norm_solution == 0 or norm_direct == 0
            else float(abs(np.vdot(solution / norm_solution, direct / norm_direct)) ** 2)
        )
        return LinearSolveResult(
            solution=np.asarray(solution, dtype=float),
            residual_norm=residual,
            metadata={
                "backend": "qiskit_finite_qpe_hhl",
                "execution": self.configuration.execution,
                **execution_metadata,
                "solution_state_fidelity_vs_direct": fidelity,
                "relative_solution_error_vs_direct": float(
                    np.linalg.norm(solution - direct) / max(norm_direct, 1e-15)
                ),
                "calibration": self.calibration.to_dict(),
            },
        )

    def solve_many(self, right_hand_sides: FloatArray) -> FloatArray:
        rhs = np.asarray(right_hand_sides, dtype=float)
        if rhs.ndim != 2 or rhs.shape[1] != self.dimension:
            raise ValueError("right_hand_sides has an incompatible shape")
        if self.configuration.execution == "statevector":
            return np.vstack([self.execute_statevector(row)[0] for row in rhs])
        return rhs @ self._inverse_map.T

    def circuit_statistics(self) -> dict[str, Any]:
        if self._stats is None:
            circuit = self.circuit_for_rhs(np.arange(1, self.dimension + 1, dtype=float))
            backend = AerSimulator(method="statevector")
            compiled = transpile(
                circuit, backend, optimization_level=self.configuration.transpile_optimization_level
            )
            self._stats = {
                "logical_qubits": self.logical_qubits,
                "phase_qubits": self.configuration.phase_qubits,
                "system_qubits": self.system_qubits,
                "ancilla_qubits": 1,
                "raw_depth": circuit.depth(),
                "raw_operation_counts": {str(k): int(v) for k, v in circuit.count_ops().items()},
                "transpiled_depth": compiled.depth(),
                "transpiled_operation_counts": {
                    str(k): int(v) for k, v in compiled.count_ops().items()
                },
                "calibration": self.calibration.to_dict(),
                "reciprocal_constant": self.reciprocal_constant,
            }
        return dict(self._stats)


@dataclass(slots=True)
class QiskitHHLFactory:
    configuration: HHLConfiguration = field(default_factory=HHLConfiguration)
    name: str = "qiskit_hhl"
    _cache: dict[bytes, QiskitHHLPreparedSystem] = field(
        default_factory=dict, init=False, repr=False
    )
    cache_hits: int = 0
    cache_misses: int = 0

    def prepare(self, matrix: FloatArray, *, label: str = "A") -> QiskitHHLPreparedSystem:
        array = np.ascontiguousarray(matrix, dtype=float)
        key = array.shape[0].to_bytes(4, "little") + array.tobytes()
        if key in self._cache:
            self.cache_hits += 1
            return self._cache[key]
        prepared = QiskitHHLPreparedSystem.build(
            array, label=label, configuration=self.configuration
        )
        self._cache[key] = prepared
        self.cache_misses += 1
        return prepared
