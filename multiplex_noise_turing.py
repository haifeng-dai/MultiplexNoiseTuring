"""Simulation tools for stochastic Turing patterns on multiplex networks.

The model is the Schnakenberg activator-inhibitor system coupled through the
aggregate Laplacian of one or more network layers.  Multiplicative noise is
integrated with the Euler-Maruyama scheme.
"""

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class ModelParameters:
    """Parameters for the Schnakenberg activator-inhibitor model."""

    a: float = 0.1
    b: float = 0.9
    activator_diffusion: float = 0.1
    inhibitor_diffusion: float = 1.0
    activator_noise: float = 0.0
    inhibitor_noise: float = 0.0


def multiplex_laplacian(
    layers: Sequence[np.ndarray], weights: Iterable[float] | None = None
) -> np.ndarray:
    """Return the weighted aggregate graph Laplacian for ``layers``.

    Layers are adjacency matrices with non-negative edge weights.  A
    symmetric matrix is required so that the returned Laplacian represents an
    undirected network.  The convention is ``L = A - diag(A.sum(axis=1))``,
    which is convenient for diffusion because it makes spatially uniform
    states stationary.
    """

    if not layers:
        raise ValueError("at least one network layer is required")
    arrays = [np.asarray(layer, dtype=float) for layer in layers]
    n = arrays[0].shape[0] if arrays[0].ndim == 2 else -1
    if n < 1:
        raise ValueError("each layer must be a non-empty square matrix")
    if weights is None:
        layer_weights = np.ones(len(arrays))
    else:
        layer_weights = np.asarray(list(weights), dtype=float)
        if layer_weights.shape != (len(arrays),):
            raise ValueError("weights must contain one value per layer")
        if np.any(layer_weights < 0):
            raise ValueError("layer weights must be non-negative")

    aggregate = np.zeros((n, n), dtype=float)
    for adjacency, weight in zip(arrays, layer_weights):
        if adjacency.shape != (n, n):
            raise ValueError("all layers must have the same square shape")
        if not np.all(np.isfinite(adjacency)):
            raise ValueError("network layers must contain finite values")
        if np.any(adjacency < 0) or not np.allclose(adjacency, adjacency.T):
            raise ValueError("layers must be symmetric and non-negative")
        if not np.allclose(np.diag(adjacency), 0):
            raise ValueError("self-loops are not supported")
        aggregate += weight * adjacency
    return aggregate - np.diag(aggregate.sum(axis=1))


def simulate(
    initial_state: np.ndarray,
    layers: Sequence[np.ndarray],
    *,
    steps: int,
    dt: float,
    parameters: ModelParameters | None = None,
    layer_weights: Iterable[float] | None = None,
    seed: int | np.random.Generator | None = None,
) -> np.ndarray:
    """Simulate the multiplex stochastic reaction-diffusion system.

    ``initial_state`` has shape ``(nodes, 2)`` and stores activator and
    inhibitor concentrations.  The returned array has shape
    ``(steps + 1, nodes, 2)`` and includes the initial state.  A supplied
    integer seed makes the trajectory reproducible.
    """

    if not isinstance(steps, (int, np.integer)) or steps < 0:
        raise ValueError("steps must be a non-negative integer")
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("dt must be positive and finite")
    state = np.asarray(initial_state, dtype=float)
    if state.ndim != 2 or state.shape[1] != 2 or state.shape[0] == 0:
        raise ValueError("initial_state must have shape (nodes, 2)")
    if not np.all(np.isfinite(state)):
        raise ValueError("initial_state must contain finite values")

    params = parameters or ModelParameters()
    if any(
        not np.isfinite(value) or value < 0
        for value in (
            params.activator_diffusion,
            params.inhibitor_diffusion,
            params.activator_noise,
            params.inhibitor_noise,
        )
    ):
        raise ValueError("diffusion and noise strengths must be finite and non-negative")
    laplacian = multiplex_laplacian(layers, layer_weights)
    if laplacian.shape[0] != state.shape[0]:
        raise ValueError("network layers and initial_state must have the same node count")

    rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
    trajectory = np.empty((steps + 1, *state.shape), dtype=float)
    trajectory[0] = state
    for index in range(steps):
        activator, inhibitor = state.T
        reaction = activator * activator / np.maximum(inhibitor, np.finfo(float).eps)
        drift_u = params.a - activator + reaction
        drift_v = params.b - reaction
        drift_u += params.activator_diffusion * (laplacian @ activator)
        drift_v += params.inhibitor_diffusion * (laplacian @ inhibitor)
        noise = rng.normal(size=state.shape) * np.sqrt(dt)
        state = state + dt * np.column_stack((drift_u, drift_v))
        state += noise * np.array(
            [params.activator_noise * activator, params.inhibitor_noise * inhibitor]
        )
        trajectory[index + 1] = state
    return trajectory

