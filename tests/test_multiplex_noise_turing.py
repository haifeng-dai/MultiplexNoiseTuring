import numpy as np
import pytest

from multiplex_noise_turing import ModelParameters, multiplex_laplacian, simulate


def test_multiplex_laplacian_combines_layers():
    first = np.array([[0, 1], [1, 0]], dtype=float)
    second = np.array([[0, 2], [2, 0]], dtype=float)

    np.testing.assert_allclose(
        multiplex_laplacian([first, second], weights=[1, 0.5]),
        [[-2, 2], [2, -2]],
    )


def test_zero_noise_simulation_is_reproducible_and_returns_initial_state():
    state = np.array([[1.0, 1.0], [1.1, 0.9]])
    layer = np.array([[0.0, 1.0], [1.0, 0.0]])
    parameters = ModelParameters(activator_noise=0, inhibitor_noise=0)

    first = simulate(state, [layer], steps=3, dt=0.01, parameters=parameters, seed=1)
    second = simulate(state, [layer], steps=3, dt=0.01, parameters=parameters, seed=2)

    np.testing.assert_array_equal(first, second)
    np.testing.assert_array_equal(first[0], state)


def test_seed_controls_multiplicative_noise():
    state = np.ones((2, 2))
    layer = np.zeros((2, 2))
    parameters = ModelParameters(activator_noise=0.2, inhibitor_noise=0.3)

    first = simulate(state, [layer], steps=2, dt=0.01, parameters=parameters, seed=7)
    second = simulate(state, [layer], steps=2, dt=0.01, parameters=parameters, seed=7)

    np.testing.assert_array_equal(first, second)


def test_invalid_layer_is_rejected():
    with pytest.raises(ValueError, match="symmetric"):
        multiplex_laplacian([[[0, 1], [0, 0]]])
