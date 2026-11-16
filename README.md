# MultiplexNoiseTuring

Simulation code for Turing patterns in a multiplex-networked
activator-inhibitor system with multiplicative noise.

## Usage

Install the runtime dependency with:

```shell
python -m pip install -r requirements.txt
```

The simulator uses a sequence of symmetric adjacency matrices, one per
network layer:

```python
import numpy as np
from multiplex_noise_turing import ModelParameters, simulate

layer = np.array([[0, 1], [1, 0]], dtype=float)
initial = np.ones((2, 2))  # columns are activator and inhibitor
trajectory = simulate(
    initial,
    [layer],
    steps=100,
    dt=0.01,
    parameters=ModelParameters(activator_noise=0.1, inhibitor_noise=0.05),
    seed=42,
)
```

`trajectory` has shape `(steps + 1, nodes, 2)` and includes the initial
state. The integer seed makes stochastic runs reproducible.