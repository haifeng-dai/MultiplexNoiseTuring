# Impact of Noise on Turing Patterns of Multiplex Networks

MATLAB simulations and figure scripts for the paper *Impact of Noise on Turing Patterns of Multiplex Networks*.

## Contents

- `fig2.m`–`fig6.m`: scripts for Figures 2–6.
- `Figure1.pdf`: network schematic for Figure 1.
- `networks/` and `generate_standard_networks.m`: network generation.
- `simulations/`: simulation and threshold calculation functions.

## Running the code

MATLAB R2021a or later is required. Scripts read network topologies, initial states, and simulation results from `results/` and write figures to `fig/`. MAT-files are not included and must be generated before plotting.

The network generator and simulation parameters require configuration: ER and WS generation are disabled by default, and some figure scripts default to reading cached results. The complete workflow from generated networks to all figures has not been validated.
