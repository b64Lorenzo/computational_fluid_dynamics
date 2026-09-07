# Computational Fluid Dynamics

A collection of MATLAB implementations developed for Computational Fluid Dynamics (CFD) exercises.

## Project Structure

```text
src/
├── mesh/          Grid generation utilities
├── schemes/       Discretization schemes
├── solvers/       Numerical solvers
├── utils/         Shared utility functions
└── problems/      Individual problem implementations

tests/             Verification and validation scripts

docs/              Theory notes, assignments, and reports

results/           Generated figures, data, and log files
```

## Features

This projects contains methods shared by the course professor. Additional components for testing and running are:

- Logging utilities
- Automated result storage
- Modular and extensible architecture

## Testing

The `tests` directory contains scripts used to verify the different methods.

## Usage

1. Open the MATLAB project.
2. Run `startup.m` to add the source directories to the MATLAB path.
3. Navigate to the desired problem in `src/problems`.
4. Execute the corresponding `main.m` script.
5. Results are written to the `results` directory.

## Author

Lorenzo Zambelli