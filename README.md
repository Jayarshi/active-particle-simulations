# Active Particle Simulations

A collection of Fortran and Python implementations of standard active matter models, including Active Brownian Particles (ABP), Run-and-Tumble Particles (RTP), and the Vicsek model. The repository is designed for both educational and research purposes, providing reference implementations, high-performance Fortran codes, Python notebooks for visualization, and F2PY interfaces for seamless integration between Fortran and Python.

## Models Included

### Active Brownian Particles (ABP)

Active particles undergoing self-propelled motion with continuous rotational diffusion.

Key topics:

* Persistent random walks
* Rotational diffusion
* Active transport
* Nonequilibrium steady states

### Run-and-Tumble Particles (RTP)

Particles that alternate between ballistic motion and stochastic tumbling events.

Key topics:

* Bacterial motility
* Persistence and reorientation
* Effective diffusion
* Nonequilibrium dynamics

### Vicsek Model

A paradigmatic model of collective motion in active matter.

Key topics:

* Flocking transitions
* Collective alignment
* Emergent order
* Nonequilibrium phase transitions

---

## Repository Structure

```text
active-particle-simulations/
├── benchmarks/
├── figures/
├── fortran/
│   ├── abpcode.f90
│   ├── rtpcode.f90
│   └── vicsekcode.f90
├── python/
│   ├── abp.ipynb
│   ├── rtp.ipynb
│   └── vicsek.ipynb
├── f2py_interfaces/
│   ├── abp.f90
│   ├── rtp.f90
│   ├── vicsek.f90
│   ├── abpf2py.ipynb
│   ├── rtpf2py.ipynb
│   └── vicsekf2py.ipynb
└── README.md
```

### `fortran/`

Standalone high-performance Fortran implementations suitable for large-scale simulations.

### `python/`

Interactive Jupyter notebooks demonstrating model implementation, simulation, and visualization.

### `f2py_interfaces/`

Examples of wrapping Fortran routines using NumPy F2PY for accelerated Python workflows.

### `benchmarks/`

Performance comparisons between pure Python and Fortran-accelerated implementations.

### `figures/`

Generated plots and simulation snapshots.

---

## Requirements

### Python

```bash
numpy
matplotlib
scipy
jupyter
```

Optional:

```bash
numba
```

### Fortran

```bash
gfortran
```

---

## Building F2PY Modules

Example:

```bash
f2py -c -m abp abp.f90
```

Similarly,

```bash
f2py -c -m rtp rtp.f90
f2py -c -m vicsek vicsek.f90
```

---

## Example Workflow

1. Develop and validate the model in Python.
2. Run large simulations using Fortran.
3. Expose computational kernels through F2PY.
4. Analyze and visualize results in Jupyter notebooks.

---

## Scientific Applications

* Active matter
* Statistical physics
* Nonequilibrium systems
* Collective behavior
* Biological locomotion
* Soft condensed matter

---

## References

1. M. C. Marchetti et al., *Hydrodynamics of Soft Active Matter*, Rev. Mod. Phys. (2013).
2. T. Vicsek et al., *Novel Type of Phase Transition in a System of Self-Driven Particles*, Phys. Rev. Lett. (1995).
3. C. Bechinger et al., *Active Particles in Complex and Crowded Environments*, Rev. Mod. Phys. (2016).

---

## License

This repository is intended for academic and educational use. Please cite the relevant literature when using these implementations in research.

