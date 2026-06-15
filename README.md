# Hydrogen and Water Energy Simulation Using XBK Transformation and Coherent Ising Machine Variations (CAC, CFC, SFC) and Simulated Bifurcation variant (dSB)

![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![License](https://img.shields.io/badge/License-See_Credits-orange)
![Status](https://img.shields.io/badge/Status-Research-green)

This repository contains a complete hybrid quantum-classical workflow to calculate the ground-state energy of the Hydrogen molecule (H2) and the Water (H2O) Molecules across dissociation curves. It integrates PySCF and OpenFermion for quantum chemistry with CIM variants and dSB solver for Ising optimization, utilizing the Extended Bravyi-Kitaev (XBK) transformation in STO-3G basis.

# TABLE OF CONTENTS
1. Overview
2. Key Features
3. Dependencies & Installation
4. Usage
5. Methodology
6. Outputs
7. Licenses, Credits & References


# 1. OVERVIEW
This project implements a simulation to compute the ground-state energy of H2 and H2O by mapping the electronic-structure problem to a Quadratic Unconstrained Binary Optimization (QUBO) problem.

The workflow consists of:
1. Hamiltonian Generation: Using PySCF and OpenFermion.
2. Symmetry Reduction: Applying Parity Tapering to reduce qubit count.
3. XBK Transformation: Transforming the problem using auxiliary qubits (r parameter) to map spectral properties to an Ising model.
4. CIM and SB Simulation: Solving the resulting Ising model using a digital simulation of CIM and SB variants.

# 2. KEY FEATURES
* Gauss-Jordan Elimination Engine: A robust implementation of linear algebra over finite fields (Rational, Prime, Binary) used for analyzing Hamiltonian symmetries.
* Coherent Ising Machine (CIM): A dynamical optimization approach inspired by networks of optical parametric oscillators, designed to solve Ising spin optimization problems by evolving toward low-energy spin configurations.
* Simulated Bifurcation (SB): A physics-inspired optimizer that uses bifurcation dynamics of coupled oscillators to search for low-energy Ising states.
* Chaotic Amplitude Control (CAC): A dynamical systems algorithm that mimics optical pulses in a Coherent Ising Machine to find global minima of spin glass problems.
* * Separated Feedback Control (SFC): A Coherent Ising Machine variant that employs separate feedback channels to stabilize dynamics and improve optimization performance.
* Chaotic Feedback Control (CFC): A Coherent Ising Machine variant that incorporates nonlinear chaotic feedback to enhance exploration of complex energy landscapes.
* Discrete Simulated Bifurcation (dSB): A discrete-time implementation of the Simulated Bifurcation algorithm for solving Ising optimization problems through oscillator-inspired dynamics.
* XBK Resolution Sweeping: Iterates through resolution parameters (r=2 to r=6) to demonstrate convergence accuracy.
* Automated Visualization: Plots the Potential Energy Surface (PES) of the molecule immediately after calculation.

# 3. DEPENDENCIES & INSTALLATION
This project requires a Python environment (3.8+) with specific scientific libraries.

Recommended Operating System: Linux or macOS (due to 'pyscf' availability).

REQUIRED LIBRARIES:
* Quantum Chemistry: qiskit, qiskit-nature, openfermion, openfermionpyscf, pyscf
* Optimization: dimod, dwave-ocean-sdk
* Math & Physics: numpy, scipy, symengine, numba, torch
* Data & Plotting: pandas, matplotlib

# INSTALLATION COMMANDS:
--------------------------------------------------------------------------------
pip install numpy scipy matplotlib pandas symengine numba torch
pip install qiskit
pip install qiskit-algorithms
pip install qiskit-nature
pip install dwave-ocean-sdk
pip install openfermion openfermionpyscf pyscf
--------------------------------------------------------------------------------

# 4. USAGE

For each variant, run the relevant file

PROCESS:
* The script "Benchmark for H2" calculates the energy for bond lengths ranging from 0.2 Å to 4.0 Å for four variants of CIM and SB.
* First of all, please install all packages and libraries in "quantum_utils.py", and after that, the linear_algebra scripts in the "linear-algebra.py" file.
* For every variant, please first run the "Creating QAIA Class.py", then "Helper Functions.py", and after that, "Applying XBK transformation and variant function.py", then "Achieving Energy.py", and at the end "Plotting and Saving the Result.py".
* It loops through XBK resolutions r=2, 3, 4, 5, 6.
* NOTE: This is computationally intensive and may take time to complete.

# 5. METHODOLOGY

[XBK Transformation]
The Extended Bravyi-Kitaev (XBK) transformation is a method to map the eigenvalue problem of a Hamiltonian into a classical optimization problem. It uses auxiliary qubits (controlled here by the parameter r) to approximate the spectral decomposition, allowing Ising machines to estimate quantum ground states.

[Coherent Ising Machine (CIM)]
CIM is a quantum-inspired algorithm for solving combinatorial optimization problems. It introduces chaotic feedback loops into the dynamics of the spins (represented as continuous variables), allowing the system to escape local minima more effectively than standard gradient descent.

# 6. OUTPUTS

The script generates the following artifacts in the working directory:


1. Visualization:
   - A Matplotlib window displaying "Energy vs. Bond Length for ".
   - Curves for all r values are plotted for comparison.


# 7. LICENSES, CREDITS & REFERENCES
This project incorporates code from several open-source and research projects.

A. FIELD/MATRIX CLASS (Gauss-Jordan Elimination)
   - Source: Project Nayuki (https://www.nayuki.io/page/gauss-jordan-elimination-over-any-field)
   - Copyright (c) 2022 Project Nayuki. All rights reserved.
   - Note: Included for functional demonstration. Contact Project Nayuki for licensing details.

B. QUANTUM CHEMISTRY HELPER FUNCTIONS
   - Source: Quantum-Chemistry-with-Annealers (GitHub)
   - Ref: https://github.com/jcopenh/Quantum-Chemistry-with-Annealers.git
   - Usage: Functions for Pauli operator construction and molecular geometry.

C. CAC ALGORITHM (Coherent Ising Machine)
   - Reference: "Coherent Ising machines with optical error correction circuits"
   - Link: https://onlinelibrary.wiley.com/doi/full/10.1002/qute.202100077

D. LIBRARIES
   - OpenFermion (Apache 2.0)
   - Qiskit (Apache 2.0)
   - PySCF (Apache 2.0 / MIT)
