# Hydrogen Energy Simulation using XBK Transformation and CAC

![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![License](https://img.shields.io/badge/License-See_Credits-orange)
![Status](https://img.shields.io/badge/Status-Research-green)

This repository contains a complete hybrid quantum-classical workflow to calculate the ground state energy of the Hydrogen molecule (H2) across dissociation curves. It integrates PySCF and OpenFermion for quantum chemistry with a custom Chaotic Amplitude Control (CAC) solver for Ising optimization, utilizing the Extended Bravyi-Kitaev (XBK) transformation.

================================================================================
TABLE OF CONTENTS
================================================================================
1. Overview
2. Key Features
3. Dependencies & Installation
4. Usage
5. Methodology
6. Outputs
7. Licenses, Credits & References

================================================================================
1. OVERVIEW
================================================================================
This project implements a simulation to find the ground state energy of H2 by mapping the electronic structure problem to a Quadratic Unconstrained Binary Optimization (QUBO) problem.

The workflow consists of:
1. Hamiltonian Generation: Using PySCF and OpenFermion.
2. Symmetry Reduction: Applying Parity Tapering to reduce qubit count.
3. XBK Transformation: Transforming the problem using auxiliary qubits (r parameter) to map spectral properties to an Ising model.
4. CAC Simulation: Solving the resulting Ising model using a digital simulation of a Coherent Ising Machine (CIM).

================================================================================
2. KEY FEATURES
================================================================================
* Gauss-Jordan Elimination Engine: A robust implementation of linear algebra over finite fields (Rational, Prime, Binary) used for analyzing Hamiltonian symmetries.
* Chaotic Amplitude Control (CAC): A dynamical systems algorithm that mimics optical pulses in a Coherent Ising Machine to find global minima of spin glass problems.
* XBK Resolution Sweeping: Iterates through resolution parameters (r=2 to r=6) to demonstrate convergence accuracy.
* Automated Visualization: Plots the Potential Energy Surface (PES) of the molecule immediately after calculation.

================================================================================
3. DEPENDENCIES & INSTALLATION
================================================================================
This project requires a Python environment (3.8+) with specific scientific libraries.

Recommended Operating System: Linux or macOS (due to 'pyscf' availability).

REQUIRED LIBRARIES:
* Quantum Chemistry: qiskit, qiskit-nature, openfermion, openfermionpyscf, pyscf
* Optimization: dimod, dwave-ocean-sdk
* Math & Physics: numpy, scipy, symengine, numba, torch
* Data & Plotting: pandas, matplotlib

INSTALLATION COMMANDS:
--------------------------------------------------------------------------------
pip install numpy scipy matplotlib pandas symengine numba torch
pip install qiskit qiskit-nature qiskit-algorithms
pip install openfermion openfermionpyscf pyscf
pip install dimod dwave-ocean-sdk
--------------------------------------------------------------------------------

================================================================================
4. USAGE
================================================================================
1. Save the provided Python code as 'main.py'.
2. Run the script:
   
   python main.py

PROCESS:
* The script calculates the energy for bond lengths ranging from 0.2A to 4.0A.
* It loops through XBK resolutions r=2, 3, 4, 5, 6.
* NOTE: This is computationally intensive and may take time to complete.

================================================================================
5. METHODOLOGY
================================================================================

[XBK Transformation]
The Extended Bravyi-Kitaev (XBK) transformation is a method to map the eigenvalue problem of a Hamiltonian into a classical optimization problem. It uses auxiliary qubits (controlled here by the parameter r) to approximate the spectral decomposition, allowing Ising machines to estimate quantum ground states.

[Chaotic Amplitude Control (CAC)]
CAC is an algorithm for solving combinatorial optimization problems. It introduces chaotic feedback loops into the dynamics of the spins (represented as continuous variables), allowing the system to escape local minima more effectively than standard gradient descent.

================================================================================
6. OUTPUTS
================================================================================
The script generates the following artifacts in the working directory:

1. Data Files (.npy):
   - SFCH2r=2.npy
   - SFCH2r=3.npy
   - SFCH2r=4.npy
   - SFCH2r=5.npy
   - SFCH2r=6.npy
   (These contain the computed energy arrays).

2. Visualization:
   - A Matplotlib window displaying "Energy vs. Bond Length for CAC".
   - Curves for all r values are plotted for comparison.

================================================================================
7. LICENSES, CREDITS & REFERENCES
================================================================================
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
