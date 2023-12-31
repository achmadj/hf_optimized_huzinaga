# Hartree-Fock with Optimized Huzinaga Gaussian Orbital Functions

Hartree-Fock Roothaan Calculations Using Optimized Huzinaga Orbitals on Small Molecules

This repository contains methods to optimize the Huzinaga Gaussian orbital functions and the code used to perform Hartree-Fock calculations using the optimized Huzinaga Gaussian orbital functions. The goal is to achieve the ground state energy of small diatomic molecules using optimized STO-6G Huzinaga basis set to be on par (or better) with larger basis set (e.g. 6-31G, cc-pVQZ). The PySCF library is utilized to perform the Hartree-Fock calculations.

## Installation

To install the required dependencies, including the PySCF library, you can use the following command:

```bash
pip install pyscf
```

## Publication
The code is developed to support and reproduce the results presented in the research paper titled "Hartree-Fock Roothaan Calculations Using Optimized Huzinaga Orbitals on Small Molecules".