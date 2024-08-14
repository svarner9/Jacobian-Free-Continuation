# Jacobian-Free Pseudo-arclength Continuation (PAC)

This repository contains numerical routines and example problems designed to solve various nonlinear equations and models encountered in statistical mechanics and applied mathematics. The examples provided demonstrate our Jacobian-free pseudo-arclength continuation technique.

If you find, this repository helpful, please cite out paper! [Varner, Balzer, Wang, _JCP_ (2024)](https://doi.org/10.1063/5.0220849)

## Folder Structure

The repository is organized into three main folders, each containing specific examples and their respective codes:

### 1. Bratu

The `Bratu` folder contains code to solve the Bratu Equation (Liouville–Bratu–Gelfand) in 2D using Dirichlet boundary conditions. This equation models the temperature distribution in a reactor with exothermic reactions. The folder includes:

- **Anderson**: Implementation of an Anderson-based Pseudo-Arclength Continuation (PAC) algorithm.
- **Newton**: Implementation of a traditional Newton-based PAC algorithm.

### 2. Ising

The `Ising` folder houses the implementation of the Ising model with an external field under the mean-field approximation. This model is used to study phase transitions and magnetic properties of materials.

### 3. PolymerAdsorption

The `PolymerAdsorption` folder contains code to calculate the polymer density profile of a Flory solution near a solid surface. This model considers the interaction between polymers and solvents via the Flory-Huggins parameter ($\chi$).

## How to Run the Examples

Each example folder includes a `README.md` file with detailed instructions on how to run the code. Generally, to execute the code, navigate to the respective folder and run the main Python script:

```sh
$ cd [ExampleFolder]
$ python main.py
