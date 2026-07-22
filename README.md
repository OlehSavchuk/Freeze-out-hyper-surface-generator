# Freeze-out-hyper-surface-generator
This code reads in coarse-graining conserved charges from transport simulations and constructs a freeze-out ypersurfaces (such as freeze-out surfaces defined by a specific density or temperature) from 4D space-time grids in relativistic heavy-ion collision simulations. 

## 📖 Overview
This tool reads macroscopic space-time evolution data (e.g., coarse-grained transport models or hydrodynamics) and extracts the 3D surface elements $d^3\sigma_\mu$ and the interpolated fluid properties (energy density, baryon density, flow velocity) across that surface. 

It utilizes the **Cornelius algorithm** to divide 4D hypercubes into continuous polyhedrons, avoiding the ambiguities of standard marching cubes algorithms in higher dimensions.

## ⚙️ Prerequisites
* A C++ compiler supporting **C++17** (e.g., GCC or Clang)
* **CMake** 3.15 or higher

## 🚀 Build Instructions
This project uses CMake for out-of-source builds. From the root directory of the repository, run:

```bash
mkdir build
cd build
cmake ..
make -j4
