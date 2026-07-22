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

./hypersurface_builder ../config.txt
```

🙏 Acknowledgments & Citations

This repository relies heavily on previously developed algorithms and framework realizations within the heavy-ion physics community. If you use this code in a publication, you must cite the following works:

The Cornelius Algorithm
The core 4D hypersurface extraction algorithm (cornelius.cpp and cornelius.h) was developed by Pasi Huovinen and Hannu Holopainen.
Citation: P. Huovinen and H. Petersen, "Particlization in hybrid models", Eur. Phys. J. A 48, 171 (2012).
ArXiv: arXiv:1206.3371 [nucl-th]

MUSIC Hydrodynamics
The realization of the continuous 4D interpolation and hypercube grid mapping (hypersurface_builder implementation) is adapted from the MUSIC 3+1D relativistic viscous hydrodynamic framework.
Citation: B. Schenke, S. Jeon, and C. Gale, "3+1D hydrodynamic simulation of relativistic heavy-ion collisions", Phys. Rev. C 82, 014903 (2010).
GitHub Repository: MUSIC Collaboration
