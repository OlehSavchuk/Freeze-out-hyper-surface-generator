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

## Local momentum anisotropy (branches `local-momentum-anisotropy`, `full-momentum-anisotropy`)
The surface output can also carry the local momentum anisotropy of the coarse-grained cells,
`exx = 3Txx/(Txx+Tyy+Tzz)`, `eyy`, `ezz` (each is 1 for an isotropic cell and `exx + eyy + ezz = 3`), extracted on the same grid
as the other fields. Put the binary files `exx.dat`, `eyy.dat`, `ezz.dat` (same layout as `e.dat`, `n.dat`: doubles, `[t][x][y][z]`)
into the data `path` and add to `config.txt`:

```
anisotropy 1        # 0 (default): output unchanged; 1: diagonal exx, eyy, ezz; 2: full tensor (also exy, exz, eyz)
aniso_fill 1.0      # value for cells without matter (NaN in the files), default 1 = isotropic
```

With `anisotropy 2` the files `exy.dat`, `exz.dat`, `eyz.dat` (`e_ij = 3<p_i p_j>/sum_k <p_k^2>`, 0 for an isotropic cell) are read as well. An optional `outfile name` sets the name of the surface file (default `surface.dat`). The time slices are streamed from the files, so only two time slices are kept in memory.

The three fields are linearly interpolated to the centroid of every surface element exactly like `e, n, gamma, vx, vy, vz`
(because the interpolation is linear, `exx + eyy + ezz = 3` is preserved). Each row of `surface.dat` then has 18 columns instead of 15 (21 columns with `anisotropy 2`):

| columns | content |
|---|---|
| 0-3 | surface element `d^3 sigma_mu` |
| 4-7 | `t, x, y, z` of the element |
| 8-14 | `e, T, n, gamma, vx, vy, vz` |
| 15-17 | `exx, eyy, ezz` (with `anisotropy 1` or `2`) |
| 18-20 | `exy, exz, eyz` (only with `anisotropy 2`) |

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
