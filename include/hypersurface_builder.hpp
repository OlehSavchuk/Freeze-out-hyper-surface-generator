#pragma once
#include <cornelius.h>
#include <cmath>
#include <vector>

class HyperSurfaceBuilder {
public:
    // Constructor initializes grid sizes and the Cornelius algorithm parameters[cite: 5]
    HyperSurfaceBuilder(int nt, int nx, int ny, int nz, double freezeout_N, double dt, double dx);
    
    // Performs a 4-dimensional linear interpolation for the fluid quantities to find precise surface values[cite: 5]
    void fD_Interpolate(std::vector<double> shift, double to_interpolate[7][2][2][2][2], std::vector<double>& interpolated);
    
    // Maps local grid data into a hypercube format and triggers the Cornelius surface extraction[cite: 5]
    void find_elements(int i, int j, int k, int l, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& Tnv);

private:
    std::ofstream s_file;           // Output file stream for surface data[cite: 5]
    enum columns{ e,T,n,gamma,vx,vy,vz }; // Mapping for the 7 fluid variables[cite: 5]
    
    // Grid dimensions and freeze-out criteria[cite: 5]
    int nt_,nx_, ny_, nz_;
    double Nfreeze_;
    double DT, DX, DY, DZ;          // Lattice spacing in time and space[cite: 5]
    
    int dim;                        // Dimensionality (4D)[cite: 5]
    Cornelius* cornelius_ptr;       // Pointer to the Cornelius surface extraction algorithm[cite: 5]
};
