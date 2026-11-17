#include "hypersurface_builder.hpp"
#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <map>
#include <sstream>

// --- HELPER FUNCTIONS ---

// Parses a simple key-value configuration file
std::map<std::string, std::string> load_parameters(const std::string& filename) {
    std::map<std::string, std::string> params;
    std::ifstream file(filename);
    std::string line, key, value;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue; // Skip comments and empty lines
        std::istringstream iss(line);
        if (iss >> key >> value) {
            params[key] = value;
        }
    }
    return params;
}

// Writes a 1D vector row to a binary output file stream[cite: 6]
void write_to_file(std::vector<double>& a, std::ofstream& f)
{
    f.write(reinterpret_cast<char*>(&a[0]),sizeof(double)*a.size());
}

// Reads binary data into a 1D vector and accumulates it[cite: 6]
void read_from_file(std::vector<double>& b, std::ifstream& f,int size)
{
    std::vector<double> read(size,0);
    f.read(reinterpret_cast<char*>(&read[0]),sizeof(double)*size);
    for(int i = 0;i<size;i++)
    {
        b[i]+=read[i];
    }
}

// Reads 4D spatial-temporal fields from binary files into the 5D storage vector[cite: 6]
void read_cells_element(std::string element,std::vector<std::vector<std::vector<std::vector<double>>>>& cells, int size)
{
    std::ifstream file(element,std::ios::binary);
    if (!file) {
        std::cerr << "cannot open " << element << std::endl;
        std::exit(1);
    }
    for(int t=0;t<size;t++)
    {
        for(int x=0;x<size;x++)
        {
            for(int y=0;y<size;y++)
            {
                read_from_file(cells[t][x][y],file,size);
            }
        }
    }
    file.close();
}

// Cells without matter have no momentum anisotropy (NaN in the extracted files). They are replaced
// by `fill` (1 = isotropic) so that the linear interpolation near the surface is not poisoned.
long long fill_nan(std::vector<std::vector<std::vector<std::vector<double>>>>& cells, double fill)
{
    long long count = 0;
    for (auto& a : cells)
        for (auto& b : a)
            for (auto& c : b)
                for (auto& v : c)
                    if (!std::isfinite(v)) { v = fill; count++; }
    return count;
}

// --- MAIN EXECUTION ---

int main(int argc, char** argv)
{
    // Require the parameter file as a command-line argument
    if (argc < 2) {
        std::cout << "Usage: ./executable config.txt\n";
        return 1;
    }

    // Load parameters from the file
    auto params = load_parameters(argv[1]);
    double dx = std::stod(params["dx"]);
    double dt = std::stod(params["dt"]);
    double xmax = std::stod(params["xmax"]);
    double n0 = std::stod(params["n0"]);
    int tstart = std::stoi(params["tstart"]);
    int tend = std::stoi(params["tend"]);
    int nx = std::stoi(params["nx"]);
    int ny = std::stoi(params["ny"]);
    int nz = std::stoi(params["nz"]);
    std::string path = params["path"];

    // Optional local momentum anisotropy: files exx.dat, eyy.dat, ezz.dat (= 3Txx/(Txx+Tyy+Tzz) etc.)
    // on the same grid as e.dat, n.dat, ...; they are interpolated and appended to the surface output.
    // config: "anisotropy 1" to enable, "aniso_fill 1.0" = value used for cells without matter (NaN)
    bool use_aniso = params.count("anisotropy") && std::stoi(params["anisotropy"]) != 0;
    double aniso_fill = params.count("aniso_fill") ? std::stod(params["aniso_fill"]) : 1.0;
    const int nq = use_aniso ? 10 : 7;

    // Grid calculations
    double dV = dx*dx*dx;
    int nbins=xmax/dx;
    double d4x=dV*dx;

    // Define indices for the 7 fluid quantities[cite: 6]
    enum columns{ e,T,n,gamma,vx,vy,vz,exx,eyy,ezz };
    
    // Allocate a 5D vector to store [quantity][time][x][y][z] initialized to zero[cite: 6]
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> Tnv(nq,std::vector<std::vector<std::vector<std::vector<double>>>>(nbins,std::vector<std::vector<std::vector<double>>>(nbins,std::vector<std::vector<double>>(nbins,std::vector<double>(nbins,0)))));
    
    std::cout<<"Begin reading temperatures densities etc"<<std::endl;
    // Load physical properties from the specified path[cite: 6]
    read_cells_element(path+"/e.dat",Tnv[e],nbins);
    read_cells_element(path+"/n.dat",Tnv[n],nbins);
    std::cout<<"read n"<<std::endl;
    read_cells_element(path+"/gamma.dat",Tnv[gamma],nbins);
    std::cout<<"read gamma"<<std::endl;
    read_cells_element(path+"/vx.dat",Tnv[vx],nbins);
    read_cells_element(path+"/vy.dat",Tnv[vy],nbins);
    read_cells_element(path+"/vz.dat",Tnv[vz],nbins);
    std::cout<<"read v"<<std::endl;
    if (use_aniso)
    {
        read_cells_element(path+"/exx.dat",Tnv[exx],nbins);
        read_cells_element(path+"/eyy.dat",Tnv[eyy],nbins);
        read_cells_element(path+"/ezz.dat",Tnv[ezz],nbins);
        long long nfill = fill_nan(Tnv[exx],aniso_fill) + fill_nan(Tnv[eyy],aniso_fill) + fill_nan(Tnv[ezz],aniso_fill);
        std::cout<<"read local momentum anisotropy exx, eyy, ezz (" << nfill << " empty-cell values set to " << aniso_fill << ")"<<std::endl;
    }
    std::cout<<"Finished reading temperatures densities etc"<<std::endl;

    int nt = tend; // Using 'tend' as the total time steps to process
    double Nfreeze = n0;

    // Initialize the surface builder with the loaded parameters[cite: 6]
    HyperSurfaceBuilder builder(nt, nx, ny, nz, Nfreeze, dt, dx, nq);
    
    // Iterate over the 4D spacetime grid to find the freeze-out surface[cite: 6]
    for(int i=tstart;i<tend-1;i++)
    {
        std::cout<< "Processing timestep: " << i << std::endl;
        for(int j=0;j<nx-1;j++)
        {
            for(int k=0;k<ny-1;k++)
            {
                for(int l=0;l<nz-1;l++)
                {
                    int test_prev=0;
                    int test_curr=0;
                    
                    // Check a 2x2x2 cube of cells to see if the density crosses the freeze-out threshold n0[cite: 6]
                    for (int jj=0;jj<2;jj++)
                    {
                        for (int kk=0;kk<2;kk++)
                        {
                            for (int ll=0;ll<2;ll++)
                            {
                                test_prev += Tnv[n][i][j+jj][k+kk][l+ll]>n0;
                                test_curr += Tnv[n][i+1][j+jj][k+kk][l+ll]<n0;
                            }
                        }
                    }
                    // If the surface crosses this local hypercube, extract the hypersurface elements[cite: 6]
                    // (Note: The if-statement condition was commented out in the original source, uncomment if needed)
                    //if(test_curr && test_prev)
                    {
                        builder.find_elements(i,j,k,l,Tnv);
                    }
                }
            }
        }
    }
    return 0;
}
