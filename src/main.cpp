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

// Reads one time slice [x][y][z] (nbins^3 doubles) of a 4D field file into cells[x][y][z]
void read_slice(std::ifstream& file, int t, std::vector<std::vector<std::vector<double>>>& cells, int size)
{
    file.clear();
    file.seekg((std::streamoff)t*size*size*size*sizeof(double), std::ios::beg);
    for(int x=0;x<size;x++)
        for(int y=0;y<size;y++)
        {
            file.read(reinterpret_cast<char*>(cells[x][y].data()),sizeof(double)*size);
            if(!file) { std::cerr << "read error at t=" << t << std::endl; std::exit(1); }
        }
}

// Cells without matter have no momentum anisotropy (NaN in the extracted files). They are replaced
// by `fill` (1 = isotropic) so that the linear interpolation near the surface is not poisoned.
long long fill_nan(std::vector<std::vector<std::vector<double>>>& cells, double fill)
{
    long long count = 0;
    for (auto& b : cells)
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
    int aniso_mode = params.count("anisotropy") ? std::stoi(params["anisotropy"]) : 0;   // 0 off, 1 diagonal (exx,eyy,ezz), 2 full tensor (+exy,exz,eyz)
    bool use_aniso = aniso_mode != 0;
    double aniso_fill = params.count("aniso_fill") ? std::stod(params["aniso_fill"]) : 1.0;
    const int nq = aniso_mode >= 2 ? 13 : (aniso_mode == 1 ? 10 : 7);

    // Grid calculations
    int nbins=xmax/dx;

    // Define indices for the 7 fluid quantities[cite: 6]
    enum columns{ e,T,n,gamma,vx,vy,vz,exx,eyy,ezz,exy,exz,eyz };
    
    // Rolling window: only two time slices [slot][x][y][z] per quantity are kept in memory
    typedef std::vector<std::vector<std::vector<double>>> Cube3;
    std::vector<std::vector<Cube3>> Tnv(nq,std::vector<Cube3>(2,Cube3(nbins,std::vector<std::vector<double>>(nbins,std::vector<double>(nbins,0)))));
    const char* fname[13]={"e.dat",nullptr,"n.dat","gamma.dat","vx.dat","vy.dat","vz.dat","exx.dat","eyy.dat","ezz.dat","exy.dat","exz.dat","eyz.dat"};
    std::vector<std::ifstream> files(nq);
    for(int q=0;q<nq;q++)
    {
        if(!fname[q]) continue;
        files[q].open(path+"/"+fname[q],std::ios::binary);
        if(!files[q]) { std::cerr << "cannot open " << path+"/"+fname[q] << std::endl; return 1; }
    }
    auto load = [&](int slot,int t){
        for(int q=0;q<nq;q++){
            if(!fname[q]) continue;
            read_slice(files[q],t,Tnv[q][slot],nbins);
            if(q>=exx) fill_nan(Tnv[q][slot], q>=exy ? 0.0 : aniso_fill);   // empty cells: isotropic (e_ii = 1, e_ij = 0)
        }
    };
    std::string outfile = params.count("outfile") ? params["outfile"] : std::string("surface.dat");
    double Nfreeze = n0;
    int nt = tend;
    HyperSurfaceBuilder builder(nt, nx, ny, nz, Nfreeze, dt, dx, nq, outfile);

    load(0,tstart); load(1,tstart+1);
    for(int i=tstart;i<tend-1;i++)
    {
        if(i>tstart){ for(int q=0;q<nq;q++) std::swap(Tnv[q][0],Tnv[q][1]); load(1,i+1); }
        std::cout<< "Processing timestep: " << i << std::endl;
        for(int j=0;j<nx-1;j++)
            for(int k=0;k<ny-1;k++)
                for(int l=0;l<nz-1;l++)
                {
                    // skip hypercubes in which n is on the same side of n0 at all 16 corners (no surface element possible)
                    int above=0;
                    for(int tt=0;tt<2;tt++) for(int jj=0;jj<2;jj++) for(int kk=0;kk<2;kk++) for(int ll=0;ll<2;ll++)
                        above += Tnv[n][tt][j+jj][k+kk][l+ll] > n0;
                    if(above==0 || above==16) continue;
                    builder.find_elements(i,0,j,k,l,Tnv);
                }
    }
    return 0;
}
