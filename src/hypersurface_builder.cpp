#include "hypersurface_builder.hpp"

// Performs sequential linear interpolation across 4 dimensions[cite: 8]
void HyperSurfaceBuilder::fD_Interpolate(std::vector<double> shift,double to_interpolate[7][2][2][2][2], std::vector<double>& interpolated)
{
    // Loop through all 7 physical quantities (energy, density, velocity, etc.)[cite: 8]
    for(int ii=0;ii<7;ii++)
    {
        // 1. Interpolate across time (DT) to collapse from 4D down to 3D[cite: 8]
        double interpolated_3D[2][2][2];
        for(int jj=0;jj<2;jj++)
        {
            for(int kk=0;kk<2;kk++)
            {
                for(int ll=0;ll<2;ll++)
                {
                    interpolated_3D[jj][kk][ll]=to_interpolate[ii][0][jj][kk][ll]+(to_interpolate[ii][1][jj][kk][ll]-to_interpolate[ii][0][jj][kk][ll])*shift[0]/DT;
                }
            }
        }
        
        // 2. Interpolate across the X dimension (DX) to collapse from 3D down to 2D[cite: 8]
        double interpolated_2D[2][2];
        for(int jj=0;jj<2;jj++)
        {
            for(int kk=0;kk<2;kk++)
            {
                    interpolated_2D[jj][kk]=interpolated_3D[0][jj][kk]+(interpolated_3D[1][jj][kk]-interpolated_3D[0][jj][kk])*shift[1]/DX;
            }
        }
        
        // 3. Interpolate across the Y dimension (DY) to collapse from 2D down to 1D[cite: 8]
        double interpolated_1D[2];
        for(int jj=0;jj<2;jj++)
        {
                    interpolated_1D[jj]=interpolated_2D[0][jj]+(interpolated_2D[1][jj]-interpolated_2D[0][jj])*shift[2]/DY;
        }
        
        // 4. Interpolate across the Z dimension (DZ) to get the final point value[cite: 8]
        interpolated[ii]=interpolated_1D[0]+(interpolated_1D[1]-interpolated_1D[0])*shift[3]/DZ;
    }
}

// Constructor: Initializes the grid, output file, and the Cornelius algorithm[cite: 8]
HyperSurfaceBuilder::HyperSurfaceBuilder(int nt, int nx, int ny, int nz, double freezeout_N, double dt, double dx)
    : nt_(nt), nx_(nx), ny_(ny), nz_(nz), Nfreeze_(freezeout_N)
    {
    this->s_file.open("surface.dat"); // File where hypersurface elements are written[cite: 8]
    this->dim = 4;
    this->DT = dt;
    this->DX = dx;
    this->DY = dx;
    this->DZ = dx;

    // Pass the 4D lattice spacing to the Cornelius algorithm[cite: 8]
    double lattice_spacing[4] = {DT, DX, DY, DZ};
    cornelius_ptr=new Cornelius();
    cornelius_ptr->init(dim, freezeout_N, lattice_spacing);
}

// Extracts hypersurface elements from a localized 2x2x2x2 grid[cite: 8]
void HyperSurfaceBuilder::find_elements(int i,int j, int k, int l, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& Tnv)
{
    double lattice_spacing[4] = {DT, DX, DY, DZ};

    // Dynamically allocate a 4D array to hold the density variable for the Cornelius algorithm[cite: 8]
    double ****cube;
    cube = new double***[2];
    for (int ii=0; ii<2; ++ii)
    {
        cube[ii] = new double**[2];
        for (int jj=0; jj<2; ++jj)
        {
            cube[ii][jj] = new double*[2];
            for (int kk=0; kk<2; ++kk)
            {
                cube[ii][jj][kk] = new double[2];
            }
        }
    }

    // Array to store all 7 physical quantities across the localized hypercube[cite: 8]
    double hydro_cube[7][2][2][2][2];

    // Populate the 4D cube with density values (used by Cornelius to locate the surface)[cite: 8]
    cube[0][0][0][0] = Tnv[n][i][j][k][l];
    cube[0][0][1][0] = Tnv[n][i][j][k+1][l];
    cube[0][1][0][0] = Tnv[n][i][j+1][k][l];
    cube[0][1][1][0] = Tnv[n][i][j+1][k+1][l];
    cube[1][0][0][0] = Tnv[n][i+1][j][k][l];
    cube[1][0][1][0] = Tnv[n][i+1][j][k+1][l];
    cube[1][1][0][0] = Tnv[n][i+1][j+1][k][l];
    cube[1][1][1][0] = Tnv[n][i+1][j+1][k+1][l];
    cube[0][0][0][1] = Tnv[n][i][j][k][l+1];
    cube[0][0][1][1] = Tnv[n][i][j][k+1][l+1];
    cube[0][1][0][1] = Tnv[n][i][j+1][k][l+1];
    cube[0][1][1][1] = Tnv[n][i][j+1][k+1][l+1];
    cube[1][0][0][1] = Tnv[n][i+1][j][k][l+1];
    cube[1][0][1][1] = Tnv[n][i+1][j][k+1][l+1];
    cube[1][1][0][1] = Tnv[n][i+1][j+1][k][l+1];
    cube[1][1][1][1] = Tnv[n][i+1][j+1][k+1][l+1];

    // Populate the full hydro_cube with all 7 variables for interpolation[cite: 8]
    for(int ii=0;ii<7;ii++)
    {
        hydro_cube[ii][0][0][0][0] = Tnv[ii][i][j][k][l];
        hydro_cube[ii][0][0][1][0] = Tnv[ii][i][j][k+1][l];
        hydro_cube[ii][0][1][0][0] = Tnv[ii][i][j+1][k][l];
        hydro_cube[ii][0][1][1][0] = Tnv[ii][i][j+1][k+1][l];
        hydro_cube[ii][1][0][0][0] = Tnv[ii][i+1][j][k][l];
        hydro_cube[ii][1][0][1][0] = Tnv[ii][i+1][j][k+1][l];
        hydro_cube[ii][1][1][0][0] = Tnv[ii][i+1][j+1][k][l];
        hydro_cube[ii][1][1][1][0] = Tnv[ii][i+1][j+1][k+1][l];
        hydro_cube[ii][0][0][0][1] = Tnv[ii][i][j][k][l+1];
        hydro_cube[ii][0][0][1][1] = Tnv[ii][i][j][k+1][l+1];
        hydro_cube[ii][0][1][0][1] = Tnv[ii][i][j+1][k][l+1];
        hydro_cube[ii][0][1][1][1] = Tnv[ii][i][j+1][k+1][l+1];
        hydro_cube[ii][1][0][0][1] = Tnv[ii][i+1][j][k][l+1];
        hydro_cube[ii][1][0][1][1] = Tnv[ii][i+1][j][k+1][l+1];
        hydro_cube[ii][1][1][0][1] = Tnv[ii][i+1][j+1][k][l+1];
        hydro_cube[ii][1][1][1][1] = Tnv[ii][i+1][j+1][k+1][l+1];
    }
    
    std::vector<std::vector<double>> surface_elements;
    
    // Execute the Cornelius algorithm to extract normal vectors and centroids[cite: 8]
    cornelius_ptr->find_surface_4d(cube);

    // Retrieve properties of the calculated freeze-out surface elements[cite: 8]
    for (int isurf = 0; isurf < cornelius_ptr->get_Nelements(); isurf++)
    {
        // 1. Get the surface normal vector d^3 \sigma_\mu[cite: 8]
        std::vector<double> surface_element;
        for (int ii = 0; ii < 4; ii++)
        {
            surface_element.push_back(cornelius_ptr->get_normal_elem(isurf, ii));
        }

        // 2. Calculate the centroid shifts within the current lattice[cite: 8]
        std::vector<double> x_fraction(4,0);
        for (int ii = 0; ii < 4; ii++) {
            double t_shift = cornelius_ptr->get_centroid_elem(isurf, ii);
            x_fraction[ii] =  t_shift;
        }
        
        // 3. Compute absolute physical coordinates[cite: 8]
        double t_center = i*DT + DT/2  +  x_fraction[0];
        surface_element.push_back(t_center);
        double x_center = j*DX - 25   +  x_fraction[1];
        surface_element.push_back(x_center);
        double y_center = k*DY - 25   +  x_fraction[2];
        surface_element.push_back(y_center);
        double z_center = l*DZ - 25   +  x_fraction[3];
        surface_element.push_back(z_center);

        // 4. Perform 4D linear interpolation for all fluid quantities at the centroid[cite: 8]
        std::vector<double> shift={x_fraction[0],x_fraction[1],x_fraction[2],x_fraction[3]};
        std::vector<double> interpolated(7,0);
        HyperSurfaceBuilder::fD_Interpolate(shift, hydro_cube, interpolated);
        
        // Append the interpolated physical quantities (T, e, v, etc.) to the element data[cite: 8]
        for(int ii=0;ii<7;ii++)
        {
            surface_element.push_back(interpolated[ii]);
        }
        surface_elements.push_back(surface_element);
    
        // 5. Write the extracted element to the output file[cite: 8]
        for(int ii=0;ii<surface_elements.size();ii++)
        {
            for(int jj=0;jj<surface_elements[ii].size();jj++)
            {
                s_file<<surface_elements[ii][jj] << "    ";
            }
            s_file<<std::endl;
        }
    }
    
    // Free the dynamically allocated 4D cube memory[cite: 8]
    for (int ia = 0; ia < 2; ++ia) {
        for (int ib = 0; ib < 2; ++ib) {
            for (int ic = 0; ic < 2; ++ic) {
                delete [] cube[ia][ib][ic];
            }
            delete [] cube[ia][ib];
        }
        delete [] cube[ia];
    }
    delete [] cube;
}
