//------------------------------------------------------------------------------
//  The @fixed_width, @begin_table, @item2 and @end_table commands are custom
//  defined commands in Doxygen.in. They are defined under ALIASES. For the page
//  created here, the 80 column limit is exceeded. Arguments of aliases are
//  separated by ','. If you intended ',' to be a string you must use an escaped
//  comma '\,'.
//
///  @page lgrid_sec Limiter Gird File
///
///  @tableofcontents
///  @section lgrid_intro_sec Introduction
///  LGRID is a code to precompute the distance from a grid point to the closest
///  point on a limiter. Limiters are specified as a series of line segments in
///  a counter clockwise fasion. That is a polygon with vertices ordered in a
///  counter clockwise will have the inside of the polygon represent the
///  interior. Clockwise ordering will have the outside of the polygon represent
///  the interior.
///
///  @section lgrid_input_sec Input Configuation
///  All limiters are defined using the same format. A key word, the number of
///  phi planes the limiter surface is associated with, a list of phi values and
///  then a list of r z segment positions terminated by a
///  @fixed_width{end_limiter} keyword.
///
///  @fixed_width{keyword\n}
///  @fixed_width{num_phi\n}
///  @fixed_width{phi1(REAL) phi2(REAL) ... phiN(REAL)\n}
///  @fixed_width{r1(REAL) z1(REAL)\n}
///  @fixed_width{r2(REAL) z2(REAL)\n}
///  @fixed_width{...\n}
///  @fixed_width{rM(REAL) zM(REAL)\n}
///  @fixed_width{end_limiter\n}
///
///  @subsection lgrid_input_keyword_sec Keywords
///  @begin_table
///  @item2{@fixed_width{new_limiter_deg}, A limiter at phi angles specified in degrees.}
///  @item2{@fixed_width{new_limiter_rad}, A limiter at phi angles specified in radians.}
///  @item2{@fixed_width{end_limiter},     End of the limiter.}
///  @end_table
///
///  @subsection lgrid_input_exam_sec Example File
///  @code
///  new_limiter_deg
///  3
///  10.0 50.0 175.0
///  0.5 -0.25
///  1.0 -0.25
///  1.0  0.25
///  0.5  0.25
///  end_limiter
///  @endcode
///
///  @section lgrid_command_line_sec Command Line Arguments
///  LGRID is run with the following command.
///
///  @fixed_width{lgrid filename r0 dr numr z0 dz numz}
///
///  @begin_table
///  @item2{@fixed_width{filename}, Path to file containing the limiter configuation.}
///  @item2{@fixed_width{r0},       Starting radial grid position.}
///  @item2{@fixed_width{dr},       Radial grid width.}
///  @item2{@fixed_width{numr},     Number of radial grid points.}
///  @item2{@fixed_width{z0},       Starting vertical grid position.}
///  @item2{@fixed_width{dz},       Vertical grid width.}
///  @item2{@fixed_width{numz},     Number of vertical grid points.}
///  @end_table
//------------------------------------------------------------------------------
//******************************************************************************
///  @file main.cpp
///  @brief Contains the main rotutines for LGRID.
//
///  @author Mark Cianciosa on 10/12/15.
///
///  LGRID is code for generating responce files for any arbitary shaped limiter
///  structure. This creates a grid a distances to the nearest plasma facing
///  component in a specified phi plane.
//******************************************************************************

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <netcdf.h>
#include <vector>

#include "vertex.hpp"

///  Limiter plane keyword specified in degrees.
const std::string new_limiter_deg = "new_limiter_deg";
///  Limiter plane keyword specified in radians.
const std::string new_limiter_rad = "new_limiter_rad";

///  Keyword token types.
///  * degrees Limiter defined in degrees.
///  * radians Limiter defined in radians.
///  * eof     End of file reached.
enum limiter_type {
    degrees,
    radians,
    eof
};

//------------------------------------------------------------------------------
///  @brief Finds the next keyword.
///
///  Reads the input file one file at a time searching for the next keyword
///  until the end of the file is reached.
///
///  @param[inout] file Input file stream to seach for next keyword.
///  @return An @ref limiter_type token.
//------------------------------------------------------------------------------
const limiter_type get_next_keyword(std::ifstream &file) {
    std::string line;
    
    for (file >> line; !file.eof() && !file.fail() && !file.bad(); file >> line) {
        if (line.compare(new_limiter_deg) == 0) {
            return degrees;
        } else if (line.compare(new_limiter_rad) == 0) {
            return radians;
        }
    }
    return eof;
}

//------------------------------------------------------------------------------
///  @brief Parses a limiter specification.
///
///  Limiters are defined by the number of phi planes this plasma facing
///  component is asscoiated with. This is followed by a list of each phi. A
///  collection of R, Z values defing the crossection is read until two doubles
///  could not be read or the file ends.
///
///  @param[inout] file     Input file stream read limiter from.
///  @param[inout] limiters Limiter vertex data.
///  @param[in]    type     Limiter type.
//------------------------------------------------------------------------------
void parse_limiter(std::ifstream &file, std::map<const double, std::vector<vertex *> > &limiters, const limiter_type type) {
    size_t num_phi;
    file >> num_phi;
    
    std::vector<double> phi_angles;
    for (size_t i = 0; i < num_phi; i++) {
        double phi;
        file >> phi;
        
        if (type == degrees) {
            phi *= M_PI/180.0;
        }
        
        phi_angles.push_back(phi);
    }
    
    double r;
    double z;
    
    vertex *limiter = nullptr;
    
    for (file >> r >> z; !file.eof() && !file.fail() && !file.bad(); file >> r >> z) {
        if (limiter == nullptr) {
            limiter = new vertex(vector3d(r, z, 0.0));
        } else {
            limiter->insert(new vertex(vector3d(r, z, 0.0)));
        }
    }
    
//  Clear the bad and fail bits so the parsing can continue. This bits get set
//  when the file fails to read two doubles. Clear doesn't actually clear error
//  bits instead see if the eof bit is set then set only that bit.
    file.eof() ? file.clear(std::iostream::eofbit) : file.clear();
    
    for (double phi : phi_angles) {
        limiters[phi].push_back(limiter);
    }
}

//------------------------------------------------------------------------------
///  @brief Parses a limiter specification.
///
///  Opens a file, loops through all the keywords and parses the limiter file.
///
///  @param[in] file_name Limiter type.
///  @return An @ref limiter_type token.
//------------------------------------------------------------------------------
std::map<const double, std::vector<vertex *> > parse_limiter_file(const std::string &file_name) {

    std::ifstream file(file_name);
    
    std::map<const double, std::vector<vertex *> > limiters;
    
    for (limiter_type type = get_next_keyword(file); type != eof; type = get_next_keyword(file)) {
        parse_limiter(file, limiters, type);
    }
    
    file.close();
    
    return limiters;
}

//------------------------------------------------------------------------------
///  @brief Main program.
///
///  Parses the limter file. Computes distances to the nearest limiter element.
///  Creates the lgrid netcdf file. This program expects 7 arguments in the
///  following order.
///
///  * Filename to parse define the limiter.
///  * Starting radial position to define the radial grid.
///  * Radial grid spacing.
///  * Number of radial grid points.
///  * Starting vertical position to define the vertical grid.
///  * Vertical grid spacing.
///  * Number of vertical grid points.
///
///  @param[in] argc Number of command line arguments.
///  @param[in] argv Arrays of argument strings.
///  @return 0 if completed.
//------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
    
    if (argc != 8) {
        std::cout << "Expected 7 arguments." << std::endl;
        std::cout << "lgrid filename r0 dr numr z0 dz numz" << std::endl;
    }
    
    std::map<const double, std::vector<vertex *> > limiters = parse_limiter_file(argv[1]);
    
    const double r0 = atof(argv[2]);
    const double dr = atof(argv[3]);
    const size_t num_r = atoll(argv[4]);

    const double z0 = atof(argv[5]);
    const double dz = atof(argv[6]);
    const size_t num_z = atoll(argv[7]);
    
    std::vector<double> grids;
    
    int ncid;
    nc_create((std::string(argv[1]) + std::string(".nc")).c_str(), NC_WRITE, &ncid);
    
    int num_phi_dimid;
    nc_def_dim(ncid, "num_phi", limiters.size(), &num_phi_dimid);
    int num_r_dimid;
    nc_def_dim(ncid, "num_r", num_r, &num_r_dimid);
    int num_z_dimid;
    nc_def_dim(ncid, "num_z", num_z, &num_z_dimid);
    
    std::vector<double> phi_angles;
    
    for (std::pair<const double, std::vector<vertex *> > &phi : limiters) {
        phi_angles.push_back(phi.first);
        for (size_t j = 0; j < num_z; j++) {
            const double z = z0 + j*dz;
            for (size_t i = 0; i < num_r; i++) {
                const double r = r0 + i*dr;
                std::vector<double> d;
                for (vertex *limiter : phi.second) {
                    vertex::distance(limiter, vector3d(r, z, 0.0), d);
                }
                std::sort(d.begin(), d.end(), [](double a, double b) {
                    return fabs(a) < fabs(b);
                });
                grids.push_back(d[0]);
            }
        }
    }
    
    int phi_angles_varid;
    nc_def_var(ncid, "phi_angles", NC_DOUBLE, 1, &num_phi_dimid, &phi_angles_varid);
    int r0_varid;
    nc_def_var(ncid, "r0", NC_DOUBLE, 0, nullptr, &r0_varid);
    int dr_varid;
    nc_def_var(ncid, "dr", NC_DOUBLE, 0, nullptr, &dr_varid);
    int z0_varid;
    nc_def_var(ncid, "z0", NC_DOUBLE, 0, nullptr, &z0_varid);
    int dz_varid;
    nc_def_var(ncid, "dz", NC_DOUBLE, 0, nullptr, &dz_varid);
    
    int grids_varid;
    std::array<int, 3> grid_dims;
    grid_dims[0] = num_phi_dimid;
    grid_dims[1] = num_z_dimid;
    grid_dims[2] = num_r_dimid;
    nc_def_var(ncid, "iso_grids", NC_DOUBLE, 3, grid_dims.data(), &grids_varid);
    
    nc_enddef(ncid);
    
    nc_put_var(ncid, phi_angles_varid, phi_angles.data());
    nc_put_var(ncid, r0_varid, &r0);
    nc_put_var(ncid, dr_varid, &dr);
    nc_put_var(ncid, z0_varid, &z0);
    nc_put_var(ncid, dz_varid, &dz);
    nc_put_var(ncid, grids_varid, grids.data());
    
    nc_close(ncid);
    
    return 0;
}
