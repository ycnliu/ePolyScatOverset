/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  Modified by Loren Greenman from LIBINT hartree-fock example
 *      Tools and objects to assist interface to LIBINT
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

#ifndef _integration_tools_cpp_header_
#define _integration_tools_cpp_header_

// standard C++ headers
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <chrono>

// Eigen matrix algebra library
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// Libint Gaussian integrals library
#include <libint2.hpp>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        Matrix;  // import dense, dynamically sized Matrix type from Eigen;
                 // this is a matrix with row-major storage (http://en.wikipedia.org/wiki/Row-major_order)
                 // to meet the layout of the integrals returned by the Libint integral library

struct Atom {
    int atomic_number;
    double x, y, z;
};

using libint2::Shell;

// These geometry/basis reading routines are not really necessary to the interface
std::vector<Atom> read_geometry(const std::string& filename);
std::vector<Shell> make_sto3g_basis(const std::vector<Atom>& atoms);

// These routines count basis size and map Shell objects to Basis objects (slightly decorated Shells)
size_t nbasis(const std::vector<Shell>& shells);
std::vector<size_t> map_shell_to_basis_function(const std::vector<Shell>& shells);
size_t max_nprim(const std::vector<Shell>& shells);
int max_l(const std::vector<Shell>& shells);

// This routine will compute a set of 1-body integrals between shells
//  These integrals can be of any 1-body type, including 
//  kinetic, nuclear, etc.
// I create a similar computation function just for nuclear integrals that 
// doesn't re-initialize the 1-body engine every time in ElectrostaticPointIntegrals.cc
Matrix compute_1body_ints(const std::vector<libint2::Shell>& shells,
                          libint2::OneBodyEngine::operator_type t,
                          const std::vector<Atom>& atoms = std::vector<Atom>());
//
Matrix compute_2body_fock(const std::vector<Shell>& shells,
			  const Matrix& D);
Matrix compute_2body_fock_simple(const std::vector<Shell>& shells,
			  const Matrix& D);

#endif

