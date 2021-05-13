/*************************   Copyright and Licenses   *************************/
/*                                                                            */
/*  Copyright(C) 2021 Xudong He                                               */
/*                                                                            */
/*  This file is part of isparallel.                                          */
/*                                                                            */
/*  This program is free software: you can redistribute it and/or modify      */
/*  it under the terms of the GNU General Public License as published by      */
/*  the Free Software Foundation, either version 3 of the License, or         */
/*  (at your option) any later version.                                       */
/*                                                                            */
/*  Isparallel is distributed in the hope that it will be useful,             */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/*  General Public License for more details.                                  */
/*                                                                            */
/*  You should have received a copy of the GNU General Public License         */
/*  along with this program.  If not, see <https://www.gnu.org/licenses/>.    */
/*                                                                            */
/******************************************************************************/

/******************************************************************************/
/*                                                                            */
/*  Author: hxd_yi (Xudong He)                                                */
/*  E-mail: hxd_yi@163.com                                                    */
/*  Description: This file define the parts function and other related        */
/*  functions, including distance calculation functions.                      */
/*                                                                            */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <set>
#include <algorithm>
#include <functional>
#include <numeric>
#include "Point.h"
#include "Atom.h"
#include "Options.h"

#ifndef PARTS_H
#define PARTS_H

using namespace std;

//const double ABSZERO = 1e-6;

typedef struct {
    Point a, b, c;
    size_t natoms;
    vector<Atom> atoms;
} Cell;

typedef struct {
    int h, k, l;
    Point n;
} xxx; 

typedef struct {
    double value;
    size_t atom_i, atom_j;
} Part_Dist;

typedef struct {
    size_t count;
    Part_Dist vdw_dist;
    Part_Dist vdw_scalar;
    Part_Dist cov_dist;
    Part_Dist cov_scalar;
} Part;

int gcd(int m, int n);
Part parts(const Cell& cif, function<double(const Point&)> fun, const IsParallelOptions& options);
double calc_dist(const Point& p, function<double(const Point& xyz)> fun, const vector<xxx>& cell_seq);
double dist_vector(const Point& vec, const Point& n);
double dist_plane(const Point& vec, const Point& n1, const Point& n2);
bool isclose(const Atom& a, const Atom& b, double dist, const IsParallelOptions& options);
vector<xxx>& construct_cell_seq(const Cell& cif, vector<xxx>& cell_seq);
vector<xxx>& construct_cell_seq(const Cell& cif, const xxx& n1, vector<xxx>& cell_seq);
vector<xxx>& construct_cell_seq(const Cell& cif, const xxx& n1, const xxx& n2, vector<xxx>& cell_seq);
vector<xxx>& clear_cell_seq(vector<xxx>& cell_seq);

#endif

