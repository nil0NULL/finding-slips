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
/*  Description: Defination of IsParallelOptions type and the                 */
/*  get_options function, which is used for providing the arguments main      */
/*  program isparallel.cpp.                                                   */
/*                                                                            */
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <string>
#include <getopt.h>

#ifndef OPTIONS_H
#define OPTIONS_H

using namespace std;

typedef struct {
    bool flag_vdw_scalar;
    bool flag_vdw_dist;
    bool flag_cov_scalar;
    bool flag_cov_dist;
    double vdw_scalar_min;
    double vdw_dist_minus;
    double cov_scalar_max;
    double cov_dist_plus;
    int maxhkl;
    string filename;
    string output;
    bool flag_details;
    bool flag_expand;
    bool flag_mol;
} IsParallelOptions;

IsParallelOptions get_options(int argc, char* argv[]);

#endif
