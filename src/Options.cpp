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

#include "Options.h"

IsParallelOptions get_options(int argc, char* argv[])
{
    static struct option long_options[] = {
        {"vdw-scalar", required_argument, 0, 'v'},
        {"vdw-dist",   required_argument, 0, 'V'},
        {"cov-scalar", required_argument, 0, 'c'},
        {"cov-dist",   required_argument, 0, 'C'},
        {"maxhkl",     required_argument, 0, 'm'},
        {"filename",   required_argument, 0, 'f'},
        {"output",     required_argument, 0, 'o'},
        {"details",    no_argument,       0, 'p'},
        {"expand",     no_argument,       0, 'x'},
        {"molecule",   no_argument,       0, 'M'},
        {0, 0, 0, 0}};

    string filename;

    IsParallelOptions ans {false, false, false, false, 0.80, 1.00, 1.20, 1.00, 3, "input.txt", "", false, false, false};

    while( true ){
        int option_index = 0;
        char c = getopt_long(argc, argv, "v:V:c:C:m:f:o:pxM", long_options, &option_index);
        if( c == -1 ) break;
        switch(c){
            case 'v':
                ans.flag_vdw_scalar = true;
                ans.vdw_scalar_min = atof(optarg);
                break;
            case 'V':
                ans.flag_vdw_dist = true;
                ans.vdw_dist_minus = atof(optarg);
                break;
            case 'c':
                ans.flag_cov_scalar = true;
                ans.cov_scalar_max = atof(optarg);
                break;
            case 'C':
                ans.flag_cov_dist = true;
                ans.cov_dist_plus = atof(optarg);
                break;
            case 'm':
                ans.maxhkl = atoi(optarg);
                break;
            case 'f':
                ans.filename = optarg;
                break;
            case 'o':
                ans.output = optarg;
                break;
            case 'p':
                ans.flag_details = true;
                break;
            case 'x':
                ans.flag_expand = true;
                break;
            case 'M':
                ans.flag_mol = true;
                break;
            case '?':
                cout << "unknow argument: " << optopt << endl;
                break;
            default:
                cout << "wrong arguments!" << endl;
                exit(-1);
        }
    }

    if( optind < argc ) ans.filename = argv[optind];

    if( ! (ans.flag_vdw_scalar || ans.flag_vdw_dist || ans.flag_cov_scalar || ans.flag_cov_dist) ) {
        ans.flag_vdw_dist = true;
        ans.vdw_dist_minus = 1.0;
    }
    if( ans.output.length() == 0 )
        if( ans.filename.rfind(".txt") == string::npos )
            ans.output = ans.filename+".out";
        else
            ans.output = ans.filename.substr(0, ans.filename.rfind(".txt")) + ".out";

    return ans;
}

