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
/*  Description: This file is the main program of the isparallel,             */
/*  including input and output subroutines. It also include expand_cell       */
/*  function which makes supercell.                                           */
/*                                                                            */
/******************************************************************************/

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <set>
#include <algorithm>
#include <functional>
#include "Point.h"
#include "Atom.h"
#include "Options.h"
#include "Parts.h"

using namespace std;

//const double ABSZERO = 1e-6;
//const double vdw_minus = 1.0;
//const double scalar_min = 0.8;

vector<xxx> lines;
vector<xxx> planes;
vector<Part> lines_result;
vector<Part> planes_result;

Cell& expand_cell222(Cell& cif);
Cell& expand_cell(Cell& cif, char x, unsigned int n);
int processing(const Cell& cif, const IsParallelOptions& options);
ofstream& output_details(ofstream& fout, xxx& x, Part& y, char type, const Cell& cif, const IsParallelOptions& options);

int main(int argc, char* argv[])
{

//  input coordinations

    IsParallelOptions options = get_options(argc, argv);

    if( options.flag_details ) {
        cout << "filename:   " << options.filename << endl;
        cout << "maxhkl:     " << setw(6) << options.maxhkl << endl;
        if( options.flag_vdw_scalar ) cout << "vdw scalar: " << fixed << setw(6) << setprecision(2) << options.vdw_scalar_min << defaultfloat << endl;
        if( options.flag_vdw_dist )   cout << "vdw dist:   " << fixed << setw(6) << setprecision(2) << options.vdw_dist_minus << defaultfloat << endl;
        if( options.flag_cov_scalar ) cout << "cov scalar: " << fixed << setw(6) << setprecision(2) << options.cov_scalar_max << defaultfloat << endl;
        if( options.flag_cov_dist )   cout << "cov dist:   " << fixed << setw(6) << setprecision(2) << options.cov_dist_plus  << defaultfloat << endl;
    }

    ifstream fin(options.filename);
    //double x, y, z;
    Cell cif;

    fin >> cif.a >> cif.b >> cif.c;

    fin >> cif.natoms;
    cif.atoms.resize(cif.natoms);
    string ln;
    getline(fin, ln);
    for(size_t i = 0; i < cif.natoms; ++i){
        int n;
        string s;
        double x, y, z, m, cov, vdw;
        string line;
        getline(fin, line);
        istringstream lin(line);
        //cout << line << endl;
        lin >> n >> s >> x >> y >> z >> m >> cov >> vdw;
        if( lin.eof() ) {
            cif.atoms.at(i) = Atom(n, s, x, y, z, m, cov, vdw);
        }else{
            int id;
            lin >> id;
            cif.atoms.at(i) = Atom(n, s, x, y, z, m, cov, vdw, id);
            options.flag_mol = true;
        }
    }

    fin.close();
    if( options.flag_details ) cout << cif.natoms << " atoms had been read." << endl;
    
    //for_each(cif.atoms.begin(), cif.atoms.end(), [](const Atom& x) { cout << x.frac_xyz() << endl; });
    
// preprocess
    if( options.flag_expand ) expand_cell222(cif);
    vector<xxx> cell_seq;
    construct_cell_seq(cif, cell_seq);
    //for_each(cell_seq.begin(), cell_seq.end(), [](const xxx& x) { cout << x.h << ", " << x.k << ", " << x.l << " " << x.n << endl; });
    clear_cell_seq(cell_seq);
    auto fun = [&cell_seq](const Point& xyz) -> double { return calc_dist(xyz, [](const Point& x) -> double { return x.length(); }, cell_seq); };
    //cout << parts(cif, fun, options) << endl;
    Part result = parts(cif, fun, options);
    //cout << result.count << endl;
    if( result.count == 1 ) expand_cell222(cif);
    //for_each(cif.atoms.begin(), cif.atoms.end(), [](const Atom& x) { cout << x.mol_id() << ": " << x.atomic_number() << endl; });

// processing

    if( options.flag_details ) cout << "Beginning" << endl;
    int type = processing(cif, options);
    if( options.flag_details ) cout << "Finished" << endl;

// output

    if( options.flag_details ) cout << "Write detailed results to the output file." << endl;

    ofstream fout(options.output);
    //fout << "12345678901234567890123456789012345678901234567890123456789012345678901234567890" << endl;
    fout << options.filename << " has " 
	 << lines.size() << " shearing slip crystal orientation(s) and " 
	 << planes.size() << " shearing slip crystal plane(s)." << endl
         << endl;
    for(size_t i = 0; i < lines.size(); ++i)
	output_details(fout, lines[i], lines_result[i], 'o', cif, options);
    for(size_t i = 0; i < planes.size(); ++i)
	output_details(fout, planes[i], planes_result[i], 'p', cif, options);
    fout.close();

    //cout << "1234567890123456789012345678901234567890" << endl;
    if( type == 0 )
        cout << options.filename << " does\'t have any shearing slip orientation." << endl;
    else if( type > 0)
        cout << options.filename << " has " << type << " shearing slip crystal orientation(s)." << endl;
    else
        cout << options.filename << " has " << -type << " shearing slip crystal plane(s)." << endl;

    return 0;
}

int processing(const Cell& cif, const IsParallelOptions& options)
{
    for(int hh = 0; hh <= options.maxhkl; ++hh)
        for(int kk = 0; kk <= options.maxhkl; ++kk)
            for(int ll = 0; ll <= options.maxhkl; ++ll) {
                if( gcd(gcd(hh, kk), ll) != 1 )continue;

                int h, k, l;
                h = hh;
                k = kk;
                l = ll;

                do{
                    do{
                        //cout << h << k << l << endl;
                        Point n = h * cif.a + k * cif.b + l * cif.c;
			n.normalize();
                        if( lines.size() > 0 && count_if(lines.begin(), lines.end(),
                            [&n](const xxx& x) -> bool { return abs(cross(x.n, n).length() / x.n.length() / n.length()) < Point::ABSZERO; }) > 0)
                            continue;

                        if( options.flag_details )
                            cout << "checking orientation: < " << setw(5) << h << " " << setw(5) << k << " " << setw(5) << l << " >" << endl;
                        //if( planes.size() > 0 && count_if(planes.begin(), planes.end(),
                            //[&n](xxx x) -> bool { return abs(dot(x.n, n) / x.n.length() / n.length()) < Point::ABSZERO; }) > 0)
                            //continue;

                        vector<xxx> cell_seq;
                        xxx new_vec = {h, k, l, n};
                        construct_cell_seq(cif, new_vec, cell_seq);
                        //for_each(cell_seq.begin(), cell_seq.end(), [](const xxx& x) { cout << x.h << ", " << x.k << ", " << x.l << " " << x.n << endl; });
                        clear_cell_seq(cell_seq);
                        //for_each(cell_seq.begin(), cell_seq.end(), [](const xxx& x) { cout << x.h << ", " << x.k << ", " << x.l << " " << x.n << endl; });

                        auto fun = [&n, &cell_seq](const Point& xyz) -> double { 
                            return calc_dist(xyz, [&n](const Point& vec) -> double { return dist_vector(vec, n); }, cell_seq); };

                        Part result = parts(cif, fun, options);
                        if( result.count > 1 ){
                            if( options.flag_details ) 
                                cout << "crystal orientation : < " << setw(5) << h << " " << setw(5) << k << " " << setw(5) << l << " > found" << endl;
                            lines.insert(lines.end(), new_vec);
			    lines_result.insert(lines_result.end(), result);
                        }
                    }while( (l *= -1) < 0 );
                }while( (k *= -1) < 0 );
            }

    if( lines.size() > 1 ) {
        for(size_t i = 0; i < lines.size() - 1; ++i)
            for(size_t j = i + 1; j < lines.size(); ++j) {

                vector<xxx> cell_seq;
                cell_seq.clear();
                construct_cell_seq(cif, lines[i], lines[j], cell_seq);
                clear_cell_seq(cell_seq);

                xxx nn {lines[i].k * lines[j].l - lines[j].k * lines[i].l,
                        lines[i].l * lines[j].h - lines[j].l * lines[i].h,
                        lines[i].h * lines[j].k - lines[j].h * lines[i].k,
                        cross(lines[i].n, lines[j].n)};
                int d = gcd(gcd(nn.h, nn.k), nn.l);
                nn.h /= d;
                nn.k /= d;
                nn.l /= d;
                    nn.n.normalize();

                if( planes.size() > 0 && count_if(planes.begin(), planes.end(),
                    [&nn](const xxx& x) -> bool { return abs(cross(x.n, nn.n) / x.n.length() / nn.n.length()) < Point::ABSZERO; }) > 0)
                    continue;
                if( options.flag_details )
                    cout << "checking plane      : ( " << setw(5) << nn.h << " " << setw(5) << nn.k << " " << setw(5) << nn.l << " )" << endl;

                Point v1 = lines[i].n;
                Point v2 = lines[j].n;

                auto fun = [&v1, &v2, &cell_seq](const Point& xyz) -> double { 
                    return calc_dist(xyz, [&v1, &v2](const Point& vec) -> double { return dist_plane(vec, v1, v2); }, cell_seq); };

                Part result = parts(cif, fun, options);
                //cout << result.count << endl;
                if( result.count > 1 ){
                    if( options.flag_details )
                        cout << "crystal plane       : ( " << setw(5) << nn.h << " " << setw(5) << nn.k << " " << setw(5) << nn.l << " ) found" << endl;
                    planes.insert(planes.end(), nn);
                    planes_result.insert(planes_result.end(), result);
                }
            }
    }
    return planes.size() ? -planes.size() : ( lines.size() ? lines.size() : 0 );
}

Cell& expand_cell222(Cell& cif)
{
    expand_cell(cif, 'a', 2);
    expand_cell(cif, 'b', 2);
    expand_cell(cif, 'c', 2);
    return cif;
}

Cell& expand_cell(Cell& cif, char x, unsigned int n)
{
    if( n <= 1 ) return cif;

    Point d;

    switch(x){
        case 'a':
            d = Point(1.0 / n, 0.0, 0.0);
            for_each(cif.atoms.begin(), cif.atoms.end(), [&n](Atom& x) { Point p = x.frac_xyz(); p.x( p.x() / n ); x.frac_xyz(p); });
            cif.a *= n;
            break;
        case 'b':
            d = Point(0.0, 1.0 / n, 0.0);
            for_each(cif.atoms.begin(), cif.atoms.end(), [&n](Atom& x) { Point p = x.frac_xyz(); p.y( p.y() / n ); x.frac_xyz(p); });
            cif.b *= n;
            break;
        case 'c':
            d = Point(0.0, 0.0, 1.0 / n);
            for_each(cif.atoms.begin(), cif.atoms.end(), [&n](Atom& x) { Point p = x.frac_xyz(); p.z( p.z() / n ); x.frac_xyz(p); });
            cif.c *= n;
            break;
    }

    vector<Atom> addlist(cif.atoms);
    vector<int> mol(cif.natoms);
    transform(cif.atoms.begin(), cif.atoms.end(), mol.begin(), [](const Atom& x) -> int { return x.mol_id(); });
    sort(mol.begin(), mol.end());
    auto last = unique(mol.begin(), mol.end());
    int nmols = last - mol.begin();
    cif.atoms.reserve(cif.natoms);

    for(int i = 1; i < n; ++i){
        for_each(addlist.begin(), addlist.end(),
            [&d, &nmols](Atom& x) {
                 x.frac_xyz( x.frac_xyz() + d );
                 x.mol_id( x.mol_id() + nmols ); });
        cif.atoms.insert(cif.atoms.end(), addlist.begin(), addlist.end());
    }
    cif.natoms *= n;

    return cif;
}

ofstream& output_details(ofstream& fout, xxx& x, Part& y, char type, const Cell& cif, const IsParallelOptions& options)
{
    char lp, rp;
    //cout << "1234567890123456789012345678901234567890" << endl;
    //cout << parts(2*a-5*b-c, 1) << endl;
    //fout << "12345678901234567890123456789012345678901234567890123456789012345678901234567890" << endl;
    if( type == 'o' ) {
        fout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	fout << " Crystal Orientation uvw: ";
	lp = '<';
	rp = '>';
    }
    if( type == 'p' ) {
        fout << "============================================================" << endl;
	fout << " Crystal Plane hkl: ";
	lp = '(';
	rp = ')';
    }

    fout << lp << " "
	 << right
	 << setw(5) << x.h << " "
	 << setw(5) << x.k << " "
	 << setw(5) << x.l << " "
	 << rp << " "
         << "Parts: " << setw(3) << y.count << endl;

    fout << " Normal Vector: " << x.n << endl;

    if( options.flag_vdw_dist   ) {
	fout << " Minimum van der Waals Contact Distance: "
	     << cif.atoms[y.vdw_dist.atom_i].element_symbol()
	     << " ... "
	     << cif.atoms[y.vdw_dist.atom_j].element_symbol() << " "
	     << setw(6) << fixed << setprecision(3) << y.vdw_dist.value << endl;
    }
    if( options.flag_vdw_scalar ) {
	fout << " Minimum van der Waals Contact Scalar: "
	     << cif.atoms[y.vdw_scalar.atom_i].element_symbol()
	     << " ... "
	     << cif.atoms[y.vdw_scalar.atom_j].element_symbol() << " "
	     << setw(6) << fixed << setprecision(3) << y.vdw_scalar.value << endl;
    }
    if( options.flag_cov_dist   ) {
	fout << " Minimum Covalent Contact Distance: "
	     << cif.atoms[y.cov_dist.atom_i].element_symbol()
	     << " ... "
	     << cif.atoms[y.cov_dist.atom_j].element_symbol() << " "
	     << setw(6) << fixed << setprecision(3) << y.cov_dist.value << endl;
    }
    if( options.flag_cov_scalar ) {
	fout << " Minimum Covalent Contact Scalar: "
	     << cif.atoms[y.cov_scalar.atom_i].element_symbol()
	     << " ... "
	     << cif.atoms[y.cov_scalar.atom_j].element_symbol() << " "
	     << setw(6) << fixed << setprecision(3) << y.cov_scalar.value << endl;
    }
    
    if( type == 'o' ) {
        fout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    }
    if( type == 'p' ) {
        fout << "============================================================" << endl;
    }
    fout << endl;
    
    return fout;
}
