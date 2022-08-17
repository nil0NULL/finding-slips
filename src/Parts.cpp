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

#include "Parts.h"

using namespace std;

int gcd(int m, int n)
{
    m = abs(m);
    n = abs(n);
    if( n == 0 ) return m;
    for(int r = m % n; r !=0 ; r = m % n){
        m = n;
        n = r;
    }
    return n;
}

Part parts(const Cell& cif, function<double(const Point&)> fun, const IsParallelOptions& options)
{
    int findmol(int k, vector<int>& mol);
    vector<vector<double>> dist(cif.natoms + 1, vector<double>(cif.natoms + 1, 0));
    vector<int> mol(cif.natoms + 1);
    vector<int> root(cif.natoms + 1);
    if( options.flag_mol )
        transform(cif.atoms.begin(), cif.atoms.end(), mol.begin(), [](const Atom& x) -> int { return x.mol_id(); });
    else
        iota(mol.begin(), mol.end(), 0);
    for_each(mol.begin(), mol.end(), [&root](const int& x) { root.at(x) = x; });
    
    //for_each(mol.begin(), mol.end(), [](const int& x) { cout << x << " "; });
    //cout << endl;

    //cout << findmol(0, mol) << findmol(18, mol) << endl;

    //cout << mol.id[0] << mol.id[18] << endl;
    for(size_t i = 0; i < cif.natoms - 1 ; ++i)
        for(size_t j = i + 1; j < cif.natoms; ++j) {
            if( findmol(mol.at(i), root) == findmol(mol.at(j), root) ) continue;
            //cout << i << ", " << j << endl;
	    Point v = cif.atoms.at(i).frac_xyz() - cif.atoms.at(j).frac_xyz();
            v -= floor(v);
	    Point xyz = cif.a * v.x() + cif.b * v.y() + cif.c * v.z();
            //cout << cif.atoms.at(i).frac_xyz() << cif.atoms.at(j).frac_xyz() << v << xyz << endl;
            dist[i][j] = dist[j][i] = fun(xyz);
            //cout << i << ", " << j << ": " << dist[i][j] << endl;
            if( isclose(cif.atoms.at(i), cif.atoms.at(j), dist[i][j], options) ) 
                root.at(findmol(mol.at(j), root)) = findmol(mol.at(i), root);
        }

    /*for(int i = 0; i < cif.natoms ; ++i) {
        for(int j = 0 ; j < cif.natoms; ++j)
            cout << fixed << setw(6) << setprecision(2) << dist[i][j] << " " ;
        cout << endl;
    }
    cout << defaultfloat;*/
    //for_each(mol.begin(), mol.end(), [](int& x) { cout << x << " " ; });
            //cout << endl;
    set<int> mols;
    for(int i = 0; i < cif.natoms; ++i) {
        mol.at(i) = findmol(mol.at(i), root);
        mols.insert(mol.at(i));
    }
    //for_each(mol.begin(), mol.end(), [](int& x) { cout << x << " " ; });
            //cout << endl;
    Part ans{mols.size(), {10.0, 0, 0}, {10.0, 0, 0}, {10.0, 0, 0}, {10.0, 0, 0}};

    for(int i = 0; i < cif.natoms - 1; ++i)
        for(int j = i + 1; j < cif.natoms; ++j) {
            if( mol.at(i) == mol.at(j) ) continue;
            if( options.flag_vdw_dist   && ans.vdw_dist.value   > dist[i][j] - (cif.atoms[i].vdw_radius() + cif.atoms[j].vdw_radius()) ) {
                ans.vdw_dist.value    = dist[i][j] - (cif.atoms[i].vdw_radius() + cif.atoms[j].vdw_radius());
                ans.vdw_dist.atom_i   = i;
                ans.vdw_dist.atom_j   = j;
            }
            if( options.flag_vdw_scalar && ans.vdw_scalar.value > dist[i][j] / (cif.atoms[i].vdw_radius() + cif.atoms[j].vdw_radius()) ) {
                ans.vdw_scalar.value  = dist[i][j] / (cif.atoms[i].vdw_radius() + cif.atoms[j].vdw_radius());
                ans.vdw_scalar.atom_i = i;
                ans.vdw_scalar.atom_j = j;
            }
            if( options.flag_cov_dist   && ans.cov_dist.value   > dist[i][j] - (cif.atoms[i].cov_radius() + cif.atoms[j].cov_radius()) ) {
                ans.cov_dist.value    = dist[i][j] - (cif.atoms[i].cov_radius() + cif.atoms[j].cov_radius());
                ans.cov_dist.atom_i   = i;
                ans.cov_dist.atom_j   = j;
            }
            if( options.flag_cov_scalar && ans.cov_scalar.value > dist[i][j] / (cif.atoms[i].cov_radius() + cif.atoms[j].cov_radius()) ) {
                ans.cov_scalar.value  = dist[i][j] / (cif.atoms[i].cov_radius() + cif.atoms[j].cov_radius());
                ans.cov_scalar.atom_i = i;
                ans.cov_scalar.atom_j = j;
            }
        }

    return ans;
}

double calc_dist(const Point& p, function<double(const Point&)> fun, const vector<xxx>& cell_seq)
{
    double d = p.length();
    vector<double> len(cell_seq.size());
    transform(cell_seq.begin(), cell_seq.end(), len.begin(),
	[&p, fun](const xxx& x) -> double { return fun(p + x.n); });
    //for_each(mol.begin(), mol.end(), [](int& x) { cout << x << " " ; });
    //for_each(len.begin(), len.end(), [](double& x) { cout << x << " " ; });
    d = min(*min_element(len.begin(), len.end()), d);
    //cout << "d: " << d << endl;
    return d;
}

double dist_vector(const Point& vec, const Point& n)
{
    return abs(vec - dot(vec, n) / (abs(n) * abs(n)) * n); 
    //return sqrt(abs(vec) * abs(vec) - abs(dot(vec, n) / abs(n)) * abs(dot(vec, n) / abs(n)));
}

double dist_plane(const Point& vec, const Point& n1, const Point& n2)
{
    Point n = cross(n1, n2);
    return abs(dot(vec, n) / abs(n));
}

int findmol(int k, vector<int>& root)
{
    if( root.at(k) == k )
        return k;
    else
        return root.at(k) = findmol(root.at(k), root);
}

bool isclose(const Atom& a, const Atom& b, double dist, const IsParallelOptions& options)
{
    return (options.flag_vdw_dist   && dist < (a.vdw_radius() + b.vdw_radius()) - options.vdw_dist_minus )
        || (options.flag_vdw_scalar && dist < (a.vdw_radius() + b.vdw_radius()) * options.vdw_scalar_min )
        || (options.flag_cov_dist   && dist < (a.cov_radius() + b.cov_radius()) + options.cov_dist_plus  )
        || (options.flag_cov_scalar && dist < (a.cov_radius() + b.cov_radius()) * options.cov_scalar_max );
}

vector<xxx>& construct_cell_seq(const Cell& cif, vector<xxx>& cell_seq)
{
    cell_seq.clear();
    cell_seq.reserve(30);
    for(int i = -1; i <= 1; ++i)
        for(int j = -1; j <= 1; ++j)
            for(int k = -1; k <= 1; ++k)
                cell_seq.push_back({i, j, k, i * cif.a + j * cif.b + k * cif.c});

    return cell_seq;
}

vector<xxx>& construct_cell_seq(const Cell& cif, const xxx& n1, vector<xxx>& cell_seq)
{
    int m1 = max({abs(n1.h), abs(n1.k), abs(n1.l)});

    vector<xxx> delta_xyz(30);
    construct_cell_seq(cif, delta_xyz);

    for(int i = 0; i <= m1; ++i) {
        int h = n1.h * i / m1;
        int k = n1.k * i / m1;
        int l = n1.l * i / m1; 
	xxx nn = {h, k, l, h * cif.a + k * cif.b + l * cif.c};
	for_each(delta_xyz.begin(), delta_xyz.end(),  
	    [&nn, &cell_seq](const xxx& x) { cell_seq.push_back({nn.h + x.h, nn.k + x.k, nn.l + x.l, nn.n + x.n}); });
    }

    return cell_seq;
}

vector<xxx>& construct_cell_seq(const Cell& cif, const xxx& n1, const xxx& n2, vector<xxx>& cell_seq)
{
    int m1 = max({abs(n1.h), abs(n1.k), abs(n1.l)});
    int m2 = max({abs(n2.h), abs(n2.k), abs(n2.l)});

    vector<xxx> delta_xyz(30);
    construct_cell_seq(cif, delta_xyz);

    for(int i = 0; i <= m1; ++i)
	for(int j = 0; j <= m2; ++j) {
	    int h = n1.h * i / m1 + n2.h * j / m2;
	    int k = n1.k * i / m1 + n2.k * j / m2;
	    int l = n1.l * i / m1 + n2.l * j / m2; 
	    xxx nn = {h, k, l, h * cif.a + k * cif.b + l * cif.c};
	    for_each(delta_xyz.begin(), delta_xyz.end(),  
		[&nn, &cell_seq](const xxx& x) { cell_seq.push_back({nn.h + x.h, nn.k + x.k, nn.l + x.l, nn.n + x.n}); });
	}

    return cell_seq;
}

vector<xxx>& clear_cell_seq(vector<xxx>& cell_seq)
{
    auto xxx_lt = [](const xxx& a, const xxx& b) -> bool {
        if( a.h < b.h ) return true;
	else if( a.h == b.h ) {
	    if( a.k < b.k ) return true;
	    else if( a.k == b.k ){
		if( a.l < b.l ) return true;
		else return false;
	    }else return false;
	} else return false; };

    auto xxx_eq = [](const xxx& a, const xxx& b) -> bool {
        return (a.h == b.h && a.k == b.k && a.l == b.l); };

    sort(cell_seq.begin(), cell_seq.end(), xxx_lt);
    auto last = unique(cell_seq.begin(), cell_seq.end(), xxx_eq);
    cell_seq.erase(last, cell_seq.end());

    return cell_seq;
}

