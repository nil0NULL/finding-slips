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
/*  Description: Defination of class Atom.                                    */
/*                                                                            */
/******************************************************************************/

#include <cmath>
#include <string>

#ifndef POINT_H
#include "Point.h"
#endif

#ifndef ATOM_H
#define ATOM_H

using namespace std;

class Atom
{
    private:
        int AtomicNumber;
        string ElementSymbol;
        double Mass;
        double CovRadius;
        double VDWRadius;
        int MolID;
        Point FracXYZ;

    public:
        Atom();
        Atom(const Atom& at);
        Atom(int an, const string& es, double x, double y, double z, double m, double cov, double vdw);
        Atom(int an, const string& es, double x, double y, double z, double m, double cov, double vdw, int id);
        Atom& atomic_number(int x);
        Atom& element_symbol(string x);
        Atom& frac_xyz(Point x);
        Atom& mass(double x);
        Atom& cov_radius(double x);
        Atom& vdw_radius(double x);
        Atom& mol_id(int x);

        int atomic_number() const; 
        string element_symbol() const;
        Point frac_xyz() const;
        double mass() const;
        double cov_radius() const;
        double vdw_radius() const;
        int mol_id() const;
};

#endif
