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

#include "Atom.h"

using namespace std;

Atom::Atom()
{
    AtomicNumber = 0;
    ElementSymbol = "";
    FracXYZ;
    Mass = 0.0;
    CovRadius = 0.0;
    VDWRadius = 0.0;
    MolID = 0;
}

Atom::Atom(const Atom& at)
{
    AtomicNumber = at.AtomicNumber;
    ElementSymbol = at.ElementSymbol;
    FracXYZ = at.FracXYZ;
    Mass = at.Mass;
    CovRadius = at.CovRadius;
    VDWRadius = at.VDWRadius;
    MolID = at.MolID;
}

Atom::Atom(int an, const string& es, double x, double y, double z, double m, double cov, double vdw)
{
    AtomicNumber = an;
    ElementSymbol = es;
    FracXYZ = Point(x, y, z);
    //FracXYZ = Point(x, y, z) - floor(Point(x, y, z));
    Mass = m;
    CovRadius = cov;
    VDWRadius = vdw;
    MolID = 0;
}

Atom::Atom(int an, const string& es, double x, double y, double z, double m, double cov, double vdw, int id)
{
    AtomicNumber = an;
    ElementSymbol = es;
    FracXYZ = Point(x, y, z);
    //FracXYZ = Point(x, y, z) - floor(Point(x, y, z));
    Mass = m;
    CovRadius = cov;
    VDWRadius = vdw;
    MolID = id;
}

Atom& Atom::atomic_number(int x){AtomicNumber = x; return *this;}
Atom& Atom::element_symbol(string x){ElementSymbol = x; return *this;}
Atom& Atom::frac_xyz(Point x){FracXYZ = x; return *this;}
Atom& Atom::mass(double x){Mass = x; return *this;}
Atom& Atom::cov_radius(double x){CovRadius = x; return *this;}
Atom& Atom::vdw_radius(double x){VDWRadius = x; return *this;}
Atom& Atom::mol_id(int x){MolID = x; return *this;}

int Atom::atomic_number() const {return AtomicNumber;} 
string Atom::element_symbol() const {return ElementSymbol;}
Point Atom::frac_xyz() const {return FracXYZ;}
double Atom::mass() const {return Mass;}
double Atom::cov_radius() const {return CovRadius;}
double Atom::vdw_radius() const {return VDWRadius;}
int Atom::mol_id() const {return MolID;}

