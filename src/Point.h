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
/*  Description: Defination of class Point.                                   */
/*                                                                            */
/******************************************************************************/

#include <cmath>
#include <iostream>
#include <iomanip>

#ifndef POINT_H

#define POINT_H

using namespace std;

class Point
{
    private:
        double X, Y, Z;

    public:
        const static double ABSZERO;

        Point();
        Point(const Point& v);
        Point(const double x, const double y, const double z);
        double length() const;
        Point& normalize();

        friend Point operator + (const Point& a, const Point& b);
        friend Point operator - (const Point& a, const Point& b);
        friend Point operator * (const double l, const Point& v);
        friend Point operator * (const Point& v, const double l);
        friend Point operator / (const Point& v, const double l);
        friend double dot(const Point& a, const Point& b);
        friend Point cross(const Point& a, const Point& b);
        friend Point floor(const Point& v);
        friend double abs(const Point& v);
        friend double det(const Point& a, const Point& b, const Point& c);
        friend istream& operator >> (istream& input, Point& x);
        friend ostream& operator << (ostream& output, const Point& x);

        Point& operator = (const Point& v);
        Point& operator += (const Point& v);
        Point& operator -= (const Point& v);
        Point& operator *= (const double l);

        Point& x(double value);
        Point& y(double value);
        Point& z(double value);

        double x();
        double y();
        double z();

};

#endif
/*
// class Atom tests
    Atom H(1, "H", 0, 0, 0, 1.008, 0.32, 1.10);
    Atom C(6, "C", 0, 0, 1, 12.01, 0.77, 1.77);
    cout << H.mass() << endl;
    H.mass(2);
    cout << H.mass() << endl;
    cout << C.atomic_number() << endl;
    C.atomic_number(7);
    cout << C.atomic_number() << endl;
*/
