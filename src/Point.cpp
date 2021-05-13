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

#include "Point.h"

using namespace std;

Point::Point()
{
    X = Y = Z = 0.0l;
}

Point::Point(const Point& v)
{
    X = v.X;
    Y = v.Y;
    Z = v.Z;
}

Point::Point(const double x, const double y, const double z)
{
    X = x;
    Y = y;
    Z = z;
}

double Point::length() const
{
    return sqrt( X * X + Y * Y + Z * Z );
}

Point& Point::normalize()
{
    double l = length();
    X /= l;
    Y /= l;
    Z /= l;
    return *this;
}

Point& Point::operator = (const Point& v)
{
    if( this == &v)return *this;
    X = v.X;
    Y = v.Y;
    Z = v.Z;
    return *this;
}

Point& Point::operator += (const Point& v)
{
    X += v.X;
    Y += v.Y;
    Z += v.Z;
    return *this;
}

Point& Point::operator -= (const Point& v)
{
    X -= v.X;
    Y -= v.Y;
    Z -= v.Z;
    return *this;
}

Point& Point::operator *= (const double l)
{
    X *= l;
    Y *= l;
    Z *= l;
    return *this;
}

Point& Point::x(double value){X = value; return *this;}
Point& Point::y(double value){Y = value; return *this;}
Point& Point::z(double value){Z = value; return *this;}

double Point::x(){return X;}
double Point::y(){return Y;}
double Point::z(){return Z;}

const double Point::ABSZERO = 1e-6;

Point operator + (const Point& a, const Point& b)
{
    return Point(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
}

Point operator - (const Point& a, const Point& b)
{
    return Point(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
}

Point operator * (const double l, const Point& v)
{
    return Point(v.X * l, v.Y * l, v.Z * l);
}

Point operator * (const Point& v, const double l)
{
    return Point(v.X * l, v.Y * l, v.Z * l);
}

Point operator / (const Point& v, const double l)
{
    return Point(v.X / l, v.Y / l, v.Z / l);
}
double dot(const Point& a, const Point& b)
{
    return a.X * b.X + a.Y * b.Y + a.Z * b.Z;
}

Point cross(const Point& a, const Point& b)
{
    return Point(a.Y * b.Z - b.Y * a.Z,
                 a.Z * b.X - b.Z * a.X,
                 a.X * b.Y - b.X * a.Y);
}

Point floor(const Point& v)
{
    return Point(floor(v.X), floor(v.Y), floor(v.Z));
}

double abs(const Point& v)
{
    return sqrt( v.X * v.X + v.Y * v.Y + v.Z * v.Z );
}

//typedef vector<vector<double>> matrix;

double det(const Point& a, const Point& b, const Point& c)
{
    return a.X * b.Y * c.Z + b.X * c.Y * a.Z + c.X * a.Y * b.Z - c.X * b.Y * a.Z - b.X * a.Y * c.Z - a.X * c.Y * b.Z;
}

Point solve(const Point& a, const Point& b, const Point& c, const Point& v)
{
    double d = det(a, b, c);
    return Point(det(v, b, c)/d, det(a, v, c)/d, det(a, b, v)/d);
}

istream& operator >> (istream& input, Point& x)
{
    input >> x.X >> x.Y >> x.Z;
    return input;
}

ostream& operator << (ostream& output, const Point& x)
{
    output << fixed << "("
           << setw(12) << setprecision(5) << x.X << ", "
           << setw(12) << setprecision(5) << x.Y << ", "
           << setw(12) << setprecision(5) << x.Z << ")"
           << defaultfloat; 
    //output << "(" << x.X << ", " << x.Y << ", " << x.Z << ")";
    return output;
}

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
