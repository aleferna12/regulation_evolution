/*

Copyright 1996-2006 Roeland Merks

This file is part of Tissue Simulation Toolkit.

Tissue Simulation Toolkit is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

Tissue Simulation Toolkit is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tissue Simulation Toolkit; if not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
02110-1301 USA

*/
#include <vector>
#include <sstream>
#define FALSE    0
#define TRUE     1
#define OK 1
#define REMARK 56


int ReadNumber(FILE *file, int *number);

int YesNoP(const char *message);

char *GetFileName(const char *message, const char *ftype);

int FileExists(FILE **fp, const char *fname, const char *ftype);

int ReadDouble(FILE *file, double *number);

//! Finds positive root of quadratic equation using Bhaskara's
double solveQuadradic(double a, double b, double c);


double dist(double x1, double y1, double x2, double y2);


//! Functions to interpret a string delimited by a character as a vector of values and vice-versa
//! When delim == \0 just concatenates the values
template<class T>
std::string vectorToString(const std::vector<T> &values, char delim = '\0') {
    std::stringstream out;
    for (auto it = values.begin(); it != values.end(); ++it) {
        out << *it;
        if (next(it) != values.end() and delim != '\0')
            out << delim;
    }
    return out.str();
}


template<class T> std::vector<T> stringToVector(const std::string &values, char delim) {
    std::vector<T> value_vec;
    std::stringstream values_ss(values);
    T value;
    std::string value_str;
    while (getline(values_ss, value_str, delim)) {
        std::stringstream (value_str) >> value;
        value_vec.push_back(value);
    }
    return value_vec;
}
template<> std::vector<std::string> stringToVector(const std::string &values, char delim);
