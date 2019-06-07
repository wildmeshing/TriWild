// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef TRIWILD_POINT_2F_H
#define TRIWILD_POINT_2F_H

#include <iostream>
#include <algorithm>
#include <math.h>
#include <cmath>

namespace triwild {

    class Point_2f {
    public:
        double x;
        double y;

        Point_2f() {}

        Point_2f(double x_, double y_) {
            x = x_;
            y = y_;
        }

        Point_2f(const Point_2f& pf) {
            x = pf.x;
            y = pf.y;
        }

        void set(double x_, double y_){
            x = x_;
            y = y_;
        }

        //+, - another point
        Point_2f operator+(const Point_2f &p) const {
            return Point_2f(x + p.x, y + p.y);
        }

        Point_2f operator-(const Point_2f &p) const {
            return Point_2f(x - p.x, y - p.y);
        }

//    friend const Point_2f operator-(const Point_2f& p1,const Point_2f& p2){
//        return Point_2f(p1.x-p2.x, p1.y-p2.y);
//    }
        //*, / double/rational
        Point_2f operator*(const double a) const {
            return Point_2f(x * a, y * a);
        }

        Point_2f operator/(const double a) const {
            return Point_2f(x / a, y / a);
        }

        //[]
        double &operator[](int i) {
            if (i == 0)
                return x;
            else
                return y;
        }

        double operator[](int i) const {
            if (i == 0)
                return x;
            else
                return y;
        }

        //=
        void operator=(const Point_2f &p) {
            x = p.x;
            y = p.y;
        }

        //<<
        friend std::ostream &operator<<(std::ostream &os, const Point_2f &p) {
            os << p.x << ", " << p.y;
            return os;
        }

        //==, !=
        friend bool operator==(const Point_2f &p1, const Point_2f &p2) {
            return p1.x == p2.x && p1.y == p2.y;
        }

        friend bool operator!=(const Point_2f &p1, const Point_2f &p2) {
            return p1.x != p2.x || p1.y != p2.y;
        }

        //
        double length_2() const {
            return x * x + y * y;
        }

        double length() const {
            return std::sqrt(x * x + y * y);
        }

        double dot(const Point_2f &p) const {
            return x * p.x + y * p.y;
        }

        double cross(const Point_2f &p) const {
            return x * p.y - y * p.x;
        }

        bool isnan() {
            return std::isnan(x) || std::isnan(y);
        }

        bool isinf() {
            return std::isinf(x) || std::isinf(y);
        }

        void normalize() {
            double l = std::sqrt(x * x + y * y);
            x /= l;
            y /= l;
        }
    };

//    class Vector_2f:public Point_2f {
//    public:
//        Vector_2f() {}
//
//        Vector_2f(double x_, double y_) {
//            x = x_;
//            y = y_;
//        }
//    };
    typedef Point_2f Vector_2f;
}

#endif //TRIWILD_POINT_2F_H
