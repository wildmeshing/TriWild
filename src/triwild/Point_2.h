// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef TRIWILD_POINT_2_H
#define TRIWILD_POINT_2_H

#include "Rational.h"

namespace triwild {

    class Point_2 {
    public:
        Rational x;
        Rational y;

        Point_2() {}

        Point_2(double x_, double y_) {
            x = x_;
            y = y_;
        }

        Point_2(Rational x_, Rational y_) {
            x = x_;
            y = y_;
        }

        Point_2(const Point_2& p) {
            x = p.x;
            y = p.y;
        }

//        ~Point_2() {}

        //+, - another point
        Point_2 operator+(const Point_2 &p) const {
            return Point_2(x + p.x, y + p.y);
        }

        Point_2 operator-(const Point_2 &p) const {
            return Point_2(x - p.x, y - p.y);
        }

        //*, / double/rational
        Point_2 operator*(double a) const {
            return Point_2(x * a, y * a);
        }

        Point_2 operator/(double a) const {
            return Point_2(x / a, y / a);
        }

        //[]
        Rational &operator[](int i) {
            if (i == 0)
                return x;
            else
                return y;
        }

        //=
//        void operator=(const Point_2 &p) {
//            x = p.x;
//            y = p.y;
//        }
        Point_2& operator=(const Point_2& p) {
            if (this == &p) return *this;
            x = p.x;
            y = p.y;
            return *this;
        }

        //<<
        friend std::ostream &operator<<(std::ostream &os, const Point_2 &p) {
            os << p.x << ", " << p.y;
            return os;
        }

        //==, !=
        friend bool operator==(const Point_2 &p1, const Point_2 &p2) {
            return p1.x == p2.x && p1.y == p2.y;
        }

        friend bool operator!=(const Point_2 &p1, const Point_2 &p2) {
            return p1.x != p2.x || p1.y != p2.y;
        }

        //
        Rational length_2() const {
            return x * x + y * y;
        }

        Rational dot(const Point_2 &p) const {
            return x * p.x + y * p.y;
        }
    };

    typedef Point_2 Vector_2;
}

#endif //TRIWILD_POINT_2_H
