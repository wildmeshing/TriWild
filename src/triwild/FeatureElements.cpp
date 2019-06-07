// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//


#include "FeatureElements.h"
#include "optimization.h"
#include <cmath>
#include <sstream>

void triwild::feature::FeatureElement::merge_after(const triwild::feature::FeatureElement &other){
    assert(curve_id == other.curve_id);
    assert(type == other.type);
    assert(degree == other.degree);

    assert(fabs(paras.back() - other.paras.front()) < 1e-15);
    assert(v_ids.back() == other.v_ids.front());


    auto other_v_ids_pos = other.v_ids.begin();
    //skip the first one since it the same
    other_v_ids_pos++;
    v_ids.insert(v_ids.end(), other_v_ids_pos, other.v_ids.end());


    auto other_paras_pos = other.paras.begin();
    //skip the first one since it the same
    other_paras_pos++;
    paras.insert(paras.end(), other_paras_pos, other.paras.end());
}

std::string triwild::feature::FeatureElement::to_eps() const
{
    auto p0 = eval(paras.front());
    auto p1 = eval(paras.back());

    std::stringstream ss;
    if(is_inflection[0])
        ss<<"0 0 0"<<" setrgbcolor\n";
    else
        ss<<"255 0 0"<<" setrgbcolor\n";
    ss << p0[0] << " " << p0[1]  << " moveto\n";
    ss << p0[0] << " " << p0[1] << " " << 3 << " 0 360 arc\n";


    if(is_inflection[1])
        ss<<"0 0 0"<<" setrgbcolor\n";
    else
        ss<<"1 0 0"<<" setrgbcolor\n";
    ss << p1[0] << " " << p1[1]  << " moveto\n";
    ss << p1[0] << " " << p1[1] << " " << 3 << " 0 360 arc\n";

    ss << "fill\n";

    return ss.str();
}

/////////////////////Line
triwild::Point_2f triwild::feature::Line_Feature::eval(double t) const {
    if(t<0 || t>1) {
        cout<<"t<0 || t>1"<<endl;
        optimization::pausee();
    }
    assert(t>=0 && t<=1);
    return start * (1 - t) + end * t;
}

triwild::Point_2f triwild::feature::Line_Feature::eval_first_derivative(double t) const {
    return end - start;
}

triwild::Point_2f triwild::feature::Line_Feature::eval_second_derivative(double t) const {
    return Point_2f(0, 0);
}

double triwild::feature::Line_Feature::inv_eval(Point_2f& p, double t, double t0, double t1) const {
    Vector_2f v1 = p - start;
    Vector_2f v2 = end - start;

    double t_tmp = v1.dot(v2) / v2.length_2();
//    double t_tmp = v1.dot(v2) / (v2.length() * v1.length()) * v1.length() / v2.length();

    if (t_tmp < t0 || t_tmp > t1)
        return t;

    return t_tmp;
//    return std::sqrt((p - start).length_2()) / std::sqrt((start - end).length_2());
}

// double triwild::feature::Line_Feature::distance(Point_2f& p, double t) const {
//     return 0;
// }

bool triwild::feature::Line_Feature::is_point() const {
    if ((start - end).length_2() < 1e-12)
        return true;
    return false;
}

bool triwild::feature::Line_Feature::is_too_short(const Eigen::MatrixXd& V, double eps) const {
    if ((start - end).length_2() < eps * eps)
        return true;
    return false;
}

std::string triwild::feature::Line_Feature::to_eps() const
{
    auto p0 = eval(paras.front());
    auto p1 = eval(paras.back());

    std::stringstream ss;

    ss<<"0 0 0"<<" setrgbcolor\n";
    ss<<p0[0]<<" "<<p0[1]<<" moveto\n";
    ss<<p1[0]<<" "<<p1[1]<<" lineto\n";
    ss<<"stroke\n";

    ss<<FeatureElement::to_eps();

    return ss.str();
}

/////////////////////RationalBezier
triwild::Point_2f triwild::feature::RationalBezierCurve_Feature::eval(double t) const {
    std::array<double, 2> p = RationalBezier::interpolate(poles, weights, t);
    return Point_2f(p[0], p[1]);
}

triwild::Point_2f triwild::feature::RationalBezierCurve_Feature::eval_first_derivative(double t) const {
    std::array<double, 2> p = RationalBezier::first_derivative(poles, weights, t);
    return Point_2f(p[0], p[1]);
}

triwild::Point_2f triwild::feature::RationalBezierCurve_Feature::eval_second_derivative(double t) const {
    std::array<double, 2> p = RationalBezier::second_derivative(poles, weights, t);
    return Point_2f(p[0], p[1]);
}

double triwild::feature::RationalBezierCurve_Feature::inv_eval(Point_2f& p, double t, double t0, double t1) const {
    double tmp_t = RationalBezier::inverse_interpolation(poles, weights, std::array<double, 2>({p[0], p[1]}), t);
    if (tmp_t > t0 && tmp_t < t1)
        return tmp_t;
    tmp_t = RationalBezier::inverse_interpolation(poles, weights, std::array<double, 2>({p[0], p[1]}), t0);
    if (tmp_t > t0 && tmp_t < t1)
        return tmp_t;
    tmp_t = RationalBezier::inverse_interpolation(poles, weights, std::array<double, 2>({p[0], p[1]}), t1);
    if (tmp_t > t0 && tmp_t < t1)
        return tmp_t;
    assert(false);
    return t;
}

// double triwild::feature::RationalBezierCurve_Feature::distance(Point_2f& p, double t) const {
//     return 0;
// }

double triwild::feature::RationalBezierCurve_Feature::how_curve(double t0, double t1) const {
    if (t1 - t0 < 1e-7)
        return 0;

    Point_2f n0 = eval_first_derivative(t0); n0 = n0 / n0.length();
    Point_2f n1 = eval_first_derivative(t1); n1 = n1 / n1.length();

    double cos_a = n0.dot(n1);
    if (cos_a > 1)
        cos_a = 1;
    if (cos_a < -1)
        cos_a = -1;

    const double angle = std::acos(cos_a);

    return angle / M_PI * 180;
}

double triwild::feature::RationalBezierCurve_Feature::get_cut_t(const double t0, const double t1, bool to_flip) const
{
    assert(!to_flip);

    auto tan0 = eval_first_derivative(t0);
    auto tan1 = eval_first_derivative(t1);
    if(tan0.length_2() < 1e-15 || tan1.length_2() < 1e-15)
        return (t0+t1)/2;

    tan0 = tan0/sqrt(tan0.length_2());
    tan1 = tan1/sqrt(tan1.length_2());

    const double tan0x = tan0.x;
    const double tan0y = tan0.y;

    const double tan1x = tan1.x;
    const double tan1y = tan1.y;

    static const std::array<int, 3> x = {{0, 2, 4}};
    static const std::array<int, 3> y = {{1, 3, 5}};

    const double cx0 = poles[x[0]];
    const double cy0 = poles[y[0]];

    const double cx1 = poles[x[1]];
    const double cy1 = poles[y[1]];

    const double cx2 = poles[x[2]];
    const double cy2 = poles[y[2]];

    const double w0 = weights[0];
    const double w1 = weights[1];
    const double w2 = weights[2];

    const double a = cx0*tan0y*w0*w1-cx0*tan0y*w0*w2+cx0*tan1y*w0*w1-cx0*tan1y*w0*w2-cx1*tan0y*w0*w1+cx1*tan0y*w1*w2-cx1*tan1y*w0*w1+cx1*tan1y*w1*w2+cx2*tan0y*w0*w2-cx2*tan0y*w1*w2+cx2*tan1y*w0*w2-cx2*tan1y*w1*w2-cy0*tan0x*w0*w1+cy0*tan0x*w0*w2-cy0*tan1x*w0*w1+cy0*tan1x*w0*w2+cy1*tan0x*w0*w1-cy1*tan0x*w1*w2+cy1*tan1x*w0*w1-cy1*tan1x*w1*w2-cy2*tan0x*w0*w2+cy2*tan0x*w1*w2-cy2*tan1x*w0*w2+cy2*tan1x*w1*w2;
    const double b = 2*w0*(((cy0-cy1)*tan0x+(-cx0+cx1)*tan0y+(cy0-cy1)*tan1x-tan1y*(cx0-cx1))*w1-(1./2.)*w2*((-cy2+cy0)*tan0x+(-cx0+cx2)*tan0y+(-cy2+cy0)*tan1x-tan1y*(cx0-cx2)));
    const double c = w1*((-cy0+cy1)*tan0x+(cx0-cx1)*tan0y+(-cy0+cy1)*tan1x+tan1y*(cx0-cx1))*w0;

    const double delta = b*b - 4*a*c;
    if (delta < 0)
        return (t0 + t1) / 2;

    double res = -1;

    if(fabs(a) < 1e-8) {
        res = fabs(b) < 1e-8 ? ((t0 + t1)/2): (-c/b);
    }
    else
    {
//        if(fabs(delta) < 1e-16) {
//            if(fabs(a) < 1e-8)
//                res = (t0 +t1) / 2.;
//            else
//                res = -b / (2 * a);
//        }
//        else
//        {
            assert(delta >= 0);
            const double delta_sqrt = sqrt(delta);

            const double sol1 = (-b - delta_sqrt)/(2*a);
            const double sol2 = (-b + delta_sqrt)/(2*a);

            if(sol1 >= t0 && sol1 <= t1)
                res = sol1;
            else
                res = sol2;
//        }
    }

//    assert(res >= t0 && res <= t1);
    if(res < t0 || res > t1){
        std::cout<<"Warning: res < t0 && res > t1"<<std::endl;
        return (t0+t1)/2;
    }

    return res;
}


void triwild::feature::RationalBezierCurve_Feature::print_info() const
{
    std::cout<<"nodes=[";
    for(int i = 0; i < poles.size(); i += 2)
    {
        std::cout << poles[i] <<", "<< poles[i+1] <<"\n";
    }

    std::cout<<"]\nweights=[";
    for(int i = 0; i < weights.size(); ++i)
    {
        std::cout << weights[i] <<"\n";
    }
    std::cout<<"]\n"<<std::endl;
}

bool triwild::feature::RationalBezierCurve_Feature::is_point() const {
    for (int i = 2; i < poles.size(); i += 2) {
        if (Point_2(poles[i] - poles[i - 2], poles[i + 1] - poles[i - 2 + 1]).length_2() > 1e-12)
            return false;
    }
    return true;
}

bool triwild::feature::RationalBezierCurve_Feature::is_too_short(const Eigen::MatrixXd& V, double eps) const {
    double l = 0;
    for (int i = 1; i < v_ids.size(); i++) {
        l += (V.row(v_ids[i]) - V.row(v_ids[i - 1])).norm();
        if (l > eps)
            return false;
    }
    return true;
}

/////////////////////BezierCurve
triwild::Point_2f triwild::feature::BezierCurve_Feature::eval(double t) const {
    if(t<0 || t>1) {
        cout<<"t<0 || t>1"<<endl;
        optimization::pausee();
    }
    assert(t>=0 && t<=1);
    std::array<double, 2> p = Bezier::interpolate(poles, t);
    return Point_2f(p[0], p[1]);
}

triwild::Point_2f triwild::feature::BezierCurve_Feature::eval_first_derivative(double t) const {
    std::array<double, 2> p = Bezier::first_derivative(poles, t);
    return Point_2f(p[0], p[1]);
}

triwild::Point_2f triwild::feature::BezierCurve_Feature::eval_second_derivative(double t) const {
    std::array<double, 2> p = Bezier::second_derivative(poles, t);
    return Point_2f(p[0], p[1]);
}

double triwild::feature::BezierCurve_Feature::inv_eval(Point_2f& p, double t, double t0, double t1) const {
    return Bezier::point_curve_distance(poles, std::array<double, 2>({p[0], p[1]}), t0, t1);
    // return Bezier::inverse_interpolation(poles, std::array<double, 2>({p[0], p[1]}), t);
}

// double triwild::feature::BezierCurve_Feature::distance(Point_2f& p, double t) const {
//     return 0;
// }

double triwild::feature::BezierCurve_Feature::how_curve(double t0, double t1) const {
    if (t1 - t0 < 1e-7)
        return 0;

    ControlVector new_poles;
    Bezier::resample(poles, t0, t1, new_poles);

    std::vector<Vector_2f> vs;
    for (int i = 2; i < new_poles.size(); i += 2) {
        Vector_2f v(new_poles[i] - new_poles[i - 2], new_poles[i + 1] - new_poles[i - 2 + 1]);
//        double l = v.length();
//        if (l == 0)
        if (v.length_2() < 1e-16)
            continue;
//        std::cout<<l<<std::endl;
        vs.push_back(v / v.length());
    }

    double angle = 0;
    for (int i = 0; i < int(vs.size()) - 1; i++) {
        double cos_a = vs[i].dot(vs[i + 1]);
        if (cos_a > 1)
            cos_a = 1;
        if (cos_a < -1)
            cos_a = -1;

        angle += std::acos(cos_a);
    }

    return angle / M_PI * 180;
}

std::vector<double> triwild::feature::BezierCurve_Feature::inflection_points(const double t0, const double t1) const {
    if(poles.size() == 6)
    {
        //Quadratic, no inflection points
        return std::vector<double>();
    }

    std::vector<double> res;

    static const std::array<int, 4> x = {{0, 2, 4, 6}};
    static const std::array<int, 4> y = {{1, 3, 5, 7}};

    const double cx0 = poles[x[0]];
    const double cy0 = poles[y[0]];

    const double cx1 = poles[x[1]];
    const double cy1 = poles[y[1]];

    const double cx2 = poles[x[2]];
    const double cy2 = poles[y[2]];

    const double cx3 = poles[x[3]];
    const double cy3 = poles[y[3]];


    const double a = (-18*cy1+36*cy2-18*cy3)*cx0+(18*cy0-54*cy2+36*cy3)*cx1+(-36*cx2+18*cx3)*cy0+(54*cx2-36*cx3)*cy1-18*cx2*cy3+18*cx3*cy2;
    const double b = (36*cy1-54*cy2+18*cy3)*cx0+(-36*cy0+54*cy2-18*cy3)*cx1+(54*(cy0-cy1))*(cx2-(1./3.)*cx3);
    const double c = (-18*cy1+18*cy2)*cx0+(18*cy0-18*cy2)*cx1-18*cx2*(cy0-cy1);

    const double delta = b*b - 4*a*c;

    //no real solutions
    if(delta < 0)
        return res;

    if(fabs(a)< 1e-15)
    {
        if(fabs(b)<1e-15)
            return  res;


        double sol = -c/b;

        if(fabs(sol-t0) < 1e-15)
            sol=t0;
        if(fabs(sol-t1) < 1e-15)
            sol=t1;

        if(sol > t0 && sol < t1)
            res.push_back(sol);
        return res;
    }

    //only one solution
    // if(delta < 1e-8)
    // {
    //     double sol = -b/(2*a);

    //     if(fabs(sol-t0) < 1e-15)
    //         sol=t0;
    //     if(fabs(sol-t1) < 1e-15)
    //         sol=t1;

    //     //sol in in range
    //     if(sol > t0 && sol < t1)
    //         res.push_back(sol);

    //     return res;
    // }



    const double delta_sqrt = sqrt(delta);

    double sol1 = (-b - delta_sqrt)/(2*a);
    double sol2 = (-b + delta_sqrt)/(2*a);

    if(fabs(sol1-t0) < 1e-15)
        sol1=t0;

    if(fabs(sol1-t1) < 1e-15)
        sol1=t1;

    if(fabs(sol2-t1) < 1e-15)
        sol2=t1;

    if(fabs(sol2-t0) < 1e-15)
        sol2=t0;

    if(sol1 > t0 && sol1 < t1)
        res.push_back(sol1);

    if(fabs(sol1-sol2) > 1e-8)
    {
        if(sol2 > t0 && sol2 < t1)
            res.push_back(sol2);
    }

    return res;
}

double triwild::feature::BezierCurve_Feature::get_cut_t(const double t0, const double t1, bool to_flip) const {
    auto tan0 = eval_first_derivative(t0);
    auto tan1 = eval_first_derivative(t1);
    if (tan0.length_2() < 1e-15 || tan1.length_2() < 1e-15)
        return (t0 + t1) / 2;

    tan0 = tan0 / sqrt(tan0.length_2());
    tan1 = tan1 / sqrt(tan1.length_2());

    if (to_flip) {
        tan0 = tan0 * (-1);
        tan1 = tan1 * (-1);
    }

    const double tan0x = tan0.x;
    const double tan0y = tan0.y;

    const double tan1x = tan1.x;
    const double tan1y = tan1.y;

    if (poles.size() == 6) {
        static const std::array<int, 3> x = {{0, 2, 4}};
        static const std::array<int, 3> y = {{1, 3, 5}};

        const double cx0 = poles[x[0]];
        const double cy0 = poles[y[0]];

        const double cx1 = poles[x[1]];
        const double cy1 = poles[y[1]];

        const double cx2 = poles[x[2]];
        const double cy2 = poles[y[2]];

        const double m = (2 * cy0 - 4 * cy1 + 2 * cy2) * tan0x + (-2 * cx0 + 4 * cx1 - 2 * cx2) * tan0y +
                         (2 * cy0 - 4 * cy1 + 2 * cy2) * tan1x - 2 * tan1y * (cx0 - 2 * cx1 + cx2);
        const double b = (-2 * cy0 + 2 * cy1) * tan0x + (2 * cx0 - 2 * cx1) * tan0y + (-2 * cy0 + 2 * cy1) * tan1x +
                         2 * tan1y * (cx0 - cx1);

        const double res = fabs(m) < 1e-8 ? ((t0 + t1) / 2) : (-b / m);

//        assert(res >= t0 && res <= t1);
        if(res < t0 || res > t1){
            std::cout<<"Warning: res < t0 && res > t1"<<std::endl;
            return (t0+t1)/2;
        }

        return res;
    } else {
        static const std::array<int, 4> x = {{0, 2, 4, 6}};
        static const std::array<int, 4> y = {{1, 3, 5, 7}};

        const double cx0 = poles[x[0]];
        const double cy0 = poles[y[0]];

        const double cx1 = poles[x[1]];
        const double cy1 = poles[y[1]];

        const double cx2 = poles[x[2]];
        const double cy2 = poles[y[2]];

        const double cx3 = poles[x[3]];
        const double cy3 = poles[y[3]];


        double a = 3 * cx0 * tan0y + 3 * cx0 * tan1y - 9 * cx1 * tan0y - 9 * cx1 * tan1y + 9 * cx2 * tan0y +
                         9 * cx2 * tan1y - 3 * cx3 * tan0y - 3 * cx3 * tan1y - 3 * cy0 * tan0x - 3 * cy0 * tan1x +
                         9 * cy1 * tan0x + 9 * cy1 * tan1x - 9 * cy2 * tan0x - 9 * cy2 * tan1x + 3 * cy3 * tan0x +
                         3 * cy3 * tan1x;
        double b = (6 * cy0 - 12 * cy1 + 6 * cy2) * tan0x + (-6 * cx0 + 12 * cx1 - 6 * cx2) * tan0y +
                         (6 * cy0 - 12 * cy1 + 6 * cy2) * tan1x - 6 * tan1y * (cx0 - 2 * cx1 + cx2);
        double c = (-3 * cy0 + 3 * cy1) * tan0x + (3 * cx0 - 3 * cx1) * tan0y + (-3 * cy0 + 3 * cy1) * tan1x +
                         3 * tan1y * (cx0 - cx1);


        const double nnn = (fabs(a)+fabs(b)+fabs(c))/3.0;

        if(nnn > 1e-10)
        {
            a/=nnn;
            b/=nnn;
            c/=nnn;
        }

        const double delta = b * b - 4 * a * c;

        if (delta < 0)
            return (t0 + t1) / 2;

        double res = -1;

        if (fabs(a) < 1e-8) {
            res = fabs(b) < 1e-8 ? ((t0 + t1) / 2) : (-c / b);
        } else {
//            if (fabs(delta) < 1e-16) {
//                if (fabs(a) < 1e-8)
//                    res = (t0 + t1) / 2.;
//                else
//                    res = -b / (2 * a);
//            } else {
            assert(delta >= 0);
            const double delta_sqrt = sqrt(delta);

            const double sol1 = (-b - delta_sqrt) / (2 * a);
            const double sol2 = (-b + delta_sqrt) / (2 * a);

            if (sol1 >= t0 && sol1 <= t1)
                res = sol1;
            else
                res = sol2;
//            }
        }

//        assert(res >= t0 && res <= t1);
        if(res < t0 || res > t1){
            std::cout<<"Warning: res < t0 && res > t1"<<std::endl;
            return (t0+t1)/2;
        }

        return res;

    }
}

#include <iomanip>
std::string triwild::feature::BezierCurve_Feature::to_maple() const
{
    std::stringstream res;
    res<<std::setprecision(16);

    for(int i = 0; i < poles.size(); i += 2)
    {
        res << "cx" << (i/2) << " := " << poles[i] <<"; ";
        res << "cy" << (i/2) << " := " << poles[i+1] <<"; ";
    }

    return res.str();
}

void triwild::feature::BezierCurve_Feature::print_info() const
{
    std::cout<<"nodes=[";
    for(int i = 0; i < poles.size(); i += 2)
    {
        std::cout << poles[i] <<", "<< poles[i+1] <<"\n";
    }
    std::cout<<"]\n"<<std::endl;
}

bool triwild::feature::BezierCurve_Feature::is_point() const {
    for (int i = 2; i < poles.size(); i += 2) {
        if (Point_2(poles[i] - poles[i - 2], poles[i + 1] - poles[i - 2 + 1]).length_2() > 1e-12)
            return false;
    }
    return true;
}

bool triwild::feature::BezierCurve_Feature::is_too_short(const Eigen::MatrixXd& V, double eps) const {
    double l = 0;
    for (int i = 1; i < v_ids.size(); i++) {
        l += (V.row(v_ids[i]) - V.row(v_ids[i - 1])).norm();
        if (l > eps)
            return false;
    }
    return true;
}


std::shared_ptr<triwild::feature::FeatureElement> triwild::feature::BezierCurve_Feature::simplify(
        const Point_2f& start, const Point_2f& end) const {
    //collinear
    std::vector<Vector_2f> vs;

    ControlVector new_poles;
    Bezier::resample(poles, paras.front(), paras.back(), new_poles);

    Point_2f p(new_poles[0], new_poles[1]);
    for (int i = 2; i < new_poles.size(); i += 2) {
        Vector_2f v = p - Point_2f(new_poles[i], new_poles[i + 1]);
        double l = v.length();
//        if (l == 0)
        if (l < 1e-8)
            continue;
        vs.push_back(v / l);
    }
    bool is_parallel = true;
    for (int i = 1; i < vs.size(); i++) {
        double abs_cos_a = std::abs(vs[0].dot(vs[i]));
        if (!(abs_cos_a > 1 - 1e-7 && abs_cos_a < 1 + 1e-7)) {
            is_parallel = false;
            break;
        }
    }
    if (!is_parallel)
        return NULL;

//    return std::make_shared<Line_Feature>(v_ids, paras, Point_2f(new_poles[0], new_poles[1]), Point_2f(new_poles[new_poles.size()-2],new_poles[new_poles.size()-1]), curve_id);
    return std::make_shared<Line_Feature>(v_ids, paras, start, end, curve_id);

    //sample
//    const int N = v_ids.size() * 50;
//    double t1 = paras.front();
//    double t2 = paras.back();
//    std::vector<Point_2f> ps;
//    ps.reserve(N);
//    Point_2f min, max;
//    for (double n = 0; n <= N; n++) {
//        ps.push_back(eval((1 - n / N) * t1 + (n / N) * t2));
//        if (n == 0 || ps.back().x < min.x)
//            min.x = ps.back().x;
//        if (n == 0 || ps.back().y < min.y)
//            min.y = ps.back().y;
//        if (n == 0 || ps.back().x > max.x)
//            max.x = ps.back().x;
//        if (n == 0 || ps.back().y > max.y)
//            max.y = ps.back().y;
//    }
//
//    Point_2f min_x_min_y(min.x, min.y), min_x_max_y(min.x, max.y);
//    double l_min_x_min_y;
//    double l_min_x_max_y;
//    for (int i = 0; i < ps.size(); i++) {
//        double l = (min_x_min_y - ps[i]).length_2();
//        if (i == 0 || l < l_min_x_min_y)
//            l_min_x_min_y = l;
//
//        l = (min_x_max_y - ps[i]).length_2();
//        if (i == 0 || l < l_min_x_max_y)
//            l_min_x_max_y = l;
//    }
//    if (l_min_x_min_y < l_min_x_max_y)
//        return std::make_shared<Line_Feature>(v_ids, paras, min, max, curve_id);
//    else
//        return std::make_shared<Line_Feature>(v_ids, paras, Point_2f(min.x, max.y), Point_2f(max.x, min.y), curve_id);
}

std::string triwild::feature::BezierCurve_Feature::to_eps() const
{
    ControlVector new_poles;
    Bezier::resample(poles, paras.front(), paras.back(), new_poles);

    std::stringstream ss;

    ss<<"0 0 0"<<" setrgbcolor\n";
    ss<<new_poles(0)<<" "<<new_poles(1)<<" moveto\n";

    for(int i = 2; i < new_poles.size(); i+=2)
        ss<<new_poles(i)<<" "<<new_poles(i+1)<<" ";

    ss<<"curveto\n";
    ss<<"stroke\n";

    ss<<FeatureElement::to_eps();

    return ss.str();
}







