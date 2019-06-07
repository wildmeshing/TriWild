// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "AMIPS.h"
#include "Curves.h"

double triwild::optimization::AMIPS_energy(const std::array<double, 6>& T){
    double helper_0[6];
    helper_0[0] = T[0];
    helper_0[1] = T[1];
    helper_0[2] = T[2];
    helper_0[3] = T[3];
    helper_0[4] = T[4];
    helper_0[5] = T[5];
    double helper_1 = helper_0[0];
    double helper_2 = helper_0[2];
    double helper_3 = helper_0[1];
    double helper_4 = helper_0[3];
    double helper_5 = helper_0[5];
    double helper_6 = helper_0[4];
    double helper_7 = 0.666666666666667*helper_6;
    double helper_8 = 0.666666666666667*helper_5;
    return -(helper_1*(-1.33333333333333*helper_1 + 0.666666666666667*helper_2 + helper_7) + helper_2*(0.666666666666667*helper_1 - 1.33333333333333*helper_2 + helper_7) + helper_3*(-1.33333333333333*helper_3 + 0.666666666666667*helper_4 + helper_8) + helper_4*(0.666666666666667*helper_3 - 1.33333333333333*helper_4 + helper_8) + helper_5*(0.666666666666667*helper_3 + 0.666666666666667*helper_4 - 1.33333333333333*helper_5) + helper_6*(0.666666666666667*helper_1 + 0.666666666666667*helper_2 - 1.33333333333333*helper_6))/((helper_1 - helper_2)*(0.577350269189626*helper_3 + 0.577350269189626*helper_4 - 1.15470053837925*helper_5) - (helper_3 - helper_4)*(0.577350269189626*helper_1 + 0.577350269189626*helper_2 - 1.15470053837925*helper_6));
}

void triwild::optimization::AMIPS_jacobian(const std::array<double, 6>& T, Eigen::Vector2d& result_0){
    double helper_0[6];
    helper_0[0] = T[0];
    helper_0[1] = T[1];
    helper_0[2] = T[2];
    helper_0[3] = T[3];
    helper_0[4] = T[4];
    helper_0[5] = T[5];
    double helper_1 = helper_0[0];
    double helper_2 = helper_0[2];
    double helper_3 = helper_0[1];
    double helper_4 = helper_0[3];
    double helper_5 = helper_0[5];
    double helper_6 = helper_0[4];
    double helper_7 = 1.0/((helper_1 - helper_2)*(0.577350269189626*helper_3 + 0.577350269189626*helper_4 - 1.15470053837925*helper_5) - (helper_3 - helper_4)*(0.577350269189626*helper_1 + 0.577350269189626*helper_2 - 1.15470053837925*helper_6));
    double helper_8 = -1.33333333333333*helper_6;
    double helper_9 = 0.666666666666667*helper_6;
    double helper_10 = 0.666666666666667*helper_5;
    double helper_11 = 1.33333333333333*helper_5;
    double helper_12 = 1.15470053837925*helper_7*(helper_1*(-1.33333333333333*helper_1 + 0.666666666666667*helper_2 + helper_9) + helper_2*(0.666666666666667*helper_1 - 1.33333333333333*helper_2 + helper_9) + helper_3*(helper_10 - 1.33333333333333*helper_3 + 0.666666666666667*helper_4) + helper_4*(helper_10 + 0.666666666666667*helper_3 - 1.33333333333333*helper_4) + helper_5*(-helper_11 + 0.666666666666667*helper_3 + 0.666666666666667*helper_4) + helper_6*(0.666666666666667*helper_1 + 0.666666666666667*helper_2 + helper_8));
    result_0(0) = helper_7*(2.66666666666667*helper_1 + helper_12*(helper_4 - helper_5) - 1.33333333333333*helper_2 + helper_8);
    result_0(1) = -helper_7*(helper_11 + helper_12*(helper_2 - helper_6) - 2.66666666666667*helper_3 + 1.33333333333333*helper_4);
}

void triwild::optimization::AMIPS_hessian(const std::array<double, 6>& T, Eigen::Matrix2d& result_0){
    double helper_0[6];
    helper_0[0] = T[0];
    helper_0[1] = T[1];
    helper_0[2] = T[2];
    helper_0[3] = T[3];
    helper_0[4] = T[4];
    helper_0[5] = T[5];
    double helper_1 = helper_0[0];
    double helper_2 = helper_0[2];
    double helper_3 = helper_0[1];
    double helper_4 = helper_0[3];
    double helper_5 = helper_0[5];
    double helper_6 = helper_0[4];
    double helper_7 = (helper_1 - helper_2)*(0.577350269189626*helper_3 + 0.577350269189626*helper_4 - 1.15470053837925*helper_5) - (helper_3 - helper_4)*(0.577350269189626*helper_1 + 0.577350269189626*helper_2 - 1.15470053837925*helper_6);
    double helper_8 = 1.0/helper_7;
    double helper_9 = helper_4 - helper_5;
    double helper_10 = 1.33333333333333*helper_6;
    double helper_11 = -2.66666666666667*helper_1 + helper_10 + 1.33333333333333*helper_2;
    double helper_12 = 2.3094010767585*helper_8;
    double helper_13 = pow(helper_7, -2);
    double helper_14 = 0.666666666666667*helper_6;
    double helper_15 = 0.666666666666667*helper_5;
    double helper_16 = 1.33333333333333*helper_5;
    double helper_17 = helper_1*(-1.33333333333333*helper_1 + helper_14 + 0.666666666666667*helper_2) + helper_2*(0.666666666666667*helper_1 + helper_14 - 1.33333333333333*helper_2) + helper_3*(helper_15 - 1.33333333333333*helper_3 + 0.666666666666667*helper_4) + helper_4*(helper_15 + 0.666666666666667*helper_3 - 1.33333333333333*helper_4) + helper_5*(-helper_16 + 0.666666666666667*helper_3 + 0.666666666666667*helper_4) + helper_6*(0.666666666666667*helper_1 - helper_10 + 0.666666666666667*helper_2);
    double helper_18 = 2.66666666666667*helper_13*helper_17;
    double helper_19 = helper_2 - helper_6;
    double helper_20 = helper_16 - 2.66666666666667*helper_3 + 1.33333333333333*helper_4;
    double helper_21 = helper_13*(-1.15470053837925*helper_11*helper_19 + 2.66666666666667*helper_17*helper_19*helper_8*helper_9 + 1.15470053837925*helper_20*helper_9);
    result_0(0) = helper_8*(helper_11*helper_12*helper_9 - helper_18*pow(helper_9, 2) + 2.66666666666667);
    result_0(1) = helper_21;
    result_0(2) = helper_21;
    result_0(3) = helper_8*(-helper_12*helper_19*helper_20 - helper_18*pow(helper_19, 2) + 2.66666666666667);
}


namespace
{

}


double triwild::optimization::AMIPS_energy(const feature::FeatureElement &feature, const double t, const std::array<double, 6>& T)
{
   const auto pt = feature.eval(t);
   const std::array<double, 6> tmp = {{
        pt.x, pt.y,
        T[2], T[3],
        T[4], T[5]
    }};

    return triwild::optimization::AMIPS_energy(tmp);
}

double triwild::optimization::AMIPS_jacobian(const feature::FeatureElement &feature, const double t, const std::array<double, 6>& T)
{
    const auto pt = feature.eval(t);
    const std::array<double, 6> tmp = {{
        pt.x, pt.y,
        T[2], T[3],
        T[4], T[5]
    }};

    Eigen::Vector2d JE;
    AMIPS_jacobian(tmp, JE);
    const auto gamma_prime = feature.eval_first_derivative(t);

    return JE(0)*gamma_prime.x + JE(1)*gamma_prime.y;
}

double triwild::optimization::AMIPS_hessian(const feature::FeatureElement &feature, const double t, const std::array<double, 6>& T)
{
    const auto pt = feature.eval(t);
    const std::array<double, 6> tmp = {{
        pt.x, pt.y,
        T[2], T[3],
        T[4], T[5]
    }};

    Eigen::Vector2d JE;
    Eigen::Matrix2d HE;
    AMIPS_jacobian(tmp, JE);
    AMIPS_hessian(tmp, HE);

    const auto gamma_prime = feature.eval_first_derivative(t);
    const auto gamma_second = feature.eval_second_derivative(t);

    const double gp_H_gp = (HE(0,0)*gamma_prime.x+HE(1,0)*gamma_prime.y)*gamma_prime.x+(HE(0,1)*gamma_prime.x+HE(1,1)*gamma_prime.y)*gamma_prime.y;
    const double J_gs = JE(0)*gamma_second.x + JE(1)*gamma_second.y;
    return gp_H_gp + J_gs;
}
