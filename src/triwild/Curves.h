// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include "TrimeshElements.h"

#include <Eigen/Dense>
#include <vector>
#include <array>

namespace triwild {
	typedef Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 8, 1> ControlVector;

	class Bezier
	{
	public:
		static std::array<double, 2> interpolate(const ControlVector &ctrls, const double t);

		static std::array<double, 2> first_derivative(const ControlVector &ctrls, const double t);
		static std::array<double, 2> second_derivative(const ControlVector &ctrls, const double t);

		static void resample(const ControlVector &ctrls, const double t0, const double t1, ControlVector &new_ctrl);

		static double inverse_interpolation(const ControlVector &ctrls, const std::array<double, 2> &p, const double start = 0);

		static double point_curve_distance(const ControlVector &ctrls, const std::array<double, 2> &p, const double t0, const double t1);

	private:
		static std::array<double, 2> interpolate_1(const ControlVector &ctrls, const double t);
		static std::array<double, 2> first_derivative_1(const ControlVector &ctrls, const double t);
		static std::array<double, 2> second_derivative_1(const ControlVector &ctrls, const double t);
		static void resample_1(const ControlVector &ctrls, const double t0, const double t1, ControlVector &new_ctrl);

		static std::array<double, 2> interpolate_2(const ControlVector &ctrls, const double t);
		static std::array<double, 2> first_derivative_2(const ControlVector &ctrls, const double t);
		static std::array<double, 2> second_derivative_2(const ControlVector &ctrls, const double t);
		static void resample_2(const ControlVector &ctrls, const double t0, const double t1, ControlVector &new_ctrl);

		static std::array<double, 2> interpolate_3(const ControlVector &ctrls, const double t);
		static std::array<double, 2> first_derivative_3(const ControlVector &ctrls, const double t);
		static std::array<double, 2> second_derivative_3(const ControlVector &ctrls, const double t);
		static void resample_3(const ControlVector &ctrls, const double t0, const double t1, ControlVector &new_ctrl);

		static double point_curve_distance_3(const ControlVector &ctrls, const std::array<double, 2> &p, const double t0, const double t1);
		static double point_curve_distance_2(const ControlVector &ctrls, const std::array<double, 2> &p, const double t0, const double t1);

	};

	class RationalBezier
	{
	public:
		static std::array<double, 2> interpolate(const ControlVector &ctrls, const ControlVector &weights, const double t);

		static std::array<double, 2> first_derivative(const ControlVector &ctrls, const ControlVector &weights, const double t);
		static std::array<double, 2> second_derivative(const ControlVector &ctrls, const ControlVector &weights, const double t);

		static double inverse_interpolation(const ControlVector &ctrls, const ControlVector &weights, const std::array<double, 2> &p, const double start = 0);
	};

	class DeterminantChecker
	{
	public:
		static bool is_positive(const std::array<Point_2f, 3>& vertices, const std::vector<Point_2f>& nodes);
	};
}
