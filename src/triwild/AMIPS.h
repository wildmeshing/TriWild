// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "FeatureElements.h"
#include <Eigen/Dense>
#include <array>

namespace triwild {
	namespace optimization{
		double AMIPS_energy(const std::array<double, 6>& T);
		void AMIPS_jacobian(const std::array<double, 6>& T, Eigen::Vector2d& J);
		void AMIPS_hessian(const std::array<double, 6>& T, Eigen::Matrix2d& H);

		double AMIPS_energy(const feature::FeatureElement &feature, const double t, const std::array<double, 6>& T);
		double AMIPS_jacobian(const feature::FeatureElement &feature, const double t, const std::array<double, 6>& T);
		double AMIPS_hessian(const feature::FeatureElement &feature, const double t, const std::array<double, 6>& T);
	}
}