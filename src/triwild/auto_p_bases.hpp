// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include <Eigen/Dense>

namespace triwild {
	namespace autogen {
		template<typename Matrix>
		void p_nodes_2d(const int p, Matrix &val);

		template<typename InMatrix, typename OutMatrix>
		void p_basis_value_2d(const int p, const int local_index, const InMatrix &uv, OutMatrix &val);

		template<typename Matrix>
		void p_grad_basis_value_2d(const int p, const int local_index, const Matrix &uv, Matrix &val);

		static const int MAX_P_BASES = 4;
	}
}
