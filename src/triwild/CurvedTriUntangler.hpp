// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include "TrimeshElements.h"

#include <vector>
#include <array>

namespace triwild {
	class CurvedTriUntangler
	{
	public:
		static bool untangle(const std::array<Point_2f, 3>& vertices, const int edge,  const Point_2f &fixed1, const Point_2f &fixed2, const Point_2f &center, std::vector<Point_2f> &new_nodes);

		static bool untangle_center(const std::array<Point_2f, 3>& vertices, const std::vector<Point_2f> &fixed, const Point_2f &center, Point_2f &new_center);

		static double ls_fit(const std::array<Point_2f, 3>& vertices, const std::vector<Point_2f> &nodes, const int edge, std::vector<Point_2f> &new_nodes, const bool all = true);
	};
}
