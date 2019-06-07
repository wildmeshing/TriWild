// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <Eigen/Dense>
#include<vector>
#include<array>

namespace triwild{
namespace autogen{
class AutoDetChecker{
public:
std::array<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 15, 15>,2> L2B;
std::array<std::array<std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2, 0, 15, 2>>, 6>, 2> loc_nodes;

static const int MAX_LEVEL=6;

static const AutoDetChecker &instance();

int get_index(const int parent, const int q) const;
private:
AutoDetChecker();
};
}
}
