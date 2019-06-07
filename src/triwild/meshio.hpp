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

namespace triwild
{
    void write_OBJ( Eigen::MatrixXd V, Eigen::MatrixXi F, const std::string & path);

    void write_msh(const MeshData& mesh, const std::string &path, const bool export_edge_tag = true);
    void load_msh(const std::string &path, MeshData& mesh);
    void write_msh_DiffusionCurve(MeshData& mesh, const std::string &path);
    bool unordered_set_intersection_own(const std::unordered_set<int> &A, const std::unordered_set<int> &B, std::vector<int> &C);

    void export_eps(MeshData& mesh,
    	const double line_width, const std::string &col, const double point_size, const std::string &point_col,
        const double feature_line_width, const std::string &feature_col, const double feature_point_size, const std::string &feature_point_col,
    	const double secondary_feature_line_width, const std::string &secondary_feature_col, const double secondary_feature_point_size, const std::string &secondary_feature_point_col,
    	const bool draw_points, const std::string &path, const double t=1);
}