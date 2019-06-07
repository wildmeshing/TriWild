// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef TRIWILD_TRIANGULATION_H
#define TRIWILD_TRIANGULATION_H

//#include <eigen/Eigen/Eigen>
#include "Point_2.h"
#include "Point_2f.h"
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>
#include "TrimeshElements.h"
#include "Args.h"

namespace triwild {
    namespace triangulation {
        bool load_input(const std::string& input, Eigen::MatrixXd& V, std::vector<std::array<int, 2>>& edges);
        void preprocessing(Eigen::MatrixXd& V, std::vector<std::array<int, 2>>& edges, GEO::Mesh& b_mesh);
        void simplify_input(Eigen::MatrixXd& V, std::vector<std::array<int, 2>>& edges, GEO::MeshFacetsAABB &b_tree);
        void BSP_subdivision(const Eigen::MatrixXd& V, const std::vector<std::array<int, 2>>& edges,
                MeshData& mesh, std::vector<std::vector<int>>& tag_boundary_es);
        bool segment_intersection(MeshData& mesh, int v_id, int v2_id, int p1_id, int p2_id,
                bool& is_cross_p1, bool& is_cross_p2, TriVertex& intersection_v, bool is_check_bbox = true);
    }
}


#endif //TRIWILD_TRIANGULATION_H
