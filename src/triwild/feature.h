// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef TRIWILD_FEATURE_H
#define TRIWILD_FEATURE_H

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>
#include <igl/writeOBJ.h>

#include "FeatureElements.h"
#include "TrimeshElements.h"

#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace triwild {
    namespace feature {
        extern std::vector<std::shared_ptr<FeatureElement>> features;
        extern std::vector<std::shared_ptr<FeatureElement>> secondary_features;
        extern double feature_eps;
        extern double feature_eps_2;

        bool init(const std::string& feature_file);
        bool init(json& feature_info);
        void map_feature2mesh(MeshData& mesh);
        bool is_on_segment(const Point_2& p, const Point_2& p1, const Point_2& p2);

        void snap_vertices(MeshData& mesh);

        void merge_inflection(MeshData& mesh);
        void curving(MeshData& mesh, GEO::MeshFacetsAABB &b_tree);
        void add_nodes(MeshData& mesh);
        void subdivide_into_2(MeshData& mesh);
        void subdivide_into_3(MeshData& mesh);
        void get_new_nodes(const std::array<int, 3>& tag_feature_e, const std::array<TriVertex, 3>& vs,
                std::vector<Point_2f>& new_nodes);
        bool is_valid_inversion(const std::array<Point_2f, 3>& ps, const std::vector<Point_2f>& ns);
        void fix_inversion(MeshData& mesh);

        void check_inversion(MeshData& mesh, bool is_output_objs = false);

        void visualize_features(MeshData& mesh);
        Point_2f json2point(const json& a);
        std::array<double, 2> json2d2stdarray(const json& a);
        std::vector<double> json2d2stdvector(const json& a);
        std::vector<double> json1d2stdvector(const json& a);

        ControlVector json2d2ctrlvector(const json& a);
        ControlVector json1d2ctrlvector(const json& a);

        void output_input_features(const std::vector<std::shared_ptr<FeatureElement>>& features,
                Eigen::MatrixXd& V, std::vector<std::array<int, 2>>& edges, const std::string& postfix);
        void output_stats(MeshData& mesh, std::ofstream& f);

        inline void reset(){
            features.clear();
            secondary_features.clear();
        }
    }
}


#endif //TRIWILD_FEATURE_H
