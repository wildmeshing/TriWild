// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef TRIWILD_TRIMESHELEMENTS_H
#define TRIWILD_TRIMESHELEMENTS_H

#include <vector>
#include <array>
#include <unordered_set>
#include <cassert>

#include "Point_2.h"
#include "Point_2f.h"

namespace triwild {
    class TriVertex {
    public:
        Point_2 pos;
        Point_2f posf;
        std::unordered_set<int> conn_tris;
        bool is_rounded = false;
        bool is_on_boundary = false;
        bool is_on_point = false;
        Point_2f input_posf;
        bool is_on_bbox = false;
        double max_scale = 1;
        double scale = 1;

        std::vector<std::array<double, 2>> feature_infos; //feature_id, t
//        bool is_on_feature = false;
//        double t;
        bool is_freezed = false;//feature endpoints or feature points

        double get_t(int feature_id) const {
            assert(feature_id >= 0);
            assert(feature_infos.size() > 0);
            if (feature_infos.size() == 1)
                return feature_infos[0][1];
            else {
                for (auto &info:feature_infos)
                    if (info[0] == feature_id)
                        return info[1];
            }
            assert(false);//should not reach here
            return std::numeric_limits<double>::max();
        }
        bool has_feature(int feature_id) const {
            for (auto &info: feature_infos)
                if (info[0] == feature_id)
                    return true;
            return false;
        }
    };

    struct MeshData {
        std::vector<TriVertex> tri_vertices;
        std::vector<bool> v_is_removed;
        std::vector<double> v_scalars;

        std::vector<std::array<int, 3>> tris;
        std::vector<std::vector<int>> tri_nodes;
        std::vector<std::array<bool, 3>> is_boundary_es;
        std::vector<std::array<bool, 3>> is_bbox_es;
        std::vector<std::array<int, 3>> tag_feature_es;
        std::vector<std::array<int, 3>> tag_secondary_feature_es;
        std::vector<double> t_quality;
        std::vector<bool> t_is_removed;
        std::vector<double> t_scalars;

        std::vector<Point_2f> nodes;
        std::vector<bool> n_is_removed;

        bool is_curved = false;

        double ideal_edge_length;
        double epsilon;
        double original_epsilon;
        double epsilon_2;
        double dd;
        int max_its = 80;
        double stop_energy = 10;

        double feature_epsilon;
        double min_scalar = 0;

        int v_empty_slot_start;
        int t_empty_slot_start;
//        bool is_preserving_feature = false;
        bool is_limit_length = true;
        bool is_edge_length_achieved = false;
    };
}

#endif //TRIWILD_TRIMESHELEMENTS_H
