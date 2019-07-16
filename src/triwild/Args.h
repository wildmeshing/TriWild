// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef TRIWILD_ARGS_H
#define TRIWILD_ARGS_H

#include <iostream>
#include <fstream>
#include <array>
#include <queue>
#include <string>
#include <limits>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <iterator>
#include <algorithm>

#include <Eigen/Core>
#include <Eigen/Dense>

namespace triwild {
    using std::cout;
    using std::cin;
    using std::endl;

    struct Arguments {
        //arguments
        double epsilon = -1;
        double feature_epsilon = 1e-3;
        double target_edge_len = -1;
        std::string input = "";
        std::string output = "";
        std::string feature_input = "";
//        std::string quality_output = "";
        std::string postfix = "";
        std::string log_file = "";
        int stage = 1;

        double stop_quality = -1;
        int max_its = 80;
        int i_dd = -1;
        float edge_length_r = 1/20.0;
        bool mute_log = false;
        double flat_feature_angle = 10;
        double min_edge_length;

        //variables
        std::string first_triangulation = "delaunay";
        Eigen::Vector2d box_min, box_max;
        double diagonal_len = -1;
        double envelope;
        double sampling_density;
        double v_bad_const = 0.5;
        double v_good_const = 1.5;
        double t_energy_threshold = 10;
        const double MAX_ENERGY = 1e50;
        bool is_preserving_feature = false;

        bool enable_debug_mesh = false;
        bool output_linear = false;

        inline void reset(){
            epsilon = -1;
            feature_epsilon = 1e-3;
            target_edge_len = -1;
            input = "";
            output = "";
            feature_input = "";
            postfix = "";
            log_file = "";
            stage = 1;

            stop_quality = -1;
            max_its = 80;
            i_dd = -1;
            edge_length_r = 1/20.0;
            mute_log = false;
            flat_feature_angle = 10;

            //variables
            first_triangulation = "delaunay";
            diagonal_len = -1;
            v_bad_const = 0.5;
            v_good_const = 1.5;
            t_energy_threshold = 10;
            is_preserving_feature = false;

            enable_debug_mesh = false;
            output_linear = false;
        }
    };

    extern Arguments args;
}

#endif //TRIWILD_ARGS_H
