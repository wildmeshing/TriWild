// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef TRIWILD_FEATURE_PREPROCESSING_H
#define TRIWILD_FEATURE_PREPROCESSING_H

#include "FeatureElements.h"
#include "TrimeshElements.h"

namespace triwild {
    namespace feature {
        void preprocessing(Eigen::MatrixXd& V, std::vector<std::array<int, 2>>& edges);

        void simplify();

        void get_inflections(std::vector<std::vector<double>>& inflections);
        void cut_reflex(std::vector<std::vector<double>>& inflections);

        void sample_features(std::vector<std::vector<Point_2f>>& samples, std::vector<std::vector<double>>& ts);

        void mu_separation(std::vector<std::vector<double>>& inflections, std::vector<std::vector<Point_2f>>& samples,
                           std::vector<std::vector<double>>& ts);

        void remove_high_curvature(std::vector<std::vector<double>>& inflections, std::vector<std::vector<Point_2f>>& samples,
                                   std::vector<std::vector<double>>& ts);

        void remove_short_features(std::vector<std::vector<double>>& inflections, std::vector<std::vector<Point_2f>>& samples,
                                   std::vector<std::vector<double>>& ts);

        void cut_inflections(std::vector<std::vector<double>>& inflections, std::vector<std::vector<Point_2f>>& samples,
                             std::vector<std::vector<double>>& ts);

        void gen_segments(Eigen::MatrixXd& V, std::vector<std::array<int, 2>>& edges);

        int push_back_new_feature(std::shared_ptr<FeatureElement> old_feature);
        int push_back_new_secondary_feature(std::shared_ptr<FeatureElement> old_feature);

        void output_features(const std::string& name);
        void check(std::vector<std::vector<double>>& inflections, std::vector<std::vector<Point_2f>>& samples,
                   std::vector<std::vector<double>>& ts);
    }
}

#endif //TRIWILD_FEATURE_PREPROCESSING_H
