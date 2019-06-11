// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
//
// Created by Yixin Hu on 2019-06-11.
//

#ifndef TRIWILD_REFERENCE_TRIANGLE_H
#define TRIWILD_REFERENCE_TRIANGLE_H

#include "Args.h"
#include <Eigen/Dense>

namespace triwild {
    Eigen::MatrixXd get_reference_triangle_vertices();
    Eigen::MatrixXi get_reference_triangle_faces();
}



#endif //TRIWILD_REFERENCE_TRIANGLE_H
