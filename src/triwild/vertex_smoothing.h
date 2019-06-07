// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef TRIWILD_VERTEX_SMOOTHING_H
#define TRIWILD_VERTEX_SMOOTHING_H

#include "TrimeshElements.h"

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>

namespace triwild {
    namespace optimization {
        void vertex_smoothing(MeshData& mesh, GEO::MeshFacetsAABB &b_tree);
        bool smooth_a_vertex(MeshData& mesh, const int v_id, Point_2f& pf, double& t);
    }
}

#endif //TRIWILD_VERTEX_SMOOTHING_H
