// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef TRIWILD_EDGE_SWAPPING_H
#define TRIWILD_EDGE_SWAPPING_H

#include "TrimeshElements.h"

namespace triwild {
    namespace optimization {
        void edge_swapping(MeshData& mesh);
        bool swap_an_edge(MeshData& mesh, const int v1_id, const int v2_id, std::vector<std::array<int, 2>>& new_es);
    }
}


#endif //TRIWILD_EDGE_SWAPPING_H
