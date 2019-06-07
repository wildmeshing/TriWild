// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

#include <vector>

namespace triwild
{
    class Logger
    {
    public:
        static Logger &instance()
        {
            static Logger tmp;

            return tmp;
        }

        long const_calls;
        long unicessary_const_calls;




        long ls_fitted;
        long ls_fitted_fail;
        long unsnapped_vertices;
        long invalid_element;
        std::vector<double> ls_fitting_distances;

    private:
        Logger()
        :const_calls(0)
        { }

    };
}
