// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "edge_swapping.h"
#include "optimization.h"
#include "feature.h"

void triwild::optimization::edge_swapping(MeshData &mesh)
{
    ////init
    std::vector<std::array<int, 2>> edges;
    get_all_edges(mesh, edges);

    std::priority_queue<ElementInQueue_l, std::vector<ElementInQueue_l>, cmp_l> queue_l;
    for (const auto &e : edges)
    {
        if (is_bbox_edge(mesh, e[0], e[1]) || is_boundary_edge(mesh, e[0], e[1]) || is_feature_edge(mesh, e[0], e[1]))
            continue;
        queue_l.push(ElementInQueue_l(e, edge_length_2(mesh, e[0], e[1])));
    }

    ////swap
    int cnt = 0;
    while (!queue_l.empty())
    {
        std::array<int, 2> v_ids = queue_l.top().v_ids;
        queue_l.pop();

        if (!is_valid_edge(mesh, v_ids[0], v_ids[1]))
            continue;
        if (!queue_l.empty() && v_ids == queue_l.top().v_ids)
            continue;

        std::vector<std::array<int, 2>> new_es;
        if (swap_an_edge(mesh, v_ids[0], v_ids[1], new_es))
        {
            cnt++;
            for (auto &new_e : new_es)
            {
                if (is_bbox_edge(mesh, new_e[0], new_e[1]) || is_boundary_edge(mesh, new_e[0], new_e[1]) || is_feature_edge(mesh, new_e[0], new_e[1]))
                    continue;
                queue_l.push(ElementInQueue_l(new_e, edge_length_2(mesh, new_e[0], new_e[1])));
            }
        }
    }
    cout << "success " << cnt << endl;
}

bool triwild::optimization::swap_an_edge(MeshData &mesh, const int v1_id, const int v2_id,
                                         std::vector<std::array<int, 2>> &new_es)
{
    std::vector<int> n_t_ids = set_intersection(mesh.tri_vertices[v1_id].conn_tris, mesh.tri_vertices[v2_id].conn_tris);
    if (n_t_ids.size() == 1)
        return false;

    std::vector<std::array<int, 3>> new_ts;
    std::vector<int> opp_v_ids;
    for (int t_id : n_t_ids)
    {
        new_ts.push_back(mesh.tris[t_id]);
        for (int j = 0; j < 3; j++)
        {
            if (mesh.tris[t_id][j] != v1_id && mesh.tris[t_id][j] != v2_id)
                opp_v_ids.push_back(mesh.tris[t_id][j]);
        }
    }
    for (int j = 0; j < 3; j++)
    {
        if (new_ts[0][j] == v2_id)
            new_ts[0][j] = opp_v_ids[1];
        if (new_ts[1][j] == v1_id)
            new_ts[1][j] = opp_v_ids[0];
    }

    ////check inversion/quality
    for (auto &t : new_ts)
    {
        if (!is_valid_inversion(mesh.tri_vertices[t[0]].posf, mesh.tri_vertices[t[1]].posf, mesh.tri_vertices[t[2]].posf,
                                mesh.tri_vertices[t[0]].pos, mesh.tri_vertices[t[1]].pos, mesh.tri_vertices[t[2]].pos))
            return false;
    }
    std::vector<double> new_qs;
    double old_max_energy = mesh.t_quality[n_t_ids[0]] > mesh.t_quality[n_t_ids[1]] ? mesh.t_quality[n_t_ids[0]]
                                                                                    : mesh.t_quality[n_t_ids[1]];
    for (auto &t : new_ts)
    {
        double new_energy;
        if (!is_valid_quality(mesh.tri_vertices[t[0]].posf, mesh.tri_vertices[t[1]].posf, mesh.tri_vertices[t[2]].posf,
                              old_max_energy, new_energy))
            return false;
        new_qs.push_back(new_energy);
    }

    assert(n_t_ids.size() == 2);
    if (mesh.is_curved && (!mesh.tri_nodes[n_t_ids[0]].empty() || !mesh.tri_nodes[n_t_ids[1]].empty()))
    {
        std::vector<std::array<int, 2>> opp_tag_feature_e;
        for (int t_id : n_t_ids)
        {
            opp_tag_feature_e.emplace_back();
            auto &tag_feature_e = opp_tag_feature_e.back();
            for (int j = 0; j < 3; j++)
            {
                if (mesh.tris[t_id][j] == v1_id)
                    tag_feature_e[0] = mesh.tag_feature_es[t_id][j];
                if (mesh.tris[t_id][j] == v2_id)
                    tag_feature_e[1] = mesh.tag_feature_es[t_id][j];
            }
        }

        std::vector<std::vector<Point_2f>> vec_new_nodes;
        for (int i = 0; i < new_ts.size(); i++)
        {
            vec_new_nodes.emplace_back();
            auto &new_nodes = vec_new_nodes.back();

            std::array<int, 3> tag_feature_e = {{-1, -1, -1}};
            for (int j = 0; j < 3; j++)
            {
                if (new_ts[i][j] == v1_id || new_ts[i][j] == v2_id)
                    tag_feature_e[j] = -1;
                else if (new_ts[i][j] == opp_v_ids[0])
                    tag_feature_e[j] = opp_tag_feature_e[1][(i + 1) % 2];
                else if (new_ts[i][j] == opp_v_ids[1])
                    tag_feature_e[j] = opp_tag_feature_e[0][(i + 1) % 2];
            }
            feature::get_new_nodes(tag_feature_e,
                                   {{mesh.tri_vertices[new_ts[i][0]], mesh.tri_vertices[new_ts[i][1]], mesh.tri_vertices[new_ts[i][2]]}},
                                   new_nodes);
            if (new_nodes.empty())
                continue;
            if (!feature::is_valid_inversion(
                    {{mesh.tri_vertices[new_ts[i][0]].posf, mesh.tri_vertices[new_ts[i][1]].posf, mesh.tri_vertices[new_ts[i][2]].posf}},
                    new_nodes))
            {
                return false;
            }
        }

        ////real update
        int cnt = 0;
        for (int t_id : n_t_ids)
        {
            auto &nodes = vec_new_nodes[cnt++];
            if (nodes.size() == mesh.tri_nodes[t_id].size())
            {
                for (int i = 0; i < mesh.tri_nodes[t_id].size(); i++)
                {
                    mesh.nodes[mesh.tri_nodes[t_id][i]] = nodes[i];
                }
            }
            else
            {
                for (int i = 0; i < mesh.tri_nodes[t_id].size(); i++)
                    mesh.n_is_removed[mesh.tri_nodes[t_id][i]] = true;
                mesh.tri_nodes[t_id].clear();
                for (auto &n : nodes)
                {
                    mesh.nodes.push_back(n);
                    mesh.n_is_removed.push_back(false);
                    mesh.tri_nodes[t_id].push_back(mesh.nodes.size() - 1);
                }
            }
        }
    }

    ////real update
    for (int j = 0; j < 3; j++)
    {
        if (mesh.tris[n_t_ids[0]][j] == v1_id)
        {
            for (int k = 0; k < 3; k++)
                if (new_ts[1][k] == opp_v_ids[1])
                {
                    mesh.is_boundary_es[n_t_ids[1]][k] = mesh.is_boundary_es[n_t_ids[0]][j];
                    mesh.is_bbox_es[n_t_ids[1]][k] = mesh.is_bbox_es[n_t_ids[0]][j];
                    mesh.tag_feature_es[n_t_ids[1]][k] = mesh.tag_feature_es[n_t_ids[0]][j];
                    mesh.tag_secondary_feature_es[n_t_ids[1]][k] = mesh.tag_secondary_feature_es[n_t_ids[0]][j];
                }
        }
        if (mesh.tris[n_t_ids[1]][j] == v2_id)
        {
            for (int k = 0; k < 3; k++)
                if (new_ts[0][k] == opp_v_ids[0])
                {
                    mesh.is_boundary_es[n_t_ids[0]][k] = mesh.is_boundary_es[n_t_ids[1]][j];
                    mesh.is_bbox_es[n_t_ids[0]][k] = mesh.is_bbox_es[n_t_ids[1]][j];
                    mesh.tag_feature_es[n_t_ids[0]][k] = mesh.tag_feature_es[n_t_ids[1]][j];
                    mesh.tag_secondary_feature_es[n_t_ids[0]][k] = mesh.tag_secondary_feature_es[n_t_ids[1]][j];
                }
        }
    }
    mesh.tris[n_t_ids[0]] = new_ts[0];
    mesh.tris[n_t_ids[1]] = new_ts[1];
    for (int j = 0; j < 3; j++)
    {
        if (mesh.tris[n_t_ids[0]][j] == v1_id)
        {
            mesh.is_boundary_es[n_t_ids[0]][j] = false;
            mesh.is_bbox_es[n_t_ids[0]][j] = false;
            mesh.tag_feature_es[n_t_ids[0]][j] = -1;
            mesh.tag_secondary_feature_es[n_t_ids[0]][j] = -1;
        }
        if (mesh.tris[n_t_ids[1]][j] == v2_id)
        {
            mesh.is_boundary_es[n_t_ids[1]][j] = false;
            mesh.is_bbox_es[n_t_ids[1]][j] = false;
            mesh.tag_feature_es[n_t_ids[1]][j] = -1;
            mesh.tag_secondary_feature_es[n_t_ids[1]][j] = -1;
        }
    }

    mesh.t_quality[n_t_ids[0]] = new_qs[0];
    mesh.t_quality[n_t_ids[1]] = new_qs[1];

    mesh.tri_vertices[v1_id].conn_tris.erase(n_t_ids[1]);
    mesh.tri_vertices[v2_id].conn_tris.erase(n_t_ids[0]);
    mesh.tri_vertices[opp_v_ids[0]].conn_tris.insert(n_t_ids[1]);
    mesh.tri_vertices[opp_v_ids[1]].conn_tris.insert(n_t_ids[0]);

    //return new_es
    new_es.push_back(std::array<int, 2>({opp_v_ids[0], opp_v_ids[1]}));
    for (int i = 0; i < 2; i++)
    {
        new_es.push_back(std::array<int, 2>({opp_v_ids[i], v1_id}));
        new_es.push_back(std::array<int, 2>({opp_v_ids[i], v2_id}));
    }

    return true;
}