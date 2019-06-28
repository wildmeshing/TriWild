// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "edge_splitting.h"
#include "optimization.h"
#include "feature.h"

void triwild::optimization::edge_splitting(MeshData &mesh)
{
    //update v_empty_slot_start and t_empty_slot_start
    mesh.v_empty_slot_start =
        std::find(mesh.v_is_removed.begin(), mesh.v_is_removed.end(), true) - mesh.v_is_removed.begin();
    mesh.t_empty_slot_start =
        std::find(mesh.t_is_removed.begin(), mesh.t_is_removed.end(), true) - mesh.t_is_removed.begin();

    ////init
    std::vector<std::array<int, 2>> edges;
    std::priority_queue<ElementInQueue_l, std::vector<ElementInQueue_l>, cmp_l> queue_l;
    get_all_edges(mesh, edges);

    for (const auto &e : edges)
    {
        double l = std::sqrt(edge_length_2(mesh, e[0], e[1]));
        if (l >
            (mesh.tri_vertices[e[0]].scale + mesh.tri_vertices[e[1]].scale) / 2 * mesh.ideal_edge_length *
                (4.0 / 3.0))
            queue_l.push(ElementInQueue_l(e, l));
    }

    ////split
    int v_slot_size = std::count(mesh.v_is_removed.begin(), mesh.v_is_removed.end(), true);
    int t_slot_size = std::count(mesh.t_is_removed.begin(), mesh.t_is_removed.end(), true);
    if (v_slot_size < queue_l.size() * 2)
        mesh.tri_vertices.reserve(queue_l.size() * 2 - v_slot_size);
    if (t_slot_size < queue_l.size() * 3 * 2)
        mesh.tris.reserve(queue_l.size() * 3 * 2 - t_slot_size + 1);

    if (mesh.is_edge_length_achieved || queue_l.top().weight >= mesh.ideal_edge_length)
        mesh.is_edge_length_achieved = true;

    bool is_push_new_edges = true;
    while (!queue_l.empty())
    {
        const ElementInQueue_l &ele = queue_l.top();
        ////split an edge
        std::vector<std::array<int, 2>> new_es;
        //        split_an_edge(mesh, ele.v_ids[0], ele.v_ids[1], new_es);
        if (!split_an_edge(mesh, ele.v_ids[0], ele.v_ids[1], new_es))
        {
            is_push_new_edges = false;
        }
        queue_l.pop();

        if (is_push_new_edges)
        {
            //        if(!args.is_preserving_feature){
            for (const auto e : new_es)
            {
                double l = std::sqrt(edge_length_2(mesh, e[0], e[1]));
                if (l > (mesh.tri_vertices[e[0]].scale + mesh.tri_vertices[e[1]].scale) / 2 * mesh.ideal_edge_length *
                            (4.0 / 3.0))
                    queue_l.push(ElementInQueue_l(e, l));
            }
        }
    }

    ////update qualities
    for (size_t i = 0; i < mesh.tris.size(); i++)
    {
        if (mesh.t_is_removed[i])
            continue;
        if (mesh.t_quality[i] < 0)
            mesh.t_quality[i] = tri_energy(mesh.tri_vertices[mesh.tris[i][0]].posf,
                                           mesh.tri_vertices[mesh.tris[i][1]].posf,
                                           mesh.tri_vertices[mesh.tris[i][2]].posf);
    }
}

bool triwild::optimization::split_an_edge(MeshData &mesh, const int v1_id, const int v2_id,
                                          std::vector<std::array<int, 2>> &new_es)
{
    //    if(mesh.is_curved)
    //        cout<<v1_id<<" "<<v2_id<<endl;
    ///compute
    int feature_id = get_feature_edge_tag(mesh, v1_id, v2_id);
    double t = 0, t1, t2;
    if (feature_id >= 0)
    {
        t1 = mesh.tri_vertices[v1_id].get_t(feature_id);
        t2 = mesh.tri_vertices[v2_id].get_t(feature_id);
        if (t1 > t2)
            std::swap(t1, t2);
        t = (t1 + t2) / 2;
    }

    std::vector<int> n_t_ids = set_intersection(mesh.tri_vertices[v1_id].conn_tris,
                                                mesh.tri_vertices[v2_id].conn_tris);

    Point_2f pf;
    Point_2 p;
    bool is_rounded = true;

    auto split_normal_edge = [&]() {
        pf = mid_point(mesh.tri_vertices[v1_id].posf, mesh.tri_vertices[v2_id].posf);
        p = Point_2(pf[0], pf[1]);
        for (int t_id : n_t_ids)
        {
            for (int j = 0; j < 3; j++)
            {
                if ((mesh.tris[t_id][j] == v1_id || mesh.tris[t_id][j] == v2_id) && !is_valid_inversion(pf, mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].posf,
                                                                                                        mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].posf,
                                                                                                        p, mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].pos,
                                                                                                        mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].pos))
                {
                    is_rounded = false;
                    break;
                }
            }
            if (!is_rounded)
                break;
        }

        if (!is_rounded)
        {
            p = mid_point(mesh.tri_vertices[v1_id].pos, mesh.tri_vertices[v2_id].pos);
            pf[0] = p[0].to_double();
            pf[1] = p[1].to_double();
        }
    };

    if (feature_id >= 0 && feature::features[feature_id]->type != "Line")
    {
        pf = mid_point(mesh.tri_vertices[v1_id].posf, mesh.tri_vertices[v2_id].posf);
        t = feature::features[feature_id]->inv_eval(pf, t, t1, t2);
        pf = feature::features[feature_id]->eval(t);
        p = Point_2(pf[0], pf[1]);

        bool is_valid = true;
        for (int i = 0; i < n_t_ids.size(); i++)
        {
            int t_id = n_t_ids[i];
            for (int j = 0; j < 3; j++)
            {
                if ((mesh.tris[t_id][j] == v1_id || mesh.tris[t_id][j] == v2_id) && !is_valid_inversion(pf, mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].posf,
                                                                                                        mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].posf,
                                                                                                        p, mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].pos,
                                                                                                        mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].pos))
                {
                    is_valid = false;
                    break;
                }
            }
            if (!is_valid)
                break;
        }

        if (!is_valid)
        {
            //split
            split_normal_edge();

            //check pocket
            std::vector<bool> is_visited(mesh.tris.size(), false);
            std::vector<bool> is_inside(mesh.tri_vertices.size(), false); //on or inside
            std::queue<int> v_queue;
            v_queue.push(n_t_ids.front());
            is_visited[n_t_ids.front()] = true;
            is_inside[v1_id] = true;
            is_inside[v2_id] = true;
            while (!v_queue.empty())
            {
                int t_id = v_queue.front();
                v_queue.pop();

                for (int j = 0; j < 3; j++)
                {
                    if (!is_inside[mesh.tris[t_id][j]] && !is_inside[mesh.tris[t_id][(j + 1) % 3]])
                        continue;
                    auto tmp = set_intersection(mesh.tri_vertices[mesh.tris[t_id][j]].conn_tris,
                                                mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].conn_tris);
                    for (int tmp_t_id : tmp)
                    {
                        if (tmp_t_id == t_id || is_visited[tmp_t_id])
                            continue;
                        int k = get_opp_r_v_id(mesh, tmp_t_id, mesh.tris[t_id][j], mesh.tris[t_id][(j + 1) % 3]);
                        if (is_inside_pocket(mesh, mesh.tri_vertices[v1_id].posf, mesh.tri_vertices[v2_id].posf,
                                             feature_id, t1, t2, mesh.tri_vertices[mesh.tris[t_id][k]].posf))
                            is_inside[mesh.tris[t_id][k]] = true;
                    }
                }
            }

            for (int i = 0; i < is_inside.size(); i++)
            {
                if (!is_inside[i])
                    continue;
                if (!mesh.tri_vertices[i].feature_infos.empty())
                    continue;
                if (!mesh.tri_vertices[i].is_on_boundary)
                    continue;
                mesh.tri_vertices[i].is_on_boundary = false;
                for (int n_t_id : mesh.tri_vertices[i].conn_tris)
                {
                    for (int j = 0; j < 3; j++)
                        if (mesh.tris[n_t_id][j] != i)
                            mesh.is_boundary_es[n_t_id][j] = false;
                }
            }

            cout << ">>>>UNTANGLED<<<<" << endl;
        }

        //        for (int i = 0; i < n_t_ids.size(); i++) {
        //            int t_id = n_t_ids[i];
        //            bool is_valid = true;
        //            for (int j = 0; j < 3; j++) {
        //                if ((mesh.tris[t_id][j] == v1_id || mesh.tris[t_id][j] == v2_id)
        //                    && !is_valid_inversion(pf, mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].posf,
        //                                           mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].posf,
        //                                           p, mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].pos,
        //                                           mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].pos)) {
        //                    is_valid = false;
        //                    break;
        //                }
        //            }
        //            if (is_valid)
        //                continue;
        //
        //            //change the feature tags
        //            assert(n_t_ids.size() == 2);
        //            int other_t_id = n_t_ids[(i + 1) % n_t_ids.size()];
        //            int j_t_id;
        //            int j_other_t_id;
        //            for (int j = 0; j < 3; j++) {
        //                if (mesh.tris[t_id][j] != v1_id && mesh.tris[t_id][j] != v2_id) {
        //                    j_t_id = j;
        //                }
        //                if (mesh.tris[other_t_id][j] != v1_id && mesh.tris[other_t_id][j] != v2_id) {
        //                    j_other_t_id = j;
        //                }
        //            }
        //
        //            cout<<"v1_id "<<v1_id<<endl;
        //            cout<<"v2_id "<<v2_id<<endl;
        //            cout<<"j_t_id "<<j_t_id<<endl;
        //            cout<<"j_other_t_id "<<j_other_t_id<<endl;
        //
        //            cout<<"t_id "<<t_id<<" "<<mesh.tris[t_id][0]<<" "<<mesh.tris[t_id][1]<<" "<<mesh.tris[t_id][2]<<endl;
        //            cout<<mesh.tag_feature_es[t_id][0]<<" "<<mesh.tag_feature_es[t_id][1]<<" "<<mesh.tag_feature_es[t_id][2]<<endl;
        //            cout<<"other_t_id"<<other_t_id<<" "<<mesh.tris[other_t_id][0]<<" "<<mesh.tris[other_t_id][1]<<" "<<mesh.tris[other_t_id][2]<<endl;
        //            cout<<mesh.tag_feature_es[other_t_id][0]<<" "<<mesh.tag_feature_es[other_t_id][1]<<" "<<mesh.tag_feature_es[other_t_id][2]<<endl;
        //
        //            mesh.tag_feature_es[t_id][j_t_id] = -1;
        //            mesh.tag_feature_es[other_t_id][j_other_t_id] = -1;
        //
        //            mesh.tag_feature_es[t_id][(j_t_id + 1) % 3] = feature_id;
        //            mesh.tag_feature_es[t_id][(j_t_id + 2) % 3] = feature_id;
        //
        //            cout<<"t_id "<<t_id<<": "<<mesh.tris[t_id][0]<<" "<<mesh.tris[t_id][1]<<" "<<mesh.tris[t_id][2]<<endl;
        //            cout<<mesh.tag_feature_es[t_id][0]<<" "<<mesh.tag_feature_es[t_id][1]<<" "<<mesh.tag_feature_es[t_id][2]<<endl;
        //            cout<<"other_t_id"<<other_t_id<<": "<<mesh.tris[other_t_id][0]<<" "<<mesh.tris[other_t_id][1]<<" "<<mesh.tris[other_t_id][2]<<endl;
        //            cout<<mesh.tag_feature_es[other_t_id][0]<<" "<<mesh.tag_feature_es[other_t_id][1]<<" "<<mesh.tag_feature_es[other_t_id][2]<<endl;
        //
        //            int v_id = mesh.tris[t_id][j_t_id];
        //            auto tmp = set_intersection(mesh.tri_vertices[v1_id].conn_tris, mesh.tri_vertices[v_id].conn_tris);
        //            if (tmp.size() == 2) {
        //                int other_t1_id = tmp[0] == t_id ? tmp[1] : tmp[0];
        //                cout<<"other_t1_id "<<other_t1_id<<": "<<mesh.tris[other_t1_id][0]<<" "<<mesh.tris[other_t1_id][1]<<" "<<mesh.tris[other_t1_id][2]<<endl;
        //                cout<<mesh.tag_feature_es[other_t1_id][0]<<" "<<mesh.tag_feature_es[other_t1_id][1]<<" "<<mesh.tag_feature_es[other_t1_id][2]<<endl;
        //                for (int j = 0; j < 3; j++) {
        //                    if (mesh.tris[other_t1_id][j] != v1_id && mesh.tris[other_t1_id][j] != v_id) {
        //                        mesh.tag_feature_es[other_t1_id][j] = feature_id;
        //                        break;
        //                    }
        //                }
        //                cout<<"other_t1_id "<<other_t1_id<<": "<<mesh.tris[other_t1_id][0]<<" "<<mesh.tris[other_t1_id][1]<<" "<<mesh.tris[other_t1_id][2]<<endl;
        //                cout<<mesh.tag_feature_es[other_t1_id][0]<<" "<<mesh.tag_feature_es[other_t1_id][1]<<" "<<mesh.tag_feature_es[other_t1_id][2]<<endl;
        //            }
        //            tmp = set_intersection(mesh.tri_vertices[v2_id].conn_tris, mesh.tri_vertices[v_id].conn_tris);
        //            if (tmp.size() == 2) {
        //                int other_t2_id = tmp[0] == t_id ? tmp[1] : tmp[0];
        //                cout<<"other_t2_id "<<other_t2_id<<": "<<mesh.tris[other_t2_id][0]<<" "<<mesh.tris[other_t2_id][1]<<" "<<mesh.tris[other_t2_id][2]<<endl;
        //                cout<<mesh.tag_feature_es[other_t2_id][0]<<" "<<mesh.tag_feature_es[other_t2_id][1]<<" "<<mesh.tag_feature_es[other_t2_id][2]<<endl;
        //                for (int j = 0; j < 3; j++) {
        //                    if (mesh.tris[other_t2_id][j] != v1_id && mesh.tris[other_t2_id][j] != v_id) {
        //                        mesh.tag_feature_es[other_t2_id][j] = feature_id;
        //                        break;
        //                    }
        //                }
        //                cout<<"other_t2_id "<<other_t2_id<<": "<<mesh.tris[other_t2_id][0]<<" "<<mesh.tris[other_t2_id][1]<<" "<<mesh.tris[other_t2_id][2]<<endl;
        //                cout<<mesh.tag_feature_es[other_t2_id][0]<<" "<<mesh.tag_feature_es[other_t2_id][1]<<" "<<mesh.tag_feature_es[other_t2_id][2]<<endl;
        //            }
        //
        //            //snap v_id to curve
        //            t = feature::features[feature_id]->inv_eval(mesh.tri_vertices[v_id].posf, t, t1, t2);
        //            mesh.tri_vertices[v_id].feature_infos.push_back({{(double) feature_id, t}});
        //
        //            check(mesh);
        //            //do normal edge splitting on this edge and return true
        ////            split_normal_edge();
        //            return false;
        //        }

        //        if (is_edge_out_of_envelop(mesh, b_tree, pf, mesh.tri_vertices[v1_id].posf)
        //            || is_edge_out_of_envelop(mesh, b_tree, pf, mesh.tri_vertices[v2_id].posf))
        //            return false;
    }
    else
    {
        split_normal_edge();
        //        pf = mid_point(mesh.tri_vertices[v1_id].posf, mesh.tri_vertices[v2_id].posf);
        //        p = Point_2(pf[0], pf[1]);
        //        for (int t_id:n_t_ids) {
        //            for (int j = 0; j < 3; j++) {
        //                if ((mesh.tris[t_id][j] == v1_id || mesh.tris[t_id][j] == v2_id)
        //                    && !is_valid_inversion(pf, mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].posf,
        //                                           mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].posf,
        //                                           p, mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].pos,
        //                                           mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].pos)) {
        //                    is_rounded = false;
        //                    break;
        //                }
        //            }
        //            if (!is_rounded)
        //                break;
        //        }
        //
        //        if (!is_rounded) {
        //            p = mid_point(mesh.tri_vertices[v1_id].pos, mesh.tri_vertices[v2_id].pos);
        //            pf[0] = p[0].to_double();
        //            pf[1] = p[1].to_double();
        //        }
    }

    std::vector<std::vector<Point_2f>> vec_new_nodes;
    if (mesh.is_curved)
    {
        for (int t_id : n_t_ids)
        {
            if (mesh.tri_nodes[t_id].empty())
                continue;

            int k1 = std::find(mesh.tris[t_id].begin(), mesh.tris[t_id].end(), v1_id) - mesh.tris[t_id].begin();
            int k2 = std::find(mesh.tris[t_id].begin(), mesh.tris[t_id].end(), v2_id) - mesh.tris[t_id].begin();
            for (int j = 0; j < 2; j++)
            {
                vec_new_nodes.emplace_back();
                auto &new_nodes = vec_new_nodes.back();

                auto tag_feature_e = mesh.tag_feature_es[t_id];
                std::array<TriVertex, 3> new_t = {{mesh.tri_vertices[mesh.tris[t_id][0]],
                                                   mesh.tri_vertices[mesh.tris[t_id][1]],
                                                   mesh.tri_vertices[mesh.tris[t_id][2]]}};
                if (j == 0)
                {
                    new_t[k1].posf = pf;
                    new_t[k1].pos = p;
                    new_t[k1].feature_infos.clear();
                    new_t[k1].feature_infos.push_back({{(double)feature_id, t}});
                    tag_feature_e[k2] = -1;
                    if (tag_feature_e[0] < 0 && tag_feature_e[1] < 0 && tag_feature_e[2] < 0)
                        continue;
                    feature::get_new_nodes(tag_feature_e, new_t, new_nodes);
                    if (new_nodes.empty())
                        continue;
                    if (!feature::is_valid_inversion({{new_t[0].posf, new_t[1].posf, new_t[2].posf}}, new_nodes))
                    {
                        return false;
                    }
                }
                else
                {
                    new_t[k2].posf = pf;
                    new_t[k2].pos = p;
                    new_t[k2].feature_infos.clear();
                    new_t[k2].feature_infos.push_back({{(double)feature_id, t}});
                    tag_feature_e[k1] = -1;
                    if (tag_feature_e[0] < 0 && tag_feature_e[1] < 0 && tag_feature_e[2] < 0)
                        continue;
                    feature::get_new_nodes(tag_feature_e, new_t, new_nodes);
                    if (new_nodes.empty())
                        continue;
                    if (!feature::is_valid_inversion({{new_t[0].posf, new_t[1].posf, new_t[2].posf}}, new_nodes))
                    {
                        return false;
                    }
                }
            }
        }
    }

    ///real update
    int v_id = mesh.v_empty_slot_start;
    if (mesh.v_empty_slot_start == mesh.tri_vertices.size())
    {
        mesh.tri_vertices.resize(mesh.tri_vertices.size() + 1);
        mesh.v_is_removed.resize(mesh.v_is_removed.size() + 1);
        mesh.v_empty_slot_start = mesh.tri_vertices.size();
    }
    else
    {
        mesh.v_empty_slot_start =
            std::find(mesh.v_is_removed.begin() + mesh.v_empty_slot_start + 1, mesh.v_is_removed.end(), true) - mesh.v_is_removed.begin();
    }
    mesh.v_is_removed[v_id] = false;
    mesh.tri_vertices[v_id] = TriVertex(); //must construct a new one!!!
    mesh.tri_vertices[v_id].pos = p;
    mesh.tri_vertices[v_id].posf = pf;
    mesh.tri_vertices[v_id].is_rounded = is_rounded;
    bool is_boudnary_e = is_boundary_edge(mesh, v1_id, v2_id);
    mesh.tri_vertices[v_id].is_on_boundary = is_boudnary_e;
    bool is_bbox_e = is_bbox_edge(mesh, v1_id, v2_id);
    mesh.tri_vertices[v_id].is_on_bbox = is_bbox_e;
    if (feature_id >= 0)
        mesh.tri_vertices[v_id].feature_infos.push_back(std::array<double, 2>({(double)feature_id, t}));
    mesh.tri_vertices[v_id].scale = (mesh.tri_vertices[v1_id].scale + mesh.tri_vertices[v2_id].scale) / 2;

    int old_tris_size = mesh.tris.size();
    std::array<int, 2> new_t_ids = {mesh.t_empty_slot_start, -1};
    if (new_t_ids[0] == mesh.tris.size())
    {
        mesh.tris.resize(mesh.tris.size() + n_t_ids.size());
        new_t_ids[0] = mesh.tris.size() - 1;
        new_t_ids[1] = mesh.tris.size() - 2;
        mesh.t_empty_slot_start = mesh.tris.size();
    }
    else
    {
        mesh.t_empty_slot_start =
            std::find(mesh.t_is_removed.begin() + mesh.t_empty_slot_start + 1, mesh.t_is_removed.end(), true) - mesh.t_is_removed.begin();
        for (int n = 0; n < int(n_t_ids.size()) - 1; n++)
        {
            if (mesh.t_empty_slot_start == mesh.tris.size())
            {
                mesh.tris.resize(mesh.tris.size() + 1);
                new_t_ids[1] = mesh.tris.size() - 1;
                mesh.t_empty_slot_start = mesh.tris.size();
            }
            else
            {
                new_t_ids[1] = mesh.t_empty_slot_start;
                mesh.t_empty_slot_start =
                    std::find(mesh.t_is_removed.begin() + mesh.t_empty_slot_start + 1, mesh.t_is_removed.end(),
                              true) -
                    mesh.t_is_removed.begin();
            }
        }
    }
    if (mesh.tris.size() > old_tris_size)
    {
        mesh.t_is_removed.resize(mesh.tris.size());
        mesh.t_quality.resize(mesh.tris.size());
        mesh.is_boundary_es.resize(mesh.tris.size());
        mesh.is_bbox_es.resize(mesh.tris.size());
        mesh.tag_feature_es.resize(mesh.tris.size());
        mesh.tag_secondary_feature_es.resize(mesh.tris.size());

        mesh.tri_nodes.resize(mesh.tris.size());
    }

    //update conn_tris
    mesh.tri_vertices[v_id].conn_tris.clear();
    for (int i = 0; i < n_t_ids.size(); i++)
    {
        mesh.tris[new_t_ids[i]] = mesh.tris[n_t_ids[i]];
        mesh.t_is_removed[new_t_ids[i]] = false;
        mesh.t_quality[new_t_ids[i]] = -1;
        mesh.t_quality[n_t_ids[i]] = -1; //need to update later
        for (int j = 0; j < 3; j++)
        {
            if (mesh.tris[n_t_ids[i]][j] == v1_id)
                mesh.tris[n_t_ids[i]][j] = v_id;
            if (mesh.tris[new_t_ids[i]][j] == v2_id)
                mesh.tris[new_t_ids[i]][j] = v_id;
            if (mesh.tris[n_t_ids[i]][j] != v1_id && mesh.tris[n_t_ids[i]][j] != v2_id)
            {
                mesh.tri_vertices[mesh.tris[n_t_ids[i]][j]].conn_tris.insert(new_t_ids[i]);
                new_es.push_back(std::array<int, 2>({v_id, mesh.tris[n_t_ids[i]][j]}));
            }
        }
        mesh.tri_vertices[v_id].conn_tris.insert(new_t_ids[i]);
        mesh.tri_vertices[v_id].conn_tris.insert(n_t_ids[i]);

        mesh.tri_vertices[v1_id].conn_tris.erase(n_t_ids[i]);
        mesh.tri_vertices[v1_id].conn_tris.insert(new_t_ids[i]);
    }
    new_es.push_back(std::array<int, 2>({v_id, v1_id}));
    new_es.push_back(std::array<int, 2>({v_id, v2_id}));

    for (int i = 0; i < n_t_ids.size(); i++)
    {
        if (n_t_ids[i] < 0)
            break;

        mesh.is_boundary_es[new_t_ids[i]] = mesh.is_boundary_es[n_t_ids[i]];
        mesh.is_bbox_es[new_t_ids[i]] = mesh.is_bbox_es[n_t_ids[i]];
        mesh.tag_feature_es[new_t_ids[i]] = mesh.tag_feature_es[n_t_ids[i]];
        mesh.tag_secondary_feature_es[new_t_ids[i]] = mesh.tag_secondary_feature_es[n_t_ids[i]];

        for (int j = 0; j < 3; j++)
        {
            if (mesh.tris[n_t_ids[i]][j] == v1_id || mesh.tris[n_t_ids[i]][j] == v2_id)
            {
                mesh.is_boundary_es[n_t_ids[i]][j] = false;
                mesh.is_bbox_es[n_t_ids[i]][j] = false;
                mesh.tag_feature_es[n_t_ids[i]][j] = -1;
                mesh.tag_secondary_feature_es[n_t_ids[i]][j] = -1;
            }
            if (mesh.tris[new_t_ids[i]][j] == v1_id || mesh.tris[new_t_ids[i]][j] == v2_id)
            {
                mesh.is_boundary_es[new_t_ids[i]][j] = false;
                mesh.is_bbox_es[new_t_ids[i]][j] = false;
                mesh.tag_feature_es[new_t_ids[i]][j] = -1;
                mesh.tag_secondary_feature_es[new_t_ids[i]][j] = -1;
            }
        }
    }

    //update nodes
    if (!vec_new_nodes.empty())
    {
        int cnt = 0;
        for (int i = 0; i < n_t_ids.size(); i++)
        {
            if (!mesh.tri_nodes[n_t_ids[i]].empty())
            {
                auto &nodes = vec_new_nodes[cnt];
                int t_id = n_t_ids[i];
                if (nodes.size() == mesh.tri_nodes[t_id].size())
                {
                    for (int j = 0; j < mesh.tri_nodes[t_id].size(); j++)
                        mesh.nodes[mesh.tri_nodes[t_id][j]] = nodes[j];
                }
                else
                {
                    for (int j = 0; j < mesh.tri_nodes[t_id].size(); j++)
                        mesh.n_is_removed[mesh.tri_nodes[t_id][j]] = true;
                    mesh.tri_nodes[t_id].clear();
                    for (auto &n : nodes)
                    {
                        mesh.nodes.push_back(n);
                        mesh.n_is_removed.push_back(false);
                        mesh.tri_nodes[t_id].push_back(mesh.nodes.size() - 1);
                    }
                }

                auto &nodes1 = vec_new_nodes[cnt + 1];
                t_id = new_t_ids[i];
                for (auto &n : nodes1)
                {
                    mesh.nodes.push_back(n);
                    mesh.n_is_removed.push_back(false);
                    mesh.tri_nodes[t_id].push_back(mesh.nodes.size() - 1);
                }

                cnt += 2;
            }
        }
    }

    return true;
}

bool triwild::optimization::is_inside_pocket(MeshData &mesh, const Point_2f &p1, const Point_2f &p2,
                                             const int feature_id, const double t1, const double t2, const Point_2f &p)
{
    const int N = (p1 - p2).length() / mesh.dd * 1.2 + 1;

    auto side_of_segment = [&](const Vector_2f &v, const Point_2f &vp, const Point_2f &p) -> int {
        Vector_2f n(-v.y, v.x);
        Vector_2f vpp = p - vp;
        double side = n.dot(vpp);

        if (side < 0)
            return -1;
        else if (side > 0)
            return 1;
        else
            return 0;
    };

    Vector_2f v = p2 - p1;
    int side = side_of_segment(v, p1, p);

    Point_2f prev_p = p2;
    Point_2f cur_p;
    bool is_inside = true;
    for (double n = 0; n <= N; n++)
    {
        cur_p = feature::features[feature_id]->eval(t1 * (n / N) + t2 * (1 - n / N));
        if (n == 0 && (cur_p - p2).length_2() > 1e-16)
        {
            v = cur_p - p2;
            if (side_of_segment(v, p2, p) != side)
            {
                is_inside = false;
                break;
            }
            prev_p = cur_p;
            continue;
        }

        v = cur_p - prev_p;
        if (side_of_segment(v, prev_p, p) != side)
        {
            is_inside = false;
            break;
        }
        prev_p = cur_p;

        if (n == N && (cur_p - p1).length_2() > 1e-16)
        {
            v = p1 - prev_p;
            if (side_of_segment(v, prev_p, p) != side)
            {
                is_inside = false;
                break;
            }
        }
    }

    return is_inside;
}