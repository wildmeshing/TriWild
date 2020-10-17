// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "edge_collapsing.h"
#include "optimization.h"
#include "feature.h"

void triwild::optimization::edge_collapsing(MeshData &mesh, GEO::MeshFacetsAABB &b_tree) {
    ////init
    std::vector<std::array<int, 2>> edges;
    get_all_edges(mesh, edges);

    std::priority_queue<ElementInQueue_s, std::vector<ElementInQueue_s>, cmp_s> queue_s;
    auto push_edge = [&](std::array<int, 2> &e) -> bool {
        bool is_pushed = false;
        double l = std::sqrt(edge_length_2(mesh, e[0], e[1]));
        if (!mesh.is_limit_length || l <
                                     (mesh.tri_vertices[e[0]].scale + mesh.tri_vertices[e[1]].scale) / 2 *
                                     mesh.ideal_edge_length * (4.0 / 5.0)) {
            for (int j = 0; j < 2; j++) {
                if (mesh.tri_vertices[e[0]].is_freezed) {
                    std::swap(e[0], e[1]);
                    continue;
                }
                if (mesh.tri_vertices[e[0]].is_on_bbox && !is_bbox_edge(mesh, e[0], e[1])) {
                    std::swap(e[0], e[1]);
                    continue;
                }
                queue_s.push(ElementInQueue_s(e, l));
                is_pushed = true;
                std::swap(e[0], e[1]);
            }
        }
        return is_pushed;
    };

    for (auto &e : edges) {
        push_edge(e);
    }

    ////collapse
    int g_timestamp = 0;
    std::vector<int> t_timestamp(mesh.tris.size(), 0);
    std::vector<std::array<int, 2>> old_es;
    std::vector<int> e_timestamp;
    auto loop_edges = [&]() -> void {
        if (queue_s.size() == 0)
            return;
        cout << "queue_s.size() = " << queue_s.size() << endl;
        int cnt = 0;
        int tmp_cnt = 0;
        while (!queue_s.empty()) {
            std::array<int, 2> v_ids = queue_s.top().v_ids;
            double old_l = queue_s.top().weight;
            queue_s.pop();

            if (!is_valid_edge(mesh, v_ids[0], v_ids[1]))
                continue;
            double l = std::sqrt(edge_length_2(mesh, v_ids[0], v_ids[1]));
            if (l != old_l)
                continue;
            if (!queue_s.empty() && v_ids == queue_s.top().v_ids)
                continue;

            std::vector<std::array<int, 2>> new_es;
            if (!collapse_an_edge(mesh, b_tree, v_ids[0], v_ids[1], l, g_timestamp, t_timestamp, new_es)) {
                old_es.push_back(v_ids);
                e_timestamp.push_back(g_timestamp);
            } else {
                for (auto &e : new_es) {
                    push_edge(e);
                    tmp_cnt++;
                }

                cnt++;
            }
        }
        cout << "success " << cnt << endl;
        cout << "re-push " << tmp_cnt << endl;
    };

    ////post processing
    do {
        old_es.clear();
        e_timestamp.clear();
        loop_edges();
        for (int i = 0; i < old_es.size(); i++) {
            auto &e = old_es[i];
            if (!is_valid_edge(mesh, e[0], e[1]))
                continue;
            double l = std::sqrt(edge_length_2(mesh, e[0], e[1]));
            if (mesh.is_limit_length && l >=
                                        (mesh.tri_vertices[e[0]].scale + mesh.tri_vertices[e[1]].scale) / 2 *
                                        mesh.ideal_edge_length * (4.0 / 5.0))
                continue;

            for (int t_id : mesh.tri_vertices[e[0]].conn_tris) {
                if (t_timestamp[t_id] > e_timestamp[i]) {
                    queue_s.push(ElementInQueue_s(e, l));
                    break;
                }
            }

            //todo: repush unpushed edges and update timestamps
        }
    } while (old_es.size() > 0);

    if (mesh.is_curved) { ; //todo: clean up dirty nodes after simplification
    }
}

bool triwild::optimization::collapse_an_edge(MeshData &mesh, GEO::MeshFacetsAABB &b_tree, int v1_id, int v2_id, double l,
                                             int &g_timestamp, std::vector<int> &t_timestamp,
                                             std::vector<std::array<int, 2>> &new_es) {
    if (mesh.tri_vertices[v1_id].is_on_boundary && mesh.tri_vertices[v2_id].is_on_bbox)
        return false; //bondary and bbox are not allowed to connect

    if (mesh.tri_vertices[v1_id].is_on_point) { //circle envelop
        double dis = (mesh.tri_vertices[v1_id].input_posf - mesh.tri_vertices[v2_id].posf).length_2();
        if (dis > mesh.epsilon_2)
            return false;
    }

    //remove isolated surface vertex (but it could still be a feature vertex
    if (mesh.tri_vertices[v1_id].is_on_boundary && is_isolated_vertex(mesh, v1_id)) {
        mesh.tri_vertices[v1_id].is_on_boundary = false;
        mesh.tri_vertices[v1_id].is_on_point = false;
    }

    ////check
    std::vector<int> changed_t_ids;
    std::vector<int> removed_t_ids;
    std::vector<double> new_qs;
    double old_max_energy = 0;
    for (int t_id : mesh.tri_vertices[v1_id].conn_tris) {
        if (mesh.t_quality[t_id] > old_max_energy)
            old_max_energy = mesh.t_quality[t_id];
    }
    bool is_snap_features = false;

    //    if(mesh.is_curved)
    //        old_max_energy = mesh.stop_energy > old_max_energy ? mesh.stop_energy : old_max_energy;

    if (mesh.tri_vertices[v1_id].feature_infos.size() > 0) {
        int feature_id = get_feature_edge_tag(mesh, v1_id, v2_id);

        if (feature_id >= 0) {
            if (!mesh.is_curved && !is_valid_feature_edge_close(mesh, v1_id, mesh.tri_vertices[v2_id].posf,
                                                                mesh.tri_vertices[v2_id].get_t(feature_id)))
                return false;
        } else {
            return false;
        }
    }

    for (int t_id : mesh.tri_vertices[v1_id].conn_tris) {
        if (std::find(mesh.tris[t_id].begin(), mesh.tris[t_id].end(), v2_id) == mesh.tris[t_id].end()) {
            int j = std::find(mesh.tris[t_id].begin(), mesh.tris[t_id].end(), v1_id) - mesh.tris[t_id].begin();
            if (!mesh.is_curved || (mesh.is_curved && mesh.tri_nodes[t_id].size() == 0)) {
                if (!is_valid_inversion(mesh.tri_vertices[v2_id].posf,
                                        mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].posf,
                                        mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].posf,
                                        mesh.tri_vertices[v2_id].pos,
                                        mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].pos,
                                        mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].pos)) {
                    //                cout << "I";
                    return false;
                }
            }

            double new_q;
            if (!is_valid_quality(mesh.tri_vertices[v2_id].posf,
                                  mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].posf,
                                  mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].posf,
                                  old_max_energy, new_q, true) //new_q need to be computed anyway
                && l != 0) {
                //                cout << "Q";
                return false;
            }

            if (mesh.tri_vertices[v1_id].is_on_boundary
                //                && mesh.tri_vertices[v1_id].feature_infos.size() == 0
                //                && !mesh.is_curved //todo: should be "skip checking feature edges"
                && l != 0 && !is_valid_envelop(mesh, b_tree, v1_id, mesh.tri_vertices[v2_id].posf)) {
                //                cout << "E";
                return false;
            }

//            if (mesh.tri_vertices[v1_id].feature_infos.size() > 0) {
//                //                if(!is_valid_feature_edge_length(mesh, v1_id, mesh.tri_vertices[v2_id].posf))
//                int feature_id = get_feature_edge_tag(mesh, v1_id, v2_id);
//
//                if (feature_id >= 0) {
//                    if (!mesh.is_curved && !is_valid_feature_edge_close(mesh, v1_id, mesh.tri_vertices[v2_id].posf,
//                                                                        mesh.tri_vertices[v2_id].get_t(feature_id)))
//                        return false;
//                } else {
//                    return false;
//
//                    // if (mesh.tri_vertices[v2_id].feature_infos.size() == 0)
//                    //     return false;
//
//                    // //snap two features --> check feature epsilon
//                    // int v2_feature_id = -1;
//                    // double v2_t;
//                    // for (auto &info: mesh.tri_vertices[v2_id].feature_infos) {
//                    //     double dis = feature::features[info[0]]->distance(mesh.tri_vertices[v1_id].posf, info[1]);//todo
//                    //     if (dis < mesh.feature_epsilon) {//todo: feature epsilon
//                    //         v2_feature_id = info[0];
//                    //         v2_t = info[1];
//                    //         break;
//                    //     }
//                    // }
//                    // if (v2_feature_id <
//                    //     mesh.tri_vertices[v1_id].feature_infos[0][0])//always snap to vertex with larger feature_id
//                    //     return false;
//
//                    // double dd = mesh.dd;//todo: sampling density
//                    // //sample feature edges of v1
//                    // Point_2f &p0 = mesh.tri_vertices[v1_id].posf;
//                    // std::unordered_set<int> n_v_ids;
//                    // for (int t_id:mesh.tri_vertices[v1_id].conn_tris) {
//                    //     for (int j = 0; j < 3; j++)
//                    //         if (mesh.tag_feature_es[t_id][j] >= 0) {
//                    //             n_v_ids.insert(mesh.tag_feature_es[t_id][(j + 1) % 3] == v1_id ?
//                    //                            mesh.tag_feature_es[t_id][(j + 2) % 3] :
//                    //                            mesh.tag_feature_es[t_id][(j + 1) % 3]);
//                    //             break;
//                    //         }
//                    // }
//                    // assert(n_v_ids.size() == 2);
//                    // //compute dis for samples
//                    // for (int v_id:n_v_ids) {
//                    //     Point_2f &p = mesh.tri_vertices[v_id].posf;
//                    //     int N = std::sqrt(edge_length_2(p, p0)) / mesh.dd + 0.5;
//                    //     for (double n = 1; n <= N; n++) {
//                    //         Point_2f tmp = p0 * (n / N) + p * (N - n) / N;
//                    //         if (feature::features[v2_feature_id]->distance(mesh.tri_vertices[v1_id].posf, v2_t)
//                    //             > mesh.feature_epsilon)
//                    //             return false;
//                    //     }
//                    // }
//                    // is_snap_features = true;
//                }
//            }

            new_qs.push_back(new_q);
            changed_t_ids.push_back(t_id);
        } else
            removed_t_ids.push_back(t_id);
    }

    //move upward for cerved elements
    std::vector<int> opp_v_ids;
    std::vector<std::array<bool, 2>> opp_is_boundary_es;
    std::vector<std::array<bool, 2>> opp_is_bbox_es;
    std::vector<std::array<int, 2>> opp_tag_feature_es;
    std::vector<std::array<int, 2>> opp_tag_secondary_feature_es;
    //    mesh.v_is_removed[v1_id] = true;
    for (int t_id : removed_t_ids) {
        //        mesh.t_is_removed[t_id] = true;
        //        mesh.tri_vertices[v2_id].conn_tris.erase(t_id);

        //feature
        std::array<bool, 2> is_boundary_e = {false, false};
        std::array<bool, 2> is_bbox_e = {false, false};
        std::array<int, 2> tag_feature_e = {-1, -1};
        std::array<int, 2> tag_secondary_feature_e = {-1, -1};
        for (int j = 0; j < 3; j++) {
            if (mesh.tris[t_id][j] != v1_id && mesh.tris[t_id][j] != v2_id) {
                //                mesh.tri_vertices[mesh.tris[t_id][j]].conn_tris.erase(t_id);
                opp_v_ids.push_back(mesh.tris[t_id][j]);
            } else if (mesh.tris[t_id][j] == v1_id) {
                is_boundary_e[0] = mesh.is_boundary_es[t_id][j];
                is_bbox_e[0] = mesh.is_bbox_es[t_id][j];
                tag_feature_e[0] = mesh.tag_feature_es[t_id][j];
                tag_secondary_feature_e[0] = mesh.tag_secondary_feature_es[t_id][j];
            } else {
                is_boundary_e[1] = mesh.is_boundary_es[t_id][j];
                is_bbox_e[1] = mesh.is_bbox_es[t_id][j];
                tag_feature_e[1] = mesh.tag_feature_es[t_id][j];
                tag_secondary_feature_e[1] = mesh.tag_secondary_feature_es[t_id][j];
            }
        }

        //        if(mesh.is_curved){
        if (tag_feature_e[0] >= 0 && tag_feature_e[1] >= 0)
            return false;
        //        }
        //        if(mesh.is_curved)
        //            cout<<tag_feature_e[0]<<" "<<tag_feature_e[1]<<endl;

        opp_is_boundary_es.push_back(is_boundary_e);
        opp_is_bbox_es.push_back(is_bbox_e);
        opp_tag_feature_es.push_back(tag_feature_e);
        opp_tag_secondary_feature_es.push_back(tag_secondary_feature_e);
    }

    //check for curved mesh
    std::vector<std::vector<Point_2f>> vec_new_nodes;
    if (mesh.is_curved) {
        std::vector<int> curved_changed_t_ids = changed_t_ids;
        for (int i = 0; i < opp_v_ids.size(); i++) {
            auto tmp = set_intersection(mesh.tri_vertices[v2_id].conn_tris,
                                        mesh.tri_vertices[opp_v_ids[i]].conn_tris);
            if (tmp.size() > 1)
                curved_changed_t_ids.push_back(tmp[0] == removed_t_ids[i] ? tmp[1] : tmp[0]);
        }
        for (int t_id : curved_changed_t_ids) {
            vec_new_nodes.emplace_back(std::vector<Point_2f>());
            auto &new_nodes = vec_new_nodes.back();
            auto tag_feature_e = mesh.tag_feature_es[t_id];
            auto new_t = mesh.tris[t_id];

            int j = std::find(mesh.tris[t_id].begin(), mesh.tris[t_id].end(), v1_id) - mesh.tris[t_id].begin();
            if (j < 3) {
                if (mesh.tris[t_id][(j + 1) % 3] == opp_v_ids[0] && opp_tag_feature_es[0][0] >= 0)
                    tag_feature_e[(j + 2) % 3] = opp_tag_feature_es[0][0];
                else if (mesh.tris[t_id][(j + 2) % 3] == opp_v_ids[0] && opp_tag_feature_es[0][0] >= 0)
                    tag_feature_e[(j + 1) % 3] = opp_tag_feature_es[0][0];
                if (opp_v_ids.size() == 2) {
                    if (mesh.tris[t_id][(j + 1) % 3] == opp_v_ids[1] && opp_tag_feature_es[1][0] >= 0)
                        tag_feature_e[(j + 2) % 3] = opp_tag_feature_es[1][0];
                    else if (mesh.tris[t_id][(j + 2) % 3] == opp_v_ids[1] && opp_tag_feature_es[1][0] >= 0)
                        tag_feature_e[(j + 1) % 3] = opp_tag_feature_es[1][0];
                }
                new_t[j] = v2_id;
            } else {
                j = std::find(mesh.tris[t_id].begin(), mesh.tris[t_id].end(), v2_id) - mesh.tris[t_id].begin();
                if (mesh.tris[t_id][(j + 1) % 3] == opp_v_ids[0] && opp_tag_feature_es[0][1] >= 0)
                    tag_feature_e[(j + 2) % 3] = opp_tag_feature_es[0][1];
                else if (mesh.tris[t_id][(j + 2) % 3] == opp_v_ids[0] && opp_tag_feature_es[0][1] >= 0)
                    tag_feature_e[(j + 1) % 3] = opp_tag_feature_es[0][1];
                if (opp_v_ids.size() == 2) {
                    if (mesh.tris[t_id][(j + 1) % 3] == opp_v_ids[1] && opp_tag_feature_es[1][1] >= 0)
                        tag_feature_e[(j + 2) % 3] = opp_tag_feature_es[1][1];
                    else if (mesh.tris[t_id][(j + 2) % 3] == opp_v_ids[1] && opp_tag_feature_es[1][1] >= 0)
                        tag_feature_e[(j + 1) % 3] = opp_tag_feature_es[1][1];
                }
            }
            if (tag_feature_e[0] < 0 && tag_feature_e[1] < 0 && tag_feature_e[2] < 0)
                continue;
            feature::get_new_nodes(tag_feature_e,
                                   {{mesh.tri_vertices[new_t[0]], mesh.tri_vertices[new_t[1]], mesh.tri_vertices[new_t[2]]}},
                                   new_nodes);
            if (new_nodes.empty())
                continue;
            if (!feature::is_valid_inversion({{mesh.tri_vertices[new_t[0]].posf,
                                                      mesh.tri_vertices[new_t[1]].posf,
                                                      mesh.tri_vertices[new_t[2]].posf}},
                                             new_nodes)) {
                //                cout << "I" << endl;
                return false;
            }
        }

        ////real update
        //delete/add nodes
        int cnt = 0;
        for (int t_id : removed_t_ids) {
            if (mesh.tri_nodes[t_id].size() > 0) {
                for (int i = 0; i < mesh.tri_nodes[t_id].size(); i++)
                    mesh.n_is_removed[mesh.tri_nodes[t_id][i]] = true;
                mesh.tri_nodes[t_id].clear();
            }
        }
        for (int t_id : curved_changed_t_ids) {
            auto &nodes = vec_new_nodes[cnt++];
            if (nodes.size() == mesh.tri_nodes[t_id].size()) {
                for (int i = 0; i < mesh.tri_nodes[t_id].size(); i++) {
                    mesh.nodes[mesh.tri_nodes[t_id][i]] = nodes[i];
                    assert(mesh.n_is_removed[mesh.tri_nodes[t_id][i]] == false);
                }
            } else {
                for (int i = 0; i < mesh.tri_nodes[t_id].size(); i++)
                    mesh.n_is_removed[mesh.tri_nodes[t_id][i]] = true;
                mesh.tri_nodes[t_id].clear();
                for (auto &n : nodes) {
                    mesh.nodes.push_back(n);
                    mesh.n_is_removed.push_back(false);
                    mesh.tri_nodes[t_id].push_back(mesh.nodes.size() - 1);
                }
            }
        }
    }

    ////real update
    //boundary & feature
    mesh.tri_vertices[v2_id].is_on_boundary =
            mesh.tri_vertices[v1_id].is_on_boundary || mesh.tri_vertices[v2_id].is_on_boundary;
    mesh.tri_vertices[v2_id].is_on_bbox = mesh.tri_vertices[v1_id].is_on_bbox || mesh.tri_vertices[v2_id].is_on_bbox;
    if (is_snap_features) {
        //        mesh.tri_vertices[v2_id].feature_infos.push_back(
        //                std::array<double, 2>({mesh.tri_vertices[v1_id].feature_infos[0][0],
        //                                       feature::features[mesh.tri_vertices[v1_id].feature_infos[0][0]]->inv_eval(
        //                                               mesh.tri_vertices[v2_id].posf,
        //                                               mesh.tri_vertices[v1_id].feature_infos[0][1])}));
        //        mesh.tri_vertices[v2_id].is_freezed = true;
        //todo: when to unfreeze?
    }
    if (mesh.tri_vertices[v2_id].is_on_point && mesh.tri_vertices[v1_id].is_on_point)
        mesh.tri_vertices[v2_id].is_on_point = false;

    mesh.v_is_removed[v1_id] = true;
    for (int i = 0; i < removed_t_ids.size(); i++) {
        mesh.t_is_removed[removed_t_ids[i]] = true;
        mesh.tri_vertices[v2_id].conn_tris.erase(removed_t_ids[i]);
        mesh.tri_vertices[opp_v_ids[i]].conn_tris.erase(removed_t_ids[i]);
    }

    //    std::vector<int> opp_v_ids;
    //    std::vector<std::array<bool, 2>> opp_is_boundary_es;
    //    std::vector<std::array<bool, 2>> opp_is_bbox_es;
    //    std::vector<std::array<int, 2>> opp_tag_feature_es;
    //    mesh.v_is_removed[v1_id] = true;
    //    for (int t_id:removed_t_ids) {
    //        mesh.t_is_removed[t_id] = true;
    //        mesh.tri_vertices[v2_id].conn_tris.erase(t_id);
    //
    //        //feature
    //        std::array<bool, 2> is_boundary_e = {false, false};
    //        std::array<bool, 2> is_bbox_e = {false, false};
    //        std::array<int, 2> tag_feature_e = {-1, -1};
    //        for (int j = 0; j < 3; j++) {
    //            if (mesh.tris[t_id][j] != v1_id && mesh.tris[t_id][j] != v2_id) {
    //                mesh.tri_vertices[mesh.tris[t_id][j]].conn_tris.erase(t_id);
    //                opp_v_ids.push_back(mesh.tris[t_id][j]);
    //            } else if (mesh.tris[t_id][j] == v1_id) {
    //                is_boundary_e[0] = mesh.is_boundary_es[t_id][j];
    //                is_bbox_e[0] = mesh.is_bbox_es[t_id][j];
    //                tag_feature_e[0] = mesh.tag_feature_es[t_id][j];
    //            } else {
    //                is_boundary_e[1] = mesh.is_boundary_es[t_id][j];
    //                is_bbox_e[1] = mesh.is_bbox_es[t_id][j];
    //                tag_feature_e[1] = mesh.tag_feature_es[t_id][j];
    //            }
    //        }
    //        opp_is_boundary_es.push_back(is_boundary_e);
    //        opp_is_bbox_es.push_back(is_bbox_e);
    //        opp_tag_feature_es.push_back(tag_feature_e);
    //    }

    //update is_boundary_es, is_bbox_es
    //update tag_feature_es --> keep the feature with larger tag
    for (int i = 0; i < opp_v_ids.size(); i++) {
        int tmp_t1_id = -1;
        int tmp_t2_id = -1;

        //is_boundary_es
        if (opp_is_boundary_es[i][0] && !opp_is_boundary_es[i][1]) {
            //            std::vector<int> tmp = set_intersection(mesh.tri_vertices[v1_id].conn_tris,
            //                                                    mesh.tri_vertices[opp_v_ids[i]].conn_tris);
            //            if (tmp.size() == 0)
            //                continue;
            //            tmp_t1_id = tmp[0];
            tmp_t1_id = set_intersection(mesh.tri_vertices[v1_id].conn_tris,
                                         mesh.tri_vertices[opp_v_ids[i]].conn_tris)[0];
            for (int j = 0; j < 3; j++) {
                if (mesh.tris[tmp_t1_id][j] != opp_v_ids[i] && mesh.tris[tmp_t1_id][j] != v1_id) {
                    mesh.is_boundary_es[tmp_t1_id][j] = true;
                    break;
                }
            }
        } else if (!opp_is_boundary_es[i][0] && opp_is_boundary_es[i][1]) {
            //            std::vector<int> tmp = set_intersection(mesh.tri_vertices[v2_id].conn_tris,
            //                                                    mesh.tri_vertices[opp_v_ids[i]].conn_tris);
            //            if (tmp.size() == 0)
            //                continue;
            //            tmp_t2_id = tmp[0];
            tmp_t2_id = set_intersection(mesh.tri_vertices[v2_id].conn_tris,
                                         mesh.tri_vertices[opp_v_ids[i]].conn_tris)[0];
            //v2 could on bbox!
            for (int j = 0; j < 3; j++) {
                if (mesh.tris[tmp_t2_id][j] != opp_v_ids[i] && mesh.tris[tmp_t2_id][j] != v2_id) {
                    mesh.is_boundary_es[tmp_t2_id][j] = true;
                    break;
                }
            }
        }

        //is_bbox_es
        if (opp_is_bbox_es[i][0] && !opp_is_bbox_es[i][1]) {
            if (tmp_t1_id < 0) {
                tmp_t1_id = set_intersection(mesh.tri_vertices[v1_id].conn_tris,
                                             mesh.tri_vertices[opp_v_ids[i]].conn_tris)[0];
            }
            for (int j = 0; j < 3; j++) {
                if (mesh.tris[tmp_t1_id][j] != opp_v_ids[i] && mesh.tris[tmp_t1_id][j] != v1_id) {
                    mesh.is_bbox_es[tmp_t1_id][j] = true;
                    break;
                }
            }
        } else if (!opp_is_bbox_es[i][0] && opp_is_bbox_es[i][1]) {
            if (tmp_t2_id < 0) {
                tmp_t2_id = set_intersection(mesh.tri_vertices[v2_id].conn_tris,
                                             mesh.tri_vertices[opp_v_ids[i]].conn_tris)[0];
            }
            for (int j = 0; j < 3; j++) {
                if (mesh.tris[tmp_t2_id][j] != opp_v_ids[i] && mesh.tris[tmp_t2_id][j] != v2_id) {
                    mesh.is_bbox_es[tmp_t2_id][j] = true;
                    break;
                }
            }
        }

        //tag_feature_es
        if (opp_tag_feature_es[i][0] > opp_tag_feature_es[i][1]) {
            if (tmp_t1_id < 0) {
                tmp_t1_id = set_intersection(mesh.tri_vertices[v1_id].conn_tris,
                                             mesh.tri_vertices[opp_v_ids[i]].conn_tris)[0];
            }
            for (int j = 0; j < 3; j++) {
                if (mesh.tris[tmp_t1_id][j] != opp_v_ids[i] && mesh.tris[tmp_t1_id][j] != v1_id) {
                    mesh.tag_feature_es[tmp_t1_id][j] = opp_tag_feature_es[i][0];
                    break;
                }
            }
        } else if (opp_tag_feature_es[i][0] < opp_tag_feature_es[i][1]) {
            if (tmp_t2_id < 0) {
                tmp_t2_id = set_intersection(mesh.tri_vertices[v2_id].conn_tris,
                                             mesh.tri_vertices[opp_v_ids[i]].conn_tris)[0];
            }
            for (int j = 0; j < 3; j++) {
                if (mesh.tris[tmp_t2_id][j] != opp_v_ids[i] && mesh.tris[tmp_t2_id][j] != v2_id) {
                    mesh.tag_feature_es[tmp_t2_id][j] = opp_tag_feature_es[i][1];
                    break;
                }
            }
        }

        //tag_secondary_feature_es
        if (opp_tag_secondary_feature_es[i][0] > opp_tag_secondary_feature_es[i][1]) {
            if (tmp_t1_id < 0) {
                tmp_t1_id = set_intersection(mesh.tri_vertices[v1_id].conn_tris,
                                             mesh.tri_vertices[opp_v_ids[i]].conn_tris)[0];
            }
            for (int j = 0; j < 3; j++) {
                if (mesh.tris[tmp_t1_id][j] != opp_v_ids[i] && mesh.tris[tmp_t1_id][j] != v1_id) {
                    mesh.tag_secondary_feature_es[tmp_t1_id][j] = opp_tag_secondary_feature_es[i][0];
                    break;
                }
            }
        } else if (opp_tag_secondary_feature_es[i][0] < opp_tag_secondary_feature_es[i][1]) {
            if (tmp_t2_id < 0) {
                tmp_t2_id = set_intersection(mesh.tri_vertices[v2_id].conn_tris,
                                             mesh.tri_vertices[opp_v_ids[i]].conn_tris)[0];
            }
            for (int j = 0; j < 3; j++) {
                if (mesh.tris[tmp_t2_id][j] != opp_v_ids[i] && mesh.tris[tmp_t2_id][j] != v2_id) {
                    mesh.tag_secondary_feature_es[tmp_t2_id][j] = opp_tag_secondary_feature_es[i][1];
                    break;
                }
            }
        }
    }
    int cnt = 0;
    std::unordered_set<int> tmp_v_ids;
    for (int t_id : changed_t_ids) {
        mesh.t_quality[t_id] = new_qs[cnt++];
        int j = std::find(mesh.tris[t_id].begin(), mesh.tris[t_id].end(), v1_id) - mesh.tris[t_id].begin();
        mesh.tris[t_id][j] = v2_id;
        mesh.tri_vertices[v2_id].conn_tris.insert(t_id);
        tmp_v_ids.insert(mesh.tris[t_id][(j + 1) % 3]);
        tmp_v_ids.insert(mesh.tris[t_id][(j + 2) % 3]);
    }
    for (int v_id : opp_v_ids)
        tmp_v_ids.erase(v_id);
    for (int v_id : tmp_v_ids)
        new_es.push_back(std::array<int, 2>({v2_id, v_id}));

    //update timestamps
    g_timestamp++;
    for (int t_id : changed_t_ids) {
        t_timestamp[t_id] = g_timestamp;
    }

    return true;
}