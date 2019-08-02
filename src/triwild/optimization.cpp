// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "optimization.h"
#include "AMIPS.h"
#include "meshio.hpp"
#include "feature.h"

#include "edge_splitting.h"
#include "edge_collapsing.h"
#include "edge_swapping.h"
#include "vertex_smoothing.h"
#include "Point_2f.h"

#include "Logger.h"

#include<igl/Timer.h>

using namespace triwild;

void triwild::optimization::init(const Eigen::MatrixXd& V, const std::vector<std::array<int, 2>>& edges,
                                  const std::vector<std::vector<int>>& tag_boundary_es,
                                  MeshData& mesh, GEO::MeshFacetsAABB &b_tree) {
    ////vertices
    for (int i = 0; i < mesh.tri_vertices.size(); i++) {
        if (tag_boundary_es[i].size() > 0)
            mesh.tri_vertices[i].is_on_boundary = true;
    }

    ////features
    std::vector<int> map_b_to_f(edges.size(), -1);
    std::vector<int> map_b_to_2nd_f(edges.size(), -1);
    if (args.is_preserving_feature) {
        std::vector<std::unordered_set<int>> conn_es(V.rows());
        for (int i = 0; i < edges.size(); i++) {
            conn_es[edges[i][0]].insert(i);
            conn_es[edges[i][1]].insert(i);
        }

//        for (int feature_id; feature_id < feature::features.size(); feature_id++) {
//            auto &v_ids = feature::features[feature_id]->v_ids;
//            auto &paras = feature::features[feature_id]->paras;
//            mesh.tri_vertices[v_ids.front()].is_freezed = true;
//            mesh.tri_vertices[v_ids.back()].is_freezed = true;
//
//            for (int i = 0; i < v_ids.size() - 1; i++) {
//                int v1_id = v_ids[i];
//                int v2_id = v_ids[i + 1];
//                auto tmp = set_intersection(conn_es[v1_id], conn_es[v2_id]);
//                assert(tmp.size() == 1);
//
//                if (map_b_to_f[tmp[0]] >= 0) {
//                    int old_feature_id = map_b_to_f[tmp[0]];
//                    if (feature::features[old_feature_id]->degree <= feature::features[feature_id]->degree)
//                        map_b_to_f[tmp[0]] = feature_id;
//                } else
//                    map_b_to_f[tmp[0]] = feature_id;
//                mesh.tri_vertices[v1_id].feature_infos.push_back({{(double) feature_id, paras[i]}});
//                mesh.tri_vertices[v2_id].feature_infos.push_back({{(double) feature_id, paras[i + 1]}});
//            }
//        }
//        //todo: cleanup mesh.tri_vertices[v_id].feature_infos?
//        for (int i = 0; i < mesh.tri_vertices.size(); i++) {
//            if (mesh.tri_vertices[i].feature_infos.size() > 1)
//                mesh.tri_vertices[i].is_freezed = true;
//        }

        for (int feature_id = 0; feature_id < feature::secondary_features.size(); feature_id++) {
            auto &v_ids = feature::secondary_features[feature_id]->v_ids;
            for (int i = 0; i < v_ids.size() - 1; i++) {
                int v1_id = v_ids[i];
                int v2_id = v_ids[i + 1];
                auto tmp = set_intersection(conn_es[v1_id], conn_es[v2_id]);
//                assert(tmp.size() == 1);//duplicate features
                if(tmp.empty())
                    continue;
                if (map_b_to_2nd_f[tmp[0]] >= 0) {
                    int old_feature_id = map_b_to_2nd_f[tmp[0]];
                    if (feature::secondary_features[old_feature_id]->degree <= feature::secondary_features[feature_id]->degree)
                        map_b_to_2nd_f[tmp[0]] = feature_id;
                } else
                    map_b_to_2nd_f[tmp[0]] = feature_id;
            }
        }
    }

    ////faces
    const int tris_size = mesh.tris.size();
    mesh.t_quality.resize(tris_size);
    mesh.is_boundary_es.resize(tris_size);
    mesh.tag_feature_es.resize(tris_size);
    mesh.tag_secondary_feature_es.resize(tris_size);
    for (int t_id = 0; t_id < tris_size; t_id++) {
        auto &t = mesh.tris[t_id];
        mesh.t_quality[t_id] = (
                tri_energy(mesh.tri_vertices[t[0]].posf, mesh.tri_vertices[t[1]].posf, mesh.tri_vertices[t[2]].posf));

        mesh.is_boundary_es[t_id] = {{false, false, false}};
        mesh.tag_feature_es[t_id] = {{-1, -1, -1}};
        mesh.tag_secondary_feature_es[t_id] = {{-1, -1, -1}};
        for (int j = 0; j < 3; j++) {
            mesh.tri_vertices[t[j]].conn_tris.insert(t_id);

            int v1_id = t[(j + 1) % 3];
            int v2_id = t[(j + 2) % 3];
            std::vector<int> b_tags = optimization::set_intersection(tag_boundary_es[v1_id], tag_boundary_es[v2_id]);
            if (!b_tags.empty()) {
                mesh.is_boundary_es[t_id][j] = true;

                //assert(b_tags.size() == 1);//wrong, the input edges could have overlapping
                if (args.is_preserving_feature) {
//                    //major feature
//                    int feature_id = map_b_to_f[b_tags[0]];
//                    mesh.tag_feature_es[t_id][j] = feature_id;
//                    if (!mesh.tri_vertices[v1_id].has_feature(feature_id))
//                        mesh.tri_vertices[v1_id].feature_infos.push_back({{(double) feature_id, -1}});
//                    if (!mesh.tri_vertices[v2_id].has_feature(feature_id))
//                        mesh.tri_vertices[v2_id].feature_infos.push_back({{(double) feature_id, -1}});
                    //secondary feature
                    int b_tag = -1;
                    for (int tag:b_tags) {
                        if (map_b_to_2nd_f[tag] < 0)
                            continue;
                        if (b_tag < 0) {
                            b_tag = tag;
                            continue;
                        }
                        if ((feature::secondary_features[map_b_to_2nd_f[tag]]->degree >
                             feature::secondary_features[map_b_to_2nd_f[b_tag]]->degree)
                            || (feature::secondary_features[map_b_to_2nd_f[tag]]->degree ==
                                feature::secondary_features[map_b_to_2nd_f[b_tag]]->degree
                                && map_b_to_2nd_f[tag] > map_b_to_2nd_f[b_tag]))
                            b_tag = tag;
                    }
                    if (b_tag < 0 || map_b_to_2nd_f[b_tag] < 0)
                        continue;
                    mesh.tag_secondary_feature_es[t_id][j] = map_b_to_2nd_f[b_tag];
                }
            }
        }
    }
    mesh.v_is_removed = std::vector<bool>(mesh.tri_vertices.size(), false);
    mesh.t_is_removed = std::vector<bool>(mesh.tris.size(), false);

    //tag bbox
    mesh.is_bbox_es = std::vector<std::array<bool, 3>>(tris_size, std::array<bool, 3>({false, false, false}));
    std::vector<std::array<int, 2>> tri_edges;
    get_all_edges(mesh, tri_edges);
    for (size_t i = 0; i < tri_edges.size(); i++) {
        auto &e = tri_edges[i];
        std::vector<int> tmp = set_intersection(mesh.tri_vertices[e[0]].conn_tris, mesh.tri_vertices[e[1]].conn_tris);
        if (tmp.size() == 1) {
            for (int j = 0; j < 3; j++)
                if (mesh.tris[tmp[0]][j] != e[0] && mesh.tris[tmp[0]][j] != e[1]) {
                    mesh.is_bbox_es[tmp[0]][j] = true;
                    mesh.tri_vertices[e[0]].is_on_bbox = true;
                    mesh.tri_vertices[e[1]].is_on_bbox = true;
                    break;
                }
        }
    }

    ////parameters
    mesh.ideal_edge_length = args.target_edge_len;
    double epsilon_user = args.diagonal_len * args.epsilon;
    if (args.i_dd < 0)
        mesh.dd = epsilon_user;
    else
        mesh.dd = args.diagonal_len / args.i_dd;
    double dd_error = mesh.dd / 2;
    mesh.epsilon = epsilon_user - dd_error;
    mesh.epsilon_2 = mesh.epsilon * mesh.epsilon;

//    if(args.stage > 1){
//        mesh.epsilon /= 2;
//        mesh.epsilon_2 /= 4;
//    }

    mesh.original_epsilon = mesh.epsilon;
    mesh.epsilon /= args.stage;
    mesh.epsilon_2 /= (args.stage * args.stage);

    mesh.max_its = args.max_its;
    mesh.stop_energy = args.stop_quality;

    mesh.v_empty_slot_start = mesh.tri_vertices.size();
    mesh.t_empty_slot_start = mesh.tris.size();

    for (auto &v:mesh.tri_vertices) {
        v.max_scale = 1;
        v.scale = 1;
    }
//    mesh.min_scalar = mesh.epsilon / mesh.ideal_edge_length;

    ////map features
    if (args.is_preserving_feature) {
        feature::map_feature2mesh(mesh);
//        feature::complete_mesh_paras(mesh);
//        feature::visualize_features(mesh);
//        mesh.min_scalar = (args.diagonal_len * 1e-4) / mesh.ideal_edge_length;
//        mesh.min_scalar = 1e-4;//give an absolute value
    }
    mesh.min_scalar = args.min_edge_length / mesh.ideal_edge_length;

    ////rounding
    for (int i = 0; i < mesh.tri_vertices.size(); i++)
        round_a_vertex(mesh, i);

    ////check & log
    output_info(mesh);
    double max_dis = check_envelope(mesh, b_tree, 1.6);
    if (max_dis > 0) {
        cout << "!check_envelope " << max_dis << " " << mesh.epsilon_2 << endl;
        output_boundary(mesh);
        pausee();
    }

//    output_boundary(mesh);
//    pausee();
}

void triwild::optimization::refine(MeshData& mesh, GEO::MeshFacetsAABB &b_tree, const std::array<int, 4>& ops) {
    //preprocessing
    if(!args.is_preserving_feature)
        mesh.is_limit_length = false;
    cout<<"//////////////// preprocessing ////////////////"<<endl;
    double flat_feature_angle = args.flat_feature_angle;
    args.flat_feature_angle *= 0.5;
    operation(mesh, b_tree, std::array<int, 4>({0, 1, 0, 0}));
    mesh.is_limit_length = true;

    bool is_just_after_update = false;
    bool is_hit_min_edge_length = false;
    std::vector<std::array<double, 2>> quality_queue;
//    int cnt_increase_epsilon = 0;
//    if(args.stage > 1)
//        cnt_increase_epsilon = 1;
    int cnt_increase_epsilon = args.stage - 1;
    for (int it = 0; it < mesh.max_its; it++) {
        double max_energy, avg_energy;
        get_max_avg_energy(mesh, max_energy, avg_energy);
        if (max_energy <= mesh.stop_energy)
            break;

        cout << "//////////////// pass " << it << " ////////////////" << endl;
        operation(mesh, b_tree, ops);
        if(mesh.is_edge_length_achieved) {
            args.flat_feature_angle = flat_feature_angle;
            if (cnt_increase_epsilon > 0 && cnt_increase_epsilon == args.stage - 1) {
//                cout<<mesh.original_epsilon<<" "<<mesh.epsilon<< endl;
                mesh.epsilon += mesh.original_epsilon/args.stage;
                mesh.epsilon_2 = mesh.epsilon * mesh.epsilon;
                cnt_increase_epsilon--;
//                cout << "increase epsilon" << endl;
//                cout<<mesh.original_epsilon<<" "<<mesh.epsilon<< endl;
//                pausee();
            }
        }

        double new_max_energy, new_avg_energy;
        get_max_avg_energy(mesh, new_max_energy, new_avg_energy);
        if (!is_just_after_update) {
            if (max_energy - new_max_energy < 5e-1 && avg_energy - new_avg_energy < 5e-2) {
                is_hit_min_edge_length = update_scaling_field(mesh, new_max_energy) || is_hit_min_edge_length;
                is_just_after_update = true;
                if(cnt_increase_epsilon > 0){
//                    cout<<mesh.original_epsilon<<" "<<mesh.epsilon<< endl;
                    mesh.epsilon += mesh.original_epsilon/args.stage;
                    mesh.epsilon_2 = mesh.epsilon * mesh.epsilon;
                    cnt_increase_epsilon--;
//                    cout << "increase epsilon" << endl;
//                    cout<<mesh.original_epsilon<<" "<<mesh.epsilon<< endl;
//                    pausee();
                }
            }
        } else
            is_just_after_update = false;

        quality_queue.push_back(std::array<double, 2>({new_max_energy, new_avg_energy}));
        if (is_hit_min_edge_length && it > 10) {
//            for(auto& q:quality_queue)
//                cout<<q[0]<<" "<<q[1]<<endl;
            bool is_loop = true;
            for (int i = 0; i < 10; i++) {
//                cout << quality_queue[it - i][0] - quality_queue[it - i - 1][0] << " "
//                     << quality_queue[it - i][1] - quality_queue[it - i - 1][1] << endl;
                if (quality_queue[it - i][0] - quality_queue[it - i - 1][0] < -1e-4
                    || quality_queue[it - i][1] - quality_queue[it - i - 1][1] < -1e-4) {
                    is_loop = false;
                    break;
                }
            }
            if (is_loop)
                break;
        }
    }

    //postprocessing
    cout<<"//////////////// postprocessing ////////////////"<<endl;
    if(args.is_preserving_feature)
        mesh.is_limit_length = false;
    for(int i=0;i<mesh.tri_vertices.size();i++) {
        if (mesh.v_is_removed[i])
            continue;
        mesh.tri_vertices[i].scale = 1;
    }
    operation(mesh, b_tree, std::array<int, 4>({0, 1, 0, 0}));

//    if(args.is_preserving_feature)
//        feature::visualize_features(mesh);
}

void triwild::optimization::operation(MeshData& mesh, GEO::MeshFacetsAABB &b_tree, const std::array<int, 4>& ops) {
    igl::Timer igl_timer;

//    int max_i;
//    double max = 0;
//    for(int i=0;i<mesh.tris.size();i++){
//        if(mesh.t_is_removed[i])
//            continue;
//        if(mesh.t_quality[i]>max){
//            max = mesh.t_quality[i];
//            max_i = i;
//        }
//    }
//
//    cout<<std::setprecision(16)<<"face "<<max_i<<endl;
//    cout<<mesh.tris[max_i][0]<<" "<<mesh.tris[max_i][1]<<" "<<mesh.tris[max_i][2]<<endl;
//    for(int j=0;j<3;j++)
//        cout<<mesh.tri_vertices[mesh.tris[max_i][j]].posf<<endl;
//    for(int j=0;j<3;j++)
//        cout<<mesh.tri_vertices[mesh.tris[max_i][j]].is_freezed<<" ";
//    cout<<endl;
//    for(int j=0;j<3;j++) {
//        if(mesh.tri_vertices[mesh.tris[max_i][j]].feature_infos.size()>0) {
//            int f_id = mesh.tri_vertices[mesh.tris[max_i][j]].feature_infos[0][0];
//            cout << "feature " << f_id << ": curve " << feature::features[f_id]->curve_id << " "
//                 << mesh.tri_vertices[mesh.tris[max_i][j]].feature_infos[0][1] << endl;
//        }
//    }
//    cout<<mesh.tag_feature_es[max_i][0]<<" "<<mesh.tag_feature_es[max_i][1]<<" "<<mesh.tag_feature_es[max_i][2]<<endl;
//    cout<<mesh.is_bbox_es[max_i][0]<<" "<<mesh.is_bbox_es[max_i][1]<<" "<<mesh.is_bbox_es[max_i][2]<<endl;
//    pausee();

    if(args.is_preserving_feature) {
        igl_timer.start();
        cout << "feature vertex snapping..." << endl;
        feature::snap_vertices(mesh);
        cout << "feature vertex snapping done!" << endl;
        cout << "time = " << igl_timer.getElapsedTime() << "s" << endl;
        output_info(mesh);
    }

    for (int i = 0; i < ops[0]; i++) {
        igl_timer.start();
        cout << "edge splitting..." << endl;
        edge_splitting(mesh);
        cout << "edge splitting done!" << endl;
        cout << "time = " << igl_timer.getElapsedTime() << "s" << endl;
        output_info(mesh);
    }

    for (int i = 0; i < ops[1]; i++) {
        igl_timer.start();
        cout << "edge collapsing..." << endl;
        edge_collapsing(mesh, b_tree);
        cout << "edge collapsing done!" << endl;
        cout << "time = " << igl_timer.getElapsedTime() << "s" << endl;
        output_info(mesh);
    }

    for (int i = 0; i < ops[2]; i++) {
        igl_timer.start();
        cout << "edge swapping..." << endl;
        edge_swapping(mesh);
        cout << "edge swapping done!" << endl;
        cout << "time = " << igl_timer.getElapsedTime() << "s" << endl;
        output_info(mesh);
    }

    for (int i = 0; i < ops[3]; i++) {
        igl_timer.start();
        cout << "vertex smoothing..." << endl;
        vertex_smoothing(mesh, b_tree);
        cout << "vertex smoothing done!" << endl;
        cout << "time = " << igl_timer.getElapsedTime() << "s" << endl;
        output_info(mesh);
    }

    double max_dis = check_envelope(mesh, b_tree);
    if (max_dis > 0) {
        cout << "!check_envelope " << max_dis << " " << mesh.epsilon_2 << endl;
        output_boundary(mesh);
        pausee();
    }

//    output_mesh(mesh);
//    pausee();

//    if(args.is_preserving_feature)
//        feature::visualize_features(mesh);
}

#include <geogram/points/kd_tree.h>
bool triwild::optimization::update_scaling_field(MeshData& mesh, double max_energy) {
    cout << "updating sclaing field ..." << endl;
    bool is_hit_min_edge_length = false;

    double radius0 = mesh.ideal_edge_length * 1.8;//increasing the radius would increase the #v in output
//    if(is_hit_min)
//        radius0 *= 2;

    double filter_energy = max_energy / 100 > mesh.stop_energy * 0.9 ? max_energy / 100 : mesh.stop_energy * 0.9;
//    if(filter_energy > 1e3)
//        filter_energy = 1e3;

    if(filter_energy > 5e2) {
        filter_energy = get_mid_energy(mesh);
        if (filter_energy < mesh.stop_energy * 0.9)
            filter_energy = mesh.stop_energy * 0.9;
    }

    cout << "filter_energy = " << filter_energy << endl;
    double recover = 1.5;
    if(args.is_preserving_feature)
        recover = 1.1;
    std::vector<double> scale_multipliers(mesh.tri_vertices.size(), recover);
    double refine_scale = 0.5;
//    double min_refine_scale = mesh.epsilon / mesh.ideal_edge_length;
    double min_refine_scale = mesh.min_scalar;

    const int N = -int(std::log2(min_refine_scale) - 1);
    std::vector<std::vector<int>> v_ids(N, std::vector<int>());
    for (size_t i = 0; i < mesh.tri_vertices.size(); i++) {
        if (mesh.v_is_removed[i])
            continue;

        bool is_refine = false;
        for (int t_id: mesh.tri_vertices[i].conn_tris) {
            if (mesh.t_quality[t_id] > filter_energy)
                is_refine = true;
        }
        if (!is_refine)
            continue;

        int n = -int(std::log2(mesh.tri_vertices[i].scale) - 0.5);
        if (n >= N)
            n = N - 1;
        v_ids[n].push_back(i);
    }

    for (int n = 0; n < N; n++) {
        if (v_ids[n].size() == 0)
            continue;

        double radius = radius0 / std::pow(2, n);

        std::unordered_set<int> is_visited;
        std::queue<int> v_queue;

        std::vector<double> pts;
        pts.reserve(v_ids[n].size() * 2);
        for (int i = 0; i < v_ids[n].size(); i++) {
            pts.push_back(mesh.tri_vertices[v_ids[n][i]].posf[0]);
            pts.push_back(mesh.tri_vertices[v_ids[n][i]].posf[1]);

            v_queue.push(v_ids[n][i]);
            is_visited.insert(v_ids[n][i]);
            scale_multipliers[v_ids[n][i]] = refine_scale;
        }
        // construct the kdtree
        GEO::NearestNeighborSearch_var nnsearch = GEO::NearestNeighborSearch::create(2, "BNN");
        nnsearch->set_points(int(v_ids[n].size()), pts.data());

        while (!v_queue.empty()) {
            int v_id = v_queue.front();
            v_queue.pop();

            for (int t_id:mesh.tri_vertices[v_id].conn_tris) {
                for (int j = 0; j < 3; j++) {
                    if (is_visited.find(mesh.tris[t_id][j]) != is_visited.end())
                        continue;
                    GEO::index_t _;
                    double sq_dist;
                    const double p[2] = {mesh.tri_vertices[mesh.tris[t_id][j]].posf[0],
                                         mesh.tri_vertices[mesh.tris[t_id][j]].posf[1]};
                    nnsearch->get_nearest_neighbors(1, p, &_, &sq_dist);
                    double dis = sqrt(sq_dist);

                    if (dis < radius) {
                        v_queue.push(mesh.tris[t_id][j]);
                        double new_ss = (dis / radius) * (1 - refine_scale) + refine_scale;
                        if (new_ss < scale_multipliers[mesh.tris[t_id][j]])
                            scale_multipliers[mesh.tris[t_id][j]] = new_ss;
                    }
                    is_visited.insert(mesh.tris[t_id][j]);
                }
            }
        }
    }

    // update scalars
    for (int i = 0; i < mesh.tri_vertices.size(); i++) {
        if (mesh.v_is_removed[i])
            continue;
        double new_scale = mesh.tri_vertices[i].scale * scale_multipliers[i];
        if (new_scale > 1)
            mesh.tri_vertices[i].scale = 1;
//        if (new_scale > mesh.tri_vertices[i].max_scale)
//            mesh.tri_vertices[i].scale = mesh.tri_vertices[i].max_scale;
        else if (new_scale < min_refine_scale) {
            is_hit_min_edge_length = true;
            mesh.tri_vertices[i].scale = min_refine_scale;
        } else
            mesh.tri_vertices[i].scale = new_scale;
    }

    cout<<"is_hit_min_edge_length = "<<is_hit_min_edge_length<<endl;
    return is_hit_min_edge_length;
}

void triwild::optimization::erase_outside(MeshData& mesh) {
    std::queue<int> t_queue;
    for (int i = 0; i < mesh.tris.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;
        if (mesh.is_bbox_es[i][0] || mesh.is_bbox_es[i][1] || mesh.is_bbox_es[i][2]) {
            t_queue.push(i);
            mesh.t_is_removed[i] = true;
        }
    }

    while (!t_queue.empty()) {
        int t_id = t_queue.front();
        t_queue.pop();

//        mesh.t_is_removed[t_id] = true;
        for (int j = 0; j < 3; j++) {
            mesh.tri_vertices[mesh.tris[t_id][j]].conn_tris.erase(t_id);
        }
        if (mesh.is_curved) {
            for (int n_id:mesh.tri_nodes[t_id])
                mesh.n_is_removed[n_id] = true;
            mesh.tri_nodes[t_id].clear();
        }

        for (int j = 0; j < 3; j++) {
            if (mesh.is_bbox_es[t_id][j] || mesh.is_boundary_es[t_id][j])
                continue;
            int v1_id = mesh.tris[t_id][(j + 1) % 3];
            int v2_id = mesh.tris[t_id][(j + 2) % 3];
            std::vector<int> tmp = set_intersection(mesh.tri_vertices[v1_id].conn_tris,
                                                    mesh.tri_vertices[v2_id].conn_tris);
            if (tmp.empty())
                continue;
            int n_t_id = set_intersection(mesh.tri_vertices[v1_id].conn_tris,
                                          mesh.tri_vertices[v2_id].conn_tris)[0];
            if(mesh.t_is_removed[n_t_id])
                continue;
            mesh.t_is_removed[n_t_id] = true;
            t_queue.push(n_t_id);
        }
    }

    for (int i = 0; i < mesh.tri_vertices.size(); i++) {
        if (mesh.v_is_removed[i])
            continue;
        if (mesh.tri_vertices[i].conn_tris.empty())
            mesh.v_is_removed[i] = true;
    }
}

void triwild::optimization::erase_holes(MeshData& mesh, const std::string& hole_file, bool is_erase) {
    //get hole points
    std::vector<GEO::vec3> ps;
    std::ifstream f(hole_file);
    if (f.is_open()) {
        double x, y;
        while (f >> x >> y) {
            ps.push_back(GEO::vec3(x, y, 0));
        }
        f.close();
    } else
        return;

    erase_holes(mesh, ps, is_erase);
}

void triwild::optimization::erase_holes(MeshData& mesh, const Eigen::MatrixXd& hole_pts, bool is_erase) {
    if (hole_pts.size() <= 0)
        return;

    assert(hole_pts.cols() == 2);
    std::vector<GEO::vec3> ps; ps.reserve(hole_pts.rows());
    for(int i = 0; i < hole_pts.rows(); ++i)
        ps.emplace_back(hole_pts(i, 0), hole_pts(i, 1), 0);

    erase_holes(mesh, ps, is_erase);
}

void triwild::optimization::erase_holes(MeshData &mesh, const std::vector<GEO::vec3> &ps, bool is_erase){
    //construct aabb tree
    GEO::Mesh f_mesh;
    f_mesh.vertices.clear();
    int cnt = std::count(mesh.v_is_removed.begin(), mesh.v_is_removed.end(), false);
    f_mesh.vertices.create_vertices(cnt);
    std::unordered_map<int, int> map_v_ids;
    cnt = 0;
    for (int i = 0; i < mesh.tri_vertices.size(); i++) {
        if (mesh.v_is_removed[i])
            continue;
        GEO::vec3 &p = f_mesh.vertices.point(cnt);
        p[0] = mesh.tri_vertices[i].posf[0];
        p[1] = mesh.tri_vertices[i].posf[1];
        p[2] = 0;
        map_v_ids[i] = cnt;
        cnt++;
    }

    f_mesh.facets.clear();
    cnt = std::count(mesh.t_is_removed.begin(), mesh.t_is_removed.end(), false);
    f_mesh.facets.create_triangles(cnt);
    cnt = 0;
    std::unordered_map<int, int> inv_map_f_ids;
    for (int i = 0; i < mesh.tris.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;
        f_mesh.facets.set_vertex(cnt, 0, map_v_ids[mesh.tris[i][0]]);
        f_mesh.facets.set_vertex(cnt, 1, map_v_ids[mesh.tris[i][1]]);
        f_mesh.facets.set_vertex(cnt, 2, map_v_ids[mesh.tris[i][2]]);
        inv_map_f_ids[cnt] = i;
        cnt++;
    }
    f_mesh.facets.compute_borders();
    GEO::MeshFacetsAABB f_tree(f_mesh, false);

    //remove
    if(is_erase) {
        GEO::vec3 _p;
        double _d;
        std::queue<int> t_queue;
        for (auto &p:ps) {
            int t_id = inv_map_f_ids[f_tree.nearest_facet(p, _p, _d)];
            mesh.t_is_removed[t_id] = true;
            t_queue.push(t_id);
        }

        while (!t_queue.empty()) {
            int t_id = t_queue.front();
            t_queue.pop();

//        mesh.t_is_removed[t_id] = true;
            for (int j = 0; j < 3; j++) {
                mesh.tri_vertices[mesh.tris[t_id][j]].conn_tris.erase(t_id);
            }
            if (mesh.is_curved) {
                for (int n_id:mesh.tri_nodes[t_id])
                    mesh.n_is_removed[n_id] = true;
                mesh.tri_nodes[t_id].clear();
            }

            for (int j = 0; j < 3; j++) {
                if (mesh.is_bbox_es[t_id][j] || mesh.is_boundary_es[t_id][j])
                    continue;
                int v1_id = mesh.tris[t_id][(j + 1) % 3];
                int v2_id = mesh.tris[t_id][(j + 2) % 3];
                std::vector<int> tmp = set_intersection(mesh.tri_vertices[v1_id].conn_tris,
                                                        mesh.tri_vertices[v2_id].conn_tris);
                if (tmp.empty())
                    continue;
                int n_t_id = set_intersection(mesh.tri_vertices[v1_id].conn_tris,
                                              mesh.tri_vertices[v2_id].conn_tris)[0];
                if (mesh.t_is_removed[n_t_id])
                    continue;
                mesh.t_is_removed[n_t_id] = true;
                t_queue.push(n_t_id);
            }
        }
    }
    else {
        std::vector<bool> is_kept(mesh.tris.size(), false);

        GEO::vec3 _p;
        double _d;
        std::queue<int> t_queue;
        for (auto &p:ps) {
            int t_id = inv_map_f_ids[f_tree.nearest_facet(p, _p, _d)];
            t_queue.push(t_id);
            is_kept[t_id] = true;
        }

        while (!t_queue.empty()) {
            int t_id = t_queue.front();
            t_queue.pop();

            for (int j = 0; j < 3; j++) {
                if (mesh.is_bbox_es[t_id][j] || mesh.is_boundary_es[t_id][j])
                    continue;
                int v1_id = mesh.tris[t_id][(j + 1) % 3];
                int v2_id = mesh.tris[t_id][(j + 2) % 3];
                std::vector<int> tmp = set_intersection(mesh.tri_vertices[v1_id].conn_tris,
                                                        mesh.tri_vertices[v2_id].conn_tris);
                if (tmp.size() < 2)
                    continue;

                int n_t_id = tmp[0] == t_id ? tmp[1] : tmp[0];
                if(is_kept[n_t_id])
                    continue;
                is_kept[n_t_id] = true;
                t_queue.push(n_t_id);
            }
        }

        for (int i = 0; i < mesh.tris.size(); i++) {
            if (mesh.t_is_removed[i])
                continue;
            if (!is_kept[i]) {
                mesh.t_is_removed[i] = true;
                for (int j = 0; j < 3; j++) {
                    mesh.tri_vertices[mesh.tris[i][j]].conn_tris.erase(i);
                }
                if (mesh.is_curved) {
                    for (int n_id:mesh.tri_nodes[i])
                        mesh.n_is_removed[n_id] = true;
                    mesh.tri_nodes[i].clear();
                }
            }
        }
    }

    for (int i = 0; i < mesh.tri_vertices.size(); i++) {
        if (mesh.v_is_removed[i])
            continue;
        if (mesh.tri_vertices[i].conn_tris.empty())
            mesh.v_is_removed[i] = true;
    }
}

#include <igl/readSTL.h>
#include <igl/writeSTL.h>
#include "../../extern/pymesh/MshSaver.h"
void triwild::optimization::output_mesh(MeshData& mesh) {
    cout << "Writing output mesh to file " << args.output << "..." << endl;
//    Eigen::MatrixXd V;
//    Eigen::MatrixXi F;
//
//    std::unordered_map<int, int> map_v_ids;
//    int cnt = 0;
//    int cnt_unrounded = 0;
//    for (size_t i = 0; i < mesh.tri_vertices.size(); i++) {
//        if (mesh.v_is_removed[i])
//            continue;
//        map_v_ids[i] = cnt++;
//
//        if (!mesh.tri_vertices[i].is_rounded)
//            cnt_unrounded++;
//    }
//    V.resize(cnt, 3);
//    cnt = 0;
//    for (size_t i = 0; i < mesh.tri_vertices.size(); i++) {
//        if (mesh.v_is_removed[i])
//            continue;
//        V(cnt, 0) = mesh.tri_vertices[i].posf[0];
//        V(cnt, 1) = mesh.tri_vertices[i].posf[1];
//        V(cnt, 2) = 0;
//        cnt++;
//    }
//
//    cnt = std::count(mesh.t_is_removed.begin(), mesh.t_is_removed.end(), false);
//    F.resize(cnt, 3);
//    cnt = 0;
//    for (size_t i = 0; i < mesh.tris.size(); i++) {
//        if (mesh.t_is_removed[i])
//            continue;
//        for (int j = 0; j < 3; j++)
//            F(cnt, j) = map_v_ids[mesh.tris[i][j]];
//        cnt++;
//    }
//    igl::writeSTL(args.output, V, F);

    if(args.is_preserving_feature)
        write_msh(mesh, args.output + "_curved.msh");
    if(!args.is_preserving_feature || args.output_linear) {
        Eigen::VectorXd V;
        Eigen::VectorXi F;

        std::unordered_map<int, int> map_v_ids;
        int cnt = 0;
        int cnt_unrounded = 0;
        for (size_t i = 0; i < mesh.tri_vertices.size(); i++) {
            if (mesh.v_is_removed[i])
                continue;
            map_v_ids[i] = cnt++;

            if (!mesh.tri_vertices[i].is_rounded)
                cnt_unrounded++;
        }
        V.resize(cnt * 3);
        cnt = 0;
        for (size_t i = 0; i < mesh.tri_vertices.size(); i++) {
            if (mesh.v_is_removed[i])
                continue;
            V(cnt * 3 + 0) = mesh.tri_vertices[i].posf[0];
            V(cnt * 3 + 1) = mesh.tri_vertices[i].posf[1];
            V(cnt * 3 + 2) = 0;
            cnt++;
        }
        cnt = std::count(mesh.t_is_removed.begin(), mesh.t_is_removed.end(), false);
        F.resize(cnt * 3);
        cnt = 0;
        for (size_t i = 0; i < mesh.tris.size(); i++) {
            if (mesh.t_is_removed[i])
                continue;
            for (int j = 0; j < 3; j++)
                F(cnt * 3 + j) = map_v_ids[mesh.tris[i][j]];
            cnt++;
        }
        PyMesh::MshSaver mSaver(args.output + "_linear.msh", true);
        mSaver.save_mesh(V, F, 3, mSaver.TRI);
    }

//    output_boundary(mesh);
}

void triwild::optimization::output_stats(const MeshData& mesh, std::ofstream& f) {
    ////#v
    int v_cnt = 0;
    int v_unrounded_cnt = 0;
    for (int i = 0; i < mesh.tri_vertices.size(); i++) {
        if (mesh.v_is_removed[i])
            continue;
        v_cnt++;
        if (!mesh.tri_vertices[i].is_rounded)
            v_unrounded_cnt++;
    }
    f << "v_cnt " << v_cnt << endl;
    f << "v_unrounded_cnt " << v_unrounded_cnt << endl;

    ////#f
    int f_cnt = 0;
    int f_inverted_cnt = 0;
    double min_area, max_area;
    double avg_area = 0;
    for (int i = 0; i < mesh.tris.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;

        double area = t_area(mesh, i);
        if (f_cnt == 0) {
            min_area = area;
            max_area = area;
        } else {
            if (area < min_area)
                min_area = area;
            if (area > max_area)
                max_area = area;
        }
        avg_area += area;

        if (!is_valid_inversion(mesh.tri_vertices[mesh.tris[i][0]].posf, mesh.tri_vertices[mesh.tris[i][1]].posf,
                                mesh.tri_vertices[mesh.tris[i][2]].posf,
                                mesh.tri_vertices[mesh.tris[i][0]].pos, mesh.tri_vertices[mesh.tris[i][1]].pos,
                                mesh.tri_vertices[mesh.tris[i][2]].pos))
            f_inverted_cnt++;

        f_cnt++;
    }
    if (f_cnt > 0)
        avg_area /= f_cnt;
    f << "f_cnt " << f_cnt << endl;
    f << "f_inverted_cnt " << f_inverted_cnt << endl;

    ////element size
    f << "min_area " << min_area << endl;
    f << "max_area " << max_area << endl;
    f << "avg_area " << avg_area << endl;

    ////energy
    double avg, max;
    get_max_avg_energy(mesh, max, avg);
    f << "max_energy " << max << endl;
    f << "avg_energy " << avg << endl;

    ////bbox
    f << "bbox_min " << Eigen::RowVector2d(args.box_min) << endl;
    f << "bbox_max " << Eigen::RowVector2d(args.box_max) << endl;
}

double triwild::optimization::t_area(const MeshData& mesh, int t_id) {
    auto &pf1 = mesh.tri_vertices[mesh.tris[t_id][0]].posf;
    auto &pf2 = mesh.tri_vertices[mesh.tris[t_id][1]].posf;
    auto &pf3 = mesh.tri_vertices[mesh.tris[t_id][2]].posf;

    return ((pf2.x - pf1.x) * (pf3.y - pf2.y) - (pf3.x - pf2.x) * (pf2.y - pf1.y)) * 0.5;
}

bool triwild::optimization::round_a_vertex(MeshData& mesh, int v_id) {
    Point_2 p(mesh.tri_vertices[v_id].posf[0], mesh.tri_vertices[v_id].posf[1]);
    for (int t_id: mesh.tri_vertices[v_id].conn_tris) {
        int j = std::find(mesh.tris[t_id].begin(), mesh.tris[t_id].end(), v_id) - mesh.tris[t_id].begin();
        if (!is_valid_inversion(mesh.tri_vertices[v_id].posf,
                                mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].posf,
                                mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].posf,
                                p,
                                mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].pos,
                                mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].pos))
            return false;
    }
    mesh.tri_vertices[v_id].is_rounded = true;
    mesh.tri_vertices[v_id].pos = p;
    return true;
}

bool triwild::optimization::is_valid_quality(Point_2f& p1, Point_2f& p2, Point_2f& p3, double old_max_energy, double& new_energy,
        bool accept_equal) {
    new_energy = tri_energy(p1, p2, p3);
    if(accept_equal) {
        if (new_energy > old_max_energy)
            return false;
        return true;
    } else {
        if (new_energy >= old_max_energy)
            return false;
        return true;
    }
}

bool triwild::optimization::is_valid_inversion(const Point_2f& pf1, const Point_2f& pf2, const Point_2f& pf3,
        const Point_2& p1, const Point_2& p2, const Point_2& p3) {
    double detf = ((pf2.x - pf1.x) * (pf3.y - pf2.y) - (pf3.x - pf2.x) * (pf2.y - pf1.y));
    if (std::abs(detf) > 1e-4) {
        if (detf > 0)
            return true;
        else
            return false;
    }

    if (((p2.x - p1.x) * (p3.y - p2.y) - (p3.x - p2.x) * (p2.y - p1.y)) > 0)
        return true;
    return false;
}

bool triwild::optimization::is_valid_envelop(MeshData& mesh, GEO::MeshFacetsAABB &b_tree, int v_id, Point_2f& p){
    std::unordered_set<int> n_v_ids;
    for(int t_id:mesh.tri_vertices[v_id].conn_tris) {
        for (int j = 0; j < 3; j++)
            if (mesh.tris[t_id][j] != v_id && mesh.is_boundary_es[t_id][j] && mesh.tag_feature_es[t_id][j] < 0) {
                int tmp_v_id = mesh.tris[t_id][(j + 1) % 3] == v_id ?
                               mesh.tris[t_id][(j + 2) % 3] : mesh.tris[t_id][(j + 1) % 3];
                n_v_ids.insert(tmp_v_id);
            }
    }
    if(n_v_ids.empty())
        return true;

    GEO::vec3 p0(p[0], p[1], 0);
    GEO::vec3 nearest_point;
    double sq_dist;
    GEO::index_t prev_facet = b_tree.nearest_facet(p0, nearest_point, sq_dist);
    if (sq_dist > mesh.epsilon_2)
        return false;
    for(int n_v_id:n_v_ids) {
        int N = (mesh.tri_vertices[n_v_id].posf - p).length() / mesh.dd + 1;
        for (double n = 1; n <= N; n++) {
            double k = n / N;
            GEO::vec3 v(mesh.tri_vertices[n_v_id].posf[0] * k + p[0] * (1 - k),
                        mesh.tri_vertices[n_v_id].posf[1] * k + p[1] * (1 - k),
                        0);
            sq_dist = v.distance2(nearest_point);
            if(sq_dist <= mesh.epsilon_2)
                continue;
            b_tree.nearest_facet_with_hint(v, prev_facet, nearest_point, sq_dist);
            if (sq_dist > mesh.epsilon_2)
                return false;
        }
    }

//    //sampling
//    std::vector<GEO::vec3> vs;
//    vs.push_back(GEO::vec3(p[0], p[1], 0));
//    for(int n_v_id:n_v_ids) {
//        Point_2f &pn = mesh.tri_vertices[n_v_id].posf;
//        int N = std::sqrt(edge_length_2(pn, p)) / mesh.dd + 1;
//        vs.reserve(vs.size() + N);
//        for (double n = 1; n <= N; n++) {
//            Point_2f tmp = pn * (n / N) + p * (N - n) / N;
//            vs.push_back(GEO::vec3(tmp[0], tmp[1], 0));
//        }
//    }
//
//    GEO::vec3 current_point = vs[0];
//    GEO::vec3 nearest_point;
//    double sq_dist;
//    GEO::index_t prev_facet = b_tree.nearest_facet(current_point, nearest_point, sq_dist);
//    if (sq_dist > mesh.epsilon_2)
//        return false;
//    for (auto &v : vs) {
//        sq_dist = v.distance2(nearest_point);
//        b_tree.nearest_facet_with_hint(v, prev_facet, nearest_point, sq_dist);
//        if (sq_dist > mesh.epsilon_2)
//            return false;
//    }

    return true;
}

bool triwild::optimization::is_edge_out_of_envelop(MeshData& mesh, GEO::MeshFacetsAABB &b_tree, Point_2f& p1, Point_2f& p2) {
    //sampling
    std::vector<GEO::vec3> vs;
    int N = std::sqrt(edge_length_2(p1, p2)) / mesh.dd + 1;
    vs.reserve(vs.size() + N);
    for (double n = 0; n <= N; n++) {
        Point_2f tmp = p1 * (n / N) + p2 * (N - n) / N;
        vs.push_back(GEO::vec3(tmp[0], tmp[1], 0));
    }

    //checking
    GEO::vec3 current_point = vs[0];
    GEO::vec3 nearest_point;
    double sq_dist;
    GEO::index_t prev_facet = b_tree.nearest_facet(current_point, nearest_point, sq_dist);
    if (sq_dist > mesh.epsilon_2)
        return false;
    for (auto &v : vs) {
        sq_dist = v.distance2(nearest_point);
        b_tree.nearest_facet_with_hint(v, prev_facet, nearest_point, sq_dist);
        if (sq_dist > mesh.epsilon_2)
            return true;
    }

    return false;
}

bool triwild::optimization::is_valid_feature_edge_length(MeshData& mesh, int v_id, Point_2f& p) {
    for (int t_id:mesh.tri_vertices[v_id].conn_tris) {
        int j = std::find(mesh.tris[t_id].begin(), mesh.tris[t_id].end(), v_id) - mesh.tris[t_id].begin();
        if (mesh.tag_feature_es[t_id][(j + 1) % 3] >= 0) {
            double l = std::sqrt(edge_length_2(mesh, v_id, mesh.tris[t_id][(j + 2) % 3]));
            if (l > mesh.ideal_edge_length *
                    (mesh.tri_vertices[v_id].max_scale + mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].max_scale) / 2)
                return false;
        }
        if (mesh.tag_feature_es[t_id][(j + 2) % 3] >= 0) {
            double l = std::sqrt(edge_length_2(mesh, v_id, mesh.tris[t_id][(j + 1) % 3]));
            if (l > mesh.ideal_edge_length *
                    (mesh.tri_vertices[v_id].max_scale + mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].max_scale) / 2)
                return false;
        }
    }
    return true;
}

bool triwild::optimization::is_valid_feature_edge_close(MeshData& mesh, int v_id, const Point_2f& p, double t) {
    std::vector<std::array<int, 2>> feature_n_v_ids;
    for (int t_id:mesh.tri_vertices[v_id].conn_tris) {
        int j = std::find(mesh.tris[t_id].begin(), mesh.tris[t_id].end(), v_id) - mesh.tris[t_id].begin();
        if (mesh.tag_feature_es[t_id][(j + 1) % 3] >= 0)
            feature_n_v_ids.push_back({{mesh.tris[t_id][(j + 2) % 3], mesh.tag_feature_es[t_id][(j + 1) % 3]}});
        if (mesh.tag_feature_es[t_id][(j + 2) % 3] >= 0)
            feature_n_v_ids.push_back({{mesh.tris[t_id][(j + 1) % 3], mesh.tag_feature_es[t_id][(j + 2) % 3]}});
    }
    vector_unique(feature_n_v_ids);

    if (feature_n_v_ids.size() != 2) {
        cout << "feature_n_v_ids.size() != 2" << endl;
        cout << feature_n_v_ids.size() << endl;
        pausee();
        return false;
    }

    for (int i = 0; i < 2; i++) {
        Point_2f &p1 = mesh.tri_vertices[feature_n_v_ids[i][0]].posf;
        if (p == p1)
            continue;

        int feature_id = feature_n_v_ids[i][1];
        double t1 = mesh.tri_vertices[feature_n_v_ids[i][0]].get_t(feature_id);
        double angle;
        if (t < t1)
            angle = feature::features[feature_id]->how_curve(t, t1);
        else
            angle = feature::features[feature_id]->how_curve(t1, t);

        if (angle > args.flat_feature_angle)//not flat enough
            return false;

//        Vector_2f v1 = feature::features[feature_id]->eval_first_derivative(t);
//        double t1 = mesh.tri_vertices[feature_n_v_ids[i][0]].get_t(feature_id);
//        Vector_2f v2 = feature::features[feature_id]->eval_first_derivative(t1);
//        double cos_a = (v1.dot(v2)) / (v1.length() * v2.length());
////        if (cos_a < -1)
////            cos_a = -1;
////        if (cos_a > 1)
////            cos_a = 1;
////        cos_a = std::abs(cos_a);
//
//        if (cos_a < 0.99)//~60 degree
//            return false;

    }

    return true;
}

bool triwild::optimization::is_bbox_edge(const MeshData& mesh, int v1_id, int v2_id) {
    if(!mesh.tri_vertices[v1_id].is_on_bbox || !mesh.tri_vertices[v2_id].is_on_bbox)
        return false;

    std::vector<int> n_t_ids = set_intersection(mesh.tri_vertices[v1_id].conn_tris,
                                                mesh.tri_vertices[v2_id].conn_tris);
    for (int t_id:n_t_ids) {
        for (int j = 0; j < 3; j++) {
            if (mesh.tris[t_id][j] != v1_id && mesh.tris[t_id][j] != v2_id && mesh.is_bbox_es[t_id][j])
                return true;
        }
    }
    return false;
}

bool triwild::optimization::is_boundary_edge(MeshData& mesh, int v1_id, int v2_id) {
    if(!mesh.tri_vertices[v1_id].is_on_boundary || !mesh.tri_vertices[v2_id].is_on_boundary)
        return false;

    std::vector<int> n_t_ids = set_intersection(mesh.tri_vertices[v1_id].conn_tris,
                                                mesh.tri_vertices[v2_id].conn_tris);
    for (int t_id:n_t_ids) {
        for (int j = 0; j < 3; j++) {
            if (mesh.tris[t_id][j] != v1_id && mesh.tris[t_id][j] != v2_id && mesh.is_boundary_es[t_id][j])
                return true;
        }
    }
    return false;
}

bool triwild::optimization::is_feature_edge(MeshData& mesh, int v1_id, int v2_id) {
    if(mesh.tri_vertices[v1_id].feature_infos.size()==0 || mesh.tri_vertices[v2_id].feature_infos.size()==0)
//    if(!mesh.tri_vertices[v1_id].is_on_feature || !mesh.tri_vertices[v2_id].is_on_feature)
        return false;

    std::vector<int> n_t_ids = set_intersection(mesh.tri_vertices[v1_id].conn_tris,
                                                mesh.tri_vertices[v2_id].conn_tris);
    for (int t_id:n_t_ids) {
        for (int j = 0; j < 3; j++) {
            if (mesh.tris[t_id][j] != v1_id && mesh.tris[t_id][j] != v2_id && mesh.tag_feature_es[t_id][j] >= 0)
                return true;
        }
    }
    return false;
}

int triwild::optimization::get_feature_edge_tag(const MeshData& mesh, int v1_id, int v2_id) {
    if (mesh.tri_vertices[v1_id].feature_infos.size() == 0 || mesh.tri_vertices[v2_id].feature_infos.size() == 0)
        return -1;

    for (int t_id:mesh.tri_vertices[v1_id].conn_tris) {
        int j2 = std::find(mesh.tris[t_id].begin(), mesh.tris[t_id].end(), v2_id) - mesh.tris[t_id].begin();
        if (j2 < 3) {
            int j = mesh.tris[t_id][(j2 + 1) % 3] == v1_id ? (j2 + 2) % 3 : (j2 + 1) % 3;
            if (mesh.tag_feature_es[t_id][j] >= 0)
                return mesh.tag_feature_es[t_id][j];
        }
    }

//    std::vector<int> n_t_ids = set_intersection(mesh.tri_vertices[v1_id].conn_tris,
//                                                mesh.tri_vertices[v2_id].conn_tris);
//    for (int t_id:n_t_ids) {
//        for (int j = 0; j < 3; j++) {
//            if (mesh.tris[t_id][j] != v1_id && mesh.tris[t_id][j] != v2_id && mesh.tag_feature_es[t_id][j] >= 0)
//                return mesh.tag_feature_es[t_id][j];
//        }
//    }
    return -1;
}

int triwild::optimization::get_secondary_feature_edge_tag(const MeshData& mesh, int v1_id, int v2_id)
{
    for (int t_id:mesh.tri_vertices[v1_id].conn_tris) {
        int j2 = std::find(mesh.tris[t_id].begin(), mesh.tris[t_id].end(), v2_id) - mesh.tris[t_id].begin();
        if (j2 < 3) {
            int j = mesh.tris[t_id][(j2 + 1) % 3] == v1_id ? (j2 + 2) % 3 : (j2 + 1) % 3;
            if (mesh.tag_secondary_feature_es[t_id][j] >= 0)
                return mesh.tag_secondary_feature_es[t_id][j];
        }
    }

    return -1;
}

bool triwild::optimization::is_valid_edge(MeshData& mesh, int v1_id, int v2_id) {
    if (mesh.v_is_removed[v1_id] || mesh.v_is_removed[v2_id])
        return false;

    std::vector<int> n_t_ids = set_intersection(mesh.tri_vertices[v1_id].conn_tris,
                                                mesh.tri_vertices[v2_id].conn_tris);
    if (n_t_ids.size() > 0)
        return true;
    return false;
}

bool triwild::optimization::is_isolated_vertex(MeshData& mesh, int v_id) {
    for (int t_id:mesh.tri_vertices[v_id].conn_tris) {
        for (int j = 0; j < 3; j++)
            if (mesh.is_boundary_es[t_id][j])
                return false;
    }
    return true;
}

void triwild::optimization::get_all_edges(const MeshData& mesh, std::vector<std::array<int, 2>>& edges) {
    edges.reserve(mesh.tris.size() * 3);
    for (size_t i = 0; i < mesh.tris.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;
        const auto &t = mesh.tris[i];
        for (int j = 0; j < 3; j++) {
            if (t[j] < t[(j + 1) % 3])
                edges.push_back(std::array<int, 2>({t[j], t[(j + 1) % 3]}));
            else
                edges.push_back(std::array<int, 2>({t[(j + 1) % 3], t[j]}));
        }
    }
    vector_unique(edges);
}

double triwild::optimization::edge_length_2(const MeshData& mesh, int v1_id, int v2_id) {
    return (mesh.tri_vertices[v1_id].posf - mesh.tri_vertices[v2_id].posf).length_2();
}

double triwild::optimization::edge_length_2(const Point_2f& p1, const Point_2f& p2) {
    return (p1 - p2).length_2();
}

Point_2 triwild::optimization::mid_point(const Point_2& p1, const Point_2& p2) {
    return (p1 + p2) / 2;
}

Point_2f triwild::optimization::mid_point(const Point_2f& p1, const Point_2f& p2) {
    return (p1 + p2) / 2;
}

double triwild::optimization::tri_energy(Point_2f& p1, Point_2f& p2, Point_2f& p3) {
    std::array<double, 6> T = {p1[0], p1[1], p2[0], p2[1], p3[0], p3[1]};
    double energy = AMIPS_energy(T);
    if (energy > args.MAX_ENERGY || std::isnan(energy) || std::isinf(energy) || energy <= 0)
        energy = args.MAX_ENERGY;

//    if (energy < 2 - 1e-8) {
//        cout << energy << endl;
//        cout<<"J_det = "<< (-T[0] + T[2])*(-0.577350269189626*T[1] - 0.577350269189626*T[3] + 1.15470053837925*T[5]) - (-T[1] + T[3])*(-0.577350269189626*T[0] - 0.577350269189626*T[2] + 1.15470053837925*T[4])<<endl;
//        cout<<"J_trace = "<<(-0.666666666666667*T[0] - 0.666666666666667*T[2] + 1.33333333333333*T[4])*T[4] + (-0.666666666666667*T[0] + 1.33333333333333*T[2] - 0.666666666666667*T[4])*T[2] + (1.33333333333333*T[0] - 0.666666666666667*T[2] - 0.666666666666667*T[4])*T[0] + (-0.666666666666667*T[1] - 0.666666666666667*T[3] + 1.33333333333333*T[5])*T[5] + (-0.666666666666667*T[1] + 1.33333333333333*T[3] - 0.666666666666667*T[5])*T[3] + (1.33333333333333*T[1] - 0.666666666666667*T[3] - 0.666666666666667*T[5])*T[1]<<endl;
//        pausee();
//    }

    return energy;
}

void triwild::optimization::get_max_avg_energy(const MeshData& mesh, double& max_energy, double& avg_energy) {
    max_energy = avg_energy = 0;
    int cnt = 0;
    for (int i = 0; i < mesh.t_quality.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;
        if (mesh.t_quality[i] > max_energy)
            max_energy = mesh.t_quality[i];
        avg_energy += mesh.t_quality[i];
        cnt++;
    }
    avg_energy /= cnt;
}

double triwild::optimization::get_mid_energy(const MeshData& mesh) {
    std::vector<double> tmp;
    for (int i = 0; i < mesh.t_quality.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;
        tmp.push_back(mesh.t_quality[i]);
    }
    std::sort(tmp.begin(), tmp.end());
    return tmp[tmp.size() / 2];
}

bool triwild::optimization::get_opp_r_v_id(const MeshData& mesh, int t_id, int v1_id, int v2_id) {
    for (int i = 0; i < 3; i++) {
        if (mesh.tris[t_id][i] != v1_id && mesh.tris[t_id][i] != v2_id)
            return i;
    }
    return -1;
}

double triwild::optimization::check_envelope(MeshData& mesh, GEO::MeshFacetsAABB &b_tree, double s) {
    return -1;

    if(args.mute_log)
        return -1;

    if(mesh.is_curved)
        return -1;

    double max_dis = -1;

    std::vector<std::array<int, 2>> b_es;
    for (int i = 0; i < mesh.tris.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;
        for (int j = 0; j < 3; j++) {
            if (mesh.is_boundary_es[i][j] && mesh.tag_feature_es[i][j] < 0) {
                std::array<int, 2> e = {mesh.tris[i][(j + 1) % 3], mesh.tris[i][(j + 2) % 3]};
                if (e[0] < e[1])
                    b_es.push_back(e);
                else {
                    std::swap(e[0], e[1]);
                    b_es.push_back(e);
                }
            }
        }
    }
    vector_unique(b_es);
    if (b_es.size() == 0)
        return max_dis;

    //sampling
    std::vector<GEO::vec3> vs;
//    cout << b_es.size() << endl;
    for (auto &e: b_es) {
        int N = std::sqrt(edge_length_2(mesh.tri_vertices[e[0]].posf, mesh.tri_vertices[e[1]].posf)) / mesh.dd + 1;
        vs.reserve(vs.size() + N);
        for (double n = 0; n <= N; n++) {
            Point_2f tmp = mesh.tri_vertices[e[0]].posf * (n / N) + mesh.tri_vertices[e[1]].posf * (N - n) / N;
            vs.push_back(GEO::vec3(tmp[0], tmp[1], 0));
        }
    }

    //check
    double real_envelope_2 = mesh.epsilon_2 * s * s;
    GEO::vec3 current_point = vs[0];
    GEO::vec3 nearest_point;
    double sq_dist;
    GEO::index_t prev_facet = b_tree.nearest_facet(current_point, nearest_point, sq_dist);
    if (sq_dist > mesh.epsilon_2)
        return false;
    for (auto &v : vs) {
        sq_dist = v.distance2(nearest_point);
        b_tree.nearest_facet_with_hint(v, prev_facet, nearest_point, sq_dist);
        if (sq_dist > real_envelope_2) {
            if (sq_dist > max_dis)
                max_dis = sq_dist;
        }
    }

    return max_dis;
}

void triwild::optimization::check(MeshData& mesh) {
    if(args.mute_log)
        return;

    //check size consistency
    if (mesh.tri_vertices.size() != mesh.v_is_removed.size()) {
        cout << "mesh.tri_vertices.size()!=mesh.v_is_removed.size()" << endl;
        pausee();
    }
    if (mesh.tris.size() != mesh.t_quality.size() || mesh.tris.size() != mesh.t_is_removed.size()
        || mesh.tris.size() != mesh.is_boundary_es.size() || mesh.tris.size() != mesh.is_bbox_es.size()
        || mesh.tris.size() != mesh.tag_feature_es.size()) {
        cout << "mesh.tris.size()!=xxx" << endl;
        pausee();
    }

    //check tags
    for (int i = 0; i < mesh.tris.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;
        for (int j = 0; j < 3; j++) {
            if (mesh.is_boundary_es[i][j]) {
                if (!(mesh.tri_vertices[mesh.tris[i][(j + 1) % 3]].is_on_boundary
                      && mesh.tri_vertices[mesh.tris[i][(j + 2) % 3]].is_on_boundary)) {
                    cout << "boundary tag error" << endl;
                    cout << "tri " << i << ": " << mesh.tris[i][0] << " " << mesh.tris[i][1] << " " << mesh.tris[i][2]
                         << endl;
                    cout << "tri " << i << ": " << mesh.is_boundary_es[i][0] << " " << mesh.is_boundary_es[i][1] << " "
                         << mesh.is_boundary_es[i][2] << endl;
                    cout << "tri " << i << ": " << mesh.tri_vertices[mesh.tris[i][0]].is_on_boundary << " "
                         << mesh.tri_vertices[mesh.tris[i][1]].is_on_boundary << " "
                         << mesh.tri_vertices[mesh.tris[i][2]].is_on_boundary << endl;
                    pausee();
                }
            }
            if (mesh.is_bbox_es[i][j]) {
                if (!(mesh.tri_vertices[mesh.tris[i][(j + 1) % 3]].is_on_bbox
                      && mesh.tri_vertices[mesh.tris[i][(j + 2) % 3]].is_on_bbox)) {
                    cout << "bbox tag error" << endl;
                    cout << mesh.is_bbox_es[i][j] << endl;
                    cout << mesh.tri_vertices[mesh.tris[i][(j + 1) % 3]].is_on_bbox << endl;
                    cout << mesh.tri_vertices[mesh.tris[i][(j + 2) % 3]].is_on_bbox << endl;
                    cout << "tri " << i << ": " << mesh.tris[i][0] << " " << mesh.tris[i][1] << " "
                         << mesh.tris[i][2] << endl;
                    cout << "tri " << i << ": " << mesh.is_bbox_es[i][0] << " " << mesh.is_bbox_es[i][1] << " "
                         << mesh.is_bbox_es[i][2] << endl;
                    cout << mesh.tri_vertices[mesh.tris[i][0]].is_on_bbox
                         << mesh.tri_vertices[mesh.tris[i][1]].is_on_bbox
                         << mesh.tri_vertices[mesh.tris[i][2]].is_on_bbox
                         << endl;
                    pausee();
                }
            }
            if (mesh.tag_feature_es[i][j] >= 0 &&
//                (!mesh.tri_vertices[mesh.tris[i][(j + 1) % 3]].is_on_feature
//                 || !mesh.tri_vertices[mesh.tris[i][(j + 2) % 3]].is_on_feature)) {
                (mesh.tri_vertices[mesh.tris[i][(j + 1) % 3]].feature_infos.size() == 0
                 || mesh.tri_vertices[mesh.tris[i][(j + 2) % 3]].feature_infos.size() == 0)) {
                cout << "feature tag error" << endl;
                cout<<"mesh.tag_feature_es["<<i<<"]["<<j<<"] = "<<mesh.tag_feature_es[i][j]<<endl;
                cout<<"mesh.tris["<<i<<"] = "<<mesh.tris[i][0]<<" "<<mesh.tris[i][1]<<" "<<mesh.tris[i][2]<<endl;
                for(int j=0;j<3;j++){
                    cout<<mesh.tris[i][j]<<endl;
                    for(auto& info:mesh.tri_vertices[mesh.tris[i][j]].feature_infos)
                        cout<<info[0]<<" "<<info[1]<<endl;
                }
                pausee();
            }
        }
    }

    //check rounding
    for (int i = 0; i < mesh.tri_vertices.size(); i++) {
        if (mesh.v_is_removed[i])
            continue;
        if (mesh.tri_vertices[i].posf.isnan()) {
            cout << "nan v" << endl;
            pausee();
        }
        if (mesh.tri_vertices[i].posf.isinf()) {
            cout << "inf v" << endl;
            pausee();
        }
        if (mesh.tri_vertices[i].is_rounded) {
            Point_2 p(mesh.tri_vertices[i].posf[0], mesh.tri_vertices[i].posf[1]);
            if (p != mesh.tri_vertices[i].pos) {
                cout << "rounding inconsistant! p!=mesh.tri_vertices[i].pos" << endl;
                cout<<"pos = "<<mesh.tri_vertices[i].pos<<endl;
                cout<<"posf = "<<mesh.tri_vertices[i].posf<<endl;
                pausee();
            }
        } else {
            Point_2f pf(mesh.tri_vertices[i].pos[0].to_double(), mesh.tri_vertices[i].pos[1].to_double());
            if (pf != mesh.tri_vertices[i].posf) {
                cout << "rounding inconsistant! pf!=mesh.tri_vertices[i].posf" << endl;
                pausee();
            }
        }

        //check conn_tris
        if(mesh.tri_vertices[i].conn_tris.size()==0){
            cout<<"v "<<i<<" conn_tris.size()==0"<<endl;
            cout<<mesh.tri_vertices[i].posf<<endl;
            pausee();
        }
        for(int t_id:mesh.tri_vertices[i].conn_tris){
            if(std::find(mesh.tris[t_id].begin(), mesh.tris[t_id].end(), i) == mesh.tris[t_id].end()){
                cout<<"v_id = "<<i<<endl;
                for(int t_id:mesh.tri_vertices[i].conn_tris)
                    cout<<t_id<<" ";
                cout<<endl;
                cout<<"tri "<<t_id<<": "<<mesh.tris[t_id][0]<<" "<<mesh.tris[t_id][1]<<" "<<mesh.tris[t_id][2]<<endl;
                pausee();
            }
        }

        //cannot on both bbox and boundary
        if(mesh.tri_vertices[i].is_on_bbox && mesh.tri_vertices[i].is_on_boundary){
            cout<<"v "<<i<<" is_on_bbox && is_on_boundary"<<endl;
            pausee();
        }
        if(mesh.tri_vertices[i].is_on_bbox && mesh.tri_vertices[i].feature_infos.size() > 0){
            cout<<"v "<<i<<" is_on_bbox && is_on_feature"<<endl;
            pausee();
        }
        if(mesh.tri_vertices[i].is_on_point && !mesh.tri_vertices[i].is_on_boundary){
            cout<<"v "<<i<<" is_on_point && !is_on_boundary"<<endl;
            pausee();
        }
    }

    for (int i = 0; i < mesh.tris.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;

        //check inversion
        if (!is_valid_inversion(mesh.tri_vertices[mesh.tris[i][0]].posf, mesh.tri_vertices[mesh.tris[i][1]].posf,
                                mesh.tri_vertices[mesh.tris[i][2]].posf,
                                mesh.tri_vertices[mesh.tris[i][0]].pos, mesh.tri_vertices[mesh.tris[i][1]].pos,
                                mesh.tri_vertices[mesh.tris[i][2]].pos)) {
            cout << "!is_valid_inversion" << endl;
            cout << "tri " << i << ": " << mesh.tris[i][0] << " " << mesh.tris[i][1] << " " << mesh.tris[i][2] << endl;
            pausee();
		exit(0);
        }

        //check quality computation
        double q = tri_energy(mesh.tri_vertices[mesh.tris[i][0]].posf, mesh.tri_vertices[mesh.tris[i][1]].posf,
                              mesh.tri_vertices[mesh.tris[i][2]].posf);
        if (false && std::abs(q - mesh.t_quality[i]) > 1e-6) {
            cout << "q!=mesh.t_quality[i]" << endl;
            cout << q << " " << mesh.t_quality[i] << endl;
            cout << q - mesh.t_quality[i] << endl;
            cout << "tri " << i << ": " << mesh.tris[i][0] << " " << mesh.tris[i][1] << " " << mesh.tris[i][2]
                 << endl;
            pausee();
        }

        for (int j = 0; j < 3; j++) {
            std::vector<int> tmp = set_intersection(mesh.tri_vertices[mesh.tris[i][(j + 1) % 3]].conn_tris,
                                                    mesh.tri_vertices[mesh.tris[i][(j + 2) % 3]].conn_tris);
            if (mesh.is_bbox_es[i][j]) {
                if (tmp.size() != 1) {
                    cout << "mesh.is_bbox_es[i][j] && tmp.size()!=1" << endl;
                    pausee();
                }
            } else {
                if (tmp.size() != 2) {
                    cout << "!mesh.is_bbox_es[i][j] && tmp.size()!=2" << endl;
                    pausee();
                }
            }

            if (mesh.is_boundary_es[i][j]) {
                int n_t_id = i == tmp[0] ? tmp[1] : tmp[0];
                for (int k = 0; k < 3; k++) {
                    if (mesh.tris[n_t_id][k] != mesh.tris[i][(j + 1) % 3]
                        && mesh.tris[n_t_id][k] != mesh.tris[i][(j + 2) % 3]) {
                        if (!mesh.is_boundary_es[n_t_id][k]) {
                            for(auto &nfid: tmp)
                                cout<<"f: "<<nfid<<endl;
                            cout<<"v0 "<<mesh.tris[i][(j + 1) % 3]<<" v1 "<<mesh.tris[i][(j + 2) % 3]<<endl;
                            cout << "is_boundary_es are not paired" << endl;
                            cout<<"t_id = "<<i<<endl;
                            cout<<mesh.tris[i][0]<<" "<<mesh.tris[i][1]<<" "<<mesh.tris[i][2]<<endl;
                            cout<<mesh.is_boundary_es[i][0]<<" "<<mesh.is_boundary_es[i][1]<<" "<<mesh.is_boundary_es[i][2]<<endl;
                            cout<<"n_t_id = "<<n_t_id<<endl;
                            cout<<mesh.tris[n_t_id][0]<<" "<<mesh.tris[n_t_id][1]<<" "<<mesh.tris[n_t_id][2]<<endl;
                            cout<<mesh.is_boundary_es[n_t_id][0]<<" "<<mesh.is_boundary_es[n_t_id][1]<<" "<<mesh.is_boundary_es[n_t_id][2]<<endl;
                            pausee();
                        }
                        break;
                    }
                }
            }
        }
    }

    //extra feature checks
    for (int i = 0; i < mesh.tris.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;

        for (int j = 0; j < 3; j++) {
            std::vector<int> tmp = set_intersection(mesh.tri_vertices[mesh.tris[i][(j + 1) % 3]].conn_tris,
                                                    mesh.tri_vertices[mesh.tris[i][(j + 2) % 3]].conn_tris);

            if (mesh.tag_feature_es[i][j] >=0 ) {
                int n_t_id = i == tmp[0] ? tmp[1] : tmp[0];
                for (int k = 0; k < 3; k++) {
                    if (mesh.tris[n_t_id][k] != mesh.tris[i][(j + 1) % 3]
                        && mesh.tris[n_t_id][k] != mesh.tris[i][(j + 2) % 3]) {
                        if (mesh.tag_feature_es[i][j] != mesh.tag_feature_es[n_t_id][k]) {
                            cout << "tag_feature_es are not paired" << endl;
                            pausee();
                        }
                        break;
                    }
                }

                if (!mesh.tri_vertices[mesh.tris[i][(j + 1) % 3]].has_feature(mesh.tag_feature_es[i][j])) {
                    cout << "v " << mesh.tris[i][(j + 1) % 3] << "get_t() failed" << endl;
                    pausee();
                }
                if (!mesh.tri_vertices[mesh.tris[i][(j + 2) % 3]].has_feature(mesh.tag_feature_es[i][j])) {
                    cout << "v " << mesh.tris[i][(j + 2) % 3] << "get_t() failed" << endl;
                    pausee();
                }
            }


            if (mesh.tag_secondary_feature_es[i][j] >=0 ) {
                int n_t_id = i == tmp[0] ? tmp[1] : tmp[0];
                for (int k = 0; k < 3; k++) {
                    if (mesh.tris[n_t_id][k] != mesh.tris[i][(j + 1) % 3]
                        && mesh.tris[n_t_id][k] != mesh.tris[i][(j + 2) % 3]) {
                        if (mesh.tag_secondary_feature_es[i][j] != mesh.tag_secondary_feature_es[n_t_id][k]) {
                            cout << "tag_secondary_feature_es are not paired" << endl;
                            pausee();
                        }
                        break;
                    }
                }
            }
        }
    }

    for(int i=0;i<mesh.tri_vertices.size();i++) {
        if (mesh.v_is_removed[i])
            continue;

        std::vector<int> tmp;
        for (auto &info:mesh.tri_vertices[i].feature_infos) {
            tmp.push_back((int)info[0]);
//            Point_2f p = feature::features[info[0]]->eval(info[1]);
//            double l = (p - mesh.tri_vertices[i].posf).length_2();
//            if (l > 1e-3) {
//                cout << "(p-mesh.tri_vertices[i].posf).length_2()>1e-3" << endl;
//                cout << "p: " << p << endl;
//                cout << "posf: " << mesh.tri_vertices[i].posf << endl;
//                cout << "v " << i << ", feature " << info[0] << " " << feature::features[info[0]]->type << endl;
//                pausee();
//            }
        }
        vector_unique(tmp);
        if(tmp.size()!=mesh.tri_vertices[i].feature_infos.size()){
            cout<<"duplicate feature id!"<<endl;
            pausee();
        }

        if(mesh.tri_vertices[i].feature_infos.size()>1 && !mesh.tri_vertices[i].is_freezed){
            cout<<"feature_infos.size()>1 && !is_freezed"<<endl;
            pausee();
        }
    }
}

void triwild::optimization::output_info(MeshData& mesh) {
    if(args.mute_log)
        return;

    //output infos
    //#
    int v_n = std::count(mesh.v_is_removed.begin(), mesh.v_is_removed.end(), false);
    int t_n = std::count(mesh.t_is_removed.begin(), mesh.t_is_removed.end(), false);
    int cnt = 0;
    for (size_t i = 0; i < mesh.tri_vertices.size(); i++) {
        if (mesh.v_is_removed[i])
            continue;
        if (!mesh.tri_vertices[i].is_rounded)
            cnt++;
    }
    cout << "#v = " << v_n << "(" << mesh.tri_vertices.size() << ")" << "    unrounded = " << cnt << endl;
    cout << "#t = " << t_n << "(" << mesh.tris.size() << ")" << endl;
    //energy
    double max_energy, avg_energy;
    get_max_avg_energy(mesh, max_energy, avg_energy);
    cout << "max_energy = " << max_energy<<endl;
    cout << "avg_energy = " << avg_energy<<endl;
    //angle
    double g_min_angle = 180;
    double g_max_angle = 0;
    std::array<double, 6> cnt5 = {0, 0, 0, 0, 0, 0};
    for (int i = 0; i < mesh.tris.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;
        double min_angle = 180;
        double max_angle = 0;
        for (int j = 0; j < 3; j++) {
            Point_2f t1 = mesh.tri_vertices[mesh.tris[i][1]].posf - mesh.tri_vertices[mesh.tris[i][0]].posf;
            Point_2f t2 = mesh.tri_vertices[mesh.tris[i][2]].posf - mesh.tri_vertices[mesh.tris[i][0]].posf;
            double l1 = std::sqrt(t1.length_2());
            double l2 = std::sqrt(t2.length_2());
            if (l1 == 0 || l2 == 0) {
                min_angle = 0;
                max_angle = 180;
                break;
            }
            double cos_a = t1.dot(t2) / (l1 * l2);
            if (std::isnan(cos_a) || std::isinf(cos_a)) {
                cout << "std::isnan(cos_a)||std::isinf(cos_a)" << endl;
                pausee();
            }
            if (cos_a > 1)
                cos_a = 1;
            else if (cos_a < -1)
                cos_a = -1;
            double a = std::acos(cos_a) / 3.1415926 * 180;
            if (a < min_angle)
                min_angle = a;
            if (a > max_angle)
                max_angle = a;
        }
        if (min_angle < g_min_angle)
            g_min_angle = min_angle;
        if (max_angle > g_max_angle)
            g_max_angle = max_angle;

        for (int j = 0; j < 3; j++) {
            if (min_angle < 5 * (j + 1))
                cnt5[j]++;
            if (max_angle > 180 - 5 * (j + 1))
                cnt5[5 - j]++;
        }
    }
    cout << "min_angle = " << g_min_angle << ", max_angle = " << g_max_angle << endl;
    cout << "< 5: " << cnt5[0] / t_n << ", " << "< 10: " << cnt5[1] / t_n << ", " << "< 15: " << cnt5[2] / t_n << endl;
    cout << "> 175: " << cnt5[5] / t_n << ", " << "> 170: " << cnt5[4] / t_n << ", " << "> 165: " << cnt5[3] / t_n
         << endl;
    cout << endl;

    check(mesh);
}

void triwild::optimization::output_boundary(MeshData& mesh, bool is_feature) {
    std::vector<std::array<int, 2>> b_es;
    for (int i = 0; i < mesh.tris.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;
        for (int j = 0; j < 3; j++) {
            if ((!is_feature && mesh.is_boundary_es[i][j]) || (is_feature && mesh.tag_feature_es[i][j] >= 0)) {
                std::array<int, 2> e = {mesh.tris[i][(j + 1) % 3], mesh.tris[i][(j + 2) % 3]};
                if (e[0] < e[1])
                    b_es.push_back(e);
                else {
                    std::swap(e[0], e[1]);
                    b_es.push_back(e);
                }
            }
        }
    }
    vector_unique(b_es);

    Eigen::MatrixXd V(b_es.size() * 3, 3);
    Eigen::MatrixXi F(b_es.size(), 3);
    cout << b_es.size() << endl;
    for (int i = 0; i < b_es.size(); i++) {
        V(i * 3 + 0, 0) = mesh.tri_vertices[b_es[i][0]].posf[0];
        V(i * 3 + 0, 1) = mesh.tri_vertices[b_es[i][0]].posf[1];
        V(i * 3 + 0, 2) = 0;

        V(i * 3 + 1, 0) = mesh.tri_vertices[b_es[i][1]].posf[0];
        V(i * 3 + 1, 1) = mesh.tri_vertices[b_es[i][1]].posf[1];
        V(i * 3 + 1, 2) = 0;

        V(i * 3 + 2, 0) = mesh.tri_vertices[b_es[i][1]].posf[0];
        V(i * 3 + 2, 1) = mesh.tri_vertices[b_es[i][1]].posf[1];
        V(i * 3 + 2, 2) = 0;

        F(i, 0) = i * 3 + 0;
        F(i, 1) = i * 3 + 1;
        F(i, 2) = i * 3 + 2;
    }

    igl::writeSTL(args.output + "_boundary.stl", V, F);
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    for (auto& el : vec) {
        os << el << ' ';
    }
    return os;
}

void triwild::optimization::pausee() {
    cout << "Is pausing... (Enter '0' to exit and other characters to continue.)" << endl;
    char c;
    cin >> c;
    if (c == '0')
        exit(0);
}

std::vector<int> triwild::optimization::set_intersection(const std::vector<int>& iv1, const std::vector<int>& iv2){
    std::vector<int> tmp;
    if(iv1.empty() || iv2.empty())
        return tmp;
    std::vector<int> v1 = iv1;
    std::vector<int> v2 = iv2;
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(tmp));
    return tmp;
}

std::vector<int> triwild::optimization::set_intersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2) {
    std::vector<int> tmp;
    if (s1.empty() || s2.empty())
        return tmp;

    if (s2.size() < s1.size()) { return set_intersection(s2, s1); }
    tmp.clear();
    tmp.reserve(std::min(s1.size(), s2.size()));
    for (int x : s1) {
        if (s2.count(x)) {
            tmp.push_back(x);
        }
    }
    std::sort(tmp.begin(), tmp.end());

    return tmp;

//    std::vector<int> v1;
//    std::vector<int> v2;
//    v1.reserve(s1.size());
//    v2.reserve(s2.size());
//    for(const auto& i: s1)
//        v1.push_back(i);
//    for(const auto& i: s2)
//        v2.push_back(i);
//    std::sort(v1.begin(), v1.end());
//    std::sort(v2.begin(), v2.end());
//    std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(tmp));
//    return tmp;
}

