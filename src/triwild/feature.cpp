// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "feature.h"
#include "optimization.h"
#include "Curves.h"

#include "edge_splitting.h"
#include "edge_collapsing.h"
#include "edge_swapping.h"
#include "vertex_smoothing.h"

#include "Logger.h"

#include "../../extern/aabbcc/src/AABB.h"
#include <igl/unique_rows.h>

namespace triwild {
    namespace feature {
        std::vector<std::shared_ptr<FeatureElement>> features;
        std::vector<std::shared_ptr<FeatureElement>> secondary_features;

        double feature_eps;
        double feature_eps_2;
    }
}

using namespace triwild;

bool triwild::feature::init(const std::string& feature_file){
    json feature_info = json({});
    if (feature_file != "") {
        std::ifstream file(feature_file);
        if (file.is_open())
            file >> feature_info;
        else {
            std::cerr << "unable to open " << feature_file << " file" << std::endl;
            return false;
        }
        file.close();
    } else
        return false;

    std::vector<int> tmp = {4179, 4189};//2321, 2322,1295, 1346,  4179, 4189
    for (int i = 0; i < feature_info.size(); i++) {
//        if(std::find(tmp.begin(), tmp.end(), i) != tmp.end())
//            continue;
        std::vector<int> v_ids;
        v_ids.reserve(feature_info[i]["v_ids"].size());
        for (int j = 0; j < feature_info[i]["v_ids"].size(); j++)
            v_ids.push_back((int) feature_info[i]["v_ids"][j]);
        std::vector<double> paras;
        paras.reserve(feature_info[i]["paras"].size());
        for (int j = 0; j < feature_info[i]["paras"].size(); j++)
            paras.push_back((double) feature_info[i]["paras"][j]);

        if (feature_info[i]["type"] == "Line") {
            features.push_back(std::make_shared<Line_Feature>(v_ids, paras,
                                                              json2point(feature_info[i]["start"]),
                                                              json2point(feature_info[i]["end"]),
                                                              feature_info[i]["curve_id"]));
        } else if (feature_info[i]["type"] == "RationalBezier") {
            features.push_back(std::make_shared<RationalBezierCurve_Feature>(v_ids, paras,
                                                                             json2d2ctrlvector(feature_info[i]["poles"]),
                                                                             json1d2ctrlvector(feature_info[i]["weigths"]),
                                                                             feature_info[i]["curve_id"]));
        } else if (feature_info[i]["type"] == "BezierCurve") {
            features.push_back(std::make_shared<BezierCurve_Feature>(v_ids, paras,
                                                                     feature_info[i]["degree"],
                                                                     json2d2ctrlvector(feature_info[i]["poles"]),
                                                                     feature_info[i]["curve_id"]));
        } else {
            cout << "Unknown curve typle " << feature_info[i]["type"] << " will not be preserved" << endl;
        }
    }

    return true;
}

void triwild::feature::map_feature2mesh(MeshData& mesh) {
    std::vector<std::vector<int>> v_feature_ids(mesh.tri_vertices.size(), std::vector<int>());

    for (int feature_id = 0; feature_id < features.size(); feature_id++) {
        auto feature = features[feature_id];
        std::vector<int> tmp_feature_v_ids;
        tmp_feature_v_ids.reserve(feature->v_ids.size() * 2);

        for (int i = 0; i < feature->v_ids.size() - 1; i++) {
            int next_v_id = feature->v_ids[i + 1];

            std::vector<bool> is_visited(mesh.tri_vertices.size(), false);
            std::queue<int> v_queue;
            v_queue.push(feature->v_ids[i]);
            while (!v_queue.empty()) {
                int v_id = v_queue.front();
                tmp_feature_v_ids.push_back(v_id);
                v_queue.pop();

                std::unordered_set<int> n_v_ids;
                for (int t_id:mesh.tri_vertices[v_id].conn_tris) {
                    for (int tmp_v_id:mesh.tris[t_id])
                        if (tmp_v_id != v_id && !is_visited[tmp_v_id])
                            n_v_ids.insert(tmp_v_id);
                }
//                for (int t_id:mesh.tri_vertices[v_id].conn_tris) {
//                    int j = std::find(mesh.tris[t_id].begin(), mesh.tris[t_id].end(), v_id) - mesh.tris[t_id].begin();
//                    if (!is_visited[mesh.tris[t_id][(j + 1) % 3]] && mesh.is_boundary_es[t_id][(j + 2) % 3]) {
//                        n_v_ids.insert(mesh.tris[t_id][(j + 1) % 3]);
//                    }
//                    if (!is_visited[mesh.tris[t_id][(j + 2) % 3]] && mesh.is_boundary_es[t_id][(j + 1) % 3]) {
//                        n_v_ids.insert(mesh.tris[t_id][(j + 2) % 3]);
//                    }
//                }

                //mark is_visited
                is_visited[v_id] = true;
                for (int n_v_id:n_v_ids)
                    is_visited[n_v_id] = true;

                //find next_v_id
                bool is_found = false;
                for (int n_v_id:n_v_ids) {
                    if (n_v_id == next_v_id) {
                        is_found = true;
                        break;
                    }
                }
                if (is_found && next_v_id == feature->v_ids.back())
                    break;

                //find v in the middle
                if (!is_found) {
                    for (int n_v_id:n_v_ids) {
                        if (is_on_segment(mesh.tri_vertices[n_v_id].pos, mesh.tri_vertices[v_id].pos,
                                          mesh.tri_vertices[next_v_id].pos)) {
                            v_queue.push(n_v_id);
                            break;
                        }
                    }
                }
            }
        }
        tmp_feature_v_ids.push_back(feature->v_ids.back());

        //update paras
        std::vector<double> old_feature_paras = feature->paras;
        int cnt = 0;
        for (int i = 0; i < tmp_feature_v_ids.size(); i++) {
            int v_id = tmp_feature_v_ids[i];
            if (v_id != feature->v_ids[cnt]) {
                double l = std::sqrt((mesh.tri_vertices[feature->v_ids[cnt - 1]].posf
                                      - mesh.tri_vertices[feature->v_ids[cnt]].posf).length_2());
                double lv = std::sqrt((mesh.tri_vertices[v_id].posf
                                       - mesh.tri_vertices[feature->v_ids[cnt]].posf).length_2());
                double start_t = old_feature_paras[cnt];
                if (l > 0)
                    start_t = (1 - lv / l) * old_feature_paras[cnt - 1] + (lv / l) * old_feature_paras[cnt];
                double tmp_t = feature->inv_eval(mesh.tri_vertices[v_id].posf,
                                                 start_t, old_feature_paras[cnt - 1], old_feature_paras[cnt]);
                if (tmp_t < old_feature_paras[cnt - 1] || tmp_t > old_feature_paras[cnt])
                    tmp_t = start_t;

                feature->paras.insert(feature->paras.begin() + i, tmp_t);

//                if (feature->type == "Line") {
//                    double tmp_t = feature->inv_eval(mesh.tri_vertices[v_id].posf,
//                                                     start_t, old_feature_paras[cnt - 1], old_feature_paras[cnt]);
//                    feature->paras.insert(feature->paras.begin() + i, tmp_t);
//                } else {
//                    double tmp_t = feature->inv_eval(mesh.tri_vertices[v_id].posf,
//                                                     start_t, old_feature_paras[cnt - 1], old_feature_paras[cnt]);
//                    if (tmp_t < old_feature_paras[cnt - 1] || tmp_t > old_feature_paras[cnt]) {
//                        std::cout << "starting point " << start_t << std::endl;
//                        std::cout << "[cnt-1] " << mesh.tri_vertices[feature->v_ids[cnt - 1]].posf << std::endl;
//                        std::cout << "[cnt] " << mesh.tri_vertices[feature->v_ids[cnt]].posf << std::endl;
//                        std::cout << old_feature_paras[cnt - 1] << std::endl;
//                        std::cout << old_feature_paras[cnt] << std::endl;
//                        cout << "feature " << feature->curve_id << endl;
//                        cout << "feature " << feature->type << endl;
//                        cout << "\n--------\n" << feature->to_maple() << "\n---------" << endl;
//                        cout << "t " << tmp_t << endl;
//                        cout << "p " << mesh.tri_vertices[v_id].posf << endl;
//                        optimization::pausee();
//                        tmp_t = start_t;
//                    }
//                    feature->paras.insert(feature->paras.begin() + i, tmp_t);
//                }
            } else
                cnt++;

            mesh.tri_vertices[v_id].feature_infos.push_back(
                    std::array<double, 2>({(double) feature_id, feature->paras[i]}));
        }
        //update v_ids
        feature->v_ids = tmp_feature_v_ids;
    }

    for (int feature_id = 0; feature_id < features.size(); feature_id++) {
        for (int i = 0; i < features[feature_id]->v_ids.size() - 1; i++) {
            int v1_id = features[feature_id]->v_ids[i];
            int v2_id = features[feature_id]->v_ids[i + 1];
            for (int t_id:mesh.tri_vertices[v1_id].conn_tris) {
                int j2 = std::find(mesh.tris[t_id].begin(), mesh.tris[t_id].end(), v2_id) - mesh.tris[t_id].begin();
                if (j2 < 3) {
                    int j = mesh.tris[t_id][(j2 + 1) % 3] == v1_id ? (j2 + 2) % 3 : (j2 + 1) % 3;
                    mesh.tag_feature_es[t_id][j] = feature_id;
                }
            }
        }
    }

    for (auto feature:features) {
        mesh.tri_vertices[feature->v_ids.front()].is_freezed = true;
        mesh.tri_vertices[feature->v_ids.back()].is_freezed = true;

        for (int i = 1; i < feature->v_ids.size() - 1; i++) {//deal with (partially) overlapped features
            std::unordered_set<int> tmp_feature_ids;
            for (int t_id:mesh.tri_vertices[feature->v_ids[i]].conn_tris) {
                for (int j = 0; j < 3; j++)
                    if (mesh.tag_feature_es[t_id][j] >= 0)
                        tmp_feature_ids.insert(mesh.tag_feature_es[t_id][j]);
            }
            if (tmp_feature_ids.size() == mesh.tri_vertices[feature->v_ids[i]].feature_infos.size())
                continue;

            auto &infos = mesh.tri_vertices[feature->v_ids[i]].feature_infos;
            for (int j = 0; j < infos.size(); j++) {
                if (tmp_feature_ids.find(infos[j][0]) == tmp_feature_ids.end()) {//feature_id info[0] is invalid
                    cout << "remove feature " << infos[j][0] << " from vertex " << feature->v_ids[i] << endl;
                    infos.erase(infos.begin() + j);
                    j--;
                }
            }
        }
    }

    for (int i = 0; i < mesh.tri_vertices.size(); i++) {
        if (mesh.tri_vertices[i].feature_infos.size() > 1)
            mesh.tri_vertices[i].is_freezed = true;
    }

//    for (int i = 0; i < mesh.tri_vertices.size(); i++) {
//        if (mesh.v_is_removed[i])
//            continue;
//        if (mesh.tri_vertices[i].feature_infos.size() != 1)//==0?
//            continue;
//
//        auto &info = mesh.tri_vertices[i].feature_infos[0];
//        Point_2f pt = feature::features[info[0]]->eval(info[1]);
//        if (pt != mesh.tri_vertices[i].posf) {
//            double l = (pt - mesh.tri_vertices[i].posf).length();
//            if (l > 1e-1) {
//                cout << "bbox" << endl;
//                cout << args.box_min(0) << " " << args.box_min(1) << endl;
//                cout << args.box_max(0) << " " << args.box_max(1) << endl;
//                cout << "l = " << l << endl;
//                cout << "v_id = " << i << endl;
//                cout << feature::features[info[0]]->curve_id << endl;
//                cout << feature::features[info[0]]->type << endl;
//                cout << mesh.tri_vertices[i].is_freezed << endl;
//                optimization::pausee();
//            }
//        }
//    }
}

bool triwild::feature::is_on_segment(const Point_2& p, const Point_2& p1, const Point_2& p2) {
    if ((p.x > p1.x && p.x > p2.x) || (p.x < p1.x && p.x < p2.x)
        || (p.y > p1.y && p.y > p2.y) || (p.y < p1.y && p.y < p2.y))
        return false;

    if ((p.y - p1.y) * (p2.x - p1.x) == (p.x - p1.x) * (p2.y - p1.y))
        return true;

    return false;
}

void triwild::feature::output_stats(MeshData& mesh, std::ofstream& f) {
    ////unsnapped
    int cnt = 0;
    for (int i = 0; i < mesh.tri_vertices.size(); i++) {
        if (mesh.v_is_removed[i])
            continue;
        if (mesh.tri_vertices[i].feature_infos.size() != 1)//==0?
            continue;

        auto &info = mesh.tri_vertices[i].feature_infos[0];
        Point_2f pt = feature::features[info[0]]->eval(info[1]);
        if (pt != mesh.tri_vertices[i].posf) {
            cnt++;
        }
    }
    Logger::instance().unsnapped_vertices = cnt;

    ////invalid
    cnt = 0;
    for (int i = 0; i < mesh.tris.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;
        if (!mesh.tri_nodes[i].empty()) {
            std::vector<Point_2f> nodes;
            for (int n_id:mesh.tri_nodes[i])
                nodes.push_back(mesh.nodes[n_id]);
            if (!is_valid_inversion({{mesh.tri_vertices[mesh.tris[i][0]].posf, mesh.tri_vertices[mesh.tris[i][1]].posf,
                                             mesh.tri_vertices[mesh.tris[i][2]].posf}}, nodes))
                cnt++;
        }
    }
    Logger::instance().invalid_element = cnt;

    f << "curved_cnt_ls_fitted " << Logger::instance().ls_fitted << "\n";
    f << "curved_cnt_ls_fitted_fail " << Logger::instance().ls_fitted_fail << "\n";
    f << "curved_cnt_unsnapped_v " << Logger::instance().unsnapped_vertices << "\n";
    f << "curved_cnt_invalid_t " << Logger::instance().invalid_element << "\n";
    f << "curved_ls_fitting_distances ";
    for (const auto d : Logger::instance().ls_fitting_distances)
        f << d << " ";
    f << "\n";
}

void triwild::feature::snap_vertices(MeshData& mesh) {
    int cnt = 0;
    int cnt_s = 0;

    for (int i = 0; i < mesh.tri_vertices.size(); i++) {
        if (mesh.v_is_removed[i])
            continue;
        if (mesh.tri_vertices[i].is_freezed)
            continue;
        if (mesh.tri_vertices[i].feature_infos.size() != 1)//==0?
            continue;

        auto &info = mesh.tri_vertices[i].feature_infos[0];
        Point_2f pt = feature::features[info[0]]->eval(info[1]);
        if (pt != mesh.tri_vertices[i].posf) {
            cnt++;

            double l = (pt - mesh.tri_vertices[i].posf).length();
            Vector_2f d = pt - mesh.tri_vertices[i].posf;
            for (int n = 0; n < 10; n++) {
                double a = l / std::pow(2, n);
                Point_2f pf;
                if(n == 0)
                    pf = pt;
                else
                    pf = mesh.tri_vertices[i].posf + d * a;
                Point_2 p(pf[0], pf[1]);
                bool is_valid = true;
                for (int t_id:mesh.tri_vertices[i].conn_tris) {
                    int j = std::find(mesh.tris[t_id].begin(), mesh.tris[t_id].end(), i) - mesh.tris[t_id].begin();
                    if (!optimization::is_valid_inversion(pf, mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].posf,
                                                          mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].posf,
                                                          p, mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].pos,
                                                          mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].pos)) {

                        is_valid = false;
                        break;
                    }
                }
                if (is_valid) {
                    mesh.tri_vertices[i].posf = pf;
                    mesh.tri_vertices[i].pos = p;

                    //update energy
                    for (int t_id:mesh.tri_vertices[i].conn_tris)
                        mesh.t_quality[t_id] = optimization::tri_energy(mesh.tri_vertices[mesh.tris[t_id][0]].posf,
                                                                        mesh.tri_vertices[mesh.tris[t_id][1]].posf,
                                                                        mesh.tri_vertices[mesh.tris[t_id][2]].posf);

                    if (mesh.is_curved) {
                        for (int t_id:mesh.tri_vertices[i].conn_tris) {
                            if (mesh.tri_nodes[t_id].empty())
                                continue;
                            std::vector<Point_2f> new_nodes;
                            get_new_nodes(mesh.tag_feature_es[t_id], {{mesh.tri_vertices[mesh.tris[t_id][0]],
                                                                              mesh.tri_vertices[mesh.tris[t_id][1]],
                                                                              mesh.tri_vertices[mesh.tris[t_id][2]]}},
                                          new_nodes);
                            for (int j = 0; j < mesh.tri_nodes[t_id].size(); j++) {
                                mesh.nodes[mesh.tri_nodes[t_id][j]] = new_nodes[j];
                            }
                        }
                    }

                    cnt_s++;
                    break;
                }
            }
        }
    }

    cout << "success " << cnt_s << "/" << cnt << endl;
}

void triwild::feature::merge_inflection(MeshData& mesh) {
    std::vector<int> indices(features.size()), map(features.size());
    for (int i = 0; i < features.size(); ++i) {
        indices[i] = i;
        map[i] = i;
    }


    std::sort(indices.begin(), indices.end(), [&](const int &ia, const int &ib) -> bool {
        const std::shared_ptr<FeatureElement> a = features[ia];
        const std::shared_ptr<FeatureElement> b = features[ib];
        if (a->curve_id == b->curve_id) {
            return a->paras.front() < b->paras.front();
        }

        return a->curve_id > b->curve_id;
    });

    int current_f_id = -1;
    double current_end = -1;
    std::shared_ptr<FeatureElement> current_feature;
    int current_feature_index = -1;

    for (int index : indices) {
        const auto feature_ptr = features[index];
        if (feature_ptr->curve_id != current_f_id) {
            current_f_id = feature_ptr->curve_id;
            current_end = feature_ptr->paras.back();
            current_feature = feature_ptr;
            current_feature_index = index;

            //end point of this feature is inflection
            // if(!feature_ptr->is_inflection[1])
            // current_f_id=-1;
            continue;
        }

        // if(feature_ptr->is_inflection[1])
        // current_f_id=-1;

//        std::cout << "same id" << current_end << " vs " << feature_ptr->paras.front() << std::endl;

        //previous id is equal current id

        //previous end point matches current start
        if (feature_ptr->paras.front() == current_end) {
            current_feature->merge_after(*feature_ptr);
            current_end = feature_ptr->paras.back();

//            std::cout << "merging" << std::endl;

            map[index] = current_feature_index;
        } else {
            //there was a gap, current feature gets updated
            current_feature = feature_ptr;
            current_feature_index = index;
            current_end = feature_ptr->paras.back();
        }
    }

    for (int i = 0; i < mesh.tris.size(); ++i) {
        if (mesh.t_is_removed[i])
            continue;

        for (int t = 0; t < 3; ++t) {
            if (mesh.tag_feature_es[i][t] >= 0)
                mesh.tag_feature_es[i][t] = map[mesh.tag_feature_es[i][t]];
        }
    }

    for (int i = 0; i < mesh.tri_vertices.size(); ++i) {
        if (mesh.v_is_removed[i])
            continue;

        for (auto &infos: mesh.tri_vertices[i].feature_infos)
            infos[0] = map[infos[0]];

        auto &infos = mesh.tri_vertices[i].feature_infos;
        std::sort(infos.begin(), infos.end(), [](const std::array<double, 2> &a, const std::array<double, 2> &b) {
            return a[0] < b[0];
        });
        int old_size = infos.size();
        for (int j = 1; j < infos.size(); j++) {
            if (infos[j][0] == infos[j - 1][0]) {
                infos.erase(infos.begin() + j);
                j--;
            }
        }
        if (old_size > 1 && infos.size() == 1) {
            mesh.tri_vertices[i].is_freezed = false;
//            cout<<"unfreezed "<<mesh.tri_vertices[i].feature_infos[0][1]<<endl;
        }
    }
}

#include <igl/triangle/triangulate.h>
#include "auto_p_bases.hpp"
#include "meshio.hpp"
void triwild::feature::curving(MeshData& mesh, GEO::MeshFacetsAABB &b_tree) {
    ////subdivide elements with >1 curved edges
    subdivide_into_2(mesh);
    subdivide_into_3(mesh);

    optimization::check(mesh);//tmp

    ////curving
    add_nodes(mesh);
    mesh.is_curved = true;

//    check_inversion(mesh);
//    write_msh(mesh, args.output + "before" + ".msh");

    merge_inflection(mesh);

    //simplification
    //done in the end of optimization: reset scales, and free the limit on edge length
    auto check = [&]() -> void {
        if(args.mute_log)
            return;
        for (int i = 0; i < mesh.tag_feature_es.size(); i++) {
            if (mesh.t_is_removed[i])
                continue;
            bool is_curved = false;
            for (int j = 0; j < 3; j++) {
                if (mesh.tag_feature_es[i][j] >= 0 && features[mesh.tag_feature_es[i][j]]->type != "Line") {
                    is_curved = true;
                    assert(!mesh.tri_nodes[i].empty());
//                    break;
                }
            }
            if(!is_curved)
                assert(mesh.tri_nodes[i].empty());
        }
    };

    mesh.is_limit_length = true;
    optimization::operation(mesh, b_tree, {{0, 1, 1, 1}});
    for (int i = 0; i < 10; i++) {
        cout << "//////////////// CURVED pass " << i << " ////////////////" << endl;

        double max_energy, avg_energy;
        optimization::get_max_avg_energy(mesh, max_energy, avg_energy);

        optimization::operation(mesh, b_tree, {{1, 1, 1, 1}});
        check();

        double new_max_energy, new_avg_energy;
        optimization::get_max_avg_energy(mesh, new_max_energy, new_avg_energy);

        if (max_energy - new_max_energy < 1e-2 && avg_energy - new_avg_energy < 1e-3)
            break;
    }

    Logger::instance().ls_fitted = 0;
    Logger::instance().ls_fitted_fail = 0;
    Logger::instance().ls_fitting_distances.clear();

    fix_inversion(mesh);
    fix_inversion(mesh);
    if(args.mute_log)
        return;
    check_inversion(mesh, true);
}

void triwild::feature::check_inversion(MeshData& mesh, bool is_output_objs){
    if(args.mute_log)
        return;

    ////check/resolve inversions
    std::vector<int> invalid_t_ids;
    for (int i = 0; i < mesh.tris.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;
        if (!mesh.tri_nodes[i].empty()) {
            assert(!mesh.tri_nodes[i].empty());
            std::vector<Point_2f> nodes;
            for (int n_id:mesh.tri_nodes[i])
                nodes.push_back(mesh.nodes[n_id]);
            if (!is_valid_inversion({{mesh.tri_vertices[mesh.tris[i][0]].posf, mesh.tri_vertices[mesh.tris[i][1]].posf,
                                             mesh.tri_vertices[mesh.tris[i][2]].posf}}, nodes)) {
//                cout << "----------\n"<<i<<" inverted curved element\n";
//
//                std::cout<<mesh.tri_vertices[mesh.tris[i][0]].posf<<"\n";
//                std::cout<<mesh.tri_vertices[mesh.tris[i][1]].posf<<"\n";
//                std::cout<<mesh.tri_vertices[mesh.tris[i][2]].posf<<"\n";
//                for(auto &n : nodes)
//                    std::cout<<n<<"\n";
//                cout<<"----------\n"<<endl;
//
//                int id = -1;
//                int index = -1;
//                for(int ii = 0; ii < 3; ++ii){
//                    auto asd = mesh.tag_feature_es[i][ii];
//                    if(asd > id)
//                    {
//                        id = asd;
//                        index = ii;
//                    }
//                }
//                assert(mesh.tag_feature_es[i][index] == id);
//                const TriVertex &v1 = mesh.tri_vertices[mesh.tris[i][(index + 1) % 3]];
//                const TriVertex &v2 = mesh.tri_vertices[mesh.tris[i][(index + 2) % 3]];
//
//                cout<<"----------\n"<<endl;
//                for(int lv = 0; lv < 3; ++lv)
//                {
//                    const auto &vvv = mesh.tri_vertices[mesh.tris[i][lv]];
//
//                    for(int fid : vvv.conn_tris)
//                    {
//                        std::cout<<mesh.tris[fid][0]<<" "<<mesh.tris[fid][1]<< " "<<mesh.tris[fid][2]<<std::endl;
//                    }
//                }
//
//                cout<<"----------\n"<<endl;
//
//                cout<<"----------\n"<<endl;
//                features[id]->print_info();
//                cout<<"----------\n"<<endl;
//                cout<<v1.get_t(id) << " "<< v2.get_t(id)<<endl;
//                cout<<"----------\n"<<endl;

                invalid_t_ids.push_back(i);
            }
        }
    }

    if (invalid_t_ids.size() > 0) {
        cout<<invalid_t_ids.size()<<" invalid found"<<endl;
        mesh.t_scalars = std::vector<double>(mesh.tris.size(), 0);

        std::ofstream f;
        f.open(args.output + ".invalid");
        for (int t_id:invalid_t_ids) {
            f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
            f << mesh.tri_vertices[mesh.tris[t_id][0]].posf << endl;
            f << mesh.tri_vertices[mesh.tris[t_id][1]].posf << endl;
            f << mesh.tri_vertices[mesh.tris[t_id][2]].posf << endl;
            for (int n_id:mesh.tri_nodes[t_id])
                f << mesh.nodes[n_id] << endl;
            f << endl;
            mesh.t_scalars[t_id] = 1;

            if (is_output_objs) {
                Eigen::MatrixXd V(3, 2), V_out, _, CN(1, 3), V_waped, tmp;
                Eigen::MatrixXi E(3, 2), F, FN;
                V << 0, 0,
                        1, 0,
                        0, 1;
                E << 0, 1,
                        1, 2,
                        2, 0;
                CN << 0, 0, 1;
                igl::triangle::triangulate(V, E, _, "Qq10a" + std::to_string(0.0001), V_out, F);
                FN = Eigen::MatrixXi::Ones(F.rows(), 3);

                V_waped.resizeLike(V_out);
                V_waped.setZero();

                int base_order = 3;
                if (mesh.tri_nodes[t_id].size() == 3)
                    base_order = 2;
                for (int j = 0; j < 3; ++j) {
                    autogen::p_basis_value_2d(base_order, j, V_out, tmp);
                    V_waped.col(0) += tmp * mesh.tri_vertices[mesh.tris[t_id][j]].posf[0];
                    V_waped.col(1) += tmp * mesh.tri_vertices[mesh.tris[t_id][j]].posf[1];
                }
                for (int j = 0; j < mesh.tri_nodes[t_id].size(); ++j) {
                    autogen::p_basis_value_2d(base_order, 3 + j, V_out, tmp);
                    V_waped.col(0) += tmp * mesh.nodes[mesh.tri_nodes[t_id][j]][0];
                    V_waped.col(1) += tmp * mesh.nodes[mesh.tri_nodes[t_id][j]][1];
                }

                V_waped.conservativeResize(V_waped.rows(), 3);
                V_waped.col(2).setZero();

                igl::writeOBJ(args.output + "_invalid" + std::to_string(t_id) + ".obj", V_waped, F, CN, FN, V_out, F);
            }
        }
        f.close();
    } else
        mesh.t_scalars.clear();
}

#include "CurvedTriUntangler.hpp"
void triwild::feature::fix_inversion(MeshData& mesh) {
    std::vector<bool> is_skip(mesh.tris.size(), false);
    for (int i = 0; i < mesh.tris.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;
        if(is_skip[i])
            continue;
        if (mesh.tri_nodes[i].empty())
            continue;

        std::vector<Point_2f> nodes;
        for (int n_id:mesh.tri_nodes[i])
            nodes.push_back(mesh.nodes[n_id]);
        std::array<Point_2f, 3> vs = {{mesh.tri_vertices[mesh.tris[i][0]].posf, mesh.tri_vertices[mesh.tris[i][1]].posf,
                                              mesh.tri_vertices[mesh.tris[i][2]].posf}};
        if (is_valid_inversion(vs, nodes))
            continue;

        ///move center nodes
        Point_2f new_p = mesh.nodes[mesh.tri_nodes[i].back()];
        if (CurvedTriUntangler::untangle_center(vs, nodes, mesh.nodes[mesh.tri_nodes[i].back()], new_p)) {
            nodes.back() = new_p;
            if (is_valid_inversion(vs, nodes)) {
                mesh.nodes[mesh.tri_nodes[i].back()] = new_p;
                continue;
            }
        }

        ///ls fitting
//        nodes.back() = mesh.nodes[mesh.tri_nodes[i].back()];
        nodes.back() = (mesh.tri_vertices[mesh.tris[i][0]].posf + mesh.tri_vertices[mesh.tris[i][1]].posf
                + mesh.tri_vertices[mesh.tris[i][2]].posf) / 3.0;
        mesh.nodes[mesh.tri_nodes[i].back()] = nodes.back();

        int e_id;
        for (e_id = 0; e_id < 3; ++e_id) {
            if (mesh.tag_feature_es[i][e_id] >= 0 && features[mesh.tag_feature_es[i][e_id]]->type != "Line") {
                e_id = (e_id + 1) % 3;
                break;
            }
        }
        bool is_fitted = true;
        std::vector<Point_2f> fit_nodes;

//        cout<<"t_id "<<i<<endl;
//        cout<<vs[0]<<endl;
//        cout<<vs[1]<<endl;
//        cout<<vs[2]<<endl;
//        for(auto&n:nodes)
//            cout<<n<<endl;
//        cout<<e_id<<endl;
//        cout<<"////////////////////"<<endl;

        const double ls_dist = CurvedTriUntangler::ls_fit(vs, nodes, e_id, fit_nodes, false);
        if (ls_dist < 0) {
            is_fitted = false;
            cout << "fail!!!!" << endl;
        }
//        else
        Logger::instance().ls_fitting_distances.push_back(ls_dist);
        //check neighbor
        int v1_id = mesh.tris[i][e_id];
        int v2_id = mesh.tris[i][(e_id + 1) % 3];
        std::vector<int> tmp = optimization::set_intersection(mesh.tri_vertices[v1_id].conn_tris,
                                                              mesh.tri_vertices[v2_id].conn_tris);
        assert(tmp.size() == 2);
        bool found = false;
        int n_t_id = tmp[0] == i ? tmp[1] : tmp[0];
        is_skip[n_t_id] = true;
        for (int j = 0; j < 3; j++) {
            if (mesh.tris[n_t_id][j] == v2_id && mesh.tris[n_t_id][(j + 1) % 3] == v1_id) {
                std::vector<Point_2f> n_nodes;
                for (int n_id: mesh.tri_nodes[n_t_id])
                    n_nodes.push_back(mesh.nodes[n_id]);
                if (is_fitted) {
                    n_nodes[2 * j] = fit_nodes[1];
                    n_nodes[2 * j + 1] = fit_nodes[0];
                }
                if (!is_fitted || !is_valid_inversion(
                        {{mesh.tri_vertices[mesh.tris[n_t_id][0]].posf, mesh.tri_vertices[mesh.tris[n_t_id][1]].posf,
                                 mesh.tri_vertices[mesh.tris[n_t_id][2]].posf}}, n_nodes)) {

                    Logger::instance().ls_fitted_fail++;
                    //set to linear
                    Point_2f p1 = mesh.tri_vertices[v1_id].posf * 2 / 3.0
                                  + mesh.tri_vertices[v2_id].posf * 1 / 3.0;
                    Point_2f p2 = mesh.tri_vertices[v1_id].posf * 1 / 3.0
                                  + mesh.tri_vertices[v2_id].posf * 2 / 3.0;
                    mesh.nodes[mesh.tri_nodes[i][e_id * 2]] = p1;
                    mesh.nodes[mesh.tri_nodes[i][e_id * 2 + 1]] = p2;
                    mesh.nodes[mesh.tri_nodes[n_t_id][j * 2]] = p2;
                    mesh.nodes[mesh.tri_nodes[n_t_id][j * 2 + 1]] = p1;

                } else {
                    mesh.nodes[mesh.tri_nodes[i][e_id * 2]] = fit_nodes[0];
                    mesh.nodes[mesh.tri_nodes[i][e_id * 2 + 1]] = fit_nodes[1];
                    mesh.nodes[mesh.tri_nodes[n_t_id][j * 2]] = fit_nodes[1];
                    mesh.nodes[mesh.tri_nodes[n_t_id][j * 2 + 1]] = fit_nodes[0];
                    Logger::instance().ls_fitted++;
                }
                found = true;
                break;
            }
        }
        assert(found);
    }

    return;
}

void triwild::feature::subdivide_into_2(MeshData& mesh) {
    std::vector<std::array<int, 2>> edges;
    for (int i = 0; i < mesh.tris.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;
        for (int j = 0; j < 3; j++) {
            if ((mesh.tag_feature_es[i][j] >= 0 && feature::features[mesh.tag_feature_es[i][j]]->type != "Line")
                && (mesh.tag_feature_es[i][(j + 1) % 3] >= 0
                    && feature::features[mesh.tag_feature_es[i][(j + 1) % 3]]->type != "Line")
                && (mesh.tag_feature_es[i][(j + 2) % 3] < 0
                    || (mesh.tag_feature_es[i][(j + 2) % 3] >= 0
                        && feature::features[mesh.tag_feature_es[i][(j + 2) % 3]]->type == "Line"))) {
                std::array<int, 2> e = {mesh.tris[i][j], mesh.tris[i][(j + 1) % 3]};
                if (e[0] > e[1])
                    std::swap(e[0], e[1]);
                edges.push_back(e);
                break;
            }
        }
    }
    optimization::vector_unique(edges);
    cout << "edge.size() = " << edges.size() << endl;

    int v_slot_size = std::count(mesh.v_is_removed.begin(), mesh.v_is_removed.end(), true);
    int t_slot_size = std::count(mesh.t_is_removed.begin(), mesh.t_is_removed.end(), true);
    for (auto &e:edges) {
        std::vector<std::array<int, 2>> _;
        optimization::split_an_edge(mesh, e[0], e[1], _);
    }
    for (size_t i = 0; i < mesh.tris.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;
        if (mesh.t_quality[i] < 0)
            mesh.t_quality[i] = optimization::tri_energy(mesh.tri_vertices[mesh.tris[i][0]].posf,
                                                         mesh.tri_vertices[mesh.tris[i][1]].posf,
                                                         mesh.tri_vertices[mesh.tris[i][2]].posf);
    }
}

void triwild::feature::subdivide_into_3(MeshData& mesh) {
    for (int i = 0; i < mesh.tris.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;
        if (!(mesh.tag_feature_es[i][0] >= 0 && feature::features[mesh.tag_feature_es[i][0]]->type != "Line")
            || !(mesh.tag_feature_es[i][1] >= 0 && feature::features[mesh.tag_feature_es[i][1]]->type != "Line")
            || !(mesh.tag_feature_es[i][2] >= 0 && feature::features[mesh.tag_feature_es[i][2]]->type != "Line"))
            continue;

        auto &t = mesh.tris[i];
        Point_2f pf = (mesh.tri_vertices[t[0]].posf + mesh.tri_vertices[t[1]].posf + mesh.tri_vertices[t[2]].posf) / 3;
        Point_2 p(pf[0], pf[1]);
        for (int j = 0; j < 3; j++) {
            if (!optimization::is_valid_inversion(pf, mesh.tri_vertices[t[(j + 1) % 3]].posf,
                                                  mesh.tri_vertices[t[(j + 2) % 3]].posf,
                                                  p, mesh.tri_vertices[t[(j + 1) % 3]].pos,
                                                  mesh.tri_vertices[t[(j + 2) % 3]].pos)) {
                p = (mesh.tri_vertices[t[0]].pos + mesh.tri_vertices[t[1]].pos + mesh.tri_vertices[t[2]].pos) / 3;
                pf[0] = p[0].to_double();
                pf[1] = p[1].to_double();
                break;
            }
        }

        mesh.tri_vertices.push_back(TriVertex());
        mesh.tri_vertices.back().pos = p;
        mesh.tri_vertices.back().posf = pf;
        mesh.v_is_removed.push_back(false);
        int v_id = mesh.tri_vertices.size()-1;
        for (int j = 0; j < 3; j++) {
            mesh.tris.push_back({{v_id, t[(j + 1) % 3], t[(j + 2) % 3]}});
            mesh.t_is_removed.push_back(false);
            mesh.is_boundary_es.push_back({{mesh.is_boundary_es[i][j], false, false}});
            mesh.is_bbox_es.push_back({{mesh.is_bbox_es[i][j], false, false}});
            mesh.tag_feature_es.push_back({{mesh.tag_feature_es[i][j], -1, -1}});
            mesh.tag_secondary_feature_es.push_back({{mesh.tag_secondary_feature_es[i][j], -1, -1}});
            mesh.t_quality.push_back(optimization::tri_energy(pf, mesh.tri_vertices[t[(j + 1) % 3]].posf,
                                                              mesh.tri_vertices[t[(j + 2) % 3]].posf));
            mesh.tri_vertices[t[j]].conn_tris.erase(i);
            mesh.tri_vertices[t[(j + 1) % 3]].conn_tris.insert(mesh.tris.size() - 1);
            mesh.tri_vertices[t[(j + 2) % 3]].conn_tris.insert(mesh.tris.size() - 1);
            mesh.tri_vertices.back().conn_tris.insert(mesh.tris.size() - 1);
        }
        mesh.t_is_removed[i]=true;
    }
}

#include <cassert>
void triwild::feature::add_nodes(MeshData& mesh) {
//    for(int i=0;i<mesh.tri_vertices.size();i++) {
//        if (mesh.v_is_removed[i])
//            continue;
//        for (auto &info:mesh.tri_vertices[i].feature_infos) {
//            assert(info[1]>=0 && info[1]<=1);
//            Point_2f p = feature::features[info[0]]->eval(info[1]);
//            double l = (p - mesh.tri_vertices[i].posf).length_2();
//            if (l > 1e-8) {
//                cout << "(p-mesh.tri_vertices[i].posf).length_2()>1e-8" << endl;
//                cout << "p: " << p << endl;
//                cout << "posf: " << mesh.tri_vertices[i].posf << endl;
//                cout << "v " << i << ", feature " << info[0] << " " << feature::features[info[0]]->type << endl;
//
////                feature::features[info[0]]->print_info();
////                cout<<feature::features[info[0]]->to_maple()<<endl;
//
////                mesh.tri_vertices[i].posf = p;
//                optimization::pausee();
//            }
//        }
//    }

    mesh.tri_nodes.resize(mesh.tris.size());

    for (int i = 0; i < mesh.tris.size(); i++) {
        if (mesh.t_is_removed[i])
            continue;

        //check degree
        int degree = 1;
        for (int j = 0; j < 3; j++) {
            if (mesh.tag_feature_es[i][j] >= 0 && features[mesh.tag_feature_es[i][j]]->degree > degree)
                degree = features[mesh.tag_feature_es[i][j]]->degree;
            if (degree == 3)
                break;
        }
        if (degree == 1)
            continue;

        //push boundary nodes
        for (int k = 0; k < 3; k++) {
            int j = (2 + k) % 3;
            int feature_id = mesh.tag_feature_es[i][j];
            TriVertex &v1 = mesh.tri_vertices[mesh.tris[i][(j + 1) % 3]];
            TriVertex &v2 = mesh.tri_vertices[mesh.tris[i][(j + 2) % 3]];
            for (double d = 1; d < degree; d++) {
                if (feature_id < 0) {
                    mesh.nodes.push_back(v1.posf * (1 - d / degree) + v2.posf * (d / degree));
                } else {
                    double t1 = v1.get_t(feature_id);
                    double t2 = v2.get_t(feature_id);
//                    Point_2f p = v1.posf * (1 - d / degree) + v2.posf * (d / degree);
//                    double t = features[feature_id]->inv_eval(p, t1 * (1 - d / degree) + t2 * (d / degree));
//                    mesh.nodes.push_back(features[feature_id]->eval(t));
                    double t = t1 * (1 - d / degree) + t2 * (d / degree);
                    assert(t >= 0 && t <= 1);
                    mesh.nodes.push_back(features[feature_id]->eval(t1 * (1 - d / degree) + t2 * (d / degree)));
                }
                mesh.tri_nodes[i].push_back(mesh.nodes.size() - 1);
            }
        }

        //push center node
        if (degree == 3) {
            Point_2f p(0, 0);
            for (int r = 0; r < 3; r++)
                p = p + mesh.tri_vertices[mesh.tris[i][r]].posf;
//            std::vector<Point_2f> nodes(6);
            for (int r = 0; r < 6; r++) {
                p = p + mesh.nodes[mesh.nodes.size() - 1 - r];
//                nodes[5 - r] = mesh.nodes[mesh.nodes.size() - 1 - r];
            }
            p = p / 9;
//            Point_2f new_p = p;
//            CurvedTriUntangler::untangle_center({{mesh.tri_vertices[mesh.tris[i][0]].posf,
//                                                         mesh.tri_vertices[mesh.tris[i][1]].posf,
//                                                         mesh.tri_vertices[mesh.tris[i][2]].posf}}, nodes, p, new_p);
//            mesh.nodes.push_back(new_p);
            mesh.nodes.push_back(p);
            mesh.tri_nodes[i].push_back(mesh.nodes.size() - 1);
        }
    }

    mesh.n_is_removed = std::vector<bool>(mesh.nodes.size(), false);
}


bool triwild::feature::is_valid_inversion(const std::array<Point_2f, 3>& ps, const std::vector<Point_2f>& ns){
    return DeterminantChecker::is_positive(ps, ns);
}

void triwild::feature::get_new_nodes(const std::array<int, 3>& tag_feature_e, const std::array<TriVertex, 3>& vs,
        std::vector<Point_2f>& new_nodes) {
    int degree = 1;
    for (int j = 0; j < 3; j++) {
        if (tag_feature_e[j] >= 0 && features[tag_feature_e[j]]->degree > degree)
            degree = features[tag_feature_e[j]]->degree;
        if (degree == 3)
            break;
    }
    if (degree == 1)
        return;

    //push boundary nodes
    for (int k = 0; k < 3; k++) {
        int j = (2 + k) % 3;
        int feature_id = tag_feature_e[j];
        const TriVertex &v1 = vs[(j + 1) % 3];
        const TriVertex &v2 = vs[(j + 2) % 3];
        for (double d = 1; d < degree; d++) {
            if (feature_id < 0) {
                new_nodes.push_back(v1.posf * (1 - d / degree) + v2.posf * (d / degree));
            } else {
                double t1 = v1.get_t(feature_id);
                double t2 = v2.get_t(feature_id);
                new_nodes.push_back(features[feature_id]->eval(t1 * (1 - d / degree) + t2 * (d / degree)));
            }
        }
    }

    //push center node
    if (degree == 3) {
        Point_2f p(0, 0);
        for (int r = 0; r < 3; r++)
            p = p + vs[r].posf;
        for (int r = 0; r < 6; r++)
            p = p + new_nodes[new_nodes.size() - 1 - r];
        p = p / 9;
        new_nodes.push_back(p);
    }
}

std::array<double, 2> triwild::feature::json2d2stdarray(const json& a){
    return std::array<double, 2>({(double)a[0], (double)a[1]});
}

Point_2f triwild::feature::json2point(const json& a){
    return Point_2f((double)a[0], (double)a[1]);
}

std::vector<double> triwild::feature::json2d2stdvector(const json& a){
    std::vector<double> r;
    if (a.size() == 0)
        return r;
    r.reserve(a.size()*a[0].size());
    for (int i = 0; i < a.size(); i++)
        for (int j = 0; j < a[i].size(); j++)
            r.push_back((double)a[i][j]);
    return r;
}

std::vector<double> triwild::feature::json1d2stdvector(const json& a){
    std::vector<double> r;
    if (a.size() == 0)
        return r;
    r.reserve(a.size());
    for (int i = 0; i < a.size(); i++)
        r.push_back((double)a[i]);
    return r;
}

ControlVector triwild::feature::json2d2ctrlvector(const json& a){
    ControlVector r;
    if (a.size() == 0)
        return r;
    r.resize(a.size()*a[0].size());
    int cnt = 0;
    for (int i = 0; i < a.size(); i++)
        for (int j = 0; j < a[i].size(); j++)
            r(cnt++) = a[i][j];
    r.conservativeResize(cnt);
    return r;
}

ControlVector triwild::feature::json1d2ctrlvector(const json& a){
    ControlVector r;
    if (a.size() == 0)
        return r;
    r.resize(a.size());
    for (int i = 0; i < a.size(); i++)
        r(i) = a[i];
    return r;
}

//#include <igl/opengl/glfw/Viewer.h>
// void triwild::feature::visualize_features(MeshData& mesh) {
//     igl::opengl::glfw::Viewer viewer;

//     cout<< viewer.data().point_size<<endl;
//     viewer.data().point_size = 10;

//     int feature_id = 0;
//     bool is_show_input = false;
//     bool is_all = true;

//     auto set_edge = [&]() -> void {
//         Eigen::MatrixXd VE;
//         Eigen::MatrixXi E;
//         if (is_show_input) {
//             if (is_all) {
//                 int cnt_e = 0;
//                 for (auto &feature:features)
//                     cnt_e += feature->v_ids.size() - 1;
//                 VE.resize(cnt_e * 2, 3);
//                 E.resize(cnt_e, 2);
//                 int cnt = 0;
//                 for (auto &feature:features) {
//                     for (int i = 0; i < feature->v_ids.size() - 1; i++) {
//                         VE(cnt * 2, 0) = mesh.tri_vertices[feature->v_ids[i]].posf[0];
//                         VE(cnt * 2, 1) = mesh.tri_vertices[feature->v_ids[i]].posf[1];
//                         VE(cnt * 2, 2) = 0;
//                         VE(cnt * 2 + 1, 0) = mesh.tri_vertices[feature->v_ids[i + 1]].posf[0];
//                         VE(cnt * 2 + 1, 1) = mesh.tri_vertices[feature->v_ids[i + 1]].posf[1];
//                         VE(cnt * 2 + 1, 2) = 0;
//                         E(cnt, 0) = cnt * 2;
//                         E(cnt, 1) = cnt * 2 + 1;
//                         cnt++;
//                     }

//                     const auto &front = mesh.tri_vertices[feature->v_ids.front()].posf;
//                     const auto &back = mesh.tri_vertices[feature->v_ids.back()].posf;
//                     viewer.data().add_points(Eigen::RowVector3d(front[0], front[1], 0), Eigen::RowVector3d(feature->is_inflection[0]?0:1,feature->is_inflection[0]?1:0,0));
//                     viewer.data().add_points(Eigen::RowVector3d(back[0], back[1], 0), Eigen::RowVector3d(feature->is_inflection[1]?0:1,feature->is_inflection[1]?1:0,0));
//                 }
//             } else {
//                 auto feature = features[feature_id];
//                 VE.resize((feature->v_ids.size() - 1) * 2, 3);
//                 E.resize(feature->v_ids.size() - 1, 2);
//                 int cnt = 0;
//                 for (int i = 0; i < feature->v_ids.size() - 1; i++) {
//                     VE(cnt * 2, 0) = mesh.tri_vertices[feature->v_ids[i]].posf[0];
//                     VE(cnt * 2, 1) = mesh.tri_vertices[feature->v_ids[i]].posf[1];
//                     VE(cnt * 2, 2) = 0;
//                     VE(cnt * 2 + 1, 0) = mesh.tri_vertices[feature->v_ids[i + 1]].posf[0];
//                     VE(cnt * 2 + 1, 1) = mesh.tri_vertices[feature->v_ids[i + 1]].posf[1];
//                     VE(cnt * 2 + 1, 2) = 0;
//                     E(cnt, 0) = cnt * 2;
//                     E(cnt, 1) = cnt * 2 + 1;
//                     cnt++;
//                 }

//                 const auto &front = mesh.tri_vertices[feature->v_ids.front()].posf;
//                 const auto &back = mesh.tri_vertices[feature->v_ids.back()].posf;
//                 viewer.data().add_points(Eigen::RowVector3d(front[0], front[1], 0), Eigen::RowVector3d(1,0,0));
//                 viewer.data().add_points(Eigen::RowVector3d(back[0], back[1], 0), Eigen::RowVector3d(1,0,0));
//             }
//         } else {
//             std::vector<std::array<int, 2>> edges;
//             for (int i = 0; i < mesh.tris.size(); i++) {
//                 if (mesh.t_is_removed[i])
//                     continue;
//                 for (int j = 0; j < 3; j++) {
//                     if ((is_all && mesh.tag_feature_es[i][j] >= 0)
//                         || (!is_all && mesh.tag_feature_es[i][j] == feature_id)) {
//                         std::array<int, 2> e = {mesh.tris[i][(j + 1) % 3], mesh.tris[i][(j + 2) % 3]};
//                         if (e[0] > e[1])
//                             std::swap(e[0], e[1]);
//                         edges.push_back(e);
//                     }
//                 }
//             }
//             std::sort(edges.begin(), edges.end());
//             edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

//             VE.resize(edges.size() * 2, 3);
//             E.resize(edges.size(), 2);
//             for (int i = 0; i < edges.size(); i++) {
//                 for (int j = 0; j < 2; j++) {
//                     VE.row(i * 2 + j)
//                             << mesh.tri_vertices[edges[i][j]].posf[0], mesh.tri_vertices[edges[i][j]].posf[1], 0;
//                     E(i, j) = i * 2 + j;
//                 }
//             }
//         }
//         Eigen::MatrixXd C(1, 3);
//         C << 1, 0, 0;
//         viewer.data().set_edges(VE, E, C);
//     };

//     viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
//         if (key == '=') {
//             feature_id++;
//             if (feature_id == features.size())
//                 feature_id = 0;
//             cout << "showing feature " << feature_id << endl;
//             set_edge();
//         } else if (key == '-') {
//             feature_id--;
//             if (feature_id == -1)
//                 feature_id = features.size() - 1;
//             cout << "showing feature " << feature_id << endl;
//             set_edge();
//         } else if (key == 'C') {
//             is_show_input = !is_show_input;
//             cout << "switch to is_show_input = "<<is_show_input<< endl;
//             set_edge();
//         } else if (key == 'A'){
//             is_all = !is_all;
//             cout << "switch to is_all = "<<is_all<< endl;
//             set_edge();
//         }
//         return false;
//     };

//     Eigen::MatrixXd V;
//     Eigen::MatrixXi F;

//     std::unordered_map<int, int> map_v_ids;
//     int cnt = 0;
//     for (size_t i = 0; i < mesh.tri_vertices.size(); i++) {
//         if (mesh.v_is_removed[i])
//             continue;
//         map_v_ids[i] = cnt++;
//     }
//     V.resize(cnt, 3);
//     cnt = 0;
//     for (size_t i = 0; i < mesh.tri_vertices.size(); i++) {
//         if (mesh.v_is_removed[i])
//             continue;
//         V(cnt, 0) = mesh.tri_vertices[i].posf[0];
//         V(cnt, 1) = mesh.tri_vertices[i].posf[1];
//         V(cnt, 2) = 0;
//         cnt++;
//     }

//     cnt = std::count(mesh.t_is_removed.begin(), mesh.t_is_removed.end(), false);
//     F.resize(cnt, 3);
//     cnt = 0;
//     for (size_t i = 0; i < mesh.tris.size(); i++) {
//         if (mesh.t_is_removed[i])
//             continue;
//         for (int j = 0; j < 3; j++)
//             F(cnt, j) = map_v_ids[mesh.tris[i][j]];
//         cnt++;
//     }

//     set_edge();

//     viewer.core.background_color << 1.0f, 1.0f, 1.0f, 0.0f;
//     viewer.data().show_overlay = true;
//     viewer.core.set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_TRACKBALL);
//     viewer.data().set_face_based(true);
//     viewer.data().set_mesh(V, F);
//     viewer.launch();
// }

#include <igl/writeSTL.h>
void triwild::feature::output_input_features(const std::vector<std::shared_ptr<FeatureElement>>& features,
        Eigen::MatrixXd& V, std::vector<std::array<int, 2>>& edges, const std::string& postfix){
    Eigen::MatrixXd V_out(V.rows(), 3);
    for (int i = 0; i < V.rows(); i++)
        V_out.row(i) << V(i, 0), V(i, 1), 0;
    int cnt = 0;
    for (int i = 0; i < features.size(); i++)
        cnt += features[i]->v_ids.size()-1;
    Eigen::MatrixXi E_out(cnt, 3);
    cnt = 0;
    for (int i = 0; i < features.size(); i++) {
        for (int j = 0; j < features[i]->v_ids.size() - 1; j++)
            E_out.row(cnt++) << features[i]->v_ids[j], features[i]->v_ids[j + 1], features[i]->v_ids[j + 1];
    }
    igl::writeSTL(args.output + "_"+postfix+".stl", V_out, E_out);









        std::ofstream os(args.output + "_"+postfix+".eps");
        if(os.fail())
            return;

        os << std::setprecision(16);

        os <<"%!PS-Adobe-3.0 EPSF-3.0\n";
        os << "%%BoundingBox: "<< (args.box_min(0)-2*args.target_edge_len) << " " << (args.box_min(1)-args.target_edge_len) << " " << (args.box_max(0)+2*args.target_edge_len) <<" "<< (args.box_max(1)+args.target_edge_len) <<"\n\n";
        os << "%%Pages: 1\n";
        os << "%%Page: 1 1\n";
        os << "/show-ctr {\ndup stringwidth pop\n -2 div 0\n rmoveto show\n} def\n\n 2 setlinejoin\n\n";


        os<<1<<" setlinewidth\n\n";

        for (const auto &f : features)
            os << f->to_eps();

        os.close();
}