// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "triangulation.h"
#include "optimization.h"
#include "feature_preprocessing.h"
#include "feature.h"
#include <igl/unique_rows.h>
#include <igl/writeSTL.h>

bool triwild::triangulation::load_input(const std::string& input, Eigen::MatrixXd& V, std::vector<std::array<int, 2>>& edges) {
    std::ifstream is(input, std::ios::in);
    if (is.fail())
        throw std::runtime_error("Unable to open input file \"" + input + "\"!");

    std::vector<std::array<double, 2>> vertices;
    int i = 0;
    std::string line_str;
    double tmp;
    while (std::getline(is, line_str)) {
        std::istringstream line(line_str);

        std::string prefix;
        line >> prefix;

        if (prefix == "v") {
            vertices.emplace_back();
            auto &v = vertices.back();
            line >> v[0] >> v[1] >> tmp;
        } else if (prefix == "l") {
            edges.emplace_back();
            auto &e = edges.back();
            line >> e[0] >> e[1];
            e[0]--;
            e[1]--;
        }
    }

    V.resize(vertices.size(), 2);
    for (int i = 0; i < vertices.size(); ++i) {
        V.row(i) << vertices[i][0], vertices[i][1];
    }

    if (V.rows() == 0)
        return false;
    else {
        for (int i = 0; i < edges.size(); i++) {
            if (edges[i][0] < 0 || edges[i][0] >= V.rows()
                || edges[i][1] < 0 || edges[i][1] >= V.rows())
                return false;
        }
    }
    return true;
}

void triwild::triangulation::preprocessing(Eigen::MatrixXd& V, std::vector<std::array<int, 2>>& edges, GEO::Mesh& b_mesh) {
    auto build_b_mesh = [](const Eigen::MatrixXd &V, const std::vector<std::array<int, 2>> &edges,
                           GEO::Mesh &b_mesh) -> void {
        b_mesh.vertices.clear();
        b_mesh.vertices.create_vertices(V.rows());
        for (int i = 0; i < V.rows(); i++) {
            GEO::vec3 &p = b_mesh.vertices.point(i);
            p[0] = V(i, 0);
            p[1] = V(i, 1);
            p[2] = 0;
        }
        b_mesh.facets.clear();
        b_mesh.facets.create_triangles(edges.size());
        for (int i = 0; i < edges.size(); i++) {
            b_mesh.facets.set_vertex(i, 0, edges[i][0]);
            b_mesh.facets.set_vertex(i, 1, edges[i][1]);
            b_mesh.facets.set_vertex(i, 2, edges[i][1]);
        }
        b_mesh.facets.compute_borders();
    };

    ////global parameters
    auto set_global_parameters = [&]() {
        args.box_max = V.colwise().maxCoeff();
        args.box_min = V.colwise().minCoeff();
        args.diagonal_len = (args.box_max - args.box_min).norm();
        if (args.target_edge_len == -1)
            args.target_edge_len = args.diagonal_len * args.edge_length_r;//args.i_edge_length == 20 in init
        double input_envelop = args.epsilon * args.diagonal_len;
        double dd = input_envelop;
        args.envelope = input_envelop - dd / 2;//sampling error = dd / 2
        args.sampling_density = dd;

        if(!args.is_preserving_feature)
            args.min_edge_length = args.envelope / args.target_edge_len;
        else
            args.min_edge_length = (args.diagonal_len * 1e-4 * 2) / args.target_edge_len;

//        if(args.stage > 1)
//            args.envelope /= 2;
        args.envelope /= args.stage;
    };

    ////remove duplicates & degenerates
    int old_cnt_v = V.rows();
    int old_cnt_e = edges.size();

//    for (auto &feature:feature::features) {
//        for (int i = 0; i < feature->paras.size(); i++) {
//            Point_2f p = feature->eval(feature->paras[i]);
//            V.row(feature->v_ids[i]) << p[0], p[1];
//        }
//    }

    //duplicate vertices
    Eigen::VectorXi VI, _;
    Eigen::MatrixXd V_tmp;
    igl::unique_rows(V, V_tmp, _, VI);
    V = V_tmp;

    //duplicate/degenerate edges
    for (int i = 0; i < edges.size(); i++) {
        auto &e = edges[i];
        edges[i][0] = VI(edges[i][0]);
        edges[i][1] = VI(edges[i][1]);
        if (e[0] == e[1]) {
            edges.erase(edges.begin() + i);
            i--;
        } else if (e[0] > e[1])
            std::swap(e[0], e[1]);
    }
    optimization::vector_unique(edges);

    cout << "remove duplicates:" << endl;
    cout << "#v " << old_cnt_v << "->" << V.rows() << endl;
    cout << "#e " << old_cnt_e << "->" << edges.size() << endl;
    old_cnt_v = V.rows();
    old_cnt_e = edges.size();

    //update feature
    if (args.is_preserving_feature) {
        //remove degenerate feature edges
        for (int i = 0; i < feature::features.size(); i++) {
            for (int j = 0; j < feature::features[i]->v_ids.size(); j++) {
                feature::features[i]->v_ids[j] = VI(feature::features[i]->v_ids[j]);
                if (j > 0 && feature::features[i]->v_ids[j] == feature::features[i]->v_ids[j - 1]) {
                    feature::features[i]->v_ids.erase(feature::features[i]->v_ids.begin() + j);
                    feature::features[i]->paras.erase(feature::features[i]->paras.begin() + j);
                    j--;
                    continue;
                }
            }
            if (feature::features[i]->v_ids.size() < 2) {
                feature::features.erase(feature::features.begin() + i);
                i--;
                continue;
            }
        }

        set_global_parameters();
        feature::preprocessing(V, edges);
        cout << "refine:" << endl;
        cout << "#v " << old_cnt_v << "->" << V.rows() << endl;
        cout << "#e " << old_cnt_e << "->" << edges.size() << endl;
        old_cnt_v = V.rows();
        old_cnt_e = edges.size();

        set_global_parameters();
        build_b_mesh(V, edges, b_mesh);

        GEO::MeshFacetsAABB b_tree(b_mesh);
        triangulation::simplify_input(V, edges, b_tree);
        cout << "simplify:" << endl;
        cout << "#v " << old_cnt_v << "->" << V.rows() << endl;
        cout << "#e " << old_cnt_e << "->" << edges.size() << endl;
        cout << "#secondary_features = " << feature::secondary_features.size() << endl;
    } else {
        set_global_parameters();
        build_b_mesh(V, edges, b_mesh);
        
        GEO::MeshFacetsAABB b_tree(b_mesh);
        triangulation::simplify_input(V, edges, b_tree);
        cout << "simplify:" << endl;
        cout << "#v " << old_cnt_v << "->" << V.rows() << endl;
        cout << "#e " << old_cnt_e << "->" << edges.size() << endl;
    }
}

void triwild::triangulation::simplify_input(Eigen::MatrixXd& V, std::vector<std::array<int, 2>>& edges, GEO::MeshFacetsAABB &b_tree) {
//    Eigen::MatrixXd oV1(V.rows(), 3), _1;
//    Eigen::MatrixXi oE1(edges.size(), 3);
//    for (int i = 0; i < V.rows(); i++)
//        oV1.row(i) << V(i, 0), V(i, 1), 0;
//    for (int i = 0; i < edges.size(); i++)
//        oE1.row(i) << edges[i][0], edges[i][1], edges[i][1];
//    igl::writeSTL(args.output+"_before_edges.stl", oV1, oE1, _1);

    using namespace optimization;

    std::vector<std::unordered_set<int>> conn_es(V.rows());
    for (int i = 0; i < edges.size(); i++) {
        conn_es[edges[i][0]].insert(i);
        conn_es[edges[i][1]].insert(i);
    }

    ////feature
    std::vector<bool> is_freezed(V.rows(), false);
    for (int i = 0; i < feature::features.size(); i++) {
        for (int j = 0; j < feature::features[i]->v_ids.size(); j++)
            is_freezed[feature::features[i]->v_ids[j]] = true;
    }
//    std::vector<int> tag_secondary_feature_vs(V.rows(), -1);
//    for (int i = 0; i < feature::secondary_features.size(); i++) {
//        for (int j = 0; j < feature::secondary_features[i]->v_ids.size(); j++)
//            tag_secondary_feature_vs[feature::secondary_features[i]->v_ids[j]] = i;
////        is_freezed[feature::secondary_features[i]->v_ids.front()] = true;
////        is_freezed[feature::secondary_features[i]->v_ids.back()] = true;
//    }
    std::vector<std::unordered_set<int>> tag_secondary_feature_vs(V.rows());
    for (int i = 0; i < feature::secondary_features.size(); i++) {
        for (int j = 0; j < feature::secondary_features[i]->v_ids.size(); j++)
            tag_secondary_feature_vs[feature::secondary_features[i]->v_ids[j]].insert(i);
//        is_freezed[feature::secondary_features[i]->v_ids.front()] = true;
//        is_freezed[feature::secondary_features[i]->v_ids.back()] = true;
    }
//    for(int i=0;i<tag_secondary_feature_vs.size();i++) {
//        if (tag_secondary_feature_vs[i].size() > 1)
//            is_freezed[i] = true;
//    }

    auto is_clean_vertex = [&](int v_id) -> bool {
        if (conn_es[v_id].size() != 2)
            return false;
        if (is_freezed[v_id])
            return false;
//        if (tag_secondary_feature_vs[v_id].size() != 1)
//            return false;
        return true;
    };

    std::priority_queue<ElementInQueue_s, std::vector<ElementInQueue_s>, cmp_s> queue_s;
    for (auto &e:edges) {
        double l = (V.row(e[0]) - V.row(e[1])).squaredNorm();
        if (is_clean_vertex(e[0]))
            queue_s.push(ElementInQueue_s(e, l));
        if (is_clean_vertex(e[1])) {
            queue_s.push(ElementInQueue_s({{e[1], e[0]}}, l));
        }
    }

    double epsilon_2 = args.envelope * args.envelope * 0.8 * 0.8;
    double dd = args.sampling_density;
    std::vector<bool> v_is_removed(V.rows(), false);
    std::vector<bool> e_is_removed(edges.size(), false);

    auto update_secondary_feature = [&](int v1_id, int v2_id) -> void {
        if (tag_secondary_feature_vs[v1_id].empty())
            return;

//        int feature_id = tag_secondary_feature_vs[v1_id][0];
        auto tmp_tag_secondary_feature_vs_v1_id = tag_secondary_feature_vs[v1_id];
        for (int feature_id:tmp_tag_secondary_feature_vs_v1_id) {
            auto &v_ids = feature::secondary_features[feature_id]->v_ids;
            auto &paras = feature::secondary_features[feature_id]->paras;
            int j = std::find(v_ids.begin(), v_ids.end(), v1_id) - v_ids.begin();
//            if (feature_id == 790) {
//                cout << "v1v2_id " << v1_id << " " << v2_id << endl;
//                cout << "feature " << feature_id << ": ";
//                for (int i:v_ids)
//                    cout << i << " ";
//                cout << endl;
//            }

            if(j == v_ids.size())
                continue;

            if ((j + 1 < v_ids.size() && v_ids[j + 1] == v2_id) || (j - 1 >= 0 && v_ids[j - 1] == v2_id)) {
                v_ids.erase(v_ids.begin() + j);
                paras.erase(paras.begin() + j);
                continue;
            }

            //deal with loop
            int j2 = std::find(v_ids.begin(), v_ids.end(), v2_id) - v_ids.begin();
            if (j2 < v_ids.size()) {
                if (j2 - 0 > v_ids.size() - j2) {
                    for (int k = j2; k < v_ids.size(); k++)
                        tag_secondary_feature_vs[v_ids[k]].erase(feature_id);
                    v_ids.erase(v_ids.begin() + j2, v_ids.end());
                    paras.erase(paras.begin() + j2, paras.end());
                } else {
                    for (int k = 0; k < j2 + 1; k++)
                        tag_secondary_feature_vs[v_ids[k]].erase(feature_id);
                    v_ids.erase(v_ids.begin(), v_ids.begin() + j2 + 1);
                    paras.erase(paras.begin(), paras.begin() + j2 + 1);
                }
                j = std::find(v_ids.begin(), v_ids.end(), v1_id) - v_ids.begin();
//                cout << "kill a loop" << endl;
            }

            v_ids[j] = v2_id;
            tag_secondary_feature_vs[v2_id].insert(feature_id);
        }

//        for (int feature_id = 0; feature_id < feature::secondary_features.size(); feature_id++) {
//            auto &feature = feature::secondary_features[feature_id];
//            for (int v_id:feature->v_ids) {
//                if (v_is_removed[v_id]) {
//                    cout << "feature " << feature_id << ": ";
//                    for (int i:feature->v_ids)
//                        cout << i << " ";
//                    cout << endl;
//                    cout << "v_id " << v_id << endl;
//                    cout << "v1_id " << v1_id << ": ";
//                    for (int i:tag_secondary_feature_vs[v1_id])
//                        cout << i << " ";
//                    cout << endl;
//                    cout << "v2_id " << v2_id << ": ";
//                    for (int i:tag_secondary_feature_vs[v2_id])
//                        cout << i << " ";
//                    cout << endl;
//                }
//                assert(v_is_removed[v_id] == false);
//            }
//        }
    };


    auto check = [&]() -> void{
        std::unordered_map<int, int> map_v_ids;
        int cnt = 0;
        for (int i = 0; i < V.rows(); i++) {
            if (v_is_removed[i])
                continue;
            map_v_ids[i] = cnt;
            cnt++;

            for(int e_id:conn_es[i]){
                if(i!=edges[e_id][0] && i!=edges[e_id][1]){
                    cout<<"invalid conn_es"<<endl;
                    cout<<"v "<<i<<endl;
                    optimization::pausee();
                }
            }
        }

        int old_edge_size = std::count(e_is_removed.begin(), e_is_removed.end(), false);
        std::vector<std::array<int, 2>> new_edges(old_edge_size);
        cnt = 0;
        for (int i = 0; i < edges.size(); i++) {
            if (e_is_removed[i])
                continue;
            new_edges[cnt] = {{map_v_ids[edges[i][0]], map_v_ids[edges[i][1]]}};
            if(new_edges[cnt][0] == new_edges[cnt][1]){
                cout<<"new_edges[cnt][0] == new_edges[cnt][1]"<<endl;
                optimization::pausee();
            }
            if (new_edges[cnt][0] > new_edges[cnt][1])
                std::swap(new_edges[cnt][0], new_edges[cnt][1]);
            cnt++;
        }
        optimization::vector_unique(new_edges);
        if(new_edges.size()!=old_edge_size){
            cout<<"new_edges.size()!=old_edge_size"<<endl;
            cout<<new_edges.size()<<" != "<<old_edge_size<<endl;
            optimization::pausee();
        }
    };

    while (!queue_s.empty()) {
        std::array<int, 2> v_ids = queue_s.top().v_ids;
        double old_l = queue_s.top().weight;
        queue_s.pop();

        if (v_is_removed[v_ids[0]] || v_is_removed[v_ids[1]])
            continue;
        if (!is_clean_vertex(v_ids[0]))
            continue;
        double l = (V.row(v_ids[0]) - V.row(v_ids[1])).squaredNorm();
        if (l != old_l)
            continue;
        if (!queue_s.empty() && v_ids == queue_s.top().v_ids)
            continue;

        ////collapse an edge
        int v1_id = v_ids[0];
        int v2_id = v_ids[1];
        int v_id = -1;
        for (int e_id:conn_es[v1_id]) {
            if (edges[e_id][0] == v2_id || edges[e_id][1] == v2_id)
                continue;
            v_id = edges[e_id][0] == v1_id ? edges[e_id][1] : edges[e_id][0];
            break;
        }
        if (v_id < 0) {
            int e_id = *(conn_es[v1_id].begin());
            conn_es[edges[e_id][0]].erase(e_id);
            conn_es[edges[e_id][1]].erase(e_id);
            e_is_removed[e_id] = true;

            //push new edge
            e_id = *(conn_es[v1_id].begin());
            if (is_clean_vertex(edges[e_id][0]))
                queue_s.push(ElementInQueue_s(edges[e_id], l));
            if (is_clean_vertex(edges[e_id][1]))
                queue_s.push(ElementInQueue_s({{edges[e_id][1], edges[e_id][0]}}, l));

//            cout<<"find a loop"<<endl;
//            check();
            continue;
        }
        double new_l = (V.row(v_id) - V.row(v2_id)).squaredNorm();

        //check envelope
        bool is_valid = true;
//        GEO::vec3 p0(V(v_id, 0), V(v_id, 1), 0);
        GEO::vec3 nearest_point;
        double sq_dist;
        GEO::index_t prev_facet;
        int N = std::sqrt(new_l) / dd + 1;
        for (double n = 1; n < N; n++) {
            double k = n / N;
            GEO::vec3 v(V(v_id, 0) * k + V(v2_id, 0) * (1 - k),
                        V(v_id, 1) * k + V(v2_id, 1) * (1 - k), 0);
            if (n == 1) {
                prev_facet = b_tree.nearest_facet(v, nearest_point, sq_dist);
                if (sq_dist > epsilon_2) {
                    is_valid = false;
                    break;
                }
            } else {
                sq_dist = v.distance2(nearest_point);
                if (sq_dist <= epsilon_2)
                    continue;
            }
            b_tree.nearest_facet_with_hint(v, prev_facet, nearest_point, sq_dist);
            if (sq_dist > epsilon_2) {
                is_valid = false;
                break;
            }
        }
        if (!is_valid)
            continue;

        //collapse
        v_is_removed[v1_id] = true;
        bool is_loop = false;
        for(int e_id:conn_es[v_id]){
            if(edges[e_id][0]==v2_id || edges[e_id][1]==v2_id){
                is_loop = true;
                break;
            }
        }
        if(!is_loop) {
            int changed_e_id;
            for (int e_id:conn_es[v1_id]) {
                if (edges[e_id][0] == v2_id || edges[e_id][1] == v2_id) {//remove
                    e_is_removed[e_id] = true;
                    conn_es[v2_id].erase(e_id);
                } else {//change
                    int j = edges[e_id][0] == v1_id ? 0 : 1;
                    edges[e_id][j] = v2_id;
                    conn_es[v2_id].insert(e_id);
                    changed_e_id = e_id;
                }
            }
            //push new edge
            if (is_clean_vertex(edges[changed_e_id][0]))
                queue_s.push(ElementInQueue_s(edges[changed_e_id], new_l));
            if (is_clean_vertex(edges[changed_e_id][1]))
                queue_s.push(ElementInQueue_s({{edges[changed_e_id][1], edges[changed_e_id][0]}}, new_l));
        } else {
            for (int e_id:conn_es[v1_id]) {
                e_is_removed[e_id] = true;
                if (edges[e_id][0] == v2_id || edges[e_id][1] == v2_id)
                    conn_es[v2_id].erase(e_id);
                if (edges[e_id][0] == v_id || edges[e_id][1] == v_id)
                    conn_es[v_id].erase(e_id);
            }
        }

        update_secondary_feature(v1_id, v2_id);

//        cout<<v_ids[0]<<" "<<v_ids[1]<<endl;
//        check();
    }
//    check();

    //clean up V and edges
    int cnt = std::count(v_is_removed.begin(), v_is_removed.end(), false);
    Eigen::MatrixXd new_V(cnt, 2);
    std::unordered_map<int, int> map_v_ids;
    cnt = 0;
    for (int i = 0; i < V.rows(); i++) {
        if (v_is_removed[i])
            continue;
        new_V.row(cnt) = V.row(i);
        map_v_ids[i] = cnt;
        cnt++;
    }
    V = new_V;

    cnt = std::count(e_is_removed.begin(), e_is_removed.end(), false);
    std::vector<std::array<int, 2>> new_edges(cnt);
    cnt = 0;
    for (int i = 0; i < edges.size(); i++) {
        if (e_is_removed[i])
            continue;
        new_edges[cnt] = {{map_v_ids[edges[i][0]], map_v_ids[edges[i][1]]}};
        cnt++;
    }
    edges = new_edges;

    for (int i = 0; i < feature::features.size(); i++) {
        for (int j = 0; j < feature::features[i]->v_ids.size(); j++)
            feature::features[i]->v_ids[j] = map_v_ids[feature::features[i]->v_ids[j]];
    }
    for (int i = 0; i < feature::secondary_features.size(); i++) {
        if (feature::secondary_features[i]->v_ids.size() < 2) {
            feature::secondary_features.erase(feature::secondary_features.begin() + i);
            i--;
            continue;
        }
        for (int j = 0; j < feature::secondary_features[i]->v_ids.size(); j++) {
//            if (v_is_removed[feature::secondary_features[i]->v_ids[j]]) {
//                feature::secondary_features[i]->v_ids.erase(feature::secondary_features[i]->v_ids.begin() + j);
//                j--;
//                continue;
//            }
            assert(v_is_removed[feature::secondary_features[i]->v_ids[j]] == false);
            feature::secondary_features[i]->v_ids[j] = map_v_ids[feature::secondary_features[i]->v_ids[j]];
        }
//        if (feature::secondary_features[i]->v_ids.size() < 2) {
//            feature::secondary_features.erase(feature::secondary_features.begin() + i);
//            i--;
//            continue;
//        }
    }

//    Eigen::MatrixXd oV(V.rows(), 3), _;
//    Eigen::MatrixXi oE(edges.size(), 3);
//    for (int i = 0; i < V.rows(); i++)
//        oV.row(i) << V(i, 0), V(i, 1), 0;
//    for (int i = 0; i < edges.size(); i++)
//        oE.row(i) << edges[i][0], edges[i][1], edges[i][1];
//    igl::writeSTL(args.output+"_edges.stl", oV, oE, _);
}

#include <geogram/delaunay/delaunay.h>
void triwild::triangulation::BSP_subdivision(const Eigen::MatrixXd& V, const std::vector<std::array<int, 2>>& edges,
        MeshData& mesh, std::vector<std::vector<int>>& tag_boundary_es) {
    ////delaunay
    Eigen::MatrixXd new_V = V;
    const int n = new_V.rows();
    new_V.conservativeResize(n + 4, 2);
    Eigen::RowVector2d min(args.box_min(0) - args.target_edge_len, args.box_min(1) - args.target_edge_len);
    Eigen::RowVector2d max(args.box_max(0) + args.target_edge_len, args.box_max(1) + args.target_edge_len);
    new_V.row(n) = min;
    new_V.row(n + 1) << min(0), max(1);
    new_V.row(n + 2) = max;
    new_V.row(n + 3) << max(0), min(1);

    GEO::Delaunay::initialize();
    GEO::Delaunay_var delaunay2d = GEO::Delaunay::create(2, "BDEL2d");
//    delaunay2d->set_vertices(new_V.rows(), new_V.data());
    auto *flatten = new double[new_V.rows() * 2];
    for (int i = 0; i < new_V.rows(); i++) {
        flatten[i * 2] = new_V(i, 0);
        flatten[i * 2 + 1] = new_V(i, 1);
    }
    delaunay2d->set_vertices(new_V.rows(), flatten);
    delete[] flatten;
    auto cell_adj = delaunay2d->cell_to_v();

    ////init bsp tree
    std::vector<std::vector<int>> bsp_faces;
    bsp_faces.resize(delaunay2d->nb_cells());
    for (int i = 0; i < delaunay2d->nb_cells(); i++)
        bsp_faces[i] = {{cell_adj[i * 3], cell_adj[i * 3 + 1], cell_adj[i * 3 + 2]}};

    std::vector<int> cnt_conn_es(V.rows(), 0);
    for (int i = 0; i < edges.size(); i++) {
        cnt_conn_es[edges[i][0]]++;
        cnt_conn_es[edges[i][1]]++;
    }
    auto &bsp_vertices = mesh.tri_vertices;
    bsp_vertices.resize(n + 4);
    for (int i = 0; i < new_V.rows(); i++) {
        bsp_vertices[i].pos[0] = new_V(i, 0);
        bsp_vertices[i].pos[1] = new_V(i, 1);
        if (i >= n)
            bsp_vertices[i].is_freezed = true;
        else if (cnt_conn_es[i] == 1) {
            bsp_vertices[i].is_on_point = true;
            bsp_vertices[i].input_posf.set(new_V(i, 0), new_V(i, 1));
        }
    }

    std::vector<std::unordered_set<int>> conn_fs(bsp_vertices.size());
    for (int i = 0; i < bsp_faces.size(); i++) {
        for (int j = 0; j < 3; j++)
            conn_fs[bsp_faces[i][j]].insert(i);
    }

    tag_boundary_es.resize(bsp_vertices.size());

//    Eigen::MatrixXd dV(bsp_vertices.size(), 3), _1;
//    Eigen::MatrixXi dF(bsp_faces.size(), 3);
//    for (int i = 0; i < bsp_vertices.size(); i++)
//        dV.row(i) << bsp_vertices[i].pos[0].to_double(), bsp_vertices[i].pos[1].to_double(), 0;
//    for (int i = 0; i < bsp_faces.size(); i++)
//        dF.row(i) << bsp_faces[i][0], bsp_faces[i][1], bsp_faces[i][2];
//    igl::writeSTL(args.output+"_delaunay.stl", dV, dF, _1);

    int old_cnt_v = bsp_vertices.size();
    int old_cnt_f = bsp_faces.size();

    ////bsp cutting
    auto split_a_face = [&](int r_v_id, int cnt, int n_size, int f_id) -> int {
        std::vector<int> new_face;
        for (int k = 0; k <= cnt; k++)//p1->
            new_face.push_back(bsp_faces[f_id][(r_v_id + k) % n_size]);
        bsp_faces.push_back(new_face);
        int new_f_id = bsp_faces.size() - 1;

        new_face.clear();
        for (int k = cnt; k <= n_size; k++)//->p2
            new_face.push_back(bsp_faces[f_id][(r_v_id + k) % n_size]);
        bsp_faces[f_id] = new_face;

        for (int v_id:bsp_faces[new_f_id]) {
            conn_fs[v_id].erase(f_id);
            conn_fs[v_id].insert(new_f_id);
        }
        for (int v_id:bsp_faces[f_id]) {
            conn_fs[v_id].insert(f_id);
        }

        return new_f_id;
    };

    TriVertex intersection_v;
    for (int e_id = 0; e_id < edges.size(); e_id++) {
//        cout<<"e_id = "<<e_id<<endl;
        int v1_id = edges[e_id][0];
        int v2_id = edges[e_id][1];
        tag_boundary_es[v1_id].push_back(e_id);

        std::vector<int> tmp = optimization::set_intersection(conn_fs[v1_id], conn_fs[v2_id]);
        if (tmp.size() > 0) { //edge exist in the init mesh
            tag_boundary_es[v2_id].push_back(e_id);
            continue;
        }

//        cout<<v1_id<<" "<<v2_id<<endl;

        assert(v1_id != v2_id);
        assert(bsp_vertices[v1_id].pos != bsp_vertices[v2_id].pos);

        int v_id = v1_id;
        std::unordered_set<int> tmp_conn_fs = conn_fs[v1_id];
        while (true) {
            bool is_finished = false;
            bool is_intersected = false;
            for (int f_id: tmp_conn_fs) {
                is_intersected = false;
                int n_size = bsp_faces[f_id].size();
                int r_v_id = std::find(bsp_faces[f_id].begin(), bsp_faces[f_id].end(), v_id) - bsp_faces[f_id].begin();
                for (int i = 1; i < n_size - 1; i++) {
                    const int j = (r_v_id + i) % n_size;

                    int p1_id = bsp_faces[f_id][j];
                    int p2_id = bsp_faces[f_id][(j + 1) % n_size];
                    if (p1_id == v_id || p2_id == v_id)
                        continue;
                    bool is_cross_p1 = false;
                    bool is_cross_p2 = false;
                    if (!segment_intersection(mesh, v_id, v2_id, p1_id, p2_id, is_cross_p1, is_cross_p2,
                                              intersection_v))
                        continue;

//                    cout<<"old v_id = "<<v_id<<endl;
//                    cout<<"p_ids = "<<p1_id<<" "<<p2_id<<endl;
//                    cout<<"intersect face "<<f_id<<endl;
//                    for(int f:bsp_faces[f_id])
//                        cout<<f<<" ";
//                    cout<<endl;

                    if (is_cross_p1) {
//                        cout<<">>>is_cross_p1"<<endl;
                        if (p1_id == v2_id) {
                            is_finished = true;
                        }
                        tag_boundary_es[p1_id].push_back(e_id);

                        tmp_conn_fs = conn_fs[p1_id];
//                        tmp_conn_fs.erase(f_id);
                        std::vector<int> n_f_ids = optimization::set_intersection(conn_fs[v_id], conn_fs[p1_id]);
                        for (int n_f_id:n_f_ids)
                            tmp_conn_fs.erase(n_f_id);
                        if (i != 1) {
                            int new_f_id = split_a_face(r_v_id, i, n_size, f_id);
//                            std::vector<int> n_f_ids = optimization::set_intersection(conn_fs[v_id], conn_fs[p1_id]);
//                            assert(n_f_ids.size()==2);
//                            for (int n_f_id:n_f_ids)
//                                tmp_conn_fs.erase(n_f_id);

//                            cout<<"f_id: ";
//                            for(int f:bsp_faces[new_f_id])
//                                cout<<f<<" ";
//                            cout<<endl;
//                            cout<<"new_f_id: ";
//                            for(int f:bsp_faces[f_id])
//                                cout<<f<<" ";
//                            cout<<endl;
                        }
                        v_id = p1_id;
//                        optimization::pausee();
                    } else if (is_cross_p2) {
//                        cout<<">>>is_cross_p2"<<endl;
                        if (p2_id == v2_id) {
                            is_finished = true;
                        }
                        tag_boundary_es[p2_id].push_back(e_id);

                        tmp_conn_fs = conn_fs[p2_id];
//                        tmp_conn_fs.erase(f_id);
                        std::vector<int> n_f_ids = optimization::set_intersection(conn_fs[v_id], conn_fs[p2_id]);
                        for (int n_f_id:n_f_ids)
                            tmp_conn_fs.erase(n_f_id);
                        if (i != n_size - 2) {
                            int new_f_id = split_a_face(r_v_id, i + 1, n_size, f_id);
//                            std::vector<int> n_f_ids = optimization::set_intersection(conn_fs[v_id], conn_fs[p2_id]);
//                            assert(n_f_ids.size()==2);
//                            for (int n_f_id:n_f_ids)
//                                tmp_conn_fs.erase(n_f_id);

//                            cout<<r_v_id<<" "<< i + 1<<" "<< n_size<<endl;
//                            cout<<"f_id: ";
//                            for(int f:bsp_faces[new_f_id])
//                                cout<<f<<" ";
//                            cout<<endl;
//                            cout<<"new_f_id: ";
//                            for(int f:bsp_faces[f_id])
//                                cout<<f<<" ";
//                            cout<<endl;
                        }
                        v_id = p2_id;
                    } else {
//                        cout<<">>>cut the edge into 2"<<endl;
                        bsp_vertices.push_back(intersection_v);
                        int new_v_id = bsp_vertices.size() - 1;
                        conn_fs.resize(conn_fs.size() + 1);
                        tag_boundary_es.resize(tag_boundary_es.size() + 1);

                        std::vector<int> p1p2_e_ids = optimization::set_intersection(tag_boundary_es[p1_id], tag_boundary_es[p2_id]);
                        for(int p1p2_e_id:p1p2_e_ids)
                            tag_boundary_es[new_v_id].push_back(p1p2_e_id);

                        bsp_faces[f_id].insert(bsp_faces[f_id].begin() + j + 1, new_v_id);
                        std::vector<int> n_f_ids = optimization::set_intersection(conn_fs[p1_id], conn_fs[p2_id]);
                        int n_f_id = n_f_ids[0] == f_id ? n_f_ids[1] : n_f_ids[0];
                        for (int k = 0; k < bsp_faces[n_f_id].size(); k++) {
                            if (bsp_faces[n_f_id][k] == p2_id
                                && bsp_faces[n_f_id][(k + 1) % bsp_faces[n_f_id].size()] == p1_id) {
                                bsp_faces[n_f_id].insert(bsp_faces[n_f_id].begin() + k + 1, new_v_id);
                                break;
                            }
                        }

                        if (j + 1 <= r_v_id)
                            r_v_id++;
                        int new_f_id = split_a_face(r_v_id, i + 1, n_size + 1, f_id);

//                        cout<<"f_id: ";
//                        for(int f:bsp_faces[f_id])
//                            cout<<f<<" ";
//                        cout<<endl;
//                        cout<<"new_f_id: ";
//                        for(int f:bsp_faces[new_f_id])
//                            cout<<f<<" ";
//                        cout<<endl;

                        conn_fs[new_v_id].insert(n_f_id);
                        tag_boundary_es[new_v_id].push_back(e_id);

                        v_id = new_v_id;
                        tmp_conn_fs.clear();
                        tmp_conn_fs.insert(n_f_id);
                    }

//                    cout<<"v_id = "<<v_id<<endl<<endl;
                    is_intersected = true;
                    break;
                }
//                for(int i=0;i<bsp_faces.size();i++){
//                    if(bsp_faces[i].size()<3){
//                        cout<<"bsp_faces[i].size()<3, face "<<i<<endl;
//                        optimization::pausee();
//                    }
//                    for(int v_id:bsp_faces[i])
//                        if(conn_fs[v_id].find(i) == conn_fs[v_id].end()){
//                            cout<<"!conn_fs["<<v_id<<"].contains("<<i<<")"<<endl;
//                        }
//                }
//                optimization::pausee();

                if (is_intersected)
                    break;
            }
//            if (!is_intersected) {
//                cout << "cannot find intersection" << endl;
//                cout << "old v_id = " << v_id << endl;
//                for (int f_id:tmp_conn_fs) {
//                    cout << "face " << f_id << ": ";
//                    for (int f:bsp_faces[f_id])
//                        cout << f << " ";
//                    cout << endl;
//                }
//                cout<<"///"<<endl;
//                for (int f_id:conn_fs[v_id]) {
//                    cout << "face " << f_id << ": ";
//                    for (int f:bsp_faces[f_id])
//                        cout << f << " ";
//                    cout << endl;
//                }
//                cout<<"///"<<endl;
//                for (int f_id:conn_fs[v2_id]) {
//                    cout << "face " << f_id << ": ";
//                    for (int f:bsp_faces[f_id])
//                        cout << f << " ";
//                    cout << endl;
//                }
//                optimization::pausee();
//            }
            if (is_finished)
                break;
        }
    }

    ////construct triangle mesh
    for (auto &v: mesh.tri_vertices) {
        v.posf[0] = v.pos[0].to_double();
        v.posf[1] = v.pos[1].to_double();
    }

    mesh.tris.reserve(bsp_faces.size());
    for (int f_id = 0; f_id < bsp_faces.size(); f_id++) {
        assert(bsp_faces[f_id].size() >= 3);
        if (bsp_faces[f_id].size() == 3) {
            mesh.tris.push_back({{bsp_faces[f_id][0], bsp_faces[f_id][1], bsp_faces[f_id][2]}});
            continue;
        }
        for (int j = 1; j < int(bsp_faces[f_id].size()) - 1; j++) {
            mesh.tris.push_back({{bsp_faces[f_id][0], bsp_faces[f_id][j], bsp_faces[f_id][j + 1]}});
        }
    }

    //output and check
//    Eigen::MatrixXd oV(mesh.tri_vertices.size(), 3), _;
//    Eigen::MatrixXi oF(mesh.tris.size(), 3);
//    for (int i = 0; i < mesh.tri_vertices.size(); i++)
//        oV.row(i) << mesh.tri_vertices[i].posf[0], mesh.tri_vertices[i].posf[1], 0;
//    for (int i = 0; i < mesh.tris.size(); i++)
//        oF.row(i) << mesh.tris[i][0], mesh.tris[i][1], mesh.tris[i][2];
//    igl::writeSTL(args.output+"_bsp.stl", oV, oF, _);

    cout << "#v " << old_cnt_v << "->" << mesh.tri_vertices.size() << endl;
    cout << "#f " << old_cnt_f << "->" << mesh.tris.size() << endl;
}

bool triwild::triangulation::segment_intersection(MeshData& mesh, int v_id, int v2_id, int p1_id, int p2_id,
                                                    bool& is_cross_p1, bool& is_cross_p2, TriVertex& intersection_v, bool is_check_bbox) {
    auto &bsp_vertices = mesh.tri_vertices;
    is_cross_p1 = false;
    is_cross_p2 = false;

    ////check bbox
    if (is_check_bbox) {
        for (int j = 0; j < 2; j++) {
            if (bsp_vertices[v_id].pos[j] < bsp_vertices[v2_id].pos[j]) {
                if (bsp_vertices[p1_id].pos[j] < bsp_vertices[p2_id].pos[j]) {
                    if (bsp_vertices[p1_id].pos[j] > bsp_vertices[v2_id].pos[j])
                        return false;
                } else {
                    if (bsp_vertices[p2_id].pos[j] > bsp_vertices[v2_id].pos[j])
                        return false;
                }
            } else {
                if (bsp_vertices[p1_id].pos[j] < bsp_vertices[p2_id].pos[j]) {
                    if (bsp_vertices[p1_id].pos[j] > bsp_vertices[v_id].pos[j])
                        return false;
                } else {
                    if (bsp_vertices[p2_id].pos[j] > bsp_vertices[v_id].pos[j])
                        return false;
                }
            }
        }
    }

//    ////check side
//    Vector_2 n(bsp_vertices[v_id].pos[1] - bsp_vertices[v2_id].pos[1],
//               bsp_vertices[v2_id].pos[0] - bsp_vertices[v_id].pos[0]);
//    Vector_2 vp1 = bsp_vertices[p1_id].pos - bsp_vertices[v_id].pos;
//    Vector_2 vp2 = bsp_vertices[p2_id].pos - bsp_vertices[v_id].pos;
//    int sign1 = (n.dot(vp1)).get_sign();
//    int sign2 = (n.dot(vp2)).get_sign();
////    cout << sign1 << " " << sign2 << endl;
////    if (sign1 == 0) {
////        is_cross_p1 = true;
////        return true;
////    }
////    if (sign2 == 0) {
////        is_cross_p2 = true;
////        return true;
////    }
//    if (sign1 == sign2)
//        return false;

    ////compute
    auto &x1 = bsp_vertices[v_id].pos[0];
    auto &y1 = bsp_vertices[v_id].pos[1];
    auto &x2 = bsp_vertices[v2_id].pos[0];
    auto &y2 = bsp_vertices[v2_id].pos[1];
    auto &x3 = bsp_vertices[p1_id].pos[0];
    auto &y3 = bsp_vertices[p1_id].pos[1];
    auto &x4 = bsp_vertices[p2_id].pos[0];
    auto &y4 = bsp_vertices[p2_id].pos[1];

    Rational td = (x4 - x3) * (y1 - y2) - (x1 - x2) * (y4 - y3);
    int td_sign = td.get_sign();
    if (td_sign == 0)
        return false;
    Rational tn = (y3 - y4) * (x1 - x3) + (x4 - x3) * (y1 - y3);
    if (tn.get_sign() * td_sign < 0 || tn > td)
        return false;
    tn = (y1 - y2) * (x1 - x3) + (x2 - x1) * (y1 - y3);
    if (tn.get_sign() * td_sign < 0 || tn > td)
        return false;

    if (tn.get_sign() == 0) {
        is_cross_p1 = true;
        return true;
    } else if (tn == td) {
        is_cross_p2 = true;
        return true;
    }

    auto &t = tn;
    t = tn / td;

    intersection_v.pos[0] = x3 + t * (x4 - x3);
    intersection_v.pos[1] = y3 + t * (y4 - y3);

    intersection_v.pos[0].canonicalize();
    intersection_v.pos[1].canonicalize();

    return true;
}