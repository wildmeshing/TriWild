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

#include "do_triwild.h"

void triwild::do_triwild(const Eigen::MatrixXd &V_in, const Eigen::MatrixXi &E_in, json &feature_info,
                         Eigen::MatrixXd &V_out, Eigen::MatrixXi &F_out, Eigen::MatrixXd &nodes, std::vector<std::vector<int>> &F_nodes,
                         double stop_quality, int max_its, int stage,
                         double epsilon, double feature_epsilon,
                         double target_edge_len, double edge_length_r,
                         double flat_feature_angle,
                         bool cut_outside,
                         const Eigen::MatrixXd hole_pts,
                         bool mute_log)
{

    args.stop_quality = stop_quality;
    args.max_its = max_its;
    args.stage = stage;
    args.epsilon = epsilon;
    args.feature_epsilon = feature_epsilon;
    args.target_edge_len = target_edge_len;
    args.edge_length_r = edge_length_r;
    args.flat_feature_angle = flat_feature_angle;
    args.mute_log = mute_log;

    auto get_outputs = [&](const MeshData &mesh) {
        std::unordered_map<int, int> map_v_ids;
        int cnt = 0;
        for (size_t i = 0; i < mesh.tri_vertices.size(); i++)
        {
            if (mesh.v_is_removed[i])
                continue;
            map_v_ids[i] = cnt++;
        }
        V_out.resize(cnt, 2);
        cnt = 0;
        for (size_t i = 0; i < mesh.tri_vertices.size(); i++)
        {
            if (mesh.v_is_removed[i])
                continue;
            V_out(cnt, 0) = mesh.tri_vertices[i].posf[0];
            V_out(cnt, 1) = mesh.tri_vertices[i].posf[1];
            cnt++;
        }
        cnt = std::count(mesh.t_is_removed.begin(), mesh.t_is_removed.end(), false);
        F_out.resize(cnt, 3);
        cnt = 0;
        for (size_t i = 0; i < mesh.tris.size(); i++)
        {
            if (mesh.t_is_removed[i])
                continue;
            for (int j = 0; j < 3; j++)
                F_out(cnt, j) = map_v_ids[mesh.tris[i][j]];
            cnt++;
        }

        map_v_ids.clear();
        cnt = 0;
        for (size_t i = 0; i < mesh.nodes.size(); i++)
        {
            if (mesh.n_is_removed[i])
                continue;
            map_v_ids[i] = cnt++;
        }
        nodes.resize(cnt, 2);
        cnt = 0;
        for (size_t i = 0; i < mesh.nodes.size(); i++)
        {
            if (mesh.n_is_removed[i])
                continue;
            nodes(cnt, 0) = mesh.nodes[i][0];
            nodes(cnt, 1) = mesh.nodes[i][1];
            cnt++;
        }
        cnt = std::count(mesh.t_is_removed.begin(), mesh.t_is_removed.end(), false);
        F_nodes.clear();
        F_nodes.resize(cnt);
        cnt = 0;
        for (size_t i = 0; i < mesh.tri_nodes.size(); i++)
        {
            if (mesh.t_is_removed[i])
                continue;
            for (int j = 0; j < mesh.tri_nodes[i].size(); j++)
                F_nodes[cnt].push_back(map_v_ids[mesh.tri_nodes[i][j]]);
            cnt++;
        }
    };

    GEO::initialize();

    igl::Timer main_timer;
    double load_mesh_time = 0;
    double preprocess_time = 0;
    double bsp_time = 0;
    double optimization_time = 0;
    double curving_time = 0;
    double cut_and_hole_time = 0;

#ifndef WIN32
    if (args.mute_log)
    {
        std::streambuf *orig_buf = cout.rdbuf();
        cout.rdbuf(NULL);
    }
#endif

    igl::Timer igl_timer;
    double t;
    ///////////////////
    cout << "Loading and preprocessing ..." << endl;
    igl_timer.start();

    Eigen::MatrixXd V = V_in;
    std::vector<std::array<int, 2>> edges(E_in.rows());
    for (int i = 0; i < E_in.rows(); i++)
    {
        for (int j = 0; j < 2; j++)
            edges[i][j] = E_in(i, j);
    }

    //    Eigen::MatrixXd V;
    //    std::vector<std::array<int, 2>> edges;
    //    triangulation::load_input(args.input, V, edges);
    int input_v_cnt = V.rows();
    int input_e_cnt = edges.size();
    load_mesh_time = igl_timer.getElapsedTime();
    igl_timer.start();

    //    if (feature::init(args.feature_input))
    //        args.is_preserving_feature = true;
    if (feature::init(feature_info))
        args.is_preserving_feature = true;

    if (args.stop_quality < 0)
        args.stop_quality = args.is_preserving_feature ? 20 : 10;
    if (args.epsilon < 0)
        //        args.epsilon = args.is_preserving_feature ? 5e-3 : 1e-3;
        args.epsilon = args.is_preserving_feature ? 2e-3 : 1e-3;

    double line_width = 0.3 * args.diagonal_len;
    double f_line_width = 0.5 * args.diagonal_len;
    double s_f_line_width = 0.4 * args.diagonal_len;
    //        double draw_points = true;
    double draw_points = false;

    double point_size = 0.0005 * args.diagonal_len;
    double f_point_size = 0.001 * args.diagonal_len;
    double s_f_point_size = 0.003 * args.diagonal_len;

    std::string line_col = "0 0 0";
    std::string f_line_col = "0.9 0.3 0.2";
    std::string s_f_line_col = "0.1 0.7 0.6";

    std::string point_col = "0 0 0";
    std::string f_point_col = "0.9 0.3 0.2";
    std::string s_f_point_col = "0.1 0.7 0.6";

    GEO::Mesh b_mesh;
    triangulation::preprocessing(V, edges, b_mesh);
    GEO::MeshFacetsAABB b_tree(b_mesh);
    t = igl_timer.getElapsedTime();
    preprocess_time = t;
    cout << "Loaded and preprocessed." << endl;
    cout << "time = " << t << "s" << endl
         << endl;

    ///////////////////
    cout << "BSP subdivision..." << endl;
    igl_timer.start();
    MeshData mesh;
    std::vector<std::vector<int>> tag_boundary_es;
    triangulation::BSP_subdivision(V, edges, mesh, tag_boundary_es);
    t = igl_timer.getElapsedTime();
    bsp_time = t;
    cout << "BSP subdivision done." << endl;
    cout << "time = " << t << "s" << endl
         << endl;

    ///////////////////
    cout << "Mesh optimization..." << endl;
    igl_timer.start();
    optimization::init(V, edges, tag_boundary_es, mesh, b_tree);

    optimization::refine(mesh, b_tree, std::array<int, 4>({1, 1, 1, 1}));
    t = igl_timer.getElapsedTime();
    optimization_time = t;
    cout << "Mesh optimization done." << endl;
    cout << "time = " << t << "s" << endl;

    //    if (!skip_eps) {
    //         export_eps(mesh,
    //             line_width, line_col, point_size, point_col,
    //             f_line_width, f_line_col, f_point_size, f_point_col,
    //             s_f_line_width, s_f_line_col, s_f_point_size, s_f_point_col,
    //             draw_points, args.output + "_lin.eps");
    //     }

    ///////////////////
    if (!args.is_preserving_feature)
    {
        //        optimization::output_mesh(mesh);//todo
        get_outputs(mesh);
        return;
    }

    cout << "Curving..." << endl;
    igl_timer.start();
    feature::curving(mesh, b_tree);
    t = igl_timer.getElapsedTime();
    curving_time = t;
    cout << "Curving done." << endl;
    cout << "time = " << t << "s" << endl;

    igl_timer.start();
    if (cut_outside)
        optimization::erase_outside(mesh);
    if (hole_pts.size() > 0)
        optimization::erase_holes(mesh, hole_pts);
    cut_and_hole_time = igl_timer.getElapsedTime();

    //    optimization::output_mesh(mesh);//todo
    get_outputs(mesh);
}