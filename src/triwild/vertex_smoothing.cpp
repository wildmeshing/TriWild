// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "vertex_smoothing.h"
#include "optimization.h"
#include "feature.h"
#include "AMIPS.h"

#include "meshio.hpp"

void triwild::optimization::vertex_smoothing(MeshData& mesh, GEO::MeshFacetsAABB &b_tree) {
    int cnt = 0;
    int cnt_b = 0;

    for (size_t v_id = 0; v_id < mesh.tri_vertices.size(); v_id++) {
        if (mesh.v_is_removed[v_id])
            continue;

        if (!mesh.tri_vertices[v_id].is_rounded && !round_a_vertex(mesh, v_id))
            continue;

        if (mesh.tri_vertices[v_id].is_freezed)
            continue;
        if (mesh.tri_vertices[v_id].is_on_bbox)
            continue;
        if (mesh.tri_vertices[v_id].is_on_point)
            continue;

//        cout<<v_id<<endl;

        Point_2f pf;
        double t = -1;
        if (!smooth_a_vertex(mesh, v_id, pf, t))
            continue;
        Point_2 p(pf[0], pf[1]);

        std::vector<double> new_qs;
        if (mesh.tri_vertices[v_id].is_on_boundary && mesh.tri_vertices[v_id].feature_infos.size() == 0) {
            GEO::vec3 geo_pf(pf[0], pf[1], 0);
            GEO::vec3 nearest_pf;
            double _;
            b_tree.nearest_facet(geo_pf, nearest_pf, _);
            pf = Point_2f(nearest_pf[0], nearest_pf[1]);
            p = Point_2(pf[0], pf[1]);
        }

        ////check inversion, quality, envelop
        bool is_valid = true;
        double old_max_energy = 0;
        for (int t_id:mesh.tri_vertices[v_id].conn_tris) {
            if (mesh.t_quality[t_id] > old_max_energy)
                old_max_energy = mesh.t_quality[t_id];
        }
        std::vector<std::vector<Point_2f>> vec_new_nodes;
        for (int t_id:mesh.tri_vertices[v_id].conn_tris) {
            int j = std::find(mesh.tris[t_id].begin(), mesh.tris[t_id].end(), v_id) - mesh.tris[t_id].begin();
            if (!mesh.is_curved || (mesh.is_curved && mesh.tri_nodes[t_id].empty())) {
                is_valid = is_valid_inversion(pf, mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].posf,
                                              mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].posf,
                                              p, mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].pos,
                                              mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].pos);
            } else {
                Point_2f old_pf = mesh.tri_vertices[mesh.tris[t_id][j]].posf;
                Point_2 old_p = mesh.tri_vertices[mesh.tris[t_id][j]].pos;
                mesh.tri_vertices[mesh.tris[t_id][j]].posf = pf;
                mesh.tri_vertices[mesh.tris[t_id][j]].pos = p;
                double old_t;
                if (t != -1) {
                    old_t = mesh.tri_vertices[mesh.tris[t_id][j]].feature_infos[0][1];
                    mesh.tri_vertices[mesh.tris[t_id][j]].feature_infos[0][1] = t;
                }

                vec_new_nodes.emplace_back();
                auto &new_nodes = vec_new_nodes.back();
                feature::get_new_nodes(mesh.tag_feature_es[t_id],
                                       {{mesh.tri_vertices[mesh.tris[t_id][0]], mesh.tri_vertices[mesh.tris[t_id][1]],
                                                mesh.tri_vertices[mesh.tris[t_id][2]]}}, new_nodes);
                assert(!new_nodes.empty());
                is_valid = feature::is_valid_inversion(
                        {{mesh.tri_vertices[mesh.tris[t_id][0]].posf, mesh.tri_vertices[mesh.tris[t_id][1]].posf,
                                 mesh.tri_vertices[mesh.tris[t_id][2]].posf}}, new_nodes);

                mesh.tri_vertices[mesh.tris[t_id][j]].posf = old_pf;
                mesh.tri_vertices[mesh.tris[t_id][j]].pos = old_p;
                if (t != -1)
                    mesh.tri_vertices[mesh.tris[t_id][j]].feature_infos[0][1] = old_t;
            }
            if (!is_valid)
                break;

            double new_energy;
            if(!is_valid_quality(pf, mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].posf,
                    mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].posf,
                    old_max_energy, new_energy)){
                is_valid = false;
                break;
            }

//            if (mesh.tri_vertices[v_id].is_on_boundary) {
//                if (mesh.tri_vertices[v_id].feature_infos.size() == 0)
//                    is_valid = is_valid && is_valid_envelop(mesh, b_tree, v_id, pf);
//                else {
//                    if (!mesh.is_curved)
//                        is_valid = is_valid && is_valid_feature_edge_close(mesh, v_id, pf, t);
//                }
//            }
//            if (!is_valid)
//                break;

            if (mesh.tri_vertices[v_id].is_on_boundary && !is_valid_envelop(mesh, b_tree, v_id, pf)){
                is_valid = false;
                break;
            }

            if (!mesh.is_curved && !mesh.tri_vertices[v_id].feature_infos.empty()
                && !is_valid_feature_edge_close(mesh, v_id, pf, t)) {
                is_valid = false;
                break;
            }

            new_qs.push_back(new_energy);
        }
        if (!is_valid)
            continue;

        ////real update
        if (!vec_new_nodes.empty()) {
            int cnt = 0;
            for (int t_id:mesh.tri_vertices[v_id].conn_tris) {
                if (!mesh.tri_nodes[t_id].empty()) {
                    auto &nodes = vec_new_nodes[cnt++];
                    assert(nodes.size() == mesh.tri_nodes[t_id].size());
                    for (size_t i = 0; i < mesh.tri_nodes[t_id].size(); i++)
                        mesh.nodes[mesh.tri_nodes[t_id][i]] = nodes[i];
                }
            }
        }

        mesh.tri_vertices[v_id].posf = pf;
        mesh.tri_vertices[v_id].pos = p;
        if (mesh.tri_vertices[v_id].feature_infos.size() > 0)
            mesh.tri_vertices[v_id].feature_infos[0][1] = t;

        if (new_qs.size() == 0) {
            for (int t_id:mesh.tri_vertices[v_id].conn_tris)
                mesh.t_quality[t_id] = tri_energy(mesh.tri_vertices[mesh.tris[t_id][0]].posf,
                                                  mesh.tri_vertices[mesh.tris[t_id][1]].posf,
                                                  mesh.tri_vertices[mesh.tris[t_id][2]].posf);
        } else {
            int i = 0;
            for (int t_id:mesh.tri_vertices[v_id].conn_tris)
                mesh.t_quality[t_id] = new_qs[i++];
        }

        if (mesh.tri_vertices[v_id].is_on_boundary)
            cnt_b++;
        cnt++;

//        if (mesh.is_curved && mesh.tri_vertices[v_id].feature_infos.size() > 0) {
//            write_msh(mesh, args.output + std::to_string(v_id) + ".msh");
//            pausee();
//        }
    }
    cout << "success all " << cnt << endl;
    cout << "success boundary " << cnt_b << endl;
}

bool triwild::optimization::smooth_a_vertex(MeshData& mesh, const int v_id, Point_2f& pf, double& t) {
    int feature_id = -1;
    if (mesh.tri_vertices[v_id].feature_infos.size() > 0)
        feature_id = mesh.tri_vertices[v_id].feature_infos[0][0];

    std::vector<std::array<double, 6>> Ts;
    for (int t_id:mesh.tri_vertices[v_id].conn_tris) {
        int j = std::find(mesh.tris[t_id].begin(), mesh.tris[t_id].end(), v_id) - mesh.tris[t_id].begin();
        Ts.push_back(std::array<double, 6>({mesh.tri_vertices[v_id].posf[0], mesh.tri_vertices[v_id].posf[1],
                                            mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].posf[0],
                                            mesh.tri_vertices[mesh.tris[t_id][(j + 1) % 3]].posf[1],
                                            mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].posf[0],
                                            mesh.tri_vertices[mesh.tris[t_id][(j + 2) % 3]].posf[1]}));
    }

    ////newton
    const int max_newton_it = 15;
    const int max_search_it = 10;
    const double f_delta = 1e-8;
    const double J_delta = 1e-8;

    Eigen::Vector2d x(mesh.tri_vertices[v_id].posf[0], mesh.tri_vertices[v_id].posf[1]);
    Eigen::Vector2d J;
    Eigen::Matrix2d H;
    double J_feature, H_feature;
    if (feature_id >= 0)
        t = mesh.tri_vertices[v_id].feature_infos[0][1];

//    double f_old, f_new;
//    int it;
    Point_2f pf_next, pf1, pf2;
    Point_2 p_next, p1, p2;
    for (int newton_it = 0; newton_it < max_newton_it; newton_it++) {
        if (newton_it > 0 && feature_id < 0) {
            for (auto &T:Ts) {
                T[0] = x(0);
                T[1] = x(1);
            }
        }

        //f
        double f = 0;
        for (auto &T:Ts) {
            if (feature_id >= 0)
                f += AMIPS_energy(*(feature::features[feature_id]), t, T);
            else
                f += AMIPS_energy(T);
        }

//        if(newton_it == 0)
//            f_old = f;
//        f_new = f;
//        it = newton_it;

        //J
        J << 0, 0;
        J_feature = 0;
        for (auto &T:Ts) {
            if (feature_id >= 0) {
                J_feature += AMIPS_jacobian(*(feature::features[feature_id]), t, T);
            } else {
                Eigen::Vector2d tmp_J;
                AMIPS_jacobian(T, tmp_J);
                J += tmp_J;
            }
        }
        if (feature_id >= 0) {
            if(std::isnan(J_feature) || std::isinf(J_feature))
                break;
        } else {
            if (!J.allFinite() || (std::abs(J(0)) < J_delta && std::abs(J(1)) < J_delta))//gradient is also zero
                break;
        }

        //H
        H << 0, 0, 0, 0;
        H_feature = 0;
        for (auto &T:Ts) {
            if (feature_id >= 0) {
                H_feature += AMIPS_hessian(*(feature::features[feature_id]), t, T);
            } else {
                Eigen::Matrix2d tmp_H;
                AMIPS_hessian(T, tmp_H);
                H += tmp_H;
            }
        }
        if(feature_id >= 0) {
            if(std::isnan(H_feature) || std::isinf(H_feature) || H_feature == 0)
                break;
        } else {
            if (!H.allFinite())
                break;
        }

        //x
        bool found_step = false;
        double a = 1;
        Eigen::Vector2d x_next;
        double t_next;
        for (int i = 0; i < max_search_it; i++) {
            if (feature_id >= 0) {
                t_next = t - a * J_feature / H_feature;
                if(std::isnan(t_next) || std::isinf(t_next) || t_next > 1)
                    break;
                if(t_next < 0){
                    a /= 2;
                    continue;
                }
            } else {
                x_next = H.colPivHouseholderQr().solve(H * x - a * J);
                //https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
                //JacobiSVD
                if(!x_next.allFinite())
                    break;
                for (auto &T:Ts) {
                    T[0] = x_next(0);
                    T[1] = x_next(1);
                }
            }


            //check energy
            double f_next = 0;
            for (auto &T:Ts) {
                if (feature_id >= 0)
                    f_next += AMIPS_energy(*(feature::features[feature_id]), t_next, T);
                else
                    f_next += AMIPS_energy(T);
            }
            if (f_next >= f) {
                a /= 2;
                continue;
            }
            //check inversion
            bool is_valid = true;
//            Point_2f pf_next(x_next(0), x_next(1));
//            Point_2 p_next(x_next(0), x_next(1));
            if(feature_id >= 0){
                pf_next = feature::features[feature_id]->eval(t_next);
                p_next[0] = pf_next[0];
                p_next[1] = pf_next[1];
            } else {
                pf_next[0] = x_next(0);
                pf_next[1] = x_next(1);
                p_next[0] = x_next(0);
                p_next[1] = x_next(1);
            }

            for (auto &T:Ts) {
                pf1[0] = T[2];
                pf1[1] = T[3];
                pf2[0] = T[4];
                pf2[1] = T[5];
                p1[0] = T[2];
                p1[1] = T[3];
                p2[0] = T[4];
                p2[1] = T[5];
//                if (!is_valid_inversion(pf_next, Point_2f(T[2], T[3]), Point_2f(T[4], T[5]),
//                                        p_next, Point_2(T[2], T[3]), Point_2(T[4], T[5]))) {
                if (!is_valid_inversion(pf_next, pf1, pf2, p_next, p1, p2)) {
                    is_valid = false;
                    break;
                }
            }
            if (!is_valid) {
                a /= 2;
                continue;
            }

            found_step = true;
            break;
        }

        if (feature_id >= 0) {
            if (!found_step || std::isnan(t_next) || std::isinf(t_next) || t>=1 || t<=0)
                break;
            t = t_next;
        } else {
            if (!found_step || !x_next.allFinite())
                break;
            x = x_next;
        }
    }

    if (feature_id >= 0) {
        if (t != mesh.tri_vertices[v_id].feature_infos[0][1]) {
            pf = feature::features[feature_id]->eval(t);
            return true;
        }
    } else {
        if (!(x(0) == mesh.tri_vertices[v_id].posf[0] && x(1) == mesh.tri_vertices[v_id].posf[1])) {
            pf[0] = x(0);
            pf[1] = x(1);
//            cout<<" "<<f_old<<"->"<<f_new<<"("<<it<<") ";
            return true;
        }
    }
    return false;
}
