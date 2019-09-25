// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef TRIWILD_OPTIMIZATION_H
#define TRIWILD_OPTIMIZATION_H

#include "Args.h"
#include "TrimeshElements.h"

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>

namespace triwild {
    namespace optimization{
        void init(const Eigen::MatrixXd& V, const std::vector<std::array<int, 2>>& edges,
                const std::vector<std::vector<int>>& tag_boundary_es,
                MeshData& mesh, GEO::MeshFacetsAABB &b_tree);
        void refine(MeshData& mesh, GEO::MeshFacetsAABB &b_tree, const std::array<int, 4>& ops= {{1,1,1,1}});
        void operation(MeshData& mesh, GEO::MeshFacetsAABB &b_tree, const std::array<int, 4>& ops= {{1,1,1,1}});
        bool update_scaling_field(MeshData& mesh, double max_energy);
        void output_mesh(MeshData& mesh);

        void erase_outside(MeshData& mesh);
        void erase_holes(MeshData& mesh, const std::string& hole_file, bool is_erase = true);
        void erase_holes(MeshData& mesh, const Eigen::MatrixXd& hole_pts, bool is_erase = true);
        void erase_holes(MeshData &mesh, const std::vector<GEO::vec3> &ps, bool is_erase = true);

        bool round_a_vertex(MeshData& mesh, int v_id);

        bool is_valid_quality(Point_2f& p1, Point_2f& p2, Point_2f& p3, double old_max_energy,
                double& new_energy, bool accept_equal = false);
        bool is_valid_inversion(const Point_2f& pf1, const Point_2f& pf2, const Point_2f& pf3,
                const Point_2& p1, const Point_2& p2, const Point_2& p3);
        bool is_valid_envelop(MeshData& mesh, GEO::MeshFacetsAABB &b_tree, int v_id, Point_2f& p);
        bool is_valid_feature_edge_length(MeshData& mesh, int v_id, Point_2f& p);
        bool is_valid_feature_edge_close(MeshData& mesh, int v_id, const Point_2f& p, double t);

        bool is_bbox_edge(const MeshData& mesh, int v1_id, int v2_id);
        bool is_boundary_edge(MeshData& mesh, int v1_id, int v2_id);
        bool is_feature_edge(MeshData& mesh, int v1_id, int v2_id);
        bool is_valid_edge(MeshData& mesh, int v1_id, int v2_id);

        int get_feature_edge_tag(const MeshData& mesh, int v1_id, int v2_id);
        int get_secondary_feature_edge_tag(const MeshData& mesh, int v1_id, int v2_id);
        bool is_edge_out_of_envelop(MeshData& mesh, GEO::MeshFacetsAABB &b_tree, Point_2f& p1, Point_2f& p2);

        bool is_isolated_vertex(MeshData& mesh, int v_id);

        void get_all_edges(const MeshData& mesh, std::vector<std::array<int, 2>>& edges);
        double edge_length_2(const MeshData& mesh, int v1_id, int v2_id);
        double edge_length_2(const Point_2f& p1, const Point_2f& p2);
        Point_2 mid_point(const Point_2& p1, const Point_2& p2);
        Point_2f mid_point(const Point_2f& p1, const Point_2f& p2);
        double tri_energy(Point_2f& p1, Point_2f& p2, Point_2f& p3);
        void get_max_avg_energy(const MeshData& mesh, double& max_energy, double& avg_energy);
        double get_mid_energy(const MeshData& mesh);

        bool get_opp_r_v_id(const MeshData& mesh, int t_id, int v1_id, int v2_id);

        //utilities
        void output_stats(const MeshData& mesh, std::ofstream& f);
        double t_area(const MeshData& mesh, int t_id);
        void check(MeshData& mesh);
        double check_envelope(MeshData& mesh, GEO::MeshFacetsAABB &b_tree, double s = 2);
        void output_info(MeshData& mesh);
        void output_boundary(MeshData& mesh, bool is_feature = false);
        template<typename T>
        std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec);
        template<typename T>
        void vector_unique(std::vector<T>& v){
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
        };
        std::vector<int> set_intersection(const std::vector<int>& v1, const std::vector<int>& v2);
        std::vector<int> set_intersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2);
        void pausee();

        class ElementInQueue_s{
        public:
            std::array<int, 2> v_ids;
            double weight;

            ElementInQueue_s(){}
            ElementInQueue_s(const std::array<int, 2>& ids, double w): v_ids(ids), weight(w){}
        };
        struct cmp_s {
            bool operator()(const ElementInQueue_s &e1, const ElementInQueue_s &e2) {
                if (e1.weight == e2.weight)
                    return e1.v_ids > e2.v_ids;
                return e1.weight > e2.weight;///choose larger edge for removal
            }
        };
        class ElementInQueue_l{
        public:
            std::array<int, 2> v_ids;
            double weight;

            ElementInQueue_l(){}
            ElementInQueue_l(const std::array<int, 2>& ids, double w): v_ids(ids), weight(w){}
        };
        struct cmp_l {
            bool operator()(const ElementInQueue_l &e1, const ElementInQueue_l &e2) {
                if (e1.weight == e2.weight)
                    return e1.v_ids < e2.v_ids;
                return e1.weight < e2.weight;///choose larger edge for removal
            }
        };
    }
}


#endif //TRIWILD_OPTIMIZATION_H
