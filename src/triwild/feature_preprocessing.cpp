// This file is part of TriWild, a software for generating linear/curved
// triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "feature_preprocessing.h"
#include "../../extern/aabbcc/src/AABB.h"
#include "Curves.h"
#include "feature.h"
#include "optimization.h"
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>
#include <igl/Timer.h>
#include <igl/unique_rows.h>
#include <igl/writeSTL.h>

void triwild::feature::preprocessing(Eigen::MatrixXd &V,
                                     std::vector<std::array<int, 2>> &edges) {

  ////throw away input segments
  //    Eigen::MatrixXd V_new(features.size() * 2, 2);
  //    std::vector<std::array<int, 2>> edges_new(features.size());
  //    for(int i=0;i<features.size();i++) {
  //        V_new.row(i * 2) = V.row(features[i]->v_ids.front());
  //        V_new.row(i * 2 + 1) = V.row(features[i]->v_ids.back());
  //        features[i]->v_ids = {i * 2, i * 2 + 1};
  //        features[i]->paras = {0, 1};
  //        edges_new[i] = {{i * 2, i * 2 + 1}};
  //    }
  //    //remove deplicates
  //    Eigen::VectorXi VI, _;
  //    igl::unique_rows(V_new, V, _, VI);
  //    for(auto& e:edges){
  //        e = {VI(e[0]), VI(e[1])};
  //    }

  gen_segments(V, edges);

  for (int i = 0; i < features.size(); i++) {
    // only for check
    {
      Point_2f p = features[i]->eval(0);
      Point_2f p0(V(features[i]->v_ids.front(), 0),
                  V(features[i]->v_ids.front(), 1));
      if (p != p0) {
        cout << "feature " << i << " ";
        cout << features[i]->type << endl;
        cout << std::setprecision(16) << p << endl << p0 << endl;
      }
    }
    {
      Point_2f p = features[i]->eval(1);
      Point_2f p0(V(features[i]->v_ids.back(), 0),
                  V(features[i]->v_ids.back(), 1));
      if (p != p0) {
        cout << "feature " << i << " ";
        cout << features[i]->type << endl;
        cout << std::setprecision(16) << p << endl << p0 << endl;
      }
    }

    features[i]->v_ids.clear();
    features[i]->paras = {0, 1};
  }

  ////simplification
  cout << "simplify" << endl;
  simplify();
  cout << "#features = " << features.size() << endl << endl;

  ////cut
  cout << "get_inflections & cut_reflex" << endl;
  std::vector<std::vector<double>> inflections;
  get_inflections(inflections);
  for (int i = 0; i < 2; i++)
    cut_reflex(inflections);
  cout << "#features = " << features.size() << endl << endl;

  ////samples
  cout << "sample_features" << endl << endl;
  igl::Timer igl_timer;
  std::vector<std::vector<Point_2f>> samples;
  std::vector<std::vector<double>> ts;
  sample_features(samples, ts);
  if (args.enable_debug_mesh)
    output_features("init");

  //    check(inflections, samples, ts);

  ////remove high curvature segments
  cout << "remove_high_curvature" << endl;
  remove_high_curvature(inflections, samples, ts);
  cout << "#features = " << features.size() << endl << endl;
  if (args.enable_debug_mesh)
    output_features("remove_high_curvature");

  //    check(inflections, samples, ts);

  ////remove short edges
  cout << "remove_short_features" << endl;
  remove_short_features(inflections, samples, ts);
  cout << "#features = " << features.size() << endl << endl;
  if (args.enable_debug_mesh)
    output_features("remove_short_features_1");

  //    check(inflections, samples, ts);

  ////mu separation
  cout << "mu_separation" << endl;
  igl_timer.start();
  mu_separation(inflections, samples, ts);
  cout << igl_timer.getElapsedTime() << "s" << endl;
  cout << "#features = " << features.size() << endl << endl;
  if (args.enable_debug_mesh)
    output_features("mu_separation");

  //    check(inflections, samples, ts);

  ////cut inflections
  cout << "cut_inflections" << endl;
  cut_inflections(inflections, samples, ts);
  cout << "#features = " << features.size() << endl << endl;

  //    check(inflections, samples, ts);

  ////remove short edges again
  cout << "remove_short_features again" << endl;
  remove_short_features(inflections, samples, ts);
  cout << "#features = " << features.size() << endl << endl;
  if (args.enable_debug_mesh)
    output_features("remove_short_features_2");

  ////construct input segments
  cout << "gen_segments" << endl;
  gen_segments(V, edges);
  cout << "#v = " << V.rows() << endl;
  cout << "#e = " << edges.size() << endl << endl;

  //    check(inflections, samples, ts);

  if (args.enable_debug_mesh) {
    output_input_features(features, V, edges, "final");
    //        output_input_features(secondary_features, V, edges,
    //        "final_secondary");
    Eigen::MatrixXd V_out(V.rows(), 3);
    for (int i = 0; i < V.rows(); i++)
      V_out.row(i) << V(i, 0), V(i, 1), 0;
    Eigen::MatrixXi E_out(edges.size(), 3);
    for (int i = 0; i < edges.size(); i++) {
      E_out.row(i) << edges[i][0], edges[i][1], edges[i][1];
    }
    igl::writeSTL(args.output + "_final_all.stl", V_out, E_out);
  }
  //    optimization::pausee();
}

void triwild::feature::simplify() {
  for (int i = 0; i < features.size(); i++) {
    if (features[i]->type != "BezierCurve")
      continue;

    if (features[i]->is_point()) {
      features.erase(features.begin() + i);
      i--;
    }
  }

  int simplified_cnt = 0;
  for (auto &feature : features) {
    if (feature->type != "BezierCurve")
      continue;

    Point_2f start = feature->eval(feature->paras.front());
    Point_2f end =
        feature->eval(feature->paras.back()); // TODO: would have gaps!!!
    auto tmp_feature = feature->simplify(start, end);
    if (tmp_feature == NULL)
      continue;
    feature = tmp_feature;
    simplified_cnt++;
  }
  cout << simplified_cnt << "/" << features.size() << " features are simplified"
       << endl;
}

void triwild::feature::get_inflections(
    std::vector<std::vector<double>> &inflections) {
  inflections.clear();
  inflections.resize(features.size());

  for (int i = 0; i < features.size(); i++) {
    inflections[i] = features[i]->inflection_points(features[i]->paras.front(),
                                                    features[i]->paras.back());
    if (inflections[i].empty())
      continue;
    std::sort(inflections[i].begin(), inflections[i].end());
    if (std::abs(inflections[i].front() - 0) < 1e-8)
      inflections[i].erase(inflections[i].begin());
    if (std::abs(inflections[i].back() - 1) < 1e-8)
      inflections[i].erase(inflections[i].begin() + inflections[i].size() - 1);
  }
}

void triwild::feature::cut_reflex(
    std::vector<std::vector<double>> &inflections) { // cut loop
  const int features_size = features.size();
  for (int feature_id = 0; feature_id < features_size; feature_id++) {
    if (features[feature_id]->type == "Line")
      continue;

    std::vector<double> ts;
    //        std::vector<bool> is_inflection;//not used

    std::vector<double> ranges = inflections[feature_id];
    ranges.insert(ranges.begin(), 0);
    ranges.push_back(1);
    for (int i = 0; i < ranges.size() - 1; i++) {
      double turning_angle =
          features[feature_id]->how_curve(ranges[i], ranges[i + 1]);
      if (turning_angle < 179)
        continue;
      double t = features[feature_id]->get_cut_t(ranges[i], ranges[i + 1],
                                                 (turning_angle >= 180));
      Point_2f p0 = features[feature_id]->eval(ranges[i]);
      Point_2f p1 = features[feature_id]->eval(ranges[i + 1]);
      Point_2f p = features[feature_id]->eval(t);
      if ((p0 - p).length_2() <= 1e-14) {
        ts.push_back(ranges[i]);
        //                if (i != 0)
        //                    is_inflection.push_back(true);
      } else if ((p1 - p).length_2() <= 1e-14) {
        ts.push_back(ranges[i + 1]);
        //                if (i + 1 != ranges.size() - 1)
        //                    is_inflection.push_back(true);
      } else {
        ts.push_back(t);
        //                is_inflection.push_back(false);
      }
    }

    // chop
    if (ts.empty())
      continue;
    ts.insert(ts.begin(), 0);
    ts.push_back(1);
    for (int i = 0; i < ts.size() - 1; i++) {
      if (ts[i + 1] == ts[i]) {
        ts.erase(ts.begin() + i + 1);
        i--;
        continue;
      }
      if (i == ts.size() - 2) {
        features[feature_id]->paras = {ts[i], ts[i + 1]};
        std::vector<double> new_inflections;
        for (double t : inflections[feature_id])
          if (t < ts[i + 1] && t > ts[i])
            new_inflections.push_back(t);
        inflections[feature_id] = new_inflections;
      } else {
        int new_feature_id = push_back_new_feature(features[feature_id]);
        features[new_feature_id]->paras = {ts[i], ts[i + 1]};
        inflections.emplace_back();
        for (double t : inflections[feature_id])
          if (t < ts[i + 1] && t > ts[i])
            inflections.back().push_back(t);
      }
    }
  }
}

void triwild::feature::sample_features(
    std::vector<std::vector<Point_2f>> &samples,
    std::vector<std::vector<double>> &ts) {
  /// sampling
  double dd = args.min_edge_length * args.target_edge_len;
  double dd_2 = dd * dd;

  samples.resize(features.size());
  ts.resize(features.size());
  for (int i = 0; i < features.size(); i++) {
    //        if (features[i]->type == "Line") {
    //            samples[i] = {features[i]->eval(features[i]->paras.front()),
    //            features[i]->eval(features[i]->paras.back())}; continue;
    //        }

    auto &feature_samples = samples[i];
    feature_samples = {features[i]->eval(features[i]->paras.front()),
                       features[i]->eval(features[i]->paras.back())};
    auto &feature_ts = ts[i];
    feature_ts = {features[i]->paras.front(), features[i]->paras.back()};
    std::vector<bool> is_short_enough = {false};
    while (true) {
      // check length
      int cnt = 0;
      for (int j = 0; j < is_short_enough.size(); j++) {
        if (is_short_enough[j])
          continue;
        if ((feature_samples[j] - feature_samples[j + 1]).length_2() <= dd_2)
          is_short_enough[j] = true;
        else
          cnt++;
      }
      if (cnt == 0)
        break;

      // insert samples
      for (int j = 0; j < feature_samples.size() - 1; j++) {
        if (is_short_enough[j])
          continue;
        double t = (feature_ts[j] + feature_ts[j + 1]) / 2;
        feature_ts.insert(feature_ts.begin() + j + 1, t);
        feature_samples.insert(feature_samples.begin() + j + 1,
                               features[i]->eval(t));
        is_short_enough.insert(is_short_enough.begin() + j + 1, false);
        j++;
      }
    }
  }
}

void triwild::feature::remove_short_features(
    std::vector<std::vector<double>> &inflections,
    std::vector<std::vector<Point_2f>> &samples,
    std::vector<std::vector<double>> &ts) {
  double min_edge_length = args.min_edge_length * args.target_edge_len;
  int cnt = 0;
  for (int i = 0; i < features.size(); i++) {
    if (samples[i].size() <= 4) {
      if (samples[i].size() > 2) {
        double l = 0;
        for (int j = 0; j < samples[i].size() - 1; j++) {
          l += (samples[i][j] - samples[i][j + 1]).length();
        }
        if (l >= min_edge_length)
          continue;
      }

      cnt++;
      push_back_new_secondary_feature(features[i]);

      features.erase(features.begin() + i);
      inflections.erase(inflections.begin() + i);
      samples.erase(samples.begin() + i);
      ts.erase(ts.begin() + i);
      i--;
    }
  }
  cout << cnt << " short edges are removed" << endl;
}

void triwild::feature::mu_separation(
    std::vector<std::vector<double>> &inflections,
    std::vector<std::vector<Point_2f>> &samples,
    std::vector<std::vector<double>> &ts) {
  double cos_angle_threshold = std::cos(20 / 180.0 * M_PI);
  double feature_eps = args.feature_epsilon * args.diagonal_len;
  double feature_eps_2 = feature_eps * feature_eps;

  /// get pairs
  aabb::Tree aabbcc_tree(2, 0, features.size(), true);
  std::vector<double> min_corner(2), max_corner(2);
  //    std::vector<std::array<Point_2f, 2>> boxes;
  //    boxes.resize(features.size());
  for (int i = 0; i < features.size(); i++) {
    Point_2f min, max;
    for (int j = 0; j < samples[i].size(); j++) {
      if (j == 0 || samples[i][j].x < min.x)
        min.x = samples[i][j].x;
      if (j == 0 || samples[i][j].x > max.x)
        max.x = samples[i][j].x;
      if (j == 0 || samples[i][j].y < min.y)
        min.y = samples[i][j].y;
      if (j == 0 || samples[i][j].y > max.y)
        max.y = samples[i][j].y;
    }
    min.x -= feature_eps;
    min.y -= feature_eps;
    max.x += feature_eps;
    max.y += feature_eps;

    //        boxes[i] = {min, max};

    min_corner[0] = min.x;
    min_corner[1] = min.y;
    max_corner[0] = max.x;
    max_corner[1] = max.y;
    aabbcc_tree.insertParticle(i, min_corner, max_corner);
  }
  std::vector<std::array<int, 2>> feature_pairs;
  for (int i = 0; i < features.size(); i++) {
    auto boxes = aabbcc_tree.query(i);
    for (auto &id : boxes) {
      //            if (features[i]->degree > features[id]->degree
      //                || features[i]->degree == features[id]->degree && i >
      //                id)
      if (samples[i].size() > samples[id].size() ||
          (samples[i].size() == samples[id].size() && i > id))
        feature_pairs.push_back({{i, (int)id}});
    }
  }
  //     for (int i = 0; i < features.size() - 1; i++) {
  //         for (int j = i + 1; j < features.size(); j++) {
  //             std::array<int, 2> pair;
  //             if (samples[i].size() > samples[j].size())
  //                 pair = {i, j};
  //             else
  //                 pair = {j, i};
  //
  //             if(i ==0 && j==2){
  //                 cout<<boxes[0][0]<<"; "<<boxes[0][1] <<endl;
  //                 cout<<boxes[2][0]<<"; "<<boxes[2][1] <<endl;
  //                 cout<<samples[0].front()<<endl;
  //                 cout<<samples[2].back()<<endl;
  //             }
  //             if (boxes[pair[0]][1].x < boxes[pair[1]][0].x ||
  //             boxes[pair[0]][0].x > boxes[pair[1]][1].x
  //                 || boxes[pair[0]][1].y < boxes[pair[1]][0].y ||
  //                 boxes[pair[0]][0].y > boxes[pair[1]][1].y) continue;
  //
  //             feature_pairs.push_back(pair);
  //         }
  //     }
  cout << "#pairs = " << feature_pairs.size() << endl;

  /// build aabb trees
  std::vector<bool> is_used(features.size(), false);
  for (auto &p : feature_pairs) {
    is_used[p[0]] = true;
    is_used[p[1]] = true;
  }
  cout << "#check_feature = "
       << std::count(is_used.begin(), is_used.end(), true) << "("
       << features.size() << ")" << endl;
  std::vector<std::shared_ptr<GEO::Mesh>> f_meshes;
  std::vector<std::shared_ptr<GEO::MeshFacetsAABB>> f_trees;
  for (int feature_id = 0; feature_id < features.size(); feature_id++) {
    f_meshes.push_back(std::make_shared<GEO::Mesh>());
    GEO::Mesh &f_mesh = *(f_meshes.back());
    f_mesh.vertices.clear();
    f_mesh.facets.clear();
    if (!is_used[feature_id]) {
      f_mesh.vertices.create_vertices(1);
      f_mesh.facets.create_triangles(1);
      f_mesh.facets.set_vertex(0, 0, 0);
      f_mesh.facets.set_vertex(0, 1, 0);
      f_mesh.facets.set_vertex(0, 2, 0);
      f_trees.push_back(std::make_shared<GEO::MeshFacetsAABB>(f_mesh, false));
      continue;
    }

    f_mesh.vertices.create_vertices(samples[feature_id].size());
    for (int i = 0; i < samples[feature_id].size(); i++) {
      GEO::vec3 &p = f_mesh.vertices.point(i);
      p[0] = samples[feature_id][i].x;
      p[1] = samples[feature_id][i].y;
      p[2] = 0;
    }
    f_mesh.facets.create_triangles(samples[feature_id].size() - 1);
    for (int i = 0; i < samples[feature_id].size() - 1; i++) {
      f_mesh.facets.set_vertex(i, 0, i);
      f_mesh.facets.set_vertex(i, 1, i + 1);
      f_mesh.facets.set_vertex(i, 2, i + 1);
    }
    f_mesh.facets.compute_borders();
    f_trees.push_back(std::make_shared<GEO::MeshFacetsAABB>(f_mesh));
  }

  std::vector<bool> f_is_removed(features.size(), false);
  std::vector<std::vector<bool>> v_is_removed(features.size(),
                                              std::vector<bool>());
  for (int i = 0; i < features.size(); i++)
    v_is_removed[i] = std::vector<bool>(samples[i].size(), false);
  for (auto &pair : feature_pairs) {
    if (f_is_removed[pair[0]] || f_is_removed[pair[1]])
      continue;

    GEO::MeshFacetsAABB &f_tree = *f_trees[pair[0]];
    int feature_id = pair[1];
    int other_feature_id = pair[0];
    auto feature = features[feature_id];

    GEO::vec3 nearest_point;
    double sq_dist;
    GEO::index_t prev_facet;
    int cnt = 0;
    std::vector<bool> cur_v_is_removed(v_is_removed[feature_id].size(), false);
    std::vector<double> cur_dist(v_is_removed[feature_id].size(), -1);
    for (int i = 0; i < samples[feature_id].size(); i++) {
      if (v_is_removed[feature_id][i])
        continue;

      GEO::vec3 v(samples[feature_id][i].x, samples[feature_id][i].y, 0);
      if (cnt == 0) {
        prev_facet = f_tree.nearest_facet(v, nearest_point, sq_dist);
        if (sq_dist <= feature_eps_2) {
          //                    v_is_removed[feature_id][i] = true;
          cur_v_is_removed[i] = true;
        }
        cur_dist[i] = sq_dist;
        cnt++;
        continue;
      }

      cnt++;
      sq_dist = v.distance2(nearest_point);
      //            if (sq_dist <= feature_eps_2) {
      ////                v_is_removed[feature_id][i] = true;
      //                cur_v_is_removed[i] = true;
      //                continue;
      //            }
      f_tree.nearest_facet_with_hint(v, prev_facet, nearest_point, sq_dist);
      if (sq_dist <= feature_eps_2) {
        //                v_is_removed[feature_id][i] = true;
        cur_v_is_removed[i] = true;
      }
      cur_dist[i] = sq_dist;
    }

    cnt = std::count(cur_v_is_removed.begin(), cur_v_is_removed.end(), true);
    if (cnt == 0)
      continue;
    if (cnt == cur_v_is_removed.size()) {
      f_is_removed[feature_id] = true;
      continue;
    }

    // preserve corners
    bool is_corner = true;
    //        if(cur_v_is_removed.front() && !v_is_removed[feature_id].front())
    //        {
    if (cur_v_is_removed.front()) {
      //            cout<<"is_corner 1"<<endl;
      Vector_2f v1 = samples[feature_id][1] - samples[feature_id][0];
      Vector_2f v2;
      if (samples[feature_id].front() == samples[other_feature_id].front())
        v2 = samples[other_feature_id][1] - samples[other_feature_id][0];
      else if (samples[feature_id].front() == samples[other_feature_id].back())
        v2 = samples[other_feature_id][samples[other_feature_id].size() - 2] -
             samples[other_feature_id].back();
      else
        is_corner = false;

      if (is_corner) {
        if (v1.dot(v2) / (v1.length() * v2.length()) < cos_angle_threshold) {
          //                    cout << "is_corner" << endl;
          //                    for (int j = 0; j < cur_v_is_removed.size();
          //                    j++)
          //                        cout << cur_v_is_removed[j];
          //                    cout << "---------" << endl;
          bool is_valid = true;
          for (int j = 1; j < cur_v_is_removed.size(); j++) {
            if (!cur_v_is_removed[j])
              break;
            if (cur_dist[j] < cur_dist[j - 1])
              is_valid = false;
          }
          if (is_valid) {
            for (int j = 0; j < cur_v_is_removed.size(); j++) {
              //                        cout<<j<<" "<<cur_v_is_removed[j]<<"
              //                        "<<v_is_removed[feature_id][j]<<endl;
              if (!cur_v_is_removed[j] /* || v_is_removed[feature_id][j]*/)
                break;
              cur_v_is_removed[j] = false;
              //                        v_is_removed[feature_id][j] = false;
            }
          }
          //                    cout<<endl;
          //                    for (int j = 0; j < cur_v_is_removed.size();
          //                    j++)
          //                        cout << cur_v_is_removed[j];
          //                    cout << "---------" << endl;
          //                    optimization::pausee();
        }
      }
      //        } else if(cur_v_is_removed.back() &&
      //        !v_is_removed[feature_id].back()) {
    } else if (cur_v_is_removed.back()) {
      //            cout<<"is_corner 2"<<endl;
      Vector_2f v1 = samples[feature_id][samples[feature_id].size() - 2] -
                     samples[feature_id].back();
      Vector_2f v2;
      if (samples[feature_id].back() == samples[other_feature_id].front())
        v2 = samples[other_feature_id][1] - samples[other_feature_id][0];
      else if (samples[feature_id].back() == samples[other_feature_id].back())
        v2 = samples[other_feature_id][samples[other_feature_id].size() - 2] -
             samples[other_feature_id].back();
      else
        is_corner = false;

      if (is_corner) {
        if (v1.dot(v2) / (v1.length() * v2.length()) < cos_angle_threshold) {
          //                    cout<<"is_corner"<<endl;
          //                    for (int j = 0; j < cur_v_is_removed.size();
          //                    j++)
          //                        cout << cur_v_is_removed[j];
          //                    cout << "---------" << endl;
          bool is_valid = true;
          for (int j = cur_v_is_removed.size() - 2; j >= 0; j--) {
            if (!cur_v_is_removed[j])
              break;
            if (cur_dist[j] < cur_dist[j + 1])
              is_valid = false;
          }
          if (is_valid) {
            for (int j = cur_v_is_removed.size() - 1; j >= 0; j--) {
              //                        cout<<j<<" "<<cur_v_is_removed[j]<<"
              //                        "<<v_is_removed[feature_id][j]<<endl;
              if (!cur_v_is_removed[j] /* || v_is_removed[feature_id][j]*/)
                break;
              cur_v_is_removed[j] = false;
              //                        v_is_removed[feature_id][j] = false;
            }
          }
          //                    cout<<endl;
          //                    for (int j = 0; j < cur_v_is_removed.size();
          //                    j++)
          //                        cout << cur_v_is_removed[j];
          //                    cout << "---------" << endl;
          //                    optimization::pausee();
        }
      }
    }

    if (std::count(cur_v_is_removed.begin(), cur_v_is_removed.end(), true) == 0)
      continue;

    // add current result
    for (int i = 0; i < cur_v_is_removed.size(); i++) {
      v_is_removed[feature_id][i] =
          v_is_removed[feature_id][i] || cur_v_is_removed[i];
    }

    // remove single vertices
    for (int i = 0; i < samples[feature_id].size(); i++) {
      if (i == 0) {
        if (!v_is_removed[feature_id][i] && v_is_removed[feature_id][i + 1])
          v_is_removed[feature_id][i] = true;
      } else if (i == samples[feature_id].size() - 1) {
        if (!v_is_removed[feature_id][i] && v_is_removed[feature_id][i - 1])
          v_is_removed[feature_id][i] = true;
      } else {
        if (!v_is_removed[feature_id][i] && v_is_removed[feature_id][i + 1] &&
            v_is_removed[feature_id][i - 1])
          v_is_removed[feature_id][i] = true;
      }
    }
    if (std::count(v_is_removed[feature_id].begin(),
                   v_is_removed[feature_id].end(),
                   true) == v_is_removed[feature_id].size()) {
      f_is_removed[feature_id] = true;
      continue;
    }

    // update aabb tree
    GEO::Mesh &f_mesh = *(f_meshes[feature_id]);
    f_mesh.vertices.clear();
    f_mesh.vertices.create_vertices(
        samples[feature_id]
            .size()); // the old vertices has been re-ordered!!!!!!!!!
    for (int i = 0; i < samples[feature_id].size(); i++) {
      GEO::vec3 &p = f_mesh.vertices.point(i);
      p[0] = samples[feature_id][i].x;
      p[1] = samples[feature_id][i].y;
      p[2] = 0;
    }
    f_mesh.facets.clear();
    cnt = 0;
    for (int i = 0; i < samples[feature_id].size() - 1; i++) {
      if (v_is_removed[feature_id][i] || v_is_removed[feature_id][i + 1])
        continue;
      cnt++;
    }
    f_mesh.facets.create_triangles(cnt);
    cnt = 0;
    for (int i = 0; i < samples[feature_id].size() - 1; i++) {
      if (v_is_removed[feature_id][i] || v_is_removed[feature_id][i + 1])
        continue;
      f_mesh.facets.set_vertex(cnt, 0, i);
      f_mesh.facets.set_vertex(cnt, 1, i + 1);
      f_mesh.facets.set_vertex(cnt, 2, i + 1);
      cnt++;
    }
    f_mesh.facets.compute_borders();
    f_trees[feature_id] = std::make_shared<GEO::MeshFacetsAABB>(f_mesh);
  }

  for (int i = 0; i < features.size(); i++) {
    if (f_is_removed[i]) {
      push_back_new_secondary_feature(features[i]);

      features.erase(features.begin() + i);
      inflections.erase(inflections.begin() + i);
      samples.erase(samples.begin() + i);
      ts.erase(ts.begin() + i);
      f_is_removed.erase(f_is_removed.begin() + i);
      v_is_removed.erase(v_is_removed.begin() + i);
      i--;
    }
  }

  /// enlarge small gaps
  for (int i = 0; i < features.size(); i++) {
    double l = 0;
    for (int j = 0; j < v_is_removed[i].size() - 1; j++) {
      if (!v_is_removed[i][j] && v_is_removed[i][j + 1]) {
        l = (samples[i][j] - samples[i][j + 1]).length();
      } else if (v_is_removed[i][j] && v_is_removed[i][j + 1]) {
        if (l > 0)
          l += (samples[i][j] - samples[i][j + 1]).length();
      } else if (v_is_removed[i][j] && !v_is_removed[i][j + 1]) {
        if (l > 0) {
          l += (samples[i][j] - samples[i][j + 1]).length();
          if (l < feature_eps) {
            for (int k = j; k < v_is_removed.size() - 1; k++) {
              l += (samples[i][k] - samples[i][k + 1]).length();
              v_is_removed[i][k + 1] = true;
              if (l >= feature_eps) {
                cout << "enlarged " << k - j + 1 << endl;
                break;
              }
            }
          }
        }
        l = 0;
      }
    }
    int cnt = std::count(v_is_removed[i].begin(), v_is_removed[i].end(), false);
    if (cnt == 0) {
      cout << "cnt == 0" << endl;
      push_back_new_secondary_feature(features[i]);

      features.erase(features.begin() + i);
      inflections.erase(inflections.begin() + i);
      samples.erase(samples.begin() + i);
      ts.erase(ts.begin() + i);
      f_is_removed.erase(f_is_removed.begin() + i);
      v_is_removed.erase(v_is_removed.begin() + i);
      i--;
    }
  }

  /// split features
  const int features_size = features.size();
  for (int i = 0; i < features_size; i++) {
    std::vector<int> ranges; // local ids of samples
    if (!v_is_removed[i].front())
      ranges.push_back(0);

    for (int j = 0; j < v_is_removed[i].size() - 1; j++) {
      if (v_is_removed[i][j] && !v_is_removed[i][j + 1])
        ranges.push_back(j + 1);
      else if (!v_is_removed[i][j] && v_is_removed[i][j + 1])
        ranges.push_back(j);
    }
    if (!v_is_removed[i].back())
      ranges.push_back(v_is_removed[i].size() - 1);

    if (ranges.size() % 2 != 0) {
      for (bool j : v_is_removed[i])
        cout << j;
      cout << endl;
      for (auto j : ranges)
        cout << j << " ";
      cout << endl;
    }
    assert(ranges.size() % 2 == 0);

    // add secondary features
    std::vector<int> removed_ranges = ranges;
    if (removed_ranges.front() == 0)
      removed_ranges.erase(removed_ranges.begin());
    else
      removed_ranges.insert(removed_ranges.begin(), 0);
    if (removed_ranges.back() == v_is_removed[i].size() - 1)
      removed_ranges.erase(removed_ranges.begin() + removed_ranges.size() - 1);
    else
      removed_ranges.push_back(v_is_removed[i].size() - 1);
    for (int j = 0; j < removed_ranges.size(); j += 2) {
      int new_feature_id = push_back_new_secondary_feature(features[i]);
      secondary_features[new_feature_id]->paras = {
          ts[i][removed_ranges[j]], ts[i][removed_ranges[j + 1]]};
    }

    // for (int j = 0; j < ranges.size(); j += 2) {
    //     if (j == ranges.size() - 2) {
    //         features[i]->paras = {ts[i][ranges[j]], ts[i][ranges[j + 1]]};
    //         samples[i] = std::vector<Point_2f>(samples[i].begin() +
    //         ranges[j],
    //                                            samples[i].begin() + ranges[j
    //                                            + 1] + 1);
    //         ts[i] = std::vector<double>(ts[i].begin() + ranges[j],
    //         ts[i].begin() + ranges[j + 1] + 1); std::vector<double>
    //         new_inflection; for (int infl: inflections[i]) {
    //             if (infl > ranges[j] && infl < ranges[j + 1])
    //                 new_inflection.push_back(infl);
    //         }
    //         inflections[i] = new_inflection;
    //     } else {
    //         int new_feature_id = push_back_new_feature(features[i]);
    //         features[new_feature_id]->paras = {ts[i][ranges[j]],
    //         ts[i][ranges[j + 1]]}; samples.push_back(
    //                 std::vector<Point_2f>(samples[i].begin() + ranges[j],
    //                 samples[i].begin() + ranges[j + 1] + 1));
    //         ts.push_back(std::vector<double>(ts[i].begin() + ranges[j],
    //         ts[i].begin() + ranges[j + 1] + 1)); inflections.emplace_back();
    //         for (int infl: inflections[i]) {
    //             if (infl > ranges[j] && infl < ranges[j + 1])
    //                 inflections.back().push_back(infl);
    //         }
    //     }
    // }
    auto old_inflection = inflections[i];
    for (int j = 0; j < ranges.size(); j += 2) {
      if (j == ranges.size() - 2) {
        features[i]->paras = {ts[i][ranges[j]], ts[i][ranges[j + 1]]};
        samples[i] =
            std::vector<Point_2f>(samples[i].begin() + ranges[j],
                                  samples[i].begin() + ranges[j + 1] + 1);
        ts[i] = std::vector<double>(ts[i].begin() + ranges[j],
                                    ts[i].begin() + ranges[j + 1] + 1);
        std::vector<double> new_inflection;
        for (double infl : old_inflection) {
          if (infl > ts[i][ranges[j]] && infl < ts[i][ranges[j + 1]])
            new_inflection.push_back(infl);
        }
        inflections[i] = new_inflection;
      } else {
        int new_feature_id = push_back_new_feature(features[i]);
        features[new_feature_id]->paras = {ts[i][ranges[j]],
                                           ts[i][ranges[j + 1]]};
        samples.push_back(
            std::vector<Point_2f>(samples[i].begin() + ranges[j],
                                  samples[i].begin() + ranges[j + 1] + 1));
        ts.push_back(std::vector<double>(ts[i].begin() + ranges[j],
                                         ts[i].begin() + ranges[j + 1] + 1));
        inflections.emplace_back();
        for (double infl : old_inflection) {
          if (infl > ts[i][ranges[j]] && infl < ts[i][ranges[j + 1]])
            inflections.back().push_back(infl);
        }
      }
    }
  }

  //    //check
  //    std::vector<int> all;
  //    for(int i=0;i<features.size();i++){
  ////        if(features[i]->type=="Line")
  //            all.push_back(i);
  //    }
  //    for(int i=0;i<features.size();i++){
  ////        if(features[i]->type!="Line")
  ////            continue;
  //        Point_2f p_1 = features[i]->eval(features[i]->paras.front());
  //        Point_2f p_2 = features[i]->eval(features[i]->paras.back());
  //        for(int j=0;j<all.size();j++) {
  //            if(all[j]==i)
  //                continue;
  //            Point_2f p1 =
  //            features[all[j]]->eval(features[all[j]]->paras.front());
  //            Point_2f p2 =
  //            features[all[j]]->eval(features[all[j]]->paras.back()); double
  //            l11 = (p_1-p1).length(); double l12 = (p_1-p2).length(); double
  //            l21 = (p_2-p1).length(); double l22 = (p_2-p2).length();
  //            if(l11<feature_eps || l12<feature_eps || l21<feature_eps ||
  //            l22<feature_eps){
  //                cout<<"feature "<<features[i]->curve_id<<"
  //                "<<features[all[j]]->curve_id<<endl;
  //                cout<<std::setprecision(16)<<endl;
  //                cout<<p_1<<"; "<<p_2<<endl;
  //                cout<<p1<<"; "<<p2<<endl;
  //                cout<<l11<<" "<<l12<<" "<<l21<<" "<<l22<<endl;
  //                optimization::pausee();
  //            }
  //        }
  //    }
}

void triwild::feature::remove_high_curvature(
    std::vector<std::vector<double>> &inflections,
    std::vector<std::vector<Point_2f>> &samples,
    std::vector<std::vector<double>> &ts) {
  int cnt_removed = 0;
  for (int feature_id = 0; feature_id < features.size(); feature_id++) {
    if (features[feature_id]->type == "Line")
      continue;

    std::vector<bool> seg_is_removed(samples[feature_id].size() - 1, false);
    for (int i = 0; i < samples[feature_id].size() - 1; i++) {
      for (int j = 0; j < inflections[feature_id].size(); j++) {
        if (inflections[feature_id][j] > ts[feature_id][i] &&
            inflections[feature_id][j] < ts[feature_id][i + 1])
          continue;
      }
      double angle = features[feature_id]->how_curve(ts[feature_id][i],
                                                     ts[feature_id][i + 1]);
      if (angle > args.flat_feature_angle) {
        seg_is_removed[i] = true;
        cnt_removed++;
      }
    }

    int cnt = std::count(seg_is_removed.begin(), seg_is_removed.end(), true);
    if (cnt == 0)
      continue;
    else if (cnt == seg_is_removed.size()) {
      push_back_new_secondary_feature(features[feature_id]);
      features.erase(features.begin() + feature_id);
      inflections.erase(inflections.begin() + feature_id);
      samples.erase(samples.begin() + feature_id);
      ts.erase(ts.begin() + feature_id);
      feature_id--;
      continue;
    }

    std::vector<int> ranges; // ids of samples
    if (!seg_is_removed.front())
      ranges.push_back(0);
    for (int j = 0; j < seg_is_removed.size() - 1; j++) {
      if (seg_is_removed[j] && !seg_is_removed[j + 1])
        ranges.push_back(j + 1);
      else if (!seg_is_removed[j] && seg_is_removed[j + 1])
        //                ranges.push_back(j);
        ranges.push_back(j + 1);
    }
    if (!seg_is_removed.back())
      ranges.push_back(samples[feature_id].size() - 1);

    if (ranges.size() % 2 != 0) {
      for (bool j : seg_is_removed)
        cout << j;
      cout << endl;
      for (auto j : ranges)
        cout << j << " ";
      cout << endl;
    }
    assert(ranges.size() % 2 == 0);

    // add secondary features
    std::vector<int> removed_ranges = ranges;
    if (removed_ranges.front() == 0)
      removed_ranges.erase(removed_ranges.begin());
    else
      removed_ranges.insert(removed_ranges.begin(), 0);
    if (removed_ranges.back() == samples[feature_id].size() - 1)
      removed_ranges.erase(removed_ranges.begin() + removed_ranges.size() - 1);
    else
      removed_ranges.push_back(samples[feature_id].size() - 1);
    for (int j = 0; j < removed_ranges.size(); j += 2) {
      int new_feature_id =
          push_back_new_secondary_feature(features[feature_id]);
      secondary_features[new_feature_id]->paras = {
          ts[feature_id][removed_ranges[j]],
          ts[feature_id][removed_ranges[j + 1]]};
    }

    // for (int j = 0; j < ranges.size(); j += 2) {
    //     if (j == ranges.size() - 2) {
    //         features[feature_id]->paras = {ts[feature_id][ranges[j]],
    //         ts[feature_id][ranges[j + 1]]}; samples[feature_id] =
    //         std::vector<Point_2f>(samples[feature_id].begin() + ranges[j],
    //                                                     samples[feature_id].begin()
    //                                                     + ranges[j + 1] + 1);
    //         ts[feature_id] = std::vector<double>(ts[feature_id].begin() +
    //         ranges[j],
    //                                              ts[feature_id].begin() +
    //                                              ranges[j + 1] + 1);
    //         std::vector<double> new_inflection;
    //         for (int infl: inflections[feature_id]) {
    //             if (infl > ranges[j] && infl < ranges[j + 1])
    //                 new_inflection.push_back(infl);
    //         }
    //         inflections[feature_id] = new_inflection;
    //     } else {
    //         int new_feature_id = push_back_new_feature(features[feature_id]);
    //         features[new_feature_id]->paras = {ts[feature_id][ranges[j]],
    //         ts[feature_id][ranges[j + 1]]}; samples.push_back(
    //                 std::vector<Point_2f>(samples[feature_id].begin() +
    //                 ranges[j],
    //                                       samples[feature_id].begin() +
    //                                       ranges[j + 1] + 1));
    //         ts.push_back(std::vector<double>(ts[feature_id].begin() +
    //         ranges[j],
    //                                          ts[feature_id].begin() +
    //                                          ranges[j + 1] + 1));
    //         inflections.emplace_back();
    //         for (int infl: inflections[feature_id]) {
    //             if (infl > ranges[j] && infl < ranges[j + 1])
    //                 inflections.back().push_back(infl);
    //             //todo: infl == ranges[j]
    //         }
    //     }
    // }
    int i = feature_id;
    auto old_inflection = inflections[i];
    for (int j = 0; j < ranges.size(); j += 2) {
      if (j == ranges.size() - 2) {
        features[i]->paras = {ts[i][ranges[j]], ts[i][ranges[j + 1]]};
        samples[i] =
            std::vector<Point_2f>(samples[i].begin() + ranges[j],
                                  samples[i].begin() + ranges[j + 1] + 1);
        ts[i] = std::vector<double>(ts[i].begin() + ranges[j],
                                    ts[i].begin() + ranges[j + 1] + 1);
        std::vector<double> new_inflection;
        for (double infl : old_inflection) {
          if (infl > ts[i][ranges[j]] && infl < ts[i][ranges[j + 1]])
            new_inflection.push_back(infl);
        }
        inflections[i] = new_inflection;
      } else {
        int new_feature_id = push_back_new_feature(features[i]);
        features[new_feature_id]->paras = {ts[i][ranges[j]],
                                           ts[i][ranges[j + 1]]};
        samples.push_back(
            std::vector<Point_2f>(samples[i].begin() + ranges[j],
                                  samples[i].begin() + ranges[j + 1] + 1));
        ts.push_back(std::vector<double>(ts[i].begin() + ranges[j],
                                         ts[i].begin() + ranges[j + 1] + 1));
        inflections.emplace_back();
        for (double infl : old_inflection) {
          if (infl > ts[i][ranges[j]] && infl < ts[i][ranges[j + 1]])
            inflections.back().push_back(infl);
        }
      }
    }
  }

  cout << cnt_removed << " segments are removed" << endl;
}

void triwild::feature::cut_inflections(
    std::vector<std::vector<double>> &inflections,
    std::vector<std::vector<Point_2f>> &samples,
    std::vector<std::vector<double>> &ts) {
  const int features_size = features.size();
  for (int feature_id = 0; feature_id < features_size; feature_id++) {
    if (inflections[feature_id].empty())
      continue;

    int cnt = 0;
    std::vector<int> ranges;
    int j = 0;
    for (double infl : inflections[feature_id]) {
      for (; j < ts[feature_id].size() - 1; j++) {
        if (ts[feature_id][j] == infl) {
          if (!ranges.empty() && j - ranges.back() < 3)
            continue;
          if (j - 0 < 3)
            continue;
          if (j - (ts[feature_id].size() - 1) < 3)
            continue;

          ranges.push_back(j);
          if (j != 0)
            ranges.push_back(j);
          break;
        } else if (ts[feature_id][j] < infl && ts[feature_id][j + 1] > infl) {
          if (!ranges.empty() && (j + 1) - ranges.back() < 3)
            continue;
          if ((j + 1) - 0 < 3)
            continue;
          if ((j + 1) - (ts[feature_id].size() - 1) < 3)
            continue;

          ts[feature_id].insert(ts[feature_id].begin() + j + 1, infl);
          samples[feature_id].insert(samples[feature_id].begin() + j + 1,
                                     features[feature_id]->eval(infl));
          ranges.push_back(j + 1);
          if (j + 1 != ts[feature_id].size() - 1)
            ranges.push_back(j + 1);
          break;
        }
      }
    }
    if (ranges.empty())
      continue;

    if (ranges.front() != 0)
      ranges.insert(ranges.begin(), 0);
    if (ranges.back() != ts[feature_id].size() - 1)
      ranges.push_back(ts[feature_id].size() - 1);

    //        for (int j = 0; j < ranges.size() - 1; j++) {
    for (int j = 0; j < ranges.size(); j += 2) {
      if (j == ranges.size() - 2) {
        features[feature_id]->paras = {ts[feature_id][ranges[j]],
                                       ts[feature_id][ranges[j + 1]]};
        samples[feature_id] = std::vector<Point_2f>(
            samples[feature_id].begin() + ranges[j],
            samples[feature_id].begin() + ranges[j + 1] + 1);
        ts[feature_id] =
            std::vector<double>(ts[feature_id].begin() + ranges[j],
                                ts[feature_id].begin() + ranges[j + 1] + 1);

        inflections[feature_id].clear();
        features[feature_id]->is_inflection = {true, false};
      } else {
        int new_feature_id = push_back_new_feature(features[feature_id]);
        features[new_feature_id]->paras = {ts[feature_id][ranges[j]],
                                           ts[feature_id][ranges[j + 1]]};
        samples.push_back(std::vector<Point_2f>(
            samples[feature_id].begin() + ranges[j],
            samples[feature_id].begin() + ranges[j + 1] + 1));
        ts.push_back(
            std::vector<double>(ts[feature_id].begin() + ranges[j],
                                ts[feature_id].begin() + ranges[j + 1] + 1));
        inflections.emplace_back();
        if (j == 0)
          features[new_feature_id]->is_inflection = {false, true};
        else
          features[new_feature_id]->is_inflection = {true, true};
      }
    }

    //        std::vector<double> ranges = inflections[feature_id];
    //        ranges.insert(ranges.begin(),
    //        features[feature_id]->paras.front());
    //        ranges.push_back(features[feature_id]->paras.back());
    //
    //        for (int j = 0; j < ranges.size() - 1; j++) {
    //            if (j == ranges.size() - 2) {
    //                features[feature_id]->paras = {ranges[j], ranges[j + 1]};
    ////                samples[feature_id] =
    /// std::vector<Point_2f>(samples[feature_id].begin() + ranges[j], /
    /// samples[feature_id].begin() + ranges[j + 1] + 1); / ts[feature_id] =
    /// std::vector<double>(ts[feature_id].begin() + ranges[j], /
    /// ts[feature_id].begin() + ranges[j + 1] + 1);
    //
    //                inflections.clear();
    //                features[feature_id]->is_inflection = {true, false};
    //            } else {
    //                int new_feature_id =
    //                push_back_new_feature(features[feature_id]);
    //                features[new_feature_id]->paras = {ranges[j], ranges[j +
    //                1]};
    ////                samples.push_back(
    //// std::vector<Point_2f>(samples[feature_id].begin() + ranges[j], /
    /// samples[feature_id].begin() + ranges[j + 1] + 1)); /
    /// ts.push_back(std::vector<double>(ts[feature_id].begin() + ranges[j], /
    /// ts[feature_id].begin() + ranges[j + 1] + 1));
    //                inflections.emplace_back();
    //                if (j == 0)
    //                    features[new_feature_id]->is_inflection = {false,
    //                    true};
    //                else
    //                    features[new_feature_id]->is_inflection = {true,
    //                    true};
    //            }
    //        }
  }
}

void triwild::feature::gen_segments(Eigen::MatrixXd &V,
                                    std::vector<std::array<int, 2>> &edges) {
  double dd = args.feature_epsilon * args.diagonal_len * 10;
  double dd_2 = dd * dd;

  double min_len_2 = (args.min_edge_length * args.target_edge_len) *
                     (args.min_edge_length * args.target_edge_len);

  std::vector<std::vector<Point_2f>> samples(features.size());
  std::vector<std::vector<Point_2f>> secondary_samples(
      secondary_features.size());
  int cnt_v = 0;
  auto sampling = [&](std::vector<std::shared_ptr<FeatureElement>> &features,
                      std::vector<std::vector<Point_2f>> &samples) {
    for (int i = 0; i < features.size(); i++) {
      if (features[i]->paras.empty())
        features[i]->paras = {0.0, 1.0};
      samples[i] = {features[i]->eval(features[i]->paras.front()),
                    features[i]->eval(features[i]->paras.back())};
      if (features[i]->type == "Line") {
        cnt_v += samples[i].size();
        continue;
      }

      std::vector<bool> is_qualified = {false};
      while (true) {
        // check length
        int cnt = 0;
        for (int j = 0; j < is_qualified.size(); j++) {
          if (is_qualified[j])
            continue;
          double l_2 = (samples[i][j] - samples[i][j + 1]).length_2();
          if (l_2 <= dd_2 &&
              features[i]->how_curve(features[i]->paras[j],
                                     features[i]->paras[j + 1]) <=
                  args.flat_feature_angle)
            is_qualified[j] = true;
          //                    else if (l_2 < min_len_2) {//should never happen
          //                        is_qualified[j] = true;
          //                    }
          else
            cnt++;
        }
        if (cnt == 0)
          break;

        // insert samples
        for (int j = 0; j < samples[i].size() - 1; j++) {
          if (is_qualified[j])
            continue;
          double t = (features[i]->paras[j] + features[i]->paras[j + 1]) / 2;
          features[i]->paras.insert(features[i]->paras.begin() + j + 1, t);
          samples[i].insert(samples[i].begin() + j + 1, features[i]->eval(t));
          is_qualified.insert(is_qualified.begin() + j + 1, false);
          j++;
        }
      }

      cnt_v += samples[i].size();
    }
  };

  sampling(features, samples);
  sampling(secondary_features, secondary_samples);

  V.resize(cnt_v, 2);
  edges.clear();
  edges.resize(cnt_v - features.size() - secondary_features.size());
  cnt_v = 0;
  int cnt_e = 0;
  for (int i = 0; i < features.size(); i++) {
    for (int j = 0; j < samples[i].size(); j++) {
      features[i]->v_ids.push_back(cnt_v);
      V.row(cnt_v) << samples[i][j].x, samples[i][j].y;
      if (j > 0) {
        edges[cnt_e] = {{cnt_v - 1, cnt_v}};
        cnt_e++;
      }
      cnt_v++;
    }
  }
  for (int i = 0; i < secondary_features.size(); i++) {
    for (int j = 0; j < secondary_samples[i].size(); j++) {
      secondary_features[i]->v_ids.push_back(cnt_v);
      V.row(cnt_v) << secondary_samples[i][j].x, secondary_samples[i][j].y;
      if (j > 0) {
        edges[cnt_e] = {{cnt_v - 1, cnt_v}};
        cnt_e++;
      }
      cnt_v++;
    }
  }

  Eigen::VectorXi VI, _;
  Eigen::MatrixXd V_tmp;
  igl::unique_rows(V, V_tmp, _, VI);
  V = V_tmp;

  // duplicate/degenerate edges
  for (int i = 0; i < edges.size(); i++) {
    edges[i][0] = VI(edges[i][0]);
    edges[i][1] = VI(edges[i][1]);
    if (edges[i][0] == edges[i][1]) {
      edges.erase(edges.begin() + i);
      i--;
    } else if (edges[i][0] > edges[i][1])
      std::swap(edges[i][0], edges[i][1]);
  }
  optimization::vector_unique(edges);

  for (int i = 0; i < features.size(); i++) {
    for (int j = 0; j < features[i]->v_ids.size(); j++) {
      features[i]->v_ids[j] = VI(features[i]->v_ids[j]);
    }
  }
  for (int i = 0; i < secondary_features.size(); i++) {
    for (int j = 0; j < secondary_features[i]->v_ids.size(); j++) {
      secondary_features[i]->v_ids[j] = VI(secondary_features[i]->v_ids[j]);
    }
  }
}

int triwild::feature::push_back_new_feature(
    std::shared_ptr<FeatureElement> old_feature) {
  if (old_feature->type == "Line")
    features.push_back(std::make_shared<Line_Feature>(
        *dynamic_cast<Line_Feature *>(old_feature.get())));
  else if (old_feature->type == "RationalBezierCurve")
    features.push_back(std::make_shared<RationalBezierCurve_Feature>(
        *dynamic_cast<RationalBezierCurve_Feature *>(old_feature.get())));
  else if (old_feature->type == "BezierCurve")
    features.push_back(std::make_shared<BezierCurve_Feature>(
        *dynamic_cast<BezierCurve_Feature *>(old_feature.get())));

  return features.size() - 1;
}

int triwild::feature::push_back_new_secondary_feature(
    std::shared_ptr<FeatureElement> old_feature) {
  if (old_feature->type == "Line")
    secondary_features.push_back(std::make_shared<Line_Feature>(
        *dynamic_cast<Line_Feature *>(old_feature.get())));
  else if (old_feature->type == "RationalBezierCurve")
    secondary_features.push_back(std::make_shared<RationalBezierCurve_Feature>(
        *dynamic_cast<RationalBezierCurve_Feature *>(old_feature.get())));
  else if (old_feature->type == "BezierCurve")
    secondary_features.push_back(std::make_shared<BezierCurve_Feature>(
        *dynamic_cast<BezierCurve_Feature *>(old_feature.get())));

  return secondary_features.size() - 1;
}

#include <igl/writeSTL.h>
void triwild::feature::output_features(const std::string &name) {
  std::vector<std::vector<Point_2f>> samples;
  std::vector<std::vector<double>> ts;
  sample_features(samples, ts);

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  int cnt = 0;
  for (int i = 0; i < features.size(); i++) {
    cnt += samples[i].size();
  }
  V.resize(cnt, 3);
  F.resize(cnt - features.size(), 3);
  int cnt_v = 0;
  int cnt_f = 0;
  for (int i = 0; i < samples.size(); i++) {
    for (int j = 0; j < samples[i].size(); j++) {
      V.row(cnt_v) << samples[i][j].x, samples[i][j].y, 0;
      if (j < samples[i].size() - 1) {
        F.row(cnt_f) << cnt_v, cnt_v, cnt_v + 1;
        cnt_f++;
      }
      cnt_v++;
    }
  }
  igl::writeSTL(args.output + name + ".stl", V, F);
}

void triwild::feature::check(std::vector<std::vector<double>> &inflections,
                             std::vector<std::vector<Point_2f>> &samples,
                             std::vector<std::vector<double>> &ts) {
  for (int i = 0; i < features.size(); i++) {
    if (samples[i].size() < 2) {
      cout << "feature " << i << " samples[i].size()<2 " << endl;
      optimization::pausee();
    }
    if (features[i]->paras.front() != ts[i].front()) {
      cout << "features[i]->paras.front() != ts[i].front()" << endl;
      optimization::pausee();
    }
    if (features[i]->paras.back() != ts[i].back()) {
      cout << "features[i]->paras.back() != ts[i].back()" << endl;
      optimization::pausee();
    }
    if (features[i]->eval(ts[i].front()) != samples[i].front()) {
      cout << "features[i]->eval(ts[i].front())!=samples[i].front()" << endl;
      optimization::pausee();
    }
    if (features[i]->eval(ts[i].back()) != samples[i].back()) {
      cout << "features[i]->eval(ts[i].back())!=samples[i].back()" << endl;
      optimization::pausee();
    }
    if (features[i]->paras.back() <= features[i]->paras.front()) {
      cout << "features[i]->paras.back()<=features[i]->paras.front()" << endl;
      cout << "feature " << i << endl;
      cout << features[i]->paras.front() << " " << features[i]->paras.back()
           << endl;
      optimization::pausee();
    }
    for (int j = 0; j < ts[i].size() - 1; j++) {
      if (ts[i][j] >= ts[i][j + 1]) {
        cout << "ts[i][j]>=ts[i][j+1]" << endl;
        cout << "feature " << i << endl;
        for (auto k : ts[i])
          cout << k << " ";
        cout << endl;
        optimization::pausee();
      }
    }
  }

  for (int i = 0; i < secondary_features.size(); i++) {
    if (secondary_features[i]->paras.size() < 2) {
      cout << "secondary_features[i]->paras.size()<2" << endl;
      optimization::pausee();
    }
    if (secondary_features[i]->paras.back() <=
        secondary_features[i]->paras.front()) {
      cout << "secondary_features[i]->paras.back()<=secondary_features[i]->"
              "paras.front()"
           << endl;
      optimization::pausee();
    }
  }
}
