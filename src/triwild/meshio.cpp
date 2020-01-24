// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "meshio.hpp"

#include "optimization.h"
#include "edge_splitting.h"
#include "feature.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <unordered_map>

#include <igl/readOBJ.h>
#include <geogram/basic/geometry.h>
#include <geogram/mesh/mesh_topology.h>
#include <geogram/mesh/mesh_geometry.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/delaunay/delaunay.h>

#include <nlohmann/json.hpp>
using json = nlohmann::json;


namespace triwild
{

	bool unordered_set_intersection_own(const std::unordered_set<int> &A, const std::unordered_set<int> &B, std::vector<int> &C)
	{
		C.clear();
		std::vector<int> A_, B_;
		for(auto &i: A)
			A_.push_back(i);
		for(auto &i: B)
			B_.push_back(i);
		std::sort(A_.begin(), A_.end());
		std::sort(B_.begin(), B_.end());
		std::set_intersection(A_.begin(), A_.end(), B_.begin(), B_.end(), std::back_inserter(C));

		return true;
	}

	void write_OBJ( Eigen::MatrixXd V, Eigen::MatrixXi F, const std::string & path)
	{
		std::cout << "Writing \"" << path << "\" (V=" << V.cols()
			<< ", F=" << F.cols() << ") .. ";
		std::cout.flush();
		std::ofstream os(path);
		if (os.fail())
			throw std::runtime_error("Unable to open OBJ file \"" + path + "\"!");

//		std::cout.precision(dbl::max_digits10);
		os<<std::setprecision(16);

		for (int i = 0; i<V.cols(); i++){
			os <<"v ";
			for (int j = 0; j<V.rows(); j++)
				os<<std::fixed<< V(j, i) << " ";
			os<< 0 << " ";
			os << std::endl;
		}

		for (int i = 0; i < F.cols(); i++) {
			os << "f ";
			for (int j = 0; j < F.rows(); j++)
				os << F(j, i) +1 << " ";
			os << std::endl;
		}
		os.close();
	}

	void write_msh(const MeshData& mesh, const std::string &path, const bool export_edge_tag)
	{
		std::ofstream os(path);
		if (os.fail())
			throw std::runtime_error("Unable to open MSH file \"" + path + "\"!");

		Eigen::MatrixXd V;
		Eigen::MatrixXi F;

		std::unordered_map<int, int> map_v_ids, map_n_ids;
		std::unordered_map<int, int> small_to_big_v;
		int cnt = 0;
		for (size_t i = 0; i < mesh.tri_vertices.size(); i++) {
			if (mesh.v_is_removed[i])
				continue;
			small_to_big_v[cnt] = i;
			map_v_ids[i] = cnt++;

		}

		for (size_t i = 0; i < mesh.nodes.size(); i++) {
			if (mesh.n_is_removed[i])
				continue;
			map_n_ids[i] = cnt++;
		}

		V.resize(cnt, 3);
		cnt = 0;
		for (size_t i = 0; i < mesh.tri_vertices.size(); i++) {
			if (mesh.v_is_removed[i])
				continue;
			V(cnt, 0) = mesh.tri_vertices[i].posf[0];
			V(cnt, 1) = mesh.tri_vertices[i].posf[1];
			V(cnt, 2) = 0;
			cnt++;
		}

		for (size_t i = 0; i < mesh.nodes.size(); i++) {
			if (mesh.n_is_removed[i])
				continue;
			V(cnt, 0) = mesh.nodes[i][0];
			V(cnt, 1) = mesh.nodes[i][1];
			V(cnt, 2) = 0;
			cnt++;
		}

		cnt = std::count(mesh.t_is_removed.begin(), mesh.t_is_removed.end(), false);

		F.resize(cnt, 3);
		std::vector<std::vector<int>> tri_nodes(cnt);

		cnt = 0;
		for (size_t i = 0; i < mesh.tris.size(); i++) {
			// std::cout<<mesh.tri_nodes[i].size()<<std::endl;
			if (mesh.t_is_removed[i])
				continue;

			for (int j = 0; j < 3; j++)
				F(cnt, j) = map_v_ids[mesh.tris[i][j]];

			if(!mesh.tri_nodes.empty())
			{
				const size_t n_nodes = mesh.tri_nodes[i].size();

				tri_nodes[cnt].reserve(n_nodes);
				for (int j = 0; j < n_nodes; j++)
					tri_nodes[cnt].push_back(map_n_ids[mesh.tri_nodes[i][j]]);
			}
			cnt++;
		}



		if(export_edge_tag)
		{

			typedef std::pair<GEO::index_t, GEO::index_t> Edge;


			GEO::Mesh M;
			// Setup vertices
			M.vertices.create_vertices((int) V.rows());
			for (int i = 0; i < (int) M.vertices.nb(); ++i) {
				GEO::vec3 &p = M.vertices.point(i);
				p[0] = V(i, 0);
				p[1] = V(i, 1);
				p[2] = V(i, 2);
			}

			// Setup faces
			assert(F.cols() == 3);
			M.facets.create_triangles((int) F.rows());

			for (int c = 0; c < (int) M.facets.nb(); ++c) {
				for (int lv = 0; lv < F.cols(); ++lv) {
					M.facets.set_vertex(c, lv, F(c, lv));
				}
			}

			M.facets.connect();
			M.cells.connect();

			// Compute a list of all the edges, and store edge index as a corner attribute
			std::vector<std::pair<Edge, GEO::index_t>> e2c; // edge to corner id
			for (GEO::index_t f = 0; f < M.facets.nb(); ++f) {
				for (GEO::index_t c = M.facets.corners_begin(f); c < M.facets.corners_end(f); ++c) {
					GEO::index_t v = M.facet_corners.vertex(c);
					GEO::index_t c2 = M.facets.next_corner_around_facet(f, c);
					GEO::index_t v2 = M.facet_corners.vertex(c2);
					e2c.emplace_back(std::make_pair(std::min(v, v2), std::max(v, v2)), c);
				}
			}
			std::sort(e2c.begin(), e2c.end());

			// Assign unique id to edges
			// GEO::Attribute<GEO::index_t> c2e(M.facet_corners.attributes(), "edge_id");
			M.edges.clear();
			Edge prev_e(-1, -1);
			GEO::index_t current_id = -1;
			std::vector<bool> boundary_edges;
			for (const auto &kv : e2c) {
				Edge e = kv.first;
				GEO::index_t c = kv.second;
				if (e != prev_e) {
					M.edges.create_edge(e.first, e.second);
					boundary_edges.push_back(true);
					++current_id;
					prev_e = e;
				} else {
					boundary_edges.back() = false;
				}
				// c2e[c] = current_id;
			}

//			std::ofstream tags(path + ".txt");
//			if(!tags.fail())
//			{
//				for(int e = 0; e < M.edges.nb(); ++e)
//				{
//
//					const int v0 = small_to_big_v[M.edges.vertex(e, 0)];
//					const int v1 = small_to_big_v[M.edges.vertex(e, 1)];
//					const int edge_tag = triwild::optimization::get_feature_edge_tag(mesh, v0, v1);
//
//					if(edge_tag >= 0)
//						tags<<feature::features[edge_tag]->curve_id + 1<<"\n";
//					else if(triwild::optimization::is_bbox_edge(mesh, v0, v1))
//						tags<<"9999999999 \n";
//					else{
//						const int secondary_edge_tag = triwild::optimization::get_secondary_feature_edge_tag(mesh, v0, v1);
//						if(secondary_edge_tag >= 0)
//							tags<<feature::secondary_features[secondary_edge_tag]->curve_id + 1<<"\n";
//						else
//							tags<<edge_tag<<"\n";
//					}
//				}
//			}
//
//			tags.close();
		}





		const auto n_nodes = V.rows();
		const auto n_elemens = F.rows();

		os << std::setprecision(16)<< "$MeshFormat\n";
		os << "2.2 0 8\n";
		os << "$EndMeshFormat\n";

		os << "$Nodes\n";
		os << n_nodes <<"\n";
		for (size_t i = 0; i<n_nodes; ++i)
		{
			os << (i+1);
			for (int j = 0; j < V.cols(); j++)
				os << " " << V(i, j);
			os << "\n";
		}
		os << "$EndNodes\n";

		os << "$Elements\n";
		os << n_elemens <<"\n";
		for (size_t i = 0; i < n_elemens; ++i) {
			const auto &loc_nodes = tri_nodes[i];
			// el_type:
			//		2: 3-node triangle
			//		9: 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges)
			//		21: 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
			//		23: 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)

			//default: 2: 3-node triangle
			int el_type = 2;
			switch(loc_nodes.size())
			{
				//2: 3-node triangle
				case 0: el_type = 2; break;

				//9: 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges)
				case 3: el_type = 9; break;

				//21: 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
				case 7: el_type = 21; break;

				//23: 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)
				case 12: el_type = 23; break;

				default: std::cerr<<"unspported number of nodes "<<tri_nodes[i].size()<<std::endl; assert(false);
			}

			os << (i+1) << " " << el_type << " 0"; //el_id el_type n_tags

			//first 3 are always the vertices
			for (int j = 0; j < F.cols(); j++)
				os << " " << F(i, j) + 1;

			// follow for node order http://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_65.php#SEC65
			for (size_t j = 0; j < loc_nodes.size(); j++)
			{
				os << " " << loc_nodes[j] + 1;
			}


			os<<"\n";

		}
		os << "$EndElements\n";


		if(!mesh.v_scalars.empty())
		{
			os << "$NodeData\n";
			os <<"1\n";
		    os << "\"Node-data\"\n";
			os << "1\n"; // num real tags.
		    os << "0.0\n"; // time value.
		    os << "3\n"; // num int tags.
		    os << "0\n"; // the time step
		    os << "1\n"; // 1-component scalar field.
		    os << n_nodes << "\n";

		    cnt = 1;
			for (size_t i = 0; i < mesh.v_scalars.size(); i++) {
				if (mesh.v_is_removed[i])
					continue;

				os << cnt << " " << mesh.v_scalars[i]<< "\n";
				++cnt;
			}

			for (size_t i = 0; i < mesh.nodes.size(); i++) {
				if (mesh.n_is_removed[i])
					continue;
				os << cnt << " " << 0 << "\n";
				cnt++;
			}


			os << "$EndNodeData\n";
		}

		if(!mesh.t_scalars.empty())
		{
			os << "$ElementData\n";
			os <<"1\n";
		    os << "\"Element-data\"\n";
			os << "1\n"; // num real tags.
		    os << "0.0\n"; // time value.
		    os << "3\n"; // num int tags.
		    os << "0\n"; // the time step
		    os << "1\n"; // 1-component scalar field.
		    os << n_elemens << "\n";

		    cnt = 1;
			for (size_t i = 0; i < mesh.t_scalars.size(); i++) {
				if (mesh.t_is_removed[i])
					continue;

				os << cnt << " " << mesh.t_scalars[i] << "\n";
				++cnt;
			}



			os << "$EndElementData\n";
		}

		os.close();
	}



	void write_msh_DiffusionCurve(MeshData& mesh, const std::string &path)
	{
		std::ofstream os(path+".msh");
		if (os.fail())
			throw std::runtime_error("Unable to open MSH file \"" + path + "\"!");

//split one edge feature piece
typedef std::pair<GEO::index_t, GEO::index_t> Edge;
std::vector<Edge> edges;
for (size_t i = 0; i < mesh.tris.size(); i++) {
	if (mesh.t_is_removed[i])
		continue;

	for (int j = 0; j < 3; j++)
	{
		auto v0 = mesh.tris[i][j];
		auto v1 = mesh.tris[i][(j+1)%3];
		edges.emplace_back(std::make_pair(std::min(v0, v1), std::max(v0, v1)));
	}
}
std::sort(edges.begin(), edges.end());
edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

std::unordered_map<int, std::vector<int>> feature_edges;
for(int i=0;i< edges.size(); i++)
{
	const int v0 = std::get<0>(edges[i]);
	const int v1 = std::get<1>(edges[i]);
	const int edge_tag = triwild::optimization::get_feature_edge_tag(mesh, v0, v1);
	if(edge_tag>=0)
		feature_edges[edge_tag].push_back(i);
}

mesh.v_empty_slot_start = std::find(mesh.v_is_removed.begin(), mesh.v_is_removed.end(), true) - mesh.v_is_removed.begin();
mesh.t_empty_slot_start = std::find(mesh.t_is_removed.begin(), mesh.t_is_removed.end(), true) - mesh.t_is_removed.begin();

std::vector <std::array<int, 2>> new_es;
for(auto &a: feature_edges)
{
	if(a.second.size() <=1)
	{
		std::cout<<"one edge feature, need to split"<<std::endl;
		int v0 = std::get<0>(edges[a.second[0]]);
		int v1 = std::get<1>(edges[a.second[0]]);
		if(!optimization::split_an_edge(mesh, v0, v1, new_es)) //no need b_tree anymore
		{
			std::cout<<"split failed, bad luck!"<<std::endl;
			;//return;
		}
	}
}
//map to local
std::cout<<"map to local"<<std::endl;
Eigen::MatrixXd V;
Eigen::MatrixXi F;

std::unordered_map<int, int> map_v_ids, map_n_ids;
std::unordered_map<int, int> small_to_big_v;
int cnt = 0;
for (size_t i = 0; i < mesh.tri_vertices.size(); i++) {
	if (mesh.v_is_removed[i])
		continue;
	small_to_big_v[cnt] = i;
	map_v_ids[i] = cnt++;

}

for (size_t i = 0; i < mesh.nodes.size(); i++) {
	if (mesh.n_is_removed[i])
		continue;
	map_n_ids[i] = cnt++;
}

V.resize(cnt, 3);
cnt = 0;
for (size_t i = 0; i < mesh.tri_vertices.size(); i++) {
	if (mesh.v_is_removed[i])
		continue;
	V(cnt, 0) = mesh.tri_vertices[i].posf[0];
	V(cnt, 1) = mesh.tri_vertices[i].posf[1];
	V(cnt, 2) = 0;
	cnt++;
}

for (size_t i = 0; i < mesh.nodes.size(); i++) {
	if (mesh.n_is_removed[i])
		continue;
	V(cnt, 0) = mesh.nodes[i][0];
	V(cnt, 1) = mesh.nodes[i][1];
	V(cnt, 2) = 0;
	cnt++;
}

cnt = std::count(mesh.t_is_removed.begin(), mesh.t_is_removed.end(), false);

F.resize(cnt, 3);
std::vector<std::vector<int>> tri_nodes(cnt);

cnt = 0;
for (size_t i = 0; i < mesh.tris.size(); i++) {
	// std::cout<<mesh.tri_nodes[i].size()<<std::endl;
	if (mesh.t_is_removed[i])
		continue;

	for (int j = 0; j < 3; j++)
		F(cnt, j) = map_v_ids[mesh.tris[i][j]];

	if(!mesh.tri_nodes.empty())
	{
		const size_t n_nodes = mesh.tri_nodes[i].size();

		tri_nodes[cnt].reserve(n_nodes);
		for (int j = 0; j < n_nodes; j++)
			tri_nodes[cnt].push_back(map_n_ids[mesh.tri_nodes[i][j]]);
	}
	cnt++;
}
//mesh,
struct HV{ int id; bool boundary = false; Eigen::Vector3d v;std::unordered_set<int> nts; std::unordered_set<int> nvs; std::unordered_set<int> nes;};
struct HE{
	int id;
	int tag = -1;
	std::vector<int> vs;
	std::vector<double> colors;//include blur
	bool boundary = false;
	std::vector<int> nts;

};
struct HF{ int id; std::vector<int> vs; std::vector<int> es;};
struct HM{std::vector<HV> Vs; std::vector<HE> Es; std::vector<HF> Fs;};

auto build_connectivity = [&](HM &m){
	m.Es.clear();
	std::vector<std::tuple<int, int, int, int>> temp;
	temp.reserve(m.Fs.size() * 3);
	for (auto & f: m.Fs) {
		for (uint32_t e = 0; e < 3; ++e) {
			uint32_t v0 = f.vs[e], v1 = f.vs[(e + 1) % 3];
			if (v0 > v1) std::swap(v0, v1);
			temp.push_back(std::make_tuple(v0, v1, f.id, e));
		}
		f.es.resize(3);
	}
	std::sort(temp.begin(), temp.end());
	m.Es.reserve(temp.size() / 2);
	int E_num = -1;
	for (uint32_t i = 0; i < temp.size(); ++i) {
		if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) ||
			std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
			E_num++;
			HE e;
			e.id = E_num;
			e.vs.push_back(std::get<0>(temp[i]));
			e.vs.push_back(std::get<1>(temp[i]));
			e.boundary = true;
			m.Es.push_back(e);
		}
		else if (i != 0 && (std::get<0>(temp[i]) == std::get<0>(temp[i - 1]) &&
			std::get<1>(temp[i]) == std::get<1>(temp[i - 1])))
			m.Es[E_num].boundary = false;

		m.Fs[std::get<2>(temp[i])].es[std::get<3>(temp[i])] = E_num;

		m.Es[E_num].nts.push_back(std::get<2>(temp[i]));
	}
	//v.nvs, v.nfs
	for(auto &v:m.Vs)
	{
		v.boundary = false;
		v.nvs.clear();
		v.nes.clear();
		v.nts.clear();
	}
	for(auto &e:m.Es)
	{
		if(e.boundary)
		{
			m.Vs[e.vs[0]].boundary = true;
			m.Vs[e.vs[1]].boundary = true;
		}
		m.Vs[e.vs[0]].nvs.insert(e.vs[1]);
		m.Vs[e.vs[1]].nvs.insert(e.vs[0]);

		m.Vs[e.vs[0]].nes.insert(e.id);
		m.Vs[e.vs[1]].nes.insert(e.id);
	}

	for(auto &f:m.Fs)
	{
		for(auto vid: f.vs)
			m.Vs[vid].nts.insert(f.id);
		// m.Vs[f.vs[0]].nts.insert(f.id);
		// m.Vs[f.vs[1]].nts.insert(f.id);
		// m.Vs[f.vs[2]].nts.insert(f.id);
	}
};

std::cout<<"build mesh"<<std::endl;
HM hm;
for(int i=0;i<V.rows();i++)
{
	HV hv;
	hv.id = i;
	hv.v = V.row(i);
	hm.Vs.push_back(hv);
}
for(int i=0;i<F.rows();i++)
{
	HF ht;
	ht.id = i;
	ht.vs.push_back(F(i,0));
	ht.vs.push_back(F(i,1));
	ht.vs.push_back(F(i,2));
	for(auto vid: tri_nodes[i])
		ht.vs.push_back(vid);

	hm.Fs.push_back(ht);
}
build_connectivity(hm);

write_OBJ(V.transpose(), F.transpose(), args.input+"_remeshed.obj");

struct feature_arch{int id;
	std::vector<int> es; std::vector<std::pair<double, int>> v_ts;
	std::vector<std::vector<int>> ess;
	std::vector<std::vector<std::pair<double, int>>> vss;
	std::vector<double> leftColor; std::vector<double> rightColor; std::vector<double> blur;
};
auto read_feature =[&](std::string path, std::vector<feature_arch> &arcs){
    json feature_info = json({});
    if (path != "") {
        std::ifstream file(path);
        if (file.is_open())
            file >> feature_info;
        else {
            std::cerr << "unable to open " << path << " file" << std::endl;
            return false;
        }
        file.close();
    } else
        return false;

    for (int i = 0; i < feature_info.size(); i++) {

		feature_arch fa;
		fa.id = (int) feature_info[i]["curve_id"];

		fa.leftColor = feature::json2d2stdvector(feature_info[i]["left_color"]);

		fa.rightColor = feature::json2d2stdvector(feature_info[i]["right_color"]);

		fa.blur = feature::json2d2stdvector(feature_info[i]["blur"]);
		arcs.push_back(fa);
    }
};
auto get_color = [&](const std::vector<feature_arch> &arcs, int curve_id, double t, bool lor, std::vector<double> &color){

	color.clear();
	std::vector<double> color_chain;
	if(lor) color_chain = arcs[curve_id].leftColor;
	else color_chain = arcs[curve_id].rightColor;
	for(int i=0;i<color_chain.size();)
	{
		if(color_chain[i+3]<= t && color_chain[i+ 4+ 3] >= t)
		{
			std::vector<double> c_l (color_chain.begin()+i, color_chain.begin()+i+4);
			std::vector<double> c_h (color_chain.begin()+i+4, color_chain.begin()+i+8);

			double t_ = 0;
			if(color_chain[i+3] != color_chain[i+ 4+ 3])
				t_ = (c_h[3] - t)/(c_h[3] - c_l[3]);
			color.push_back(c_h[0] - t_* (c_h[0]-c_l[0]));
			color.push_back(c_h[1] - t_* (c_h[1]-c_l[1]));
			color.push_back(c_h[2] - t_* (c_h[2]-c_l[2]));
			break;
		}
		i+=4;
	}

	for(int i=0;i<arcs[curve_id].blur.size();)
	{
		if(arcs[curve_id].blur[i+1]<= t && arcs[curve_id].blur[i+ 2+ 1] >= t)
		{
			std::vector<double> b_l (arcs[curve_id].blur.begin()+i, arcs[curve_id].blur.begin()+i+2);
			std::vector<double> b_h (arcs[curve_id].blur.begin()+i+2, arcs[curve_id].blur.begin()+i+4);

			double t_ = 0;
			if(arcs[curve_id].blur[i+1] != arcs[curve_id].blur[i+ 2+ 1])
				t_ = (b_h[1] - t)/(b_h[1] - b_l[1]);
			color.push_back(b_h[0] - t_* (b_h[0] - b_l[0]));
			break;
		}
		i+=2;
	}
};

std::cout<<"read Json feature file"<<std::endl;
std::vector<feature_arch> arches;
read_feature(args.feature_input, arches);

std::cout<<"build feature arches"<<std::endl;
std::vector<int> E_label(hm.Es.size(), -1);
for(auto &e: hm.Es)
{
	const int v0 = small_to_big_v[e.vs[0]];
	const int v1 = small_to_big_v[e.vs[1]];
	const int edge_tag = triwild::optimization::get_feature_edge_tag(mesh, v0, v1);
	E_label[e.id] = edge_tag;
	if(edge_tag>=0){
		arches[feature::features[edge_tag]->curve_id].es.push_back(e.id);

		for(int i =0;i<2;i++)
		{
			double t = -1;
			auto v = small_to_big_v[e.vs[i]];
			for(auto ft: mesh.tri_vertices[v].feature_infos)
			{
				if(ft[0] == edge_tag) {t = ft[1];break;}
			}
			if(t==-1) {
				std::cout<<"error"<<std::endl;
				std::cin.get();
			}

			if(feature::features[edge_tag]->curve_id == 130)
			{
				std::cout<<"t: "<<t<<std::endl;
				std::cout<<"evid: "<<e.vs[i]<<std::endl;
			}

			arches[feature::features[edge_tag]->curve_id].v_ts.push_back(std::make_pair(t, e.vs[i]));
		}
	}
}
//into independent triangles.
std::cout<<"into independent triangles"<<std::endl;
std::vector<std::vector<int>> V2Vs(hm.Vs.size());
std::vector<int> Vs2V;
HM hm_;
for(int i=0;i<F.rows();i++)
{
	HF ht;
	ht.id = i;
	ht.vs.push_back(F(i,0));
	ht.vs.push_back(F(i,1));
	ht.vs.push_back(F(i,2));
	for(auto vid: tri_nodes[i])
		ht.vs.push_back(vid);

	hm_.Fs.push_back(ht);
}
for(auto &v: hm.Vs)
{
	HV v_;
	v_.id = hm_.Vs.size();
	v_.v = v.v;
	hm_.Vs.push_back(v_);
	V2Vs[v.id].push_back(v_.id);
	Vs2V.push_back(v.id);
}
for(auto &v: hm.Vs)
{
	assert(v.nts.size());
	int i = -1;
	for ( auto it = v.nts.begin(); it != v.nts.end(); ++it)
	{
		i++;
		if(i==0) continue;

		HV v_;
		v_.id = hm_.Vs.size();
		v_.v = v.v;
		for(auto &tv: hm_.Fs[*it].vs)
			if(tv == v.id) tv = v_.id;
		hm_.Vs.push_back(v_);
		V2Vs[v.id].push_back(v_.id);
		Vs2V.push_back(v.id);
	}
}
build_connectivity(hm_);
V.resize(hm_.Vs.size(), 3);
for(auto &v: hm_.Vs)
{
	V(v.id, 0) = v.v[0];
	V(v.id, 1) = v.v[1];
	V(v.id, 2) = v.v[2];
}
F.resize(hm_.Fs.size(), 3);
for(auto &f: hm_.Fs)
{
	F(f.id, 0) = f.vs[0];
	F(f.id, 1) = f.vs[1];
	F(f.id, 2) = f.vs[2];
}
write_OBJ(V.transpose(), F.transpose(), args.input+"_separated.obj");
//merge back edges
std::cout<<"merge back edges"<<std::endl;
std::vector<std::vector<int>> E2Es(hm.Es.size());
std::vector<int> Es2E;
for(auto &e: hm_.Es)
{
	auto v0 = Vs2V[e.vs[0]], v1 = Vs2V[e.vs[1]];

	const auto &nes0 = hm.Vs[v0].nes;
	const auto &nes1 = hm.Vs[v1].nes;
	std::vector<int> sharedes;
	unordered_set_intersection_own(nes0, nes1, sharedes);
	if(!sharedes.size())
	{
		std::cout<<"ev0, ev1 "<<e.vs[0]<<" "<<e.vs[1]<<std::endl;
		std::cout<<"v0, v1 "<<v0<<" "<<v1<<std::endl;
	}
	E2Es[sharedes[0]].push_back(e.id);
	Es2E.push_back(sharedes[0]);
}

std::vector<bool> E_tag(hm_.Es.size(), true);
std::vector<std::unordered_set<int>> Vs_2Vs_reverse(hm_.Vs.size());
std::vector<int> Vs2Vs_(hm_.Vs.size());
for(auto &v: hm_.Vs)
{
	Vs2Vs_[v.id] = v.id;
	Vs_2Vs_reverse[v.id].insert(v.id);
}
for(auto &e: hm_.Es)
{
	auto e_map = Es2E[e.id];
	if(E_tag[e.id] && E_label[e_map] == -1 && E2Es[e_map].size() == 2)
	{
		auto e2 = E2Es[e_map][0] == e.id? E2Es[e_map][1]:E2Es[e_map][0];
		E_tag[e2] = false;

		auto v0 = Vs2Vs_[e.vs[0]];
		auto v1 = Vs2Vs_[e.vs[1]];

		auto v2 = Vs2Vs_[hm_.Es[e2].vs[0]];
		auto v3 = Vs2Vs_[hm_.Es[e2].vs[1]];
		//find corresponding pairs
		auto v0_map = Vs2V[v0];
		auto v1_map = Vs2V[v1];

		auto v2_map = Vs2V[v2];
		auto v3_map = Vs2V[v3];

		if((v0_map != v2_map && v0_map !=v3_map)||
		(v1_map != v2_map && v1_map !=v3_map))
		{
			std::cout<<"error "<<std::endl;
		}

		if(v0_map == v3_map && v1_map ==v2_map)
		{
			std::swap(v2, v3);
		}

		if(v0 != v2)
		{
			Vs_2Vs_reverse[v0].insert(Vs_2Vs_reverse[v2].begin(), Vs_2Vs_reverse[v2].end());

			for(auto &vid: Vs_2Vs_reverse[v2])
				Vs2Vs_[vid] = v0;
			Vs_2Vs_reverse[v2].clear();
		}
		if(v1 != v3)
		{
			Vs_2Vs_reverse[v1].insert(Vs_2Vs_reverse[v3].begin(), Vs_2Vs_reverse[v3].end());

			for(auto &vid: Vs_2Vs_reverse[v3])
				Vs2Vs_[vid] = v1;
			Vs_2Vs_reverse[v3].clear();
		}
	}
}
//reindexing
std::cout<<"reindexing"<<std::endl;
HM hm_2;
for(auto &f: hm_.Fs)
{
	HF ht;
	ht.id = f.id;
	for(auto &vid : f.vs) ht.vs.push_back(Vs2Vs_[vid]);
	hm_2.Fs.push_back(ht);
}
std::vector<int> V22V(hm_.Vs.size());
std::vector<int> V22V_reverse;
for(auto &v: hm_.Vs)
{
	if(Vs_2Vs_reverse[v.id].size())
	{
		HV v_;
		v_.id = hm_2.Vs.size();
		v_.v = v.v;
		hm_2.Vs.push_back(v_);
		V22V[v.id] = v_.id;
		V22V_reverse.push_back(v.id);
	}
}
for(auto &f: hm_2.Fs)
{
	for(auto &vid : f.vs) vid = V22V[vid];
}
build_connectivity(hm_2);

V.resize(hm_2.Vs.size(), 3);
for(auto &v: hm_2.Vs)
{
	V(v.id, 0) = v.v[0];
	V(v.id, 1) = v.v[1];
	V(v.id, 2) = v.v[2];
}
F.resize(hm_2.Fs.size(), 3);
for(auto &f: hm_2.Fs)
{
	F(f.id, 0) = f.vs[0];
	F(f.id, 1) = f.vs[1];
	F(f.id, 2) = f.vs[2];
}
write_OBJ(V.transpose(), F.transpose(), args.input+"_merged.obj");
//color edges
std::cout<<"color edges"<<std::endl;
for(auto &e: hm_2.Es)
{
	auto ev0 = e.vs[0];
	auto ev1 = e.vs[1];

	auto v0 = V22V_reverse[e.vs[0]];
	auto v1 = V22V_reverse[e.vs[1]];
	auto v0_map = Vs2V[v0];
	auto v1_map = Vs2V[v1];

	//curve id
	const auto &nes0 = hm.Vs[v0_map].nes;
	const auto &nes1 = hm.Vs[v1_map].nes;
	std::vector<int> sharedes;
	unordered_set_intersection_own(nes0, nes1, sharedes);
	assert(sharedes.size());
	auto edge_tag = E_label[sharedes[0]];

	if(!e.boundary) {
		if(edge_tag != -1)
		{
			std::cout<<"one edge feature not cleaned"<<std::endl;
			std::cout<<"ev0 ev1: "<<ev0 << " "<<ev1<<std::endl;
			std::cout<<"v0_map v1_map: "<<v0_map << " "<<v1_map<<std::endl;
		}
		continue;
	}
	if( edge_tag != -1)
	{
		//e color
		auto a = arches[feature::features[edge_tag]->curve_id];
		e.tag = a.id;
		double t0 = -1, t1= -1;
		for(auto &v: a.v_ts)
		{
			if(v0_map == std::get<1> (v))
				t0 = std::get<0> (v);
			if(v1_map == std::get<1> (v))
				t1 = std::get<0> (v);
		}
		assert(t0 !=-1);
		assert(t1 !=-1);

		if(t0 > t1){
			std::swap(v0_map, v1_map);
			std::swap(v0, v1);
			std::swap(t0, t1);
			std::swap(ev0, ev1);
		}
		int tid = e.nts[0];

		int v0id = std::find(hm_2.Fs[tid].vs.begin(),hm_2.Fs[tid].vs.end(),ev0) -hm_2.Fs[tid].vs.begin();
		int v1id = std::find(hm_2.Fs[tid].vs.begin(),hm_2.Fs[tid].vs.end(),ev1) -hm_2.Fs[tid].vs.begin();

		bool left = true;
		if((v0id+1)%3 != v1id)
			left = false;

		if(std::find(hm_2.Fs[tid].vs.begin(),hm_2.Fs[tid].vs.end(),14039) != hm_2.Fs[tid].vs.end())
		{
			std::cout<<"t0 t1: "<<t0 << " "<<t1<<std::endl;
			std::cout<<"ev0 ev1: "<<ev0 << " "<<ev1<<std::endl;
			std::cout<<"v0_map v1_map: "<<v0_map << " "<<v1_map<<std::endl;
			std::cout<<"a.id: "<<a.id <<std::endl;

			std::cout<<"left: "<<left <<std::endl;

			for(auto &tvid: hm_2.Fs[tid].vs)
				std::cout<<"vid: "<<tvid <<std::endl;
		}
		if(std::find(hm_2.Fs[tid].vs.begin(),hm_2.Fs[tid].vs.end(),14095) != hm_2.Fs[tid].vs.end())
		{
			std::cout<<"t0 t1: "<<t0 << " "<<t1<<std::endl;
			std::cout<<"ev0 ev1: "<<ev0 << " "<<ev1<<std::endl;
			std::cout<<"v0_map v1_map: "<<v0_map << " "<<v1_map<<std::endl;
			std::cout<<"a.id: "<<a.id <<std::endl;

			std::cout<<"left: "<<left <<std::endl;
			for(auto &tvid: hm_2.Fs[tid].vs)
				std::cout<<"vid: "<<tvid <<std::endl;
		}

		if(std::find(hm_2.Fs[tid].vs.begin(),hm_2.Fs[tid].vs.end(),14109) != hm_2.Fs[tid].vs.end())
		{
			std::cout<<"t0 t1: "<<t0 << " "<<t1<<std::endl;
			std::cout<<"ev0 ev1: "<<ev0 << " "<<ev1<<std::endl;
			std::cout<<"v0_map v1_map: "<<v0_map << " "<<v1_map<<std::endl;
			std::cout<<"a.id: "<<a.id <<std::endl;

			std::cout<<"left: "<<left <<std::endl;
			for(auto &tvid: hm_2.Fs[tid].vs)
				std::cout<<"vid: "<<tvid <<std::endl;
		}
		if(std::find(hm_2.Fs[tid].vs.begin(),hm_2.Fs[tid].vs.end(),2090) != hm_2.Fs[tid].vs.end())
		{
			std::cout<<"t0 t1: "<<t0 << " "<<t1<<std::endl;
			std::cout<<"ev0 ev1: "<<ev0 << " "<<ev1<<std::endl;
			std::cout<<"v0_map v1_map: "<<v0_map << " "<<v1_map<<std::endl;
			std::cout<<"a.id: "<<a.id <<std::endl;

			std::cout<<"left: "<<left <<std::endl;
			for(auto &tvid: hm_2.Fs[tid].vs)
				std::cout<<"vid: "<<tvid <<std::endl;
		}


		std::vector<double> color;
		get_color(arches, a.id, t0, left, color);
		e.colors = color;
		get_color(arches, a.id, t1, left, color);
		e.colors.insert(e.colors.end(), color.begin(), color.end());
	}
}

hm = hm_2;

V.resize(3, hm.Vs.size());
F.resize(3, hm.Fs.size());
for(auto &v: hm.Vs)
{
	V(0, v.id) = v.v[0];
	V(1, v.id) = v.v[1];
	V(2, v.id) = v.v[2];
}
for(auto &f: hm.Fs)
{
	F(0, f.id) = f.vs[0];
	F(1, f.id) = f.vs[1];
	F(2, f.id) = f.vs[2];
}
write_OBJ(V, F, args.output+"_cut.obj");

GEO::Mesh M;
// Setup vertices
M.vertices.create_vertices((int) hm.Vs.size());
for (int i = 0; i < (int) M.vertices.nb(); ++i) {
	GEO::vec3 &p = M.vertices.point(i);
	p[0] = hm.Vs[i].v[0];
	p[1] = hm.Vs[i].v[1];
	p[2] = hm.Vs[i].v[2];
}
// Setup faces
M.facets.create_triangles((int) hm.Fs.size());

for (int c = 0; c < (int) M.facets.nb(); ++c) {
	for (int lv = 0; lv < 3; ++lv) {
		M.facets.set_vertex(c, lv, hm.Fs[c].vs[lv]);
	}
}
M.facets.connect();
M.cells.connect();

// Compute a list of all the edges, and store edge index as a corner attribute
std::vector<std::pair<Edge, GEO::index_t>> e2c; // edge to corner id
for (GEO::index_t f = 0; f < M.facets.nb(); ++f) {
	for (GEO::index_t c = M.facets.corners_begin(f); c < M.facets.corners_end(f); ++c) {
		GEO::index_t v = M.facet_corners.vertex(c);
		GEO::index_t c2 = M.facets.next_corner_around_facet(f, c);
		GEO::index_t v2 = M.facet_corners.vertex(c2);
		e2c.emplace_back(std::make_pair(std::min(v, v2), std::max(v, v2)), c);
	}
}
std::sort(e2c.begin(), e2c.end());

// Assign unique id to edges
// GEO::Attribute<GEO::index_t> c2e(M.facet_corners.attributes(), "edge_id");
M.edges.clear();
Edge prev_e(-1, -1);
GEO::index_t current_id = -1;
std::vector<bool> boundary_edges;
for (const auto &kv : e2c) {
	Edge e = kv.first;
	GEO::index_t c = kv.second;
	if (e != prev_e) {
		M.edges.create_edge(e.first, e.second);
		boundary_edges.push_back(true);
		++current_id;
		prev_e = e;
	} else {
		boundary_edges.back() = false;
	}
}

std::ofstream tags(path + ".txt");
std::ofstream tag_r(path + "_b.txt");
std::ofstream tag_g(path + "_g.txt");
std::ofstream tag_b(path + "_r.txt");

if(!tags.fail())
{
	for(int e = 0; e < M.edges.nb(); ++e)
	{
		const int v0 = M.edges.vertex(e, 0);
		const int v1 = M.edges.vertex(e, 1);


		const auto &nes0 = hm.Vs[v0].nes;
		const auto &nes1 = hm.Vs[v1].nes;
		std::vector<int> sharedes;
		unordered_set_intersection_own(nes0, nes1, sharedes);
		if(!sharedes.size())
		{
			std::cout<<v0<<" "<<v1<<std::endl;
			for(auto &vid: hm.Vs[v0].nvs)
				std::cout<<vid<<std::endl;
			std::cout<<"v1 nvs: ";
			for(auto &vid: hm.Vs[v1].nvs)
				std::cout<<vid<<std::endl;
			std::cout<<"v0 nes: ";
			for(auto &eid: hm.Vs[v0].nes)
				std::cout<<eid<<std::endl;
			std::cout<<"v1 nes: ";
			for(auto &eid: hm.Vs[v1].nes)
				std::cout<<eid<<std::endl;

		}
		assert(sharedes.size());

		tags<<(hm.Es[sharedes[0]].tag+1)<<"\n";

		auto color = hm.Es[sharedes[0]].colors;
		if(color.size())
		{
			tag_r<<e<<" "<<"1 "<<color[0]/255.0<<" "<<color[4]/255.0<<"\n";
			tag_g<<e<<" "<<"1 "<<color[1]/255.0<<" "<<color[5]/255.0<<"\n";
			tag_b<<e<<" "<<"1 "<<color[2]/255.0<<" "<<color[6]/255.0<<"\n";
		}
		else{
			continue;
			// tag_r<<e<<" "<<"1 "<<0/255.0<<" "<<0/255.0<<"\n";
			// tag_g<<e<<" "<<"1 "<<0/255.0<<" "<<0/255.0<<"\n";
			// tag_b<<e<<" "<<"1 "<<0/255.0<<" "<<0/255.0<<"\n";
		}
	}
}
tags.close();
tag_r.close();
tag_g.close();
tag_b.close();


	const auto n_nodes = hm.Vs.size();
	const auto n_elemens = hm.Fs.size();

	os << std::setprecision(16)<< "$MeshFormat\n";
	os << "2.2 0 8\n";
	os << "$EndMeshFormat\n";

	os << "$Nodes\n";
	os << n_nodes <<"\n";
	for (size_t i = 0; i<n_nodes; ++i)
	{
		os << (i+1);
		for (int j = 0; j < 3; j++)
			os << " " << hm.Vs[i].v[j];
		os << "\n";
	}
	os << "$EndNodes\n";

	os << "$Elements\n";
	os << n_elemens <<"\n";
	for (size_t i = 0; i < n_elemens; ++i) {
		int loc_nodes = hm.Fs[i].vs.size() - 3;
		// el_type:
		//		2: 3-node triangle
		//		9: 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges)
		//		21: 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
		//		23: 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)

		//default: 2: 3-node triangle
		int el_type = 2;
		switch(loc_nodes)
		{
			//2: 3-node triangle
			case 0: el_type = 2; break;

			//9: 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges)
			case 3: el_type = 9; break;

			//21: 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
			case 7: el_type = 21; break;

			//23: 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)
			case 12: el_type = 23; break;

			default: std::cerr<<"unspported number of nodes "<<hm.Fs[i].vs.size()<<std::endl; assert(false);
		}

		os << (i+1) << " " << el_type << " 0"; //el_id el_type n_tags

		//first 3 are always the vertices
		for (int j = 0; j < 3; j++)
			os << " " << hm.Fs[i].vs[j] + 1;

		// follow for node order http://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_65.php#SEC65
		for (size_t j = 0; j < loc_nodes; j++)
		{
			os << " " << hm.Fs[i].vs[j+3] + 1;
		}
		os<<"\n";

	}
	os << "$EndElements\n";
	os.close();
	}


	void export_eps(MeshData& mesh,
    	const double line_width, const std::string &col, const double point_size, const std::string &point_col,
    	const double feature_line_width, const std::string &feature_col, const double feature_point_size, const std::string &feature_point_col,
    	const double secondary_feature_line_width, const std::string &secondary_feature_col, const double secondary_feature_point_size, const std::string &secondary_feature_point_col,
    	const bool draw_points, const std::string &path, const double t)
	{
		Eigen::MatrixXd V;
		Eigen::MatrixXi F;

		std::unordered_map<int, int> map_v_ids, map_n_ids;
		std::unordered_map<int, int> small_to_big_v;
		int cnt = 0;
		for (size_t i = 0; i < mesh.tri_vertices.size(); i++) {
			if (mesh.v_is_removed[i])
				continue;
			small_to_big_v[cnt] = i;
			map_v_ids[i] = cnt++;

		}

		// for (size_t i = 0; i < mesh.nodes.size(); i++) {
		// 	if (mesh.n_is_removed[i])
		// 		continue;
		// 	map_n_ids[i] = cnt++;
		// }

		V.resize(cnt, 3);
		cnt = 0;
		for (size_t i = 0; i < mesh.tri_vertices.size(); i++) {
			if (mesh.v_is_removed[i])
				continue;
			V(cnt, 0) = mesh.tri_vertices[i].posf[0];
			V(cnt, 1) = mesh.tri_vertices[i].posf[1];
			V(cnt, 2) = 0;
			cnt++;
		}

		// for (size_t i = 0; i < mesh.nodes.size(); i++) {
		// 	if (mesh.n_is_removed[i])
		// 		continue;
		// 	V(cnt, 0) = mesh.nodes[i][0];
		// 	V(cnt, 1) = mesh.nodes[i][1];
		// 	V(cnt, 2) = 0;
		// 	cnt++;
		// }

		cnt = std::count(mesh.t_is_removed.begin(), mesh.t_is_removed.end(), false);

		F.resize(cnt, 3);
		// std::vector<std::vector<int>> tri_nodes(cnt);

		cnt = 0;
		for (size_t i = 0; i < mesh.tris.size(); i++) {
			// std::cout<<mesh.tri_nodes[i].size()<<std::endl;
			if (mesh.t_is_removed[i])
				continue;

			for (int j = 0; j < 3; j++)
				F(cnt, j) = map_v_ids[mesh.tris[i][j]];

			// if(!mesh.tri_nodes.empty())
			// {
			// 	const size_t n_nodes = mesh.tri_nodes[i].size();

			// 	tri_nodes[cnt].reserve(n_nodes);
			// 	for (int j = 0; j < n_nodes; j++)
			// 		tri_nodes[cnt].push_back(map_n_ids[mesh.tri_nodes[i][j]]);
			// }
			cnt++;
		}

		typedef std::pair<GEO::index_t, GEO::index_t> Edge;


		GEO::Mesh M;
		// Setup vertices
		M.vertices.create_vertices((int) V.rows());
		for (int i = 0; i < (int) M.vertices.nb(); ++i) {
			GEO::vec3 &p = M.vertices.point(i);
			p[0] = V(i, 0);
			p[1] = V(i, 1);
			p[2] = V(i, 2);
		}

		// Setup faces
		assert(F.cols() == 3);
		M.facets.create_triangles((int) F.rows());

		for (int c = 0; c < (int) M.facets.nb(); ++c) {
			for (int lv = 0; lv < F.cols(); ++lv) {
				M.facets.set_vertex(c, lv, F(c, lv));
			}
		}

		M.facets.connect();
		M.cells.connect();

		// Compute a list of all the edges, and store edge index as a corner attribute
		std::vector<std::pair<Edge, GEO::index_t>> e2c; // edge to corner id
		for (GEO::index_t f = 0; f < M.facets.nb(); ++f) {
			for (GEO::index_t c = M.facets.corners_begin(f); c < M.facets.corners_end(f); ++c) {
				GEO::index_t v = M.facet_corners.vertex(c);
				GEO::index_t c2 = M.facets.next_corner_around_facet(f, c);
				GEO::index_t v2 = M.facet_corners.vertex(c2);
				e2c.emplace_back(std::make_pair(std::min(v, v2), std::max(v, v2)), c);
			}
		}
		std::sort(e2c.begin(), e2c.end());

		// Assign unique id to edges
		// GEO::Attribute<GEO::index_t> c2e(M.facet_corners.attributes(), "edge_id");
		M.edges.clear();
		Edge prev_e(-1, -1);
		GEO::index_t current_id = -1;
		std::vector<bool> boundary_edges;
		for (const auto &kv : e2c) {
			Edge e = kv.first;
			GEO::index_t c = kv.second;
			if (e != prev_e) {
				M.edges.create_edge(e.first, e.second);
				boundary_edges.push_back(true);
				++current_id;
				prev_e = e;
			} else {
				boundary_edges.back() = false;
			}
			// c2e[c] = current_id;
		}

		std::ofstream os(path);
		if(os.fail())
			return;

		os << std::setprecision(16);
		Eigen::Matrix<double, 4, 4> L2B; L2B <<
			1, 0, 0, 0,
  			-5./6., 3,  -3./2, 1./3.,
   			1./3.,  -3./2., 3, -5./6.,
            0, 0, 0, 1;

		os <<"%!PS-Adobe-3.0 EPSF-3.0\n";
		os << "%%BoundingBox: "<< (args.box_min(0)-2*args.target_edge_len) << " " << (args.box_min(1)-args.target_edge_len) << " " << (args.box_max(0)+2*args.target_edge_len) <<" "<< (args.box_max(1)+args.target_edge_len) <<"\n\n";
		os << "%%Pages: 1\n";
		os << "%%Page: 1 1\n";
		os << "/show-ctr {\ndup stringwidth pop\n -2 div 0\n rmoveto show\n} def\n\n 2 setlinejoin\n\n";

		std::vector<bool> primary_vertices(mesh.tri_vertices.size(), false);

		os<<col<<" setrgbcolor\n";
		os<<line_width<<" setlinewidth\n\n";
		for(int e = 0; e < M.edges.nb(); ++e)
		{
			const int v0 = M.edges.vertex(e, 0);
			const int v1 = M.edges.vertex(e, 1);
			const int edge_tag = triwild::optimization::get_feature_edge_tag(mesh, small_to_big_v[v0], small_to_big_v[v1]);
			if(edge_tag >= 0)
				continue;

			const int secondary_edge_tag = triwild::optimization::get_secondary_feature_edge_tag(mesh, small_to_big_v[v0], small_to_big_v[v1]);
			if(secondary_edge_tag >= 0)
				continue;

			const double *p0 = M.vertices.point_ptr(v0);
			const double *p1 = M.vertices.point_ptr(v1);

			primary_vertices[small_to_big_v[v0]] = true;
			primary_vertices[small_to_big_v[v1]] = true;

			os<<p0[0]<<" "<<p0[1]<<" moveto\n";
			os<<p1[0]<<" "<<p1[1]<<" lineto\n";
		}
		os<<"stroke\n\n\n";

		if(draw_points)
		{
			os<<point_col<<" setrgbcolor\n";
			for (int i = 0; i < V.rows(); ++i)
			{
				if(!mesh.tri_vertices[small_to_big_v[i]].feature_infos.empty())
					continue;

				if(!primary_vertices[small_to_big_v[i]])
					continue;

				//TODO skip secondary vertices
				os << V(i, 0) << " " << V(i, 1)  << " moveto\n";
				os << V(i, 0) << " " << V(i, 1) << " " << point_size << " 0 360 arc\n";
			}
			os<<"fill\n\n";
		}




		os<<secondary_feature_col<<" setrgbcolor\n";
		os<<secondary_feature_line_width<<" setlinewidth\n\n";
		for(int e = 0; e < M.edges.nb(); ++e)
		{
			const int v0 = M.edges.vertex(e, 0);
			const int v1 = M.edges.vertex(e, 1);
			const int edge_tag = triwild::optimization::get_feature_edge_tag(mesh, small_to_big_v[v0], small_to_big_v[v1]);
			if(edge_tag >= 0)
				continue;

			const int secondary_edge_tag = triwild::optimization::get_secondary_feature_edge_tag(mesh, small_to_big_v[v0], small_to_big_v[v1]);
			if(secondary_edge_tag < 0)
				continue;

			const double *p0 = M.vertices.point_ptr(v0);
			const double *p1 = M.vertices.point_ptr(v1);

			os<<p0[0]<<" "<<p0[1]<<" moveto\n";
			os<<p1[0]<<" "<<p1[1]<<" lineto\n";
		}
		os<<"stroke\n\n\n";

		if(draw_points)
		{
			os<<secondary_feature_point_col<<" setrgbcolor\n";
			for (int i = 0; i < V.rows(); ++i)
			{
				if(!mesh.tri_vertices[small_to_big_v[i]].feature_infos.empty())
					continue;

				if(primary_vertices[small_to_big_v[i]])
					continue;

				//TODO skip secondary vertices
				os << V(i, 0) << " " << V(i, 1)  << " moveto\n";
				os << V(i, 0) << " " << V(i, 1) << " " << secondary_feature_point_size << " 0 360 arc\n";
			}
			os<<"fill\n\n";
		}



		// os<<"["<<line_width*2<<"] 0 setdash\n";
		// os<<line_width/2<<" setlinewidth\n\n";
		// for(int e = 0; e < M.edges.nb(); ++e)
		// {
		// 	const int v0 = M.edges.vertex(e, 0);
		// 	const int v1 = M.edges.vertex(e, 1);
		// 	const int edge_tag = triwild::optimization::get_feature_edge_tag(mesh, small_to_big_v[v0], small_to_big_v[v1]);
		// 	if(edge_tag < 0)
		// 		continue;

		// 	const double *p0 = M.vertices.point_ptr(v0);
		// 	const double *p1 = M.vertices.point_ptr(v1);

		// 	os<<p0[0]<<" "<<p0[1]<<" moveto\n";
		// 	os<<p1[0]<<" "<<p1[1]<<" lineto\n";
		// }
		// os<<"stroke\n\n\n";

		os<<"[] 0 setdash\n";
		os<<feature_col<<" setrgbcolor\n";
		os<<feature_line_width<<" setlinewidth\n\n";
		for(int e = 0; e < M.edges.nb(); ++e)
		{
			const int v0 = M.edges.vertex(e, 0);
			const int v1 = M.edges.vertex(e, 1);
			const int v0b = small_to_big_v[v0];
			const int v1b = small_to_big_v[v1];

			const int edge_tag = triwild::optimization::get_feature_edge_tag(mesh, v0b, v1b);
			if(edge_tag < 0)
				continue;

			if(!feature::features.empty() && (feature::features[edge_tag]->type == "Line" || mesh.tri_nodes.empty()))
			{
				os << mesh.tri_vertices[v0b].posf[0] << " " << mesh.tri_vertices[v0b].posf[1] << " moveto\n";
				os << mesh.tri_vertices[v1b].posf[0] << " " << mesh.tri_vertices[v1b].posf[1] << " lineto\n";

				continue;
			}

			const int tri_id = optimization::set_intersection(mesh.tri_vertices[v0b].conn_tris, mesh.tri_vertices[v1b].conn_tris)[0];

			if(mesh.tri_nodes[tri_id].empty())
			{
				os << mesh.tri_vertices[v0b].posf[0] << " " << mesh.tri_vertices[v0b].posf[1] << " moveto\n";
				os << mesh.tri_vertices[v1b].posf[0] << " " << mesh.tri_vertices[v1b].posf[1] << " lineto\n";

				continue;
			}


			int lv0 = -1;
			for(int i = 0; i < 3; ++i)
			{
				if(mesh.tris[tri_id][i] == v0b){
					lv0=i;
					break;
				}
			}

			assert(lv0 >= 0);

			Eigen::Matrix<double, 4, 2> nodes;
			nodes.row(0) << mesh.tri_vertices[v0b].posf[0], mesh.tri_vertices[v0b].posf[1];


			if(mesh.tris[tri_id][(lv0+1)%3] == v1b)
			{
				//forward
				nodes.row(1) << mesh.nodes[mesh.tri_nodes[tri_id][2*lv0]][0], mesh.nodes[mesh.tri_nodes[tri_id][2*lv0]][1];
				nodes.row(2) << mesh.nodes[mesh.tri_nodes[tri_id][2*lv0+1]][0], mesh.nodes[mesh.tri_nodes[tri_id][2*lv0+1]][1];
			}
			else
			{
				const int lv1 = (lv0+2)%3;
				assert(mesh.tris[tri_id][lv1] == v1b);
				//backward
				nodes.row(1) << mesh.nodes[mesh.tri_nodes[tri_id][2*lv1+1]][0], mesh.nodes[mesh.tri_nodes[tri_id][2*lv1+1]][1];
				nodes.row(2) << mesh.nodes[mesh.tri_nodes[tri_id][2*lv1]][0], mesh.nodes[mesh.tri_nodes[tri_id][2*lv1]][1];
			}

			nodes.row(3) << mesh.tri_vertices[v1b].posf[0], mesh.tri_vertices[v1b].posf[1];

			if(t != 1)
			{
				nodes.row(1) = (1-t)*(1./3*nodes.row(0)+2./3.*nodes.row(3)) + t*nodes.row(1);
				nodes.row(2) = (1-t)*(2./3*nodes.row(0)+1./3.*nodes.row(3)) + t*nodes.row(2);
			}

			assert(lv0 >= 0);

#ifndef NDEBUG
			const double *p0 = M.vertices.point_ptr(v0);
			assert(fabs(p0[0] - nodes(0, 0)) < 1e-10);
			assert(fabs(p0[1] - nodes(0, 1)) < 1e-10);
			const double *p1 = M.vertices.point_ptr(v1);

			assert(fabs(p1[0] - nodes(3, 0)) < 1e-10);
			assert(fabs(p1[1] - nodes(3, 1)) < 1e-10);
#endif
			const Eigen::Matrix<double, 4, 2> bezier = L2B*nodes;

			os<<bezier(0,0)<<" "<<bezier(0,1)<<" moveto\n";
			os<<bezier(1,0)<<" "<<bezier(1,1)<<" ";
			os<<bezier(2,0)<<" "<<bezier(2,1)<<" ";
			os<<bezier(3,0)<<" "<<bezier(3,1)<<" curveto\n";
		}
		os<<"stroke\n\n\n";

		if(draw_points)
		{
			os<<feature_point_col<<" setrgbcolor\n";
			for (int i = 0; i < V.rows(); ++i)
			{
				if(mesh.tri_vertices[small_to_big_v[i]].feature_infos.empty())
					continue;
				os << V(i, 0) << " " << V(i, 1)  << " moveto\n";
				os << V(i, 0) << " " << V(i, 1) << " " << feature_point_size << " 0 360 arc\n";
			}
			os<<"fill\n\n";
		}

		os.close();
	}


	void load_msh(const std::string &path, MeshData& mesh)
	{
		std::ifstream infile(path.c_str());

		std::string line;

		int phase = -1;
		int line_number = -1;
		bool size_read = false;

		int n_triangles = 0;

		std::vector<std::vector<int>> all_elements;

		Eigen::MatrixXd vertices;

		while (std::getline(infile, line))
		{
			++line_number;

			if(line.empty())
				continue;

			if(line[0] == '$')
			{
				if(line.substr(1,3) == "End")
					phase = -1;
				else
				{
					const auto header = line.substr(1);

					if(header.find("MeshFormat") == 0)
						phase = 0;
					else if(header.find("Nodes") == 0)
						phase = 1;
					else if(header.find("Elements") == 0)
						phase = 2;
					else
					{
						phase = -1;
					}
				}

				size_read = false;

				continue;
			}


			if(phase == -1)
				continue;

			std::istringstream iss(line);
			//header
			if(phase == 0)
			{
				double version_number;
				int file_type;
				int data_size;

				iss >> version_number >> file_type >> data_size;

				assert(version_number == 2.2);
				assert(file_type == 0);
				assert(data_size == 8);
			}
			//coordiantes
			else if(phase == 1)
			{
				if(!size_read)
				{
					int n_vertices;
					iss >> n_vertices;
					vertices.resize(n_vertices, 2);
					size_read = true;
				}
				else
				{
					int node_number;
					double x_coord, y_coord, z_coord;

					iss >> node_number >> x_coord >> y_coord >> z_coord;
					//node_numbers starts with 1
					vertices.row(node_number-1) << x_coord, y_coord; //, z_coord;
				}
			}
			//elements
			else if(phase == 2)
			{
				if(!size_read)
				{
					int number_of_elements;
					iss >> number_of_elements;
					all_elements.resize(number_of_elements);
					size_read = true;
				}
				else
				{
					int elm_number, elm_type, number_of_tags;

					iss >> elm_number >> elm_type >> number_of_tags;

					//point
					if(elm_type == 15)
						continue;

					//edges
					if(elm_type == 1 || elm_type == 8 || elm_type == 26 || elm_type == 27 || elm_type == 28)
						continue;


					if(elm_type == 2 || elm_type == 9 || elm_type == 21 || elm_type == 23)
						++n_triangles;
					else
						assert(false);

					//skipping tags
					for(int i = 0; i < number_of_tags; ++i)
					{
						int tmp;
						iss >> tmp;
					}

					auto &node_list = all_elements[elm_number-1];
					node_list.push_back(elm_type);

					while(iss.good())
					{
						int tmp;
						iss >> tmp;
						node_list.push_back(tmp);
					}
				}
			}
			else
			{
				assert(false);
			}
		}

		args.box_max = vertices.colwise().maxCoeff();
   		args.box_min = vertices.colwise().minCoeff();

		mesh.tri_vertices.resize(vertices.rows());
		mesh.v_is_removed.resize(vertices.rows(), true);
		mesh.v_scalars.resize(vertices.rows());

		mesh.nodes.resize(vertices.rows());
		mesh.n_is_removed.resize(vertices.rows(), true);


		for(int i = 0; i < vertices.rows(); ++i)
		{
			mesh.tri_vertices[i].posf[0]=vertices(i,0);
			mesh.tri_vertices[i].posf[1]=vertices(i,1);

			mesh.tri_vertices[i].feature_infos.push_back({{0, 0}});

			mesh.nodes[i][0]=vertices(i,0);
			mesh.nodes[i][1]=vertices(i,1);
		}

		int index = 0;


		mesh.tris.resize(n_triangles);
		mesh.tri_nodes.resize(n_triangles);
		mesh.is_boundary_es.resize(n_triangles);
		mesh.is_bbox_es.resize(n_triangles);
		mesh.tag_feature_es.resize(n_triangles);
		mesh.t_quality.resize(n_triangles);
		mesh.t_is_removed.resize(n_triangles, false);
		mesh.t_scalars.resize(n_triangles);
		mesh.tri_nodes.resize(n_triangles);


		for(const auto &els : all_elements)
		{
			if(els.empty())
				continue;

			const int elm_type = els[0];
			if(elm_type != 2 && elm_type != 9 && elm_type != 21 && elm_type != 23)
				continue;

			auto &t = mesh.tris[index];
			for(size_t i = 1; i <= 3; ++i){
				const int v_id = els[i]-1;
				mesh.v_is_removed[v_id]=false;
				t[i-1] = v_id;

				mesh.tri_vertices[v_id].conn_tris.insert(index);
			}



			auto &el = mesh.tri_nodes[index];
			for(size_t i = 4; i < els.size(); ++i){
				mesh.n_is_removed[els[i] - 1] = false;
				el.push_back(els[i] - 1);
			}

			mesh.tag_feature_es[index] = {{0, 0, 0}};


			++index;
		}
	}

}
