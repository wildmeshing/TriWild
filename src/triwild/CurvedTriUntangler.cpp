// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "CurvedTriUntangler.hpp"

#include "Curves.h"
#include "auto_p_bases.hpp"

//#include <igl/triangle/triangulate.h>
#include <Eigen/Dense>
#include <nlopt.hpp>
#include <cassert>

#include "reference_triangle.h"

namespace triwild {
	namespace
	{

		typedef struct {
			std::array<Point_2f, 3> vertices;
			int edge;

			Point_2f fixed_node1;
			Point_2f fixed_node2;

			std::vector<Eigen::MatrixXd> grads;

			Eigen::VectorXd f1dx;
			Eigen::VectorXd f2dx;

			Eigen::VectorXd f1dy;
			Eigen::VectorXd f2dy;

			std::vector<double> nodes;
		} tri_data;


		typedef struct {
			Eigen::VectorXd f1dx_pre;
			Eigen::VectorXd f2dx_pre;

			Eigen::VectorXd f1dy_pre;
			Eigen::VectorXd f2dy_pre;


			Eigen::VectorXd f1dx;
			Eigen::VectorXd f2dx;

			Eigen::VectorXd f1dy;
			Eigen::VectorXd f2dy;

			Eigen::MatrixXd grad_center;
		} center_tri_data;


		typedef struct {
			std::vector<Point_2f> curve_samples;

			Point_2f vertex1, vertex2;
			int n1_index, n2_index;


			std::vector<Eigen::MatrixXd> grads;

			// Eigen::VectorXd f1dx;
			// Eigen::VectorXd f2dx;

			// Eigen::VectorXd f1dy;
			// Eigen::VectorXd f2dy;

			std::array<Point_2f, 3> vertices;
			std::vector<Point_2f> nodes;

			int point_index;
		} ls_data_global;

	 	struct ls_data_main{
			ls_data_global &data;
			int point_index;

			ls_data_main(ls_data_global &data, int point_index)
			: data(data), point_index(point_index)
			{ }

		};


		void lagr3(const std::array<Point_2f, 4> &c, const double t, Point_2f &res)
		{
			const double l0 = (t-1./3.)/(0./3.-1./3.) * (t-2./3.)/(0./3.-2./3.) * (t-3./3.)/(0./3.-3./3.);
			const double l1 = (t-0./3.)/(1./3.-0./3.) * (t-2./3.)/(1./3.-2./3.) * (t-3./3.)/(1./3.-3./3.);
			const double l2 = (t-0./3.)/(2./3.-0./3.) * (t-1./3.)/(2./3.-1./3.) * (t-3./3.)/(2./3.-3./3.);
			const double l3 = (t-0./3.)/(3./3.-0./3.) * (t-1./3.)/(3./3.-1./3.) * (t-2./3.)/(3./3.-2./3.);

			res = c[0]*l0 + c[1]*l1 + c[2]*l2 + c[3]*l3;
		}


		void get_tri(const std::vector<double> &params, const tri_data &tmp, std::vector<double> &res)
		{
			assert(params.size() == 6);

			const double t1 = params[0];
			const double t2 = params[1];

			const double t3 = params[2];
			const double t4 = params[3];

			const double px = params[4];
			const double py = params[5];


			const int i1 = (tmp.edge+1) % 3;
			const int i2 = (i1+1) % 3;
			const int i3 = (i2+1) % 3;

			const double p1x = (1-t1)*tmp.vertices[i1].x + t1*tmp.vertices[i2].x;
			const double p1y = (1-t1)*tmp.vertices[i1].y + t1*tmp.vertices[i2].y;

			const double p2x = (1-t2)*tmp.vertices[i1].x + t2*tmp.vertices[i2].x;
			const double p2y = (1-t2)*tmp.vertices[i1].y + t2*tmp.vertices[i2].y;

			const double p3x = (1-t3)*tmp.vertices[i2].x + t3*tmp.vertices[i3].x;
			const double p3y = (1-t3)*tmp.vertices[i2].y + t3*tmp.vertices[i3].y;

			const double p4x = (1-t4)*tmp.vertices[i2].x + t4*tmp.vertices[i3].x;
			const double p4y = (1-t4)*tmp.vertices[i2].y + t4*tmp.vertices[i3].y;

			const auto &bp1 = tmp.fixed_node1;
			const auto &bp2 = tmp.fixed_node2;

			res.resize(20);

			int index = 0;

			res[index++] = tmp.vertices[0].x; res[index++] = tmp.vertices[0].y;
			res[index++] = tmp.vertices[1].x; res[index++] = tmp.vertices[1].y;
			res[index++] = tmp.vertices[2].x; res[index++] = tmp.vertices[2].y;

			if(tmp.edge == 0)
			{
				res[index++] = bp1[0]; res[index++] = bp1[1];
				res[index++] = bp2[0]; res[index++] = bp2[1];

				res[index++] = p1x; res[index++] = p1y;
				res[index++] = p2x; res[index++] = p2y;

				res[index++] = p3x; res[index++] = p3y;
				res[index++] = p4x; res[index++] = p4y;
			}
			else if(tmp.edge == 1)
			{
				res[index++] = p1x; res[index++] = p1y;
				res[index++] = p2x; res[index++] = p2y;

				res[index++] = bp1[0]; res[index++] = bp1[1];
				res[index++] = bp2[0]; res[index++] = bp2[1];

				res[index++] = p3x; res[index++] = p3y;
				res[index++] = p4x; res[index++] = p4y;
			}
			else
			{
				res[index++] = p1x; res[index++] = p1y;
				res[index++] = p2x; res[index++] = p2y;

				res[index++] = p3x; res[index++] = p3y;
				res[index++] = p4x; res[index++] = p4y;

				res[index++] = bp1[0]; res[index++] = bp1[1];
				res[index++] = bp2[0]; res[index++] = bp2[1];
			}

			res[index++] = px; res[index++] = py;
			assert(index == 20);
		}

		double objective_fun(const std::vector<double> &params, std::vector<double> &grad_obj, void *data)
		{
			//never use the gradient
			assert(grad_obj.empty());
			tri_data &tmp = *reinterpret_cast<tri_data *>(data);

			std::vector<double> &nodes = tmp.nodes;
			get_tri(params, tmp, nodes);

			Eigen::VectorXd &f1dx = tmp.f1dx; f1dx.setZero();
			Eigen::VectorXd &f2dx = tmp.f2dx; f2dx.setZero();

			Eigen::VectorXd &f1dy = tmp.f1dy; f1dy.setZero();
			Eigen::VectorXd &f2dy = tmp.f2dy; f2dy.setZero();

			for(size_t i = 0; i < nodes.size(); i+=2)
			{
				// autogen::p_grad_basis_value_2d(3, i/2, tmp.points, grad);
				const Eigen::MatrixXd &grad = tmp.grads[i/2];

				f1dx += grad.col(0)*nodes[i];
				f2dx += grad.col(0)*nodes[i+1];

				f1dy += grad.col(1)*nodes[i];
				f2dy += grad.col(1)*nodes[i+1];
			}

			const Eigen::VectorXd ddet = f1dx.array()*f2dy.array() - f1dy.array()*f2dx.array();

			double min_det = 0;

			// min_det = -sum(ddet(ddet<0));
			for(int i = 0; i < ddet.size(); ++i)
			{
				if(ddet(i) < 0)
					min_det += ddet(i);
			}

			if(fabs(min_det)<1e-10)
				min_det = ddet.minCoeff();

			return min_det;
		}


		double center_objective_fun(const std::vector<double> &params, std::vector<double> &grad_obj, void *data)
		{
			//never use the gradient
			assert(grad_obj.empty());
			center_tri_data &tmp = *reinterpret_cast<center_tri_data *>(data);

			Eigen::VectorXd &f1dx = tmp.f1dx; f1dx = tmp.f1dx_pre;
			Eigen::VectorXd &f2dx = tmp.f2dx; f2dx = tmp.f2dx_pre;

			Eigen::VectorXd &f1dy = tmp.f1dy; f1dy = tmp.f1dy_pre;
			Eigen::VectorXd &f2dy = tmp.f2dy; f2dy = tmp.f2dy_pre;

			f1dx += tmp.grad_center.col(0)*params[0];
			f2dx += tmp.grad_center.col(0)*params[1];

			f1dy += tmp.grad_center.col(1)*params[0];
			f2dy += tmp.grad_center.col(1)*params[1];


			const Eigen::VectorXd ddet = f1dx.array()*f2dy.array() - f1dy.array()*f2dx.array();

			double min_det = 0;

			// min_det = -sum(ddet(ddet<0));
			for(int i = 0; i < ddet.size(); ++i)
			{
				if(ddet(i) < 0)
					min_det += ddet(i);
			}

			if(fabs(min_det)<1e-10)
				min_det = ddet.minCoeff();

			return min_det;
		}


		double ls_energy(const std::vector<double> &params, std::vector<double> &grad_obj, void *data)
		{
			//never use the gradient
			assert(grad_obj.empty());

			static const int n_pts = 10;

			const ls_data_global &tmp = reinterpret_cast<ls_data_main *>(data)->data;

			std::array<Point_2f, 4> ctrl = {{tmp.vertex1, Point_2f(params[0], params[1]), Point_2f(params[2], params[3]), tmp.vertex2}};


			double dist = 0;
			Point_2f p;

			for(int i = 0; i <= n_pts; ++i)
			{
				double min_dist = std::numeric_limits<double>::max();

				lagr3(ctrl, (1.0*i)/n_pts, p);

				for(size_t j=0; j < tmp.curve_samples.size(); ++j)
				{
					const double dd = (p - tmp.curve_samples[j]).length_2();
					min_dist=std::min(dd, min_dist);
					if(min_dist < 1e-10)
						break;
				}

				dist += min_dist;
			}

			return dist;
		}


		double ls_constraints(const std::vector<double> &params, std::vector<double> &grad_constraint, void *data)
		{
			//never use the gradient
			assert(grad_constraint.empty());
			ls_data_main &tmp1 = *reinterpret_cast<ls_data_main *>(data);
			ls_data_global &tmp = tmp1.data;

			const std::array<Point_2f, 3> &vertices = tmp.vertices;

			std::vector<Point_2f> &nodes = tmp.nodes;
			nodes[tmp.n1_index][0] = params[0];
			nodes[tmp.n1_index][1] = params[1];

			nodes[tmp.n2_index][0] = params[2];
			nodes[tmp.n2_index][1] = params[3];

			// std::cout<<"--------------------"<<std::endl;
			// for(const auto &p : vertices)
			// 	std::cout<<p<<std::endl;

			// for(const auto &p : nodes)
			// 	std::cout<<p<<std::endl;


			// Eigen::VectorXd &f1dx = tmp.f1dx; f1dx.setZero();
			// Eigen::VectorXd &f2dx = tmp.f2dx; f2dx.setZero();

			// Eigen::VectorXd &f1dy = tmp.f1dy; f1dy.setZero();
			// Eigen::VectorXd &f2dy = tmp.f2dy; f2dy.setZero();

			Eigen::Matrix2d J; J.setZero();
			const int point_index = tmp1.point_index;

			for(size_t i = 0; i < vertices.size(); ++i)
			{
				// autogen::p_grad_basis_value_2d(3, i/2, tmp.points, grad);
				const Eigen::MatrixXd &grad = tmp.grads[i];

				J(0, 0) += grad(point_index, 0)*vertices[i][0];
				J(0, 1) += grad(point_index, 0)*vertices[i][1];

				J(1, 0) += grad(point_index, 1)*vertices[i][0];
				J(1, 1) += grad(point_index, 1)*vertices[i][1];
			}

			// 	f1dx += grad.col(0)*vertices[i][0];
			// 	f2dx += grad.col(0)*vertices[i][1];

			// 	f1dy += grad.col(1)*vertices[i][0];
			// 	f2dy += grad.col(1)*vertices[i][1];
			// }

			assert(nodes.size() == 7);
			for(size_t i = 0; i < nodes.size(); ++i)
			{
				// autogen::p_grad_basis_value_2d(3, i/2, tmp.points, grad);
				const Eigen::MatrixXd &grad = tmp.grads[3+i];

				J(0, 0) += grad(point_index, 0)*nodes[i][0];
				J(0, 1) += grad(point_index, 0)*nodes[i][1];

				J(1, 0) += grad(point_index, 1)*nodes[i][0];
				J(1, 1) += grad(point_index, 1)*nodes[i][1];

				// f1dx += grad.col(0)*nodes[i][0];
				// f2dx += grad.col(0)*nodes[i][1];

				// f1dy += grad.col(1)*nodes[i][0];
				// f2dy += grad.col(1)*nodes[i][1];
			}

			// const Eigen::VectorXd ddet = f1dx.array()*f2dy.array() - f1dy.array()*f2dx.array();
			const double ddet = J.determinant();

			// return (-ddet.minCoeff() + 1e-6)*10;
			return (-ddet + 1e-6);
		}
	}

	bool CurvedTriUntangler::untangle(const std::array<Point_2f, 3>& vertices, const int edge,  const Point_2f &fixed1, const Point_2f &fixed2, const Point_2f &center, std::vector<Point_2f> &new_nodes)
	{
		//only cubic

		tri_data data;
		data.fixed_node1 = fixed1;
		data.fixed_node2 = fixed2;
		data.vertices = vertices;
		data.edge = edge;


		Eigen::MatrixXd V(3, 2), V_out, _, tmp;
		Eigen::MatrixXi E(3, 2), F;
		V << 0, 0,
		1, 0,
		0, 1;
		E << 0, 1,
		1, 2,
		2, 0;
//		igl::triangle::triangulate(V, E, _, "Qq10a" + std::to_string(0.0001), V_out, F);
        V_out = get_reference_triangle_vertices();
        F = get_reference_triangle_faces();

		// std::cout<<V_out.rows()<<std::endl;
		// std::cout<<V_out<<std::endl;

		// data.points = V_out;
		data.grads.resize(10);
		for(size_t i = 0; i < 10; ++i)
		{
			autogen::p_grad_basis_value_2d(3, i, V_out, data.grads[i]);
		}

		data.f1dx.resize(data.grads.front().rows());
		data.f2dx.resize(data.grads.front().rows());

		data.f1dy.resize(data.grads.front().rows());
		data.f2dy.resize(data.grads.front().rows());

		 /* algorithm and dimensionality */
		nlopt::opt opt(nlopt::LN_SBPLX, 6);
		opt.set_max_objective(objective_fun, &data);

		std::vector<double> params = {1./3., 2./3., 1./3., 2./3., center.x, center.y};

		double min_j;


		// std::vector<double> asd;
		// get_tri(params, data, asd);
		// for(auto p : asd)
		// 	std::cout<<p<<std::endl;
		// asd.clear();
		// std::cout<<"objective_fun "<<objective_fun(params, asd, &data)<<std::endl;
		// return false;

		try{
			nlopt::result error_code = opt.optimize(params, min_j);

			// std::cout<<"min_j "<<min_j<<std::endl;

			if (error_code < 0) {
				std::cerr<<"nlopt failed! "<< error_code <<std::endl;
				return false;
			}

			if(min_j < 0)
			{
				return false;
			}
		}
		catch(std::exception &e) {
			std::cout << "nlopt failed: " << e.what() << std::endl;
		}

		//		 for(auto p : params)
		//		 	std::cout<<p<<" ";
		//		 std::cout<<std::endl;

		 for(int i = 0; i < 4; ++i)
		 {
		 	if(params[i] <= 0 || params[i] >= 1 )
		 		return false;
		 }

		get_tri(params, data, data.nodes);
		assert(data.nodes.size() == 20);
		new_nodes.resize(7);

		for(int i = 3*2; i < data.nodes.size(); i+=2)
		{
			const int index = (i-3*2)/2;
			new_nodes[index].x = data.nodes[i];
			new_nodes[index].y = data.nodes[i+1];
		}

		return true;
	}



	bool CurvedTriUntangler::untangle_center(const std::array<Point_2f, 3>& vertices, const std::vector<Point_2f> &fixed, const Point_2f &center, Point_2f &new_center)
	{
		//only cubic
		center_tri_data data;

		Eigen::MatrixXd V(3, 2), V_out, _, tmp;
		Eigen::MatrixXi E(3, 2), F;
		V << 0, 0,
		1, 0,
		0, 1;
		E << 0, 1,
		1, 2,
		2, 0;
//		igl::triangle::triangulate(V, E, _, "Qq10a" + std::to_string(0.0001), V_out, F);
        V_out = get_reference_triangle_vertices();
        F = get_reference_triangle_faces();

		// std::cout<<V_out.rows()<<std::endl;
		// std::cout<<V_out<<std::endl;

		// data.points = V_out;
		autogen::p_grad_basis_value_2d(3, 9, V_out, data.grad_center);


		Eigen::VectorXd &f1dx_pre = data.f1dx_pre; f1dx_pre.resize(V_out.rows()); f1dx_pre.setZero();
		Eigen::VectorXd &f2dx_pre = data.f2dx_pre; f2dx_pre.resize(V_out.rows()); f2dx_pre.setZero();

		Eigen::VectorXd &f1dy_pre = data.f1dy_pre; f1dy_pre.resize(V_out.rows()); f1dy_pre.setZero();
		Eigen::VectorXd &f2dy_pre = data.f2dy_pre; f2dy_pre.resize(V_out.rows()); f2dy_pre.setZero();

		Eigen::MatrixXd grad;
		for(size_t i = 0; i < vertices.size(); i++)
		{
			autogen::p_grad_basis_value_2d(3, i, V_out, grad);

			f1dx_pre += grad.col(0)*vertices[i][0];
			f2dx_pre += grad.col(0)*vertices[i][1];

			f1dy_pre += grad.col(1)*vertices[i][0];
			f2dy_pre += grad.col(1)*vertices[i][1];
		}

		for(size_t i = 0; i < 6; i++)
		{
			autogen::p_grad_basis_value_2d(3, 3+i, V_out, grad);

			f1dx_pre += grad.col(0)*fixed[i][0];
			f2dx_pre += grad.col(0)*fixed[i][1];

			f1dy_pre += grad.col(1)*fixed[i][0];
			f2dy_pre += grad.col(1)*fixed[i][1];
		}

		data.f1dx.resize(data.f1dx_pre.rows());
		data.f2dx.resize(data.f2dx_pre.rows());

		data.f1dy.resize(data.f1dy_pre.rows());
		data.f2dy.resize(data.f2dy_pre.rows());

		 /* algorithm and dimensionality */
		nlopt::opt opt(nlopt::LN_SBPLX, 2);
		opt.set_max_objective(center_objective_fun, &data);
		//		opt.set_maxeval(100);
		std::vector<double> params = {center.x, center.y};

		double min_j;


		// std::vector<double> asd;
		// get_tri(params, data, asd);
		// for(auto p : asd)
		// 	std::cout<<p<<std::endl;
		// asd.clear();
		// std::cout<<"objective_fun "<<objective_fun(params, asd, &data)<<std::endl;
		// return false;

		try{
			nlopt::result error_code = opt.optimize(params, min_j);

			// std::cout<<"min_j "<<min_j<<std::endl;

			if (error_code < 0) {
				std::cerr<<"nlopt failed! "<< error_code <<std::endl;
				return false;
			}

			if(min_j < 0)
			{
				return false;
			}
		}
		catch(std::exception &e) {
			std::cout << "nlopt failed: " << e.what() << std::endl;
		}

		// for(auto p : params)
		// 	std::cout<<p<<" ";
		// std::cout<<std::endl;


		new_center[0] = params[0];
		new_center[1] = params[1];

		return true;
	}

	double CurvedTriUntangler::ls_fit(const std::array<Point_2f, 3>& vertices, const std::vector<Point_2f> &nodes, const int edge, std::vector<Point_2f> &new_nodes, const bool all) {
		static const int n_samples = 100;

		const int v1 = edge;
		const int v2 = (edge + 1) % 3;

		const int n1 = v1 * 2;
		const int n2 = n1 + 1;

		std::array<Point_2f, 4> ctrl = {{vertices[v1], nodes[n1], nodes[n2], vertices[v2]}};

		// std::cout<<ctrl[0]<<" "<<ctrl[1]<<" "<<ctrl[2]<<" "<<ctrl[3]<<std::endl;

		ls_data_global data;
		data.curve_samples.resize(n_samples + 1);
		data.vertex1 = vertices[v1];
		data.vertex2 = vertices[v2];

		data.vertices = vertices;

		for (int i = 0; i <= n_samples; ++i) {
			const double t = (1.0 * i) / n_samples;
			lagr3(ctrl, t, data.curve_samples[i]);
			// std::cout<<data.curve_samples[i]<<std::endl;
		}

		std::vector<double> initial_pos = {
				(1 - 1. / 3.) * vertices[v1][0] + 1. / 3. * vertices[v2][0],
				(1 - 1. / 3.) * vertices[v1][1] + 1. / 3. * vertices[v2][1],
				(1 - 2. / 3.) * vertices[v1][0] + 2. / 3. * vertices[v2][0],
				(1 - 2. / 3.) * vertices[v1][1] + 2. / 3. * vertices[v2][1]
		};


		Eigen::MatrixXd V(3, 2), V_out, _, tmp;
		Eigen::MatrixXi E(3, 2), F;
		V << 0, 0,
				1, 0,
				0, 1;
		E << 0, 1,
				1, 2,
				2, 0;
//		igl::triangle::triangulate(V, E, _, "Qq10a" + std::to_string(0.0001), V_out, F);
        V_out = get_reference_triangle_vertices();
        F = get_reference_triangle_faces();

		// std::cout<<V_out.rows()<<std::endl;
		// std::cout<<V_out<<std::endl;

		// data.points = V_out;
		data.grads.resize(10);
		for (size_t i = 0; i < 10; ++i) {
			autogen::p_grad_basis_value_2d(3, i, V_out, data.grads[i]);
		}

		// data.f1dx.resize(data.grads.front().rows());
		// data.f2dx.resize(data.grads.front().rows());

		// data.f1dy.resize(data.grads.front().rows());
		// data.f2dy.resize(data.grads.front().rows());

		data.n1_index = n1;
		data.n2_index = n2;
		data.nodes = nodes;
		// data.nodes.back() = (vertices[0]+vertices[1]+vertices[2])/3.;


		nlopt::opt opt(nlopt::LN_COBYLA, 4);
		ls_data_main tmp_data(data, -1);
		opt.set_min_objective(ls_energy, &tmp_data);

		std::vector<ls_data_main> constraints;
		for(int i = 0; i < V_out.rows(); ++i)
			constraints.emplace_back(data, i);
		for(size_t i = 0; i < constraints.size(); ++i)
			opt.add_inequality_constraint(ls_constraints, &constraints[i], 0);
		opt.set_maxeval(100);


		std::vector<double> grad;
//#ifndef NDEBUG
		std::cout << "initial energy " << ls_energy(initial_pos, grad, &tmp_data) << std::endl;
		double max = -1000000;
		for(size_t i = 0; i < constraints.size(); ++i)
			max = std::max(ls_constraints(initial_pos, grad, &constraints[i]), max);
		std::cout << "initial constraint " << max << std::endl;
//#endif

		double min_dist;

		try {
			nlopt::result error_code = opt.optimize(initial_pos, min_dist);

			// std::cout<<"min_j "<<min_j<<std::endl;

			if (error_code < 0) {
				std::cerr << "nlopt failed! " << error_code << std::endl;
				return -999999;
			}
		}
		catch (std::exception &e) {
			std::cout << "nlopt failed: " << e.what() << std::endl;
			return -999999;
		}

//#ifndef NDEBUG
		std::cout << "avg_min_dist = " << min_dist / 11 << std::endl;
		max = -1000000;
		for(size_t i = 0; i < constraints.size(); ++i)
			max = std::max(ls_constraints(initial_pos, grad, &constraints[i]), max);
		std::cout << "final constraint " << max << std::endl;
//#endif

		for(size_t i = 0; i < constraints.size(); ++i)
		{
			const double ci = ls_constraints(initial_pos, grad, &constraints[i]);
			if(ci > 0)
				return -min_dist/11;
		}


		// for(auto p : params)
		// 	std::cout<<p<<" ";
		// std::cout<<std::endl;

		if (all) {
			new_nodes = data.nodes;

			new_nodes[n1][0] = initial_pos[0];
			new_nodes[n1][1] = initial_pos[1];

			new_nodes[n2][0] = initial_pos[2];
			new_nodes[n2][1] = initial_pos[3];
		} else {
			new_nodes.resize(2);
			new_nodes[0][0] = initial_pos[0];
			new_nodes[0][1] = initial_pos[1];

			new_nodes[1][0] = initial_pos[2];
			new_nodes[1][1] = initial_pos[3];
		}

		return min_dist / 11;

	}
}
