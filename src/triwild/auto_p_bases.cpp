// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "auto_p_bases.hpp"


namespace triwild {
	namespace autogen {
		namespace {
			template<typename InMatrix, typename OutMatrix>
			void p_1_basis_value_2d(const int local_index, const InMatrix &uv, OutMatrix &result_0){

				auto &x=uv.col(0).array();
				auto &y=uv.col(1).array();

				switch(local_index){
					case 0: {result_0 = -x - y + 1;} break;
					case 1: {result_0 = x;} break;
					case 2: {result_0 = y;} break;
					default: assert(false);
				}
			}

			template<typename Matrix>
			void p_1_basis_grad_value_2d(const int local_index, const Matrix &uv, Matrix &val){

				auto &x=uv.col(0).array();
				auto &y=uv.col(1).array();

				val.resize(uv.rows(), uv.cols());
				switch(local_index){
					case 0: {{val.col(0).setConstant(-1); }{val.col(1).setConstant(-1); }} break;
					case 1: {{val.col(0).setOnes(); }{val.col(1).setZero(); }} break;
					case 2: {{val.col(0).setZero(); }{val.col(1).setOnes(); }} break;
					default: assert(false);
				}
			}

			template<typename Matrix>
			void p_1_nodes_2d(Matrix &res) {
				res.resize(3, 2); res <<
				0, 0,
				1, 0,
				0, 1;
			}

			template<typename InMatrix, typename OutMatrix>
			void p_2_basis_value_2d(const int local_index, const InMatrix &uv, OutMatrix &result_0){

				auto &x=uv.col(0).array();
				auto &y=uv.col(1).array();

				switch(local_index){
					case 0: {result_0 = (x + y - 1)*(2*x + 2*y - 1);} break;
					case 1: {result_0 = x*(2*x - 1);} break;
					case 2: {result_0 = y*(2*y - 1);} break;
					case 3: {result_0 = -4*x*(x + y - 1);} break;
					case 4: {result_0 = 4*x*y;} break;
					case 5: {result_0 = -4*y*(x + y - 1);} break;
					default: assert(false);
				}
			}

			template<typename Matrix>
			void p_2_basis_grad_value_2d(const int local_index, const Matrix &uv, Matrix &val){

				auto &x=uv.col(0).array();
				auto &y=uv.col(1).array();

				val.resize(uv.rows(), uv.cols());
				Eigen::ArrayXd result_0(uv.rows());
				switch(local_index){
					case 0: {{result_0 = 4*x + 4*y - 3;val.col(0) = result_0; }{result_0 = 4*x + 4*y - 3;val.col(1) = result_0; }} break;
					case 1: {{result_0 = 4*x - 1;val.col(0) = result_0; }{result_0.setZero();val.col(1) = result_0; }} break;
					case 2: {{result_0.setZero();val.col(0) = result_0; }{result_0 = 4*y - 1;val.col(1) = result_0; }} break;
					case 3: {{result_0 = 4*(-2*x - y + 1);val.col(0) = result_0; }{result_0 = -4*x;val.col(1) = result_0; }} break;
					case 4: {{result_0 = 4*y;val.col(0) = result_0; }{result_0 = 4*x;val.col(1) = result_0; }} break;
					case 5: {{result_0 = -4*y;val.col(0) = result_0; }{result_0 = 4*(-x - 2*y + 1);val.col(1) = result_0; }} break;
					default: assert(false);
				}
			}

			template<typename Matrix>
			void p_2_nodes_2d(Matrix &res) {
				res.resize(6, 2); res <<
				0, 0,
				1, 0,
				0, 1,
				1.0/2.0, 0,
				1.0/2.0, 1.0/2.0,
				0, 1.0/2.0;
			}

			template<typename InMatrix, typename OutMatrix>
			void p_3_basis_value_2d(const int local_index, const InMatrix &uv, OutMatrix &result_0){

				auto &x=uv.col(0).array();
				auto &y=uv.col(1).array();

				switch(local_index){
					case 0: {
						const auto &helper_0 = x*x;
						const auto &helper_1 = y*y;
						result_0 = -27.0/2.0*helper_0*y + 9*helper_0 - 27.0/2.0*helper_1*x + 9*helper_1 - 9.0/2.0*x*x*x + 18*x*y - 11.0/2.0*x - 9.0/2.0*y*y*y - 11.0/2.0*y + 1;
					} break;
					case 1: {result_0 = (1.0/2.0)*x*(9*x*x - 9*x + 2);} break;
					case 2: {result_0 = (1.0/2.0)*y*(9*y*y - 9*y + 2);} break;
					case 3: {result_0 = (9.0/2.0)*x*(x + y - 1)*(3*x + 3*y - 2);} break;
					case 4: {result_0 = -9.0/2.0*x*(3*x*x + 3*x*y - 4*x - y + 1);} break;
					case 5: {result_0 = (9.0/2.0)*x*y*(3*x - 1);} break;
					case 6: {result_0 = (9.0/2.0)*x*y*(3*y - 1);} break;
					case 7: {result_0 = -9.0/2.0*y*(3*x*y - x + 3*y*y - 4*y + 1);} break;
					case 8: {result_0 = (9.0/2.0)*y*(x + y - 1)*(3*x + 3*y - 2);} break;
					case 9: {result_0 = -27*x*y*(x + y - 1);} break;
					default: assert(false);
				}
			}

			template<typename Matrix>
			void p_3_basis_grad_value_2d(const int local_index, const Matrix &uv, Matrix &val){

				auto &x=uv.col(0).array();
				auto &y=uv.col(1).array();

				val.resize(uv.rows(), uv.cols());
				Eigen::ArrayXd result_0(uv.rows());
				switch(local_index){
					case 0: {{result_0 = -27.0/2.0*x*x - 27*x*y + 18*x - 27.0/2.0*y*y + 18*y - 11.0/2.0;val.col(0) = result_0; }{result_0 = -27.0/2.0*x*x - 27*x*y + 18*x - 27.0/2.0*y*y + 18*y - 11.0/2.0;val.col(1) = result_0; }} break;
					case 1: {{result_0 = (27.0/2.0)*x*x - 9*x + 1;val.col(0) = result_0; }{result_0.setZero();val.col(1) = result_0; }} break;
					case 2: {{result_0.setZero();val.col(0) = result_0; }{result_0 = (27.0/2.0)*y*y - 9*y + 1;val.col(1) = result_0; }} break;
					case 3: {{result_0 = 9*((9.0/2.0)*x*x + 6*x*y - 5*x + (3.0/2.0)*y*y - 5.0/2.0*y + 1);val.col(0) = result_0; }{result_0 = (9.0/2.0)*x*(6*x + 6*y - 5);val.col(1) = result_0; }} break;
					case 4: {{result_0 = 9*(-9.0/2.0*x*x - 3*x*y + 4*x + (1.0/2.0)*y - 1.0/2.0);val.col(0) = result_0; }{result_0 = -9.0/2.0*x*(3*x - 1);val.col(1) = result_0; }} break;
					case 5: {{result_0 = (9.0/2.0)*y*(6*x - 1);val.col(0) = result_0; }{result_0 = (9.0/2.0)*x*(3*x - 1);val.col(1) = result_0; }} break;
					case 6: {{result_0 = (9.0/2.0)*y*(3*y - 1);val.col(0) = result_0; }{result_0 = (9.0/2.0)*x*(6*y - 1);val.col(1) = result_0; }} break;
					case 7: {{result_0 = -9.0/2.0*y*(3*y - 1);val.col(0) = result_0; }{result_0 = 9*(-3*x*y + (1.0/2.0)*x - 9.0/2.0*y*y + 4*y - 1.0/2.0);val.col(1) = result_0; }} break;
					case 8: {{result_0 = (9.0/2.0)*y*(6*x + 6*y - 5);val.col(0) = result_0; }{result_0 = 9*((3.0/2.0)*x*x + 6*x*y - 5.0/2.0*x + (9.0/2.0)*y*y - 5*y + 1);val.col(1) = result_0; }} break;
					case 9: {{result_0 = -27*y*(2*x + y - 1);val.col(0) = result_0; }{result_0 = -27*x*(x + 2*y - 1);val.col(1) = result_0; }} break;
					default: assert(false);
				}
			}

			template<typename Matrix>
			void p_3_nodes_2d(Matrix &res) {
				res.resize(10, 2); res <<
				0, 0,
				1, 0,
				0, 1,
				1.0/3.0, 0,
				2.0/3.0, 0,
				2.0/3.0, 1.0/3.0,
				1.0/3.0, 2.0/3.0,
				0, 2.0/3.0,
				0, 1.0/3.0,
				1.0/3.0, 1.0/3.0;
			}

			template<typename InMatrix, typename OutMatrix>
			void p_4_basis_value_2d(const int local_index, const InMatrix &uv, OutMatrix &result_0){

				auto &x=uv.col(0).array();
				auto &y=uv.col(1).array();

				switch(local_index){
					case 0: {
						const auto &helper_0 = x*x;
						const auto &helper_1 = x*x*x;
						const auto &helper_2 = y*y;
						const auto &helper_3 = y*y*y;
						result_0 = 64*helper_0*helper_2 - 80*helper_0*y + (70.0/3.0)*helper_0 + (128.0/3.0)*helper_1*y - 80.0/3.0*helper_1 - 80*helper_2*x + (70.0/3.0)*helper_2 + (128.0/3.0)*helper_3*x - 80.0/3.0*helper_3 + (32.0/3.0)*x*x*x*x + (140.0/3.0)*x*y - 25.0/3.0*x + (32.0/3.0)*y*y*y*y - 25.0/3.0*y + 1;
					} break;
					case 1: {result_0 = (1.0/3.0)*x*(32*x*x*x - 48*x*x + 22*x - 3);} break;
					case 2: {result_0 = (1.0/3.0)*y*(32*y*y*y - 48*y*y + 22*y - 3);} break;
					case 3: {
						const auto &helper_0 = x*x;
						const auto &helper_1 = y*y;
						result_0 = -16.0/3.0*x*(24*helper_0*y - 18*helper_0 + 24*helper_1*x - 18*helper_1 + 8*x*x*x - 36*x*y + 13*x + 8*y*y*y + 13*y - 3);
					} break;
					case 4: {
						const auto &helper_0 = 32*x*x;
						const auto &helper_1 = y*y;
						result_0 = 4*x*(helper_0*y - helper_0 + 16*helper_1*x - 4*helper_1 + 16*x*x*x - 36*x*y + 19*x + 7*y - 3);
					} break;
					case 5: {
						const auto &helper_0 = x*x;
						result_0 = -16.0/3.0*x*(8*helper_0*y - 14*helper_0 + 8*x*x*x - 6*x*y + 7*x + y - 1);
					} break;
					case 6: {result_0 = (16.0/3.0)*x*y*(8*x*x - 6*x + 1);} break;
					case 7: {
						const auto &helper_0 = 4*x;
						result_0 = helper_0*y*(-helper_0 + 16*x*y - 4*y + 1);
					} break;
					case 8: {result_0 = (16.0/3.0)*x*y*(8*y*y - 6*y + 1);} break;
					case 9: {
						const auto &helper_0 = y*y;
						result_0 = -16.0/3.0*y*(8*helper_0*x - 14*helper_0 - 6*x*y + x + 8*y*y*y + 7*y - 1);
					} break;
					case 10: {
						const auto &helper_0 = x*x;
						const auto &helper_1 = 32*y*y;
						result_0 = 4*y*(16*helper_0*y - 4*helper_0 + helper_1*x - helper_1 - 36*x*y + 7*x + 16*y*y*y + 19*y - 3);
					} break;
					case 11: {
						const auto &helper_0 = x*x;
						const auto &helper_1 = y*y;
						result_0 = -16.0/3.0*y*(24*helper_0*y - 18*helper_0 + 24*helper_1*x - 18*helper_1 + 8*x*x*x - 36*x*y + 13*x + 8*y*y*y + 13*y - 3);
					} break;
					case 12: {result_0 = 32*x*y*(x + y - 1)*(4*x + 4*y - 3);} break;
					case 13: {result_0 = -32*x*y*(4*y - 1)*(x + y - 1);} break;
					case 14: {result_0 = -32*x*y*(4*x - 1)*(x + y - 1);} break;
					default: assert(false);
				}
			}

			template<typename Matrix>
			void p_4_basis_grad_value_2d(const int local_index, const Matrix &uv, Matrix &val){

				auto &x=uv.col(0).array();
				auto &y=uv.col(1).array();

				val.resize(uv.rows(), uv.cols());
				Eigen::ArrayXd result_0(uv.rows());
				switch(local_index){
					case 0: {
						{
							const auto &helper_0 = x*x;
							const auto &helper_1 = y*y;
							result_0 = 128*helper_0*y - 80*helper_0 + 128*helper_1*x - 80*helper_1 + (128.0/3.0)*x*x*x - 160*x*y + (140.0/3.0)*x + (128.0/3.0)*y*y*y + (140.0/3.0)*y - 25.0/3.0;val.col(0) = result_0;
						}{
							const auto &helper_0 = x*x;
							const auto &helper_1 = y*y;
							result_0 = 128*helper_0*y - 80*helper_0 + 128*helper_1*x - 80*helper_1 + (128.0/3.0)*x*x*x - 160*x*y + (140.0/3.0)*x + (128.0/3.0)*y*y*y + (140.0/3.0)*y - 25.0/3.0;val.col(1) = result_0;
						}
					} break;
					case 1: {{result_0 = (128.0/3.0)*x*x*x - 48*x*x + (44.0/3.0)*x - 1;val.col(0) = result_0; }{result_0.setZero();val.col(1) = result_0; }} break;
					case 2: {{result_0.setZero();val.col(0) = result_0; }{result_0 = (128.0/3.0)*y*y*y - 48*y*y + (44.0/3.0)*y - 1;val.col(1) = result_0; }} break;
					case 3: {
						{
							const auto &helper_0 = 24*y;
							const auto &helper_1 = x*x;
							const auto &helper_2 = y*y;
							result_0 = -16*helper_0*helper_1 + 16*helper_0*x + 288*helper_1 - 256*helper_2*x + 96*helper_2 - 512.0/3.0*x*x*x - 416.0/3.0*x - 128.0/3.0*y*y*y - 208.0/3.0*y + 16;val.col(0) = result_0;
						}
						{
							result_0 = -16.0/3.0*x*(24*x*x + 48*x*y - 36*x + 24*y*y - 36*y + 13);val.col(1) = result_0;
						}
					} break;
					case 4: {
						{
							const auto &helper_0 = 96*x*x;
							const auto &helper_1 = y*y;
							result_0 = 4*helper_0*y - 4*helper_0 + 128*helper_1*x - 16*helper_1 + 256*x*x*x - 288*x*y + 152*x + 28*y - 12;val.col(0) = result_0;
						}
						{
							result_0 = 4*x*(32*x*x + 32*x*y - 36*x - 8*y + 7);val.col(1) = result_0;
						}
					} break;
					case 5: {
						{
							const auto &helper_0 = x*x;
							result_0 = -128*helper_0*y + 224*helper_0 - 512.0/3.0*x*x*x + 64*x*y - 224.0/3.0*x - 16.0/3.0*y + 16.0/3.0;val.col(0) = result_0;
						}
						{
							result_0 = -16.0/3.0*x*(8*x*x - 6*x + 1);val.col(1) = result_0;
						}
					} break;
					case 6: {
						{
							result_0 = (16.0/3.0)*y*(24*x*x - 12*x + 1);val.col(0) = result_0;
						}
						{
							result_0 = (16.0/3.0)*x*(8*x*x - 6*x + 1);val.col(1) = result_0;
						}
					} break;
					case 7: {
						{
							const auto &helper_0 = 4*y;
							result_0 = helper_0*(-helper_0 + 32*x*y - 8*x + 1);val.col(0) = result_0;
						}
						{
							const auto &helper_0 = 4*x;
							result_0 = helper_0*(-helper_0 + 32*x*y - 8*y + 1);val.col(1) = result_0;
						}
					} break;
					case 8: {
						{
							result_0 = (16.0/3.0)*y*(8*y*y - 6*y + 1);val.col(0) = result_0;
						}
						{
							result_0 = (16.0/3.0)*x*(24*y*y - 12*y + 1);val.col(1) = result_0;
						}
					} break;
					case 9: {
						{
							result_0 = -16.0/3.0*y*(8*y*y - 6*y + 1);val.col(0) = result_0;
						}
						{
							const auto &helper_0 = y*y;
							result_0 = -128*helper_0*x + 224*helper_0 + 64*x*y - 16.0/3.0*x - 512.0/3.0*y*y*y - 224.0/3.0*y + 16.0/3.0;val.col(1) = result_0;
						}
					} break;
					case 10: {
						{
							result_0 = 4*y*(32*x*y - 8*x + 32*y*y - 36*y + 7);val.col(0) = result_0;
						}
						{
							const auto &helper_0 = x*x;
							const auto &helper_1 = 96*y*y;
							result_0 = 128*helper_0*y - 16*helper_0 + 4*helper_1*x - 4*helper_1 - 288*x*y + 28*x + 256*y*y*y + 152*y - 12;val.col(1) = result_0;
						}
					} break;
					case 11: {
						{
							result_0 = -16.0/3.0*y*(24*x*x + 48*x*y - 36*x + 24*y*y - 36*y + 13);val.col(0) = result_0;
						}
						{
							const auto &helper_0 = 24*x;
							const auto &helper_1 = x*x;
							const auto &helper_2 = y*y;
							result_0 = -16*helper_0*helper_2 + 16*helper_0*y - 256*helper_1*y + 96*helper_1 + 288*helper_2 - 128.0/3.0*x*x*x - 208.0/3.0*x - 512.0/3.0*y*y*y - 416.0/3.0*y + 16;val.col(1) = result_0;
						}
					} break;
					case 12: {{result_0 = 32*y*(12*x*x + 16*x*y - 14*x + 4*y*y - 7*y + 3);val.col(0) = result_0; }{result_0 = 32*x*(4*x*x + 16*x*y - 7*x + 12*y*y - 14*y + 3);val.col(1) = result_0; }} break;
					case 13: {{result_0 = -32*y*(8*x*y - 2*x + 4*y*y - 5*y + 1);val.col(0) = result_0; }{result_0 = -32*x*(8*x*y - x + 12*y*y - 10*y + 1);val.col(1) = result_0; }} break;
					case 14: {
						{
							result_0 = -32*y*(12*x*x + 8*x*y - 10*x - y + 1);val.col(0) = result_0;
						}
						{
							result_0 = -32*x*(4*x*x + 8*x*y - 5*x - 2*y + 1);val.col(1) = result_0;
						}
					} break;
					default: assert(false);
				}
			}

			template<typename Matrix>
			void p_4_nodes_2d(Matrix &res) {
				res.resize(15, 2); res <<
				0, 0,
				1, 0,
				0, 1,
				1.0/4.0, 0,
				1.0/2.0, 0,
				3.0/4.0, 0,
				3.0/4.0, 1.0/4.0,
				1.0/2.0, 1.0/2.0,
				1.0/4.0, 3.0/4.0,
				0, 3.0/4.0,
				0, 1.0/2.0,
				0, 1.0/4.0,
				1.0/4.0, 1.0/4.0,
				1.0/4.0, 1.0/2.0,
				1.0/2.0, 1.0/4.0;
			}
		}

		template<typename Matrix>
		void p_nodes_2d(const int p, Matrix &val){
			switch(p){
				case 1: p_1_nodes_2d(val); break;
				case 2: p_2_nodes_2d(val); break;
				case 3: p_3_nodes_2d(val); break;
				case 4: p_4_nodes_2d(val); break;
				default: assert(false);
			}
		}

		template<typename InMatrix, typename OutMatrix>
		void p_basis_value_2d(const int p, const int local_index, const InMatrix &uv, OutMatrix &val){
			switch(p){
				case 1: p_1_basis_value_2d(local_index, uv, val); break;
				case 2: p_2_basis_value_2d(local_index, uv, val); break;
				case 3: p_3_basis_value_2d(local_index, uv, val); break;
				case 4: p_4_basis_value_2d(local_index, uv, val); break;
				default: assert(false);
			}
		}

		template<typename Matrix>
		void p_grad_basis_value_2d(const int p, const int local_index, const Matrix &uv, Matrix &val){
			switch(p){
				case 1: p_1_basis_grad_value_2d(local_index, uv, val); break;
				case 2: p_2_basis_grad_value_2d(local_index, uv, val); break;
				case 3: p_3_basis_grad_value_2d(local_index, uv, val); break;
				case 4: p_4_basis_grad_value_2d(local_index, uv, val); break;
				default: assert(false);
			}
		}


		//template instantiation
		template void p_nodes_2d<Eigen::MatrixXd>(const int, Eigen::MatrixXd&);
		template void p_nodes_2d<Eigen::Matrix<double, Eigen::Dynamic, 2, 0, 15, 2>>(const int, Eigen::Matrix<double, Eigen::Dynamic, 2, 0, 15, 2>&);

		template void p_grad_basis_value_2d<Eigen::MatrixXd>(const int, const int, const Eigen::MatrixXd&, Eigen::MatrixXd&);
		template void p_grad_basis_value_2d<Eigen::Matrix<double, 1, 2, 1, 1, 2>>(const int, const int, const Eigen::Matrix<double, 1, 2, 1, 1, 2>&, Eigen::Matrix<double, 1, 2, 1, 1, 2>&);

		template void p_basis_value_2d<Eigen::MatrixXd, Eigen::MatrixXd>(const int, const int, const Eigen::MatrixXd&, Eigen::MatrixXd&);
		template void p_basis_value_2d<Eigen::Matrix<double, Eigen::Dynamic, 2, 0, 15, 2>, Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 15, 1> >(int, int, Eigen::Matrix<double, Eigen::Dynamic, 2, 0, 15, 2> const&, Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 15, 1>&);
	}
}

