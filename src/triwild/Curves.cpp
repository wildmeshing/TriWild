// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "Curves.h"

#include "optimization.h"
#include "auto_p_bases.hpp"
#include "auto_det_checker.hpp"

#include <functional>
#include <cassert>
#include <cmath>
#include <iostream>

namespace triwild {
	namespace
	{
		double radians(const double angle)
		{
			return angle/180*M_PI;
		}

		template<typename Gamma>
		double compute_distance(const std::array<double, 2> &p, const double t, const Gamma &gamma)
		{
			// const auto gamma_val = Bezier::interpolate(ctrls, t);
			const auto gamma_val = gamma(t);
			return (gamma_val[0]-p[0])*(gamma_val[0]-p[0])+(gamma_val[1]-p[1])*(gamma_val[1]-p[1]);
		}

		template<typename Gamma>
		double compute_distance_grad(const std::array<double, 2> &p, const double t, const Gamma &gamma, const Gamma &gamma_prime)
		{
			// const auto gamma_val = Bezier::interpolate(ctrls, t);
			// const auto gamma_prime_val = Bezier::first_derivative(ctrls, t);

			const auto gamma_val = gamma(t);
			const auto gamma_prime_val = gamma_prime(t);

			return 2*(gamma_val[0]-p[0])*gamma_prime_val[0] + 2*(gamma_val[1]-p[1])*gamma_prime_val[1];
		}

		template<typename Gamma>
		double compute_distance_hessian(const std::array<double, 2> &p, const double t, const Gamma &gamma, const Gamma &gamma_prime, const Gamma &gamma_second)
		{
			// const auto gamma = Bezier::interpolate(ctrls, t);
			// const auto gamma_prime_val = Bezier::first_derivative(ctrls, t);
			// const auto gamma_second_val = Bezier::second_derivative(ctrls, t);

			const auto gamma_val = gamma(t);
			const auto gamma_prime_val = gamma_prime(t);
			const auto gamma_second_val = gamma_second(t);

			return 2*(gamma_val[0]-p[0])*gamma_second_val[0] + 2*(gamma_val[1]-p[1])*gamma_second_val[1] + 2*gamma_prime_val[0]*gamma_prime_val[0] + 2*gamma_prime_val[1]*gamma_prime_val[1];
		}

		template<typename Gamma>
		double armijo_linesearch(const std::array<double, 2> &p, const double t, const double searchDir, const Gamma &gamma, const Gamma &gamma_prime, const double alpha_init = 1.0)
		{
			static const int MAX_STEP_SIZE_ITER = 12;

			const double c = 0.5;
			const double tau = 0.5;

			double alpha = alpha_init;
			double f = compute_distance(p, t + alpha * searchDir, gamma);
			const double f_in = compute_distance(p, t, gamma);

			double grad = compute_distance_grad(p, t, gamma, gamma_prime);
			const double Cache = c * grad * searchDir;

			int cur_iter = 0;
			while(f > f_in + alpha * Cache && cur_iter < MAX_STEP_SIZE_ITER) {
				alpha *= tau;
				f = compute_distance(p, t + alpha * searchDir, gamma);

				cur_iter++;
			}

			return alpha;
		}

		template<typename Gamma>
		double newton_inverse_interpolation(const std::array<double, 2> &p, const double start, const Gamma &gamma, const Gamma &gamma_prime, const Gamma &gamma_second)
		{
#if 0
			//Check derivatives with FD
			const double t = ((double)std::rand())/RAND_MAX;
			const double dt = 1e-8;

			const double f = compute_distance(p, t, gamma);
			const double fdt = compute_distance(p, t+dt, gamma);
			const double ffd = (fdt - f)/dt;

			const double fp = compute_distance_grad(p, t, gamma, gamma_prime);
//			cout<<"fabs((fp - ffd)/fp) = "<<fabs((fp - ffd)/fp)<<endl;
//			assert(fabs((fp - ffd)/fp) < 1e-4);
			assert(fabs((fp - ffd)/(fabs(fp) > 1 ? fp : 1)) < 1e-3);

			const double fpdt = compute_distance_grad(p, t+dt, gamma, gamma_prime);
			const double fpfd = (fpdt - fp)/dt;
			const double fpp = compute_distance_hessian(p, t, gamma, gamma_prime, gamma_second);
			assert(fabs((fpfd - fpp)/fpp) < 1e-3);
#endif
			static const int MAX_ITER = 100;
			static const double GRAD_TOLERANCE = 1e-10;

			bool status = true;
			size_t iter = 0;

			double current_t = start;
			do
			{
				const double hessian = compute_distance_hessian(p, current_t, gamma, gamma_prime, gamma_second);
				const double gradient = compute_distance_grad(p, current_t, gamma, gamma_prime);
				const double delta_t = - gradient / hessian;

				const double rate = armijo_linesearch(p, current_t, delta_t, gamma, gamma_prime);
				current_t += rate * delta_t;

				++iter;

				const double energy = compute_distance(p, current_t, gamma);
				const double step = fabs(rate * delta_t);

				const double abs_grad = fabs(gradient);


				if (abs_grad < GRAD_TOLERANCE)
					status = false;
				else if (iter > MAX_ITER){
					status = false;
					// std::cerr<<"Max iters"<<std::endl;
				}
				else if (step < 1e-10){
					status = false;
					// std::cerr<<"Small step"<<std::endl;
				}

				// std::cout<<iter<<" "<<step<<" "<<rate<<" "<<current_t<<" "<<abs_grad<<" "<<energy<<std::endl;
			}
			while (status);

			return current_t;
		}

	// 	template<int N, typename Gamma>
	// 	double find_roots_in_interval(std::array<double, N> &intervals, const std::array<double, 2> &p, const Gamma &gamma, const Gamma &gamma_prime, const Gamma &gamma_second, const double t0, const double t1)
	// 	{
	// 		static const int MAX_ITER = 20;
	// 		assert(t0 >= 0);
	// 		assert(t1 <= 1);
	// 		assert(t0 < t1);

	// 		std::sort(intervals.begin(), intervals.end());

	// 		double min_t = 0;
	// 		double min_dist = std::numeric_limits<double>::max();


	// 		for(int i = 0; i < N-1; ++i)
	// 		{
	// 			// std::cout<<"considering "<<intervals[i]<<" "<<intervals[i+1]<<std::endl;

	// 			//ignoring non real roots
	// 			if(std::isnan(intervals[i]) || std::isnan(intervals[i+1]))
	// 				continue;

	// 			//both are negative
	// 			if(intervals[i+1] < 0)
	// 				continue;

	// 			//both are larger than 1
	// 			if(intervals[i] > 1)
	// 				break;

	// 			const double start = std::max(t0, intervals[i]);
	// 			const double end = std::min(t1, intervals[i+1]);

	// 			if(end - start < 1e-12)
	// 				continue;
	// 			// std::cout<<"checking "<<start<<", "<<end<<" -> "<<mid<<std::endl;


	// 			double current_t = -1;

	// 			const double dps = compute_distance_grad(p, start, gamma, gamma_prime);
	// 			const double dpe = compute_distance_grad(p, end, gamma, gamma_prime);

	// 			//if both derivative are positive or negative, no min in the interval
	// 			if(dps*dpe > 0)
	// 			{
	// 				const double ds = compute_distance(p, start, gamma);
	// 				const double de = compute_distance(p, end, gamma);

	// 				if(ds < de)
	// 					current_t = start;
	// 				else
	// 					current_t = end;
	// 			}
	// 			//otherwise use newton to check
	// 			else
	// 			{
	// 				double current_start = start;
	// 				double current_end = end;
	// 				for(int j = 0; j < MAX_ITER; ++j)
	// 				{
	// 					const double mid = (current_start + current_end)/2;
	// 					current_t = newton_inverse_interpolation(p, mid, gamma, gamma_prime, gamma_second);

	// 					if(current_t >= start && current_t <= end){
	// 						break;
	// 					}

	// 					const double dps = compute_distance_grad(p, current_start, gamma, gamma_prime);
	// 					const double dpe = compute_distance_grad(p, current_end, gamma, gamma_prime);
	// 					const double dp = compute_distance_grad(p, mid, gamma, gamma_prime);

	// 					if(fabs(dps)<1e-15)
	// 					{
	// 						current_t = current_start;
	// 						break;
	// 					}


	// 					if(fabs(dpe)<1e-15)
	// 					{
	// 						current_t = current_end;
	// 						break;
	// 					}

	// 					//both have the same sign, the min is in [mid, current_end]
	// 					if(dps*dp > 0)
	// 						current_start = mid;
	// 					else
	// 					{
	// 						//both have the same sign, the min is in [current_start, mid]
	// 						assert(dp*dpe > 0);
	// 						current_end = mid;
	// 					}
	// 				}

	// 				//we should have converged in the interval
	// 				assert(current_t >= start);
	// 				assert(current_t <= end);
	// 			}

	// 			if(current_t < start || current_t > end)
	// 			{
	// 				// std::cerr<<"Converged to "<<current_t<<" out or range  "<<start<<", "<<end<<std::endl;
	// 				continue;
	// 			}

	// 			const double dist = compute_distance(p, current_t, gamma);

	// 			if(dist < min_dist)
	// 			{
	// 				min_dist = dist;
	// 				min_t = current_t;
	// 			}
	// 		}

	// 		return min_t;
	// 	}
	}


	double Bezier::point_curve_distance(const ControlVector &ctrls, const std::array<double, 2> &p, const double t0, const double t1)
	{
		if(ctrls.size() == 6)
			return point_curve_distance_2(ctrls, p, t0, t1);
		else
			return point_curve_distance_3(ctrls, p, t0, t1);
	}

	double Bezier::point_curve_distance_2(const ControlVector &ctrls, const std::array<double, 2> &p, const double t0, const double t1)
	{
		assert(ctrls.size() == 6);

		static const std::array<int, 3> x = {{0, 2, 4}};
		static const std::array<int, 3> y = {{1, 3, 5}};

		const double cx0 = ctrls[x[0]];
		const double cy0 = ctrls[y[0]];

		const double cx1 = ctrls[x[1]];
		const double cy1 = ctrls[y[1]];

		const double cx2 = ctrls[x[2]];
		const double cy2 = ctrls[y[2]];

		const double px = p[0];
		const double py = p[1];

		double a = 4 * cx0 * cx0 + (-16 * cx1 + 8 * cx2) * cx0 + 4 * cy0 * cy0 + (-16 * cy1 + 8 * cy2) * cy0 + 16 * cx1 * cx1 - 16 * cx1 * cx2 + 4 * cx2 * cx2 + 16 * cy1 * cy1 - 16 * cy1 * cy2 + 4 * cy2 * cy2;
		double b = -12 * cx0 * cx0 + (36 * cx1 - 12 * cx2) * cx0 - 12 * cy0 * cy0 + (36 * cy1 - 12 * cy2) * cy0 - 24 * cx1 * cx1 + 12 * cx1 * cx2 - 24 * cy1 * cy1 + 12 * cy1 * cy2;
		double c = 12 * cx0 * cx0 + (-4 * px - 24 * cx1 + 4 * cx2) * cx0 + 12 * cy0 * cy0 + (-4 * py - 24 * cy1 + 4 * cy2) * cy0 + 8 * cx1 * px - 4 * px * cx2 + 8 * cy1 * py - 4 * py * cy2 + 8 * cx1 * cx1 + 8 * cy1 * cy1;
		double d = -4 * cx0 * cx0 + (4 * px + 4 * cx1) * cx0 - 4 * cy0 * cy0 + (4 * py + 4 * cy1) * cy0 - 4 * cx1 * px - 4 * cy1 * py;

		b/=a;
		c/=a;
		d/=a;

		typedef Eigen::Matrix<double, 3, 3> MatType;
		const std::function<std::array<double, 2>(const double)> gamma = [&ctrls](const double t){ return Bezier::interpolate_2(ctrls, t); };
		using std::abs;


		MatType companion; companion.setZero();
		for(int i = 0; i < 2; ++i)
			companion(i+1,i)=1;

		companion(0,2) = -d;
		companion(1,2) = -c;
		companion(2,2) = -b;

		Eigen::EigenSolver<MatType> es(companion, false);

		double min_dist = std::numeric_limits<double>::max();
		double min_t = -1;

		const auto &vals = es.eigenvalues();

		for(int i = 0; i < vals.size(); ++i)
		{
			const auto lambda = vals(i);
			const double current_t = lambda.real();

			if( abs(abs(lambda)-abs(current_t)) > 1e-8)
				continue;

			if(current_t >= t0 && current_t <= t1){
				double dist = compute_distance(p, current_t, gamma);
				if(dist < min_dist)
				{
					min_dist = dist;
					min_t = current_t;
				}
			}
		}

		if(min_t == -1){
			// std::cerr<<"Uable to find root in interval, checking endmpoints"<<std::endl;


			double dist0 = compute_distance(p, t0, gamma);
			double dist1 = compute_distance(p, t1, gamma);

			return dist0 < dist1 ? t0 : t1;
		}
		return min_t;

		// std::array<double, 4> intervals;
		// intervals[0] = -p2/3;
		// intervals[1] = (9*p0-p1*p2)/(2*p2*p2-6*p1);

		// intervals[2] = t0;
		// intervals[3] = t1;

		// const std::function<std::array<double, 2>(const double)> gamma = [&ctrls](const double t){ return Bezier::interpolate_2(ctrls, t); };
		// const std::function<std::array<double, 2>(const double)> gamma_prime = [&ctrls](const double t){ return Bezier::first_derivative_2(ctrls, t); };
		// const std::function<std::array<double, 2>(const double)> gamma_second = [&ctrls](const double t){ return Bezier::second_derivative_2(ctrls, t); };

		// return find_roots_in_interval<4>(intervals, p, gamma, gamma_prime, gamma_second, t0, t1);
	}

	double Bezier::point_curve_distance_3(const ControlVector &ctrls, const std::array<double, 2> &p, const double t0, const double t1)
	{
		assert(ctrls.size() == 8);

		static const std::array<int, 4> x = {{0, 2, 4, 6}};
		static const std::array<int, 4> y = {{1, 3, 5, 7}};

		const double cx0 = ctrls[x[0]];
		const double cy0 = ctrls[y[0]];

		const double cx1 = ctrls[x[1]];
		const double cy1 = ctrls[y[1]];

		const double cx2 = ctrls[x[2]];
		const double cy2 = ctrls[y[2]];

		const double cx3 = ctrls[x[3]];
		const double cy3 = ctrls[y[3]];

		const double px = p[0];
		const double py = p[1];

		// std::cout<<
		// "cx0 := " << cx0 <<"; "<<
		// "cy0 := " << cy0 <<"; "<<

		// "cx1 := " << cx1 <<"; "<<
		// "cy1 := " << cy1 <<"; "<<

		// "cx2 := " << cx2 <<"; "<<
		// "cy2 := " << cy2 <<"; "<<

		// "cx3 := " << cx3 <<"; "<<
		// "cy3 := " << cy3 <<"; "<<

		// "px := " << px <<"; "<<
		// "py := " << py <<"; "<<std::endl;


		double a = 6 * cx0 * cx0 + (-36 * cx1 + 36 * cx2 - 12 * cx3) * cx0 + 6 * cy0 * cy0 + (36 * cy2 - 12 * cy3 - 36 * cy1) * cy0 + 54 * cx1 * cx1 + (-108 * cx2 + 36 * cx3) * cx1 + 54 * cy1 * cy1 + (-108 * cy2 + 36 * cy3) * cy1 + 54 * cy2 * cy2 - 36 * cy2 * cy3 + 6 * cy3 * cy3 + 54 * cx2 * cx2 - 36 * cx2 * cx3 + 6 * cx3 * cx3;
		double b = -(30 * cx0 * cx0) + ((150 * cx1 - 120 * cx2 + 30 * cx3) * cx0) - (30 * cy0 * cy0) + ((-120 * cy2 + 30 * cy3 + 150 * cy1) * cy0) - (180 * cx1 * cx1) + ((270 * cx2 - 60 * cx3) * cx1) - (180 * cy1 * cy1) + ((270 * cy2 - 60 * cy3) * cy1) - (90 * cx2 * cx2) + (30 * cx2 * cx3) - 0.90e2 * cy2 * (cy2 - cy3 / 3.);
		double c = 60 * cx0 * cx0 + (-240 * cx1 + 144 * cx2 - 24 * cx3) * cx0 + 60 * cy0 * cy0 + (144 * cy2 - 24 * cy3 - 240 * cy1) * cy0 + 216 * cx1 * cx1 + (-216 * cx2 + 24 * cx3) * cx1 + 216 * cy1 * cy1 + (-216 * cy2 + 24 * cy3) * cy1 + 36 * cy2 * cy2 + 36 * cx2 * cx2;
		double d = -60 * cx0 * cx0 + (6 * px + 180 * cx1 - 72 * cx2 + 6 * cx3) * cx0 - 60 * cy0 * cy0 + (-72 * cy2 + 6 * cy3 + 6 * py + 180 * cy1) * cy0 - 108 * cx1 * cx1 + (-18 * px + 54 * cx2) * cx1 - 108 * cy1 * cy1 + (54 * cy2 - 18 * py) * cy1 + 18 * cy2 * py - 6 * cy3 * py + 18 * px * cx2 - 6 * px * cx3;
		double e = 30 * cx0 * cx0 + (-12 * px - 60 * cx1 + 12 * cx2) * cx0 + 30 * cy0 * cy0 + (12 * cy2 - 12 * py - 60 * cy1) * cy0 - 12 * cy2 * py + 24 * px * cx1 - 12 * px * cx2 + 24 * py * cy1 + 18 * cx1 * cx1 + 18 * cy1 * cy1;
		double f = -6 * cx0 * cx0 + (6 * px + 6 * cx1) * cx0 - 6 * cy0 * cy0 + (6 * py + 6 * cy1) * cy0 - 6 * px * cx1 - 6 * py * cy1;

		// const double k = std::max(std::abs(a),std::max(std::abs(b), std::max(std::abs(c), std::max(std::abs(d), std::max(std::abs(e),std::abs(f))))));
		// a/=k;
		b/=a;
		c/=a;
		d/=a;
		e/=a;
		f/=a;

		typedef Eigen::Matrix<double, 5, 5> MatType;
		const std::function<std::array<double, 2>(const double)> gamma = [&ctrls](const double t){ return Bezier::interpolate_3(ctrls, t); };
		using std::abs;


		MatType companion; companion.setZero();
		for(int i = 0; i < 4; ++i)
			companion(i+1,i)=1;

		companion(0,4) = -f;
		companion(1,4) = -e;
		companion(2,4) = -d;
		companion(3,4) = -c;
		companion(4,4) = -b;

		Eigen::EigenSolver<MatType> es(companion, false);

		double min_dist = std::numeric_limits<double>::max();
		double min_t = -1;

		const auto &vals = es.eigenvalues();

		for(int i = 0; i < vals.size(); ++i)
		{
			const auto lambda = vals(i);
			const double current_t = lambda.real();

			if( abs(abs(lambda)-abs(current_t)) > 1e-8)
				continue;

			if(current_t >= t0 && current_t <= t1){
				double dist = compute_distance(p, current_t, gamma);
				if(dist < min_dist)
				{
					min_dist = dist;
					min_t = current_t;
				}
			}
		}

		if(min_t == -1){
			// std::cerr<<"Uable to find root in interval, checking endmpoints"<<std::endl;


			double dist0 = compute_distance(p, t0, gamma);
			double dist1 = compute_distance(p, t1, gamma);

			return dist0 < dist1 ? t0 : t1;
		}
		return min_t;


		// // https://en.wikipedia.org/wiki/Quintic_function
		// const double p3=(5*a*c-2*b*b)/(5*a*a);
		// const double p2=(25*a*a*d - 15*a*b*c+4*b*b*b)/(25*a*a*a);
		// const double p1=(125*a*a*a*e - 50*a*a*b*d + 15*a*b*b*c - 3*b*b*b*b)/(125*a*a*a*a);
		// const double p0=(3125*a*a*a*a*f - 625*a*a*a*b*e + 125*a*a*b*b*d - 25*a*b*b*b*c+4*b*b*b*b*b)/(3125*a*a*a*a*a);

		// // std::cout<<"pp: "<<p0<<" "<<p1<<" "<<p2<<" "<<p3<<std::endl;


		// const double aa = 12*p3*p3*p3 + 45*p2*p2 - 40*p3*p1;
		// const double ba = 8*p3*p3*p2 + 60 *p2*p1 - 50*p3*p0;
		// const double ca = 4*p3*p3*p1 + 75*p2*p0;

		// const double ab = 10*p3;
		// const double bb = -15*p2;
		// const double cb = 4*p3*p3;

		// std::array<double, 6> intervals;
		// const double adelta = ba*ba - 4*aa*ca;
		// const double bdelta = bb*bb - 4*ab*cb;

		// // std::cout<<aa<<" "<<ab<<std::endl;
		// // std::cout<<p3<<" "<<p2<<" "<<p1<<" "<<p0<<std::endl;

		// intervals[0] = (-ba - sqrt(adelta))/(2*aa);
		// intervals[1] = (-ba + sqrt(adelta))/(2*aa);

		// intervals[2] = (-bb - sqrt(bdelta))/(2*ab);
		// intervals[3] = (-bb + sqrt(bdelta))/(2*ab);

		// // std::cout<<"i b: "<<intervals[0]<<" "<<intervals[1]<<" "<<intervals[2]<<" "<<intervals[3]<<std::endl;

		// for(int i = 0; i < 4; ++i)
		// 	intervals[i] = intervals[i] - b/(5.*a);

		// // std::cout<<"i a: "<<intervals[0]<<" "<<intervals[1]<<" "<<intervals[2]<<" "<<intervals[3]<<std::endl;

		// intervals[4] = t0;
		// intervals[5] = t1;

		// const std::function<std::array<double, 2>(const double)> gamma_prime = [&ctrls](const double t){ return Bezier::first_derivative_3(ctrls, t); };
		// const std::function<std::array<double, 2>(const double)> gamma_second = [&ctrls](const double t){ return Bezier::second_derivative_3(ctrls, t); };

		// return find_roots_in_interval<6>(intervals, p, gamma, gamma_prime, gamma_second, t0, t1);

	}

	double Bezier::inverse_interpolation(const ControlVector &ctrls, const std::array<double, 2> &p, const double start)
	{

		const std::function<std::array<double, 2>(const double)> gamma = [&ctrls](const double t){ return Bezier::interpolate(ctrls, t); };
		const std::function<std::array<double, 2>(const double)> gamma_prime = [&ctrls](const double t){ return Bezier::first_derivative(ctrls, t); };
		const std::function<std::array<double, 2>(const double)> gamma_second = [&ctrls](const double t){ return Bezier::second_derivative(ctrls, t); };

		return newton_inverse_interpolation(p, start, gamma, gamma_prime, gamma_second);
	}

	std::array<double, 2> Bezier::interpolate(const ControlVector &ctrls, const double t)
	{
		assert(ctrls.size() == 4 || ctrls.size() == 6 || ctrls.size() == 8);

		if(ctrls.size() == 4)
			return interpolate_1(ctrls, t);
		else if(ctrls.size() == 6)
			return interpolate_2(ctrls, t);
		else //if(ctrls.size() == 8)
			return interpolate_3(ctrls, t);
	}

	std::array<double, 2> Bezier::first_derivative(const ControlVector &ctrls, const double t)
	{
		assert(ctrls.size() == 4 || ctrls.size() == 6 || ctrls.size() == 8);

		if(ctrls.size() == 4)
			return first_derivative_1(ctrls, t);
		else if(ctrls.size() == 6)
			return first_derivative_2(ctrls, t);
		else //if(ctrls.size() == 8)
			return first_derivative_3(ctrls, t);
	}

	std::array<double, 2> Bezier::second_derivative(const ControlVector &ctrls, const double t)
	{
		assert(ctrls.size() == 4 || ctrls.size() == 6 || ctrls.size() == 8);

		if(ctrls.size() == 4)
			return second_derivative_1(ctrls, t);
		else if(ctrls.size() == 6)
			return second_derivative_2(ctrls, t);
		else //if(ctrls.size() == 8)
			return second_derivative_3(ctrls, t);
	}

	void Bezier::resample(const ControlVector &ctrls, const double t0, const double t1, ControlVector &new_ctrl)
	{
		assert(ctrls.size() == 4 || ctrls.size() == 6 || ctrls.size() == 8);

		if(ctrls.size() == 4)
			resample_1(ctrls, t0, t1, new_ctrl);
		else if(ctrls.size() == 6)
			resample_2(ctrls, t0, t1, new_ctrl);
		else //if(ctrls.size() == 8)
			resample_3(ctrls, t0, t1, new_ctrl);
	}

	void Bezier::resample_1(const ControlVector &ctrls, const double t0, const double t1, ControlVector &new_ctrl)
	{
		assert(ctrls.size() == 4);

		static const std::array<int, 2> x = {{0, 2}};
		static const std::array<int, 2> y = {{1, 3}};

		new_ctrl.resize(ctrls.size());
		const auto p0 = interpolate_1(ctrls, t0);
		const auto p1 = interpolate_1(ctrls, t1);
		new_ctrl[0] = p0[0];
		new_ctrl[1] = p0[1];

		new_ctrl[2] = p1[0];
		new_ctrl[3] = p1[1];
	}

	void Bezier::resample_2(const ControlVector &ctrls, const double t0, const double t1, ControlVector &new_ctrl)
	{
		assert(ctrls.size() == 6);

		static const std::array<int, 3> x = {{0, 2, 4}};
		static const std::array<int, 3> y = {{1, 3, 5}};

		const auto p0 = interpolate_2(ctrls, t0);
		const auto p1 = interpolate_2(ctrls, t1);

		const auto tt0 = first_derivative_2(ctrls, t0);

		new_ctrl.resize(ctrls.size());

		new_ctrl[0] = p0[0];
		new_ctrl[1] = p0[1];

		new_ctrl[2] = p0[0] + 1./2.*tt0[0]*(t1-t0);
		new_ctrl[3] = p0[1] + 1./2.*tt0[1]*(t1-t0);

		new_ctrl[4] = p1[0];
		new_ctrl[5] = p1[1];

#ifndef NDEBUG
		const auto tt1 = first_derivative_2(ctrls, t1);
		assert(fabs(new_ctrl[2] - (p1[0] - 1./2.*tt1[0]*(t1-t0))) < 1e-10);
		assert(fabs(new_ctrl[3] - (p1[1] - 1./2.*tt1[1]*(t1-t0))) < 1e-10);
#endif
	}

	void Bezier::resample_3(const ControlVector &ctrls, const double t0, const double t1, ControlVector &new_ctrl)
	{
		assert(ctrls.size() == 8);

		static const std::array<int, 4> x = {{0, 2, 4, 6}};
		static const std::array<int, 4> y = {{1, 3, 5, 7}};

		const auto p0 = interpolate_3(ctrls, t0);
		const auto p1 = interpolate_3(ctrls, t1);

		const auto tt0 = first_derivative_3(ctrls, t0);
		const auto tt1 = first_derivative_3(ctrls, t1);

		new_ctrl.resize(ctrls.size());

		new_ctrl[0] = p0[0];
		new_ctrl[1] = p0[1];

		new_ctrl[2] = p0[0] + 1./3.*tt0[0]*(t1-t0);
		new_ctrl[3] = p0[1] + 1./3.*tt0[1]*(t1-t0);

		new_ctrl[4] = p1[0] - 1./3.*tt1[0]*(t1-t0);
		new_ctrl[5] = p1[1] - 1./3.*tt1[1]*(t1-t0);

		new_ctrl[6] = p1[0];
		new_ctrl[7] = p1[1];
	}



	std::array<double, 2> Bezier::interpolate_1(const ControlVector &ctrls, const double t)
	{
		assert(ctrls.size() == 4);

		static const std::array<int, 2> x = {{0, 2}};
		static const std::array<int, 2> y = {{1, 3}};

		std::array<double, 2> out;
		out[0] = (1-t)*ctrls[x[0]] + t*ctrls[x[1]];
		out[1] = (1-t)*ctrls[y[0]] + t*ctrls[y[1]];

		return out;
	}

	std::array<double, 2> Bezier::first_derivative_1(const ControlVector &ctrls, const double t)
	{
		assert(ctrls.size() == 4);

		static const std::array<int, 2> x = {{0, 2}};
		static const std::array<int, 2> y = {{1, 3}};

		std::array<double, 2> out;
		out[0] = ctrls[x[1]]-ctrls[x[0]];
		out[1] = ctrls[y[1]]-ctrls[y[0]];

		return out;
	}

	std::array<double, 2> Bezier::second_derivative_1(const ControlVector &ctrls, const double t)
	{
		std::array<double, 2> out;
		out[0]=out[1]=0;

		return out;
	}


	std::array<double, 2> Bezier::interpolate_2(const ControlVector &ctrls, const double t)
	{
		assert(ctrls.size() == 6);

		static const std::array<int, 3> x = {{0, 2, 4}};
		static const std::array<int, 3> y = {{1, 3, 5}};

		std::array<double, 2> out;
		out[0] = (1-t)*(1-t)*ctrls[x[0]] + 2*(1-t)*t*ctrls[x[1]] + t*t*ctrls[x[2]];
		out[1] = (1-t)*(1-t)*ctrls[y[0]] + 2*(1-t)*t*ctrls[y[1]] + t*t*ctrls[y[2]];

		return out;
	}

	std::array<double, 2> Bezier::first_derivative_2(const ControlVector &ctrls, const double t)
	{
		assert(ctrls.size() == 6);

		static const std::array<int, 3> x = {{0, 2, 4}};
		static const std::array<int, 3> y = {{1, 3, 5}};

		std::array<double, 2> out;
		out[0] = 2*(1-t)*(ctrls[x[1]]-ctrls[x[0]]) + 2*t*(ctrls[x[2]]-ctrls[x[1]]);
		out[1] = 2*(1-t)*(ctrls[y[1]]-ctrls[y[0]]) + 2*t*(ctrls[y[2]]-ctrls[y[1]]);

		return out;
	}

	std::array<double, 2> Bezier::second_derivative_2(const ControlVector &ctrls, const double t)
	{
		assert(ctrls.size() == 6);

		static const std::array<int, 3> x = {{0, 2, 4}};
		static const std::array<int, 3> y = {{1, 3, 5}};

		std::array<double, 2> out;
		out[0] = 2*(ctrls[x[2]] - 2*ctrls[x[1]] + ctrls[x[0]]);
		out[1] = 2*(ctrls[y[2]] - 2*ctrls[y[1]] + ctrls[y[0]]);

		return out;
	}


	std::array<double, 2> Bezier::interpolate_3(const ControlVector &ctrls, const double t)
	{
		assert(ctrls.size() == 8);

		static const std::array<int, 4> x = {{0, 2, 4, 6}};
		static const std::array<int, 4> y = {{1, 3, 5, 7}};

		std::array<double, 2> out;
		out[0] = (1-t)*(1-t)*(1-t)*ctrls[x[0]] + 3*(1-t)*(1-t)*t*ctrls[x[1]] + 3*(1-t)*t*t*ctrls[x[2]] + t*t*t*ctrls[x[3]];
		out[1] = (1-t)*(1-t)*(1-t)*ctrls[y[0]] + 3*(1-t)*(1-t)*t*ctrls[y[1]] + 3*(1-t)*t*t*ctrls[y[2]] + t*t*t*ctrls[y[3]];

		return out;
	}

	std::array<double, 2> Bezier::first_derivative_3(const ControlVector &ctrls, const double t)
	{
		assert(ctrls.size() == 8);

		static const std::array<int, 4> x = {{0, 2, 4, 6}};
		static const std::array<int, 4> y = {{1, 3, 5, 7}};

		std::array<double, 2> out;
		out[0] = 3*(1-t)*(1-t)*(ctrls[x[1]] - ctrls[x[0]]) + 6*(1-t)*t*(ctrls[x[2]] - ctrls[x[1]]) + 3*t*t*(ctrls[x[3]]-ctrls[x[2]]);
		out[1] = 3*(1-t)*(1-t)*(ctrls[y[1]] - ctrls[y[0]]) + 6*(1-t)*t*(ctrls[y[2]] - ctrls[y[1]]) + 3*t*t*(ctrls[y[3]]-ctrls[y[2]]);

		return out;
	}

	std::array<double, 2> Bezier::second_derivative_3(const ControlVector &ctrls, const double t)
	{
		assert(ctrls.size() == 8);

		static const std::array<int, 4> x = {{0, 2, 4, 6}};
		static const std::array<int, 4> y = {{1, 3, 5, 7}};

		std::array<double, 2> out;
		out[0] = 6*(1-t)*(ctrls[x[2]] - 2*ctrls[x[1]] + ctrls[x[0]]) + 6*t*(ctrls[x[3]] - 2*ctrls[x[2]] + ctrls[x[1]]);
		out[1] = 6*(1-t)*(ctrls[y[2]] - 2*ctrls[y[1]] + ctrls[y[0]]) + 6*t*(ctrls[y[3]] - 2*ctrls[y[2]] + ctrls[y[1]]);

		return out;
	}


	std::array<double, 2> RationalBezier::interpolate(const ControlVector &ctrls, const ControlVector &weights, const double t)
	{
		assert(ctrls.size() == 6);
		assert(weights.size() == 3);

		const double b0=(1-t)*(1-t);
		const double b1=2*(1-t)*t;
		const double b2=t*t;

		static const std::array<int, 3> x = {{0, 2, 4}};
		static const std::array<int, 3> y = {{1, 3, 5}};

		const double c0x = ctrls[x[0]];
		const double c0y = ctrls[y[0]];

		const double c1x = ctrls[x[1]];
		const double c1y = ctrls[y[1]];

		const double c2x = ctrls[x[2]];
		const double c2y = ctrls[y[2]];

		const double denom = b0*weights[0] + b1*weights[1] + b2*weights[2];

		std::array<double, 2> out;
		out[0] = (b0*c0x*weights[0] + b1*c1x*weights[1] + b2*c2x*weights[2])/denom;
		out[1] = (b0*c0y*weights[0] + b1*c1y*weights[1] + b2*c2y*weights[2])/denom;

		return out;
	}


	std::array<double, 2> RationalBezier::first_derivative(const ControlVector &ctrls, const ControlVector &weights, const double t)
	{
		assert(ctrls.size() == 6);
		assert(weights.size() == 3);

		const double b0=(1-t)*(1-t);
		const double b1=2*(1-t)*t;
		const double b2=t*t;

		static const std::array<int, 3> x = {{0, 2, 4}};
		static const std::array<int, 3> y = {{1, 3, 5}};

		const double c0x = ctrls[x[0]];
		const double c0y = ctrls[y[0]];

		const double c1x = ctrls[x[1]];
		const double c1y = ctrls[y[1]];

		const double c2x = ctrls[x[2]];
		const double c2y = ctrls[y[2]];

		const double weights0 = weights[0];
		const double weights1 = weights[1];
		const double weights2 = weights[2];

		std::array<double, 2> out;
		out[0] = (((((-2 * c0x + 2 * c1x) * weights1 + 2 * weights2 * (c0x - c2x)) * weights0 - 2 * weights1 * weights2 * (c1x - c2x)) * t * t) + 4. * weights0 * (((c0x - c1x) * weights1) - (weights2 * (c0x - c2x)) / 2.) * t - (2 * weights0 * weights1 * (c0x - c1x))) * pow(((weights0 - 2 * weights1 + weights2) * t * t + (-2 * weights0 + 2 * weights1) * t + weights0), (-2));
		out[1] = (((((-2 * c0y + 2 * c1y) * weights1 + 2 * weights2 * (c0y - c2y)) * weights0 - 2 * weights1 * weights2 * (c1y - c2y)) * t * t) + 4. * weights0 * (((c0y - c1y) * weights1) - (weights2 * (c0y - c2y)) / 2.) * t - (2 * weights0 * weights1 * (c0y - c1y))) * pow(((weights0 - 2 * weights1 + weights2) * t * t + (-2 * weights0 + 2 * weights1) * t + weights0), (-2));

		return out;
	}

	std::array<double, 2> RationalBezier::second_derivative(const ControlVector &ctrls, const ControlVector &weights, const double t)
	{
		assert(ctrls.size() == 6);
		assert(weights.size() == 3);

		const double b0=(1-t)*(1-t);
		const double b1=2*(1-t)*t;
		const double b2=t*t;

		static const std::array<int, 3> x = {{0, 2, 4}};
		static const std::array<int, 3> y = {{1, 3, 5}};

		const double c0x = ctrls[x[0]];
		const double c0y = ctrls[y[0]];

		const double c1x = ctrls[x[1]];
		const double c1y = ctrls[y[1]];

		const double c2x = ctrls[x[2]];
		const double c2y = ctrls[y[2]];

		const double weights0 = weights[0];
		const double weights1 = weights[1];
		const double weights2 = weights[2];

		std::array<double, 2> out;
		out[0] = (4. * (((c0x - c1x) * weights1 - weights2 * (c0x - c2x)) * t + (-c0x + c1x) * weights1 - weights2 * (c0x - c2x) / 2.) * pow(t - 1., 2.) * weights0 * weights0 + (((-8. * c0x + 8. * c1x) * weights1 * weights1 + 12. * weights2 * (c0x - c2x) * weights1 - 4. * weights2 * weights2 * (c0x - c2x)) * pow(t, 3.) + 24. * (weights1 - weights2 / 2.) * ((c0x - c1x) * weights1 - weights2 * (c0x - c2x) / 2.) * t * t - 24. * (weights1 - weights2 / 2.) * weights1 * (c0x - c1x) * t + 8. * weights1 * weights1 * (c0x - c1x)) * weights0 - 8. * (weights1 - weights2 / 2.) * (c1x - c2x) * pow(t, 3.) * weights1 * weights2) * pow(pow(t - 1., 2.) * weights0 - 2. * t * ((weights1 - weights2 / 2.) * t - weights1), -3.);
		out[1] = (4. * (((c0y - c1y) * weights1 - weights2 * (c0y - c2y)) * t + (-c0y + c1y) * weights1 - weights2 * (c0y - c2y) / 2.) * pow(t - 1., 2.) * weights0 * weights0 + (((-8. * c0y + 8. * c1y) * weights1 * weights1 + 12. * weights2 * (c0y - c2y) * weights1 - 4. * weights2 * weights2 * (c0y - c2y)) * pow(t, 3.) + 24. * (weights1 - weights2 / 2.) * ((c0y - c1y) * weights1 - weights2 * (c0y - c2y) / 2.) * t * t - 24. * (weights1 - weights2 / 2.) * weights1 * (c0y - c1y) * t + 8. * weights1 * weights1 * (c0y - c1y)) * weights0 - 8. * (weights1 - weights2 / 2.) * (c1y - c2y) * pow(t, 3.) * weights1 * weights2) * pow(pow(t - 1., 2.) * weights0 - 2. * t * ((weights1 - weights2 / 2.) * t - weights1), -3.);

		return out;
	}

	double RationalBezier::inverse_interpolation(const ControlVector &ctrls, const ControlVector &weights, const std::array<double, 2> &p, const double start)
	{
		const std::function<std::array<double, 2>(const double)> gamma = [&ctrls, &weights](const double t){ return RationalBezier::interpolate(ctrls, weights, t); };
		const std::function<std::array<double, 2>(const double)> gamma_prime = [&ctrls, &weights](const double t){ return RationalBezier::first_derivative(ctrls, weights, t); };
		const std::function<std::array<double, 2>(const double)> gamma_second = [&ctrls, &weights](const double t){ return RationalBezier::second_derivative(ctrls, weights, t); };

		return newton_inverse_interpolation(p, start, gamma, gamma_prime, gamma_second);
	}

	namespace
	{
		void compute_initial_j(const std::array<Point_2f, 3>& vertices, const std::vector<Point_2f>& nodes, const int base_order, Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 15, 1> &lagr_nodes)
		{
			const int j_order = 2*(base_order-1);

			Eigen::Matrix<double, Eigen::Dynamic, 2, 0, 15, 2> ref_nodes;
			Eigen::Matrix<double, 1, 2> grad;
			Eigen::Matrix<double, 1, 2> uv;

			autogen::p_nodes_2d(j_order, ref_nodes);
			lagr_nodes.resize(ref_nodes.rows(), 1);

			Eigen::Matrix2d jacobian;

			for(int i = 0; i < lagr_nodes.size(); ++i)
			{
				jacobian.setZero();
				uv = ref_nodes.row(i);

				for(int j = 0; j < 3; ++j)
				{
					autogen::p_grad_basis_value_2d(base_order, j, uv, grad);
					jacobian.row(0) += grad * vertices[j][0];
					jacobian.row(1) += grad * vertices[j][1];
				}


				for(int j = 0; j < nodes.size(); ++j)
				{
					autogen::p_grad_basis_value_2d(base_order, 3+j, uv, grad);
					jacobian.row(0) += grad * nodes[j].x;
					jacobian.row(1) += grad * nodes[j].y;
				}

				lagr_nodes(i) = jacobian.determinant();
			}
		}

		Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 15, 1> interpolate(const int order, const Eigen::Matrix<double, Eigen::Dynamic, 2, 0, 15, 2> &ref_nodes, const Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 15, 1> &vals)
		{
			Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 15, 1> res(ref_nodes.rows(), 1);
			res.setZero();
			Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 15, 1> val;

			for(int i = 0; i < vals.rows(); ++i)
			{
				autogen::p_basis_value_2d(order, i, ref_nodes, val);

				res += vals[i]*val;
			}

			return res;
		}

		bool recursive_check(const int level, const int prev_index, const int piece, const int order, const Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 15, 1> &orig_lagr_nodes)
		{
			static const double TOL = 1e-8;
			const auto &autogen = autogen::AutoDetChecker::instance();

			// std::cout<<"level:"<<level<<" prev_index:"<<prev_index;
			// std::cout<<" piece:"<<piece<<" order:"<<order<<std::endl;

			if(level >= autogen.MAX_LEVEL){
				std::cerr<<"Max iter reached"<<std::endl;
				return false;
			}

			const int iindex = (order == 2) ? 0 : 1;

			const int current_index = autogen.get_index(prev_index, piece);

			//get the eval pts
			const Eigen::Matrix<double, Eigen::Dynamic, 2, 0, 15, 2> &loc_nodes = autogen.loc_nodes[iindex][level][current_index];

			//evaluate the function there
			const Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 15, 1> lagr_nodes = interpolate(order, loc_nodes, orig_lagr_nodes);


			// std::cout<<"-------------\n"<<lagr_nodes<<"\n----------------"<<std::endl;
			// std::cout<<"-------------\n"<<loc_nodes<<"\n----------------\n\n\n"<<std::endl;


			const int n = lagr_nodes.size();

			for(int i = 0; i < n; ++i)
			{
				//one sample point is negative
				if(lagr_nodes(i) < TOL) {
					// std::cout<<level<<" negative"<<lagr_nodes(i)<<std::endl;
					return false;
				}
			}

			//convert into bezier form
			const Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 15, 1> &bezier_nodes = autogen.L2B[iindex]*lagr_nodes;
			bool some_negative = false;
			for(int i = 0; i < n; ++i)
			{
				if(bezier_nodes(i) < TOL)
					some_negative = true;
			}

			//convex hull is above zero
			if(!some_negative)
				return true;

			for(int p = 0; p < 4; ++p)
			{
				bool check = recursive_check(level+1, current_index, p, order, orig_lagr_nodes);
				if(!check)
					return false;
			}

			return true;
		}
	}

	bool DeterminantChecker::is_positive(const std::array<Point_2f, 3>& vertices, const std::vector<Point_2f>& nodes)
	{
		//NOT Linear triangle
		assert(!nodes.empty());


		Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 15, 1> lagr_nodes;
		int base_order = -1;

		switch(nodes.size())
		{
			//quadratic
			case 3: base_order = 2; break;

			//cubic
			case 7: base_order = 3; break;

			default: assert(false);
		}

		compute_initial_j(vertices, nodes, base_order, lagr_nodes);
		const int j_order = 2*(base_order-1);
		return recursive_check(0, 0, 0, j_order, lagr_nodes);
	}
}