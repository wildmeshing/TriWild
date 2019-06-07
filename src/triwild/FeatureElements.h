// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef TRIWILD_FEATUREELEMENTS_H
#define TRIWILD_FEATUREELEMENTS_H

#include "Curves.h"

#include <array>
#include <string>
#include <vector>
#include <iostream>
#include <memory>

#include <Eigen/Dense>

#include "Point_2f.h"

namespace triwild {
    namespace feature {
        class FeatureElement {
        public:
            std::vector<int> v_ids;
            std::vector<double> paras;
            std::string type;
            int degree = 1;
            int curve_id;
            std::array<bool, 2> is_inflection = {{false, false}};

            virtual void print_info() const  = 0;

            virtual Point_2f eval(double t)  const = 0;
            virtual Point_2f eval_first_derivative(double t)  const = 0;
            virtual Point_2f eval_second_derivative(double t)  const = 0;
            virtual double inv_eval(Point_2f& p, double t, double t0, double t1) const = 0;
            // virtual double distance(Point_2f& p, double t = 0) const = 0;

            virtual std::shared_ptr<FeatureElement> simplify(const Point_2f& start, const Point_2f& end) const{ return NULL; }

            virtual std::vector<double> inflection_points(const double t0, const double t1) const { return std::vector<double>(); }
            virtual double get_cut_t(const double t0, const double t1, bool to_flip = false) const = 0;

            virtual std::string to_maple() const { return ""; }
            virtual std::string to_eps() const;

            virtual double how_curve(double t1, double t2) const = 0;
            virtual bool is_point() const = 0;
            virtual bool is_too_short(const Eigen::MatrixXd& V, double eps) const = 0;

            virtual ~FeatureElement(){}
            FeatureElement(const FeatureElement& other){
                v_ids = other.v_ids;
                paras = other.paras;
                type = other.type;
                degree = other.degree;
                curve_id = other.curve_id;
            }
            FeatureElement(const std::vector<int> &v_ids_, const std::vector<double> &paras_,
                    std::string type_, int degree_, int curve_id_):
                    v_ids(v_ids_), paras(paras_), type(type_), degree(degree_), curve_id(curve_id_){}

            void merge_after(const FeatureElement &other);
        };

        class Line_Feature: public FeatureElement {
        public:
            Line_Feature(const std::vector<int> &v_ids_, const std::vector<double> &paras_,
                    const Point_2f &start_, const Point_2f &end_, int curve_id_) :
            FeatureElement(v_ids_, paras_, "Line", 1, curve_id_),
            start(start_), end(end_) {}

            Line_Feature(const Line_Feature& other): FeatureElement(other){
                start = other.start;
                end = other.end;
            }

            void print_info() const override { std::cout<<"start="<<start<<" end="<<end<<std::endl; }

            std::string to_eps() const override;

            Point_2f eval(double t) const override;
            Point_2f eval_first_derivative(double t) const override;
            Point_2f eval_second_derivative(double t) const override;
            double inv_eval(Point_2f& p, double t, double t0, double t1) const override;
            // double distance(Point_2f& p, double t = 0) const override;

            double get_cut_t(const double t0, const double t1, bool to_flip = false) const override { return (t0 + t1)/2; }
            double how_curve(double t1, double t2) const override { return 0; }
            bool is_point() const override;
            bool is_too_short(const Eigen::MatrixXd& V, double eps) const override;

        protected:
            Point_2f start;
            Point_2f end;
        };

        class RationalBezierCurve_Feature: public FeatureElement {
        public:
            RationalBezierCurve_Feature(const std::vector<int> &v_ids_, const std::vector<double> &paras_,
                        const ControlVector &poles_, const ControlVector &weights_, int curve_id_) :
                    FeatureElement(v_ids_, paras_, "RationalBezierCurve", 3, curve_id_),
                    poles(poles_), weights(weights_) {}

            RationalBezierCurve_Feature(const RationalBezierCurve_Feature& other): FeatureElement(other){
                poles = other.poles;
                weights = other.weights;
            }

            Point_2f eval(double t) const override;
            Point_2f eval_first_derivative(double t)  const override;
            Point_2f eval_second_derivative(double t)  const override;
            double inv_eval(Point_2f& p, double t, double t0, double t1) const override;
            // double distance(Point_2f& p, double t = 0) const override;

            double get_cut_t(const double t0, const double t1, bool to_flip = false) const override;
            double how_curve(double t1, double t2) const override;
            bool is_point() const override;
            bool is_too_short(const Eigen::MatrixXd& V, double eps) const override;

            void print_info() const override;

        protected:
            ControlVector poles;
            ControlVector weights;
        };

        class BezierCurve_Feature: public FeatureElement {
        public:
            BezierCurve_Feature(const std::vector<int> &v_ids_, const std::vector<double> &paras_,
                            int degree_, const ControlVector &poles_, int curve_id_) :
                    FeatureElement(v_ids_, paras_, "BezierCurve", 3, curve_id_),
                    poles(poles_) {}

            BezierCurve_Feature(const BezierCurve_Feature& other): FeatureElement(other){
                poles = other.poles;
            }

            Point_2f eval(double t) const override;
            Point_2f eval_first_derivative(double t)  const override;
            Point_2f eval_second_derivative(double t)  const override;
            double inv_eval(Point_2f& p, double t, double t0, double t1) const override;
            // double distance(Point_2f& p, double t = 0) const override;

            std::vector<double> inflection_points(const double t0, const double t1) const  override;

            double get_cut_t(const double t0, const double t1, bool to_flip = false) const override;
            double how_curve(double t1, double t2) const override;
            bool is_point() const override;
            bool is_too_short(const Eigen::MatrixXd& V, double eps) const override;

            std::string to_maple() const override;
            std::string to_eps() const override;
            void print_info() const override;

            // void compute_roots_intervals();
            std::shared_ptr<FeatureElement> simplify(const Point_2f& start, const Point_2f& end) const override;

        protected:
            ControlVector poles;
        };

    }
}


#endif //TRIWILD_FEATUREELEMENTS_H
