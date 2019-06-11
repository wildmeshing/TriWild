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

#include "reference_triangle.h"
#include <igl/list_to_matrix.h>


namespace triwild{
    Eigen::MatrixXd get_reference_triangle_vertices() {
        static bool init = false;
        static Eigen::MatrixXd V;
        if(!init) {
            std::fstream file;
//            file.open("../src/triwild/reference_triangle_vertices.txt");
            file.open(_REF_VS);
//            std::cout<<_REF_VS<<std::endl;

            if (!file.good()) {
                std::cout << "Failed to open file" << std::endl;
                file.close();
                return V;
            }

            std::string s;
            std::vector<std::vector<double>> matrix;
            while (getline(file, s)) {
                std::stringstream input(s);
                double temp;
                matrix.emplace_back();

                std::vector<double> &currentLine = matrix.back();

                while (input >> temp)
                    currentLine.push_back(temp);
            }
            file.close();

            if (!igl::list_to_matrix(matrix, V)) {
                std::cout << matrix.size() << endl;
                std::cout << "list to matrix error" << std::endl;
            }

            init = false;
        }

//        std::cout << V.rows() << std::endl;
//        std::cout << V.row(0) << std::endl;
//        std::cout << V.row(V.rows() - 1) << std::endl;
//
//        std::cout<<"pausing"<<std::endl;
//        char c;
//        std::cin>>c;

        return V;
    }

    Eigen::MatrixXi get_reference_triangle_faces(){
        static bool init = false;
        static Eigen::MatrixXi F;
        if(!init) {
            std::fstream file;
//            file.open("../src/triwild/reference_triangle_faces.txt");
            file.open(_REF_FS);
//            std::cout<<_REF_VS<<std::endl;

            if (!file.good()) {
                std::cout << "Failed to open file" << std::endl;
                file.close();
                return F;
            }

            std::string s;
            std::vector<std::vector<double>> matrix;
            while (getline(file, s)) {
                std::stringstream input(s);
                double temp;
                matrix.emplace_back();

                std::vector<double> &currentLine = matrix.back();

                while (input >> temp)
                    currentLine.push_back(temp);
            }
            file.close();

            if (!igl::list_to_matrix(matrix, F))
                std::cout << "list to matrix error" << std::endl;
            init = false;
        }

//        std::cout << F.rows() << std::endl;
//        std::cout << F.row(0) << std::endl;
//        std::cout << F.row(F.rows() - 1) << std::endl;
//
//        std::cout<<"pausing"<<std::endl;
//        char c;
//        std::cin>>c;

        return F;
    }
}