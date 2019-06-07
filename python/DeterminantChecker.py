import os
import numpy as np
from math import factorial
from sympy import *

import matplotlib.pyplot as plt

p2_nodes=np.array([
    [0, 0],
    [1, 0],
    [0, 1],
    [1/2, 0],
    [1/2, 1/2],
    [0, 1/2]
])

bezier_indices_2 = np.array([
    [0, 0],
    [1, 0], [1, 1], [0, 1],
    [2, 0], [0, 2]
])

p4_nodes=np.array([
    [0, 0],
    [1, 0],
    [0, 1],
    [1/4, 0],
    [1/2, 0],
    [3/4, 0],
    [3/4, 1/4],
    [1/2, 1/2],
    [1/4, 3/4],
    [0, 3/4],
    [0, 1/2],
    [0, 1/4],
    [1/4, 1/4],
    [1/4, 1/2],
    [1/2, 1/4]
])


bezier_indices_4 = np.array([
    [0, 0],
    [1, 0], [1, 1], [0, 1],
    [2, 0], [2, 1], [2, 2], [0, 2], [1, 2],
    [3, 0], [3, 1], [0, 3], [1, 3],
    [4, 0], [0, 4]
])


def bezier(p, i, j, xi, eta):
    k = (p-i-j)
    return factorial(p) / factorial(i) / factorial(j) / factorial(k) * xi**i * eta**j * (1 - xi - eta)**k


def create_matrix(deg, nodes, bezier_indices):
    n = len(bezier_indices)
    assert(n == len(nodes))

    B2L = zeros(n)

    for i in range(n):
        xi = nodes[i,0]
        eta = nodes[i,1]

        for j in range(n):
            B2L[i,j] = bezier(deg, bezier_indices[j,0], bezier_indices[j,1], xi, eta)

    return [B2L, B2L.inv()]


def subdivide_one(p1, p2, p3, nodes):
    n = len(nodes)
    # print(n)
    new_nodes = np.zeros([n, 2])

    for i in range(n):
        b2 = nodes[i, 0]
        b3 = nodes[i, 1]
        b1 = 1 - b2 - b3
        assert(b3 >= 0)

        new_p = p1*b1 + p2*b2 + p3*b3
        # print(new_p)
        new_nodes[i, :] = new_p

    return new_nodes


def subdivide(nodes, main_nodes):
    new_nodes_1 = subdivide_one(
        nodes[0,:],
        (nodes[0,:]+nodes[1,:])/2,
        (nodes[0,:]+nodes[2,:])/2, main_nodes)

    new_nodes_2 = subdivide_one(
        (nodes[0,:]+nodes[1,:])/2,
        nodes[1,:],
        (nodes[1,:]+nodes[2,:])/2, main_nodes)

    new_nodes_3 = subdivide_one(
        (nodes[0,:]+nodes[2,:])/2,
        (nodes[0,:]+nodes[1,:])/2,
        (nodes[1,:]+nodes[2,:])/2, main_nodes)

    new_nodes_4 = subdivide_one(
        (nodes[2,:]+nodes[0,:])/2,
        (nodes[2,:]+nodes[1,:])/2,
        nodes[2,:], main_nodes)

    return [new_nodes_1, new_nodes_2, new_nodes_3, new_nodes_4]


def get_index(parent, q):
    return parent*4 + q


def create_level(iindex, level, parent_index, q, deg, nodes, main_nodes, bezier_indices, max_levels):
    if level > max_levels:
        return ""

    # [B2L, L2B] = create_matrix(deg, main_nodes, bezier_indices)
    # size = L2B.rows
    # # print(size)

    current_index = get_index(parent_index, q)

    res = ""
    if current_index == 0:
        # res += "L2B[{}].resize({});\n".format(level, 4**level)
        res += "loc_nodes[{}][{}].resize({});\n".format(iindex, level, 4**level)

        print("level " + str(level))


    # current_mat = "L2B[{}][{}]".format(level, current_index)
    # res += current_mat + ".resize({}, {});\n".format(size, size)
    # res += current_mat + "<<\n"

    # for i in range(size):
    #     for j in range(size):
    #         res += "{},".format(L2B[i, j])
    #     res+="\n"

    # res = res[:-2]
    # res += ";\n\n"

    size = nodes.shape[0]

    current_mat = "loc_nodes[{}][{}][{}]".format(iindex, level, current_index)
    # res += current_mat + ".resize({}, 2);\n".format(size)
    # res += current_mat + "<<\n"
    # for i in range(size):
    #     for j in range(2):
    #         res += "{},".format(nodes[i, j])
    #     res+="\n"
    # res = res[:-2]
    # res += ";\n\n"


    current_mat_arr = current_mat.replace('[', '_').replace(']', '')
    res += 'static double {}[] = '.format(current_mat_arr)
    res += "{"

    for i in range(size):
        for j in range(2):
            res += "{},".format(nodes[i, j])
        res+="\n"

    res = res[:-2]
    res += "};\n"
    res += "{} = Eigen::Map<Eigen::Matrix<double, 2, Eigen::Dynamic, 0, 2, 15>>({},2,{}).transpose();\n\n".format(current_mat, current_mat_arr, size)



    new_nodes = subdivide(nodes, main_nodes)

    q = 0
    for nn in new_nodes:
        res += create_level(iindex, level+1, current_index, q, deg, nn, main_nodes, bezier_indices, max_levels)
        q += 1

    return res


def compute_one(iindex, order, max_levels, nodes, bezier_indices):
    matrices = ""

    [B2L, L2B] = create_matrix(order, nodes, bezier_indices)
    size = L2B.rows
    current_mat = "L2B[{}]".format(iindex)
    # matrices += current_mat + ".resize({}, {});\n".format(size, size)
    # matrices += current_mat + "<<\n"

    current_mat_arr = current_mat.replace('[', '_').replace(']', '')
    matrices += 'static double {}[] = '.format(current_mat_arr)
    matrices += "{"

    for i in range(size):
        for j in range(size):
            matrices += "{},".format(L2B[i, j])
        matrices+="\n"

    matrices = matrices[:-2]
    matrices += "};\n"
    matrices += "{} = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 15, 15>>({},{},{}).transpose();\n\n".format(current_mat, current_mat_arr, size, size)

    matrices += create_level(iindex, 0, 0, 0, order, nodes, nodes, bezier_indices, max_levels)

    return matrices


def main():
    path = "../"
    max_levels = 5

    hpp = "#include <Eigen/Dense>\n"
    hpp += "#include<vector>\n"
    hpp += "#include<array>\n\n"
    hpp += "namespace tri_wild{\nnamespace autogen{\n"
    hpp += "class AutoDetChecker{\npublic:\n"
    hpp += "std::array<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 15, 15>,2> L2B;\n"
    hpp += "std::array<std::array<std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2, 0, 15, 2>>, {}>, 2> loc_nodes;\n\n".format(max_levels+1)
    hpp += "static const int MAX_LEVEL={};\n\n".format(max_levels+1)
    hpp += "static const AutoDetChecker &instance();\n\n"

    cpp = "#include \"auto_det_checker.hpp\"\n\n"
    cpp += "namespace tri_wild{\nnamespace autogen{\n"
    cpp += "const AutoDetChecker& AutoDetChecker::instance(){static AutoDetChecker tmp; return tmp; }\n\n"

    hpp += "int get_index(const int parent, const int q) const;\n"
    cpp += "int AutoDetChecker::get_index(const int parent, const int q) const { return parent*4 + q; }\n"

    hpp += "private:\nAutoDetChecker();\n"
    cpp += "AutoDetChecker::AutoDetChecker(){\n"

    order2 = compute_one(0, 2, max_levels, p2_nodes, bezier_indices_2)
    order4 = compute_one(1, 4, max_levels, p4_nodes, bezier_indices_4)

    cpp += order2 + "\n\n" + order4 + "}"




    cpp += "}\n}\n"
    hpp += "};\n}\n}\n"



    print("saving...")
    with open(os.path.join(path, "auto_det_checker.cpp"), "w") as file:
        file.write(cpp)

    with open(os.path.join(path, "auto_det_checker.hpp"), "w") as file:
        file.write(hpp)


    # [B2L4, L2B4] = create_matrix(4, p4_nodes, bezier_indices_4)

    # test = 3
    # new_nodes = subdivide(p4_nodes, p4_nodes)
    # new_nodes_1 = subdivide(new_nodes[test], p4_nodes)

    # plt.plot(p4_nodes[:,0], p4_nodes[:,1], 'k*')
    # plt.plot(new_nodes[test][:,0], new_nodes[test][:,1], '*')

    # for n in new_nodes_1:
    #     plt.plot(n[:,0], n[:,1], '.')


    # for i in range(len(new_nodes_1)):
    #     plt.text(new_nodes_1[i,0], new_nodes_1[i,1], str(i+1))

    # for i in range(len(p4_nodes)):
    #     plt.text(p4_nodes[i,0], p4_nodes[i,1], str(i+1))
    plt.show()



if __name__ == '__main__':
    main()
