#include "../extern/CLI11.hpp"
#include "triwild/meshio.hpp"

#include <igl/Timer.h>
#include <igl/writeSTL.h>
#include "triwild/optimization.h"
#include "triwild/feature.h"
#include "triwild/triangulation.h"

using namespace triwild;

//#include "triwild/do_triwild.h"

int main(int argc, char* argv[]) {
    CLI::App app{"Robust Triangulation"};
    app.add_option("--input", args.input, "Input segments in .obj format.")->required();
    app.add_option("--output", args.output, "Output path.");
    app.add_option("--postfix", args.postfix, "Add postfix into outputs' file name.");
    app.add_option("--feature-input", args.feature_input, "Input feature json file.");
//    app.add_option("--quality-output", args.quality_output, "Quality_output");
    app.add_option("--stop-quality", args.stop_quality, "Specify max AMIPS energy for stopping mesh optimization.");
    app.add_option("--max-its", args.max_its, "Max number of mesh optimization iterations.");
    app.add_option("--stage", args.stage, "Specify envelope stage");
    app.add_option("--envelope-r", args.epsilon, "relative envelope epsilon_r. Absolute epsilonn = epsilon_r * diagonal_of_bbox");
    app.add_option("--feature-envelope-r", args.feature_epsilon, "Relative feature envelope mu_r. Absolute mu = mu_r * diagonal_of_bbox");
    app.add_option("--target-edge-length", args.target_edge_len, "Absolute target edge length l.");//absolute value fore edge length
    app.add_option("--target-edge-length-r", args.edge_length_r, "Relative target edge length l_r. Absolute l = l_r * diagonal_of_bbox");//relative value for edge length
//    app.add_option("--dd", args.i_dd, "");
    app.add_option("--log-file", args.log_file, "Output a log file.");
    app.add_option("--min-angle", args.flat_feature_angle, "Desired minimal angle.");

    app.add_flag("--mute-log", args.mute_log, "Mute prints.");
//    app.add_flag("--enable-debug-mesh", args.enable_debug_mesh, "Enable output debug mesh.");

    bool cut_outside = false;
    app.add_flag("--cut-outside", cut_outside, "Remove \"outside part\".");

    bool diffusion_curves = false;
//    app.add_flag("--diffusion-curves", diffusion_curves, "Export for diffusion curves");

    bool skip_eps = false;
    app.add_flag("--skip-eps", skip_eps, "Skip saving eps.");

    std::string hole_file;
    app.add_option("--cut-holes", hole_file, "Input a .xyz file for specifying points inside holes you want to remove.");

    app.add_flag("--output-linear-mesh", args.output_linear, "Output linear mesh for curved pipeline.");

    try {
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    GEO::initialize();

    igl::Timer main_timer;
    double load_mesh_time = 0;
    double preprocess_time = 0;
    double bsp_time = 0;
    double optimization_time = 0;
    double curving_time = 0;
    double cut_and_hole_time = 0;


    if (args.mute_log) {
        std::streambuf *orig_buf = cout.rdbuf();
        cout.rdbuf(NULL);
    }

    if (args.output == "")
        args.output = args.input + "_" + args.postfix;

    igl::Timer igl_timer;
    double t;
    ///////////////////
    cout << "Loading and preprocessing ..." << endl;
    igl_timer.start();

    Eigen::MatrixXd V;
    std::vector<std::array<int, 2>> edges;
    triangulation::load_input(args.input, V, edges);

//    Eigen::MatrixXi E(edges.size(), 2);
//    for(int i=0;i<edges.size();i++){
//        for(int j=0;j<2;j++)
//            E(i, j) = edges[i][j];
//    }
//
//    using json = nlohmann::json;
//    json feature_info = json({});
//    if (args.feature_input != "") {
//        std::ifstream file(args.feature_input);
//        if (file.is_open())
//            file >> feature_info;
//        file.close();
//    }
//
//    Eigen::MatrixXd V_out;
//    Eigen::MatrixXi F_out;
//    Eigen::MatrixXd nodes;
//    std::vector<std::vector<int>> F_nodes;
//    do_triwild(V, E, feature_info, V_out, F_out, nodes, F_nodes);
//
//    return 0;

    int input_v_cnt = V.rows();
    int input_e_cnt = edges.size();
    load_mesh_time = igl_timer.getElapsedTime();
    igl_timer.start();

    if (feature::init(args.feature_input))
        args.is_preserving_feature = true;

    if (args.stop_quality < 0)
        args.stop_quality = args.is_preserving_feature ? 20 : 10;
    if (args.epsilon < 0)
//        args.epsilon = args.is_preserving_feature ? 5e-3 : 1e-3;
        args.epsilon = args.is_preserving_feature ? 2e-3 : 1e-3;

    double line_width = 0.3 * args.diagonal_len;
    double f_line_width = 0.5 * args.diagonal_len;
    double s_f_line_width = 0.4 * args.diagonal_len;
//        double draw_points = true;
    double draw_points = false;

    double point_size = 0.0005 * args.diagonal_len;
    double f_point_size = 0.001 * args.diagonal_len;
    double s_f_point_size = 0.003 * args.diagonal_len;

    std::string line_col = "0 0 0";
    std::string f_line_col = "0.9 0.3 0.2";
    std::string s_f_line_col = "0.1 0.7 0.6";

    std::string point_col = "0 0 0";
    std::string f_point_col = "0.9 0.3 0.2";
    std::string s_f_point_col = "0.1 0.7 0.6";

    GEO::Mesh b_mesh;
    triangulation::preprocessing(V, edges, b_mesh);
    GEO::MeshFacetsAABB b_tree(b_mesh);
    t = igl_timer.getElapsedTime();
    preprocess_time = t;
    cout << "Loaded and preprocessed." << endl;
    cout << "time = " << t << "s" << endl << endl;

    ///////////////////
    cout << "BSP subdivision..." << endl;
    igl_timer.start();
    MeshData mesh;
    std::vector<std::vector<int>> tag_boundary_es;
    triangulation::BSP_subdivision(V, edges, mesh, tag_boundary_es);
    t = igl_timer.getElapsedTime();
    bsp_time = t;
    cout << "BSP subdivision done." << endl;
    cout << "time = " << t << "s" << endl << endl;

    ///////////////////
    cout << "Mesh optimization..." << endl;
    igl_timer.start();
    optimization::init(V, edges, tag_boundary_es, mesh, b_tree);

    optimization::refine(mesh, b_tree, std::array<int, 4>({1, 1, 1, 1}));
    t = igl_timer.getElapsedTime();
    optimization_time = t;
    cout << "Mesh optimization done." << endl;
    cout << "time = " << t << "s" << endl;

//    if (!skip_eps) {
//         export_eps(mesh,
//             line_width, line_col, point_size, point_col,
//             f_line_width, f_line_col, f_point_size, f_point_col,
//             s_f_line_width, s_f_line_col, s_f_point_size, s_f_point_col,
//             draw_points, args.output + "_lin.eps");
//     }

    ///////////////////
    if (!args.is_preserving_feature) {
        optimization::output_mesh(mesh);

        if (args.log_file != "") {
            std::ofstream file(args.log_file);
            if (!file.fail()) {
                file << "load_mesh_time " << load_mesh_time << "\n";
                file << "preprocess_time " << preprocess_time << "\n";
                file << "bsp_time " << bsp_time << "\n";
                file << "optimization_time " << optimization_time << "\n";
                file << "curving_time " << curving_time << "\n";
                file << "cut_and_hole_time " << cut_and_hole_time << "\n";
                file << "input_v_cnt " << input_v_cnt  << "\n";
                file << "input_e_cnt " << input_e_cnt  << "\n";
                optimization::output_stats(mesh, file);
            }
            file.close();
        }

        return 0;
    }

    cout << "Curving..." << endl;
    igl_timer.start();
    feature::curving(mesh, b_tree);
    t = igl_timer.getElapsedTime();
    curving_time = t;
    cout << "Curving done." << endl;
    cout << "time = " << t << "s" << endl;

    igl_timer.start();
    if (cut_outside)
        optimization::erase_outside(mesh);
    if (hole_file != "")
        optimization::erase_holes(mesh, hole_file);
    cut_and_hole_time = igl_timer.getElapsedTime();

    optimization::output_mesh(mesh);

    if (diffusion_curves)
        write_msh_DiffusionCurve(mesh, args.output + "DC");

    if (!skip_eps) {
       export_eps(mesh,
                   line_width, line_col, point_size, point_col,
                   f_line_width, f_line_col, f_point_size, f_point_col,
                   s_f_line_width, s_f_line_col, s_f_point_size, s_f_point_col,
                   draw_points, args.output + ".eps");

        // const double N = 10;
        // for(int i = 0; i <= N; ++i)
        // {
        // const double t = i/N;
        //     const double t = i/N;
        //     export_eps(mesh,
        //      line_width, line_col, point_size, point_col,
        //      f_line_width, f_line_col, f_point_size, f_point_col,
        //      s_f_line_width, s_f_line_col, s_f_point_size, s_f_point_col,
        //      draw_points, args.output + "_" + std::to_string(i) + ".eps", t);
        // }
    }

    if (args.log_file != "") {
        std::ofstream file(args.log_file);
        if (!file.fail()) {
            file << "load_mesh_time " << load_mesh_time << "\n";
            file << "preprocess_time " << preprocess_time << "\n";
            file << "bsp_time " << bsp_time << "\n";
            file << "optimization_time " << optimization_time << "\n";
            file << "curving_time " << curving_time << "\n";
            file << "cut_and_hole_time " << cut_and_hole_time << "\n";
            file << "input_v_cnt " << input_v_cnt  << "\n";
            file << "input_e_cnt " << input_e_cnt  << "\n";
            optimization::output_stats(mesh, file);
            feature::output_stats(mesh, file);
        }
        file.close();
    }

    return 0;
}