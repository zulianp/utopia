#include "utopia_Testing.hpp"

#include "utopia_ui.hpp"

#include <utility>

#include "utopia_libmesh_FunctionSpace_new.hpp"
#include "utopia_libmesh_Mesh.hpp"

using namespace utopia;

void lm_create_mesh() {
    using Mesh_t = utopia::libmesh::Mesh;

    InputParameters params;
    params.set("type", "cube");
    params.set("elem_type", "TET4");
    params.set("show", true);

    Mesh_t mesh;
    mesh.read(params);
}

// template <class Arg>
// inline std::pair<std::string, Arg> param(const std::string &key, Arg &&value) {
//     return {key, std::forward<Arg>(value)};
// }

// template <class First>
// void param_set(InputParameters &params, std::pair<std::string, First> &&p) {
//     params.set(std::forward<std::string>(p.first), std::forward<First>(p.second));
// }

// void param_set(InputParameters &params, std::pair<std::string, InputParameters> &&p) {
//     params.set(std::move(p.first), std::make_shared<InputParameters>(std::move(p.second)));
// }

// template <class First, class... Args>
// void param_append(InputParameters &params, std::pair<std::string, First> &&first) {
//     params.set(std::forward<std::string>(first.first), std::forward<First>(first.second));
// }

// template <class First, class... Args>
// void param_append(InputParameters &params,
//                   std::pair<std::string, First> &&first,
//                   std::pair<std::string, Args> &&... list) {
//     params.set(std::forward<std::string>(first.first), std::forward<First>(first.second));
//     param_append(params, list...);
// }

// template <class... Args>
// inline InputParameters param_list(std::pair<std::string, Args> &&... list) {
//     InputParameters ret;
//     param_append(ret, std::forward<Args>(list)...);
//     return ret;
// }

// template <class Arg>
// inline InputParameters param_list(std::pair<std::string, Arg> &&p) {
//     InputParameters ret;
//     param_set(ret, std::forward<Arg>(p));
//     return ret;
// }

void lm_create_function_space() {
    using FunctionSpace_t = utopia::libmesh::FunctionSpace;

    InputParameters params;
    params.set("system_type", "linear_implicit");
    // params.set("show", true);

    // auto params = param_list(param("system_type", "linear_implicit"),
    //                          param("show", true),
    //                          param("mesh", param_list(param("type", "sphere"))));

    FunctionSpace_t space;
    space.read(params);

    auto s0 = space[0];

    std::cout << s0.n_local_dofs() << " " << s0.n_dofs() << std::endl;

    UVector v;
    USparseMatrix m;

    space.create_vector(v);
    space.create_matrix(m);
}

void lm() {
    UTOPIA_RUN_TEST(lm_create_mesh);
    UTOPIA_RUN_TEST(lm_create_function_space);
}

UTOPIA_REGISTER_TEST_FUNCTION(lm);
