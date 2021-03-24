#include "utopia_GeometryTest.hpp"
#include "utopia_libmesh.hpp"

#include "utopia_NormalTangentialCoordinateSystem.hpp"

#include <libmesh/const_function.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/petsc_matrix.h>
#include <libmesh/petsc_vector.h>
#include "libmesh/nemesis_io.h"

#include <memory>

namespace utopia {

    void GeometryTest::run(Input &in) {
        using namespace libMesh;
        using namespace std;

        const auto elem_order = FIRST;

        auto mesh = make_shared<libMesh::Mesh>(this->comm());
        mesh->read("../data/stent2_101.e");
        auto dim = mesh->mesh_dimension();

        auto equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);
        auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("geo-test");

        ////////////////////////////////////////////

        auto Vx = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "n_x");
        auto Vy = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "n_y");
        auto V = Vx * Vy;

        if (dim == 3) {
            V *= LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "n_z");
        }

        Vx.initialize();

        UVector is_normal_component;
        UVector normals;
        USparseMatrix mat;
        assemble_normal_tangential_transformation(*mesh, Vx.dof_map(), {101}, is_normal_component, normals, mat);
        // mat.implementation().set_name("t");
        // write("O.m", mat);

        // normals.implementation().set_name("n");
        // write("vn.m", normals);

        UVector t_normals = mat * normals;

        libMesh::Nemesis_IO io(*mesh);
        convert(t_normals, *sys.solution);
        sys.solution->close();
        io.write_equation_systems("geo-test.e", *equation_systems);
    }
}  // namespace utopia
