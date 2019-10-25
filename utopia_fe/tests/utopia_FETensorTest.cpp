#include "utopia_FETensorTest.hpp"
#include "utopia_MultiTensor.hpp"

#include "utopia_libmesh.hpp"
#include "libmesh/mesh_generation.h"

namespace utopia {

    void FETensorTest::run(Input &in)
    {
        std::cout << "[FETensorTest]" << std::endl;
        typedef utopia::LibMeshFunctionSpace FunctionSpaceT;

        auto mesh = std::make_shared<libMesh::DistributedMesh>(comm());

        const unsigned int n = 2;
        libMesh::MeshTools::Generation::build_cube(*mesh,
            n, n, n,
            0, 1,
            0, 1.,
            0, 1.,
            libMesh::TET4
        );

        auto es = std::make_shared<libMesh::EquationSystems>(*mesh);
        es->add_system<libMesh::LinearImplicitSystem>("lapl");

        auto V = FunctionSpaceT(es);

        MultiMatrixd t2;
        MultiVectord t1;
        MultiScalard t0;

        t0.resize(1);
        t0[0] = 0.5;

        t1.resize(1);
        t1[0] = values(3, 1.0);
        // t1 = t1 + t1;

        // MultiVectord t1_2 = t0 * t1;

        // MultiScalard dot_t1 = dot(t1, t1);

        // disp(dot_t1);
    }
}