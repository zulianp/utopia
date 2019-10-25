#include "utopia_FETensorTest.hpp"
#include "utopia_MultiTensor.hpp"
#include "utopia_FormTensor.hpp"
#include "utopia_FiniteElement.hpp"
#include "utopia_libmesh_FiniteElement.hpp"

#include "utopia_libmesh.hpp"
#include "libmesh/mesh_generation.h"

namespace utopia {
    template<class Traits, int Backend>
    class Eval<
        Gradient<
            TrialFunction<
                FiniteElement<LibMeshFunctionSpace>
                >
            >, Traits, Backend> {
    public:

        using Expr = Gradient<
            TrialFunction<
                FiniteElement<LibMeshFunctionSpace>
                >
            >;

        template<class T>
        inline static void apply(const Expr &expr, FormTensor<T, 1> &result)
        {
            std::cout << expr.get_class() << std::endl;
        }   

    };

    // template<class T>
    // void assemble(
    //     const FunctionSpace<T> &V, )

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

        FormVectord g;
        FiniteElement<FunctionSpaceT> element(V);

        //Which basis-function
        auto u = trial(element);

        for(auto e_it = mesh->active_local_elements_begin(); e_it != mesh->active_local_elements_end(); ++e_it) {
            //Change element
            element.set((*e_it)->id());
            //What we need
            element.init(grad(u) + u);

            g = grad(u);
        }
    }
}