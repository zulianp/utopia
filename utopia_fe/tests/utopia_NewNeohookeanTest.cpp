#include "utopia_NewNeohookeanTest.hpp"

#include "utopia_ui.hpp"
#include "utopia_SymbolicFunction.hpp"

#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIScalarSampler.hpp"
#include "utopia_MeshTransferOperator.hpp"
#include "utopia_Flow.hpp"

#include "utopia_FEEval_Local.hpp"
#include "utopia_libmesh_AssembleLocal.hpp"

#include "utopia_libmesh.hpp"
#include "libmesh/mesh_generation.h"

#include "libmesh/mesh_refinement.h"
#include "libmesh/boundary_mesh.h"

namespace utopia {

    template<class Vector, class Traits, int Backend>
    class Eval<Interpolate<
                Vector,
                TrialFunction<
                    FiniteElement<ProductFunctionSpace<LibMeshFunctionSpace>, LIBMESH_TAG> 
                    >
                >, Traits, Backend> 
    {
    public:

        using Expr = utopia::Interpolate<
                Vector,
                TrialFunction<
                    FiniteElement<ProductFunctionSpace<LibMeshFunctionSpace>, LIBMESH_TAG> 
                    >
                >;

        template<class Out>
        inline static void apply(const Expr &expr, Out &out)
        {
            assert(false && "IMPLEMENT ME");
        }

    };


    template<class Vector, class Traits, int Backend>
    class Eval<Gradient<Interpolate<
                Vector,
                TrialFunction<
                    FiniteElement<ProductFunctionSpace<LibMeshFunctionSpace>, LIBMESH_TAG> 
                    >
                >>, Traits, Backend> 
    {
    public:

        using Expr = utopia::Gradient<Interpolate<
                Vector,
                TrialFunction<
                    FiniteElement<ProductFunctionSpace<LibMeshFunctionSpace>, LIBMESH_TAG> 
                    >
                >>;

        template<class Out>
        inline static void apply(const Expr &expr, Out &out)
        {
            assert(false && "IMPLEMENT ME");
        }

    };

    

    void NewNeohookeanTest::run(Input &in)
    {
        std::cout << "[NewNeohookeanTest]" << std::endl;
        using Space = utopia::LibMeshFunctionSpace;
        using FE = utopia::FiniteElement<Space>;
        using TFE = utopia::FiniteElement<ProductFunctionSpace<Space>>;

        UIMesh<libMesh::DistributedMesh> mesh(comm());
        UIFunctionSpace<Space> space(make_ref(mesh));
        in.get("mesh", mesh);
        in.get("space", space);

        UIForcingFunction<Space, UVector> forcing_function(space.subspace(0));
        in.get("forcing-function", forcing_function);


        auto &V =space.space();

        //FIXME use ghost
        UVector x = local_zeros(V[0].dof_map().n_local_dofs());

        MultiMatrixd uk;
        
        USparseMatrix A;
        assemble(
            V,
            inner(grad(trial(V)), grad(trial(V))) * dX,
            A,
            [&](TFE &element, USerialMatrix &mat)
            {
                    //Symbolic
                    auto u = trial(element);
                    uk = grad(interpolate(x, u));

                    disp(uk);
                    exit(0);

                    //Symbolic to Numeric conversion
                    // f = u;

                    // f *= 0.0;
                    // dx  = measure(element);
                    // vec = form1(f, dx);
            }
        );

    }

}