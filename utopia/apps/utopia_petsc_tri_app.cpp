#include "utopia_Base.hpp"

// include edsl components
#include "utopia_AppRunner.hpp"
#include "utopia_AssemblyView.hpp"
#include "utopia_BratuFE.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_Core.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_LaplacianView.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_MassMatrixView.hpp"
#include "utopia_PetscDM.hpp"
#include "utopia_PhaseField.hpp"
#include "utopia_PoissonFE.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_SampleView.hpp"
#include "utopia_Tri3.hpp"
#include "utopia_TrivialPreconditioners.hpp"
#include "utopia_petsc.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_dma_FunctionSpace.hpp"
#include "utopia_petsc_impl.hpp"

#include "utopia_LinearElasticityFE.hpp"

#include "utopia_app_utils.hpp"

#include "utopia_petsc_NeumannBoundaryConditions.hpp"

#include "utopia_petsc_DMDA.hpp"
#include "utopia_petsc_DMDA_FunctionSpace.hpp"

#include <cmath>

namespace utopia {

    void petsc_tri(Input &in) {
        std::cout << "petsc_tri" << std::endl;
        static const int Dim = 2;

        using Mesh = utopia::PetscDM<Dim>;
        // using Comm = Mesh::Comm;
        using SizeType = Mesh::SizeType;
        using Scalar = Mesh::Scalar;
        using Point = Mesh::Point;

        using FunctionSpace = utopia::FunctionSpace<Mesh, 1, Tri3<Scalar, 2>>;

        FunctionSpace space;
        space.read(in);

        auto &mesh = space.mesh();

        mesh.describe();

        ArrayView<const SizeType> nodes;

        for (SizeType i = 0; i < mesh.n_elements(); ++i) {
            mesh.nodes(i, nodes);
            disp(nodes);
        }

        PetscVector v;

        space.create_vector(v);

        space.sample(v, UTOPIA_LAMBDA(const Point &p)->Scalar { return p[0] * p[1]; });

        rename("f", v);
        space.write("trifun.vtk", v);
    }

    UTOPIA_REGISTER_APP(petsc_tri);

}  // namespace utopia
