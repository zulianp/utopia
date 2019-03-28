#ifndef UTOPIA_HOMEMADE_FUNCTION_SPACE_HPP
#define UTOPIA_HOMEMADE_FUNCTION_SPACE_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_homemade_FEForwardDeclarations.hpp"

#include "utopia_FunctionalTraits.hpp"
#include "utopia_FunctionSpace.hpp"
#include "utopia_homemade_Mesh.hpp"
#include "utopia_homemade_AssemblyContext.hpp"

namespace utopia {

    class HMFESpace : public FunctionSpace<HMFESpace> {
    public:
        HMFESpace(const std::shared_ptr<Mesh> &mesh, const int order = 1)
        : mesh_(mesh), order_(order)
        {
            make_dof_map();
        }

        HMFESpace()
        : mesh_(std::make_shared<Mesh>())
        {}

        inline Mesh &mesh()
        {
            return *mesh_;
        }

        inline const Mesh &mesh() const
        {
            return *mesh_;
        }

        void dof_indices(const int element_index, std::vector<int> &indices);

        void make_dof_map();

        std::shared_ptr<Mesh> mesh_;
        std::vector<int> dof_ptr;
        std::vector<int> dof_index;
        int order_;
    };

    template<>
    class Traits<HMFESpace> {
    public:
        static const int Backend = HOMEMADE;
        static const int Order = 1;
        static const int FILL_TYPE = FillType::DENSE;

        typedef double Scalar;
        typedef utopia::Vectord Vector;
        typedef utopia::Matrixd Matrix;

        typedef utopia::HMFESpace Implementation;
        typedef utopia::HMDerivative GradientType;
        typedef utopia::HMFun DivergenceType;
        typedef utopia::HMJacobian JacobianType;
        typedef utopia::HMDerivative CurlType;

        typedef utopia::Mesh MeshType;
    };

    template<>
    class FunctionalTraits<HMFESpace, AssemblyContext<HOMEMADE> > {
    public:
        inline static int type(const HMFESpace &space,  const AssemblyContext<HOMEMADE> &ctx)
        {
            return utopia::POLYNOMIAL_FUNCTION;
        }

        inline static int order(const HMFESpace &space, const AssemblyContext<HOMEMADE> &ctx)
        {
            return space.mesh().element_order(ctx.current_element);
        }
    };
}

#endif //UTOPIA_HOMEMADE_FUNCTION_SPACE_HPP
