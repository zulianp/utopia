#ifndef UTOPIA_FRACTURE_FLOW_HPP
#define UTOPIA_FRACTURE_FLOW_HPP

#include "utopia_Describable.hpp"
#include "utopia_FluxPostProcessor.hpp"
#include "utopia_Input.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIScalarSampler.hpp"
#include "utopia_WeakDirichletBoundaryConditions.hpp"

#include "libmesh/parallel_mesh.h"

namespace utopia {

    class FractureFlow final : public Configurable, public Describable {
    public:
        typedef utopia::LibMeshFunctionSpace FunctionSpaceT;
        typedef utopia::Traits<FunctionSpaceT> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;

        FractureFlow(libMesh::Parallel::Communicator &comm);
        void read(Input &is) override;
        bool empty() const;
        void describe(std::ostream &os = std::cout) const override;

        UIMesh<libMesh::DistributedMesh> mesh;
        UIFunctionSpace<FunctionSpaceT> space;
        std::shared_ptr<UIForcingFunction<FunctionSpaceT, UVector>> forcing_function;
        std::shared_ptr<UIFunction<double>> sampler;
        ElementMatrix diffusion_tensor;

        std::vector<std::shared_ptr<PostProcessor<FunctionSpaceT, UVector>>> post_processors_;

        void post_process(const UVector &sol);
        void export_post_process();

        void post_process(LibMeshFunctionSpace &space, const UVector &pressure, const UVector &concentration);

        void apply_weak_BC(USparseMatrix &A, UVector &b) const;

    private:
        std::shared_ptr<WeakDirichletBoundaryConditions<FunctionSpaceT, USparseMatrix, UVector>> weak_BC_;
    };

    class FractureFlowAuxSystem {
    public:
        FractureFlowAuxSystem(LibMeshFunctionSpace &V, const std::string &name = "aperture");
        void sample(const std::shared_ptr<UIFunction<double>> &sampler);

        UVector sampled;

    private:
        libMesh::LinearImplicitSystem &aux_;
        std::vector<int> var_nums_;
    };

}  // namespace utopia

#endif  // UTOPIA_FRACTURE_FLOW_HPP
