#ifndef UTOPIA_FRACTURED_POUROUS_MEDIA_HPP
#define UTOPIA_FRACTURED_POUROUS_MEDIA_HPP

#include "utopia_BackgroundModel.hpp"
#include "utopia_EmbeddedModel.hpp"

#include "utopia_fe_base.hpp"

namespace utopia {

    template<class Matrix, class Vector>
    class LagrangeMultiplierTerms {
    public:
        using FunctionSpace = LibMeshFunctionSpace;

        std::shared_ptr<FunctionSpace> space_;
    };

    template<class Matrix, class Vector>
    class FracturedPourousMedia : public Model<Matrix, Vector> {
    public:
        using BackgroundModel = utopia::BackgroundModel<Matrix, Vector>;
        using EmbeddedModel   = utopia::EmbeddedModel<Matrix, Vector>;
        using LagrangeMultiplierTerms = utopia::LagrangeMultiplierTerms<Matrix, Vector>;


        void read(Input &in) override
        {
            in.get("pourous-matrix",   pourous_matrix_);
            
            in.get("fracture-networks", [this](Input &in) {
                in.get_all([this](Input &in) {
                    auto dfn = std::make_shared<EmbeddedModel>(this->comm_);
                    dfn->read(in);
                    fracture_network_.push_back(dfn);
                });;
            });
        }

        inline bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override
        {
            return false;
        }

        FracturedPourousMedia(libMesh::Parallel::Communicator &comm)
        : comm_(comm), pourous_matrix_(comm)
        {}

    private:
        libMesh::Parallel::Communicator &comm_;
        BackgroundModel pourous_matrix_;
        std::vector<std::shared_ptr<EmbeddedModel>> fracture_network_;
        std::vector<std::shared_ptr<LagrangeMultiplierTerms>> lagrange_multiplier_;
    };
    
}

#endif //UTOPIA_FRACTURED_POUROUS_MEDIA_HPP
