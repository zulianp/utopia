#ifndef UTOPIA_FRACTURED_POUROUS_MEDIA_HPP
#define UTOPIA_FRACTURED_POUROUS_MEDIA_HPP

#include "utopia_BackgroundModel.hpp"
#include "utopia_EmbeddedModel.hpp"

#include "utopia_fe_base.hpp"

namespace utopia {

    template<class Matrix, class Vector>
    class FracturedPourousMedia : public Model<Matrix, Vector> {
    public:
        using BackgroundModel = utopia::BackgroundModel<Matrix, Vector>;
        using EmbeddedModel   = utopia::EmbeddedModel<Matrix, Vector>;


        void read(Input &in) override
        {
            in.get("pourous-matrix",   pourous_matrix_);
            in.get("fracture-network", fracture_network_);
        }

        inline bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override
        {
            return false;
        }

    private:
        BackgroundModel pourous_matrix_;
        EmbeddedModel   fracture_network_;
    };
    
}

#endif //UTOPIA_FRACTURED_POUROUS_MEDIA_HPP
