#ifndef UTOPIA_FLOW_WITH_FRACTURES_HPP
#define UTOPIA_FLOW_WITH_FRACTURES_HPP

#include "utopia_Model.hpp"
#include "utopia_UIScalarSampler.hpp"
#include "utopia_FEForwardDeclarations.hpp"

#include <iostream>

namespace utopia {

    template<class FunctionSpace, class Matrix, class Vector>
    class FlowWithFractures final : public Model<Matrix, Vector> {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);
        typedef utopia::Traits<FunctionSpace> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;
        typedef typename TraitsT::Vector ElementVector;

        FlowWithFractures(FunctionSpace &space);
        ~FlowWithFractures();

        class Impl;

        bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override;

        bool is_linear() const override;

        void clear() override;

        void read(Input &in) override;

        void rescale(const Scalar rescale);

        std::shared_ptr<UIFunction<Scalar>> sampler();
        std::shared_ptr<UIFunction<USerialMatrix>> tangent_sampler();

    private:
        std::unique_ptr<Impl> impl_;
    };
}

#endif //UTOPIA_FLOW_WITH_FRACTURES_HPP
