#ifndef UTOPIA_NEW_FLOW_HPP
#define UTOPIA_NEW_FLOW_HPP

#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_Flow.hpp"
#include "utopia_Integral.hpp"
#include "utopia_Model.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIScalarSampler.hpp"

#include <iostream>

namespace utopia {

    template <class FunctionSpace, class Matrix, class Vector>
    class NewFlow final : public Model<Matrix, Vector> {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);
        typedef utopia::Traits<FunctionSpace> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;

        NewFlow(FunctionSpace &space);

        bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override;
        bool assemble_lower_dimensional_features(const Vector &x, Matrix &hessian, Vector &gradient);
        void read(Input &in) override;

        inline bool is_linear() const override { return true; }
        inline void clear() override {}
        inline void rescale(const Scalar rescale) { rescale_ = rescale; }

    private:
        FunctionSpace &space_;
        UIForcingFunction<FunctionSpace, Vector> forcing_function_;

        UIScalarFunction<Scalar> permeability_;
        ElementMatrix diffusion_tensor_;

        std::vector<std::shared_ptr<UIFunction<Scalar>>> lower_dimensional_permeability_;
        std::vector<int> lower_dimensional_tags_;
        Scalar rescale_;

        void read_permeability_tensor(Input &in);
    };

}  // namespace utopia

#endif  // UTOPIA_NEW_FLOW_HPP
