#ifndef UTOPIA_FE_FUNCTION_HPP
#define UTOPIA_FE_FUNCTION_HPP

#include <memory>
#include "utopia_DeviceView.hpp"
#include "utopia_Coefficient.hpp"
#include "utopia_AssemblyView.hpp"

namespace utopia {

    template<class FunctionSpaceView, class DataView, int Order>
    class FEFunctionView {};

    template<class FunctionSpaceView, class DataView>
    class FEFunctionView<FunctionSpaceView, DataView, 1> {
    public:
        using Elem = typename FunctionSpaceView::Elem;

        FEFunctionView(const FunctionSpaceView &space, const DataView &data)
        : space_(space), data_(data)
        {}

    private:
        FunctionSpaceView space_;
        DataView data_;
    };

    template<class FunctionSpace>
    class FEFunction {
    public:
        using Vector                  = typename FunctionSpace::Vector;
        using FunctionSpaceViewDevice = typename FunctionSpace::ViewDevice;

        using VectorViewDevice        = utopia::DeviceView<Vector, 1>;
        using ViewDevice              = utopia::FEFunctionView<FunctionSpaceViewDevice, VectorViewDevice, 1>;
        using Coefficient             = utopia::Coefficient<FunctionSpace>;

        FEFunction(
            const std::shared_ptr<Coefficient> &coeff)
        : space_(make_ref(coeff->space())), coeff_(coeff)
        {}

        FEFunction(
            const std::shared_ptr<FunctionSpace> &space,
            Vector &global_vector)
        : space_(space), coeff_(std::make_shared<Coefficient>(*space_, global_vector))
        {}

        FEFunction(
            const std::shared_ptr<FunctionSpace> &space
        )
        : space_(space), coeff_(std::make_shared<Coefficient>(*space_))
        {}

        FEFunction(
            const FunctionSpace &space,
            Vector &global_vector)
        : space_(make_ref(space)), coeff_(std::make_shared<Coefficient>(*space_, global_vector))
        {}

        std::shared_ptr<Coefficient> coefficient()
        {
            return coeff_;
        }

        void update(const Vector &x)
        {
            coeff_->update(x);
        }

        template<class Quadrature>
        NodalInterpolate<FunctionSpace, Quadrature> value(const Quadrature &q)
        {
            return NodalInterpolate<FunctionSpace, Quadrature>(coeff_, q);
        }

        template<class Quadrature>
        GradInterpolate<FunctionSpace, Quadrature> gradient(const Quadrature &q)
        {
            return GradInterpolate<FunctionSpace, Quadrature>(coeff_, q);
        }

    private:
        std::shared_ptr<const FunctionSpace> space_;
        std::shared_ptr<Coefficient> coeff_;
    };
}

#endif //UTOPIA_FE_FUNCTION_HPP
