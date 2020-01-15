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



        FEFunction(
            const std::shared_ptr<FunctionSpace> &space,
            const std::shared_ptr<Vector> &data)
        : space_(space), data_(data)
        {}

        FEFunction(
            const std::shared_ptr<FunctionSpace> &space)
        : space_(space), data_(std::make_shared<Vector>())
        {
            space_->create_vector(*data_);
        }

        FEFunction(
            FunctionSpace &space)
        : space_(make_ref(space)), data_(std::make_shared<Vector>())
        {
            space_->create_vector(*data_);
        }

        inline DeviceView<Vector, 1> view_device()
        {
            return space_->assembly_view_device(*data_);
        }

        Coefficient<FunctionSpace> coefficient()
        {
            return Coefficient<FunctionSpace>(*space_);
        }

        template<class Quadrature>
        NodalInterpolate<FunctionSpace, Quadrature> value(const Quadrature &q)
        {
            return NodalInterpolate<FunctionSpace, Quadrature>(*space_, q);
        }

        template<class Quadrature>
        GradInterpolate<FunctionSpace, Quadrature> gradient(const Quadrature &q)
        {
            return GradInterpolate<FunctionSpace, Quadrature>(*space_, q);
        }

        // template<class Quadrature>
        // ShapeFunction<FunctionSpace, Quadrature> shape_function(const Quadrature &q)
        // {
        //     return ShapeFunction<FunctionSpace, Quadrature>(*space_, q);
        // }

    private:
        std::shared_ptr<FunctionSpace> space_;
        std::shared_ptr<Vector> data_;
    };
}

#endif //UTOPIA_FE_FUNCTION_HPP
