#ifndef UTOPIA_DEVICE_TENSOR_REDUCE_HPP
#define UTOPIA_DEVICE_TENSOR_REDUCE_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_ViewForwardDeclarations.hpp"
#include "utopia_StoreAs.hpp"
#include <string>

namespace utopia {

    template<class Expr, class Op, int Dim>
    class DeviceTensorReduce : public DeviceExpression<DeviceTensorReduce<Expr, Op, Dim>> {
    public:
        using SizeType = typename Traits<Expr>::SizeType;
        using Scalar   = typename Traits<Expr>::Scalar;

        UTOPIA_INLINE_FUNCTION DeviceTensorReduce(const Expr &expr)
        : expr_(expr)
        {}

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i) const
        {
            if(Dim == 1) {
                const SizeType cols = expr_.cols();
                Scalar ret = 0.0;
                for(SizeType j = 0; j < cols; ++j) {
                    ret += expr_(i, j);
                }

                return ret;

            } else {
                const SizeType rows = expr_.rows();
                Scalar ret = 0.0;
                for(SizeType j = 0; j < rows; ++j) {
                    ret += expr_(j, i);
                }

                return ret;
            }
        }

        inline std::string get_class() const override
        {
            return "DeviceTensorReduce<" + expr_.get_class() + "," + GetClass<Op>() + "," + std::to_string(Dim) + ">";
        }

    private:
        UTOPIA_STORE_CONST(Expr)  expr_;
    };

    template<class Expr, class Op, int Dim>
    class Traits< DeviceTensorReduce<Expr, Op, Dim> > : public Traits<typename Traits<Expr>::Vector> {};

}

#endif //UTOPIA_DEVICE_TENSOR_REDUCE_HPP
