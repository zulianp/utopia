#ifndef UTOPIA_DEVICE_TRANSPOSE_HPP
#define UTOPIA_DEVICE_TRANSPOSE_HPP

namespace utopia {

    template<class Expr>
    class DeviceTranspose :  public DeviceExpression<DeviceTranspose<Expr>> {
    public:

        using Scalar   = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        UTOPIA_FUNCTION DeviceTranspose(const Expr &expr)
        : expr_(expr) 
        {}

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i, const SizeType &j) const
        {
            return DeviceOp<Scalar, Conjugate>::apply(expr_(j, i));
        }


    private:
        UTOPIA_STORE_CONST(Expr) expr_;

    };
}

#endif //UTOPIA_DEVICE_TRANSPOSE_HPP
