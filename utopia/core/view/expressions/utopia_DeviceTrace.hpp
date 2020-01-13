#ifndef UTOPIA_DEVICE_TRACE_HPP
#define UTOPIA_DEVICE_TRACE_HPP

#include "utopia_Base.hpp"
#include "utopia_Algorithms.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template<class Expr>
    class DeviceTrace {
    public:
        using Scalar   = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        UTOPIA_INLINE_FUNCTION static Scalar apply(const Expr &expr)
        {
            Scalar ret = 0.0;

            const SizeType n = device::min(rows(expr), cols(expr));
            for(SizeType i = 0; i < n; ++i) {
                ret += expr(i, i);
            }

            return ret;
        }
    };

    template<class Expr>
    class Traits< DeviceTrace<Expr> > : public Traits<Expr> {
    public:
        static const int Order = 0;
    };



}

#endif //UTOPIA_DEVICE_TRACE_HPP
