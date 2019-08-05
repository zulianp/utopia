#ifndef UTOPIA_FE_KERNEL_HPP
#define UTOPIA_FE_KERNEL_HPP

#include "utopia_FEForwardDeclarations.hpp"
#include <utility>
#include <functional>

namespace utopia {
    template<typename T, int Order_, class Fun>
    class FEKernel : public Expression< FEKernel<T, Order_, Fun> > {
    public:
        static const int Order = Order_;
        typedef T Scalar;


        FEKernel(Fun fun, const unsigned int integration_order)
        : fun_(fun), integration_order_(integration_order)
        {}
        template<class... Args>
        void apply(Args &&...args)
        {
            fun_(std::forward<Args...>(args...));
        }

        std::string getClass() const override
        {
            return "FEKernel";
        }

        unsigned int integration_order() const
        {
            return integration_order_;
        }

    private:
        Fun fun_;
        unsigned int integration_order_;
    };


    template<typename T, class Fun>
    FEKernel<T, 2, Fun> bilinear_kernel(Fun fun, const unsigned int integration_order)
    {
        return FEKernel<T, 2, Fun>(fun, integration_order);
    }

    template<typename T, class Fun>
    FEKernel<T, 1, Fun> linear_kernel(Fun fun, const unsigned int integration_order)
    {
        return FEKernel<T, 1, Fun>(fun, integration_order);
    }
}

#endif //UTOPIA_FE_KERNEL_HPP
