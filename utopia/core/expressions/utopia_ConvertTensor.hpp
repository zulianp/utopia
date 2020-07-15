#ifndef UTOPIA_CONVERT_TENSOR_HPP
#define UTOPIA_CONVERT_TENSOR_HPP

#include "utopia_Base.hpp"

#include "utopia_Range.hpp"
#include "utopia_Tensor.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <class From,
              class To,
              int Order = Traits<From>::Order,
              int BackendFrom = Traits<From>::Backend,
              int BackendTo = Traits<To>::Backend>
    class ConvertTensor {};

    template <class T, int Order, int BE>
    class ConvertTensor<T, T, Order, BE, BE> {
    public:
        static void apply(const T &in, T &out) { out.copy(in); }
    };

    /**
     * @ingroup    interoperability
     * @brief      Converts the From tensor to the To tensor type
     *
     * @param      from       Source tensor
     * @param      to         Destination tensor
     */
    template <class From, class To, int Order>
    inline void convert(const Tensor<From, Order> &from, Tensor<To, Order> &to) {
        ConvertTensor<From, To>::apply(from.derived(), to.derived());
    }

    template <class T1, class T2, int Order>
    bool cross_backend_approxeq(const Tensor<T1, Order> &l, const Tensor<T2, Order> &r) {
        T1 r_copy;
        convert(r, r_copy);
        return approxeq(l, r_copy, 1e-10);
    }

}  // namespace utopia

#endif  // UTOPIA_CONVERT_TENSOR_HPP
