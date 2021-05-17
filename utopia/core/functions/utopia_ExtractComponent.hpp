#ifndef UTOPIA_EXTRACT_COMPONENT_HPP
#define UTOPIA_EXTRACT_COMPONENT_HPP

#include "utopia_Traits.hpp"
#include "utopia_Layout.hpp"
#include "utopia_Tensor.hpp"

namespace utopia {
    template <class Vector>
    class VectorComponent {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        static void apply_extract(const Vector &from, const int block_size, const int c, Vector &to)  {
            auto from_lo = layout(from);

            auto to_size = from_lo.size() / block_size;
            auto to_local_size = from_lo.local_size() / block_size;

            if(to.empty()) {
                to.zeros(layout(from_lo.comm(), to_local_size, to_size));
            } else {
                assert(to_local_size == to.local_size());
                assert(to_size == to.size());
            }

            auto from_view = local_view_device(from);
            auto to_view = local_view_device(to);

            parallel_for(
                local_range_device(to), UTOPIA_LAMBDA(const SizeType i) {
                    to_view.set(i, from_view.get(i * block_size + c));
                });

        }


        static void apply_set(const Vector &from, const int block_size, const int c, Vector &to)  {
            auto from_lo = layout(from);

            auto to_size = from_lo.size() * block_size;
            auto to_local_size = from_lo.local_size() * block_size;

            if(to.empty()) {
                to.zeros(layout(from_lo.comm(), to_local_size, to_size));
            } else {
                assert(to_local_size == to.local_size());
                assert(to_size == to.size());
            }

            auto from_view = local_view_device(from);
            auto to_view = local_view_device(to);

            parallel_for(
                local_range_device(from), UTOPIA_LAMBDA(const SizeType i) {
                    to_view.set(i * block_size + c, from_view.get(i));
                });

        }
    };


    template<typename T>
    void extract_component(const Tensor<T, 1> &from, const int block_size, const int c, Tensor<T, 1> &to)
    {
        VectorComponent<T>::apply_extract(from.derived(), block_size, c, to.derived());
    }


    template<typename T>
    void set_component(const Tensor<T, 1> &from, const int block_size, const int c, Tensor<T, 1> &to)
    {
        VectorComponent<T>::apply_set(from.derived(), block_size, c, to.derived());
    }
}

#endif //UTOPIA_EXTRACT_COMPONENT_HPP
