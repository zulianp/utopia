#ifndef UTOPIA_PARALLEL_EACH_HPP
#define UTOPIA_PARALLEL_EACH_HPP

#include <string>
#include "utopia_Each.hpp"

namespace utopia {
    template <class Tensor, int Order = Tensor::Order, int FILL_TYPE = Traits<Tensor>::FILL_TYPE>
    class ParallelEach : public Each<Tensor, Order, FILL_TYPE> {
    public:
        using Super = utopia::Each<Tensor, Order, FILL_TYPE>;

        using Super::apply_read;
        using Super::apply_transform;
        using Super::apply_write;

        template <class Fun>
        inline static void apply_read(const Tensor &v, Fun fun, const std::string &) {
            Super::apply_read(v, fun);
        }

        template <class Fun>
        inline static void apply_write(Tensor &v, Fun fun, const std::string &) {
            Super::apply_write(v, fun);
        }

        template <class Fun>
        inline static void apply_transform(Tensor &v, Fun fun, const std::string &) {
            Super::apply_transform(v, fun);
        }
    };

    template <class Tensor, class Fun>
    inline void parallel_each_read(const Tensor &v, Fun fun, const std::string &name = "utopia::parallel_each_read") {
        ParallelEach<Tensor>::apply_read(v, fun, name);
    }

    template <class Tensor, class Fun>
    inline void parallel_each_write(Tensor &v, Fun fun, const std::string &name = "utopia::parallel_each_write") {
        ParallelEach<Tensor>::apply_write(v, fun, name);
    }

    template <class Tensor, class Fun>
    inline void parallel_transform(Tensor &v, Fun fun, const std::string &name = "utopia::parallel_transform") {
        ParallelEach<Tensor>::apply_transform(v, fun, name);
    }
}  // namespace utopia

#endif  // UTOPIA_PARALLEL_EACH_HPP
