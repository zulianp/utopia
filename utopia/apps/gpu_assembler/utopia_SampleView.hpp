#ifndef UTOPIA_SAMPLE_VIEW_HPP
#define UTOPIA_SAMPLE_VIEW_HPP

#include <utility>

#include "utopia_AssemblyView.hpp"
#include "utopia_FunctionSpace.hpp"
#include "utopia_LaplacianView.hpp"

namespace utopia {

    template <class Elem, class Function, class MemType = typename Elem::MemType, typename...>
    class Sampler {};

    template <class Mesh, class Function, typename... Args>
    class Sampler<FunctionSpace<Mesh, 1, Args...>, Function> {
    public:
        static const int NComponents = 1;
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Vector = typename FunctionSpace::Vector;
        using Scalar = typename FunctionSpace::Scalar;
        using Point = typename FunctionSpace::Point;
        using SizeType = typename FunctionSpace::SizeType;
        using Elem = typename FunctionSpace::ViewDevice::Elem;

        class ViewDevice {
        public:
            template <class Elem, class Accumulator>
            UTOPIA_INLINE_FUNCTION void assemble(const Elem &e, Accumulator &acc) const {
                Point p;
                const auto n = e.n_nodes();
                // only for scalars
                for (SizeType j = 0; j < n; ++j) {
                    e.node(j, p);
                    acc(j) = fun_(p);
                }
            }

            ViewDevice(Function fun) : fun_(std::move(fun)) {}

            Function fun_;
        };

        Sampler(const FunctionSpace &, Function fun) : fun_(std::move(fun)) {}

        ViewDevice view_device() const { return ViewDevice(fun_); }

    private:
        Function fun_;
    };

    template <class FunctionSpace, class Function>
    Sampler<FunctionSpace, Function> sampler(const FunctionSpace &space, Function fun) {
        return Sampler<FunctionSpace, Function>(space, fun);
    }
}  // namespace utopia

#endif  // UTOPIA_SAMPLE_VIEW_HPP
