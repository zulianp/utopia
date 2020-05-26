#ifndef UTOPIA_FE_FUNCTION_HPP
#define UTOPIA_FE_FUNCTION_HPP

#include <functional>
#include "utopia_FEForwardDeclarations.hpp"

namespace utopia {

    template <class Out, class Fun>
    class ContextFunction : public Expression<ContextFunction<Out, Fun> > {
    public:
        virtual ~ContextFunction() {}
        // FIXME
        static const int Order = 0;
        typedef double Scalar;

        template <int Backend>
        auto eval(const AssemblyContext<Backend> &ctx) const -> Out {
            return fun(ctx);
        }

        ContextFunction(Fun fun) : fun(fun) {}

        std::string get_class() const { return "ContextFunction"; }

    private:
        Fun fun;
    };

    // FIXME
    template <class Out, class Fun>
    class Traits<ContextFunction<Out, Fun> > : public Traits<double> {};

    template <class Out, class Fun>
    inline ContextFunction<Out, Fun> ctx_fun(Fun f) {
        return ContextFunction<Out, Fun>(f);
    }
}  // namespace utopia

#endif  // UTOPIA_FE_FUNCTION_HPP
