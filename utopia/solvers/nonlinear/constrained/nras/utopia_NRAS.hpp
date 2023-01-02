#ifndef UTOPIA_NRAS_HPP
#define UTOPIA_NRAS_HPP

#include "utopia_Traits.hpp"

#include <memory>

namespace utopia {
    template <class Matrix, class Vector = typename Traits<Matrix>::Vector>
    class NRAS {
    public:
        using Function = utopia::Function<Matrix, Vector>;

        NRAS();
        ~NRAS();

        void set_local_function(const std::shared_ptr<Function> &lfun);
        void set_global_function(const std::shared_ptr<Function> &gfun);
        void set_cutoff_operator(const std::shared_ptr<Operator<Vector>> &op);
        void set_projection_operator(const std::shared_ptr<Operator<Vector>> &op);
        bool solve(Vector &x);

        class Impl;

    private:
        std::unique_ptr<Impl> impl_;
    };
}  // namespace utopia

#endif  // UTOPIA_NRAS_HPP
