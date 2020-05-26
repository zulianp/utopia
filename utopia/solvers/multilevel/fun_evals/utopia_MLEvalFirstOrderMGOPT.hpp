#ifndef UTOPIA_ML_EVAL_FIRST_ORDER_MGOPT_HPP
#define UTOPIA_ML_EVAL_FIRST_ORDER_MGOPT_HPP

#include "utopia_Core.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_Function.hpp"
#include "utopia_LevelMemory.hpp"
#include "utopia_MultiLevelEvaluations.hpp"

namespace utopia {
    template <typename Matrix, typename Vector>
    class MultilevelDerivEval<Matrix, Vector, FIRST_ORDER_MGOPT> final {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

    public:
        MultilevelDerivEval(const SizeType &nl_levels) : n_levels_(nl_levels), initialized_(false) {}

        inline Scalar compute_energy(const SizeType &level,
                                     const ExtendedFunction<Matrix, Vector> &fun,
                                     const Vector &x,
                                     const Vector & /*s_global*/) {
            Scalar energy = 0.0;
            fun.value(x, energy);

            if (level < n_levels_ - 1) {
                energy += dot(g_diff[level], x);
            }

            return energy;
        }

        // s_global is assummed to be zero
        inline Scalar compute_energy(const SizeType &level,
                                     const ExtendedFunction<Matrix, Vector> &fun,
                                     const Vector &x) {
            Scalar energy = 0.0;
            fun.value(x, energy);

            if (level < n_levels_ - 1) {
                energy += dot(g_diff[level], x);
            }

            return energy;
        }

        inline bool compute_gradient(const SizeType &level,
                                     const ExtendedFunction<Matrix, Vector> &fun,
                                     const Vector &x,
                                     const Vector & /* s_global*/) {
            fun.gradient(x, g[level]);

            if (level < n_levels_ - 1) {
                g[level] += g_diff[level];
            }

            return true;
        }

        // s_global is assummed to be zero
        inline bool compute_gradient(const SizeType &level,
                                     const ExtendedFunction<Matrix, Vector> &fun,
                                     const Vector &x) {
            fun.gradient(x, g[level]);

            if (level < n_levels_ - 1) {
                g[level] += g_diff[level];
            }

            return true;
        }

        inline Scalar compute_gradient_energy(const SizeType &level,
                                              const ExtendedFunction<Matrix, Vector> &fun,
                                              const Vector &x,
                                              const Vector & /*s_global*/) {
            Scalar energy = 0.0;
            fun.value(x, energy);
            fun.gradient(x, g[level]);

            if (level < n_levels_ - 1) {
                energy += dot(g_diff[level], x);
                g[level] += g_diff[level];
            }

            return energy;
        }

        inline bool compute_hessian(const SizeType &level,
                                    const ExtendedFunction<Matrix, Vector> &fun,
                                    const Vector &x) {
            return fun.hessian(x, H[level]);
        }

        void init_memory(const std::vector<Layout> &layouts,
                         const std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector> > > &level_functions) {
            g_diff.resize(n_levels_);
            g.resize(n_levels_);
            H.resize(n_levels_);

            for (auto l = 0; l < n_levels_; l++) {
                g_diff[l].zeros(layouts[l]);
                g[l].zeros(layouts[l]);
                level_functions[l]->initialize_hessian(H[l], H[l]);
            }

            initialized_ = true;
        }

        bool initialized() const { return initialized_; }

    private:
        SizeType n_levels_;
        bool initialized_;

    public:
        std::vector<Vector> g, g_diff;
        std::vector<Matrix> H;
    };

}  // namespace utopia

#endif  // UTOPIA_ML_EVAL_FIRST_ORDER_MGOPT_HPP