#ifndef UTOPIA_ML_EVAL_SECOND_ORDER_HPP
#define UTOPIA_ML_EVAL_SECOND_ORDER_HPP

#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_LevelMemory.hpp"
#include "utopia_MultiLevelEvaluations.hpp"

namespace utopia
{
    template<typename Matrix, typename Vector>
    class MultilevelDerivEval<Matrix, Vector, SECOND_ORDER> final
    {
        using Scalar   = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout   = typename Traits<Vector>::Layout;

    public:

        MultilevelDerivEval(const SizeType & nl_levels): n_levels_(nl_levels), initialized_(false)
        {

        }

        inline Scalar compute_energy(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & s_global)
        {
            Scalar energy = 0.0;
            fun.value(x, energy);

            if(level < n_levels_-1){
                help_[level] = H_diff[level] * s_global;
                energy += (0.5 * dot(help_[level], s_global)) + dot(g_diff[level], s_global);
            }
            return energy;
        }

        // s_global is assummed to be zero
        inline Scalar compute_energy(const SizeType & /*level*/,
                                     const ExtendedFunction<Matrix, Vector> &fun,
                                     const Vector &x) {
            Scalar energy = 0.0;
            fun.value(x, energy);
            return energy;
        }

        inline bool compute_gradient(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & s_global)
        {
            UTOPIA_NO_ALLOC_BEGIN("RMTR::compute_gradient1");
            fun.gradient(x, g[level]);
            UTOPIA_NO_ALLOC_END();

            if(level < n_levels_-1){
                UTOPIA_NO_ALLOC_BEGIN("RMTR::compute_gradient2");
                help_[level] = H_diff[level] * s_global;
                UTOPIA_NO_ALLOC_END();

                UTOPIA_NO_ALLOC_BEGIN("RMTR::compute_gradient3");
                // g[level] += g_diff[level] + help_[level];
                g[level] = g[level] + g_diff[level] + help_[level];
                UTOPIA_NO_ALLOC_END();
            }
            return true;
        }

        // s_global is assummed to be zero
        inline bool compute_gradient(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x)
        {
            fun.gradient(x, g[level]);

            if(level < n_levels_-1){
                g[level] += g_diff[level];
            }
            return true;
        }

        inline Scalar compute_gradient_energy(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & s_global)
        {
            Scalar energy = 0.0;
            fun.value(x, energy);
            fun.gradient(x, g[level]);

            if(level < n_levels_-1){
                help_[level] = H_diff[level] * s_global;

                energy += (0.5 * dot(help_[level], s_global)) + dot(g_diff[level], s_global);
                g[level] += g_diff[level] + help_[level];
            }

            return energy;
        }

       inline bool compute_hessian(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x)
        {
            fun.hessian(x, H[level]);

            if(level < n_levels_-1){
                H[level] += H_diff[level];
            }
            return true;
        }

        void init_memory(const std::vector<Layout> &layouts, const std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector> > > & level_functions)
        {
            help_.resize(n_levels_);
            g_diff.resize(n_levels_);
            g.resize(n_levels_);
            H_diff.resize(n_levels_);
            H.resize(n_levels_);

            for(auto l=0; l < n_levels_; l++){
                help_[l].zeros(layouts[l]);
                g_diff[l].zeros(layouts[l]);
                g[l].zeros(layouts[l]);

                level_functions[l]->initialize_hessian(H[l], H[l]);
                H_diff[l] = H[l];
            }

            initialized_ = true;
        }

        bool initialized() const
        {
            return initialized_;
        }

        private:
            SizeType n_levels_;
            bool initialized_;
        public:
            std::vector<Vector> g, g_diff,  help_;
            std::vector<Matrix> H, H_diff;

    };

}

#endif //UTOPIA_ML_EVAL_SECOND_ORDER_HPP