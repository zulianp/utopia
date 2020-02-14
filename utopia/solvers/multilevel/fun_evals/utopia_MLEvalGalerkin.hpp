#ifndef UTOPIA_ML_EVAL_GALERKIN_HPP
#define UTOPIA_ML_EVAL_GALERKIN_HPP

#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_LevelMemory.hpp"
#include "utopia_MultiLevelEvaluations.hpp"

namespace utopia
{
    template<typename Matrix, typename Vector>
    class MultilevelDerivEval<Matrix, Vector, GALERKIN> final
    {
        typedef UTOPIA_SCALAR(Vector)           Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)        SizeType;

    public:
        MultilevelDerivEval(const SizeType & nl_levels): n_levels_(nl_levels), initialized_(false)
        {

        }

        inline Scalar compute_energy(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & s_global)
        {
            if(level < n_levels_-1){
                help_[level] = H_diff[level] * s_global; 
                return (dot(g_diff[level], s_global) + 0.5 * dot(help_[level], s_global));
            }
            else
            {
                Scalar energy = 0.0; 
                fun.value(x, energy); 
                return energy; 
            }
        }

        // s_global is assummed to be zero 
        inline Scalar compute_energy(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x)
        {
            if(level < n_levels_-1){
                return 0.0;
            }
            else{
                Scalar energy = 0.0; 
                fun.value(x, energy); 
                return energy; 
            }
        }                

        inline bool compute_gradient(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & s_global)
        {
            if(level < n_levels_-1){
                g[level] = g_diff[level] + (H_diff[level] * s_global);
            }
            else
            {
                fun.gradient(x, g[level]); 
            }

            return true;
        }

        // s_global is assummed to be zero 
        inline bool compute_gradient(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x)
        {
            if(level < n_levels_-1){
                g[level] = g_diff[level];
            }
            else{
                fun.gradient(x, g[level]); 
            }

            return true;
        }        

        inline Scalar compute_gradient_energy(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & s_global)
        {
            if(level < n_levels_-1)
            {
                help_[level]    = H_diff[level] * s_global; 
                g[level]        = g_diff[level] + help_[level];
                return (dot(g_diff[level], s_global) + 0.5 * dot(help_[level], s_global));
            }
            else
            {
                Scalar energy = 0.0; 
                fun.value(x, energy); 
                fun.gradient(x, g[level]); 
                return energy; 
            }
        }

        inline bool compute_hessian(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x)
        {
            if(level < n_levels_-1){
                H[level] = H_diff[level];
            }
            else
            {
                fun.hessian(x, H[level]); 
            }
            return true;
        }

        void init_memory(const std::vector<SizeType> & n_dofs_, const std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector> > > & level_functions)
        {
            help_.resize(n_levels_); 
            g_diff.resize(n_levels_);
            g.resize(n_levels_);
            H_diff.resize(n_levels_);
            H.resize(n_levels_);

            for(auto l=0; l < n_levels_; l++){
                help_[l]    = local_zeros(n_dofs_[l]); 
                g_diff[l]   = local_zeros(n_dofs_[l]); 
                g[l]        = local_zeros(n_dofs_[l]); 

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


#endif //UTOPIA_ML_EVAL_GALERKIN_HPP