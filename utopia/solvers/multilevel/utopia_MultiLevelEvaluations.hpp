#ifndef UTOPIA_ML_EVAL_HPP
#define UTOPIA_ML_EVAL_HPP

#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_LevelMemory.hpp"


namespace utopia
{
    enum LocalSolveType{    PRE_SMOOTHING  = 1,
                            POST_SMOOTHING = 2,
                            COARSE_SOLVE   = 0};


    template<class Matrix, class Vector, MultiLevelCoherence CONSISTENCY_TYPE>
    class MultilevelDerivEval  { }; 

    template<typename Matrix, typename Vector>
    class MultilevelDerivEval<Matrix, Vector, FIRST_ORDER> final
    {

        typedef UTOPIA_SCALAR(Vector)           Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)        SizeType;

    public:

        MultilevelDerivEval(const SizeType & nl_levels): n_levels_(nl_levels), initialized_(false)
        {

        }

        inline Scalar compute_energy(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & s_global)
        {
            Scalar energy = 0.0;
            fun.value(x, energy);

            if(level < n_levels_-1)
            {
                energy += dot(g_diff[level], s_global);
            }

            return energy;
        }

        inline bool compute_gradient(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & /* s_global*/)
        {
            fun.gradient(x, g[level]);

            if(level < n_levels_-1)
            {
                g[level] += g_diff[level];
            }

            return true;
        }

        inline Scalar compute_gradient_energy(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & s_global)
        {
            Scalar energy = 0.0;
            fun.value(x, energy);
            fun.gradient(x, g[level]);

            if(level < n_levels_-1)
            {
                energy += dot(g_diff[level], s_global);
                g[level] += g_diff[level];
            }            

            return energy;
        }


        inline bool compute_hessian(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x)
        {
            fun.hessian(x, H[level]);
            return true;
        }


        inline bool init_deriv_loc_solve(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & s_global, const LocalSolveType & solve_type)
        {
            bool make_hess_updates = true;

            if(solve_type==PRE_SMOOTHING && level==this->n_levels()-1)
            {
                return make_hess_updates; 
            }
            else
            {
                if( (solve_type==PRE_SMOOTHING && level < this->n_levels()-1) || (solve_type == COARSE_SOLVE))
                {
                    this->ml_derivs_.g[level] += this->ml_derivs_.g_diff[level]; 
                }
                else
                {
                    this->compute_gradient(level, fun, x, s_global);
                }

                // energy computations ... 
                if(solve_type != POST_SMOOTHING)
                {
                    this->memory_.energy[level] = compute_energy(level, fun,x, s_global);
                }
            }

            return make_hess_updates;             
        }


        void init_memory(const std::vector<SizeType> & n_dofs_)
        {
            g_diff.resize(n_levels_);
            g.resize(n_levels_);

            // H_diff should not be necessary 
            H_diff.resize(n_levels_);
            H.resize(n_levels_);

            for(auto l=0; l < n_levels_; l++){
                g_diff[l]   = local_zeros(n_dofs_[l]); 
                g[l]        = local_zeros(n_dofs_[l]); 
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
            std::vector<Vector> g, g_diff; 
            std::vector<Matrix> H, H_diff;
    }; 


    template<typename Matrix, typename Vector>
    class MultilevelDerivEval<Matrix, Vector, SECOND_ORDER> final
    {

        typedef UTOPIA_SCALAR(Vector)           Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)        SizeType;

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

        inline bool compute_gradient(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & s_global)
        {
            fun.gradient(x, g[level]);

            if(level < n_levels_-1){
                g[level] += g_diff[level] + (H_diff[level] * s_global); 
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

        void init_memory(const std::vector<SizeType> & n_dofs_)
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

        void init_memory(const std::vector<SizeType> & n_dofs_)
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



    // // Matrix free first order 
    // template<typename Matrix, typename Vector>
    // class MultilevelDerivEval<Matrix, Vector, FIRST_ORDER_DF> final
    // {

    //     typedef UTOPIA_SCALAR(Vector)           Scalar;
    //     typedef UTOPIA_SIZE_TYPE(Vector)        SizeType;

    // public:

    //     MultilevelDerivEval(const SizeType & nl_levels): n_levels_(nl_levels), initialized_(false)
    //     {

    //     }

    //     inline Scalar compute_energy(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & s_global)
    //     {
    //         Scalar energy = 0.0;
    //         fun.value(x, energy);

    //         if(level < n_levels_-1)
    //         {
    //             energy += dot(g_diff[level], s_global);
    //         }

    //         return energy;
    //     }

    //     inline bool compute_gradient(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & /* s_global*/)
    //     {
    //         fun.gradient(x, g[level]);

    //         if(level < n_levels_-1)
    //         {
    //             g[level] += g_diff[level];
    //         }

    //         return true;
    //     }

    //    inline bool compute_hessian(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x)
    //     {
    //         // fun.hessian(x, H[level]);
    //         return true;
    //     }

    //     void init_memory(const std::vector<SizeType> & n_dofs_)
    //     {
    //         g_diff.resize(n_levels_);
    //         g.resize(n_levels_);

    //         // H_diff should not be necessary 
    //         // H_diff.resize(n_levels_);
    //         // H.resize(n_levels_);

    //         for(auto l=0; l < n_levels_; l++){
    //             g_diff[l]   = local_zeros(n_dofs_[l]); 
    //             g[l]        = local_zeros(n_dofs_[l]); 
    //         }

    //         initialized_ = true; 
    //     }

    //     bool initialized() const 
    //     {
    //         return initialized_; 
    //     }

    //     private:
    //         SizeType n_levels_; 
    //         bool initialized_; 

    //     // public:            
    //     //     std::vector<Vector> g, g_diff; 
    //     //     std::vector<Matrix> H, H_diff;
    // }; 


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////                                

    template <MultiLevelCoherence T>
    struct is_first_order
    {
      static const bool value = false;
    };

    template <>
    struct is_first_order<FIRST_ORDER> {
      static const bool value = true;
    };

}




#endif //UTOPIA_RMTR_INFTY_HPP