#ifndef UTOPIA_ML_EVAL_HPP
#define UTOPIA_ML_EVAL_HPP

#include "utopia_Core.hpp"
#include "utopia_Function.hpp"


namespace utopia 
{

    enum MultiLevelCoherence{   FIRST_ORDER  = 1, 
                                SECOND_ORDER = 2, 
                                GALERKIN     = 0};

//  --------------------------------------------- Hessians ---------------------------------------------------------------------------

    template<typename Matrix, typename Vector, typename FunctionType, MultiLevelCoherence MC>  
    class MultilevelHessianEval;

    template<typename Matrix, typename Vector, typename FunctionType>
    class MultilevelHessianEval<Matrix, Vector, FunctionType, GALERKIN>
    {
        public:
           inline static bool compute_hessian(const FunctionType & fun, const Vector & x, Matrix & H, const Matrix & H_diff)
            {
                H = H_diff; 
                return true; 
            }
    }; 


    template<typename Matrix, typename Vector, typename FunctionType>
    class MultilevelHessianEval<Matrix, Vector, FunctionType, SECOND_ORDER>
    {
        public: 
           inline static bool compute_hessian(const FunctionType & fun, const Vector & x,  Matrix & H, const Matrix & H_diff)
            {
                fun.hessian(x, H); 
                H += H_diff; 
                return true; 
            }
    }; 

    template<typename Matrix, typename Vector, typename FunctionType>
    class MultilevelHessianEval<Matrix, Vector, FunctionType, FIRST_ORDER>
    {
        public: 
           inline static bool compute_hessian(const FunctionType & fun, const Vector & x,  Matrix & H, const Matrix & H_diff)
            {
                fun.hessian(x, H); 
                return true; 
            }
    }; 


//  --------------------------------------------- Gradients ---------------------------------------------------------------------------

    template<typename Matrix, typename Vector, typename FunctionType, MultiLevelCoherence MC>  
    class MultilevelGradientEval;

    template<typename Matrix, typename Vector, typename FunctionType>
    class MultilevelGradientEval<Matrix, Vector, FunctionType, GALERKIN>
    {
        public:
            inline static bool compute_gradient(const FunctionType & fun, const Vector & x,  Vector & g, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global)
            {
                g = g_diff + H_diff * s_global; 
                return true; 
            }
    }; 


    template<typename Matrix, typename Vector, typename FunctionType>
    class MultilevelGradientEval<Matrix, Vector, FunctionType, SECOND_ORDER>
    {
        public: 
            inline static bool compute_gradient(const FunctionType & fun, const Vector & x,  Vector & g, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global)
            {
                fun.gradient(x, g);
                g += g_diff + H_diff * s_global; 
                return true; 
            }
    }; 

    template<typename Matrix, typename Vector, typename FunctionType>
    class MultilevelGradientEval<Matrix, Vector, FunctionType, FIRST_ORDER>
    {
        public: 
            inline static bool compute_gradient(const FunctionType & fun, const Vector & x,  Vector & g, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global)
            {
                fun.gradient(x, g);
                g += g_diff; 
                return true; 
            }
    }; 



//  --------------------------------------------- Energies ---------------------------------------------------------------------------

    template<typename Matrix, typename Vector, typename FunctionType, MultiLevelCoherence MC>  
    class MultilevelEnergyEval;

    template<typename Matrix, typename Vector, typename FunctionType>
    class MultilevelEnergyEval<Matrix, Vector, FunctionType, GALERKIN>
    {
        public:
            inline static typename Traits<Vector>::Scalar compute_energy(const FunctionType & fun, const Vector & x, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global)
            {
                return (dot(g_diff, s_global) + 0.5 * dot(s_global, H_diff * s_global)); 
            }
    }; 


    template<typename Matrix, typename Vector, typename FunctionType>
    class MultilevelEnergyEval<Matrix, Vector, FunctionType, SECOND_ORDER>
    {
        public: 
            inline static typename Traits<Vector>::Scalar compute_energy(const FunctionType & fun, const Vector & x, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global)
            {
                typename Traits<Vector>::Scalar energy = 0.0; 
                fun.value(x, energy); 
                energy += dot(g_diff, s_global) + 0.5 * dot(s_global, H_diff * s_global); 
                return energy; 
            }
    }; 

    template<typename Matrix, typename Vector, typename FunctionType>
    class MultilevelEnergyEval<Matrix, Vector, FunctionType, FIRST_ORDER>
    {
        public: 
            inline static typename Traits<Vector>::Scalar compute_energy(const FunctionType & fun, const Vector & x, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global)
            {
                typename Traits<Vector>::Scalar energy = 0.0; 
                fun.value(x, energy); 
                energy += dot(g_diff, s_global); 
                return energy; 
            }
    }; 



}




#endif //UTOPIA_RMTR_INFTY_HPP