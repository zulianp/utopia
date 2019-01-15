#ifndef UTOPIA_ML_EVAL_HPP
#define UTOPIA_ML_EVAL_HPP

#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_ExtendedFunction.hpp"


namespace utopia 
{

    enum MultiLevelCoherence{   FIRST_ORDER  = 1, 
                                SECOND_ORDER = 2, 
                                GALERKIN     = 0};

    template <MultiLevelCoherence T>
    struct is_first_order 
    {
      static const bool value = false;
    };

    template <>
    struct is_first_order<FIRST_ORDER> {
      static const bool value = true;
    };                                

//  --------------------------------------------- Hessians ---------------------------------------------------------------------------

    template<typename Matrix, typename Vector, MultiLevelCoherence MC>  
    class MultilevelHessianEval;

    template<typename Matrix, typename Vector>
    class MultilevelHessianEval<Matrix, Vector, GALERKIN>
    {
        public:
           inline static bool compute_hessian(const ExtendedFunction<Matrix, Vector> & /*fun*/, const Vector & /*x*/, Matrix & H, const Matrix & H_diff)
            {
                H = H_diff; 
                return true; 
            }
    }; 


    template<typename Matrix, typename Vector>
    class MultilevelHessianEval<Matrix, Vector, SECOND_ORDER>
    {
        public: 
           inline static bool compute_hessian(const ExtendedFunction<Matrix, Vector> & fun, const Vector & x,  Matrix & H, const Matrix & H_diff)
            {
                fun.hessian(x, H); 
                H += H_diff; 
                return true; 
            }
    }; 

    template<typename Matrix, typename Vector>
    class MultilevelHessianEval<Matrix, Vector, FIRST_ORDER>
    {
        public: 
           inline static bool compute_hessian(const ExtendedFunction<Matrix, Vector> & fun, const Vector & x,  Matrix & H, const Matrix & /*H_diff*/)
            {
                fun.hessian(x, H); 
                return true; 
            }
    }; 


//  --------------------------------------------- Gradients ---------------------------------------------------------------------------

    template<typename Matrix, typename Vector, MultiLevelCoherence MC>  
    class MultilevelGradientEval;

    template<typename Matrix, typename Vector>
    class MultilevelGradientEval<Matrix, Vector, GALERKIN>
    {
        public:
            inline static bool compute_gradient(const ExtendedFunction<Matrix, Vector> & /*fun*/, const Vector & /*x*/,  Vector & g, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global)
            {   
                g = g_diff + H_diff * s_global; 
                return true; 
            }
    }; 


    template<typename Matrix, typename Vector>
    class MultilevelGradientEval<Matrix, Vector, SECOND_ORDER>
    {
        public: 
            inline static bool compute_gradient(const ExtendedFunction<Matrix, Vector> & fun, const Vector & x,  Vector & g, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global)
            {
                fun.gradient(x, g);
                g += g_diff + H_diff * s_global; 
                return true; 
            }
    }; 

    template<typename Matrix, typename Vector>
    class MultilevelGradientEval<Matrix, Vector, FIRST_ORDER>
    {
        public: 
            inline static bool compute_gradient(const ExtendedFunction<Matrix, Vector> & fun, const Vector & x,  Vector & g, const Vector & g_diff, const Matrix & /*H_diff*/, const Vector & /*s_global*/)
            {
                fun.gradient(x, g);
                g += g_diff; 
                return true; 
            }
    }; 



//  --------------------------------------------- Energies ---------------------------------------------------------------------------

    template<typename Matrix, typename Vector, MultiLevelCoherence MC>  
    class MultilevelEnergyEval;

    template<typename Matrix, typename Vector>
    class MultilevelEnergyEval<Matrix, Vector, GALERKIN>
    {
        public:
            inline static typename Traits<Vector>::Scalar compute_energy(const ExtendedFunction<Matrix, Vector> & /*fun*/, const Vector & /*x*/, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global)
            {
                return (dot(g_diff, s_global) + 0.5 * dot(H_diff * s_global, s_global)); 
            }
    }; 


    template<typename Matrix, typename Vector>
    class MultilevelEnergyEval<Matrix, Vector, SECOND_ORDER>
    {
        public: 
            inline static typename Traits<Vector>::Scalar compute_energy(const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global)
            {
                typename Traits<Vector>::Scalar energy = 0.0; 
                fun.value(x, energy); 
                energy += dot(g_diff, s_global) + 0.5 * dot(H_diff * s_global, s_global); 
                return energy; 
            }
    }; 

    template<typename Matrix, typename Vector>
    class MultilevelEnergyEval<Matrix, Vector, FIRST_ORDER>
    {
        public: 
            inline static typename Traits<Vector>::Scalar compute_energy(const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & g_diff, const Matrix & /*H_diff*/, const Vector & s_global)
            {
                typename Traits<Vector>::Scalar energy = 0.0; 
                fun.value(x, energy); 
                energy += dot(g_diff, s_global); 
                return energy; 
            }
    }; 



}




#endif //UTOPIA_RMTR_INFTY_HPP