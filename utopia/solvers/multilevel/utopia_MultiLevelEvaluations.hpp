#ifndef UTOPIA_ML_EVAL_HPP
#define UTOPIA_ML_EVAL_HPP

#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_LevelMemory.hpp"


namespace utopia
{

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

        inline Scalar compute_energy(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & g_diff, const Matrix & /*H_diff*/, const Vector & s_global)
        {
            Scalar energy = 0.0;
            fun.value(x, energy);

            if(level < n_levels_-1)
            {
                energy += dot(g_diff, s_global);
            }

            return energy;
        }

        inline bool compute_gradient(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x,  Vector & g, const Vector & g_diff, const Matrix & /*H_diff*/, const Vector & /* s_global*/)
        {
            fun.gradient(x, g);

            if(level < n_levels_-1)
            {
                g += g_diff;
            }

            return true;
        }

       inline bool compute_hessian(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x,  Matrix & H, const Matrix & /*H_diff*/)
        {
            fun.hessian(x, H);
            return true;
        }

        void init_memory(const std::vector<SizeType> & /*n_dofs_*/)
        {
            initialized_ = true; 
        }

        bool initialized() const 
        {
            return initialized_; 
        }

        private:
            SizeType n_levels_; 
            bool initialized_; 

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

        inline Scalar compute_energy(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global)
        {
            Scalar energy = 0.0;
            fun.value(x, energy);

            if(level < n_levels_-1){
                help_[level] = H_diff * s_global; 
                energy += (0.5 * dot(help_[level], s_global)) + dot(g_diff, s_global); 
            }
            return energy;
        }

        inline bool compute_gradient(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x,  Vector & g, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global)
        {
            fun.gradient(x, g);
            if(level < n_levels_-1){
                g += g_diff + (H_diff * s_global); 
            }
            return true;
        }

       inline bool compute_hessian(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x,  Matrix & H, const Matrix & H_diff)
        {
            fun.hessian(x, H);
            if(level < n_levels_-1){
                H += H_diff;
            }
            return true;
        }

        void init_memory(const std::vector<SizeType> & n_dofs_)
        {   
            help_.resize(n_levels_); 

            for(auto l=0; l < n_levels_; l++){
                help_[l] = local_zeros(n_dofs_[l]); 
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
            std::vector<Vector> help_; 

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

        inline Scalar compute_energy(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global)
        {
            if(level < n_levels_-1){
                help_[level] = H_diff * s_global; 
                return (dot(g_diff, s_global) + 0.5 * dot(help_[level], s_global));
            }
            else
            {
                Scalar energy = 0.0; 
                fun.value(x, energy); 
                return energy; 
            }
        }

        inline bool compute_gradient(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x,  Vector & g, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global)
        {
            if(level < n_levels_-1){
                g = g_diff + (H_diff * s_global);
            }
            else
            {
                fun.gradient(x, g); 
            }

            return true;
        }

        inline bool compute_hessian(const SizeType & level, const ExtendedFunction<Matrix, Vector> & fun, const Vector & x,  Matrix & H, const Matrix & H_diff)
        {
            if(level < n_levels_-1){
                H = H_diff;
            }
            else
            {
                fun.hessian(x, H); 
            }
            return true;
        }

        void init_memory(const std::vector<SizeType> & n_dofs_)
        {
            help_.resize(n_levels_); 

            for(auto l=0; l < n_levels_; l++){
                help_[l] = local_zeros(n_dofs_[l]); 
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
            std::vector<Vector> help_; 

    }; 



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


    // template<typename Matrix, typename Vector, MultiLevelCoherence MC>
    // class MultilevelHessianEval;

    // template<typename Matrix, typename Vector, MultiLevelCoherence MC>
    // class MultilevelGradientEval;


    // template<typename Matrix, typename Vector, MultiLevelCoherence MC>
    // class MultilevelEnergyEval;


    // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ////////////////////////////////////////////////// First order///////////////////////////////////////////////////
    // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // template<typename Matrix, typename Vector>
    // class MultilevelEnergyEval<Matrix, Vector, FIRST_ORDER>
    // {
    //     public:
    //         inline static typename Traits<Vector>::Scalar compute_energy(const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & g_diff, const Matrix & /*H_diff*/, const Vector & s_global)
    //         {
    //             typename Traits<Vector>::Scalar energy = 0.0;
    //             fun.value(x, energy);
    //             energy += dot(g_diff, s_global);
    //             return energy;
    //         }
    // };


    // template<typename Matrix, typename Vector>
    // class MultilevelGradientEval<Matrix, Vector, FIRST_ORDER>
    // {
    //     public:
    //         inline static bool compute_gradient(const ExtendedFunction<Matrix, Vector> & fun, const Vector & x,  Vector & g, const Vector & g_diff, const Matrix & /*H_diff*/, const Vector & /* s_global*/)
    //         {
    //             fun.gradient(x, g);
    //             g += g_diff;
    //             return true;
    //         }
    // };


    // template<typename Matrix, typename Vector>
    // class MultilevelHessianEval<Matrix, Vector, FIRST_ORDER>
    // {
    //     public:
    //        inline static bool compute_hessian(const ExtendedFunction<Matrix, Vector> & fun, const Vector & x,  Matrix & H, const Matrix & /*H_diff*/)
    //         {
    //             fun.hessian(x, H);
    //             return true;
    //         }
    // };




    // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ////////////////////////////////////////////////// Second order///////////////////////////////////////////////////
    // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////    

    // template<typename Matrix, typename Vector>
    // class MultilevelEnergyEval<Matrix, Vector, SECOND_ORDER>
    // {
    //     public:
    //         inline static typename Traits<Vector>::Scalar compute_energy(const ExtendedFunction<Matrix, Vector> & fun, const Vector & x, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global)
    //         {
    //             typename Traits<Vector>::Scalar energy = 0.0;
    //             fun.value(x, energy);
    //             energy += (0.5 * dot(H_diff * s_global, s_global)) + dot(g_diff, s_global); 

    //             return energy;
    //         }
    // };

    // template<typename Matrix, typename Vector>
    // class MultilevelGradientEval<Matrix, Vector, SECOND_ORDER>
    // {
    //     public:
    //         inline static bool compute_gradient(const ExtendedFunction<Matrix, Vector> & fun, const Vector & x,  Vector & g, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global)
    //         {
    //             fun.gradient(x, g);
    //             g += g_diff + (H_diff * s_global); 
                
    //             return true;
    //         }
    // };


    // template<typename Matrix, typename Vector>
    // class MultilevelHessianEval<Matrix, Vector, SECOND_ORDER>
    // {
    //     public:
    //        inline static bool compute_hessian(const ExtendedFunction<Matrix, Vector> & fun, const Vector & x,  Matrix & H, const Matrix & H_diff)
    //         {
    //             fun.hessian(x, H);
    //             H += H_diff;
                                
    //             return true;
    //         }
    // };



    // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ////////////////////////////////////////////////// Galerkin   ///////////////////////////////////////////////////
    // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////    

    // template<typename Matrix, typename Vector>
    // class MultilevelEnergyEval<Matrix, Vector, GALERKIN>
    // {
    //     public:
    //         inline static typename Traits<Vector>::Scalar compute_energy(const ExtendedFunction<Matrix, Vector> & /*fun*/, const Vector & /*x*/, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global)
    //         {
    //             return (dot(g_diff, s_global) + 0.5 * dot(H_diff * s_global, s_global));
    //         }
    // };


    // template<typename Matrix, typename Vector>
    // class MultilevelGradientEval<Matrix, Vector, GALERKIN>
    // {
    //     public:
    //         inline static bool compute_gradient(const ExtendedFunction<Matrix, Vector> & /*fun*/, const Vector & /*x*/,  Vector & g, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global)
    //         {
    //             g = g_diff + (H_diff * s_global);
    //             return true;
    //         }
    // };


    // template<typename Matrix, typename Vector>
    // class MultilevelHessianEval<Matrix, Vector, GALERKIN>
    // {
    //     public:
    //        inline static bool compute_hessian(const ExtendedFunction<Matrix, Vector> & /*fun*/, const Vector & /*x*/, Matrix & H, const Matrix & H_diff)
    //         {
    //             H = H_diff;
    //             return true;
    //         }
    // };



}




#endif //UTOPIA_RMTR_INFTY_HPP