// #ifndef UTOPIA_NON_LINEAR_FE_FUNCTION_HPP
// #define UTOPIA_NON_LINEAR_FE_FUNCTION_HPP

// #include "utopia_Base.hpp"

// namespace utopia 
// {
//     /**
//      * @brief      Base class for Nonlinear NonLinearFEFunction. All application context needed by solver is usually provided inside of this functions.
//      *             In optimization settings, user needs to supply value(energy), gradient, hessian.
//      *
//      * @todo       Intorduce approximate Hessian updates strategies, e.g. BFGS, ... 
//      * @tparam     Matrix  
//      * @tparam     Vector  
//      */
// //     template<class Matrix, class Vector, class LHS, class RHS, int Backend>
// //     class NonLinearFEFunction  {
// //     // public:
// //     //     DEF_UTOPIA_SCALAR(Matrix)


// //     //     NonLinearFEFunction(const LHS &lhs, const RHS &rhs)
// //     //     : lhs_(lhs), rhs_(rhs)
// //     //     {}

// //     //     virtual ~NonLinearFEFunction() { }

// //     //     virtual bool value(const Vector &/*point*/, Scalar &/*value*/) const 
// //     //     {
// //     //         return 1.;
// //     //     }

// //     //     virtual bool gradient(const Vector &/*point*/, Vector &/*result*/) const
// //     //     {
// //     //         //assemble rhs
// //     //         return false;
// //     //     }


// //     //     virtual bool hessian(const Vector &x, Matrix &H) const
// //     //     {
// //     //         //assemble lhs
// //     //         return false;
// //     //     }

// //     //     virtual bool update(const Vector &/*point*/) { return true; };

// //     //     LHS lhs_;
// //     //     RHS rhs_;
// //     };
// // }
// #endif //UTOPIA_NON_LINEAR_FE_FUNCTION_HPP
