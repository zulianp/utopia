// /*
// * @Author: alenakopanicakova
// * @Date:   2017-04-17
// * @Last Modified by:   Alena Kopanicakova
// * @Last Modified time: 2017-05-05
// */

// #ifndef UTOPIA_TR_SMOOTHER_HPP
// #define UTOPIA_TR_SMOOTHER_HPP

// #include "utopia_NonLinearSmoother.hpp"




// namespace utopia 
// {


//     template<class Matrix, class Vector>
//     class NonLinearTRSmoother<Matrix, Vector> : public NonLinearSmoother<Matrix, Vector> 
//     {
//         typedef UTOPIA_SCALAR(Vector)                           Scalar;
//         typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;
//         typedef utopia::NonLinearSmoother<Matrix, Vector>       Smoother;
//         typedef utopia::TRSubproblem<Matrix, Vector>            TRSubproblem; 
//         typedef utopia::Function<Matrix, Vector>                Function;

//         public:
//         NonLinearTRSmoother(const std::shared_ptr<TRSubproblem> &tr_subproblem = std::shared_ptr<TRSubproblem>(), 
//                             const Parameters params = Parameters()):
//         _tr_subproblem(tr_subproblem), 
//         _eps_smooth()


//         { 
//             set_parameters(params); 
//         }

//         virtual void set_parameters(const Parameters params) override
//         {
//             Smoother::set_parameters(params); 
//         }

//         virtual bool nonlinear_smooth(Function & fun,  Vector &x, const Vector &rhs) override
//         {





//             return true; 
//         }



//     private: 
//         Scalar _eps_smooth; 
//         Vector _Rg; 

//         Scalar _delta0;                     // initial tr radius - restr from fine level 
//         Scalar _delta_fine;                 // fine level TR radius
//         Scalar _delta;                      // current TR radius
//         Scalar _intermediate_delta;         // current TR radius

//         SizeType _l;                        // level 



//         std::shared_ptr<TRSubproblem> _tr_subproblem;  // solving TR subproblems on different levels   






//     };

// }

// #endif //UTOPIA_TR_SMOOTHER_HPP

