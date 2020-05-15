// #ifndef UTOPIA_NOX_NEWTON_HPP
// #define UTOPIA_NOX_NEWTON_HPP

// #include "utopia_Core.hpp"
// #include "utopia_LinearSolver.hpp"
// #include "utopia_Function.hpp"
// #include "utopia_NonLinearSolver.hpp"
// #include "utopia_Newton.hpp"

// //// NOX headers

// ////

// #include <iomanip>
// #include <limits>

// namespace utopia
// {

//     template<class Matrix, class Vector>
//     class Newton<Matrix, Vector, TRILINOS> : public SNESSolver<Matrix, Vector>
//     {
//         typedef UTOPIA_SCALAR(Vector)                           Scalar;
//         typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;

//         typedef utopia::SNESSolver<Matrix, Vector>                  SNESSolver;
//         typedef utopia::LSStrategy<Matrix, Vector>                  LSStrategy;
//         typedef typename NonLinearSolver<Matrix, Vector>::Solver    LinearSolver;

//         public:
//         Newton( const std::shared_ptr <LinearSolver> &linear_solver = std::shared_ptr<LinearSolver>(),
//                 const Parameters params = Parameters(),
//                 const Scalar & alpha = 1.0,
//                 const SizeType & order = 3,
//                 const std::vector<std::string> ls_types    = {"basic", "bt", "l2", "cp", "nleqerr"}) :

//                     SNESSolver(linear_solver, params), alpha_(alpha), line_search_types_(ls_types)
//         {
//             line_search_type_ = line_search_types_.at(0);
//             set_parameters(params);
//             this->set_snes_type("newtonls");

//             line_search_order_ = (order < 4 && order > 0) ?  order : 3;
//         }

//         virtual void set_parameters(const Parameters params) override
//         {
//             SNESSolver::set_parameters(params);
//         }

//         virtual void set_line_search_type(const std::string & ls_type)
//         {
//             line_search_type_ = in_array(ls_type, line_search_types_) ? ls_type : line_search_types_.at(0);;
//         }

//         virtual void set_line_search_order(const SizeType & order)
//         {
//             line_search_order_ = (order < 4 && order > 0) ?  order : 3;
//         }

//     protected:
//         virtual void set_snes_options(SNES & snes,  const Scalar & atol     = SNESSolver::atol(),
//                                                     const Scalar & rtol     = SNESSolver::rtol(),
//                                                     const Scalar & stol     = SNESSolver::stol(),
//                                                     const SizeType & max_it = SNESSolver::max_it()) override
//         {
//             SNESSolver::set_snes_options(snes, atol, rtol, stol, max_it);

//             SNESLineSearch linesearch;
//             SNESGetLineSearch(snes, &linesearch);
//             SNESLineSearchSetType(linesearch,   line_search_type_.c_str() );
//             SNESLineSearchSetOrder(linesearch, line_search_order_);
//             SNESLineSearchSetDamping(linesearch, alpha_);
//         }

//     private:
//         Scalar alpha_;
//         const std::vector<std::string> line_search_types_;
//         std::string line_search_type_;
//         SizeType line_search_order_;

//     };

// }
// #endif //UTOPIA_NOX_NEWTON_HPP
