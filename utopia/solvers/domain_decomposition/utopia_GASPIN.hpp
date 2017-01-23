// /*
// * @Author: alenakopanicakova
// * @Date:   2016-05-11
// * @Last Modified by:   alenakopanicakova
// * @Last Modified time: 2016-06-06
// */
// //    ..........................DO NOT USE THIS FUNCTION RIGHT NOW...................
// #ifndef UTOPIA_GASPIN_HPP
// #define UTOPIA_GASPIN_HPP

// // #include "tr_params.hpp"
// #include "utopia_TRSubproblem.hpp"
// #include "utopia_CauchyPoint.hpp"
// #include "utopia_TR_base.hpp"
// #include "utopia_TR_LocalSolve.hpp" 
// #include "utopia_GLFunction.hpp"

// #include <iomanip>
// #include <iostream>
// #include <limits>
// #include <algorithm>


// namespace utopia 
// {
//     template<class GlobalMatrix, class GlobalVector, class LocalMatrix, class LocalVector>
//     //! GASPIN solver 
//     /*!
//         Globalized aspin solver 
//         NOTE: - linear solvers are required
//     */
//     class GASPIN :  public TrustRegionBase<GlobalMatrix, GlobalVector> 
//     {
//       typedef typename utopia::Traits<GlobalVector>::Scalar Scalar;
//       typedef utopia::LinearSolver<LocalMatrix, LocalVector> Local_linear_solver;
//       typedef utopia::LinearSolver<GlobalMatrix, GlobalVector> Global_linear_solver;
//       typedef utopia::TRSubproblem<GlobalMatrix, GlobalVector> TRSubproblem; 
//       typedef typename utopia::Traits<LocalVector>::SizeType SizeType;


//     public:
//       GASPIN( Local_TR<GlobalMatrix, GlobalVector, LocalMatrix, LocalVector> local_TR, 
//               const std::shared_ptr<Global_linear_solver> &global_linear_solver = std::shared_ptr<Global_linear_solver>(), 
//               const std::shared_ptr<Local_linear_solver> &local_linear_solver = std::shared_ptr<Local_linear_solver>()):   
      
//       TrustRegionBase<GlobalMatrix, GlobalVector>(global_linear_solver),
//       global_linear_solver(global_linear_solver), 
//       local_linear_solver(local_linear_solver),
//       local_TR(local_TR) 
//       {}


//       bool solve(GLFunction<GlobalMatrix, GlobalVector, LocalMatrix, LocalVector> &fun, GlobalVector &x_k) 
//       {
//         using namespace utopia;

// /////////////////////////////////////////////////////////////////////////////////////////////////////////////
// /// --------------------------------------------- GASPIN  ---------------------------------------------------
// ///   symbols without subscript are meant to represent global objects, 
// ///   while variables with subscript k represent objects on k-th iterate 
// ///   e.g.  s = \sum I_k s_k 
// ///   subscript s stands for iteration               
// /////////////////////////////////////////////////////////////////////////////////////////////////////////////

//     //     Scalar it = 0, rad_flg, local_red;
//     //     Scalar g_norm = std::numeric_limits<Scalar>::infinity(); 

//     //     GlobalVector g, p_CP = x_k, dx = x_k, p_N = x_k, p_k = x_k, x_k1 = x_k, count_domains, s, g_precond;
//     //     GlobalMatrix H, C; 



//     //     LocalVector s_k; 

//     //     // passing solver and parameters into subproblem 
//     //   //  this->tr_subproblem->init(this->linear_solver(), params); 
      


//     //      fun.hessian(x_k, H); 
//     //      fun.gradient(x_k, g);


//     //     // #define DEBUG_mode
//     //     // #define LS_check


//     //     // TR delta initialization
//     //     delta = params.delta0(); 
//     //     rad_flg = this->delta_init(x_k ,delta);
//     //     radius_L = delta;  


//     //     // test with which kind of matrix we are working 
//     //     //  mat_tests.test(H, x_k); 


//     //     // just to start  - do this as const int - bla-bla -bla - count subdomains 
//     //     count_domains = localValues(1,1.0);
//     //     Scalar domains = sum(count_domains);
//     //     this->info().GASPIN_init_message(domains); 


//     //     Scalar dummy = 0; 
//     //     g_norm = norm2(g);

//     //     // print out - just to have idea how are we starting 
//     //     #ifdef DEBUG_mode
//     //       this->info().TR_iter_status_debug(it, g_norm, dummy, dummy,  dummy , dummy, dummy, dummy,  dummy, dummy, dummy, dummy);
//     //     #else                 
//     //       this->info().TR_iter_status(it, g_norm, dummy, dummy, dummy, dummy); 
//     //     #endif



//     //     // found out if there is a linear solution - or start with the newton step 
//     //     #ifdef LS_check
//     //       Scalar ls_flag = this->linear_solution_check(fun, g, H, p_N, x_k);
//     //       if(ls_flag)
//     //       {
//     //         this->info().exitSolver(1, g_norm, SolverInfo::CONVERGED_FNORM_RELATIVE);
//     //         return true; 
//     //       }
//     //     #endif

//     //     it++; 
//     //     fun.value(x_k, E); 

//     //     while(it < this->info().max_it())
//     //     {
//     //       x_k = x_k1; 
//     //       E_k = E; 

//     //       fun.hessian(x_k, H); 

//     //       local_red = local_TR.solve(fun,H, x_k, g, s_k, radius_L); 
//     //       fun.interpolation(s_k, s);                    // assign local corrections into glob. vector 

//     //       fun.build_preconditioner(x_k, C);             // build inverse of additive schwarz preconditioner - no overlap
//     //       recombination(g, C, s, g_precond, radius_L);


//     // //----------------------------------------------------------------------------
//     // //     new step p_k w.r. ||p_k|| <= delta
//     // //----------------------------------------------------------------------------
//     //       this->tr_subproblem->get_pk(g_precond, H, delta, p_k, pred); 
//     // //----------------------------------------------------------------------------
//     // //----------------------------------------------------------------------------

//     //       if(it == 1 && rad_flg)
//     //       {
//     //         delta = norm2(p_k);
//     //         delta *= 0.2; 
//     //       }


//     //       // value of the objective function with correction 
//     //       fun.value(x_k + p_k, E_k1);
//     //       product = dot(g, p_k);  // just to do tests 


//     //       // decrease ratio 
//     //       ared = E_k - E_k1;        // reduction observed on objective function
//     //       rho = ared/ (2 * pred);   // decrease ratio 

//     // //----------------------------------------------------------------------------
//     // //     acceptance of trial point 
//     // //     MISSING SDC2 CHECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
//     // //     and both CAUCHY POINTs comparisons !!!!!!!!!!!!!!!!!
//     // //----------------------------------------------------------------------------
//     //       if(ared < 0 || pred < 0 || ared == pred)
//     //       {
//     //         rho = 0; 
//     //       } 

//     //       this->trial_point_acceptance(rho, E, E_k, E_k1, p_k, x_k, x_k1);

//     // //----------------------------------------------------------------------------
//     // //    convergence check 
//     // //----------------------------------------------------------------------------
//     //       fun.gradient(x_k1, g); 
//     //       g_norm = norm2(g); 

//     //       #ifdef DEBUG_mode
//     //         this->info().TR_iter_status_debug(it, g_norm, product, E_k,  E_k1, E, ared, pred,  rho, delta, dummy, dummy);
//     //       #else                 
//     //         this->info().TR_iter_status(it, g_norm, E_k, E, rho, delta); 
//     //       #endif
                

//     //       // termination because norm of grad is down
//     //       if(g_norm < this->info().tol())
//     //       {
//     //         x_k = x_k1; 
//     //         this->info().exitSolver(it, g_norm, SolverInfo::CONVERGED_FNORM_RELATIVE);
//     //         return true; 
//     //       }

//     //       // step size so small that we rather exit than wait for nan's
//     //       g_norm = norm2(p_k);
//     //       if(g_norm < params.eps())
//     //       {
//     //         x_k = x_k1; 
//     //         this->info().exitSolver(it, g_norm, SolverInfo::CONVERGED_SNORM_RELATIVE);
//     //         return true; 
//     //       }

//     //       // radius is super small - there is no hope for TR to converge ever
//     //       g_norm = norm2(x_k1); 
//     //       if(delta < params.eps())
//     //       {
//     //         x_k = x_k1; 
//     //         this->info().exitSolver(it, g_norm, SolverInfo::CONVERGED_TR_DELTA);
//     //         return true; 
//     //       }


//     // //----------------------------------------------------------------------------
//     // //      tr. radius update 
//     // //----------------------------------------------------------------------------

//     //       // local radius update 
//     //       radius_L_hat = params.gamma2() * radius_L; 
//     //       radius_L = std::min(delta, radius_L_hat); 


//     //       if(rho < params.eta1())
//     //       {
//     //         delta = params.gamma1() * delta;
//     //       }
          
//     //       else if (rho > params.eta2()) 
//     //       {
//     //         delta = std::min(params.gamma2() * delta, params.delta_max()); 
//     //       }      

//     //       it++; 

//     //     }

//     //       x_k = x_k1; 
//     //       this->info().exitSolver(it, g_norm, SolverInfo::DIVERGED_MAX_IT);
//           return false;


//       }


// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////        

//     /*!
//     \details
//     set up strategy

//     @note
//     \param subproblem            - choice of subproblem 
//       */
//     void strategy(const std::shared_ptr<TRSubproblem> &subproblem = std::shared_ptr<TRSubproblem>())
//     {
//       tr_subproblem = subproblem; 
//     }




//     /*!
//     \details
//     Parametersinit 

//     @note
//     \param Parameters                - initialization Parameterscomming from moose 
//       */
//     // bool params_init(const TR_Parameters& parameters)
//     // {
//     //   params.set_allparams(parameters); 
//     //   return true; 
//     // }



// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////        
//     private:

//       /*!
//       \ details
//       recombination approach from G-ASPIN paper 
//       - check out minuses - if correct          

//       \f$  \alpha_i = \min \{ 1, \max \{ 0, 1 - \frac{\Delta_i^L}{||g_i|| + || C_i   * \sum_k I_k s_k|| \} \} \} \\
//       \tilde{g}_i = \alpha_i * g_i - (1 - \alpha_i) * C_i \sum_k I_k s_k \f$  

//       @note
//       \param g                  - gradient
//       \param C                  - preconditioner
//       \param s                  - local corrections
//       \param g_precond          - preconditioned gradient 
//       \param delta_L            - local TR radius 
//       */
//       template<class Matrix, class Vector>
//       bool recombination( const Vector& g, const Matrix& C, const Vector& s, Vector & g_precond, const Scalar delta_L)
//       {

//         Scalar g_norm = norm2(g); 
//         Vector C_corr =  C * s; 
//         Scalar corr_length = norm2(C_corr); 
//         Scalar term2 = 1 - delta_L/(g_norm + corr_length);


//         Scalar max_term =  (0 > term2)? 0 : term2; 
//         Scalar alpha = (1 < max_term) ? 1 : max_term;  


//         std::cout<<"alpha: "<< alpha << "   -->>   \n"; 
//         g_precond = alpha * g - (1 - alpha) * C_corr; 

//         return true; 
//       }




//       /*!
//       \ details
//       - combination approach to get precond gradient 
//       - used for the first approach from G-ASPIN paper 

//       @note
//       \param C                  - preconditioner
//       \param s                  - local corrections
//       \param g_precond          - preconditioned gradient 
//       */
//       template<class Matrix, class Vector>
//       bool combination(const Matrix& C, const Vector& s, Vector & g_precond)
//       {

//         Vector C_corr =  C * s; 
//         g_precond = -1 *  C_corr; 

//         return true; 
//       }

// //----------------------------------------------------------------------------------------------------------

//       // LS solvers 
//       std::shared_ptr<Global_linear_solver> global_linear_solver;
//       std::shared_ptr<Local_linear_solver> local_linear_solver;


//       Scalar radius_L = 999999999999; 
//       Scalar radius_L_hat = 999999999; 

//       // check this out 
//      // TRParametersparams; 

//       Scalar delta = 0, product = 0;                   
//       Scalar ared, pred; 
//       Scalar rho; 
//       Scalar E, E_k, E_k1; 

//       //Mat_tests<GlobalMatrix, GlobalVector > mat_tests;  
//       std::shared_ptr<TRSubproblem> tr_subproblem;

//       Local_TR<GlobalMatrix, GlobalVector, LocalMatrix, LocalVector> local_TR;       /*!< local TR   */  




//   };

//     }
// #endif //UTOPIA_GASPIN_HPP
