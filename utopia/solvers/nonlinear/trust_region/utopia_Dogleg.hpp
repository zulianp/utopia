// /*
// * @Author: alenakopanicakova
// * @Date:   2016-05-11
// * @Last Modified by:   alenakopanicakova
// * @Last Modified time: 2016-10-14
// */

// #ifndef UTOPIA_TR_SUBPROBLEM_DOGLEG_HPP
// #define UTOPIA_TR_SUBPROBLEM_DOGLEG_HPP
// #include "utopia_TRSubproblem.hpp"
// #include "utopia_CauchyPoint.hpp"
// #include "utopia_Parameters.hpp"    


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// NEEEEEDS PROPER INTERFACE ///////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// namespace utopia 
// {

// 	template<class Matrix, class Vector>
//     class Dogleg : public TRSubproblem<Matrix, Vector>
//     {
//         typedef UTOPIA_SCALAR(Vector) Scalar;

//     public:

//         Dogleg(const Parameters & params = Parameters()): 
//                                                           TRSubproblem<Matrix, Vector>(params)
//         {
//             // set_parameters(); 
//         };
		
//         // THIS FUNCTION GOT A BIT BROKEN BY MOVING FROM PASSO TO UTOPIA ----- > TODO:: fix it .... 
		
//         /*!
//     	\details
		
// 		Dogleg alg. to obtain \f$ p_k \f$ step, wrt. \f$ ||p_k|| \leq \Delta_k \f$

//     	@note
//     	\param g  		gradient  
//     	\param B 		Hessian 
//     	\param delta 	tr. radius on current iterate 
//     	\param p_k   	new step
//       	*/
//       //   bool get_pk(const Vector &g, const Matrix &B, const Scalar &delta, Vector &p_k, Scalar &pred) override 
//       //   {

//       //       // check if LS was initialized 
//       //       LS_check(); 
//       //       Vector p_N = g, p_SD = g, p_CP = g;
//       //       Scalar pred_N, g_B_g = dot(g, B * g), pNlen; 

//       //       // newton direction 
//       //       this->linear_solver_->solve(B, -1 * g, p_N);


//       //       // model reduction for newton point
//       //       this->get_pred(g, B, p_N, pred_N); 


//       //       // good situation: p_N reduces m_k(p)
// 	    	// if(pred_N >= 0)
//       //       {   
//       //           pNlen = norm2(p_N);
//       //           // newton point inside of TR
//       //           if(pNlen <= delta)
//       //           {
//       //               p_k = p_N;
//       //               this->get_pred(g, B, p_k, pred); 
//       //           }
//       //           else
//       //           {
//       //               // check if CP or SD 
//       //               //steepest descent step 
//       //               p_k = (dot(g, g)/g_B_g) * g;
//       //               p_k *= -1.0;
//       //               Scalar SD_norm = norm2(p_k); 

//       //               if(SD_norm >= delta)
//       //               {
//       //                   // p_CP lies on boundary => take it 
//       //                   // p_k = -1 * p_CP;
                        
//       //               }
//       //               else
//       //               {
//       //                   // p_SD inside, p_N outside => dogleg 
//       //                   // p_k = p_CP + \tau (p_N - P_CP), 
//       //           		// where \tau is the largest value in [0,1]
//       //                   // such that || p_k || <= \delta                   	
//       //                   Vector d = p_N - p_SD; 
//       //                   this->quad_solver(p_CP, d, delta,  p_k);
//       //                   this->get_pred(g, B, p_k, pred); 
//       //               }
//       //           }


//       //       }

//       //       //  bad situation: p_N increases m_k(p)
//       //       else
//       //       {
//       //           if(g_B_g > 5 * std::sqrt(2e-16))
//       //           {
//       //               //  1D function is convex enough
//       //               //  find unconstrained SD minimizer for model fun
//       //               p_SD = -(dot(g, g)/g_B_g) * g;
//       //               Scalar pSDlen = norm2(p_SD);

//       //               if(pSDlen <= delta)
//       //               {
//       //                   p_k = p_SD;
//       //                   this->get_pred(g, B, p_k, pred); 
//       //               }
//       //               else
//       //               {
//       //                   p_k = delta * (p_SD* (1/pSDlen));
//       //                   this->get_pred(g, B, p_k, pred);  
//       //               }
//       //           }
//       //           else
//       //           {
//       //               //  1D slice of model function in gradient direction is concave:
//       //               //  SD minimizer lies on TR boundary, in direction of -g'
//       //               //cauchyPoint(J, J_T, r, delta, p_CP);
//       //               cauchy_point.get_pk(g, B, delta, p_k, pred); 
//       //               p_k = p_CP; 
//       //           }
//       //       }

// 	    	// return true; 
//       //   }




//         /*!
//         \details
        
//         Dogleg alg. to obtain \f$ p_k \f$ step, wrt. \f$ ||p_k|| \leq \Delta_k \f$

//         @note
//         \param g        gradient  
//         \param B        Hessian 
//         \param delta    tr. radius on current iterate 
//         \param p_k      new step
//         */


//         /////////////////////////////////////////////////////////////////////////////////////////////////////
//         bool get_pk(const Vector &g, const Matrix &B, const Scalar &delta, Vector &p_k, Scalar &pred) override 
//         {

//             // // check if LS was initialized 
//             // LS_check(); 
//             // Vector p_N = g, p_SD = g, p_CP = g;
//             // Scalar pred_N, g_B_g = dot(g, B * g), p_N_norm; 

//             // // newton direction 
//             // this->linear_solver_->solve(B, -1 * g, p_N);
//             // p_N_norm = norm2(p_N);

//             // // newton point inside of TR
//             // if(p_N_norm <= delta)
//             // {
//             //     p_k = p_N;
//             //     // model reduction for newton point
//             //     this->get_pred(g, B, p_k, pred); 
//             // }
//             // else
//             // {
//             //     // dogleg part ... 
//             //     p_SD = (dot(g, g)/g_B_g) * g;
//             //     p_SD *= -1.0;
//             //     Scalar SD_norm = norm2(p_SD); 

                        
//             //     Scalar p_NN = dot(p_N, p_N);
//             //     Scalar p_ND = dot(p_N, p_SD);
//             //     Scalar p_DD = dot(p_SD, p_SD);

//             //     Scalar a = p_NN - 2.0 * p_ND + p_DD;
//             //     Scalar b = -2.0 * p_DD + 2.0 * p_ND;
//             //     Scalar c = p_DD - delta * delta;


//             //     Scalar discriminant = b*b - 4.0*a*c;

//             //     if(discriminant < 0.0 && fabs(discriminant) < 1e-16)
//             //         discriminant = 0.0;

//             //     if (discriminant < 0.0) //TODO::FIXME
//             //     {
//             //         std::cout << "Trust Region: negative discriminant of dogleg = " << discriminant << ", using Newton direction."<< std::endl;
//             //         Scalar pn_norm = norm2(p_N); 
//             //         p_k = (delta/pn_norm) * p_N;

//             //         this->get_pred(g, B, p_k, pred); 
//             //         return true ;
//             //     }

//             //     Scalar x1 = (-b + std::sqrt(discriminant))/(2 * a);
//             //     Scalar x2 = (-b - std::sqrt(discriminant))/(2 * a);

//             //     Scalar tau = std::max(x1, x2);
//             //     p_k = (1 - tau) * p_SD + tau * p_N; 
//             //     this->get_pred(g, B, p_k, pred); 
//             // }

//             return true; 
//         }
//         /////////////////////////////////////////////////////////////////////////////////////////////////////




//         // /*!
//         // \details
        
//         // Dogleg alg. to obtain \f$ p_k \f$ step, wrt. \f$ ||p_k|| \leq \Delta_k \f$
//         // - imlpementation for normal equation

//         // @note
//         // \param g        gradient  
//         // \param B        Hessian 
//         // \param delta    tr. radius on current iterate 
//         // \param p_k      new step
//         // */
//         // bool get_pk(const Vector &r, const Vector &g, const Matrix &J, const Matrix &B, const Scalar &delta, Vector &p_k, Scalar &pred) override
//         // {

//         //     // check if LS was initialized 
//         //     LS_check(); 
            
//         //     Vector p_N = g, p_SD = g, p_CP = g;

//         //     // newton direction 
//         //     this->linear_solver_->solve(J, -1 * r, p_N);

//         //     Scalar pNlen; 

//         //     pNlen = norm2(p_N);
//         //     // newton point inside of TR
//         //     if(pNlen <= delta)
//         //     {
//         //         p_k = p_N;
//         //         this->get_pred(g, B, p_k, pred); 
//         //     }
//         //     else
//         //     {
//         //         // check if CP or SD 
//         //         cauchy_point.get_pk(g, B, delta, p_k, pred); 
//         //         Scalar CP_norm = norm2(p_CP); 

//         //         if(CP_norm >= delta)
//         //         {
//         //             // p_CP lies on boundary => take it 
//         //             p_k = p_CP;
//         //         }
//         //         else
//         //         {
//         //             // p_SD inside, p_N outside => dogleg 
//         //             // p_k = p_CP + \tau (p_N - P_CP), 
//         //             // where \tau is the largest value in [0,1]
//         //             // such that || p_k || <= \delta                    
//         //             Vector d = p_N - p_SD; 
//         //             this->quad_solver(p_CP, d, delta,  p_k);
//         //             this->get_pred(g, B, p_k, pred); 
//         //         }
//         //     }

//         //     return true; 
//         // }







//         /*!
//         \details
        
//         Dogleg alg. to obtain \f$ p_k \f$ step, wrt. \f$ ||p_k|| \leq \Delta_k \f$
//         - imlpementation for normal equation

//         @note
//         \param g        gradient  
//         \param B        Hessian 
//         \param delta    tr. radius on current iterate 
//         \param p_k      new step
//         */
//         /////////////////////////////////////////////////////////////////////////////////////////////////////
//         bool get_pk(const Vector &r, const Vector &g, const Matrix &J, const Matrix &B, const Scalar &delta, Vector &p_k, Scalar &pred) override
//         {

//             // // check if LS was initialized 
//             // LS_check(); 
            
//             // Vector p_N = g, p_SD = g, p_CP = g;
//             // Scalar pred_N, g_B_g = dot(g, B * g), p_N_norm; 

//             // // newton direction 
//             // this->linear_solver_->solve(J, -1 * r, p_N);

//             // // newton direction 
//             // p_N_norm = norm2(p_N);

//             // // newton point inside of TR
//             // if(p_N_norm <= delta)
//             // {
//             //     p_k = p_N;
//             //     // model reduction for newton point
//             //     this->get_pred(g, B, p_k, pred); 
//             // }
//             // //else
//             // // {
//             // //   p_k = p_N * (delta/ p_N_norm ); 
//             // //   this->get_pred(g, B, p_k, pred); 

//             // // }
//             // else
//             // {
//             //     // dogleg part ... 
//             //     p_SD = (dot(g, g)/g_B_g) * g;
//             //     p_SD *= -1.0;
//             //     Scalar SD_norm = norm2(p_SD); 

                        
//             //     Scalar p_NN = dot(p_N, p_N);
//             //     Scalar p_ND = dot(p_N, p_SD);
//             //     Scalar p_DD = dot(p_SD, p_SD);

//             //     Scalar a = p_NN - 2.0 * p_ND + p_DD;
//             //     Scalar b = -2.0 * p_DD + 2.0 * p_ND;
//             //     Scalar c = p_DD - delta * delta;


//             //     Scalar discriminant = b*b - 4.0*a*c;

//             //     if(discriminant < 0.0 && fabs(discriminant) < 1e-16)
//             //         discriminant = 0.0;

//             //     if (discriminant < 0.0) //TODO::FIXME
//             //     {
//             //         std::cout << "Trust Region: negative discriminant of dogleg = " << discriminant << ", using Newton direction."<< std::endl;
//             //         Scalar pn_norm = norm2(p_N); 
//             //         p_k = (delta/pn_norm) * p_N;

//             //         this->get_pred(g, B, p_k, pred); 
//             //         return true ;
//             //     }

//             //     Scalar x1 = (-b + std::sqrt(discriminant))/(2 * a);
//             //     Scalar x2 = (-b - std::sqrt(discriminant))/(2 * a);

//             //     Scalar tau = std::max(x1, x2);
//             //     p_k = (1 - tau) * p_SD + tau * p_N; 
//             //     this->get_pred(g, B, p_k, pred); 
//             // }

//             return true; 
//         }
// /////////////////////////////////////////////////////////////////////////////////////////////////////





//         /*!
//         \details
//             Quick check if LS was initialized - needed in order to obtain p_N
//         */
//         bool LS_check()
//         {
//             // quite brute-force, figure out what's wrong in factory ... 
            
// ////////////////
//             ////////////////if(this->linear_solver_ == NULL)
//             ////////////////{
//               ////////////////auto cg = std::make_shared< ConjugateGradient<Matrix, Vector> >();
//               ////////////////cg->verbose(false); 
//               ////////////////this->linear_solver_ = cg; 
//             ////////////////}
//             return true; 
//         }




//         /*!
//         \details
//             Quick check if LS was initialized - needed in order to obtain p_N
//         */
//         bool set_LS( const std::shared_ptr <utopia::LinearSolver<Matrix, Vector>> &linear_solver = std::shared_ptr<utopia::LinearSolver<Matrix, Vector>>())
//         {
//               //////////////// this->linear_solver_ = linear_solver;
//             return true; 
//         }


//     private: 
//         CauchyPoint<Matrix, Vector> cauchy_point;       /*!< Cauchy Point  */  


//     };
// }


// #endif //UTOPIA_TR_SUBPROBLEM_DOGLEG_HPP

