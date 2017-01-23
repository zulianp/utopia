// /*
// * @Author: alenakopanicakova
// * @Date:   2016-10-12
// * @Last Modified by:   alenakopanicakova
// * @Last Modified time: 2016-10-04
// */

// #ifndef UTOPIA_TR_SUBPROBLEM_MG_INFTY_HPP
// #define UTOPIA_TR_SUBPROBLEM_MG_INFTY_HPP
// #include "utopia_TRSubproblem.hpp"
// #include "utopia_IterativeSolver.hpp"



// namespace utopia 
// {

// 	template<class Matrix, class Vector>
//     class MGSubproblem : public TRSubproblem<Matrix, Vector>, public MultigridConstrained<Matrix, Vector>
//     {
// 		typedef UTOPIA_SCALAR(Vector) 						Scalar;
//         typedef UTOPIA_SIZE_TYPE(Vector) 					SizeType;

//         typedef utopia::LinearSolver<Matrix, Vector>        Solver;
//         typedef utopia::Smoother<Matrix, Vector>            Smoother;
        

//     public:

//     	MGSubproblem(	const std::shared_ptr<Smoother> &smoother = std::shared_ptr<Smoother>(), 
//                         const std::shared_ptr<Solver> &direct_solver = std::shared_ptr<Solver>()) :
//                         MultigridConstrained<Matrix, Vector>(smoother, direct_solver)                                                             
//         {

//         }

//         ~MGSubproblem(){} 


// 		/*!
//     	\details
// 		Implementation of CG-Steihaug to obtain step p_k 

// 		- please note, that choice of tol_Steihaug and initial radius 
// 		can be crucial and should be properly set up

//     	@note
//     	\param g  		gradient  
//     	\param B 		Hessian 
//     	\param delta 	tr. radius on current iterate 
//     	\param p_k   	new step
//       	*/
//         bool get_pk(const Vector &g, const Matrix &B, const Scalar &delta, Vector &p_k, Scalar &pred) override 
//         {
// 			// Vector lb = -1 * delta * local_values(local_size(g).get(0), 1);
// 			// Vector ub = delta * local_values(local_size(g).get(0), 1);

// 			// this->galerkin_assembly(B);
// 			// this->solve(g, p_k, ub, lb); 
// 			// p_k *= -1; 


// 	  //   	TRSubproblem<Matrix, Vector>::get_pred(g, B, p_k, pred);
// 	    	return true; 
//         }


//     };
// }


// #endif //UTOPIA_TR_SUBPROBLEM_MG_INFTY_HPP