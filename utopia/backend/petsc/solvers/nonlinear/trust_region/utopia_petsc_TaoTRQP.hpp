#ifndef TR_BOX_BASED_ON_TAO_SUBPROBLEM
#define TR_BOX_BASED_ON_TAO_SUBPROBLEM

#include "utopia_TRBoxSubproblem.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Traits.hpp"
#include "utopia_TRQuadraticFunction.hpp"

#include <string>
#include <memory>

namespace utopia {

    template<class Matrix, class Vector>
    class TaoTRSubproblem : public TRBoxSubproblem<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector) Scalar;

        typedef utopia::LinearSolver<Matrix, Vector>            LinearSolver;
        typedef utopia::TRBoxSubproblem<Matrix, Vector>         TRBoxSubproblem;

        public:
            TaoTRSubproblem(const std::shared_ptr<LinearSolver> &linear_solver = std::make_shared<KSPSolver<Matrix, Vector>>("gmres"),
                            const Parameters params = Parameters()):
                            TRBoxSubproblem(params)
            {
              
            }

            TaoTRSubproblem * clone() const override
            {
                return new TaoTRSubproblem(std::shared_ptr<LinearSolver>(linear_solver_->clone()));
            }
            
            virtual ~TaoTRSubproblem( ){}

            virtual bool tr_constrained_solve(const Matrix &H, const Vector &g, Vector &p_k, const BoxConstraints<Vector> & box) override
            {
                // TODO:: if we store and re-initialize solver, there is a lot of memory leaks comming 
                utopia::TaoSolver<Matrix, Vector> tao_solver_(linear_solver_); 
                tao_solver_.set_box_constraints(box);

                TRQuadraticFunction<Matrix, Vector, Traits<Vector>::Backend> fun(make_ref(H) , make_ref(g));
                
                //  TODO:: investigate suitable options: gpcg, tron, ...
                tao_solver_.set_type("gpcg");
                
                // this is already way much better than even neccessary... 
                tao_solver_.atol(1e-11);
                tao_solver_.rtol(1e-11); 
                tao_solver_.stol(1e-11);

                tao_solver_.solve(fun, p_k);
                
                return true;
            }

            virtual void set_linear_solver(const std::shared_ptr<LinearSolver > &ls) override
            {
                linear_solver_ = ls; 
            }    

        protected: 
            std::shared_ptr<LinearSolver> linear_solver_; 

    };
}

#endif //TR_BOX_BASED_ON_TAO_SUBPROBLEM