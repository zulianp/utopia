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
            TaoTRSubproblem(const std::shared_ptr<LinearSolver> &linear_solver = std::make_shared<KSPSolver<Matrix, Vector>>(),
                            const Parameters params = Parameters()):
                            TRBoxSubproblem(params)
            {
              tao_solver_ = std::make_shared<utopia::TaoSolver<Matrix, Vector> >(); 
            }

            TaoTRSubproblem * clone() const override
            {
                return new TaoTRSubproblem(std::shared_ptr<LinearSolver>(linear_solver_->clone()));
            }
            
            virtual ~TaoTRSubproblem( ){}

            virtual bool tr_constrained_solve(const Matrix &H, const Vector &g, Vector &p_k, const BoxConstraints<Vector> & box) override
            {
                tao_solver_->set_box_constraints(box);

                // note, we do not need to switch sign ... 
                TRQuadraticFunction<Matrix, Vector, Traits<Vector>::Backend> fun(make_ref(H) , make_ref(g));
                
                //  suitable options: gpcg, BQPIP, pgd, bncg, bqpip - to be checked from rmtr ... 
                tao_solver_->set_type("gpcg");
                tao_solver_->solve(fun, p_k);
                
                return true;
            }

            virtual void set_linear_solver(const std::shared_ptr<LinearSolver > &ls) override
            {
                linear_solver_ = ls; 
                tao_solver_->set_linear_solver(ls); 
            }    

        protected: 
            std::shared_ptr<utopia::TaoSolver<Matrix, Vector> > tao_solver_; 
            std::shared_ptr<LinearSolver> linear_solver_; 

    };
}

#endif //TR_BOX_BASED_ON_TAO_SUBPROBLEM