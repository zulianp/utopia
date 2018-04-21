#ifndef TR_ACTIVE_SET_TR_BOX_SUBPROBLEM
#define TR_ACTIVE_SET_TR_BOX_SUBPROBLEM
#include <string>
#include "utopia_TRSubproblem.hpp"
#include "utopia_BoxConstraints.hpp"

namespace  utopia 
{

    /**
     * @brief      Wrapper for Active set method to solve TR subproblems 
     */
    template<class Matrix, class Vector>
    class ActiveSetTRSubproblem : public TRBoxSubproblem<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector) Scalar;

        typedef utopia::LinearSolver<Matrix, Vector>            LinearSolver;
        typedef utopia::TRBoxSubproblem<Matrix, Vector>         TRBoxSubproblem;
        typedef utopia::SemismoothNewton<Matrix, Vector>        SemismoothNewton;
        typedef utopia::Preconditioner<Vector>                  Preconditioner;

        public:
            ActiveSetTRSubproblem(  const std::shared_ptr <LinearSolver> & linear_solver,
                                    const Parameters params = Parameters()):
                TRBoxSubproblem(params), ls_solver_(linear_solver)
            {
                active_set_solver_ = std::make_shared<SemismoothNewton>(ls_solver_); 
            }
            
            virtual ~ActiveSetTRSubproblem( ){}

            virtual bool tr_constrained_solve(const Matrix &H, const Vector &g, Vector &p_k, const BoxConstraints<Vector> & up_constrain) override
            {
                active_set_solver_->set_box_constraints(up_constrain);
                Vector g_minus = -1.0 * g; 
                active_set_solver_->solve(H, g_minus, p_k);

                return true;
            };

        public:
            virtual void set_linear_solver(const std::shared_ptr<LinearSolver > &ls) override
            {
                ls_solver_ = ls; 
                active_set_solver_->set_linear_solver(ls); 
            }                

        private:  
            std::shared_ptr<LinearSolver> ls_solver_;               
            std::shared_ptr<SemismoothNewton> active_set_solver_;           
        
    };
}

#endif //TR_ACTIVE_SET_TR_BOX_SUBPROBLEM
