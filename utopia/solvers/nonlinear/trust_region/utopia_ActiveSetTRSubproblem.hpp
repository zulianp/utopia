/*
* @Author: alenakopanicakova
* @Date:   2017-06-15
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-02
*/
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
            ActiveSetTRSubproblem(  const std::shared_ptr <LinearSolver> &linear_solver = std::shared_ptr<LinearSolver>(),
                                    const Parameters params = Parameters()):
                TRBoxSubproblem(params)
            {
                _active_set_solver = std::make_shared<SemismoothNewton>(linear_solver); 
            }
            
            virtual ~ActiveSetTRSubproblem( ){}


            virtual bool tr_constrained_solve(const Matrix &H, const Vector &g, Vector &p_k, const BoxConstraints<Vector> & up_constrain) override
            {
                _active_set_solver->set_box_constraints(up_constrain);
                
                Vector g_minus = -1.0 * g; 

                // just for debugging
                _active_set_solver->verbose(true);
                _active_set_solver->solve(H, g_minus, p_k);

                return true;
            };



        protected: 
            std::shared_ptr<SemismoothNewton> _active_set_solver; 



        
    };
}

#endif //TR_ACTIVE_SET_TR_BOX_SUBPROBLEM
