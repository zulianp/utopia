#ifndef TR_BOX_BASED_ON_TAO_SUBPROBLEM
#define TR_BOX_BASED_ON_TAO_SUBPROBLEM

#include <string>

// a lot of problems with includes :( 
#include "utopia_TRBoxSubproblem.hpp"
#include "utopia_QuadraticFunction.hpp"
#include "utopia_petsc_TaoSolver.hpp"

namespace  utopia 
{

    template<class Matrix, class Vector>
    class TaoTRSubproblem : public TRBoxSubproblem<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector) Scalar;

        typedef utopia::LinearSolver<Matrix, Vector>            LinearSolver;
        typedef utopia::TRBoxSubproblem<Matrix, Vector>         TRBoxSubproblem;

        public:
            TaoTRSubproblem(    const std::shared_ptr <LinearSolver> &linear_solver = std::shared_ptr<LinearSolver>(),
                                const Parameters params = Parameters()):
                                TRBoxSubproblem(params)
            {
              solver_ = std::make_shared<utopia::TaoSolver<Matrix, Vector> >( ); 
            }
            
            virtual ~TaoTRSubproblem( ){}


            virtual bool tr_constrained_solve(const Matrix &H, const Vector &g, Vector &p_k, const BoxConstraints<Vector> & box) override
            {
                solver_->set_box_constraints(box);


                // to be fixed ... 
                Matrix * m3 = const_cast<Matrix*> (&H); 
                Vector * v3 = const_cast<Vector*> (&g); 

                *v3  *= -1.0; 


                QuadraticFunction<Matrix, Vector> fun(make_ref( *m3) , make_ref( *v3 ));

                
                //  just for debugging
                
                //  options working nicely: gpcg, BQPIP, pgd,  ----- BNCG - to be checked from rmtr ... 
                solver_->verbose(false);
                solver_->set_type("gpcg");
                solver_->solve(fun, p_k);

                return true;
            };



        protected: 
            std::shared_ptr<utopia::TaoSolver<Matrix, Vector> > solver_; 

    };
}

#endif //TR_BOX_BASED_ON_TAO_SUBPROBLEM
