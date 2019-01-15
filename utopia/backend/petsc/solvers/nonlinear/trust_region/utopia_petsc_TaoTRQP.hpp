#ifndef TR_BOX_BASED_ON_TAO_SUBPROBLEM
#define TR_BOX_BASED_ON_TAO_SUBPROBLEM

#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Traits.hpp"
#include "utopia_TRQuadraticFunction.hpp"

#include <string>
#include <memory>

namespace utopia 
{

    template<class Matrix, class Vector>
    class TaoQPSolver final: public QPSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector) Scalar;

        typedef utopia::LinearSolver<Matrix, Vector>            LinearSolver;

        public:
            TaoQPSolver(const std::shared_ptr<LinearSolver> &linear_solver = std::make_shared<KSPSolver<Matrix, Vector>>("gmres"))
            {
              tao_solver_ = TaoSolver<Matrix, Vector>(linear_solver); 
            }

            TaoQPSolver * clone() const override
            {
                return new TaoQPSolver(*this);
            }
            

            bool solve(const Matrix &H, const Vector &g, Vector &p_k) override
            {
                TRQuadraticFunction<Matrix, Vector, Traits<Vector>::Backend> fun(make_ref(H) , make_ref(g));
                
                configure_tao(); 

                tao_solver_.smooth(fun, p_k);
                
                return true;
            }

            void set_linear_solver(const std::shared_ptr<LinearSolver > &ls)
            {
                tao_solver_->set_linear_solver(ls); 
            }    

            void pc_type(const std::string & pc_type)
            {
                tao_solver_.set_pc_type(pc_type); 
            }


        private:
            void configure_tao()
            {
                // TODO:: investigate reasonable options... 
                tao_solver_.set_type("gpcg");
                
                // default in tao is hudge overshooot.... 
                tao_solver_.atol(this->atol());
                tao_solver_.rtol(this->rtol()); 
                tao_solver_.stol(this->stol());

                // counts + 1 ... 
                tao_solver_.max_it(this->max_it());
                tao_solver_.verbose(this->verbose());


                auto &box = this->get_box_constraints();
                tao_solver_.set_box_constraints(box);                
            }



        protected: 
            TaoSolver<Matrix, Vector> tao_solver_; 

    };
}

#endif //TR_BOX_BASED_ON_TAO_SUBPROBLEM