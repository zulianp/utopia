/*
* @Author: alenakopanicakova
* @Date:   2016-04-07
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-02
*/
#ifndef TR_SUBPROBLEM
#define TR_SUBPROBLEM
#include <string>
#include "utopia_IterativeSolver.hpp"

namespace  utopia 
{

    /**
     * @brief      Wrapper for TR subproblems 
     */
    template<class Matrix, class Vector>
    class TRSubproblem : public IterativeSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector) Scalar;
        typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;
        typedef utopia::Preconditioner<Vector> Preconditioner;

        public:
            TRSubproblem(const Parameters params = Parameters())
            {
                set_parameters(params); 
            };
            
            virtual ~TRSubproblem( ){}        

            /**
             * @brief      Sets the parameters.
             *
             * @param[in]  params  The parameters
             */
            virtual void set_parameters(const Parameters params) override
            {
                IterativeSolver::set_parameters(params); 
                current_radius(params.delta0()); 
            } 

            /**
             * @brief      Setter for current radius. 
             *
             * @param[in]  radius  The radius
             */
            virtual void current_radius(const Scalar &radius)
            {
                current_radius_ = radius; 
            }; 
            
            /**
             * @brief      Getter for current radius. 
             */
            virtual Scalar current_radius()
            {
                return current_radius_;
            }; 


        protected:
            /**
             * @brief      Solver solves easy problem of finding minimum, 
             *             such that \f$ p_k = s + \tau d \f$ minimizes \f$ m_k(p_k) \f$
             *             and satisfies \f$ ||p_k|| = \Delta_k  \f$
             *
             * @param[in]  s      
             * @param[in]  d      
             * @param[in]  delta  The TR radius. 
             * @param      p_k    The current iterate. 
             *
             * @return     tau
             */
            Scalar quad_solver(const Vector &s, const Vector &d, const Scalar & delta,  Vector &p_k)
            {
                Scalar a, b, c, x1, x2, nom, denom,tau;
                
                a = dot(d, d);
                b = dot(2 * s, d); 
                c = dot(s, s) - delta * delta;

                nom = b * b - 4 * a * c; 
                nom =  std::sqrt(nom);
                denom = 2 * a; 

                x1 = (- b + nom)/denom; 
                x2 = (- b - nom)/denom; 

                tau = std::max(x1, x2); 

                if(tau != tau)
                    tau = 0; 

                p_k = s + tau * d; 
                return tau; 
            }


        public: 
            virtual bool tr_constrained_solve(const Matrix &H, const Vector &g, Vector &p_k)
            {
                update(make_ref(H));
                apply(g, p_k); 
                return true; 
            }



    private:

        virtual bool unpreconditioned_solve(const Matrix &/*B*/, const Vector &/*g*/, Vector &/*p_k*/){ return false; };
        virtual bool preconditioned_solve(const Matrix &/*B*/, const Vector &/*g*/, Vector &/*p_k*/){ return false; };

        /**
         * @brief                Solution routine for CG. 
         *
         * @param[in]  b         The right hand side. 
         * @param      x         The initial guess/solution.
         *
         * @return true if the linear system has been solved up to required tollerance. False otherwise
         */
         virtual bool apply(const Vector &b, Vector &x) override
         {
            if(precond_) {
                return preconditioned_solve(*this->get_operator(), b, x);
            } else {
                return unpreconditioned_solve(*this->get_operator(), b, x);
            }
         }

         /**
          * @brief      Sets the preconditioner.
          *
          * @param[in]  precond  The preconditioner.
          */
         virtual void set_preconditioner(const std::shared_ptr<Preconditioner> &precond)
         {
             precond_ = precond;
         }

         /*! @brief if overriden the subclass has to also call this one first.
          */
         virtual void update(const std::shared_ptr<const Matrix> &op) override
         {
             IterativeSolver::update(op);
             if(precond_) {
                auto ls_ptr = dynamic_cast<LinearSolver<Matrix, Vector> *>(precond_.get());
                if(ls_ptr) {
                    ls_ptr->update(op);    
                }
             }
         }

        std::shared_ptr<Preconditioner> precond_;   /*!< Preconditioner to be used. */  



    protected: 
        Scalar current_radius_;                     /*!< Radius on current iterate - used to solve constrained QP wrt TR bound. */  
        
    };
}

#endif //TR_SUBPROBLEM
