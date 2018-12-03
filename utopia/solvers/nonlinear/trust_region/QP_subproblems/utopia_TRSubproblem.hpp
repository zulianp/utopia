#ifndef TR_SUBPROBLEM_L2NORM_HPP
#define TR_SUBPROBLEM_L2NORM_HPP
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

            virtual void set_linear_solver(const std::shared_ptr<LinearSolver<Matrix, Vector> > &ls)
            {
                if(this->verbose())
                    std::cout<<"current TR strategy does not need linear solver \n";
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


            virtual TRSubproblem * clone() const override = 0;

        protected:
            /**
             * @brief      Solver solves easy problem of finding minimum,
             *             such that \f$ result = s + \tau p_k \f$ minimizes \f$ m_k(result) \f$
             *             and satisfies \f$ ||result|| = \Delta_k  \f$
             *
             * @param[in]  s
             * @param[in]  p_k
             * @param[in]  delta  The TR radius.
             * @param      result    The new iterate.
             *
             * @return     tau
             */
            Scalar quad_solver(const Vector &s, const Vector &p_k, const Scalar & delta,  Vector &result)
            {
                Scalar a, b, c, x1, x2, nom, denom,tau;

                a = dot(p_k, p_k);
                b = dot(2 * s, p_k);
                c = dot(s, s) - delta * delta;

                nom = b * b - 4 * a * c;
                nom =  std::sqrt(nom);
                denom = 2 * a;

                x1 = (- b + nom)/denom;
                x2 = (- b - nom)/denom;

                tau = std::max(x1, x2);

                if(tau != tau)
                    tau = 0;

                result = s + tau * p_k;
                return tau;
            }



        Scalar quadratic_function(const Scalar & a,  const Scalar & b, const Scalar &c)
        {
            Scalar sqrt_discriminant = std::sqrt( b * b - 4.0 * a * c);

            Scalar lower = (-b + sqrt_discriminant)/ (2.0 * a);
            Scalar upper = (-b - sqrt_discriminant)/ (2.0 * a);

            return std::max(upper, lower);
        }


    public:
        virtual bool unpreconditioned_solve(const Matrix &/*B*/, const Vector &/*g*/, Vector &/*p_k*/){ return false; };
        virtual bool preconditioned_solve(const Matrix &/*B*/, const Vector &/*g*/, Vector &/*p_k*/){ return false; };

        virtual bool unpreconditioned_solve(const Operator<Vector> &/*B*/, const Vector &/*g*/, Vector &/*p_k*/){ return false; };
        virtual bool preconditioned_solve(const Operator<Vector> &/*B*/, const Vector &/*g*/, Vector &/*p_k*/){ return false; };


        virtual bool tr_constrained_solve(const Operator<Vector> &H, const Vector &g, Vector &s, const Scalar & tr_radius)
        {
            this->current_radius(tr_radius);

            if(this->precond_) 
            {
                return preconditioned_solve(H, g, s);
            } else {
                return unpreconditioned_solve(H, g, s);
            }
            return true;
        };

        virtual bool tr_constrained_solve(const Matrix &H, const Vector &g, Vector &p_k, const Scalar & tr_radius)
        {
            current_radius(tr_radius);
            update(make_ref(H));
            apply(g, p_k);
            return true;
        }

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



    protected:
        std::shared_ptr<Preconditioner> precond_;   /*!< Preconditioner to be used. */
        Scalar current_radius_;                     /*!< Radius on current iterate - used to solve constrained QP wrt TR bound. */



    };
}

#endif //TR_SUBPROBLEM_L2NORM_HPP
