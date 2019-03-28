#ifndef TR_SUBPROBLEM_L2NORM_HPP
#define TR_SUBPROBLEM_L2NORM_HPP
#include <string>
#include "utopia_IterativeSolver.hpp"

namespace  utopia
{


    template<class Vector>
    class TRSubproblemBase : public virtual Configurable
    {
        public:
            typedef UTOPIA_SCALAR(Vector) Scalar;

            virtual ~TRSubproblemBase( ){}

            TRSubproblemBase(): current_radius_(1e14)
            {

            }

            void current_radius(const Scalar &radius)
            {
                current_radius_ = radius;
            }

            Scalar current_radius()
            {
                return current_radius_;
            }

            virtual void read(Input &in) override
            {
                in.get("current_radius", current_radius_);
            }

            virtual void print_usage(std::ostream &os) const override
            {
                this->print_param_usage(os, "current_radius", "real", "Value of trust region radius.", "1e14");
            }

        protected:
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


        protected:
            Scalar current_radius_;
    };


    template<class Matrix, class Vector>
    class TRSubproblem : public IterativeSolver<Matrix, Vector>, public virtual TRSubproblemBase<Vector>
    {
        typedef UTOPIA_SCALAR(Vector) Scalar;

        public:
            TRSubproblem()
            {

            }

            virtual ~TRSubproblem( ){}
            virtual TRSubproblem * clone() const override = 0;

            virtual void read(Input &in) override
            {
                IterativeSolver<Matrix, Vector>::read(in);
                TRSubproblemBase<Vector>::read(in);
            }

            virtual void print_usage(std::ostream &os) const override
            {
                IterativeSolver<Matrix, Vector>::print_usage(os);
                TRSubproblemBase<Vector>::print_usage(os);
            }
    };


    template<class Vector>
    class MatrixFreeTRSubproblem : public MatrixFreeLinearSolver<Vector>, public virtual TRSubproblemBase<Vector>
    {
        typedef UTOPIA_SCALAR(Vector) Scalar;

        public:
            MatrixFreeTRSubproblem()
            {

            }

            virtual ~MatrixFreeTRSubproblem( ){}
            virtual MatrixFreeTRSubproblem * clone() const override= 0;

            virtual void read(Input &in) override
            {
                MatrixFreeLinearSolver<Vector>::read(in);
                TRSubproblemBase<Vector>::read(in);
            }

            virtual void print_usage(std::ostream &os) const override
            {
                MatrixFreeLinearSolver<Vector>::print_usage(os);
                TRSubproblemBase<Vector>::print_usage(os);
            }
    };

}

#endif //TR_SUBPROBLEM_L2NORM_HPP
