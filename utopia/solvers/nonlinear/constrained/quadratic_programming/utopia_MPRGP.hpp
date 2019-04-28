#ifndef UTOPIA_CONSTRAINT_QP_MPRGP
#define UTOPIA_CONSTRAINT_QP_MPRGP

#include <string>
#include "utopia_BoxConstraints.hpp"
#include "utopia_QPSolver.hpp"

namespace  utopia
{

    template<class Matrix, class Vector>
    class MPGRP final:  public MatrixFreeQPSolver<Vector>, public QPSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                   Scalar;
        typedef utopia::LinearSolver<Matrix, Vector>    Solver;

        public:
            MPGRP()
            {

            }

            void read(Input &in) override
            {
                MatrixFreeQPSolver<Vector>::read(in);
                QPSolver<Matrix, Vector>::read(in);

                // in.get("Cauchy-point", cp_);

                // if(precond_) {
                //     in.get("precond", *precond_);
                // }
                // if(linear_solver_) {
                //     in.get("linear-solver", *linear_solver_);
                // }
            }


            void print_usage(std::ostream &os) const override
            {
                MatrixFreeQPSolver<Vector>::print_usage(os);
                QPSolver<Matrix, Vector>::print_usage(os);

                // this->print_param_usage(os, "Cauchy-point", "GeneralizedCauchyPoint", "Input parameters for Generalized Cauchy point solver.", "-");
                // this->print_param_usage(os, "precond", "Preconditioner", "Input parameters for Preconditioner.", "-");
                // this->print_param_usage(os, "linear-solver", "LinearSolver", "Input parameters for LinearSolver.", "-");
            }



            MPGRP * clone() const override
            {
                return new MPGRP(*this);
            }

            void set_preconditioner(const std::shared_ptr<Preconditioner<Vector> > &precond)
            {
                precond_ = precond;
            }

            bool solve(const Operator<Vector> &A, const Vector &rhs, Vector &sol) override
            {
                auto &box = this->get_box_constraints();
                // init(local_size(rhs).get(0));
                return aux_solve(A, rhs, sol, box);
            }

            bool solve(const Matrix &A, const Vector &rhs, Vector &sol) override
            {
                auto A_op_ptr = utopia::op_ref(A);
                auto &box = this->get_box_constraints();
                // init(local_size(rhs).get(0));
                return aux_solve(*A_op_ptr, rhs, sol, box);
            }


        private:
            bool aux_solve(const Operator<Vector> &H, const Vector &g, Vector &s, const BoxConstraints<Vector> & constraints)
            {
                const auto &ub = constraints.upper_bound();
                const auto &lb = constraints.lower_bound();

                if(this->verbose()){
                    this->init_solver("ProjectedGradientActiveSet", {" "});
                }

    

                return true;
            }


          

        private:

            void get_fi(const Vector &x, const Vector &g, const Vector &lb, const Vector &ub, Vector & fi) const
            {
                fi = local_values(local_size(x).get(0), 0); 

                {
                    Read<Vector> r_ub(ub), r_lb(lb), r_x(x), r_g(g);
                    
                    each_write(fi, [&ub, &lb, &x, &g](const SizeType i) -> double 
                    {
                        Scalar li =  lb.get(i); Scalar ui =  ub.get(i); Scalar xi =  x.get(i);
                        if(li < xi && xi < ui)
                            return g.get(i);
                        else
                            return 0.0; 
                    }   );
                }
            }


            Scalar get_alpha_f(const Vector &x, const Vector &p, const Vector &lb, const Vector &ub) const
            {
                Vector alpha_f1 = local_values(local_size(x).get(0), 1e15); 
                Vector alpha_f2 = local_values(local_size(x).get(0), 1e15); 

                {
                    Read<Vector> r_ub(ub), r_lb(lb), r_x(x), r_g(p);
                    
                    each_write(alpha_f1, [&lb, &x, &p](const SizeType i) -> double 
                    {
                        Scalar li = lb.get(i); 
                        Scalar xi = x.get(i); 
                        Scalar pi = p.get(i);

                        if(pi > 0)
                        {
                            return (xi-li)/pi; 
                        }
                        else
                        {
                            return 1e15; 
                        }

                    }   );



                    each_write(alpha_f2, [&ub, &x, &p](const SizeType i) -> double 
                    {
                        Scalar ui = ub.get(i); 
                        Scalar xi = x.get(i); 
                        Scalar pi = p.get(i);

                        if(pi < 0)
                        {
                            return (xi-ui)/pi; 
                        }
                        else
                        {
                            return 1e15; 
                        }

                    }   );

                }

                Scalar alpha_f1_min = min(alpha_f1); 
                Scalar alpha_f2_min = min(alpha_f2); 

                return std::min(alpha_f1_min, alpha_f2_min); 
            }



        private:
            std::shared_ptr<Preconditioner<Vector> > precond_;

            Vector r, q, d, Hd;

    };
}

#endif //UTOPIA_CONSTRAINT_QP_MPRGP
