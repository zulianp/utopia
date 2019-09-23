#ifndef UTOPIA_POISSON_HPP
#define UTOPIA_POISSON_HPP

#include "utopia_Core.hpp"
#include "utopia_ExtendedFunction.hpp"

#include <vector>

namespace utopia {

    template<class FunctionSpace, class Matrix, class Vector>
    class FormPoisson final : public ExtendedFunction<Matrix, Vector> {
    public:
        using Scalar   = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using IndexSet = typename utopia::Traits<Vector>::IndexSet;

        FormPoisson(FunctionSpace &V) : V_(V)
        {
            initialize();
        }

        bool value(const Vector &x, typename Vector::Scalar &energy) const override
        {
            IndexSet ghost_nodes;
            convert(V_.dof_map().get_send_list(), ghost_nodes);
            Vector x_ =  ghosted(V_.dof_map().n_local_dofs(), V_.dof_map().n_dofs(), ghost_nodes);
            x_ = x;
            synchronize(x_);

            auto u  = trial(V_);

            auto uk = interpolate(x_, u);
            auto r  = interpolate(rhs_, u);

            utopia::assemble(0.5 * inner(grad(uk), grad(uk)) * dX - inner(uk, coeff(1.)) * dX, energy);

            if(mpi_world_rank() == 0) {
                std::cout << "energy:        " << energy << std::endl;
            }

            return true;
        }

        bool gradient_no_rhs(const Vector &x, Vector &gradient) const override
        {
            IndexSet ghost_nodes;
            convert(V_.dof_map().get_send_list(), ghost_nodes);
            Vector x_ =  ghosted(V_.dof_map().n_local_dofs(), V_.dof_map().n_dofs(), ghost_nodes);
            
            x_ = x;
            synchronize(x_);

            auto u  = trial(V_);
            auto v  = test(V_);
            auto uk = interpolate(x_, u);
            auto r  = interpolate(rhs_, u);

            auto l_form = inner(grad(uk), grad(v)) * dX;

            utopia::assemble(l_form, gradient);

            gradient -= rhs_;
            apply_zero_boundary_conditions(V_.dof_map(), gradient);


            Scalar n_grad = norm2(gradient);
            Scalar n_x = norm2(x);

            // static int n_out = 0;
            // write("g" + std::to_string(n_out++) + ".m", gradient);

            if(mpi_world_rank() == 0) {
                std::cout << "norm(x):        " << n_x << std::endl;
                std::cout << "norm(gradient): " << n_grad << std::endl;
            }

            return true;
        }

        bool hessian(const Vector &x, Matrix &hessian) const override
        {
            hessian = H_;

            Scalar n_hess = norm2(hessian);

            if(mpi_world_rank() == 0) {
                std::cout << "norm(hessian): " << n_hess << std::endl;
            }

            return true;
        }

    private:
        FunctionSpace &V_;
        Vector rhs_;
        Matrix H_;

        void initialize()
        {
            auto u  = trial(V_);
            auto v  = test(V_);
            auto b_form = inner(grad(u), grad(v)) * dX;
            utopia::assemble(b_form, H_);
            utopia::assemble(inner(coeff(1.), v) * dX, rhs_);

            apply_boundary_conditions(V_.dof_map(), H_, rhs_);


            /// stuff for extended function
             Vector x = local_zeros(V_.dof_map().n_local_dofs());
             apply_boundary_conditions(V_.dof_map(), x);

             Vector marked;
             mark_constrained_dofs(V_.dof_map(), marked);
             this->set_equality_constrains(marked, x);


        }

    };


    template<class FunctionSpace, class Matrix, class Vector>
    class Poisson final : public ExtendedFunction<Matrix, Vector> {
    public:
        typedef typename utopia::Traits<Vector>::Scalar Scalar;
        typedef typename utopia::Traits<Vector>::SizeType SizeType;

        Poisson(FunctionSpace &V) : V_(V), rhs_value_(10000.)
        {
            initialize();
        }

        bool value(const Vector &x, typename Vector::Scalar &energy) const override
        {
            energy = 0.5 * dot(H_ * x, x) - dot(x, rhs_);
            return true;
        }

        bool gradient_no_rhs(const Vector &x, Vector &gradient) const override
        {
            gradient = H_ * x - rhs_;
            apply_zero_boundary_conditions(V_.dof_map(), gradient);
            return true;
        }

        bool hessian(const Vector &x, Matrix &hessian) const override
        {
            hessian = H_;
            return true;
        }

        bool update(const Vector &) override
        {
            return true;
        }

    private:
        FunctionSpace &V_;
        Vector rhs_;
        Matrix H_;
        Scalar rhs_value_;

        void initialize()
        {
            auto u  = trial(V_);
            auto v  = test(V_);

            utopia::assemble(inner(grad(u), grad(v)) * dX, H_);
            utopia::assemble(inner(coeff(rhs_value_), v) * dX, rhs_);


            apply_boundary_conditions(V_.dof_map(), H_, rhs_);


            ////////////////////////////////////////////////////////////
            /// stuff for extended function

            Vector x = local_zeros(V_.dof_map().n_local_dofs());
            apply_boundary_conditions(V_.dof_map(), x);

            Vector marked;
            mark_constrained_dofs(V_.dof_map(), marked);
            this->set_equality_constrains(marked, x);
        }

    };

}

#endif // UTOPIA_POISSON_HPP
