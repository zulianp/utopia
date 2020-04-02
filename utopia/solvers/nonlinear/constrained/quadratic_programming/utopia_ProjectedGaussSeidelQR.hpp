#ifndef UTOPIA_PROJECTED_GAUSS_SEIDEL_QR_HPP
#define UTOPIA_PROJECTED_GAUSS_SEIDEL_QR_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_QPSolver.hpp"
#include "utopia_RowView.hpp"
#include "utopia_ProjectedGaussSeidel.hpp"
#include <cmath>

namespace utopia {
    template<class Matrix, class Vector>
    class ProjectedGaussSeidelQR final: public ProjectedGaussSeidel<Matrix, Vector>
    {
    public:

        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        ProjectedGaussSeidelQR()
        {
            ProjectedGaussSeidel<Matrix, Vector>::use_line_search(false);

            if(mpi_world_size() > 1 ){
                utopia_error("ProjectedGaussSeidelQR does not support parallel computations. \n");
            }
        }

        ProjectedGaussSeidelQR(const ProjectedGaussSeidelQR &) = default;

        inline ProjectedGaussSeidelQR * clone() const override
        {
            auto ptr = new ProjectedGaussSeidelQR(*this);
            ptr->set_box_constraints(this->get_box_constraints());
            // ptr->set_R(this->get_R());
            return ptr;
        }

        void read(Input &in) override
        {
            ProjectedGaussSeidel<Matrix, Vector>::read(in);
        }

        void set_R(const Matrix &R)
        {
            R_ = R;
        }

        const Matrix & get_R()
        {
            return R_;
        }

        bool smooth(const Vector &b, Vector &x) override
        {
            const Matrix &A = *this->get_operator();

            SizeType it = 0;
            SizeType n_sweeps = this->sweeps();
            if(this->has_bound()) {
                while(it++ < n_sweeps){
                    this->step(A, b, x);
                }
            } else {
                while(it++ < n_sweeps) {
                    ProjectedGaussSeidel<Matrix, Vector>::unconstrained_step(A, b, x);
                }
            }
            return it == SizeType(this->sweeps() - 1);
        }

        bool apply(const Vector &b, Vector &x) override
        {
            if(this->verbose())
                this->init_solver("utopia ProjectedGaussSeidel", {" it. ", "|| u - u_old ||"});

            const Matrix &A = *this->get_operator();

            x_old = x;
            bool converged = false;
            const SizeType check_s_norm_each = 1;

            int iteration = 0;
            while(!converged) {
                if(this->has_bound()) {
                    step(A, b, x);
                } else {
                    this->unconstrained_step(A, b, x);
                }

                if(iteration % check_s_norm_each == 0) {
                    const Scalar diff = norm2(x_old - x);

                    if(this->verbose()) {
                        PrintInfo::print_iter_status({static_cast<Scalar>(iteration), diff});
                    }

                    converged = this->check_convergence(iteration, 1, 1, diff);
                }

                ++iteration;

                if(converged)
                    break;

                x_old = x;
            }
            return converged;
        }


        void print_usage(std::ostream &os) const override
        {
            ProjectedGaussSeidel<Matrix, Vector>::print_usage(os);

            this->print_param_usage(os, "use_line_search", "bool", "Determines if line-search should be used.", "true");
            this->print_param_usage(os, "n_local_sweeps", "int", "Number of local sweeps.", "3");
        }

        bool step(const Matrix &A, const Vector &b, Vector &x) override
        {
            active_set_.set(0.0);

            const Vector & g = this->get_upper_bound();
            const Vector & l = this->get_lower_bound();

            Scalar g_i, l_i;

            Range rr = row_range(A);
            {
                Write<Vector> w_a(active_set_);
                ReadAndWrite<Vector> rw_x(x);
                Read<Vector> r_d_inv(d_inv), r_g(g), r_l(l);
                Read<Matrix> r_A(A);
                Read<Vector> r_b(b);
                Read<Matrix> r_R(R_);
                SizeType n_rows = local_size(R_).get(0);


                for(auto i = rr.begin(); i != rr.end(); ++i)
                {
                    RowView<const Matrix> row_view(A, i);
                    decltype(i) n_values = row_view.n_values();

                    Scalar s = 0;//x.get(i);

                    for(auto index = 0; index < n_values; ++index)
                    {
                        const decltype(i) j = row_view.col(index);
                        const auto a_ij = row_view.get(index);

                        if(rr.inside(j))// && j != i)
                        {
                            s += a_ij * x.get(j);
                        }
                    }

                    //update correction
                    x.set(i, d_inv.get(i)*(b.get(i) - s) + x.get(i)) ;

                    if  (i < n_rows)
                    {
                        RowView<const Matrix> row_viewR(R_, i);
                        decltype(i) nnz_R = row_viewR.n_values();

                        Scalar r_sum_c = 0.0;
                        Scalar r_ii = 0.0;

                        for(auto index = 0; index < nnz_R; ++index)
                        {
                            // r_ii = 0.0;
                            const decltype(i) j = row_viewR.col(index);
                            const auto r_ij = row_viewR.get(index);

                            if(j < i)
                            {
                                r_sum_c += r_ij * x.get(j);
                            }
                            else if(i == j)
                            {
                                r_ii = r_ij;
                            }
                        }

                        if (r_ii > 0)
                        {
                            g_i = (g.get(i) - r_sum_c)/r_ii;
                            l_i = (l.get(i) - r_sum_c)/r_ii;
                        }
                        else if (r_ii < 0)
                        {
                            l_i = (g.get(i) - r_sum_c)/r_ii;
                            g_i = (l.get(i) - r_sum_c)/r_ii;
                        }
                        else
                        {
                            std::cerr<<"--------- ProjectedGaussSeidelQR::PGSQR......... \n";
                        }

                        //update correction
                        if (( g_i <= x.get(i) ) || (x.get(i) <= l_i))
                        {
                            x.set(i, std::max(std::min( x.get(i), g_i), l_i));
                            active_set_.set(i, 1.0);
                        }
                        else
                        {
                            active_set_.set(i, 0.0);
                        }
                    }
                }
            }

            return true;
        }

        void init(const Matrix &A)
        {
            d_inv = diag(A);
            d_inv = 1./d_inv;
            active_set_.zeros(layout(d_inv));
        }


        virtual void update(const std::shared_ptr<const Matrix> &op) override
        {
            ProjectedGaussSeidel<Matrix, Vector>::update(op);
            init(*op);
        }

        const Vector& get_active_set()
        {
            return active_set_;
        }


    private:
        Vector d_inv, x_old;
        Vector active_set_;
        Matrix R_;
    };
}

#endif //UTOPIA_PROJECTED_GAUSS_SEIDEL_QR_HPP
