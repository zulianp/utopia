
#include "utopia_ProjectedGaussSeidelNew.hpp"
#include "utopia_ProjectedGaussSeidelSweep.hpp"
#include "utopia_Algorithms.hpp"
#include "utopia_make_unique.hpp"

namespace utopia {

    template<class Matrix, class Vector>
    ProjectedGaussSeidel<Matrix, Vector, PETSC>::~ProjectedGaussSeidel() {}

    template<class Matrix, class Vector>
    ProjectedGaussSeidel<Matrix, Vector, PETSC>::ProjectedGaussSeidel()
    : use_line_search_(false), use_symmetric_sweep_(true), l1_(false), n_local_sweeps_(3), use_sweeper_(true)
    {}

    template<class Matrix, class Vector>
    ProjectedGaussSeidel<Matrix, Vector, PETSC>::ProjectedGaussSeidel(const ProjectedGaussSeidel &other)
    :   use_line_search_(other.use_line_search_),
        use_symmetric_sweep_(other.use_symmetric_sweep_),
        l1_(other.l1_),
        n_local_sweeps_(other.n_local_sweeps_),
        use_sweeper_(other.use_sweeper_)
    {}

    template<class Matrix, class Vector>
    ProjectedGaussSeidel<Matrix, Vector, PETSC> * ProjectedGaussSeidel<Matrix, Vector, PETSC>::clone() const
    {
        auto ptr = new ProjectedGaussSeidel(*this);
        ptr->set_box_constraints(this->get_box_constraints());
        return ptr;
    }

    template<class Matrix, class Vector>
    void ProjectedGaussSeidel<Matrix, Vector, PETSC>::read(Input &in)
    {
        QPSolver<Matrix, Vector>::read(in);

        in.get("use_line_search", use_line_search_);
        in.get("use_symmetric_sweep", use_symmetric_sweep_);
        in.get("n_local_sweeps", n_local_sweeps_);
        in.get("l1", l1_);
        in.get("use_sweeper", use_sweeper_);
    }

    template<class Matrix, class Vector>
    void ProjectedGaussSeidel<Matrix, Vector, PETSC>::print_usage(std::ostream &os) const
    {
        QPSolver<Matrix, Vector>::print_usage(os);

        this->print_param_usage(os, "use_line_search", "bool", "Determines if line-search should be used.", "true");
        this->print_param_usage(os, "use_symmetric_sweep", "bool", "Determines if symmetric local should be used.", "true");
        this->print_param_usage(os, "n_local_sweeps", "int", "Number of local sweeps.", "3");
    }

    template<class Matrix, class Vector>
    bool ProjectedGaussSeidel<Matrix, Vector, PETSC>::smooth(const Vector &b, Vector &x)
    {
        UTOPIA_TRACE_REGION_BEGIN("ProjectedGaussSeidel::smooth");

        const Matrix &A = *this->get_operator();

        // init(A);
        SizeType it = 0;
        SizeType n_sweeps = this->sweeps();
        if(this->has_bound()) {
            while(it++ < n_sweeps) {
                step(A, b, x);
                std::cout<<"--------------------------------------------- sweep --------------------------------- \n";
            }
        } else {
            while(unconstrained_step(A, b, x) && it++ < n_sweeps) {}
        }

        UTOPIA_TRACE_REGION_END("ProjectedGaussSeidel::smooth");
        return it == SizeType(this->sweeps() - 1);
    }

    template<class Matrix, class Vector>
    bool ProjectedGaussSeidel<Matrix, Vector, PETSC>::apply(const Vector &b, Vector &x)
    {
        UTOPIA_TRACE_REGION_BEGIN("ProjectedGaussSeidel::apply");

        if(this->verbose()) {
            if(l1_) {
                this->init_solver("utopia ProjectedL1GaussSeidelNew", {" it. ", "|| u - u_old ||"});
            } else {
                this->init_solver("utopia ProjectedGaussSeidel", {" it. ", "|| u - u_old ||"});
            }
        }


        const Matrix &A = *this->get_operator();

        x_old = x;
        bool converged = false;
        const SizeType check_s_norm_each = 1;

        int iteration = 0;
        while(!converged) {
            if(this->has_bound()) {
                step(A, b, x);
            } else {
                unconstrained_step(A, b, x);
            }

            if(iteration % check_s_norm_each == 0) {
                c = x - x_old;
                const Scalar diff = norm2(c);

                if(this->verbose()) {
                    PrintInfo::print_iter_status(iteration, {diff});
                }

                converged = this->check_convergence(iteration, 1, 1, diff);
            }

            ++iteration;

            if(converged) break;

            x_old = x;
        }

        UTOPIA_TRACE_REGION_END("ProjectedGaussSeidel::apply");
        return converged;
    }

    template<class Matrix, class Vector>
    void ProjectedGaussSeidel<Matrix, Vector, PETSC>::non_linear_jacobi_step(const Matrix &A, const Vector &b, Vector &x)
    {
        r = b - A * x;
        x = min(x + e_mul(d_inv, r), this->get_upper_bound());
    }

    ///residual must be computed outside
    template<class Matrix, class Vector>
    void ProjectedGaussSeidel<Matrix, Vector, PETSC>::apply_local_sweeps_unconstrained(const Matrix &A, const Vector &r, Vector &c) const
    {
        //reset correction
        c.set(0.0);

        if(use_sweeper_) {

            auto &&r_view     = const_local_view_device(r);
            auto &&c_view     = local_view_device(c);

            sweeper_->set_residual_view(r_view.array());
            sweeper_->set_correction_view(c_view.array());
            sweeper_->apply_unconstrained(this->n_local_sweeps());

        } else {

            //FIXME (use sytax as for the vectors)
            Matrix A_view;
            local_block_view(A, A_view);

            auto &&r_view     = const_local_view_device(r);
            auto &&d_inv_view = const_local_view_device(d_inv);

            auto &&c_view     = local_view_device(c);

            //WARNING THIS DOES NOT WORK IN PARALLEL (i.e., no OpenMP no Cuda)
            SizeType prev_i = 0;
            Scalar val = 0.0;

            /////////////////////////////////////////////
            auto sweeper = [&](const SizeType &i, const SizeType &j, const Scalar &a_ij) {
                if(i != prev_i) {
                    //next
                    const Scalar temp = d_inv_view.get(prev_i) * val;

                    c_view.set(prev_i, temp );
                    val = r_view.get(i);

                    prev_i = i;
                }

                val -= Scalar(i != j) * a_ij * c_view.get(j);
            };

            /////////////////////////////////////////////

            for(SizeType il = 0; il < this->n_local_sweeps(); ++il) {
                prev_i = 0;
                val = r_view.get(0);

                A_view.read(sweeper);

                //complete for last entry
                c_view.set(prev_i, d_inv_view.get(prev_i) * val );

                if(use_symmetric_sweep_) {
                    val = r_view.get(prev_i);

                    A_view.read_reverse(sweeper);

                    //complete for last entry
                    c_view.set(prev_i, d_inv_view.get(prev_i) * val );
                }
            }

        }
    }

    ///residual must be computed outside
    template<class Matrix, class Vector>
    void ProjectedGaussSeidel<Matrix, Vector, PETSC>::apply_local_sweeps(
        const Matrix &A,
        const Vector &r,
        const Vector &lb,
        const Vector &ub,
        Vector &c) const
    {
        // UTOPIA_TRACE_REGION_BEGIN("ProjectedGaussSeidel::apply_local_sweeps");

        //reset correction
        c.set(0.0);

        //FIXME (use sytax as for the vectors)
        Matrix A_view;
        local_block_view(A, A_view);


        if(use_sweeper_) {

            auto &&r_view     = const_local_view_device(r);
            auto &&lb_view    = const_local_view_device(lb);
            auto &&ub_view    = const_local_view_device(ub);
            auto &&c_view     = local_view_device(c);

            sweeper_->set_residual_view(r_view.array());
            sweeper_->set_bounds(lb_view.array(), ub_view.array());
            sweeper_->set_correction_view(c_view.array());
            sweeper_->apply(this->n_local_sweeps());

        } else {
            auto &&r_view     = const_local_view_device(r);
            auto &&lb_view    = const_local_view_device(lb);
            auto &&ub_view    = const_local_view_device(ub);
            auto &&d_inv_view = const_local_view_device(d_inv);

            auto &&c_view     = local_view_device(c);
            //WARNING THIS DOES NOT WORK IN PARALLEL
            SizeType prev_i = 0;
            Scalar val = 0.0;

            /////////////////////////////////////////////
            auto sweeper = [&](const SizeType &i, const SizeType &j, const Scalar &a_ij) {
                if(i != prev_i) {
                    //next
                    const Scalar temp = device::max(
                        lb_view.get(prev_i),
                        device::min(
                            d_inv_view.get(prev_i) * val,
                            ub_view.get(prev_i)
                            )
                        );


                    c_view.set(prev_i, temp );
                    val = r_view.get(i);

                    prev_i = i;
                }

                val -= Scalar(i != j) * a_ij * c_view.get(j);
            };

            /////////////////////////////////////////////

            for(SizeType il = 0; il < this->n_local_sweeps(); ++il) {
                prev_i = 0;
                val = r_view.get(0);

                A_view.read(sweeper);

                //complete for last entry
                c_view.set(prev_i, d_inv_view.get(prev_i) * val );

                if(use_symmetric_sweep_) {
                    val = r_view.get(prev_i);

                    A_view.read_reverse(sweeper);

                    //complete for last entry
                    c_view.set(prev_i, d_inv_view.get(prev_i) * val );
                }
            }
        }

        // UTOPIA_TRACE_REGION_END("ProjectedGaussSeidel::apply_local_sweeps");
    }

    template<class Matrix, class Vector>
    bool ProjectedGaussSeidel<Matrix, Vector, PETSC>::unconstrained_step(const Matrix &A, const Vector &b, Vector &x)
    {
        r = A * x;
        r = b - r;

        apply_local_sweeps_unconstrained(A, r, c);

        Scalar alpha = 1.;

        if(use_line_search_) {

            Ac = A * c;
            alpha = dot(c, r)/dot(Ac, c);

            if(std::isinf(alpha)) {
                return true;
            }

            if(std::isnan(alpha)) {
                return false;
            }

            if(alpha <= 0) {
                std::cerr << "[Warning] negative alpha" << std::endl;
                alpha = 1.;
                c = r;
            }
        }

        x += alpha * c;
        return true;
    }

    template<class Matrix, class Vector>
    bool ProjectedGaussSeidel<Matrix, Vector, PETSC>::step(const Matrix &A, const Vector &b, Vector &x)
    {
        r = A * x;
        r = b - r;

        //localize gap function for correction
        this->fill_empty_bounds(layout(x));
        ub_loc = this->get_upper_bound() - x;
        lb_loc = this->get_lower_bound() - x;

        apply_local_sweeps(A, r, lb_loc, ub_loc, c);

        if(use_line_search_)
        {
            Range rr = row_range(A);

            UTOPIA_NO_ALLOC_BEGIN("ProjectedGaussSeidel21");
            inactive_set_ *= 0.;
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("ProjectedGaussSeidel22");
            {
                Read<Vector> r_c(c), r_g(ub_loc), r_l(lb_loc);
                Write<Vector> w_a(inactive_set_);

                for(auto i = rr.begin(); i != rr.end(); ++i) {
                    if((c.get(i) < ub_loc.get(i)) && (c.get(i) > lb_loc.get(i))) {
                        inactive_set_.set(i, 1.);
                    }
                }
            }
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("ProjectedGaussSeidel23");
            is_c_ = e_mul(c, inactive_set_);

            Ac = A * is_c_;
            Scalar alpha = dot(is_c_, r)/dot(Ac, is_c_);
            UTOPIA_NO_ALLOC_END();

            if(std::isinf(alpha)) {
                return true;
            }

            if(std::isnan(alpha)) {
                return false;
            }

            assert(alpha > 0);

            if(alpha <= 0) {
                std::cerr << "[Warning] negative alpha" << std::endl;
                alpha = 1.;
                UTOPIA_NO_ALLOC_BEGIN("ProjectedGaussSeidel24");
                descent_dir = utopia::max(utopia::min(r, ub_loc), ub_loc);
                UTOPIA_NO_ALLOC_END();
            } else if(alpha <= 1.) {
                UTOPIA_NO_ALLOC_BEGIN("ProjectedGaussSeidel25");
                descent_dir = alpha * c;
                UTOPIA_NO_ALLOC_END();
            } else {
                UTOPIA_NO_ALLOC_BEGIN("ProjectedGaussSeidel26");
                c *= alpha;
                descent_dir = utopia::max(utopia::min(c, ub_loc), ub_loc);
                UTOPIA_NO_ALLOC_END();
            }
            UTOPIA_NO_ALLOC_BEGIN("ProjectedGaussSeidel27");
            x += descent_dir;
            UTOPIA_NO_ALLOC_END();
        }
        else
        {
            UTOPIA_NO_ALLOC_BEGIN("ProjectedGaussSeidel3");
            x += c;
            UTOPIA_NO_ALLOC_END();
        }

        return true;
    }

    template<class Matrix, class Vector>
    void ProjectedGaussSeidel<Matrix, Vector, PETSC>::init_memory(const Layout &layout)
    {
        Super::init_memory(layout);

        if(!use_sweeper_) {
            d.zeros(layout);
        }

        c.zeros(layout);

        if(use_line_search_) {
            inactive_set_.zeros(layout);
            Ac.zeros(layout);
            is_c_.zeros(layout);
            descent_dir.zeros(layout);
        }
    }

    template<class Matrix, class Vector>
    void ProjectedGaussSeidel<Matrix, Vector, PETSC>::init(const Matrix &A)
    {
        auto A_layout = row_layout(A);
        const bool reset = empty(c) || !A_layout.same(layout(c));

        if(reset) {
            init_memory(A_layout);
        } else {
            c.set(0);
            if(use_line_search_) {
                inactive_set_.set(0);
            }
        }

        if(use_sweeper_) {
            if(!sweeper_) {
                sweeper_ = utopia::make_unique< ProjectedGaussSeidelSweep<Scalar, SizeType> >();
            }

            sweeper_->symmetric(use_symmetric_sweep_);
            sweeper_->l1(l1_);

            Matrix A_local;
            local_block_view(A, A_local);

            if(reset) {
                sweeper_->init_from_local_matrix(A_local);
            } else {
                sweeper_->update_from_local_matrix(A_local);
            }

        } else {
            d = diag(A);

            if(l1_) {
                auto &&d_view = local_view_device(d);

                A.read(
                    UTOPIA_LAMBDA(const SizeType &i, const SizeType &, const Scalar &value) {
                         d_view.atomic_add(i, device::abs(value));
                    }
                );
            }

            d_inv = 1./d;
        }
    }

    template<class Matrix, class Vector>
    void ProjectedGaussSeidel<Matrix, Vector, PETSC>::update(const std::shared_ptr<const Matrix> &op)
    {
        UTOPIA_TRACE_REGION_BEGIN("ProjectedGaussSeidel::update");
        IterativeSolver<Matrix, Vector>::update(op);
        init(*op);
        UTOPIA_TRACE_REGION_END("ProjectedGaussSeidel::update");
    }

}
