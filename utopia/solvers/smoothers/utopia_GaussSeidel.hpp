#ifndef UTOPIA_BLOCK_GAUSS_SEIDEL_HPP
#define UTOPIA_BLOCK_GAUSS_SEIDEL_HPP

#include "utopia_Smoother.hpp"
#include "utopia_Core.hpp"
#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_Allocations.hpp"


namespace utopia {

    /**
     * @brief      Homemade implementation of SOR.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
    class GaussSeidel final : public IterativeSolver<Matrix, Vector>
    {
        using Scalar   = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout   = typename Traits<Vector>::Layout;

        typedef utopia::IterativeSolver<Matrix, Vector> Solver;

    public:
        GaussSeidel(): use_line_search_(false), use_symmetric_sweep_(true), l1_(false), n_local_sweeps_(1), check_convergence_each_(10)
        {}

        void read(Input &in) override
        {
            Solver::read(in);
            in.get("l1", l1_);
        }

        void print_usage(std::ostream &os) const override
        {
            Solver::print_usage(os);
        }

        void check_convergence_each(const SizeType &n)
        {
            check_convergence_each_ = n;
        }


        /**
         * @brief      Smoothing of GS from Petsc. Currently we are using symmetric block GS (builds block jacobi and on blocks calls GS).
         *
         * @param[in]  A     The stiffness matrix.
         * @param[in]  rhs   The right hand side.
         * @param      x     The solution.
         */
        bool smooth(const Vector &rhs, Vector &x) override
        {
            const Matrix &A = *this->get_operator();
            SizeType it = 0;
            SizeType n_sweeps = this->sweeps();

            bool success = true;
            while(success && it++ < n_sweeps) {
                r = rhs - A * x;
                success = local_sweeps(A, r, c);
                x += c;
            }

            return it == SizeType(this->sweeps() - 1);
        }

        /**
         * @brief      Solving system with Gauss-Seidel method.
         *
         * @param[in]  rhs   The right hand side.
         * @param      x     The solution.
         */
        bool apply(const Vector &rhs, Vector &x) override
        {
            UTOPIA_NO_ALLOC_BEGIN("GaussSeidel::apply");

            if(this->verbose()) {
                if(l1_) {
                    this->init_solver("utopia L1GaussSeidel", {" it. ", "|| r ||"});
                } else {
                    this->init_solver("utopia GaussSeidel", {" it. ", "|| r ||"});
                }
            }

            const Matrix &A = *this->get_operator();
            bool converged = false;
            int iteration = 0;
            Scalar r_norm;

            r = rhs - A * x;
            r_norm = norm2(r);

            if(this->verbose()) {
                PrintInfo::print_iter_status(iteration, {r_norm});
            }

            while(!converged) {
                local_sweeps(A, r, c);
                x += c;
                r = rhs - A * x;

                ++iteration;

                if(iteration % check_convergence_each_ == 0) {
                    r_norm = norm2(r);

                    if(this->verbose()) {
                        PrintInfo::print_iter_status(iteration, {r_norm});
                    }

                    converged = this->check_convergence(iteration, 1, 1, r_norm);

                    if(converged) break;
                }
            }

            UTOPIA_NO_ALLOC_END();
            return converged;
        }


        inline GaussSeidel * clone() const override
        {
            return new GaussSeidel(*this);
        }

        void update(const std::shared_ptr<const Matrix> &op) override {
            Solver::update(op);
            init(*op);
        }

        void use_line_search(const bool val)
        {
            use_line_search_ = val;
        }

        inline SizeType n_local_sweeps() const
        {
            return n_local_sweeps_;
        }

        inline void n_local_sweeps(const SizeType n_local_sweeps)
        {
            n_local_sweeps_ = n_local_sweeps;
        }

        inline void use_symmetric_sweep(const bool use_symmetric_sweep)
        {
            use_symmetric_sweep_ = use_symmetric_sweep;
        }

        inline void l1(const bool val)
        {
            l1_ = val;
        }

  private:
        bool local_sweeps(const Matrix &A, const Vector &r, Vector &c)
        {
            c *= 0.;

            Range rr = row_range(A);
            {
                ReadAndWrite<Vector> rw_c(c);
                Read<Vector> r_r(r), r_d_inv(d_inv);
                Read<Matrix> r_A(A);

                for(SizeType il = 0; il < this->n_local_sweeps(); ++il) {
                    for(auto i = rr.begin(); i != rr.end(); ++i) {
                        RowView<const Matrix> row_view(A, i);
                        decltype(i) n_values = row_view.n_values();

                        auto s = r.get(i);

                        for(auto index = 0; index < n_values; ++index) {
                            const decltype(i) j = row_view.col(index);
                            const auto a_ij = row_view.get(index);

                            if(rr.inside(j) && i != j) {
                                s -= a_ij * c.get(j);
                            }
                        }

                        //update correction
                        c.set(i, d_inv.get(i) * s );
                    }

                    if(use_symmetric_sweep_) {
                        for(auto i = rr.end()-1; i >= rr.begin(); --i) {
                            RowView<const Matrix> row_view(A, i);
                            decltype(i) n_values = row_view.n_values();

                            auto s = r.get(i);

                            for(auto index = 0; index < n_values; ++index) {
                                const decltype(i) j = row_view.col(index);
                                const auto a_ij = row_view.get(index);

                                if(rr.inside(j) && i != j) {
                                    s -= a_ij * c.get(j);
                                }
                            }

                            //update correction
                            c.set(i, d_inv.get(i) * s);
                        }
                    }
                }
            }

            Scalar alpha = 1.;

            if(use_line_search_) {
                Ac = A*c;
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

            c *= alpha;
            return true;
        }

        void init(const Matrix &A)
        {
            d = diag(A);

            if(l1_) {
                Write<Vector> w(d);
                each_read(A, [this](const SizeType &i, const SizeType &/*j*/, const Scalar &value) {
                    d.add(i, std::abs(value));
                });
            }

            d_inv = 1./d;
            init_memory(row_layout(A));
        }

    public:
        void init_memory(const Layout &layout) override
        {
            c.zeros(layout);
            r.zeros(layout);

            if(use_line_search_) {
                Ac.zeros(layout);
            }
        }



    private:
        bool use_line_search_;
        bool use_symmetric_sweep_;
        bool l1_;
        SizeType n_local_sweeps_;
        SizeType check_convergence_each_;

        Vector r, d, c, d_inv, Ac;

    };

}

#endif //UTOPIA_BLOCK_GAUSS_SEIDEL_HPP

