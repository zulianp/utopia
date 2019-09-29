#ifndef UTOPIA_BLOCK_GAUSS_SEIDEL_HPP
#define UTOPIA_BLOCK_GAUSS_SEIDEL_HPP

#include "utopia_Smoother.hpp"
#include "utopia_Core.hpp"
#include "utopia_LinearSolverInterfaces.hpp"


namespace utopia {

    /**
     * @brief      Wrapper for PETSC implementation of SOR.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
    class GaussSeidel: public IterativeSolver<Matrix, Vector>, public Smoother<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                   Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                SizeType;
        typedef utopia::IterativeSolver<Matrix, Vector> Solver;
        typedef utopia::Smoother<Matrix, Vector>        Smoother;

    public:
        GaussSeidel(): use_line_search_(true), use_symmetric_sweep_(true), n_local_sweeps_(1)
        {

        }

        void read(Input &in) override
        {
            Solver::read(in);
            Smoother::read(in);
        }

        void print_usage(std::ostream &os) const override
        {
            Solver::print_usage(os);
            Smoother::print_usage(os);
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

            // init(A);
            SizeType it = 0;
            SizeType n_sweeps = this->sweeps();

            while(unconstrained_step(A, rhs, x) && it++ < n_sweeps) {}
            
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
            if(this->verbose())
                this->init_solver("utopia GaussSeidel", {" it. ", "|| r ||"});

            const Matrix &A = *this->get_operator();
            bool converged = false;

            int iteration = 0;

            Scalar r_norm = norm2(A * x - rhs); 
            if(this->verbose()) {
                PrintInfo::print_iter_status(iteration, {r_norm});
            }            

            while(!converged) {
                unconstrained_step(A, rhs, x);
                r_norm = norm2(A * x - rhs); 
                
                ++iteration;
                if(this->verbose()) {
                    PrintInfo::print_iter_status(iteration, {r_norm});
                }

                converged = this->check_convergence(iteration, 1, 1, r_norm);

                if(converged) break;
            }

            return converged;
        }


        inline GaussSeidel * clone() const override
        {
            return new GaussSeidel(*this);
        }

        virtual void update(const std::shared_ptr<const Matrix> &op) override
        {
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



  protected:
        bool unconstrained_step(const Matrix &A, const Vector &b, Vector &x)
        {
            r = b - A * x;
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

                alpha = dot(c, r)/dot(A * c, c);

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


        void init(const Matrix &A)
        {
            d = diag(A);
            d_inv = 1./d;
            c = local_zeros(local_size(A).get(0));
        }        





      private:
        bool use_line_search_;
        bool use_symmetric_sweep_;
        SizeType n_local_sweeps_;

        Vector r, d, c, d_inv, x_old, descent_dir;

    };

}

#endif //UTOPIA_BLOCK_GAUSS_SEIDEL_HPP

