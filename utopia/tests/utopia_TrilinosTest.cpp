#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS
#include "utopia.hpp"
#include "utopia_Jacobi.hpp"
#include "utopia_Testing.hpp"
#include "utopia_trilinos.hpp"
#include "utopia_trilinos_Each_impl.hpp"
#include "utopia_trilinos_Utils.hpp"
#include "utopia_trilinos_solvers.hpp"

#include <algorithm>
#include "utopia_MultilevelTestProblem1D.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

#ifdef WITH_PETSC
#include "utopia_petsc_trilinos.hpp"
#endif

#include <cmath>
#include "utopia_Bratu1D.hpp"
#include "utopia_Eval_Structure.hpp"
#include "utopia_IPTransfer.hpp"
#include "utopia_Poisson1D.hpp"
#include "utopia_Structure.hpp"
#include "utopia_TestProblems.hpp"

namespace utopia {

    class TrilinosTest {
    public:
        using Traits = utopia::Traits<TpetraVectord>;
        using Comm = Traits::Communicator;
        using SizeType = Traits::SizeType;
        using LocalSizeType = Traits::LocalSizeType;
        using Scalar = Traits::Scalar;

        Comm comm_;

        template <class Matrix>
        void build_rectangular_matrix(const SizeType &n, const SizeType &m, Matrix &mat) const {
            using TraitsT = utopia::Traits<Matrix>;
            // mat  = local_sparse(n, m, 2);
            mat.sparse(layout(comm_, n, m, TraitsT::determine(), TraitsT::determine()), 2, 2);

            Write<Matrix> w_(mat);
            auto r = row_range(mat);
            auto cols = size(mat).get(1);
            for (auto i = r.begin(); i < r.end(); ++i) {
                if (i >= cols) {
                    break;
                }

                mat.set(i, i, 1.);
            }
        }

        template <class Matrix>
        void build_rectangular_matrix_2(const SizeType &n, const SizeType &m, Matrix &mat) const {
            using TraitsT = utopia::Traits<Matrix>;
            // mat  = local_sparse(n, m, 2);
            mat.sparse(layout(comm_, n, m, TraitsT::determine(), TraitsT::determine()), 2, 2);

            Write<Matrix> w_(mat);
            auto r = row_range(mat);
            auto cols = size(mat).get(1);
            for (auto i = r.begin(); i < r.end(); ++i) {
                if (i >= cols) {
                    break;
                }

                mat.set(i, i, 1.);
            }

            mat.set(0, m - 1, 1.);
        }

        void trilinos_build() {
            auto n = 10;
            TpetraMatrixd m;
            m.sparse(layout(comm_, n, n, Traits::determine(), Traits::determine()), 3, 2);
            assemble_laplacian_1D(m);

            TpetraVectord v(row_layout(m), 1.);

            TpetraVectord z = m * v;
            double nz = norm2(z);
            utopia_test_assert(approxeq(nz, 0.));
        }

        void trilinos_build_identity() {
            auto n = 10;
            TpetraMatrixd id;
            id.identity(layout(comm_, n, n, Traits::determine(), Traits::determine()), 1.0);

            TpetraVectord v(row_layout(id), 2.);
            double actual = norm1(id * v);

            utopia_test_assert(approxeq(size(v).get(0) * 2., actual));

            TpetraMatrixd id_t = transpose(id);
            actual = norm1(id * v);

            utopia_test_assert(approxeq(size(v).get(0) * 2., actual));
        }

        void trilinos_rect_matrix() {
            TpetraMatrixd P;
            // build_rectangular_matrix(5, 10, P);
            build_rectangular_matrix(10, 5, P);

            // auto rm = P.implementation().implementation().getRangeMap();
            // auto dm = P.implementation().implementation().getDomainMap();

            // auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));

            // rm->describe(*out);
            // dm->describe(*out);

            utopia_test_assert(P.is_valid(true));
        }

        void trilinos_accessors() {
            auto n = 10;
            TpetraVectord v(layout(comm_, n, Traits::determine()), 10.);

            {
                Range r = range(v);
                Write<TpetraVectord> w_v(v);

                // set first and last entries of each process, are to be 0
                v.set(r.begin(), 0.);
                v.add(r.end() - 1, -10.);
            }

            Size s = size(v);
            Size ls = local_size(v);

            utopia_test_assert(s.get(0) == n * mpi_world_size());
            utopia_test_assert(ls.get(0) == n);

            // FIXME replace this with an actual test
            // disp(v);
        }

        void trilinos_matrix_access() {
            auto n = 10;
            TpetraMatrixd Y;
            Y.sparse(layout(comm_, n, n, Traits::determine(), Traits::determine()), 3, 2);

            assemble_laplacian_1D(Y);

            auto rr = row_range(Y);
            auto i = rr.begin();

            {
                Read<TpetraMatrixd> r_(Y);
                if (i == 0 || i + 1 == size(Y).get(0)) {
                    utopia_test_assert(approxeq(Y.get(i, i), 1.));
                } else {
                    utopia_test_assert(approxeq(Y.get(i, i), 2.));
                }
            }

            {
                Write<TpetraMatrixd> w_(Y);
                Y.set(i, i, 4.);
            }

            {
                Read<TpetraMatrixd> r_(Y);
                utopia_test_assert(approxeq(Y.get(i, i), 4.));
            }
        }

        void trilinos_set() {
            auto n = 10;
            TpetraVectord v(layout(comm_, n, Traits::determine()), 0.);
            double sum_v = sum(v);

            utopia_test_assert(approxeq(sum_v, 0.));
        }

        void trilinos_vec_minus() {
            auto n = 10;
            TpetraVectord y(layout(comm_, n, Traits::determine()), 1.);
            TpetraVectord x(layout(y), 5.);

            TpetraVectord z;
            z = y - x;

            TpetraVectord expected(layout(y), -4.);

            double sum_z = double(sum(z));

            utopia_test_assert(approxeq(sum_z, size(z).get(0) * (-4.)));
            utopia_test_assert(approxeq(z, expected));
        }

        void trilinos_vec_axpy() {
            auto n = 10;
            TpetraVectord y(layout(comm_, n, Traits::determine()), 1.);
            TpetraVectord x(layout(y), 5.);
            auto alpha = 0.1;
            y += alpha * x;

            double val = norm1(y);
            utopia_test_assert(approxeq(val, n * mpi_world_size() * 1.5));
        }

        void trilinos_residual() {
            auto n = 10;
            TpetraVectord y(layout(comm_, n, Traits::determine()), 2.);
            TpetraVectord x(layout(y), 1.);
            TpetraMatrixd Id;
            Id.identity(layout(comm_, n, n, Traits::determine(), Traits::determine()), 1.0);

            TpetraVectord z = x - Id * y;

            double val = norm1(z);
            utopia_test_assert(approxeq(val, size(y).get(0)));
        }

        void trilinos_vec_scale() {
            auto n = 10;
            TpetraVectord y(layout(comm_, n, Traits::determine()), 1.);

            y *= 2.;

            double val = norm1(y);
            utopia_test_assert(approxeq(val, size(y).get(0) * 2.));

            TpetraVectord y2 = y * 2.;

            val = norm1(y2);
            utopia_test_assert(approxeq(val, size(y2).get(0) * 4.));
        }

        void trilinos_mat_scale() {
            auto n = 10;
            TpetraMatrixd Y;
            Y.sparse(layout(comm_, n, n, Traits::determine(), Traits::determine()), 3, 2);
            assemble_laplacian_1D(Y);

            Y *= 2.;

            TpetraVectord x(layout(comm_, n, Traits::determine()), 1.);
            double val = norm1(Y * x);
            utopia_test_assert(approxeq(val, 0.));

            TpetraMatrixd Y2 = Y * 2.;

            val = norm1(Y2 * x);
            utopia_test_assert(approxeq(val, 0.));
        }

        void trilinos_mat_axpy() {
            auto n = 10;
            TpetraMatrixd X;
            X.sparse(layout(comm_, n, n, Traits::determine(), Traits::determine()), 3, 2);
            assemble_laplacian_1D(X);

            TpetraMatrixd Y = X;

            auto alpha = 0.1;
            Y += alpha * X;

            TpetraVectord v(layout(comm_, n, Traits::determine()), 5.);

            double val = norm1(Y * v);
            double tolerance = 30. * std::numeric_limits<double>::epsilon();
            // std::cout << "val " << val <<std::endl;
            utopia_test_assert(approxeq(val, 0., tolerance));

            TpetraMatrixd Id;
            Id.identity(layout(comm_, n, n, Traits::determine(), Traits::determine()), 1.0);
            Id += 2. * Id;

            v.set(1.);
            val = norm1(Id * v);
            utopia_test_assert(approxeq(val, size(v).get(0) * 3., 1e-14));
        }

        void trilinos_mv() {
            auto n = 10;
            TpetraVectord x(layout(comm_, n, Traits::determine()), 5.);
            TpetraMatrixd m;
            m.sparse(layout(comm_, n, n, Traits::determine(), Traits::determine()), 3, 2);
            assemble_laplacian_1D(m);

            TpetraVectord y = m * x;

            const double val = norm2(y);
            utopia_test_assert(approxeq(val, 0.));
        }

        void trilinos_apply_transpose() {
            auto rows = 5;
            auto cols = 6;
            TpetraMatrixd A;
            A.sparse(layout(comm_, rows, cols, Traits::determine(), Traits::determine()), 2, 2);

            {
                Write<TpetraMatrixd> w_A(A);
                Range r = row_range(A);

                for (auto i = r.begin(); i < r.end(); ++i) {
                    A.set(i, i, 1.);
                    A.set(i, i + 1, 1.);
                }
            }

            TpetraVectord v(layout(comm_, rows, Traits::determine()), 1.);
            TpetraVectord At_v = transpose(A) * v;

            each_read(At_v, [](const SizeType i, const double val) { utopia_test_assert(val <= 2. + 1e-16); });

            double s_At_v = sum(At_v);
            utopia_test_assert(approxeq(s_At_v, size(A).get(0) * 2.));
        }

        void trilinos_apply_transpose_explicit() {
            auto rows = 5;
            auto cols = 6;
            TpetraMatrixd A;
            A.sparse(layout(comm_, rows, cols, Traits::determine(), Traits::determine()), 2, 2);

            {
                Write<TpetraMatrixd> w_A(A);
                Range r = row_range(A);

                for (auto i = r.begin(); i < r.end(); ++i) {
                    A.set(i, i, 1.);
                    A.set(i, i + 1, 1.);
                }
            }

            TpetraVectord v(layout(comm_, rows, Traits::determine()), 1.);
            // Explicit transpose
            TpetraMatrixd At = transpose(A);
            TpetraVectord At_v = At * v;

            // disp(At);

            each_read(At_v, [](const SizeType i, const double val) { utopia_test_assert(val <= 2. + 1e-16); });

            double s_At_v = sum(At_v);

            // disp(s_At_v);

            utopia_test_assert(approxeq(s_At_v, size(A).get(0) * 2.));
        }

        void trilinos_transpose() {
            auto rows = 5;
            auto cols = 5;
            TpetraMatrixd A;
            A.sparse(layout(comm_, rows, cols, Traits::determine(), Traits::determine()), 2, 2);
            auto gs = size(A);

            {
                Write<TpetraMatrixd> w_A(A);
                Range r = row_range(A);

                for (auto i = r.begin(); i < r.end(); ++i) {
                    A.set(i, 0, 1.);
                    A.set(i, gs.get(1) - 1, 2.);
                }
            }

            auto s = size(A);

            TpetraMatrixd At = transpose(A);
            auto s_t = size(At);
            utopia_test_assert(s_t.get(0) == s.get(1));

            // disp(A);
            // std::cout << "-----------------------" << std::endl;
            // disp(At);
            // std::cout << "-----------------------" << std::endl;

            TpetraMatrixd id;
            id.identity(layout(comm_, rows, cols, Traits::determine(), Traits::determine()), 1.0);
            TpetraVectord v(layout(comm_, cols, Traits::determine()), 2.);

            TpetraVectord actual = id * v;
            double norm_actual = norm1(actual);

            utopia_test_assert(approxeq(size(actual).get(0) * 2, norm_actual));

            // Does not work in parallel
            TpetraMatrixd id_t = transpose(id);
            TpetraVectord v2(layout(comm_, rows, Traits::determine()), 2.);
            actual = id_t * v2;
            norm_actual = norm1(actual);
            double norm_expected = size(v2).get(0) * 2.;

            // disp(id);
            // std::cout << "-----------------------" << std::endl;
            // disp(id_t);
            // std::cout << "-----------------------" << std::endl;

            // std::cout << norm_expected << " == " << norm_actual << std::endl;
            utopia_test_assert(approxeq(norm_expected, norm_actual));
        }

        void trilinos_mm() {
            auto n = 10;
            TpetraMatrixd A;
            A.sparse(layout(comm_, n, n, Traits::determine(), Traits::determine()), 3, 2);
            assemble_laplacian_1D(A);
            TpetraMatrixd B = A;
            TpetraMatrixd C = transpose(A) * B;

            TpetraVectord x(layout(comm_, n, Traits::determine()), 5.);
            TpetraVectord y = C * x;

            const double val = norm2(y);
            utopia_test_assert(approxeq(val, 0.));
        }

        void trilinos_m_tm() {
            auto n = 10;
            auto m = 3;
            TpetraMatrixd A;
            A.sparse(layout(comm_, n, n, Traits::determine(), Traits::determine()), 3, 2);
            assemble_laplacian_1D(A);
            TpetraMatrixd P;
            build_rectangular_matrix(n, m, P);
            TpetraMatrixd B = A * P;
            TpetraMatrixd P_t = transpose(P);
            TpetraMatrixd C_1 = P_t * A;
            TpetraMatrixd C_2 = transpose(P) * A;

            // FIXME write test here
        }

        void trilinos_diag() {
            auto n = 10;
            TpetraMatrixd A;
            A.sparse(layout(comm_, n, n, Traits::determine(), Traits::determine()), 3, 2);
            assemble_laplacian_1D(A);
            TpetraVectord d = diag(A);
            //        disp(A);
            //        disp(d);

            const double val = norm1(d);
            utopia_test_assert(approxeq(val, size(d).get(0) * 2. - 2.));

            TpetraMatrixd D = diag(d);
            TpetraVectord x(layout(comm_, n, Traits::determine()), 1.);
            //        disp(D);
            //        disp(x);
            utopia_test_assert(approxeq(d, D * x));
        }

        void trilinos_diag_rect_matrix() {
            auto n = 10;
            auto m = 3;
            TpetraMatrixd A;
            A.identity(layout(comm_, n, m, Traits::determine(), Traits::determine()), 1.0);
            TpetraVectord d;
            d = diag(A);
            const double val = norm1(d);

            if (!approxeq(val, size(A).get(1) * 1.)) {
                m_utopia_error("diag does not work on for tpetra rectangular matrices in parallel (nor serial)");
            }
        }

        void test_ptap(const int n, const int m) {
            TpetraMatrixd A;
            A.sparse(layout(comm_, n, n, Traits::determine(), Traits::determine()), 3, 2);
            assemble_laplacian_1D(A);

            TpetraMatrixd P;
            build_rectangular_matrix(n, m, P);

            TpetraMatrixd R_2 = transpose(P) * A;
            utopia_test_assert(R_2.is_valid(true));

            // For the moment this is computing (transpose(P) * A) * P
            TpetraMatrixd R = utopia::ptap(A, P);  // equiv: transpose(P) * A * P;

            utopia_test_assert(R.is_valid(true));

#ifdef WITH_PETSC
            // using petsc to test trilinos

            PetscMatrix A_petsc;
            A_petsc.sparse(layout(comm_, n, n, PetscTraits::determine(), PetscTraits::determine()), 3, 2);
            assemble_laplacian_1D(A_petsc);

            PetscMatrix P_petsc;
            build_rectangular_matrix(n, m, P_petsc);

            PetscMatrix R_2_petsc = transpose(P_petsc) * A_petsc;
            PetscMatrix R_petsc = utopia::ptap(A_petsc, P_petsc);

            PetscMatrix R_tpetra;
            PetscMatrix R_2_tpetra;

            backend_convert_sparse(R_2, R_2_tpetra);
            backend_convert_sparse(R, R_tpetra);

            // disp(R_2_tpetra);
            // disp("-----------------------------");
            // disp(R_2_petsc);

            // write("R_t.mm", R);
            // write("R_p.m", R_petsc);

            // write("R2_t.mm", R_2);
            // write("R2_p.m", R_2_petsc);

            double diff_2 = norm2(R_2_petsc - R_2_tpetra);
            double diff = norm2(R_petsc - R_tpetra);

            utopia_test_assert(approxeq(diff_2, 0.));
            utopia_test_assert(approxeq(diff, 0.));
#endif  // WITH_PETSC
        }

        void test_rap(const int n, const int m) {
            TpetraMatrixd A;
            A.sparse(layout(comm_, n, n, Traits::determine(), Traits::determine()), 3, 2);
            assemble_laplacian_1D(A);

            TpetraMatrixd P;
            build_rectangular_matrix(n, m, P);

            TpetraMatrixd R = transpose(P);

            TpetraMatrixd res = R * A * P;

            // disp(res);
        }

        void trilinos_ptap_square_mat() { test_ptap(10, 10); }

        void trilinos_ptap() {
            // does not work in parallel
            test_ptap(10, 3);
        }

        void trilinos_rap_square_mat() {
            // does not work in parallel
            test_rap(10, 10);
        }

        void trilinos_rap() {
            // does not work in parallel
            test_rap(10, 10);
        }

        void trilinos_cg() {
            SizeType n = 20;
            Poisson1D<TpetraMatrixd, TpetraVectord> problem(n);

            TpetraMatrixd H;
            TpetraVectord g, x;

            x.zeros(layout(comm_, Traits::decide(), n));

            problem.hessian(x, H);
            problem.gradient(x, g);
            g *= 0.0001;

            ConjugateGradient<TpetraMatrixd, TpetraVectord> cg;
            cg.rtol(1e-6);
            cg.atol(1e-6);
            cg.max_it(800);
            cg.solve(H, g, x);

            double diff = norm2(g - H * x);
            ;
            utopia_test_assert(approxeq(diff, 0., 1e-6));
        }

        template <class Matrix, class Vector>
        void test_mg() {
            using TransferT = utopia::Transfer<Matrix, Vector>;
            using IPTransferT = utopia::IPTransfer<Matrix, Vector>;
            using MatrixTransferT = utopia::MatrixTransfer<Matrix, Vector>;

            const static bool verbose = false;
            const static bool use_masks = false;

            using ProblemType = utopia::Poisson1D<Matrix, Vector>;
            MultiLevelTestProblem1D<Matrix, Vector, ProblemType> ml_problem(4, 10, !use_masks);
            // ml_problem.describe();
            // ml_problem.write_matlab("./");

            // auto smoother      = std::make_shared<Jacobi<Matrix, Vector>>();
            auto smoother = std::make_shared<ConjugateGradient<Matrix, Vector, HOMEMADE>>();
            auto coarse_solver = std::make_shared<ConjugateGradient<Matrix, Vector, HOMEMADE>>();

            // smoother->set_preconditioner(std::make_shared< InvDiagPreconditioner<Matrix, Vector> >());
            coarse_solver->set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector>>());
            coarse_solver->max_it(1000);
            // coarse_solver->verbose(true);

            Multigrid<Matrix, Vector> multigrid(smoother, coarse_solver);

            multigrid.max_it(10);
            multigrid.atol(1e-13);
            multigrid.stol(1e-13);
            multigrid.rtol(1e-9);
            multigrid.pre_smoothing_steps(3);
            multigrid.post_smoothing_steps(3);
            multigrid.fix_semidefinite_operators(true);
            multigrid.must_generate_masks(use_masks);
            ;
            multigrid.verbose(verbose);

            multigrid.set_transfer_operators(ml_problem.get_transfer());

            auto funs = ml_problem.get_functions();

            Vector x, g;
            Matrix A;
            funs.back()->get_eq_constrains_values(x);
            funs.back()->gradient(x, g);
            funs.back()->hessian(x, A);

            multigrid.update(make_ref(A));

            if (verbose) {
                multigrid.describe();
            }

            multigrid.apply(g, x);

            double diff0 = norm2(A * x);
            double diff = norm2(g - A * x);
            double rel_diff = diff / diff0;

            utopia_test_assert(rel_diff < 1e-8);
        }

        void trilinos_local_row_view() {
            auto rows = 3;
            auto cols = 4;
            TpetraMatrixd A;
            A.sparse(layout(comm_, rows, cols, Traits::determine(), Traits::determine()), 2, 2);

            {
                Write<TpetraMatrixd> w_A(A);
                Range r = row_range(A);

                for (auto i = r.begin(); i < r.end(); ++i) {
                    A.set(i, i, i);
                    A.set(i, i + 1, i + 1);
                }
            }

            TpetraMatrixd At = transpose(A);

            auto &M = A;

            auto rr = row_range(M);
            for (auto i = rr.begin(); i < rr.end(); ++i) {
                RowView<TpetraMatrixd> row(M, i, true);
                for (auto j = 0; j < row.n_values(); ++j) {
                    int col = row.col(j);
                    int val = row.get(j);
                    utopia_test_assert(col == val);
                }
            }
        }

        void trilinos_range() {
            SizeType n = 10, m = 5;
            SizeType gn = mpi_world_size() * n, gm = mpi_world_size() * m;

            TpetraMatrixd P;
            build_rectangular_matrix(n, m, P);

            Range rr = row_range(P);

            utopia_test_assert(rr.begin() == n * mpi_world_rank());
            utopia_test_assert(rr.end() == n * (mpi_world_rank() + 1));
            utopia_test_assert(gn == size(P).get(0));
            utopia_test_assert(gm == size(P).get(1));
        }

        void trilinos_e_mul() {
            int n = 10;
            TpetraVectord v(layout(comm_, n, Traits::determine()), 1.);
            TpetraVectord ones(layout(v), 1.0);

            TpetraVectord ones_mul_v = e_mul(ones, v);
            v = e_mul(ones, v);

            double sv = sum(v);
            utopia_test_assert(approxeq(sv, n * mpi_world_size()));
        }

        template <class Matrix, class Vector>
        void st_cg_test() {
            typename Vector::SizeType _n = 10;

            SteihaugToint<Matrix, Vector, HOMEMADE> cg;
            cg.set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector>>());
            // cg.set_preconditioner(std::make_shared<IdentityPreconditioner<Vector> >());

            cg.rtol(1e-7);
            cg.atol(1e-6);
            cg.max_it(_n);
            cg.verbose(false);

            Matrix A;
            A.sparse(layout(comm_, Traits::decide(), Traits::decide(), _n, _n), 3, 2);
            assemble_symmetric_laplacian_1D(A, true);

            Vector rhs(row_layout(A), 975.9);

            {
                auto r = range(rhs);

                Write<Vector> w(rhs);
                if (r.inside(0)) {
                    rhs.set(0, 0.0);
                }
                if (r.inside(_n - 1)) {
                    rhs.set(_n - 1, 0.0);
                }
            }

            Vector x(layout(rhs), 0.0);
            cg.solve(A, rhs, x);
            utopia_test_assert(approxeq(rhs, A * x, 1e-5));
        }

        void stcg_pt_test() {
#ifdef WITH_PETSC
            // petsc version
            st_cg_test<PetscMatrix, PetscVector>();
#endif  // WITH_PETSC
            st_cg_test<TpetraMatrixd, TpetraVectord>();
        }

        void trilinos_mg_1D() {
            // if(mpi_world_size() > 1) return;
            // petsc version
#ifdef WITH_PETSC
            test_mg<PetscMatrix, PetscVector>();
#endif  // WITH_PETSC
        // trilinos version
            test_mg<TpetraMatrixd, TpetraVectord>();
        }

        void trilinos_mg() {
            // if(mpi_world_size() > 1) return;

            using MatrixT = utopia::TpetraMatrixd;
            using VectorT = utopia::TpetraVectord;

            // using MatrixT = utopia::PetscMatrix;
            // using VectorT = utopia::PetscVector;

            VectorT rhs;
            MatrixT A, I;

            Multigrid<MatrixT, VectorT> multigrid(std::make_shared<ConjugateGradient<MatrixT, VectorT, HOMEMADE>>(),
                                                  std::make_shared<ConjugateGradient<MatrixT, VectorT, HOMEMADE>>()
                                                  // std::make_shared<SOR<MatrixT, VectorT>>(),
                                                  // std::make_shared<Factorization<MatrixT, VectorT>>()
            );

#ifdef WITH_PETSC

            bool ok = true;
            // FIXME needs trilinos formats but for the moment lets use petsc's
            {
                PetscMatrix petsc_A, petsc_I;
                PetscVector petsc_rhs;

                const std::string folder = Utopia::instance().get("data_path") + "/laplace/matrices_for_petsc";

                ok = read(folder + "/f_rhs", petsc_rhs);
                utopia_test_assert(ok);
                ok = read(folder + "/f_A", petsc_A);
                utopia_test_assert(ok);
                ok = read(folder + "/I_3", petsc_I);
                utopia_test_assert(ok);

                backend_convert_sparse(petsc_I, I);
                backend_convert_sparse(petsc_A, A);
                backend_convert(petsc_rhs, rhs);
            }

            // write("A.mm", A);
            // write("I.mm", I);

            std::vector<std::shared_ptr<MatrixT>> interpolation_operators;
            interpolation_operators.push_back(make_ref(I));

            multigrid.set_transfer_operators(interpolation_operators);
            multigrid.max_it(20);
            multigrid.atol(1e-15);
            multigrid.stol(1e-15);
            multigrid.rtol(1e-15);
            // multigrid.verbose(true);
            multigrid.fix_semidefinite_operators(true);
            multigrid.must_generate_masks(true);
            VectorT x(layout(rhs), 0.0);

            try {
                multigrid.update(make_ref(A));
                ok = multigrid.apply(rhs, x);
                utopia_test_assert(ok);
            } catch (const std::exception &ex) {
                std::cout << ex.what() << std::endl;
                utopia_test_assert(false);
            }

            std::cout << std::flush;

            double diff = norm2(rhs - A * x);
            utopia_test_assert(approxeq(diff, 0., 1e-6));

#endif  // WITH_PETSC
        }

        void trilinos_row_view() {
            TpetraMatrixd A;
            A.sparse(layout(comm_, 4, 4, Traits::determine(), Traits::determine()), 3, 2);
            assemble_laplacian_1D(A);

            auto rr = row_range(A);

            for (auto i = rr.begin(); i != rr.end(); ++i) {
                RowView<TpetraMatrixd> row(A, i);
                utopia_test_assert(row.n_values() >= 2);
                auto col = row.col(0);
                // auto val = row.get(0);

                utopia_test_assert(col == i || col == i - 1 || col == i + 1);
            }

            TpetraMatrixd B;
            B.sparse(layout(A), 3, 2);

            {
                Write<TpetraMatrixd> w_(B, GLOBAL_ADD);

                auto r = row_range(B);
                auto s = size(B);

                for (auto i = r.begin(); i < r.end(); ++i) {
                    B.set(i, i, i);

                    if (i + 1 < s.get(1)) {
                        B.set(i, i + 1, i + 1);
                    }

                    if (i - 1 >= 0) {
                        B.set(i, i - 1, i - 1);
                    }
                }
            }

            each_read(B, [](const SizeType i, const SizeType j, const double val) {
                SizeType j_val = val;
                utopia_test_assert(j_val == val);
            });
        }

        void trilinos_row_view_and_loops() {
            int n = 10;
            int m = 3;

            TpetraMatrixd P;
            build_rectangular_matrix(n, m, P);

            auto rr = row_range(P);

            int nnz = std::min(n, int(size(P).get(1)));
            if (rr.begin() >= nnz) {
                nnz = 0;
            }

            SizeType count = 0;
            each_read(P, [&count](const SizeType i, const SizeType j, const double val) {
                utopia_test_assert(val == 1.);
                ++count;
            });

            utopia_test_assert(nnz == count);

            TpetraMatrixd P_t = transpose(P);

            each_read(P_t, [&count](const SizeType i, const SizeType j, const double val) {
                utopia_test_assert(val == 1.);
                --count;
            });

            if (mpi_world_size() == 1) {
                utopia_test_assert(count == 0);
            }

            // disp(P);
            // disp(P_t);
        }

        void trilinos_each_read_transpose() {
            int n = 10;
            int m = 3;

            TpetraMatrixd P;
            build_rectangular_matrix(n, m, P);

            TpetraMatrixd R = transpose(P);
            TpetraMatrixd R_copy = R;
            R_copy *= 0.;

            TpetraVectord v(col_layout(R), 10.);
            TpetraVectord Rv(row_layout(R), 0.0);

            Rv = R * v;

            double nrv = norm2(Rv);
            utopia_test_assert(nrv > 10.);
        }

        void trilinos_read() {
            TpetraMatrixd m;
            auto path = Utopia::instance().get("data_path") + "/matrixmarket/gre_343_343_crg.mm";
            bool ok = read(path, m);
            utopia_test_assert(ok);
        }

#ifdef WITH_PETSC
        void trilinos_petsc_interop() {
            KSPSolver<TpetraMatrixd, TpetraVectord> solver;

            SizeType n = 20;
            Poisson1D<TpetraMatrixd, TpetraVectord> problem(n);

            TpetraVectord x(layout(comm_, Traits::decide(), n), 0.0), g;
            TpetraMatrixd H;

            problem.gradient(x, g);
            problem.hessian(x, H);

            g *= 0.0001;

            solver.solve(H, g, x);

            PetscMatrix p_mat;
            backend_convert_sparse(H, p_mat);

            // disp(p_mat);
            double diff = norm2(g - H * x);
            utopia_test_assert(approxeq(diff, 0., 1e-8));
        }
#endif  // WITH_PETSC

        void trilinos_structure() {
            auto n = 10;
            TpetraMatrixd A;
            A.sparse(layout(comm_, n, n, Traits::determine(), Traits::determine()), 3, 2);
            assemble_laplacian_1D(A);

            auto expr = structure(A);

            TpetraMatrixd B(expr);
        }

        void trilinos_exp() {
#ifdef WITH_PETSC
            TpetraVectord x(layout(comm_, 10, Traits::determine()), 2.);
            PetscVector y(layout(comm_, 10, PetscTraits::determine()), 2.);

            TpetraVectord ex = exp(x);
            PetscVector ey = exp(y);

            utopia_test_assert(cross_backend_approxeq(ey, ex));

#endif  // WITH_PETSC
        }

        void trilinos_diag_ops() {
            int n = 10;

            TpetraMatrixd m;
            m.sparse(layout(comm_, n, n, Traits::determine(), Traits::determine()), 3, 2);
            assemble_laplacian_1D(m);

            TpetraMatrixd m_copy = m;

            TpetraVectord d(layout(comm_, n, Traits::determine()), -1.);
            TpetraMatrixd D = diag(d);
            m += D;

            TpetraVectord ones(layout(comm_, n, Traits::determine()), 1);

            double sum_m = sum(m * ones);
            double sum_d = sum(d);
            double sum_D = sum(D * ones);
            double sum_m_copy = sum(m_copy * ones);

            utopia_test_assert(approxeq(sum_d, sum_D));
            utopia_test_assert(approxeq(sum_m, sum_d + sum_m_copy));

            TpetraMatrixd m_new = m_copy + D;
            double sum_m_new = sum(m * ones);

            utopia_test_assert(approxeq(sum_m_new, sum_d + sum_m_copy));
        }

        void trilinos_bratu_1D() {
#ifdef WITH_PETSC
            int n = 10;

            auto fun_tpetra = std::make_shared<Bratu1D<TpetraMatrixd, TpetraVectord>>(n);
            auto fun_petsc = std::make_shared<Bratu1D<PetscMatrix, PetscVector>>(n);

            TpetraVectord x_tpetra(layout(comm_, Traits::decide(), n), 1.0);
            PetscVector x_petsc(layout(comm_, PetscTraits::decide(), n), 1.0);

            double val_tpetra = 0.;
            double val_petsc = 0.;

            bool ok = true;
            ok = fun_tpetra->value(x_tpetra, val_tpetra);
            assert(ok);
            ok = fun_petsc->value(x_petsc, val_petsc);
            assert(ok);

            utopia_test_assert(cross_backend_approxeq(x_petsc, x_tpetra));

            utopia_test_assert(approxeq(val_tpetra, val_petsc, 1e-8));

            TpetraVectord grad_tpetra;
            PetscVector grad_petsc;

            ok = fun_tpetra->gradient(x_tpetra, grad_tpetra);
            assert(ok);
            ok = fun_petsc->gradient(x_petsc, grad_petsc);
            assert(ok);

            utopia_test_assert(cross_backend_approxeq(grad_petsc, grad_tpetra));

            // last part fails

            TpetraMatrixd H_tpetra;
            PetscMatrix H_petsc;

            ok = fun_tpetra->hessian(x_tpetra, H_tpetra);
            assert(ok);
            ok = fun_petsc->hessian(x_petsc, H_petsc);
            assert(ok);

            // write("H_p.m", H_petsc);
            // write("H_t.m", H_tpetra);

            PetscMatrix H_converted;
            backend_convert_sparse(H_tpetra, H_converted);

            // write("H_c.m", H_converted);

            utopia_test_assert(cross_backend_approxeq(H_petsc, H_tpetra));

#endif  // WITH_PETSC
        }

        template <class Matrix, class Vector>
        void rmtr_test() {
            using IPTransferT = utopia::IPTransfer<Matrix, Vector>;
            using ProblemType = utopia::Bratu1D<Matrix, Vector>;

            auto problem = MultiLevelTestProblem1D<Matrix, Vector, ProblemType>(2, 10, true);

            auto funs = problem.get_functions();

            Vector x;
            funs.back()->get_eq_constrains_values(x);

            auto tr_strategy_coarse = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE>>();
            tr_strategy_coarse->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector>>());
            tr_strategy_coarse->atol(1e-12);
            tr_strategy_coarse->rtol(1e-12);

            auto tr_strategy_fine = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE>>();
            tr_strategy_fine->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector>>());
            tr_strategy_fine->atol(1e-12);
            tr_strategy_fine->rtol(1e-12);

            auto rmtr = std::make_shared<RMTR_l2<Matrix, Vector, GALERKIN>>(problem.n_levels());
            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
            rmtr->set_fine_tr_strategy(tr_strategy_fine);

            rmtr->set_transfer_operators(problem.get_transfer());

            rmtr->max_it(1000);
            rmtr->max_coarse_it(1);
            rmtr->max_QP_smoothing_it(1);
            rmtr->delta0(1);
            rmtr->atol(1e-6);
            rmtr->rtol(1e-10);
            rmtr->grad_smoothess_termination(0.000001);

            rmtr->verbose(false);
            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);
            rmtr->set_functions(funs);
            bool ok = rmtr->solve(x);

            utopia_test_assert(ok);
        }

        void trilinos_replace_value() {
            TpetraMatrixd A;
            A.sparse(layout(comm_, 2, 3, Traits::determine(), Traits::determine()), 4, 4);

            auto gs = size(A);
            auto r = row_range(A);

            {
                Write<TpetraMatrixd> w_(A);

                for (auto i = r.begin(); i < r.end(); ++i) {
                    // A.set(i, i, 1.);
                    A.set(i, i + 1, 2.);
                    A.set(i, 0, 3.);
                    A.set(i, device::max(SizeType(gs.get(1) - i - 1), SizeType(0)), 4.);
                }
            }

            {
                Write<TpetraMatrixd> w_(A);

                for (auto i = r.begin(); i < r.end(); ++i) {
                    A.set(i, 0, -3.);
                    A.set(i, device::max(SizeType(gs.get(1) - i - 1), SizeType(0)), -4.);
                }
            }
        }

        void trilinos_copy_write() {
            TpetraMatrixd P;
            build_rectangular_matrix(10, 5, P);

            P = transpose(P);

            TpetraMatrixd P2 = P;
            P2 *= 0.;

            {
                Write<TpetraMatrixd> w_(P2);
                each_read(P,
                          [&P2](const SizeType i, const SizeType j, const double value) { P2.set(i, j, value * 2.); });
            }

            TpetraMatrixd A;
            MultiLevelTestProblem1D<TpetraMatrixd, TpetraVectord, Poisson1D<TpetraMatrixd, TpetraVectord>> ml_problem(
                10, 2);
            // A = *ml_problem.interpolators[0];
            auto transfer = ml_problem.get_transfer();

            if (dynamic_cast<MatrixTransfer<TpetraMatrixd, TpetraVectord> *>(transfer.back().get()) != nullptr) {
                auto *mat_transfer =
                    dynamic_cast<MatrixTransfer<TpetraMatrixd, TpetraVectord> *>(transfer.back().get());
                A = mat_transfer->I();
            } else {
                utopia_error("trilinos_copy_write::error...");
            }

            A = transpose(A);
            TpetraMatrixd A2 = A;
            A2 *= 0.;

            {
                Write<TpetraMatrixd> w_(A2);
                each_read(A, [&A2](const SizeType i, const SizeType j, const double value) { A2.set(i, j, 1.); });
            }
        }

#ifdef WITH_PETSC
        void trilinos_copy_write_big() {
            PetscMatrix petsc_P;

            const std::string folder = Utopia::instance().get("data_path") + "/laplace/matrices_for_petsc";
            bool ok = read(folder + "/I_2", petsc_P);
            utopia_test_assert(ok);

            TpetraMatrixd P;
            backend_convert_sparse(petsc_P, P);

            P = transpose(P);

            TpetraMatrixd P2 = P;
            P2 *= 0.;

            {
                Write<TpetraMatrixd> w_(P2);
                each_read(P,
                          [&P2](const SizeType i, const SizeType j, const double value) { P2.set(i, j, value * 2.); });
            }
        }
#endif  // WITH_PETSC

        void trilinos_ghosted() {
            const int n = mpi_world_size() * 2;
            const int off = mpi_world_rank() * 2;

            std::vector<TpetraVectord::SizeType> ghosts{(off + 3) % n};
            TpetraVectord v = ghosted(2, n, ghosts);

            auto r = range(v);

            {
                Write<TpetraVectord> w_v(v);
                for (auto i = r.begin(); i != r.end(); ++i) {
                    v.set(i, i);
                }
            }

            // synchronize(v);

            {
                Read<TpetraVectord> r_v(v);
                std::vector<SizeType> index{(off + 3) % n};
                std::vector<double> values;
                v.get(index, values);
                utopia_test_assert(index[0] == SizeType(values[0]));
            }

            // disp(v);
        }

        void trilinos_rmtr() {
#ifdef WITH_PETSC
            // petsc version
            rmtr_test<PetscMatrix, PetscVector>();
#endif  // WITH_PETSC

            rmtr_test<TpetraMatrixd, TpetraVectord>();
        }

        void trilinos_matrix_norm() {
            TpetraMatrixd m;
            m.identity(layout(comm_, 10, 10, Traits::determine(), Traits::determine()), 1.0);

            double nm = norm2(m);
            utopia_test_assert(approxeq(nm, std::sqrt(1. * size(m).get(0))));
        }

#ifdef HAVE_BELOS_TPETRA

        void trilinos_belos() {
            {
                m_utopia_warning(
                    "TrilinosTest::trilinos_belos commented out because of excpetion. Fix and remove this fallback.");
                return;
            }

            std::string xml_file = Utopia::instance().get("data_path") + "/xml/UTOPIA_belos.xml";

            BelosSolver<TpetraMatrixd, TpetraVectord> solver;
            // solver.read_xml(xml_file);
            solver.import("linear-solver", Utopia::instance().get("data_path") + "/json/belos.json");

            Poisson1D<TpetraMatrixd, TpetraVectord> fun(10);
            TpetraVectord x = fun.initial_guess();
            TpetraVectord g;
            TpetraMatrixd A;

            fun.gradient(x, g);
            fun.hessian(x, A);

            g *= 0.0001;

            double diff0 = norm2(g - A * x);
            solver.solve(A, g, x);
            double diff = norm2(g - A * x);

            utopia_test_assert(approxeq(diff / diff0, 0., 1e-6));
        }

#endif  // HAVE_BELOS_TPETRA

#ifdef HAVE_AMESOS2_KOKKOS

        void trilinos_amesos2() {
            std::string xml_file = Utopia::instance().get("data_path") + "/xml/UTOPIA_amesos.xml";

            Amesos2Solver<TpetraMatrixd, TpetraVectord> solver;
            solver.read_xml(xml_file);

            Poisson1D<TpetraMatrixd, TpetraVectord> fun(10);
            TpetraMatrixd A;
            TpetraVectord x, g;
            x = fun.initial_guess();
            fun.get_rhs(g);
            fun.hessian(x, A);

            g *= 0.0001;

            double diff0 = norm2(g - A * x);

            solver.solve(A, g, x);
            double diff = norm2(g - A * x);

            utopia_test_assert(approxeq(diff / diff0, 0., 1e-6));
        }

#endif  // HAVE_AMESOS2_KOKKOS

#ifdef WITH_PETSC
        void trilinos_transform() {
            PetscMatrix petsc_P;

            const std::string folder = Utopia::instance().get("data_path") + "/laplace/matrices_for_petsc";
            bool ok = read(folder + "/I_2", petsc_P);
            utopia_test_assert(ok);

            TpetraMatrixd P;
            backend_convert_sparse(petsc_P, P);

            double sum_P = sum(P);
            P = transpose(P);

            each_apply(P, [](const double value) -> double { return value * 2.; });

            double sum_P_2 = sum(P);

            utopia_test_assert(approxeq(sum_P * 2., sum_P_2));

            each_transform(P, [](const SizeType i, const SizeType j, const double) -> double { return j; });

            each_read(P, [](const SizeType i, const SizeType j, const double value) {
                if (j != SizeType(value)) {
                    std::cout << i << " " << j << " " << value << std::endl;
                }

                utopia_test_assert(j == SizeType(value));
            });
        }
#endif  // WITH_PETSC

        void trilinos_set_zeros() {
            TpetraMatrixd m;
            m.identity(layout(comm_, 10, 10, Traits::determine(), Traits::determine()), 1.0);
            using SizeT = UTOPIA_SIZE_TYPE(TpetraMatrixd);
            auto rr = row_range(m);

            std::vector<SizeT> index;
            index.push_back(rr.begin());
            set_zero_rows(m, index, 2.);

            // disp(m);

            Read<TpetraMatrixd> r(m);
            auto val = m.get(rr.begin(), rr.begin());
            utopia_test_assert(val == 2.);
        }

        void trilinos_decompose() {
            TrilinosCommunicator comm;

            SizeType n = 17;
            SizeType n_local = decompose(comm, n);

            SizeType n_sum = comm.sum(n_local);
            utopia_test_assert(n_sum == n);
        }

        void test_global_matrix() {
            SizeType n = 17;
            TpetraMatrixd A;
            A.sparse(layout(comm_, n, n, Traits::determine(), Traits::determine()), 3, 2);
            assemble_laplacian_1D(A);

            TpetraVectord x(col_layout(A), 2.0);
            TpetraVectord b = A * x;

            const Scalar zero = sum(b);

            utopia_test_assert(approxeq(zero, 0.0));
        }

        template <typename T>
        using ArrayT = Teuchos::ArrayRCP<T>;
        static void trilinos_crs_construct() {
            const SizeType n_rows = 2;
            const SizeType n_cols = 3;

            ArrayT<std::size_t> row_ptr(n_rows + 1);
            row_ptr[0] = 0;
            row_ptr[1] = 2;
            row_ptr[2] = 4;

            ArrayT<LocalSizeType> columns(row_ptr[n_rows]);
            columns[0] = 0;
            columns[1] = 1;
            columns[2] = 1;
            columns[3] = 2;

            ArrayT<Scalar> values(row_ptr[n_rows]);
            values[0] = 1;
            values[1] = -1;
            values[2] = -1;
            values[3] = 1;

            TpetraMatrixd A = utopia::crs(n_rows, n_cols, row_ptr, columns, values);
        }

        void run() {
            UTOPIA_RUN_TEST(stcg_pt_test);
            UTOPIA_RUN_TEST(trilinos_structure);
            UTOPIA_RUN_TEST(trilinos_build);
            UTOPIA_RUN_TEST(trilinos_build_identity);
            UTOPIA_RUN_TEST(trilinos_accessors);
            UTOPIA_RUN_TEST(trilinos_vec_scale);
            UTOPIA_RUN_TEST(trilinos_mat_scale);
            UTOPIA_RUN_TEST(trilinos_vec_axpy);
            UTOPIA_RUN_TEST(trilinos_vec_minus);
            UTOPIA_RUN_TEST(trilinos_mv);
            UTOPIA_RUN_TEST(trilinos_mm);
            UTOPIA_RUN_TEST(trilinos_m_tm);
            UTOPIA_RUN_TEST(trilinos_diag);
            UTOPIA_RUN_TEST(trilinos_read);

            UTOPIA_RUN_TEST(trilinos_rect_matrix);
            UTOPIA_RUN_TEST(trilinos_e_mul);
            UTOPIA_RUN_TEST(trilinos_row_view);
            UTOPIA_RUN_TEST(trilinos_apply_transpose);
            UTOPIA_RUN_TEST(trilinos_set);
            UTOPIA_RUN_TEST(trilinos_residual);

            UTOPIA_RUN_TEST(trilinos_matrix_access);
            UTOPIA_RUN_TEST(trilinos_matrix_norm);
            UTOPIA_RUN_TEST(trilinos_diag_ops);

            UTOPIA_RUN_TEST(trilinos_rmtr);

            UTOPIA_RUN_TEST(trilinos_transpose);

            UTOPIA_RUN_TEST(trilinos_apply_transpose_explicit);
            UTOPIA_RUN_TEST(trilinos_each_read_transpose);
            UTOPIA_RUN_TEST(trilinos_local_row_view);
            UTOPIA_RUN_TEST(trilinos_ptap_square_mat);

            UTOPIA_RUN_TEST(trilinos_range);
            UTOPIA_RUN_TEST(trilinos_mg_1D);
            UTOPIA_RUN_TEST(trilinos_bratu_1D);
            UTOPIA_RUN_TEST(trilinos_replace_value);
            UTOPIA_RUN_TEST(trilinos_ghosted);
            UTOPIA_RUN_TEST(trilinos_set_zeros);

            // FIXME fails from mpi_world_size() > 3
            UTOPIA_RUN_TEST(trilinos_copy_write);

            ////////////////////////////////////////////
            // test that fail on GPU if the env variables are not set correctly for cuda
            UTOPIA_RUN_TEST(trilinos_exp);
            UTOPIA_RUN_TEST(trilinos_mg);
            UTOPIA_RUN_TEST(trilinos_cg);
            UTOPIA_RUN_TEST(trilinos_ptap);
            ////////////////////////////////////////////

            UTOPIA_RUN_TEST(trilinos_rap);
            UTOPIA_RUN_TEST(trilinos_rap_square_mat);
            UTOPIA_RUN_TEST(trilinos_decompose);
            UTOPIA_RUN_TEST(test_global_matrix);

#ifdef HAVE_BELOS_TPETRA
            UTOPIA_RUN_TEST(trilinos_belos);
#endif  // HAVE_BELOS_TPETRA

            //#ifdef HAVE_AMESOS2_TPETRA
            UTOPIA_RUN_TEST(trilinos_amesos2);
            //#endif //HAVE_AMESOS2_TPETRA

#ifdef WITH_PETSC
            UTOPIA_RUN_TEST(trilinos_transform);
            UTOPIA_RUN_TEST(trilinos_petsc_interop);
            UTOPIA_RUN_TEST(trilinos_copy_write_big);
#endif  // WITH_PETSC

            // Fails on multinode GPU
            UTOPIA_RUN_TEST(trilinos_mat_axpy);
            ////////////////////////////////////////////

            if (mpi_world_size() <= 3) {
                // working up to 3 processes
                UTOPIA_RUN_TEST(trilinos_row_view_and_loops);
            }

            // tests that fail in parallel
            if (mpi_world_size() == 1) {
                UTOPIA_RUN_TEST(trilinos_crs_construct);
            }
            // else {
            //     m_utopia_warning_once("several tests left out for parallel execution");
            // }

            // tests that always fail
            // UTOPIA_RUN_TEST(trilinos_diag_rect_matrix);
        }
    };

    void trilinos_specific() { TrilinosTest().run(); }

    UTOPIA_REGISTER_TEST_FUNCTION(trilinos_specific);
}  // namespace utopia

#endif  // WITH_TRILINOS
