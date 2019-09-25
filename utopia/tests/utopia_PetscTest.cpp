#include "utopia.hpp"
#include "utopia_Testing.hpp"
#include "test_problems/utopia_TestFunctionsND.hpp"
#include "utopia_QuadraticFunction.hpp"
#include "utopia_Device.hpp"
#include "utopia_assemble_laplacian_1D.hpp"
#include "utopia_ZeroRowsToIdentity.hpp"

namespace utopia {

#ifdef WITH_PETSC

    //FIXME
    void petc_optional()
    {
       //  MPI_Comm sub_comm;
       //  MPI_Comm_split(
       //      PETSC_COMM_WORLD,
       //      mpi_world_rank() % 2 == 1,
       //      mpi_world_rank(),
       //      &sub_comm);

       //  int rank;
       //  MPI_Comm_rank(sub_comm, &rank);

       //  {
       //      //optionals only work for this builder at the moment
       //      PetscMatrix m = local_sparse(10, 10, 1, sub_comm, str("my_mat"));

       //      auto r = row_range(m);
       //      {
       //          Write<PetscMatrix> w_m(m);
       //          for(auto i = r.begin(); i != r.end(); ++i) {
       //              m.set(i, i, rank);
       //          }
       //      }


       //     PetscVector vec;
       //     std::string path = Utopia::instance().get("data_path");
       //     read(path + "/RHS_10x10x10_hexa_3D", vec, sub_comm, str("my_vec"));
       //     }

       // MPI_Comm_free(&sub_comm);
    }


    void petsc_reciprocal() {
        // test also  diag
        PetscMatrix A = values(4, 4, 1.0);
        auto rr = row_range(A);

        {
            Write<PetscMatrix> w(A);

            if(rr.inside(1)) { A.set(1, 1, 99); }
            if(rr.inside(2)) { A.set(2, 2, 77); }
        }

        PetscVector diag_A = diag(A);
        PetscVector v_expected = values(4, 1.0);
        {
            Write<PetscVector> w(v_expected);
            if(rr.inside(1)) { v_expected.set(1, 99); }
            if(rr.inside(2)) { v_expected.set(2, 77); }
        }

        utopia_test_assert(approxeq(v_expected, diag_A));

        A = diag(diag_A);

        diag_A = values(4, 9.0);
        PetscVector result_petsc = 5 / diag_A;
        v_expected = values(4, 5.0/9.0);
        utopia_test_assert(approxeq(v_expected, result_petsc));

        PetscVector v_expected_test= values(4.0, 1.0);

        {
            Write<PetscVector> w(v_expected_test);
            Range rr=range(v_expected_test);
            for (auto ii=rr.begin(); ii<rr.end(); ++ii){
                 v_expected_test.set(ii, 2 * ii);
        }

        }

        PetscVector p = power(v_expected_test, 3.0);
    }

    void petsc_axpy() {

        {
            //! [axpy (petsc)]
            int n = 10;
            PetscVector x = values(n, 1.0);
            PetscVector y = values(n, 2.0);
            PetscVector expected = values(n, 1.0 * 0.1 + 2.0);
            const double alpha = 0.1;

            PetscVector result = alpha * x + y;
            utopia_test_assert(norm1(result - expected) < 1e-5);
            //! [axpy (petsc)]
        }

        ///////////////////////////////////

        {
            const PetscInt n = 10;
            const PetscInt m = 20;

            PetscMatrix X = values(m, n, 1.0);
            PetscMatrix Y = dense_identity(m, n);
            const double alpha = 4;
            PetscMatrix result = alpha * X + Y;
            double actual = sum(result);
            utopia_test_assert(approxeq(810., actual));
        }

        ///////////////////////////////////

        {
            const int m = mpi_world_size() * 3;
            const int n = mpi_world_size() * 2;

            PetscMatrix X = 2.  * identity(m, n);
            PetscMatrix Y = 0.1 * identity(m, n);
            double alpha = 4.;
            PetscMatrix res = alpha * X + Y;
        }

        {
            const int n = mpi_world_size() * 10;

            PetscMatrix expected = sparse(n, n, 2);
            PetscMatrix X = identity(n, n);
            PetscMatrix Y = sparse(n, n, 1);

            {
                Write<PetscMatrix> w_e(expected), w_Y(Y);
                Range r = row_range(expected);

                for(auto i = r.begin(); i < r.end(); ++i) {
                    expected.set(i, i, 0.1);
                    expected.set(i, n-i-1, 0.1);
                }

                for(auto i = r.begin(); i < r.end(); ++i) {
                    Y.set(i, n-i-1, 0.1);
                }
            }

            PetscMatrix actual = 0.1 * X + Y;

            utopia_test_assert(approxeq(expected, actual));
        }
    }

    void petsc_vector_accessors() {
        int mult;
        MPI_Comm_size(PETSC_COMM_WORLD, &mult);
        const PetscInt n = 10 * mult;
        PetscVector x = values(n, 1, 1);
        PetscVector y = values(n, 1, 2);

        //for parallel data structures (works also for serial ones)
        Range xr = range(x);
        const PetscInt xb = xr.begin();


        //The use of Read/Write/ReadAndWrite locks is important for universal compatibility with all (parallel) backends
        { //scoped write lock
            Write<PetscVector> writing(x, GLOBAL_INSERT);
            //the petsc way
            std::vector<PetscScalar> values{-1., 10., 10., 10, -1.};
            std::vector<PetscInt> index{xb, xb + 2, xb + 3, xb + 4, xb + 9};
            x.set(index, values);
        }

        { //scoped read lock
            Read<PetscVector> reading(x);
            PetscScalar v = x.get(xb + 1);
            utopia_test_assert(approxeq(1, v));
        }

        { //scoped read-write lock
            ReadAndWrite<PetscVector> readingAndWriting(x);
            x.set(xb + 1, x.get(xb + 0) * x.get(xb + 1));

        }

        each_read(x, [](const SizeType i, const double value) {
            switch (i % 10) {
                case 0:
                case 1:
                case 9:
                utopia_test_assert(approxeq(-1, value));
                break;
                case 2:
                case 3:
                case 4:
                utopia_test_assert(approxeq(10, value));
                break;
                case 5:
                case 6:
                case 7:
                case 8:
                utopia_test_assert(approxeq(1, value));
            }
        });


        x.set(0.);
    }


    void petsc_matrix_accessors() {
        int mult;
        MPI_Comm_size(PETSC_COMM_WORLD, &mult);
        const PetscInt n = 10 * mult;
        PetscMatrix x = values(n, 2, 1);
        PetscMatrix y = values(n, 2, 2);

        //for parallel data structures (works also for serial ones)
        Range xr = row_range(x);
        const PetscInt xb = xr.begin();

    //  Petsc matrix does not support ReadAndWrite
    //    //The use of Read/Write locks is important for universal compatibility with all (parallel) backends
        { //scoped write lock
            Write<PetscMatrix> writing(x);
            x.set(xb, 0, -1);
            x.set(xb + 9, 0, -1);
        }

        { //scoped read lock
            Read<PetscMatrix> reading(x);
            PetscScalar v = x.get(xb + 1, 0);
            utopia_test_assert(approxeq(1, v));
        }
    }

    void petsc_sparse_matrix_accessors() {

        //! [Read write matrix]
        int mult = mpi_world_size();
        const PetscInt n = 10 * mult;
        PetscMatrix x = sparse(n, 10, 7);
        PetscMatrix y = sparse(n, 10, 2);

        // For parallel data structures (works also for serial ones. Adopt paridgm for code portability)
        Range xr = row_range(x);
        const PetscInt xb = xr.begin();

        // Petsc matrix does not support ReadAndWrite only separate Read and Write
        // The use of Read/Write locks is important for universal compatibility with all (parallel) backends
        {
            Write<PetscMatrix> w(y);
            y.set(xb, 0, 0);
            y.set(xb + 1, 1, 1);
            y.set(xb + 2, 0, 4);
            y.set(xb + 3, 0, 6);
            y.set(xb + 4, 0, 8);
            y.set(xb + 9, 0, 18);
            y.set(xb + 9, 1, 19);
        }

        {
            Read<PetscMatrix> r(x);
            PetscScalar v = x.get(xb + 1, 1);
            utopia_test_assert(approxeq(3, v));
        }

        x = sparse(n, n, 1);
        {
            Write<PetscMatrix> w(x);

            x.set(xb, xb, 1);
            x.set(xb + 1, xb + 1, 2);
        }

        {
            Read<PetscMatrix> r(x);
            PetscScalar v = x.get(xb + 1, xb +1);
            utopia_test_assert(approxeq(2, v));
            v = x.get(xb, xb);
            utopia_test_assert(approxeq(1, v));
        }

        //! [Read write matrix]
    }

    void petsc_mv()
    {
        const SizeType n = 10;
        const SizeType m = 20;

        //creates a dense matrix
        PetscVector v   = values(n, 1, 1);
        PetscMatrix mat = identity(m, n);
        PetscVector result = mat * v;

        PetscVector expected = values(m, 1, 1);
        {
            auto r = range(expected);
            Write<PetscVector> w(expected);
            SizeType e_begin = std::max(n, r.begin());
            SizeType e_end   = std::min(m, r.end());
            for (SizeType i = e_begin; i < e_end; i++) {
                expected.set(i, 0);
            }
        }

        utopia_test_assert(approxeq(expected, result));
    }

    void petsc_copy()
    {
        const PetscInt n = 10;
        const PetscInt m = 20;

        PetscVector v1 = values(n, 1, 1);
        PetscVector v2 = v1;
        utopia_test_assert(approxeq(v1, v2));

        PetscMatrix m1 = identity(m, n);
        PetscMatrix m2 = m1;
        utopia_test_assert(approxeq(m1, m2));
    }

    void petsc_wrapper()
    {
        PetscMatrix m = identity(2, 2);
        auto expected_ptr = raw_type(m);
        Mat pmat = raw_type(m);
        PetscMatrix wmat;
        wrap(pmat, wmat);
        utopia_test_assert( raw_type(wmat) == expected_ptr );
    }

    void petsc_vector_composite() {
        const PetscInt n = 10;
        PetscVector v1 = values(n, 1, 1.0);
        PetscVector v2 = values(n, 1, 2.0);
        PetscVector v3 = values(n, 1, 0.2);

        auto expr = v1 * 0.1 + v2 + v1 - v3;
        PetscVector vresult = expr;
        PetscVector expected = values(n, 1, 2.9);
        utopia_test_assert(approxeq(expected, vresult));

        PetscScalar value = norm2(expr);
        utopia_test_assert(approxeq(std::sqrt(10*2.9*2.9), value));

        const double valueSum = sum(expr);
        utopia_test_assert(std::abs(valueSum - 29) < 1e-14);
    }

    void petsc_matlab_connection() {
        PetscVector rhs, sol, M_rhs;
        PetscMatrix K, M;
        std::string path = Utopia::instance().get("data_path");

        read(path + "/RHS_10x10x10_hexa_3D", rhs);
        read(path + "/K_hexa_10x10x10_3D", K);
        read(path + "/M_hexa_10x10x10_3D", M);

        M_rhs = M * rhs;
        sol = local_zeros(local_size(M_rhs));

        // Linear solver
        // auto lsolver = std::make_shared<Factorization<PetscMatrix, PetscVector>>();
        auto lsolver = std::make_shared<BiCGStab<PetscMatrix, PetscVector>>();
        lsolver->solve(K, M_rhs, sol);

        double diff = norm2(M_rhs - K * sol);
        utopia_test_assert(diff < 1e-6);

        // Non-linear solver
        Newton<PetscMatrix, PetscVector> newton(lsolver);
        newton.enable_differentiation_control(false);

        QuadraticFunction<PetscMatrix, PetscVector> fun(make_ref(K), make_ref(M_rhs));
        sol = local_zeros(local_size(M_rhs));
        newton.solve(fun, sol);

        diff = norm2(M_rhs - K * sol);
        utopia_test_assert(diff < 1e-6);
    }


    void petsc_matrix_composite() {
        // std::cout << "begin: petsc_matrix_composite_test" << std::endl;

        const PetscInt n = 10;
        PetscMatrix m1 = values(n, n, 1.0);
        PetscMatrix m2 = values(n, n, 2.0);
        PetscMatrix m3 = values(n, n, 0.2);

        PetscMatrix result = m1 * 0.1 + m2 + m1 - m3;
        PetscMatrix expected = values(n, n, 2.9);
        utopia_test_assert(approxeq(expected, result));

        // std::cout << "end: petsc_matrix_composite_test" << std::endl;
    }

    // void petsc_view()
    // {
    //     //! [Global views]
    //     const PetscInt offset = mpi_world_size();
    //     const PetscInt n = offset + 5;
    //     const PetscInt m = offset + 10;
    //     PetscMatrix m1 = identity(m, n);
    //     PetscMatrix m2 = m1.range(
    //         1, offset + 6,
    //         0, offset + 5);


    //     each_read(m2, [](const SizeType i, const SizeType j, const double value) {
    //         if (i + 1 == j) {
    //             if(!approxeq(1, value)) {
    //                 std::cout << i << ", " << j << " -> " << value << std::endl;
    //             }
    //             utopia_test_assert(approxeq(1, value));
    //         } else {
    //             utopia_test_assert(approxeq(0, value));
    //         }
    //     });

    //     //! [Global views]
    // }

    void petsc_mat_tests()
    {
        PetscBool assembled;
        PetscMatrix M = zeros(5, 5);
        MatAssembled(raw_type(M), &assembled);
        utopia_test_assert(!empty(M));

        PetscMatrix A;
        utopia_test_assert(empty(A));
    }

    void petsc_vec_tests()
    {
        PetscVector v = zeros(5);
        utopia_test_assert(!empty(v));

        PetscVector e;
        utopia_test_assert(empty(e));
    }

    void petsc_read_and_write()
    {
        //! [Input and output (petsc)]

        PetscVector x = values(3, 1, 3.0);

        // Write vector to disk
        write("test_vec.txt", x);

        // Display values on the console
        // disp(x);

        PetscVector y;
        // Read vector from disk
        read("test_vec.txt", y);

        utopia_test_assert(approxeq(x, y));

        PetscMatrix m = identity(3, 3);

        PetscMatrix w;
        // Write matrix to disk
        write("test_mat.txt", m);

        // Read matrix from disk
        read("test_mat.txt", w);

        utopia_test_assert(approxeq(m, w));

        //! [Input and output (petsc)]
    }



    void petsc_to_blas()
    {
#ifdef WITH_BLAS

    PetscVector x = zeros(16);

    Range xr = range(x);
    const PetscInt xb = xr.begin();

    BlasVectord y = values(xr.extent(), 2.0);

    {
        Read<BlasVectord> r_y(y);
        Write<PetscVector> w_x(x);
        for (SizeType i = 0; i < xr.extent(); ++i)
            x.set(xb + i, y.get(i));
    }

    PetscVector expected = values(16, 2.0);
    utopia_test_assert(approxeq(expected, x));

#endif //WITH_BLAS
    }



    //FIXME does not work
    void petsc_local_entities() {
        std::cout << "Begin: petsc_local_entities_test." << std::endl;
        PetscMatrix matrix;

        int rank = matrix.comm().rank();
        int size = matrix.comm().size();

        PetscVector ones = local_values(2, rank);

        matrix = local_identity((rank && rank < size - 1) ? 4 : 3, 2);

        // disp(matrix);

        if (size > 1) {
            Write<PetscMatrix> write(matrix);
            Range mr = row_range(matrix);
            Range vr = range(ones);

            if (mr.begin() != 0 && mr.end() != matrix.size().get(0)) {
                matrix.set(mr.begin() + 2, vr.begin() - 1, 1.0);
                matrix.set(mr.begin() + 3, vr.end(), 1.0);
            }

            if (mr.begin() == 0) {
                matrix.set(mr.begin() + 2, vr.end(), 1.0);
            }
            if (mr.end() == matrix.size().get(0)) {
                matrix.set(mr.begin() + 2, vr.begin() - 1, 1.0);
            }
        }

        PetscVector result = matrix * ones;

        // disp(matrix);
        // disp(result);
        std::cout << "End: petsc_local_entities_test." << std::endl;
    }



    void petsc_conversion() {
        // std::cout << "Begin: petsc_conversion_test." << std::endl;
        PetscVector vec = values(10, 1.0);


        PetscVector dvec;
        convert(raw_type(vec), dvec);
        utopia_test_assert(approxeq(vec, dvec));
        convert(dvec, raw_type(vec));
        utopia_test_assert(approxeq(vec, dvec));


        PetscMatrix mat_dense = values(2, 3, 1);
        PetscMatrix dmat_dense;
        convert(raw_type(mat_dense), dmat_dense);
        utopia_test_assert(approxeq(mat_dense, dmat_dense));
        convert(dmat_dense, raw_type(mat_dense));

        PetscMatrix mat = sparse(3, 3, 2);
        {
            Write<PetscMatrix> w(mat);
            mat.set(0, 0, 1);
            mat.set(0, 1, 2);
            mat.set(1, 1, 3);
            mat.set(2, 2, 4);
        }

        PetscMatrix dmat;
        convert(raw_type(mat), dmat);
        convert(dmat, raw_type(mat));

        // std::cout << "End: petsc_conversion_test." << std::endl;
    }


    void petsc_factory_and_operations()
    {
        const int n = mpi_world_size() * 3;
        PetscMatrix m = sparse(n, n, 3);
        {
            auto r = row_range(m);

            Write<PetscMatrix> w(m);
            for(auto i = r.begin(); i < r.end(); ++i) {
                if(i == n - 1 || i == 0) {
                    m.set(i, i, 1.);
                    continue;
                }

                m.set(i, i,  2.);
                m.set(i, i - 1, -1.);
                m.set(i, i + 1, -1.);
            }
        }

        PetscVector c = (m + 0.1 * identity(n, n)) * values(n, 0.5);
    }

    void maria()
    {
        if(mpi_world_size() != 2) return;

        int master_sizes[2] = {
            15, 10
        };

        int slave_sizes[2] = {
            25, 20
        };

        const std::string data_path = Utopia::instance().get("data_path") + "/master_and_slave";

        const auto r = mpi_world_rank();

        PetscMatrix B = local_sparse(slave_sizes[r], master_sizes[r], 2);
        PetscMatrix D = local_sparse(slave_sizes[r], slave_sizes[r],  2);

        PetscViewer fd;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, (data_path + "/master_matrix").c_str(), FILE_MODE_READ, &fd);
        MatLoad(raw_type(B),fd);
        PetscViewerDestroy(&fd);

        PetscViewerBinaryOpen(PETSC_COMM_WORLD, (data_path + "/slave_matrix").c_str(), FILE_MODE_READ, &fd);
        MatLoad(raw_type(D), fd);
        PetscViewerDestroy(&fd);
        PetscMatrix Dinv = diag(1./sum(D, 1));
        PetscMatrix T = Dinv * B;

        PetscMatrix Tt = transpose(Dinv * B);

        // disp(Tt.size());
        // disp(T.size());

        PetscMatrix m = transpose(Dinv + D * Dinv);
    }

    void local_diag_block()
    {
        PetscMatrix a = sparse(4, 4, 3);

        {
            Write<PetscMatrix> write(a);
            a.set(0, 0, 1);
            a.set(1, 1, 2);
            a.set(2, 2, 3);
            a.set(3, 3, 4);
            a.set(1, 3, 4);
            a.set(3, 1, 4);
        }

        each_read(a, [](const SizeType i, const SizeType j, const double value) {
            if (i == j)
                utopia_test_assert(approxeq(i + 1, value));
            else
                utopia_test_assert(approxeq(4, value));
        });

        PetscMatrix b(PETSC_COMM_SELF);
        b = local_diag_block(a);
        each_read(b, [](const SizeType i, const SizeType j, const double value) {
            if (i == j)
                utopia_test_assert(approxeq(i + 1, value));
            else
                utopia_test_assert(approxeq(4, value));
        });
    }

    void petsc_each_sparse_matrix()
    {
        const SizeType n = mpi_world_size();
        PetscMatrix a = sparse(n, n, 3);

        {
            auto r = row_range(a);
            Write<PetscMatrix> write(a);

            for(auto i = r.begin(); i < r.end(); ++i) {
                a.set(i, i, i);
            }
        }

        each_read(a, [](const SizeType i, const SizeType /*j*/, const double value) {
            utopia_test_assert(approxeq(i, value));
        });
    }

    void petsc_test_ptap()
    {
        const int n = mpi_world_size() * 3;
        PetscMatrix P = identity(n, n);
        PetscMatrix A = identity(n, n);
        PetscMatrix C = identity(n, n);
        PetscMatrix PtAP, ABC, PAPt;

        //The next line is equivalent to this:  PtAP = ptap(P, A); since it is pattern matched
        PtAP = transpose(P) * A * P;
        ABC = A * P * C;

        PetscMatrix expected = identity(n, n);
        utopia_test_assert(approxeq(expected, PtAP));
        utopia_test_assert(approxeq(expected, ABC));
    }

    void petsc_test_rart()
    {
        const int n = mpi_world_size() * 4;
        const int m = mpi_world_size() * 3;

        PetscMatrix R = identity(m, n);
        PetscMatrix A = identity(n, n);
        PetscMatrix C = identity(n, m);
        PetscMatrix RARt, ABC;

        RARt = R * A * transpose(R);
        ABC = R * A * C;

        PetscMatrix expected = identity(m, m);
        utopia_test_assert(approxeq(expected, RARt));
        utopia_test_assert(approxeq(expected, ABC));
    }

    void petsc_matrix_composition()
    {
        const SizeType n = mpi_world_size() * 2;
        PetscMatrix m = identity(n, n);
        PetscVector v = values(n, 2.);
        auto expr = abs(transpose(0.1 * (m * m) - m) * (m) * v);

        PetscVector res = expr;
        PetscVector expected = values(n, 0.9 * 2.);
        utopia_test_assert(approxeq(expected, res));
    }



    // template<class Left, class Right>
    // utopia::Assign<Left, Right> assign(utopia::Expression<Left> &left,
    //                                    const utopia::Expression<Right> &right) {
    //     return utopia::Assign<Left, Right>(left.derived(), right.derived());
    // }

    template<class Expr>
    void eval(const Expr &expr)
    {
        utopia::Eval<Expr>::apply(expr);
    }

    void petsc_new_eval()
    {
        const int n = mpi_world_size() * 2;
        double alpha = 0.5;
        PetscVector x = values(n, 1.0);
        PetscVector y = values(n, 1.0);

        eval(construct(y, alpha * x + y));
        PetscMatrix m;
        eval(construct(m, outer(x, y)));
        // disp(m);

    //    PetscVector v2 = values(11, 0.0);
    //    auto view = v2.range(0, 10);
    //    eval(assign(view, to_wrapper(x)));
    //    disp(v2);

    }

    void petsc_precond()
    {
        int n = mpi_world_size() * 10;
        PetscVector d   = values(n, 2.0);
        PetscMatrix dm = diag(d);
        PetscVector rhs = values(n, 1.);
        PetscVector sol = zeros(n);

        PetscVector v = diag(d) * rhs;

        auto m = std::make_shared<PetscMatrix>(2.*identity(n, n));

        auto precond = std::make_shared< InvDiagPreconditioner<PetscMatrix, PetscVector> >();
        precond->update(m);

        auto cg = std::make_shared< ConjugateGradient<PetscMatrix, PetscVector> >();
        cg->set_preconditioner(precond);
        cg->solve(*m, rhs, sol);

        PetscVector expected = values(n, 0.5);
        utopia_test_assert(approxeq(expected, sol));


        // PetscViewer   viewer;
        // PetscViewerDrawOpen(PETSC_COMM_WORLD,0,"",300,0,300,300,&viewer);
        // VecView(raw_type(sol),viewer);
        // PetscViewerDestroy(&viewer);

    }


    void petsc_tensor_reduction()
    {
        int n = mpi_world_size() * 10;
        PetscMatrix mat = identity(n, n);

        //summing columns
        PetscVector v = sum(mat, 1);
        PetscVector expected = values(n, 1);
        utopia_test_assert(approxeq(expected, v));
    }

    void petsc_inverse()
    {
        if(mpi_world_size() == 1) {
            PetscMatrix mat = dense_identity(3, 3);

            {
                Write<PetscMatrix> w(mat);
                mat.set(0, 1, 2);
                mat.set(1, 0, 2);
            }

            PetscMatrix inv_mat  = inv(mat);
            PetscMatrix actual   = inv_mat * mat;
            PetscMatrix expected = dense_identity(3, 3);

            utopia_test_assert( approxeq(actual, expected) );
        }
    }

    void petsc_hardcoded_cg()
    {
        const int n = mpi_world_size() * 2;
        const int i_max = 2;
        const double eps = 1e-8;

        PetscVector x  = zeros(n);
        PetscVector b  = values(n, 1.0);
        PetscMatrix A = identity(n, n);
        PetscMatrix M = diag(diag(A));

        /////////////////////
        PetscVector r = zeros(n), d = zeros(n), q = zeros(n), s = zeros(n);
        double delta_new = 0, delta_0 = 0, delta_old = 0, alpha = 0, beta = 0;

        int i = 0;
        r = b - A * x;
        solve(M, r, d);

        delta_new = dot(r, d);
        delta_0   = delta_new;

        while(i < i_max && delta_new > eps * eps * delta_0) {
            q = A * d;
            alpha = delta_new/dot(d, q);
            x += alpha * d;

            if(i % 50 == 0) {
                r = b - A * x;
            } else {
                r -= alpha * q;
            }

            solve(M, r, s);
            delta_old = delta_new;
            delta_new = dot(r, s);
            beta      = delta_new/delta_old;

            d = s + beta * d;
            ++i;
        }

        // disp(x);
    }

    void petsc_ksp_precond_delegate()
    {
        const int n = 10;
        if(mpi_world_size() > n) return;

        TestFunctionND_1<PetscMatrix, PetscVector> fun(n);

        ConjugateGradient<PetscMatrix, PetscVector> cg;
        cg.set_preconditioner(std::make_shared<DelegatePreconditioner<PetscMatrix, PetscVector> >());
        Newton<PetscMatrix, PetscVector> newton(make_ref(cg));

        PetscVector x = zeros(n);
        newton.solve(fun, x);
        // disp(x);

        PetscVector expected = values(n, 0.468919);
        utopia_test_assert(approxeq(expected, x));

        cg.set_preconditioner(std::make_shared<PointJacobi<PetscMatrix, PetscVector> >());
        x = zeros(n);

        newton.solve(fun, x);
        utopia_test_assert(approxeq(expected, x));
    }

    void petsc_is_nan_or_inf()
    {
        const int n     = 10;
        if(mpi_world_size() > n) return;

        PetscVector denom  = zeros(n);
        PetscVector nom    = values(n, 1.0);

        PetscVector sol    = nom/denom;

        {
            Write<PetscVector> w(sol);
            if(range(sol).inside(2)) { sol.set(2, 1); }
        }

        //FIXME
        if(!has_nan_or_inf(sol)) {
            utopia_error("petsc_is_nan_or_inf: failed. Known problem related to petsc-3.9. They disallow divisions by 0 and instead set the value to 0");
        }

        utopia_test_assert(!has_nan_or_inf(denom));

    }

    void petsc_mat_mul_add()
    {
        int n = mpi_world_size() * 2;
        for (size_t i = 0; i < 50; i++) {
            double d1 = rand()/(double)RAND_MAX,
            d2 = rand()/(double)RAND_MAX,
            d3 = rand()/(double)RAND_MAX;

            PetscMatrix m = values(n, n, d1);
            PetscVector v1 = values(n, d2);
            PetscVector v2 = values(n, d3);
            PetscVector r1, r2, r3;

            r1 = m * v1 + v2;
            r2 = v2 + m * v1;

            r3 = transpose(m) * v2 + v1;

            PetscVector expected = values(n, n * d1 * d2 + d3);

            utopia_test_assert(approxeq(expected, r1));
            utopia_test_assert(approxeq(expected, r2));

            expected = values(n, n * d1 * d3 + d2);
            utopia_test_assert(approxeq(expected, r3));
        }
    }

    void petsc_min()
    {
        const int n = mpi_world_size() * 2;
        PetscVector v  = values(n, 1.0);
        PetscMatrix A = identity(n, n);

        double min_v = min(v);
        utopia_test_assert(approxeq(1.0, min_v));

        double min_A = min(A);
        utopia_test_assert(approxeq(0.0, min_A));

        PetscVector min_row_A = min(A, 1);
        PetscVector expected  = values(n, 0.0);
        utopia_test_assert(approxeq(expected, min_row_A));
    }

    void petsc_max()
    {
        const int n = mpi_world_size() * 2;
        PetscVector v  = values(n, 1.0);
        PetscMatrix A = identity(n, n);

        double max_v = max(v);
        utopia_test_assert(approxeq(1.0, max_v));

        double max_A = max(A);
        utopia_test_assert(approxeq(1.0, max_A));

        PetscVector max_row_A = max(A, 1);
        PetscVector expected  = values(n, 1.0);
        utopia_test_assert(approxeq(expected, max_row_A));
    }

    void petsc_binary_min_max()
    {
        const int n = mpi_world_size() * 2;
        PetscVector one = values(n, 1.);
        PetscVector two = values(n, 2.);

        PetscVector actual_min = min(one, two);
        PetscVector actual_max = max(one, two);

        utopia_test_assert(approxeq(one, actual_min));
        utopia_test_assert(approxeq(two, actual_max));

        actual_min = min(two, values(n, 1.));
        actual_max = max(values(n, 2.), one);

        utopia_test_assert(approxeq(one, actual_min));
        utopia_test_assert(approxeq(two, actual_max));
    }

    void petsc_ghosted()
    {
        const int n = mpi_world_size() * 2;
        const int off = mpi_world_rank() * 2;

        std::vector<PetscInt> ghosts{ (off + 3) % n };
        PetscVector v = ghosted(2, n, ghosts);

        auto r = range(v);

        {
            Write<PetscVector> w_v(v, GLOBAL_INSERT);
            for(auto i = r.begin(); i != r.end(); ++i) {
                std::vector<PetscInt> index(1); index[0] = i;
                std::vector<PetscScalar> values(1); values[0] = static_cast<double>(i);
                v.set(index, values);
            }
        }

        {
            Read<PetscVector> r_v(v);
            std::vector<PetscInt> index{(off + 3) % n};
            std::vector<PetscScalar> values;
            v.get(index, values);
            utopia_test_assert(index[0] == PetscInt(values[0]));
        }

    }


    void petsc_block_mat()
    {
        const SizeType n = mpi_world_size() * 2;
        PetscMatrix mat = identity(n, n);
        mat.convert_to_mat_baij(2);
    }


    void petsc_line_search()
    {
        auto n = 10;
        PetscVector v = local_values(n, 1.);
        PetscMatrix m = local_identity(n, n);

        auto expr = dot(v, v)/dot(m * v, v);
        // std::cout << tree_format(expr.get_class()) << std::endl;

        double s = expr;
        utopia_test_assert(approxeq(1., s));
    }

    void petsc_residual()
    {
        auto n = 10;
        PetscVector x = local_values(n, 1.);
        PetscMatrix A = local_identity(n, n);
        PetscVector  b = local_values(n, 2.);

        PetscVector res = b - A * x;
        // disp(res);
        double s = sum(res);
        utopia_test_assert(approxeq(n * mpi_world_size(), s));
    }

    void petsc_transform()
    {
        PetscMatrix P;

        const std::string folder =  Utopia::instance().get("data_path") + "/laplace/matrices_for_petsc";
        bool ok = read(folder + "/I_2", P); utopia_test_assert(ok);


        double sum_P = sum(P);
        P = transpose(P);

        each_apply(P, [](const double value) -> double {
            return value * 2.;
        });

        double sum_P_2 = sum(P);

        utopia_test_assert(approxeq(sum_P * 2., sum_P_2));
    }

    void petsc_dot_test()
    {
        auto n = 10;

        PetscVector x1 = values(n, 1.);
        PetscVector x2 = values(n, 2.);
        PetscVector x3 = values(n, 3.);
        PetscVector x4 = values(n, 4.);

        std::vector<std::shared_ptr<PetscVector> > vectors_x;
        std::vector<PetscScalar> result_x;

        vectors_x.push_back(make_ref(x2));
        vectors_x.push_back(make_ref(x3));
        vectors_x.push_back(make_ref(x4));

        dots(x1, vectors_x, result_x);

        utopia_test_assert(approxeq(result_x[0], 2.*n));
        utopia_test_assert(approxeq(result_x[1], 3.*n));
        utopia_test_assert(approxeq(result_x[2], 4.*n));

        PetscScalar r1 = dot(x1, x2);
        PetscScalar r2 = dot(x1, x3);

        PetscScalar result_sum1 = r1 + r2;
        PetscScalar result_sum2 = dot(x1, x2) + dot(x1, x3);
        utopia_test_assert(approxeq(result_sum1, result_sum2));

        PetscScalar result_div1 = r1/r2;
        PetscScalar result_div2 = dot(x1, x2)/dot(x1, x3);
        utopia_test_assert(approxeq(result_div1, result_div2));

        PetscScalar result_mul1 = r1*r2;
        PetscScalar result_multv2 = dot(x1, x2)*dot(x1, x3);
        utopia_test_assert(approxeq(result_mul1, result_multv2));

        PetscScalar result_min1 = r1-r2;
        PetscScalar result_min2 = dot(x1, x2)-dot(x1, x3);
        utopia_test_assert(approxeq(result_min1, result_min2));


        PetscMatrix B = diag(PetscVector(values(n,1.0)));
        PetscScalar pred = -1.0 * dot(x1, x2) - 0.5 * dot(B * x3, x4);

        PetscScalar pred1 = dot(x1, x2);
        PetscScalar pred2 = dot(B * x3, x4);
        PetscScalar pred_sum = -1.0 * pred1 - 0.5 * pred2;
        utopia_test_assert(approxeq(pred, pred_sum));
    }


    void petsc_norm_test()
    {
        auto n = 10;

        PetscVector x1 = values(n, 1.);
        PetscVector x2 = values(n, 2.);
        PetscVector x3 = values(n, 3.);

        PetscScalar x1_norm_original = norm2(x1);
        PetscScalar x2_norm_original = norm2(x2);
        PetscScalar x3_norm_original = norm2(x3);

        PetscScalar r1, r2;
        PetscScalar r11, r12, r13;

        norms2(x1, x2, r1, r2);

        utopia_test_assert(approxeq(x1_norm_original, r1));
        utopia_test_assert(approxeq(x2_norm_original, r2));

        norms2(x1, x2, x3, r11, r12, r13);

        utopia_test_assert(approxeq(x1_norm_original, r11));
        utopia_test_assert(approxeq(x2_norm_original, r12));
        utopia_test_assert(approxeq(x3_norm_original, r13));

    }


    void petsc_get_col_test()
    {
        auto n = 10;
        auto m = 5;
        auto col_id = 2;

        PetscMatrix M = values(n, m, 0.0);
        {
            Write<PetscMatrix> w_m(M);
            auto r = row_range(M);
            auto c = col_range(M);

            for(auto i = r.begin(); i != r.end(); ++i)
            {
                for(auto j = c.begin(); j != c.end(); ++j)
                {
                    M.set(i, j, j);
                }
            }
        }

        PetscVector col_result = zeros(n);
        // mat_get_col(M, col_result, col_id);
        M.col(col_id, col_result);

        PetscVector col_expected = local_values(local_size(col_result).get(0), col_id);
        utopia_test_assert(approxeq(col_result, col_expected));

    }

    void petsc_memcheck()
    {
        UTOPIA_PETSC_MEMCHECK_BEGIN();
        UTOPIA_PETSC_MEMCHECK_END();
    }


    void petsc_dense_mat_mult_test()
    {
        if(mpi_world_size()>5)
            return;

        PetscMatrix A = values(5, 5, 2.0);
        PetscMatrix B = values(5, 5, 10.0);
        PetscMatrix C = A*B;

        utopia_test_assert(approxeq(norm_infty(C), 500));
    }

    void petsc_chop_test()
    {
        auto n = 10;

        PetscMatrix M = local_identity(n, n);
        {
            Write<PetscMatrix> w_m(M);
            auto r = row_range(M);
            // auto c = col_range(M);

            for(auto i = r.begin(); i != r.end(); ++i)
            {
                if(i < 5)
                {
                    M.set(i, i, 1.0);
                }
                else
                {
                    M.set(i,i, -1.0); 
                }
            }
        }

        PetscMatrix M_p = M; 
        PetscMatrix M_n = M;  

        chop_smaller_than(M_p, 1e-15); 
        chop_greater_than(M_n, 1e-15); 

        M_p += M_n; 
        utopia_test_assert(approxeq(M_p, M));        
    }

    void petsc_zero_rows_to_id()
    {   
        SizeType n = 4;

        PetscMatrix m = 2.*local_identity(n, n);
        m *= 0.;

        zero_rows_to_identity(m, 1e-10);
        // disp(m);
    }


    #endif //WITH_PETSC;

    void petsc_specific() {
#ifdef WITH_PETSC

        UTOPIA_RUN_TEST(petsc_memcheck);
        UTOPIA_RUN_TEST(petsc_line_search);
        UTOPIA_RUN_TEST(petsc_residual);
        UTOPIA_RUN_TEST(petsc_block_mat);
        UTOPIA_RUN_TEST(petsc_ghosted);

        // UTOPIA_RUN_TEST(petc_optional); // fails to compile with gpu

        // UTOPIA_RUN_TEST(petsc_view); //FIXME
        UTOPIA_RUN_TEST(petsc_ksp_precond_delegate);
        UTOPIA_RUN_TEST(petsc_hardcoded_cg);
        UTOPIA_RUN_TEST(petsc_reciprocal);
        UTOPIA_RUN_TEST(petsc_axpy);
        UTOPIA_RUN_TEST(petsc_vector_accessors);
        UTOPIA_RUN_TEST(petsc_matrix_accessors);
        UTOPIA_RUN_TEST(petsc_mv);
        UTOPIA_RUN_TEST(petsc_copy);
        UTOPIA_RUN_TEST(petsc_wrapper);
        UTOPIA_RUN_TEST(petsc_vector_composite);
        UTOPIA_RUN_TEST(petsc_matlab_connection);
        UTOPIA_RUN_TEST(petsc_matrix_composite);
        UTOPIA_RUN_TEST(petsc_mat_tests);
        UTOPIA_RUN_TEST(petsc_vec_tests);
        UTOPIA_RUN_TEST(petsc_read_and_write);
        UTOPIA_RUN_TEST(petsc_to_blas);
        UTOPIA_RUN_TEST(petsc_is_nan_or_inf);
        UTOPIA_RUN_TEST(petsc_mat_mul_add);
        UTOPIA_RUN_TEST(petsc_min);
        UTOPIA_RUN_TEST(petsc_max);
        UTOPIA_RUN_TEST(petsc_factory_and_operations);
        UTOPIA_RUN_TEST(petsc_each_sparse_matrix);
        UTOPIA_RUN_TEST(petsc_matrix_composition);
        UTOPIA_RUN_TEST(petsc_test_ptap);
        UTOPIA_RUN_TEST(petsc_test_rart);
        UTOPIA_RUN_TEST(petsc_new_eval);
        UTOPIA_RUN_TEST(petsc_tensor_reduction);
        UTOPIA_RUN_TEST(petsc_precond);
        UTOPIA_RUN_TEST(petsc_binary_min_max);
        UTOPIA_RUN_TEST(petsc_dot_test);
        UTOPIA_RUN_TEST(petsc_transform);
        UTOPIA_RUN_TEST(petsc_get_col_test);
        UTOPIA_RUN_TEST(petsc_dense_mat_mult_test);
        UTOPIA_RUN_TEST(petsc_norm_test);
        UTOPIA_RUN_TEST(petsc_chop_test); 
        UTOPIA_RUN_TEST(petsc_zero_rows_to_id);


        //serial tests
#ifdef PETSC_HAVE_MUMPS
        UTOPIA_RUN_TEST(petsc_inverse);
#endif //PETSC_HAVE_MUMPS

        // petsc_sparse_matrix_accessors();  // TODO:: here something doesnt work in parallel !

        // // petsc_local_entities(); //FIXME does not work
        //  petsc_conversion();
        //  maria();
        //  //local_diag_block();              // TODO:: utopia_test_assert fails in parallel

        #endif // WITH_PETSC
    }

    UTOPIA_REGISTER_TEST_FUNCTION(petsc_specific);
}
