#include "utopia.hpp"
#include "utopia_PetscTest.hpp"
#include "test_problems/utopia_TestFunctionsND.hpp"
#include "utopia_QuadraticFunction.hpp"
#include "utopia_Device.hpp"

namespace utopia {

#ifdef WITH_PETSC

    void petc_optional()
    {
        MPI_Comm sub_comm;
        MPI_Comm_split(
            PETSC_COMM_WORLD,
            mpi_world_rank() % 2 == 1,
            mpi_world_rank(),
            &sub_comm);

        int rank;
        MPI_Comm_rank(sub_comm, &rank);

        {
	        //optionals only work for this builder at the moment
	        DSMatrixd m = local_sparse(10, 10, 1, sub_comm, str("my_mat"));

	        auto r = row_range(m);
	        {
	            Write<DSMatrixd> w_m(m);
	            for(auto i = r.begin(); i != r.end(); ++i) {
	                m.set(i, i, rank);
	            }
	        }


	       DVectord vec;
	       std::string path = Utopia::instance().get("data_path");
	       read(path + "/RHS_10x10x10_hexa_3D", vec, sub_comm, str("my_vec"));
   		}

       MPI_Comm_free(&sub_comm);
    }


    void petsc_reciprocal() {
        // test also  diag
        DSMatrixd A = values(4, 4, 1.0);
        {
            Write<DSMatrixd> w(A);
            A.set(1, 1, 99);
            A.set(2, 2, 77);
        }

        DVectord diag_A = diag(A);
        DVectord v_expected = values(4, 1.0);
        {
            Write<DVectord> w(v_expected);
            v_expected.set(1, 99);
            v_expected.set(2, 77);
        }
        utopia_test_assert(approxeq(v_expected, diag_A));

        A = diag(diag_A);

        diag_A = values(4, 9.0);
        DVectord result_petsc = 5 / diag_A;
        v_expected = values(4, 5.0/9.0);
        utopia_test_assert(approxeq(v_expected, result_petsc));

        DVectord v_expected_test= values(4.0, 1.0);

        {
            Write<DVectord> w(v_expected_test);
            Range rr=range(v_expected_test);
            for (auto ii=rr.begin(); ii<rr.end(); ++ii){
                 v_expected_test.set(ii,2*ii);
        }

        }

        DVectord p = power(v_expected_test, 3.0);
    }

    void petsc_axpy() {

        {
            //! [axpy (petsc)]
            int n = 10;
            DVectord x = values(n, 1.0);
            DVectord y = values(n, 2.0);
            DVectord expected = values(n, 1.0 * 0.1 + 2.0);
            const double alpha = 0.1;

            DVectord result = alpha * x + y;
            utopia_test_assert(norm1(result - expected) < 1e-5);
            //! [axpy (petsc)]
        }

        ///////////////////////////////////

        {
            const PetscInt n = 10;
            const PetscInt m = 20;

            DMatrixd X = values(m, n, 1.0);
            DMatrixd Y = identity(m, n);
            const double alpha = 4;
            DMatrixd result = alpha * X + Y;
            double actual = sum(result);
            utopia_test_assert(approxeq(810., actual));
        }

        ///////////////////////////////////

        {
            const int m = mpi_world_size() * 3;
            const int n = mpi_world_size() * 2;

            DSMatrixd X = 2.  * identity(m, n);
            DSMatrixd Y = 0.1 * identity(m, n);
            double alpha = 4.;
            DSMatrixd res = alpha * X + Y;
        }

        {
            const int n = mpi_world_size() * 10;

            DSMatrixd expected = sparse(n, n, 2);
            DSMatrixd X = identity(n, n);
            DSMatrixd Y = sparse(n, n, 1);

            {
                Write<DSMatrixd> w_e(expected), w_Y(Y);
                Range r = row_range(expected);

                for(auto i = r.begin(); i < r.end(); ++i) {
                    expected.set(i, i, 0.1);
                    expected.set(i, n-i-1, 0.1);
                }

                for(auto i = r.begin(); i < r.end(); ++i) {
                    Y.set(i, n-i-1, 0.1);
                }
            }

            DSMatrixd actual = 0.1 * X + Y;

            utopia_test_assert(approxeq(expected, actual));
        }
    }

    void petsc_vector_accessors() {
        int mult;
        MPI_Comm_size(PETSC_COMM_WORLD, &mult);
        const PetscInt n = 10 * mult;
        DVectord x = values(n, 1, 1);
        DVectord y = values(n, 1, 2);

        //for parallel data structures (works also for serial ones)
        Range xr = range(x);
        const PetscInt xb = xr.begin();


        //The use of Read/Write/ReadAndWrite locks is important for universal compatibility with all (parallel) backends
        { //scoped write lock
            Write<DVectord> writing(x);
            //usual way
            x.set(xb, -1);
            x.set(xb + 9, -1);

            //the petsc way
            std::vector <PetscScalar> values{10., 10., 10};
            std::vector <PetscInt> index{xb + 2, xb + 3, xb + 4};
            x.set(index, values);
        }

        { //scoped read lock
            Read<DVectord> reading(x);
            PetscScalar v = x.get(xb + 1);
            utopia_test_assert(approxeq(1, v));
        }

        { //scoped read-write lock
            ReadAndWrite<DVectord> readingAndWriting(x);
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
        DMatrixd x = values(n, 2, 1);
        DMatrixd y = values(n, 2, 2);

        //for parallel data structures (works also for serial ones)
        Range xr = row_range(x);
        const PetscInt xb = xr.begin();

    //  Petsc matrix does not support ReadAndWrite
    //    //The use of Read/Write locks is important for universal compatibility with all (parallel) backends
        { //scoped write lock
            Write<DMatrixd> writing(x);
            x.set(xb, 0, -1);
            x.set(xb + 9, 0, -1);
        }

        { //scoped read lock
            Read<DMatrixd> reading(x);
            PetscScalar v = x.get(xb + 1, 0);
            utopia_test_assert(approxeq(1, v));
        }
    }

    void petsc_sparse_matrix_accessors() {

        //! [Read write matrix]
        int mult = mpi_world_size();
        const PetscInt n = 10 * mult;
        DSMatrixd x = sparse(n, 10, 7);
        DSMatrixd y = sparse(n, 10, 2);

        // For parallel data structures (works also for serial ones. Adopt paridgm for code portability)
        Range xr = row_range(x);
        const PetscInt xb = xr.begin();

        // Petsc matrix does not support ReadAndWrite only separate Read and Write
        // The use of Read/Write locks is important for universal compatibility with all (parallel) backends
        {
            Write<DSMatrixd> w(y);
            y.set(xb, 0, 0);
            y.set(xb + 1, 1, 1);
            y.set(xb + 2, 0, 4);
            y.set(xb + 3, 0, 6);
            y.set(xb + 4, 0, 8);
            y.set(xb + 9, 0, 18);
            y.set(xb + 9, 1, 19);
        }

        {
            Read<DSMatrixd> r(x);
            PetscScalar v = x.get(xb + 1, 1);
            utopia_test_assert(approxeq(3, v));
        }

        x = sparse(n, n, 1);
        {
            Write<DSMatrixd> w(x);

            x.set(xb, xb, 1);
            x.set(xb + 1, xb + 1, 2);
        }

        {
            Read<DSMatrixd> r(x);
            PetscScalar v = x.get(xb + 1, xb +1);
            utopia_test_assert(approxeq(2, v));
            v = x.get(xb, xb);
            utopia_test_assert(approxeq(1, v));
        }

        //! [Read write matrix]
    }

    void petsc_mv()
    {
        const PetscInt n = 10;
        const PetscInt m = 20;

        //creates a dense matrix
        DVectord v   = values(n, 1, 1);
        DMatrixd mat = identity(m, n);
        DVectord result = mat * v;

        DVectord expected = values(m, 1, 1);
        {
            Write<DVectord> w(expected);
            for (size_t i = n; i < m; i++) {
                expected.set(i, 0);
            }
        }

        utopia_test_assert(approxeq(expected, result));
    }

    void petsc_copy()
    {
        const PetscInt n = 10;
        const PetscInt m = 20;

        DVectord v1 = values(n, 1, 1);
        DVectord v2 = v1;
        utopia_test_assert(approxeq(v1, v2));

        DMatrixd m1 = identity(m, n);
        DMatrixd m2 = m1;
        utopia_test_assert(approxeq(m1, m2));
    }

    void petsc_wrapper()
    {
        DSMatrixd m = identity(2, 2);
        auto expected_ptr = raw_type(m);
        Mat pmat = raw_type(m);
        DSMatrixd wmat;
        wrap(pmat, wmat);
        utopia_test_assert( raw_type(wmat) == expected_ptr );
    }

    void petsc_vector_composite() {
        const PetscInt n = 10;
        DVectord v1 = values(n, 1, 1.0);
        DVectord v2 = values(n, 1, 2.0);
        DVectord v3 = values(n, 1, 0.2);

        auto expr = v1 * 0.1 + v2 + v1 - v3;
        DVectord vresult = expr;
        DVectord expected = values(n, 1, 2.9);
        utopia_test_assert(approxeq(expected, vresult));

        PetscScalar value = norm2(expr);
        utopia_test_assert(approxeq(std::sqrt(10*2.9*2.9), value));

        const double valueSum = sum(expr);
        utopia_test_assert(std::abs(valueSum - 29) < 1e-14);
    }

    void petsc_matlab_connection() {
        DVectord rhs, sol, M_rhs;
        DSMatrixd K, M;
        std::string path = Utopia::instance().get("data_path");

        read(path + "/RHS_10x10x10_hexa_3D", rhs);
        read(path + "/K_hexa_10x10x10_3D", K);
        read(path + "/M_hexa_10x10x10_3D", M);

        M_rhs = M * rhs;
        sol = local_zeros(local_size(M_rhs));

        // Linear solver
        // auto lsolver = std::make_shared<Factorization<DSMatrixd, DVectord>>();
        auto lsolver = std::make_shared<BiCGStab<DSMatrixd, DVectord>>();
        lsolver->solve(K, M_rhs, sol);

        double diff = norm2(M_rhs - K * sol);
        utopia_test_assert(diff < 1e-6);

        // Non-linear solver
        Newton<DSMatrixd, DVectord> newton(lsolver);
        newton.enable_differentiation_control(false);

        QuadraticFunction<DSMatrixd, DVectord> fun(make_ref(K), make_ref(M_rhs));
        sol = local_zeros(local_size(M_rhs));
        newton.solve(fun, sol);

        diff = norm2(M_rhs - K * sol);
        utopia_test_assert(diff < 1e-6);
    }


    void petsc_matrix_composite() {
        // std::cout << "begin: petsc_matrix_composite_test" << std::endl;

        const PetscInt n = 10;
        DMatrixd m1 = values(n, n, 1.0);
        DMatrixd m2 = values(n, n, 2.0);
        DMatrixd m3 = values(n, n, 0.2);

        DMatrixd result = m1 * 0.1 + m2 + m1 - m3;
        DMatrixd expected = values(n, n, 2.9);
        utopia_test_assert(approxeq(expected, result));

        // std::cout << "end: petsc_matrix_composite_test" << std::endl;
    }

    void petsc_view()
    {
        //! [Global views]
        const PetscInt offset = mpi_world_size();
        const PetscInt n = offset + 5;
        const PetscInt m = offset + 10;
        DSMatrixd m1 = identity(m, n);
        DSMatrixd m2 = m1.range(
            1, offset + 6,
            0, offset + 5);


        each_read(m2, [](const SizeType i, const SizeType j, const double value) {
            if (i + 1 == j) {
                if(!approxeq(1, value)) {
                    std::cout << i << ", " << j << " -> " << value << std::endl;
                }
                utopia_test_assert(approxeq(1, value));
            } else {
                utopia_test_assert(approxeq(0, value));
            }
        });

        //! [Global views]
    }

    void petsc_mat_tests()
    {
        PetscBool assembled;
        DSMatrixd M = zeros(5, 5);
        MatAssembled(raw_type(M), &assembled);
        utopia_test_assert(!empty(M));

        DSMatrixd A;
        utopia_test_assert(empty(A));
    }

    void petsc_vec_tests()
    {
        DVectord v = zeros(5);
        utopia_test_assert(!empty(v));

        DVectord e;
        utopia_test_assert(empty(e));
    }

    void petsc_read_and_write()
    {
        //! [Input and output (petsc)]

        DVectord x = values(3, 1, 3.0);

        // Write vector to disk
        write("test_vec.txt", x);

        // Display values on the console
        // disp(x);

        DVectord y;
        // Read vector from disk
        read("test_vec.txt", y);

        utopia_test_assert(approxeq(x, y));

        DSMatrixd m = sparse(3, 3, 2);
        {
            Write<DSMatrixd> w(m);
            m.set(0, 0, 1);
            m.set(0, 1, 2);
            m.set(1, 1, 3);
            m.set(2, 2, 4);
        }

        DSMatrixd w;
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

    DVectord x = zeros(16);

    Range xr = range(x);
    const PetscInt xb = xr.begin();

    Vectord y = values(xr.extent(), 2.0);

    {
        Read<Vectord> r_y(y);
        Write<DVectord> w_x(x);
        for (SizeType i = 0; i < xr.extent(); ++i)
            x.set(xb + i, y.get(i));
    }

    DVectord expected = values(16, 2.0);
    utopia_test_assert(approxeq(expected, x));

#endif //WITH_BLAS
    }



    //FIXME does not work
    void petsc_local_entities() {
        std::cout << "Begin: petsc_local_entities_test." << std::endl;
        DSMatrixd matrix;

        PetscInt rank, size;
        MPI_Comm_rank(matrix.implementation().communicator(), &rank);
        MPI_Comm_size(matrix.implementation().communicator(), &size);
        DVectord ones = local_values(2, rank);

        matrix = local_identity((rank && rank < size - 1) ? 4 : 3, 2);

        // disp(matrix);

        if (size > 1) {
            Write<DSMatrixd> write(matrix);
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

        DVectord result = matrix * ones;

        // disp(matrix);
        // disp(result);
        std::cout << "End: petsc_local_entities_test." << std::endl;
    }



    void petsc_conversion() {
        // std::cout << "Begin: petsc_conversion_test." << std::endl;
        DVectord vec = values(10, 1.0);


        DVectord dvec;
        convert(raw_type(vec), dvec);
        utopia_test_assert(approxeq(vec, dvec));
        convert(dvec, raw_type(vec));
        utopia_test_assert(approxeq(vec, dvec));


        DMatrixd mat_dense = values(2, 3, 1);
        DMatrixd dmat_dense;
        convert(raw_type(mat_dense), dmat_dense);
        utopia_test_assert(approxeq(mat_dense, dmat_dense));
        convert(dmat_dense, raw_type(mat_dense));

        DSMatrixd mat = sparse(3, 3, 2);
        {
            Write<DSMatrixd> w(mat);
            mat.set(0, 0, 1);
            mat.set(0, 1, 2);
            mat.set(1, 1, 3);
            mat.set(2, 2, 4);
        }

        DSMatrixd dmat;
        convert(raw_type(mat), dmat);
        convert(dmat, raw_type(mat));

        // std::cout << "End: petsc_conversion_test." << std::endl;
    }


    void petsc_factory_and_operations()
    {
        const int n = mpi_world_size() * 3;
        DSMatrixd m = sparse(n, n, 3);
        {
            auto r = row_range(m);

            Write<DSMatrixd> w(m);
            for(auto i = r.begin(); i < r.end(); ++i) {
                if(i == 0 || i == n - 1) {
                    m.set(1, 1, 1.0);
                    continue;
                }

                m.set(i, i, 2.);
                m.set(i, i-1, -1.);
                m.set(i, i+1, -1.);
            }
        }

        DVectord c = (m + 0.1 * identity(n, n)) * values(n, 0.5);
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

        DSMatrixd B = local_sparse(slave_sizes[r], master_sizes[r], 2);
        DSMatrixd D = local_sparse(slave_sizes[r], slave_sizes[r],  2);

        PetscViewer fd;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, (data_path + "/master_matrix").c_str(), FILE_MODE_READ, &fd);
        MatLoad(raw_type(B),fd);
        PetscViewerDestroy(&fd);

        PetscViewerBinaryOpen(PETSC_COMM_WORLD, (data_path + "/slave_matrix").c_str(), FILE_MODE_READ, &fd);
        MatLoad(raw_type(D), fd);
        PetscViewerDestroy(&fd);
        DSMatrixd Dinv = diag(1./sum(D, 1));
        DSMatrixd T = Dinv * B;

        DSMatrixd Tt = transpose(Dinv * B);

        // disp(Tt.size());
        // disp(T.size());

        DSMatrixd m = transpose(Dinv + D * Dinv);
    }

    void local_diag_block()
    {
        DSMatrixd a = sparse(4, 4, 3);

        {
            Write<DSMatrixd> write(a);
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

        SSMatrixd b = local_diag_block(a);
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
        DSMatrixd a = sparse(n, n, 3);

        {
            auto r = row_range(a);
            Write<DSMatrixd> write(a);

            for(auto i = r.begin(); i < r.end(); ++i) {
                a.set(i, i, i);
            }
        }

        each_read(a, [](const SizeType i, const SizeType j, const double value) {
            utopia_test_assert(approxeq(i, value));
        });
    }

    void petsc_test_ptap()
    {
        const int n = mpi_world_size() * 3;
        DSMatrixd P = identity(n, n);
        DSMatrixd A = identity(n, n);
        DSMatrixd C = identity(n, n);
        DSMatrixd PtAP, ABC, PAPt;

        //The next line is equivalent to this:  PtAP = ptap(P, A); since it is pattern matched
        PtAP = transpose(P) * A * P;
        ABC = A * P * C;

        DSMatrixd expected = identity(n, n);
        utopia_test_assert(approxeq(expected, PtAP));
        utopia_test_assert(approxeq(expected, ABC));
    }


    void petsc_matrix_composition()
    {
        const SizeType n = mpi_world_size() * 2;
        DSMatrixd m = identity(n, n);
        DVectord v = values(n, 2.);
        auto expr = abs(transpose(0.1 * (m * m) - m) * (m) * v);

        DVectord res = expr;
        DVectord expected = values(n, 0.9 * 2.);
        utopia_test_assert(approxeq(expected, res));
    }


    template<class Tensor, int Order>
    utopia::Wrapper<Tensor, Order> &to_wrapper(utopia::Wrapper<Tensor, Order> &t)
    {
        return t;
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
        DVectord x = values(n, 1.0);
        DVectord y = values(n, 1.0);

        eval(construct(y, alpha * x + y));
        DMatrixd m;
        eval(construct(m, outer(x, y)));
        // disp(m);

    //    DVectord v2 = values(11, 0.0);
    //    auto view = v2.range(0, 10);
    //    eval(assign(view, to_wrapper(x)));
    //    disp(v2);

    }

    void petsc_precond()
    {
        int n = mpi_world_size() * 10;
        DVectord d   = values(n, 2.0);
        DSMatrixd dm = diag(d);
        DVectord rhs = values(n, 1.);
        DVectord sol = zeros(n);

        DVectord v = diag(d) * rhs;

        auto m = std::make_shared<DMatrixd>(2.*identity(n, n));

        auto precond = std::make_shared< InvDiagPreconditioner<DMatrixd, DVectord> >();
        precond->update(m);

        auto cg = std::make_shared< ConjugateGradient<DMatrixd, DVectord> >();
        cg->set_preconditioner(precond);
        cg->solve(*m, rhs, sol);

        DVectord expected = values(n, 0.5);
        utopia_test_assert(approxeq(expected, sol));


        // PetscViewer   viewer;
        // PetscViewerDrawOpen(PETSC_COMM_WORLD,0,"",300,0,300,300,&viewer);
        // VecView(raw_type(sol),viewer);
        // PetscViewerDestroy(&viewer);

    }


    void petsc_tensor_reduction()
    {
        int n = mpi_world_size() * 10;
        DSMatrixd mat = identity(n, n);

        //summing columns
        DVectord v = sum(mat, 1);
        DVectord expected = values(n, 1);
        utopia_test_assert(approxeq(expected, v));
    }

    void petsc_inverse()
    {
        if(mpi_world_size() == 1) {
            DMatrixd mat = identity(3, 3);

            {
                Write<DMatrixd> w(mat);
                mat.set(0, 1, 2);
                mat.set(1, 0, 2);
            }

            DMatrixd inv_mat  = inv(mat);
            DMatrixd actual   = inv_mat * mat;
            DMatrixd expected = identity(3, 3);

            utopia_test_assert( approxeq(actual, expected) );
        }
    }

    void petsc_hardcoded_cg()
    {
        const int n = mpi_world_size() * 2;
        const int i_max = 2;
        const double eps = 1e-8;

        DVectord x  = zeros(n);
        DVectord b  = values(n, 1.0);
        DSMatrixd A = identity(n, n);
        DSMatrixd M = diag(diag(A));

        /////////////////////
        DVectord r = zeros(n), d = zeros(n), q = zeros(n), s = zeros(n);
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

        TestFunctionND_1<DMatrixd, DVectord> fun(n);

        ConjugateGradient<DMatrixd, DVectord> cg;
        cg.set_preconditioner(std::make_shared<DelegatePreconditioner<DMatrixd, DVectord> >());
        Newton<DMatrixd, DVectord> newton(make_ref(cg));

        DVectord x = zeros(n);
        newton.solve(fun, x);
        // disp(x);

        DVectord expected = values(n, 0.468919);
        utopia_test_assert(approxeq(expected, x));

        cg.set_preconditioner(std::make_shared<PointJacobi<DMatrixd, DVectord> >());
        x = zeros(n);

        newton.solve(fun, x);
        utopia_test_assert(approxeq(expected, x));
    }

    void petsc_is_nan_or_inf()
    {
        const int n     = 10;
        if(mpi_world_size() > n) return;

        DVectord denom  = zeros(n);
        DVectord nom    = values(n, 1.0);

        DVectord sol    = nom/denom;

        {
            Write<DVectord> w(sol);
            sol.set(2, 1);
        }

        //FIXME
        if(!has_nan_or_inf(sol)) {
            utopia_error("petsc_is_nan_or_inf: failed. Known problem on alpine-linux");
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

            DSMatrixd m = values(n, n, d1);
            DVectord v1 = values(n, d2);
            DVectord v2 = values(n, d3);
            DVectord r1, r2, r3;

            r1 = m * v1 + v2;
            r2 = v2 + m * v1;

            r3 = transpose(m) * v2 + v1;

            DVectord expected = values(n, n * d1 * d2 + d3);

            utopia_test_assert(approxeq(expected, r1));
            utopia_test_assert(approxeq(expected, r2));

            expected = values(n, n * d1 * d3 + d2);
            utopia_test_assert(approxeq(expected, r3));
        }
    }

    void petsc_min()
    {
        const int n = mpi_world_size() * 2;
        DVectord v  = values(n, 1.0);
        DSMatrixd A = identity(n, n);

        double min_v = min(v);
        utopia_test_assert(approxeq(1.0, min_v));

        double min_A = min(A);
        utopia_test_assert(approxeq(0.0, min_A));

        DVectord min_row_A = min(A, 1);
        DVectord expected  = values(n, 0.0);
        utopia_test_assert(approxeq(expected, min_row_A));
    }

    void petsc_max()
    {
        const int n = mpi_world_size() * 2;
        DVectord v  = values(n, 1.0);
        DSMatrixd A = identity(n, n);

        double max_v = max(v);
        utopia_test_assert(approxeq(1.0, max_v));

        double max_A = max(A);
        utopia_test_assert(approxeq(1.0, max_A));

        DVectord max_row_A = max(A, 1);
        DVectord expected  = values(n, 1.0);
        utopia_test_assert(approxeq(expected, max_row_A));
    }

    void petsc_binary_min_max()
    {
        const int n = mpi_world_size() * 2;
        DVectord one = values(n, 1.);
        DVectord two = values(n, 2.);

        DVectord actual_min = min(one, two);
        DVectord actual_max = max(one, two);

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
        DVectord v = ghosted(2, n, ghosts);

        auto r = range(v);

        {
            Write<DVectord> w_v(v);
            for(auto i = r.begin(); i != r.end(); ++i) {
                v.set(i, i);
            }
        }

        {
            Read<DVectord> r_v(v);
            std::vector<PetscInt> index{(off + 3) % n};
            std::vector<PetscScalar> values;
            v.get(index, values);
            utopia_test_assert(index[0] == PetscInt(values[0]));
        }

    }


    void petsc_block_mat()
    {
        const SizeType n = mpi_world_size() * 2;
        DSMatrixd mat = sparse(n, n, 3);
        {
            Write<DSMatrixd> w_m(mat);
            mat.set(0, 0, 1.);
            mat.set(1, 1, 1.);
        }

        mat.implementation().convert_to_mat_baij(2);
    }


    void petsc_line_search()
    {
        auto n = 10;
        DVectord v = local_values(n, 1.);
        DSMatrixd m = local_identity(n, n);

        auto expr = dot(v, v)/dot(m * v, v);
        // std::cout << tree_format(expr.getClass()) << std::endl;

        double s = expr;
        utopia_test_assert(approxeq(1., s));
    }

    void petsc_residual()
    {
        auto n = 10;
        DVectord x = local_values(n, 1.);
        DSMatrixd A = local_identity(n, n);
        DVectord  b = local_values(n, 2.);

        DVectord res = b - A * x;
        // disp(res);
        double s = sum(res);
        utopia_test_assert(approxeq(n * mpi_world_size(), s));
    }

    void petsc_transform()
    {
        DSMatrixd P;

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

        DVectord x1 = values(n, 1.);
        DVectord x2 = values(n, 2.);
        DVectord x3 = values(n, 3.);
        DVectord x4 = values(n, 4.);

        std::vector<std::shared_ptr<DVectord> > vectors_x; 
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


        DSMatrixd B = diag(DVectord(values(n,1.0)));
        PetscScalar pred = -1.0 * dot(x1, x2) - 0.5 * dot(B * x3, x4);

        PetscScalar pred1 = dot(x1, x2); 
        PetscScalar pred2 = dot(B * x3, x4);
        PetscScalar pred_sum = -1.0 * pred1 - 0.5 * pred2; 
        utopia_test_assert(approxeq(pred, pred_sum));  
    }

    void petsc_get_col_test()
    {
        auto n = 10; 
        auto m = 5; 
        auto col_id = 2; 

        DMatrixd M = values(n, m, 0.0); 
        {
            Write<DMatrixd> w_m(M);
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

        DVectord col_result = zeros(n); 
        mat_get_col(M, col_result, col_id); 
        
        DVectord col_expected = local_values(local_size(col_result).get(0), col_id); 
        utopia_test_assert(approxeq(col_result, col_expected));  

    }


    #endif //WITH_PETSC;

    void runPetscTest() {
#ifdef WITH_PETSC

        UTOPIA_UNIT_TEST_BEGIN("PetscTest");
        UTOPIA_RUN_TEST(petsc_line_search);
        UTOPIA_RUN_TEST(petsc_residual);
        UTOPIA_RUN_TEST(petsc_block_mat);
        UTOPIA_RUN_TEST(petsc_ghosted);

        // UTOPIA_RUN_TEST(petc_optional); // fails to compile with gpu
        
        UTOPIA_RUN_TEST(petsc_view);
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
        UTOPIA_RUN_TEST(petsc_new_eval);
        UTOPIA_RUN_TEST(petsc_tensor_reduction);
        UTOPIA_RUN_TEST(petsc_precond);
        UTOPIA_RUN_TEST(petsc_binary_min_max);
        UTOPIA_RUN_TEST(petsc_dot_test); 
        UTOPIA_RUN_TEST(petsc_transform);
        UTOPIA_RUN_TEST(petsc_get_col_test); 

        //serial tests
#ifdef PETSC_HAVE_MUMPS
        UTOPIA_RUN_TEST(petsc_inverse);
#endif //PETSC_HAVE_MUMPS

        // petsc_sparse_matrix_accessors();  // TODO:: here something doesnt work in parallel !

        // // petsc_local_entities(); //FIXME does not work
        //  petsc_conversion();
        //  maria();
        //  //local_diag_block();              // TODO:: utopia_test_assert fails in parallel

        UTOPIA_UNIT_TEST_END("PetscTest");
        #endif // WITH_PETSC
    }
}
