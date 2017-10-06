/*
* @Author: Eric Botter
* @Date:   2016-11-15
*/
#include "utopia.hpp"
#include "utopia_PETScTest.hpp"
#include "test_problems/utopia_TestFunctionsND.hpp"

namespace utopia {

#ifdef WITH_PETSC

    void petsc_reciprocal_test() {
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
        assert(approxeq(v_expected, diag_A));

        A = diag(diag_A);

        diag_A = values(4, 9.0);
        DVectord result_petsc = 5 / diag_A;
        v_expected = values(4, 5.0/9.0);
        assert(approxeq(v_expected, result_petsc));
    }

    void petsc_axpy_test() {
       
        {
            //! [axpy (petsc)]
            int n = 10;
            DVectord x = values(n, 1.0);
            DVectord y = values(n, 2.0);
            DVectord expected = values(n, 1.0 * 0.1 + 2.0);
            const double alpha = 0.1;

            DVectord result = alpha * x + y;
            assert(norm1(result - expected) < 1e-5);
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
            assert(approxeq(810., actual));
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
    }

    void petsc_vector_accessors_test() {
        // std::cout << "begin: petsc_vector_accessors_test" << std::endl;

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
            assert(approxeq(1, v));
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
                assert(approxeq(-1, value));
                break;
                case 2:
                case 3:
                case 4:
                assert(approxeq(10, value));
                break;
                case 5:
                case 6:
                case 7:
                case 8:
                assert(approxeq(1, value));
            }
        });

        // std::cout << "end: petsc_vector_accessors_test" << std::endl;
    }


    void petsc_matrix_accessors_test() {
        // std::cout << "begin: petsc_matrix_accessors_test" << std::endl;

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
            //usual way
            x.set(xb, 0, -1);
            x.set(xb + 9, 0, -1);
    //
    //        //the petsc way
            std::vector <PetscScalar> values{10., 10., 10};
            std::vector <PetscInt> rowIndex{xb + 2, xb + 3, xb + 4};
            std::vector <PetscInt> colIndex{0, 0, 0};
            x.set(rowIndex, colIndex, values);
        }

        { //scoped read lock
            Read<DMatrixd> reading(x);
            PetscScalar v = x.get(xb + 1, 0);
            assert(approxeq(1, v));
        }

        // std::cout << "end: petsc_matrix_accessors_test" << std::endl;
    }

    void petsc_sparse_matrix_accessors_test() {
        // std::cout << "begin: petsc_sparse_matrix_accessors_test" << std::endl;

        //! [Read write matrix]
        int mult = mpi_world_size();
        const PetscInt n = 10 * mult;
        DSMatrixd x = sparse(n, 10, 2);
        DSMatrixd y = sparse(n, 10, 2);

        // For parallel data structures (works also for serial ones. Adopt paridgm for code portability)
        Range xr = row_range(x);
        const PetscInt xb = xr.begin();

        // Petsc matrix does not support ReadAndWrite only separate Read and Write
        // The use of Read/Write locks is important for universal compatibility with all (parallel) backends
        {
            Write<DSMatrixd> w(x);
            x.set({xb, xb + 1, xb + 2, xb + 3, xb + 4, xb + 9, xb + 9},
              {0, 1, 0, 0, 0, 0, 1},
              {0, 3, 4, 6, 8, 18, 19});
        }

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
            assert(approxeq(3, v));
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
            assert(approxeq(2, v));
            v = x.get(xb, xb);
            assert(approxeq(1, v));
        }

        //! [Read write matrix]

        // std::cout << "end: petsc_sparse_matrix_accessors_test" << std::endl;
    }

    void petsc_mv_test()
    {
        // std::cout << "begin: petsc_mv_test" << std::endl;

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
        assert(approxeq(expected, result));

        // std::cout << "end: petsc_mv_test" << std::endl;
    }

    void petsc_copy_test()
    {
        // std::cout << "begin: petsc_copy_test" << std::endl;

        const PetscInt n = 10;
        const PetscInt m = 20;

        DVectord v1 = values(n, 1, 1);
        DVectord v2 = v1;
        assert(approxeq(v1, v2));

        DMatrixd m1 = identity(m, n);
        DMatrixd m2 = m1;
        assert(approxeq(m1, m2));

        // std::cout << "end: petsc_copy_test" << std::endl;
    }


    void petsc_wrapper_test()
    {
        // std::cout << "begin: petsc_wrapper_test" << std::endl;

        DSMatrixd m = identity(2, 2);
        auto expected_ptr = raw_type(m);
        Mat pmat = raw_type(m);
        DSMatrixd wmat = sparse_mref(pmat);
        assert( raw_type(wmat) == expected_ptr );

        // std::cout << "end: petsc_wrapper_test" << std::endl;
    }


    void petsc_vector_composite_test() {
        // std::cout << "begin: petsc_vector_composite_test" << std::endl;

        const PetscInt n = 10;
        DVectord v1 = values(n, 1, 1.0);
        DVectord v2 = values(n, 1, 2.0);
        DVectord v3 = values(n, 1, 0.2);

        auto expr = v1 * 0.1 + v2 + v1 - v3;
        DVectord vresult = expr;
        DVectord expected = values(n, 1, 2.9);
        assert(approxeq(expected, vresult));

        PetscScalar value = norm2(expr);
        assert(approxeq(std::sqrt(10*2.9*2.9), value));

        const double valueSum = sum(expr);
        assert(std::abs(valueSum - 29) < 1e-14);

        // std::cout << "end: petsc_vector_composite_test" << std::endl;
    }

    void petsc_matlab_connection_test() {
        // test of hexahedrons 3D FEM - made by matlab assembly
        DVectord rhs, sol;
        DSMatrixd K, M;

        // with 2 double-dots works
        std::string path = Utopia::Instance().get("data_path");

        read(path + "/RHS_10x10x10_hexa_3D", rhs);
        read(path + "/K_hexa_10x10x10_3D", K);
        read(path + "/M_hexa_10x10x10_3D", M);

        // test of direct solver
        auto lsolver = std::make_shared<BiCGStab<DSMatrixd, DVectord>>();
        // lsolver->solve(K, rhs, sol );
        // std::cout<< "LS - done \n";


        // test newton
        Newton<DSMatrixd, DVectord> newton(lsolver);
        newton.enable_differentiation_control(false);
        //newton.enableVerbose(false);

        // since case is lin. - QP minimization
        // std::cout << "Running FEM 3D poisson example using Newton solver \n";
        QuadraticFunction<DSMatrixd, DVectord> funn(rhs, K);
        newton.solve(funn, rhs);

        // write sol to the file
        // write(path+ "SOLUTION_hexa_3D", sol);
        // std::cout<< "Solution written into file ... \n";

    }


    void petsc_matrix_composite_test() {
        // std::cout << "begin: petsc_matrix_composite_test" << std::endl;

        const PetscInt n = 10;
        DMatrixd m1 = values(n, n, 1.0);
        DMatrixd m2 = values(n, n, 2.0);
        DMatrixd m3 = values(n, n, 0.2);

        DMatrixd result = m1 * 0.1 + m2 + m1 - m3;
        DMatrixd expected = values(n, n, 2.9);
        assert(approxeq(expected, result));

        // std::cout << "end: petsc_matrix_composite_test" << std::endl;
    }

    void petsc_view_test()
    {
        // std::cout << "begin: petsc_view_test" << std::endl;

        //! [Global views]
        const PetscInt n = 5;
        const PetscInt m = 10;
        DMatrixd m1 = identity(m, n);
        DMatrixd m2 = m1.range(
            1, 6, 
            0, 5);

        disp(m2);

        each_read(m2, [](const SizeType i, const SizeType j, const double value) {
            if (i + 1 == j) {
                if(!approxeq(1, value)) {
                    std::cout << i << ", " << j << " -> " << value << std::endl;
                }
                assert(approxeq(1, value));
            } else {
                assert(approxeq(0, value));
            }
        });

        //! [Global views]
        // std::cout << "end: petsc_view_test" << std::endl;
    }


    void petsc_mat_tests()
    {
        PetscBool assembled;
        DSMatrixd M = zeros(5,5);
        MatAssembled(raw_type(M), &assembled);
        assert(!empty(M));

        DSMatrixd A;
        assert(empty(A));
    }



    void petsc_read_and_write_test()
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

        assert(approxeq(x, y));

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

        assert(approxeq(m, w));

        //! [Input and output (petsc)]
    }



    void petsc_to_blas_test()
    {
    #ifdef WITH_BLAS
        DVectord x = zeros(16);

        Range xr = range(x);
        const PetscInt xb = xr.begin();

        Vectord y = values(xr.extent(), 2.0);

        {
            Read<Vectord> r_y(y);
            Write<DVectord> w_x(x);
            for (SizeType i = xb; i < xb + xr.extent(); ++i)
                x.set(i, y.get(i - xb));
        }

        DVectord expected = values(16, 2.0);
        assert(approxeq(expected, x));
    #endif //WITH_BLAS
    }



    //FIXME does not work
    void petsc_local_entities_test() {
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



    void petsc_conversion_test() {
        // std::cout << "Begin: petsc_conversion_test." << std::endl;
        DVectord vec = values(10, 1.0);


        DVectord dvec;
        convert(raw_type(vec), dvec);
        assert(approxeq(vec, dvec));
        convert(dvec, raw_type(vec));
        assert(approxeq(vec, dvec));


        DMatrixd mat_dense = values(2, 3, 1);
        DMatrixd dmat_dense;
        convert(raw_type(mat_dense), dmat_dense);
        assert(approxeq(mat_dense, dmat_dense));
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


    void petsc_factory_and_operations_test()
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

    void maria_test()
    {
        if(mpi_world_size() != 2) return;

        int master_sizes[2] = {
            15, 10
        };

        int slave_sizes[2] = {
            25, 20
        };

        const std::string data_path = Utopia::Instance().get("data_path") + "/master_and_slave";

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

    void local_diag_block_test()
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
                assert(approxeq(i + 1, value));
            else
                assert(approxeq(4, value));
        });

        SSMatrixd b = local_diag_block(a);
        each_read(b, [](const SizeType i, const SizeType j, const double value) {
            if (i == j)
                assert(approxeq(i + 1, value));
            else
                assert(approxeq(4, value));
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
            assert(approxeq(i, value));
        });
    }

    void petsc_test_mat_ptap_product()
    {
        const int n = mpi_world_size() * 3;
        DSMatrixd P = identity(n, n);
        DSMatrixd A = identity(n, n);
        DSMatrixd C = identity(n, n);
        DSMatrixd PtAP, ABC, PAPt;

        //The next line is equivalent to this:  PtAP = mat_PtAP_product(P, A); since it is pattern matched
        PtAP = transpose(P) * A * P;
        ABC = A * P * C;

        DSMatrixd expected = identity(n, n);
        assert(approxeq(expected, PtAP));
        assert(approxeq(expected, ABC));
    }


    void petsc_matrix_composition_test()
    {
        const SizeType n = mpi_world_size() * 2;
        DSMatrixd m = identity(n, n);
        DVectord v = values(n, 2.);
        auto expr = abs(transpose(0.1 * (m * m) - m) * (m) * v);
        
        DVectord res = expr;
        DVectord expected = values(n, 0.9 * 2.);
        assert(approxeq(expected, res));
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

    void petsc_new_eval_test()
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

    void petsc_precond_test()
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
        assert(approxeq(expected, sol));


        // PetscViewer   viewer;
        // PetscViewerDrawOpen(PETSC_COMM_WORLD,0,"",300,0,300,300,&viewer);
        // VecView(raw_type(sol),viewer);
        // PetscViewerDestroy(&viewer);

    }


    void petsc_tensor_reduction_test()
    {
        int n = mpi_world_size() * 10;
        DSMatrixd mat = identity(n, n);

        //summing columns
        DVectord v = sum(mat, 1);
        DVectord expected = values(n, 1);
        assert(approxeq(expected, v));
    }

    void petsc_inverse_test()
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

            assert( approxeq(actual, expected) );
        }   
    }

    void petsc_harcoded_cg_test()
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

    void petsc_ksp_precond_delegate_test()
    {
        const int n = 10;
        TestFunctionND_1<DMatrixd, DVectord> fun(n);

        ConjugateGradient<DMatrixd, DVectord> cg;
        cg.set_preconditioner(std::make_shared<DelegatePreconditioner<DMatrixd, DVectord> >());
        Newton<DMatrixd, DVectord> newton(make_ref(cg));

        DVectord x = zeros(n);
        newton.solve(fun, x);
        // disp(x);

        DVectord expected = values(n, 0.468919);
        assert(approxeq(expected, x));

        cg.set_preconditioner(std::make_shared<PointJacobi<DMatrixd, DVectord> >());
        x = zeros(n);

        newton.solve(fun, x);
        assert(approxeq(expected, x));
    }

    void petsc_is_nan_or_inf_test()
    {
        const int n     = 10; 
        DVectord denom  = zeros(n);
        DVectord nom    = values(n, 1.0);

        DVectord sol    = nom/denom; 
        
        {
            Write<DVectord> w(sol);
            sol.set(2, 1); 
        }

        assert(has_nan_or_inf(sol) == 1);
        assert(has_nan_or_inf(denom) == 0);

    }


    void petsc_mat_mul_add_test()
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

            assert(approxeq(expected, r1));
            assert(approxeq(expected, r2));

            expected = values(n, n * d1 * d3 + d2);
            assert(approxeq(expected, r3));
        }
    }


    void petsc_min_test()
    {
        const int n = mpi_world_size() * 2;
        DVectord v  = values(n, 1.0);
        DSMatrixd A = identity(n, n);

        double min_v = min(v);
        assert(approxeq(1.0, min_v));

        double min_A = min(A);
        assert(approxeq(0.0, min_A));

        DVectord min_row_A = min(A, 1);
        DVectord expected  = values(n, 0.0);
        assert(approxeq(expected, min_row_A));
    }

    void petsc_max_test()
    {
        const int n = mpi_world_size() * 2;
        DVectord v  = values(n, 1.0);
        DSMatrixd A = identity(n, n);

        double max_v = max(v);
        assert(approxeq(1.0, max_v));

        double max_A = max(A);
        assert(approxeq(1.0, max_A));

        DVectord max_row_A = max(A, 1);
        DVectord expected  = values(n, 1.0);
        assert(approxeq(expected, max_row_A));
    }

    #endif //WITH_PETSC;

    void runPETScTest() {
#ifdef WITH_PETSC
        UTOPIA_UNIT_TEST_BEGIN("PETScTest");
        
        UTOPIA_RUN_TEST(petsc_ksp_precond_delegate_test);
        UTOPIA_RUN_TEST(petsc_harcoded_cg_test);
        UTOPIA_RUN_TEST(petsc_reciprocal_test);
        UTOPIA_RUN_TEST(petsc_axpy_test);
        UTOPIA_RUN_TEST(petsc_vector_accessors_test);
        UTOPIA_RUN_TEST(petsc_matrix_accessors_test);
        UTOPIA_RUN_TEST(petsc_mv_test);
        UTOPIA_RUN_TEST(petsc_copy_test);
        UTOPIA_RUN_TEST(petsc_wrapper_test);
        UTOPIA_RUN_TEST(petsc_vector_composite_test);
        UTOPIA_RUN_TEST(petsc_matlab_connection_test);
        UTOPIA_RUN_TEST(petsc_matrix_composite_test);
        UTOPIA_RUN_TEST(petsc_mat_tests);
        UTOPIA_RUN_TEST(petsc_read_and_write_test);
        UTOPIA_RUN_TEST(petsc_to_blas_test);
        UTOPIA_RUN_TEST(petsc_is_nan_or_inf_test); 
        UTOPIA_RUN_TEST(petsc_mat_mul_add_test);
        UTOPIA_RUN_TEST(petsc_min_test);
        UTOPIA_RUN_TEST(petsc_max_test);
        UTOPIA_RUN_TEST(petsc_factory_and_operations_test);
        UTOPIA_RUN_TEST(petsc_each_sparse_matrix);
        UTOPIA_RUN_TEST(petsc_matrix_composition_test);
        UTOPIA_RUN_TEST(petsc_test_mat_ptap_product);
        UTOPIA_RUN_TEST(petsc_new_eval_test);
        UTOPIA_RUN_TEST(petsc_tensor_reduction_test);
        UTOPIA_RUN_TEST(petsc_precond_test);
        
        //serial tests
        UTOPIA_RUN_TEST(petsc_inverse_test);

        // petsc_sparse_matrix_accessors_test();  // TODO:: here something doesnt work in parallel !
        // petsc_view_test();                   // TODO:: assert fails in parallel MatGetSubMatrix does not work if the distro changes
        // // petsc_local_entities_test(); //FIXME does not work
        //  petsc_conversion_test();
        //  maria_test();
        //  //local_diag_block_test();              // TODO:: assert fails in parallel
        
        UTOPIA_UNIT_TEST_END("PETScTest");
        #endif // WITH_PETSC
    }
}
