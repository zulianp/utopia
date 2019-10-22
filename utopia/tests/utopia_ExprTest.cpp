#include "utopia.hpp"
#include "utopia_Testing.hpp"
#include "test_problems/utopia_TestProblems.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_For.hpp"
#include "utopia_Eval_Residual.hpp"
#include <utility>

namespace utopia {

    template<typename Matrix, typename Vector>
    class ExpressionTests
    {

    public:
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef UTOPIA_SCALAR(Vector) Scalar;

        ExpressionTests(const SizeType & n): n_(n)
        {}

        void run()
        {
            UTOPIA_RUN_TEST(axpy_test);
            UTOPIA_RUN_TEST(e_div_test);
            UTOPIA_RUN_TEST(e_mul_test);
            UTOPIA_RUN_TEST(for_each_loop_test);
            UTOPIA_RUN_TEST(multi_reduce_test);
            UTOPIA_RUN_TEST(mv_test);
            UTOPIA_RUN_TEST(negate_alpha_test);
            UTOPIA_RUN_TEST(negate_test);
            UTOPIA_RUN_TEST(parallel_each_write_test);
            UTOPIA_RUN_TEST(quad_form_test);
            UTOPIA_RUN_TEST(residual_test);
            UTOPIA_RUN_TEST(transform_test);
            UTOPIA_RUN_TEST(mat_copy); 
            UTOPIA_RUN_TEST(reciprocal_test); 
            UTOPIA_RUN_TEST(max_min_test); 
            UTOPIA_RUN_TEST(multi_axpy);
            UTOPIA_RUN_TEST(inv_diag);
            UTOPIA_RUN_TEST(comp_mat);
            UTOPIA_RUN_TEST(bratu_grad);
            
            // FIXME (mem allocs)
            

            // Seq. fault 
            UTOPIA_RUN_TEST(mat_transp_mult_test); 
        }

        void mat_transp_mult_test()
        {
            Matrix H = sparse(n_, n_, 3); 
            assemble_laplacian_1D(H);
            Matrix D = diag(diag(H));
            H = H + transpose(H) - D;
        }


        void negate_test()
        {
            Vector x = values(n_, 1.0);
            Vector y = values(n_, 10.0);

            UTOPIA_NO_ALLOC_BEGIN("negate_test");
            x = -x;
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("copy_test");
            x = y;
            UTOPIA_NO_ALLOC_END();            
        }

        void mv_test()
        {
            Vector x = values(n_, 1.0);
            Vector b = values(n_, 1.0);
            Vector p = values(n_, 1.0);

            Matrix A = sparse(n_, n_, 3);
            assemble_laplacian_1D(A);

            UTOPIA_NO_ALLOC_BEGIN("mv_test1");
            p = A * x + b;
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("mv_test2");
            p = A * x - b;
            UTOPIA_NO_ALLOC_END();     

            UTOPIA_NO_ALLOC_BEGIN("mv_test3");
            p = b - A * x;
            UTOPIA_NO_ALLOC_END();   

            UTOPIA_NO_ALLOC_BEGIN("mv_test4");
            p = b + A * x;
            UTOPIA_NO_ALLOC_END();                 

            UTOPIA_NO_ALLOC_BEGIN("mv_test5");
            p = 0.5 * b + A * x;
            UTOPIA_NO_ALLOC_END();         

            UTOPIA_NO_ALLOC_BEGIN("mv_test6");
            p = 0.5 * b - A * x;
            UTOPIA_NO_ALLOC_END();   

            //This will always create a copy (A * x needs to put the result somewhere, unless the backend specialization is used and does not allocate anything)
            //a better way would be to use a third vector p = b - A * x;
            b = b - A * x;

            UTOPIA_NO_ALLOC_BEGIN("mv_test7");
            p = b - A * x;
            UTOPIA_NO_ALLOC_END();      
        }

        void mat_copy()
        {

            Matrix I = identity(n_, n_); 
            Matrix D = identity(n_, n_); 
            Vector v = diag(I); 

            // I do not know how relevant are these tests, as sparsity pattern might be different... 
            //create same_sparsity_copy()
            // UTOPIA_NO_ALLOC_BEGIN("mat_copy1");
            //FIME
            // D = I; 
            // UTOPIA_NO_ALLOC_END();     

            // UTOPIA_NO_ALLOC_BEGIN("mat_copy2");
            //FIXME this it is equivalent to making a copy (but it could be treated as an AXPY)
            D = -1.0 * I;
            D = -I;
            // UTOPIA_NO_ALLOC_END();     

            UTOPIA_NO_ALLOC_BEGIN("mat_copy3");
            D -= I; 
            UTOPIA_NO_ALLOC_END();     

            UTOPIA_NO_ALLOC_BEGIN("mat_copy4");
            D += I; 
            UTOPIA_NO_ALLOC_END();        

            UTOPIA_NO_ALLOC_BEGIN("mat_copy5");
            D += 5.0*D; 
            UTOPIA_NO_ALLOC_END();                             

            D += Matrix(diag(v)); 
            
            UTOPIA_NO_ALLOC_BEGIN("mat_copy5");//
            //https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatDiagonalSet.html#MatDiagonalSet
            D += diag(v); //bad way D += Matrix(diag(v)); 
            UTOPIA_NO_ALLOC_END();  

            // UTOPIA_NO_ALLOC_BEGIN("mat_copy6");
            //FIME still creates a temporary (but now it is just a vector)
            D -= diag(v); //bad way D -= Matrix(diag(v)); 
            // UTOPIA_NO_ALLOC_END();  
        }


        void e_mul_test()
        {
            Vector x = values(n_, 1.0);
            Vector y = values(n_, 2.0);
            Vector z = values(n_, 0.0);

            UTOPIA_NO_ALLOC_BEGIN("e_mul_test");
            z = e_mul(x, y);
            UTOPIA_NO_ALLOC_END();
            utopia_test_assert(approxeq(sum(z), 2*n_));
        }


        void max_min_test()
        {
            Vector x = values(n_, 1.0);
            Vector y = values(n_, 2.0);
            Vector z = values(n_, 0.0);

            UTOPIA_NO_ALLOC_BEGIN("max_min_test1");
            z = utopia::max(x,y); 
            UTOPIA_NO_ALLOC_END();
            utopia_test_assert(approxeq(sum(z), 2*n_));

            UTOPIA_NO_ALLOC_BEGIN("max_min_test2");
            z = utopia::min(x,y); 
            UTOPIA_NO_ALLOC_END();
            utopia_test_assert(approxeq(sum(z), n_));

            UTOPIA_NO_ALLOC_BEGIN("max_min_test3");
            z = utopia::max(utopia::min(x,y), y); 
            UTOPIA_NO_ALLOC_END();

            utopia_test_assert(approxeq(sum(z), 2*n_));

        }


        void reciprocal_test()
        {
            Vector x = values(n_, 2.0);

            UTOPIA_NO_ALLOC_BEGIN("reciprocal_test");
            x = 1./x; 
            UTOPIA_NO_ALLOC_END();
            utopia_test_assert(approxeq(sum(x), 0.5*n_));
        }        

        void e_div_test()
        {
            Vector x = values(n_, 6.0);
            Vector z = values(n_, 3.0);

            UTOPIA_NO_ALLOC_BEGIN("e_div_test");
            z = x / z;
            
            Scalar sum_z = sum(z);
            utopia_test_assert(approxeq(sum_z, 2.0*n_));
            z = x / x;
            sum_z = sum(z);
            utopia_test_assert(approxeq(sum_z, 1.0*n_));

            z = z / x;
            sum_z = sum(z);
            utopia_test_assert(approxeq(sum_z, 1.0/6.0*n_));

            UTOPIA_NO_ALLOC_END();
        }

        void transform_test()
        {
            Vector x = values(n_, 1.0);
            
            parallel_transform(
                x,
                UTOPIA_LAMBDA(const SizeType &i, const Scalar &v) -> Scalar {
                    return (i+1)*v;
            });
            
            Scalar expected = ((n_ + 1) * n_)/2.0;
            Scalar sum_x = sum(x);

            utopia_test_assert(approxeq(sum_x, expected, 1e-10));
        }

        void negate_alpha_test()
        {
            Vector x = values(n_, 1.0);
            Vector y = values(n_, 1.0);
            
            UTOPIA_NO_ALLOC_BEGIN("negate_alpha_test");
            y = -0.5 * x;
            UTOPIA_NO_ALLOC_END();
        }

        void quad_form_test()
        {
            Vector x = values(n_, 1.0);
            Vector y = values(n_, 2.0);

            auto expr = 0.5 * dot(x, y) - 0.5 * dot(x, y);

            Scalar val = expr;

            utopia_test_assert(approxeq(val, 0.0));
        }

        void residual_test()
        {
            Vector x = values(n_, 1.0);
            Vector b = values(n_, 1.0);
            Matrix A = sparse(n_, n_, 3);
            Vector r = values(n_, 0.0);

            UTOPIA_NO_ALLOC_BEGIN("residual_test");
            r = x - b;
            UTOPIA_NO_ALLOC_END();
        }

        void axpy_test()
        {
            Vector x = values(n_, 1.0);
            Vector y = values(n_, 1.0);
            Vector p = values(n_, 2.0);

            UTOPIA_NO_ALLOC_BEGIN("axpy_test0");
            x = x - 0.5 * p;
            y = x - 0.5 * p;
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("axpy_test1");
            x.set(1.0); 
            x = x + 0.5 * p;
            UTOPIA_NO_ALLOC_END();      
            Scalar val = sum(x);      
            utopia_test_assert(approxeq(val, n_*2.0));

            UTOPIA_NO_ALLOC_BEGIN("axpy_test2");
            x.set(2.0); 
            y.set(1.0); 
            x = y + 0.5*x; 
            UTOPIA_NO_ALLOC_END();           
            val = sum(x);          
            utopia_test_assert(approxeq(val, n_*2.0));
        }

        void for_each_loop_test()
        {
            Vector x = values(n_, 2.);
            Vector y = values(n_, 1.);
            Vector z = values(n_, 0.);

            using ForLoop = utopia::ParallelFor<Traits<Vector>::Backend>;

           {
                auto d_x = const_device_view(x);
                auto d_y = const_device_view(y);
                auto d_z = device_view(z);

                ForLoop::apply(range(z), UTOPIA_LAMBDA(const SizeType i)
                {
                    const Scalar xi = d_x.get(i);
                    const Scalar yi = d_y.get(i);
                    d_z.set(i, xi - yi);
                });
            }

            Scalar sum_z = sum(z);
            utopia_test_assert(approxeq(Scalar(n_), sum_z));
        }

        void parallel_each_write_test()
        {
            Vector x = values(n_, 2.);
            Vector y = values(n_, 1.);
            Vector z = values(n_, 0.);

            {
                auto d_x = const_device_view(x);
                auto d_y = const_device_view(y);
                auto d_z = device_view(z);

                parallel_each_write(z, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                {
                    const Scalar xi = d_x.get(i);
                    const Scalar yi = d_y.get(i);
                    return xi - yi;
                });
            }

            Scalar sum_z = sum(z);
            utopia_test_assert(approxeq(Scalar(n_), sum_z));
        }
        
        void multi_reduce_test()
        {
            Vector x = values(n_, 2.);
            Vector y = values(n_, 1.);

            each_write(x, [](const SizeType &i) -> Scalar {
                return -(i + 1.0);
            });

            each_write(y, [](const SizeType &i) -> Scalar {
                return (i + 1.0);
            });

            const Scalar m = multi_min(x, y);
            utopia_test_assert(approxeq(m, Scalar(-n_)));
        }

        void multi_axpy()
        {
            //Assign<Vec, Minus<Plus<Vec, Multiplies<Number, Vec>>, Multiplies<Number, Vec>>>
            Vector a = values(n_, 2.);
            Vector b = values(n_, 2.);
            Vector c = values(n_, 2.);
            Vector result = zeros(n_);

            UTOPIA_NO_ALLOC_BEGIN("multi_axpy");
            Scalar alpha = 1.0, beta = 2.0; 
            result = a + (alpha * b) - (beta * c);
            UTOPIA_NO_ALLOC_END();
        }

        void inv_diag()
        {   
            Matrix H = sparse(n_, n_, 3);
            Vector d = zeros(n_);

            assemble_laplacian_1D(H);

            UTOPIA_NO_ALLOC_BEGIN("inv_diag");
            d = 1./diag(H);
            UTOPIA_NO_ALLOC_END();
        }

        void comp_mat()
        {
            Matrix H = identity(n_, n_);
            Matrix R = H;

            UTOPIA_NO_ALLOC_BEGIN("comp_mat");
            R = H + H * H;
            UTOPIA_NO_ALLOC_END();

            Matrix Id = 2.0 * identity(n_, n_);
            utopia_test_assert(approxeq(Id, R));
        }

        void bratu_grad()
        {
            Vector result = zeros(n_);
            Vector x = zeros(n_);

            UTOPIA_NO_ALLOC_BEGIN("bratu_grad");
            result = x - (0.5 * exp(x)); 
            UTOPIA_NO_ALLOC_END();
        }

    private:
        SizeType n_;
    };

    void expr()
    {
#ifdef WITH_PETSC
        auto n_dofs     = 10;

        ExpressionTests<PetscMatrix, PetscVector>(n_dofs).run();

#ifdef WITH_TRILINOS
        ExpressionTests<TpetraMatrixd, TpetraVectord>(n_dofs).run();
#endif
#endif

    }

    UTOPIA_REGISTER_TEST_FUNCTION(expr);
}
