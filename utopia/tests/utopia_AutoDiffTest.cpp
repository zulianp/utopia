#include "utopia_AutoDiffTest.hpp"
#include "utopia.hpp"
#include "test_problems/utopia_TestProblems.hpp"

class ExampleDiffFun {
public:
    template<class X>
    inline auto operator()(const X &x) const -> decltype( dot(x, x) ) {
        return dot(x, x);
    }
};

template<class Matrix, class Vector>
void diff_test()
{
    using namespace utopia;

    const int n = 10;
    //some example vectors
   Vector x = values(n, 25.0);
   Vector b = values(n, 1.0);

   //some example matrices 
   Matrix A = zeros(n, n);
   Matrix B = 0.5 * identity(n, n);
   Matrix C = identity(n, n);

   {
       Write<Matrix> w(A);
       Range r = row_range(A);
       for(SizeType i = r.begin(); i != r.end(); ++i) {
           if(i > 0) {    
               A.add(i, i - 1, -1.0);    
           }

           if(i < n-1) {
               A.add(i, i + 1, -1.0);
           }

           A.add(i, i, 2.0);
       }
   }

    //Create the independent variable vector (VERY IMPORTANT)
    auto d_x = independent_variable(x);

    //Extended example
    {
        //create valued expression 
        auto expr = 0.5 * dot(B * d_x, A * pow2(d_x)) + dot(C * d_x,  A * b);

        //create derivative 
        auto d_expr = derivative(expr);

        //evaluate derivative
        Vector df = d_expr;
    }

    //Short version
    {
        Vector df = derivative( 0.5 * dot(B * d_x, A * pow2(d_x)) + dot(C * d_x,  A * b) );
        // disp(df);   
    }

    //2nd order derivative
    {   
        auto expr     =  0.5 * dot(B * d_x, A * pow2(d_x)) + dot(C * d_x,  A * b) ;
        auto d_expr   = derivative(expr);
        auto d_expr_2 = derivative(d_expr);


        // std::cout << tree_format(d_expr_2.getClass()) << std::endl;
     
        Number<double> f = expr;    //or double f = scalar_cast<double>(expr);
        Vector g = d_expr;
        Matrix H = d_expr_2;
        // disp(H);
    }

    //Using automatic diff in the context of newton solvers
    {
        //To be worked on
        auto f = auto_diff_fun<Matrix, Vector>(ExampleDiffFun());
        Vector sol = values(n, 1.0);
    
        Newton<Matrix, Vector> newton(std::make_shared<ConjugateGradient<Matrix, Vector>>(A));
        newton.solve(f, sol);
        // disp(sol);
    }
}

template<class Matrix, class Vector>
void sparse_diff_test()
{
    using namespace utopia;

    const int n = 100;
    //some example vectors
    Vector x = values(n, 25.0);
    Vector b = values(n, 1.0);

    //some example matrices 
    Matrix A = sparse(n, n, 3);
    Matrix B = 0.5 * identity(n, n);
    Matrix C = identity(n, n);

    {
        Write<Matrix> w(A);
        Range r = row_range(A);
        for(SizeType i = r.begin(); i != r.end(); ++i) {
            if(i > 0) {    
                A.add(i, i - 1, -1.0);    
            }

            if(i < n-1) {
                A.add(i, i + 1, -1.0);
            }

            A.add(i, i, 2.0);
        }
    }

    //Create the independent variable vector (VERY IMPORTANT)
    auto d_x = independent_variable(x);

    //Extended example
    {
        //create valued expression 
        auto expr = 0.5 * dot(B * d_x, A * pow2(d_x)) + dot(C * d_x,  A * b);

        //create derivative 
        auto d_expr = derivative(expr);

        Number<double> f = expr; //or double f = scalar_cast<double>(expr);
        
        //evaluate derivative
        Vector df = d_expr;
    }

    //Short version
    {
        Vector df = derivative( 0.5 * dot(B * d_x, A * pow2(d_x)) + dot(C * d_x,  A * b) );
        // disp(df);   
    }

    //2nd order derivative
    {   
        auto expr     = 0.5 * dot(B * d_x, A * pow2(d_x)) + dot(C * d_x,  A * b) ;
        auto d_expr   = derivative(expr);
        auto d_expr_2 = derivative(d_expr);

        // std::cout << tree_format(d_expr_2.getClass()) << std::endl;

        Number<double> f = expr;
        Vector g = d_expr;
        Matrix H = d_expr_2;
    }
}

template<class Matrix, class Vector>
void trace_test()
{
    using namespace utopia;

    int n = 10;
    Matrix x = identity(n, n);

    // double r = trace(x);
    //disp(r);

    auto d_x  = independent_variable(x);
    auto expr = trace(d_x);
    auto d_expr = derivative(expr);

   // std::cout << tree_format(d_expr.getClass()) << std::endl;
    Matrix f = d_expr;
    //disp(f);

}

void run_autodiff_test()
{
	using namespace utopia;

	std::cout << "Begin: AutoDiffTest" << std::endl;
#ifdef WITH_BLAS
	// diff_test<Matrixd, Vectord>();

	//FIXME does not work -> incomplete backend
	// sparse_diff_test<CRSMatrixd, Vectord>();

#endif //WITH_BLAS

#ifdef WITH_PETSC
	if(mpi_world_size() == 1) {
		//petsc does not support parallel dense operations
	    // diff_test<DMatrixd, DVectord>();
	}

    // sparse_diff_test<DSMatrixd, DVectord>();
    trace_test<DSMatrixd, DVectord>();
#endif //WITH_PETSC	

    std::cout << "End: AutoDiffTest" << std::endl;
}
