#include <utopia.hpp>


void assemble_laplacian_1D(const utopia::SizeType n, utopia::DSMatrixd &A)
{
    using namespace utopia;

    // n x n matrix with maximum 3 entries x row
    A = sparse(n, n, 3);
        
    {
        Write<DSMatrixd> w(A);
        Range r = row_range(A);

        // You can use add instead of set. [Warning] Petsc does not allow to mix add and set.
        for(SizeType i = r.begin(); i != r.end(); ++i) 
        {
            if(i > 0) 
            {    
                A.set(i, i - 1, -1.0);    
            }

            if(i < n-1) 
            {
                A.set(i, i + 1, -1.0);
            }

            A.set(i, i, 2.0);
        }
    }

    // final 1D FEM Laplacian
    const double h = 1.0/n; 
    A = h * A; 
}

//use "make run_assemble" to compile and run this example 
int main(int argc, char** argv)
{
    using namespace utopia;

    Utopia::Init(argc, argv);
    
    // size of problem
    const SizeType n   = 10;
    const double solution = 1.0;

    { //Only in the main file: put this scope so that the petsc objects will be destroyed before the call to finalize

        DSMatrixd A;
        assemble_laplacian_1D(n, A);

        // exact solution to our problem
        const DVectord u_exact  = values(n, solution);

        // constructing initial guess
        DVectord u = zeros(n); 

        // constructing rhs
        const DVectord rhs   = A * u_exact;

        // setting up parameters of solver 
        Parameters params; 
        params.tol(1e-9); 
        params.lin_solver_type("UTOPIA_CG"); 
        params.linear_solver_verbose(true); 
        
        // solve 
        solve(A, rhs, u, params); 

        // display solution 
        // disp(u); 

        // comparing obtained solution with 
        std::cout << "Correct solution: " <<  (approxeq(u, u_exact)? "true." : "false." ) << std::endl;

    }

    return Utopia::Finalize();
}
