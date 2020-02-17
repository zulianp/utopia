#include <utopia.hpp>


void assemble_laplacian_1D(const utopia::SizeType n, utopia::PetscMatrix &A)
{
    using namespace utopia;

    // n x n matrix with maximum 3 entries x row
    A = sparse(n, n, 3);
        
    {
        Write<PetscMatrix> w(A);
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

        PetscMatrix A;
        assemble_laplacian_1D(n, A);

        // exact solution to our problem
        const PetscVector u_exact  = values(n, solution);

        // constructing initial guess
        PetscVector u = zeros(n); 

        // constructing rhs
        const PetscVector rhs   = A * u_exact;
        
        // solve 
        solve(A, rhs, u); 

        // display solution 
        // disp(u); 

        // comparing obtained solution with 
        std::cout << "Correct solution: " <<  (approxeq(u, u_exact)? "true." : "false." ) << std::endl;

    }

    return Utopia::Finalize();
}
