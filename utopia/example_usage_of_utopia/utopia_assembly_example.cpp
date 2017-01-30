#include <utopia.hpp>



//only serial
template<class Matrix, class Vector>
void nl_gauss_seidel_step(const Matrix &mat, const Vector &diag, const Vector &rhs, const Vector &upper_bound, Vector &solution)
{
    using namespace utopia;
    typedef double Scalar;

    Read<Vector> r_rhs(rhs);
    Read<Vector> r_diag(diag);
    
    Range rr = row_range(mat);
    SizeType offset = rr.begin();

    std::vector<double> sol(rr.extent());

    each_read(solution, [&](const SizeType i, const Scalar value) {
        assert(i - offset >= 0);
        sol[i - offset] = value;
    });

    Scalar sum = 0.;
    Scalar sol_value = 0.0;
    Scalar upper_bound_i = 0.0;
    SizeType prev_I = rr.begin();

    each_read(mat, [&](const SizeType i, const SizeType j, const Scalar value) {
        // std::cout << i << ", " << j << " -> " << value << std::endl;
    
        if(i != prev_I) {
            sol_value = (rhs.get(prev_I) - sum)/diag.get(prev_I);

            upper_bound_i = upper_bound.get(prev_I);
            
            if(sol_value > upper_bound_i) {
                sol_value = upper_bound_i;
            }

            sol[prev_I - offset] = sol_value;

            prev_I = i;
            sum = 0.;
        }

        if(i != j) {
            sum += value * sol[j - offset] ;
        }

    });

    sol_value =  (rhs.get(prev_I) - sum)/diag.get(prev_I);;
    upper_bound_i = upper_bound.get(prev_I);
    
    if(sol_value > upper_bound_i) {
        sol_value = upper_bound_i;
    }

    sol[prev_I - rr.begin()] = sol_value;
    solution.set(prev_I, sol_value);

    // std::cout << "HERE:\n";
    each_write(solution, [&](const SizeType i) -> Scalar {
        assert(i - offset >= 0);
        // std::cout << sol[i - offset] << "\t";
        return sol[i - offset];
    });

    // std::cout << "\nTHERE:\n";
}

//only serial
template<class Matrix, class Vector>
bool nl_gauss_seidel_solve(const Matrix &mat, const Vector &rhs, const Vector &upper_bound, Vector &solution, const int max_iterations = 1000)
{   
    using namespace utopia;
    int check_residual_each = 1;

    Vector diag          = utopia::diag(mat);
    Vector prev_solution = zeros(size(solution));

    disp(diag);

    for(int k = 0; k < max_iterations; ++k) {
        nl_gauss_seidel_step(mat, diag, rhs, upper_bound, solution);
    
        if((k + 1) % check_residual_each == 0) {
            double diff  = norm2(prev_solution - solution);
            double res   = norm2(mat * solution - rhs);

            std::cout << "iter: " << k << " change in the solution: " << diff <<  " residual: " << res << std::endl;
            if(diff < 1e-14) {
                return true;
            }

            prev_solution = solution;
        }  

        if(k == 0) {
            prev_solution = solution;
        }
    }

    return false;
}

template<class Matrix>
void assemble_laplacian_1D(const utopia::SizeType n, Matrix &m)
{
    using namespace utopia;

    // n x n matrix with maximum 3 entries x row        
    {
        Write<Matrix> w(m);
        Range r = row_range(m);

        //You can use set instead of add. [Warning] Petsc does not allow to mix add and set.
        for(SizeType i = r.begin(); i != r.end(); ++i) {
            if(i > 0) {    
                m.add(i, i - 1, -1.0);    
            }

            if(i < n-1) {
                m.add(i, i + 1, -1.0);
            }

            m.add(i, i, 2.0);
        }
    }
}


template<class Matrix, class Vector>
void solve_constrained_poisson_problem()
{
    using namespace utopia;

    const SizeType n   = 10;

    Matrix m = zeros(n, n);
    assemble_laplacian_1D(n, m);
    {
        Range r = row_range(m);
        Write<Matrix> w(m);
        if(r.begin() == 0) {
            m.set(0, 0, 1.);
            m.set(0, 1, 0);
        }

        if(r.end() == n) {
            m.set(n-1, n-1, 1.);
            m.set(n-1, n-2, 0);
        }
    }

    Vector rhs = values(n, 1.);
    { 
        //Creating test vector (alternative way see [assemble vector alternative], which might be easier for beginners)
        Range r = range(rhs);
        Write<Vector> w(rhs);

        if(r.begin() == 0) {
            rhs.set(0, 0);
        }

        if(r.end() == n) {
            rhs.set(n-1, 0.);
        }
    }

    Vector upper_bound = values(n, 9.0);
    Vector solution    = zeros(n);
    nl_gauss_seidel_solve(m, rhs, upper_bound, solution);
    disp(solution);
    disp(m);
    disp(rhs);

}

//use "make run_assemble" to compile and run this example 
int main(int argc, char** argv)
{
    using namespace utopia;

    Utopia::Init(argc, argv);

    const SizeType n   = 3;
    const double value = 6.0;

    { //Only in the main file: put this scope so that the petsc objects will be destroyed before the call to finalize

        DSMatrixd m = sparse(n, n, 3);
        assemble_laplacian_1D(n, m);

        //Made const since it is immutable
        const DVectord v  = values(n, value);
        DVectord actual   = m * v;

        DVectord expected = zeros(n);
        
        { 
            //Creating test vector (alternative way see [assemble vector alternative], which might be easier for beginners)
            Range r = range(expected);
            Write<DVectord> w(expected);

            if(r.begin() == 0) {
                expected.set(0, value);
            }

            if(r.end() == n) {
                expected.set(n-1, value);
            }
        }

        //expecting a vector of zeros except the 1st and last entries which will be "value"
        std::cout << "correct: " <<  (approxeq(actual, expected)? "true" : "false" ) << std::endl;

        //to display the ouput in the terminal (alternative way see [print vector alternative])
        disp(actual);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////  ALTERNATIVES  ////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //[assemble vector alternative]:
        DVectord rhs = zeros(n);
        each_write(rhs, [value](const SizeType i) -> double {
            //The returned value will be written in the vector
            if(i == 0 || i == n-1) {
                return 0;
            }
            return value;
        });

        //Another way to print it is to iterate over the vector and read 
        //[print vector alternative]:
        each_read(rhs, [](const SizeType i, const double entry) {
            std::cout << "rhs(" << i << ") = " << entry << std::endl;
        });

        {
            SizeType local_n = 3;
            SizeType block_size = 2; 
            DVectord x       = local_zeros (local_n * block_size);
            DVectord f       = local_values(local_n * block_size, 0.5);
            DSMatrixd mat    = local_identity(local_n, local_n);

            // kron_prod(mat x identity) * f
            {
                Read<DVectord>  r(x);
                Write<DVectord> w(f);

                each_read(mat, [block_size, &x, &f](const SizeType i, const SizeType j, const double entry) 
                {
                    std::cout << "m(" << i << ", " << j << ") = " << entry << std::endl;

                    for(SizeType k = 0; k < block_size; ++k) 
                    {
                        f.add(i * block_size + k, x.get(j * block_size + k) * entry);
                    }
                });
            }

            disp(f);
        }
    }

    solve_constrained_poisson_problem<Matrixd, Vectord>();
    return Utopia::Finalize();
}
