#include <utopia.hpp>

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

    const SizeType n = 40;

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

    Vector upper_bound = values(n, 100.0);
    Vector solution    = zeros(n);

    ProjectedGaussSeidel<Matrix, Vector> pgs;
    pgs.set_box_constraints(make_upper_bound_constraints(make_ref(upper_bound)));
    pgs.solve(m, rhs, solution);

    write("U.m", solution);
    write("M.m", m);
    write("R.m", rhs);
}

//use "make run_assemble" to compile and run this example 
int main(int argc, char** argv)
{
    using namespace utopia;

    Utopia::Init(argc, argv);

    const SizeType n   = 3;
    const double value = 6.0;

    { //Only in the main file: put this scope so that the petsc objects will be destroyed before the call to finalize

        PetscMatrix m = sparse(n, n, 3);
        assemble_laplacian_1D(n, m);

        //Made const since it is immutable
        const PetscVector v  = values(n, value);
        PetscVector actual   = m * v;

        PetscVector expected = zeros(n);
        
        { 
            //Creating test vector (alternative way see [assemble vector alternative], which might be easier for beginners)
            Range r = range(expected);
            Write<PetscVector> w(expected);

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
        PetscVector rhs = zeros(n);
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
            PetscVector x       = local_zeros (local_n * block_size);
            PetscVector f       = local_values(local_n * block_size, 0.5);
            PetscMatrix mat    = local_identity(local_n, local_n);

            // kron_prod(mat x identity) * f
            {
                Read<PetscVector>  r(x);
                Write<PetscVector> w(f);

                each_read(mat, [block_size, &x, &f](const SizeType i, const SizeType j, const double entry) 
                {
                    std::cout << "m(" << i << ", " << j << ") = " << entry << std::endl;

                    for(SizeType k = 0; k < block_size; ++k) 
                    {
                        f.add(i * block_size + k, x.get(j * block_size + k) * entry);
                    }
                });
            }

            // disp(f);
        }
    }

    solve_constrained_poisson_problem<PetscMatrix, PetscVector>();
    return Utopia::Finalize();
}
