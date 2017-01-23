#include <utopia.hpp>

template<class MatrixT, class VectorT>
void run()
{
    using namespace utopia;

    const SizeType n = 3;

    MatrixT m = sparse(n, n, 2);
    {
        Write<MatrixT> w(m);
        auto r = row_range(m);

        for(auto i = r.begin(); i != r.end(); ++i) {
            m.set(i, i, 0.5);

            if(i > 0) {
                m.set(i, i-1, 1.0);
            }

            if(i < n-1) {
                m.set(i, i+1, 1.0);
            }
        }
    }

    { 
        VectorT a = zeros(n);
        VectorT b = zeros(n);
        VectorT c = m * values(n, 0.5);  //multipliting a vector of 0.5s

        //writing i/n in the vector
        each_write(a, [](const SizeType i) -> double  { return i/double(n); }   );

        {
            //if another vector is needed, just provide a lock and pass it to the lambda functor
            Read<VectorT> r(c);

            //applying a filter to the content of a and writing it into b. a cannot be equal to b (for the moment)
            each_transform(a, b, [&c](const SizeType i, const double entry) -> double  { return exp(-entry) * c.get(i); }    );
        }

        // disp(a);
        // disp(b);    
    }
}

//use "make run_rwt" to compile and run this example. This example shows how to create simple implementation 
//independent code and how to use the macros to instantiate the specific implementations. 
int main(int argc, char** argv)
{
    using namespace utopia;

    Utopia::Init(argc, argv);

//if it has compiled with blas or petsc WITH_BLAS or WITH_PETSC macros are available (if you want to make it compile no matter the utopia installation)
#ifdef WITH_PETSC    
    //run with petsc types 
    run<DSMatrixd, DVectord>();
#endif //WITH_PETSC 

#ifdef WITH_BLAS    
    //run with blas types 
    // run<CRSMatrixd, Vectord>();
#endif //WITH_BLAS    

    return Utopia::Finalize();
}
