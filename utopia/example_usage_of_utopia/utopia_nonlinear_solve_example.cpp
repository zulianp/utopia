#include <utopia.hpp>

template<class Matrix, class Vector>
class Rosenbrock2DFunction : public utopia::Function<Matrix, Vector> 
{
    public:
        typedef UTOPIA_SCALAR(Matrix) Scalar;

        Rosenbrock2DFunction() 
        {
            assert(!utopia::is_parallel<Matrix>::value || utopia::mpi_world_size() == 1 && "does not work for parallel matrices");
        }

        bool value(const Vector &point, Scalar &result) const override 
        {
            using namespace utopia; 
            assert(point.size().get(0) == 2);

            const Read<Vector> read(point);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            result = 100.0 * pow(y - x * x, 2.0) + pow(1.0 - x, 2.0);
            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override 
        {
            using namespace utopia; 
            assert(point.size().get(0) == 2);
            result = zeros(2);

            const Read<Vector> read(point);
            const Write<Vector> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);

            result.set(0, (400.0 * x * x * x - 400 * x * y + 2.0 * x - 2.0));
            result.set(1, 200.0 * (y - x * x));
            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override 
        {
            using namespace utopia; 
            assert(point.size().get(0) == 2);
            result = zeros(2, 2);

            const Read<Vector> read(point);
            const Write<Matrix> write(result);

            const Scalar x = point.get(0);
            const Scalar y = point.get(1);
            const Scalar mixed = -400.0 * x;

            result.set(0, 0, 1200 * x * x - 400 * y + 2);
            result.set(0, 1, mixed);
            result.set(1, 0, mixed);
            result.set(1, 1, 200.0);
            return true;
        }
    };

int main(int argc, char** argv)
{
    using namespace utopia;
    Utopia::Init(argc, argv);
    
    { //Only in the main file: put this scope so that the petsc objects will be destroyed before the call to finalize
    
        // instatiating Rosenbrock 2D banana function
        Rosenbrock2DFunction<utopia::DMatrixd, utopia::DVectord> rosenbrock_fun;

        // exact solution to our problem
        DVectord rosenbrock_exact = values(2, 1);

        // constructing initial guess
        DVectord x = values(2, 2); 

        // setting up parameters of solver 
        // Parameters params; 
        // params.tol(1e-9); 
        // params.solver_type("TRUST_REGION"); 
        // params.lin_solver_type(BICGSTAB_TAG); 
        // params.trust_region_alg(DOGLEG_TAG);
        // params.verbose(true);  

        // nonlinear solve 
        solve(rosenbrock_fun, x); 

        // comparing obtained solution with exact one 
        std::cout << "Correct solution: " <<  (approxeq(x, rosenbrock_exact)? "true." : "false." ) << std::endl;
    }

    return Utopia::Finalize();
}
