
#include <iostream>
#include "utopia_petsc.hpp"
#include "utopia_Instance.hpp"

extern "C" {

    void utopia_initialize(int argc, char *argv[])
    {
        using namespace utopia;
        std::cout << "utopia_initialize" << std::endl;
        Utopia::Init(argc, argv);
    }

    void utopia_finalize()
    {
        std::cout << "utopia_finalize" << std::endl;

        using namespace utopia;
        Utopia::Finalize();
    }

    void utopia_create_solver(const char name[], void **ptr)
    {
        using namespace utopia;
        std::cout << "utopia_create_solver: " << name << std::endl;

        auto solver_ptr = reinterpret_cast<LinearSolver<PetscMatrix, PetscVector>**>(ptr);

        static const std::string ksp = "ksp";
        static const std::string fact = "fact";

        if(fact == name) {
            *solver_ptr = new Factorization<PetscMatrix, PetscVector>();
        } else {
            *solver_ptr = new KSPSolver<PetscMatrix, PetscVector>();
        }
    }

    void utopia_print_solver_info(void *ptr)
    {
        using namespace utopia;

        auto solver_ptr = reinterpret_cast<LinearSolver<PetscMatrix, PetscVector>*>(ptr);
        auto ksp = dynamic_cast<KSPSolver<PetscMatrix, PetscVector>*>(solver_ptr);

        if(ksp) {
            std::cout << "KSP" << std::endl;
            ksp->describe(std::cout);
        } else {
            auto fact = dynamic_cast<Factorization<PetscMatrix, PetscVector>*>(solver_ptr);
            fact->describe(std::cout);
        }

    }

    void utopia_destroy_solver(void **ptr)
    {
        using namespace utopia;
        std::cout << "utopia_destroy_solver" << std::endl;

        auto solver_ptr = reinterpret_cast<LinearSolver<PetscMatrix, PetscVector>**>(ptr);

        delete *solver_ptr;
        *solver_ptr = nullptr;
    }

    void utopia_solve(void *solver, void*A, void*b, void *x)
    {
        std::cout << "(void" << std::endl;
    }

}
