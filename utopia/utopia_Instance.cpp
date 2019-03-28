#include "utopia_Instance.hpp"
#include "utopia_Base.hpp"
#include "utopia_Tracer.hpp"
#include "utopia_Config.hpp"

#ifdef WITH_TRILINOS
#include <Tpetra_Core.hpp>
#endif //WITH_TRILINOS

#ifdef WITH_PETSC
#include "petscsys.h"
#include "utopia_petsc_build_ksp.hpp"
#else
#ifdef WITH_MPI
#include <mpi.h>
#endif //WITH_MPI
#endif //WITH_PETSC



#include <cassert>

namespace utopia {

    // #define DISABLE_LOGGER

    void Utopia::Init(int argc, char *argv[]) {

#ifdef WITH_PETSC
        static char help[] = "initializing utopia environment through petsc";

    #ifdef WITH_SLEPC
        SlepcInitialize(&argc,&argv,(char*)0, help); // calls PetscInitialize inside
    #else
        PetscInitialize(&argc, &argv, (char *) 0, help);
    #endif    //WITH_SLEPC

        // is this proper place for doing this ???
        KSPRegister("utopia", KSPCreate_UTOPIA);

        //        PetscInitializeNoArguments();
#else
#ifdef WITH_MPI
        MPI_Init(&argc, &argv);
#endif //WITH_MPI
#endif //WITH_PETSC

#ifdef WITH_TRILINOS
        Tpetra::initialize(&argc, &argv);
#endif //WITH_TRILINOS
    }


    int Utopia::Finalize() {
#ifdef UTOPIA_TRACE_ENABLED
        Tracer::instance().save_collected_log();
#endif

        instance().logger().flush();

        if(mpi_world_rank() == 0) {
            instance().maintenance_logger().flush();
        }

#ifdef WITH_PETSC
        #ifdef WITH_SLEPC
            SlepcFinalize(); // calls PetscFinalize inside
        #else
            PetscFinalize();
        #endif    //WITH_SLEPC

#else
#ifdef WITH_MPI
        return MPI_Finalize();
#endif //WITH_MPI
#endif //WITH_PETSC

#ifdef WITH_TRILINOS
        Tpetra::finalize();
#endif //WITH_TRILINOS

        return  instance().exit_code_;
    }

    Utopia &Utopia::instance()
    {
        static Utopia instance;
        return instance;
    }

    Utopia::Utopia()
    : exit_code_(EXIT_SUCCESS)
    {
        set("data_path", "../data");
        set("opencl_templates_path", "../backend/opencl/templates");

#ifdef DISABLE_LOGGER
        logger_ = std::make_shared<NullLogger>();
        auto temp = std::make_shared<StandardLogger>();
        maintenance_logger_ = temp;
#else
        logger_ = std::make_shared<StandardLogger>();
        auto temp = std::make_shared<StandardLogger>();
        temp->set_direct_output(false, false, false);
        maintenance_logger_ = temp;
#endif

    }

    bool Utopia::verbose() const
    {
        return get("verbose") == "true";
    }

    void Utopia::Abort()
    {
        int error_code = -1;
#ifdef WITH_MPI
        MPI_Abort(MPI_COMM_WORLD, error_code);
#else
        exit(error_code);
#endif
    }
}
